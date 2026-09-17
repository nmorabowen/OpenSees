"""WP-109 -- the Python module must carry exactly ONE libstdc++ runtime.

This file PINS the root cause of the gcc + `-fopenmp` Zone-A segfault banked by
PR #843 (run 35164371356) and root-caused by WP-109 (run 35171085324):

    test_adr30_projection_p0.py::test_massless_dof_is_not_policeable_by_the_soe_layer
    Fatal Python error: Segmentation fault      -> pytest exit 139

gdb on the ubuntu-latest runner put the fault at the first `opserr << int` on the
zero-mass `system Diagonal` refusal:

    std::codecvt<char16_t,char>::do_unshift     <- wrong facet, wrong vtable
    std::ostream::_M_insert<long>
    PythonStream::err_out<int>                  SRC/interpreter/PythonStream.h
    DiagonalDirectSolver::solve                 `opserr << i`

The std::stringstream inside PythonStream asked its locale for the num_put facet
and got a codecvt facet: the process held TWO libstdc++ runtimes. The fork sets
`CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++"`; FindOpenMP's probe is
an executable try_compile, so with LADRUNO_OPENMP=ON it reported the static
`libstdc++.a` as an OpenMP "implicit library" (`OpenMP_CXX_LIB_NAMES =
libstdc++;gomp;pthread`), and CMakeLists.txt appended `${OpenMP_CXX_LIBRARIES}` to
the SHARED OpenSeesPy module -- which also DT_NEEDs `libstdc++.so.6`. The loader
then bound `std::num_put<char>::id` to one copy and `_M_insert<long>` to the other.

WHAT THIS FILE ASSERTS (Linux/ELF only; a PE `.pyd` is skipped with a reason --
MSVC has no libstdc++ and `/MT` is a different, deliberate choice):

  1. the loaded module DEFINES no `std::locale::*`, `std::ios_base::*` or
     `std::codecvt<char16_t,...>` symbol in its dynamic symbol table. Those are
     non-template classes that live ONLY in libstdc++.so; a definition inside
     opensees.so means a second libstdc++ was linked in.
  2. `libstdc++.so.6` IS in DT_NEEDED -- the one true runtime.
  3. sanity: `.dynsym` was really parsed (PyInit_opensees is defined) so an
     empty-table false pass is impossible.

MUTATION GATE: reverting the WP-109 filter in CMakeLists.txt (linking
`${OpenMP_CXX_LIBRARIES}` again) turns (1) red on every gcc build with
LADRUNO_OPENMP=ON -- measured on esmeralda (gcc 11.4, 151 libstdc++ internals
exported) BEFORE the fix and 0 after. Pure-Python ELF reader: no binutils needed,
so the test can never skip for a missing tool (BUILD_GOTCHAS section 15's
"biased to RUN" rule).

Theory: Ladruno_implementation/75b_ladruno_threaded_assembly_adr.md section 14.5;
Ladruno_internal/BUILD_GOTCHAS.md section 16.
"""
import os
import struct

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

# std::locale, std::ios_base (Itanium-mangled nested names) and the C++11
# char16_t/char32_t codecvt facets -- exactly the family the crash dispatched into.
_LIBSTDCXX_ONLY_PREFIXES = (
    "_ZNSt6locale",       # std::locale::*   (classic(), _Impl, id::_M_id ...)
    "_ZNKSt6locale",
    "_ZNSt8ios_base",     # std::ios_base::* (Init, _M_init, failure ...)
    "_ZNKSt8ios_base",
    "_ZNSt7codecvtIDs",   # std::codecvt<char16_t, ...>
    "_ZNKSt7codecvtIDs",
    "_ZNSt7codecvtIDi",   # std::codecvt<char32_t, ...>
    "_ZNKSt7codecvtIDi",
)

SHT_DYNSYM = 11
SHT_DYNAMIC = 6
SHN_UNDEF = 0
DT_NEEDED = 1


def _module_path():
    p = getattr(ops, "__file__", None)
    assert p and os.path.isfile(p), f"cannot locate the loaded module file ({p!r})"
    return p


def _read_elf64(path):
    """Minimal ELF64 little-endian reader: returns (defined_dynsym_names, needed)."""
    with open(path, "rb") as f:
        data = f.read()
    if data[:4] != b"\x7fELF":
        return None
    ei_class, ei_data = data[4], data[5]
    if ei_class != 2 or ei_data != 1:
        pytest.skip(f"ELF but not 64-bit little-endian (class={ei_class}, data={ei_data})")
    e_shoff, = struct.unpack_from("<Q", data, 0x28)
    e_shentsize, e_shnum = struct.unpack_from("<HH", data, 0x3A)
    assert e_shentsize == 64, e_shentsize
    sections = []
    for i in range(e_shnum):
        off = e_shoff + i * e_shentsize
        (sh_name, sh_type, sh_flags, sh_addr, sh_offset, sh_size,
         sh_link, sh_info, sh_addralign, sh_entsize) = struct.unpack_from("<IIQQQQIIQQ", data, off)
        sections.append((sh_type, sh_offset, sh_size, sh_link, sh_entsize))

    def _cstr(strtab_off, idx):
        end = data.index(b"\x00", strtab_off + idx)
        return data[strtab_off + idx:end].decode("ascii", "replace")

    defined = set()
    needed = []
    for sh_type, sh_offset, sh_size, sh_link, sh_entsize in sections:
        if sh_type == SHT_DYNSYM:
            str_off = sections[sh_link][1]
            ent = sh_entsize or 24
            for off in range(sh_offset, sh_offset + sh_size, ent):
                st_name, st_info, st_other, st_shndx = struct.unpack_from("<IBBH", data, off)
                if st_shndx != SHN_UNDEF and st_name:
                    defined.add(_cstr(str_off, st_name))
        elif sh_type == SHT_DYNAMIC:
            str_off = sections[sh_link][1]
            ent = sh_entsize or 16
            for off in range(sh_offset, sh_offset + sh_size, ent):
                d_tag, d_val = struct.unpack_from("<qQ", data, off)
                if d_tag == 0:
                    break
                if d_tag == DT_NEEDED:
                    needed.append(_cstr(str_off, d_val))
    return defined, needed


def test_module_defines_no_second_libstdcxx_and_needs_the_shared_one():
    path = _module_path()
    parsed = _read_elf64(path)
    if parsed is None:
        pytest.skip(f"{os.path.basename(path)} is not ELF (PE/Mach-O): the duplicate-"
                    "libstdc++ hazard is gcc/ELF-specific")
    defined, needed = parsed

    # (3) the parser saw a real table -- guards against a vacuous pass.
    assert any(s.startswith("PyInit_") for s in defined), (
        f"no PyInit_* symbol found in .dynsym of {path}: ELF parse is broken, not the module")

    # (1) no private libstdc++ inside the module.
    leaked = sorted(s for s in defined if s.startswith(_LIBSTDCXX_ONLY_PREFIXES))
    assert not leaked, (
        f"{os.path.basename(path)} DEFINES {len(leaked)} libstdc++-internal symbol(s) -- a second, "
        "static libstdc++ is linked into the shared module (FindOpenMP + -static-libstdc++, "
        "WP-109). Two C++ runtimes in one process => wrong locale facet => segfault at the "
        f"first `opserr << int` under PythonStream. First few: {leaked[:6]}")

    # (2) the one true runtime is the shared one.
    assert any(n.startswith("libstdc++.so") for n in needed), (
        f"libstdc++.so.6 is not DT_NEEDED by {os.path.basename(path)} (NEEDED = {needed}); "
        "the module must use the system C++ runtime, never a private copy")
