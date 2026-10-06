"""`wire_venv_pth.py`'s generated boot module must let `LADRUNO_OPENSEES_BIN`
override the baked-in install dir at RUNTIME, not just at wire-time.

Born from the same trap as [[LEDGER_quirks]] "An INSTALLED Ladruno hijacks
`import opensees` in every venv it has wired": the generated
`_ladruno_opensees_boot.py` eagerly `import opensees`s from a path baked in
when `wire_venv_pth.py` last ran, before any user code executes -- so a
worktree build's own `sys.path.insert` bootstrap is a no-op in a venv the
Ladruno installer has wired. `LADRUNO_OPENSEES_BIN` / `LADRUNO_OPENSEESMP_BIN`
are the escape hatch. This test renders the actual `BOOT_TEMPLATE` (the
literal string a real venv gets) and executes it in an isolated namespace,
so a future edit that silently drops the override cannot pass unnoticed.

No built engine needed -- this only checks which DIRECTORY the boot module
picks, not that `opensees.pyd` actually imports from it.
"""
import ast
import importlib
import importlib.util
import os
import subprocess
import sys
import types
from pathlib import Path

import pytest

_WIRE_VENV_PTH = (
    Path(__file__).resolve().parents[1] / "Ladruno_scripts" / "wire_venv_pth.py"
)


def _load_wire_venv_pth():
    spec = importlib.util.spec_from_file_location("wire_venv_pth", _WIRE_VENV_PTH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _run_boot_template(rendered: str) -> dict:
    """Execute the rendered boot module body in an isolated namespace,
    restoring this test process's sys.path / PATH / meta_path afterward (the
    template mutates all three as a side effect of the real boot sequence).

    sys.meta_path matters since the alias went lazy: the boot module installs a
    finder there, and a test that left it behind would arm every LATER test in
    the session with a hook that resolves `openseespy` to whatever `opensees`
    happens to be importable."""
    saved_path = list(sys.path)
    saved_meta = list(sys.meta_path)
    saved_environ_path = os.environ.get("PATH", "")
    ns: dict = {}
    try:
        exec(compile(rendered, "<boot>", "exec"), ns)
    finally:
        sys.path[:] = saved_path
        sys.meta_path[:] = saved_meta
        os.environ["PATH"] = saved_environ_path
    return ns


def _engine_imports_outside_functions(tree: ast.AST) -> list:
    """Every `import opensees[...]` that is NOT inside a function body.

    Inside a function is the whole point (lazy); at module level it is the bug —
    that single statement is what pinned the install's DLLs in every interpreter
    the venv started, and what made `sys.path.insert` a no-op.
    """
    found: list = []

    class _V(ast.NodeVisitor):
        def visit_FunctionDef(self, node):      # deliberately do NOT descend
            return
        visit_AsyncFunctionDef = visit_FunctionDef

        def visit_Import(self, node):
            for alias in node.names:
                if alias.name.split(".")[0] in ("opensees", "openseesmp"):
                    found.append(alias.name)

        def visit_ImportFrom(self, node):
            if (node.module or "").split(".")[0] in ("opensees", "openseesmp"):
                found.append(node.module)

    _V().visit(tree)
    return found


@pytest.mark.skipif(sys.platform != "win32", reason="os.add_dll_directory is Windows-only")
def test_env_override_wins_over_baked_in_install_dir(tmp_path, monkeypatch):
    wire_venv_pth = _load_wire_venv_pth()
    installed = tmp_path / "installed"
    override = tmp_path / "override"
    installed.mkdir()
    override.mkdir()

    rendered = wire_venv_pth.BOOT_TEMPLATE.format(
        bin_repr=repr(str(installed)), mp_repr=repr(""), alias_repr=repr(False)
    )
    monkeypatch.setenv("LADRUNO_OPENSEES_BIN", str(override))
    ns = _run_boot_template(rendered)

    assert ns["_dirs"][0] == str(override), (
        "the boot module must prefer LADRUNO_OPENSEES_BIN over the dir "
        "wire_venv_pth.py baked in -- without this, an installer-wired venv "
        "silently hijacks `import opensees` to the installed build"
    )


def test_no_override_falls_back_to_baked_in_install_dir(tmp_path, monkeypatch):
    wire_venv_pth = _load_wire_venv_pth()
    installed = tmp_path / "installed"
    installed.mkdir()

    rendered = wire_venv_pth.BOOT_TEMPLATE.format(
        bin_repr=repr(str(installed)), mp_repr=repr(""), alias_repr=repr(False)
    )
    monkeypatch.delenv("LADRUNO_OPENSEES_BIN", raising=False)
    ns = _run_boot_template(rendered)

    assert ns["_dirs"][0] == str(installed)


# --------------------------------------------------------------------------
# PASSIVE wiring (#735) — the generated boot module must not load the engine at
# interpreter startup. One module-level `import opensees` cost us twice:
#   * it pinned the install's MKL DLLs in every venv interpreter (VS Code's
#     linters, Jupyter), so the next installer upgrade died with
#     "DeleteFile failed; code 5" (BUILD_GOTCHAS §5);
#   * it pinned sys.modules before user code ran, which is the hijack the two
#     tests above exist to work around.
# --------------------------------------------------------------------------
def test_boot_template_has_no_module_level_engine_import():
    """The invariant, checked structurally so it cannot regress unnoticed."""
    wire_venv_pth = _load_wire_venv_pth()
    rendered = wire_venv_pth.BOOT_TEMPLATE.format(
        bin_repr=repr(r"C:\some\bin"), mp_repr=repr(""), alias_repr=repr(False)
    )
    offenders = _engine_imports_outside_functions(ast.parse(rendered))
    assert offenders == [], (
        f"the boot module imports {offenders} outside a function, i.e. at "
        "interpreter startup. That pins the install's DLLs (installer upgrades "
        "then fail with DeleteFile code 5) and hijacks sys.modules before user "
        "code runs. Move it inside the lazy meta_path finder."
    )


@pytest.mark.skipif(sys.platform != "win32", reason="os.add_dll_directory is Windows-only")
def test_openseespy_alias_resolves_lazily(tmp_path, monkeypatch):
    """...and behaviourally: nothing aliased at exec, but the name still works.

    Uses a stand-in `opensees` module, so this needs no built engine and cannot
    be fooled by a real one already sitting in sys.modules.
    """
    wire_venv_pth = _load_wire_venv_pth()
    rendered = wire_venv_pth.BOOT_TEMPLATE.format(
        bin_repr=repr(str(tmp_path)), mp_repr=repr(""), alias_repr=repr(True)
    )
    monkeypatch.delenv("PMI_RANK", raising=False)
    monkeypatch.delenv("PMI_SIZE", raising=False)

    fake = types.ModuleType("opensees")
    monkeypatch.setitem(sys.modules, "opensees", fake)
    for name in ("openseespy", "openseespy.opensees"):
        monkeypatch.delitem(sys.modules, name, raising=False)

    saved_path = list(sys.path)
    saved_meta = list(sys.meta_path)
    saved_environ_path = os.environ.get("PATH", "")
    try:
        exec(compile(rendered, "<boot>", "exec"), {})

        # PASSIVE: executing the boot module aliases nothing yet.
        assert "openseespy" not in sys.modules, (
            "the boot module resolved the openseespy alias at startup; it must "
            "wait until something actually asks for the name"
        )

        # ...but the name resolves on demand, to the same object, by both routes.
        assert importlib.import_module("openseespy") is fake
        assert importlib.import_module("openseespy.opensees") is fake
        assert fake.opensees is fake, "`openseespy.opensees` attribute access"
    finally:
        sys.path[:] = saved_path
        sys.meta_path[:] = saved_meta
        os.environ["PATH"] = saved_environ_path
        for name in ("openseespy", "openseespy.opensees"):
            sys.modules.pop(name, None)


@pytest.mark.skipif(sys.platform != "win32", reason="os.add_dll_directory is Windows-only")
def test_lazy_alias_is_skipped_under_mpi(tmp_path, monkeypatch):
    """Under an MPI rank the alias must not be installed at all — aliasing the
    sequential runtime into an `import openseesmp` rank would put two OpenSees
    runtimes in one process. Intel MPI / Hydra set PMI_RANK per rank."""
    wire_venv_pth = _load_wire_venv_pth()
    rendered = wire_venv_pth.BOOT_TEMPLATE.format(
        bin_repr=repr(str(tmp_path)), mp_repr=repr(""), alias_repr=repr(True)
    )
    monkeypatch.setenv("PMI_RANK", "0")
    before = list(sys.meta_path)
    _run_boot_template(rendered)
    assert list(sys.meta_path) == before, (
        "an MPI rank installed the openseespy alias finder anyway"
    )


# --------------------------------------------------------------------------
# F2-c: the openseespy alias is OPT-IN. Default wiring must leave `openseespy`
# alone (apeGmsh #1055: tests silently ran against the installed fork binary).
# These run in a fresh subprocess with a stand-in `opensees` module on the
# generated boot module's bin dir, so no built engine is needed.
# --------------------------------------------------------------------------
_PROBE = (
    "import site, sys; site.addsitedir(sys.argv[1]); "
    "import _ladruno_opensees_boot; "
    "import importlib; "
    "m = importlib.import_module('openseespy.opensees'); "
    "print('FORK' if getattr(m, 'IS_FAKE_FORK', False) else 'OTHER')"
)


def _wire_into(tmp_path, alias):
    site_dir = tmp_path / "site"
    bin_dir = tmp_path / "bin"
    site_dir.mkdir()
    bin_dir.mkdir()
    (bin_dir / "opensees.py").write_text("IS_FAKE_FORK = True\n")
    # a stand-in "real" openseespy so the non-alias import has something to find
    pkg = tmp_path / "real" / "openseespy"
    pkg.mkdir(parents=True)
    (pkg / "__init__.py").write_text("")
    (pkg / "opensees.py").write_text("IS_FAKE_FORK = False\n")
    wire = _load_wire_venv_pth()
    (site_dir / "_ladruno_opensees_boot.py").write_text(
        wire.BOOT_TEMPLATE.format(
            bin_repr=repr(str(bin_dir)), mp_repr=repr(""), alias_repr=repr(alias),
        ),
        encoding="ascii",
    )
    return site_dir, tmp_path / "real"


def _probe(tmp_path, env_extra, alias=False):
    site_dir, real = _wire_into(tmp_path, alias)
    env = {k: v for k, v in os.environ.items()
           if k not in ("PMI_RANK", "PMI_SIZE", "LADRUNO_OPENSEESPY_ALIAS",
                        "LADRUNO_OPENSEES_BIN", "PYTHONPATH")}
    env["PYTHONPATH"] = str(real)
    env.update(env_extra)
    out = subprocess.run(
        [sys.executable, "-c", _PROBE, str(site_dir)],
        capture_output=True, text=True, env=env, timeout=60,
    )
    assert out.returncode == 0, out.stderr
    return out.stdout.strip()


def test_alias_off_by_default(tmp_path):
    assert _probe(tmp_path, {}) == "OTHER"


def test_alias_on_with_baked_flag(tmp_path):
    assert _probe(tmp_path, {}, alias=True) == "FORK"


def test_alias_on_with_env_var(tmp_path):
    assert _probe(tmp_path, {"LADRUNO_OPENSEESPY_ALIAS": "1"}) == "FORK"


def test_main_flag_is_baked_into_generated_boot(tmp_path, monkeypatch):
    import sysconfig
    wire = _load_wire_venv_pth()
    monkeypatch.setattr(sysconfig, "get_paths",
                        lambda *a, **k: {"purelib": str(tmp_path)})
    monkeypatch.setattr(sys, "argv",
                        ["wire_venv_pth.py", "--alias-openseespy", "BINDIR"])
    wire.main()
    boot = (tmp_path / "_ladruno_opensees_boot.py").read_text()
    assert "_alias = True or" in boot and "'BINDIR'" in boot
    monkeypatch.setattr(sys, "argv", ["wire_venv_pth.py", "BINDIR"])
    wire.main()
    assert "_alias = False or" in (tmp_path / "_ladruno_opensees_boot.py").read_text()
