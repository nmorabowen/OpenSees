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
import importlib.util
import os
import sys
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
    restoring this test process's sys.path / PATH afterward (the template
    mutates both as a side effect of the real boot sequence)."""
    saved_path = list(sys.path)
    saved_environ_path = os.environ.get("PATH", "")
    ns: dict = {}
    try:
        exec(compile(rendered, "<boot>", "exec"), ns)
    finally:
        sys.path[:] = saved_path
        os.environ["PATH"] = saved_environ_path
    return ns


@pytest.mark.skipif(sys.platform != "win32", reason="os.add_dll_directory is Windows-only")
def test_env_override_wins_over_baked_in_install_dir(tmp_path, monkeypatch):
    wire_venv_pth = _load_wire_venv_pth()
    installed = tmp_path / "installed"
    override = tmp_path / "override"
    installed.mkdir()
    override.mkdir()

    rendered = wire_venv_pth.BOOT_TEMPLATE.format(
        bin_repr=repr(str(installed)), mp_repr=repr("")
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
        bin_repr=repr(str(installed)), mp_repr=repr("")
    )
    monkeypatch.delenv("LADRUNO_OPENSEES_BIN", raising=False)
    ns = _run_boot_template(rendered)

    assert ns["_dirs"][0] == str(installed)
