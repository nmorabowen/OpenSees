---
wp: WP-179
title: "A test of the boot .pth gets the INSTALLED _ladruno_opensees_boot, not the one it wrote — any child probing site/.pth behaviour must start with -S (WP-179)"
pr: "#965"
date: 2026-10-07
---
### A test of the boot `.pth` imports the INSTALLED `_ladruno_opensees_boot`, not the one it wrote (WP-179)
- **Bites:** `tests/test_wire_venv_pth_override.py::test_alias_on_with_baked_flag` and `::test_alias_on_with_env_var` failed on the build desktop (`'OTHER' == 'FORK'`) and passed on a clean interpreter. The file is not `zone_a`, so no CI ran it, and the failure stayed hidden until WP-179's Windows test set picked it up.
- **Why:** the probe child ran `[sys.executable, "-c", PROBE]` without `-S`. `site` therefore ran the machine's own `ladruno_opensees.pth`, which imports `_ladruno_opensees_boot` from site-packages at startup. Here that was an older template with no `_alias` flag. The probe then called `site.addsitedir(<tmp>)` and `import _ladruno_opensees_boot`, which returned the **cached** installed module, so the template under test never ran. Any box where `wire_venv_pth.py` was ever run (every dev machine, the Windows runner) reads the installed boot's flag instead.
- **Workaround/status:** fixed. The probe starts with `-S`; it adds its own site dir and needs nothing from `site`. The same family as WP-175/176 (a boot `.pth` decides what a child imports): **a child process that probes import, `site` or `.pth` behaviour must start with `-S`**, or it tests the developer's machine instead of the code.
