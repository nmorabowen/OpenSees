# CMS profiling deck (ADR-1000 P4 §1)

`cms_profile_chain.tcl` — the physical smoke scaled to a size where the hierarchy
phase timings are resolvable. The smoke is 4 elements per rank: correct, but far
too small to profile.

```bash
mpiexec -n 4 OpenSeesMP.exe cms_profile_chain.tcl
```

- Requires an `OpenSeesMP` built with `LADRUNO_CMS=ON`, and exactly 4 ranks.
- Size: `LADRUNO_CMS_PROFILE_ELEMENTS` (default 2000 elements per rank).
- Running it straight out of the build tree needs `TCL_LIBRARY` pointed at the
  conan `lib/tcl8.6`, else the interpreter starts with no commands
  (`invalid command name "wipe"`). The curated `dist/` resolves this.

The deck checks its first eigenvalue against the analytic chain spectrum and
errors out on mismatch, so a profile can never be taken from a run that silently
produced garbage.

Results and interpretation live in [[1000_ladruno_cms_adr]] §27. Two cautions
recorded there and repeated here:

- it is a **1-D chain**, the friendliest topology (1-DOF interfaces), so `S2`
  and `T1` have almost nothing to do — do not generalise the phase mix to models
  with real interfaces;
- at 8000 elements per rank it **fails** rather than profiles: the local
  fixed-interface Lanczos exhausts `maximumRestarts`, which is hard-coded to 20
  and exposed by no command option.
