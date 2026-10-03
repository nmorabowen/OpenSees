---
wp: LEGACY
title: "Parallel build OOM (cl.exe C1060 \"out of heap\") on the giant template TUs under RAM pressure"
legacy_seq: 41
---
### Parallel build OOM (`cl.exe` C1060 "out of heap") on the giant template TUs under RAM pressure
- **Bites:** with low free RAM (~1–2 GB of 28), `cmake --build ... -j8`/`-j16` dies
  with `fatal error C1060: compiler is out of heap space` on the huge template TUs
  (`OPS_AllASDPlasticMaterial3Ds.cpp`, `MPCORecorder.cpp`, `LadrunoRecorder.cpp`),
  and even ordinary TUs get OS-killed (ninja `FAILED: [code=2]` with no compiler
  diagnostic = the process was killed, not a code error).
- **Workaround:** compile the monsters **serially first** —
  `ninja -j1 CMakeFiles/OPS_Material.dir/SRC/material/nD/ASDPlasticMaterial3D/OPS_AllASDPlasticMaterial3Ds.cpp.obj`
  (and the two MPCO recorder objs) — then `cmake --build build\build\Release
  --target OpenSeesPy -j2` for the rest (ninja resumes cached objs). Don't assume a
  `code=2` with no error text is a code bug; check free RAM first. Learned 2026-05-31.
