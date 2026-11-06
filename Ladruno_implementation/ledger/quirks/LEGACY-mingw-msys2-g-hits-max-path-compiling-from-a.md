---
wp: LEGACY
title: "MinGW/MSYS2 g++ hits MAX_PATH compiling from a long Windows path -- compile syntax-check probes from a short path (C:\\tmp\\...), never from the Claude scratchpa…"
legacy_seq: 420
---
### MinGW/MSYS2 g++ hits `MAX_PATH` compiling from a long Windows path -- compile syntax-check probes from a short path (`C:\tmp\...`), never from the Claude scratchpad path
- **Bites:** the per-session scratchpad directory (`C:\Users\<user>\AppData\Local\Temp\claude\<long-encoded-project-path>\<uuid>\scratchpad\...`) routinely exceeds Windows' legacy 260-character `MAX_PATH`. MinGW/MSYS2's `g++` (used for the header-only `-fsyntax-only` pre-flight, since it needs no OpenSees link) can fail to open its own intermediate files, or fail more confusingly deep inside a header include chain, when invoked with a working directory or `-I` path near that limit -- and the resulting error looks like a header/include problem, not a path-length problem.
- **Workaround/status:** compile every g++ pre-flight probe from a short path (`C:\tmp\p5preflight\...` or similar), never directly inside the scratchpad tree. Only the STAGING scripts (the `.py` fix appliers) need to live in the scratchpad; the actual g++ invocation's cwd and `-o` target should not.
