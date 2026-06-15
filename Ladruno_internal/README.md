---
title: Ladruno internal docs index
project: Ladruno
tags:
  - index
  - internal
---

# Ladruno internal docs

Internal-only notes about the Ladruno OpenSees fork: how it was built, why specific patches exist, and the running history of changes we made on top of upstream OpenSees. Not meant for end users.

> End user looking for the build procedure? See [BUILDING.md](../BUILDING.md) at the fork root — that's the *what to run, in what order* doc. The notes here are the *why behind the decisions*.

## Contents

- [[01_compilation_journal]] — full record of the Windows build (toolchain, the eight CMakeLists patches, MUMPS LP64 saga, splash banner). Read this when something breaks and you need to remember *why* a patch is there.
- [[02_esmeralda_linux_build_guide]] — Linux build on the `esmeralda` SSH host: machine inventory, where the binaries live, the `opensees.sh` run wrapper (`TCL_LIBRARY` gotcha), proven conan+cmake recipe, and how to rebuild when `ladruno` moves.
- [[BUILD_GOTCHAS]] — the env/runtime layer around the journal: Python 3.12 ABI target, the four Windows batch traps in `setup_env.bat`, the CMake-4.3-shadow fix, MUMPS `intsize64` stickiness, the pytest bootstrap, and the installer DLL-lock.
- [[WORKFLOW_GOTCHAS]] — PR/CI/git process lessons unique to this fast-auto-merging fork: stranded commits, stale-PR "no checks", fast-check-only CI gates, the `gh --body-file -` trap, the mandatory header stamp, and the vanilla-file footprint rule.

## Conventions

- One file per topic, numbered to suggest reading order.
- Obsidian-flavored Markdown: wikilinks (`[[note]]`) for cross-references, YAML frontmatter for metadata, fenced code blocks with language tags. Plain Markdown also renders fine in any other viewer (GitHub, VS Code).
- Long-lived. When upstream OpenSees changes and we re-apply patches, **update the journal**, don't fork it.

## Companion folder

- [[../Ladruno_implementation/README|Ladruno_implementation]] — forward-looking: plans for new operators, elements, materials, recorders we want to add to this fork.
