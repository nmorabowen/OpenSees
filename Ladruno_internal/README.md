---
title: Ladruno internal docs index
project: Ladruno
tags:
  - index
  - internal
---

# Ladruno internal docs

Internal-only notes about the Ladruno OpenSees fork: how it was built, why specific patches exist, and the running history of changes we made on top of upstream OpenSees. Not meant for end users.

## Contents

- [[01_compilation_journal]] — full record of the Windows build (toolchain, the eight CMakeLists patches, MUMPS LP64 saga, splash banner). Read this when something breaks and you need to remember *why* a patch is there.

## Conventions

- One file per topic, numbered to suggest reading order.
- Obsidian-flavored Markdown: wikilinks (`[[note]]`) for cross-references, YAML frontmatter for metadata, fenced code blocks with language tags. Plain Markdown also renders fine in any other viewer (GitHub, VS Code).
- Long-lived. When upstream OpenSees changes and we re-apply patches, **update the journal**, don't fork it.

## Companion folder

- [[../Ladruno_implementation/README|Ladruno_implementation]] — forward-looking: plans for new operators, elements, materials, recorders we want to add to this fork.
