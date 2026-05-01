---
title: <Feature name>
project: Ladruno
status: draft
priority: medium
owner: nmora
tags:
  - implementation
  - <category, e.g. element / material / recorder / solver>
---

# <Feature name>

## What

One-paragraph description. What does this feature do? What's the scope? What's *not* in scope (will live in a follow-up plan)?

## Why

Motivation. What user problem does this solve? Is there a paper / spec / project deadline driving it? Without this section, future-you won't remember why this was prioritized.

## Where

OpenSees source files this will touch. Reference similar existing implementations to copy patterns from.

- New code: `OpenSees/SRC/...`
- Modify: `OpenSees/SRC/...`
- Reference (similar existing impl): `OpenSees/SRC/...`
- Build: does this need a new CMake target, library, or external dep? If yes, also touch [[../Ladruno_internal/01_compilation_journal]].

## How

Design sketch. Include:

- Public API (Tcl command signature, Python binding, recorder syntax, etc.)
- Internal class hierarchy / data flow
- Integration points (what calls into the new code, what does the new code call)
- Testing: what's the minimum viable smoke test, and what's the analytical / reference solution we'll compare to?

```cpp
// Optional pseudo-code or API sketch
```

## Risks / open questions

> [!question]
> Question 1?

> [!question]
> Question 2?

- Dependencies (Conan packages, external libraries, header-only deps)
- Numerical concerns (stability, conditioning, matching reference solutions)
- ABI / build-system risks (Fortran integer size, MPI integer size, MKL threading model)
- Backwards compatibility with existing models

## Implementation log

*(filled in once the plan is being executed; move file to `Ladruno_internal/implemented_<name>.md` when done)*
