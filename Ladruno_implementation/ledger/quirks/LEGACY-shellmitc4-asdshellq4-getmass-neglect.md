---
wp: LEGACY
title: "ShellMITC4 / ASDShellQ4 getMass() NEGLECT rotational inertia (zero on DOFs 4,5,6)"
legacy_seq: 142
---
### `ShellMITC4` / `ASDShellQ4` `getMass()` NEGLECT rotational inertia (zero on DOFs 4,5,6)
- **Bites:** any code that reasons about a shell node's mass PER DOF — e.g. the ADR-62 P3 shell mesh-tie,
  which keeps tied slave DOFs in the explicit equation set and so needs nonzero mass on **every** tied DOF
  (a massless tied DOF ⇒ singular projection). A shell tie that defaults to tying rotations (DOFs 4,5,6)
  will find them massless.
- **Why:** both stock shells lump **translational** mass only. `ShellMITC4::formInertiaTerms` says verbatim
  "translational mass only // rotational inertia terms are neglected"; `ASDShellQ4::getMass()` says
  "Rotational mass neglected for the moment" and only fills the translational diagonals (`index+q`, q=0..2).
  So `getMass()` returns exactly 0.0 on the rotational diagonals.
- **How to apply:** don't assume a shell's rotational DOFs carry mass. LadrunoTie's BLOCKER-2 was made
  **per-DOF** for exactly this: it names-refuses a tied rotational DOF with no mass and tells the user to add
  nodal rotary mass (`mass $node 0 0 0 mrx mry mrz`) or drop the rotations (`-dof 3 1 2 3`). The per-DOF
  element scan maps global DOF `d` of the node at external position `k` to mass diagonal
  `(Σ preceding-node getNumberDOF()) + d − 1` (standard consecutive-per-node layout). Related:
  [[project_explicit_constraint_projection]] (the handler requires lumped mass on every tied DOF).
