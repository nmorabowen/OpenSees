---
wp: PR-426
title: "LadrunoRigidBody"
pr: "#426, #427, #432, #435"
status: "shipped"
section: "table"
legacy_seq: 96
---
| **LadrunoRigidBody** ([[58_ladruno_rigid_body_adr]]) — 6-DOF rigid-body Element (LS-DYNA `*CONSTRAINED_NODAL_RIGID_BODY` / Abaqus `*RIGID BODY`): a zero-stiffness Element owning a private internal 6-DOF CoM `Node` (condensed mass + body-frame inertia tensor) with rigid-link `MP_Constraint`s to each slave; body-frame momentum-conserving SO(3) integrator in `commitState` (integrated OFF the global solve); finite-rotation slave-following (`u_i=u_R+(R−I)d_i⁰` imposed in `update()`) + CoM moment gather + Housner rocking MVP. Explicit-only, serial-only v1. Full build history (P1 ballistic → P2 SO(3) → P2-S2 slaving/gather/rocking) in the section below. | Element | **ELE_TAG 33015** | `SRC/element/ladrunoRigidBody/`, `tests/test_ladrunoRigidBody_element.py` | shipped | [#426](https://github.com/nmorabowen/OpenSees/pull/426), [#427](https://github.com/nmorabowen/OpenSees/pull/427), [#432](https://github.com/nmorabowen/OpenSees/pull/432), [#435](https://github.com/nmorabowen/OpenSees/pull/435) |
