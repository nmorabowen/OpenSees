---
wp: ADR-96
title: "FE_Element::setID() is greedy: it copies EVERY equation of every DOF_Group (ADR-96)"
legacy_seq: 404
---
## `FE_Element::setID()` is greedy: it copies EVERY equation of every DOF_Group (ADR-96)

`SRC/analysis/fe_ele/FE_Element.cpp` `setID()` walks `myDOF_Groups`, copies each
group's full `getID()` into `myID` and returns `-3` the moment it runs past `numDOF`
— leaving a half-filled map behind. A handler-level adapter with a fixed
per-node slot count (the contact `LadrunoContactFE`, `3·(1+n_ps)`) therefore
cannot connect an ndf-4 (`LadrunoUP`) or ndf-6 node through the base method: the
pressure/rotation equation lands in a translation slot. `numDOF` and `theModel`
are **private** in `FE_Element`, so an override must size itself from
`myID.Size()` and reach the groups through `Node::getDOF_GroupPtr()`. Fix pattern:
`LadrunoContactFE::setID()` (ADR-96). Domain elements do not hit this because
`FE_Element(ele)` sizes `numDOF` from `ele->getNumDOF()` — which is why a
mixed-ndf `ZeroLength` has to REPORT the element size (`dofNd1 + dofNd2`) and
scatter its core, not just relax its count check.
