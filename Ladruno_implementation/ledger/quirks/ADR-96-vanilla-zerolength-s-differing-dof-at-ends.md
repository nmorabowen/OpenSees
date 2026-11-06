---
wp: ADR-96
title: "Vanilla ZeroLength's \"differing dof at ends\" refusal CRASHED at the element command (ADR-96)"
legacy_seq: 406
---
## Vanilla `ZeroLength`'s "differing dof at ends" refusal CRASHED at the `element` command (ADR-96)

`ZeroLength::setDomain()` refused a mixed-ndf pair with a warning and a bare
`return`, leaving `t1d` NULL. `Domain::addElement()` calls `element->update()`
right after `element->setDomain()` (`Domain.cpp:493-494`), and
`ZeroLength::update()` dereferences `t1d` through `computeCurrentStrain1d()`:
access violation, upstream, on every such deck — measured on a `(2,3)` pair
while building ADR-96's G3 gate (the guide's "warn + bail" row was wrong: it
was warn + crash). The passenger-mode rotational-`-dir` refusal inherited the
same path. Both now call `ladrunoDisable()` (a zero `t1d` of the default 2-slot
width, `update()` a no-op); the warning text is unchanged and the element stays
in the domain contributing nothing. `tests/test_adr96_passenger_dof.py`
(`test_g3_rotational_dir_is_refused_on_a_passenger_pair`,
`test_g3_ndf_2_pair_is_still_refused_as_vanilla`) are the regression guards.
