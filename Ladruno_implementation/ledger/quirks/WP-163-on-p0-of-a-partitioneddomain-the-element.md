---
wp: WP-163
title: "On P0 of a PartitionedDomain the ELEMENT iterator also yields the ShadowSubdomains — any per-element sweep must skip isSubdomain() (WP-163)"
date: 2026-10-03
---
### On P0 of a PartitionedDomain the ELEMENT iterator also yields the ShadowSubdomains — any per-element sweep must skip `isSubdomain()` (WP-163)
- **Bites:** `PartitionedDomain::getElements()` returns the domain's own elements and then each `ShadowSubdomain`
  (`PartitionedDomainEleIter`). A Subdomain reports DOFs and external nodes but `getNodePtrs()` returns 0
  (`Subdomain.cpp`), `getMass`/`getDamp` print "DOES NOT DO ANYTHING", and `getResistingForce` round-trips to the
  remote actor. The energy kernel gathered nodal velocities through `getNodePtrs()` → null deref on P0 at the first
  record (OpenSeesSP + `recorder ladruno ... -G energy`, and the standalone EnergyBalance recorder).
  `PartitionedDomain::getElement(tag)` does NOT search subdomains, so tag-based region sweeps silently see only P0's
  own elements instead.
- **Workaround/status:** fixed in `ebkernel::addElementEnergy` (skip `isSubdomain()` / null node pointers) and the
  sizing loops (WP-163 M1). `mapElements` (`Ladruno_ElementResults.h`) and `Domain::calculateNodalReactions`
  already skip `ELE_TAG_Subdomain`; copy that guard into any new per-element loop.
