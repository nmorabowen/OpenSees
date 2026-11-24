---
wp: WP-253
title: "ElasticIsotropicPlaneStrain2D is Update-safe but not Form-safe: getStress writes a class-wide static"
pr: "#1147"
date: 2026-11-23
---
### A material on the Update allowlist can still race in a host that calls getStress / getInitialTangent

`ElasticIsotropicPlaneStrain2D::ladrunoThreadSafeUpdate()` is `true`, and `ladrunoThreadSafeForm()` is
`false`, for a reason: `getStress()` and `getInitialTangent()` return the class-wide statics `sigma` and
`D`. The Update contract (WP-107) covers `setTrialStrain` alone. A host that reads `getStress()` right
after `setTrialStrain()` on worker threads, as the MPM particle update does, races on that static
under the Update answer alone. Mode F's threaded material loop (WP-253) therefore requires
`ladrunoThreadSafeUpdate() && ladrunoThreadSafeForm()` on every particle material; the 2D elastic
particles run their material calls serially in particle order while the transfer loops thread.
Check the Form answer too whenever a threaded caller does more than `setTrialStrain`.

#### Also learned (WP-253): a threaded commit cannot name its refuser from the shared counter

The WP-99 refusal counter is process-wide (atomic since WP-224). The serial MPM commit names the
first refusing particle from a before/after read of that counter (WP-232k D13); under threads,
another worker can bump it between the two reads. WP-253 adds a `thread_local` tally next to it
(`ladrunoThreadCommitRefusals()`, bumped in the same place, never reset), whose before/after
difference on the committing thread is exact. The global counter's meaning is unchanged.
