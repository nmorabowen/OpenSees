---
wp: LEGACY
title: "A settle / convergence gate written as \"the output stopped changing over a chunk\" measures the chunk length, not the convergence"
legacy_seq: 287
---
### A settle / convergence gate written as "the output stopped changing over a chunk" measures the chunk length, not the convergence
- **Bites:** you gate a relaxation loop on `|R − R_prev|/R < tol` evaluated once per `analyze(chunk, dt)` call. It passes. You shorten `chunk` to get finer reporting and it passes *sooner*, on a state that is *less* converged — because less happens inside a shorter chunk. The verdict tracks a reporting parameter.
- **Why:** a per-chunk *change* is an increment, and increments scale with the interval they are measured over. A *residual* does not.
- **Rule:** gate on a chunk-free quantity — the true static unbalance `‖f_ext − f_int‖_∞` (`ladrunoDR residualNorm`, which DR computes from the pre-damping solved acceleration and is exactly `‖M*·a‖_∞`) — and demote the change measure to a reported diagnostic. The same rule applies to any "it stopped moving" stall detector whose window is a tunable. *Learned 2026-08-11 (note 83 §1.1).*
