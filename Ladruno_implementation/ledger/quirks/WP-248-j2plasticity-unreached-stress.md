---
wp: WP-248
title: "J2Plasticity returns a stress after its local iteration fails"
date: 2026-11-23
pr: "#1142"
---
### J2Plasticity prints "More than 25 iterations" and still returns a stress (a stress target beyond yield looks reached)

Stock `J2Plasticity` (`SRC/material/nD/J2Plasticity.cpp`) caps its return-mapping
iteration at 25; past that it prints `More than 25 iterations in constituive
subroutine J2-plasticity` and returns the unconverged state with a success code.
Under stress control beyond the perfectly plastic plateau (`ndTest mixed ssssss
500 ...` with σy = 420) the outer Newton therefore sees a stress it accepts, and
nothing is refused. Its converged plateau is also only good to about 1e-8
relative (420.0000048 for σy = 420). Found writing the WP-248 deck, whose oracles
measure unloading from the plateau actually reached and test non-convergence
with `-maxIter 1 -tol 1e-300` instead of relying on J2 to refuse. A caller that
must not trust such a state needs a material that refuses, or its own check of
the material's warnings.
