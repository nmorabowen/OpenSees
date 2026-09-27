# WP-133: PDMY03 critical-state constants, and what the PDMY "dilation brake" does (TIMs F23)

Revision 1, 2026-09-26. Not adversarially reviewed. Branch `wp/133-pdmy03-cs-params`, PR #866.
Intake: `_tims_2d_model_requests_2026-09-25.md` §F23 (on `wp/127-sanisand-replay-counters`).

## (a) New optional flags on `PressureDependMultiYield03`

```tcl
nDMaterial PressureDependMultiYield03 $tag $nd $rho $G $B $phi $gammaPeak $refP $d $PTAng \
    $mType $ca $cb $cc $cd $ce $da $db $dc \
    <$NYS=20 <$r1 $Gs1 ...>> <$liq1=1 $liq2=0 $pa=101 <$c=1.73>> \
    <-ei $e0=0.6> <-cs1 $v=0.9> <-cs2 $v=0.02> <-cs3 $v=0.7>
```

```python
ops.nDMaterial('PressureDependMultiYield03', tag, *positional, '-cs1', 0.62, '-cs2', 0.03)
```

- The four constants used to be hard-coded in the constructor (`ei = 0.6, cs1 = 0.9, cs2 = 0.02,
  cs3 = 0.7`). PDMY01 and PDMY02 already take them as positional arguments (`e`, `volLimit1..3`);
  PDMY03 did not.
- **Flags, not trailing positionals.** PDMY03's optional tail is already positional and
  variable-length: a negative `NYS` inserts `2|NYS|` backbone values in the middle of it. A fifth,
  sixth... positional after `c` would be readable only when every earlier optional is given, and
  both parsers write into a fixed `param[23]` without a bounds check, so an extra positional is an
  unchecked write past the array. Flags must come after every positional argument, in any order;
  the first recognised flag ends the positional list. An unknown option or a missing value refuses
  the command (no material is created).
- Both command surfaces: the Tcl `nDMaterial` ladder (`TclModelBuilderNDMaterialCommand.cpp`) and
  `OPS_PressureDependMultiYield03` (openseespy / openseesmp). `SRC/runtime/commands/modeling/material/
  nDMaterial.cpp` (the OpenSeesRT/xara surface) also parses PDMY03 but is not compiled by this fork's
  build; it is not changed.
- `sendSelf`/`recvSelf` already carried `einit` and `volLimit1..3` (`data(1)`, `data(13..15)`);
  `getCopy` shares `matN`, so copies read the same slot. Nothing to change there.

**Fixed along the way (same file):** the per-material constants live in static arrays indexed by
`matN`, reallocated every 20 materials (`matCount % 20 == 0`). The reallocation loop wrote
`einitx[i] = ei` (and `volLimit*x[i] = cs*`) for every EXISTING material `i`, i.e. the NEW
material's constants, and leaked the old arrays. Invisible while the constants were hard-coded (all
equal); a real cross-material leak once they are user-set. It now copies each material's own
values and frees the old arrays. PDMY02 always copied them correctly.

### Evidence

| Claim | How | Status |
|---|---|---|
| Omitted flags are byte-identical to the pre-WP-133 engine | Baselines captured from a full build of untouched `origin/ladruno` (`84fedcf13`) before any source edit: `tests/wp133_pdmy03_byteid_baseline.json` (Python, `float.hex` of every stress/strain component per step, three decks incl. user-defined backbone surfaces) and `tests/wp133_pdmy03_tcl_baseline.txt` (Tcl, `%.17g`). Compared after the change. Explicit flags equal to the defaults are identical too. | verified (G1) |
| The constants reach the brake | moving `cs1` so the path crosses the line changes the response from the crossing step on; higher `cs1` crosses later, higher `cs2`, `cs3` or `ei` earlier | verified (G2) |
| No cross-material leakage | two materials with different constants in one model reproduce their single runs bitwise, either creation order, and after 25 more PDMY03 materials force the reallocation | verified (G3) |
| The G3 reallocation gate can see the old leak | mutation build with the vanilla copy loop (`einitx[i] = ei`, ...) restored: both bitwise G3 gates FAIL, the other 12 pass; fix restored and rebuilt: 14/14 | verified (mutation) |

Gate: `tests/test_wp133_pdmy03_cs_params.py` (14 cases, ~45 s).

## (b) The dilation brake: the act's evidence and a checked reading

**The act's evidence (claimed, not re-run here).** In drained plane-strain compression, ten dense-sand
candidates across PDMY01/02/03 (with and without a retuned critical-state curve, `dilationParam3`
0.6–2.0) failed a saturation gate: the volumetric strain rate at the end never fell below 27 % of
its peak, usually above 95 %. On the strip footing the dense set drove p' under the footing from
19.7 to 1 652 kPa without a plateau.

**The act's reading, to be checked:** "PDMY's dilation brake is a switch on void ratio
(`PressureDependMultiYield.cpp:2189-2214`) that a genuinely dense sand does not reach."

### What the code does (verified by reading; PDMY01 `:2189-2214`, PDMY02 `:2195-2219`, PDMY03 `isCriticalState`)

```
e_trial = ei + eps_v,trial (1 + ei)        eps_v = tr(total strain since the material was created)
e_curr  = ei + eps_v,curr  (1 + ei)
e_cr(p) = cs1 - cs2 (p/pa)^cs3             (cs3 = 0: cs1 - cs2 ln(p/pa))
isCriticalState = 0  if e_curr and e_trial are on the SAME side of the line
                  1  otherwise (the increment crosses it)
```

Callers set the plastic potential (the dilatancy) to zero when it returns 1: PDMY01 in the dilation
branch (`:2168`) and at exit (`:2185`); PDMY02/03 at the end of `getPlasticPotential`. Three things
follow:

1. **It is a crossing detector, not a state switch.** It returns 1 only for an increment whose start
   and end sit on opposite sides of the line. Once the material is past the line, both states are on
   the same side again, the function returns 0, and dilation continues at the full rule. There is no
   branch that holds the material at critical state. *Verified numerically* on PDMY03 (G2c): with the
   line placed across the path, the volumetric rate dips only at the crossing step and is back to the
   reference rate (within 0.5 %) ten steps later.
2. **The void ratio is not the user's.** `e` starts at `ei` for every material (0.6 in PDMY03 before
   this WP, whatever the relative density the other parameters represent) and moves with the total
   strain recorded since the material was created, gravity stage included.
3. **With the default line a dense sand is far from it.** At the defaults `e_cr` = 0.894 at 19.7 kPa (18.4 % to reach it),
   0.880 at 100 kPa, 0.759 at 1 652 kPa. From `ei = 0.6` the sand has to dilate by
   `(e_cr − 0.6)/1.6`: 17.5 % volumetric strain at 100 kPa, 9.9 % at 1 652 kPa. In the WP-133 deck
   (1 element, φ = 40°, 100 kPa lateral, 5 % axial strain) `e` goes from 0.598 (minimum, early contraction) to about 0.63 at the end; the line is
   never approached, and the response is unchanged by moving `cs1` far away (G2a).

**Verdict.** The reading is right that the brake is keyed to void ratio and that a dense sand at
default constants never reaches it. It is incomplete in a way that matters for the next step:
**reaching it would not help.** Because the check fires only on the crossing increment, retuning
`ei`/`cs1..3` (now possible on PDMY03 too) moves *when* one increment loses its dilatancy; it cannot
produce a volumetric-rate plateau. That is consistent with the act's evidence that the candidates
with a retuned critical-state curve failed the saturation gate as well. Within PDMY01/02/03 no
constant produces saturation; a model with a genuine state-dependent dilatancy that vanishes at
critical state (e.g. SANISAND's `D ∝ (M_d − η)` with `M_d` a function of the state parameter ψ, or
PM4Sand) is the route.

Verified: items 1–3 on PDMY03 (source reading + gates); the PDMY01/02 functions are the same code by
reading (PDMY01 evaluates `e_cr` at `currentStress`, PDMY02/03 at `updatedTrialStress`), not run
here. Not verified: the act's ten-candidate and strip numbers (not re-run); whether a crossing-only
brake was intended by the model's authors.

## Open items

- **A two-element model stalls after the brake fires.** Two decoupled quads in one Newton system
  (`KrylovNewton`, one material with `cs1` set to cross mid-run): `analyze(1)` did not return within
  15 minutes on the step after the crossing, where each element alone runs to the end. Not
  diagnosed. Hypothesis (unverified): the zeroed plastic potential makes the tangent jump, a wild
  iterate follows, and `setSubStrainRate` sizes the substep count as `|Δε|/1e-5` with no cap, so the
  material takes millions of substeps. G3's two-element gate stops before the crossing for this
  reason.
- `pAtm` is a static member of PDMY01/02/03: the last material created sets it for all (quirk row).
- PDMY03's per-material arrays and `matCount` are process-wide and never reset on `wipe` (vanilla;
  outside this WP).
