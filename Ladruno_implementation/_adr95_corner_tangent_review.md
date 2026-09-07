# ADR-95 -- UW DruckerPrager corner (Jact=(1,1)) tangent review

READ-ONLY review of `SRC/material/nD/UWmaterials/DruckerPrager.cpp`. `git blame` on the return map: `^8dc15f69d ...
fmckenna 2011` -- **vanilla upstream**, not a Ladruno edit; a fix needs a `LEDGER_vanilla_files.md` row.

## 1. Trace of the corner branch

`Jact` holds **flags 0/1 only** (474-488, 585-586, 599-605, 613-614), but the residual/Jacobian assembly (516-525,
copy at 551-560) switches on the **value**, not the index:

    516  for (int i = 0; i < 2; i++) {
    517      if (Jact(i) == 1)      { R(0) = <f1 residual>; g(0,0) = <df1/dg0>; }
    521      else if (Jact(i) == 2) { R(1) = <f2 residual>; g(1,1) = <df2/dg1>; }

`Jact(i) == 2` is **unreachable**. At the corner both i=0 and i=1 take the `==1` branch and write row 0 twice: **R(1)
stays 0** (505/550) and **g(1,1) stays at the dummy `1`** from the "initialize such that det(g)=1" lines 509-512 /
546-549. Only the off-diagonals 526-529 are right. With rho_bar = 0 (this campaign) g(1,0) = 0, so

    g = [ g00 , -9*K*rho ]   R = [ R0 , 0 ]  =>  dgamma = -inv(g)*R = [ -R0/g00 , 0 ]
        [  0  ,     1    ]                       ==> gamma(1) == 0 EXACTLY, every iterate.

Newton (536) converges on |R0| having enforced f1 only; f2 recomputed at 573 is still > 0. The corner Kuhn-Tucker test
597-609 accepts it anyway, because gTOL = -1e-10 (374) makes `gamma(1) > gTOL` TRUE for gamma(1)==0: 607-608 sets
okay=true on the first pass, `count` never grows, and the `count > 3` bailout (612-617) **cannot** fire. Matches the
measurement.

**Blow-up.** g(1,1)=1 makes the inverse dimensionally wrong: g_contra(1,1) = **+1** where it should be -1/(9K), and
g_contra(0,1) = -9*K*rho/(2G). Through 688-689:

    688  temp1 = -n - (3*K*rho/(2G))*I1 - (27*K*K*rho/(2G))*I1   <-- last term spurious
    689  temp2 = 3*K*I1                                          <-- should be -(1/3)*I1

into the assembly 694-697:

    695   + 3*K*I1(i)*temp2(j) = +9*K^2 * I1(x)I1        (correct: -K * I1(x)I1)
    696   + 2*G*n(i)*temp1(j) -> -27*rho*K^2 * n(x)I1    (correct: -2G*n(x)n only)

Two spurious rank-one terms in **stress^2 units where a stress belongs**; relative to 2G they scale as 9K^2/(2G) and
27*rho*K^2/(2G) -- the numeric value of K in the model's stress unit, i.e. the 5-8 decades seen in det(A)/(2G)^3
(-3e5..-6.6e7 vs -0.06). They are non-symmetric with different left vectors (I1 vs n), so det(A) gets a rank-2
indefinite perturbation: large, negative, orientation-dependent (the ~200x spread over the 4 GPs). `NormCep < 1e-10`
(702) is a floor, never a ceiling. Correct corner algebra (theta=0,H=0,delta2=0 => g11 = -9K): det g = 18GK, temp1 =
-n, temp2 = -(1/3)*I1, giving `Cep = 2G*(1 - 2G*gamma0/||eta_trial||)*(IIdev - n(x)n)` -- volumetric part cancels
exactly (both surfaces constrain I1), tangent O(2G) = healthy-GP scale.

## 2. Verdict: (a) mathematically wrong

A coding defect, not vertex conditioning. Two riders. (i) The **f2-only** branch Jact=(0,1) is broken the same way:
i=1 hits `==1` and assembles the **f1** residual into R(0)/g(0,0), so a pure cutoff step does a cone return with
gamma(1)=0 and the same `+9K^2 I1(x)I1` tangent. (ii) There **is** a 1/norm_eta at 697, `-
4*G*G/norm_eta*gamma(0)*(IIdev - n(x)n)`, and it uses the **final** norm (recomputed 673) while `n` is the **trial**
normal (497-502; eta is never updated inside the loop) -- radial return requires ||eta_trial|| there. Today it is
merely large; once (a) is fixed the corner return lands exactly on the apex (||eta_final||=0 when Kiso=sigma_y) and
697 becomes a **divide by zero**. Fix both or you trade a big number for a NaN. Line 498 guards `n` against
norm_eta->0; 697 has no such guard.

## 3. Can the cone branch leave a GP at f2 > 0?

Yes -- it is the observed path. With rho_bar = 0 a cone return changes I1 by `-9*K*rho_bar*gamma(0) = 0` (665): the
f1-only map **cannot move I1 at all**, so I1_tr > T stays > T. 581-587 detects `f2 >= fTOL`, promotes to the corner,
count += 1 -- the corner returns gamma(1)=0, leaves I1 untouched, and is accepted at 607-608. The GP commits at f2 > 0
with the cutoff violated, reporting branch 3 (`mLadBranch = Jact(0)+2*Jact(1)`, 716) and gamma1 == 0.

## 4. Options for P4 (ranked)

1. **(ii) Fix the corner/apex tangent** -- rewrite 516-525 and 551-560 index-driven:
   `if (Jact(0)==1){R(0)=..;g(0,0)=..;} if (Jact(1)==1){R(1)=..; g(1,1) = -9*mK +
   mdelta2*T(alpha2);}`, and divide 697 by the **trial** norm_eta (latch it before 673).
   *Risk:* vanilla algebra used by every DP model, and 697 changes f1-only runs too -- the gate
   must expect a real (correct) delta, not bit-identity.
2. **(iv) Elastic fallback at corner GPs** (`if (Jact(0)&&Jact(1)) mCep = mCe;`). *Risk:*
   Newton-safe and cheap, but the return map stays wrong -- it converges quietly onto a stress
   that violates the cutoff.
3. **(i) `-noTensionCutoff` / user-settable T** -- the `mTo = 1e10` path exists (151-153, 965,
   976). *Risk:* deletes the physics, and for rho>0 the cone apex is still at I1 = T, so 697 is
   singular there anyway.
4. **(iii) Apex smoothing (hyperbolic DP / Abbo-Sloan)** -- the fork already ships smooth pieces
   (`ASDPlasticMaterial3D/YieldFunctions/DruckerPrager_YF.h`,
   `PlasticFlowDirections/MohrCoulombTensionCutoff_PF.h`). *Risk:* new material + new
   calibration; right long-term, far too much for P4.

**Recommend 1, with 4 as the eventual replacement.** Only 1 makes the reported stress admissible; 2 merely makes
Newton quiet.

## 5. Falsifier (one GP driven to I1 > T)

    1. DruckerPrager3D, phi=20 (rho=0.2), rho_bar=0, SY=0.2 (T=0.816), H=0; hydrostatic
       tension strain ramp until f2_trial > 0.  Read ladrunoBranch.
    2. assert branch  == 3          # [0]  corner
    3. assert gamma1  > 1e-12       # [2]  today exactly 0.0
    4. assert I1_ret  <= T + 1e-8   # [6]  cutoff actually enforced
    5. assert abs(detAmin) < 10.0   # [7]  today 3e5 .. 6.6e7
    6. assert ||Cep - 2G*(1-2G*g0/||eta_tr||)*(IIdev-n(x)n)|| / ||Ce|| < 1e-10
