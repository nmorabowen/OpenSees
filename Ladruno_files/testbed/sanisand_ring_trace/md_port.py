"""WP-128: a line-by-line Python port of the IntScheme-1 (ModifiedEuler) path
of ManzariDafalias / LadrunoSANISAND, as the WP-127 build runs it.

WHY.  F18(a) asks what an error floor `err / max(2||sigma||, sigma_ref)`
would admit.  The C++ has no such flag yet (that is WP-129), and changing the
norm changes the substep sequence, so it cannot be re-derived from a trace.
This port reproduces the C++ replay (validated by validate_port.py against
`ladrunoSANISANDReplay` on the ring states and the probe set) and then lets
the norm, TolE, dT_min and the forced-accept policy be varied OFFLINE.

Ported, in call order, with the source lines they mirror (WP-127 tree):
  integrate()                ManzariDafalias.cpp:1019  (alpha_in reversal test)
  ladrunoGuardReversalNoise  LadrunoSANISAND.cpp:2897  (P2-5 guard)
  explicit_integrator        :1144
  IntersectionFactor(_Unloading) :2827 / :2880
  ModifiedEuler              :1547
  Stress_Correction          :2924
  GetStateDependent          :5365, GetF :5049, GetNormalToYield :5308
  commitState moduli         :508 + GetElasticModuli :5143 (e_init, ElastFlag 1)
NOT ported: tangents (they do not feed the state), MaxStrain/MaxEnergy, RK,
CPPM, -implex, the stage flip, -pRe (0 in the campaign).

Conventions exactly as the C++: Voigt xx yy zz xy yz zx; stress-like
(contravariant) vectors carry tensor shear, strain-like (covariant) vectors
carry ENGINEERING shear; sigma and strain compression-positive (internal).
"""
import math

ONE3 = 1.0 / 3.0
TWO3 = 2.0 / 3.0
ROOT23 = math.sqrt(2.0 / 3.0)
SMALL = 1e-10
I1 = (1.0, 1.0, 1.0, 0.0, 0.0, 0.0)

# trace codes as ManzariDafalias::ladrunoTraceSubstep
ACCEPT, REJECT, FORCED, FORCED_CLAMP, LOWP1, LOWP2, ABANDON, CAP = range(8)


# ------------------------------------------------------------------ tensors
def tr(v):
    return v[0] + v[1] + v[2]


def dev(v):
    t = tr(v)
    return [v[0] - ONE3 * t, v[1] - ONE3 * t, v[2] - ONE3 * t, v[3], v[4], v[5]]


def add(a, b):
    return [a[i] + b[i] for i in range(6)]


def sub(a, b):
    return [a[i] - b[i] for i in range(6)]


def scal(s, a):
    return [s * x for x in a]


def contr(a, b):   # DoubleDot2_2_Contr
    r = 0.0
    for i in range(6):
        r += a[i] * b[i] + (i > 2) * a[i] * b[i]
    return r


def mixed(a, b):   # DoubleDot2_2_Mixed
    return sum(a[i] * b[i] for i in range(6))


def cov_dot(a, b):   # DoubleDot2_2_Cov
    r = 0.0
    for i in range(6):
        r += a[i] * b[i] - (i > 2) * 0.5 * a[i] * b[i]
    return r


def ncontr(a):
    return math.sqrt(contr(a, a))


def ncov(a):
    return math.sqrt(cov_dot(a, a))


def to_contra(v):
    return [v[0], v[1], v[2], 0.5 * v[3], 0.5 * v[4], 0.5 * v[5]]


def to_cov(v):
    return [v[0], v[1], v[2], 2.0 * v[3], 2.0 * v[4], 2.0 * v[5]]


def single_dot(v1, v2):
    return [v1[0] * v2[0] + v1[3] * v2[3] + v1[5] * v2[5],
            v1[3] * v2[3] + v1[1] * v2[1] + v1[4] * v2[4],
            v1[5] * v2[5] + v1[4] * v2[4] + v1[2] * v2[2],
            0.5 * (v1[0] * v2[3] + v1[3] * v2[0] + v1[3] * v2[1] + v1[1] * v2[3] + v1[5] * v2[4] + v1[4] * v2[5]),
            0.5 * (v1[3] * v2[5] + v1[5] * v2[3] + v1[1] * v2[4] + v1[4] * v2[1] + v1[4] * v2[2] + v1[2] * v2[4]),
            0.5 * (v1[0] * v2[5] + v1[5] * v2[0] + v1[3] * v2[4] + v1[4] * v2[3] + v1[5] * v2[2] + v1[2] * v2[5])]


def stiff_apply(K, G, eps_cov):   # GetStiffness(K,G) * v  (covariant -> contravariant)
    a = K + 4.0 * ONE3 * G
    b = K - 2.0 * ONE3 * G
    e = eps_cov
    return [a * e[0] + b * e[1] + b * e[2],
            b * e[0] + a * e[1] + b * e[2],
            b * e[0] + b * e[1] + a * e[2],
            G * e[3], G * e[4], G * e[5]]


class Material:
    """One LadrunoSANISAND point, campaign flags, IntScheme 1."""

    def __init__(self, params, TolF=1e-7, TolE=1e-4, Pres=0.0, Pmin=0.0101,
                 maxSubsteps=20000, reversalTol=1e-10, reversalRel=0.05,
                 err_floor=None, dT_min=1e-6, forced_policy="vanilla",
                 drag="vanilla", correction=True, alpha_err=False, fabric_err=False):
        (self.G0, self.nu, self.e_init, self.Mc, self.c, self.lamc, self.e0,
         self.ksi, self.Patm, self.m, self.h0, self.ch, self.nb, self.A0,
         self.nd, self.zmax, self.cz, self.den) = params
        self.TolF, self.TolE, self.Pres, self.Pmin = TolF, TolE, Pres, Pmin
        self.maxSub = maxSubsteps
        self.revTol, self.revRel = reversalTol, reversalRel
        # err_floor None = TODAY's norm (abs below ||sigma|| 0.5, /2||sigma|| above);
        # a number s_ref = the F18(a) proposal err/max(2||sigma||, s_ref)
        self.err_floor = err_floor
        self.dT_min = dT_min
        self.forced_policy = forced_policy   # "vanilla" | "refuse"
        # COUNTERFACTUALS (not the C++; for attribution only):
        #   drag "vanilla": a stage with dgamma < 0 moves alpha by the stress-ratio
        #                   change (C++); "frozen": that stage leaves alpha alone
        #   correction False: Stress_Correction skipped
        #   alpha_err True: the substep error also measures the alpha stages
        #                   (as RungeKutta45 does), relative to max(||alpha||, 0.5)
        self.drag = drag
        self.correction = correction
        self.alpha_err = alpha_err
        #   fabric_err True: ... and on the fabric stages (same form)
        self.fabric_err = fabric_err

    # ------------------------------------------------------------ model
    def g(self, c3):
        return 2 * self.c / ((1 + self.c) - (1 - self.c) * c3)

    def F(self, s, a):
        p = ONE3 * tr(s) + self.Pres
        x = sub(dev(s), scal(p, a))
        return ncontr(x) - ROOT23 * self.m * p

    def normal(self, s, a):
        p = ONE3 * tr(s) + self.Pres
        if abs(p) < SMALL:
            return [0.0] * 6
        x = scal(-p, a)
        x = add(x, dev(s))
        nn = ncontr(x)
        nn = 1.0 if nn < SMALL else nn
        return [xi / nn for xi in x]

    def moduli(self, s):
        pn = ONE3 * tr(s)
        pn = self.Pmin if pn <= self.Pmin else pn
        eG = self.e_init
        G = self.G0 * self.Patm * (2.97 - eG) ** 2 / (1 + eG) * math.sqrt(pn / self.Patm)
        K = TWO3 * (1 + self.nu) / (1 - 2 * self.nu) * G
        return K, G

    def state_dep(self, s, a, z, e, a_in):
        p = ONE3 * tr(s) + self.Pres
        p = SMALL if p < SMALL else p
        n = self.normal(s, a)
        aain = contr(sub(a, a_in), n)
        psi = e - (self.e0 - self.lamc * (p / self.Patm) ** self.ksi)
        c3 = math.sqrt(6.0) * tr(single_dot(n, single_dot(n, n)))
        c3 = 1.0 if c3 > 1 else c3
        c3 = -1.0 if c3 < -1 else c3
        gg = self.g(c3)
        ab = gg * self.Mc * math.exp(-self.nb * psi) - self.m
        ad = gg * self.Mc * math.exp(self.nd * psi) - self.m
        b0 = self.G0 * self.h0 * (1.0 - self.ch * e) / math.sqrt(p / self.Patm)
        d = sub(scal(ROOT23 * ad, n), a)
        b = sub(scal(ROOT23 * ab, n), a)
        h = 1.0e10 if abs(aain) < SMALL else b0 / aain
        zn = contr(z, n)
        A = self.A0 * (1 + (zn if zn > 0 else 0.0))
        D = A * contr(d, n)
        if p < 0.05 * self.Patm:
            D *= 1.0 / (1.0 + math.exp(7.6349 - 7.2713 * 101.0 / self.Patm * p))
        B = 1.0 + 1.5 * (1 - self.c) / self.c * gg * c3
        C = 3.0 * math.sqrt(1.5) * (1 - self.c) / self.c * gg
        nn = single_dot(n, n)
        R = [B * n[i] - C * (nn[i] - ONE3 * I1[i]) + ONE3 * D * I1[i] for i in range(6)]
        return dict(n=n, d=d, b=b, h=h, psi=psi, ab=ab, ad=ad, b0=b0, A=A, D=D,
                    B=B, C=C, R=R, c3=c3, aain=aain)

    # ------------------------------------------------------------ stage
    def _stage(self, s, a, z, e, a_in, dvol, ddev, K, G, strict):
        """one Heun stage: returns (dsig, dalpha, dfabric, dgamma, kind)"""
        st = self.state_dep(s, a, z, e, a_in)
        p = ONE3 * tr(s) + self.Pres
        n, b, h, B, C, D = st["n"], st["b"], st["h"], st["B"], st["C"], st["D"]
        r = [x / p for x in dev(s)]
        Kp = TWO3 * p * h * contr(b, n)
        nnn = tr(single_dot(n, single_dot(n, n)))
        temp4 = Kp + 2.0 * G * (B - C * nnn) - K * D * contr(n, r)
        if abs(temp4) < SMALL:
            return [0.0] * 6, [0.0] * 6, [0.0] * 6, None, "neutral", st, Kp
        dg = (2.0 * G * mixed(n, ddev) - K * dvol * contr(n, r)) / temp4
        neg = dg < -SMALL if strict else dg < 0.0
        if neg:
            ds = add(scal(2.0 * G, to_contra(ddev)), scal(K * dvol, I1))
            s2 = add(s, ds)
            da = scal(3.0, sub(scal(1.0 / tr(s2), dev(s2)), scal(1.0 / tr(s), dev(s))))
            if self.drag == "frozen":
                da = [0.0] * 6
            return ds, da, [0.0] * 6, 0.0, "elasticDrag", st, Kp
        mdg = dg if dg > 0 else 0.0
        nn = single_dot(n, n)
        corr = [2.0 * G * (B * n[i] - C * (nn[i] - ONE3 * I1[i])) + K * D * I1[i] for i in range(6)]
        ds = [2.0 * G * to_contra(ddev)[i] + K * dvol * I1[i] - mdg * corr[i] for i in range(6)]
        da = scal(mdg * TWO3 * h, b)
        mD = -D if -D > 0 else 0.0
        return ds, da, None, dg, "plastic", st, Kp   # fabric done by caller (needs z+dz1)

    def _dfab(self, n, zbase, dg, D):
        mdg = dg if dg > 0 else 0.0
        mD = -D if -D > 0 else 0.0
        return scal(-1.0 * mdg * self.cz * mD, add(scal(self.zmax, n), zbase))

    # ------------------------------------------------------------ ME
    def modified_euler(self, S0, E0, A0, Z0, a_in, Enext, K, G, out):
        dStrain = sub(Enext, E0)
        T, dT, dT_min, TolE = 0.0, 1.0, self.dT_min, self.TolE
        S, A, Z = list(S0), list(A0), list(Z0)
        out["meCalls"] += 1
        p = ONE3 * tr(S) + self.Pres
        if p < self.Pmin + self.Pres:
            out["entryPmin"] += 1
            S = add(dev(S), scal(self.Pmin, I1))
        while T < 1.0:
            out["nsub"] += 1
            out["substeps"] += 1
            if self.maxSub > 0 and out["nsub"] > self.maxSub:
                out["cap"] += 1
                out["trace"].append((T, dT, float("nan"), CAP, dT == dT_min, None))
                return S, A, Z, False
            e = self.e_init - (1 + self.e_init) * tr(add(scal(T, dStrain), E0))
            dvol = dT * tr(dStrain)
            ddev = scal(dT, dev(dStrain))
            ds1, da1, dz1, dg1, k1, st1, Kp1 = self._stage(S, A, Z, e, a_in, dvol, ddev, K, G, True)
            if k1 == "plastic":
                dz1 = self._dfab(st1["n"], Z, dg1, st1["D"])
            S1 = add(S, ds1)
            p = ONE3 * tr(S1) + self.Pres
            if p < self.Pres:
                if dT == dT_min:
                    out["abandon"] += 1
                    out["trace"].append((T, dT, float("nan"), ABANDON, True, None))
                    return S, A, Z, True
                out["lowp"] += 1
                out["trace"].append((T, dT, float("nan"), LOWP1, False, None))
                dT = max(0.1 * dT, dT_min)
                continue
            A1 = add(A, da1)
            Z1 = add(Z, dz1)
            ds2, da2, dz2, dg2, k2, st2, Kp2 = self._stage(S1, A1, Z1, e, a_in, dvol, ddev, K, G, False)
            # stage 2 elastic-drag alpha uses NextStress + dSigma2 vs NextStress
            if k2 == "elasticDrag" and self.drag != "frozen":
                s2 = add(S, ds2)
                da2 = scal(3.0, sub(scal(1.0 / tr(s2), dev(s2)), scal(1.0 / tr(S), dev(S))))
            if k2 == "plastic":
                dz2 = self._dfab(st2["n"], add(Z, dz1), dg2, st2["D"])
            nS = add(S, scal(0.5, add(ds1, ds2)))
            nZ = add(Z, scal(0.5, add(dz1, dz2)))
            nA = add(A, scal(0.5, add(da1, da2)))
            p = ONE3 * tr(nS) + self.Pres
            if p < self.Pres:
                if dT == dT_min:
                    out["abandon"] += 1
                    out["trace"].append((T, dT, float("nan"), ABANDON, True, None))
                    return S, A, Z, True
                out["lowp"] += 1
                out["trace"].append((T, dT, float("nan"), LOWP2, False, None))
                dT = max(0.1 * dT, dT_min)
                continue
            sn = ncontr(S)
            dd = ncontr(sub(ds2, ds1))
            if self.err_floor is None:
                err = dd if sn < 0.5 else dd / (2 * sn)
            else:
                den = max(2 * sn, self.err_floor)
                err = dd / den if den > 0 else dd
            if self.alpha_err:
                an = ncontr(A)
                ea = ncontr(sub(da2, da1))
                ea = ea if an < 0.5 else ea / (2 * an)
                err = max(err, ea)
            if self.fabric_err:
                zn = ncontr(Z)
                ez = ncontr(sub(dz2, dz1))
                ez = ez if zn < 0.5 else ez / (2 * zn)
                err = max(err, ez)
            info = dict(kinds=(k1, k2), Kp=(Kp1, Kp2), h=(st1["h"], st2["h"]),
                        aain=(st1["aain"], st2["aain"]), bn=(contr(st1["b"], st1["n"]),),
                        errabs=dd, snorm=sn, nS=nS, nA=nA, e=e,
                        dA=(da1, da2), dS=(ds1, ds2), S=S, A=A, dT=dT)
            if err > TolE:
                q = max(0.8 * math.sqrt(TolE / err), 0.1)
                if dT == dT_min:
                    if self.forced_policy == "refuse":
                        out["forced"] += 1
                        out["trace"].append((T, dT, err, FORCED, True, info))
                        return S, A, Z, False
                    S = nS
                    eta = math.sqrt(13.5) * ncontr(dev(S)) / tr(S)
                    clamped = False
                    if eta > self.Mc:
                        S = add(scal(ONE3 * tr(S), I1), scal(self.Mc / eta, dev(S)))
                        clamped = True
                    A = add(A0, scal(3.0, sub(scal(1.0 / tr(S), dev(S)), scal(1.0 / tr(S0), dev(S0)))))
                    out["forced"] += 1
                    out["clamp"] += clamped
                    out["trace"].append((T, dT, err, FORCED_CLAMP if clamped else FORCED, True, info))
                    T += dT
                else:
                    out["rej"] += 1
                    out["trace"].append((T, dT, err, REJECT, False, info))
                dT = max(q * dT, dT_min)
            else:
                out["acc"] += 1
                out["trace"].append((T, dT, err, ACCEPT, dT == dT_min, info))
                S, A, Z = nS, nA, nZ
                if self.correction:
                    S, A, Z = self.stress_correction(S0, A0, Z0, a_in, S, A, Z, e, K, G, out)
                T += dT
                q = max(0.8 * math.sqrt(TolE / err), 0.5) if err > 0 else float("inf")
                dT = max(q * dT, dT_min)
                dT = min(dT, 1 - T)
        return S, A, Z, True

    # ------------------------------------------------------------ correction
    def stress_correction(self, S0, A0, Z0, a_in, S, A, Z, e, K, G, out):
        maxIter = 50
        p = ONE3 * tr(S) + self.Pres
        if p < self.Pmin + self.Pres:
            out["corrLowP"] += 1
            p = self.Pmin + self.Pres
            return scal(p, I1), [0.0] * 6, Z
        fr = self.F(S, A)
        if abs(fr) < self.TolF:
            return S, A, Z
        out["corrCalls"] += 1
        nS, nA = list(S), list(A)
        NS, NA, NZ = S, A, Z
        for i in range(1, maxIter + 1):
            devS = dev(nS)
            st = self.state_dep(nS, nA, NZ, e, a_in)
            n, b, h, R = st["n"], st["b"], st["h"], st["R"]
            dSP = stiff_apply(K, G, to_cov(R))
            aBar = scal(TWO3 * h, b)
            r = [x / p for x in devS]
            nr = contr(n, r)
            dfs = [n[k] - ONE3 * nr * I1[k] for k in range(6)]
            dfa = scal(-p, n)
            lam = fr / (contr(dfs, dSP) - contr(dfa, aBar))
            if abs(self.F(sub(nS, scal(lam, dSP)), add(nA, scal(lam, aBar)))) < abs(fr):
                nS = sub(nS, scal(lam, dSP))
                nA = add(nA, scal(lam, aBar))
            else:
                lam = fr / contr(dfs, dfs)
                if abs(self.F(sub(nS, scal(lam, dfs)), nA)) < abs(fr):
                    nS = sub(nS, scal(lam, dfs))
                else:
                    out["corrGiveUp"] += 1          # silent return, f unchanged
                    return NS, NA, NZ
            fr = self.F(nS, nA)
            if abs(fr) < self.TolF:
                NS, NA = nS, nA
                break
            if i == maxIter:
                out["corrMaxIter"] += 1
                if self.F(S0, NA) < self.TolF:
                    dS = sub(NS, S0)
                    up, mid, down = 1.0, 0.5, 0.0
                    fo = self.F(add(S0, scal(mid, dS)), NA)
                    for jj in range(maxIter):
                        if fo < 0.0:
                            down = mid
                            mid = 0.5 * (up + mid)
                        else:
                            up = mid
                            mid = 0.5 * (down + mid)
                        fo = self.F(add(S0, scal(mid, dS)), NA)
                        if abs(fo) < self.TolF:
                            NS = add(S0, scal(mid, dS))
                            break
                else:
                    NS, NA, NZ = S0, A0, Z0
            p = ONE3 * tr(NS) + self.Pres
        return NS, NA, NZ

    # ------------------------------------------------------------ intersections
    def isect(self, S0, E0, E1, A0, a0, a1):
        inc = sub(E1, E0)
        K, G = self.moduli(S0)
        f0 = self.F(add(S0, scal(a0, stiff_apply(K, G, inc))), A0)
        f1 = self.F(add(S0, scal(a1, stiff_apply(K, G, inc))), A0)
        a = a0
        for i in range(1, 11):
            a = a1 - f1 * (a1 - a0) / (f1 - f0)
            f = self.F(add(S0, scal(a, stiff_apply(K, G, inc))), A0)
            if abs(f) < self.TolF:
                break
            if f * f0 < 0:
                a1, f1 = a, f
            else:
                f1 = f1 * f0 / (f0 + f)
                a0, f0 = a, f
            if i == 10:
                a = 0.0
                break
        if a > 1 - SMALL:
            a = 1.0
        if a < SMALL:
            a = 0.0
        return a

    def isect_unload(self, S0, E0, E1, A0):
        a0, a1 = 0.0, 1.0
        K, G = self.moduli(S0)
        dS = stiff_apply(K, G, sub(E1, E0))
        for i in range(1, 20):
            da = (a1 - a0) / 2.0
            a = a1 - da
            f = self.F(add(S0, scal(a, dS)), A0)
            if f > self.TolF:
                a1 = a
            elif f < -self.TolF:
                a0 = a
                break
            else:
                return a
        return self.isect(S0, E0, E1, A0, a0, a1)

    # ------------------------------------------------------------ one update
    def update(self, sig, alpha, alpha_in, z, e, dstrain, primed=True,
               prev_incr_norm=0.0, dt=1.0):
        """== ladrunoSANISANDReplay (compressionPositive, 3D)."""
        out = dict(meCalls=0, nsub=0, substeps=0, acc=0, rej=0, forced=0, clamp=0,
                   lowp=0, abandon=0, cap=0, entryPmin=0, pnReset=0, corrCalls=0,
                   corrGiveUp=0, corrMaxIter=0, corrLowP=0, trace=[], path=-1, rc=0)
        # committed state (replay): eps_n carries only the void ratio
        vol = (self.e_init - e) / (1.0 + self.e_init)
        En = [0.5 * vol, 0.5 * vol, 0.0, 0.0, 0.0, 0.0]
        E1 = add(En, dstrain)
        K, G = self.moduli(sig)
        # integrate(): alpha_in reversal test with mCe = Ce(K, G)
        td = stiff_apply(K, G, dstrain)
        a_in = list(alpha) if contr(sub(alpha, alpha_in), td) < 0.0 else list(alpha_in)
        # explicit_integrator
        dS = stiff_apply(K, G, dstrain)
        St = add(sig, dS)
        f = self.F(St, alpha)
        p = ONE3 * tr(St) + self.Pres
        ok = True
        if p >= self.Pres and f <= self.TolF:
            S, A, Z = St, list(alpha), list(z)
            out["path"] = 0
        else:
            fn = self.F(sig, alpha)
            pn = ONE3 * tr(sig) + self.Pres
            if pn < self.Pres:
                out["path"] = 5
                out["pnReset"] += 1
                S, A, Z = scal(self.Pmin, I1), [0.0] * 6, list(z)
            elif fn > self.TolF:
                out["path"] = 1
                S, A, Z, ok = self.modified_euler(sig, En, alpha, z, a_in, E1, K, G, out)
            elif fn < -self.TolF:
                out["path"] = 2
                ratio = self.isect(sig, En, E1, alpha, 0.0, 1.0)
                out["ratio"] = ratio
                de = scal(ratio, dstrain)
                S, A, Z, ok = self.modified_euler(add(sig, stiff_apply(K, G, de)), add(En, de),
                                                  alpha, z, a_in, E1, K, G, out)
            elif abs(fn) < self.TolF:
                n = self.normal(sig, alpha)
                nd = ncontr(dS)
                if contr(n, dS) / (1.0 if nd == 0 else nd) > -math.sqrt(self.TolF):
                    out["path"] = 3
                    S, A, Z, ok = self.modified_euler(sig, En, alpha, z, a_in, E1, K, G, out)
                else:
                    out["path"] = 4
                    ratio = self.isect_unload(sig, En, E1, alpha)
                    out["ratio"] = ratio
                    de = scal(ratio, dstrain)
                    S, A, Z, ok = self.modified_euler(add(sig, stiff_apply(K, G, de)), add(En, de),
                                                      alpha, z, a_in, E1, K, G, out)
            else:
                S, A, Z = St, list(alpha), list(z)
        # P2-5 reversal-noise guard
        if primed:
            if dt == 0.0 or ncov(dstrain) < max(self.revRel * prev_incr_norm, self.revTol):
                a_in = list(alpha_in)
        out["rc"] = 0 if ok else -3
        e1 = self.e_init - (1 + self.e_init) * tr(E1)
        out.update(sigma=S, alpha=A, alpha_in=a_in, z=Z, e=e1,
                   f_after=self.F(S, A), p=ONE3 * tr(S))
        return out
