"""The ENERGY plug (task item 2): the only place where the harness knows which elastic energy the oracle uses.

The rest of the harness (drivers, objective, fit, smoke) asks for an energy BY NAME and gets back the oracle
parameters of the elastic part as common names (tests/conftest.py names; model.py renames per oracle). The
energy choice (BA06 vs HAR, plan §2.5) is made (owner decision (a) 2026-10-02: HAR, n = 1/2, for TIMs/Toyoura; BA06
stays the default and the paper mode), so nothing outside this module depends on it.

Contract of a plug (class with these members):
  name                           registry key
  available(oracle) -> (ok, why) does the oracle module ('O1'/'O2') implement this energy?
  params(sand, p_init, e_init, policy) -> dict
                                 the elastic parameters (common names) for one element test, from the sand's
                                 ElasticTargets (G(p, e), nu), the test's initial mean stress magnitude p_init
                                 (kPa > 0) and void ratio e_init, and the ElasticPolicy
  oracle_extra(oracle) -> dict   the energy switch of the oracle's Params ({'energy': ...})
  describe(sand, p_init, e_init, policy) -> dict
                                 the same numbers plus the stiffness they imply, for the report
  tcl_flags(sand, policy, e_ref) -> list
                                 the same energy as `nDMaterial LadrunoNorSand` flags for a BVP deck (one material:
                                 policy 'global'); [flag, value, ...] pairs, p_a EXCLUDED (it is the CSL's -p_a)

Sheet references: §2.1 (interface), §2.2 (BA06: K = -p/kappa_hat, constant mu0 for alpha0 = 0), §2.3 (HAR, (S.4h)-
(S.5h'')), §2.4 (the option bookkeeping: p_a is ONE flag shared by HAR and the fork CSL; the BA06 values are REFUSED
under HAR, never ignored), §15 (parameter transfer: BA06 mu0 = G(p_rep, e), kappa_hat = p_rep / K(p_rep);
HAR g, k from G0, nu, e_ref with p_a = p_atm).

HAR mapping (axis q = 0, sheet §2.3 "K = k p_a (|p|/p_a)^n, G = g p_a (|p|/p_a)^n"): the DM04 targets have G, K both
proportional to p^(1/2) and K/G constant, so n = the targets' exponent (0.5), g = G(p_a, e_ref)/p_a and
k = K(p_a, e_ref)/p_a reproduce G(p, e_ref) and K(p, e_ref) at EVERY p with no representative pressure, and the axis
Poisson ratio (3k - 2g)/(6k + 2g) equals nu. p_a is the shared `-p_a` flag: the sand's p_a (101 kPa for TIMs, sand.TIMS_P_A).
"""
from __future__ import annotations

import dataclasses
from dataclasses import dataclass

from .sand import Sand


class EnergyUnavailable(RuntimeError):
    pass


@dataclass(frozen=True)
class ElasticPolicy:
    """Where the constant-coefficient energies (BA06) are matched to the pressure-dependent targets.
    mode 'per_test': p_rep = factor * p_init and e_rep = e_init of each test (the element fits; every test
                     gets the stiffness of its own confinement, which a single BVP material cannot have).
    mode 'global'  : one p_rep / e_rep for every test (what a BVP needs: one material). e_rep None = test e.
    HAR ignores p_rep (its stiffness scales with p^n itself) and uses e_rep (or the test's e) as e_ref;
    a 'global' policy for HAR therefore needs only e_rep. Build the BVP policy with ElasticPolicy.bvp()."""
    mode: str = "per_test"
    factor: float = 1.0
    p_rep: float | None = None
    e_rep: float | None = None

    @classmethod
    def bvp(cls, p_rep: float | None = None, e_rep: float | None = None) -> "ElasticPolicy":
        """The 'global' policy for a BVP (one material for the whole soil body). BA06 needs p_rep (kPa > 0: the
        pressure at which the constant mu0, kappa_hat match the targets); e_rep (the void ratio of the body) is
        needed by both energies for a BVP, because the targets' G depends on e. HAR ignores p_rep."""
        return cls(mode="global", p_rep=p_rep, e_rep=e_rep)

    def rep(self, p_init: float, e_init: float):
        if self.mode == "per_test":
            return self.factor * p_init, e_init
        if self.mode == "global":
            if self.p_rep is None:
                raise ValueError("ElasticPolicy('global') needs p_rep")
            return float(self.p_rep), (e_init if self.e_rep is None else float(self.e_rep))
        raise ValueError(f"unknown ElasticPolicy mode {self.mode!r}")

    def e_ref(self, e_init: float | None) -> float:
        """The void ratio the targets are evaluated at (HAR: no p_rep needed)."""
        if self.mode not in ("per_test", "global"):
            raise ValueError(f"unknown ElasticPolicy mode {self.mode!r}")
        e = e_init if self.e_rep is None else self.e_rep
        if e is None:
            raise ValueError("ElasticPolicy needs e_rep (or the test's e)")
        return float(e)


def _oracle_params_cls(oracle: str):
    from . import model
    return model.oracle_module(oracle).Params


def _fields(oracle: str):
    return {f.name for f in dataclasses.fields(_oracle_params_cls(oracle))}


class BA06:
    """Borja & Andrade 2006 energy (S.4), alpha0 = 0: K = -p/kappa_hat, G = mu0 (constant).
    Matched at (p_rep, e_rep): mu0 = G(p_rep, e_rep), kappa_hat = p_rep / K(p_rep, e_rep) (sheet §15).
    p0 = -p_rep, eps_v0 = 0: the energy's reference point is the representative state (any other choice
    only shifts eps^e, sheet S.4)."""
    name = "BA06"

    def available(self, oracle: str):
        fields = _fields(oracle)
        if "energy" in fields:
            return True, "oracle has an energy switch; BA06 selected explicitly"
        return True, "oracle energy is BA06 (no switch)"

    def params(self, sand: Sand, p_init: float, e_init: float, policy: ElasticPolicy) -> dict:
        p_rep, e_rep = policy.rep(p_init, e_init)
        G = sand.elastic.G(p_rep, e_rep)
        K = sand.elastic.K(p_rep, e_rep)
        out = dict(p0=-p_rep, kappa_hat=p_rep / K, eps_v0=0.0, mu0=G, alpha0=0.0)
        return out

    def oracle_extra(self, oracle: str) -> dict:
        return {"energy": "BA06"} if "energy" in _fields(oracle) else {}

    def describe(self, sand: Sand, p_init: float, e_init: float, policy: ElasticPolicy) -> dict:
        p_rep, e_rep = policy.rep(p_init, e_init)
        d = self.params(sand, p_init, e_init, policy)
        d.update(energy=self.name, p_rep=p_rep, e_rep=e_rep, G_target_at_rep=sand.elastic.G(p_rep, e_rep),
                 K_target_at_rep=sand.elastic.K(p_rep, e_rep),
                 note="K ∝ p (exponent 1), G constant: matches the sqrt(p) targets only at p_rep")
        return d

    def tcl_flags(self, sand: Sand, policy: ElasticPolicy, e_ref: float | None = None) -> list:
        """-p0 -kappa_hat -mu0 at the policy's (p_rep, e_rep) (no -energy flag: BA06 is the material's default,
        so the list also works on a build that predates -energy)."""
        if policy.mode != "global":
            raise ValueError("tcl_flags: a deck material needs ElasticPolicy.bvp(p_rep, e_rep)")
        d = self.params(sand, 1.0, e_ref if e_ref is not None else 0.0, policy)
        return [["-p0", d["p0"]], ["-kappa_hat", d["kappa_hat"]], ["-mu0", d["mu0"]]]


class HAR:
    """Houlsby, Amorosi & Rojas (2005) energy (sheet §2.3, option `-energy HAR`): Psi (S.4h), stiffness ∝ p^n with
    the same exponent for K and G. n = the targets' exponent (0.5); k = K(p_a, e_ref)/p_a, g = G(p_a, e_ref)/p_a,
    so G(p, e_ref), K(p, e_ref) follow the DM04 targets at every p with no representative pressure. p_a is the shared
    `-p_a` flag = sand.p_a (TIMs 101 kPa). Under HAR the five BA06 values are REFUSED by both oracles (sheet §2.4),
    so params() returns none of them; p0 := -p_a and the p_min default 5e-3 p_a come from the oracle."""
    name = "HAR"
    # HAR constant -> oracle Params field name (O1 and O2 use the same names, round 3b).
    O2_FIELDS = {"n": "n_e", "G_r": "g", "K_r": "k"}

    def available(self, oracle: str):
        fields = _fields(oracle)
        if "energy" not in fields:
            return False, f"{oracle} has no energy switch (BA06 only)"
        missing = [v for v in self.O2_FIELDS.values() if v not in fields]
        if missing:
            return False, f"{oracle} has an energy switch but not the HAR fields {missing} (bind HAR.O2_FIELDS)"
        return True, "HAR fields present"

    def constants(self, sand: Sand, e_ref: float) -> dict:
        """(n, p_a, G_r, K_r, g, k) at the sand's shared p_a: G_r = G(p_a, e_ref) kPa, g = G_r / p_a."""
        n = sand.elastic.n_exp
        if not (0.0 <= n < 1.0):
            raise ValueError(f"HAR needs 0 <= n < 1, the targets' exponent is {n}")
        p_a = sand.p_a
        G_r, K_r = sand.elastic.G(p_a, e_ref), sand.elastic.K(p_a, e_ref)
        return dict(n=n, p_a=p_a, G_r=G_r, K_r=K_r, g=G_r / p_a, k=K_r / p_a)

    def params(self, sand: Sand, p_init: float, e_init: float, policy: ElasticPolicy) -> dict:
        c = self.constants(sand, policy.e_ref(e_init))
        return {self.O2_FIELDS["n"]: c["n"], self.O2_FIELDS["G_r"]: c["g"], self.O2_FIELDS["K_r"]: c["k"]}

    def oracle_extra(self, oracle: str) -> dict:
        ok, why = self.available(oracle)
        if not ok:
            raise EnergyUnavailable(f"HAR: {why}")
        return {"energy": "HAR"}

    def describe(self, sand: Sand, p_init: float, e_init: float, policy: ElasticPolicy) -> dict:
        e_ref = policy.e_ref(e_init)
        c = self.constants(sand, e_ref)
        nu_axis = (3.0 * c["k"] - 2.0 * c["g"]) / (6.0 * c["k"] + 2.0 * c["g"])
        return dict(energy=self.name, e_ref=e_ref, nu_axis=nu_axis, nu_target=sand.elastic.nu,
                    note="G, K follow the targets' p^n at every p; p_a is the shared -p_a flag", **c)

    def tcl_flags(self, sand: Sand, policy: ElasticPolicy, e_ref: float | None = None) -> list:
        """-energy HAR -k -g -n (sheet §2.4 names); -p_a is the CSL's flag (the same number)."""
        if policy.mode != "global":
            raise ValueError("tcl_flags: a deck material needs ElasticPolicy.bvp(e_rep=...)")
        c = self.constants(sand, policy.e_ref(e_ref))
        return [["-energy", "HAR"], ["-k", c["k"]], ["-g", c["g"]], ["-n", c["n"]]]


ENERGIES = {"BA06": BA06(), "HAR": HAR()}


def get(name: str):
    try:
        return ENERGIES[name]
    except KeyError:
        raise KeyError(f"unknown energy {name!r}; registered: {sorted(ENERGIES)}") from None
