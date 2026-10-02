"""The ENERGY plug (task item 2): the only place where the harness knows which elastic energy the oracle uses.

The rest of the harness (drivers, objective, fit, smoke) asks for an energy BY NAME and gets back the oracle
parameters of the elastic part as common names (tests/conftest.py names; model.py renames per oracle). The
energy choice (BA06 vs HAR, plan §2.5 convexity gate) is being decided in parallel, so nothing outside this
module may depend on it.

Contract of a plug (class with these members):
  name                           registry key
  available(oracle) -> (ok, why) does the oracle module ('O1'/'O2') implement this energy?
  params(sand, p_init, e_init, policy) -> dict
                                 the elastic parameters (common names) for one element test, from the sand's
                                 ElasticTargets (G(p, e), nu), the test's initial mean stress magnitude p_init
                                 (kPa > 0) and void ratio e_init, and the ElasticPolicy
  describe(sand, p_init, e_init, policy) -> dict
                                 the same numbers plus the stiffness they imply, for the report

Adding HAR when O2 gets it: fill HAR.O2_FIELDS with the O2 field names (and HAR.constants if the sheet's
§2.3 constants are defined differently from (n, p_r, G_r, K_r)); no other harness file changes.

Sheet references: §2.1 (interface), §2.2 (BA06: K = -p/kappa_hat, constant mu0 for alpha0 = 0), §2.3 (HAR
slot), §15 (BA06 refit rule: mu0 = G(p_rep, e), kappa_hat = p_rep / K(p_rep)).
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
    HAR ignores p_rep (its stiffness scales with p^n itself)."""
    mode: str = "per_test"
    factor: float = 1.0
    p_rep: float | None = None
    e_rep: float | None = None

    def rep(self, p_init: float, e_init: float):
        if self.mode == "per_test":
            return self.factor * p_init, e_init
        if self.mode == "global":
            if self.p_rep is None:
                raise ValueError("ElasticPolicy('global') needs p_rep")
            return float(self.p_rep), (e_init if self.e_rep is None else float(self.e_rep))
        raise ValueError(f"unknown ElasticPolicy mode {self.mode!r}")


def _oracle_params_cls(oracle: str):
    from . import model
    return model.oracle_module(oracle).Params


class BA06:
    """Borja & Andrade 2006 energy (S.4), alpha0 = 0: K = -p/kappa_hat, G = mu0 (constant).
    Matched at (p_rep, e_rep): mu0 = G(p_rep, e_rep), kappa_hat = p_rep / K(p_rep, e_rep) (sheet §15).
    p0 = -p_rep, eps_v0 = 0: the energy's reference point is the representative state (any other choice
    only shifts eps^e, sheet S.4)."""
    name = "BA06"

    def available(self, oracle: str):
        # Both oracles implement BA06 and only BA06 (O2 kernel.elastic, O1 model.energy) as of 2026-10-02.
        fields = {f.name for f in dataclasses.fields(_oracle_params_cls(oracle))}
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
        fields = {f.name for f in dataclasses.fields(_oracle_params_cls(oracle))}
        return {"energy": "BA06"} if "energy" in fields else {}

    def describe(self, sand: Sand, p_init: float, e_init: float, policy: ElasticPolicy) -> dict:
        p_rep, e_rep = policy.rep(p_init, e_init)
        d = self.params(sand, p_init, e_init, policy)
        d.update(energy=self.name, p_rep=p_rep, e_rep=e_rep, G_target_at_rep=sand.elastic.G(p_rep, e_rep),
                 K_target_at_rep=sand.elastic.K(p_rep, e_rep),
                 note="K ∝ p (exponent 1), G constant: matches the sqrt(p) targets only at p_rep")
        return d


class HAR:
    """Houlsby, Amorosi & Rojas (2005) slot (sheet §2.3; plan §2.5). Stiffness ∝ p^n: matched to the DM04
    targets at p_r = p_a with n = the targets' exponent (0.5), so G and K follow sqrt(p) without a p_rep.
    NOT available until an oracle gets the energy: available() checks the oracle's Params for an 'energy'
    field and for every field of O2_FIELDS; params() raises EnergyUnavailable otherwise."""
    name = "HAR"
    # HAR constant -> oracle Params field name. Bind to O2's names when O2 gets HAR (only this dict changes).
    O2_FIELDS = {"n": "har_n", "p_r": "har_p_r", "G_r": "har_G_r", "K_r": "har_K_r"}

    def available(self, oracle: str):
        fields = {f.name for f in dataclasses.fields(_oracle_params_cls(oracle))}
        if "energy" not in fields:
            return False, f"{oracle} has no energy switch (BA06 only, sheet §2.3 slot empty as of 2026-10-02)"
        missing = [v for v in self.O2_FIELDS.values() if v not in fields]
        if missing:
            return False, f"{oracle} has an energy switch but not the HAR fields {missing} (bind HAR.O2_FIELDS)"
        return True, "HAR fields present"

    def constants(self, sand: Sand, e_ref: float) -> dict:
        p_r = sand.elastic.p_a
        return dict(n=sand.elastic.n_exp, p_r=p_r, G_r=sand.elastic.G(p_r, e_ref), K_r=sand.elastic.K(p_r, e_ref))

    def params(self, sand: Sand, p_init: float, e_init: float, policy: ElasticPolicy) -> dict:
        e_ref = e_init if policy.e_rep is None else policy.e_rep
        c = self.constants(sand, e_ref)
        # p0 stays the reference pressure of the energy (sheet §2.3 "a reference state where p = p0").
        out = {self.O2_FIELDS[k]: v for k, v in c.items()}
        out.update(p0=-p_init, eps_v0=0.0)
        return out

    def oracle_extra(self, oracle: str) -> dict:
        ok, why = self.available(oracle)
        if not ok:
            raise EnergyUnavailable(f"HAR: {why}")
        return {"energy": "HAR"}

    def describe(self, sand: Sand, p_init: float, e_init: float, policy: ElasticPolicy) -> dict:
        e_ref = e_init if policy.e_rep is None else policy.e_rep
        d = dict(energy=self.name, **self.constants(sand, e_ref), e_ref=e_ref)
        return d


ENERGIES = {"BA06": BA06(), "HAR": HAR()}


def get(name: str):
    try:
        return ENERGIES[name]
    except KeyError:
        raise KeyError(f"unknown energy {name!r}; registered: {sorted(ENERGIES)}") from None
