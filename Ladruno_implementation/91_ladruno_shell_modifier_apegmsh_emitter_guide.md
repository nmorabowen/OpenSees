---
title: LadrunoShellModifier — implementation guide for apeGmsh
project: Ladruno
status: ready to implement — fork dependency shipped in PR #793
priority: high
adr: ADR 91
tags:
  - apegmsh
  - emitter
  - section
  - shell
  - etabs
  - cracked-section
---

# LadrunoShellModifier — implementation guide for apeGmsh

**What we are asking for:** a `LadrunoShellModifier` section primitive in
`apeGmsh.opensees.section.plate`, emitted to both the Tcl and openseespy backends.

That part is small — it is a decorator dataclass with one dependency and a nine-flag emit. Most
of this note is about the **two things around it that are wrong today**, one of which is
silently producing incorrect models right now.

Background: [[91_ladruno_shell_stiffness_modifiers_adr]] (spec),
[[LadrunoShellModifier_guide]] (consumer-language guide).

---

## 1. The headline: `etabs_import` is silently dropping cracked-section stiffness

`src/apeGmsh/interop/etabs_import.py` (571 lines) contains **no** handling of area section
property modifiers. Grep it for `modifier`, `f11`, `f22`, `m11` — nothing.

That means: when someone imports a real ETABS building — where the shear walls are almost
certainly assigned `f11 = f22 = f12 = 0.35` per ACI 318-25 §6.6.3.1.1, because that is what
every practising engineer does — apeGmsh builds those walls at **gross** stiffness. The walls
come out roughly 3× too stiff in plane. Nothing warns. The model runs, converges, and produces
confidently wrong drifts, wrong periods, and a wrong distribution of base shear between walls
and frames.

This is not a missing feature so much as an active correctness bug in the import path, and it
predates WP-91. It is worth fixing even if the emitter work below is deferred.

**Minimum acceptable fix, in priority order:**

1. **Refuse loudly.** If the imported model carries any area modifier not equal to 1.0 and the
   importer cannot represent it, raise — do not warn-and-continue. A silently-3×-stiff wall is
   worse than a failed import.
2. **Then** represent it, per §2–§4 below.

Whoever picks this up should check whether `apeETABS`'s reader surfaces the modifiers at all —
`cAreaObj.GetModifiers` / `cPropArea.GetModifiers`, a 10-entry array — since `etabs_import`
consumes a `StructuralModel` and the data may be getting lost upstream of apeGmsh rather than
inside it.

## 2. The primitive

Follows the existing `Section` shape in `src/apeGmsh/opensees/section/plate.py` exactly. It is a
decorator, so unlike `ElasticMembranePlateSection` it has a real `dependencies()`:

```python
@dataclass(frozen=True)
class LadrunoShellModifier(Section):
    """ETABS-style stiffness modifiers on any order-8 plate section (ADR 91)."""

    inner: Section
    f11: float = 1.0
    f22: float = 1.0
    f12: float = 1.0
    m11: float = 1.0
    m22: float = 1.0
    m12: float = 1.0
    v13: float = 1.0
    v23: float = 1.0
    mass: float = 1.0

    def __post_init__(self) -> None:
        for name in ("f11", "f22", "f12", "m11", "m22", "m12", "v13", "v23", "mass"):
            v = getattr(self, name)
            if v < 0.0:
                raise ValueError(f"{name} must be >= 0.0, got {v}")

    def _emit(self, emitter: "Emitter", tag: int) -> None:
        inner_tag = resolve_section_tag(emitter, self.inner)
        args: list = []
        for name in ("f11", "f22", "f12", "m11", "m22", "m12", "v13", "v23", "mass"):
            v = getattr(self, name)
            if v != 1.0:                      # emit only what differs from default
                args += [f"-{name}", v]
        emitter.section("LadrunoShellModifier", tag, inner_tag, *args)

    def dependencies(self) -> tuple[Primitive, ...]:
        return (self.inner,)
```

Tag resolution for `inner` should reuse whatever `LayeredShell` already does for its layer
materials — `plate.py:179` uses a closure-captured `resolve_mat_tag` attached to the emitter,
and `plate.py:165-168` flags this as an open coordinator question. Use the same mechanism, do
not invent a second one.

**Emit only non-default flags.** Every flag defaults to `1.0` on the fork side and an
all-defaults wrap is byte-identical to the inner section, so a sparse emit is both shorter and
self-documenting in the deck. It also means an all-1.0 wrapper is a harmless no-op, which lets
a generator wrap unconditionally.

## 3. The trap: area bucketing is keyed on section name alone

`etabs_import.py:189-198` buckets areas by section:

```python
area_buckets.setdefault(ar.section, []).append(surf)
```

ETABS modifiers can be assigned **per area object**, overriding the section property. Two walls
sharing section `"W30"` can legitimately carry different modifiers — that is a normal thing to
do (a coupling-beam-adjacent pier cracked harder than its neighbour).

So the moment modifiers enter the picture, `ar.section` is **no longer a sufficient bucket
key**. It must become something like `(ar.section, modifier_tuple)`, with one emitted
`LadrunoShellModifier` section per distinct combination and a stable, human-readable PG name.

Get this wrong and the failure is silent and severe: every wall in the group inherits whichever
area happened to be first, so a subset of walls gets another wall's cracking. This is the same
class of bug as §1 and will not show up as an error.

## 4. The `Ep_mod` gap, while you are in there

`plate.py:94-99` emits `ElasticMembranePlateSection` as `E, nu, h, rho` — four arguments. The
upstream section takes a **fifth**, `Ep_mod`, an out-of-plane stiffness modifier (upstream
contribution, Degenkolb). apeGmsh cannot currently express it at all.

It is optional and defaults to `1.0`, so nothing is broken today. But once `LadrunoShellModifier`
exists there are two overlapping ways to express an out-of-plane reduction, and users will find
both. Recommendation: **add the optional `Ep_mod` field for round-trip fidelity** (so a deck
read from elsewhere survives), but have the ETABS import path and any new authoring API route
through `LadrunoShellModifier` instead, which is strictly more expressive. Document that
`Ep_mod=r` is exactly equivalent to `m11=m22=m12=v13=v23=r` — the fork's G5 gate pins that
equivalence to round-off, so it is a fact you can rely on rather than an approximation.

## 5. Mapping ETABS → the emitter

The ETABS OAPI area modifier array is 10 entries, in this order:

```
[f11, f22, f12, m11, m22, m12, v13, v23, mass, weight]
```

The first nine map straight through. **The tenth does not exist on the fork side** and must not
be silently discarded:

In OpenSees the shell self-weight body force is derived from the same `getRho()` that builds the
mass matrix (`ShellMITC4.cpp:1725-1755`), so a weight modifier could only alias the mass
modifier. ADR 91 §5 declines to ship an argument that would quietly do something other than its
name.

So the importer must **refuse, or warn explicitly, on `weight != 1.0`** and tell the user to
scale self-weight at the load level instead (the `eleLoad -type -selfWeight` factors, or the
gravity pattern factor). Dropping it silently reintroduces exactly the class of bug in §1.

## 6. One modelling note for whoever writes the docs

The modifiers are applied as a congruence, `D' = S·D·S` with `S = diag(√f11 … √v23)`, so the
Poisson coupling term moves as `√(f11·f22)` rather than being left alone. When `f11 == f22` —
every standard cracked-wall recipe — this is indistinguishable from scaling the whole membrane
block, and the distinction never surfaces. It only matters for strongly unequal `f11`/`f22`.

The reason it is not a diagonal-only rescale: that destroys positive definiteness at exactly the
cracked-wall modifiers this feature exists for (ADR 91 §4). Whether ETABS itself uses the same
convention is genuinely unknown — CSI does not document it — and is tracked as ADR 91 OQ-1. If
apeGmsh gains an ETABS round-trip test with unequal `f11`/`f22`, **that test would settle the
open question**, and the result is worth reporting back to the fork.

## 7. Acceptance

- [ ] `LadrunoShellModifier` primitive, both backends, sparse-flag emit
- [ ] All-defaults wrap emits and round-trips as a no-op
- [ ] `etabs_import` either represents area modifiers or refuses loudly (§1)
- [ ] Area bucket key includes the modifier tuple (§3)
- [ ] `weight != 1.0` refused with an actionable message (§5)
- [ ] Optional `Ep_mod` on `ElasticMembranePlateSection` for round-trip fidelity (§4)

Fork dependency: PR #793 (draft at time of writing). The verb is frozen —
`section LadrunoShellModifier $tag $innerSecTag [-f11 v] ... [-mass v]` — and covered by 15
gates plus an ADR-87 D2 mutation gate, so it is safe to build against before merge.
