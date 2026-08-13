"""ADR-78 P2b review follow-up: the ALL-LANES definitions round trip.

The shipped db_roundtrip.py proves save -> wipe -> restore for ONE lane
(a single NTS `-kn auto` contact). packDefinitions/unpackDefinitions also
carry the mortar lane (37 slots, friction + tie + edge-edge fields), the
rigid-plane lane (12 slots) and NTS option fields (kt/mu/cell/visc/...),
none of which had ever EXECUTED a round trip. This closes that gap with
two instruments:

1. BEHAVIORAL PARITY, one column pair per lane in one model:
     A  NTS contact, explicit kn + kt/mu + -cell + -geomtan
     B  mortar -tie (-epsTie), loaded in TENSION (tie is load-bearing;
        a contact could not masquerade as it)
     C  contactPlane + -visc (statics-inert, field packed)
     D  mortar friction (-epsN/-mu/-epsT/-cohesion/-tauMax) + the
        -edgeedge option block (inert on parallel facets; fields packed)
   build -> analyze = reference; build -> save -> wipe -> restore ->
   analyze must reproduce every lane's tip EXACTLY.

2. FIELD SENSITIVITY, byte-level: the packed definitions Vector is the
   distinctive large .VECs record FileDatastore writes. Saving the same
   model twice must produce IDENTICAL defs bytes (determinism control);
   saving a variant with ONE field changed must produce DIFFERENT bytes
   (the field is really packed). One variant per lane's slot block:
   kt (NTS), mu (mortar friction), edgeKn (edge block), epsTie (tie),
   kn + normal (plane). A field the pack silently dropped would fail its
   sensitivity case here without needing the physics to depend on it.

Plus: the live-restore VERIFY path stays silent on a matching model and
WARNS when the live engine's definitions no longer match the stream
(declaring one more contact after the save changes the packed size).

Usage:  python db_roundtrip_all_lanes.py <dist-bin-dir>
Prints one `P2DBALL <json>` line; run.py asserts on it.
"""
import glob
import json
import os
import sys
import tempfile

PAIR_DX = 10.0     # x-offset between the four column pairs
GAP = 1e-3


def _pair_nodes(ops, off, x0):
    """Two stacked unit bricks, P0 geometry, tags offset by `off`."""
    lo = {1: (0, 0, 0), 2: (1, 0, 0), 3: (1, 1, 0), 4: (0, 1, 0),
          5: (0, 0, 1), 6: (1, 0, 1), 7: (1, 1, 1), 8: (0, 1, 1)}
    hi = {11: (0, 0, 1 - GAP), 12: (1, 0, 1 - GAP), 13: (1, 1, 1 - GAP),
          14: (0, 1, 1 - GAP), 15: (0, 0, 2 - GAP), 16: (1, 0, 2 - GAP),
          17: (1, 1, 2 - GAP), 18: (0, 1, 2 - GAP)}
    for t, (x, y, z) in {**lo, **hi}.items():
        ops.node(t + off, float(x) + x0, float(y), float(z))


def build(ops, variant=None):
    """The four-lane model. `variant` changes exactly one packed field."""
    v = variant or {}
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("ElasticIsotropic", 1, 2.0e7, 0.0)

    # --- pair A (off 0): NTS contact, explicit kn, friction fields ---
    _pair_nodes(ops, 0, 0.0)
    ops.element("stdBrick", 101, 1, 2, 3, 4, 5, 6, 7, 8, 1)
    ops.element("stdBrick", 102, 11, 12, 13, 14, 15, 16, 17, 18, 1)
    for n in (1, 2, 3, 4):
        ops.fix(n, 1, 1, 1)
    for n in (5, 6, 7, 8, 11, 12, 13, 14, 15, 16, 17, 18):
        ops.fix(n, 1, 1, 0)
    ops.contactSurface(1, "-master", 4, 5, 6, 7, 8)
    ops.contactSurface(2, "-slave", 11, 12, 13, 14)
    ops.contact(1, 1, 2, 2.0e7, v.get("ktA", 1.0e6), 0.3,
                "-outward", 0.0, 0.0, 1.0, "-cell", 1.5, "-geomtan")

    # --- pair B (off 100): mortar mesh-tie, loaded in TENSION ---
    _pair_nodes(ops, 100, PAIR_DX)
    ops.element("stdBrick", 111, 101, 102, 103, 104, 105, 106, 107, 108, 1)
    ops.element("stdBrick", 112, 111, 112, 113, 114, 115, 116, 117, 118, 1)
    for n in (101, 102, 103, 104):
        ops.fix(n, 1, 1, 1)
    for n in (105, 106, 107, 108, 111, 112, 113, 114, 115, 116, 117, 118):
        ops.fix(n, 1, 1, 0)
    ops.contactSurface(3, "-master", 4, 105, 106, 107, 108)
    ops.contactSurface(4, "-slave-segments", 4, 111, 112, 113, 114)
    ops.contact(2, 3, 4, "-mortar", "-tie",
                "-epsTie", v.get("epsTie", 1.0e7),
                "-outward", 0.0, 0.0, 1.0)

    # --- pair C (off 200): single brick on a rigid plane, -visc packed ---
    # The declared normal is deliberately TILTED (norm sqrt(1+49) -- nowhere
    # near a power of two), so the stored unit normal has bits that the F2
    # renormalize-on-unpack defect would flip. An axis-aligned (0,0,1) normal
    # is exact under sqrt/divide and can never see F2 (why v1 missed it).
    hi = {211: (0, 0, GAP), 212: (1, 0, GAP), 213: (1, 1, GAP),
          214: (0, 1, GAP), 215: (0, 0, 1 + GAP), 216: (1, 0, 1 + GAP),
          217: (1, 1, 1 + GAP), 218: (0, 1, 1 + GAP)}
    for t, (x, y, z) in hi.items():
        ops.node(t, float(x) + 2 * PAIR_DX, float(y), float(z))
    ops.element("stdBrick", 122, 211, 212, 213, 214, 215, 216, 217, 218, 1)
    for n in hi:
        ops.fix(n, 1, 1, 0)
    ops.contactSurface(5, "-slave", 211, 212, 213, 214)
    ops.contactPlane(3, 5,
                     v.get("planeNx", 1.0), 0.0, v.get("planeNz", 7.0),
                     2 * PAIR_DX, 0.0, 0.0,
                     v.get("planeKn", 2.0e7), "-visc", 5.0)

    # --- pair D (off 300): mortar friction + the -edgeedge option block ---
    _pair_nodes(ops, 300, 3 * PAIR_DX)
    ops.element("stdBrick", 131, 301, 302, 303, 304, 305, 306, 307, 308, 1)
    ops.element("stdBrick", 132, 311, 312, 313, 314, 315, 316, 317, 318, 1)
    for n in (301, 302, 303, 304):
        ops.fix(n, 1, 1, 1)
    for n in (305, 306, 307, 308, 311, 312, 313, 314, 315, 316, 317, 318):
        ops.fix(n, 1, 1, 0)
    ops.contactSurface(6, "-master", 4, 305, 306, 307, 308)
    ops.contactSurface(7, "-slave-segments", 4, 311, 312, 313, 314)
    ops.contact(4, 6, 7, "-mortar", "-epsN", 2.0e7,
                "-mu", v.get("muD", 0.3), "-epsT", 1.0e6,
                "-cohesion", 10.0, "-tauMax", 1.0e5,
                "-outward", 0.0, 0.0, 1.0,
                "-edgeedge", "-edgeKn", v.get("edgeKn", 5.0e5),
                "-edgeBand", 0.01, "-edgeMu", 0.2, "-edgeKt", 1.0e5)

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n in (15, 16, 17, 18):                     # A: compression
        ops.load(n, 0.0, 0.0, -2500.0)
    for n in (115, 116, 117, 118):                 # B: TENSION (tie carries it)
        ops.load(n, 0.0, 0.0, +1000.0)
    for n in (215, 216, 217, 218):                 # C: onto the plane
        ops.load(n, 0.0, 0.0, -2500.0)
    for n in (315, 316, 317, 318):                 # D: compression
        ops.load(n, 0.0, 0.0, -2500.0)


def analyse(ops):
    ops.constraints("LadrunoContact")
    ops.numberer("RCM")
    ops.system("UmfPack")
    ops.test("NormDispIncr", 1.0e-10, 100, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    return ops.analyze(1)


def tips(ops):
    return {"A": ops.nodeDisp(15, 3), "B": ops.nodeDisp(115, 3),
            "C": ops.nodeDisp(215, 3), "D": ops.nodeDisp(315, 3)}


def defs_bytes(dbdir):
    """The packed-definitions Vector record: the distinctive LARGE .VECs
    file (node-state vectors are ndf-sized; domainTime is 1)."""
    cands = []
    for f in glob.glob(os.path.join(dbdir, "*.VECs.*")):
        # FileDatastore names records `<db>.VECs.<size>.<commitTag>`
        parts = f.split(".")
        try:
            size = int(parts[-2])
        except (ValueError, IndexError):
            continue
        if size > 50:
            cands.append((size, f))
    if len(cands) != 1:
        raise RuntimeError("expected exactly one large .VECs record, got %r"
                           % cands)
    size, path = cands[0]
    with open(path, "rb") as fh:
        raw = fh.read()
    # The FileDatastore record = a small header (its dbTag, which depends on
    # how many databases this PROCESS opened before -- nondeterministic here)
    # + the vector payload. Compare only the trailing size*8 payload bytes.
    payload = raw[-size * 8:]
    if len(payload) != size * 8:
        raise RuntimeError("record shorter than its own payload?")
    return payload


def save_to(ops, dbdir, variant=None):
    build(ops, variant)
    ops.database("File", os.path.join(dbdir, "db"))
    ops.save(1)


def main():
    distbin = os.path.abspath(sys.argv[1])
    sys.path.insert(0, distbin)
    if hasattr(os, "add_dll_directory"):
        os.add_dll_directory(distbin)
    import opensees as ops
    got = os.path.normcase(os.path.abspath(ops.__file__))
    if not got.startswith(os.path.normcase(distbin)):
        raise RuntimeError("opensees resolved to %s, expected under %s"
                           % (ops.__file__, distbin))
    out = {"build": ops.ladrunoBuild()}

    root = tempfile.mkdtemp(prefix="p2dball_")

    # 1. reference
    build(ops)
    out["ref_ok"] = analyse(ops)
    ref = tips(ops)
    out["tips_ref"] = ref

    # 2. save -> wipe -> restore -> analyze (behavioral parity, all lanes)
    rt_dir = os.path.join(root, "rt"); os.makedirs(rt_dir)
    save_to(ops, rt_dir)
    ops.wipe()
    ops.database("File", os.path.join(rt_dir, "db"))
    ops.restore(1)
    # 2b. BIT-EXACT round trip (review F2): re-save the just-restored engine
    #     and byte-compare the packed payload against the original save. The
    #     tilted plane normal makes this probe live -- pre-fix the unpack
    #     renormalized the already-unit normal (1 +/- 1 ulp) and this compare
    #     failed on the plane slots.
    rs_dir = os.path.join(root, "resave"); os.makedirs(rs_dir)
    ops.database("File", os.path.join(rs_dir, "db"))
    ops.save(1)
    out["resave_bitexact"] = (defs_bytes(rs_dir) == defs_bytes(rt_dir))
    ops.database("File", os.path.join(rt_dir, "db"))
    out["rt_ok"] = analyse(ops)
    out["tips_rt"] = tips(ops)

    # 3. live-restore verify path: matching model -> must stay silent
    #    (run.py greps the captured output for the warning text); then
    #    declare ONE MORE definition and restore again -> size mismatch
    #    -> the verify instrument must WARN.
    lv_dir = os.path.join(root, "lv"); os.makedirs(lv_dir)
    save_to(ops, lv_dir)
    ops.restore(1)                       # populated + matching: silent
    print("P2DBALL-MARK pre-mutation-restore-done")
    ops.contactSurface(8, "-slave", 15, 16, 17, 18)
    ops.contactPlane(9, 8, 0.0, 0.0, 1.0, 0.0, 0.0, 5.0, 1.0e6)
    ops.restore(1)                       # live defs != stream: must warn

    # 4. field sensitivity: baseline bytes, determinism twin, one variant
    #    per lane slot-block. Different bytes == the field is packed.
    base_dir = os.path.join(root, "base"); os.makedirs(base_dir)
    save_to(ops, base_dir)
    base = defs_bytes(base_dir)
    twin_dir = os.path.join(root, "twin"); os.makedirs(twin_dir)
    save_to(ops, twin_dir)
    out["determinism_same"] = (defs_bytes(twin_dir) == base)
    sens = {}
    for name, variant in (
            ("ktA_nts", {"ktA": 2.0e6}),
            ("muD_mortar", {"muD": 0.35}),
            ("edgeKn_edge", {"edgeKn": 6.0e5}),
            ("epsTie_tie", {"epsTie": 2.0e7}),
            ("planeKn_plane", {"planeKn": 3.0e7}),
            ("planeN_plane", {"planeNx": 0.1}),   # re-normalized on add
    ):
        vd = os.path.join(root, "v_" + name); os.makedirs(vd)
        save_to(ops, vd, variant)
        sens[name] = (defs_bytes(vd) != base)
    out["sensitivity"] = sens

    # 5. F3 rollback probe: corrupt the stream's nSurf count (slot 1) so the
    #    unpack dies mid-record, then restore TWICE. Both attempts must fail
    #    LOUDLY. Pre-fix, the first failure left the added prefix in the
    #    engine, so the retry routed to the verify path and "succeeded" with
    #    a warning -- laundering a corrupt restore into a truncated model.
    import struct
    cr_dir = os.path.join(root, "corrupt"); os.makedirs(cr_dir)
    save_to(ops, cr_dir)
    f = glob.glob(os.path.join(cr_dir, "*.VECs.*"))
    f = [x for x in f if int(x.split(".")[-2]) > 50][0]
    with open(f, "r+b") as fh:
        raw = fh.read()
        size = int(f.split(".")[-2])
        hdr = len(raw) - size * 8
        fh.seek(hdr + 8)                       # slot 1 = nSurf
        fh.write(struct.pack("<d", 9999.0))
    ops.wipe()
    ops.database("File", os.path.join(cr_dir, "db"))
    print("P2DBALL-MARK corrupt-restores-begin")
    n_raised = 0
    for _ in range(2):                         # both attempts must fail LOUDLY
        try:
            ops.restore(1)
        except Exception:
            n_raised += 1                      # the loud path raises here
    out["corrupt_restores_raised"] = n_raised

    ops.wipe()
    print("P2DBALL " + json.dumps(out))


if __name__ == "__main__":
    main()
