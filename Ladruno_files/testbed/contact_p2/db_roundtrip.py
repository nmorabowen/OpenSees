"""ADR-78 P2b gate: contact definitions survive Domain::sendSelf/recvSelf.

Serial, through the `database File` path (FE_Datastore IS a Channel, and
`restore` after a wipe exercises Domain::recvSelf's geometry-REBUILD branch,
which must reconstruct the contact engine from the stream).

Three cases:
  ref       -- uninterrupted build + analyze: the reference w15 WITH contact.
  nocontact -- same model, contact never declared: the discriminator value a
               silently-dropped definitions block would reproduce (top block
               falls through the gap; |w15| far larger).
  roundtrip -- build, save(1), WIPE, re-open the database, restore(1),
               fresh analysis objects, analyze. Must equal `ref` exactly.
  liveverify-- build, save(1), restore(1) with NO wipe (the recv-again branch):
               the engine is populated, so recvSelf takes the VERIFY path.
               Must not warn, must still analyze to `ref`.

Usage:  python db_roundtrip.py <dist-bin-dir>
Prints one `P2DB <json>` line; the run.py driver asserts on it.
"""
import json
import os
import sys
import tempfile


def build(ops, with_contact=True):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    coords = {
        1: (0, 0, 0), 2: (1, 0, 0), 3: (1, 1, 0), 4: (0, 1, 0),
        5: (0, 0, 1), 6: (1, 0, 1), 7: (1, 1, 1), 8: (0, 1, 1),
        11: (0, 0, 0.999), 12: (1, 0, 0.999), 13: (1, 1, 0.999), 14: (0, 1, 0.999),
        15: (0, 0, 1.999), 16: (1, 0, 1.999), 17: (1, 1, 1.999), 18: (0, 1, 1.999),
    }
    for t, (x, y, z) in coords.items():
        ops.node(t, float(x), float(y), float(z))
    ops.nDMaterial("ElasticIsotropic", 1, 2.0e7, 0.0)
    ops.element("stdBrick", 101, 1, 2, 3, 4, 5, 6, 7, 8, 1)
    ops.element("stdBrick", 102, 11, 12, 13, 14, 15, 16, 17, 18, 1)
    for n in (1, 2, 3, 4):
        ops.fix(n, 1, 1, 1)
    for n in (5, 6, 7, 8, 11, 12, 13, 14, 15, 16, 17, 18):
        ops.fix(n, 1, 1, 0)
    if with_contact:
        ops.contactSurface(1, "-master", 4, 5, 6, 7, 8)
        ops.contactSurface(2, "-slave", 11, 12, 13, 14)
        ops.contact(1, 1, 2, "auto", "-outward", 0.0, 0.0, 1.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n in (15, 16, 17, 18):
        ops.load(n, 0.0, 0.0, -2500.0)


def analyse(ops):
    ops.constraints("LadrunoContact")
    ops.numberer("RCM")
    ops.system("UmfPack")
    ops.test("NormDispIncr", 1.0e-10, 50, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    return ops.analyze(1)


def main():
    distbin = os.path.abspath(sys.argv[1])
    sys.path.insert(0, distbin)
    import opensees as ops
    got = os.path.normcase(os.path.abspath(ops.__file__))
    if not got.startswith(os.path.normcase(distbin)):
        raise RuntimeError("opensees resolved to %s, expected under %s" % (ops.__file__, distbin))

    out = {"build": ops.ladrunoBuild(), "pyd": ops.__file__}

    # reference (uninterrupted, with contact)
    build(ops, with_contact=True)
    out["ref_ok"] = analyse(ops)
    out["w15_ref"] = ops.nodeDisp(15, 3)

    # discriminator: what a dropped definitions block would look like
    build(ops, with_contact=False)
    out["nc_ok"] = analyse(ops)
    out["w15_nocontact"] = ops.nodeDisp(15, 3)

    dbdir = tempfile.mkdtemp(prefix="p2db_")
    dbpath = os.path.join(dbdir, "p2")

    # roundtrip: save -> wipe -> restore (geometry-rebuild branch) -> analyze
    build(ops, with_contact=True)
    ops.database("File", dbpath)
    ops.save(1)
    ops.wipe()
    ops.database("File", dbpath)
    ops.restore(1)
    out["rt_ok"] = analyse(ops)
    out["w15_roundtrip"] = ops.nodeDisp(15, 3)

    # vanilla control: a CONTACT-FREE model through the same save -> wipe ->
    # restore. The 17->19 domainData extension rides EVERY model's stream, so
    # the vanilla path must round-trip too (contactSize=0, no Vector sent).
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, (x, y, z) in {1: (0, 0, 0), 2: (1, 0, 0), 3: (1, 1, 0), 4: (0, 1, 0),
                         5: (0, 0, 1), 6: (1, 0, 1), 7: (1, 1, 1), 8: (0, 1, 1)}.items():
        ops.node(t, float(x), float(y), float(z))
    ops.nDMaterial("ElasticIsotropic", 1, 2.0e7, 0.0)
    ops.element("stdBrick", 101, 1, 2, 3, 4, 5, 6, 7, 8, 1)
    for n in (1, 2, 3, 4):
        ops.fix(n, 1, 1, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n in (5, 6, 7, 8):
        ops.load(n, 0.0, 0.0, -2500.0)
    out["van_ref_ok"] = analyse(ops)
    out["w5_van_ref"] = ops.nodeDisp(5, 3)
    # rebuild, save, wipe, restore, analyze
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, (x, y, z) in {1: (0, 0, 0), 2: (1, 0, 0), 3: (1, 1, 0), 4: (0, 1, 0),
                         5: (0, 0, 1), 6: (1, 0, 1), 7: (1, 1, 1), 8: (0, 1, 1)}.items():
        ops.node(t, float(x), float(y), float(z))
    ops.nDMaterial("ElasticIsotropic", 1, 2.0e7, 0.0)
    ops.element("stdBrick", 101, 1, 2, 3, 4, 5, 6, 7, 8, 1)
    for n in (1, 2, 3, 4):
        ops.fix(n, 1, 1, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n in (5, 6, 7, 8):
        ops.load(n, 0.0, 0.0, -2500.0)
    ops.database("File", dbpath + "_van")
    ops.save(1)
    ops.wipe()
    ops.database("File", dbpath + "_van")
    ops.restore(1)
    out["van_rt_ok"] = analyse(ops)
    out["w5_van_rt"] = ops.nodeDisp(5, 3)

    # live verify: restore onto the still-populated model (recv-again branch)
    build(ops, with_contact=True)
    ops.database("File", dbpath + "_live")
    ops.save(1)
    ops.restore(1)
    out["lv_ok"] = analyse(ops)
    out["w15_liveverify"] = ops.nodeDisp(15, 3)

    ops.wipe()
    print("P2DB " + json.dumps(out))


if __name__ == "__main__":
    main()
