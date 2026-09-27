"""WP-128 side note: the two pre-existing failures in
tests/test_ladruno_sanisand_flip_determinism.py.  Runs the test's own child
deck (default -flipAlphaIn, 2 holds, 10 push steps) in a subprocess at
MKL_NUM_THREADS = 1 and prints the step list and the tail of stderr, then the
same with the push under `system FullGeneral` instead of Pardiso (is the -3 the
linear solver or the material?)."""
import json, os, subprocess, sys
import _boot as B
sys.path.insert(0, os.path.join(B.ROOT, "tests"))
import test_ladruno_sanisand_flip_determinism as T

def child(src, threads=1):
    env = dict(os.environ)
    env["MKL_NUM_THREADS"] = env["OMP_NUM_THREADS"] = str(threads)
    p = subprocess.run([sys.executable, "-S", "-u", "-c",
                        "import os,sys; os.add_dll_directory(%r); sys.path.insert(0,%r); sys.path.append(%r)\n" % (B.BIN, B.BIN, B.SITE) + src,
                        os.path.join(B.ROOT, "tests"), json.dumps([["default", 2]]), "10", "1"],
                       env=env, stdin=subprocess.DEVNULL, capture_output=True, text=True, timeout=900)
    out = [l for l in p.stdout.splitlines() if l.startswith("RESULT ")]
    return (json.loads(out[-1][7:]) if out else None), p.stderr

for label, src in (("Pardiso", T._CHILD),
                   ("FullGeneral", T._CHILD.replace("ops.system('Pardiso', *(('-stats',) if stats else ()))", "ops.system('FullGeneral')")),
                   ("Pardiso,NormDispIncr1e-6", T._CHILD.replace("ops.test('NormDispIncr', 1.0e-8, 100, 0)", "ops.test('NormDispIncr', 1.0e-6, 100, 0)"))):
    res, err = child(src)
    print("==", label)
    if res:
        r = res["runs"][0]
        print("   build", res["build"], "roundoff", r["roundoff"], "/", r["census"])
        print("   steps", [(s[0], round(s[2], 6)) for s in r["steps"]])
    tail = [l for l in err.splitlines() if "WARNING" in l or "failed" in l.lower() or "Error" in l]
    print("   stderr (WARNING/failed lines, first 12):")
    for l in tail[:12]:
        print("     ", l[:300])
