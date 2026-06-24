"""Probe v3: Hertz via a FIXED rigid spherical indenter + the ladrunoContactForce query.

The textbook Hertz problem (rigid sphere indenting an elastic half-space), set up to
CONVERGE robustly and read pressures EXACTLY:
  * half-space block = MASTER (flat top facets), bottom fixed, sides roller.
  * rigid sphere = slave nodes FIXED at the lowered sphere surface (NO free body ⇒ no
    penalty-contact bootstrap singularity; the block is the only deformable body).
  * stiff penalty kn ⇒ ~rigid contact; each slave's EXACT normal force tn=kn·<−gap>₊ is
    read by the new `ladrunoContactForce` query (the adapter traction is not in nodeReaction).
Hertz (approach δ): a=sqrt(R δ), p0=(2 E*/π)(a/R), p(r)=p0 sqrt(1−(r/a)²), E*=E/(1−ν²),
P=(4/3)E* sqrt(R) δ^1.5.
"""
import os, sys, math
D = os.path.join(os.path.dirname(__file__), "..", "..", "dist", "bin")
sys.path.insert(0, os.path.abspath(D))
import opensees as ops


def hertz(R, delta, E, nu):
    Estar = E / (1 - nu * nu)
    a = math.sqrt(R * delta)
    p0 = (2.0 * Estar / math.pi) * (a / R)
    P = (4.0 / 3.0) * Estar * math.sqrt(R) * delta ** 1.5
    return Estar, a, p0, P


def run(nx=30, nz=10, W=3.0, H=1.5, R=20.0, delta=0.008, E=1000.0, nu=0.3,
        kn=1.0e6, ns=37, geom=False, nstep=4, verbose=True, tol=1e-8):
    Estar, aH, p0H, PH = hertz(R, delta, E, nu)
    ops.wipe(); ops.model("basic", "-ndm", 3, "-ndf", 3)
    dx = W / nx; dz = H / nz
    nid = {}; t = 1
    for k in range(nz + 1):
        for j in range(nx + 1):
            for i in range(nx + 1):
                ops.node(t, -W/2 + i*dx, -W/2 + j*dx, -H + k*dz); nid[(i, j, k)] = t; t += 1
    for k in range(nz + 1):
        for j in range(nx + 1):
            for i in range(nx + 1):
                if k == 0:
                    ops.fix(nid[(i, j, 0)], 1, 1, 1)
                elif i in (0, nx) or j in (0, nx):
                    ops.fix(nid[(i, j, k)], 1, 1, 0)
    ops.nDMaterial("ElasticIsotropic", 1, E, nu)
    eid = 1
    for k in range(nz):
        for j in range(nx):
            for i in range(nx):
                ops.element("LadrunoBrick", eid,
                            nid[(i,j,k)], nid[(i+1,j,k)], nid[(i+1,j+1,k)], nid[(i,j+1,k)],
                            nid[(i,j,k+1)], nid[(i+1,j,k+1)], nid[(i+1,j+1,k+1)], nid[(i,j+1,k+1)], 1)
                eid += 1
    mtags = []
    for j in range(nx):
        for i in range(nx):
            mtags += [nid[(i,j,nz)], nid[(i+1,j,nz)], nid[(i+1,j+1,nz)], nid[(i,j+1,nz)]]
    ops.contactSurface(10, "-master", 4, *mtags)
    span = 1.8 * aH
    sdx = 2 * span / (ns - 1)
    slaves = []; scoord = {}; sid = 500000
    for jj in range(ns):
        for ii in range(ns):
            x = -span + ii*sdx; y = -span + jj*sdx; r = math.hypot(x, y)
            if r > span:
                continue
            zsl = (R - math.sqrt(R*R - r*r)) - delta
            ops.node(sid, x, y, zsl); ops.fix(sid, 1, 1, 1)
            slaves.append(sid); scoord[sid] = (r, sdx*sdx); sid += 1
    ops.contactSurface(20, "-slave", *slaves)
    args = [1, 10, 20, kn, 0.0, 0.0, "-outward", 0.0, 0.0, 1.0]
    if geom:
        args.append("-geomtan")
    ops.contact(*args)
    ops.constraints("LadrunoContact"); ops.numberer("Plain"); ops.system("ProfileSPD")
    ops.integrator("LoadControl", 1.0 / nstep)
    ops.test("NormDispIncr", tol, 100, 0)
    ops.algorithm("Newton"); ops.analysis("Static")
    rc = 0; itlast = 0
    for _ in range(nstep):
        rc = ops.analyze(1); itlast = ops.testIter()
        if rc != 0:
            break
    data = []
    for s in slaves:
        r, area = scoord[s]
        f = ops.ladrunoContactForce(s)
        if f > 1e-9:
            data.append((r, f / area, f))
    data.sort()
    a_fe = max((r for r, p, f in data), default=0.0)
    Psum = sum(f for r, p, f in data)
    p0_fe = max((p for r, p, f in data), default=0.0)
    if verbose:
        print(f"geom={geom} rc={rc} it={itlast}  nC={len(data)}")
        print(f"   a:  FE={a_fe:.4f} Hertz={aH:.4f}  ratio={a_fe/aH:.3f}")
        print(f"   p0: FE={p0_fe:.2f} Hertz={p0H:.2f}  ratio={p0_fe/p0H:.3f}")
        print(f"   P:  FE={Psum:.3f} Hertz={PH:.3f}  ratio={Psum/PH:.3f}")
        # elliptic shape check at r≈0.5a
        if a_fe > 0:
            mid = [p for r, p, f in data if 0.4*a_fe < r < 0.6*a_fe]
            if mid:
                print(f"   p(0.5a)/p0: FE={sum(mid)/len(mid)/p0_fe:.3f}  Hertz={math.sqrt(1-0.25):.3f}")
    return dict(rc=rc, it=itlast, a=a_fe, p0=p0_fe, P=Psum, aH=aH, p0H=p0H, PH=PH, data=data)


if __name__ == "__main__":
    print("=== coarse ==="); run(nx=24, nz=8, ns=31)
    print("=== fine   ==="); run(nx=36, nz=12, ns=45)
