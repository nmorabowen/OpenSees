"""Parameterized 3D elastic solid box for FEAST L2 profiling.

A structured hex (stdBrick) grid, fixed at the z=0 base, elastic isotropic with
mass density -> a solve-dominated 3D generalized eigenproblem K phi = w^2 M phi
(3D fill-in makes the distributed factorization the cost centre L3 parallelizes
and L2 would multiply). Replicated: every openseesmp rank builds the FULL model.

Size knob: `build_box(ops, ne)` makes an ne x ne x ne element grid.
  DOF ~ 3*((ne+1)^3 - (ne+1)^2)   [base nodes fully fixed]
    ne=20 ->  ~26.5k    ne=30 ->  ~86.4k    ne=40 ->  ~206k
    ne=55 ->  ~527k     ne=70 ->  ~1.06M
"""
import math

E, NU, RHO = 30.0e9, 0.2, 2400.0        # concrete-ish elastic solid
LX = LY = LZ = 10.0                      # box dimensions (m); aspect 1:1:1
TWO_PI = 2.0 * math.pi


def build_box(ops, ne):
    """Build the ne^3 stdBrick box on the active openseesmp interpreter.

    Returns (ndof_free, nnode, nele)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("ElasticIsotropic", 1, E, NU, RHO)

    nn = ne + 1
    dx, dy, dz = LX / ne, LY / ne, LZ / ne

    def nid(i, j, k):
        return 1 + i + j * nn + k * nn * nn

    for k in range(nn):
        z = k * dz
        for j in range(nn):
            y = j * dy
            for i in range(nn):
                t = nid(i, j, k)
                ops.node(t, i * dx, y, z)
                if k == 0:
                    ops.fix(t, 1, 1, 1)

    ele = 0
    for k in range(ne):
        for j in range(ne):
            for i in range(ne):
                ele += 1
                n1 = nid(i, j, k);     n2 = nid(i + 1, j, k)
                n3 = nid(i + 1, j + 1, k);   n4 = nid(i, j + 1, k)
                n5 = nid(i, j, k + 1); n6 = nid(i + 1, j, k + 1)
                n7 = nid(i + 1, j + 1, k + 1); n8 = nid(i, j + 1, k + 1)
                ops.element("stdBrick", ele, n1, n2, n3, n4, n5, n6, n7, n8, 1)

    nnode = nn * nn * nn
    ndof_free = 3 * (nnode - nn * nn)      # base layer fully fixed
    return ndof_free, nnode, ele


def hz(lam):
    return math.sqrt(abs(lam)) / TWO_PI


def calibrate_fmax(ops, ne, want_modes):
    """Fixed, deterministic band [0, fmax] enclosing ~want_modes low modes, via a
    one-shot ARPACK eigen (same on every replicated rank). Returns (fmax, lambdas)."""
    build_box(ops, ne)
    lam = list(ops.eigen(want_modes + 2))
    # midpoint between the want_modes-th and the next -> a clean band edge
    fmax = 0.5 * (hz(lam[want_modes - 1]) + hz(lam[want_modes]))
    return fmax, lam
