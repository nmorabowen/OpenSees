"""WP-128 bootstrap: load THIS worktree's opensees.pyd (WP-127 build) and the
replay helper.  Run every script here with

    C:\\Users\\nmora\\AppData\\Local\\Python\\pythoncore-3.12-64\\python.exe -S <script>

(-S: a boot .pth otherwise preloads a stale pyd from another worktree).
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
BIN = os.path.join(ROOT, "dist", "bin")
SITE = r"C:\Users\nmora\AppData\Local\Python\pythoncore-3.12-64\Lib\site-packages"

os.add_dll_directory(BIN)
sys.path.insert(0, BIN)
sys.path.insert(0, os.path.join(ROOT, "Ladruno_scripts"))
sys.path.insert(0, os.path.join(ROOT, "tests"))
if SITE not in sys.path:
    sys.path.append(SITE)
os.environ["PYTHONPATH"] = os.pathsep.join([BIN, os.path.join(ROOT, "tests")])

import opensees as ops  # noqa: E402

assert os.path.normcase(os.path.dirname(os.path.abspath(ops.__file__))) == \
    os.path.normcase(BIN), ops.__file__

import sanisand_replay as sr  # noqa: E402

OUT = os.path.join(HERE, "out")
os.makedirs(OUT, exist_ok=True)

# ---------------------------------------------------------------------------
# campaign material (attachments README; nu 0.312885)
P = list(sr.CAMPAIGN_PARAMS)
(G0, NU, E_INIT, MC, C_, LAMC, E0, KSI, PATM, M_, H0, CH, NB, A0, ND, ZMAX,
 CZ, DEN) = P

# prototype tags
TAG_ME = 1        # IntScheme 1, campaign (TolE hard-coded 1e-4 under honorTolR 0)
TAG_ME8 = 2       # IntScheme 1, TolR 1e-8, -honorTolR 1
TAG_RK45 = 3      # IntScheme 45 (RungeKutta45), TolR 1e-10 (dT_min 1e-3 hard-coded!)
TAG_CPPM = 4      # IntScheme 2 (BackwardEuler_CPPM), TolR 1e-10
PROTOS = {TAG_ME: "ME(campaign)", TAG_ME8: "ME(1e-8,honor)",
          TAG_RK45: "RK45(1e-10)", TAG_CPPM: "CPPM(1e-10)"}
SCHEME_OF = {TAG_ME: 1, TAG_ME8: 1, TAG_RK45: 45, TAG_CPPM: 2}


def _opts(scheme, tolr, honor, maxsub=20000):
    return (scheme, 0, 1, 1.0e-7, tolr,
            "-Presidual", 0.0, "-Pmin", 0.0101, "-maxSubsteps", maxsub,
            "-honorTolR", honor, "-flipAlphaIn", "init")


def define_prototypes():
    ops.wipe()
    ops.nDMaterial("LadrunoSANISAND", TAG_ME, *P, *_opts(1, 1.0e-7, 0))
    ops.nDMaterial("LadrunoSANISAND", TAG_ME8, *P, *_opts(1, 1.0e-8, 1))
    ops.nDMaterial("LadrunoSANISAND", TAG_RK45, *P, *_opts(45, 1.0e-10, 1))
    ops.nDMaterial("LadrunoSANISAND", TAG_CPPM, *P, *_opts(2, 1.0e-10, 1, 0))


# ---------------------------------------------------------------------------
# tensor helpers (Voigt xx yy zz xy yz zx, stress-like components)
import math  # noqa: E402


def tr(v):
    return v[0] + v[1] + v[2]


def dev(v):
    t = tr(v) / 3.0
    return [v[0] - t, v[1] - t, v[2] - t, v[3], v[4], v[5]]


def ddot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + 2.0 * (a[3] * b[3] + a[4] * b[4] + a[5] * b[5])


def norm(a):
    return math.sqrt(max(ddot(a, a), 0.0))


def eta_alpha(alpha):
    """sqrt(3/2)||alpha||: the stress ratio alpha alone encodes."""
    return math.sqrt(1.5) * norm(alpha)


def eta_sigma(sig):
    p = tr(sig) / 3.0
    return math.sqrt(1.5) * norm(dev(sig)) / p if p != 0 else float("nan")


def cos3theta_of(n):
    # n deviatoric unit: cos3theta = sqrt(6) tr(n.n.n)
    import numpy as np
    N = np.array([[n[0], n[3], n[5]], [n[3], n[1], n[4]], [n[5], n[4], n[2]]])
    c = math.sqrt(6.0) * float(np.trace(N @ N @ N))
    return max(-1.0, min(1.0, c))


def g_theta(c3):
    return 2 * C_ / ((1 + C_) - (1 - C_) * c3)


def psi_of(e, p):
    return e - (E0 - LAMC * (p / PATM) ** KSI)


def yield_f(sig, alpha):
    p = tr(sig) / 3.0
    s = dev(sig)
    x = [s[i] - p * alpha[i] for i in range(6)]
    return norm(x) - math.sqrt(2.0 / 3.0) * M_ * p


def bounding(sig, alpha, e):
    """The model's own alpha^b_theta at the point (Lode angle of n), the
    bounding-surface image sqrt(2/3) alpha^b n, and b:n, alpha:n."""
    p = tr(sig) / 3.0
    s = dev(sig)
    x = [s[i] - p * alpha[i] for i in range(6)]
    nx = norm(x)
    n = [xi / nx for xi in x] if nx > 1e-300 else [0.0] * 6
    c3 = cos3theta_of(n)
    psi = psi_of(e, p)
    ab = g_theta(c3) * MC * math.exp(-NB * psi) - M_
    an = ddot(alpha, n)
    bn = math.sqrt(2.0 / 3.0) * ab - an
    return dict(p=p, psi=psi, alpha_b=ab, Mb=ab + M_, alpha_dot_n=an,
                b_dot_n=bn, alpha_over_b=math.sqrt(1.5) * norm(alpha) / ab
                if ab > 0 else float("nan"), n=n, c3=c3)
