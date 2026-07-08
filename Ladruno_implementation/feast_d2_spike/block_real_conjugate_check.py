import numpy as np
rng = np.random.default_rng(0)
n = 60
A = rng.standard_normal((n, n))
K = A @ A.T + n * np.eye(n)          # SPD-ish stiffness
Bm = rng.standard_normal((n, n))
M = Bm @ Bm.T + n * np.eye(n)        # SPD mass
Y = rng.standard_normal((n, 4))
RHS = M @ Y

z = 3.7 + 1.9j   # a FEAST contour node, non-real shift
Xc = np.linalg.solve(z * M - K, RHS)   # true complex solve (baseline oracle)

a, b = z.real, z.imag

# Non-symmetric block form: [[aM-K, -bM], [bM, aM-K]] @ [Xr; Xi] = [B; 0]
top_ns = np.hstack([a*M - K, -b*M])
bot_ns = np.hstack([ b*M,   a*M - K])
Areal_ns = np.vstack([top_ns, bot_ns])

# Symmetric indefinite form (negate the imaginary-row block):
# [[aM-K, -bM], [-bM, -(aM-K)]] @ [Xr; Xi] = [B; 0] -- MUMPS SYM=2 / LDL^T target
top_sym = np.hstack([ a*M - K, -b*M])
bot_sym = np.hstack([-b*M,   -(a*M - K)])
Areal_sym = np.vstack([top_sym, bot_sym])

Bg = np.vstack([RHS, np.zeros_like(RHS)])

for label, Areal in [("non-symmetric", Areal_ns), ("symmetric (SYM=2 target)", Areal_sym)]:
    Xg = np.linalg.solve(Areal, Bg)
    Xr_from_block = Xg[:n]
    Xi_from_block = Xg[n:]
    Xblock = Xr_from_block + 1j*Xi_from_block
    err = np.max(np.abs(Xblock - Xc)) / np.max(np.abs(Xc))
    print(f"[{label}] relative error vs true complex solve: {err:.3e}"
          f" | real: {np.isrealobj(Areal)} | shape: {Areal.shape} (orig n={n})"
          f" | symmetric: {np.allclose(Areal, Areal.T)}")
    assert err < 1e-12, f"{label} block-real form diverged from the complex solve: {err:.3e}"
assert np.allclose(Areal_sym, Areal_sym.T), "SYM=2 target form lost symmetry"
print("OK: both block-real forms reproduce the complex solve (asserted < 1e-12)")
