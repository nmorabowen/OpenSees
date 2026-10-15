// LADRUNO-HEADER-START
// ==========================================================================
//
//   ▄█          ▄████████ ████████▄     ▄████████ ███    █▄  ███▄▄▄▄    ▄██████▄
//  ███         ███    ███ ███   ▀███   ███    ███ ███    ███ ███▀▀▀██▄ ███    ███
//  ███         ███    ███ ███    ███   ███    ███ ███    ███ ███   ███ ███    ███
//  ███         ███    ███ ███    ███  ▄███▄▄▄▄██▀ ███    ███ ███   ███ ███    ███
//  ███       ▀███████████ ███    ███ ▀▀███▀▀▀▀▀   ███    ███ ███   ███ ███    ███
//  ███         ███    ███ ███    ███ ▀███████████ ███    ███ ███   ███ ███    ███
//  ███▌    ▄   ███    ███ ███   ▄███   ███    ███ ███    ███ ███   ███ ███    ███
//  █████▄▄██   ███    █▀  ████████▀    ███    ███ ████████▀   ▀█   █▀   ▀██████▀
//  ▀                                   ███    ███
//
//  Ladruno — a research fork of OpenSees
//  Created by:  Nicolas Mora Bowen  ·  Patricio Palacios  ·  José Abell  ·  Guppi
//
// Header auto-stamped by Ladruno_scripts/stamp_headers.py (art: banner_ASCII.txt).
// Do not hand-edit between the markers; edit the script/art and re-run instead.
// ==========================================================================
// LADRUNO-HEADER-END

// Ladruno WP-194 -- mixed-precision MUMPS (SMUMPS factor, double residual).
// See LadrunoMumpsMixed.h and
// Ladruno_implementation/194_ladruno_mumps_mixed_precision_adr.md.

#include <LadrunoMumpsMixed.h>

#ifdef _LADRUNO_SMUMPS

#include <OPS_Globals.h>
#include <cfloat>
#include <cmath>
#include <cstring>

#ifdef _WIN32
// ---------------------------------------------------------------------------
// Link shim. SMUMPS (ifx: upper-case externals, arguments by reference) calls
// the single-precision BLAS-1 routines SAXPY and SDOT. MKL's mkl_intel_lp64.lib
// defines each in one object together with the lower-case C names `saxpy` and
// `sdot` (and mkl_blas__saxpy/_sdot, which MKL's own cblas_* entry points call,
// so going through CBLAS pulls the same objects). The fork already defines both
// lower-case names itself (`saxpy(double*, double*, double, int)` in
// sparseSYM/utility.c, `sdot` in the in-tree METIS 4 util.c), so pulling those
// MKL objects fails the link with LNK2005. These two reference-BLAS
// implementations (Fortran semantics, negative increments included) resolve
// SAXPY/SDOT before MKL is searched. Callers: SMUMPS's BLR solve (ssol_lr)
// and MKL ScaLAPACK's single-precision PBLAS helpers (the root-node dense
// factorization); all BLAS-1 on vectors, while the factorization kernels
// (SGEMM/STRSM/...) still come from MKL. Linux builds are untested here.
// The lower-case fork functions are untouched and still bind as before.
static inline int ladrunoBlasStart(int n, int inc) { return (inc < 0) ? (1 - n) * inc : 0; }

extern "C" float SDOT(const int *n, const float *x, const int *incx, const float *y, const int *incy)
{
  float s = 0.0f;
  int nn = *n, ix = ladrunoBlasStart(nn, *incx), iy = ladrunoBlasStart(nn, *incy);
  for (int i = 0; i < nn; i++, ix += *incx, iy += *incy)
    s += x[ix] * y[iy];
  return s;
}

extern "C" void SAXPY(const int *n, const float *a, const float *x, const int *incx, float *y, const int *incy)
{
  int nn = *n;
  if (nn <= 0 || *a == 0.0f) return;
  int ix = ladrunoBlasStart(nn, *incx), iy = ladrunoBlasStart(nn, *incy);
  for (int i = 0; i < nn; i++, ix += *incx, iy += *incy)
    y[iy] += (*a) * x[ix];
}
#endif

#define ICNTL(I) icntl[(I)-1]
#define CNTL(I)  cntl[(I)-1]

LadrunoMumpsMixed::LadrunoMumpsMixed(int refineSteps)
  : init(false), needsAnalysis(true), refine(refineSteps), rank(0),
    aF(0), aFsize(0), rhsF(0), rhsSize(0), xw(0), yw(0), rbuf(0), wSize(0),
    tAnalysis(-1.0), tSolveSum(0.0), nSolveSum(0), firstSolveAfterFactor(false)
{
  memset(&sid, 0, sizeof(sid));
  if (refine < 0) refine = 0;
}

LadrunoMumpsMixed::~LadrunoMumpsMixed()
{
  // MPI may already be finalized when a static solver is destroyed; MUMPS'
  // JOB=-2 is collective and needs it.
  int mpiFinal = 0;
  MPI_Finalized(&mpiFinal);
  if (!mpiFinal)
    terminate();
  delete [] aF;
  delete [] rhsF;
  delete [] xw;
  delete [] yw;
  delete [] rbuf;
}

void
LadrunoMumpsMixed::terminate(void)
{
  if (init) {
    sid.job = -2;
    smumps_c(&sid);
    init = false;
  }
  needsAnalysis = true;
}

void
LadrunoMumpsMixed::setControls(const LadrunoMumpsControls &c)
{
  // same controls, same values as the DMUMPS path (MumpsParallelSolver)
  sid.ICNTL(5) = 0; sid.ICNTL(18) = 3;            // assembled, distributed input
  sid.ICNTL(1) = -1; sid.ICNTL(2) = -1; sid.ICNTL(3) = -1; sid.ICNTL(4) = 0;
  sid.ICNTL(14) = c.icntl14;
  sid.ICNTL(7) = c.icntl7;
  sid.ICNTL(35) = c.icntl35;                       // BLR in the float factor
  if (c.cntl7 != 0.0)
    sid.CNTL(7) = (float)c.cntl7;
}

int
LadrunoMumpsMixed::initialize(MPI_Comm comm, int sym, const LadrunoMumpsControls &c)
{
  terminate();
  sid.job = -1;
  sid.par = 1;
  sid.sym = sym;
#ifdef _OPENMPI
  if (comm == MPI_COMM_WORLD)
    sid.comm_fortran = 0;
  else
    sid.comm_fortran = MPI_Comm_c2f(comm);
#else
  sid.comm_fortran = MPI_Comm_c2f(comm);
#endif
  smumps_c(&sid);
  MPI_Comm_rank(comm, &rank);
  if (sid.infog[0] < 0)
    return report("initialization");
  init = true;
  setControls(c);
  return 0;
}

int
LadrunoMumpsMixed::report(const char *phase)
{
  int info = sid.infog[0];
  if (info >= 0)
    return 0;
  if (rank == 0) {
    int info2 = sid.infog[1];
    opserr << "WARNING MumpsParallelSolver -precision single: SMUMPS " << phase
	   << " failed, INFOG(1)=" << info << " INFOG(2)=" << info2;
    switch (info) {
    case -6:  opserr << " (matrix singular in structure)"; break;
    case -10: opserr << " (matrix numerically singular in single precision)"; break;
    case -8: case -9:
      opserr << " (work array too small: raise -ICNTL14)"; break;
    case -5: case -7: case -13:
      opserr << " (out of memory)"; break;
    default: break;
    }
    opserr << "\n";
  }
  return info;
}

// A -> float. A double entry beyond FLT_MAX would become inf and poison the
// factor; refuse instead (collectively: every rank must agree before JOB=1/2).
int
LadrunoMumpsMixed::convert(int nnz, const double *A, MPI_Comm comm)
{
  int need = (nnz > 0) ? nnz : 1;
  if (aFsize < need) {
    delete [] aF;
    aF = new float[need];
    aFsize = need;
  }
  int bad = 0;
  for (int i = 0; i < nnz; i++) {
    double a = A[i];
    if (!(std::fabs(a) <= (double)FLT_MAX))   // also catches NaN
      bad = 1;
    aF[i] = (float)a;
  }
  MPI_Allreduce(MPI_IN_PLACE, &bad, 1, MPI_INT, MPI_MAX, comm);
  if (bad) {
    if (rank == 0)
      opserr << "WARNING MumpsParallelSolver -precision single: a matrix entry is NaN or "
	     << "exceeds the single-precision range (|a| > " << (double)FLT_MAX
	     << "); not factored. Use the default double precision for this model.\n";
    return -1;
  }
  return 0;
}

int
LadrunoMumpsMixed::floatSolve(int n)
{
  if (rank == 0) {
    sid.n = n;
    sid.rhs = rhsF;
  }
  sid.job = 3;
  smumps_c(&sid);
  return report("solve");
}

void
LadrunoMumpsMixed::residual(MPI_Comm comm, int sym, int n, int nnz,
			    const int *rowA, const int *colA, const double *A,
			    const double *B, const double *X)
{
  if (rank == 0)
    for (int i = 0; i < n; i++) xw[i] = X[i];
  MPI_Bcast(xw, n, MPI_DOUBLE, 0, comm);
  for (int i = 0; i < n; i++) yw[i] = 0.0;
  // local triplets, 0-based; duplicates (and the same entry held by several
  // ranks) add up, exactly as MUMPS assembles them. Symmetric storage holds
  // one triangle: mirror the off-diagonal entries.
  for (int k = 0; k < nnz; k++) {
    int i = rowA[k], j = colA[k];
    yw[i] += A[k] * xw[j];
    if (sym != 0 && i != j)
      yw[j] += A[k] * xw[i];
  }
  MPI_Reduce(yw, rbuf, n, MPI_DOUBLE, MPI_SUM, 0, comm);
  if (rank == 0)
    for (int i = 0; i < n; i++) rbuf[i] = B[i] - rbuf[i];
}

static double
infNorm(const double *v, int n)
{
  double m = 0.0;
  for (int i = 0; i < n; i++) {
    double a = std::fabs(v[i]);
    if (a > m) m = a;
  }
  return m;
}

int
LadrunoMumpsMixed::solve(MPI_Comm comm, int sym, int n, int nnz, int *rowA,
			 int *colA, const double *A, const double *B, double *X,
			 bool &factored, const LadrunoMumpsControls &c)
{
  int res = 0;
  bool stats = (c.printStats != 0);

  if (!init) {
    if ((res = initialize(comm, sym, c)) < 0)
      return res;
    needsAnalysis = true;
  }
  setControls(c);

  if (rhsSize < n || wSize < n) {
    delete [] rhsF; delete [] xw; delete [] yw; delete [] rbuf;
    int m = (n > 0) ? n : 1;
    rhsF = new float[m]; xw = new double[m]; yw = new double[m]; rbuf = new double[m];
    rhsSize = wSize = m;
  }

  bool converted = false;

  // ---- analysis (once per sparsity pattern)
  if (needsAnalysis) {
    if ((res = convert(nnz, A, comm)) < 0)
      return res;
    converted = true;
    for (int i = 0; i < nnz; i++) { rowA[i]++; colA[i]++; }
    sid.n = n;
    sid.nz_loc = nnz;
    sid.irn_loc = rowA;
    sid.jcn_loc = colA;
    sid.a_loc = aF;
    sid.job = 1;
    double t0 = stats ? MPI_Wtime() : 0.0;
    smumps_c(&sid);
    if (stats) tAnalysis = MPI_Wtime() - t0;
    for (int i = 0; i < nnz; i++) { rowA[i]--; colA[i]--; }
    if ((res = report("analysis")) < 0)
      return res;
    needsAnalysis = false;
    factored = false;
  }

  // ---- factorization in single precision (the kept factor under WP-188 reuse)
  if (!factored) {
    if (!converted && (res = convert(nnz, A, comm)) < 0)
      return res;
    for (int i = 0; i < nnz; i++) { rowA[i]++; colA[i]++; }
    sid.n = n;
    sid.nz_loc = nnz;
    sid.irn_loc = rowA;
    sid.jcn_loc = colA;
    sid.a_loc = aF;
    sid.job = 2;
    double t0 = stats ? MPI_Wtime() : 0.0;
    smumps_c(&sid);
    double tFactor = stats ? MPI_Wtime() - t0 : 0.0;
    for (int i = 0; i < nnz; i++) { rowA[i]--; colA[i]--; }
    if ((res = report("factorization")) < 0)
      return res;
    factored = true;
    firstSolveAfterFactor = true;

    if (stats) {
      double tmax[3] = { tAnalysis, tFactor, tSolveSum };
      MPI_Allreduce(MPI_IN_PLACE, tmax, 3, MPI_DOUBLE, MPI_MAX, comm);
      int nSolves = nSolveSum;
      bool hadAnalysis = (tmax[0] >= 0.0);
      tAnalysis = -1.0; tSolveSum = 0.0; nSolveSum = 0;
      if (rank == 0) {
	double fEntries = (double)sid.infog[8];
	if (fEntries < 0) fEntries = -fEntries * 1.0e6;
	opserr << "MUMPS stats: n=" << sid.n << " BLR(ICNTL35)=" << c.icntl35;
	if (c.icntl35 != 0)
	  opserr << " eps(CNTL7)=" << c.cntl7;
	opserr << "\n";
	// own line, so the header above keeps the exact text that existing log
	// parsers match with anchored regexes
	opserr << "  precision           = single (SMUMPS factor, double residual), refine "
	       << refine << "\n";
	opserr << "  factor entries INFOG(9)   = " << fEntries << "\n";
	opserr << "  factor MB/proc INFOG(21)  = " << sid.infog[20] << "  (max over procs)\n";
	opserr << "  factor MB tot  INFOG(22)  = " << sid.infog[21] << "  (sum over procs)\n";
	opserr << "  elim flops     RINFOG(3)  = " << sid.rinfog[2] << "\n";
	if (c.icntl35 != 0) {
	  opserr << "  BLR RINFOG(14)            = " << sid.rinfog[13] << "\n";
	  opserr << "  BLR RINFOG(15)            = " << sid.rinfog[14] << "\n";
	}
	if (hadAnalysis)
	  opserr << "  t_analysis s (max)  = " << tmax[0] << "\n";
	opserr << "  t_factor s (max)    = " << tmax[1] << "\n";
	opserr << "  t_solve s (max, n)  = " << tmax[2] << ", " << nSolves << "\n";
      }
    }
  }

  // ---- solve: scaled float RHS through the float factor, double correction.
  // Scaling by ||b||_inf keeps the float RHS clear of overflow/underflow.
  double t0 = stats ? MPI_Wtime() : 0.0;
  double s = 1.0;
  if (rank == 0) {
    s = infNorm(B, n);
    if (!(s > 0.0) || !std::isfinite(s)) s = 1.0;
    for (int i = 0; i < n; i++) rhsF[i] = (float)(B[i] / s);
  }
  if ((res = floatSolve(n)) < 0)
    return res;
  if (rank == 0)
    for (int i = 0; i < n; i++) X[i] = s * (double)rhsF[i];

  // residual history, only on the first solve after a factorization and only
  // with -stats (it costs an extra matvec + reduction per line)
  bool history = stats && firstSolveAfterFactor;
  double bnorm = (rank == 0) ? infNorm(B, n) : 0.0;
  if (bnorm == 0.0) bnorm = 1.0;
  if (history) {
    residual(comm, sym, n, nnz, rowA, colA, A, B, X);
    if (rank == 0)
      opserr << "  mixed: ||b-Ax||inf/||b||inf  step 0 = " << infNorm(rbuf, n) / bnorm << "\n";
  }

  // ---- -refine k: classical mixed-precision iterative refinement
  for (int it = 1; it <= refine; it++) {
    if (!history)                 // with history the residual is already in rbuf
      residual(comm, sym, n, nnz, rowA, colA, A, B, X);
    double rs = 1.0;
    if (rank == 0) {
      rs = infNorm(rbuf, n);
      if (!(rs > 0.0) || !std::isfinite(rs)) rs = 1.0;
      for (int i = 0; i < n; i++) rhsF[i] = (float)(rbuf[i] / rs);
    }
    if ((res = floatSolve(n)) < 0)
      return res;
    if (rank == 0)
      for (int i = 0; i < n; i++) X[i] += rs * (double)rhsF[i];
    if (history) {
      residual(comm, sym, n, nnz, rowA, colA, A, B, X);
      if (rank == 0)
	opserr << "  mixed: ||b-Ax||inf/||b||inf  step " << it << " = "
	       << infNorm(rbuf, n) / bnorm << "\n";
    }
  }
  firstSolveAfterFactor = false;

  if (stats) {
    tSolveSum += MPI_Wtime() - t0;
    nSolveSum++;
  }
  return 0;
}

#endif // _LADRUNO_SMUMPS
