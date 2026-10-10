#!/bin/bash
# WP-138 esmeralda build: fresh clone of ladruno at a pinned SHA, Linux MKL
# PARDISO opt-in (SEQUENTIAL 3-layer, spack oneMKL 2024.2.2), OpenSeesPy target.
# Run ON A COMPUTE NODE:  srun -N1 -n1 -c32 --exclusive bash build_wp138.sh <sha>
set -u
ROOT=$HOME/ladruno_wp138
SRC=$ROOT/OpenSees
M=/mnt/nfshare/software/spack/opt/spack/linux-zen3/intel-oneapi-mkl-2024.2.2-4olf4bnrh4yzaqbhq742doqvnim6yynl/mkl/2024.2
export PATH=$HOME/ladruno_build_test/conan_venv/bin:$PATH
export TMPDIR=$HOME/ladruno_build_test/tmp
mkdir -p $TMPDIR $ROOT/bin
echo "host=$(hostname) cpus=$(nproc) start=$(date -Is)"
# clone + patch + conan install were done on the HEAD node (compute nodes have
# no internet and no cmake); cmake comes from a pip venv on shared storage.
cd $SRC || exit 1
echo "HEAD=$(git log --oneline -1)"; git diff --stat
CM=$ROOT/cmake_venv/bin/cmake
$CM --version | head -1
TC=$PWD/build/Release/generators/conan_toolchain.cmake
BDIR=build/wp138_seq
LAP="$M/lib/libmkl_intel_lp64.so;$M/lib/libmkl_sequential.so;$M/lib/libmkl_core.so;-lm;-ldl"
$CM -S . -B $BDIR -DCMAKE_TOOLCHAIN_FILE=$TC -DCMAKE_BUILD_TYPE=Release \
  -DLADRUNO_MKL_PARDISO_LINUX=ON -DLADRUNO_MKL_PARDISO_LINUX_THREADED=OFF -DMKL_RT_HINT=$M/lib \
  "-DLAPACK_LIBRARIES=$LAP" "-DBLAS_LIBRARIES=$LAP" > $ROOT/configure.log 2>&1 || { echo CONFIGURE_FAIL; tail -30 $ROOT/configure.log; exit 1; }
grep -E "PARDISO|Ladruno: MKL|LAPACK_LIBRARIES|OPENMP|Python" $ROOT/configure.log
$CM --build $BDIR --target OpenSeesPy -j30 > $ROOT/build.log 2>&1
RC=$?
echo "BUILD_EXIT=$RC end=$(date -Is)"
if [ $RC -ne 0 ]; then grep -nE "error" $ROOT/build.log | head -30; exit 1; fi
cp $BDIR/OpenSeesPy.so $ROOT/bin/opensees.so
ldd $ROOT/bin/opensees.so | grep -iE "mkl|blas|lapack|gomp|python|not found"
# smoke: import + ladrunoBuild + trivial truss + Pardiso solve
export LD_LIBRARY_PATH=$M/lib:${LD_LIBRARY_PATH:-}
export MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_CBWR=COMPATIBLE LADRUNO_OPENSEES_QUIET=1
cd $ROOT/bin
$HOME/ladruno_build_test/conan_venv/bin/python -S - <<'EOF'
import sys; sys.path.insert(0, '.')
import opensees as ops
print("file", ops.__file__)
print("build", ops.ladrunoBuild().strip().splitlines()[0])
for sysname in ("Pardiso", "UmfPack"):
    ops.wipe(); ops.model('basic', '-ndm', 2, '-ndf', 2)
    ops.node(1, 0, 0); ops.node(2, 1, 0); ops.fix(1, 1, 1); ops.fix(2, 0, 1)
    ops.uniaxialMaterial('Elastic', 1, 1000.0)
    ops.element('Truss', 1, 1, 2, 1.0, 1)
    ops.timeSeries('Linear', 1); ops.pattern('Plain', 1, 1); ops.load(2, 1.0, 0.0)
    ops.system(sysname); ops.numberer('RCM'); ops.constraints('Plain')
    ops.integrator('LoadControl', 1.0); ops.algorithm('Linear'); ops.analysis('Static')
    rc = ops.analyze(1)
    print(sysname, "rc", rc, "disp", ops.nodeDisp(2, 1), "(expect 0.001)")
EOF
echo "SMOKE_EXIT=$?"
