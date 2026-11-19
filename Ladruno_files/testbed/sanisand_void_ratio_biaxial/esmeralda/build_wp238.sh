#!/bin/bash
# WP-238 esmeralda build: the OpenSees worktree of ~/ladruno_wp238 (branch pushed from the
# PC into ~/ladruno_wp152/OpenSees, then `git worktree add`), same configuration as
# build_wp152.sh (Linux MKL PARDISO opt-in, sequential layer), OpenSeesPy target only.
#   srun -N1 -n1 -c32 --exclusive bash build_wp238.sh
set -u
ROOT=$HOME/ladruno_wp238
SRC=$ROOT/OpenSees
M=/mnt/nfshare/software/spack/opt/spack/linux-zen3/intel-oneapi-mkl-2024.2.2-4olf4bnrh4yzaqbhq742doqvnim6yynl/mkl/2024.2
export PATH=$HOME/ladruno_build_test/conan_venv/bin:$PATH
export TMPDIR=$ROOT/tmp
mkdir -p $TMPDIR $ROOT/bin $ROOT/logs
echo "host=$(hostname) cpus=$(nproc) start=$(date -Is)"
cd $SRC || exit 1
echo "HEAD=$(git log --oneline -1)"; git status --short | head
CM=$HOME/ladruno_wp138/cmake_venv/bin/cmake
TC=$HOME/ladruno_wp152/OpenSees/build/Release/generators/conan_toolchain.cmake
BDIR=build/wp238_seq
LAP="$M/lib/libmkl_intel_lp64.so;$M/lib/libmkl_sequential.so;$M/lib/libmkl_core.so;-lm;-ldl"
$CM -S . -B $BDIR -DCMAKE_TOOLCHAIN_FILE=$TC -DCMAKE_BUILD_TYPE=Release \
  -DLADRUNO_MKL_PARDISO_LINUX=ON -DLADRUNO_MKL_PARDISO_LINUX_THREADED=OFF -DMKL_RT_HINT=$M/lib \
  "-DLAPACK_LIBRARIES=$LAP" "-DBLAS_LIBRARIES=$LAP" > $ROOT/logs/configure.log 2>&1 \
  || { echo CONFIGURE_FAIL; tail -30 $ROOT/logs/configure.log; exit 1; }
$CM --build $BDIR --target OpenSeesPy -j30 > $ROOT/logs/build.log 2>&1
RC=$?
echo "BUILD_EXIT=$RC end=$(date -Is)"
if [ $RC -ne 0 ]; then grep -nE "error" $ROOT/logs/build.log | head -30; exit 1; fi
cp $BDIR/OpenSeesPy.so $ROOT/bin/opensees.so
export LD_LIBRARY_PATH=$M/lib:${LD_LIBRARY_PATH:-}
export MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_CBWR=COMPATIBLE LADRUNO_OPENSEES_QUIET=1
PY=$HOME/ladruno_build_test/conan_venv/bin/python
cd $ROOT/bin
$PY -S -c "import sys; sys.path.insert(0,'.'); import opensees as ops; print('build', ops.ladrunoBuild().strip().splitlines()[0])"
echo "SMOKE_EXIT=$?"
# Linux check of the WP-238 tests (the generator and the per-element parameter)
cd $SRC/tests
export PYTHONPATH=$ROOT/bin:$SRC/tests:$SRC/Ladruno_scripts
$PY -m pytest -q -p no:cacheprovider test_sanisand_void_ratio_field.py test_ladruno_sanisand_void_ratio_param.py 2>&1 | tail -15
echo "PYTEST_EXIT=${PIPESTATUS[0]}"
