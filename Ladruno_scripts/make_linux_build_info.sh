#!/usr/bin/env bash
# WP-166: write BUILD_INFO.txt for the Linux opensees.so artifact / release.
# usage: make_linux_build_info.sh <build_dir_with_OpenSees_binary> <out_file>
# Best effort for `ladrunoBuild` (needs the conan Tcl runtime; TCL_LIBRARY set).
set -u
bdir="${1:?build dir}"; out="${2:?output file}"
{
  echo "git_sha: $(git rev-parse HEAD)"
  echo "git_ref: ${GITHUB_REF:-$(git rev-parse --abbrev-ref HEAD)}"
  echo "built_utc: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "python_abi: $(python -c 'import sys,sysconfig;print("cp%d%d"%sys.version_info[:2], sysconfig.get_config_var("EXT_SUFFIX"))')"
  echo "python_version: $(python -V 2>&1)"
  echo "os: $(. /etc/os-release 2>/dev/null; echo "${PRETTY_NAME:-$(uname -sr)}")"
  echo "--- ladrunoBuild ---"
  echo 'puts [ladrunoBuild]' > "$bdir/_lb.tcl"
  "$bdir/OpenSees" "$bdir/_lb.tcl" 2>&1 | head -40 || true
  rm -f "$bdir/_lb.tcl"
} > "$out"
cat "$out"
