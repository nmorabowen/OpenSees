#!/usr/bin/env bash
# WP-130 / F18(c), second round: the IntScheme-2 arms with -cppmTangent fixed
# (the sign-corrected algorithmic tangent), SEQUENTIALLY. Needs LADRUNO_DIST_BIN, PY.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
W=${WALL:-1200}
A="$HERE/run_arm.sh"
echo "=== s2t_T2: IntScheme 2, -cppmTangent fixed, vanilla ladder ===";  bash $A s2t_T2 2 $W -- -cppmTangent fixed
echo "=== s2t_ref0: + -cppmOnFail refuse -cppmHalvings 0 ===";          bash $A s2t_ref0 2 $W -- -cppmTangent fixed -cppmOnFail refuse -cppmHalvings 0
echo "=== s2t_ref0_se_ls: + -cppmStart explicit -cppmLineSearch on ==="; bash $A s2t_ref0_se_ls 2 $W -- -cppmTangent fixed -cppmOnFail refuse -cppmHalvings 0 -cppmStart explicit -cppmLineSearch on
echo "=== s2t_ref3_se_ls: -cppmHalvings 3 + start + LS ===";           bash $A s2t_ref3_se_ls 2 $W -- -cppmTangent fixed -cppmOnFail refuse -cppmHalvings 3 -cppmStart explicit -cppmLineSearch on
echo "=== DONE ==="
