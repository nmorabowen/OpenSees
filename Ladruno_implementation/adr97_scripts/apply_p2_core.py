"""ADR-97 P2 (wp/97c) -- re-runnable, idempotent staging script for the two big
files: ASDPlasticMaterial3D.h (the principal-space return + dispatch) and
OPS_AllASDPlasticMaterial3Ds.cpp (the parse-time refusal message).

Run AFTER apply_p2_cpp.py, from the worktree root:
    python3.12 Ladruno_implementation/adr97_scripts/apply_p2_core.py
"""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SRC = os.path.join(ROOT, "SRC", "material", "nD", "ASDPlasticMaterial3D")

EDITS = []


def edit(relpath, anchor, new, tag):
    EDITS.append((relpath, anchor, new, tag))


# ===========================================================================
# A. explicit Eigen eigen-solver include
# ===========================================================================
edit(
    "ASDPlasticMaterial3D.h",
    '#include "AllASDHardeningFunctions.h"\n',
    '#include "AllASDHardeningFunctions.h"\n'
    "#include <Eigen/Eigenvalues>   // Ladruno (ADR-97 wp/97c): SelfAdjointEigenSolver\n",
    "Eigen/Eigenvalues include",
)

# ===========================================================================
# B. the all-IVs-inert compile-time fold
# ===========================================================================
edit(
    "ASDPlasticMaterial3D.h",
    """template <typename... Ts>
struct asdp_all_ivs_support_cp<std::tuple<Ts...>>
{
    static constexpr bool value = (Ts::hardening_supports_cp() && ... && true);
};
""",
    """template <typename... Ts>
struct asdp_all_ivs_support_cp<std::tuple<Ts...>>
{
    static constexpr bool value = (Ts::hardening_supports_cp() && ... && true);
};

// Ladruno (ADR-97 wp/97c): is EVERY internal variable of this specialization
// FIXED (h == 0)?  The principal-stress-space Mohr-Coulomb return of P2 is a
// closed-form projection onto a surface it assumes does not move: it carries no
// q-row, so a live hardening law would silently be ignored.  Folding this at
// compile time is what keeps a hypothetical MohrCoulomb_YF<ArmstrongFrederick..>
// out of that path and REFUSED at parse time instead.
template <typename Tuple>
struct asdp_all_ivs_are_inert;

template <typename... Ts>
struct asdp_all_ivs_are_inert<std::tuple<Ts...>>
{
    static constexpr bool value = (Ts::hardening_is_perfectly_plastic() && ... && true);
};
""",
    "asdp_all_ivs_are_inert fold",
)

# ===========================================================================
# C. support predicate: smooth OR principal
# ===========================================================================
edit(
    "ASDPlasticMaterial3D.h",
    """    static constexpr bool ladruno_cp_supported =
        yf_has_cp_derivatives<YieldFunctionType>::value &&
        pf_has_cp_derivatives<PlasticFlowType>::value &&
        asdp_all_ivs_support_cp<iv_concat_types>::value;

    static constexpr bool supportsClosestPoint() { return ladruno_cp_supported; }""",
    """    static constexpr bool ladruno_cp_smooth_supported =
        yf_has_cp_derivatives<YieldFunctionType>::value &&
        pf_has_cp_derivatives<PlasticFlowType>::value &&
        asdp_all_ivs_support_cp<iv_concat_types>::value;

    // Ladruno (ADR-97 wp/97c): the PRINCIPAL-STRESS-SPACE family of this
    // specialization (0 = none).  Non-zero only when the yield function and the
    // plastic flow direction carry the SAME non-zero marker -- so the generator's
    // cross pairings (MohrCoulomb_YF x VonMises_PF, VonMises_YF x MohrCoulomb_PF,
    // HoekBrown_YF x MohrCoulomb_PF, ...) are NOT enabled by P2: the principal
    // return assumes both the surface and the potential are the piecewise-linear
    // Mohr-Coulomb ones, and no oracle covers the mixed maps.  They stay refused.
    // Perfect plasticity is required as well (the map has no q-row).
    static constexpr int ladruno_cp_principal_family =
        (yf_cp_principal_family<YieldFunctionType>::value != 0
         && yf_cp_principal_family<YieldFunctionType>::value
            == pf_cp_principal_family<PlasticFlowType>::value
         && asdp_all_ivs_are_inert<iv_concat_types>::value)
        ? yf_cp_principal_family<YieldFunctionType>::value : 0;

    static constexpr bool ladruno_cp_supported =
        ladruno_cp_smooth_supported || (ladruno_cp_principal_family != 0);

    static constexpr bool supportsClosestPoint() { return ladruno_cp_supported; }""",
    "ladruno_cp_supported = smooth || principal",
)

# ===========================================================================
# D. the principal-space return map itself, inserted before Closest_Point
# ===========================================================================
PRINCIPAL = r'''
    //==================================================================================================
    // Ladruno (ADR-97 wp/97c): PRINCIPAL-STRESS-SPACE multi-surface closest-point
    // return for the Mohr-Coulomb family (Clausen, Damkilde & Andersen 2006/2007)
    //==================================================================================================
    //
    // WHY A SECOND ALGORITHM (ADR-97 D3).  The shipped 6D Mohr-Coulomb yield
    // function is the exact Lode-angle form
    //     f = A(theta) sqrt(J2) + I1 sin(phi)/3 - c cos(phi),
    //     A(theta) = cos(theta) - sin(theta) sin(phi)/sqrt(3),
    // whose gradient carries a 1/cos(3 theta) that BLOWS UP on the corners; the
    // header dodges it by swapping in a Drucker-Prager gradient for
    // |theta| >= 29 deg and otherwise central-differencing f over the six raw
    // Voigt slots.  A coupled Newton built on that gradient cannot converge
    // quadratically at a corner -- which is precisely where Mohr-Coulomb models
    // spend their time.  In PRINCIPAL STRESS SPACE the same surface is six
    // PLANES: on the sorted sextant s1 >= s2 >= s3 (tension positive) it is the
    // single plane
    //     f = a_13 . s - k ,  a_ij = [ (1+sin phi)/2, 0, -(1-sin phi)/2 ], k = c cos(phi),
    // its two corners are LINES and its vertex is the point (k/sin phi)*[1,1,1].
    // Every return is then a LINEAR projection in the elastic metric -- closed
    // form, no Newton at all -- and the tangent is constant on each region.
    // The P0 oracle `adr97_oracle/cppm_mc.py` verifies numerically, over 2000
    // random states, that this principal-stress f IS the header's own invariant
    // expression to 1e-14 relative; every number below is pinned against it.
    //
    // FLOW DIRECTION.  The header's `m` is deviator(dg/dsigma with PHI) plus
    // sin(psi)/3 * delta -- the deviatoric shape of the phi-surface with a
    // psi-controlled volumetric part, NOT the textbook non-associated gradient.
    // In principal space that is m_ij = a_ij - (sin phi)/3 * 1 + (sin psi)/3 * 1,
    // which collapses to a_ij when psi == phi.  Mirrored here, not corrected.
    //
    // REGION SELECTION IS BY BOUNDARY PLANES, NEVER BY AN ACTIVE-SET SEARCH.
    // Each edge line L has a boundary plane spanned by the line direction ell_L
    // and the face return direction rp = D m_13, passing through the apex; the
    // side of it the trial principal point falls on decides face vs edge, and the
    // sign of the edge-line parameter of the edge return decides edge vs apex.
    // This is exact and branch-free, and it REPLACES the |theta| < 29 deg guard.
    // It is also not optional: at the apex with psi < phi the three Koiter
    // multipliers are NOT all positive (the cone of return directions no longer
    // contains the hydrostatic direction), so a "grow the active set while
    // dLambda >= 0" search classifies those states wrongly -- ADR-97 P0 header
    // finding 6, and the reason the oracle prints the multipliers with that note.
    //
    // The boundary-plane SIGNS are taken analytically here, not from the oracle's
    // dimensional reference point: any point strictly inside the open face is
    // apex + b1*ell_1 + b2*ell_2 with b1, b2 > 0, and n_1 . ell_1 == 0 by
    // construction, so sgn_1 = sign(n_1 . ell_2) and symmetrically.  Verified
    // identical to the oracle's calibration over 60 (phi, psi, c) combinations
    // and 13094 trial states, mismatches 0 (see the P2 report).
    //
    // BACK-TRANSFORM.  sigma_ret is an isotropic tensor function of sigma_tr, so
    //     d sigma / d sigma_tr = Rs * T * Rs^{-1}
    // with Rs the Voigt image of E -> Q E Q^T (Q the TRIAL eigenvectors),
    // T's normal block = dy/dx (the principal 3x3 Koiter block) and T's shear
    // slots = (y_i - y_j)/(x_i - x_j) -- the eigenprojection ROTATION term.  That
    // ratio is 0/0 on every axisymmetric path (most triaxial decks, so this is
    // the common case, not the corner case); the l'Hopital limit is
    // dy_i/dx_i - dy_i/dx_j.  The switch threshold is RELATIVE to the yield
    // function's strength scale (ADR-94 M5: an absolute one reintroduces the
    // unit dependence that decided pass/fail in kPa vs Pa).  Rs^{-1} is built
    // directly as the Voigt image of E -> Q^T E Q, not inverted numerically.
    //
    // NO NEWTON, so `cp_iterations` reads 1 on a plastic step here and 0 on an
    // elastic one -- it is a region count, not an iteration count.
    //
    // `Backward_Euler` is untouched (ADR-97 D1): it still runs the scalar Newton
    // on the smoothed 6D gradient, and for MohrCoulombTensionCutoff still calls
    // the same ADR-84 `special_return` hook this map calls.

    // Apply the material's configured tangent_type to a raw active-set (Koiter)
    // tangent.  Mirrors the policy ADR-84 P3 set for `special_return` and the one
    // P1 set for the Drucker-Prager apex: `Continuum` maps to the same raw
    // operator, because the continuum elastoplastic operator of a multi-surface
    // corner is not defined by a single n : E : m.
    void cp_apply_tangent_policy(const VoigtMatrix& raw, const VoigtMatrix& Eelastic)
    {
        using TOT = ASDPlasticMaterial3D_Tangent_Operator_Type;
        switch (INT_OPT_tangent_operator_type[ASDP_TAG])
        {
        case TOT::Elastic:
            Stiffness = Eelastic;
            break;
        case TOT::Continuum:
        case TOT::Algorithmic:
        case TOT::Numerical_Algorithmic_FirstOrder:
        case TOT::Numerical_Algorithmic_SecondOrder:
            Stiffness = raw;
            break;
        case TOT::Secant:
        default:
            Stiffness = VoigtMatrix((raw + Eelastic) / 2.0);
            break;
        }
    }

    int cp_principal_return(const VoigtVector& depsilon,
                            const VoigtVector& sigma_tr,
                            const VoigtMatrix& Eelastic,
                            double tol_f)
    {
        using namespace ASDPlasticMaterial3DGlobals;
        (void) depsilon;

        // ---- surface constants ------------------------------------------
        double sin_phi = 0.0, k_coh = 0.0, sin_phi_pf = 0.0, sin_psi = 0.0;
        if (!yf.cp_mc_face_params(iv_storage, parameters_storage, sin_phi, k_coh)
                || !pf.cp_mc_flow_params(iv_storage, parameters_storage,
                                         sin_phi_pf, sin_psi))
        {
            opserr << "ASDPlasticMaterial3D::Closest_Point (tag " << ASDP_TAG
                   << ") - this yield function / plastic flow pair is marked as a"
                   << " principal-space Mohr-Coulomb family but does not supply the"
                   << " face parameters (ADR-97 P2) -- rejecting step" << endln;
            return LADRUNO_MATERIAL_REFUSED;
        }
        {
            const double d = sin_phi - sin_phi_pf;
            if ((d > 1e-12) || (d < -1e-12))
            {
                opserr << "ASDPlasticMaterial3D::Closest_Point (tag " << ASDP_TAG
                       << ") - the yield function and the plastic flow direction"
                       << " disagree on MC_phi (" << sin_phi << " vs " << sin_phi_pf
                       << ") -- rejecting step (ADR-97 P2)" << endln;
                return LADRUNO_MATERIAL_REFUSED;
            }
        }
        if (!(sin_phi > 100 * MACHINE_EPSILON))
        {
            // phi == 0 is Tresca: the Mohr-Coulomb apex is at infinity, and
            // Clausen's boundary planes are defined THROUGH it.  Rather than
            // silently reclassify, refuse and name the alternative.
            opserr << "ASDPlasticMaterial3D::Closest_Point (tag " << ASDP_TAG
                   << ") - integration_method Closest_Point needs MC_phi > 0"
                   << " (ADR-97 P2): at phi = 0 the Mohr-Coulomb surface is a"
                   << " Tresca prism whose apex is at infinity, so the"
                   << " boundary-plane region test is undefined."
                   << " Use Backward_Euler." << endln;
            return LADRUNO_MATERIAL_REFUSED;
        }

        // ---- isotropy of the elastic tangent ----------------------------
        // The principal-space projection needs D in the principal frame, which
        // exists only for an isotropic E.  Same guard ADR-84's special_return
        // uses; LinearIsotropic3D_EL is the only elasticity registered with this
        // family, so this is a future-proofing refusal, not a live branch.
        const double lam = Eelastic(0, 1);
        const double Gmu = (Eelastic(0, 0) - Eelastic(0, 1)) / 2.0;
        {
            const double sE = ((Eelastic(0, 0) < 0) ? -Eelastic(0, 0) : Eelastic(0, 0))
                              + ((Gmu < 0) ? -Gmu : Gmu);
            const double t = 1e-8 * sE;
            const double d02 = Eelastic(0, 2) - lam;
            const double d12 = Eelastic(1, 2) - lam;
            const double d33 = Eelastic(3, 3) - Gmu;
            const double d11 = Eelastic(1, 1) - Eelastic(0, 0);
            if (!(Gmu > 0.0)
                    || (d02 > t) || (d02 < -t) || (d12 > t) || (d12 < -t)
                    || (d33 > t) || (d33 < -t) || (d11 > t) || (d11 < -t))
            {
                opserr << "ASDPlasticMaterial3D::Closest_Point (tag " << ASDP_TAG
                       << ") - the principal-space Mohr-Coulomb return needs an"
                       << " ISOTROPIC elastic tangent -- rejecting step (ADR-97 P2)"
                       << endln;
                return LADRUNO_MATERIAL_REFUSED;
            }
        }
        Eigen::Matrix3d D3;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                D3(i, j) = lam + ((i == j) ? 2.0 * Gmu : 0.0);

        // ---- spectral decomposition of the TRIAL stress, DESCENDING ------
        Eigen::Matrix3d st;
        st << sigma_tr.v11(), sigma_tr.v12(), sigma_tr.v13(),
              sigma_tr.v12(), sigma_tr.v22(), sigma_tr.v23(),
              sigma_tr.v13(), sigma_tr.v23(), sigma_tr.v33();
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(st);
        if (es.info() != Eigen::Success)
        {
            opserr << "ASDPlasticMaterial3D::Closest_Point (tag " << ASDP_TAG
                   << ") - the trial stress eigen-decomposition failed"
                   << " -- rejecting step" << endln;
            return LADRUNO_MATERIAL_REFUSED;
        }
        Eigen::Vector3d xv;
        Eigen::Matrix3d Q;
        for (int i = 0; i < 3; ++i)          // eigenvalues come out ASCENDING
        {
            xv(i) = es.eigenvalues()(2 - i);
            Q.col(i) = es.eigenvectors().col(2 - i);
        }

        // ---- the three surfaces of the sorted sextant --------------------
        // index 0 = (1,3), 1 = (2,3), 2 = (1,2), in 0-based principal slots.
        static const int PIJ[3][2] = {{0, 2}, {1, 2}, {0, 1}};
        const double Aco = 0.5 * (1.0 + sin_phi);
        const double Bco = 0.5 * (1.0 - sin_phi);
        Eigen::Vector3d av[3], mv[3];
        for (int t = 0; t < 3; ++t)
        {
            av[t].setZero();
            av[t](PIJ[t][0]) =  Aco;
            av[t](PIJ[t][1]) = -Bco;
            mv[t] = av[t];
            const double vol = (sin_psi - sin_phi) / 3.0;
            mv[t](0) += vol; mv[t](1) += vol; mv[t](2) += vol;
        }

        const Eigen::Vector3d apexv = Eigen::Vector3d::Constant(k_coh / sin_phi);
        Eigen::Vector3d Dm13 = D3 * mv[0];

        // edge directions: null(A) of the two active gradients = their cross
        // product, oriented (like the oracle) so ell(0) >= ell(2).
        Eigen::Vector3d ell1 = av[0].cross(av[1]);      // edge s1 == s2
        if (ell1(0) < ell1(2)) ell1 = -ell1;
        ell1.normalize();
        Eigen::Vector3d ell2 = av[2].cross(av[0]);      // edge s2 == s3
        if (ell2(0) < ell2(2)) ell2 = -ell2;
        ell2.normalize();
        Eigen::Vector3d nb1 = ell1.cross(Dm13);
        Eigen::Vector3d nb2 = ell2.cross(Dm13);
        const double dot12 = nb1.dot(ell2);
        const double dot21 = nb2.dot(ell1);
        if (!(dot12 > 0.0 || dot12 < 0.0) || !(dot21 > 0.0 || dot21 < 0.0))
        {
            opserr << "ASDPlasticMaterial3D::Closest_Point (tag " << ASDP_TAG
                   << ") - degenerate Mohr-Coulomb boundary-plane geometry"
                   << " -- rejecting step (ADR-97 P2)" << endln;
            return LADRUNO_MATERIAL_REFUSED;
        }
        const double sgn1 = (dot12 > 0.0) ? 1.0 : -1.0;
        const double sgn2 = (dot21 > 0.0) ? 1.0 : -1.0;

        // ---- region classification --------------------------------------
        const Eigen::Vector3d dtr = xv - apexv;
        const double pl1 = sgn1 * nb1.dot(dtr);
        const double pl2 = sgn2 * nb2.dot(dtr);

        Eigen::Vector3d yv;
        Eigen::Matrix3d dydx = Eigen::Matrix3d::Zero();
        const Eigen::Matrix3d I3 = Eigen::Matrix3d::Identity();
        int region = 0;                      // 0 face, 1 edge s1==s2, 2 edge s2==s3, 3 apex

        if (pl1 >= 0.0 && pl2 >= 0.0)
        {
            const double den = av[0].dot(Dm13);
            if (!(den > MACHINE_EPSILON))
            {
                opserr << "ASDPlasticMaterial3D::Closest_Point (tag " << ASDP_TAG
                       << ") - singular Mohr-Coulomb face return (a . E m = " << den
                       << ") -- rejecting step" << endln;
                return LADRUNO_MATERIAL_REFUSED;
            }
            const double dl = (av[0].dot(xv) - k_coh) / den;
            if (dl < 0.0)
            {
                opserr << "ASDPlasticMaterial3D::Closest_Point (tag " << ASDP_TAG
                       << ") - the Mohr-Coulomb face return produced a NEGATIVE"
                       << " plastic multiplier (" << dl << ") from a trial state"
                       << " outside the surface -- rejecting step" << endln;
                return LADRUNO_MATERIAL_REFUSED;
            }
            yv = xv - dl * Dm13;
            Eigen::Matrix3d outer = Dm13 * av[0].transpose();
            dydx = I3 - outer / den;
            region = 0;
        }
        else
        {
            const bool use1 = (pl1 < pl2);
            const int t0 = use1 ? 0 : 2;     // LINE1 = (13,23); LINE2 = (12,13)
            const int t1 = use1 ? 1 : 0;
            Eigen::Matrix<double, 2, 3> Am;
            Am.row(0) = av[t0].transpose();
            Am.row(1) = av[t1].transpose();
            Eigen::Matrix<double, 3, 2> Mm;
            Mm.col(0) = mv[t0];
            Mm.col(1) = mv[t1];
            Eigen::Matrix<double, 3, 2> DM = D3 * Mm;
            Eigen::Matrix2d ADM = Am * DM;
            const double det = ADM(0, 0) * ADM(1, 1) - ADM(0, 1) * ADM(1, 0);
            const double dscale = (ADM.cwiseAbs().maxCoeff() + MACHINE_EPSILON);
            const double adet = (det < 0) ? -det : det;
            if (!(adet > MACHINE_EPSILON * dscale * dscale))
            {
                opserr << "ASDPlasticMaterial3D::Closest_Point (tag " << ASDP_TAG
                       << ") - singular Mohr-Coulomb edge return (det A E M = "
                       << det << ") -- rejecting step" << endln;
                return LADRUNO_MATERIAL_REFUSED;
            }
            Eigen::Matrix2d ADMi;
            ADMi(0, 0) =  ADM(1, 1) / det;  ADMi(0, 1) = -ADM(0, 1) / det;
            ADMi(1, 0) = -ADM(1, 0) / det;  ADMi(1, 1) =  ADM(0, 0) / det;
            Eigen::Vector2d rhs;
            rhs(0) = av[t0].dot(xv) - k_coh;
            rhs(1) = av[t1].dot(xv) - k_coh;
            Eigen::Vector2d dl2 = ADMi * rhs;
            Eigen::Vector3d yline = xv - DM * dl2;
            const double tpar = (yline - apexv).dot(use1 ? ell1 : ell2);
            if (tpar < 0.0)
            {
                // Beyond the apex along the edge: the vertex is the return.  The
                // three Koiter multipliers there need NOT be positive for psi <
                // phi, which is exactly why this is decided geometrically.
                yv = apexv;
                dydx.setZero();
                region = 3;
            }
            else
            {
                yv = yline;
                Eigen::Matrix<double, 3, 2> DMi = DM * ADMi;
                Eigen::Matrix3d proj = DMi * Am;
                dydx = I3 - proj;
                region = use1 ? 1 : 2;
            }
        }

        // ---- back to 6D --------------------------------------------------
        Eigen::Matrix3d Sret = Q * yv.asDiagonal() * Q.transpose();
        VoigtVector sigma_ret(Sret(0, 0), Sret(1, 1), Sret(2, 2),
                              Sret(0, 1), Sret(1, 2), Sret(0, 2));
        for (int i = 0; i < 6; ++i)
        {
            if (!(sigma_ret(i) == sigma_ret(i)))
            {
                opserr << "ASDPlasticMaterial3D::Closest_Point (tag " << ASDP_TAG
                       << ") - NaN in the principal-space Mohr-Coulomb return"
                       << " -- rejecting step" << endln;
                return LADRUNO_MATERIAL_REFUSED;
            }
        }

        // Admissibility against the header's OWN composite f.  For plain
        // Mohr-Coulomb this is a self-check (the return is exact by
        // construction); for MohrCoulombTensionCutoff it is the real gate: this
        // path is only reached after `special_return` declined, and the MC return
        // it performs must still respect the cutoff plane.
        {
            const double f_ret = yf(sigma_ret, iv_storage, parameters_storage);
            if (!(f_ret <= tol_f))
            {
                opserr << "ASDPlasticMaterial3D::Closest_Point (tag " << ASDP_TAG
                       << ") - the principal-space return landed OUTSIDE this yield"
                       << " function's own surface (f = " << f_ret << " > tol = "
                       << tol_f << ", region " << region
                       << ") -- rejecting step (ADR-97 P2)" << endln;
                return LADRUNO_MATERIAL_REFUSED;
            }
        }

        // ---- the Koiter tangent, rotated to the global Voigt frame -------
        VoigtMatrix Tp = VoigtMatrix::Zero();
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) Tp(i, j) = dydx(i, j);
        {
            double xmax = 0.0;
            for (int i = 0; i < 3; ++i)
            {
                const double a = (xv(i) < 0) ? -xv(i) : xv(i);
                if (a > xmax) xmax = a;
            }
            double sref = yf.strength_scale(iv_storage, parameters_storage);
            if (sref < 0) sref = -sref;
            const double eps_deg = 1e-9 * ((xmax > sref) ? xmax : sref);
            static const int SIJ[3][2] = {{0, 1}, {1, 2}, {0, 2}};   // slots 3,4,5
            for (int s = 0; s < 3; ++s)
            {
                const int i = SIJ[s][0], j = SIJ[s][1];
                const double dx = xv(i) - xv(j);
                const double adx = (dx < 0) ? -dx : dx;
                Tp(3 + s, 3 + s) = (adx > eps_deg)
                                   ? (yv(i) - yv(j)) / dx
                                   : (dydx(i, i) - dydx(i, j));
            }
        }

        VoigtMatrix Rs = VoigtMatrix::Zero();
        VoigtMatrix Rsi = VoigtMatrix::Zero();
        {
            static const int BI[6][2] = {{0, 0}, {1, 1}, {2, 2}, {0, 1}, {1, 2}, {0, 2}};
            for (int k = 0; k < 6; ++k)
            {
                Eigen::Matrix3d Ek = Eigen::Matrix3d::Zero();
                Ek(BI[k][0], BI[k][1]) = 1.0;
                Ek(BI[k][1], BI[k][0]) = 1.0;
                Eigen::Matrix3d F  = Q * Ek * Q.transpose();
                Eigen::Matrix3d Fi = Q.transpose() * Ek * Q;
                for (int r = 0; r < 6; ++r)
                {
                    Rs(r, k)  = F(BI[r][0], BI[r][1]);
                    Rsi(r, k) = Fi(BI[r][0], BI[r][1]);
                }
            }
        }
        VoigtMatrix RT   = Rs * Tp;
        VoigtMatrix RTR  = RT * Rsi;
        VoigtMatrix Calg = RTR * Eelastic;

        // ---- plastic strain increment ------------------------------------
        // d eps^p = E^{-1} (sigma_tr - sigma_ret): convention-safe (it never
        // touches the engineering-shear factor of the principal flow vectors)
        // and exact for every region including the vertex.
        VoigtVector dep = VoigtVector::Zero();
        {
            Eigen::Matrix<double, 6, 6> Ee;
            for (int i = 0; i < 6; ++i)
                for (int j = 0; j < 6; ++j) Ee(i, j) = Eelastic(i, j);
            Eigen::FullPivLU< Eigen::Matrix<double, 6, 6> > elu(Ee);
            if (!elu.isInvertible())
            {
                opserr << "ASDPlasticMaterial3D::Closest_Point (tag " << ASDP_TAG
                       << ") - the elastic tangent is singular"
                       << " -- rejecting step" << endln;
                return LADRUNO_MATERIAL_REFUSED;
            }
            Eigen::Matrix<double, 6, 1> rhs6;
            for (int i = 0; i < 6; ++i) rhs6(i) = sigma_tr(i) - sigma_ret(i);
            Eigen::Matrix<double, 6, 1> d6 = elu.solve(rhs6);
            for (int i = 0; i < 6; ++i) dep(i) = d6(i);
        }

        TrialStress         = sigma_ret;
        TrialPlastic_Strain = CommitPlastic_Strain + dep;
        cp_last_iterations  = 1;      // closed form: no Newton on any region
        cp_apply_tangent_policy(Calg, Eelastic);

        if (ladruno_strict_rejects("Closest_Point (principal)", TrialStress))
            return LADRUNO_MATERIAL_REFUSED;
        return 0;
    }
'''

edit(
    "ASDPlasticMaterial3D.h",
    "    int Closest_Point(const VoigtVector & strain_incr)\n",
    PRINCIPAL + "\n    int Closest_Point(const VoigtVector & strain_incr)\n",
    "cp_principal_return + tangent policy helper",
)

# ===========================================================================
# E. dispatch inside Closest_Point, before the smooth apex block
# ===========================================================================
edit(
    "ASDPlasticMaterial3D.h",
    """            Stiffness = Eelastic;
            return 0;
        }

        if constexpr (yf_has_apex<YieldFunctionType>::value)
        {
            if (cp_apex_region(depsilon, sigma_tr, Eelastic, f_tr, tol_f))""",
    """            Stiffness = Eelastic;
            return 0;
        }

        // Ladruno (ADR-97 wp/97c): the Mohr-Coulomb family leaves here and never
        // reaches the smooth 6D Newton below -- its yield function has no usable
        // 6D gradient (Drucker-Prager substitution for |theta| >= 29 deg,
        // central differences otherwise).  MohrCoulombTensionCutoff first offers
        // the trial to ADR-84's `special_return`, whose closed-form cutoff
        // face / Rankine edge / MC-cutoff corner / compound-corner / apex returns
        // and RAW Koiter tangent are reused verbatim -- no geometry is
        // re-derived; when that hook declines (cutoff inactive at the trial, or
        // an MC-dominant trial none of whose cutoff features validated) the plain
        // Mohr-Coulomb principal return takes over and re-checks the COMPOSITE f.
        if constexpr (ladruno_cp_principal_family != 0)
        {
            if constexpr (yf_has_special_return<YieldFunctionType>::value)
            {
                VoigtVector sigma_sr, dep_sr;
                VoigtMatrix stiff_sr;
                int sr_quality = SR_QUALITY_EXACT;
                if (yf.special_return(sigma_tr, Eelastic, tol_f,
                                      iv_storage, parameters_storage,
                                      sigma_sr, dep_sr, stiff_sr, sr_quality))
                {
                    if (sr_quality == SR_QUALITY_FALLBACK
                            && INT_OPT_strict_convergence[ASDP_TAG] != 0)
                    {
                        opserr << "ASDPlasticMaterial3D::Closest_Point (tag " << ASDP_TAG
                               << ") - special_return could not land this trial on any"
                               << " exact cutoff feature (face/edge/corner/compound);"
                               << " the only remaining candidate is the conservative"
                               << " vertex projection, which would discard the"
                               << " deviatoric state -- rejecting step"
                               << " (strict_convergence)" << endln;
                        return LADRUNO_MATERIAL_REFUSED;
                    }
                    TrialStress         = sigma_sr;
                    TrialPlastic_Strain = CommitPlastic_Strain + dep_sr;
                    cp_last_iterations  = 1;
                    // ADR-84 P3: stiff_sr is the RAW active-set (Koiter) tangent,
                    // i.e. the consistent tangent of THAT exact return -- which is
                    // precisely what `Algorithmic` means here.
                    cp_apply_tangent_policy(stiff_sr, Eelastic);
                    if (ladruno_strict_rejects("Closest_Point (special_return)", TrialStress))
                        return LADRUNO_MATERIAL_REFUSED;
                    return 0;
                }
            }
            return cp_principal_return(depsilon, sigma_tr, Eelastic, tol_f);
        }

        if constexpr (yf_has_apex<YieldFunctionType>::value)
        {
            if (cp_apex_region(depsilon, sigma_tr, Eelastic, f_tr, tol_f))""",
    "Closest_Point principal dispatch",
)

# ===========================================================================
# F. refusal messages, in the material and in the parser
# ===========================================================================
edit(
    "ASDPlasticMaterial3D.h",
    """                   << " (ADR-97 D3). P1 ships VonMises and DruckerPrager with Null,"
                   << " Linear (scalar or tensor) and ArmstrongFrederick hardening;"
                   << " MohrCoulomb and MohrCoulombTensionCutoff are ADR-97 P2,"
                   << " HoekBrown P3, StiffSoil P5. Use Backward_Euler." << endln;""",
    """                   << " (ADR-97 D3). Supported today: VonMises and DruckerPrager"
                   << " with Null / Linear (scalar or tensor) / ArmstrongFrederick"
                   << " hardening (P1), and the perfectly plastic MohrCoulomb and"
                   << " MohrCoulombTensionCutoff families where BOTH the yield"
                   << " function and the plastic flow direction are of that family"
                   << " (P2). MIXED pairings such as MohrCoulomb_YF x VonMises_PF or"
                   << " VonMises_YF x MohrCoulomb_PF are NOT supported: the"
                   << " principal-space map assumes both are piecewise linear, and"
                   << " the smooth 6D map cannot use MohrCoulomb's Lode-angle"
                   << " gradient. HoekBrown is ADR-97 P3, StiffSoil P5."
                   << " Use Backward_Euler." << endln;""",
    "material refusal message P2",
)

edit(
    "OPS_AllASDPlasticMaterial3Ds.cpp",
    """                            opserr << "   ADR-97 D3: the closest-point map is added"
                                   << " family by family. P1 covers VonMises and"
                                   << " DruckerPrager with Null / Linear / "
                                   << "ArmstrongFrederick hardening; MohrCoulomb and"
                                   << " MohrCoulombTensionCutoff are P2, HoekBrown"
                                   << " P3, StiffSoil P5. Use Backward_Euler."
                                   << endln;""",
    """                            opserr << "   ADR-97 D3: the closest-point map is added"
                                   << " family by family. P1 covers VonMises and"
                                   << " DruckerPrager with Null / Linear / "
                                   << "ArmstrongFrederick hardening; P2 adds the"
                                   << " perfectly plastic MohrCoulomb and"
                                   << " MohrCoulombTensionCutoff families, but ONLY"
                                   << " where the yield function AND the plastic flow"
                                   << " direction are both of that family -- a mixed"
                                   << " pairing (MohrCoulomb_YF with VonMises_PF or"
                                   << " DruckerPrager_PF, VonMises_YF or"
                                   << " DruckerPrager_YF with MohrCoulomb_PF,"
                                   << " HoekBrown_PF with anything) is verified by no"
                                   << " oracle and stays refused. HoekBrown is P3,"
                                   << " StiffSoil P5. Use Backward_Euler."
                                   << endln;""",
    "parser refusal message P2",
)

# Ladruno (ADR-97 wp/97c): the P1 comment above the parse-time refusal block
edit(
    "OPS_AllASDPlasticMaterial3Ds.cpp",
    """                        // Ladruno (ADR-97 wp/97b): the closest-point return map.
                        // Available family by family (ADR-97 D3): P1 ships
                        // VonMises and DruckerPrager (flank + apex) with Null,
                        // Linear scalar/tensor and ArmstrongFrederick hardening.""",
    """                        // Ladruno (ADR-97 wp/97b, extended wp/97c): the
                        // closest-point return map.  Available family by family
                        // (ADR-97 D3): P1 ships VonMises and DruckerPrager
                        // (flank + apex) with Null, Linear scalar/tensor and
                        // ArmstrongFrederick hardening; P2 adds the principal-
                        // stress-space MohrCoulomb and MohrCoulombTensionCutoff
                        // returns for MATCHED YF/PF pairs only.""",
    "parser comment P2",
)


def main():
    applied, skipped = [], []
    for relpath, anchor, new, tag in EDITS:
        path = os.path.join(SRC, relpath)
        if not os.path.isfile(path):
            print("MISSING FILE: %s" % path)
            return 2
        with open(path, "r", encoding="utf-8", newline="") as fh:
            txt = fh.read()
        # the sources in this tree are CRLF; convert the LF-written anchors.
        crlf = "\r\n" in txt
        an = anchor.replace("\n", "\r\n") if crlf else anchor
        nw = new.replace("\n", "\r\n") if crlf else new
        if nw in txt:
            skipped.append(tag)
            continue
        n = txt.count(an)
        if n != 1:
            print("ANCHOR MISS (%d matches) in %s for %r" % (n, relpath, tag))
            return 3
        txt = txt.replace(an, nw, 1)
        with open(path, "w", encoding="utf-8", newline="") as fh:
            fh.write(txt)
        applied.append(tag)
    print("applied  : %d" % len(applied))
    for t in applied:
        print("   + %s" % t)
    print("already  : %d" % len(skipped))
    for t in skipped:
        print("   = %s" % t)
    return 0


if __name__ == "__main__":
    sys.exit(main())
