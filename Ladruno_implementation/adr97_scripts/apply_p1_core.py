"""ADR-97 WP-b (P1) -- re-runnable C++ edit script, part 2: the integrator.

Adds `Closest_Point` (fully implicit closest-point return map) and the
`Algorithmic` consistent tangent to ASDPlasticMaterial3D.h, and the parser
tokens + refusals to OPS_AllASDPlasticMaterial3Ds.cpp.

Run apply_p1_cpp.py FIRST (it adds the interface members this uses).

usage:  python3.12 apply_p1_core.py [<worktree root>]
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from apply_p1_cpp import apply_edits            # noqa: E402

ROOT = sys.argv[1] if len(sys.argv) > 1 else \
    r"C:\Users\nmb\Documents\Github\OpenSees\.claude\worktrees\asdplastic-review-plan-62585c"
ASD = os.path.join(ROOT, "SRC", "material", "nD", "ASDPlasticMaterial3D")

E = []


def edit(rel, marker, anchor, replacement, last=False):
    E.append((rel, marker, anchor, replacement, last))


# NOTE: the CRTP two-phase-lookup fix (`dm_dsigma_buffer` lives in the DEPENDENT
# base PlasticFlowBase<T>, so an unqualified reference from a derived class is
# ill-formed on GCC and accepted by MSVC -- the ADR-94 wave's red CI) is a plain
# textual substitution applied in main() below.

# ---------------------------------------------------------------------------
# 1. compile-time fold over the IV tuple + the CP support flag
# ---------------------------------------------------------------------------
edit(
    "ASDPlasticMaterial3D.h",
    "asdp_all_ivs_support_cp",
    """template <
    class ElasticityType,
    class YieldFunctionType,
    class PlasticFlowType,
    int thisClassTag >
class ASDPlasticMaterial3D : public NDMaterial
{""",
    """// Ladruno (ADR-97 wp/97b): does EVERY internal variable of this specialization
// carry a hardening law with analytic closest-point derivatives (dh/dq, dh/dm)?
// Folded over the concatenated IV tuple at compile time so the parser can refuse
// `integration_method Closest_Point` for a specialization this ADR has not
// converted, instead of silently believing the inert zero defaults.
template <typename Tuple>
struct asdp_all_ivs_support_cp;

template <typename... Ts>
struct asdp_all_ivs_support_cp<std::tuple<Ts...>>
{
    static constexpr bool value = (Ts::hardening_supports_cp() && ... && true);
};

template <
    class ElasticityType,
    class YieldFunctionType,
    class PlasticFlowType,
    int thisClassTag >
class ASDPlasticMaterial3D : public NDMaterial
{""")

edit(
    "ASDPlasticMaterial3D.h",
    "ladruno_cp_supported",
    """    using parameters_storage_t = utuple_storage<parameters_concat_types>;
""",
    """    using parameters_storage_t = utuple_storage<parameters_concat_types>;

    // Ladruno (ADR-97 wp/97b): is `integration_method Closest_Point` available for
    // THIS specialization?  All three legs must opt in: the yield function must
    // supply the uncontracted df/dq, the plastic flow direction the analytic
    // dm/dsigma and dm/dq, and every internal variable's hardening law dh/dq and
    // dh/dm.  P1 ships VonMises + DruckerPrager x {Null, Linear scalar/tensor,
    // ArmstrongFrederick} = 20 of the 46 registered specializations; the rest are
    // refused at parse time naming the ADR phase that will deliver them.
    static constexpr bool ladruno_cp_supported =
        yf_has_cp_derivatives<YieldFunctionType>::value &&
        pf_has_cp_derivatives<PlasticFlowType>::value &&
        asdp_all_ivs_support_cp<iv_concat_types>::value;

    static constexpr bool supportsClosestPoint() { return ladruno_cp_supported; }

    // Ladruno (ADR-97 wp/97b): Newton system size cap.  6 stress rows + the
    // internal variables + one consistency row; the widest registered
    // specialization has 14 IV components, so 26 leaves headroom and keeps the
    // Jacobian on the stack (Eigen fixed-max dynamic storage, no heap traffic in
    // the Gauss-point loop).  A specialization that exceeds it is REFUSED, loudly.
    static constexpr int ASDP_CP_MAXN = 26;
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
                          Eigen::ColMajor, ASDP_CP_MAXN, ASDP_CP_MAXN> cp_matrix_t;
    typedef Eigen::Matrix<double, Eigen::Dynamic, 1,
                          Eigen::ColMajor, ASDP_CP_MAXN, 1> cp_vector_t;
""")

# ---------------------------------------------------------------------------
# 2. dispatch case
# ---------------------------------------------------------------------------
edit(
    "ASDPlasticMaterial3D.h",
    "exitflag = this->Closest_Point(strain_increment);",
    """        case ASDPlasticMaterial3D_Constitutive_Integration_Method::Backward_Euler :
            exitflag = this->Backward_Euler(strain_increment);;
            break;""",
    """        case ASDPlasticMaterial3D_Constitutive_Integration_Method::Backward_Euler :
            exitflag = this->Backward_Euler(strain_increment);;
            break;
        // Ladruno (ADR-97 wp/97b): the fully implicit closest-point return map.
        // Backward_Euler above is untouched (ADR-97 D1: byte-identical).
        case ASDPlasticMaterial3D_Constitutive_Integration_Method::Closest_Point :
            exitflag = this->Closest_Point(strain_increment);
            break;""")

# ---------------------------------------------------------------------------
# 3. accept list
# ---------------------------------------------------------------------------
edit(
    "ASDPlasticMaterial3D.h",
    "::Closest_Point   // Ladruno (ADR-97",
    """                || method == (int) ASDPlasticMaterial3D_Constitutive_Integration_Method::Backward_Euler_LineSearch)""",
    """                || method == (int) ASDPlasticMaterial3D_Constitutive_Integration_Method::Closest_Point   // Ladruno (ADR-97 wp/97b)
                || method == (int) ASDPlasticMaterial3D_Constitutive_Integration_Method::Backward_Euler_LineSearch)""")

# ---------------------------------------------------------------------------
# 4. cp_iterations getter + response
# ---------------------------------------------------------------------------
edit(
    "ASDPlasticMaterial3D.h",
    "getCPIterations",
    """    const Vector &getInternalVariableByPos(int pos)
    {""",
    """    // Ladruno (ADR-97 wp/97b): how many Newton iterations the LAST Closest_Point
    // step took at this Gauss point (0 for an elastic step).  Exposed so ADR-97's
    // gate 1 can assert the quadratic <= 5 without parsing stdout -- the material
    // prints a page per construction and pytest's capfd cannot see the .pyd's
    // cout anyway (LEDGER_quirks).
    const Vector &getCPIterations(void)
    {
        static Vector result(1);
        result(0) = (double) cp_last_iterations;
        return result;
    }

    const Vector &getInternalVariableByPos(int pos)
    {""")

edit(
    "ASDPlasticMaterial3D.h",
    '"cp_iterations"',
    """        else
        {
            const char *iv_name = argv[0];""",
    """        else if (strcmp(argv[0], "cp_iterations") == 0 ) {
            // Ladruno (ADR-97 wp/97b)
            s.tag("ResponseType", "cp_iterations");
            return new MaterialResponse(this, 9, this->getCPIterations());
        }
        else
        {
            const char *iv_name = argv[0];""")

edit(
    "ASDPlasticMaterial3D.h",
    "*(matInformation.theVector) = getCPIterations();",
    """        else if (responseID == 8)
            *(matInformation.theVector) = getJ2Strain();""",
    """        else if (responseID == 8)
            *(matInformation.theVector) = getJ2Strain();
        else if (responseID == 9)   // Ladruno (ADR-97 wp/97b)
            *(matInformation.theVector) = getCPIterations();""")

# ---------------------------------------------------------------------------
# 5. the member
# ---------------------------------------------------------------------------
edit(
    "ASDPlasticMaterial3D.h",
    "int cp_last_iterations = 0;",
    """    bool first_step;
    bool stress_set_externally;""",
    """    bool first_step;
    bool stress_set_externally;

    // Ladruno (ADR-97 wp/97b): per-instance Closest_Point Newton iteration count
    // (NSDMI, so both constructors get it without touching their bodies).
    int cp_last_iterations = 0;""")

# ---------------------------------------------------------------------------
# 6. THE INTEGRATOR
# ---------------------------------------------------------------------------
CP_CODE = r'''
    //==================================================================================================
    // Ladruno (ADR-97 wp/97b): CLOSEST-POINT return map (CPPM) + consistent tangent
    //==================================================================================================
    //
    // Solves, as ONE coupled Newton system in x = (sigma_{n+1}, q_{n+1}, dlambda):
    //
    //     R_sigma = s - s_tr + dl * E(s) * m(s, q)                     (6 rows)
    //     R_q     = q - q_n  - dl * h(s, q, m(s,q))                    (n_q rows)
    //     R_f     = f(s, q)                                            (1 row)
    //
    // with s_tr = sigma_n + E(sigma_n)*depsilon.  This is NOT what
    // `Backward_Euler` does: that is an Ortiz-Simo CUTTING PLANE which
    // re-evaluates n, m and H at the running iterate and accumulates
    // sigma -= deltaLambda*(E*m), so its fixed point is
    // s_tr - sum_k dl_k E m(sigma^k), and its internal-variable update is a
    // Newton-path quadrature exact only for constant h.  The two coincide exactly
    // when the flow direction does not rotate over the step (proportional loading)
    // and diverge otherwise -- measurably so for Armstrong-Frederick, which 22 of
    // the 46 registered specializations carry (ADR-94 M3; ADR-97 P0 measured a
    // 9.15 stress spread over four equally valid cutting-plane iterate paths on a
    // stress of ~46, against 1.1e-12 for this map).
    //
    // `Backward_Euler` is untouched (ADR-97 D1).  Every result the fork has
    // published on this material was obtained on it.
    //
    // Conventions (ADR-94 wp/94c, do not re-derive): tension positive,
    // p = trace/3, Voigt [11 22 33 12 23 13], strain-like slots ENGINEERING, every
    // YF/PF gradient taken w.r.t. the STORED slot, every contraction a plain dot.
    //
    // E(sigma_{n+1}) is evaluated INSIDE the residual (ADR-97 D6), which makes the
    // converged state hyperelastically consistent; the dE/dsigma JACOBIAN block is
    // optional (identically zero for LinearIsotropic3D_EL, the only elasticity
    // registered in all 43 non-StiffSoil specializations) -- dropping it costs
    // iterations, never accuracy, because the converged x still satisfies the
    // exact residual.
    //
    // There is NO intersection / elastic-fraction split, and none is missing:
    // `Backward_Euler` zeroes `intersection_stress`/`intersection_strain` at setup
    // and never reads either again (the yield-crossing bisection lives in
    // Forward_Euler and compute_local_stress only).
    //
    // No line search: for convex f with associated flow the CPPM residual is the
    // stationarity system of a strictly convex closest-point projection in the
    // E-metric, so Newton from the elastic predictor is globally convergent on the
    // yielding branch.  A line search that "succeeds" on a non-converged system is
    // ADR-94 M7's Backward_Euler_LineSearch failure mode (2/20 vs BE's 20/20).

    int cp_n_iv()
    {
        int n = 0;
        iv_storage.apply([&n](auto & internal_variable)
        {
            n += internal_variable.size();
        });
        return n;
    }

    // Scatter x into (TrialStress, IV trial values) and form the residual; when
    // `J` is non-null, form the Jacobian at the SAME point.  Returns false on a
    // non-finite residual.
    bool cp_assemble(const cp_vector_t& x, int N,
                     const VoigtVector& sigma_tr, const VoigtVector& depsilon,
                     VoigtVector& m_out, cp_vector_t& R, cp_matrix_t* J)
    {
        using namespace ASDPlasticMaterial3DGlobals;

        const double dl = x(N - 1);
        for (int i = 0; i < 6; ++i) TrialStress(i) = x(i);
        {
            int off = 6;
            iv_storage.apply([&](auto & internal_variable)
            {
                const int nq = internal_variable.size();
                for (int i = 0; i < nq; ++i) internal_variable.trial_value(i) = x(off + i);
                off += nq;
            });
        }

        // GCC: never bind an Eigen product (or a functor's mutable return buffer,
        // which the NEXT call overwrites) to a reference -- named value locals only.
        VoigtMatrix Ecur = et(TrialStress, parameters_storage);          // ADR-97 D6
        VoigtVector m    = pf(depsilon, TrialStress, iv_storage, parameters_storage);
        m_out = m;
        VoigtVector Em   = Ecur * m;
        const double f_val = yf(TrialStress, iv_storage, parameters_storage);

        for (int i = 0; i < 6; ++i)
            R(i) = x(i) - sigma_tr(i) + dl * Em(i);
        {
            int off = 6;
            iv_storage.apply([&](auto & internal_variable)
            {
                const int nq = internal_variable.size();
                // hardening_function reads the TRIAL value, i.e. q_{n+1}: this is
                // what makes the internal-variable update implicit (ADR-97 D4).
                auto h = internal_variable.hardening_function(depsilon, m, TrialStress, parameters_storage);
                for (int i = 0; i < nq; ++i)
                    R(off + i) = x(off + i) - internal_variable.committed_value(i) - dl * h(i);
                off += nq;
            });
        }
        R(N - 1) = f_val;

        for (int i = 0; i < N; ++i)
            if (!(R(i) == R(i))) return false;      // NaN

        if (J == 0) return true;

        VoigtVector nvec = yf.df_dsigma_ij(TrialStress, iv_storage, parameters_storage);
        VoigtMatrix dmds = pf.dm_dsigma(depsilon, TrialStress, iv_storage, parameters_storage);
        VoigtMatrix Edmds = Ecur * dmds;

        cp_matrix_t& Jm = *J;
        Jm.setZero();

        //   J_ss = I + dl * E * dm/ds     J_sl = E*m     J_fs = (df/ds)^T
        for (int i = 0; i < 6; ++i)
        {
            for (int j = 0; j < 6; ++j)
                Jm(i, j) = ((i == j) ? 1.0 : 0.0) + dl * Edmds(i, j);
            Jm(i, N - 1) = Em(i);
            Jm(N - 1, i) = nvec(i);
        }
        //   + dl * (dE/ds : m)  -- ADR-97 D6; identically zero for LinearIsotropic3D_EL
        if constexpr (el_is_stress_dependent<ElasticityType>::value)
        {
            VoigtMatrix dEm = VoigtMatrix::Zero();
            et.dE_dsigma_contract(TrialStress, m, parameters_storage, dEm);
            for (int i = 0; i < 6; ++i)
                for (int j = 0; j < 6; ++j) Jm(i, j) += dl * dEm(i, j);
        }

        //   J_sq = dl * E * dm/dq        J_fq = (df/dq)^T
        {
            int offb = 6;
            iv_storage.apply([&](auto & iv_b)
            {
                const int nb = iv_b.size();
                VoigtMatrix dmdq = VoigtMatrix::Zero();
                pf.dm_dq(iv_b, depsilon, TrialStress, iv_storage, parameters_storage, dmdq);
                VoigtMatrix Edmdq = Ecur * dmdq;
                double dfdq[6] = {0., 0., 0., 0., 0., 0.};
                yf.df_dq(iv_b, TrialStress, iv_storage, parameters_storage, dfdq);
                for (int c = 0; c < nb; ++c)
                {
                    for (int i = 0; i < 6; ++i) Jm(i, offb + c) = dl * Edmdq(i, c);
                    Jm(N - 1, offb + c) = dfdq[c];
                }
                offb += nb;
            });
        }

        //   J_ql = -h
        //   J_qs = -dl * (dh/dm)(dm/ds)                 [dh/ds is zero by construction:
        //                                                no policy in this tree reads sigma]
        //   J_qq = I - dl * (dh/dq + (dh/dm)(dm/dq))    [the OFF-DIAGONAL IV blocks are
        //                                                NOT zero: every h reads m, and m
        //                                                reads the flow direction's own
        //                                                back stress]
        {
            int offa = 6;
            iv_storage.apply([&](auto & iv_a)
            {
                const int na = iv_a.size();
                auto h = iv_a.hardening_function(depsilon, m, TrialStress, parameters_storage);
                VoigtMatrix dhdm = VoigtMatrix::Zero();
                VoigtMatrix dhdq = VoigtMatrix::Zero();
                iv_a.hardening_dh_dm(depsilon, m, TrialStress, parameters_storage, dhdm);
                iv_a.hardening_dh_dq(depsilon, m, TrialStress, parameters_storage, dhdq);
                for (int rr = 0; rr < na; ++rr)
                {
                    Jm(offa + rr, N - 1) = -h(rr);
                    for (int j = 0; j < 6; ++j)
                    {
                        double v = 0.0;
                        for (int k = 0; k < 6; ++k) v += dhdm(rr, k) * dmds(k, j);
                        Jm(offa + rr, j) = -dl * v;
                    }
                }
                int offb = 6;
                iv_storage.apply([&](auto & iv_b)
                {
                    const int nb = iv_b.size();
                    VoigtMatrix dmdq_b = VoigtMatrix::Zero();
                    pf.dm_dq(iv_b, depsilon, TrialStress, iv_storage, parameters_storage, dmdq_b);
                    const bool same_iv = (offa == offb);
                    for (int rr = 0; rr < na; ++rr)
                        for (int c = 0; c < nb; ++c)
                        {
                            double v = same_iv ? dhdq(rr, c) : 0.0;
                            for (int k = 0; k < 6; ++k) v += dhdm(rr, k) * dmdq_b(k, c);
                            Jm(offa + rr, offb + c) =
                                ((same_iv && rr == c) ? 1.0 : 0.0) - dl * v;
                        }
                    offb += nb;
                });
                offa += na;
            });
        }
        return true;
    }

    // ELASTIC-METRIC apex classification.  `check_apex_region` in the yield
    // function is EUCLIDEAN (p - p_apex >= eta*q) and says so in its own comment:
    // the exact condition needs K, G and the dilatancy, none of which the YF's
    // signature can see.  `Closest_Point` classifies HERE, where E is in scope, and
    // does not call `check_apex_region` at all.  ADR-97's P0 oracle quantified the
    // cost of the Euclidean test in BOTH directions: at etabar = 0.2 (exact slope
    // 0.333 < the header's 0.4) a trial at (p-p_apex)/q = 0.36 is classified CONE
    // and the cone return then gives sqrt(J2)_{n+1} = -0.47, an INADMISSIBLE
    // negative deviatoric norm; at etabar = eta = 0.4 (exact slope 0.667) trials at
    // 0.45 and 0.60 are classified APEX although the correct return is to the cone.
    //
    // The test itself is family-agnostic: take the LINEARISED cone step at the
    // trial state, dl0 = f_tr/(n:E:m - H), and ask whether the deviatoric part of
    // the returned stress has FLIPPED sign.  For Drucker-Prager E*m has the
    // deviatoric part (G/q)*r exactly, so this reduces to q_tr - G*dl0 < 0 -- the
    // oracle's exact test, with K*etabar and G entering through E.
    bool cp_apex_region(const VoigtVector& depsilon, const VoigtVector& sigma_tr,
                        const VoigtMatrix& Eelastic, double f_tr)
    {
        using namespace ASDPlasticMaterial3DGlobals;
        VoigtVector m0 = pf(depsilon, sigma_tr, iv_storage, parameters_storage);
        VoigtVector n0 = yf.df_dsigma_ij(sigma_tr, iv_storage, parameters_storage);
        VoigtVector Em0 = Eelastic * m0;
        // The plastic modulus formed from CLOSEST_POINT's own df/dq, not from
        // yf.hardening(): for Drucker-Prager the two differ, because the shipped
        // yf_hardening carries a df/dk = -1 term for a cohesion IV that its own f
        // does not contain (ADR-97 P0 header finding 2).
        double H0 = 0.0;
        iv_storage.apply([&](auto & internal_variable)
        {
            double d[6] = {0., 0., 0., 0., 0., 0.};
            yf.df_dq(internal_variable, sigma_tr, iv_storage, parameters_storage, d);
            auto h = internal_variable.hardening_function(depsilon, m0, sigma_tr, parameters_storage);
            const int nq = internal_variable.size();
            for (int i = 0; i < nq; ++i) H0 += d[i] * h(i);
        });
        const double den = n0.dot(Em0) - H0;
        if (!(den > MACHINE_EPSILON)) return false;
        const double dl0 = f_tr / den;
        VoigtVector dev_tr  = sigma_tr.deviator();
        VoigtVector dev_Em0 = Em0.deviator();
        VoigtVector dev_ret = dev_tr - dl0 * dev_Em0;
        return tensor_dot_stress_like(dev_ret, dev_tr) < 0.0;
    }

    // Reduced apex system.  At the vertex dm/ds blows up as 1/sqrt(J2) and the
    // plastic flow direction is a whole subdifferential cone, so the 6-row
    // residual cannot produce it: sigma_{n+1} is PINNED at the yield function's
    // own apex, and the only unknowns are the internal variables and dlambda.
    //     sigma_{n+1} = sigma_apex(q_{n+1})
    //     d eps^p     = E^{-1}(sigma_tr - sigma_apex)
    //     dl          = <m_apex, d eps^p>_e / <m_apex, m_apex>_e   (>= 0)
    //     q_{n+1}     = q_n + dl * h(sigma_apex, q_{n+1}, m_apex)
    // solved by fixed point (ONE pass for perfect plasticity, where h == 0).
    // Returns false when the yield function's own apex is not on its own surface,
    // in which case the caller falls through to the generic Newton -- the same
    // admissibility gate `Backward_Euler` uses, kept verbatim.
    bool cp_apex_return(const VoigtVector& depsilon, const VoigtVector& sigma_tr,
                        const VoigtMatrix& Eelastic, double tol_f, double tol_s,
                        int max_iter, int& rc)
    {
        using namespace ASDPlasticMaterial3DGlobals;

        Eigen::Matrix<double, 6, 6> E_eig;
        for (int i = 0; i < 6; ++i)
            for (int j = 0; j < 6; ++j) E_eig(i, j) = Eelastic(i, j);
        Eigen::FullPivLU< Eigen::Matrix<double, 6, 6> > elu(E_eig);
        if (!elu.isInvertible())
        {
            opserr << "ASDPlasticMaterial3D::Closest_Point (tag " << ASDP_TAG
                   << ") - the elastic tangent is singular at the apex return"
                   << " -- rejecting step" << endln;
            rc = LADRUNO_MATERIAL_REFUSED;
            return true;
        }

        VoigtVector sigma_apex = VoigtVector::Zero();
        VoigtVector m_apex     = VoigtVector::Zero();
        VoigtVector dep        = VoigtVector::Zero();
        double dl = 0.0;
        bool ok = false;

        for (int it = 0; it < max_iter; ++it)
        {
            // Copy out of the YF's return buffer immediately: apex_stress and
            // df_dsigma_ij share one mutable member per functor.
            VoigtVector sa = yf.apex_stress(iv_storage, parameters_storage);
            sigma_apex = sa;
            const double f_apex = yf(sigma_apex, iv_storage, parameters_storage);
            const double f_apex_abs = (f_apex < 0) ? -f_apex : f_apex;
            if (!(f_apex_abs <= tol_f))
            {
                // The apex this YF names is not on its own surface: fall through
                // to the generic map rather than commit a fabricated stress.
                iv_storage.revert_all();
                return false;
            }
            VoigtVector ma = pf(depsilon, sigma_apex, iv_storage, parameters_storage);
            m_apex = ma;

            Eigen::Matrix<double, 6, 1> rhs_e;
            for (int i = 0; i < 6; ++i) rhs_e(i) = sigma_tr(i) - sigma_apex(i);
            Eigen::Matrix<double, 6, 1> dep_e = elu.solve(rhs_e);
            for (int i = 0; i < 6; ++i) dep(i) = dep_e(i);

            // Both operands are ENGINEERING-strain-like Voigt vectors, so the
            // projection uses the engineering contraction (ADR-94 wp/94c B5).
            const double mm = tensor_dot_engineering_strain_like(m_apex, m_apex);
            double dl_new = 0.0;
            if (mm > MACHINE_EPSILON)
            {
                dl_new = tensor_dot_engineering_strain_like(m_apex, dep) / mm;
                if (dl_new < 0.0) dl_new = 0.0;
            }

            double dq_max = 0.0;
            iv_storage.apply([&](auto & internal_variable)
            {
                auto h = internal_variable.hardening_function(depsilon, m_apex, sigma_apex, parameters_storage);
                const int nq = internal_variable.size();
                for (int i = 0; i < nq; ++i)
                {
                    const double nv = internal_variable.committed_value(i) + dl_new * h(i);
                    const double d = nv - internal_variable.trial_value(i);
                    const double ad = (d < 0) ? -d : d;
                    if (ad > dq_max) dq_max = ad;
                    internal_variable.trial_value(i) = nv;
                }
            });
            dl = dl_new;
            cp_last_iterations = it + 1;
            if (dq_max <= tol_s) { ok = true; break; }
        }

        if (!ok)
        {
            opserr << "ASDPlasticMaterial3D::Closest_Point (tag " << ASDP_TAG
                   << ") - the apex internal-variable fixed point did not converge in "
                   << max_iter << " iterations -- rejecting step" << endln;
            rc = LADRUNO_MATERIAL_REFUSED;
            return true;
        }

        {
            VoigtVector sa = yf.apex_stress(iv_storage, parameters_storage);
            sigma_apex = sa;
            Eigen::Matrix<double, 6, 1> rhs_e;
            for (int i = 0; i < 6; ++i) rhs_e(i) = sigma_tr(i) - sigma_apex(i);
            Eigen::Matrix<double, 6, 1> dep_e = elu.solve(rhs_e);
            for (int i = 0; i < 6; ++i) dep(i) = dep_e(i);
        }
        TrialStress         = sigma_apex;
        TrialPlastic_Strain = CommitPlastic_Strain + dep;
        (void) dl;

        // Algorithmic apex tangent.  sigma_{n+1} = sigma_apex(q_{n+1}), so
        // C_alg = (d sigma_apex/dq)(dq/d eps).  For every yield function this ADR
        // ships, apex_stress() is built from MODEL PARAMETERS only
        // (DruckerPrager_YF: p_apex = xi_c/eta), so d sigma_apex/dq is identically
        // zero and the consistent tangent is EXACTLY the zero matrix -- rank 0,
        // which is what the P0 oracle pins (cppm_dp.py, "apex, associated,
        // perfect": rank C = 0, FD error 0).  That is rank deficient by
        // construction and WILL make an element whose every Gauss point sits at
        // the apex singular, which is the true state of affairs.  The finite
        // difference below MEASURES d sigma_apex/dq instead of assuming it; a
        // yield function whose apex moves with an internal variable gets a
        // one-time warning and the same zero, which ADR-97 P5 replaces.
        VoigtMatrix apex_stiff = VoigtMatrix::Zero();
        {
            double a_max = 0.0;
            VoigtVector s0 = yf.apex_stress(iv_storage, parameters_storage);
            iv_storage.apply([&](auto & internal_variable)
            {
                const int nq = internal_variable.size();
                for (int i = 0; i < nq; ++i)
                {
                    const double v0 = internal_variable.trial_value(i);
                    const double av = (v0 < 0) ? -v0 : v0;
                    const double dv = 1e-8 * (av + 1.0);
                    internal_variable.trial_value(i) = v0 + dv;
                    VoigtVector sp = yf.apex_stress(iv_storage, parameters_storage);
                    internal_variable.trial_value(i) = v0;
                    for (int r = 0; r < 6; ++r)
                    {
                        const double d = (sp(r) - s0(r)) / dv;
                        const double ad = (d < 0) ? -d : d;
                        if (ad > a_max) a_max = ad;
                    }
                }
            });
            if (a_max > 0.0
                && INT_OPT_tangent_operator_type[ASDP_TAG]
                   == ASDPlasticMaterial3D_Tangent_Operator_Type::Algorithmic)
            {
                static bool warned_apex_tangent = false;
                if (!warned_apex_tangent)
                {
                    opserr << "ASDPlasticMaterial3D::Closest_Point (tag " << ASDP_TAG
                           << ") - this yield function's apex MOVES with an internal"
                           << " variable (|d sigma_apex/dq| = " << a_max << "), so the"
                           << " zero apex tangent reported under tangent_type"
                           << " Algorithmic is approximate (ADR-97 P5)." << endln;
                    warned_apex_tangent = true;
                }
            }
        }

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
                Stiffness = apex_stiff;
                break;
            case TOT::Secant:
            default:
                Stiffness = VoigtMatrix((apex_stiff + Eelastic) / 2.0);
                break;
            }
        }

        if (ladruno_strict_rejects("Closest_Point (apex)", TrialStress))
        {
            rc = LADRUNO_MATERIAL_REFUSED;
            return true;
        }
        rc = 0;
        return true;
    }

    int Closest_Point(const VoigtVector & strain_incr)
    {
        using namespace ASDPlasticMaterial3DGlobals;

        if (!ladruno_cp_supported)
        {
            opserr << "ASDPlasticMaterial3D::Closest_Point (tag " << ASDP_TAG
                   << ") - integration_method Closest_Point is not implemented for"
                   << " this yield function / plastic flow / hardening family"
                   << " (ADR-97 D3). P1 ships VonMises and DruckerPrager with Null,"
                   << " Linear (scalar or tensor) and ArmstrongFrederick hardening;"
                   << " MohrCoulomb and MohrCoulombTensionCutoff are ADR-97 P2,"
                   << " HoekBrown P3, StiffSoil P5. Use Backward_Euler." << endln;
            return LADRUNO_MATERIAL_REFUSED;
        }

        const int    max_iter = INT_OPT_n_max_iterations[ASDP_TAG];
        const double tol_f    = yf_tolerance();

        VoigtVector depsilon = strain_incr;
        const VoigtVector& sigma_n = CommitStress;

        iv_storage.revert_all();
        dsigma.setZero();
        // The non-Algorithmic tangent types read `depsilon_elpl` as the `depsilon`
        // argument of pf()/hardening(); Backward_Euler leaves it stale (it never
        // sets it), Closest_Point does not.
        depsilon_elpl = depsilon;

        VoigtMatrix Eelastic = et(sigma_n, parameters_storage);
        VoigtVector sigma_tr = sigma_n + Eelastic * depsilon;

        TrialStrain         = CommitStrain + depsilon;
        TrialStress         = sigma_tr;
        TrialPlastic_Strain = CommitPlastic_Strain;
        cp_last_iterations  = 0;

        // The stress- and IV-row tolerances get ADR-94 M5's relative floor too, so
        // the same problem in Pa and in kPa converges identically.  Every internal
        // variable in the P1 families is stress dimensioned (a yield stress, a
        // cohesion, a back stress), so one scaled tolerance serves all of them.
        double scale = yf.strength_scale(iv_storage, parameters_storage);
        if (scale < 0) scale = -scale;
        double tol_s = DBL_OPT_stress_absolute_tol[ASDP_TAG];
        {
            const double rel_s = DBL_OPT_f_relative_tol[ASDP_TAG] * scale;
            if (rel_s > tol_s) tol_s = rel_s;
        }

        const double f_tr = yf(sigma_tr, iv_storage, parameters_storage);
        if (!(f_tr == f_tr))
        {
            opserr << "ASDPlasticMaterial3D::Closest_Point (tag " << ASDP_TAG
                   << ") - the elastic predictor's yield value is NaN"
                   << " -- rejecting step" << endln;
            return LADRUNO_MATERIAL_REFUSED;
        }
        if (f_tr <= tol_f)
        {
            // Genuinely elastic.  Note this is NOT Backward_Euler's second
            // disjunct ("f merely DECREASED by more than tol => call it elastic"),
            // which accepts an end state that is still outside the surface; ADR-84
            // P2a had to gate that with strict_convergence.  A new integrator has
            // no compatibility reason to inherit it.
            Stiffness = Eelastic;
            return 0;
        }

        if constexpr (yf_has_apex<YieldFunctionType>::value)
        {
            if (cp_apex_region(depsilon, sigma_tr, Eelastic, f_tr))
            {
                int rc = -1;
                if (cp_apex_return(depsilon, sigma_tr, Eelastic, tol_f, tol_s, max_iter, rc))
                    return rc;
                // else: the YF's apex failed its own admissibility gate; fall
                // through to the generic Newton below.
                TrialStress = sigma_tr;
            }
        }

        const int n_iv = cp_n_iv();
        const int N = 6 + n_iv + 1;
        if (N > ASDP_CP_MAXN)
        {
            opserr << "ASDPlasticMaterial3D::Closest_Point (tag " << ASDP_TAG
                   << ") - this specialization needs a " << N << "x" << N
                   << " Newton system, over the ASDP_CP_MAXN = " << ASDP_CP_MAXN
                   << " cap -- rejecting step (raise the cap in "
                   << "ASDPlasticMaterial3D.h)" << endln;
            return LADRUNO_MATERIAL_REFUSED;
        }

        cp_vector_t x(N), R(N), dx(N);
        cp_matrix_t J(N, N);

        // Elastic-predictor start: s = s_tr, q = q_n, dl = 0.
        for (int i = 0; i < 6; ++i) x(i) = sigma_tr(i);
        {
            int off = 6;
            iv_storage.apply([&](auto & internal_variable)
            {
                const int nq = internal_variable.size();
                for (int i = 0; i < nq; ++i) x(off + i) = internal_variable.trial_value(i);
                off += nq;
            });
        }
        x(N - 1) = 0.0;

        VoigtVector m_conv = VoigtVector::Zero();
        bool converged = false;
        int iter = 0;
        double rs = 0.0, rq = 0.0, rf = 0.0;

        for (iter = 0; iter < max_iter; ++iter)
        {
            if (!cp_assemble(x, N, sigma_tr, depsilon, m_conv, R, &J))
            {
                opserr << "ASDPlasticMaterial3D::Closest_Point (tag " << ASDP_TAG
                       << ") - NaN in the closest-point residual at iteration "
                       << iter << " -- rejecting step" << endln;
                return LADRUNO_MATERIAL_REFUSED;
            }

            rs = 0.0; rq = 0.0;
            for (int i = 0; i < 6; ++i)
            {
                const double a = (R(i) < 0) ? -R(i) : R(i);
                if (a > rs) rs = a;
            }
            for (int i = 6; i < N - 1; ++i)
            {
                const double a = (R(i) < 0) ? -R(i) : R(i);
                if (a > rq) rq = a;
            }
            rf = (R(N - 1) < 0) ? -R(N - 1) : R(N - 1);

            if (rf <= tol_f && rs <= tol_s && rq <= tol_s)
            {
                converged = true;
                break;      // J is the Jacobian AT the converged point: reuse it
            }

            Eigen::FullPivLU<cp_matrix_t> lu(J);
            if (!lu.isInvertible())
            {
                opserr << "ASDPlasticMaterial3D::Closest_Point (tag " << ASDP_TAG
                       << ") - singular closest-point Jacobian at iteration " << iter
                       << " -- rejecting step" << endln;
                return LADRUNO_MATERIAL_REFUSED;
            }
            dx = lu.solve(R);
            for (int i = 0; i < N; ++i)
                if (!(dx(i) == dx(i)))
                {
                    opserr << "ASDPlasticMaterial3D::Closest_Point (tag " << ASDP_TAG
                           << ") - NaN in the Newton correction at iteration " << iter
                           << " -- rejecting step" << endln;
                    return LADRUNO_MATERIAL_REFUSED;
                }
            x -= dx;
        }

        if (!converged)
        {
            opserr << "ASDPlasticMaterial3D::Closest_Point (tag " << ASDP_TAG
                   << ") - the closest-point Newton exhausted " << max_iter
                   << " iterations: |f| = " << rf << " (tol " << tol_f
                   << "), |R_sigma| = " << rs << ", |R_q| = " << rq
                   << " (tol " << tol_s << ") -- rejecting step" << endln;
            return LADRUNO_MATERIAL_REFUSED;
        }
        cp_last_iterations = iter;

        const double dl = x(N - 1);
        if (dl < 0.0)
        {
            opserr << "ASDPlasticMaterial3D::Closest_Point (tag " << ASDP_TAG
                   << ") - converged with a NEGATIVE plastic multiplier (" << dl
                   << "), which is inadmissible -- rejecting step" << endln;
            return LADRUNO_MATERIAL_REFUSED;
        }

        TrialPlastic_Strain = CommitPlastic_Strain + dl * m_conv;

        if (INT_OPT_tangent_operator_type[ASDP_TAG]
                == ASDPlasticMaterial3D_Tangent_Operator_Type::Algorithmic)
        {
            // The converged system depends on epsilon_{n+1} only through s_tr, and
            // d(s_tr)/d(eps) = E(sigma_n).  Differentiating R(x(eps)) = 0:
            //     J * [dsigma ; dq ; ddl] = [E*deps ; 0 ; 0]
            //     C_alg = (J^{-1})_{sigma,sigma} * E
            // Solved with the ONE factorization already in hand -- the Schur form
            // in the ADR is for the doc, not for the code.  C_alg is UNSYMMETRIC
            // whenever m != n (non-associated Drucker-Prager, etabar != eta), so
            // models using it must run on an unsymmetric solver (UmfPack; NOT
            // ProfileSPD, and NOT PARDISO's symmetric -matrixType).
            if constexpr (el_is_stress_dependent<ElasticityType>::value)
            {
                static bool warned_dEds = false;
                if (!warned_dEds)
                {
                    opserr << "ASDPlasticMaterial3D::Closest_Point (tag " << ASDP_TAG
                           << ") - this elasticity model is stress dependent and has"
                           << " no ELASTICITY_STRESS_DERIVATIVE, so tangent_type"
                           << " Algorithmic is missing the E,sigma : (sigma - sigma_tr)"
                           << " term (ADR-97 D6/P5). The committed stress is"
                           << " unaffected." << endln;
                    warned_dEds = true;
                }
            }
            Eigen::FullPivLU<cp_matrix_t> lu(J);
            if (!lu.isInvertible())
            {
                opserr << "ASDPlasticMaterial3D::Closest_Point (tag " << ASDP_TAG
                       << ") - the converged closest-point Jacobian is singular, so"
                       << " no algorithmic tangent exists -- rejecting step" << endln;
                return LADRUNO_MATERIAL_REFUSED;
            }
            cp_matrix_t rhs(N, 6);
            rhs.setZero();
            for (int i = 0; i < 6; ++i)
                for (int j = 0; j < 6; ++j) rhs(i, j) = Eelastic(i, j);
            cp_matrix_t Z = lu.solve(rhs);
            for (int i = 0; i < 6; ++i)
                for (int j = 0; j < 6; ++j) Stiffness(i, j) = Z(i, j);
        }
        else
        {
            // Secant / Continuum / Elastic / Numerical_Algorithmic_* are evaluated
            // at the CONVERGED closest-point state, but they are still the same
            // operators ADR-94 M3 measured at 80 / 57 / 103 / 5 % against a central
            // difference of the material's own response: none of them is the
            // tangent of THIS map.  Only `Algorithmic` is.
            ComputeTangentStiffness();
        }

        if (ladruno_strict_rejects("Closest_Point", TrialStress))
            return LADRUNO_MATERIAL_REFUSED;

        return 0;
    }

'''

edit(
    "ASDPlasticMaterial3D.h",
    "int Closest_Point(const VoigtVector & strain_incr)",
    """    int Backward_Euler(const VoigtVector & strain_incr)
    {
        using namespace ASDPlasticMaterial3DGlobals;

        int errorcode = -1;""",
    CP_CODE + """    int Backward_Euler(const VoigtVector & strain_incr)
    {
        using namespace ASDPlasticMaterial3DGlobals;

        int errorcode = -1;""")


# ===========================================================================
# parser
# ===========================================================================
PARSER = os.path.join(ROOT, "SRC", "material", "nD", "ASDPlasticMaterial3D")

P = []


def pedit(marker, anchor, replacement, last=False):
    P.append(("OPS_AllASDPlasticMaterial3Ds.cpp", marker, anchor, replacement, last))


pedit(
    "Closest_Point (ADR-97",
    """       "    integration_method (string) : Backward_Euler | Forward_Euler |\\n"
       "        Forward_Euler_Subincrement | Modified_Euler_Error_Control |\\n"
       "        Runge_Kutta_45_Error_Control\\n"
""",
    """       "    integration_method (string) : Backward_Euler | Closest_Point (ADR-97) |\\n"
       "        Forward_Euler | Forward_Euler_Subincrement |\\n"
       "        Modified_Euler_Error_Control | Runge_Kutta_45_Error_Control\\n"
       "    (Closest_Point is the fully implicit closest-point return map; it is the\\n"
       "     ONLY integrator for which tangent_type Algorithmic -- the exact\\n"
       "     consistent tangent -- is defined, and it is currently available for\\n"
       "     VonMises and DruckerPrager only.)\\n"
""")

pedit(
    "Algorithmic (Closest_Point only)",
    """       "    tangent (string) : Elastic | Numerical_Algorithmic_FirstOrder | Numerical_Algorithmic_SecondOrder\\\\ \\n"
""",
    """       "    tangent_type (string) : Secant (default) | Continuum | Elastic |\\n"
       "        Algorithmic (Closest_Point only) | Numerical_Algorithmic_FirstOrder |\\n"
       "        Numerical_Algorithmic_SecondOrder\\\\ \\n"
""")

pedit(
    "ADR-97 wp/97b): the closest-point return map",
    """                    else if (std::strcmp(method_name, "Backward_Euler") == 0)
                        method = (int) ASDPlasticMaterial3D_Constitutive_Integration_Method::Backward_Euler;                    """,
    """                    else if (std::strcmp(method_name, "Backward_Euler") == 0)
                        method = (int) ASDPlasticMaterial3D_Constitutive_Integration_Method::Backward_Euler;
                    else if (std::strcmp(method_name, "Closest_Point") == 0)
                    {
                        // Ladruno (ADR-97 wp/97b): the closest-point return map.
                        // Available family by family (ADR-97 D3): P1 ships
                        // VonMises and DruckerPrager (flank + apex) with Null,
                        // Linear scalar/tensor and ArmstrongFrederick hardening.
                        // A specialization whose YF, PF or hardening law has not
                        // opted in is refused HERE rather than silently
                        // approximated with the inert zero/finite-difference
                        // defaults -- the class of defect ADR-94 M3 found.
                        if (!instance->supportsClosestPoint())
                        {
                            opserr << "nDMaterial ASDPlasticMaterial3D - "
                                   << "integration_method 'Closest_Point' is not"
                                   << " available for this model (YF "
                                   << instance->getYFName().c_str() << ", PF "
                                   << instance->getPFName().c_str() << ", IV "
                                   << instance->getIVName().c_str() << ")." << endln;
                            opserr << "   ADR-97 D3: the closest-point map is added"
                                   << " family by family. P1 covers VonMises and"
                                   << " DruckerPrager with Null / Linear / "
                                   << "ArmstrongFrederick hardening; MohrCoulomb and"
                                   << " MohrCoulombTensionCutoff are P2, HoekBrown"
                                   << " P3, StiffSoil P5. Use Backward_Euler."
                                   << endln;
                            asdp_parse_rejected = true;
                            return false;
                        }
                        method = (int) ASDPlasticMaterial3D_Constitutive_Integration_Method::Closest_Point;
                    }""")

pedit(
    "Backward_Euler, Closest_Point\" << endln;",
    """                        opserr << "   Valid values: Forward_Euler, Forward_Euler_Subincrement, "
                               << "Modified_Euler_Error_Control, Runge_Kutta_45_Error_Control, "
                               << "Backward_Euler" << endln;""",
    """                        opserr << "   Valid values: Forward_Euler, Forward_Euler_Subincrement, "
                               << "Modified_Euler_Error_Control, Runge_Kutta_45_Error_Control, "
                               << "Backward_Euler, Closest_Point" << endln;  // Ladruno (ADR-97 wp/97b)""")

pedit(
    'tangent_type_name, "Algorithmic"',
    """                    else if (std::strcmp(tangent_type_name, "Numerical_Algorithmic_FirstOrder") == 0)""",
    """                    else if (std::strcmp(tangent_type_name, "Algorithmic") == 0)
                        // Ladruno (ADR-97 wp/97b, D2): the EXACT consistent tangent
                        // of the Closest_Point map.  The enum value has existed
                        // since upstream with no dispatch case anywhere -- it was
                        // dead, which is why nobody has been silently getting it.
                        // The cross-check against integration_method is below,
                        // after the whole option loop (the two tokens may arrive in
                        // either order).
                        tangent = (int) ASDPlasticMaterial3D_Tangent_Operator_Type::Algorithmic;
                    else if (std::strcmp(tangent_type_name, "Numerical_Algorithmic_FirstOrder") == 0)""")

pedit(
    "Valid values: Elastic, Continuum, Secant, Algorithmic",
    """                        opserr << "   Valid values: Elastic, Continuum, Secant, "
                               << "Numerical_Algorithmic_FirstOrder, "
                               << "Numerical_Algorithmic_SecondOrder" << endln;""",
    """                        opserr << "   Valid values: Elastic, Continuum, Secant, Algorithmic, "  // Ladruno (ADR-97 wp/97b)
                               << "Numerical_Algorithmic_FirstOrder, "
                               << "Numerical_Algorithmic_SecondOrder" << endln;""")

pedit(
    "ADR-97 wp/97b, D2): a consistent tangent",
    """    instance->set_constitutive_integration_method(method, tangent,""",
    """    // Ladruno (ADR-97 wp/97b, D2): a consistent tangent is defined only relative
    // to a specific committed map.  Offering `Algorithmic` on the cutting-plane
    // Backward_Euler would ship a FOURTH almost-right tangent, which is exactly
    // the class of defect ADR-94 M3 found (Continuum 57 %, Secant 80 %, Elastic
    // 103 %, Numerical_Algorithmic 31 % against a central difference of the
    // material's own committed response).  Refuse, naming both tokens.
    if (tangent == (int) ASDPlasticMaterial3D_Tangent_Operator_Type::Algorithmic
            && method != (int) ASDPlasticMaterial3D_Constitutive_Integration_Method::Closest_Point)
    {
        opserr << "nDMaterial ASDPlasticMaterial3D - tangent_type 'Algorithmic' is"
               << " the exact consistent tangent of the 'Closest_Point' return map"
               << " and is REFUSED with any other integration_method (ADR-97 D2)."
               << endln;
        opserr << "   Set 'integration_method Closest_Point', or pick a tangent_type"
               << " that the chosen integrator defines: Secant, Continuum, Elastic,"
               << " Numerical_Algorithmic_FirstOrder, Numerical_Algorithmic_SecondOrder."
               << endln;
        asdp_parse_rejected = true;
        return false;
    }

    instance->set_constitutive_integration_method(method, tangent,""")


# ---------------------------------------------------------------------------
# 7. apex region test: the degenerate hydrostatic trial state
# ---------------------------------------------------------------------------
edit(
    "ASDPlasticMaterial3D.h",
    "a trial state ON the hydrostatic axis",
    """    bool cp_apex_region(const VoigtVector& depsilon, const VoigtVector& sigma_tr,
                        const VoigtMatrix& Eelastic, double f_tr)""",
    """    bool cp_apex_region(const VoigtVector& depsilon, const VoigtVector& sigma_tr,
                        const VoigtMatrix& Eelastic, double f_tr, double tol_f)""")

edit(
    "ASDPlasticMaterial3D.h",
    "degenerates to 0 < 0",
    """        const double den = n0.dot(Em0) - H0;
        if (!(den > MACHINE_EPSILON)) return false;
        const double dl0 = f_tr / den;
        VoigtVector dev_tr  = sigma_tr.deviator();
        VoigtVector dev_Em0 = Em0.deviator();
        VoigtVector dev_ret = dev_tr - dl0 * dev_Em0;
        return tensor_dot_stress_like(dev_ret, dev_tr) < 0.0;""",
    """        VoigtVector dev_tr = sigma_tr.deviator();
        double q2_tr = tensor_dot_stress_like(dev_tr, dev_tr);
        if (!(q2_tr > 0.0)) q2_tr = 0.0;
        const double q_tr = std::sqrt(q2_tr);
        // A trial state ON the hydrostatic axis and outside the surface can only
        // return to the vertex.  The flip test below is a SIGN test on
        // dot(dev_ret, dev_tr), which degenerates to 0 < 0 -- i.e. says CONE --
        // when the trial deviator vanishes; the cone Newton then has no flow
        // direction at all and exhausts its iterations.  That degenerate state
        // is not exotic: it is exactly ADR-94 B4's hydrostatic-tension
        // reproducer, the one that used to commit NaN.  `tol_f` is the yield
        // tolerance, so this comparison is in stress units and unit consistent
        // (ADR-94 M5).
        if (q_tr <= tol_f) return true;
        const double den = n0.dot(Em0) - H0;
        if (!(den > MACHINE_EPSILON)) return false;
        const double dl0 = f_tr / den;
        VoigtVector dev_Em0 = Em0.deviator();
        VoigtVector dev_ret = dev_tr - dl0 * dev_Em0;
        return tensor_dot_stress_like(dev_ret, dev_tr) < 0.0;""")

edit(
    "ASDPlasticMaterial3D.h",
    "cp_apex_region(depsilon, sigma_tr, Eelastic, f_tr, tol_f)",
    """            if (cp_apex_region(depsilon, sigma_tr, Eelastic, f_tr))""",
    """            if (cp_apex_region(depsilon, sigma_tr, Eelastic, f_tr, tol_f))""")


def main():
    # the CRTP fix is a plain textual substitution, applied to every occurrence
    fixed = 0
    for pf in ("VonMises_PF", "DruckerPrager_PF"):
        path = os.path.join(ASD, "PlasticFlowDirections", "%s.h" % pf)
        with open(path, "r", encoding="utf-8", errors="surrogateescape",
                  newline="") as f:
            src = f.read()
        if "this->dm_dsigma_buffer" not in src and "dm_dsigma_buffer" in src:
            src = src.replace("dm_dsigma_buffer", "this->dm_dsigma_buffer")
            with open(path, "w", encoding="utf-8", errors="surrogateescape",
                      newline="") as f:
                f.write(src)
            fixed += 1
            print("  applied  PlasticFlowDirections/%-32s CRTP this-> qualification" % (pf + ".h"))
    a1, s1 = apply_edits(E, ASD)
    a2, s2 = apply_edits(P, PARSER)
    print("apply_p1_core: %d applied (+%d CRTP), %d already present"
          % (a1 + a2, fixed, s1 + s2))


if __name__ == "__main__":
    main()
