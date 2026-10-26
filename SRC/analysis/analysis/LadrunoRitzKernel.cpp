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

// Ladruno WP-213: shared block Krylov-Ritz kernel on (K, M) (ADR 213 D1).
// Algorithm and conventions: LadrunoRitzKernel.h.

#include <LadrunoRitzKernel.h>
#include <AnalysisModel.h>
#include <LinearSOE.h>
#include <Domain.h>
#include <Node.h>
#include <FE_Element.h>
#include <FE_EleIter.h>
#include <DOF_Group.h>
#include <DOF_GrpIter.h>
#include <TransformationDOF_Group.h>
#include <MP_Constraint.h>
#include <MP_ConstraintIter.h>
#include <SP_Constraint.h>
#include <SP_ConstraintIter.h>
#include <Matrix.h>
#include <ID.h>
#include <LoadPattern.h>
#include <NodalLoad.h>
#include <NodalLoadIter.h>
#include <ElementalLoadIter.h>
#include <OPS_Globals.h>
#include <cmath>
#include <map>
#include <set>
#include <utility>
#include <sstream>
#include <algorithm>
#include <limits>

// dsyevx, not dsyev: the LAPACK subset bundled in OTHER/LAPACK (what the
// Linux build links) has dsyevx but not dsyev (WP-213 Zone-A link failure)
#ifdef _WIN32
extern "C" int DSYEVX(char *JOBZ, char *RANGE, char *UPLO, int *N, double *A,
                      int *LDA, double *VL, double *VU, int *IL, int *IU,
                      double *ABSTOL, int *M, double *W, double *Z, int *LDZ,
                      double *WORK, int *LWORK, int *IWORK, int *IFAIL,
                      int *INFO);
#define LADRUNO_RITZ_DSYEVX DSYEVX
#else
extern "C" int dsyevx_(char *JOBZ, char *RANGE, char *UPLO, int *N, double *A,
                       int *LDA, double *VL, double *VU, int *IL, int *IU,
                       double *ABSTOL, int *M, double *W, double *Z, int *LDZ,
                       double *WORK, int *LWORK, int *IWORK, int *IFAIL,
                       int *INFO);
#define LADRUNO_RITZ_DSYEVX dsyevx_
#endif

// ===========================================================================
// mode-kind registry
// ===========================================================================
namespace
{
    struct ModeKindRecord
    {
        Vector              eigenvalues;   // snapshot of the domain's set
        std::vector<int>    kinds;
        std::vector<double> saPeriods;     // 0 = no override
        std::vector<int>    dirMasks;      // residualModes directions, bit d-1
    };

    std::map<const Domain*, ModeKindRecord>& modeKindRegistry()
    {
        static std::map<const Domain*, ModeKindRecord> theRegistry;
        return theRegistry;
    }

    // the domain's eigenvalues, or an empty vector when none were ever set
    // (Domain::getEigenvalues exits in that case)
    Vector currentEigenvalues(const Domain* theDomain)
    {
        Domain* d = const_cast<Domain*>(theDomain);
        if (d == 0 || d->getNumEigenvalues() <= 0)
            return Vector();
        return d->getEigenvalues();
    }

    bool bitIdentical(const Vector& a, const Vector& b)
    {
        if (a.Size() != b.Size())
            return false;
        for (int i = 0; i < a.Size(); i++)
            if (a(i) != b(i))
                return false;
        return true;
    }
}

namespace
{
    // the record for the domain's CURRENT mode set, or 0 when absent / stale
    const ModeKindRecord* validRecord(const Domain* theDomain, int& numModes)
    {
        Vector ev = currentEigenvalues(theDomain);
        numModes = ev.Size();
        std::map<const Domain*, ModeKindRecord>& reg = modeKindRegistry();
        std::map<const Domain*, ModeKindRecord>::const_iterator it = reg.find(theDomain);
        if (it != reg.end() && bitIdentical(it->second.eigenvalues, ev) &&
            (int)it->second.kinds.size() == numModes)
            return &it->second;
        return 0;
    }
}

std::vector<int>
LadrunoModeKinds::get(const Domain* theDomain)
{
    int n = 0;
    const ModeKindRecord* rec = validRecord(theDomain, n);
    return rec ? rec->kinds : std::vector<int>(n, EIGEN);
}

std::vector<double>
LadrunoModeKinds::saPeriods(const Domain* theDomain)
{
    int n = 0;
    const ModeKindRecord* rec = validRecord(theDomain, n);
    if (rec && (int)rec->saPeriods.size() == n)
        return rec->saPeriods;
    return std::vector<double>(n, 0.0);
}

std::vector<int>
LadrunoModeKinds::dirMasks(const Domain* theDomain)
{
    int n = 0;
    const ModeKindRecord* rec = validRecord(theDomain, n);
    if (rec && (int)rec->dirMasks.size() == n)
        return rec->dirMasks;
    return std::vector<int>(n, 0);
}

double
LadrunoModeKinds::saPeriod(const Domain* theDomain, int mode0)
{
    int n = 0;
    const ModeKindRecord* rec = validRecord(theDomain, n);
    if (rec == 0 || mode0 < 0 || mode0 >= (int)rec->saPeriods.size())
        return 0.0;
    return rec->saPeriods[mode0];
}

void
LadrunoModeKinds::set(const Domain* theDomain, const std::vector<int>& kinds,
                      const std::vector<double>& saPeriods,
                      const std::vector<int>& dirMasks)
{
    ModeKindRecord rec;
    rec.eigenvalues = currentEigenvalues(theDomain);
    rec.kinds = kinds;
    rec.saPeriods = saPeriods;
    rec.saPeriods.resize(kinds.size(), 0.0);
    rec.dirMasks = dirMasks;
    rec.dirMasks.resize(kinds.size(), 0);
    modeKindRegistry()[theDomain] = rec;
}

void
LadrunoModeKinds::clear(const Domain* theDomain)
{
    modeKindRegistry().erase(theDomain);
}

char
LadrunoModeKinds::label(int kind)
{
    switch (kind) {
    case RESIDUAL: return 'R';
    case RITZ:     return 'Z';
    default:       return 'E';
    }
}

bool
LadrunoModeKinds::has(const Domain* theDomain, int kind)
{
    int n = 0;
    const ModeKindRecord* rec = validRecord(theDomain, n);
    return rec != 0 && std::find(rec->kinds.begin(), rec->kinds.end(), kind) != rec->kinds.end();
}

bool
LadrunoModeKinds::warnRitzBasis(const Domain* theDomain, const char* cmd)
{
    if (!has(theDomain, RITZ))
        return false;
    opserr << cmd << ": WARNING the mode set holds Ritz vectors (ritzVectors): they are"
           << " tuned to a load, not to a frequency band, and the higher ones are not"
           << " natural modes - use eigen for frequency-domain work\n";
    return true;
}

const char*
LadrunoModeKinds::Report::tag(int mode0) const
{
    if (mode0 < 0 || mode0 >= (int)tags.size())
        return "";
    return tags[mode0].c_str();
}

namespace
{
    // 1-based ascending mode numbers as "4-5" / "3, 6-7"
    std::string modeRanges(const std::vector<int>& modes)
    {
        std::ostringstream s;
        for (size_t i = 0; i < modes.size(); ) {
            size_t j = i;
            while (j + 1 < modes.size() && modes[j + 1] == modes[j] + 1)
                j++;
            if (i > 0)
                s << ", ";
            s << modes[i];
            if (j > i)
                s << "-" << modes[j];
            i = j + 1;
        }
        return s.str();
    }
}

LadrunoModeKinds::Report
LadrunoModeKinds::report(const Domain* theDomain, const Vector& lambdas)
{
    Report r;
    int n = 0;
    const ModeKindRecord* rec = validRecord(theDomain, n);
    if (rec == 0 || !bitIdentical(rec->eigenvalues, lambdas))
        return r;
    std::vector<int> residual, ritz;
    for (int k = 0; k < n; k++) {
        if (rec->kinds[k] == RESIDUAL)
            residual.push_back(k + 1);
        else if (rec->kinds[k] == RITZ)
            ritz.push_back(k + 1);
    }
    if (residual.empty() && ritz.empty())
        return r;

    r.kinds = rec->kinds;
    r.tags.assign(n, std::string());
    for (int k = 0; k < n; k++)
        if (rec->kinds[k] != EIGEN)
            r.tags[k] = std::string("  # ") + label(rec->kinds[k]);

    std::ostringstream note;
    if (!residual.empty()) {
        const bool one = residual.size() == 1;
        note << "# NOTE: " << (one ? "mode " : "modes ") << modeRanges(residual)
             << " (R) " << (one ? "is a residual pseudo-mode" : "are residual pseudo-modes")
             << " (static correction); " << (one ? "its period is" : "their periods are")
             << " Rayleigh, not natural.";
        // -period $T: one override shared by the whole residual set (D3b)
        const double T = rec->saPeriods[residual[0] - 1];
        bool same = T > 0.0;
        for (size_t i = 1; i < residual.size(); i++)
            if (rec->saPeriods[residual[i] - 1] != T)
                same = false;
        if (same)
            note << " responseSpectrumAnalysis reads " << (one ? "its" : "their")
                 << " Sa at T = " << T << ".";
        note << "\n";
    }
    if (!ritz.empty()) {
        const bool one = ritz.size() == 1;
        note << "# NOTE: " << (one ? "mode " : "modes ") << modeRanges(ritz)
             << " (Z) " << (one ? "is a Ritz vector" : "are Ritz vectors")
             << "; only the lowest approximate natural periods (upper bounds),"
             << " the higher ones are not modes: do not report their periods as"
             << " natural periods.\n";
    }
    r.note = note.str();
    return r;
}

// ===========================================================================
// kernel
// ===========================================================================
LadrunoRitzKernel::LadrunoRitzKernel(AnalysisModel& theModel, LinearSOE& theSOE)
    : m_model(&theModel), m_soe(&theSOE),
      m_domain(theModel.getDomainPtr()), m_neq(theModel.getNumEqn())
{
}

double
LadrunoRitzKernel::dot(const Vector& a, const Vector& b)
{
    double s = 0.0;
    const int n = a.Size();
    for (int i = 0; i < n; i++)
        s += a(i) * b(i);
    return s;
}

void
LadrunoRitzKernel::applyM(const Vector& x, Vector& y)
{
    if (y.Size() != m_neq)
        y.resize(m_neq);
    y.Zero();

    FE_EleIter& theFEs = m_model->getFEs();
    FE_Element* fe;
    while ((fe = theFEs()) != 0) {
        // an FE_Element without an Element (Penalty / Lagrange constraint FEs,
        // contact adapters) carries no mass (same rule as eigen, WP-183)
        if (fe->getElement() == 0)
            continue;
        const Vector& f = fe->getM_Force(x, 1.0);
        const ID& id = fe->getID();
        for (int i = 0; i < id.Size(); i++) {
            const int loc = id(i);
            if (loc >= 0 && loc < m_neq)
                y(loc) += f(i);
        }
    }

    DOF_GrpIter& theDOFs = m_model->getDOFs();
    DOF_Group* dof;
    while ((dof = theDOFs()) != 0) {
        // Lagrange multiplier groups have no node and no mass
        if (m_domain->getNode(dof->getNodeTag()) == 0)
            continue;
        const Vector& f = dof->getM_Force(x, 1.0);
        const ID& id = dof->getID();
        for (int i = 0; i < id.Size(); i++) {
            const int loc = id(i);
            if (loc >= 0 && loc < m_neq)
                y(loc) += f(i);
        }
    }
}

void
LadrunoRitzKernel::applyK(const std::vector<Vector>& X, std::vector<Vector>& KX,
                          std::vector<double>* absQuad,
                          std::vector<Vector>* absKX)
{
    const int n = (int)X.size();
    KX.assign(n, Vector(m_neq));
    if (absQuad != 0)
        absQuad->assign(n, 0.0);
    if (absKX != 0)
        absKX->assign(n, Vector(m_neq));
    if (n == 0)
        return;

    Vector xl, yl;
    FE_EleIter& theFEs = m_model->getFEs();
    FE_Element* fe;
    while ((fe = theFEs()) != 0) {
        // the same element contribution eigen assembles into K
        fe->zeroTangent();
        fe->addKtToTang(1.0);
        const Matrix& kt = fe->getTangent(0);
        const ID& id = fe->getID();
        const int nd = id.Size();
        if (nd == 0)
            continue;
        xl.resize(nd);
        yl.resize(nd);
        for (int j = 0; j < n; j++) {
            for (int i = 0; i < nd; i++) {
                const int loc = id(i);
                xl(i) = (loc >= 0 && loc < m_neq) ? X[j](loc) : 0.0;
            }
            yl.addMatrixVector(0.0, kt, xl, 1.0);
            if (absQuad != 0) {
                double q = 0.0;
                for (int i = 0; i < nd; i++)
                    for (int l = 0; l < nd; l++)
                        q += std::fabs(kt(i, l)) * std::fabs(xl(i)) * std::fabs(xl(l));
                (*absQuad)[j] += q;
            }
            if (absKX != 0) {
                for (int i = 0; i < nd; i++) {
                    const int loc = id(i);
                    if (loc < 0 || loc >= m_neq)
                        continue;
                    double s = 0.0;
                    for (int l = 0; l < nd; l++)
                        s += std::fabs(kt(i, l)) * std::fabs(xl(l));
                    (*absKX)[j](loc) += s;
                }
            }
            for (int i = 0; i < nd; i++) {
                const int loc = id(i);
                if (loc >= 0 && loc < m_neq)
                    KX[j](loc) += yl(i);
            }
        }
    }
}

template <class Getter>
void
LadrunoRitzKernel::nodeToEq(Getter getter, Vector& out)
{
    if (out.Size() != m_neq)
        out.resize(m_neq);
    out.Zero();

    // constrained node -> its MP: a TransformationDOF_Group with an MP orders
    // its ID as [own unconstrained DOFs..., retained node's DOFs...]; only the
    // own part belongs to this node (the retained part is written by the
    // retained node's group).
    std::map<int, MP_Constraint*> mpOf;
    MP_ConstraintIter& theMPs = m_domain->getMPs();
    MP_Constraint* mp;
    while ((mp = theMPs()) != 0)
        mpOf[mp->getNodeConstrained()] = mp;

    Vector v;
    DOF_GrpIter& theDOFs = m_model->getDOFs();
    DOF_Group* dof;
    while ((dof = theDOFs()) != 0) {
        Node* nd = m_domain->getNode(dof->getNodeTag());
        if (nd == 0)
            continue;   // Lagrange multiplier group
        const int ndf = nd->getNumberDOF();
        v.resize(ndf);
        v.Zero();
        getter(nd, v);
        const ID& id = dof->getID();

        MP_Constraint* theMP = 0;
        if (dynamic_cast<TransformationDOF_Group*>(dof) != 0) {
            std::map<int, MP_Constraint*>::iterator it = mpOf.find(nd->getTag());
            if (it != mpOf.end())
                theMP = it->second;
        }

        if (theMP == 0) {
            for (int i = 0; i < ndf && i < id.Size(); i++) {
                const int loc = id(i);
                if (loc >= 0 && loc < m_neq)
                    out(loc) = v(i);
            }
        } else {
            const ID& cdof = theMP->getConstrainedDOFs();
            int k = 0;
            for (int i = 0; i < ndf; i++) {
                if (cdof.getLocation(i) >= 0)
                    continue;
                if (k < id.Size()) {
                    const int loc = id(k);
                    if (loc >= 0 && loc < m_neq)
                        out(loc) = v(i);
                }
                k++;
            }
        }
    }
}

int
LadrunoRitzKernel::rigidInfluence(int dir, Vector& iota)
{
    if (dir < 1) {
        m_err = "direction must be >= 1";
        return -1;
    }
    std::set<std::pair<int, int> > fixed;
    SP_ConstraintIter& theSPs = m_domain->getSPs();
    SP_Constraint* sp;
    while ((sp = theSPs()) != 0)
        fixed.insert(std::make_pair(sp->getNodeTag(), sp->getDOF_Number()));
    nodeToEq([dir, &fixed](Node* nd, Vector& v) {
        const int ndm = nd->getCrds().Size();
        if (dir <= ndm && dir - 1 < v.Size() &&
            fixed.count(std::make_pair(nd->getTag(), dir - 1)) == 0)
            v(dir - 1) = 1.0;
    }, iota);
    return 0;
}

int
LadrunoRitzKernel::patternLoad(int patternTag, Vector& f, int& numElementalLoads)
{
    numElementalLoads = 0;
    if (f.Size() != m_neq)
        f.resize(m_neq);
    f.Zero();
    LoadPattern* pat = m_domain->getLoadPattern(patternTag);
    if (pat == 0) {
        std::ostringstream s;
        s << "no load pattern with tag " << patternTag;
        m_err = s.str();
        return -1;
    }
    ElementalLoadIter& theEleLoads = pat->getElementalLoads();
    while (theEleLoads() != 0)
        numElementalLoads++;

    // reference load per node (a node may carry several NodalLoads)
    std::map<int, Vector> pNode;
    NodalLoadIter& theLoads = pat->getNodalLoads();
    NodalLoad* nl;
    while ((nl = theLoads()) != 0) {
        Node* nd = m_domain->getNode(nl->getNodeTag());
        if (nd == 0)
            continue;
        int type = 0;
        const Vector& data = nl->getData(type);
        std::map<int, Vector>::iterator it = pNode.find(nd->getTag());
        if (it == pNode.end())
            it = pNode.insert(std::make_pair(nd->getTag(), Vector(nd->getNumberDOF()))).first;
        for (int i = 0; i < data.Size() && i < it->second.Size(); i++)
            it->second(i) += data(i);
    }

    Vector tp;
    DOF_GrpIter& theDOFs = m_model->getDOFs();
    DOF_Group* dof;
    while ((dof = theDOFs()) != 0) {
        std::map<int, Vector>::const_iterator it = pNode.find(dof->getNodeTag());
        if (it == pNode.end())
            continue;   // unloaded node, or a Lagrange multiplier group
        const Vector& p = it->second;
        const ID& id = dof->getID();
        TransformationDOF_Group* tdof = dynamic_cast<TransformationDOF_Group*>(dof);
        Matrix* T = tdof ? tdof->getT() : 0;
        if (T != 0) {
            // generalized load of the group: T^T p (own DOFs, then the
            // retained node's DOFs, the ID order)
            tp.resize(T->noCols());
            tp.Zero();
            for (int k = 0; k < T->noCols(); k++)
                for (int i = 0; i < T->noRows() && i < p.Size(); i++)
                    tp(k) += (*T)(i, k) * p(i);
            for (int k = 0; k < tp.Size() && k < id.Size(); k++) {
                const int loc = id(k);
                if (loc >= 0 && loc < m_neq)
                    f(loc) += tp(k);
            }
        } else {
            for (int i = 0; i < p.Size() && i < id.Size(); i++) {
                const int loc = id(i);
                if (loc >= 0 && loc < m_neq)
                    f(loc) += p(i);
            }
        }
    }
    return 0;
}

void
LadrunoRitzKernel::massDiagonal(Vector& d, bool& offDiagonal)
{
    if (d.Size() != m_neq)
        d.resize(m_neq);
    d.Zero();
    offDiagonal = false;
    auto addMass = [&](const Matrix& m, const ID& id) {
        const int nd = id.Size();
        for (int i = 0; i < nd && i < m.noRows(); i++) {
            const int loc = id(i);
            for (int j = 0; j < nd && j < m.noCols(); j++) {
                if (i == j || m(i, j) == 0.0)
                    continue;
                const int locj = id(j);
                if (loc >= 0 && loc < m_neq && locj >= 0 && locj < m_neq && locj != loc)
                    offDiagonal = true;
            }
            if (loc >= 0 && loc < m_neq)
                d(loc) += m(i, i);
        }
    };
    FE_EleIter& theFEs = m_model->getFEs();
    FE_Element* fe;
    while ((fe = theFEs()) != 0) {
        if (fe->getElement() == 0)
            continue;   // constraint / contact FE: no mass (applyM's rule)
        fe->zeroTangent();
        fe->addMtoTang(1.0);
        addMass(fe->getTangent(0), fe->getID());
    }
    DOF_GrpIter& theDOFs = m_model->getDOFs();
    DOF_Group* dof;
    while ((dof = theDOFs()) != 0) {
        if (m_domain->getNode(dof->getNodeTag()) == 0)
            continue;   // Lagrange multiplier group
        dof->zeroTangent();
        dof->addMtoTang(1.0);
        addMass(dof->getTangent(0), dof->getID());
    }
}

int
LadrunoRitzKernel::domainModes(std::vector<Vector>& phi)
{
    phi.clear();
    const int p = m_domain->getNumEigenvalues();
    if (p <= 0)
        return 0;
    phi.assign(p, Vector(m_neq));
    for (int k = 0; k < p; k++) {
        bool bad = false;
        nodeToEq([k, &bad](Node* nd, Vector& v) {
            const Matrix& E = nd->getEigenvectors();
            if (E.noCols() <= k) {
                bad = true;
                return;
            }
            for (int i = 0; i < v.Size() && i < E.noRows(); i++)
                v(i) = E(i, k);
        }, phi[k]);
        if (bad) {
            m_err = "a node carries fewer eigenvectors than the domain has "
                    "eigenvalues (node added after eigen?) - re-run eigen";
            return -1;
        }
    }
    return 0;
}

int
LadrunoRitzKernel::install(const std::vector<Vector>& modes,
                           const std::vector<double>& lambdas)
{
    const int n = (int)modes.size();
    if (n == 0 || (int)lambdas.size() != n) {
        m_err = "install: empty or inconsistent mode set";
        return -1;
    }
    m_model->setNumEigenvectors(n);
    Vector ev(n);
    for (int i = 0; i < n; i++) {
        m_model->setEigenvector(i + 1, modes[i]);
        ev(i) = lambdas[i];
    }
    m_model->setEigenvalues(ev);
    return 0;
}

int
LadrunoRitzKernel::formTangent()
{
    if (m_soe->getNumEqn() != m_neq) {
        std::ostringstream s;
        s << "the analysis LinearSOE has " << m_soe->getNumEqn()
          << " equations but the AnalysisModel has " << m_neq
          << " - run eigen (or an analysis step) on the current model first";
        m_err = s.str();
        return -1;
    }
    m_soe->zeroA();
    FE_EleIter& theFEs = m_model->getFEs();
    FE_Element* fe;
    while ((fe = theFEs()) != 0) {
        fe->zeroTangent();
        fe->addKtToTang(1.0);
        if (m_soe->addA(fe->getTangent(0), fe->getID()) < 0) {
            m_err = "LinearSOE::addA failed while forming K";
            return -1;
        }
    }
    return 0;
}

int
LadrunoRitzKernel::solve(const Vector& b, Vector& x)
{
    m_soe->setB(b);
    if (m_soe->solve() < 0) {
        m_err = "LinearSOE::solve failed on K (singular tangent, or a system "
                "that cannot factor this matrix - Lagrange needs an "
                "indefinite solver such as UmfPack, Mumps or Pardiso)";
        return -1;
    }
    x = m_soe->getX();
    return 0;
}

int
LadrunoRitzKernel::run(const std::vector<Vector>& F,
                       const std::vector<Vector>* Phi,
                       const Options& opt,
                       Result& out)
{
    out = Result();
    const int q = (int)F.size();
    if (q == 0) {
        m_err = "no starting vectors";
        return -1;
    }
    if (m_neq <= 0) {
        m_err = "the AnalysisModel has no equations";
        return -1;
    }
    for (int j = 0; j < q; j++)
        if (F[j].Size() != m_neq) {
            m_err = "a starting vector is not in equation numbering";
            return -1;
        }
    const int maxBlocks = std::max(1, opt.maxBlocks);
    out.keptPerStart.assign(q, 0);
    out.droppedPerStart.assign(q, 0);
    out.firstBlockRelNorm.assign(q, 0.0);

    // ---- M-orthonormal copy of the existing basis ------------------------
    std::vector<Vector> Q, MQ;
    Vector My(m_neq);
    if (Phi != 0) {
        for (size_t k = 0; k < Phi->size(); k++) {
            Vector y = (*Phi)[k];
            if (y.Size() != m_neq) {
                m_err = "an existing mode is not in equation numbering";
                return -1;
            }
            applyM(y, My);
            const double n0 = std::sqrt(std::max(0.0, dot(y, My)));
            if (!(n0 > 0.0))
                continue;   // a massless mode spans nothing in the M-norm
            for (int pass = 0; pass < 2; pass++)
                for (size_t i = 0; i < Q.size(); i++)
                    y.addVector(1.0, Q[i], -dot(MQ[i], y));
            applyM(y, My);
            const double n = std::sqrt(std::max(0.0, dot(y, My)));
            if (!(n > 1.0e-8 * n0))
                continue;   // dependent on earlier modes
            y /= n;
            My /= n;
            Q.push_back(y);
            MQ.push_back(My);
        }
    }

    // ---- factor K ----------------------------------------------------------
    if (formTangent() < 0)
        return -1;

    // x = K^-1 b for every b; the first batch is checked: the solves must
    // actually be K^-1 b (catches an SOE that cannot hold K - a diagonal /
    // explicit system - or one sized for another model).
    // WP-215 (#1025): the measure is the normwise backward error
    //   ||K x - b|| / (|| |K||x| || + ||b||),
    // |K||x| from the same element sweep as K x. A backward-stable solve sits
    // at ~eps whatever the scale of K; ||K x - b|| / ||b|| alone reached
    // 1.6e-4 under Penalty (alpha = 1e14, ||K|| ||x|| / ||b|| ~ 1e12) and
    // refused good solves. It is never above the old measure, so every solve
    // accepted before is still accepted.
    bool checked = false;
    auto solveAll = [&](const std::vector<Vector>& B, std::vector<Vector>& Xs) -> int {
        Xs.assign(B.size(), Vector(m_neq));
        for (size_t j = 0; j < B.size(); j++)
            if (solve(B[j], Xs[j]) < 0)
                return -1;
        if (checked || B.empty())
            return 0;
        checked = true;
        std::vector<Vector> KX0, absKX0;
        applyK(Xs, KX0, 0, &absKX0);
        for (size_t j = 0; j < B.size(); j++) {
            const double nb = B[j].Norm();
            if (!(nb > 0.0))
                continue;
            Vector res = KX0[j];
            res.addVector(1.0, B[j], -1.0);
            out.maxSolveResidual = std::max(out.maxSolveResidual,
                                            res.Norm() / (absKX0[j].Norm() + nb));
        }
        if (!(out.maxSolveResidual < 1.0e-4)) {
            std::ostringstream s;
            s << "the analysis LinearSOE did not solve K x = f (backward error "
              << out.maxSolveResidual << "); use a general sparse system "
              << "(UmfPack, Pardiso, Mumps, BandGeneral, ProfileSPD for "
              << "Transformation / Plain)";
            m_err = s.str();
            return -1;
        }
        return 0;
    };

    // ---- first block: K^-1 F (load starts) or F itself (shape starts) -----
    std::vector<Vector> Y;      // raw columns of the current block
    std::vector<int>    origin; // starting vector each column descends from
    for (int j = 0; j < q; j++) {
        if (!(F[j].Norm() > 0.0)) {
            m_err = "a starting vector is zero";
            return -1;
        }
        origin.push_back(j);
    }
    if (opt.solveStart) {
        if (solveAll(F, Y) < 0)
            return -1;
    } else {
        Y = F;
    }

    // ---- block loop -------------------------------------------------------
    std::vector<Vector> X, MX;
    bool full = false;
    for (int b = 0; b < maxBlocks && !Y.empty() && !full; b++) {
        std::vector<int> acceptedIdx, acceptedOrigin;
        for (size_t j = 0; j < Y.size(); j++) {
            Vector& y = Y[j];
            const int start = origin[j];
            applyM(y, My);
            const double n0 = std::sqrt(std::max(0.0, dot(y, My)));
            if (!(n0 > 0.0)) {
                out.droppedPerStart[start]++;
                out.numDropped++;
                continue;
            }
            // twice is enough (classical re-orthogonalization)
            for (int pass = 0; pass < 2; pass++) {
                for (size_t i = 0; i < Q.size(); i++)
                    y.addVector(1.0, Q[i], -dot(MQ[i], y));
                for (size_t i = 0; i < X.size(); i++)
                    y.addVector(1.0, X[i], -dot(MX[i], y));
            }
            applyM(y, My);
            const double n = std::sqrt(std::max(0.0, dot(y, My)));
            if (b == 0)
                out.firstBlockRelNorm[start] = (n > opt.dropTol * n0) ? n / n0 : 0.0;
            if (!(n > opt.dropTol * n0)) {
                // breakdown: already spanned by Phi and the earlier columns
                out.droppedPerStart[start]++;
                out.numDropped++;
                continue;
            }
            y /= n;
            My /= n;
            X.push_back(y);
            MX.push_back(My);
            out.keptPerStart[start]++;
            acceptedIdx.push_back((int)X.size() - 1);
            acceptedOrigin.push_back(start);
            if (opt.maxVectors > 0 && (int)X.size() >= opt.maxVectors) {
                full = true;
                break;
            }
        }
        if (full || b + 1 >= maxBlocks)
            break;

        // next block: K^-1 M x for every surviving column whose start still
        // has cycles left
        std::vector<Vector> rhs;
        std::vector<int> originNext;
        for (size_t a = 0; a < acceptedIdx.size(); a++) {
            const int start = acceptedOrigin[a];
            const int cap = (start < (int)opt.cyclesPerStart.size() &&
                             opt.cyclesPerStart[start] > 0)
                ? opt.cyclesPerStart[start] : maxBlocks;
            if (b + 1 >= cap)
                continue;
            rhs.push_back(MX[acceptedIdx[a]]);
            originNext.push_back(start);
        }
        if (solveAll(rhs, Y) < 0)
            return -1;
        origin.swap(originNext);
    }

    const int r = (int)X.size();
    if (r == 0)
        return 0;   // everything was already spanned

    // ---- Rayleigh-Ritz ----------------------------------------------------
    std::vector<Vector> KX;
    applyK(X, KX);
    std::vector<double> A((size_t)r * r), w(r);
    for (int i = 0; i < r; i++)
        for (int j = i; j < r; j++) {
            const double kij = 0.5 * (dot(X[i], KX[j]) + dot(X[j], KX[i]));
            A[(size_t)i + (size_t)j * r] = kij;
            A[(size_t)j + (size_t)i * r] = kij;
        }

    // all eigenpairs (RANGE = 'A'); Z receives the eigenvectors, column-major
    char jobz = 'V', range = 'A', uplo = 'U';
    int n = r, lda = r, ldz = r, il = 1, iu = r, m = 0, info = 0, lwork = -1;
    // abstol = 2 * underflow: the most accurate eigenvalues dsyevx gives
    double vl = 0.0, vu = 0.0, abstol = 2.0 * std::numeric_limits<double>::min(), wq = 0.0;
    std::vector<double> Z((size_t)r * r);
    std::vector<int> iwork(5 * (size_t)r), ifail(r);
    LADRUNO_RITZ_DSYEVX(&jobz, &range, &uplo, &n, A.data(), &lda, &vl, &vu, &il, &iu,
                        &abstol, &m, w.data(), Z.data(), &ldz, &wq, &lwork,
                        iwork.data(), ifail.data(), &info);
    if (info != 0) {
        m_err = "LAPACK dsyevx workspace query failed";
        return -1;
    }
    lwork = std::max(8 * r, (int)wq);
    std::vector<double> work(lwork);
    LADRUNO_RITZ_DSYEVX(&jobz, &range, &uplo, &n, A.data(), &lda, &vl, &vu, &il, &iu,
                        &abstol, &m, w.data(), Z.data(), &ldz, work.data(), &lwork,
                        iwork.data(), ifail.data(), &info);
    if (info != 0 || m != r) {
        m_err = "LAPACK dsyevx failed on the reduced stiffness";
        return -1;
    }
    A.swap(Z);   // eigenvectors, read below as A[i + k r]
    for (int k = 0; k < r; k++)
        if (!(w[k] > 0.0) || !std::isfinite(w[k])) {
            std::ostringstream s;
            s << "non-positive Ritz eigenvalue " << w[k]
              << " (the tangent is not positive definite on the Ritz subspace:"
              << " mechanism, or a P-Delta / softened state past the limit)";
            m_err = s.str();
            return -1;
        }

    out.psi.assign(r, Vector(m_neq));
    out.lambda.assign(w.begin(), w.end());
    for (int k = 0; k < r; k++)
        for (int i = 0; i < r; i++)
            out.psi[k].addVector(1.0, X[i], A[(size_t)i + (size_t)k * r]);
    return 0;
}
