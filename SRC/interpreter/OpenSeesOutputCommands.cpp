/* *****************************************************************************
Copyright (c) 2015-2017, The Regents of the University of California (Regents).
All rights reserved.

Redistribution and use in source and binary forms, with or without 
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies,
either expressed or implied, of the FreeBSD Project.

REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF ANY, PROVIDED HEREUNDER IS 
PROVIDED "AS IS". REGENTS HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, 
UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

*************************************************************************** */

// Written: Minjie

// Description: commands to output

#include <elementAPI.h>
#include <ConstraintHandler.h>                 // Ladruno: ADR-30 P3 tie-force query
#include <LadrunoProjectionHandler.h>          // Ladruno: ADR-30 P3
extern ConstraintHandler **OPS_GetHandler(void);   // Ladruno: defined in OpenSeesCommands.cpp
#include <Domain.h>
#include <Node.h>
#include <NodeIter.h>
#include <DOF_Group.h>
#include <Matrix.h>
#include <LoadPattern.h>
#include <LoadPatternIter.h>
#include <FileStream.h>
#include <ID.h>
#include <NodalLoad.h>
#include <NodalLoadIter.h>
#include <ElementalLoad.h>
#include <ElementalLoadIter.h>
#include <Element.h>
#include <ElementIter.h>
#include <CrdTransf.h>
#include <SP_Constraint.h>
#include <SP_ConstraintIter.h>
#include <MP_Constraint.h>
#include <MP_ConstraintIter.h>
#include <map>
#include <set>
#include <algorithm>
#include <Recorder.h>
#include <Pressure_Constraint.h>
#include <vector>
#include <Parameter.h>
#include <ParameterIter.h>
#include <DummyStream.h>
#include <Response.h>
#include <Mesh.h>
#include <BackgroundMesh.h>
#include <Parameter.h>
#include <ParameterIter.h>
#include <UniaxialMaterial.h>
#include <SectionForceDeformation.h>
#include <Damping.h>

void* OPS_NodeRecorder();
void* OPS_EnvelopeNodeRecorder();
void* OPS_ElementRecorder();
void* OPS_EnvelopeElementRecorder();
void* OPS_PVDRecorder();
void* OPS_AlgorithmRecorder();
void* OPS_RemoveRecorder();
#ifdef _HDF5
void* OPS_MPCORecorder();
void* OPS_LadrunoRecorder();
void* OPS_LadrunoMonitorRecorder();   // Ladruno: live analysis-monitor recorder
void* OPS_VTKHDF_Recorder();
#endif
void* OPS_GmshRecorder();
BackgroundMesh& OPS_getBgMesh();

void* OPS_DriftRecorder();
void* OPS_EnvelopeDriftRecorder();
void* OPS_EnergyBalanceRecorder();

int OPS_sectionLocation();
int OPS_sectionWeight();
int OPS_sectionTag();
int OPS_sectionDisplacement();

namespace {

    struct char_cmp {
        bool operator () (const char *a, const char *b) const
        {
            return strcmp(a, b)<0;
        }
    };

    typedef std::map<const char *, void *(*)(void), char_cmp> OPS_ParsingFunctionMap;

    static OPS_ParsingFunctionMap recordersMap;

    static int setUpRecorders(void) {
        recordersMap.insert(std::make_pair("Node", &OPS_NodeRecorder));
        recordersMap.insert(std::make_pair("EnvelopeNode", &OPS_EnvelopeNodeRecorder));
        recordersMap.insert(std::make_pair("Element", &OPS_ElementRecorder));
        recordersMap.insert(std::make_pair("EnvelopeElement", &OPS_EnvelopeElementRecorder));
	recordersMap.insert(std::make_pair("PVD", &OPS_PVDRecorder));
	recordersMap.insert(std::make_pair("BgPVD", &OPS_PVDRecorder));
	recordersMap.insert(std::make_pair("Remove", &OPS_RemoveRecorder));
	recordersMap.insert(std::make_pair("ElementRemoval", &OPS_RemoveRecorder));
	recordersMap.insert(std::make_pair("NodeRemoval", &OPS_RemoveRecorder));
	recordersMap.insert(std::make_pair("Collapse", &OPS_RemoveRecorder));
	recordersMap.insert(std::make_pair("Drift", &OPS_DriftRecorder));
	recordersMap.insert(std::make_pair("EnvelopeDrift", &OPS_EnvelopeDriftRecorder));
	recordersMap.insert(std::make_pair("gmsh", &OPS_GmshRecorder));
	recordersMap.insert(std::make_pair("EnergyBalance", &OPS_EnergyBalanceRecorder));
#ifdef _HDF5
	recordersMap.insert(std::make_pair("mpco", &OPS_MPCORecorder));
	recordersMap.insert(std::make_pair("ladruno", &OPS_LadrunoRecorder));
	recordersMap.insert(std::make_pair("Monitor", &OPS_LadrunoMonitorRecorder)); // Ladruno live monitor
    recordersMap.insert(std::make_pair("VTKHDF", &OPS_VTKHDF_Recorder));
#endif
        //recordersMap.insert(std::make_pair("Drift", &OPS_DriftRecorder));
        //recordersMap.insert(std::make_pair("Pattern", &OPS_PatternRecorder));

        return 0;
    }
}

int OPS_Recorder()
{
    static bool initDone = false;
    if (initDone == false) {
        setUpRecorders();
        initDone = true;
    }

    if (OPS_GetNumRemainingInputArgs() < 2) {
        opserr << "WARNING too few arguments: recorder type? tag? ...\n";
        return -1;
    }

    const char* type = OPS_GetString();
    OPS_ParsingFunctionMap::const_iterator iter = recordersMap.find(type);
    if (iter == recordersMap.end()) {
        opserr << "WARNING recorder type " << type << " is unknown\n";
        return -1;
    }

    Recorder* theRecorder = (Recorder*)(*iter->second)();
    if (theRecorder == 0) {
        opserr << "WARNING failed to create recorder\n";
        return -1;
    }

    if (strcmp(type,"BgPVD") == 0) {
	BackgroundMesh& bg = OPS_getBgMesh();
	bg.addRecorder(theRecorder);
    } else {

	// now add the recorder to the domain
	Domain* theDomain = OPS_GetDomain();
	if (theDomain == 0) return -1;

	if (theDomain->addRecorder(*theRecorder) < 0) {
	    opserr << "ERROR could not add to domain - recorder.\n";
	    delete theRecorder;
	    return -1;
	}
    }

    // set recorder tag as result
    int size = 1;
    int tag = theRecorder->getTag();
    if (OPS_SetIntOutput(&size, &tag, true) < 0) {
        opserr << "ERROR: failed to return recorder tag\n";
        return -1;
    }

    return 0;
}

int OPS_nodeDisp()
{
    if (OPS_GetNumRemainingInputArgs() < 1) {
	opserr << "WARNING insufficient args: nodeDisp nodeTag <dof ...>\n";
	return -1;
    }

    // tag and dof
    int data[2] = {0, -1};
    int numdata = OPS_GetNumRemainingInputArgs();
    if (numdata > 2) {
	numdata = 2;
    }

    if (OPS_GetIntInput(&numdata, data) < 0) {
	opserr<<"WARNING nodeDisp - failed to read int inputs\n";
	return -1;
    }
    data[1]--;

    // get response
    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;
    const Vector *nodalResponse = theDomain->getNodeResponse(data[0], Disp);

    if (nodalResponse == 0) {
	opserr << "WARNING no response is found\n";
	return -1;
    }

    // set outputs
    int size = nodalResponse->Size();
    if (data[1] >= 0) {
	if (data[1] >= size) {
	    opserr << "WARNING nodeDisp nodeTag? dof? - dofTag? too large\n";
	    return -1;
	}

	double value = (*nodalResponse)(data[1]);
	numdata = 1;

	if (OPS_SetDoubleOutput(&numdata, &value, true) < 0) {
	    opserr<<"WARNING nodeDisp - failed to read double inputs\n";
	    return -1;
	}


    } else {

	std::vector<double> values(size);
	for (int i=0; i<size; i++) {
	    values[i] = (*nodalResponse)(i);
	}
	if (OPS_SetDoubleOutput(&size, &values[0], false) < 0) {
	    opserr<<"WARNING nodeDisp - failed to read double inputs\n";
	    return -1;
	}
    }

    return 0;
}


// Ladruno (ADR-30 P3): ladrunoProjectionTieForce nodeTag? dof?
// Returns the constraint tie force f = M(a_raw - a_proj) at (node, dof) from the last
// projection step. Requires the active handler to be LadrunoProjection.
int OPS_LadrunoProjectionTieForce()
{
    if (OPS_GetNumRemainingInputArgs() < 2) {
	opserr << "WARNING want - ladrunoProjectionTieForce nodeTag? dof?\n";
	return -1;
    }
    int data[2];
    int numdata = 2;
    if (OPS_GetIntInput(&numdata, data) < 0) {
	opserr << "WARNING ladrunoProjectionTieForce - failed to read int inputs\n";
	return -1;
    }
    int nodeTag = data[0];
    int dof = data[1] - 1;          // 1-based (user) -> 0-based (local DOF index)

    ConstraintHandler **theHandler = OPS_GetHandler();
    LadrunoProjectionHandler *lph = 0;
    if (theHandler != 0 && *theHandler != 0)
	lph = dynamic_cast<LadrunoProjectionHandler *>(*theHandler);
    if (lph == 0) {
	opserr << "WARNING ladrunoProjectionTieForce - the active constraint handler is not "
		  "LadrunoProjection (use `constraints LadrunoProjection`).\n";
	return -1;
    }

    double f = lph->getTieForce(nodeTag, dof);
    int one = 1;
    if (OPS_SetDoubleOutput(&one, &f, true) < 0)
	return -1;
    return 0;
}


// ---------------------------------------------------------------------------
// Ladruno (ADR-39): ContactDomain commands. P1b = define surfaces + contacts +
// attach the engine; ZERO force (narrow phase is P2).
// ---------------------------------------------------------------------------
#include <LadrunoContactDomain.h>
#include <LadrunoContactSurface.h>

// lazily get/create the contact engine on the active Domain
static LadrunoContactDomain *OPS_getOrCreateContactDomain()
{
    Domain *theDomain = OPS_GetDomain();
    if (theDomain == 0) return 0;
    LadrunoContactDomain *cd = theDomain->getLadrunoContactDomain();
    if (cd == 0) {
        cd = new LadrunoContactDomain();
        theDomain->setLadrunoContactDomain(cd);
    }
    return cd;
}

// contactSurface tag (-slave | -master nodesPerSeg | -slave-segments nodesPerSeg) nodeTag...
int OPS_LadrunoContactSurface()
{
    if (OPS_GetNumRemainingInputArgs() < 3) {
        opserr << "WARNING want - contactSurface tag (-slave | -master nodesPerSeg | "
                  "-slave-segments nodesPerSeg) nodeTag...\n";
        return -1;
    }
    int tag, n = 1;
    if (OPS_GetIntInput(&n, &tag) < 0) {
        opserr << "WARNING contactSurface - could not read tag\n";
        return -1;
    }
    const char *kindStr = OPS_GetString();
    LadrunoContactSurface::Kind kind;
    int nodesPerSeg = 0;
    if (strcmp(kindStr, "-slave") == 0) {
        kind = LadrunoContactSurface::SLAVE_NODES;
    } else if (strcmp(kindStr, "-master") == 0) {
        kind = LadrunoContactSurface::MASTER_SEGMENTS;
        n = 1;
        if (OPS_GetIntInput(&n, &nodesPerSeg) < 0 || nodesPerSeg < 3) {
            opserr << "WARNING contactSurface -master - need nodesPerSeg (>=3)\n";
            return -1;
        }
    } else if (strcmp(kindStr, "-slave-segments") == 0) {
        // Ladruno ADR-41 C2: a FACETED slave surface (same layout as -master, slave side).
        kind = LadrunoContactSurface::SLAVE_SEGMENTS;
        n = 1;
        if (OPS_GetIntInput(&n, &nodesPerSeg) < 0 || nodesPerSeg < 3) {
            opserr << "WARNING contactSurface -slave-segments - need nodesPerSeg (>=3)\n";
            return -1;
        }
    } else {
        opserr << "WARNING contactSurface - kind must be -slave, -master or -slave-segments\n";
        return -1;
    }
    int nrem = OPS_GetNumRemainingInputArgs();
    if (nrem <= 0) {
        opserr << "WARNING contactSurface - no node tags given\n";
        return -1;
    }
    ID nodeTags(nrem);
    for (int i = 0; i < nrem; i++) {
        int v, one = 1;
        if (OPS_GetIntInput(&one, &v) < 0) {
            opserr << "WARNING contactSurface - could not read node tag\n";
            return -1;
        }
        nodeTags(i) = v;
    }
    LadrunoContactDomain *cd = OPS_getOrCreateContactDomain();
    if (cd == 0) return -1;
    LadrunoContactSurface *s = new LadrunoContactSurface(tag, kind, nodeTags, nodesPerSeg);
    if (cd->addSurface(s) < 0) { delete s; return -1; }
    return 0;
}

// contact tag masterSurfTag slaveSurfTag <kn kt mu> <-outward ox oy oz>
int OPS_LadrunoContact()
{
    if (OPS_GetNumRemainingInputArgs() < 3) {
        opserr << "WARNING want - contact tag masterSurfTag slaveSurfTag <kn kt mu> <-outward ox oy oz>\n";
        return -1;
    }
    int idata[3], n = 3;
    if (OPS_GetIntInput(&n, idata) < 0) {
        opserr << "WARNING contact - could not read tag/master/slave\n";
        return -1;
    }
    double kn = 0.0, kt = 0.0, mu = 0.0;   // P1b zero-force defaults; P2 parses/uses
    bool knAuto = false;                   // Ladruno ADR-39 P2b-2b: `auto` -> size kn from master ele
    // The kn slot accepts a number ("kn kt mu") OR the literal `auto`. Either may be
    // omitted entirely (then the next token is the -outward flag, or end-of-args).
    if (OPS_GetNumRemainingInputArgs() >= 1) {
        const char *peek = OPS_GetString();   // consume + classify
        // ANY option flag (-outward / -cell / -consistanttan / -mortar / -epsN / ...)
        // means the numeric kn/kt/mu slot was omitted ⇒ un-read and let the options loop
        // handle it. (Was: only `-outward` was recognized here — ADR-41 C2 generalized it.)
        bool isFlag = (peek != 0 && peek[0] == '-');
        bool isAuto = (peek != 0 && strcmp(peek, "auto") == 0);
        if (isAuto) {
            knAuto = true;                    // 'auto' consumed; kn is resolved at handle()
            // optional kt mu (friction, P3) may follow unless the next token is a flag.
            if (OPS_GetNumRemainingInputArgs() >= 2) {
                const char *p2 = OPS_GetString();
                bool flag2 = (p2 != 0 && p2[0] == '-');
                OPS_ResetCurrentInputArg(-1);
                if (!flag2) {
                    double d[2]; int m = 2;
                    if (OPS_GetDoubleInput(&m, d) < 0) {
                        opserr << "WARNING contact - could not read kt mu after 'auto'\n";
                        return -1;
                    }
                    kt = d[0]; mu = d[1];
                }
            }
        } else if (isFlag) {
            OPS_ResetCurrentInputArg(-1);     // un-read for the option-flag loop below
        } else {
            OPS_ResetCurrentInputArg(-1);     // un-read; read the numeric kn (kt mu)
            // existing form is exactly `kn kt mu`; tolerate a bare `kn` too.
            int avail = OPS_GetNumRemainingInputArgs();
            int m = (avail >= 3) ? 3 : 1;
            double d[3] = {0.0, 0.0, 0.0};
            if (OPS_GetDoubleInput(&m, d) < 0) {
                opserr << "WARNING contact - could not read kn kt mu\n";
                return -1;
            }
            kn = d[0]; kt = d[1]; mu = d[2];
        }
    }
    // optional -outward ox oy oz : orientation direction toward the allowed half-space
    // optional -cell <frac>      : P2.5 bucket-sort cell = frac * median seg diagonal
    bool hasOutward = false;
    double outward[3] = {0.0, 0.0, 0.0};
    double cellFrac = 1.0;
    // Ladruno ADR-60: finite-sliding NTS re-emit (off ⇒ byte-identical). -reemit opts in; -resortFrac
    // sets the drift fraction of the search band that triggers a re-sort; -resortEvery sets the min
    // commits between re-sorts (0 ⇒ the default floor). NTS lane only (refused with -mortar below).
    bool   enableReemit = false;
    double resortFrac   = 0.5;
    int    resortEvery  = 0;
    bool   smoothNormal = false;  // Ladruno ADR-63 #4a: averaged nodal-normal smoothing (off ⇒ identical)
    bool consistentTan = false;   // Ladruno ADR-39 P3.5: friction tangent symmetry
    // Ladruno ADR-39 B3 (P2b-2c): `-geomtan` opts the NTS SEGMENT lane into the consistent
    // ∂n/∂u geometric NORMAL tangent (kn·gN·∂²gN/∂u²) ⇒ quadratic Newton on CURVED / large-
    // sliding interfaces. SYMMETRIC ⇒ solver-safe on any system (unlike the friction Csl). Off
    // (default) ⇒ the shipped kn·BᵀB main term (byte-identical; EXACT for a flat/fixed master).
    bool consistentNormal = false;
    // Ladruno ADR-41 C2: the mortar lane. `-mortar` selects the clipped-GP mortar +
    // frictionless commit-cycle ALM formulation; -epsN/-augTol/-maxAug/-ngp tune it.
    bool isMortar = false;
    double epsN = 0.0;
    bool epsNAuto = false;
    double augTol = 1e-8;
    int maxAug = 20, ngp = 2;
    // ADR-41 C3.1 mortar friction (Coulomb/Tresca unified cone min(μN+c, τmax)); all 0 ⇒ frictionless.
    double mortarMu = 0.0, epsT = 0.0, cohesion = 0.0, tauMax = 0.0;
    bool epsTAuto = false;
    // ADR-41 C4 mesh-tying: `-tie` makes the mortar pair a PERMANENT bond (the zero-gap limit — full
    // 3-vec r→0, no clamp, no friction). `-epsTie auto|val` is the tie penalty (an alias for the
    // -epsN penalty slot — a tie has one penalty). Mutually exclusive with -mu/-cohesion/-tauMax.
    bool isTie = false;
    // ADR-41 D2 viscous stabilization: `-visc <μ_c>` adds a velocity-proportional normal contact
    // damper (p_visc = μ_c·gap_rate) that bleeds chatter/snap-through energy in the pounding/rocking/
    // uplift regime. 0 (default) ⇒ no viscous term (byte-identical). Works on NTS (D2.1) AND mortar
    // CONTACT (D2.2); refused with -tie (a bond has no chatter). Naturally inert in statics (v ≡ 0).
    double muc = 0.0;
    // Ladruno ADR-39 B1/B2 (P4/P5): `-soft <SOFSCL>` opts into the LS-DYNA §26.15 SOFT Courant-stable
    // explicit penalty. The FORMULATION sets SOFT=1 vs SOFT=2: WITHOUT -mortar it is B1 SOFT=1 (NODE-
    // to-segment); WITH -mortar it is B2 SOFT=2 (SEGMENT-based — integrates the clipped facet overlap,
    // catching the corner/edge/T-intersection cases NTS misses). Both size the contact stiffness from
    // the nodal MASS + the timestep (k_soft = SOFSCL·4·m_eff/dt²) instead of the material, so the
    // contact never throttles the explicit dt_cr (impact/pounding/recontact runs at the STRUCTURAL dt).
    // A modifier on the penalty: still needs a base penalty (-kn/auto for NTS; -epsN/auto for mortar) —
    // that base is what an IMPLICIT run uses (SOFT is explicit-only ⇒ implicit byte-identical). SOFSCL
    // optional (default 0.10, the LS-DYNA SLSFAC default); ≤0 ⇒ off; >1 warns (ω·dt = 2√SOFSCL > 2).
    double softScale = 0.0;
    // Ladruno ADR-57 E2: `-edgeedge` opts a -mortar contact into the perpendicular edge-edge
    // fallback (the cos_t→0 pairs the face-mortar clip degenerates on get a dedicated segment-
    // segment penalty). OFF by default ⇒ byte-identical. `-edgeKn auto|<val>` sets the edge penalty
    // (default = the mortar epsN); `-edgeBand <d>` sets the gap activation band (default from the
    // facet edge length). Requires -mortar (validated after the loop).
    // E3: `-edgeMu`/`-edgeKt`/`-edgeCohesion`/`-edgeTauMax` add edge-edge Coulomb/Tresca friction
    // (the unified cone min(μN+c, τmax)); `-edgeConsistentTan` opts into the non-symmetric Csl tangent.
    bool edgeEdge = false, edgeKnAuto = false;
    double edgeKn = 0.0, edgeBand = 0.0;
    double edgeMu = 0.0, edgeKt = 0.0, edgeCohesion = 0.0, edgeTauMax = 0.0;
    bool edgeConsistentTan = false;
    // Ladruno ADR-57 E5: `-edgeSoft <SOFSCL>` opts the edge-edge fallback into the explicit Courant-stable
    // SOFT penalty (the B1/B2 SOFT analogue for the edge operator). Under the explicit CentralDifferenceLadruno
    // the edge penalty is replaced by k_soft = SOFSCL·4·m_eff/dt² so edge-on impact runs at the structural
    // dt_cr. SOFSCL optional (default 0.10); ≤0 ⇒ off; inert under implicit ⇒ byte-identical. Needs -edgeedge.
    double edgeSoftScale = 0.0;
    // Ladruno ADR-57 E6: `-edgeAlm` opts the edge-edge fallback into the one-scalar commit-cycle ALM
    // (a per-pair λ_N drives penetration → an εN-independent tol; the held-load analyze_augmented proc
    // reads ladrunoEdgePenetration). OFF by default ⇒ the E2 penalty path ⇒ byte-identical. Implicit-only.
    // `-edgeAugTol <tol>` is the augmentation tolerance (metadata; the proc passes its own augTol). Needs -edgeedge.
    bool edgeAlm = false;
    double edgeAugTol = 0.0;
    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char *opt = OPS_GetString();
        if (opt != 0 && strcmp(opt, "-mortar") == 0) {
            isMortar = true;
        } else if (opt != 0 && strcmp(opt, "-soft") == 0) {
            // optional numeric SOFSCL; default 0.10 if the next token is a flag / end of args.
            softScale = 0.10;
            if (OPS_GetNumRemainingInputArgs() > 0) {
                const char *p = OPS_GetString();          // string read consumes on Tcl AND Py
                bool isFlag = (p != 0 && p[0] == '-');
                OPS_ResetCurrentInputArg(-1);
                if (!isFlag) {
                    double v[1]; int m = 1;
                    if (OPS_GetDoubleInput(&m, v) < 0) {
                        opserr << "WARNING contact -soft - need a SOFSCL value or omit it (default 0.1)\n";
                        return -1;
                    }
                    softScale = v[0];
                }
            }
            if (softScale > 1.0)
                opserr << "WARNING contact -soft SOFSCL=" << softScale << " > 1: the contact mode "
                          "ω·dt = 2√SOFSCL > 2 is UNSTABLE under central difference; use SOFSCL ≤ 1 "
                          "(default 0.1).\n";
        } else if (opt != 0 && strcmp(opt, "-tie") == 0) {
            // ADR-41 C4: a PERMANENT mesh-tie bond (requires -mortar; validated after the loop).
            isTie = true;
        } else if (opt != 0 && strcmp(opt, "-visc") == 0) {
            // ADR-41 D2: viscous normal-stabilization coefficient μ_c (p_visc = μ_c·gap_rate).
            double v[1]; int m = 1;
            if (OPS_GetDoubleInput(&m, v) < 0) {
                opserr << "WARNING contact -visc - need a coefficient value\n";
                return -1;
            }
            muc = v[0];
        } else if (opt != 0 && strcmp(opt, "-epsTie") == 0) {
            // -epsTie auto | <value> : the tie penalty (an alias for the -epsN penalty slot — a tie
            // has a single penalty; auto ⇒ sized from the owning solid, like the contact penalty).
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "WARNING contact -epsTie - need auto or a value\n";
                return -1;
            }
            const char *p = OPS_GetString();
            if (p != 0 && strcmp(p, "auto") == 0) {
                epsNAuto = true;
            } else {
                OPS_ResetCurrentInputArg(-1);
                double v[1]; int m = 1;
                if (OPS_GetDoubleInput(&m, v) < 0) {
                    opserr << "WARNING contact -epsTie - need auto or a value\n";
                    return -1;
                }
                epsN = v[0];
            }
        } else if (opt != 0 && strcmp(opt, "-epsN") == 0) {
            // -epsN auto | <value> : ALM normal penalty (auto => sized like -kn auto).
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "WARNING contact -epsN - need auto or a value\n";
                return -1;
            }
            const char *e = OPS_GetString();
            if (e != 0 && strcmp(e, "auto") == 0) {
                epsNAuto = true;
            } else {
                OPS_ResetCurrentInputArg(-1);
                double v[1]; int m = 1;
                if (OPS_GetDoubleInput(&m, v) < 0) {
                    opserr << "WARNING contact -epsN - need auto or a value\n";
                    return -1;
                }
                epsN = v[0];
            }
        } else if (opt != 0 && strcmp(opt, "-augTol") == 0) {
            double v[1]; int m = 1;
            if (OPS_GetDoubleInput(&m, v) < 0) {
                opserr << "WARNING contact -augTol - need a value\n";
                return -1;
            }
            augTol = v[0];
        } else if (opt != 0 && strcmp(opt, "-maxAug") == 0) {
            int v[1]; int m = 1;
            if (OPS_GetIntInput(&m, v) < 0) {
                opserr << "WARNING contact -maxAug - need an integer\n";
                return -1;
            }
            maxAug = v[0];
        } else if (opt != 0 && strcmp(opt, "-ngp") == 0) {
            int v[1]; int m = 1;
            if (OPS_GetIntInput(&m, v) < 0) {
                opserr << "WARNING contact -ngp - need an integer\n";
                return -1;
            }
            ngp = v[0];
        } else if (opt != 0 && strcmp(opt, "-mu") == 0) {
            // ADR-41 C3.1: Coulomb friction coefficient on the mortar interface.
            double v[1]; int m = 1;
            if (OPS_GetDoubleInput(&m, v) < 0) {
                opserr << "WARNING contact -mu - need a value\n";
                return -1;
            }
            mortarMu = v[0];
        } else if (opt != 0 && strcmp(opt, "-epsT") == 0) {
            // -epsT auto | <value> : tangential penalty (auto => = the normal penalty epsN).
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "WARNING contact -epsT - need auto or a value\n";
                return -1;
            }
            const char *e = OPS_GetString();
            if (e != 0 && strcmp(e, "auto") == 0) {
                epsTAuto = true;
            } else {
                OPS_ResetCurrentInputArg(-1);
                double v[1]; int m = 1;
                if (OPS_GetDoubleInput(&m, v) < 0) {
                    opserr << "WARNING contact -epsT - need auto or a value\n";
                    return -1;
                }
                epsT = v[0];
            }
        } else if (opt != 0 && strcmp(opt, "-cohesion") == 0) {
            double v[1]; int m = 1;
            if (OPS_GetDoubleInput(&m, v) < 0) {
                opserr << "WARNING contact -cohesion - need a value\n";
                return -1;
            }
            cohesion = v[0];
        } else if (opt != 0 && strcmp(opt, "-tauMax") == 0) {
            double v[1]; int m = 1;
            if (OPS_GetDoubleInput(&m, v) < 0) {
                opserr << "WARNING contact -tauMax - need a value\n";
                return -1;
            }
            tauMax = v[0];
        } else if (opt != 0 && strcmp(opt, "-edgeedge") == 0) {
            // ADR-57 E2: enable the perpendicular edge-edge fallback (requires -mortar).
            edgeEdge = true;
        } else if (opt != 0 && strcmp(opt, "-edgeKn") == 0) {
            // -edgeKn auto | <value> : the edge-edge penalty (auto ⇒ sized per master facet like
            // -epsN auto; a value ⇒ fixed; omitted ⇒ defaults to the resolved mortar penalty).
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "WARNING contact -edgeKn - need auto or a value\n";
                return -1;
            }
            const char *e = OPS_GetString();
            if (e != 0 && strcmp(e, "auto") == 0) {
                edgeKnAuto = true;
            } else {
                OPS_ResetCurrentInputArg(-1);
                double v[1]; int m = 1;
                if (OPS_GetDoubleInput(&m, v) < 0) {
                    opserr << "WARNING contact -edgeKn - need auto or a value\n";
                    return -1;
                }
                edgeKn = v[0];
            }
        } else if (opt != 0 && strcmp(opt, "-edgeBand") == 0) {
            // -edgeBand <d> : the gap activation band d_band for the edge-edge fallback
            // (default ⇒ sized from the facet edge length at handle() time).
            double v[1]; int m = 1;
            if (OPS_GetDoubleInput(&m, v) < 0) {
                opserr << "WARNING contact -edgeBand - need a value\n";
                return -1;
            }
            edgeBand = v[0];
        } else if (opt != 0 && strcmp(opt, "-edgeMu") == 0) {
            double v[1]; int m = 1;
            if (OPS_GetDoubleInput(&m, v) < 0) { opserr << "WARNING contact -edgeMu - need a value\n"; return -1; }
            edgeMu = v[0];
        } else if (opt != 0 && strcmp(opt, "-edgeKt") == 0) {
            double v[1]; int m = 1;
            if (OPS_GetDoubleInput(&m, v) < 0) { opserr << "WARNING contact -edgeKt - need a value\n"; return -1; }
            edgeKt = v[0];
        } else if (opt != 0 && strcmp(opt, "-edgeCohesion") == 0) {
            double v[1]; int m = 1;
            if (OPS_GetDoubleInput(&m, v) < 0) { opserr << "WARNING contact -edgeCohesion - need a value\n"; return -1; }
            edgeCohesion = v[0];
        } else if (opt != 0 && strcmp(opt, "-edgeTauMax") == 0) {
            double v[1]; int m = 1;
            if (OPS_GetDoubleInput(&m, v) < 0) { opserr << "WARNING contact -edgeTauMax - need a value\n"; return -1; }
            edgeTauMax = v[0];
        } else if (opt != 0 && strcmp(opt, "-edgeConsistentTan") == 0) {
            edgeConsistentTan = true;
        } else if (opt != 0 && strcmp(opt, "-edgeSoft") == 0) {
            // ADR-57 E5: optional numeric SOFSCL; default 0.10 if the next token is a flag / end of args
            // (mirrors -soft). Under explicit CDL the edge penalty becomes k_soft = SOFSCL·4·m_eff/dt².
            edgeSoftScale = 0.10;
            if (OPS_GetNumRemainingInputArgs() > 0) {
                const char *p = OPS_GetString();          // string read consumes on Tcl AND Py
                bool isFlag = (p != 0 && p[0] == '-');
                OPS_ResetCurrentInputArg(-1);
                if (!isFlag) {
                    double v[1]; int m = 1;
                    if (OPS_GetDoubleInput(&m, v) < 0) {
                        opserr << "WARNING contact -edgeSoft - need a SOFSCL value or omit it (default 0.1)\n";
                        return -1;
                    }
                    edgeSoftScale = v[0];
                }
            }
            if (edgeSoftScale > 1.0)
                opserr << "WARNING contact -edgeSoft SOFSCL=" << edgeSoftScale << " > 1: the edge mode "
                          "ω·dt = 2√SOFSCL > 2 is UNSTABLE under central difference; use SOFSCL ≤ 1 "
                          "(default 0.1).\n";
        } else if (opt != 0 && strcmp(opt, "-edgeAlm") == 0) {
            // ADR-57 E6: enable the one-scalar commit-cycle ALM on the edge-edge fallback (a flag).
            edgeAlm = true;
        } else if (opt != 0 && strcmp(opt, "-edgeAugTol") == 0) {
            // ADR-57 E6: the edge-edge ALM augmentation tolerance (metadata; the held-load proc passes it).
            double v[1]; int m = 1;
            if (OPS_GetDoubleInput(&m, v) < 0) { opserr << "WARNING contact -edgeAugTol - need a value\n"; return -1; }
            edgeAugTol = v[0];
        } else if (opt != 0 && strcmp(opt, "-outward") == 0) {
            double o[3]; int m = 3;
            if (OPS_GetDoubleInput(&m, o) < 0) {
                opserr << "WARNING contact -outward - need ox oy oz\n";
                return -1;
            }
            outward[0] = o[0]; outward[1] = o[1]; outward[2] = o[2];
            hasOutward = true;
        } else if (opt != 0 && strcmp(opt, "-cell") == 0) {
            // Ladruno ADR-39 P2.5: broad-phase cell-size scale (median seg diag).
            // A huge value => 1 bucket => every segment a candidate = brute force.
            double f[1]; int m = 1;
            if (OPS_GetDoubleInput(&m, f) < 0) {
                opserr << "WARNING contact -cell - need a positive frac\n";
                return -1;
            }
            cellFrac = f[0];
        } else if (opt != 0 && strcmp(opt, "-consistanttan") == 0) {
            // Ladruno ADR-39 P3.5: opt in to the NON-SYMMETRIC consistent friction
            // tangent (true quadratic Newton convergence). It REQUIRES a non-symmetric
            // linear solver — `system FullGeneral` / `UmfPack` / `BandGeneral`. A
            // symmetric solver (ProfileSPD/BandSPD/SparseSYM) reads only the upper
            // triangle and SILENTLY drops the d_TN coupling, corrupting the solve. The
            // default (symmetric) tangent is correct on any solver (design-gate Q2).
            consistentTan = true;
            opserr << "WARNING contact -consistanttan: the non-symmetric consistent "
                      "friction tangent needs a non-symmetric solver (system FullGeneral/"
                      "UmfPack/BandGeneral); symmetric solvers will silently corrupt it.\n";
        } else if (opt != 0 && strcmp(opt, "-geomtan") == 0) {
            // Ladruno ADR-39 B3 (P2b-2c): opt the NTS SEGMENT lane into the consistent ∂n/∂u
            // geometric NORMAL tangent ⇒ quadratic Newton on CURVED / large-sliding interfaces
            // (the Hertz benchmark). SYMMETRIC ⇒ correct on ANY solver (no -consistanttan needed).
            // Off ⇒ the shipped kn·BᵀB main term (byte-identical; EXACT for a flat/fixed master).
            consistentNormal = true;
        } else if (opt != 0 && strcmp(opt, "-reemit") == 0) {
            // Ladruno ADR-60: opt the NTS contact into finite-sliding re-emit (deformed-config
            // broad-phase re-sort when a slave migrates off its reference-config candidate
            // segments — the silent pass-through fix). Off ⇒ byte-identical.
            enableReemit = true;
        } else if (opt != 0 && strcmp(opt, "-resortFrac") == 0) {
            double v[1]; int m = 1;
            if (OPS_GetDoubleInput(&m, v) < 0) {
                opserr << "WARNING contact -resortFrac - need a positive fraction of the search band\n";
                return -1;
            }
            resortFrac = v[0];
        } else if (opt != 0 && strcmp(opt, "-resortEvery") == 0) {
            int v[1]; int m = 1;
            if (OPS_GetIntInput(&m, v) < 0) {
                opserr << "WARNING contact -resortEvery - need an integer (min commits between re-sorts)\n";
                return -1;
            }
            resortEvery = v[0];
        } else if (opt != 0 && strcmp(opt, "-smoothNormal") == 0) {
            // Ladruno ADR-63 #4a: opt the NTS contact into averaged nodal-normal smoothing — a smooth
            // N(X) master-normal field (coherent winding + a GLOBAL outward sign) that resolves the
            // ADR-60 R3 curved-master ridge flip + the ADR-41 Q-NORMAL junction chatter. Off ⇒ the
            // faceted normalOriented() path ⇒ byte-identical. NTS-only (refused with -mortar below).
            smoothNormal = true;
        } else {
            // Ladruno ADR-39 P2b-2b (gate MINOR-1): error on an unexpected trailing
            // token rather than silently swallowing it (e.g. a stray friction value
            // after `auto`, or a mistyped `-outwards`) — silent acceptance hid input
            // mistakes. Recognized forms: `kn kt mu`/`auto` then optional `-outward`.
            opserr << "WARNING contact - unexpected token '" << (opt ? opt : "")
                   << "' (expected -outward or end of arguments)\n";
            return -1;
        }
    }
    Domain *theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;
    LadrunoContactDomain *cd = theDomain->getLadrunoContactDomain();
    if (cd == 0) {
        opserr << "WARNING contact - define a contactSurface first\n";
        return -1;
    }
    if (isTie && !isMortar) {
        opserr << "WARNING contact -tie requires -mortar (mesh-tying is a mortar formulation)\n";
        return -1;
    }
    if (isTie && (mortarMu > 0.0 || cohesion > 0.0 || tauMax > 0.0)) {
        // A tie is an EQUALITY bond — it has no friction cone. Refuse the combination explicitly.
        opserr << "WARNING contact -tie is mutually exclusive with friction "
                  "(-mu/-cohesion/-tauMax): a mesh-tie has no friction cone\n";
        return -1;
    }
    // Ladruno contact-review P3: a τmax-only friction config (no -mu, no -cohesion) is the
    // unified cone min(μN+c, τmax) = min(0, τmax) = 0 ⇒ FREE SLIP — almost certainly not what
    // the user wanted (a shear-capped BOND is -cohesion, optionally capped by -tauMax). The old
    // kernel silently turned this config into an UNBOUNDED elastic bond; refuse it loudly.
    if (tauMax > 0.0 && mortarMu <= 0.0 && cohesion <= 0.0) {
        opserr << "WARNING contact -tauMax without -mu or -cohesion: the friction cone "
                  "min(mu*N + c, tauMax) = min(0, tauMax) = 0 (frictionless free slip). For a "
                  "shear-capped bond use -cohesion <c> (a Tresca cap), optionally with -tauMax\n";
        return -1;
    }
    if (edgeEdge && edgeTauMax > 0.0 && edgeMu <= 0.0 && edgeCohesion <= 0.0) {
        opserr << "WARNING contact -edgeTauMax without -edgeMu or -edgeCohesion: the edge friction "
                  "cone min(mu*N + c, tauMax) = min(0, tauMax) = 0 (frictionless free slip). For a "
                  "shear-capped bond use -edgeCohesion <c>, optionally with -edgeTauMax\n";
        return -1;
    }
    if (muc > 0.0 && isTie) {
        // ADR-41 D2.2: a mesh-tie is a permanent bond (no contact-status flips ⇒ no chatter to damp).
        // Refuse -visc with -tie rather than silently ignore it.
        opserr << "WARNING contact -visc is not allowed with -tie (a bond has no contact-chatter "
                  "regime); drop -tie for a viscous-stabilized mortar CONTACT\n";
        return -1;
    }
    if (edgeEdge && !isMortar) {
        // ADR-57 E2: the edge-edge fallback is a modifier on a -mortar contact (the surfaces are
        // already declared faceted master/slave; the cos_t→0 pairs are routed off the mortar lane).
        opserr << "WARNING contact -edgeedge requires -mortar (the perpendicular edge-edge fallback "
                  "is a mortar-lane modifier)\n";
        return -1;
    }
    if (!edgeEdge && (edgeKn > 0.0 || edgeKnAuto || edgeBand > 0.0 || edgeMu > 0.0 ||
                      edgeKt > 0.0 || edgeCohesion > 0.0 || edgeTauMax > 0.0 || edgeConsistentTan ||
                      edgeSoftScale > 0.0 || edgeAlm || edgeAugTol > 0.0))
        opserr << "WARNING contact -edgeKn/-edgeBand/-edgeMu/-edgeSoft/-edgeAlm/... given without "
                  "-edgeedge; ignored (enable the edge-edge fallback with -edgeedge)\n";
    if (edgeEdge && edgeConsistentTan && (edgeMu > 0.0 || edgeCohesion > 0.0 || edgeTauMax > 0.0))
        // the non-symmetric Coulomb Csl tangent needs a non-symmetric solver (FullGeneral/UmfPack/
        // BandGeneral); a symmetric SOE silently drops the lower triangle. Warn once (like -consistanttan).
        opserr << "WARNING contact -edgeConsistentTan: the non-symmetric edge-edge Coulomb friction "
                  "tangent needs a non-symmetric solver (system FullGeneral/UmfPack/BandGeneral); "
                  "symmetric solvers will silently corrupt it.\n";
    if (consistentNormal && isMortar) {
        // B3 (P2b-2c) is the NTS SEGMENT geometric tangent; the mortar lane's geometric
        // ∂{D,M,n}/∂u block is a SEPARATE deferral (C2 shipped the penalty Gram only). Refuse
        // rather than silently ignore -geomtan on a mortar contact.
        opserr << "WARNING contact -geomtan is an NTS (node-to-segment) option; it does not "
                  "apply to -mortar (the mortar geometric tangent is separately deferred)\n";
        return -1;
    }
    if (enableReemit && isMortar) {
        // ADR-60: re-emit is an NTS feature. The mortar lane is brute-force (every master facet is
        // already a candidate per slave facet), so it handles finite sliding without re-emission.
        // Refuse rather than silently ignore -reemit on a -mortar contact.
        opserr << "WARNING contact -reemit is an NTS (node-to-segment) option; -mortar is brute-force "
                  "and already finite-sliding-correct (re-emit not needed)\n";
        return -1;
    }
    if (smoothNormal && isMortar) {
        // ADR-63 #4a: nodal-normal smoothing is wired to the NTS lane only (it resolves the NTS R3
        // ridge flip + chatter). The mortar lane derives its GP normal from its own facet clip — a
        // SEPARATE Q-NORMAL instance, a future pass (the header is reusable). Refuse rather than
        // silently ignore -smoothNormal on a -mortar contact.
        opserr << "WARNING contact -smoothNormal is an NTS (node-to-segment) option; it does not apply "
                  "to -mortar (the mortar GP-normal smoothing is a separate, deferred pass)\n";
        return -1;
    }
    if (softScale > 0.0 && isMortar) {
        // B2 (P5): `-mortar -soft <SOFSCL>` selects the SOFT=2 SEGMENT-BASED explicit penalty — the
        // segment-to-segment generalization of the NTS SOFT=1 lane that catches the corner/edge/
        // T-intersection cases NTS misses (it integrates the clipped facet overlap, not slave NODES).
        // -visc μ_c (D2.2 normal damper) AND -mu/-cohesion/-tauMax (Courant-stable Coulomb/Tresca over
        // the segment overlap) are both supported on the SOFT=2 active set under explicit; under
        // implicit they fall through to the regular mortar penalty/friction (-soft is explicit-only).
        // Only -tie is REFUSED (a permanent bond is not a soft contact penalty).
        if (isTie) {
            opserr << "WARNING contact -soft (SOFT=2 segment-based explicit penalty) is not allowed "
                      "with -tie (a permanent bond is not a soft contact penalty)\n";
            return -1;
        }
        // SOFT=2 still needs a base penalty (the implicit fall-through + byte-identity contract use
        // it): -epsN auto|<val>, -epsTie (sets epsNAuto/epsN), or a positional auto|<kn>.
        if (!epsNAuto && !knAuto && epsN <= 0.0 && kn <= 0.0) {
            opserr << "WARNING contact -mortar -soft needs a base penalty (-epsN auto|<val> or a "
                      "positional auto|<kn>): SOFT=2 sizes k_soft under explicit, implicit uses it\n";
            return -1;
        }
        // B2 coupled-stability guard: SOFT=2 sizes a PER-NODE Courant penalty (ω·dt = 2√SOFSCL ≤ 2
        // for SOFSCL ≤ 1), but the ASSEMBLED segment stiffness K_c = Σ_I k_soft,I B_Iᵀ B_I couples
        // shared facet nodes, raising the true central-difference limit ~2× — the coupled step can
        // go unstable near SOFSCL ≈ 0.3 (oracle proto_b2_soft2_segment.py T4 / [[LEDGER_quirks]]).
        // The per-node `>1` warning below understates this for the segment lane, so warn earlier.
        if (softScale > 0.25)
            opserr << "WARNING contact -mortar -soft SOFSCL=" << softScale << ": SOFT=2's per-node "
                      "Courant bound is necessary-not-sufficient — inter-node coupling in the assembled "
                      "segment stiffness can make the central-difference step unstable near SOFSCL≈0.3. "
                      "Use SOFSCL ≤ 0.25 (default 0.1) unless a finer dt margin is verified.\n";
    } else if (softScale > 0.0 && !knAuto && kn <= 0.0) {
        // B1 (P4): NTS SOFT=1 is a MODIFIER on the penalty — under explicit it sizes k_soft, but an
        // implicit run (and the byte-identity contract) needs a real base kn. Refuse -soft with none.
        opserr << "WARNING contact -soft needs a base penalty: give a positional `auto` or `<kn>` "
                  "before the options (SOFT=1 sizes k_soft under explicit, implicit uses the base kn)\n";
        return -1;
    }
    if (isMortar) {
        // Ladruno ADR-41 C2.0/C2.2/C3.1: the mortar definition (normal ALM + C3.1 Coulomb/Tresca
        // friction via -mu/-epsT/-cohesion/-tauMax). friction params ≤0 ⇒ the frictionless C2 path.
        // C4: -tie ⇒ a permanent mesh-tie bond (full 3-vec r→0; friction refused above).
        // D2.2: -visc μ_c ⇒ viscous normal stabilization on the mortar contact (refused with -tie above).
        // B2: softScale>0 ⇒ the SOFT=2 segment-based explicit penalty (off ⇒ byte-identical mortar).
        // ADR-57 E2/E3: edgeEdge ⇒ the perpendicular edge-edge fallback (+ E3 friction; off ⇒ byte-identical).
        // ADR-57 E5: edgeSoftScale>0 ⇒ the explicit Courant-stable SOFT penalty on the edge fallback.
        // ADR-57 E6: edgeAlm ⇒ the one-scalar commit-cycle ALM (off ⇒ the E2 penalty path).
        return cd->addMortarContact(idata[0], idata[1], idata[2], kn, knAuto, epsN, epsNAuto,
                                    augTol, maxAug, ngp, hasOutward ? outward : 0, cellFrac,
                                    mortarMu, epsT, epsTAuto, cohesion, tauMax, consistentTan, isTie,
                                    muc, softScale, edgeEdge, edgeKn, edgeKnAuto, edgeBand,
                                    edgeMu, edgeKt, edgeCohesion, edgeTauMax, edgeConsistentTan,
                                    edgeSoftScale, edgeAlm, edgeAugTol);
    }
    // D2: -visc μ_c (NTS viscous normal stabilization; 0 ⇒ off, byte-identical).
    // B3: -geomtan ⇒ the consistent ∂n/∂u geometric normal tangent (off ⇒ byte-identical).
    // B1: -soft SOFSCL ⇒ the explicit SOFT=1 Courant-stable penalty (off ⇒ byte-identical).
    return cd->addContact(idata[0], idata[1], idata[2], kn, kt, mu,
                          hasOutward ? outward : 0, knAuto, cellFrac, consistentTan, muc,
                          consistentNormal, softScale, enableReemit, resortFrac, resortEvery,
                          smoothNormal);
}

// contactPlane tag slaveSurfTag  nx ny nz  px py pz  kn  <-visc μ_c> <-soft SOFSCL>  (P2a rigid plane)
int OPS_LadrunoContactPlane()
{
    if (OPS_GetNumRemainingInputArgs() < 9) {
        opserr << "WARNING want - contactPlane tag slaveSurfTag nx ny nz px py pz kn "
                  "<-visc muc> <-soft SOFSCL>\n";
        return -1;
    }
    int idata[2], ni = 2;
    if (OPS_GetIntInput(&ni, idata) < 0) {
        opserr << "WARNING contactPlane - could not read tag/slaveSurfTag\n";
        return -1;
    }
    double d[7];
    int nd = 7;
    if (OPS_GetDoubleInput(&nd, d) < 0) {
        opserr << "WARNING contactPlane - could not read nx ny nz px py pz kn\n";
        return -1;
    }
    double n[3] = {d[0], d[1], d[2]};
    double p0[3] = {d[3], d[4], d[5]};
    double kn = d[6];
    // ADR-41 D2: optional trailing `-visc <μ_c>` viscous normal-stabilization coefficient (0 ⇒ off).
    double muc = 0.0;
    // ADR-39 B1 (P4): optional `-soft <SOFSCL>` — the explicit SOFT=1 Courant-stable penalty (size
    // kn from the slave mass + dt under CentralDifferenceLadruno; inert/kn under implicit). The
    // rigid-plane m_eff = m_s (the plane is fixed). SOFSCL optional (default 0.10); ≤0 ⇒ off.
    double softScale = 0.0;
    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char *opt = OPS_GetString();
        if (opt != 0 && strcmp(opt, "-visc") == 0) {
            double v[1]; int m = 1;
            if (OPS_GetDoubleInput(&m, v) < 0) {
                opserr << "WARNING contactPlane -visc - need a coefficient value\n";
                return -1;
            }
            muc = v[0];
        } else if (opt != 0 && strcmp(opt, "-soft") == 0) {
            softScale = 0.10;                          // optional numeric SOFSCL; default 0.10
            if (OPS_GetNumRemainingInputArgs() > 0) {
                const char *p = OPS_GetString();
                bool isFlag = (p != 0 && p[0] == '-');
                OPS_ResetCurrentInputArg(-1);
                if (!isFlag) {
                    double v[1]; int m = 1;
                    if (OPS_GetDoubleInput(&m, v) < 0) {
                        opserr << "WARNING contactPlane -soft - need a SOFSCL value or omit it (default 0.1)\n";
                        return -1;
                    }
                    softScale = v[0];
                }
            }
            if (softScale > 1.0)
                opserr << "WARNING contactPlane -soft SOFSCL=" << softScale << " > 1 is UNSTABLE "
                          "(ω·dt = 2√SOFSCL > 2); use SOFSCL ≤ 1 (default 0.1).\n";
        } else {
            opserr << "WARNING contactPlane - unexpected token '" << (opt ? opt : "")
                   << "' (expected -visc, -soft, or end of arguments)\n";
            return -1;
        }
    }
    Domain *theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;
    LadrunoContactDomain *cd = theDomain->getLadrunoContactDomain();
    if (cd == 0) {
        opserr << "WARNING contactPlane - define the slave contactSurface first\n";
        return -1;
    }
    return cd->addRigidPlane(idata[0], idata[1], p0, n, kn, muc, softScale);
}

// ladrunoContactInfo -> [numContacts, numCommits, numReverts, numMortarContacts]
// (zeros if no engine). The 4th element (ADR-41 C2) is appended — existing callers that
// index [0..2] are unaffected.
int OPS_LadrunoContactInfo()
{
    Domain *theDomain = OPS_GetDomain();
    int out[4] = {0, 0, 0, 0};
    if (theDomain != 0) {
        LadrunoContactDomain *cd = theDomain->getLadrunoContactDomain();
        if (cd != 0) {
            out[0] = cd->getNumContacts();
            out[1] = cd->getNumCommits();
            out[2] = cd->getNumReverts();
            out[3] = cd->getNumMortarContacts();
        }
    }
    int four = 4;
    if (OPS_SetIntOutput(&four, out, false) < 0)
        return -1;
    return 0;
}

// ladrunoMortarPenetration -> the max KKT-active normal penetration ‖ḡ‖_∞ over all mortar
// slave nodes (max of max(0, −ḡ_I) on the active set; 0 if no engine / no contact). ADR-41
// C2.2: the convergence measure a held-load `analyzeAugmented` loop reads to stop augmenting
// (penetration → an epsN-INDEPENDENT augTol within maxAug = the headline ALM accuracy win).
int OPS_LadrunoMortarPenetration()
{
    Domain *theDomain = OPS_GetDomain();
    double pen = 0.0;
    if (theDomain != 0) {
        LadrunoContactDomain *cd = theDomain->getLadrunoContactDomain();
        if (cd != 0)
            pen = cd->getMaxMortarPenetration();
    }
    int one = 1;
    if (OPS_SetDoubleOutput(&one, &pen, true) < 0)
        return -1;
    return 0;
}

// ladrunoEdgePenetration -> the max KKT-active normal penetration ‖gN‖_∞ over all edge-edge ALM
// pairs (max of −gN where λ_N + εN·gN < 0; 0 if no engine / no -edgeAlm pair). ADR-57 E6: the
// point-like analogue of ladrunoMortarPenetration — the convergence measure a held-load
// `analyze_augmented` loop reads (pass query=ops.ladrunoEdgePenetration) to drive the edge-edge
// penetration → an εN-INDEPENDENT augTol via the one-scalar Uzawa.
int OPS_LadrunoEdgePenetration()
{
    Domain *theDomain = OPS_GetDomain();
    double pen = 0.0;
    if (theDomain != 0) {
        LadrunoContactDomain *cd = theDomain->getLadrunoContactDomain();
        if (cd != 0)
            pen = cd->getMaxEdgePenetration();
    }
    int one = 1;
    if (OPS_SetDoubleOutput(&one, &pen, true) < 0)
        return -1;
    return 0;
}

// ladrunoMortarTieResidual -> the max weighted relative-displacement bond ‖r̄‖_∞ over all mortar
// TIE slave nodes (max over nodes/components of |r̄_I,d| = |rtGlobal_I,d / aGlobal_I|; 0 if no
// engine / no tie). ADR-41 C4: the convergence measure a held-load `analyzeAugmented` loop reads to
// stop augmenting (‖r‖ → an epsTie-INDEPENDENT augTol within maxAug = the headline ALM tie win).
int OPS_LadrunoMortarTieResidual()
{
    Domain *theDomain = OPS_GetDomain();
    double res = 0.0;
    if (theDomain != 0) {
        LadrunoContactDomain *cd = theDomain->getLadrunoContactDomain();
        if (cd != 0)
            res = cd->getMaxMortarTieResidual();
    }
    int one = 1;
    if (OPS_SetDoubleOutput(&one, &res, true) < 0)
        return -1;
    return 0;
}

// ladrunoContactForce slaveNodeTag -> the total normal contact-force magnitude on an NTS slave
// node (Σ over its active master-segment pairs of tn = kn·<−gap>₊; 0 if no engine / not in
// contact). ADR-39 B3 (P2b-2c): direct nodal contact-pressure readout for the Hertz benchmark
// (the contact traction is computed by an injected FE adapter, so it is NOT in nodeReaction).
int OPS_LadrunoContactForce()
{
    if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING want - ladrunoContactForce slaveNodeTag\n";
        return -1;
    }
    int tag = 0, one = 1;
    if (OPS_GetIntInput(&one, &tag) < 0) {
        opserr << "WARNING ladrunoContactForce - could not read slaveNodeTag\n";
        return -1;
    }
    Domain *theDomain = OPS_GetDomain();
    double f = 0.0;
    if (theDomain != 0) {
        LadrunoContactDomain *cd = theDomain->getLadrunoContactDomain();
        if (cd != 0)
            f = cd->getNtsForce(tag);
    }
    if (OPS_SetDoubleOutput(&one, &f, true) < 0)
        return -1;
    return 0;
}

// ladrunoBeginAugment / ladrunoEndAugment -> open / close a held-load within-step contact
// augmentation sweep (ADR-41 D1). Between the two, Domain::commit() performs the contact Uzawa
// lambda update but fires NO recorders and bumps NO commitTag, so the held-load re-commits the
// analyze_augmented proc issues at a zero load increment produce no spurious recorder samples and
// no time advance (the capstone contract #3 "no recorder/load corruption" clause). Outside the
// pair the flag is OFF => Domain::commit() is byte-identical to stock. Idempotent; no args.
int OPS_LadrunoBeginAugment()
{
    Domain *theDomain = OPS_GetDomain();
    if (theDomain != 0)
        theDomain->setContactAugmenting(true);
    return 0;
}

int OPS_LadrunoEndAugment()
{
    Domain *theDomain = OPS_GetDomain();
    if (theDomain != 0)
        theDomain->setContactAugmenting(false);
    return 0;
}


int OPS_nodeCrd()
{
    if (OPS_GetNumRemainingInputArgs() < 1) {
	opserr << "WARNING insufficient args: nodeDisp nodeTag <dof ...>\n";
	return -1;
    }

    // tag and dof
    int data[2] = {0, -1};
    int numdata = OPS_GetNumRemainingInputArgs();
    if (numdata > 2) {
	numdata = 2;
    }

    if (OPS_GetIntInput(&numdata, data) < 0) {
	opserr<<"WARNING nodeDisp - failed to read int inputs\n";
	return -1;
    }
    data[1]--;

    // get Crds
    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;
    Node *theNode = theDomain->getNode(data[0]);
    if (theNode == 0) return -1;

    const Vector &crd = theNode->getCrds();


    // set outputs
    int size = crd.Size();
    if (data[1] >= 0) {
	if (data[1] >= size) {
	    opserr << "WARNING nodeDisp nodeTag? dof? - dofTag? too large\n";
	    return -1;
	}

	double value = crd(data[1]);
	numdata = 1;

	if (OPS_SetDoubleOutput(&numdata, &value, true) < 0) {
	    opserr<<"WARNING nodeDisp - failed to read double inputs\n";
	    return -1;
	}


    } else {

      int size = crd.Size();
      std::vector<double> values(size);
      for (int i=0; i<size; i++) {
	values[i] =  crd(i);
      }
      
      if (OPS_SetDoubleOutput(&size, &values[0], false) < 0) {
	opserr<<"WARNING nodeCrd - failed to set double inputs\n";
	return -1;
      }
    }
    
    return 0;
}


int OPS_nodeReaction()
{
    // make sure at least one other argument to contain type of system
    if (OPS_GetNumRemainingInputArgs() < 1) {
	opserr << "WARNING want - nodeReaction nodeTag? <dof?>\n";
	return -1;
   }

    int data[2] = {0, -1};
    int numdata = OPS_GetNumRemainingInputArgs();
    if (numdata > 2) numdata = 2;

    if (OPS_GetIntInput(&numdata, data) < 0) {
	opserr<<"WARNING nodeReaction - failed to read int inputs\n";
	return -1;
    }

    data[1]--;

    // get response
    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;
    const Vector *nodalResponse = theDomain->getNodeResponse(data[0], Reaction);

    if (nodalResponse == 0) {
	return -1;
    }

    int size = nodalResponse->Size();

    if (data[1] >= 0) {

      if (data[1] >= size) {
	  opserr << "WARNING nodeReaction nodeTag? dof? - dofTag? too large\n";
	  return -1;
      }

      double value = (*nodalResponse)(data[1]);
      numdata = 1;

      // now we copy the value to the tcl string that is returned
      if (OPS_SetDoubleOutput(&numdata, &value, true) < 0) {
	  opserr<<"WARNING nodeReaction - failed to set double output\n";
	  return -1;
      }

    } else {

	std::vector<double> values(size);
	for (int i=0; i<size; i++) {
	    values[i] = (*nodalResponse)(i);
	}
	if (OPS_SetDoubleOutput(&size, &values[0], false) < 0) {
	    opserr<<"WARNING nodeReaction - failed to set double output\n";
	    return -1;
	}
    }

    return 0;
}

int OPS_nodeEigenvector()
{
    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    int numdata = OPS_GetNumRemainingInputArgs();
    if (numdata < 2) {
	opserr << "WARNING want - nodeEigenVector nodeTag? eigenVector? <dof?>\n";
	return -1;
    }

    // tag, eigenvector, dof
    if (numdata > 3) numdata = 3;
    int data[3] = {0, 0, -1};

    if (OPS_GetIntInput(&numdata, data) < 0) {
	opserr<<"WARNING invalid int inputs\n";
	return -1;
    }
    int eigenvector = data[1];
    int dof = data[2];

    dof--;
    eigenvector--;

    // get eigen vectors
    Node* theNode = theDomain->getNode(data[0]);
    if (theNode == 0) {
	    opserr << "nodeEigenvector - node with tag " << data[0] << " not found\n";
	    return -1;
    }
    const Matrix &theEigenvectors = theNode->getEigenvectors();

    int size = theEigenvectors.noRows();
    int numEigen = theEigenvectors.noCols();

    if (eigenvector < 0 || eigenvector >= numEigen) {
	opserr << "WARNING nodeEigenvector nodeTag? dof? - eigenvecor too large\n";
	return -1;
    }

    if (dof >= 0) {
	if (dof >= size) {
	    opserr << "WARNING nodeEigenvector nodeTag? dof? - dofTag? too large\n";
	    return -1;
	}

	double value = theEigenvectors(dof, eigenvector);
	size = 1;

	// now we copy the value to the tcl string that is returned
	if (OPS_SetDoubleOutput(&size, &value, true) < 0) {
	    opserr<<"WARNING nodeEigenvector - failed to set double output\n";
	    return -1;
	}

    } else {

	Vector values(size);
	for (int i=0; i<size; i++) {
	    values(i) = theEigenvectors(i, eigenvector);
	}

	// now we copy the value to the tcl string that is returned
	if (OPS_SetDoubleOutput(&size, &values(0), false) < 0) {
	    opserr<<"WARNING nodeEigenvector - failed to set double output\n";
	    return -1;
	}
    }

    return 0;
}

int OPS_getTime()
{
    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;
    double time = theDomain->getCurrentTime();
    int numdata = 1;
    if (OPS_SetDoubleOutput(&numdata, &time, true) < 0) {
	opserr << "WARNING failed to get current time\n";
	return -1;
    }

    return 0;
}

int OPS_eleResponse()
{
    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    int numdata = OPS_GetNumRemainingInputArgs();
    if (numdata < 2) {
	opserr << "WARNING want - eleResponse eleTag? eleArgs...\n";
	return -1;
    }

    int tag;
    numdata = 1;
    if (OPS_GetIntInput(&numdata, &tag) < 0) {
	opserr << "could not read eleTag\n";
	return -1;
    }

    numdata = OPS_GetNumRemainingInputArgs();
    if (numdata > 0) {
	const char** argv = new const char*[numdata];
	for (int i=0; i<numdata; i++) {
	  argv[i] = new char[128];
	  // Turn everything in to a string for setResponse
	  OPS_GetStringFromAll((char*)argv[i], 128);
	}
	const Vector* data = theDomain->getElementResponse(tag, argv, numdata);

	for (int i=0; i<numdata; i++)
	  delete [] argv[i];
	delete [] argv;

	if (data != 0) {
	    int size = data->Size();
	    double* newdata = new double[size];
	    for (int i=0; i<size; i++) {
		newdata[i] = (*data)(i);
	    }
	    if (OPS_SetDoubleOutput(&size, newdata, false) < 0) {
		opserr << "WARNING failed to et response\n";
		delete [] newdata;
		return -1;
	    }
	    delete [] newdata;

	} else {
        int size = 0;
        double* newdata = 0;
        if (OPS_SetDoubleOutput(&size, newdata, false) < 0) {
            opserr << "WARNING failed to et response\n";
            return -1;
        }
	}
    } else {
        int size = 0;
        double* newdata = 0;
        if (OPS_SetDoubleOutput(&size, newdata, false) < 0) {
            opserr << "WARNING failed to et response\n";
            return -1;
        }
    }
    return 0;

}

int OPS_getLoadFactor()
{
    if (OPS_GetNumRemainingInputArgs() < 1) {
	opserr << "WARNING no load pattern supplied -- getLoadFactor\n";
	return -1;
    }

    int pattern;
    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &pattern) < 0) {
	opserr << "ERROR reading load pattern tag -- getLoadFactor\n";
	return -1;
    }

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;
    LoadPattern *thePattern = theDomain->getLoadPattern(pattern);
    if (thePattern == 0) {
	opserr << "ERROR load pattern with tag " << pattern << " not found in domain -- getLoadFactor\n";
	return -1;
    }

    double factor = thePattern->getLoadFactor();
    if (OPS_SetDoubleOutput(&numdata, &factor, true) < 0) {
	opserr << "WARNING failed to set load factor\n";
	return -1;
    }

    return 0;
}

// printNode():
// function to print out the nodal information contained in line
//     print <filename> node <flag int> <int int int>
// input: nodeArg: integer equal to arg count to node plus 1
//        output: output stream to which the results are sent
//
int printNode(OPS_Stream &output)
{
    int flag = 0; // default flag sent to a nodes Print() method
    int nodeArg = 0;
    int argc = OPS_GetNumRemainingInputArgs();

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    // if just 'print <filename> node' print all the nodes - no flag
    if (argc == 0) {
        NodeIter &theNodes = theDomain->getNodes();
        Node *theNode;
        while ((theNode = theNodes()) != 0) {
            theNode->Print(output);
        }
        return 0;
    }

    // if 'print <filename> node flag int <int int ..>' get the flag
    const char* flagArg = OPS_GetString();
    if ((strcmp(flagArg, "flag") == 0) || (strcmp(flagArg, "-flag") == 0)) {
        // get the specified flag
        if (argc < 2) {
            opserr << "WARNING print <filename> node <flag int> no int specified \n";
            return -1;
        }
        int numdata = 1;
        if (OPS_GetIntInput(&numdata, &flag) < 0) {
            opserr << "WARNING print node failed to get integer flag: \n";
            return -1;
        }
        nodeArg += 2;
    }
    else {
        OPS_ResetCurrentInputArg(2);
    }

    // now print the nodes with the specified flag, 0 by default

    // if 'print <filename> node flag'
    //     print out all the nodes in the domain with flag
    if (argc == nodeArg) {
        NodeIter &theNodes = theDomain->getNodes();
        Node *theNode;
        while ((theNode = theNodes()) != 0) {
            theNode->Print(output, flag);
        }
        return 0;
    }
    else {
        // otherwise print out the specified nodes i j k .. with flag
        int numNodes = argc - nodeArg;
        ID *theNodes = new ID(numNodes);
        for (int i = 0; i < numNodes; i++) {
            int nodeTag;
            int numdata = 1;
            if (OPS_GetIntInput(&numdata, &nodeTag) < 0) {
                opserr << "WARNING print node failed to get integer: " << endln;
                delete theNodes;
                return -1;
            }
            (*theNodes)(i) = nodeTag;
            nodeArg++;
        }
        theDomain->Print(output, theNodes, 0, flag);
        delete theNodes;
    }

    return 0;
}

int printElement(OPS_Stream &output)
{
    int flag = 0; // default flag sent to a nodes Print() method
    int eleArg = 0;
    int argc = OPS_GetNumRemainingInputArgs();

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    // if just 'print <filename> node' print all the nodes - no flag
    if (argc == 0) {
        ElementIter &theElements = theDomain->getElements();
        Element *theElement;
        while ((theElement = theElements()) != 0) {
            theElement->Print(output);
        }
        return 0;
    }

    // if 'print <filename> Element flag int <int int ..>' get the flag
    const char* eleflag = OPS_GetString();
    if ((strcmp(eleflag, "flag") == 0) || (strcmp(eleflag, "-flag")) == 0) {
        // get the specified flag
        if (argc < 2) {
            opserr << "WARNING print <filename> ele <flag int> no int specified \n";
            return -1;
        }
        int numdata = 1;
        if (OPS_GetIntInput(&numdata, &flag) < 0) {
            opserr << "WARNING print ele failed to get integer flag: \n";
            return -1;
        }
        eleArg += 2;
    }
    else {
        OPS_ResetCurrentInputArg(2);
    }

    // now print the Elements with the specified flag, 0 by default
    if (argc == eleArg) {
        ElementIter &theElements = theDomain->getElements();
        Element *theElement;
        while ((theElement = theElements()) != 0) {
            theElement->Print(output, flag);
        }
        return 0;
    }
    else {
        // otherwise print out the specified nodes i j k .. with flag
        int numEle = argc - eleArg;
        ID *theEle = new ID(numEle);
        for (int i = 0; i < numEle; i++) {
            int eleTag;
            int numdata = 1;
            if (OPS_GetIntInput(&numdata, &eleTag) < 0) {
                opserr << "WARNING print ele failed to get integer: " << endln;
                delete theEle;
                return -1;
            }
            (*theEle)(i) = eleTag;
        }
        theDomain->Print(output, 0, theEle, flag);
        delete theEle;
    }

    return 0;
}

int printAlgorithm(OPS_Stream &output)
{
    /*int eleArg = 0;
    int argc = OPS_GetNumRemainingInputArgs();

    EquiSolnAlgo** theAlgorithm = OPS_GetAlgorithm();
    if (theAlgorithm == 0) return -1;

    // if just 'print <filename> algorithm'- no flag
    if (argc == 0) {
        theAlgorithm->Print(output);
        return 0;
    }

    // if 'print <filename> Algorithm flag' get the flag
    int flag;
    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &flag) < 0) {
        opserr << "WARNING print algorithm failed to get integer flag: \n";
        return -1;
    }
    theAlgorithm->Print(output, flag);*/

    return 0;
}

int printIntegrator(OPS_Stream &output)
{
    /*int eleArg = 0;
    int argc = OPS_GetNumRemainingInputArgs();

    StaticIntegrator** theStaticIntegrator = OPS_GetStaticIntegrator();
    TransientIntegrator** theTransientIntegrator = OPS_GetTransientIntegrator();

    if (theStaticIntegrator == 0 && theTransientIntegrator == 0)
        return 0;

    IncrementalIntegrator *theIntegrator;
    if (theStaticIntegrator != 0)
        theIntegrator = theStaticIntegrator;
    else
        theIntegrator = *theTransientIntegrator;

    // if just 'print <filename> algorithm'- no flag
    if (argc == 0) {
        theIntegrator->Print(output);
        return 0;
    }

    // if 'print <filename> Algorithm flag' get the flag
    int flag;
    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &flag) < 0) {
        opserr << "WARNING print integrator failed to get integer flag: \n";
        return -1;
    }
    theIntegrator->Print(output, flag);*/

    return 0;
}

int OPS_printModelGID()
{
    // This function print's a file with node and elements in a format useful for GID
    int res = 0;
    bool hasLinear = 0;
    bool hasTri3  = 0;
    bool hasQuad4 = 0;
    bool hasQuad8 = 0;
    bool hasQuad9 = 0;
    bool hasBrick = 0;
    int startEle = 1;
    int endEle = 1;
    int eleRange = 0;
    int numdata = 1;

    FileStream outputFile;
    // OPS_Stream *output = &opserr;
    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    if (OPS_GetNumRemainingInputArgs() < 1) {
	opserr << "WARNING printGID fileName? - no filename supplied\n";
	return -1;
    }
    const char* filename = OPS_GetString();
    openMode mode = OVERWRITE;
    while (OPS_GetNumRemainingInputArgs() > 0) {
	const char* flag = OPS_GetString();
	if (strcmp(flag,"-append") == 0) {
	    mode = APPEND;
	}
	if (strcmp(flag,"-eleRange") == 0 && OPS_GetNumRemainingInputArgs()>1) {
	    //opserr<<"WARNING:commands: eleRange defined"<<endln;
	    eleRange = 1;
	    if (OPS_GetIntInput(&numdata, &startEle) < 0) {
		opserr << "WARNING print node failed to get integer: "  << endln;
		return -1;
	    }
	    if (OPS_GetIntInput(&numdata, &endEle) < 0) {
		opserr << "WARNING print node failed to get integer: " << endln;
		return -1;
	    }
	    //opserr<<"startEle = "<<startEle<<" endEle = "<<endEle<<endln;
	}
    }

    if (outputFile.setFile(filename, mode) < 0) {
        opserr << "WARNING printGID " << filename << " failed to set the file\n";
	return -1;
    }

    // Cycle over Elements to understand what type of elements are there
    ElementIter &theElements = theDomain->getElements();
    Element *theElement;
    while ((theElement = theElements()) != 0) {
	// int tag = theElement->getTag();

	// Check type of Element with Number of Nodes
	// if 2 Nodes print the Element
	int nNode = theElement->getNumExternalNodes();
	if (nNode == 2) {
	    hasLinear = 1;
	} else if (nNode == 4) {
	    hasQuad4 = 1;
	} else if (nNode == 3) {
	    hasTri3 = 1;
	} else if (nNode == 9) {
	    hasQuad9 = 1;
	} else if (nNode == 8) {
	    const char *name = theElement->getClassType();
	    if (strcmp(name,"Brick") == 0) {
		hasBrick = 1;
	    } else {
		hasQuad8 = 1;
	    }
	}
    }
    // **** Linear Elements - 2 Nodes
    if (hasLinear == 1) {
	// Print HEADER
	outputFile << "MESH \"2NMESH\" dimension 3 ElemType Linear Nnode 2" << endln;
	outputFile << "#color 0 0 255" << endln << endln;

	// Print node coordinates
	outputFile << "Coordinates" << endln;
	NodeIter &theNodes = theDomain->getNodes();
	Node *theNode;
	while ((theNode = theNodes()) != 0) {
	    int tag = theNode->getTag();
	    const Vector &crds = theNode->getCrds();
	    //outputFile << tag << "\t\t" << crds(0) << "\t" << crds(1) << "\t" << crds(2) << endln;
	    int l_tmp = crds.Size();
	    outputFile << tag << "\t\t";
	    for (int ii = 0; ii<l_tmp; ii++) {
		outputFile << crds(ii) << "\t";
	    };
	    for (int ii = l_tmp; ii<3; ii++) {
		outputFile << 0.0 << "\t";
	    };
	    outputFile << endln;
	}
	outputFile << "End coordinates" << endln << endln;

	// Print elements connectivity
	outputFile << "Elements" << endln;
	ElementIter &theElements = theDomain->getElements();
	Element *theElement;
	while ((theElement = theElements()) != 0) {
	    int tag = theElement->getTag();
	    // Check if element tag is inside theRange
	    if (((tag<=endEle) & (tag>=startEle)) || (eleRange == 0)) {
		// Check type of Element with Number of Nodes
		// if 2 Nodes print the Element
		int nNode = theElement->getNumExternalNodes();
		if (nNode == 2) {
		    Node **NodePtrs;
		    NodePtrs = theElement->getNodePtrs();
		    ID tagNodes(nNode);
		    for (int i = 0; i < nNode; i++) {
			tagNodes(i)=NodePtrs[i]->getTag();
		    }
		    outputFile << tag << "\t\t";
		    for (int i = 0; i < nNode; i++) {
			outputFile << tagNodes(i) << "\t";
		    }
		    outputFile << endln;
		}
	    }
	}
	outputFile << "End elements" << endln;
    }
    // **** Quadrilateral Elements - 4 Nodes
    if (hasQuad4 == 1) {
	// Print HEADER
	outputFile << "MESH \"4NMESH\" dimension 3 ElemType Quadrilateral Nnode 4" << endln;
	outputFile << "#color 0 255 0" << endln << endln;

	// Print node coordinates
	outputFile << "Coordinates" << endln;
	NodeIter &theNodes = theDomain->getNodes();
	Node *theNode;
	while ((theNode = theNodes()) != 0) {
	    int tag = theNode->getTag();
	    const Vector &crds = theNode->getCrds();
	    //outputFile << tag << "\t\t" << crds(0) << "\t" << crds(1) << "\t" << crds(2) << endln;
	    int l_tmp = crds.Size();
	    outputFile << tag << "\t\t";
	    for (int ii = 0; ii<l_tmp; ii++) {
		outputFile << crds(ii) << "\t";
	    };
	    for (int ii = l_tmp; ii<3; ii++) {
		outputFile << 0.0 << "\t";
	    };
	    outputFile << endln;
	}
	outputFile << "End coordinates" << endln << endln;

	// Print elements connectivity
	outputFile << "Elements" << endln;
	ElementIter &theElements = theDomain->getElements();
	Element *theElement;
	while ((theElement = theElements()) != 0) {
	    int tag = theElement->getTag();
	    // Check if element tag is inside theRange
	    if (((tag<=endEle) & (tag>=startEle)) || (eleRange == 0)) {

		// Check type of Element with Number of Nodes
		// if 2 Nodes print the Element
		int nNode = theElement->getNumExternalNodes();
		if (nNode == 4) {
		    Node **NodePtrs;
		    NodePtrs = theElement->getNodePtrs();
		    ID tagNodes(nNode);
		    for (int i = 0; i < nNode; i++) {
			tagNodes(i)=NodePtrs[i]->getTag();
		    }
		    outputFile << tag << "\t\t";
		    for (int i = 0; i < nNode; i++) {
			outputFile << tagNodes(i) << "\t";
		    }
		    outputFile << endln;
		}
	    }
	}
	outputFile << "End elements" << endln;
    }
    // **** Triangular Elements - 3 Nodes
    if (hasTri3 == 1) {
	// Print HEADER
	outputFile << "MESH \"3NMESH\" dimension 3 ElemType Triangle Nnode 3" << endln;
	outputFile << "#color 0 255 0" << endln << endln;

	// Print node coordinates
	outputFile << "Coordinates" << endln;
	NodeIter &theNodes = theDomain->getNodes();
	Node *theNode;
	while ((theNode = theNodes()) != 0) {
	    int tag = theNode->getTag();
	    const Vector &crds = theNode->getCrds();
	    //outputFile << tag << "\t\t" << crds(0) << "\t" << crds(1) << "\t" << crds(2) << endln;
	    int l_tmp = crds.Size();
	    outputFile << tag << "\t\t";
	    for (int ii = 0; ii<l_tmp; ii++) {
		outputFile << crds(ii) << "\t";
	    };
	    for (int ii = l_tmp; ii<3; ii++) {
		outputFile << 0.0 << "\t";
	    };
	    outputFile << endln;
	}
	outputFile << "End coordinates" << endln << endln;

	// Print elements connectivity
	outputFile << "Elements" << endln;
	ElementIter &theElements = theDomain->getElements();
	Element *theElement;
	while ((theElement = theElements()) != 0) {
	    int tag = theElement->getTag();
	    // Check if element tag is inside theRange
	    if (((tag<=endEle) & (tag>=startEle)) || (eleRange ==0)) {

		// Check type of Element with Number of Nodes
		// if 3 Nodes print the Element
		int nNode = theElement->getNumExternalNodes();
		if (nNode == 3) {
		    Node **NodePtrs;
		    NodePtrs = theElement->getNodePtrs();
		    ID tagNodes(nNode);
		    for (int i = 0; i < nNode; i++) {
			tagNodes(i)=NodePtrs[i]->getTag();
		    }
		    outputFile << tag << "\t\t";
		    for (int i = 0; i < nNode; i++) {
			outputFile << tagNodes(i) << "\t";
		    }
		    outputFile << endln;
		}
	    }
	}
	outputFile << "End elements" << endln;
    }
    // **** Quadrilateral Elements - 9 Nodes
    if (hasQuad9 == 1) {
	// Print HEADER
	outputFile << "MESH \"9NMESH\" dimension 3 ElemType Linear Nnode 9" << endln;
	outputFile << "#color 0 255 0" << endln << endln;

	// Print node coordinates
	outputFile << "Coordinates" << endln;
	NodeIter &theNodes = theDomain->getNodes();
	Node *theNode;
	while ((theNode = theNodes()) != 0) {
	    int tag = theNode->getTag();
	    const Vector &crds = theNode->getCrds();
	    //outputFile << tag << "\t\t" << crds(0) << "\t" << crds(1) << "\t" << crds(2) << endln;
	    int l_tmp = crds.Size();
	    outputFile << tag << "\t\t";
	    for (int ii = 0; ii<l_tmp; ii++) {
		outputFile << crds(ii) << "\t";
	    };
	    for (int ii = l_tmp; ii<3; ii++) {
		outputFile << 0.0 << "\t";
	    };
	    outputFile << endln;
	}
	outputFile << "End coordinates" << endln << endln;

	// Print elements connectivity
	outputFile << "Elements" << endln;
	ElementIter &theElements = theDomain->getElements();
	Element *theElement;
	while ((theElement = theElements()) != 0) {
	    int tag = theElement->getTag();
	    // Check if element tag is inside theRange
	    if (((tag<=endEle) & (tag>=startEle)) || (eleRange ==0)) {

		// Check type of Element with Number of Nodes
		// if 2 Nodes print the Element
		int nNode = theElement->getNumExternalNodes();
		if (nNode == 9) {
		    Node **NodePtrs;
		    NodePtrs = theElement->getNodePtrs();
		    ID tagNodes(nNode);
		    for (int i = 0; i < nNode; i++) {
			tagNodes(i)=NodePtrs[i]->getTag();
		    }
		    outputFile << tag << "\t\t";
		    for (int i = 0; i < nNode; i++) {
			outputFile << tagNodes(i) << "\t";
		    }
		    outputFile << endln;
		}
	    }
	}
	outputFile << "End elements" << endln;
    }
    // **** Hexahedra Elements - 8 Nodes
    if (hasBrick == 1) {
	// Print HEADER
	outputFile << "MESH \"8NMESH\" dimension 3 ElemType Hexahedra Nnode 8" << endln;
	outputFile << "#color 255 0 0" << endln << endln;

	// Print node coordinates
	outputFile << "Coordinates" << endln;
	NodeIter &theNodes = theDomain->getNodes();
	// MeshRegion *myRegion = theDomain->getRegion(0);
	Node *theNode;
	while ((theNode = theNodes()) != 0) {
	    int tag = theNode->getTag();
	    const Vector &crds = theNode->getCrds();
	    //outputFile << tag << "\t\t" << crds(0) << "\t" << crds(1) << "\t" << crds(2) << endln;
	    int l_tmp = crds.Size();
	    outputFile << tag << "\t\t";
	    for (int ii = 0; ii<l_tmp; ii++) {
		outputFile << crds(ii) << "\t";
	    };
	    for (int ii = l_tmp; ii<3; ii++) {
		outputFile << 0.0 << "\t";
	    };
	    outputFile << endln;
	}
	outputFile << "End coordinates" << endln << endln;

	// Print elements connectivity
	outputFile << "Elements" << endln;
	ElementIter &theElements = theDomain->getElements();
	Element *theElement;
	while ((theElement = theElements()) != 0) {
	    int tag = theElement->getTag();
	    // Check if element tag is inside theRange
	    if (((tag<=endEle) & (tag>=startEle)) || (eleRange == 0)) {

		// Check type of Element with Number of Nodes
		// if 2 Nodes print the Element
		int nNode = theElement->getNumExternalNodes();
		if (nNode == 8) {
		    Node **NodePtrs;
		    NodePtrs = theElement->getNodePtrs();
		    ID tagNodes(nNode);
		    for (int i = 0; i < nNode; i++) {
			tagNodes(i)=NodePtrs[i]->getTag();
		    }
		    outputFile << tag << "\t\t";
		    for (int i = 0; i < nNode; i++) {
			outputFile << tagNodes(i) << "\t";
		    }
		    outputFile << endln;
		}
	    }
	}
	outputFile << "End elements" << endln;
    }

    outputFile.close();
    return res;
}

int OPS_eleForce()
{
    // make sure at least one other argument to contain type of system
    if (OPS_GetNumRemainingInputArgs() < 1) {
	opserr << "WARNING want - eleForce eleTag? <dof?>\n";
	return -1;
    }

    int tag;
    int dof = -1;
    int numdata = 1;

    if (OPS_GetIntInput(&numdata, &tag) < 0) {
	opserr << "WARNING eleForce eleTag? dof? - could not read nodeTag? \n";
	return -1;
    }

    if (OPS_GetNumRemainingInputArgs() > 0) {
	if (OPS_GetIntInput(&numdata, &dof) < 0) {
	    opserr << "WARNING eleForce eleTag? dof? - could not read dof? \n";
	    return -1;
	}
    }

    dof--;

    /*
      Element *theEle = theDomain.getElement(tag);
      if (theEle == 0)
      return TCL_ERROR;

      const Vector &force = theEle->getResistingForce();
    */

    const char *myArgv[1];
    char myArgv0[8];
    strcpy(myArgv0,"forces");
    myArgv[0] = myArgv0;

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;
    const Vector *force = theDomain->getElementResponse(tag, &myArgv[0], 1);
    if (force != 0) {
	int size = force->Size();

	if (dof >= 0) {

	    if (size < dof) {
		opserr << "WARNING eleForce dof > size\n";
		return -1;
	    }

	    double value = (*force)(dof);

	    // now we copy the value to the tcl string that is returned
	    if (OPS_SetDoubleOutput(&numdata, &value, true) < 0) {
		opserr << "WARNING eleForce failed to set output\n";
		return -1;
	    }

	} else {

	    double* data = new double[size];
	    for (int i=0; i<size; i++) {
		data[i] = (*force)(i);
	    }
	    if (OPS_SetDoubleOutput(&size, data, false) < 0) {
		opserr << "WARNING eleForce failed to set outputs\n";
		delete [] data;
		return -1;
	    }

	    delete [] data;
	    return 0;
	}
    }

    return 0;
}

int OPS_eleDynamicalForce()
{
    // make sure at least one other argument to contain type of system
    if (OPS_GetNumRemainingInputArgs() < 1) {
	opserr << "WARNING want - eleForce eleTag? <dof?>\n";
	return -1;
    }

    int tag;
    int dof = -1;
    int numdata = 1;

    if (OPS_GetIntInput(&numdata, &tag) < 0) {
	opserr << "WARNING eleForce eleTag? dof? - could not read nodeTag? \n";
	return -1;
    }

    if (OPS_GetNumRemainingInputArgs() > 0) {
	if (OPS_GetIntInput(&numdata, &dof) < 0) {
	    opserr << "WARNING eleForce eleTag? dof? - could not read dof? \n";
	    return -1;
	}
    }

    dof--;

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    Element *theEle = theDomain->getElement(tag);
    if (theEle == 0) {
	opserr << "WARNING element "<<tag<<" does not exist\n";
	return -1;
    }

    const Vector &force = theEle->getResistingForceIncInertia();
    int size = force.Size();

    if (dof >= 0) {

	if (size < dof) {
	    opserr << "WARNING eleDyanmicalForce size < dof\n";
	    return -1;
	}

	double value = force(dof);

	// now we copy the value to the tcl string that is returned
	if (OPS_SetDoubleOutput(&numdata, &value, false) < 0) {
	    opserr << "WARNING eleDyanmicalForce failed to set output\n";
	    return -1;
	}

    } else {

	double* data = new double[size];
	for (int i=0; i<size; i++) {
	    data[i] = force(i);
	}
	if (OPS_SetDoubleOutput(&size, data, false) < 0) {
	    opserr << "WARNING eleDyanmicalForce failed to set outputs\n";
	    delete [] data;
	    return -1;
	}

	delete [] data;
	return 0;
    }

    return 0;

}

int OPS_nodeUnbalance()
{
    // make sure at least one other argument to contain type of system
    if (OPS_GetNumRemainingInputArgs() < 1) {
	opserr << "WARNING want - nodeUnbalance nodeTag? <dof?>\n";
	return -1;
    }

    int tag;
    int dof = -1;
    int numdata = 1;

    if (OPS_GetIntInput(&numdata, &tag) < 0) {
	opserr << "WARNING eleForce eleTag? dof? - could not read nodeTag? \n";
	return -1;
    }

    if (OPS_GetNumRemainingInputArgs() > 0) {
	if (OPS_GetIntInput(&numdata, &dof) < 0) {
	    opserr << "WARNING eleForce eleTag? dof? - could not read dof? \n";
	    return -1;
	}
    }

    dof--;

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;
    const Vector *nodalResponse = theDomain->getNodeResponse(tag, Unbalance);
    if (nodalResponse == 0) {
	opserr << "WARNING failed to get nodal response\n";
	return -1;
    }

    int size = nodalResponse->Size();

    if (dof >= 0) {

	if (size < dof) {
	    opserr << "WARNING nodeUnbalance size < dof\n";
	    return -1;
	}

	double value = (*nodalResponse)(dof);

	// now we copy the value to the tcl string that is returned
	if (OPS_SetDoubleOutput(&numdata, &value, true) < 0) {
	    opserr << "WARNING nodeUnbalance failed to set output\n";
	    return -1;
	}

    } else {

	double* data = new double[size];
	for (int i=0; i<size; i++) {
	    data[i] = (*nodalResponse)(i);
	}
	if (OPS_SetDoubleOutput(&size, data, false) < 0) {
	    opserr << "WARNING eleDyanmicalForce failed to set outputs\n";
	    delete [] data;
	    return -1;
	}

	delete [] data;
	return 0;
    }

    return 0;
}

int OPS_nodeVel()
{
    // make sure at least one other argument to contain type of system
    if (OPS_GetNumRemainingInputArgs() < 1) {
	opserr << "WARNING want - nodeVel nodeTag? <dof?>\n";
	return -1;
    }

    int tag;
    int dof = -1;
    int numdata = 1;

    if (OPS_GetIntInput(&numdata, &tag) < 0) {
	opserr << "WARNING nodeVel nodeTag? dof? - could not read nodeTag? \n";
	return -1;
    }

    if (OPS_GetNumRemainingInputArgs() > 0) {
	if (OPS_GetIntInput(&numdata, &dof) < 0) {
	    opserr << "WARNING nodeVel nodeTag? dof? - could not read dof? \n";
	    return -1;
	}
    }

    dof--;

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    const Vector *nodalResponse = theDomain->getNodeResponse(tag, Vel);
    if (nodalResponse == 0) {
	opserr << "WARNING failed to get nodal response\n";
	return -1;
    }

    int size = nodalResponse->Size();

    if (dof >= 0) {

	if (size < dof) {
	    opserr << "WARNING nodeVel size < dof\n";
	    return -1;
	}

	double value = (*nodalResponse)(dof);

	// now we copy the value to the tcl string that is returned
	if (OPS_SetDoubleOutput(&numdata, &value, true) < 0) {
	    opserr << "WARNING nodeVel failed to set output\n";
	    return -1;
	}

    } else {

	double* data = new double[size];
	for (int i=0; i<size; i++) {
	    data[i] = (*nodalResponse)(i);
	}
	if (OPS_SetDoubleOutput(&size, data, false) < 0) {
	    opserr << "WARNING nodeVel failed to set outputs\n";
	    delete [] data;
	    return -1;
	}

	delete [] data;
	return 0;
    }

    return 0;
}

int OPS_nodeAccel()
{
    // make sure at least one other argument to contain type of system
    if (OPS_GetNumRemainingInputArgs() < 1) {
	opserr << "WARNING want - nodeAccel nodeTag? <dof?>\n";
	return -1;
    }

    int tag;
    int dof = -1;
    int numdata = 1;

    if (OPS_GetIntInput(&numdata, &tag) < 0) {
	opserr << "WARNING nodeAccel nodeTag? dof? - could not read nodeTag? \n";
	return -1;
    }

    if (OPS_GetNumRemainingInputArgs() > 0) {
	if (OPS_GetIntInput(&numdata, &dof) < 0) {
	    opserr << "WARNING nodeAccel nodeTag? dof? - could not read dof? \n";
	    return -1;
	}
    }

    dof--;

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    const Vector *nodalResponse = theDomain->getNodeResponse(tag, Accel);
    if (nodalResponse == 0) {
	opserr << "WARNING failed to get nodal response\n";
	return -1;
    }

    int size = nodalResponse->Size();

    if (dof >= 0) {

	if (size < dof) {
	    opserr << "WARNING nodeAccel size < dof\n";
	    return -1;
	}

	double value = (*nodalResponse)(dof);

	// now we copy the value to the tcl string that is returned
	if (OPS_SetDoubleOutput(&numdata, &value, true) < 0) {
	    opserr << "WARNING nodeAccel failed to set output\n";
	    return -1;
	}

    } else {

	double* data = new double[size];
	for (int i=0; i<size; i++) {
	    data[i] = (*nodalResponse)(i);
	}
	if (OPS_SetDoubleOutput(&size, data, false) < 0) {
	    opserr << "WARNING nodeAccel failed to set outputs\n";
	    delete [] data;
	    return -1;
	}

	delete [] data;
	return 0;
    }

    return 0;
}

int OPS_nodeResponse()
{
    // make sure at least one other argument to contain type of system
    if (OPS_GetNumRemainingInputArgs() < 3) {
	opserr << "WARNING want - nodeResponse nodeTag? dof? responseID?\n";
	return -1;
    }

    int data[3];
    int numdata = 3;
    if (OPS_GetIntInput(&numdata, data) < 0) {
	opserr << "WARNING nodeResponse - could not read int inputs \n";
	return -1;
    }

    int tag=data[0], dof=data[1], responseID=data[2];

    dof--;

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    const Vector *nodalResponse = theDomain->getNodeResponse(tag, (NodeResponseType)responseID);
    if (nodalResponse == 0 || nodalResponse->Size() < dof || dof < 0) {
	opserr << "WARNING errors in read response\n";
	return -1;
    }

    double value = (*nodalResponse)(dof);
    numdata = 1;

    // now we copy the value to the tcl string that is returned
    if (OPS_SetDoubleOutput(&numdata, &value, true) < 0) {
	opserr << "WARNING failed to set output\n";
	return -1;
    }

    return 0;
}

int OPS_nodeCoord()
{
    // make sure at least one other argument to contain type of system
    if (OPS_GetNumRemainingInputArgs() < 1) {
	opserr << "WARNING want - nodeCoord nodeTag? <dim?>\n";
	return -1;
    }

    int tag;
    int numdata = 1;

    if (OPS_GetIntInput(&numdata, &tag) < 0) {
	opserr << "WARNING nodeCoord nodeTag? dim? - could not read nodeTag? \n";
	return -1;
    }

    int dim = -1;

    if (OPS_GetNumRemainingInputArgs() > 0) {
	const char* flag = OPS_GetString();
	if (strcmp(flag,"X") == 0 || strcmp(flag,"x") == 0) {
	    dim = 0;
	} else if (strcmp(flag,"Y") == 0 || strcmp(flag,"y") == 0) {
	    dim = 1;
	} else if (strcmp(flag,"Z") == 0 || strcmp(flag,"z") == 0) {
	    dim = 2;
	} else {
	    OPS_ResetCurrentInputArg(-1);
	    if (OPS_GetIntInput(&numdata, &dim) < 0) {
		opserr << "WARNING nodeCoord nodeTag? dim? - could not read dim? \n";
		return -1;
	    }
	    dim--;
	}
    }

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;
    Node *theNode = theDomain->getNode(tag);

    if (theNode == 0) {
	opserr << "WARNING node "<<tag<<" does not exist\n";
	return -1;
    }

    const Vector &coords = theNode->getCrds();

    int size = coords.Size();
    if (dim == -1) {

	double* data = new double[size];
	for (int i=0; i<size; i++) {
	    data[i] = coords(i);
	}
	if (OPS_SetDoubleOutput(&size, data, false) < 0) {
	    opserr << "WARNING failed to set output\n";
	    delete [] data;
	    return -1;
	}
	delete [] data;

    } else if (dim < size) {
	double value = coords(dim); // -1 for OpenSees vs C indexing
	if (OPS_SetDoubleOutput(&numdata, &value, true) < 0) {
	    opserr << "WARNING failed to set output\n";
	    return -1;
	}

    } else {
	opserr << "WARNING invalid dim\n";
	return -1;
    }

    return 0;
}

int OPS_setNodeCoord()
{
    // make sure at least one other argument to contain type of system
    if (OPS_GetNumRemainingInputArgs() < 3) {
	opserr << "WARNING want - setNodeCoord nodeTag? dim? value?\n";
	return -1;
    }

    int tag;
    int numdata = 1;

    if (OPS_GetIntInput(&numdata, &tag) < 0) {
	opserr << "WARNING setNodeCoord nodeTag? dim? value? - could not read nodeTag? \n";
	return -1;
    }

    int dim;
    double value;

    if (OPS_GetIntInput(&numdata, &dim) < 0) {
	opserr << "WARNING setNodeCoord nodeTag? dim? value? - could not read dim? \n";
	return -1;
    }
    if (OPS_GetDoubleInput(&numdata, &value) < 0) {
	opserr << "WARNING setNodeCoord nodeTag? dim? value? - could not read value? \n";
	return -1;
    }

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;
    Node *theNode = theDomain->getNode(tag);

    if (theNode == 0) {
	opserr << "WARNING node "<< tag<<" does not exist\n";
	return -1;
    }

    Vector coords(theNode->getCrds());
    coords(dim-1) = value;
    theNode->setCrds(coords);

    return 0;
}

int OPS_getPatterns()
{
  Domain* theDomain = OPS_GetDomain();
  if (theDomain == 0) return -1;
    
  LoadPattern *thePattern;
  LoadPatternIter &thePatterns = theDomain->getLoadPatterns();

  std::vector <int> data;
  
  while ((thePattern = thePatterns()) != 0)
    data.push_back(thePattern->getTag());

  int size = data.size();
  
  if (OPS_SetIntOutput(&size, data.data(), false) < 0) {
    opserr << "WARNING getPatterns - failed to set output\n";
    return -1;
  }
  
  return 0;
}

int OPS_getFixedNodes()
{
    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    SP_ConstraintIter &spIter = theDomain->getDomainAndLoadPatternSPs();
    SP_Constraint *theSP;

	std::vector <int> data;

    while ((theSP = spIter()) != 0) {
        data.push_back(theSP->getNodeTag());
    }

	// sort and make it unique
	sort( data.begin(), data.end() );
	data.erase( unique( data.begin(), data.end() ), data.end() );

	int size = data.size();

	if (OPS_SetIntOutput(&size, data.data(), false) < 0) {
	  opserr << "WARNING failed to set output\n";
	  return -1;
	}

    return 0;
}

int OPS_getFixedDOFs()
{

    if (OPS_GetNumRemainingInputArgs() < 1) {
	opserr << "WARNING want - getFixedDOFs fNodeTag?\n";
	return -1;
    }

    int tag;
    int numdata = 1;

    if (OPS_GetIntInput(&numdata, &tag) < 0) {
	  opserr << "WARNING getFixedDOFs fNodeTag? \n";
	  return -1;
    }

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    SP_ConstraintIter &spIter = theDomain->getDomainAndLoadPatternSPs();
    SP_Constraint *theSP;

	std::vector <int> data;

    while ((theSP = spIter()) != 0) {
        if (theSP->getNodeTag() == tag) {
		  data.push_back(theSP->getDOF_Number() + 1);
        }
    }

	int size = data.size();

	if (OPS_SetIntOutput(&size, data.data(), false) < 0) {
	  opserr << "WARNING failed to set output\n";
	  return -1;
	}

    return 0;
}

int OPS_getConstrainedNodes()
{

    bool all = 1;
    int rNodeTag;
    int numdata = 1;

    if (OPS_GetNumRemainingInputArgs() > 0) {
	  if (OPS_GetIntInput(&numdata, &rNodeTag) < 0) {
		opserr << "WARNING getConstrainedNodes <rNodeTag?> - could not read rNodeTag\n";
		return -1;
	  }
	all = 0;
    }

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    MP_ConstraintIter &mpIter = theDomain->getMPs();
    MP_Constraint *theMP;

	std::vector <int> data;

    while ((theMP = mpIter()) != 0) {
        if (all || rNodeTag == theMP->getNodeRetained()) {
			data.push_back(theMP->getNodeConstrained());
        }
    }

	// sort and make it unique
	sort( data.begin(), data.end() );
	data.erase( unique( data.begin(), data.end() ), data.end() );

	int size = data.size();

	if (OPS_SetIntOutput(&size, data.data(), false) < 0) {
	  opserr << "WARNING failed to set output\n";
	  return -1;
	}

    return 0;
}

int OPS_getConstrainedDOFs()
{

    if (OPS_GetNumRemainingInputArgs() < 1) {
	opserr << "WARNING want - getConstrainedDOFs cNode? <rNode?> <rDOF?>\n";
	return -1;
    }

    int cNode;
    int numdata = 1;

    if (OPS_GetIntInput(&numdata, &cNode) < 0) {
	  opserr << "WARNING getConstrainedDOFs cNode? <rNode?> <rDOF?> - could not read cNode? \n";
	  return -1;
    }

    int rNode;
	bool allNodes = 1;

    if (OPS_GetNumRemainingInputArgs() > 0) {
	  if (OPS_GetIntInput(&numdata, &rNode) < 0) {
		opserr << "WARNING getConstrainedDOFs cNode? <rNode?> <rDOF?> - could not read rNode? \n";
		return -1;
	  }
	  allNodes = 0;
	}

    int rDOF;
	bool allDOFs = 1;
    if (OPS_GetNumRemainingInputArgs() > 0) {
	  if (OPS_GetIntInput(&numdata, &rDOF) < 0) {
		opserr << "WARNING getConstrainedDOFs cNode? <rNode?> <rDOF?> - could not read rDOF? \n";
		return -1;
	  }
	  rDOF--;
	  allDOFs = 0;
	}

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    MP_ConstraintIter &mpIter = theDomain->getMPs();
    MP_Constraint *theMP;

    int tag;
    int i;
    int n;
	std::vector <int> data;

    while ((theMP = mpIter()) != 0) {
        tag = theMP->getNodeConstrained();
        if (tag == cNode) {
            if (allNodes || rNode == theMP->getNodeRetained()) {
                const ID &cDOFs = theMP->getConstrainedDOFs();
                n = cDOFs.Size();
                if (allDOFs) {
                    for (i = 0; i < n; i++) {
					  data.push_back(cDOFs(i) + 1);
                    }
                }
                else {
                    const ID &rDOFs = theMP->getRetainedDOFs();
                    for (i = 0; i < n; i++) {
                        if (rDOF == rDOFs(i))
						    data.push_back(cDOFs(i) + 1);
                    }
                }
            }
        }
    }

	int size = data.size();

	if (OPS_SetIntOutput(&size, data.data(), false) < 0) {
	  opserr << "WARNING failed to set output\n";
	  return -1;
	}

    return 0;
}


int OPS_getRetainedNodes()
{
    bool all = 1;
    int cNodeTag;
    int numdata = 1;

    if (OPS_GetNumRemainingInputArgs() > 0) {
	  if (OPS_GetIntInput(&numdata, &cNodeTag) < 0) {
		opserr << "WARNING getRetainedNodes <cNodeTag?> - could not read cNodeTag\n";
		return -1;
	  }
	all = 0;
    }

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    MP_ConstraintIter &mpIter = theDomain->getMPs();
    MP_Constraint *theMP;

	std::vector <int> data;

    while ((theMP = mpIter()) != 0) {
        if (all || cNodeTag == theMP->getNodeConstrained()) {
			data.push_back(theMP->getNodeRetained());
        }
    }

	// sort and make it unique
	sort( data.begin(), data.end() );
	data.erase( unique( data.begin(), data.end() ), data.end() );

	int size = data.size();

	if (OPS_SetIntOutput(&size, data.data(), false) < 0) {
	  opserr << "WARNING failed to set output\n";
	  return -1;
	}

    return 0;
}


int OPS_getRetainedDOFs()
{

    if (OPS_GetNumRemainingInputArgs() < 1) {
	opserr << "WARNING want - getRetainedDOFs rNode? <cNode?> <cDOF?>\n";
	return -1;
    }

    int rNode;
    int numdata = 1;

    if (OPS_GetIntInput(&numdata, &rNode) < 0) {
	  opserr << "WARNING getRetainedDOFs rNode? <cNode?> <cDOF?> - could not read rNode? \n";
	  return -1;
    }

    int cNode;
	bool allNodes = 1;

    if (OPS_GetNumRemainingInputArgs() > 0) {
	  if (OPS_GetIntInput(&numdata, &cNode) < 0) {
		opserr << "WARNING getRetainedDOFs rNode? <cNode?> <cDOF?> - could not read cNode? \n";
		return -1;
	  }
	  allNodes = 0;
	}

    int cDOF;
	bool allDOFs = 1;
    if (OPS_GetNumRemainingInputArgs() > 0) {
	  if (OPS_GetIntInput(&numdata, &cDOF) < 0) {
		opserr << "WARNING getRetainedDOFs rNode? <cNode?> <cDOF?> - could not read cDOF? \n";
		return -1;
	  }
	  cDOF--;
	  allDOFs = 0;
	}

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    MP_ConstraintIter &mpIter = theDomain->getMPs();
    MP_Constraint *theMP;

    int tag;
    int i;
    int n;
	std::vector <int> data;	
    while ((theMP = mpIter()) != 0) {
        tag = theMP->getNodeRetained();
        if (tag == rNode) {
            if (allNodes || cNode == theMP->getNodeConstrained()) {
                const ID &rDOFs = theMP->getRetainedDOFs();
                n = rDOFs.Size();
                if (allDOFs) {
                    for (i = 0; i < n; i++) {
					  data.push_back(rDOFs(i) + 1);
                    }
                }
                else {
                    const ID &cDOFs = theMP->getConstrainedDOFs();
                    for (i = 0; i < n; i++) {
                        if (cDOF == cDOFs(i))
						    data.push_back(rDOFs(i) + 1);
                    }
                }
            }
        }
    }

	int size = data.size();

	if (OPS_SetIntOutput(&size, data.data(), false) < 0) {
	  opserr << "WARNING failed to set output\n";
	  return -1;
	}

    return 0;
}


int OPS_updateElementDomain()
{
    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    // Need to "setDomain" to make the change take effect.
    ElementIter &theElements = theDomain->getElements();
    Element *theElement;
    while ((theElement = theElements()) != 0) {
	theElement->setDomain(theDomain);
    }

    return 0;
}

int OPS_getNDMM()
{

	int ndm;
    int numdata = 1;

    if (OPS_GetNumRemainingInputArgs() > 0) {

	  int tag;

	  if (OPS_GetIntInput(&numdata, &tag) < 0) {
		opserr << "WARNING getNDM nodeTag? \n";
		return -1;
	  }

	  Domain* theDomain = OPS_GetDomain();
	  if (theDomain == 0) return -1;
	  Node *theNode = theDomain->getNode(tag);

	  if (theNode == 0) {
		opserr << "WARNING node "<< tag <<" does not exist\n";
		return -1;
	  }

	  const Vector& crds = theNode->getCrds();
	  ndm = crds.Size();

    } else {

	  ndm = OPS_GetNDM();

	}

	int size = 1;

	if (OPS_SetIntOutput(&size, &ndm, false) < 0) {
	  opserr << "WARNING failed to set output\n";
	  return -1;
	}

    return 0;
}

int OPS_getNDFF()
{

	int ndf;
    int numdata = 1;

    if (OPS_GetNumRemainingInputArgs() > 0) {

	  int tag;

	  if (OPS_GetIntInput(&numdata, &tag) < 0) {
		opserr << "WARNING getNDF nodeTag? \n";
		return -1;
	  }

	  Domain* theDomain = OPS_GetDomain();
	  if (theDomain == 0) return -1;
	  Node *theNode = theDomain->getNode(tag);

	  if (theNode == 0) {
		opserr << "WARNING node "<< tag <<" does not exist\n";
		return -1;
	  }

	  ndf = theNode->getNumberDOF();

    } else {

	  ndf = OPS_GetNDF();

	}

	int size = 1;

	if (OPS_SetIntOutput(&size, &ndf, false) < 0) {
	  opserr << "WARNING failed to set output\n";
	  return -1;
	}

    return 0;
}

int OPS_classType()
{
  if (OPS_GetNumRemainingInputArgs() < 2) {
    opserr << "ERROR want - classType objectType tag?\n";
    return -1;
  }
  
  std::string type = OPS_GetString();
  int tag;
  int numdata = 1;
  
  if (OPS_GetIntInput(&numdata, &tag) < 0) {
    opserr << "ERROR classType objectType tag? - unable to read tag" << endln;
    return -1;
  }
  
  if (type == "uniaxialMaterial") {
    UniaxialMaterial *theMaterial = OPS_GetUniaxialMaterial(tag);
    if (theMaterial == 0) {
      opserr << "ERROR classType - uniaxialMaterial with tag " << tag << " not found" << endln;
      return -1;
    }
    
    std::string classType = theMaterial->getClassType();
    if (OPS_SetString(classType.c_str()) < 0) {
      opserr << "ERROR failed to set classType" << endln;
      return -1;
    }      
  }
  
  else if (type == "section") {
    SectionForceDeformation *theSection = OPS_getSectionForceDeformation(tag);
    if (theSection == 0) {
      opserr << "ERROR classType - section with tag " << tag << " not found" << endln;
      return -1;
    }
    
    std::string classType = theSection->getClassType();
    if (OPS_SetString(classType.c_str()) < 0) {
      opserr << "ERROR failed to set classType" << endln;
      return -1;
    }      
  }

  else if (type == "damping") {
    Damping *theDamping = OPS_getDamping(tag);
    if (theDamping == 0) {
      opserr << "ERROR classType - damping with tag " << tag << " not found" << endln;
      return -1;
    }
    
    std::string classType = theDamping->getClassType();
    if (OPS_SetString(classType.c_str()) < 0) {
      opserr << "ERROR failed to set classType" << endln;
      return -1;
    }      
  }
	  
  else {
    opserr << "WARNING classType - " << type.c_str() << " not yet supported" << endln;
    return 0;
  }

  return 0;
}

int OPS_eleType()
{
    if (OPS_GetNumRemainingInputArgs() < 1) {
	opserr << "WARNING want - eleType eleTag?\n";
	return -1;
    }

    int tag;
    int numdata = 1;

    if (OPS_GetIntInput(&numdata, &tag) < 0) {
	opserr << "WARNING eleType eleTag? \n";
	return -1;
    }

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;
    
    char buffer[80];
    Element *theElement = theDomain->getElement(tag);
    if (theElement == 0) {
        opserr << "WARNING eleType ele " << tag << " not found" << endln;
        return -1;
    }
    const char* type = theElement->getClassType();
    sprintf(buffer, "%s", type);    
    
    if (OPS_SetString(buffer) < 0) {
      opserr << "WARNING failed to set eleType\n";
      return -1;
    }

    return 0;
}

int OPS_eleNodes()
{
    if (OPS_GetNumRemainingInputArgs() < 1) {
	opserr << "WARNING want - eleNodes eleTag?\n";
	return -1;
    }

    int tag;
    int numdata = 1;

    if (OPS_GetIntInput(&numdata, &tag) < 0) {
	opserr << "WARNING eleNodes eleTag? \n";
	return -1;
    }

    const char *myArgv[1];
    char myArgv0[80];
    strcpy(myArgv0,"nodeTags");
    myArgv[0] = myArgv0;

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    const Vector *tags = theDomain->getElementResponse(tag, &myArgv[0], 1);
    //  Element *theElement = theDomain.getElement(tag);
    if (tags != 0) {
	int numTags = tags->Size();
	int* data = new int[numTags];
	for (int i = 0; i < numTags; i++) {
	    data[i] = (int)(*tags)(i);
	}

	if (OPS_SetIntOutput(&numTags, data, false) < 0) {
	    opserr << "WARNING failed to set outputs\n";
	    delete [] data;
	    return -1;
	}

	delete [] data;
    } else {
        int numTags = 0;
        int* data = 0;
        if (OPS_SetIntOutput(&numTags, data, false) < 0) {
            opserr << "WARNING failed to set outputs\n";
            return -1;
        }
    }

    return 0;
}

int OPS_nodeDOFs()
{
    if (OPS_GetNumRemainingInputArgs() < 1) {
	opserr << "WARNING want - nodeDOFs nodeTag?\n";
	return -1;
    }

    int tag;
    int numdata = 1;

    if (OPS_GetIntInput(&numdata, &tag) < 0) {
	opserr << "WARNING nodeDOFs nodeTag?\n";
	return -1;
    }

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    Node *theNode = theDomain->getNode(tag);
    if (theNode == 0) {
	opserr << "WARNING nodeDOFs node " << tag << " not found" << endln;
	return -1;
    }
    int numDOF = theNode->getNumberDOF();

    DOF_Group *theDOFgroup = theNode->getDOF_GroupPtr();
    if (theDOFgroup == 0) {
      opserr << "WARNING nodeDOFs DOF group null" << endln;
      return -1;
    }
    const ID &eqnNumbers = theDOFgroup->getID();
    int *data = new int[numDOF];
    for (int i = 0; i < numDOF; i++) {
      data[i] = eqnNumbers(i);
    }
    if (OPS_SetIntOutput(&numDOF, data, false) < 0) {
      opserr << "WARNING nodeDOFs failed to set outputs\n";
      delete [] data;
      return -1;
    }
    delete [] data;

    return 0;
}

int OPS_nodeMass() {
    if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING want - nodeMass nodeTag? <dof>\n";
        return -1;
    }

    int tag[2] = {0, -1};
    int numdata = OPS_GetNumRemainingInputArgs();
    if (numdata > 2) {
        numdata = 2;
    }

    if (OPS_GetIntInput(&numdata, &tag[0]) < 0) {
        opserr << "WARNING nodeMass nodeTag?\n";
        return -1;
    }
    tag[1]--;

    Domain *theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    Node *theNode = theDomain->getNode(tag[0]);
    if (theNode == 0) {
        opserr << "WARNING nodeMass node " << tag << " not found" << endln;
        return -1;
    }

    int numDOF = theNode->getNumberDOF();
    const Matrix &mass = theNode->getMass();
    if (tag[1] >= 0) {
        if (tag[1] >= numDOF) {
            opserr << "WARNING: nodeMass nodeTag? dof? - dof too large\n";
            return -1;
        }
        double value = mass(tag[1], tag[1]);
        numdata = 1;
        if (OPS_SetDoubleOutput(&numdata, &value, true) < 0) {
            opserr << "WARNING: nodeMass - failed to set mass output\n";
            return -1;
        }
    } else {
        std::vector<double> data(numDOF);
        for (int i = 0; i < numDOF; i++) {
            data[i] = mass(i, i);
        }

        if (OPS_SetDoubleOutput(&numDOF, &data[0], false) < 0) {
            opserr << "WARNING nodeMass failed to set mass\n";
            return -1;
        }
    }

    return 0;
}

int OPS_nodePressure()
{
    if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING: want - nodePressure nodeTag?\n";
        return -1;
    }

    int tag;
    int numdata = 1;

    if (OPS_GetIntInput(&numdata, &tag) < 0) {
        opserr << "WARNING: nodePressure invalid tag\n";
        return -1;
    }

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    double pressure = 0.0;
    Pressure_Constraint* thePC = theDomain->getPressure_Constraint(tag);
    if(thePC != 0) {
        pressure = thePC->getPressure();
    }
    if (OPS_SetDoubleOutput(&numdata, &pressure, true) < 0) {
	opserr << "WARNING failed to get presure\n";
	return -1;
    }

    return 0;
}

int OPS_setNodePressure()
{
    if (OPS_GetNumRemainingInputArgs() < 2) {
        opserr << "WARNING: want - setNodePressure nodeTag? Pressure?\n";
        return -1;
    }

    int tag;
    int numdata = 1;

    if (OPS_GetIntInput(&numdata, &tag) < 0) {
        opserr << "WARNING: setNodePressure invalid tag\n";
        return -1;
    }

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    double pressure = 0.0;

    if (OPS_GetDoubleInput(&numdata, &pressure) < 0) {
        opserr << "WARNING: setNodePressure invalid pressure\n";
        return -1;
    }

    Pressure_Constraint* thePC = theDomain->getPressure_Constraint(tag);
    if(thePC != 0) {
         thePC->setPressure(pressure);
    }

    return 0;
}


int OPS_nodeBounds()
{

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;
    const Vector &bounds = theDomain->getPhysicalBounds();
    int size = bounds.Size();
    double* data = new double[size];
    for (int i = 0; i < size; i++)
      data[i] = bounds(i);

    if (OPS_SetDoubleOutput(&size, data, false) < 0) {
	opserr << "WARNING failed to get node bounds\n";
	delete [] data;
	return -1;
    }

    delete [] data;

    return 0;
}

int OPS_setPrecision()
{
    if (OPS_GetNumRemainingInputArgs() < 1) {
	opserr << "WARNING setPrecision precision? - no precision value supplied\n";
	return -1;
    }
    int precision;
    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &precision) < 0) {
	opserr << "WARNING setPrecision precision? - error reading precision value supplied\n";
	return -1;
    }
    opserr.setPrecision(precision);

    return 0;
}

int OPS_getEleTags()
{
    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    std::vector<int> eletags;
    if (OPS_GetNumRemainingInputArgs() < 1) {
	// return all elements

	Element *theEle;
	ElementIter &eleIter = theDomain->getElements();

	while ((theEle = eleIter()) != 0) {
	    eletags.push_back(theEle->getTag());
	}
    } else if (OPS_GetNumRemainingInputArgs() == 2) {

	// return nodes in mesh
	const char* type = OPS_GetString();
	if (strcmp(type,"-mesh") == 0) {
	    int tag;
	    int num = 1;
	    if (OPS_GetIntInput(&num, &tag) < 0) {
		opserr << "WARNING: failed to get mesh tag\n";
		return -1;
	    }
	    Mesh* msh = OPS_getMesh(tag);
	    if (msh == 0) {
		opserr << "WARNING: mesh "<<tag<<" does not exist\n";
		return -1;
	    }
	    const ID& tags = msh->getEleTags();
	    for (int i=0; i<tags.Size(); ++i) {
		eletags.push_back(tags(i));
	    }
	}
    }

    int size = 0;
    int *data = 0;
    if (!eletags.empty()) {
        size = (int) eletags.size();
        data = &eletags[0];
    }

    if (OPS_SetIntOutput(&size, data, false) < 0) {
	opserr << "WARNING failed to set outputs\n";
	return -1;
    }

    return 0;
}

int OPS_getCrdTransfTags()
{
  // Defined in CrdTransf.cpp
  ID transfTags = OPS_getAllCrdTransfTags();

  int size = transfTags.Size();
  int *data = 0;
  if (size > 0) {
    data = &transfTags[0];
  }
  
  if (OPS_SetIntOutput(&size, data, false) < 0) {
    opserr << "WARNING failed to set outputs\n";
    return -1;
  }
  
  return 0;
}

int OPS_getNodeTags() {
    Domain *theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    std::vector<int> nodetags;
    if (OPS_GetNumRemainingInputArgs() < 1) {
        // return all nodes
        Node *theNode;
        NodeIter &nodeIter = theDomain->getNodes();

        while ((theNode = nodeIter()) != 0) {
            nodetags.push_back(theNode->getTag());
        }
    } else if (OPS_GetNumRemainingInputArgs() > 1) {
        // return nodes in mesh
        const char *type = OPS_GetString();
        if (strcmp(type, "-mesh") == 0) {
            int numtags = OPS_GetNumRemainingInputArgs();
            std::set<int> nodeset;
            for (int i = 0; i < numtags; ++i) {
                int tag;
                int num = 1;
                if (OPS_GetIntInput(&num, &tag) < 0) {
                    opserr << "WARNING: failed to get mesh tag\n";
                    return -1;
                }
                Mesh *msh = OPS_getMesh(tag);
                if (msh == 0) {
                    opserr << "WARNING: mesh " << tag
                           << " does not exist\n";
                    return -1;
                }
                const ID &tags = msh->getNodeTags();
                for (int i = 0; i < tags.Size(); ++i) {
                    nodeset.insert(tags(i));
                }
                const ID &newtags = msh->getNewNodeTags();
                for (int i = 0; i < newtags.Size(); ++i) {
                    nodeset.insert(newtags(i));
                }
            }
            nodetags.assign(nodeset.begin(), nodeset.end());
        }
    }

    int size = 0;
    int *data = 0;
    if (!nodetags.empty()) {
        size = (int)nodetags.size();
        data = &nodetags[0];
    }

    if (OPS_SetIntOutput(&size, data, false) < 0) {
        opserr << "WARNING failed to set outputs\n";
        return -1;
    }

    return 0;
}

int OPS_sectionForce()
{
    // make sure at least one other argument to contain type of system
    if (OPS_GetNumRemainingInputArgs() < 2) {
	opserr << "WARNING want - sectionForce eleTag? secNum? <dof?> \n";
	return -1;
    }

    //opserr << "sectionForce: ";
    //for (int i = 0; i < argc; i++)
    //  opserr << argv[i] << ' ' ;
    //opserr << endln;

    int numdata = 2;
    int data[3];

    if (OPS_GetIntInput(&numdata, data) < 0) {
	opserr << "WARNING sectionForce eleTag? secNum? <dof?> - could not read int input? \n";
	return -1;
    }

    int tag = data[0];
    int secNum = data[1];
    int dof = -1;
    
    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    Element *theElement = theDomain->getElement(tag);
    if (theElement == 0) {
	opserr << "WARNING sectionForce element with tag " << tag << " not found in domain \n";
	return -1;
    }

    int argcc = 3;
    char a[80] = "section";
    char b[80];
    sprintf(b, "%d", secNum);
    char c[80] = "force";
    const char *argvv[3];
    argvv[0] = a;
    argvv[1] = b;
    argvv[2] = c;

    DummyStream dummy;

    Response *theResponse = theElement->setResponse(argvv, argcc, dummy);
    if (theResponse == 0) {
	return 0;
    }

    theResponse->getResponse();
    Information &info = theResponse->getInformation();
    const Vector &theVec = *(info.theVector);

    if (OPS_GetNumRemainingInputArgs() > 0) {
      numdata = 1;
      if (OPS_GetIntInput(&numdata, &dof) < 0) {
	opserr << "WARNING sectionForce eleTag? secNum? dof? - could not read int input? \n";
	delete theResponse;
	return -1;
      }      
    }

    int ndof = theVec.Size();
    if (dof > 0 && dof <= ndof) {
    
      double value = theVec(dof-1);
      numdata = 1;
      
      if (OPS_SetDoubleOutput(&numdata, &value, true) < 0) {
	opserr << "WARNING failed to set output\n";
	delete theResponse;
	return -1;
      }
    } else {
      std::vector<double> values;
      values.reserve(ndof);
	for (int i = 0; i < ndof; i++) {
	  values.push_back(theVec(i));
	}
      
      if (OPS_SetDoubleOutput(&ndof, &values[0], false) < 0) {
	opserr << "WARNING failed to set output\n";
	delete theResponse;
	return -1;
      }      
    }

    delete theResponse;

    return 0;
}

int OPS_sectionDeformation()
{
    // make sure at least one other argument to contain type of system
    if (OPS_GetNumRemainingInputArgs() < 2) {
	opserr << "WARNING want - sectionDeformation eleTag? secNum? <dof?> \n";
	return -1;
    }

    //opserr << "sectionDeformation: ";
    //for (int i = 0; i < argc; i++)
    //  opserr << argv[i] << ' ' ;
    //opserr << endln;

    int numdata = 2;
    int data[3];

    if (OPS_GetIntInput(&numdata, data) < 0) {
	opserr << "WARNING sectionDeformation eleTag? secNum? <dof?> - could not read int input? \n";
	return -1;
    }

    int tag = data[0];
    int secNum = data[1];
    int dof = -1;

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    Element *theElement = theDomain->getElement(tag);
    if (theElement == 0) {
	opserr << "WARNING sectionDeformation element with tag " << tag << " not found in domain \n";
	return -1;
    }

    int argcc = 3;
    char a[80] = "section";
    char b[80];
    sprintf(b, "%d", secNum);
    char c[80] = "deformation";
    const char *argvv[3];
    argvv[0] = a;
    argvv[1] = b;
    argvv[2] = c;

    DummyStream dummy;

    Response *theResponse = theElement->setResponse(argvv, argcc, dummy);
    if (theResponse == 0) {
	return 0;
    }

    theResponse->getResponse();
    Information &info = theResponse->getInformation();
    const Vector &theVec = *(info.theVector);
    
    if (OPS_GetNumRemainingInputArgs() > 0) {
      numdata = 1;
      if (OPS_GetIntInput(&numdata, &dof) < 0) {
	opserr << "WARNING sectionForce eleTag? secNum? dof? - could not read int input? \n";
	delete theResponse;
	return -1;
      }      
    }

    int ndof = theVec.Size();
    if (dof > 0 && dof <= ndof) {
    
      double value = theVec(dof-1);
      numdata = 1;
      
      if (OPS_SetDoubleOutput(&numdata, &value, true) < 0) {
	opserr << "WARNING failed to set output\n";
	delete theResponse;
	return -1;
      }
    } else {
      std::vector<double> values;
      values.reserve(ndof);
	for (int i = 0; i < ndof; i++) {
	  values.push_back(theVec(i));
	}
      
      if (OPS_SetDoubleOutput(&ndof, &values[0], false) < 0) {
	opserr << "WARNING failed to set output\n";
	delete theResponse;
	return -1;
      }      
    }

    delete theResponse;

    return 0;
}

int OPS_sectionStiffness()
{
    // make sure at least one other argument to contain type of system
    if (OPS_GetNumRemainingInputArgs() < 2) {
	opserr << "WARNING want - sectionStiffness eleTag? secNum? \n";
	return -1;
    }

    //opserr << "sectionStiffness: ";
    //for (int i = 0; i < argc; i++)
    //  opserr << argv[i] << ' ' ;
    //opserr << endln;

    int numdata = 2;
    int data[2];

    if (OPS_GetIntInput(&numdata, data) < 0) {
	opserr << "WARNING sectionStiffness eleTag? secNum? - could not read int input? \n";
	return -1;
    }

    int tag = data[0];
    int secNum = data[1];

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    Element *theElement = theDomain->getElement(tag);
    if (theElement == 0) {
	opserr << "WARNING sectionStiffness element with tag " << tag << " not found in domain \n";
	return -1;
    }

    int argcc = 3;
    char a[80] = "section";
    char b[80];
    sprintf(b, "%d", secNum);
    char c[80] = "stiffness";
    const char *argvv[3];
    argvv[0] = a;
    argvv[1] = b;
    argvv[2] = c;

    DummyStream dummy;

    Response *theResponse = theElement->setResponse(argvv, argcc, dummy);
    if (theResponse == 0) {
	return 0;
    }

    theResponse->getResponse();
    Information &info = theResponse->getInformation();

    const Matrix &theMat = *(info.theMatrix);
    int nsdof = theMat.noCols();
    int size = nsdof*nsdof;
    if (size == 0) {
        if (OPS_SetDoubleOutput(&size, 0, false) < 0) {
            opserr << "WARNING failed to set output\n";
            delete theResponse;
            return -1;
        }
	delete theResponse;
	return 0;
    }

    std::vector<double> values;
    values.reserve(size);

    for (int i = 0; i < nsdof; i++) {
	for (int j = 0; j < nsdof; j++) {
	    values.push_back(theMat(i,j));
	}
    }

    if (OPS_SetDoubleOutput(&size, &values[0], false) < 0) {
	opserr << "WARNING failed to set output\n";
	delete theResponse;
	return -1;
    }

    delete theResponse;

    return 0;
}

int OPS_sectionFlexibility()
{
    // make sure at least one other argument to contain type of system
    if (OPS_GetNumRemainingInputArgs() < 2) {
	opserr << "WARNING want - sectionFlexibility eleTag? secNum? \n";
	return -1;
    }

    //opserr << "sectionFlexibility: ";
    //for (int i = 0; i < argc; i++)
    //  opserr << argv[i] << ' ' ;
    //opserr << endln;

    int numdata = 2;
    int data[2];

    if (OPS_GetIntInput(&numdata, data) < 0) {
	opserr << "WARNING sectionFlexibility eleTag? secNum? - could not read int input? \n";
	return -1;
    }

    int tag = data[0];
    int secNum = data[1];

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    Element *theElement = theDomain->getElement(tag);
    if (theElement == 0) {
	opserr << "WARNING sectionFlexibility element with tag " << tag << " not found in domain \n";
	return -1;
    }

    int argcc = 3;
    char a[80] = "section";
    char b[80];
    sprintf(b, "%d", secNum);
    char c[80] = "flexibility";
    const char *argvv[3];
    argvv[0] = a;
    argvv[1] = b;
    argvv[2] = c;

    DummyStream dummy;

    Response *theResponse = theElement->setResponse(argvv, argcc, dummy);
    if (theResponse == 0) {
	return 0;
    }

    theResponse->getResponse();
    Information &info = theResponse->getInformation();

    const Matrix &theMat = *(info.theMatrix);
    int nsdof = theMat.noCols();
    int size = nsdof*nsdof;
    if (size == 0) {
        if (OPS_SetDoubleOutput(&size, 0, false) < 0) {
            opserr << "WARNING failed to set output\n";
            delete theResponse;
            return -1;
        }
	delete theResponse;
	return 0;
    }

    std::vector<double> values;
    values.reserve(size);

    for (int i = 0; i < nsdof; i++) {
	for (int j = 0; j < nsdof; j++) {
	    values.push_back(theMat(i,j));
	}
    }

    if (OPS_SetDoubleOutput(&size, &values[0], false) < 0) {
	opserr << "WARNING failed to set output\n";
	delete theResponse;
	return -1;
    }

    delete theResponse;

    return 0;
}

int OPS_cbdiDisplacement()
{
    // make sure at least one other argument to contain type of system
    if (OPS_GetNumRemainingInputArgs() < 2) {
	opserr << "WARNING want - cbdiDisplacement eleTag? x/L? \n";
	return -1;
    }

    //opserr << "sectionWeight: ";
    //for (int i = 0; i < argc; i++)
    //  opserr << argv[i] << ' ' ;
    //opserr << endln;

    int numdata = 1;
    int data[1];
    double ddata[1];

    if (OPS_GetIntInput(&numdata, data) < 0) {
	opserr << "WARNING cbdiDisplacement eleTag? x/L? - could not read int input? \n";
	return -1;
    }
    if (OPS_GetDoubleInput(&numdata, ddata) < 0) {
	opserr << "WARNING cbdiDisplacement eleTag? x/L? - could not read double input? \n";
	return -1;
    }    

    int tag = data[0];
    double xOverL = ddata[0];

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    Element *theElement = theDomain->getElement(tag);
    if (theElement == 0) {
	opserr << "WARNING cbdiDisplacment element with tag " << tag << " not found in domain \n";
	return -1;
    }

    int argcc = 1;
    char a[80] = "cbdiDisplacements";
    const char *argvv[1];
    argvv[0] = a;

    DummyStream dummy;

    Response *theResponse = theElement->setResponse(argvv, argcc, dummy);
    if (theResponse == 0) {
	return 0;
    }
    Information &info = theResponse->getInformation();
    info.theDouble = xOverL;
    
    theResponse->getResponse();

    const Matrix &theMatrix = *(info.theMatrix);
    if (xOverL < 0.0 || xOverL > 1.0) {
	opserr << "WARNING invalid xOverL\n";
	delete theResponse;
	return -1;
    }

    double value[3];
    value[0] = theMatrix(0,0);
    value[1] = theMatrix(0,1);
    value[2] = theMatrix(0,2);
    
    numdata = 3;
    if (OPS_SetDoubleOutput(&numdata, &value[0], false) < 0) {
	opserr << "WARNING failed to set output\n";
	delete theResponse;
	return -1;
    }

    delete theResponse;

    return 0;
}

int OPS_basicDeformation()
{
    // make sure at least one other argument to contain type of system
    if (OPS_GetNumRemainingInputArgs() < 1) {
	opserr << "WARNING want - basicDeformation eleTag? \n";
	return -1;
    }

    //opserr << "basicDeformation: ";
    //for (int i = 0; i < argc; i++)
    //  opserr << argv[i] << ' ' ;
    //opserr << endln;

    int numdata = 1;
    int tag;

    if (OPS_GetIntInput(&numdata, &tag) < 0) {
	opserr << "WARNING basicDeformation eleTag? - could not read eleTag? \n";
	return -1;
    }

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    Element *theElement = theDomain->getElement(tag);
    if (theElement == 0) {
	opserr << "WARNING basicDeformation element with tag " << tag << " not found in domain \n";
	return -1;
    }

    int argcc = 1;
    char a[80] = "basicDeformation";
    const char *argvv[1];
    argvv[0] = a;

    DummyStream dummy;

    Response *theResponse = theElement->setResponse(argvv, argcc, dummy);

    // Try "basicDeformations"
    if (theResponse == 0) {
      char a[80] = "basicDeformations";
      argvv[0] = a;
      theResponse = theElement->setResponse(argvv, argcc, dummy);      
    }
    
    if (theResponse == 0) {
	return 0;
    }

    theResponse->getResponse();
    Information &info = theResponse->getInformation();

    // Vector
    if (info.theVector != 0) {
      const Vector &theVec = *(info.theVector);
      int nbf = theVec.Size();

      std::vector<double> data(nbf);
      for (int i=0; i<nbf; i++) {
	data[i] = theVec(i);
      }

      if (OPS_SetDoubleOutput(&nbf, &data[0], false) < 0) {
	opserr << "WARNING failed to set output\n";
	delete theResponse;
	return -1;
      }
    }
    // Scalar
    else {
      int nbf = 1;
      double data = info.theDouble;
      if (OPS_SetDoubleOutput(&nbf, &data, false) < 0) {
	opserr << "WARNING failed to set output\n";
	delete theResponse;
	return -1;
      }      
    }

    delete theResponse;

    return 0;
}

int OPS_basicForce()
{
    // make sure at least one other argument to contain type of system
    if (OPS_GetNumRemainingInputArgs() < 1) {
	opserr << "WARNING want - basicForce eleTag? \n";
	return -1;
    }

    //opserr << "basicForce: ";
    //for (int i = 0; i < argc; i++)
    //  opserr << argv[i] << ' ' ;
    //opserr << endln;

    int numdata = 1;
    int tag;

    if (OPS_GetIntInput(&numdata, &tag) < 0) {
	opserr << "WARNING basicForce eleTag? - could not read eleTag? \n";
	return -1;
    }

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    Element *theElement = theDomain->getElement(tag);
    if (theElement == 0) {
	opserr << "WARNING basicForce element with tag " << tag << " not found in domain \n";
	return -1;
    }

    int argcc = 1;
    char a[80] = "basicForce";
    const char *argvv[1];
    argvv[0] = a;

    DummyStream dummy;

    Response *theResponse = theElement->setResponse(argvv, argcc, dummy);

    // Try "basicForces"
    if (theResponse == 0) {
      char a[80] = "basicForces";
      argvv[0] = a;
      theResponse = theElement->setResponse(argvv, argcc, dummy);      
    }
    
    if (theResponse == 0) {
	return 0;
    }

    theResponse->getResponse();
    Information &info = theResponse->getInformation();

    // Vector
    if (info.theVector != 0) {
      const Vector &theVec = *(info.theVector);
      int nbf = theVec.Size();

      std::vector<double> data(nbf);
      for (int i=0; i<nbf; i++) {
	data[i] = theVec(i);
      }

      if (OPS_SetDoubleOutput(&nbf, &data[0], false) < 0) {
	opserr << "WARNING failed to set output\n";
	delete theResponse;
	return -1;
      }
    }
    // Scalar
    else {
      int nbf = 1;
      double data = info.theDouble;
      if (OPS_SetDoubleOutput(&nbf, &data, false) < 0) {
	opserr << "WARNING failed to set output\n";
	delete theResponse;
	return -1;
      }      
    }

    delete theResponse;

    return 0;
}

int OPS_basicStiffness()
{
    // make sure at least one other argument to contain type of system
    if (OPS_GetNumRemainingInputArgs() < 1) {
	opserr << "WARNING want - basicStiffness eleTag? \n";
	return -1;
    }

    //opserr << "basicStiffness: ";
    //for (int i = 0; i < argc; i++)
    //  opserr << argv[i] << ' ' ;
    //opserr << endln;

    int numdata = 1;
    int tag;

    if (OPS_GetIntInput(&numdata, &tag) < 0) {
	opserr << "WARNING basicStiffness eleTag? - could not read eleTag? \n";
	return -1;
    }

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    Element *theElement = theDomain->getElement(tag);
    if (theElement == 0) {
	opserr << "WARNING basicStiffness element with tag " << tag << " not found in domain \n";
	return -1;
    }

    int argcc = 1;
    char a[80] = "basicStiffness";
    const char *argvv[1];
    argvv[0] = a;

    DummyStream dummy;

    Response *theResponse = theElement->setResponse(argvv, argcc, dummy);
    if (theResponse == 0) {
	return 0;
    }

    theResponse->getResponse();
    Information &info = theResponse->getInformation();

    // Matrix
    if (info.theMatrix != 0) {
      const Matrix &theMatrix = *(info.theMatrix);
      int nbf = theMatrix.noCols();

      std::vector<double> values;
      int size = nbf*nbf;
      if (size == 0) {
        if (OPS_SetDoubleOutput(&size, 0, false) < 0) {
	  opserr << "WARNING failed to set output\n";
	  delete theResponse;
	  return -1;
        }
        return 0;
      }
      values.reserve(size);
      
      for (int i = 0; i < nbf; i++) {
	for (int j = 0; j < nbf; j++) {
	  values.push_back(theMatrix(i,j));
	}
      }

      if (OPS_SetDoubleOutput(&size, &values[0], false) < 0) {
	opserr << "WARNING failed to set output\n";
	delete theResponse;
	return -1;
      }
    }
    // Scalar
    else {
      int nbf = 1;
      double data = info.theDouble;
      if (OPS_SetDoubleOutput(&nbf, &data, false) < 0) {
	opserr << "WARNING failed to set output\n";
	delete theResponse;
	return -1;
      }      
    }

    delete theResponse;

    return 0;
}

int OPS_version()
{
    if (OPS_SetString(OPS_VERSION) < 0) {
	opserr << "WARNING failed to set version string\n";
	return -1;
    }

    return 0;
}

int OPS_logFile()
{
    if (OPS_GetNumRemainingInputArgs() < 1) { 
	opserr << "WARNING logFile fileName? - no filename supplied\n";
	return -1;
    }
    openMode mode = OVERWRITE;
    bool echo = true;

    const char* filename = OPS_GetString();
    if (strcmp(filename, "Invalid String Input!") == 0) {
	opserr << "WARNING: invalid string input\n";
	return -1;
    }

    while (OPS_GetNumRemainingInputArgs() > 0) {

	const char* opt = OPS_GetString();
	
	if (strcmp(opt,"-append") == 0) {
	    mode = APPEND;
	} else if (strcmp(opt,"-noEcho") == 0) {
	    echo = false;
	}
    }

    if (opserr.setFile(filename, mode, echo) < 0) {
	opserr << "WARNING logFile " << filename << " failed to set the file\n";
	return -1;
    }

    // const char *pwd = getInterpPWD(interp);  
    // simulationInfo.addOutputFile(argv[1], pwd);

    return 0;
}

// Sensitivity:BEGIN /////////////////////////////////////////////
int OPS_sensNodeDisp()
{
    // make sure at least one other argument to contain type of system
    if (OPS_GetNumRemainingInputArgs() < 3) {
	opserr << "WARNING want - sensNodeDisp nodeTag? dof? paramTag?\n";
	return -1;
    }    

    int data[3];
    int numdata = 3;
    if (OPS_GetIntInput(&numdata, &data[0]) < 0) {
	opserr << "WARNING: failed to get tag, dof or paramTag\n";
	return -1;
    }

    Domain* domain = OPS_GetDomain();
    if (domain == 0) return 0;


    Node *theNode = domain->getNode(data[0]);
    if (theNode == 0) {
	opserr << "sensNodeDisp: node " << data[0] << " not found" << "\n";
	return -1;
    }

    Parameter *theParam = domain->getParameter(data[2]);
    if (theParam == 0) {
	opserr << "sensNodeDisp: parameter " << data[2] << " not found" << "\n";
	return -1;
    }

    int gradIndex = theParam->getGradIndex();

    double value = theNode->getDispSensitivity(data[1],gradIndex);

    numdata = 1;
    if (OPS_SetDoubleOutput(&numdata, &value, true) < 0) {
	opserr<<"WARNING failed to set output\n";
	return -1;
    }

    return 0;
}

int OPS_sensNodeVel()
{
    // make sure at least one other argument to contain type of system
    if (OPS_GetNumRemainingInputArgs() < 3) {
	opserr << "WARNING want - sensNodeVel nodeTag? dof? paramTag?\n";
	return -1;
    }    

    int data[3];
    int numdata = 3;
    if (OPS_GetIntInput(&numdata, &data[0]) < 0) {
	opserr << "WARNING: failed to get tag, dof or paramTag\n";
	return -1;
    }

    Domain* domain = OPS_GetDomain();
    if (domain == 0) return 0;


    Node *theNode = domain->getNode(data[0]);
    if (theNode == 0) {
	opserr << "sensNodeVel: node " << data[0] << " not found" << "\n";
	return -1;
    }

    Parameter *theParam = domain->getParameter(data[2]);
    if (theParam == 0) {
	opserr << "sensNodeVel: parameter " << data[2] << " not found" << "\n";
	return -1;
    }

    int gradIndex = theParam->getGradIndex();

    double value = theNode->getVelSensitivity(data[1],gradIndex);

    numdata = 1;
    if (OPS_SetDoubleOutput(&numdata, &value, true) < 0) {
	opserr<<"WARNING failed to set output\n";
	return -1;
    }

    return 0;
}

int OPS_sensNodeAccel()
{
    // make sure at least one other argument to contain type of system
    if (OPS_GetNumRemainingInputArgs() < 3) {
	opserr << "WARNING want - sensNodeAccel nodeTag? dof? paramTag?\n";
	return -1;
    }    

    int data[3];
    int numdata = 3;
    if (OPS_GetIntInput(&numdata, &data[0]) < 0) {
	opserr << "WARNING: failed to get tag, dof or paramTag\n";
	return -1;
    }

    Domain* domain = OPS_GetDomain();
    if (domain == 0) return 0;


    Node *theNode = domain->getNode(data[0]);
    if (theNode == 0) {
	opserr << "sensNodeAccel: node " << data[0] << " not found" << "\n";
	return -1;
    }

    Parameter *theParam = domain->getParameter(data[2]);
    if (theParam == 0) {
	opserr << "sensNodeAccel: parameter " << data[2] << " not found" << "\n";
	return -1;
    }

    int gradIndex = theParam->getGradIndex();

    double value = theNode->getAccSensitivity(data[1],gradIndex);

    numdata = 1;
    if (OPS_SetDoubleOutput(&numdata, &value, true) < 0) {
	opserr<<"WARNING failed to set output\n";
	return -1;
    }

    return 0;
}

int OPS_sensLambda()
{
    if (OPS_GetNumRemainingInputArgs() < 2) {
	opserr << "WARNING no load pattern supplied -- getLoadFactor\n";
	return -1;
    }

    int data[2];
    int numdata = 2;
    if (OPS_GetIntInput(&numdata, &data[0]) < 0) {
	opserr << "WARNING: failed to read patternTag or paramTag\n";
	return -1;
    }

    Domain* domain = OPS_GetDomain();
    if (domain == 0) return 0;

    LoadPattern *thePattern = domain->getLoadPattern(data[0]);
    if (thePattern == 0) {
	opserr << "ERROR load pattern with tag " << data[0] << " not found in domain\n";
	return -1;
    }

    Parameter *theParam = domain->getParameter(data[1]);
    if (theParam == 0) {
	opserr << "sensLambda: parameter " << data[1] << " not found" << "\n";
	return -1;
    }
  
    int gradIndex = theParam->getGradIndex();
    double factor = thePattern->getLoadFactorSensitivity(gradIndex);

    numdata = 1;
    if (OPS_SetDoubleOutput(&numdata, &factor, true) < 0) {
	opserr<<"WARNING failed to set output\n";
	return -1;
    }

    return 0;
}

int OPS_sensSectionForce()
{
    if (OPS_GetNumRemainingInputArgs() < 3) {
	opserr << "WARNING want - sensSectionForce eleTag? <secNum?> dof? paramTag?\n";
	return -1;
    }    
  
    //opserr << "sensSectionForce: ";
    //for (int i = 0; i < argc; i++) 
    //  opserr << argv[i] << ' ' ;
    //opserr << endln;

    int numdata = OPS_GetNumRemainingInputArgs();
    std::vector<int> data(numdata);
    if (OPS_GetIntInput(&numdata, &data[0]) < 0) {
	opserr << "WARNING: failed to read input data\n";
	return -1;
    }

    int tag, dof, paramTag;
    int secNum = -1;
    if (numdata == 3) {
	tag = data[0];
	dof = data[1];
	paramTag = data[2]; 
    } else {
	tag = data[0];
	secNum = data[1];
	dof = data[2];
	paramTag = data[3]; 
    }

    Domain* domain = OPS_GetDomain();
    if (domain == 0) return 0;

    ParameterIter &pIter = domain->getParameters();
    Parameter *theParam;
    while ((theParam = pIter()) != 0) {
	theParam->activate(false);
    }

    theParam = domain->getParameter(paramTag);
    int gradIndex = theParam->getGradIndex();
    theParam->activate(true);

    Element *theElement = domain->getElement(tag);
    if (theElement == 0) {
	opserr << "WARNING sensSectionForce element with tag " << tag << " not found in domain \n";
	return -1;
    }

    char a[80] = "section";
    char b[80];
    sprintf(b, "%d", secNum);
    char c[80] = "dsdh";
    const char *argvv[3];
    int argcc = 3;
    argvv[0] = a;
    argvv[1] = b;
    argvv[2] = c;
    if (secNum < 0) { // For zeroLengthSection
	argcc = 2;
	argvv[1] = c;
    }

    DummyStream dummy;

    Response *theResponse = theElement->setResponse(argvv, argcc, dummy);
    if (theResponse == 0) {
	numdata = 1;
	double res = 0.0;
	if (OPS_SetDoubleOutput(&numdata, &res, true) < 0) {
	    opserr<<"WARNING failed to set output\n";
	    return -1;
	}
	return 0;
    }

    theResponse->getResponseSensitivity(gradIndex);
    Information &info = theResponse->getInformation();

    Vector theVec = *(info.theVector);

    numdata = theVec.Size();
    if (OPS_SetDoubleOutput(&numdata, &theVec(dof-1), false) < 0) {
	opserr<<"WARNING failed to set output\n";
	return -1;
    }

    theParam->activate(false);

    delete theResponse;

    return 0;
}

int OPS_sensNodePressure()
{
    // make sure at least one other argument to contain type of system
    if (OPS_GetNumRemainingInputArgs() < 2) {
	opserr << "WARNING want - sensNodePressure nodeTag? paramTag?\n";
	return -1;
    }    

    int data[2];
    int numdata = 2;
    if (OPS_GetIntInput(&numdata, &data[0]) < 0) {
	opserr << "WARNING: failed to get tag or paramTag\n";
	return -1;
    }

    Domain* domain = OPS_GetDomain();
    if (domain == 0) return 0;

    double dp = 0.0;
    Pressure_Constraint* thePC = domain->getPressure_Constraint(data[0]);
    if(thePC != 0) {
        // int ptag = thePC->getPressureNode();
        // Node* pNode = theDomain.getNode(ptag);
        Node* pNode = thePC->getPressureNode();
        if(pNode != 0) {

            Parameter *theParam = domain->getParameter(data[1]);
            if (theParam == 0) {
                opserr << "sensNodePressure: parameter " << data[1] << " not found" << endln;
                return -1;
            }

            int gradIndex = theParam->getGradIndex();
            dp = pNode->getVelSensitivity(1,gradIndex);
        }
    }

    numdata = 1;
    if (OPS_SetDoubleOutput(&numdata, &dp, true) < 0) {
	opserr<<"WARNING failed to set output\n";
	return -1;
    }

    return 0;
}

int OPS_getEleClassTags()
{
    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    int numdata = OPS_GetNumRemainingInputArgs();

	std::vector <int> data;

	// all element tags
    if (numdata < 1) {
	  Element *theEle;
	  ElementIter &theEles = theDomain->getElements();

	  while ((theEle = theEles()) != 0) {
		data.push_back(theEle->getClassTag());
	  }

	  // specific element tag
    } else if (numdata == 1) {
	  int eleTag;

	  if (OPS_GetIntInput(&numdata, &eleTag) < 0) {
		opserr << "could not read eleTag\n";
		return -1;
	  }

	  Element *theEle = theDomain->getElement(eleTag);
	  if (theEle == 0) {
		  opserr << "getEleClassTags - element with tag " << eleTag << " not found" << endln;
		  return -1;
	  }

	  data.push_back(theEle->getClassTag());

	} else {
	  opserr << "WARNING want - getEleClassTags <eleTag?>\n";
	  return -1;
    }

	int size = data.size();

	if (OPS_SetIntOutput(&size, data.data(), false) < 0) {
	  opserr << "WARNING failed to set output\n";
	  return -1;
	}

    return 0;
}

int OPS_getEleLoadClassTags()
{
    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    int numdata = OPS_GetNumRemainingInputArgs();

	std::vector <int> data;

    if (numdata < 1) {
	  LoadPattern *thePattern;
	  LoadPatternIter &thePatterns = theDomain->getLoadPatterns();

	  while ((thePattern = thePatterns()) != 0) {
		ElementalLoadIter theEleLoads = thePattern->getElementalLoads();
		ElementalLoad* theLoad;

		while ((theLoad = theEleLoads()) != 0) {
		  data.push_back(theLoad->getClassTag());
		}

	  }

	} else if (numdata == 1) {

	  int patternTag;
	  if (OPS_GetIntInput(&numdata, &patternTag) < 0) {
		opserr << "could not read patternTag\n";
		return -1;
	  }

	  LoadPattern *thePattern = theDomain->getLoadPattern(patternTag);
	  if (thePattern == nullptr) {
		opserr << "ERROR load pattern with tag " << patternTag << " not found in domain -- getEleLoadClassTags\n";
		return -1;
	  }
	  ElementalLoadIter theEleLoads = thePattern->getElementalLoads();
	  ElementalLoad* theLoad;

	  while ((theLoad = theEleLoads()) != 0) {
		data.push_back(theLoad->getClassTag());
	  }

	} else {
	opserr << "WARNING want - getEleLoadClassTags <patternTag?>\n";
	return -1;
    }


	int size = data.size();

	if (OPS_SetIntOutput(&size, data.data(), false) < 0) {
	  opserr << "WARNING failed to set output\n";
	  return -1;
	}

    return 0;
}

int OPS_getEleLoadTags()
{
    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    int numdata = OPS_GetNumRemainingInputArgs();

	std::vector <int> data;

    if (numdata < 1) {
	  LoadPattern *thePattern;
	  LoadPatternIter &thePatterns = theDomain->getLoadPatterns();

	  while ((thePattern = thePatterns()) != 0) {
		ElementalLoadIter theEleLoads = thePattern->getElementalLoads();
		ElementalLoad* theLoad;

		while ((theLoad = theEleLoads()) != 0) {
		  data.push_back(theLoad->getElementTag());
		}

	  }

	} else if (numdata == 1) {

	  int patternTag;
	  if (OPS_GetIntInput(&numdata, &patternTag) < 0) {
		opserr << "could not read patternTag\n";
		return -1;
	  }

	  LoadPattern* thePattern = theDomain->getLoadPattern(patternTag);
	  if (thePattern == nullptr) {
		opserr << "ERROR load pattern with tag " << patternTag << " not found in domain -- getEleLoadTags\n";
		return -1;
	  }
	  ElementalLoadIter& theEleLoads = thePattern->getElementalLoads();
	  ElementalLoad* theLoad;

	  while ((theLoad = theEleLoads()) != 0) {
		data.push_back(theLoad->getElementTag());
	  }

	} else {
	opserr << "WARNING want - getEleLoadTags <patternTag?>\n";
	return -1;
    }

	int size = data.size();

	if (OPS_SetIntOutput(&size, data.data(), false) < 0) {
	  opserr << "WARNING failed to set output\n";
	  return -1;
	}

    return 0;
}

int OPS_getEleLoadData()
{
    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    int numdata = OPS_GetNumRemainingInputArgs();

	std::vector <double> data;

    if (numdata < 1) {
	  LoadPattern *thePattern;
	  LoadPatternIter &thePatterns = theDomain->getLoadPatterns();

	  int typeEL;

	  while ((thePattern = thePatterns()) != 0) {
		ElementalLoadIter &theEleLoads = thePattern->getElementalLoads();
		ElementalLoad* theLoad;

		while ((theLoad = theEleLoads()) != 0) {
		  const Vector &eleLoadData = theLoad->getData(typeEL, 1.0);

		  int eleLoadDataSize = eleLoadData.Size();
		  for (int i = 0; i < eleLoadDataSize; i++) {
			data.push_back(eleLoadData(i));
		  }
		}
	  }

	} else if (numdata == 1) {

	  int patternTag;
	  if (OPS_GetIntInput(&numdata, &patternTag) < 0) {
		opserr << "could not read patternTag\n";
		return -1;
	  }

	  LoadPattern* thePattern = theDomain->getLoadPattern(patternTag);
	  if (thePattern == nullptr) {
		opserr << "ERROR load pattern with tag " << patternTag << " not found in domain -- getEleLoadData\n";
		return -1;
	  }
	  ElementalLoadIter& theEleLoads = thePattern->getElementalLoads();
	  ElementalLoad* theLoad;

	  int typeEL;

	  while ((theLoad = theEleLoads()) != 0) {
		const Vector &eleLoadData = theLoad->getData(typeEL, 1.0);

		int eleLoadDataSize = eleLoadData.Size();
		for (int i = 0; i < eleLoadDataSize; i++) {
		  data.push_back(eleLoadData(i));
		}
	  }

	} else {
	opserr << "WARNING want - getEleLoadData <patternTag?>\n";
	return -1;
    }

	int size = data.size();

	if (OPS_SetDoubleOutput(&size, data.data(), false) < 0) {
	  opserr << "WARNING failed to set output\n";
	  return -1;
	}

    return 0;
}

int OPS_getNodeLoadTags()
{
    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    int numdata = OPS_GetNumRemainingInputArgs();

	std::vector <int> data;

    if (numdata < 1) {
	  LoadPattern *thePattern;
	  LoadPatternIter &thePatterns = theDomain->getLoadPatterns();

	  while ((thePattern = thePatterns()) != 0) {
		NodalLoadIter theNodLoads = thePattern->getNodalLoads();
		NodalLoad* theNodLoad;

		while ((theNodLoad = theNodLoads()) != 0) {
		  data.push_back(theNodLoad->getNodeTag());
		}

	  }

	} else if (numdata == 1) {

	  int patternTag;
	  if (OPS_GetIntInput(&numdata, &patternTag) < 0) {
		opserr << "could not read patternTag\n";
		return -1;
	  }

	  LoadPattern* thePattern = theDomain->getLoadPattern(patternTag);
	  if (thePattern == nullptr) {
		opserr << "ERROR load pattern with tag " << patternTag << " not found in domain -- getEleLoadTags\n";
		return -1;
	  }
	  NodalLoadIter& theNodLoads = thePattern->getNodalLoads();
	  NodalLoad* theNodLoad;

	  while ((theNodLoad = theNodLoads()) != 0) {
		data.push_back(theNodLoad->getNodeTag());
	  }

	} else {
	opserr << "WARNING want - getNodeLoadTags <patternTag?>\n";
	return -1;
    }

	int size = data.size();

	if (OPS_SetIntOutput(&size, data.data(), false) < 0) {
	  opserr << "WARNING failed to set output\n";
	  return -1;
	}

    return 0;
}

int OPS_getNodeLoadData()
{
    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    int numdata = OPS_GetNumRemainingInputArgs();

	std::vector <double> data;

    if (numdata < 1) {
	  LoadPattern *thePattern;
	  LoadPatternIter &thePatterns = theDomain->getLoadPatterns();

	  int typeEL;

	  while ((thePattern = thePatterns()) != 0) {
		NodalLoadIter &theNodLoads = thePattern->getNodalLoads();
		NodalLoad* theNodLoad;

		while ((theNodLoad = theNodLoads()) != 0) {
		  const Vector &nodeLoadData = theNodLoad->getData(typeEL);

		  int nodeLoadDataSize = nodeLoadData.Size();
		  for (int i = 0; i < nodeLoadDataSize; i++) {
			data.push_back(nodeLoadData(i));
		  }
		}
	  }

	} else if (numdata == 1) {

	  int patternTag;
	  if (OPS_GetIntInput(&numdata, &patternTag) < 0) {
		opserr << "could not read patternTag\n";
		return -1;
	  }

	  LoadPattern* thePattern = theDomain->getLoadPattern(patternTag);
	  if (thePattern == nullptr) {
		opserr << "ERROR load pattern with tag " << patternTag << " not found in domain -- getNodeLoadData\n";
		return -1;
	  }
	  NodalLoadIter& theNodLoads = thePattern->getNodalLoads();
	  NodalLoad* theNodLoad;

	  int typeEL;

	  while ((theNodLoad = theNodLoads()) != 0) {
		const Vector &nodeLoadData = theNodLoad->getData(typeEL);

		int nodeLoadDataSize = nodeLoadData.Size();
		for (int i = 0; i < nodeLoadDataSize; i++) {
		  data.push_back(nodeLoadData(i));
		}
	  }

	} else {
	opserr << "WARNING want - getNodeLoadData <patternTag?>\n";
	return -1;
    }

	int size = data.size();

	if (OPS_SetDoubleOutput(&size, data.data(), false) < 0) {
	  opserr << "WARNING failed to set output\n";
	  return -1;
	}

    return 0;
}


int OPS_getNumElements()
{
    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

	int nEles = theDomain->getNumElements();
	int size = 1;

	if (OPS_SetIntOutput(&size, &nEles, false) < 0) {
	  opserr << "WARNING failed to set output\n";
	  return -1;
	}

    return 0;
}

// Sensitivity:END /////////////////////////////////////////////
