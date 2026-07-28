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

// Description: This file contains the class definition for Python Module
// PythonModule implements a DL_Interpreter for the python language
//

#include "PythonModule.h"
#include "PythonStream.h"
#include <OPS_Globals.h>
#include <cstring>
#include <cctype>

// Ladruno Patch 9: this translation unit is compiled once per Python target.
// The sequential target imports as `opensees`; the parallel shim translation
// units (PythonMPIModule.cpp / PythonSPModule.cpp) pre-#define these before
// #include "PythonModule.cpp" so the same source emits `openseesmp` /
// `openseessp` with distinct PyInit_* symbols. Defaults keep the sequential
// build byte-identical.
#ifndef OPS_PY_MODULE_NAME
#define OPS_PY_MODULE_NAME "opensees"
#endif
#ifndef OPS_PY_INIT_FUNC
#define OPS_PY_INIT_FUNC PyInit_opensees
#endif

// define opserr
static PythonStream sserr;
OPS_Stream *opserrPtr = &sserr;


PythonModule::PythonModule()
        : wrapper(), cmds(this) {
    // does nothing
}

PythonModule::~PythonModule() {
    // does nothing
}

int
PythonModule::run() {
    return 0;
}

int
PythonModule::addCommand(const char *, Command &) {
    return -1;
}

int
PythonModule::removeCommand(const char *) {
    return -1;
}

const char *PythonModule::trimSpaces(PyObject *o) {
  Py_ssize_t size = 0;
  const char *s = PyUnicode_AsUTF8AndSize(o, &size);

  // empty string
  if (size == 0) {
    return s;
  }

  // find first char that is not space
  Py_ssize_t firstLoc = 0;
  for (Py_ssize_t i = 0; i < size; ++i) {
    if (isspace(s[i])) {
      firstLoc = i + 1;
    } else {
      break;
    }
  }

  // find last char that is not space
  Py_ssize_t lastLoc = size - 1;
  for (Py_ssize_t i = size - 1; i >= 0; --i) {
    if (isspace(s[i])) {
      lastLoc = i - 1;
    } else {
      break;
    }
  }

  // new Py Object
  PyObject *newo = 0;

  // all spaces: make an empty string
  if (firstLoc == size || lastLoc < 0) {
    newo = PyUnicode_FromString("");
  }

  // get substring
  else if (firstLoc > 0 || lastLoc < size - 1) {
    newo = PyUnicode_Substring(o, firstLoc, lastLoc + 1);
  }

  // get string
  if (newo == 0) {
    return s;
  }

  const char* news = PyUnicode_AsUTF8(newo);
  Py_DECREF(newo);
  return news;
}

int
PythonModule::getNumRemainingInputArgs(void) {
    return wrapper.getNumberArgs() - wrapper.getCurrentArg();
}

int
PythonModule::getInt(int *data, int numArgs) {
    if ((wrapper.getNumberArgs() - wrapper.getCurrentArg()) < numArgs) {
        return -1;
    }

    for (int i = 0; i < numArgs; i++) {
        PyObject *o = PyTuple_GetItem(wrapper.getCurrentArgv(), wrapper.getCurrentArg());
        wrapper.incrCurrentArg();
        if (PyLong_Check(o) || PyFloat_Check(o) || PyBool_Check(o)) {
            PyErr_Clear();
            data[i] = PyLong_AsLong(o);
            if (PyErr_Occurred()) {
                return -1;
            }
        } else {
            return -1;
        }
    }

    return 0;
}

int
PythonModule::getDouble(double *data, int numArgs) {
    if ((wrapper.getNumberArgs() - wrapper.getCurrentArg()) < numArgs) {
        return -1;
    }

    for (int i = 0; i < numArgs; i++) {
        PyObject *o = PyTuple_GetItem(wrapper.getCurrentArgv(), wrapper.getCurrentArg());
        wrapper.incrCurrentArg();
        if (PyLong_Check(o) || PyFloat_Check(o) || PyBool_Check(o)) {
            PyErr_Clear();
            data[i] = PyFloat_AsDouble(o);
            if (PyErr_Occurred()) {
                return -1;
            }
        } else {
            return -1;
        }
    }

    return 0;
}

int PythonModule::getDoubleList(int* size, Vector* data)
{
    if (wrapper.getCurrentArg() >= wrapper.getNumberArgs()) {
        return -1;
    }

    PyObject* o = PyTuple_GetItem(wrapper.getCurrentArgv(), wrapper.getCurrentArg());
    wrapper.incrCurrentArg();

    if (PyList_Check(o)) {
        *size = PyList_Size(o);
        data->resize(*size);
        for (int i = 0; i < *size; i++) {
            PyErr_Clear();
            PyObject* item = PyList_GetItem(o, i);
            if (!(PyLong_Check(item) || PyFloat_Check(item) || PyBool_Check(item))) {
                opserr << "PythonModule::getDoubleList error: item " << i << " in list is not a float (or int or bool)\n";
                return -1;
            }
            (*data)(i) = PyFloat_AsDouble(item);
            if (PyErr_Occurred()) {
                return -1;
            }
        }
    }
    else if (PyTuple_Check(o)) {
        *size = PyTuple_Size(o);
        data->resize(*size);
        for (int i = 0; i < *size; i++) {
            PyErr_Clear();
            PyObject* item = PyTuple_GetItem(o, i);
            if (!(PyLong_Check(item) || PyFloat_Check(item) || PyBool_Check(item))) {
                opserr << "PythonModule::getDoubleList error: item " << i << " in tuple is not a float (or int or bool)\n";
                return -1;
            }
            (*data)(i) = PyFloat_AsDouble(item);
            if (PyErr_Occurred()) {
                return -1;
            }
        }
    }
    else {
      // Removing this error message for Path series list inputs -- MHS
      //opserr << "PythonModule::getDoubleList error: input is neither a list nor a tuple\n";
        return -1;
    }

    return 0;
}

const char *
PythonModule::getString() {
    if (wrapper.getCurrentArg() >= wrapper.getNumberArgs()) {
        return 0;
    }

    PyObject *o = PyTuple_GetItem(wrapper.getCurrentArgv(), wrapper.getCurrentArg());
    wrapper.incrCurrentArg();
#if PY_MAJOR_VERSION >= 3
    if (!PyUnicode_Check(o)) {
        return 0;
    }

    // PyObject* space = PyUnicode_FromString(" ");
    // PyObject* empty = PyUnicode_FromString("");
    // PyObject* newo = PyUnicode_Replace(o, space, empty, -1);
    const char* res = trimSpaces(o);
    // Py_DECREF(newo);
    // Py_DECREF(space);
    // Py_DECREF(empty);

    return res;
#else
    if (!PyString_Check(o)) {
        return 0;
    }

    return PyString_AS_STRING(o);
#endif
}

const char *PythonModule::getStringFromAll(char* buffer, int len) {
    if (wrapper.getCurrentArg() >= wrapper.getNumberArgs()) {
        return 0;
    }

    PyObject *o =
        PyTuple_GetItem(wrapper.getCurrentArgv(), wrapper.getCurrentArg());
    wrapper.incrCurrentArg();

    // check if int
    if (PyLong_Check(o) || PyBool_Check(o)) {
        PyErr_Clear();
        int data = PyLong_AsLong(o);
        if (PyErr_Occurred()) {
            return 0;
        }
        snprintf(buffer, len, "%d", data);
        return buffer;
    }
    // check if double
    else if (PyFloat_Check(o)) {
        PyErr_Clear();
        double data = PyFloat_AsDouble(o);
        if (PyErr_Occurred()) {
            return 0;
        }
        snprintf(buffer, len, "%.20f", data);
        return buffer;
    }

#if PY_MAJOR_VERSION >= 3
    if (!PyUnicode_Check(o)) {
        return 0;
    }

    // PyObject *space = PyUnicode_FromString(" ");
    // PyObject *empty = PyUnicode_FromString("");
    // PyObject *newo = PyUnicode_Replace(o, space, empty, -1);
    const char* res = trimSpaces(o);
    // Py_DECREF(newo);
    // Py_DECREF(space);
    // Py_DECREF(empty);

    int lenres = int(strlen(res)) + 1;
    if (lenres > len) {
        lenres = len;
    }

    strncpy(buffer, res, lenres);

    return buffer;
#else
    if (!PyString_Check(o)) {
        return 0;
    }

    return PyString_AS_STRING(o);
#endif
}

void*
PythonModule::getVoidPtr()
{
    if (wrapper.getCurrentArg() >= wrapper.getNumberArgs()) {
        return nullptr;
    }
    PyObject *obj =
        PyTuple_GetItem(wrapper.getCurrentArgv(), wrapper.getCurrentArg());
    wrapper.incrCurrentArg();
    return static_cast<void*>(obj);
}

int
PythonModule::getStringCopy(char **stringPtr) {
    return -1;
}

int 
PythonModule::evalDoubleStringExpression(const char* theExpression, double& current_val)
{
    if (theExpression == 0) {
        opserr << "OPS_EvalDoubleStringExpression Error: Expression not set\n";
        return -1;
    }

    // run the string and get results
    PyObject* py_main = PyImport_AddModule("__main__");
    if (py_main == NULL) {
        opserr << "OPS_EvalDoubleStringExpression Error: cannot add module  __main__\n";
        return -1;
    }
    PyObject* py_dict = PyModule_GetDict(py_main);
    if (py_main == NULL) {
        opserr << "OPS_EvalDoubleStringExpression Error: cannot get dict of module __main__\n";
        return -1;
    }
    PyObject* PyRes = PyRun_String(theExpression, Py_eval_input, py_dict, py_dict);

    if (PyRes == NULL) {
        opserr << "OPS_EvalDoubleStringExpression Error: failed to evaluate expression\n";
        return -1;
    }

    // get results
    if (!(PyLong_Check(PyRes) || PyFloat_Check(PyRes) || PyBool_Check(PyRes))) {
        opserr << "OPS_EvalDoubleStringExpression Error: the expression must return a float (or int or bool)\n";
        return -1;
    }
    current_val = PyFloat_AsDouble(PyRes);

    // done
    return 0;
}

void
PythonModule::resetInput(int cArg) {
    wrapper.resetCommandLine(cArg);
}

int
PythonModule::setInt(int *data, int numArgs, bool scalar) {
    wrapper.setOutputs(data, numArgs, scalar);

    return 0;
}

int PythonModule::setInt(std::vector<std::vector<int>> &data) {
    wrapper.setOutputs(data);
    return 0;
}

int PythonModule::setInt(std::map<const char*, int>& data) {
    wrapper.setOutputs(data);
    return 0;
}

int PythonModule::setInt(std::map<const char*, std::vector<int>>& data) {
    wrapper.setOutputs(data);
    return 0;
}

int
PythonModule::setDouble(double *data, int numArgs, bool scalar) {
    wrapper.setOutputs(data, numArgs, scalar);

    return 0;
}

int PythonModule::setDouble(std::vector<std::vector<double>> &data) {
    wrapper.setOutputs(data);
    return 0;
}

int PythonModule::setDouble(std::map<const char*, double>& data) {
    wrapper.setOutputs(data);
    return 0;
}

int PythonModule::setDouble(std::map<const char*, std::vector<double>>& data) {
    wrapper.setOutputs(data);
    return 0;
}

int
PythonModule::setString(const char *str) {
    wrapper.setOutputs(str);

    return 0;
}

int
PythonModule::setString(std::vector<const char*>& data) {
    wrapper.setOutputs(data);
    return 0;
}

int
PythonModule::setString(std::vector<std::vector<const char*>>& data) {
    wrapper.setOutputs(data);
    return 0;
}

int PythonModule::setString(std::map<const char*, const char*>& data) {
    wrapper.setOutputs(data);
    return 0;
}

int PythonModule::setString(std::map<const char*, std::vector<const char*>>& data) {
    wrapper.setOutputs(data);
    return 0;
}

int PythonModule::setGenericDict(GenericDict& data) {
    wrapper.setGenericDictOutput(data);
    return 0;
}

int
PythonModule::runCommand(const char *cmd) {
    return PyRun_SimpleString(cmd);
}

static PythonModule *module = 0;

PyMethodDef *getmethodsFunc() {
    module = new PythonModule;
    PythonWrapper *wrapper = module->getWrapper();
    wrapper->addOpenSeesCommands();

    return wrapper->getMethods();
}

void cleanupFunc() {
    module->getCmds().wipe();
    if (module != 0) {
        delete module;
    }
}

struct module_state {
    PyObject *error;
};

#if PY_MAJOR_VERSION >= 3
#define GETSTATE(m) ((struct module_state*)PyModule_GetState(m))
#else
#define GETSTATE(m) (&_state)
static struct module_state _state;
#endif

// static PyObject *
// error_out(PyObject *m)
// {
//     struct module_state *st = GETSTATE(m);
//     PyErr_SetString(st->error, "something bad happened");
//
//     return NULL;
// }

#if PY_MAJOR_VERSION >= 3

static int opensees_traverse(PyObject *m, visitproc visit, void *arg) {
    Py_VISIT(GETSTATE(m)->error);

    return 0;
}

static int opensees_clear(PyObject *m) {
    Py_CLEAR(GETSTATE(m)->error);

    return 0;
}

static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        OPS_PY_MODULE_NAME,
        NULL,
        sizeof(struct module_state),
        getmethodsFunc(),
        NULL,
        opensees_traverse,
        opensees_clear,
        NULL
};

#define INITERROR return NULL

PyMODINIT_FUNC
OPS_PY_INIT_FUNC(void)

#else
#define INITERROR return

//void
PyMODINIT_FUNC
initopensees(void)
#endif
{
#if PY_MAJOR_VERSION >= 3
    PyObject *pymodule = PyModule_Create(&moduledef);
#else
    PyObject *pymodule = Py_InitModule(OPS_PY_MODULE_NAME, getmethodsFunc());
#endif

    if (pymodule == NULL)
        INITERROR;

    // Ladruno banner. Suppress with LADRUNO_OPENSEES_QUIET=1.
    if (getenv("LADRUNO_OPENSEES_QUIET") == nullptr) {
        static const char *ladrunoBanner =
            "\n"
            " ▄█          ▄████████ ████████▄     ▄████████ ███    █▄  ███▄▄▄▄    ▄██████▄\n"
            "███         ███    ███ ███   ▀███   ███    ███ ███    ███ ███▀▀▀██▄ ███    ███\n"
            "███         ███    ███ ███    ███   ███    ███ ███    ███ ███   ███ ███    ███\n"
            "███         ███    ███ ███    ███  ▄███▄▄▄▄██▀ ███    ███ ███   ███ ███    ███\n"
            "███       ▀███████████ ███    ███ ▀▀███▀▀▀▀▀   ███    ███ ███   ███ ███    ███\n"
            "███         ███    ███ ███    ███ ▀███████████ ███    ███ ███   ███ ███    ███\n"
            "███▌    ▄   ███    ███ ███   ▄███   ███    ███ ███    ███ ███   ███ ███    ███\n"
            "█████▄▄██   ███    █▀  ████████▀    ███    ███ ████████▀   ▀█   █▀   ▀██████▀\n"
            "▀                                   ███    ███\n"
            "\n";
        PySys_FormatStdout("%s", ladrunoBanner);
#ifdef OPENSEES_VERSION
        PySys_FormatStdout("Ladruno OpenSees build: %s\n", OPENSEES_VERSION);
#endif
        // Ladruno feature list. Auto-generated between FEATURES-START/END by
        // Ladruno_scripts/patch_banner.py from Ladruno_scripts/banner_features.txt.
        // Edit that .txt then re-run the script — do not hand-edit between markers.
        // FEATURES-START
        static const char *kFeatures =
"      Ladruno fork — active features:\n"
"        • OpenSeesPyMP — import openseesmp (MPI-parallel Python)\n"
"        • EnergyBalance recorder (per-region sidecar; -v2 channels: E_inject, E_lnvd, E_modal, E_hg, KE_ms; -ownedNodes MPI gate + per-rank files)\n"
"        • ExplicitBathe — Noh–Bathe explicit (-lnvd FLAC damping, -sms [-consistent] mass scaling)\n"
"        • CentralDifferenceLadruno — robust explicit central difference\n"
"        • CentralDifferenceSMS — selective mass scaling (DT2MS-style, explicit)\n"
"        • CentralDifferenceSMSConsistent — consistent (Olovsson) mass scaling (explicit)\n"
"        • HRZ mass-conserving lumping (-lump hrz) for dt_cr\n"
"        • Queryable critical time step (dt_cr)\n"
"        • BezierTri6 — quadratic Bézier triangle (Kadapa 2018)\n"
"        • BezierTet10 — quadratic Bézier tetrahedron (Kadapa 2018)\n"
"        • LadrunoBrick — unified hex (std/bbar/uri/ssp/eas + hourglass)\n"
"        • LadrunoBrick20 — 20-node serendipity quadratic hex (std 27-pt + uri 2x2x2 Barlow + HRZ -lumped, reduce-to 20NodeBrick)\n"
"        • LadrunoQuad/LadrunoCST — 2D continuum (Quad std/bbar/ssp/eas; CST std)\n"
"        • LadrunoQuad/LadrunoCST/LadrunoLST/LadrunoCSTPair — 2D continuum (Quad std/bbar/ssp/eas; CST + LST-T6 std; -geom finite via shared 2D F-kernel; CSTPair = 2-tri F-bar-Patch macro dSNPO 15.1.9)\n"
"        • Solid geometry methods — -geom linear/corot/finite (+ F-bar)\n"
"        • LadrunoJ2 — combined iso + Chaboche AF kinematic J2\n"
"        • LadrunoJ2Finite — finite-strain J2, co-rotating backstress (+IMPL-EX)\n"
"        • LadrunoUniaxialJ2 — uniaxial Chaboche AF J2 fiber/truss\n"
"        • StagedDefGrad/InitDefGrad — finite-strain stress-free staged activation (F_rel=F·F0⁻¹)\n"
"        • StagedStrain — small-strain (2D+3D) stress-free staged activation (ε_rel=ε−ε0)\n"
"        • Lemaitre ductile damage (-damage lemaitre, +lch, +IMPL-EX)\n"
"        • LadrunoIMKBeam — concentrated-plasticity IMK beam hinges\n"
"        • LadrunoRebarBuckling — rebar buckling wrapper (Dhakal–Maekawa / Gomes–Appleton, cyclic re-straightening)\n"
"        • LadrunoRCConcrete — RC plastic-damage + MCFT compression softening (shell)\n"
"        • LadrunoRCFiniteStrain — finite-strain (Hencky) RC plastic-damage view (-geom finite)\n"
"        • LadrunoConcrete3D — CDPM2-grade solid concrete (Menétrey-Willam + dual damage; 3D/plane/axisym/shell)\n"
"        • LadrunoBondSlip — 1D bond-slip τ–s (CEB-FIP MC2010)\n"
"        • LadrunoCohesiveHinge — discrete cohesive moment–rotation hinge (Gf per hinge, ∫M dθ=Gf)\n"
"        • LadrunoCohesiveHingeBiaxial — coupled Mz–My cohesive interaction-surface hinge (elliptical onset, per-axis Gf)\n"
"        • LadrunoEmbeddedRebar — rebar in non-matching mesh (bond/penalty/AL)\n"
"        • LadrunoEmbeddedNode — stress-free node-to-host embed tie (penalty/AL, implicit+explicit)\n"
"        • LadrunoDistributingCoupling — RBE3 distributing coupling (ndf-mismatch moment transfer)\n"
"        • LadrunoKinematicCoupling — RBE2 rigid kinematic coupling (selectable -dof, transport)\n"
"        • LadrunoArcLength — adaptive arc-length radius (Ramm desired-iter)\n"
"        • LadrunoDynamicRelaxation — quasi-static DR (Gershgorin mass + kinetic damping)\n"
"        • LadrunoIndirectControl — indirect/CMOD control (snap-back tracking)\n"
"        • LadrunoHHT — sensitivity/DDM HHT (reliability/fragility/FORM on damped α)\n"
"        • LadrunoGeneralizedAlpha — sensitivity/DDM Chung-Hulbert generalized-α (DDM)\n"
"        • LadrunoStabilizedUnbalance — true-equilibrium test for -stabilize\n"
"        • LadrunoProjection — momentum-conserving explicit constraint projection (equalDOF/rigidLink/diaphragm)\n"
"        • LadrunoTie — kinematic non-conforming mesh-tie generator (collocation + integral-mortar, -dual sparse biorthogonal basis; solid + ndf-6 shell rotational DOFs, -hermite rotation-consistent w–θ edge transfer, -shellSolid plane-section shell-edge↔solid-face tie; emits EQ-constraints for LadrunoProjection; exact, dt_cr-neutral)\n"
"        • LadrunoContact — NTS + mortar/ALM contact (penalty + commit-cycle Uzawa, Coulomb/Tresca friction, mesh-tying, viscous stabilization, consistent dn/du normal tangent, SOFT=1 + SOFT=2 Courant-stable explicit penalty, finite-sliding re-emit, averaged nodal-normal smoothing -smoothNormal w/ convex-ridge facet ownership)\n"
"        • stressesPlaneStrain — exact out-of-plane σ_zz element response (plane strain)\n"
"        • Ladruno — modular HDF5 .ladruno recorder\n"
"        • LadrunoMonitor — live SWMR-HDF5 analysis monitor (recorder Monitor; -every/-hz gated streaming)\n"
"        • LadrunoRigidBody — 6-DOF rigid body (condensed mass + side-channel SO(3); finite-rotation slaving + moment gather; rocking-foundation MVP)\n"
"        • LadrunoSolidShell — 8-node ANS/EAS solid-shell (through-thickness σ33 punching/bearing host; Dvorkin-Bathe + Betsch-Stein + state-EAS)\n"
"        • LadrunoComplexEigen — complexEigen: complex/state-space modal for non-classical damping (true per-mode zeta + omega_d + phased mode shapes; assembled getDamp() projection, QZ reduced pencil)\n"
"        • FeastEigen — band-targeted eigensolver: eigen -feast fmin fmax (ALL modes in the band, MKL contour integration; -certify = Sturm/inertia completeness; -rci = dfeast_srci RCI, distributed dmumps inner solve under openseesmp)\n"
"        • LadrunoModalResponse — modalResponseHistory: exact PWL modal-superposition transient (per-mode Nigam-Jennings recurrence, base-accel or -load nodal-force pattern w/ -series; -damp/-rayleigh/-modalDamp; commits per step so recorders capture the history)\n"
"        • modal combination — responseSpectrumAnalysis -combine SRSS|CQC|ABS|TenPercent (native CQC/SRSS/ABS/Ten-Percent of the per-mode response-spectrum fields; -damp/-modalDamp for CQC zeta)\n"
"        • modal frequency domain — frequencyResponse / steadyStateDynamics / randomResponse: modal FRF H_a(Om) sweep, steady-state harmonic amplitude + stationary random response PSD->RMS (one-sided Hz input PSD; RMS/nu0/Davenport peak -stats; base-accel or -load nodal-force pattern; lin/log/resonance-biased grid; disp/vel/accel)\n"
"        • LadrunoUP — unified Biot u-p saturated-porous element (33017): honest pore-pressure DOF (p IS p — statics well-posed: drained + steady seepage), T3/Q4/H8 equal-order + McGann auto-stab, Bezier T6/Tet10 Taylor-Hood vertex-p (heterogeneous ndf), -formulation bbar, -dynSeepage; -geom corot 3D large-rotation u-p (core-frame total force through the SolidTransformation seam, H/S reference-frame, R^T seepage drive); needs a general solver (unsymmetric tangent)\n"
"        • LadrunoPorousOverlay — persistent-fluid staggered u-p overlay (pattern LadrunoPorousOverlay 33022): pore-pressure field OUTSIDE the DOF graph over any T3/Q4/H8 solid region (zero element changes, symmetric solvers legal again); fixed-stress RATE-form fs1 (default L=classic alpha^2/K_dr, -fsL oedometric opt-in), fluid survives remove element (-onRemoval keep|drain), -subcycle auto, -staticMode hold|steady; LadrunoStaggeredAnalyze iterated fixed-stress driver (-tol/-maxIter/-pScale/-verbose/-stats) + parameter-route moduli/stage transport (parameter ... loadPattern E|nu|layerE|layerNu); explicit lane -fsL zero (L=0 under CD at dt<=0.5x the undrained pencil) + overlay-aware criticalTimeStep (undrained Qe*Se^-1*Qe^T pencil); p-field recorder channels (recorder ladruno -overlay / Monitor -overlay + MODEL/OVERLAYS topology rows); -fluidUpdate explicit fully matrix-free lane (lumped-S* forward fluid step at load-application time, dual-CFL advisory with diffusion-slack report) + SMS undrained-pencil composability (lumped CentralDifferenceSMS certified; consistent Olovsson priced-but-warned)\n"
"        • system Pardiso — threaded MKL PARDISO desktop sparse-direct solver (symbolic/numeric/solve split for factorization reuse; -matrixType 1|2 symmetric upper-triangle half-storage = 1.35x faster + 42% less peak memory, exact; -krylov L factorization-preconditioned CGS reuses the L/U across a CHANGED tangent = 1.51x under full Newton at 51k DOF, 1.57x stacked with -matrixType 1, SPD/unsym only; -stats peak-memory counters; Tcl + Python)\n"
"        • system Mumps — distributed MUMPS cluster sparse-direct solver (Ladruno additions: -BLR <eps> Block Low-Rank approximate factorization ICNTL(35)+CNTL(7), APPROXIMATE so keep it OFF on oracle paths; -stats dumps INFOG/RINFOG incl. INFOG(21) peak MB/proc; -ICNTL35/-CNTL7 expert overrides; -matrixType 0|1|2 symmetric half-storage — all of these Tcl + Python since ADR-75 P2h, which fixed a Tcl ladder that DROPPED them silently; -commSplit sub-communicator groups is Python/openseesmp ONLY)\n"
"\n";
        // FEATURES-END
        PySys_FormatStdout("%s\n", kFeatures);
        PySys_FormatStdout("Suppress this banner with LADRUNO_OPENSEES_QUIET=1\n\n");
    }

    struct module_state *st = GETSTATE(pymodule);

    // add OpenSeesError
    st->error = PyErr_NewExceptionWithDoc(OPS_PY_MODULE_NAME ".OpenSeesError", "Internal OpenSees errors.", NULL, NULL);
    if (st->error == NULL) {
        Py_DECREF(pymodule);
        INITERROR;
    }
    Py_INCREF(st->error);
    PyModule_AddObject(pymodule, "OpenSeesError", st->error);

    // add OpenSeesParameter dict
    auto *par = PyDict_New();
    if (par == NULL) {
      INITERROR;
    }
    if (PyModule_AddObject(pymodule, "OpenSeesParameter", par) < 0) {
        Py_DECREF(par);
        INITERROR;
    }

    sserr.setError(st->error);

    Py_AtExit(cleanupFunc);

#if PY_MAJOR_VERSION >= 3
    return pymodule;
#endif
}
