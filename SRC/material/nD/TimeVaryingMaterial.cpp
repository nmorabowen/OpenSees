/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

// José L. Larenas & José A. Abell (UANDES)
// Massimo Petracca - ASDEA Software, Italy
//
// A Wapper material that allow's time varying stiffness and
// resistance properties for the wrapped material.
//

#include <TimeVaryingMaterial.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <OPS_Globals.h>
#include <elementAPI.h>
#include <Parameter.h>
#include <Domain.h>


void *OPS_TimeVaryingMaterial(void)
{
    const char *usage =
        "nDMaterial TimeVarying $tag $projMatTag $N t1 .. tN E1 .. EN K1 .. KN A1 .. AN\n";

    // check arguments
    if (OPS_GetNumRemainingInputArgs() < 3) {
        opserr << "nDMaterial TimeVarying Error: too few arguments.\n" << usage;
        return nullptr;
    }

    // get integer data
    int iData[2];
    int numData = 2;
    if (OPS_GetInt(&numData, iData) != 0)  {
        opserr << "nDMaterial TimeVarying Error: invalid nDMaterial tags.\n";
        return nullptr;
    }

    int Ndatapoints = 0;
    numData = 1;
    if (OPS_GetInt(&numData, &Ndatapoints) != 0 || Ndatapoints < 1)  {
        opserr << "nDMaterial TimeVarying Error: invalid number of data points N\n";
        return nullptr;
    }
    if (OPS_GetNumRemainingInputArgs() < 4 * Ndatapoints) {
        opserr << "nDMaterial TimeVarying Error: expected " << 4 * Ndatapoints
               << " values (t, E, K, A histories) for tag " << iData[0] << ".\n" << usage;
        return nullptr;
    }

    Vector t(Ndatapoints);
    Vector E(Ndatapoints);
    Vector K(Ndatapoints);
    Vector A(Ndatapoints);
    Vector *histories[4] = {&t, &E, &K, &A};
    const char *names[4] = {"t", "E", "K", "A"};
    for (int h = 0; h < 4; ++h) {
        numData = Ndatapoints;
        if (OPS_GetDouble(&numData, &(*histories[h])(0)) != 0) {
            opserr << "nDMaterial TimeVarying Error: reading " << names[h]
                   << " data for nDMaterial TimeVarying with tag " << iData[0] << ".\n";
            return nullptr;
        }
    }

    // validate histories
    for (int k = 0; k < Ndatapoints; ++k) {
        if (k > 0 && t(k) <= t(k - 1)) {
            opserr << "nDMaterial TimeVarying Error: times must be strictly increasing (tag " << iData[0] << ").\n";
            return nullptr;
        }
        if (A(k) <= 0.0 || K(k) <= 0.0 || E(k) <= 0.0 || E(k) >= 9.0 * K(k)) {
            opserr << "nDMaterial TimeVarying Error: need A > 0, K > 0 and 0 < E < 9K at every point (tag "
                   << iData[0] << ", point " << k + 1 << ").\n";
            return nullptr;
        }
    }

    // get the projected material to map
    NDMaterial *theProjMaterial = OPS_getNDMaterial(iData[1]);
    if (theProjMaterial == 0) {
        opserr << "WARNING: nDMaterial does not exist.\n";
        opserr << "nDMaterial: " << iData[1] << "\n";
        opserr << "nDMaterial TimeVarying: " << iData[0] << "\n";
        return nullptr;
    }

    {
        NDMaterial *test = theProjMaterial->getCopy("ThreeDimensional");
        if (test == 0) {
            opserr << "nDMaterial TimeVarying Error: material " << iData[1]
                   << " has no ThreeDimensional version.\n";
            return nullptr;
        }
        delete test;
    }

    // create the TimeVarying wrapper
    NDMaterial* theTimeVaryingMaterial = new TimeVaryingMaterial(
        iData[0],
        *theProjMaterial,
        t, E, K , A);
    if (theTimeVaryingMaterial == 0) {
        opserr << "nDMaterial TimeVarying Error: failed to allocate a new material.\n";
        return nullptr;
    }

    // done
    return theTimeVaryingMaterial;
}

std::map<int, Vector> TimeVaryingMaterial::time_histories;
std::map<int, Vector> TimeVaryingMaterial::E_histories;
std::map<int, Vector> TimeVaryingMaterial::K_histories;
std::map<int, Vector> TimeVaryingMaterial::A_histories;
std::map<int, bool> TimeVaryingMaterial::new_time_step ;
std::map<int, double> TimeVaryingMaterial::E;
std::map<int, double> TimeVaryingMaterial::G;
std::map<int, double> TimeVaryingMaterial::A;
std::map<int, double> TimeVaryingMaterial::nu;

TimeVaryingMaterial::TimeVaryingMaterial()
    : NDMaterial(0, ND_TAG_TimeVaryingMaterial)
{

}

TimeVaryingMaterial::TimeVaryingMaterial(
    int tag,
    NDMaterial &theProjMat,
    const Vector& t_, const Vector& E_, const Vector& K_, const Vector& A_)
    : NDMaterial(tag, ND_TAG_TimeVaryingMaterial)
{
    // copy the isotropic material
    theProjectedMaterial = theProjMat.getCopy("ThreeDimensional");
    if (theProjectedMaterial == 0) {
        opserr << "nDMaterial TimeVarying Error: failed to get a (3D) copy of the isotropic material\n";
    }

    int Ndatapoints = E_.Size();
    time_histories[tag] = Vector(Ndatapoints);
    E_histories[tag] = Vector(Ndatapoints);
    K_histories[tag] = Vector(Ndatapoints);
    A_histories[tag] = Vector(Ndatapoints);
    time_histories[tag] = t_;
    E_histories[tag] = E_;
    K_histories[tag] = K_;
    A_histories[tag] = A_;

    new_time_step[tag] = true;
}

TimeVaryingMaterial::~TimeVaryingMaterial()
{
    if (theProjectedMaterial)
    {
        delete theProjectedMaterial;
        theProjectedMaterial = 0;
    }
}

double TimeVaryingMaterial::getRho(void)
{
    return theProjectedMaterial->getRho();
}

int TimeVaryingMaterial::setTrialStrain(const Vector & strain)
{
    if (theProjectedMaterial == 0)
        return -1;

    //Compute the real strain increment
    static Vector depsilon_real(6);
    epsilon_real = strain - epsilon_internal ; // epsilon_internal comes from thermal analysis
    depsilon_real = epsilon_real - epsilon_real_n ; // depsilon_real = epsilon_real_new - epsilon_real_old

    // Get the current parameters at current time
    double current_time = OPS_GetDomain()->getCurrentTime();  //should not be used in displacement control user could set the current time as a parameter (check ASDConcrete3D)
    getParameters(current_time);

    int tag = this->getTag();
    double Ex  = E[tag];  double Ey  = E[tag];  double Ez  = E[tag];
    double Gxy = G[tag];  double Gyz = G[tag];  double Gzx = G[tag];
    double vxy = nu[tag]; double vyz = nu[tag]; double vzx = nu[tag];

    // compute the initial orthotropic constitutive tensor
    static Matrix C0(6, 6);
    C0.Zero();
    double vyx = vxy * Ey / Ex;
    double vzy = vyz * Ez / Ey;
    double vxz = vzx * Ex / Ez;
    double d = (1.0 - vxy * vyx - vyz * vzy - vzx * vxz - 2.0 * vxy * vyz * vzx) / (Ex * Ey * Ez);
    C0(0, 0) = (1.0 - vyz * vzy) / (Ey * Ez * d);
    C0(1, 1) = (1.0 - vzx * vxz) / (Ez * Ex * d);
    C0(2, 2) = (1.0 - vxy * vyx) / (Ex * Ey * d);
    C0(1, 0) = (vxy + vxz * vzy) / (Ez * Ex * d);
    C0(0, 1) = C0(1, 0);
    C0(2, 0) = (vxz + vxy * vyz) / (Ex * Ey * d);
    C0(0, 2) = C0(2, 0);
    C0(2, 1) = (vyz + vxz * vyx) / (Ex * Ey * d);
    C0(1, 2) = C0(2, 1);
    C0(3, 3) = Gxy;
    C0(4, 4) = Gyz;
    C0(5, 5) = Gzx;

    // compute the Asigma and its inverse
    if (A[tag] <= 0 ) {
        opserr << "nDMaterial TimeVarying Error: A must be greater than 0 for tag = " << tag << "\n";
        opserr << "A = " << A[tag] << endln;
        return -1;
    }

    // compute the initial projected constitutive tensor and its inverse
    static Matrix C0proj(6, 6);
    static Matrix C0proj_inv(6, 6);
    C0proj = theProjectedMaterial->getInitialTangent();
    int res = C0proj.Invert(C0proj_inv);
    if (res < 0) {
        opserr << "nDMaterial TimeVarying Error: the isotropic material gave a singular initial tangent.\n";
        return -1;
    }

    // compute the strain tensor map inv(C0_iso) * Asigma * C0_ortho
    Aepsilon.addMatrixProduct(0.0, C0proj_inv, C0, A[tag]);

    //Compute the projected strain increment
    static Vector depsilon_proj(6);
    depsilon_proj.addMatrixVector(0.0, Aepsilon, depsilon_real, 1.0); // depsilon_proj = Aepsilon(t) * depsilon_real

    //Compute the projected total strain
    static Vector epsilon_proj(6);
    epsilon_proj = epsilon_proj_n;
    epsilon_proj.addVector(1.0, depsilon_proj, 1.0); // epsilon_proj = epsilon_proj_old + depsilon_proj

    // call projected material with the projected total strain
    res = theProjectedMaterial->setTrialStrain(epsilon_proj);
    if (res != 0) {
        opserr << "TimeVaryingMaterial::setTrialStrain\n";
        opserr << "nDMaterial Projected Material Error: the isotropic material failed in setTrialStrain.\n";
        return res;
    }

    //cool! we're done!
    return 0;
}

const Vector &TimeVaryingMaterial::getStrain(void)
{
    // the strain to report back is the real strain plus the internal,
    // the internal gets removed before converting it to "real"
    static Vector epsilon_return(6);
    epsilon_return = epsilon_real + epsilon_internal;
    return epsilon_return;
}

const Vector &TimeVaryingMaterial::getStress(void)
{
    // stress in projected space
    const Vector& sigma_proj = theProjectedMaterial->getStress();

    // compute the projected stress increment
    static Vector dsigma_proj(6);
    dsigma_proj = sigma_proj - sigma_proj_n; // dsigma_proj = sigma_proj_new - sigma_proj_old

    // rescale the stress increment to back to real space
    static Vector dsigma_real(6);
    // dsigma_real = Asigma^-1 * dsigma_proj
    dsigma_real = dsigma_proj / A[this->getTag()];  // dsigma_real = Asigma^-1 * dsigma_proj

    //add real stress increment to the previous real stress
    sigma_real = sigma_real_n + dsigma_real;

    return sigma_real;
}

const Matrix &TimeVaryingMaterial::getTangent(void)
{
    // tensor in isotropic space
    const Matrix &C_proj = theProjectedMaterial->getTangent();

    // compute orthotripic tangent
    static Matrix C_real(6, 6);
    C_real.addMatrixProduct(0.0, C_proj, Aepsilon, 1 / A[this->getTag()]);

    return C_real;
}

const Matrix &TimeVaryingMaterial::getInitialTangent(void)
{
    // elasticity tensor in projected space
    const Matrix& C_proj = theProjectedMaterial->getInitialTangent();

    // compute the real tangent from the projected one
    static Matrix C_real(6, 6);
    C_real.addMatrixProduct(0.0, C_proj, Aepsilon, 1 / A[this->getTag()]);
    return C_real;
}

int TimeVaryingMaterial::commitState(void)
{
    new_time_step[this->getTag()] = true;

    const Vector& sigma_proj = theProjectedMaterial->getStress();
    const Vector& epsilon_proj = theProjectedMaterial->getStrain();

    sigma_real_n = sigma_real;
    sigma_proj_n = sigma_proj;
    epsilon_real_n = epsilon_real;
    epsilon_proj_n = epsilon_proj;
    epsilon_new_n = epsilon_new;

    return theProjectedMaterial->commitState();
}

int TimeVaryingMaterial::revertToLastCommit(void)
{
    sigma_real = sigma_real_n ;
    epsilon_real = epsilon_real_n ;
    epsilon_new = epsilon_new_n ;
    return theProjectedMaterial->revertToLastCommit();
}

int TimeVaryingMaterial::revertToStart(void)
{
    sigma_real_n.Zero();
    sigma_proj_n.Zero();
    epsilon_real_n.Zero();
    epsilon_proj_n.Zero();
    epsilon_new_n .Zero();
    return theProjectedMaterial->revertToStart();
}

NDMaterial * TimeVaryingMaterial::getCopy(void)
{
    TimeVaryingMaterial *theCopy = new TimeVaryingMaterial();
    theCopy->setTag(getTag());
    theCopy->theProjectedMaterial = theProjectedMaterial->getCopy("ThreeDimensional");
    theCopy->Aepsilon = Aepsilon;
    theCopy->epsilon_internal = epsilon_internal;
    theCopy->sigma_real = sigma_real;
    theCopy->epsilon_real = epsilon_real;
    theCopy->sigma_real_n = sigma_real_n;
    theCopy->sigma_proj_n = sigma_proj_n;
    theCopy->epsilon_real_n = epsilon_real_n;
    theCopy->epsilon_proj_n = epsilon_proj_n;
    theCopy->epsilon_new = epsilon_new;
    theCopy->epsilon_new_n = epsilon_new_n;
    return theCopy;
}

NDMaterial* TimeVaryingMaterial::getCopy(const char* code)
{
    if (strcmp(code, "ThreeDimensional") == 0)
        return getCopy();
    return NDMaterial::getCopy(code);
}

const char* TimeVaryingMaterial::getType(void) const
{
    return "ThreeDimensional";
}

int TimeVaryingMaterial::getOrder(void) const
{
    return 6;
}

void TimeVaryingMaterial::Print(OPS_Stream &s, int flag)
{
    s << "Time Varying Material, tag: " << this->getTag() << "\n";
}

int TimeVaryingMaterial::sendSelf(int commitTag, Channel &theChannel)
{

    return 0;
}

int TimeVaryingMaterial::recvSelf(int commitTag, Channel & theChannel, FEM_ObjectBroker & theBroker)
{

    return 0;
}

int TimeVaryingMaterial::setParameter(const char** argv, int argc, Parameter& param)
{
    // 4000 - init strain
    if (strcmp(argv[0], "initNormalStrain") == 0) {
        double initNormalStrain = epsilon_internal(0);
        param.setValue(initNormalStrain);
        return param.addObject(4001, this);
    }

    // forward to the adapted (isotropic) material
    return theProjectedMaterial->setParameter(argv, argc, param);
}


int TimeVaryingMaterial::updateParameter(int parameterID, Information& info)
{
    switch (parameterID) {

    case 4001:
    {
        double initNormalStrain = info.theDouble; // alpha * (Temp(t)  - Temp())
        epsilon_internal.Zero();
        epsilon_internal(0) = initNormalStrain;
        epsilon_internal(1) = initNormalStrain;
        epsilon_internal(2) = initNormalStrain;
        return 0;
    }

    // default
    default:
        return -1;
    }
}


Response* TimeVaryingMaterial::setResponse(const char** argv, int argc, OPS_Stream& s)
{
    if (argc > 0) {
        if (strcmp(argv[0], "stress") == 0 ||
                strcmp(argv[0], "stresses") == 0 ||
                strcmp(argv[0], "strain") == 0 ||
                strcmp(argv[0], "strains") == 0 ||
                strcmp(argv[0], "Tangent") == 0 ||
                strcmp(argv[0], "tangent") == 0) {
            // stresses, strain and tangent should be those of this adapter (orthotropic)
            return NDMaterial::setResponse(argv, argc, s);
        }

        else {
            // any other response should be obtained from the adapted (isotropic) material
            return theProjectedMaterial->setResponse(argv, argc, s);
        }
    }
    return NDMaterial::setResponse(argv, argc, s);
}

void TimeVaryingMaterial::getParameters(double time)
{
    int tag = this->getTag();
    if (new_time_step[tag]) {
        double K = 0;
        new_time_step[tag] = false;

        // Find the interval in which 'time' falls within time_history
        int index = 0;
        for (; index < time_histories[tag].Size(); ++index)
        {
            if (time_histories[tag](index) >= time)
            {
                break;
            }
        }

        if (index == 0)
        {
            E[tag] = E_histories[tag](0) ;
            A[tag] = A_histories[tag](0) ;
            K = K_histories[tag](0) ;

            G[tag]  = (3.0 * K * E[tag]) / (9.0 * K - E[tag]);
            nu[tag] = (3.0 * K - E[tag]) / (6.0 * K);
        }

        else if (index > 0 && index < time_histories[tag].Size())
        {
            // Perform linear interpolation for the values of E, A, and  K
            double t1    = time_histories[tag](index - 1);
            double t2    = time_histories[tag](index);
            double alpha = (time - t1) / (t2 - t1);

            E[tag] = (1.0 - alpha) * E_histories[tag](index - 1) + alpha * E_histories[tag](index);
            A[tag] = (1.0 - alpha) * A_histories[tag](index - 1) + alpha * A_histories[tag](index);
            K = (1.0 - alpha) * K_histories[tag](index - 1) + alpha * K_histories[tag](index);

            G[tag]  = (3.0 * K * E[tag]) / (9.0 * K - E[tag]);
            nu[tag] = (3.0 * K - E[tag]) / (6.0 * K);
        }

        else {
            E[tag] = E_histories[tag](time_histories[tag].Size() - 1) ;
            A[tag] = A_histories[tag](time_histories[tag].Size() - 1) ;
            K = K_histories[tag](time_histories[tag].Size() - 1) ;

            G[tag]  = (3.0 * K * E[tag]) / (9.0 * K - E[tag]);
            nu[tag] = (3.0 * K - E[tag]) / (6.0 * K);
        }
    }
}
