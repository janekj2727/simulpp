/*
 * An implementation of Verlet integration method (Verlet + TRVP + Nose–Hoover thermostat + barostat)
 * for simul++ simulation package
 * Author JJ, Date Feb 2022
 */

/*
 * In atoms R following values are stored:
 * R(0, :): x(i) y(i) z(i)             common for all integrators
 * R(1, :): h*vx(i) h*vy(i) h*vz(i)    common for all integrators (h* because of Gear)
 * R(2, :): x(i) - x(i-1), ...         this integrator specific (h*v(i-1/2)) (!!!) (not to be rescaled by lambda(i+1))
 * R(3, :): x(i+1) y(i+1) z(i+1)       Verlet specific
 *
 * R(0, :) and R(1, :) must be defined before integrator initialization
 *
 * (In case of TRVP – this IS the case) the next values of R are:
 * R(4, :): h*vx(i)^P, ...             predicted velocity
 * R(5, :): x(i-1) - x(i-2), ...       difference of positions
 * R(6, :): x(i-2) - x(i-3), ...       dtto
 * ...
 */

/*
 * MACSIMUS versions of VERLET velocity:
 * #define VERLET 0 The simplest, fastest, and least accurate O(h) version, v(t) = [r(t + h) − r(t)]/h.
 *   Essentially, velocity is shifted by h/2 from positions. Correct averaged energy of
 *   the harmonic oscillator. Good enough with friction (Nose) thermostat. (Added in
 *   V2.4a)
 * #define VERLET 1 Compromised speed, velocity accurate to O(h^2): v(t) = [r(t + h) − r(t − h)]/(2h).
 *   Averaged energy of harmonic oscillator has error O(h^2). (Added in V2.4a) This
 *   is equivalent to the velocity Verlet algorithm and to k3m2e ala Gear method.
 * #define VERLET 2 Best energy conservation (exact for harmonic oscillator but with O(h^2)
 *   error), slowest, v(t)^2 = [r(t)−r(t−h)]/h·[r(t+h)−r(t)]/h. (The off-diagonal components
 *   of the pressure tensor are approximated as average of both possible h/2-shifted terms).
 * #define VERLET 3 Compromised speed, velocity and energy conservation accurate to O(h^2)
 *   (slightly worse than for VERLET=1), but with correct averaged energy of harmonic
 *   oscillator. Energy and pressure tensor components are averages of both shifted values
 *   v(t) = [r(t + h) − r(t)]/h and v(t) = [r(t) − r(t − h)]/h. This is the recommended default
 *   (in MACSIMUS).
 */

/*
 * Nose–Hoover thermostat + barostat implemented according to our (Janek+Kolafa, 2022) article and MACSIMUS manual (p 169(equations) and 302(TRVP))
 * and TRVP article (Kolafa+Lisal 2011)
 * Inspired by the Nose–Hoover thermostat implementation in VerletNVTNoseTRVP.cpp
 */

#include <iostream>
#include <cmath>
#include <cstring>

#include "Vector.hpp"
#include "Matrix.hpp"
#include "Atom.hpp"
#include "Molecule.hpp"
#include "ConstraintBond.hpp"
#include "SimulatedSystem.hpp"
#include "VerletNPTNoseTRVP.hpp"
#include "general_utils.hpp"

#define ATOMR(i, j) ((*atomR)(i, j))
// #define DEBUG

// initialization of this integrator (actually not needed)
VerletNPTNoseTRVP::VerletNPTNoseTRVP(double step, SimulatedSystem *simsystem, double Tfinal, double tauT, double Pfinal, double tauP, bool startWithShake)
    : AbstractVerletIntegrator{step, simsystem, startWithShake}
{
    tau_T = tauT;
    T_f = Tfinal;
    tau_P = tauP;
    P_f = Pfinal;
    M_P = (simsystem->noDegrOfFred + 3) * T_f * tau_P * tau_P;
    M_T = simsystem->noDegrOfFred * T_f * tau_T * tau_T;
    recM_T = 1.0 / M_T;
    recM_P = 1.0 / M_P;

    // R initialization correction (R(2) = hv(i-1/2) instead of r(i-1))
    std::vector<Molecule>::iterator itmol;
    int i, k;
    double Ekininst;
    for (itmol = system->molecules.begin(); itmol != system->molecules.end(); itmol++)
    {
        for (i = 0; i < itmol->noAtoms; i++)
        {
            for (k = 0; k < 3; k++)
            {
                itmol->atoms[i].R->operator()(2, k) = itmol->atoms[i].R->operator()(0, k) - itmol->atoms[i].R->operator()(2, k);
            }
        }
    }
    positionsAt.clear();
    positionsAt.push_back(0);
    positionsAt.push_back(3);
    Ekininst = system->CalculateEkin(h, 1);

    Xi = new Vector(4);
    Lambda = new Vector(4);
    Xi_old = new Vector(*Xi);
    Lambda_old = new Vector(*Lambda);
    InitXi(predVel, Ekininst);
    InitLambda(predVel, Ekininst);
    simsystem->Eextra = EnergyOfExtraDOF();
}

// initialization of this integrator + TRVP (only this version is used)
VerletNPTNoseTRVP::VerletNPTNoseTRVP(double step, SimulatedSystem *simsystem, Vector trvpB, double Tfinal, double tauT, double Pfinal, double tauP, bool startWithShake)
    : AbstractVerletIntegrator{step, simsystem, trvpB, startWithShake}
{
    tau_T = tauT;
    T_f = Tfinal;
    tau_P = tauP;
    P_f = Pfinal;
    M_P = (simsystem->noDegrOfFred + 3) * T_f * tau_P * tau_P;
    M_T = simsystem->noDegrOfFred * T_f * tau_T * tau_T;
    recM_T = 1.0 / M_T;
    recM_P = 1.0 / M_P;

    // R initialization correction (R(2) = hv(i-1/2) instead of r(i-1))
    std::vector<Molecule>::iterator itmol;
    int i, k;
    double Ekininst;
    for (itmol = system->molecules.begin(); itmol != system->molecules.end(); itmol++)
    {
        for (i = 0; i < itmol->noAtoms; i++)
        {
            for (k = 0; k < 3; k++)
            {
                itmol->atoms[i].R->operator()(2, k) = itmol->atoms[i].R->operator()(0, k) - itmol->atoms[i].R->operator()(2, k);
            }
        }
    }
    positionsAt.clear();
    positionsAt.push_back(0);
    positionsAt.push_back(3);
    Ekininst = system->CalculateEkin(h, 1);

    Xi = new Vector(simsystem->GetRsize());
    Lambda = new Vector(simsystem->GetRsize());
    InitXi(predVel, Ekininst);
    VelocityVec(Xi);
    InitLambda(predVel, Ekininst);
    VelocityVec(Lambda);
    simsystem->Eextra = EnergyOfExtraDOF();
    system->Xi = Xi->operator()(0);
    system->Xivel = Xi->operator()(1) / h;
    system->Lambda = Lambda->operator()(0);
    system->Lambdavel = Lambda->operator()(1) / h;
    Xi_old = new Vector(*Xi);
    Lambda_old = new Vector(*Lambda);
}

// destructor
VerletNPTNoseTRVP::~VerletNPTNoseTRVP()
{
    delete Xi;
    delete Lambda;
    delete Xi_old;
    delete Lambda_old;
}

// perform N  steps of integration
int VerletNPTNoseTRVP::Integrate(int Nsteps)
{
    int l;
    double hvXiP = 0.0;        // predicted value of h*xi'(t)
    double hvLamP = 0.0;       // predicted value of h*lambda'(t)
    double boxscalePred = 1.0; // predicted value of box scaling for SHAKE
    double Ekininst = 0.0;     // instanteneous kinetic energy
    double Pinst = 0.0;        // instanteneous pressure
    double hvXi = 0.0;         // h*xi'(t)
    double Pkin = 0.0;         // kinetic pressure
    double Pvir = 0.0;         // virial pressure
    double boxscaling = 1.0;   // final box-scaling
    std::vector<Molecule>::iterator itmol;

    static double Ekinfintot2 = (system->noDegrOfFred + 1) * T_f;
    // static std::vector<int> positionsToRescale = {3};

    for (l = 0; l < Nsteps; l++)
    {
        Ekininst = 0.0;

        // time shift r(i-1) = r(i) from previous and r(i) = r(i+1) from previous step
        // must be before boxscaling in order not to affect by scaling the velocity v(i-1/2)
        TimeShift();

        // Box rescaling using Lambda(3)
        boxscaling = exp(Lambda->operator()(3) - Lambda->operator()(0));
        system->RescaleBox(boxscaling, positionsAt, false); // atom-based rescaling

        // time shift for Xi
        TimeShiftVec(Xi);

        // time shift for lambda
        TimeShiftVec(Lambda);

        // forces calculation, includes virconstr nullifying...
        system->CalculateForces();

        // TRVP for atomic velocities
        TRVP();
        // TRVP for Xi
        hvXiP = TRVPvec(Xi);
        // TRVP for Lambda
        hvLamP = TRVPvec(Lambda);
        // predict box scaling for SHAKE
        // lambdaP(t+h) = lambda(t) + h (2lambda'(t-h/2) - lambda'(t-3h/2))
        // boxLengthP(t+h)/boxLength(t) = exp(h (2lambda'(t-h/2) - lambda'(t-3h/2)))
        boxscalePred = exp(2.0 * Lambda->operator()(2) - Lambda->operator()(predVel + 1));

#ifdef PARALLEL
#pragma omp parallel for num_threads(thread_count) private(itmol) reduction(+: Ekininst)
#endif
        for (itmol = system->molecules.begin(); itmol != system->molecules.end(); itmol++)
        {
            // one step of Verlet (new position calculation), uses predicted value of h*xi' and h*lambda'
            VerletStep(&(*itmol), hvXiP, hvLamP);

            // SHAKE (shakeType = -3)
            Shake(&(*itmol), system->virconstr, boxscalePred); // accumulation to system->virconstr is controlled by omp atomic inside

            // new differences r(t+h) - r(t), kinetic temperature calculation (according to VERLET)
            Ekininst += VelocityCalculation(&(*itmol));
        }

        // calculate pressure (Pinst)
        // instanteneous pressure calculation
        Pkin = system->CalculatePkin(h, Ekininst);
        Pvir = system->CalculatePconf(); // includes Pcorr

        Pinst = Pkin + Pvir;

        if (system->boundaryCond > 0)
        {
            system->ApplyPeriodicBC(positionsAt);
        }

        // Verlet step for Xi
        // xi'' = 1/M_T * (2*Ekin(t) + M_P*lambda'(t)^2 - (N_f + 1) * k_B * T_f)
        Xi->operator()(3) = Xi->operator()(0) + Xi->operator()(2) +
                            h * h * recM_T * (2 * Ekininst - Ekinfintot2) + M_P * hvLamP * hvLamP * recM_T;
        // Xi velocity calculation (final velocity, not predicted)
        hvXi = VelocityVec(Xi);

        // Verlet step for Lambda
        Lambda->operator()(3) = Lambda->operator()(0) + Lambda->operator()(2) +
                                h * h * recM_P * 3.0 * (system->CalculateVolume() * (Pinst - P_f) + 2 * Ekininst / system->noDegrOfFred) - hvLamP * hvXi;
        // lambda velocity (not needed now, but to be updated regularly)
        VelocityVec(Lambda);
    }
    system->Eextra = EnergyOfExtraDOF();
    system->Xi = Xi->operator()(0);
    system->Xivel = Xi->operator()(1) / h;
    system->Lambda = Lambda->operator()(0);
    system->Lambdavel = Lambda->operator()(1) / h;

    return 0;
}

int VerletNPTNoseTRVP::VerletStep(Molecule *mol, double hvXiPred, double hvLamPred)
{
    int j, k;
    Matrix *atomR;

    for (j = 0; j < mol->noAtoms; j++)
    {
        atomR = mol->atoms[j].R;
        for (k = 0; k < 3; k++)
        {
            ATOMR(3, k) = ATOMR(0, k) + ATOMR(2, k) +
                          h * h * mol->atoms[j].force[k] / mol->atoms[j].mass -
                          ATOMR(predVel, k) * (hvXiPred + (1.0 + 3.0 / system->noDegrOfFred) * hvLamPred);
            // r(i+1) = 2r(i) - r(i-1) + h^2 * (f(i)/m - v(i)*(xi'(i) + 3*lambda'(i)/N_f + lambda'(i))
            // velocities (v(i), lambda'(i) and xi'(i)) use predicted values
        }
    }
    return 0;
}

// Energy of the extra degree of freedom
double VerletNPTNoseTRVP::EnergyOfExtraDOF() const
{
    // Eextra = Ekinextra + Epotextra =
    // 1/2*M_T*xi'^2 + 1/2*M_P*lambda'^2 + (N_f + 1)*k_B*T_f*xi + P_f * V
    return 0.5 * (M_T * (Xi->operator()(1) / h) * (Xi->operator()(1) / h) +
                  M_P * (Lambda->operator()(1) / h) * (Lambda->operator()(1) / h)) +
           (system->noDegrOfFred + 1.0) * T_f * Xi->operator()(0) +
           P_f * system->CalculateVolume();
}

// Initialize Xi
int VerletNPTNoseTRVP::InitXi(int predVel, double Ekininst)
{
    double force;
    int m;
    bool fromconfig = false;
    char aux[10];

    if ((system->extendedDOFs != nullptr) && (system->extendedDOFs->GetNumberOfRows() > 2))
    {
        fromconfig = true;
        sscanf(system->extDOFnames, "%10s", aux);
        if (strstr(aux, "Xi") == NULL)
        {
            print_warning(1, "Extended DOFs from .config don't start by Xi, cannot use info from .config file\n");
            fromconfig = false;
        }
    }

    if (fromconfig)
    {
        Xi->operator()(0) = system->extendedDOFs->operator()(0, 0);
        Xi->operator()(1) = system->extendedDOFs->operator()(1, 0) * h;
        if (system->extendedDOFs->operator()(2, 0) == 0.0) // if empty, use differential equation
        {
            force = h * h * recM_T * (2 * Ekininst + (system->noDegrOfFred) * T_f);
            Xi->operator()(2) = 0.5 * force;
        }
        else
        {
            Xi->operator()(2) = h * h * 0.5 * system->extendedDOFs->operator()(2, 0);
        }
        Xi->operator()(3) = Xi->operator()(0) + Xi->operator()(1) + Xi->operator()(2); // r(i+1) = r(i) + h*v(i-1/2) + h^2/2*f(i)/m (for NVE, more complicated for other ensembles)

        if (system->extendedDOFs->GetNumberOfRows() > 3) // expected true
        {
            Xi->operator()(4) = system->extendedDOFs->operator()(4, 0) * h;
            Xi->operator()(2) = Xi->operator()(0) - Xi->operator()(4); // if TRVP, than this value needed...
        }
        else
        {
            Xi->operator()(2) = Xi->operator()(0) - Xi->operator()(1);
        }

        // TRVP – stored past velocities
        for (m = predVel + 1; m < Xi->GetSize(); m++)
        {
            if ((m < system->extendedDOFs->GetNumberOfRows()) && (system->extendedDOFs->operator()(m, 0) != 0.0))
            {
                Xi->operator()(m) = system->extendedDOFs->operator()(m, 0) * h;
            }
            else
            {
                Xi->operator()(m) = Xi->operator()(1) - (m - predVel + 0.5) * 0.5 * h * h / tau_T / tau_T * (Ekininst / (0.5 * system->noDegrOfFred * T_f) - 1.0);
            }
        }
        // this integrator specific Xi(2) (see R arrange in the beginning of this file)
        Xi->operator()(2) = Xi->operator()(0) - Xi->operator()(2);
    }
    else
    {
        // calculate xi'' (force on xi)
        force = h * h * recM_T * (2 * Ekininst + (system->noDegrOfFred) * T_f);

        // initialize past values of Xi
        Xi->operator()(2) = Xi->operator()(1) - 0.5 * force;
        Xi->operator()(3) = Xi->operator()(0) + Xi->operator()(1) + 0.5 * force;
        // TRVP – stored past velocities
        if (predVel > 1)
        {
            Xi->operator()(predVel) = Xi->operator()(1);
            for (m = predVel + 1; m < Xi->GetSize(); m++)
            {
                if (Xi->operator()(m) != 0.0)
                {
                    Xi->operator()(m) *= h;
                }
                else
                {
                    Xi->operator()(m) = Xi->operator()(1) - (m - predVel + 1) * 0.5 * force;
                }
            }
        }
    }

    return 0;
}

// Initialization of Lambda
int VerletNPTNoseTRVP::InitLambda(int predVel, double Ekininst)
{
    double force, V, Pkin, Pvir, Pinst;
    int m;
    bool fromconfig = false;
    char aux1[10], aux2[10], aux3[10];
    int lambda_pos;

    // calculate pressure
    V = system->CalculateVolume();
    Pkin = system->CalculatePkin(h, Ekininst);
    Pvir = system->CalculatePconf();
    Pinst = Pkin + Pvir;

    if ((system->extendedDOFs != nullptr) && (system->extendedDOFs->GetNumberOfRows() > 2))
    {
        fromconfig = true;
        sscanf(system->extDOFnames, "%10s %10s %10s", aux1, aux2, aux3);
        if (strstr(aux2, "Lambda") != NULL)
        {
            lambda_pos = 1;
        }
        else if (strstr(aux3, "Lambda") != NULL)
        {
            lambda_pos = 2;
        }
        else
        {
            print_warning(1, "Extended DOFs from .config don't define Lambda, cannot use info from .config file\n");
            fromconfig = false;
            delete system->extendedDOFs;
            system->extendedDOFs = nullptr;
        }
    }

    if (fromconfig)
    {
        Lambda->operator()(0) = system->extendedDOFs->operator()(0, lambda_pos);
        Lambda->operator()(1) = system->extendedDOFs->operator()(1, lambda_pos) * h;
        if (system->extendedDOFs->operator()(2, lambda_pos) == 0.0) // if empty, use differential equation
        {
            force = 3.0 * h * h * recM_P * (V * (Pinst - P_f) + 2.0 * Ekininst / system->noDegrOfFred) - Xi->operator()(1) * Lambda->operator()(1);
            Lambda->operator()(2) = 0.5 * force;
        }
        else
        {
            Lambda->operator()(2) = h * h * 0.5 * system->extendedDOFs->operator()(2, lambda_pos);
            if (fabs(Lambda->operator()(2)) > h * h * 100) // maybe pressure was stored on this position by MTTK/Alejandre version of MTK NPT
            {
                print_warning(1, "Large value stored in Lambda(2) (" + std::to_string(system->extendedDOFs->operator()(2, lambda_pos)) + ").\n",
                "Check if you start from the configuration (.config) created by VerletNPTNoseTRVP integrator.\n",
                "Simulation may crash.\n");
            }
        }
        Lambda->operator()(3) = Lambda->operator()(0) + Lambda->operator()(1) + Lambda->operator()(2); // r(i+1) = r(i) + h*v(i-1/2) + h^2/2*f(i)/m (for NVE, more complicated for other ensembles)

        if (system->extendedDOFs->GetNumberOfRows() > 3) // expected true
        {
            Lambda->operator()(4) = system->extendedDOFs->operator()(4, lambda_pos) * h;
            Lambda->operator()(2) = Lambda->operator()(0) - Lambda->operator()(4); // if TRVP, than this value needed...
        }
        else
        {
            Lambda->operator()(2) = Lambda->operator()(0) - Lambda->operator()(1);
        }

        // TRVP – stored past velocities
        for (m = predVel + 1; m < Xi->GetSize(); m++)
        {
            if ((m < system->extendedDOFs->GetNumberOfRows()) && (system->extendedDOFs->operator()(m, lambda_pos) != 0.0))
            {
                Lambda->operator()(m) = system->extendedDOFs->operator()(m, lambda_pos) * h;
            }
            else
            {
                force = 3.0 * h * h * recM_P * (V * (Pinst - P_f) + 2.0 * Ekininst / system->noDegrOfFred) - Xi->operator()(1) * Lambda->operator()(1);
                Lambda->operator()(m) = Lambda->operator()(1) - (m - predVel + 0.5) * 0.5 * force;
            }
        }

        // this integrator specific (see the R design for this integrator)
        Lambda->operator()(2) = Lambda->operator()(0) - Lambda->operator()(2);
        delete system->extendedDOFs;
        system->extendedDOFs = nullptr;
        // }
    }
    else
    {
        // calculate lambda'' (force on lambda)
        force = 3.0 * h * h * recM_P * (V * (Pinst - P_f) + 2.0 * Ekininst / system->noDegrOfFred) - Xi->operator()(1) * Lambda->operator()(1);

        // initialize values of Lambda
        // Lambda->operator()(0) = log(V)/3.0;
        Lambda->operator()(0) = 0.0;
        Lambda->operator()(2) = Lambda->operator()(1) - 0.5 * force;
        Lambda->operator()(3) = Lambda->operator()(0) + Lambda->operator()(1) + 0.5 * force;
        // TRVP – stored past velocities
        Lambda->operator()(predVel) = Lambda->operator()(1);
        for (m = predVel + 1; m < Lambda->GetSize(); m++)
        {
            if (Lambda->operator()(m) != 0.0)
            {
                Lambda->operator()(m) *= h;
            }
            else
            {
                Lambda->operator()(m) = Lambda->operator()(1) - (m - predVel + 1) * 0.5 * force;
            }
        }
    }

    return 0;
}

// time-shift for vectors (Xi and Lambda)
int VerletNPTNoseTRVP::TimeShiftVec(Vector *vec)
{
    int i;
    // TRVP values (differences h*vec')
    for (i = predVel + trvpCoeff->GetSize() - 2; i > predVel; i--)
    {
        vec->operator()(i + 1) = vec->operator()(i);
    }
    vec->operator()(predVel + 1) = vec->operator()(2);
    // normal TimeShift
    vec->operator()(2) = vec->operator()(3) - vec->operator()(0); // hv(i-1/2) = r(i+1) - r(i) from previous
    vec->operator()(0) = vec->operator()(3);                      // r(i) = r(i+1) from previous
    return 0;
}

// TRVP for vectors (Xi and Lambda)
double VerletNPTNoseTRVP::TRVPvec(Vector *vec)
{
    int i;
    // TRVP
    vec->operator()(predVel) = trvpCoeff->operator()(0) * vec->operator()(2);
    for (i = 1; i < trvpCoeff->GetSize(); i++)
    {
        vec->operator()(predVel) += trvpCoeff->operator()(i) * vec->operator()(predVel + i);
    }
    return vec->operator()(predVel);
}

// velocity of the extended degree of freedom
double VerletNPTNoseTRVP::VelocityVec(Vector *vec)
{
    double vhm, vhp;

    vhp = vec->operator()(3) - vec->operator()(0); // hv(t) = [r(t + h) − r(t)]
    vhm = vec->operator()(2);                      // hv(t) = [r(t) − r(t - h)] (already in vec(2));

#if (VERLET == 0) // hv(t) = [r(t + h) − r(t)]
    vec->operator()(1) = vhp;
#elif (VERLET == 1) // hv(i) = (r(i+1)-r(i-1))/2 - equivalent to k3m2e, velocity Verlet
    vec->operator()(1) = (vhp + vhm) * 0.5;
#elif (VERLET == 2) // (hv(t))^2 = [r(t)−r(t−h)]·[r(t+h)−r(t)]
    vec->operator()(1) = sqrt(fabs(vhm * vhp)) * ((vhp > 0) - (vhp < 0));
#else               // VERLET == 3  energy - averages of both shifted values v(t) = [r(t + h) − r(t)]/h and v(t) = [r(t) − r(t − h)]/h.
    vec->operator()(1) = sqrt((vhm * vhm + vhp * vhp) * 0.5) * ((vhp > 0) - (vhp < 0));
#endif
    return vec->operator()(1);
}

// Calculate velocity according to VERLET #define and return kinetic energy
double VerletNPTNoseTRVP::VelocityCalculation(Molecule *mol)
{
    // double Ekininsthp = 0.0;
    double Ekininst = 0.0;
    Matrix *atomR;
    int j, k;
    double vhm, vhp;

    for (j = 0; j < mol->noAtoms; j++)
    {
        atomR = mol->atoms[j].R;
        for (k = 0; k < 3; k++)
        {
            vhp = ATOMR(3, k) - ATOMR(0, k); // hv(t) = [r(t + h) − r(t)]
            vhm = ATOMR(2, k);               // hv(t) = [r(t) − r(t - h)]; // already in ATOMR(2)
                                             // Ekininsthp += mol->atoms[j].mass * vhp * vhp / (h * h);

#if (VERLET == 0) // hv(t) = [r(t + h) − r(t)]
            ATOMR(1, k) = vhp;
#elif (VERLET == 1) // hv(i) = (r(i+1)-r(i-1))/2 - equivalent to k3m2e, velocity Verlet
            ATOMR(1, k) = (vhp + vhm) * 0.5;
#elif (VERLET == 2) // (hv(t))^2 = [r(t)−r(t−h)]·[r(t+h)−r(t)]
            ATOMR(1, k) = sqrt(fabs(vhm * vhp)) * ((vhp > 0) - (vhp < 0));
#else               // VERLET == 3  energy - averages of both shifted values v(t) = [r(t + h) − r(t)]/h and v(t) = [r(t) − r(t − h)]/h.
            ATOMR(1, k) = sqrt((vhm * vhm + vhp * vhp) * 0.5) * ((vhp > 0) - (vhp < 0));
#endif

            Ekininst += mol->atoms[j].mass * ATOMR(1, k) * ATOMR(1, k) / h / h;
        }
    }

    // return Ekininsthp * 0.5;
    return Ekininst * 0.5;
}

int VerletNPTNoseTRVP::TRVP()
{
    // implementation of time-reversible velocity prediction (Kolafa+Lisal:JCTC 2011) for Verlet
    // version which saves differences (h*v) is used with coefficients B (in original paper)
    // R(2) = r(t) - r(t-h), R(5) = r(t-h) - r(t-2h), R(6) = r(t-2h) - r(t-3h), ...
    int i, j, k, l;
    Matrix *atomR;

    for (i = 0; i < system->noMolecules; i++)
    {
        for (j = 0; j < system->molecules[i].noAtoms; j++)
        {
            atomR = system->molecules[i].atoms[j].R;
            for (k = 0; k < 3; k++)
            {
                ATOMR(predVel, k) = trvpCoeff->operator()(0) * ATOMR(2, k);
                for (l = 1; l < trvpCoeff->GetSize(); l++)
                {
                    ATOMR(predVel, k) += trvpCoeff->operator()(l) * ATOMR(predVel + l, k);
                }
            }
        }
    }
    return 0;
}

int VerletNPTNoseTRVP::TimeShift()
{
    Matrix *atomR;
    int i, j, k, l;

    if (predVel == 1)
    {
        for (i = 0; i < system->noMolecules; i++)
        {
            for (j = 0; j < system->molecules[i].noAtoms; j++)
            {
                atomR = system->molecules[i].atoms[j].R;

                for (k = 0; k < 3; k++)
                {
                    ATOMR(2, k) = ATOMR(3, k) - ATOMR(0, k); // hv(i-1/2) = r(i+1) - r(i) from previous
                    ATOMR(0, k) = ATOMR(3, k);               // r(i) = r(i+1) from previous
                }
            }
        }
    }
    else
    {
        for (i = 0; i < system->noMolecules; i++)
        {
            for (j = 0; j < system->molecules[i].noAtoms; j++)
            {
                atomR = system->molecules[i].atoms[j].R;

                for (k = 0; k < 3; k++)
                {
                    // TRVP values (differences h*v)
                    for (l = predVel + trvpCoeff->GetSize() - 2; l > predVel; l--)
                    {
                        ATOMR(l + 1, k) = ATOMR(l, k);
                    }
                    ATOMR(predVel + 1, k) = ATOMR(2, k);
                    // normal Verlet TimeShift
                    ATOMR(2, k) = ATOMR(3, k) - ATOMR(0, k); // hv(i-1/2) = r(i+1) - r(i) from previous
                    ATOMR(0, k) = ATOMR(3, k);               // r(i) = r(i+1) from previous
                }
            }
        }
    }

    return 0;
}

// Prepare Xi and Lambda for .config printing
int VerletNPTNoseTRVP::PrepareConfig()
{
    int configLevel;
    double temp = 0.0;
    std::vector<Molecule>::iterator itmol;
    int i, k;

    (*Xi_old) = (*Xi);
    (*Lambda_old) = (*Lambda);
    
    // this integrator specific rearrange
    Lambda->operator()(2) = Lambda->operator()(0) - Lambda->operator()(2);
    Xi->operator()(2) = Xi->operator()(0) - Xi->operator()(2);
    // R initialization correction (R(2) = hv(i-1/2) instead of r(i-1))
    for (itmol = system->molecules.begin(); itmol != system->molecules.end(); itmol++)
    {
        for (i = 0; i < itmol->noAtoms; i++)
        {
            for (k = 0; k < 3; k++)
            {
                itmol->atoms[i].R->operator()(2, k) = itmol->atoms[i].R->operator()(0, k) - itmol->atoms[i].R->operator()(2, k);
            }
        }
    }

    // normal PrepareConfig
    configLevel = AbstractVerletIntegrator::PrepareConfig();

    temp = Xi->operator()(4);
    Xi->operator()(4) = Xi->operator()(0) - Xi->operator()(2);

    Xi->operator()(1) = Xi->operator()(0) - Xi->operator()(2);
    Xi->operator()(2) = Xi->operator()(3) - Xi->operator()(0) - Xi->operator()(1); // to enable precise continue in whatever ensemble

    Xi->operator()(3) = temp;

    temp = Lambda->operator()(4);
    Lambda->operator()(4) = Lambda->operator()(0) - Lambda->operator()(2);

    Lambda->operator()(1) = Lambda->operator()(0) - Lambda->operator()(2);
    Lambda->operator()(2) = Lambda->operator()(3) - Lambda->operator()(0) - Lambda->operator()(1); // to enable precise continue in whatever ensemble

    Lambda->operator()(3) = temp;

    system->extendedDOFs = new Matrix(Xi->GetSize(), 2);
    for (i = 0; i < Xi->GetSize(); i++)
    {
        system->extendedDOFs->operator()(i, 0) = Xi->operator()(i);
        system->extendedDOFs->operator()(i, 1) = Lambda->operator()(i);
    }

    configLevel += 200;
    strcpy(system->extDOFnames, "Xi Lambda");

    return configLevel;
}

int VerletNPTNoseTRVP::RestoreExtendedDOFs()
{
    (*Lambda) = (*Lambda_old);
    (*Xi) = (*Xi_old);
    return 0;
}
