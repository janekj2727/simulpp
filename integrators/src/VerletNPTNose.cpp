/*
 * An implementation of Verlet integration method (Velocity Verlet + Nose–Hoover thermostat)
 * "MTTK" – Trotter decomposition scheme, without TRVP (with TRVP the algorithm is completely different – defined in VerletNPTNoseTRVP.cpp)
 * for trial simulations of contrained dynamics
 * Author JJ, Date Oct 2021
 */

/*
 * In atoms R following values are stored:
 * R(0, :): x(i) y(i) z(i)             common for all integrators
 * R(1, :): h*vx(i) h*vy(i) h*vz(i)    common for all integrators (h* because of Gear)
 * R(2, :): x(i-1) y(i-1) z(i-1)       Verlet specific
 * R(3, :): x(i-2) y(i-2) z(i-2)       meanless, this integrator specific
 *
 * R(0, :) and R(1, :) must be defined before integrator initialization
 * note velocity is `hard-coded` in this algorithm
 */

/*
 * In this integrator MACSIMUS versions of VERLET velocity does not apply.
 * Velocities are uniquely determined by the integration algorithm (and correspond to VERLET=1).
 */

/*
 *  Nosé–Hoover thermostat implemented according to the DL POLY manual p. 98–103
 *  Checked against Frenkel&Smit ???
 *  MACSIMUS notation used...
 *  See also notes in control.md
 */

#include <iostream>
#include <cmath>
#include <cassert>
#include <cstring>

#include "Vector.hpp"
#include "Matrix.hpp"
#include "math_utils.hpp"
#include "Atom.hpp"
#include "Molecule.hpp"
#include "ConstraintBond.hpp"
#include "SimulatedSystem.hpp"
#include "VerletNPTNose.hpp"
#include "general_utils.hpp"

#define ATOMR(i, j) ((*atomR)(i, j))
// #define DEBUG
// #define GRADLAMBDAUPDATE // no effect

// initialization of this integrator (main part is done in AbstractVerletIntegrator construction)
VerletNPTNose::VerletNPTNose(double step, SimulatedSystem *simsystem, double Tfinal, double tauT, double Pfinal, double tauP, bool startWithShake)
    : AbstractVerletIntegrator{step, simsystem, startWithShake, true}
{
    int i;

    T_f = Tfinal;
    tau_T = tauT;
    P_f = Pfinal;
    tau_P = tauP;
    M_T = simsystem->noDegrOfFred * T_f * tau_T * tau_T;
    M_P = (simsystem->noDegrOfFred + 3) * T_f * tau_P * tau_P;
    // predVel = 1; // already done in AbstractIntegrator

    // Extra degree of freedom – variable connected with heat bath
    Xi = new Vector(4);
    Lambda = new Vector(4);
    for (i = 0; i < 3; i++)
    {
        orig_size[i] = simsystem->GetBox(i);
    }
    Initialize(startWithShake); // initialize R here (using simpleinit version of AbstractVerletIntegrator), argument has no meaning
    InitXi();
    InitLambda();
    simsystem->Eextra = EnergyOfExtraDOF();
    Xi_old = new Vector(*Xi);
    Lambda_old = new Vector(*Lambda);
}

// empty destructor
VerletNPTNose::~VerletNPTNose()
{
    delete Xi;
    delete Lambda;
    delete Xi_old;
    delete Lambda_old;
}

// perform N  steps of integration
int VerletNPTNose::Integrate(int Nsteps)
{
    int i, l;
    Molecule *p_mol;
    double Ekininst = 0.0; // instanteneous kinetic energy
    double virialc = 0.0;  // virial of constraining forces
    double velscale = 1.0; // velocity scaling factor

    // calculate kinetic energy (initialization of Ekininst)
    // during integration the kinetic energy is updated
    Ekininst = system->CalculateEkin(h, 1);

    for (l = 0; l < Nsteps; l++)
    {
        virialc = 0.0; // virial of constraining forces
        velscale = 1.0;

        // time shift r(i-1) = r(i) from previous and r(i) = r(i+1) from previous step
        // special version of TimeShift for this integrator
        TimeShift();
        // time shift xi... as for r
        Xi->operator()(3) = Xi->operator()(2);
        Xi->operator()(2) = Xi->operator()(0);
        // time shift lambda ... as for xi
        Lambda->operator()(3) = Lambda->operator()(2);
        Lambda->operator()(2) = Lambda->operator()(0);

        // thermostat h/4
        ThermostatQuarterStep(Ekininst, velscale);

        // barostat h/2
        BarostatHalfStep(Ekininst, velscale);

        // thermostat h/4
        ThermostatQuarterStep(Ekininst, velscale);

        // scale velocities (scaling factor accumulated form thermostat and barostat)
        system->RescaleVelocities(velscale);
        velscale = 1.0;

        // rescale box (update volume), atom based rescaling...
        system->RescaleBox(exp(Lambda->operator()(1)), positionsAt, false);

#ifdef PARALLEL
#pragma omp parallel for num_threads(thread_count) private(p_mol, i)
#endif
        for (i = 0; i < system->noMolecules; i++)
        {
            p_mol = &(system->molecules[i]);

            // VV1 (velocity update with old forces)
            VerletStep(p_mol);

            // RATTLE1
            // bond lengths correction
            // virialc inside protected by omp atomic
            Shake(p_mol, virialc);
        }

        // Verlet step for lambda (in DL POLY notation: integration of eta)
        // not needed actually...
        Lambda->operator()(0) += Lambda->operator()(1);

        // force calculation f(t+h) - parallelized inside
        system->CalculateForces();

#ifdef PARALLEL
#pragma omp parallel for num_threads(thread_count) private(p_mol, i)
#endif
        for (i = 0; i < system->noMolecules; i++)
        {
            p_mol = &(system->molecules[i]);

            // VV2
            // recalculate Ekininst
            VerletStep2(p_mol);

            // RATTLE2
            // velocity perpendicular...
            Rattle(p_mol, virialc);
        }

        // add virial of constraint forces to the system virial of constr. forces
        system->virconstr = virialc / h / h;

        // recalculate Ekininst (after VerletStep2 and Rattle)
        Ekininst = system->CalculateEkin(h, 1);

        // thermostat h/4 (same as the 1st step)
        ThermostatQuarterStep(Ekininst, velscale);

        // barostat h/2
        BarostatHalfStep(Ekininst, velscale);

        // thermostat h/4
        ThermostatQuarterStep(Ekininst, velscale);

        // finally, rescale velocities to final values
        system->RescaleVelocities(velscale);
        velscale = 1.0;

        // boundary conditions
        if (system->boundaryCond > 0)
        {
            system->ApplyPeriodicBC(positionsAt);
        }
    }

    // calculate energy of extra DOF and pass it to system
    system->Eextra = EnergyOfExtraDOF();
    // pass thermostat and barostat variables to system for printing in measurement
    system->Xi = Xi->operator()(0);
    system->Xivel = Xi->operator()(1) / h;
    system->Lambda = Lambda->operator()(0);
    system->Lambdavel = Lambda->operator()(1) / h;

    return 0;
}

// calculate energy of the extra degree of freedom (xi)
double VerletNPTNose::EnergyOfExtraDOF() const
{
    // Eextra = Ekinextra + Epotextra =
    // 1/2*M_T*xi'^2 + 1/2*M_P*lambda'^2 + (N_f + 1)*k_B*T_f*xi + P_f * V
    return 0.5 * (M_T * (Xi->operator()(1) / h) * (Xi->operator()(1) / h) +
                  M_P * (Lambda->operator()(1) / h) * (Lambda->operator()(1) / h)) +
           (system->noDegrOfFred + 1.0) * T_f * Xi->operator()(0) +
           P_f * system->CalculateVolume();
}

// quarter step of thermostat
int VerletNPTNose::ThermostatQuarterStep(double &Ekin, double &velscale)
{
    double scaling;
    static double Ekinfinaltotal = 0.5 * (system->noDegrOfFred + 1) * T_f; // including Ekin of barostat

    // integration of xi here (by h/16) in analogy with NVTNose implementation in DL POLY
    Xi->operator()(0) += Xi->operator()(1) / 16.0;

    // xi'(t+h/8) = xi'(t) + h/8*(2Ekininst - 2Ekinfinal + M_P*v_lambda^2 - kB*T_final)/M_T
    // M_P*v_lambda^2 ... 2*kinetic energy connected with lambda
    // kB*T_f ... 2*desired kinetic energy of lambda
    // Ekinfinal + kB*T_f = Ekinfinaltotal
    Xi->operator()(1) += h * h / 8.0 * (2 * Ekin + M_P * Lambda->operator()(1) * Lambda->operator()(1) / h / h - 2.0 * Ekinfinaltotal) / M_T;

    // integration of xi (by h/8) in analogy with NVT Nose...
    Xi->operator()(0) += Xi->operator()(1) / 8.0;

    // v(t) = v(t) * exp(-xi'(t+h/8)*h/4)
    //     (E_kininst(new) = 0.5 * sum(m*(v*exp(-xi'(t+h/4)*h/2))^2)
    //      E_kininst(new) = E_kininst * exp(-2*xi'(t+h/4)*h/2))
    scaling = exp(-Xi->operator()(1) / 4.0);

    // system->RescaleVelocities(scaling); // velocities can be rescaled later here only save scaling
    velscale *= scaling;
    // kinetic energy rescaled immediately
    Ekin *= scaling * scaling;

    // xi'(t+h/8) = xi'(t) + h/8*(2Ekininst - 2Ekinfinal + M_P*v_lambda^2 - kB*T_final)/M_T
    // M_P*v_lambda^2 ... 2*kinetic energy connected with lambda
    // kB*T_f ... 2*desired kinetic energy of lambda
    // Ekinfinal + kB*T_f = Ekinfinaltotal
    Xi->operator()(1) += h * h / 8.0 * (2 * Ekin + M_P * Lambda->operator()(1) * Lambda->operator()(1) / h / h - 2.0 * Ekinfinaltotal) / M_T;

    // in analogy with NVT Nose, the 3rd part of xi integration (by h/16)
    Xi->operator()(0) += Xi->operator()(1) / 16.0;

    return 0;
}

// half step of barostat
int VerletNPTNose::BarostatHalfStep(double &Ekin, double &velscale)
{
    double scaling;
    double V = system->CalculateVolume();
    double Pvir, Pkin; // virial and kinetic pressure

    // #ifdef GRADLAMBDAUPDATE // switch to analogy with termostat variable integration
    //     // integration of lambda here (by h/8) in analogy with thermostat
    //     // not done in DL_POLY, integration of lambda is performed directly by V(t+h) = V(t) * exp(3*lambda_v(t+h/2)* h)
    //     Lambda->operator()(0) += Lambda->operator()(1) / 8.0;
    //     // volume update to actual value
    //     V *= exp(3.0 * Lambda->operator()(0));
    // #endif

    // lambda'(t) *= exp(-xi'(t+h/4)*h/8)
    Lambda->operator()(1) *= exp(-Xi->operator()(1) / 8.0);

    // instanteneous pressure calculation
    Pkin = system->CalculatePkin(h, Ekin);
    Pvir = system->CalculatePconf();
    // #ifdef CUTOFFCORR
    //     // Pkin = (Tkin * (3*noAtoms - noConstr)) / (3.0 * V)
    //     // Tkin = 2Ekin/N_f
    //     Pkin = (2.0 * Ekin * (3.0 * system->noAtomsTotal - system->noConstrTotal) / (3.0 * V * system->noDegrOfFred)); // best agreement with MACSIMUS
    //     Pcorr = system->EnergyCutoffCorr() / V;                                                                        // not exact with gradual lambda update... !!!
    // #else
    //     Pkin = (2.0 * Ekin) / (3.0 * V); // best agreement with DL POLY
    //     // Pkin = 2.0 * Ekin * (3.0 * system->noAtomsTotal - system->noConstrTotal) / (3.0 * V * system->noDegrOfFred);
    //     Pcorr = 0.0;
    // #endif
    //     Pvir = -vir / (3.0 * V);

    // lambda'(t+h/4) = lambda'(t) + 3*h/4*((P(t) - Pfinal)*V/M_P + 2Ekin/N_f/M_P)
    // P(t) = Pkin + Pvir + Pcorr
    Lambda->operator()(1) += 3.0 * h * h / 4.0 * ((Pkin + Pvir - P_f) * V + 2.0 * Ekin / system->noDegrOfFred) / M_P;

    // lambda'(t+h/4) *= exp(-xi'(t+h/4)*h/8)
    Lambda->operator()(1) *= exp(-Xi->operator()(1) / 8.0);

    // velocity scaling
    // v *= exp(-(1+3/N_f)*lambda'(t+h/4)*h/2)
    // velocities can be rescaled later here only save scaling and update Ekin
    scaling = exp(-(1.0 + 3.0 / system->noDegrOfFred) * Lambda->operator()(1) / 2.0);
    velscale *= scaling;
    Ekin *= scaling * scaling;

    // integration lambda here??? or later...
    // #ifdef GRADLAMBDAUPDATE // switch to analogy with termostat variable integration
    //     // integration of lambda here (by h/8) in analogy with thermostat
    //     // not done in DL_POLY, integration of lambda is performed directly by V(t+h) = V(t) * exp(3*lambda_v(t+h/2)* h)
    //     Lambda->operator()(0) += Lambda->operator()(1) / 4.0;
    //     // volume update to actual value
    //     V *= exp(3.0 * Lambda->operator()(1) / 4.0);
    //     Pvir = -vir / (3.0 * V);
    // #endif

    // lambda'(t+h/4) *= exp(-xi'(t+h/4)*h/8)
    Lambda->operator()(1) *= exp(-Xi->operator()(1) / 8.0);

    // instanteneous pressure calculation kinetic part
    Pkin = system->CalculatePkin(h, Ekin);

    // lambda'(t+h/2) = lambda'(t+h/4) + 3*h/4*((P(t) - Pfinal)*V/M_P + 2Ekin/N_f/M_P)
    // P(t) = Pkin + Pvir + Pcorr
    Lambda->operator()(1) += 3.0 * h * h / 4.0 * ((Pkin + Pvir - P_f) * V + 2 * Ekin / system->noDegrOfFred) / M_P;

    // lambda'(t+h/2) *= exp(-xi'(t+h/4)*h/8)
    Lambda->operator()(1) *= exp(-Xi->operator()(1) / 8.0);

    // #ifdef GRADLAMBDAUPDATE // switch to analogy with termostat variable integration
    //     // integration of lambda here (by h/8) in analogy with thermostat
    //     // not done in DL_POLY, integration of lambda is performed directly by V(t+h) = V(t) * exp(3*lambda_v(t+h/2)* h)
    //     Lambda->operator()(0) += Lambda->operator()(1) / 8.0;
    //     // volume update to actual value
    //     V *= exp(3.0 * Lambda->operator()(1) / 8.0);
    // #endif

    return 0;
}

// VV1 – 1st part of velocity verlet
int VerletNPTNose::VerletStep(Molecule *p_mol)
{
    int j, k;
    Matrix *atomR;
    // double expL = exp(Lambda->operator()(1));
    // double expL = 1.0;

    for (j = 0; j < p_mol->noAtoms; j++)
    {
        atomR = p_mol->atoms[j].R;
        for (k = 0; k < 3; k++)
        {
            // v(t+h/2) = v(t) + h/2 * f(t)/m
            ATOMR(1, k) += h * h * 0.5 * p_mol->atoms[j].force[k] / p_mol->atoms[j].mass;
            // r(t+h) = exp(lambda'(t+h/2) * h) * r(t) + h * v(t+h/2)
            // scaling done before with the entire box
            ATOMR(0, k) += ATOMR(1, k);
        }
    }

    return 0;
}

// VV2 – 2nd part of velocity verlet
int VerletNPTNose::VerletStep2(Molecule *p_mol)
{
    int j, k;
    Matrix *atomR;

    for (j = 0; j < p_mol->noAtoms; j++)
    {
        atomR = p_mol->atoms[j].R;
        for (k = 0; k < 3; k++)
        {
            // v(t+h) = v(t+h/2) + h/2 * f(t+h)/m
            // ATOMR(1, k) += h * h * 0.5 * p_mol->atoms[j].force[k]/p_mol->atoms[j].mass;  // original version (works without SHAKE)
            ATOMR(1, k) = ATOMR(0, k) - ATOMR(2, k) + h * h * 0.5 * p_mol->atoms[j].force[k] / p_mol->atoms[j].mass; // with SHAKE
        }
    }

    return 0;
    // return p_mol->CalculateEkin(h, 1);
}

// TimeShift version for this integrator (substantially different from the others)
int VerletNPTNose::TimeShift()
{
    Matrix *atomR;
    int i, j, k;

    for (i = 0; i < system->noMolecules; i++)
    {
        for (j = 0; j < system->molecules[i].noAtoms; j++)
        {
            atomR = system->molecules[i].atoms[j].R;

            for (k = 0; k < 3; k++)
            {
                ATOMR(3, k) = ATOMR(2, k); // r(i-2) = r(i-1) from previous
                ATOMR(2, k) = ATOMR(0, k); // r(i-1) = r(i) from previous
            }
        }
    }
    // TRVP version not needed (just for sure make assertion...)
    assert(predVel == 1);
    return 0;
}

// SHAKE call RATTLE 1 version
double VerletNPTNose::Shake(Molecule *mol, double &virial)
{
    double maxRelErr = 0.0;
    static const int iterShake = 100; // maximum number of SHAKE iterations (currently fixed number)
    int m, i, j, k;
    Matrix *atomR;
    double relErr = 0.0;
    double omega = (mol->noConstrBonds > 1) ? (system->omegaShake) : (1.0); // overrelaxation parameter
    double local_virial = 0.0;

    for (m = 0; m < iterShake; m++)
    {
        maxRelErr = 0.0;

        // for each bond calculate contribution to constr forces
        for (i = 0; i < mol->noConstrBonds; i++)
        {
            if (relErr = mol->constrBonds[i].CalculateForces(local_virial, -1, system->epsShake, 1.0, omega), relErr > maxRelErr)
            {
                maxRelErr = relErr;
            }

            j = mol->constrBonds[i].atomI;
            atomR = mol->atoms[j].R;
            for (k = 0; k < 3; k++)
            {
                ATOMR(0, k) += mol->atoms[j].constraintForces[k];
                ATOMR(1, k) = ATOMR(0, k) - ATOMR(2, k); // change appropriately to r(i+1) - r(i)
            }
            j = mol->constrBonds[i].atomJ;
            atomR = mol->atoms[j].R;
            for (k = 0; k < 3; k++)
            {
                ATOMR(0, k) += mol->atoms[j].constraintForces[k];
                ATOMR(1, k) = ATOMR(0, k) - ATOMR(2, k); // change appropriately to r(i+1) - r(i)
            }
        }
#ifdef PARALLEL
#pragma omp atomic
#endif
        virial += local_virial;
        local_virial = 0.0;

        if (maxRelErr < system->epsShake)
        {
            (m > system->maxShakeIter) ? (system->maxShakeIter = m) : (true);
            return maxRelErr;
        }
    }

    (m > system->maxShakeIter) ? (system->maxShakeIter = m) : (true);
    return maxRelErr;
}

// RATTLE to correct velocity (to be perpendicular to the constr. bond)
double VerletNPTNose::Rattle(Molecule *p_mol, double &virial)
{
    double maxRelErr = 0.0;
    static const int iterShake = 100; // maximum number of RATTLE iterations (currently fixed number)
    int m, i, j, k;
    Matrix *atomR;
    double relErr;
    double omega = (p_mol->noConstrBonds > 1) ? (system->omegaShake) : (1.0); // overrelaxation parameter
    double local_virial = 0.0;

    // what should be the value of eps and how it should be calculated???
    // for diatomics one iteration is enough, for larger molecules this could be a problem...

    for (m = 0; m < iterShake; m++)
    {
        maxRelErr = 0.0;

        // for each bond calculate contribution to constr forces
        for (i = 0; i < p_mol->noConstrBonds; i++)
        {
            if (relErr = p_mol->constrBonds[i].CalculateForces(local_virial, -2, system->epsShake, 1.0, omega), relErr > maxRelErr)
            {
                maxRelErr = relErr;
            }

            j = p_mol->constrBonds[i].atomI;
            atomR = p_mol->atoms[j].R;
            for (k = 0; k < 3; k++)
            {
                ATOMR(1, k) += p_mol->atoms[j].constraintForces[k];
                // velocity correction to be perpendicular to bond
            }
            j = p_mol->constrBonds[i].atomJ;
            atomR = p_mol->atoms[j].R;
            for (k = 0; k < 3; k++)
            {
                ATOMR(1, k) += p_mol->atoms[j].constraintForces[k];
                // velocity correction to be perpendicular to bond
            }
        }

#ifdef PARALLEL
#pragma omp atomic
#endif
        virial += local_virial;
        local_virial = 0.0;

        if (maxRelErr < system->epsShake)
        {
            (m > system->maxShakeIter) ? (system->maxShakeIter = m) : (true);
            return maxRelErr;
        }
    }

    (m > system->maxShakeIter) ? (system->maxShakeIter = m) : (true);
    return maxRelErr;
}

// initialize values of Xi
// Initialize Xi
int VerletNPTNose::InitXi()
{
    double Ekin;
    std::vector<Molecule>::iterator itmol;
    char aux[10];
    bool fromconfig = false;
    int m;
    double Ekinfinaltotal = 0.5 * (system->noDegrOfFred) * T_f; // not including Ekin of barostat (Lambda is not initialized now)

    if (system->extendedDOFs != nullptr)
    {
        fromconfig = true;
        sscanf(system->extDOFnames, "%s", aux);
        if (strstr(aux, "Xi") == NULL)
        {
            print_warning(1, "Extended DOFs from .config don't start by Xi, cannot use info from .config file\n");
            fromconfig = false;
        }
    }

    // calculate Ekin
    Ekin = 0.0;
    for (itmol = system->molecules.begin(); itmol != system->molecules.end(); itmol++)
    {
        Ekin += itmol->CalculateEkin(h, 1);
    }

    if (fromconfig)
    {
        for (m = 0; m < std::min(4, system->extendedDOFs->GetNumberOfRows()); m++)
        {
            Xi->operator()(m) = system->extendedDOFs->operator()(m, 0) * pow(h, (double)m) / fact(m);
        }
        if (Xi->operator()(2) == 0.0)
        {
            // calculate xi'' (force on xi)
            Xi->operator()(2) = h * h / 2.0 * (2 * Ekin - 2.0 * Ekinfinaltotal) / M_T;
        }
        // initialize past values of Xi (meanless but why not...)
        Xi->operator()(2) = Xi->operator()(0) - Xi->operator()(1) + Xi->operator()(2);
        // calculation of Xi(3) would be useless because the step starts with TimeShift()...
        // if different cfglevel, this can lead to incorrect initialization but nevermind...
    }
    else
    {
        // initialize past values of Xi (meanless but why not...)
        Xi->operator()(2) = Xi->operator()(0) - Xi->operator()(1) + h * h / 2.0 * (2 * Ekin - 2.0 * Ekinfinaltotal) / M_T;
        // calculation of Xi(3) would be useless because the step starts with TimeShift()...
    }
    return 0;
}

// initialize values of Lambda
int VerletNPTNose::InitLambda()
{
    double Ekininst;
    std::vector<Molecule>::iterator itmol;
    char aux1[10], aux2[10], aux3[10];
    bool fromconfig = false;
    int m, lambda_pos;
    double Pvir, Pkin;
    double V;

    // calculate Ekininst
    Ekininst = 0.0;
    for (itmol = system->molecules.begin(); itmol != system->molecules.end(); itmol++)
    {
        Ekininst += itmol->CalculateEkin(h, 1);
    }
    // calculate Pkin
    V = system->CalculateVolume();
    Pkin = system->CalculatePkin(h, Ekininst);

    if (system->extendedDOFs != nullptr)
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
        for (m = 0; m < std::min(4, system->extendedDOFs->GetNumberOfRows()); m++)
        {
            Lambda->operator()(m) = system->extendedDOFs->operator()(m, lambda_pos) * pow(h, (double)m) / fact(m);
        }
        if (Lambda->operator()(2) == 0.0) // Lambda(3) serves as a storage for Pvir from the previous run (to include also virconstr)
        {
            Pvir = system->CalculatePconf();
        }
        else
        {
            Pvir = Lambda->operator()(2);
            system->virconstr = (Pvir - system->CalculatePconf()) * (-3.0) * V; // to correct virial of constraints
        }
        // if different cfglevel, this can lead to incorrect initialization but nevermind...
        delete system->extendedDOFs;
        system->extendedDOFs = nullptr;
    }
    else
    {
        Pvir = system->CalculatePconf();
        // initialize past values of Xi (meanless but why not...)
        Lambda->operator()(2) = Lambda->operator()(0) - Lambda->operator()(1) + 3.0 * h * h / 2.0 * ((Pkin + Pvir - P_f) * V + 2.0 * Ekininst / system->noDegrOfFred) / M_P;
        // calculation of Xi(3) would be useless because the step starts with TimeShift()...
    }
    return 0;
}

// prepare Xi for .config printing
int VerletNPTNose::PrepareConfig()
{
    int i, j, k;
    int configLevel;
    double Ekin = 0.0;
    std::vector<Molecule>::iterator itmol;
    Atom *atom;

    (*Xi_old) = (*Xi);
    (*Lambda_old) = (*Lambda);

    for (i = 0; i < system->noMolecules; i++)
    {
        for (j = 0; j < system->molecules[i].noAtoms; j++)
        {
            atom = &(system->molecules[i].atoms[j]);
            for (k = 0; k < 3; k++)
            {
                atom->R->operator()(2, k) = h * h * 0.5 * atom->force[k] / atom->mass;
            }
        }
    }

    for (itmol = system->molecules.begin(); itmol != system->molecules.end(); itmol++)
    {
        Ekin += itmol->CalculateEkin(h, 1);
    }

    // calculate xi'' (force on xi)
    Xi->operator()(2) = 0.5 * h * h / tau_T / tau_T * (Ekin / (0.5 * system->noDegrOfFred * T_f) - 1.0);

    Lambda->operator()(2) = system->CalculatePconf();

    system->extendedDOFs = new Matrix(3, 2);
    for (i = 0; i < 3; i++)
    {
        system->extendedDOFs->operator()(i, 0) = Xi->operator()(i);
        system->extendedDOFs->operator()(i, 1) = Lambda->operator()(i);
    }

    configLevel = 202;
    strcpy(system->extDOFnames, "Xi Lambda");
    return configLevel;
}

// new initialize
int VerletNPTNose::Initialize(bool startWithShake)
{
    int i, j, k;
    Matrix *atomR;

    system->CalculateForces(); // returns potential energy (saved in SimulatedSystem itself)

    for (i = 0; i < system->noMolecules; i++)
    {
        for (j = 0; j < system->molecules[i].noAtoms; j++)
        {
            atomR = system->molecules[i].atoms[j].R;
            for (k = 0; k < 3; k++)
            {
                if (ATOMR(2, k) == 0.0)
                {
                    // calculate acceleration from forces
                    ATOMR(2, k) = 0.5 * h * h * system->molecules[i].atoms[j].force[k] / system->molecules[i].atoms[j].mass;
                }
                ATOMR(1, k) = ATOMR(1, k) * h;
                // initialize past values of Xi (meanless but why not...)
                ATOMR(2, k) = ATOMR(0, k) - ATOMR(1, k) + ATOMR(2, k);
            }
        }
        // giving up trying to get correct pressure for the first measurement ensured by Lambda(2) saving
    }

    return 0;
}

int VerletNPTNose::RestoreExtendedDOFs()
{
    (*Lambda) = (*Lambda_old);
    (*Xi) = (*Xi_old);
    return 0;
}