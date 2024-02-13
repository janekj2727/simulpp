/*
 * An implementation of Verlet integration method (Velocity Verlet + Nose–Hoover thermostat)
 * "MTTK" – Trotter decomposition scheme, without TRVP (with TRVP the algorithm is completely different – defined in VerletNVTNoseTRVP.cpp)
 * for trial simulations of contrained dynamics
 * Author JJ, Date Oct 2021
 * Feb 2022 – still do not work ideally, energy conservation is good but fluctuations too large (esp. when compared to the DL POLY implementation)
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
 * Velocities are uniquely determined by the integration algorithm (and corresponds to VERLET=1).
 */

/*
 *  Nosé–Hoover thermostat implemented according to the DL POLY manual p. 83
 *  Checked against Frenkel&Smit p. 539
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
#include "VerletNVTNose.hpp"
#include "general_utils.hpp"

#define ATOMR(i, j) ((*atomR)(i, j))
// #define DEBUG

// initialization of this integrator (main part is done in AbstractVerletIntegrator construction)
VerletNVTNose::VerletNVTNose(double step, SimulatedSystem *simsystem, double Tfinal, double tauT, bool startWithShake)
    : AbstractVerletIntegrator{step, simsystem, startWithShake, true}
{
    T_f = Tfinal;
    tau_T = tauT;

    Initialize(startWithShake); // this constructor uses simpleinit version of AbstractVerletIntegrator not to overwrite ATOMR, shake is not performed

    // Extra degree of freedom – variable connected with heat bath
    Xi = new Vector(4);
    InitXi();
    system->Eextra = EnergyOfExtraDOF();
    Xi_old = new Vector(*Xi);
}

// destructor – delete thermostat variable Xi
VerletNVTNose::~VerletNVTNose()
{
    delete Xi;
    delete Xi_old;
}

// perform N  steps of integration
int VerletNVTNose::Integrate(int Nsteps)
{
    int i, l;
    Molecule *p_mol;
    double Ekininst = 0.0; // instanteneous kinetic energy
    double virial = 0.0;   // virial

    // calculate kinetic energy (initialization of Ekininst)
    // during integration the kinetic energy is updated
    for (i = 0; i < system->noMolecules; i++)
    {
        p_mol = &(system->molecules[i]);
        Ekininst += p_mol->CalculateEkin(h, 1);
    }

    for (l = 0; l < Nsteps; l++)
    {
        virial = 0.0;
        // time shift r(i-1) = r(i) from previous and r(i) = r(i+1) from previous step
        // special version of TimeShift for this integrator
        TimeShift();
        // time shift xi... as for r
        Xi->operator()(3) = Xi->operator()(2);
        Xi->operator()(2) = Xi->operator()(0);

        // thermostat h/2
        ThermostatHalfStep(Ekininst);

#ifdef PARALLEL
#pragma omp parallel for num_threads(thread_count) private(p_mol, i)
#endif
        for (i = 0; i < system->noMolecules; i++)
        {
            p_mol = &(system->molecules[i]);

            // VV1 (velocity update with old forces)
            VerletStep(p_mol);

            // RATTLE1 (different from SHAKE because of the different order in R)
            // bond lengths correction
            // virial protected inside by omp atomic
            Shake(p_mol, virial);
        }

        // boundary conditions
        if (system->boundaryCond > 0)
        {
            system->ApplyPeriodicBC(positionsAt);
        }

        // force calculation f(t+h)
        system->CalculateForces();

#ifdef PARALLEL
#pragma omp parallel for num_threads(thread_count) private(p_mol, i)
#endif
        for (i = 0; i < system->noMolecules; i++)
        {
            p_mol = &(system->molecules[i]);

            // VV2
            VerletStep2(p_mol);

            // RATTLE2
            // velocity perpendicular... – for nitrogen makes very little difference
            // virial protected by omp atomic inside
            Rattle(p_mol, virial); // changes energy but do not change Ekininst!!!
        }
        // recalculate Ekin
        Ekininst = system->CalculateEkin(h, 1);

        // thermostat h/2 (same as the 1st step)
        ThermostatHalfStep(Ekininst);

        // // boundary conditions // better to do before forces calculation (?)
        // if (system->boundaryCond > 0)
        // {
        //     system->ApplyPeriodicBC(positionsAt);
        // }
    }

    // calculate energy of extra DOF and pass it to system
    system->Eextra = EnergyOfExtraDOF();
    // add virial of constraint forces to the system virial of constr. forces
    system->virconstr = 1.0 / h / h * virial;
    // thermostat variables to system (to be printed in measurements if desired)
    system->Xi = Xi->operator()(0);
    system->Xivel = Xi->operator()(1) / h;

    return 0;
}

// calculate energy of the extra degree of freedom (xi)
double VerletNVTNose::EnergyOfExtraDOF() const
{
    // M_T = N_f * k_B * T_f * tauT^2
    // Eextra = Ekinextra + Epotextra = 0.5 * M_T * xi'^2 + N_f * k_B * T_f * xi
    return system->noDegrOfFred * T_f * (0.5 * tau_T * tau_T * Xi->operator()(1) * Xi->operator()(1) / h / h + Xi->operator()(0));
}

// half step of thermostat
int VerletNVTNose::ThermostatHalfStep(double &Ekin)
{
    double scaling;
    static double recTauT2 = 1.0 / (tau_T * tau_T);
    static double Ekinfinal = 0.5 * system->noDegrOfFred * T_f;

    /*
    // Try to compute the kinetic energy afresh - no change ( tested without SHAKE)
    int i;
    Ekin = 0.0;
    Molecule *p_mol;

    for (i = 0; i < system->noMolecules; i++)
    {
        p_mol = &(system->molecules[i]);
        Ekin += p_mol->CalculateEkin(h, 1);
    }
    */

    // in DL_POLY 1st integration of xi here (by h/8) (seems not to lead to better conservation)
    // Xi->operator()(0) += Xi->operator()(1) / 8.0; // DL_POLY version, originally nothing here

    // xi'(t+h/4) = xi'(t) + h/4 * (2Ekininst-2Ekinfinal)/M_T
    // (2Ekininst(new) - 2Ekinfinal)/M_T = 2Ekinfinal*(Ekininst/Ekinfinal - 1)/M_T = (Ekininst/Ekinfinal - 1)/tau_T^2
    // (M_T = N_f * k_B * T_f * tau_T^2 = 2Ekinfinal * tau_T^2)
    Xi->operator()(1) += h * h * 0.25 * (Ekin / Ekinfinal - 1.0) * recTauT2; // works (big fluctuations)
    // Xi->operator()(1) += h * h / 4.0 * (2.0 * Ekin - 2.0 * Ekinfinal) * recTauT2 / (2.0 * Ekinfinal); // no difference

    // iLxi(F&S): xi(t+h/2) = xi(t) + xi'(t+h/4)*h/2 (+ according to the fortran code, in contrast with − in equation E.2.9)
    // '+' also in the original paper by MTTK(1996)
    // originally here integration of Xi:
    Xi->operator()(0) += Xi->operator()(1) / 2.0;
    // in DL_POLY here only h/4
    // Xi->operator()(0) += Xi->operator()(1) / 4.0; // DL_POLY

    // v(t) = v(t) * exp(-xi'(t+h/4)*h/2)
    //     (E_kininst(new) = 0.5 * sum(m*(v*exp(-xi'(t+h/4)*h/2))^2)
    //      E_kininst(new) = E_kininst * exp(-2*xi'(t+h/4)*h/2))
    scaling = exp(-Xi->operator()(1) / 2.0); // works (big fluctuations)
    // scaling = 1.0 - Xi->operator()(1) / 2.0; // slightly different, but similar fluctuations
    system->RescaleVelocities(scaling);
    Ekin *= scaling * scaling;

    // xi'(t+h/2) = xi'(t+h/4) + h/4 * (2Ekininst(new)-2Ekinfinal)/M_T
    // (2Ekininst(new) - 2Ekinfinal)/M_T = 2Ekinfinal*(Ekininst/Ekinfinal - 1)/M_T = (Ekininst/Ekinfinal - 1)/tau_T^2
    Xi->operator()(1) += h * h * 0.25 * (Ekin / Ekinfinal - 1.0) * recTauT2;

    // in DL_POLY here the 3rd part of xi integration (by h/8)
    // Xi->operator()(0) += Xi->operator()(1) / 8.0; // DL_POLY version, originally nothing here

    return 0;
}

// VV1 – 1st part of velocity verlet
int VerletNVTNose::VerletStep(Molecule *p_mol)
{
    int j, k;
    Matrix *atomR;

    for (j = 0; j < p_mol->noAtoms; j++)
    {
        atomR = p_mol->atoms[j].R;
        for (k = 0; k < 3; k++)
        {
            // v(t+h/2) = v(t) + h/2 * f(t)/m
            ATOMR(1, k) += h * h * 0.5 * p_mol->atoms[j].force[k] / p_mol->atoms[j].mass;
            // r(t+h) = r(t) + h * v(t+h/2)
            ATOMR(0, k) += ATOMR(1, k);
        }
    }

    return 0;
}

// VV2 – 2nd part of velocity verlet
int VerletNVTNose::VerletStep2(Molecule *p_mol)
{
    int j, k;
    Matrix *atomR;

    for (j = 0; j < p_mol->noAtoms; j++)
    {
        atomR = p_mol->atoms[j].R;
        for (k = 0; k < 3; k++)
        {
            // v(t+h) = v(t+h/2) + h/2 * f(t+h)/m
            ATOMR(1, k) += h * h * 0.5 * p_mol->atoms[j].force[k] / p_mol->atoms[j].mass; // original version (works without SHAKE and with SHAKE with velocity update)
            // ATOMR(1, k) = ATOMR(0, k) - ATOMR(2, k) + h * h * 0.5 * p_mol->atoms[j].force[k] / p_mol->atoms[j].mass; // with SHAKE (now ATOM(1) updated during SHAKE)
        }
    }

    return 0;
}

// initialize values of Xi
// Initialize Xi
int VerletNVTNose::InitXi()
{
    double Ekin, force;
    std::vector<Molecule>::iterator itmol;
    char aux[10];
    bool fromconfig = false;
    int m;
    // int normallev = (predVel == 1)?(size):(size - predVel);

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
        for (m = 0; m < std::max(3, system->extendedDOFs->GetNumberOfRows()); m++)
        {
            Xi->operator()(m) = system->extendedDOFs->operator()(m, 0) * pow(h, (double)m) / fact(m);
        }
        if (Xi->operator()(2) == 0.0)
        {
            // calculate xi'' (force on xi)
            Xi->operator()(2) = 0.5 * h * h / tau_T / tau_T * (Ekin / (0.5 * system->noDegrOfFred * T_f) - 1.0);
        }
        // initialize past values of Xi (meanless but why not...)
        Xi->operator()(2) = Xi->operator()(0) - Xi->operator()(1) + Xi->operator()(2);
        // calculation of Xi(3) would be useless because the step starts with TimeShift()...
        // if different cfglevel, this can lead to incorrect initialization but nevermind...
        delete system->extendedDOFs;
        system->extendedDOFs = nullptr;
    }
    else
    {
        // calculate xi'' (force on xi)
        force = h * h / tau_T / tau_T * (Ekin / (0.5 * system->noDegrOfFred * T_f) - 1.0);

        // initialize past values of Xi (meanless but why not...)
        Xi->operator()(2) = Xi->operator()(0) - Xi->operator()(1) + 0.5 * force;
        // calculation of Xi(3) would be useless because the step starts with TimeShift()...
    }
    return 0;
}

// TimeShift version for this integrator (substantially different from the others)
int VerletNVTNose::TimeShift()
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
double VerletNVTNose::Shake(Molecule *mol, double &virial)
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
double VerletNVTNose::Rattle(Molecule *p_mol, double &virial)
{
    double maxRelErr = 0.0;
    static const int iterShake = 10; // maximum number of RATTLE iterations (currently fixed number)
    int m, i, j, k;
    Matrix *atomR;
    double relErr;
    double local_virial = 0.0;

    // what should be the value of eps and how it should be calculated???
    // for diatomics one iteration is enough, for larger molecules this could be a problem...

    for (m = 0; m < iterShake; m++)
    {
        maxRelErr = 0.0;

        // for each bond calculate contribution to constr forces
        for (i = 0; i < p_mol->noConstrBonds; i++)
        {
            if (relErr = p_mol->constrBonds[i].CalculateForces(local_virial, -2, system->epsShake), relErr > maxRelErr)
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

// prepare Xi for .config printing
int VerletNVTNose::PrepareConfig()
{
    int i, j, k;
    int configLevel;
    double Ekin;
    std::vector<Molecule>::iterator itmol;
    Atom *atom;

    (*Xi_old) = (*Xi);

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

    // calculate Ekin
    Ekin = 0.0;
    for (itmol = system->molecules.begin(); itmol != system->molecules.end(); itmol++)
    {
        Ekin += itmol->CalculateEkin(h, 1);
    }
    // calculate xi'' (force on xi)
    Xi->operator()(2) = 0.5 * h * h / tau_T / tau_T * (Ekin / (0.5 * system->noDegrOfFred * T_f) - 1.0);

    system->extendedDOFs = new Matrix(3, 1);
    for (i = 0; i < 3; i++)
    {
        system->extendedDOFs->operator()(i, 0) = Xi->operator()(i);
    }

    configLevel = 102;
    strcpy(system->extDOFnames, "Xi");
    return configLevel;
}

int VerletNVTNose::Initialize(bool startWithShake)
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
        // giving it up, not trying to get the correct pressure for the first measurement...
        // if (startWithShake)
        // {
        //     // SHAKE original for Verlet
        //     Shake(&(system->molecules[i]), system->virconstr);
        // }
        // for (j = 0; j < system->molecules[i].noAtoms; j++)
        // {
        //     atomR = system->molecules[i].atoms[j].R;
        //     for (k = 0; k < 3; k++)
        //     {
        //         ATOMR(3, k) = system->molecules[i].atoms[j].force[k];
        //     }
        // }
        // VelocityCalculation(&(system->molecules[i]));
        // if (startWithShake)
        // {
        //     system->virconstr *= 1.0;
        // }
    }

    return 0;
}

int VerletNVTNose::RestoreExtendedDOFs()
{
    (*Xi) = (*Xi_old);
    return 0;
}