/*
 * An implementation of Gear integration method with Nose thermostat
 * for trial simulations of constrained dynamics
 * Author JJ, Date Sep 2021
 */

/*
 * In atoms R following values are stored:
 * R(0, :): x(i) y(i) z(i)             common for all integrators
 * R(1, :): h*vx(i) h*vy(i) h*vz(i)    common for all integrators (h* because of Gear)
 * R(2, :): h^2/2*ax(i)...             Gear specific
 * R(3, :): h^3/6*x'''(i)...           Gear specific
 * R(4, :): h^4/24*x^(iv)(i)...        Gear specific
 * ...
 *
 * R(0, :) and R(1, :) must be defined before integrator initialization
 *
 * In case of TRVP the last rows in R contain differences of previous positions:
 * R(predVel, :): h*vx(i)^P, ...             predicted velocity
 * R(predVel+1, :): x(i-1) - x(i-2), ...     difference of positions
 * R(predVel+2, :): x(i-2) - x(i-3), ...     dtto
 * ...
 */

/*
 * Nose–Hoover thermostat implemented according to differential equation in MACSIMUS manual p. 163
 * compliant with the TRVP paper (Kolafa, Lísal (2011))
 */

#include <iostream>
#include <cmath>
#include <cstring>

#include "Vector.hpp"
#include "Matrix.hpp"
#include "math_utils.hpp"
#include "Atom.hpp"
#include "Molecule.hpp"
#include "ConstraintBond.hpp"
#include "SimulatedSystem.hpp"
#include "GearNVTNose.hpp"
#include "general_utils.hpp"

#define ATOMR(i, j) ((*atomR)(i, j))

// initialization of this integrator (main part is done in AbstractGearIntegrator construction)
GearNVTNose::GearNVTNose(double step, SimulatedSystem *simsystem, int Rsize, char *integName, double Tfinal, double tauT, Vector trvpB, bool startWithShake)
    : AbstractGearIntegrator{step, simsystem, Rsize, integName, trvpB, startWithShake}
{
    // desired temperature
    Tf = Tfinal;
    // M_T = N_f * k_B * T_f * tauT^2
    // M_T stored instead of tau_T
    tau_T = system->noDegrOfFred * Tfinal * tauT * tauT;
    // temperature 'friction' factor (Xi')
    lambda_T = 1.0;
    // intaneous temperature set to final (just to have reasonable value)
    Tinst = Tf;

    // Extra degree of freedom – variable connected with heat bath (Xi)
    // allocation
    ExtDoFs = new Matrix(size, 1);
    GForExtDoFs = new Vector(1);
    // initialization of Xi
    InitXi();
    // energy of extDoF to system->Eextra
    system->Eextra = EnergyOfExtDoF();
}

// destructor
GearNVTNose::~GearNVTNose()
{
    if (ExtDoFs != nullptr)
    {
        delete ExtDoFs;
        ExtDoFs = nullptr;
    }
    if (GForExtDoFs != nullptr)
    {
        delete GForExtDoFs;
        GForExtDoFs = nullptr;
    }
}

// RHS of the differential equation
inline int GearNVTNose::CalculateRHS(Atom *at)
{
    int i;

    for (i = 0; i < 3; i++)
    {
        at->errorG->operator()(i) = 0.5 * (h * h * at->force[i] / at->mass - lambda_T * at->R->operator()(predVel, i)) - at->R->operator()(2, i);
        // G = h^2/2*(f/m - xi'*v) - h^2/2*a^P – Nose according to MACSIMUS manual p. 163
    }

    return 0;
}

// energy of extra degrees of freedom
inline double GearNVTNose::EnergyOfExtDoF()
{
    // Eextra = Ekinextra + Epotextra = 0.5 * M_T * xi'^2 + N_f * k_B * T_f * xi
    return 0.5 * tau_T * ExtDoFs->operator()(1, 0) * ExtDoFs->operator()(1, 0) / h / h + system->noDegrOfFred * Tf * ExtDoFs->operator()(0, 0);
}

// get temperature scaling ('friction') factor
inline double GearNVTNose::GetTempScaling()
{
    // simple but put here for consistence
    return ExtDoFs->operator()(predVel, 0);
}

// right-hand side for extra degrees of freedom
inline double GearNVTNose::RHSForExtDoF(int extDOFno)
{
    // correction for Xi: G = h^2/2*(1/tauT^2 * (Tkin/Tfinal - 1))  - h^2/2*a^P
    // because xi'' = 1/(tauT^2) * (Tkin/Tfinal - 1)
    // M_T = N_f * k_B * T_f * tauT^2
    // xi'' = 1 / M_T * (2 * Ekininst - N_f * k_B * Tfinal)
    double Ekininst = system->CalculateEkin(h);
    return h * h * 0.5 / tau_T * (2.0 * Ekininst - system->noDegrOfFred * Tf);
}

// Initialization of Xi
// If no values in .config file, calculate kinetic energy and start anew
// else use values from .config file (stored in system->extendedDOFs)
int GearNVTNose::InitXi()
{
    double Ekin;
    std::vector<Molecule>::iterator itmol;
    int m;
    char aux[10];
    bool fromconfig = false;

    // check if extended DOFs were loaded in SimulatedSystem::ReadConfig()
    if (system->extendedDOFs != nullptr)
    {
        // if they were, then set fromconfig to true and read extDOF names
        fromconfig = true;
        sscanf(system->extDOFnames, "%s", aux);
        // the first name must be 'Xi'...
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

    // if extDOFs present and start by 'Xi' then read Xi from them
    if (fromconfig)
    {
        // we can read only those values stored (in system->extendedDOFs)
        // and only those we have space for (shakeType = size if not TRVP and shakeType = predVel if TRVP)
        for (m = 0; m < std::min(std::min(size_Shake, system->extendedDOFs->GetNumberOfRows()), ExtDoFs->GetNumberOfRows()); m++)
        {
            // standard Nordsiek's scaling (as in Taylor expansion)
            ExtDoFs->operator()(m, 0) = system->extendedDOFs->operator()(m, 0) * pow(h, (double)m) / fact(m);
        }
        if (ExtDoFs->operator()(2, 0) == 0.0)
        {
            // initialize Xi(2) (should not be needed)
            ExtDoFs->operator()(2, 0) = 0.5 * h * h / tau_T * (2.0 * Ekin - system->noDegrOfFred * Tf);
        }
        // TRVP – stored past velocities
        if (predVel > 1)
        {
            for (m = predVel; m < std::min(size, system->extendedDOFs->GetNumberOfRows()); m++)
            {
                // we can read only those values stored (in system->extendedDOFs)
                // and only those we have space for (size = Rsize)
                ExtDoFs->operator()(m, 0) = system->extendedDOFs->operator()(m, 0) * h;
                // if not present in .config file...
                if (ExtDoFs->operator()(m, 0) == 0.0)
                {
                    ExtDoFs->operator()(m, 0) = ExtDoFs->operator()(m - 1, 0);
                }
            }
        }
        // if different normal and TRVP level, this can lead to incorrect initialization
        // the same integrator must be used when init continue requested
        delete system->extendedDOFs;
        system->extendedDOFs = nullptr;
    }
    else
    {
        // initialize Xi(2) using instanteneous value of kinetic energy
        ExtDoFs->operator()(2, 0) = 0.5 * h * h / tau_T * (2.0 * Ekin - system->noDegrOfFred * Tf);

        // initialize past TRVP values using xi''
        if (predVel > 1)
        {
            for (m = predVel + 1; m < size; m++)
            {
                ExtDoFs->operator()(m, 0) = -(m - predVel) * ExtDoFs->operator()(2, 0);
            }
        }
    }
    return 0;
}

// prepare Xi for .config printing
int GearNVTNose::PrepareConfig()
{
    int i;
    int configLevel = AbstractGearIntegrator::PrepareConfig();

    system->extendedDOFs = new Matrix(size, 1);
    for (i = 0; i < size; i++)
    {
        system->extendedDOFs->operator()(i, 0) = ExtDoFs->operator()(i, 0);
    }
    configLevel += 100;
    strcpy(system->extDOFnames, "Xi");
    return configLevel;
}

// perform N steps of integration
int GearNVTNose::Integrate(int Nsteps)
{
    #ifndef GEARTHERMO
    #define GEARTHERMO
    #endif
    #ifdef GEARBOXRESC
    #undef GEARBOXRESC
    #endif
    #ifndef GEAREXTDOFS
    #define GEAREXTDOFS
    #endif
    #include "GearIntegrateMacro.cpp"
}