/*
 * An implementation of Gear integration method with Nosé–Hoover (MTK) barostat (+ thermostat)
 * for simul++
 * Author JJ, Date Sep 2021
 *
 * Rewritten according to the template method pattern in Aug 2023
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
 * Nose–Hoover barostat implemented according to the universal integration scheme described in the documentation/Gear_formalism_for_MD...
 * Differential equations from the article by Janek+Kolafa(2022) (inspired by MACSIMUS and checked against DL_POLY)
 * Originally described by Martyna et al. (thus MTK)
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
#include "GearNPTNose.hpp"
#include "general_utils.hpp"

#define ATOMR(i, j) ((*atomR)(i, j))

// initialization of this integrator (main part is done in AbstractGearIntegrator construction) 
GearNPTNose::GearNPTNose(double step, SimulatedSystem *simsystem, int Rsize, char *integName, double Tfinal, double tauT, double Pfinal, double tauP, Vector trvpB, bool startWithShake)
    : AbstractGearIntegrator{step, simsystem, Rsize, integName, trvpB, startWithShake}
{
    double Ekininst = 0.0;
    Tf = Tfinal;
    tau_T = system->noDegrOfFred * Tfinal * tauT * tauT; // M_T instead of tauT;
    lambda_T = 1.0;
    Tinst = Tf;
    Pf = Pfinal;
    tau_P = (system->noDegrOfFred + 3) * Tfinal * tauP * tauP; // M_P instead of tauP;
    molecular_based = false;
    lambda_P = 1.0;

    // Extra degrees of freedom – variables connected with heat bath and size of the system Xi(0) and Lambda(1)
    ExtDoFs = new Matrix(size, 2);
    GForExtDoFs = new Vector(2);

    Ekininst = system->CalculateEkin(h, 1);
    InitXi(Ekininst);
    InitLambda(Ekininst);
    // (do it in Abstract... but problem with allocation and type of extDOF)

    // save old_scaling from Lambda(0)
    lambda_P = ExtDoFs->operator()(0, 1);
    system->Eextra = EnergyOfExtDoF();
}

// destructor
GearNPTNose::~GearNPTNose()
{
    delete ExtDoFs;
    ExtDoFs = nullptr;
    delete GForExtDoFs;
    GForExtDoFs = nullptr;
}

inline double GearNPTNose::GetScalingFactor(int order)
{
    double scaling;

    scaling = exp(ExtDoFs->operator()(0, 1) - lambda_P);
    lambda_P = ExtDoFs->operator()(0, 1);
    return scaling;
}

inline double GearNPTNose::GetShakeScaling()
{
    return exp(((*A) * (*ExtDoFs))(0, 1) - lambda_P);
}

// temperature scaling
inline double GearNPTNose::GetTempScaling()
{
    // constant friction factor (the same for all atoms) - h * (xi' + 3lambda'/N_f + lambda')
    return ExtDoFs->operator()(predVel, 0) + (1.0 + 3.0 / system->noDegrOfFred) * ExtDoFs->operator()(predVel, 1);
}

// RHS for ExtDoFs
inline double GearNPTNose::RHSForExtDoF(int extDOFno)
{
    double Ekin, V, Pkin, Pvir, Pinst;

    switch (extDOFno)
    {
    case 0:
        // correction for Xi: G = h^2/2*xi''  - h^2/2*xi''^P
        // xi'' = 1/M_T * (2*Ekin + M_P * lambda'^2 - (N_f + 1) * k_B * T_f)
        Ekin = system->CalculateEkin(h);
        return 0.5 / tau_T * (h * h * (2.0 * Ekin - (system->noDegrOfFred + 1) * Tf) + tau_P * ExtDoFs->operator()(predVel, 1) * ExtDoFs->operator()(predVel, 1));
    case 1:
        // pressure calculation
        Ekin = system->CalculateEkin(h);
        V = system->CalculateVolume();
        Pkin = system->CalculatePkin(h, Ekin);
        Pvir = system->CalculatePconf();
        Pinst = Pkin + Pvir;
        // correction for Lambda: G = h^2/2lambda'' - h^2/2*lambda''^P
        // lambda'' = 3/M_P * (V * (Pinst - P_f) + 2Ekin/N_f) - xi' * lambda'
        return 0.5 * (h * h * 3.0 / tau_P * (V * (Pinst - Pf) + 2.0 * Ekin / system->noDegrOfFred) - ExtDoFs->operator()(1, 0) * ExtDoFs->operator()(predVel, 1));
    }

    return 0.0;
}

// RHS for physical degrees of freedom
inline int GearNPTNose::CalculateRHS(Atom *at)
{
    int i;
    // normal correction (part 1 – before SHAKE)
    for (i = 0; i < 3; i++)
    {
        at->errorG->operator()(i) =
            0.5 * (h * h * at->force[i] / at->mass - lambda_T * at->R->operator()(predVel, i)) - at->R->operator()(2, i);
        // G = h^2/2*(f/m - (xi' + 3lambda'/N_f + lambda')*v) - h^2/2*a^P
        // h * (xi' + 3lambda'/N_f + lambda') stored in lambda_T (done in GetTempScaling())
    }
    return 0;
}

// Energy of extra degrees of freedom
inline double GearNPTNose::EnergyOfExtDoF()
{
    // Eextra = Ekinextra + Epotextra = 0.5 * M_T * xi'^2 + 0.5 * M_P * lambda'^2 + (N_f + 1) * k_B * T_f * xi + P_f * V

    return 0.5 / h / h * (tau_T * ExtDoFs->operator()(1, 0) * ExtDoFs->operator()(1, 0) + tau_P * ExtDoFs->operator()(1, 1) * ExtDoFs->operator()(1, 1)) + (system->noDegrOfFred + 1) * Tf * ExtDoFs->operator()(0, 0) + Pf * system->CalculateVolume();
}

// Initialization of Xi
int GearNPTNose::InitXi(double Ekininst)
{
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

    // if extDOFs present and start by 'Xi' then read Xi from them
    if (fromconfig)
    {
        // we can read only those values stored (in system->extendedDOFs)
        // and only those we have space for (size_Shake = size if not TRVP and size_Shake = predVel if TRVP)
        for (m = 0; m < std::min(size_Shake, system->extendedDOFs->GetNumberOfRows()); m++)
        {
            // standard Nordsiek's scaling (as in Taylor expansion)
            ExtDoFs->operator()(m, 0) = system->extendedDOFs->operator()(m, 0) * pow(h, (double)m) / fact(m);
        }
        // TRVP – stored past velocities
        if (predVel > 1) // means TRVP used
        {
            // we can read only those values stored (in system->extendedDOFs)
            // and only those we have space for (size = Rsize)
            for (m = predVel; m < std::min(size, system->extendedDOFs->GetNumberOfRows()); m++)
            {
                // scaling as for velocities (it IS a past velocity)
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
    }
    else
    {
        // initialize Xi(2) using instanteneous value of kinetic energy
        ExtDoFs->operator()(2, 0) = 1.0 / tau_T * (h * h * (Ekininst - 0.5 * (system->noDegrOfFred + 1) * Tf) + 0.5 * tau_P * ExtDoFs->operator()(predVel, 1) * ExtDoFs->operator()(predVel, 1));

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

// Initialization of Lambda
int GearNPTNose::InitLambda(double Ekininst)
{
    double V, Pkin, Pvir, Pinst;
    int m;
    bool fromconfig = false;
    char aux1[10], aux2[10], aux3[10];
    int lambda_pos = 0; // position of lambda in .config file

    // calculate instanteneous pressure
    V = system->CalculateVolume();
    Pkin = system->CalculatePkin(h, Ekininst);
    Pvir = system->CalculatePconf();
    Pinst = Pkin + Pvir;

    // check if extended DOFs were loaded in SimulatedSystem::ReadConfig()
    if (system->extendedDOFs != nullptr)
    {
        // if they were, then set fromconfig to true and read extDOF names
        fromconfig = true;
        m = sscanf(system->extDOFnames, "%s %s %s", aux1, aux2, aux3);
        // check the position of 'Lambda' ('Xi' should be first)
        if ((m > 1) && (strstr(aux2, "Lambda") != NULL))
        {
            lambda_pos = 1;
        }
        // lambda_pos should be 1, next lines are redundant
        else if ((m > 2) && (strstr(aux3, "Lambda") != NULL))
        {
            lambda_pos = 2;
        }
        // if there is no extDOF named 'Lambda', print warning and don't read Lambda from system->extendedDOFs
        else
        {
            print_warning(1, "Extended DOFs from .config don't define Lambda, cannot use info from .config file\n");
            fromconfig = false;
        }
    }

    // if Lambda is to be read from system->extendedDOFs
    if (fromconfig)
    {
        // we can read only those values stored (in system->extendedDOFs)
        // and only those we have space for (size_Shake = size if not TRVP and size_Shake = predVel if TRVP)
        for (m = 0; m < std::min(size_Shake, system->extendedDOFs->GetNumberOfRows()); m++)
        {
            // standard Nordsiek's scaling (as in Taylor expansion)
            ExtDoFs->operator()(m, 1) = system->extendedDOFs->operator()(m, lambda_pos) * pow(h, (double)m) / fact(m);
        }
        // TRVP – stored past velocities
        if (predVel > 1)
        {
            // we can read only those values stored (in system->extendedDOFs)
            // and only those we have space for (size = Rsize)
            for (m = predVel; m < std::min(size, system->extendedDOFs->GetNumberOfRows()); m++)
            {
                // scaling as for velocities (it IS a past velocity)
                ExtDoFs->operator()(m, 1) = system->extendedDOFs->operator()(m, lambda_pos) * h;
                // if not present in .config file...
                if (ExtDoFs->operator()(m, 1) == 0.0)
                {
                    ExtDoFs->operator()(m, 1) = ExtDoFs->operator()(m - 1, 1);
                }
            }
        }
        if (ExtDoFs->operator()(2, 1) == 0.0)
        {
            // initialize Lambda(2) (using Ekininst, Pinst, V and Xi(predVel) = ExtDoFs(predVel, 0))
            ExtDoFs->operator()(2, 1) = h * h * (1.5 / tau_P * (V * (Pinst - Pf) + 3.0 * Ekininst / system->noDegrOfFred)) - 0.5 * ExtDoFs->operator()(predVel, 0) * ExtDoFs->operator()(predVel, 1);
        }
        // if different normal and TRVP level, this can lead to incorrect initialization
        // the same integrator must be used when init continue requested
        delete system->extendedDOFs;
        system->extendedDOFs = nullptr;
    }
    else
    {
        // initialize Lambda(2) (using Ekininst, Pinst, V and Xi(predVel) = ExtDoFs(predVel, 0))
        ExtDoFs->operator()(2, 1) = h * h * (1.5 / tau_P * (V * (Pinst - Pf) + 3.0 * Ekininst / system->noDegrOfFred)) - 0.5 * ExtDoFs->operator()(predVel, 0) * ExtDoFs->operator()(predVel, 1);

        // initialize past TRVP values using lambda''
        if (predVel > 1)
        {
            for (m = predVel + 1; m < size; m++)
            {
                ExtDoFs->operator()(m, 1) = -(m - predVel) * ExtDoFs->operator()(2, 1);
            }
        }
    }
    return 0;
}

// prepare Xi for .config printing
int GearNPTNose::PrepareConfig()
{
    int i;
    int configLevel = AbstractGearIntegrator::PrepareConfig();

    system->extendedDOFs = new Matrix(size, 2);
    for (i = 0; i < size; i++)
    {
        system->extendedDOFs->operator()(i, 0) = ExtDoFs->operator()(i, 0);
        system->extendedDOFs->operator()(i, 1) = ExtDoFs->operator()(i, 1);
    }
    configLevel += 200;
    strcpy(system->extDOFnames, "Xi Lambda");
    return configLevel;
}

// perform N  steps of integration
int GearNPTNose::Integrate(int Nsteps)
{
    #ifndef GEARTHERMO
    #define GEARTHERMO
    #endif
    #ifndef GEARBOXRESC
    #define GEARBOXRESC
    #endif
    #ifndef GEAREXTDOFS
    #define GEAREXTDOFS
    #endif
    #include "GearIntegrateMacro.cpp"
}