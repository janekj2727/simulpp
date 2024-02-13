/*
 * An implementation of Verlet integration method (Verlet + TRVP + Nose–Hoover thermostat)
 * for trial simulations of contrained dynamics
 * Author JJ, Date Oct 2021
 */

/*
 * In atoms R following values are stored:
 * R(0, :): x(i) y(i) z(i)             common for all integrators
 * R(1, :): h*vx(i) h*vy(i) h*vz(i)    common for all integrators (h* because of Gear)
 * R(2, :): x(i-1) y(i-1) z(i-1)       Verlet specific
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
 * Nose thermostat implemented according to MACSIMUS manual (p 163(equations) and 302(TRVP))
 * and TRVP article (Kolafa+Lisal 2011)
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
#include "VerletNVTNoseTRVP.hpp"
#include "general_utils.hpp"

#define ATOMR(i, j) ((*atomR)(i, j))
// #define DEBUG

// initialization of this integrator
// version without TRVP (not needed for this integrator)
VerletNVTNoseTRVP::VerletNVTNoseTRVP(double step, SimulatedSystem *simsystem, double Tfinal, double tauT, bool startWithShake)
    : AbstractVerletIntegrator{step, simsystem, startWithShake}
{
    tau_T = tauT;
    T_f = Tfinal;
    Xi = new Vector(4);
    Xi_old = new Vector(*Xi);
}

// initialization of this integrator + TRVP
VerletNVTNoseTRVP::VerletNVTNoseTRVP(double step, SimulatedSystem *simsystem, Vector trvpB, double Tfinal, double tauT, bool startWithShake)
    : AbstractVerletIntegrator{step, simsystem, trvpB, startWithShake}
{
    tau_T = tauT;
    T_f = Tfinal;
    Xi = new Vector(simsystem->GetRsize());
    InitXi(predVel);
    system->Eextra = EnergyOfExtraDOF();
    system->Xi = Xi->operator()(0);
    system->Xivel = Xi->operator()(1) / h;
    Xi_old = new Vector(*Xi);
}

// empty destructor
VerletNVTNoseTRVP::~VerletNVTNoseTRVP()
{
    delete Xi;
    delete Xi_old;
}

// perform N  steps of integration
int VerletNVTNoseTRVP::Integrate(int Nsteps)
{
    int i, l;
    Molecule *p_mol;
    double hvXiP = 0.0;    // predicted value of xi(t)
    double Ekininst = 0.0; // instanteneous kinetic energy derived from vhp
    static double recTauT2 = 1.0 / (tau_T * tau_T);
    static double Ekinfinal = 0.5 * system->noDegrOfFred * T_f;

    for (l = 0; l < Nsteps; l++)
    {
        Ekininst = 0.0;

        // time shift r(i-1) = r(i) from previous and r(i) = r(i+1) from previous step
        TimeShift();

        // time shift for Xi
        // TRVP values (differences h*xi')
        for (i = predVel + trvpCoeff->GetSize() - 2; i > predVel; i--)
        {
            Xi->operator()(i + 1) = Xi->operator()(i);
        }
        Xi->operator()(predVel + 1) = Xi->operator()(0) - Xi->operator()(2);
        // normal TimeShift for xi
        Xi->operator()(2) = Xi->operator()(0); // r(i-1) = r(i) from previous
        Xi->operator()(0) = Xi->operator()(3); // r(i) = r(i+1) from previous
        system->CalculateForces();

        // TRVP for atomic velocities
        TRVP();
        // TRVP for Xi
        Xi->operator()(predVel) = trvpCoeff->operator()(0) * (Xi->operator()(0) - Xi->operator()(2));
        for (i = 1; i < trvpCoeff->GetSize(); i++)
        {
            Xi->operator()(predVel) += trvpCoeff->operator()(i) * Xi->operator()(predVel + i);
        }
        hvXiP = Xi->operator()(predVel);

#ifdef PARALLEL
#pragma omp parallel for num_threads(thread_count) private(p_mol, i) reduction(+ \
                                                                               : Ekininst)
#endif
        for (i = 0; i < system->noMolecules; i++)
        {
            p_mol = &(system->molecules[i]);
            // one step of Verlet (new position calculation), uses predicted value of h*xi'
            VerletStep(p_mol, hvXiP);

            // SHAKE (shakeType = 0)
            Shake(p_mol, system->virconstr); // accumulation to system->virconstr is controlled by omp atomic to prevent collisions

            // new differences r(t+h) - r(t), kinetic temperature calculation (according to VERLET)
            Ekininst += VelocityCalculation(p_mol);
        }

        if (system->boundaryCond > 0)
        {
            system->ApplyPeriodicBC(positionsAt);
        }

        // Verlet step for Xi
        Xi->operator()(3) = 2.0 * Xi->operator()(0) - Xi->operator()(2) +
                            h * h * recTauT2 * (Ekininst / Ekinfinal - 1.0);
    }
    // pass quantities to system for measuring...
    system->Eextra = EnergyOfExtraDOF(); // calculates also Xivel...
    system->Xi = Xi->operator()(0);
    system->Xivel = Xi->operator()(1) / h;

    return 0;
}

int VerletNVTNoseTRVP::VerletStep(Molecule *mol, double hvXiPred)
{
    int j, k;
    Matrix *atomR;

    for (j = 0; j < mol->noAtoms; j++)
    {
        atomR = mol->atoms[j].R;
        for (k = 0; k < 3; k++)
        {
            ATOMR(3, k) = 2.0 * ATOMR(0, k) - ATOMR(2, k) +
                          h * h * mol->atoms[j].force[k] / mol->atoms[j].mass -
                          ATOMR(predVel, k) * hvXiPred;
            // r(i+1) = 2r(i) - r(i-1) + h^2 * (f(i)/m - v(i)*xi'(i))
            // both velocities (v(i) and xi'(i)) use predicted values
        }
    }
    return 0;
}

// Energy of the extra degree of freedom
double VerletNVTNoseTRVP::EnergyOfExtraDOF() const
{
    double vhp, vhm;
    // Calculate Xi velocity according to VERLET definition
    vhp = Xi->operator()(3) - Xi->operator()(0); // hv(t) = [r(t + h) − r(t)]
    vhm = Xi->operator()(0) - Xi->operator()(2); // hv(t) = [r(t) − r(t - h)];

#if (VERLET == 0) // hv(t) = [r(t + h) − r(t)]
    Xi->operator()(1) = vhp;
#elif (VERLET == 1) // hv(i) = (r(i+1)-r(i-1))/2 - equivalent to k3m2e, velocity Verlet
    Xi->operator()(1) = (vhp + vhm) * 0.5;
#elif (VERLET == 2) // (hv(t))^2 = [r(t)−r(t−h)]·[r(t+h)−r(t)]
    Xi->operator()(1) = sqrt(fabs(vhm * vhp)) * ((vhp > 0) - (vhp < 0));
#else               // VERLET == 3  energy - averages of both shifted values v(t) = [r(t + h) − r(t)]/h and v(t) = [r(t) − r(t − h)]/h.
    Xi->operator()(1) = sqrt((vhm * vhm + vhp * vhp) * 0.5) * ((vhp > 0) - (vhp < 0));
#endif

    // M_T = N_f * k_B * T_f * tauT^2
    // Eextra = Ekinextra + Epotextra = 0.5 * M_T * xi'^2 + N_f * k_B * T_f * xi
    return system->noDegrOfFred * T_f * (0.5 * tau_T * tau_T * Xi->operator()(1) * Xi->operator()(1) / h / h + Xi->operator()(0));

    // DEBUGGING
    // return Xi->operator()(1);
}

// Initialize Xi
int VerletNVTNoseTRVP::InitXi(int predVel)
{
    double Ekin, force;
    std::vector<Molecule>::iterator itmol;
    int m;
    char aux[10];
    bool fromconfig = false;

    // calculate kinetic energy
    Ekin = 0.0;
    for (itmol = system->molecules.begin(); itmol != system->molecules.end(); itmol++)
    {
        Ekin += itmol->CalculateEkin(h, 1);
    }

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
            force = h * h / tau_T / tau_T * (Ekin / (0.5 * system->noDegrOfFred * T_f) - 1.0);
            Xi->operator()(2) = force;
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
                Xi->operator()(m) = Xi->operator()(1) - (m - predVel + 0.5) * 0.5 * h * h / tau_T / tau_T * (Ekin / (0.5 * system->noDegrOfFred * T_f) - 1.0);
            }
        }
        // }
    }
    else
    {
        // calculate xi'' (force on xi)
        force = h * h / tau_T / tau_T * (Ekin / (0.5 * system->noDegrOfFred * T_f) - 1.0);

        // initialize past values of Xi
        Xi->operator()(2) = Xi->operator()(0) - Xi->operator()(1) + 0.5 * force;
        Xi->operator()(3) = Xi->operator()(0) + Xi->operator()(1) + 0.5 * force;
        // TRVP – stored past velocities
        for (m = predVel + 1; m < Xi->GetSize(); m++)
        {
            if (Xi->operator()(m) != 0.0)
            {
                Xi->operator()(m) *= h;
            }
            else
            {
                Xi->operator()(m) = Xi->operator()(1) - (m - predVel + 0.5) * 0.5 * force;
            }
        }
    }

    return 0;
}

// Prepare Xi for .config printing
int VerletNVTNoseTRVP::PrepareConfig()
{
    int i;
    int configLevel = AbstractVerletIntegrator::PrepareConfig();
    double temp = 0.0;

    (*Xi_old) = (*Xi);

    temp = Xi->operator()(4);
    Xi->operator()(4) = Xi->operator()(0) - Xi->operator()(2);

    Xi->operator()(1) = Xi->operator()(0) - Xi->operator()(2);
    Xi->operator()(2) = Xi->operator()(3) - Xi->operator()(0) - Xi->operator()(1); // to enable precise continue in whatever ensemble

    Xi->operator()(3) = temp;

    system->extendedDOFs = new Matrix(Xi->GetSize(), 1);
    for (i = 0; i < Xi->GetSize(); i++)
    {
        system->extendedDOFs->operator()(i, 0) = Xi->operator()(i);
    }

    configLevel += 100;
    strcpy(system->extDOFnames, "Xi");

    return configLevel;
}

int VerletNVTNoseTRVP::RestoreExtendedDOFs()
{
    (*Xi) = (*Xi_old);
    return 0;
}
