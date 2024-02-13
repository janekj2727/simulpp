/*
 * An implementation of Verlet integration method (leapfrog + Berendsen thermostat)
 * for trial simulations of contrained dynamics
 * Author JJ, Date Jan 2021
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
 * In case of TRVP the next values of R are:
 * R(4, :): h*vx(i)^P, ...             predicted velocity
 * R(5, :): x(i-1) - x(i-2), ...       difference of positions
 * R(6, :): x(i-2) - x(i-3), ...       dtto
 * ...
 */

/*
 * MACSIMUS versions of VERLET velocity:
 * #define VERLET 0 The simplest, fastest, and least accurate O(h) version, v(t) = [r(t + h) − r(t)]/h.
 *   Essentially, velocity is shifted by h/2 from positions. Correct averaged energy of
 *   the harmonic oscillator. Good enough with friction (Berendsen) thermostat. (Added in
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
 * Berendsen thermostat implemented according to Frenkel&Smith (Understanding molecular simulation, p. 162)
 * Different versions of Berendsen thermostat can be found in MACSIMUS an DL_POLY (similar to this)
 */

#include <iostream>
#include <cmath>

#include "Vector.hpp"
#include "Matrix.hpp"
#include "Atom.hpp"
#include "Molecule.hpp"
#include "ConstraintBond.hpp"
#include "SimulatedSystem.hpp"
#include "VerletNVTBerendsen.hpp"

#define ATOMR(i, j) ((*atomR)(i, j))
// #define DEBUG

// initialization of this integrator
VerletNVTBerendsen::VerletNVTBerendsen(double step, SimulatedSystem *simsystem, double Tfinal, double tauT, bool startWithShake)
    : AbstractVerletIntegrator{step, simsystem, startWithShake}
{
    double instEkin;
    tau_T = tauT;
    T_f = Tfinal;
    instEkin = system->CalculateEkin(h);
    lambda = sqrt(1.0 + h / tau_T * (T_f * system->noDegrOfFred / (2.0 * instEkin) - 1.0));
}

// initialization of this integrator + TRVP
VerletNVTBerendsen::VerletNVTBerendsen(double step, SimulatedSystem *simsystem, Vector trvpB, double Tfinal, double tauT, bool startWithShake)
    : AbstractVerletIntegrator{step, simsystem, trvpB, startWithShake}
{
    double instEkin;
    tau_T = tauT;
    T_f = Tfinal;
    instEkin = system->CalculateEkin(h);
    lambda = sqrt(1.0 + h / tau_T * (T_f * system->noDegrOfFred / (2.0 * instEkin) - 1.0));
}

// empty destructor
VerletNVTBerendsen::~VerletNVTBerendsen()
{
}

// perform N  steps of integration
int VerletNVTBerendsen::Integrate(int Nsteps)
{
    int i, l;
    Molecule *p_mol;
    double Ekininst = 0.0; // instanteneous kinetic energy derived from vhp

    for (l = 0; l < Nsteps; l++)
    {
        Ekininst = 0.0;

        // time shift r(i-1) = r(i) from previous and r(i) = r(i+1) from previous step
        TimeShift();
        // TRVP (if needed)
        if (predVel != 1)
        {
            TRVP();
        }

        system->CalculateForces(); // returns potential energy which is saved elsewhere (in SimulatedSystem.cpp itself)

#ifdef PARALLEL
#pragma omp parallel for num_threads(thread_count) private(p_mol, i) reduction(+: Ekininst)
#endif
        for (i = 0; i < system->noMolecules; i++)
        {
            p_mol = &(system->molecules[i]);
            // one step of Verlet (new position calculation)
            VerletStep(p_mol, lambda);

            // SHAKE (shakeType = 0)
            Shake(p_mol, system->virconstr);

            Ekininst += VelocityCalculation(p_mol);
        }
        // calculate lambda from vhp temperature
        if (system->boundaryCond > 0)
        {
            system->ApplyPeriodicBC(positionsAt);
        }
        // calculate velocity scaling factor lambda from Ekin (from vhp – velocity in time t + h/2)
        lambda = sqrt(1.0 + h / tau_T * (T_f * system->noDegrOfFred / (2.0 * Ekininst) - 1.0));
    }
    return 0;
}

int VerletNVTBerendsen::VerletStep(Molecule *mol, double lambda)
{
    int j, k;
    Matrix *atomR;

    for (j = 0; j < mol->noAtoms; j++)
    {
        atomR = mol->atoms[j].R;
        for (k = 0; k < 3; k++)
        {
            ATOMR(3, k) = ATOMR(0, k) + lambda * (ATOMR(0, k) - ATOMR(2, k) + (h * h / mol->atoms[j].mass) * mol->atoms[j].force[k]);
            // r(i+1) = r(i) + lambda*(r(i) - r(i-1) + h^2/m * f(i))
        }
    }
    return 0;
}