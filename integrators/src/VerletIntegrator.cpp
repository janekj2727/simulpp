/*
 * An implementation of Verlet integration method for trial simulations of contrained dynamics
 * NVE ensemble – the simplest case
 * Author JJ, Date Jan 2021
 * Considerably modified Jul 2021
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

#include <iostream>
#include <cmath>

#include "Vector.hpp"
#include "Matrix.hpp"
#include "Atom.hpp"
#include "Molecule.hpp"
#include "ConstraintBond.hpp"
#include "SimulatedSystem.hpp"
#include "VerletIntegrator.hpp"

#define ATOMR(i, j) ((*atomR)(i, j))
// #define DEBUG

// initialization of this integrator
VerletIntegrator::VerletIntegrator(double step, SimulatedSystem *simsystem, bool startWithShake)
    : AbstractVerletIntegrator{step, simsystem, startWithShake}
{
    // constructor of AbstractVerletIntegrator should be sufficient
}

// initialization of this integrator + TRVP
VerletIntegrator::VerletIntegrator(double step, SimulatedSystem *simsystem, Vector trvpB, bool startWithShake)
    : AbstractVerletIntegrator{step, simsystem, trvpB, startWithShake}
{
    // constructor of AbstractVerletIntegrator should be sufficient
}

// empty destructor
VerletIntegrator::~VerletIntegrator()
{
}

// perform N  steps of integration
int VerletIntegrator::Integrate(int Nsteps)
{
    int i, l;
    Molecule *p_mol;
    double Ekininst;

    for (l = 0; l < Nsteps; l++)
    {
        Ekininst = 0.0; // zero kinetic energy cumulator
        // time shift r(i-1) = r(i) from previous and r(i) = r(i+1) from previous step
        TimeShift();
        // TRVP (if needed)
        // TRVP()

        system->CalculateForces(); // returns potential energy which is saved elsewhere (in SimulatedSystem.cpp itself)
#ifdef PARALLEL
#pragma omp parallel for num_threads(thread_count) reduction(+ \
                                                             : Ekininst) private(i, p_mol)
#endif
        for (i = 0; i < system->noMolecules; i++)
        {
            p_mol = &(system->molecules[i]);
            // one step of Verlet (new position calculation)
            VerletStep(p_mol);

            // SHAKE (shakeType = 0)
            Shake(p_mol, system->virconstr); // accumulation to systm->virconstr is controlled by omp atomic inside

            // calculate velocity and kinetic energy
            Ekininst += VelocityCalculation(p_mol);
        }

        if (system->boundaryCond > 0)
        {
            system->ApplyPeriodicBC(positionsAt);
        }
    }
    return 0;
}

int VerletIntegrator::VerletStep(Molecule *mol)
{
    int j, k;
    Matrix *atomR;
    for (j = 0; j < mol->noAtoms; j++)
    {
        atomR = mol->atoms[j].R;
        for (k = 0; k < 3; k++)
        {
            ATOMR(3, k) = 2 * ATOMR(0, k) - ATOMR(2, k) +
                          (h * h / mol->atoms[j].mass) * mol->atoms[j].force[k];
            // r(i+1) = 2*r(i) - r(i-1) + h^2/m * f(i)
            // ATOMR(1, k) = (ATOMR(3, k) - ATOMR(2, k)) * 0.5; // hv(i) = (r(i+1)-r(i-1))/2 (calculated later)
        }
    }
    return 0;
}
