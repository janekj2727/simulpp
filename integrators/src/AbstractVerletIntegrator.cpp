/*
 * An implementation of Verlet integration method for trial simulations of contrained dynamics
 * General Verlet method as a base for Verlet versions for particular ensembles
 * Author JJ, Date Jul 2021
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
 *
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
#include "math_utils.hpp"
#include "Atom.hpp"
#include "Molecule.hpp"
#include "ConstraintBond.hpp"
#include "SimulatedSystem.hpp"
#include "AbstractVerletIntegrator.hpp"

#define ATOMR(i, j) ((*atomR)(i, j))
// #define DEBUG

// initialization of this integrator
AbstractVerletIntegrator::AbstractVerletIntegrator(double step, SimulatedSystem *simsystem, bool startWithShake)
    : AbstractIntegrator{step, simsystem, 1}
{
    trvpCoeff = nullptr;

    Initialize(startWithShake);

    positionsAt.push_back(0);
    positionsAt.push_back(2);
    positionsAt.push_back(3);
}

// initialization of this integrator (simpleinit version – without initialization)
AbstractVerletIntegrator::AbstractVerletIntegrator(double step, SimulatedSystem *simsystem, bool startWithShake, bool simpleinit)
    : AbstractIntegrator{step, simsystem, 1}
{
    trvpCoeff = nullptr;

    if (!simpleinit)
    {
        Initialize(startWithShake);
    }

    positionsAt.push_back(0);
    positionsAt.push_back(2);
    positionsAt.push_back(3);
}

// initialization of this integrator + TRVP
AbstractVerletIntegrator::AbstractVerletIntegrator(double step, SimulatedSystem *simsystem, Vector trvpB, bool startWithShake)
    : AbstractIntegrator{step, simsystem, 4}
{
    trvpCoeff = new Vector(trvpB);

    Initialize(predVel, startWithShake);

    positionsAt.push_back(0);
    positionsAt.push_back(2);
    positionsAt.push_back(3);
}

// destructor deletes predictor matrix and corrector vector created during construction
AbstractVerletIntegrator::~AbstractVerletIntegrator()
{
    if (trvpCoeff != nullptr)
    {
        delete trvpCoeff;
        trvpCoeff = nullptr;
    }
}

// initialize values in R in atoms...
int AbstractVerletIntegrator::Initialize(bool startWithShake)
{
    int i, j, k;
    Matrix *atomR;
    double temp;

    // // system->instEpot = system->CalculateForces(); // returns potential energy
    // system->CalculateForces(); // returns potential energy (potential energy is saved in SimulatedSystem)

    // for (i = 0; i < system->noMolecules; i++)
    // {
    //     for (j = 0; j < system->molecules[i].noAtoms; j++)
    //     {
    //         atomR = system->molecules[i].atoms[j].R;
    //         for (k = 0; k < 3; k++)
    //         {
    //             ATOMR(1, k) = ATOMR(1, k) * h;                                                                                                       // v(i) -> h*v(i)
    //             ATOMR(2, k) = ATOMR(0, k) - ATOMR(1, k) + h * h * 0.5 * system->molecules[i].atoms[j].force[k] / system->molecules[i].atoms[j].mass; // r(i-1) = r(i) - h*v(i) + h^2/2*f(i)/m
    //             ATOMR(3, k) = ATOMR(0, k) + ATOMR(1, k) + h * h * 0.5 * system->molecules[i].atoms[j].force[k] / system->molecules[i].atoms[j].mass; // r(i+1) = r(i) + h*v(i) + h^2/2*f(i)/m
    //         }
    //     }
    //     if (startWithShake)
    //     {
    //         // SHAKE original for Verlet
    //         Shake(&(system->molecules[i]), system->virconstr);
    //     }
    // }
    // return 0;

    system->CalculateForces(); // returns potential energy (saved in SimulatedSystem itself)

    for (i = 0; i < system->noMolecules; i++)
    {
        for (j = 0; j < system->molecules[i].noAtoms; j++)
        {
            atomR = system->molecules[i].atoms[j].R;
            for (k = 0; k < 3; k++)
            {
                ATOMR(1, k) = ATOMR(1, k) * h;
                if (ATOMR(2, k) == 0.0) // if empty, use force
                {
                    ATOMR(2, k) = system->molecules[i].atoms[j].force[k] / system->molecules[i].atoms[j].mass;
                }
                ATOMR(2, k) = ATOMR(2, k) * h * h * 0.5;
                ATOMR(3, k) = ATOMR(0, k) + ATOMR(1, k) + ATOMR(2, k); // r(i+1) = r(i) + h*v(i) + h^2/2*f(i)/m (for NVE, more complicated for other ensembles)
                ATOMR(2, k) = ATOMR(0, k) - ATOMR(1, k);               // r(i-1) = r(i) - h*v(i-1/2)  (for NVE, more complicated for other ensembles)
                                                                       // if not TRVP (predVel==1) than not actually needed

                temp = ATOMR(3, k);
                ATOMR(3, k) = ATOMR(0, k) + ATOMR(1, k) + h * h * system->molecules[i].atoms[j].force[k] / system->molecules[i].atoms[j].mass;
                system->molecules[i].atoms[j].force[k] = temp;
            }
        }
        if (startWithShake)
        {
            // SHAKE original for Verlet
            Shake(&(system->molecules[i]), system->virconstr);
        }
        for (j = 0; j < system->molecules[i].noAtoms; j++)
        {
            atomR = system->molecules[i].atoms[j].R;
            for (k = 0; k < 3; k++)
            {
                ATOMR(3, k) = system->molecules[i].atoms[j].force[k];
            }
        }
        VelocityCalculation(&(system->molecules[i]));
        if (startWithShake)
        {
            system->virconstr *= 1.0;
        }
    }

    return 0;
}

// initialize values in R in atoms... + TRVP differences
int AbstractVerletIntegrator::Initialize(int predVel, bool startWithShake)
{
    int i, j, k, m;
    Matrix *atomR;
    double temp;

    // system->instEpot =
    system->CalculateForces(); // returns potential energy (saved in SimulatedSystem itself)

    for (i = 0; i < system->noMolecules; i++)
    {
        for (j = 0; j < system->molecules[i].noAtoms; j++)
        {
            atomR = system->molecules[i].atoms[j].R;
            for (k = 0; k < 3; k++)
            {
                ATOMR(1, k) = ATOMR(1, k) * h;
                if (ATOMR(2, k) == 0.0) // if empty, use force
                {
                    ATOMR(2, k) = system->molecules[i].atoms[j].force[k] / system->molecules[i].atoms[j].mass;
                }
                ATOMR(2, k) = ATOMR(2, k) * h * h * 0.5;
                ATOMR(3, k) = ATOMR(0, k) + ATOMR(1, k) + ATOMR(2, k); // r(i+1) = r(i) + h*v(i-1/2) + h^2/2*f(i)/m (for NVE, more complicated for other ensembles)
                if ((predVel == 1) || (ATOMR(4, k) == 0.0))
                {
                    ATOMR(2, k) = ATOMR(0, k) - ATOMR(1, k) + ATOMR(2, k); // r(i-1) = r(i) - h*v(i) + h^2/2*f(i)/m (for NVE, more complicated for other ensembles)
                    // if not TRVP (predVel==1) than not actually needed
                }
                else
                {
                    ATOMR(4, k) = ATOMR(4, k) * h;
                    ATOMR(2, k) = ATOMR(0, k) - ATOMR(4, k); // if TRVP, than this value needed...

                    // TRVP – stored past velocities
                    for (m = predVel + 1; m < atomR->GetNumberOfRows(); m++)
                    {
                        if (ATOMR(m, k) != 0.0)
                        {
                            ATOMR(m, k) *= h;
                        }
                        else
                        {
                            ATOMR(m, k) = ATOMR(1, k) - (m - predVel + 0.5) * h * h * 0.5 * system->molecules[i].atoms[j].force[k] / system->molecules[i].atoms[j].mass;
                        }
                    }
                }
                temp = ATOMR(3, k);
                ATOMR(3, k) = ATOMR(0, k) + ATOMR(1, k) + h * h * system->molecules[i].atoms[j].force[k] / system->molecules[i].atoms[j].mass;
                system->molecules[i].atoms[j].force[k] = temp;
            }
        }
        if (startWithShake)
        {
            // SHAKE original for Verlet
            Shake(&(system->molecules[i]), system->virconstr);
        }
        for (j = 0; j < system->molecules[i].noAtoms; j++)
        {
            atomR = system->molecules[i].atoms[j].R;
            for (k = 0; k < 3; k++)
            {
                ATOMR(3, k) = system->molecules[i].atoms[j].force[k];
            }
        }
        VelocityCalculation(&(system->molecules[i]));
        // if (startWithShake)
        // {
        //     system->virconstr *= 1.0;
        // }
    }

    return 0;
}

// prepare R for .config printing
int AbstractVerletIntegrator::PrepareConfig()
{
    int i, j, k;
    Matrix *atomR;
    int configLevel = 2;
    double temp = 0.0;

    for (i = 0; i < system->noMolecules; i++)
    {
        for (j = 0; j < system->molecules[i].noAtoms; j++)
        {
            atomR = system->molecules[i].atoms[j].R;
            for (k = 0; k < 3; k++)
            {
                if (predVel != 1)
                {
                    temp = ATOMR(4, k);
                    ATOMR(4, k) = ATOMR(0, k) - ATOMR(2, k);
                }
                ATOMR(1, k) = ATOMR(0, k) - ATOMR(2, k);
                ATOMR(2, k) = ATOMR(3, k) - ATOMR(0, k) - ATOMR(1, k); // to enable precise continue in whatever ensemble
                if (predVel != 1)
                {
                    ATOMR(3, k) = temp;
                }
            }
        }
    }
    if (predVel != 1)
    {
        configLevel += 10 * (trvpCoeff->GetSize() - 1);
    }

    return configLevel;
}

// SHAKE for 1 molecule
double AbstractVerletIntegrator::Shake(Molecule *mol, double &virial)
{
    double maxRelErr = 0.0;
    static const int iterShake = 100; // maximum number of SHAKE iterations (currently fixed number)
    int m, i, j, k;
    Matrix *atomR;
    double local_virial = 0.0;
    double relErr = 0.0;
    double omega = (mol->noConstrBonds > 1) ? (system->omegaShake) : (1.0); // overrelaxation parameter

    for (m = 0; m < iterShake; m++)
    {
        maxRelErr = 0.0;

        // for each bond calculate contribution to constr forces
        for (i = 0; i < mol->noConstrBonds; i++)
        {
            if (relErr = mol->constrBonds[i].CalculateForces(local_virial, 0, system->epsShake, 1.0, omega), relErr > maxRelErr)
            {
                maxRelErr = relErr;
            }
            j = mol->constrBonds[i].atomI;
            atomR = mol->atoms[j].R;
            for (k = 0; k < 3; k++)
            {
                ATOMR(3, k) += mol->atoms[j].constraintForces[k];
                // r(i+1)^c = r(i+1)^Verlet + h^2/m*f_c...
                mol->atoms[j].constraintForces[k] = 0.0;
            }
            j = mol->constrBonds[i].atomJ;
            atomR = mol->atoms[j].R;
            for (k = 0; k < 3; k++)
            {
                ATOMR(3, k) += mol->atoms[j].constraintForces[k];
                // r(i+1)^c = r(i+1)^Verlet + h^2/m*f_c...
                mol->atoms[j].constraintForces[k] = 0.0;
            }
        }
#ifdef PARALLEL
#pragma omp atomic
#endif
        virial += 1.0 / h / h * local_virial;
        local_virial = 0.0;
        // omega = (omega - 1.0) / 2.0 + 1.0;

        if (maxRelErr < system->epsShake)
        {
            (m > system->maxShakeIter) ? (system->maxShakeIter = m) : (true);
            return maxRelErr;
        }
    }

    (m > system->maxShakeIter) ? (system->maxShakeIter = m) : (true);
    return maxRelErr;
}

// SHAKE for 1 molecule with predicted rescaling
double AbstractVerletIntegrator::Shake(Molecule *mol, double &virial, double rescaling)
{
    double maxRelErr = 0.0, relErr = 0.0;
    static const int iterShake = 100; // maximum number of SHAKE iterations (currently fixed number)
    int m, i, j, k;
    Matrix *atomR;
    double local_virial = 0.0;
    double omega = (mol->noConstrBonds > 1) ? (system->omegaShake) : (1.0); // overrelaxation parameter

    for (m = 0; m < iterShake; m++)
    {
        maxRelErr = 0.0;
        // for each bond calculate contribution to constr forces
        for (i = 0; i < mol->noConstrBonds; i++)
        {
            if (relErr = mol->constrBonds[i].CalculateForces(local_virial, -3, system->epsShake, rescaling, omega), relErr > maxRelErr)
            {
                maxRelErr = relErr;
            }
            j = mol->constrBonds[i].atomI;
            atomR = mol->atoms[j].R;
            for (k = 0; k < 3; k++)
            {
                ATOMR(3, k) += mol->atoms[j].constraintForces[k];
                // r(i+1)^c = r(i+1)^Verlet + h^2/m*f_c...
            }
            j = mol->constrBonds[i].atomJ;
            atomR = mol->atoms[j].R;
            for (k = 0; k < 3; k++)
            {
                ATOMR(3, k) += mol->atoms[j].constraintForces[k];
                // r(i+1)^c = r(i+1)^Verlet + h^2/m*f_c...
            }
        }

#ifdef PARALLEL
#pragma omp atomic
#endif
        virial += 1.0 / h / h * local_virial;
        local_virial = 0.0;
        // omega = (omega - 1.0) / 2.0 + 1.0;

        if (maxRelErr < system->epsShake)
        {
            (m > system->maxShakeIter) ? (system->maxShakeIter = m) : (true);
            return maxRelErr;
        }
    }
    (m > system->maxShakeIter) ? (system->maxShakeIter = m) : (true);
    return maxRelErr;
}

int AbstractVerletIntegrator::TimeShift()
{
#ifdef PARALLEL // seems not to help (for 512 N2 molekules and trvp3)
#pragma omp parallel num_threads(thread_count)
    {
#endif
        int i, j, k, l;
        Matrix *atomR;
        if (predVel == 1)
        {
#ifdef PARALLEL
#pragma omp for
#endif
            for (i = 0; i < system->noMolecules; i++)
            {
                for (j = 0; j < system->molecules[i].noAtoms; j++)
                {
                    atomR = system->molecules[i].atoms[j].R;

                    for (k = 0; k < 3; k++)
                    {
                        ATOMR(2, k) = ATOMR(0, k); // r(i-1) = r(i) from previous
                        ATOMR(0, k) = ATOMR(3, k); // r(i) = r(i+1) from previous
                    }
                }
            }
        }
        else
        {
#ifdef PARALLEL
#pragma omp for
#endif
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
                        ATOMR(predVel + 1, k) = ATOMR(0, k) - ATOMR(2, k);
                        // normal Verlet TimeShift
                        ATOMR(2, k) = ATOMR(0, k); // r(i-1) = r(i) from previous
                        ATOMR(0, k) = ATOMR(3, k); // r(i) = r(i+1) from previous
                    }
                }
            }
        }
#ifdef PARALLEL
    }
#endif

    return 0;
}

// Calculate velocity according to VERLET #define and return kinetic energy
double AbstractVerletIntegrator::VelocityCalculation(Molecule *mol)
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
            vhm = ATOMR(0, k) - ATOMR(2, k); // hv(t) = [r(t) − r(t - h)];
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

int AbstractVerletIntegrator::VerletStep(Molecule *mol)
{
    // should be defined in each version of Verlet separately
    // this only to ensure existence
    return 0;
}

int AbstractVerletIntegrator::TRVP()
{
    // implementation of time-reversible velocity prediction (Kolafa+Lisal:JCTC 2011) for Verlet
    // version which saves differences (h*v) is used with coefficients B (in original paper)
    // R(0) - R(2) = r(t) - r(t-h), R(5) = r(t-h) - r(t-2h), R(6) = r(t-2h) - r(t-3h), ...

#ifdef PARALLEL // seems not to help (for 512 N2 molecules and trvp3)
#pragma omp parallel num_threads(thread_count)
#endif
    {
        int i, j, k, l;
        Matrix *atomR;
#ifdef PARALLEL
#pragma omp for
#endif
        for (i = 0; i < system->noMolecules; i++)
        {
            for (j = 0; j < system->molecules[i].noAtoms; j++)
            {
                atomR = system->molecules[i].atoms[j].R;
                for (k = 0; k < 3; k++)
                {
                    ATOMR(predVel, k) = trvpCoeff->operator()(0) * (ATOMR(0, k) - ATOMR(2, k));
                    for (l = 1; l < trvpCoeff->GetSize(); l++)
                    {
                        ATOMR(predVel, k) += trvpCoeff->operator()(l) * ATOMR(predVel + l, k);
                    }
                }
            }
        }
    }
    return 0;
}

// energy of extra DOF – default return 0.0
double AbstractVerletIntegrator::EnergyOfExtraDOF() const
{
    return 0.0;
}