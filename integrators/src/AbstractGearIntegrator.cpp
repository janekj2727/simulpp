/*
 * An implementation of Gear integration method for trial simulations of constrained dynamics
 * General Gear methods – particular versions for different ensembles should be derived from this
 * Author JJ, Date Jul 2021
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
 * R(predVel+1, :): x(i) - x(i-1), ...       difference of positions
 * R(predVel+2, :): x(i-1) - x(i-2), ...     dtto
 * ...
 */

#include <iostream>
#include <cmath>
#include <cstring>
#include <string>

#include "Vector.hpp"
#include "Matrix.hpp"
#include "math_utils.hpp"
#include "general_utils.hpp"
#include "Atom.hpp"
#include "Molecule.hpp"
#include "ConstraintBond.hpp"
#include "SimulatedSystem.hpp"
#include "AbstractGearIntegrator.hpp"

#define ATOMR(i, j) ((*atomR)(i, j))

// initialization of Gear integratorsď
AbstractGearIntegrator::AbstractGearIntegrator(double step, SimulatedSystem *simsystem, int Rsize, char *integName, Vector trvpB, bool startWithShake)
    : AbstractIntegrator{step, simsystem, (trvpB.GetSize() == 0) ? (1) : (Rsize - trvpB.GetSize())}, size(Rsize)
{
    int i, j, m;
    double X;
    char auxc;
    std::string integrator(integName);

    // if TRVP then trvpB is not nullptr and contains TRVP coeffs
    if (trvpB.GetSize() == 0)
    {
        // without TRVP, shakeType is the size of the method and predVel = 1
        size_Shake = Rsize;
    }
    else
    {
        // with TRVP, shakeType is the size of the method without TRVP part
        // predVel is the next term, so its coordinate is the same as shakeType
        size_Shake = Rsize - trvpB.GetSize(); // = predVel
    }

    // predictor and corrector construction from integratorName
    // method size was checked in Simulation::ReadControl and is equal to size_Shake
    // method order (number after m)
    auxc = integName[3];
    // convert number in auxc (ASCII) to integer (method order)
    m = (int)auxc - '0'; // convert char to int
    // create corrector
    r = new Vector(Rsize);
    // correctors (long and boring section...)
    switch (size_Shake)
    {
    case 3:
        if (m == size_Shake)
        {
            // k3m3 orig
            r->operator()(0) = 1.0 / 6.0;
        }
        else if (integrator.find("v") != std::string::npos)
        {
            // k3m2v velocity
            r->operator()(0) = 1.0 / 3.0;
        }
        else
        {
            // k3m2e energy (velocity Verlet equivalent)
            r->operator()(0) = 0.0;
        }
        r->operator()(1) = 1.0;
        r->operator()(2) = 1.0;
        break;
    case 4:
        if (m == 4)
        {
            // k4m4 orig
            r->operator()(0) = 1.0 / 6.0;
            r->operator()(1) = 5.0 / 6.0;
        }
        else if ((m == 3) && (integrator.find("v") != std::string::npos))
        {
            // k4m3v velocity
            r->operator()(0) = 1.0 / 4.0;
            r->operator()(1) = 5.0 / 6.0;
        }
        else if ((m == 2) && (integrator.find("e") != std::string::npos))
        {
            // k4m2e energy (Verlet -- Beeman equivalent)
            r->operator()(0) = 0.0;
            r->operator()(1) = 2.0 / 3.0;
        }
        else
        {
            // unknown integrator, size defined (4) --> default: k4m2e
            print_warning(1, "Unknown integrator name: " + integrator + ".\n", "Using default integrator of the same size: k4m2e.\n");
            // k4m2e energy (Verlet -- Beeman equivalent)
            r->operator()(0) = 0.0;
            r->operator()(1) = 2.0 / 3.0;
        }
        r->operator()(2) = 1.0;
        r->operator()(3) = 1.0 / 3.0;
        break;
    case 5:
        if (m == size)
        {
            // k5m5 orig
            r->operator()(0) = 1.9 / 12.0;
            r->operator()(1) = 3.0 / 4.0;
        }
        else if ((m == 2) && (integrator.find("e") != std::string::npos))
        {
            // k5m2e -- new Verlet equivalent
            r->operator()(0) = 0.0;
            r->operator()(1) = 0.0;
        }
        else if ((m == 4) && (integrator.find("v") != std::string::npos))
        {
            // k5m4v velocity
            r->operator()(0) = 1.9 / 9.0;
            r->operator()(1) = 3.0 / 4.0;
        }
        else if ((m == 4) && (integrator.find("s") != std::string::npos))
        {
            // k5m4s shake version
            r->operator()(0) = 0.0;
            r->operator()(1) = 3.0 / 4.0;
        }
        else if ((m == 4) && (integrator.find("e") != std::string::npos))
        {
            // k5m4e energy
            r->operator()(0) = 1.0 / 12.0;
            r->operator()(1) = 3.0 / 4.0;
        }
        else
        {
            // unknown integrator, size defined (5) --> default: k5m4e
            print_warning(1, "Unknown integrator name: " + integrator + ".\n", "Using default integrator of the same size: k5m4e.\n");
            // k5m4e energy
            r->operator()(0) = 1.0 / 12.0;
            r->operator()(1) = 3.0 / 4.0;
        }
        r->operator()(2) = 1.0;
        r->operator()(3) = 0.5;
        r->operator()(4) = 1.0 / 12.0;
        break;
    case 6:
        if (m == size)
        {
            // k6m6 orig
            r->operator()(0) = 3.0 / 20.0;
            r->operator()(1) = 25.1 / 36.0;
        }
        else if ((m == 5) && (integrator.find("v") != std::string::npos))
        {
            // k6m5v velocity
            r->operator()(0) = 3.0 / 16.0;
            r->operator()(1) = 25.1 / 36.0;
        }
        else if ((m == 4) && (integrator.find("s") != std::string::npos))
        {
            // k6m4s shake
            r->operator()(0) = 0.0;
            r->operator()(1) = 56.0 / 90.0;
        }
        else if ((m == 4) && (integrator.find("e") != std::string::npos))
        {
            // k6m4e energy
            r->operator()(0) = 1.0 / 30.0;
            r->operator()(1) = 2.3 / 3.6;
        }
        else
        {
            // unknown integrator, size defined (6) --> default: k6m4e
            print_warning(1, "Unknown integrator name: " + integrator + ".\n", "Using default integrator of the same size: k6m4e.\n");
            // k6m4e energy
            r->operator()(0) = 1.0 / 30.0;
            r->operator()(1) = 2.3 / 3.6;
        }
        r->operator()(2) = 1.0;
        r->operator()(3) = 1.1 / 1.8;
        r->operator()(4) = 1.0 / 6.0;
        r->operator()(5) = 1.0 / 60.0;
        break;
    case 7:
        if (m == size)
        {
            // k7m7 orig
            r->operator()(0) = 86.3 / 604.8;
        }
        else if ((m == 6) && (integrator.find("v") != std::string::npos))
        {
            // k7m6v velocity
            r->operator()(0) = 86.3 / 504.0;
        }
        else if ((m == 6) && (integrator.find("s") != std::string::npos))
        {
            // k7m6s shake
            r->operator()(0) = 0.0;
        }
        else if ((m == 6) && (integrator.find("e") != std::string::npos))
        {
            // k7m6e energy
            r->operator()(0) = 7.0 / 72.0;
        }
        else
        {
            // unknown integrator, size defined (7) --> default: k7m4e
            print_warning(1, "Unknown integrator name: " + integrator + ".\n", "Using default integrator of the same size: k7m6e.\n");
            // k7m6e energy
            r->operator()(0) = 7.0 / 72.0;
        }
        r->operator()(1) = 9.5 / 14.4;
        r->operator()(2) = 1.0;
        r->operator()(3) = 2.5 / 3.6;
        r->operator()(4) = 3.5 / 14.4;
        r->operator()(5) = 1.0 / 24.0;
        r->operator()(6) = 1.0 / 360.0;
        break;
    case 8:
        if (m == 8)
        {
            // k8m8 orig
            r->operator()(0) = 27.5 / 201.6;
            r->operator()(1) = 19.087 / 30.24;
        }
        else if ((m == 7) && (integrator.find("v") != std::string::npos))
        {
            // k8m7v velocity
            r->operator()(0) = 27.5 / 172.8;
            r->operator()(1) = 19.087 / 30.24;
        }
        else if ((m == 6) && (integrator.find("s") != std::string::npos))
        {
            // k8m6s shake
            r->operator()(0) = 0.0;
            r->operator()(1) = 12.3 / 21.0;
        }
        else if ((m == 6) && (integrator.find("e") != std::string::npos))
        {
            // k8m6e energy
            r->operator()(0) = 4.3 / 84.0;
            r->operator()(1) = 2.17 / 3.6;
        }
        else
        {
            // unknown integrator, size defined (8) --> default: k8m6e
            print_warning(1, "Unknown integrator name: " + integrator + ".\n", "Using default integrator of the same size: k8m6e.\n");
            // k8m6e energy
            r->operator()(0) = 4.3 / 84.0;
            r->operator()(1) = 2.17 / 3.6;
        }
        r->operator()(2) = 1.0;
        r->operator()(3) = 1.37 / 1.8;
        r->operator()(4) = 5.0 / 16.0;
        r->operator()(5) = 0.17 / 2.4;
        r->operator()(6) = 1.0 / 120.0;
        r->operator()(7) = 1.0 / 2520.0;
        break;
    case 9:
        if (m == 9)
        {
            // k9m9 orig
            r->operator()(0) = 33.953 / 259.2;
        }
        else if ((m == 8) && (integrator.find("v") != std::string::npos))
        {
            // k9m8v velocity
            r->operator()(0) = 33.953 / 226.8;
        }
        else if ((m == 8) && (integrator.find("s") != std::string::npos))
        {
            // k9m8s shake
            r->operator()(0) = 0.0;
        }
        else if ((m == 8) && (integrator.find("e") != std::string::npos))
        {
            // k9m8e energy
            r->operator()(0) = 0.859 / 8.64;
        }
        else
        {
            // unknown integrator, size defined (9) --> default: k9m8e
            print_warning(1, "Unknown integrator name: " + integrator + ".\n", "Using default integrator of the same size: k9m8e.\n");
            // k9m8e energy
            r->operator()(0) = 0.859 / 8.64;
        }
        r->operator()(1) = 5.257 / 8.64;
        r->operator()(2) = 1.0;
        r->operator()(3) = 4.9 / 6.0;
        r->operator()(4) = 2.03 / 5.4;
        r->operator()(5) = 0.49 / 4.8;
        r->operator()(6) = 0.7 / 43.2;
        r->operator()(7) = 1.0 / 720.0;
        r->operator()(8) = 1.0 / 20160.0;
        break;
    }
    // if TRVP, then corrector[predVel+1] must be the same as corrector[0]
    if (predVel != 1)
    {
        r->operator()(predVel + 1) = r->operator()(0);
    }
    // predictor
    A = new Matrix(Rsize, Rsize);
    // if not TRVP, binomial coeffs
    if (predVel == 1)
    {
        for (i = 0; i < Rsize; i++)
        {
            for (j = i; j < Rsize; j++)
            {
                A->operator()(i, j) = (double)(fact(j) / (fact(j - i) * fact(i)));
            }
        }
    }
    else // if TRVP...
    {
        // until row and column predVel, binomial coeffs
        for (i = 0; i < predVel; i++)
        {
            for (j = i; j < predVel; j++)
            {
                A->operator()(i, j) = (double)(fact(j) / (fact(j - i) * fact(i)));
            }
        }
        // row predVel: prediction (TRVP) (see Gear_formalism_for_MD)
        for (i = 1; i < predVel; i++)
        {
            A->operator()(predVel, i) = trvpB(0);
        }
        j = 1;
        for (i = predVel + 1; i < Rsize; i++)
        {
            A->operator()(predVel, i) = trvpB(j);
            j++;
        }
        // row predVel + 1: r(t+h)^P - r(t) (for Verlet r(t+h)^P = r(t+h)^C)
        for (i = 1; i < predVel; i++)
        {
            A->operator()(predVel + 1, i) = 1.0;
        }
        // next rows: time shift
        for (i = predVel + 2; i < Rsize; i++)
        {
            A->operator()(i, i - 1) = 1.0;
        }
    }
    // k5m2eX -- different predictor, based on value of X (parameter k = X/12)
    if ((size_Shake == 5) && (m == 2))
    {
        // if name contains value of free parameter k (X)
        if (strlen(integName) > 5)
        {
            // name should be `k5m2eX_NUMBER_`
            if (toupper(integName[5]) != 'X')
            {
                print_warning(1, "Wrong integrator name (5-sized Verlet). Expect X in the name string but obtained: " + std::string(integName) + ". Defaulted to 5-sized Verlet (1/3).");
                X = 4.0;
            }
            else
            {
                // parse X
                X = std::atof(integName + 6);
            }
            // X should be 'small number' [-12 .. 24]
            if ((X < -12) || (X > 24))
            {
                print_warning(1, "Wrong free parameter to 5-sized Verlet: " + std::to_string(X / 12.0) + "\n", "Defaulted to 1/3.\n");
                X = 1.0 / 3.0;
            }
            // 'k' (saved as X) = X/12
            else
            {
                X /= 12.0;
            }
        }
        else // default value 1/3
        {
            X = 1.0 / 3.0;
        }
        // modify predictor
        A->operator()(0, 3) = 3.0 - 6.0 * X;
        A->operator()(0, 4) = 36.0 * X - 6.0;
        A->operator()(1, 4) = 24.0 * X - 6.0;
    }
    // initialize 'sumr'
    sumr = ((*A) * (*r))[0];

    // omega and relative error for shake
    omega_Shake = (simsystem == nullptr) ? (1.0) : (simsystem->omegaShake); // overrelaxation
    eps_Shake = (simsystem == nullptr) ? (1.0) : (simsystem->epsShake);
    // max number of iterations for SHAKE
    maxiter_Shake = 100;
    // rescaling for SHAKE
    lambda_Shake = 1.0;
    // initialize configuration
    Initialize(startWithShake);
    // positions to be rescaled when box is being rescaled
    positionsAt.push_back(0);
}

// destructor deletes predictor matrix and corrector vector created during construction
AbstractGearIntegrator::~AbstractGearIntegrator()
{
    // delete predictor matrix
    delete A;
    // delete corrector vector;
    delete r;
}

// initialize values in R in atoms...
int AbstractGearIntegrator::Initialize(bool startWithShake)
{
    if (system == nullptr)
    {
        return 1;
    }
    int i, j, k, m;
    Matrix *atomR;
    double temp;
    bool accGiven = false;

    Matrix G(1, 3);

    // best agreement with continue is achieved when old accelerations are used (if present)
    // but then the pressure is off, due to the wrong virial of constraints
    // the solution is to calculate forces anew, perform shake (without correction2), multiply virconstr by sumr
    // and finally assign back the original value of acceleration

    // system->instEpot =
    system->CalculateForces(); // returns potential energy (saved in SimulatedSystem itself)

    for (i = 0; i < system->noMolecules; i++)
    {
        for (j = 0; j < system->molecules[i].noAtoms; j++)
        {
            atomR = system->molecules[i].atoms[j].R;
            for (k = 0; k < 3; k++)
            {
                ATOMR(1, k) = ATOMR(1, k) * h; // v(i) -> h*v(i)
                if ((!accGiven) && (ATOMR(2, k) == 0.0))
                {
                    ATOMR(2, k) = 0.5 * h * h * system->molecules[i].atoms[j].force[k] / system->molecules[i].atoms[j].mass;
                    // a(i) = h^2/2*f(i)/m
                }
                else
                {
                    accGiven = true;
                    temp = ATOMR(2, k);
                    ATOMR(2, k) = 0.5 * h * h * system->molecules[i].atoms[j].force[k] / system->molecules[i].atoms[j].mass;
                    // ATOMR(2, k) = 0.5 * h * h * ATOMR(2, k);
                    system->molecules[i].atoms[j].force[k] = temp; // temporarily save the original acceleration value to force
                }
                // higher derivatives
                for (m = 3; m < size_Shake; m++)
                {
                    if (ATOMR(m, k) != 0.0)
                    {
                        ATOMR(m, k) *= pow(h, (double)m) / fact(m);
                    }
                }
                // TRVP – stored past velocities
                if (predVel > 1)
                {
                    for (m = predVel; m < size; m++)
                    {
                        if (ATOMR(m, k) != 0.0)
                        {
                            ATOMR(m, k) *= h;
                        }
                        else
                        {
                            ATOMR(m, k) = ATOMR(1, k) - (m - predVel) * ATOMR(2, k);
                        }
                    }
                }
            }
        }
    }

    if (startWithShake) // default true maybe unnecessary
    {
        system->virconstr = 0.0;
        for (i = 0; i < system->noMolecules; i++)
        {
            ShakeHook(&(system->molecules[i])); // calculate Shake forces
        }
        if (accGiven) // if accelerations given, do not perform correction after Shake (not to spoil the values), but multiply the virconstr by sum of corrector coeffs
        {
            system->virconstr *= sumr;
            for (i = 0; i < system->noMolecules; i++)
            {
                for (j = 0; j < system->molecules[i].noAtoms; j++)
                {
                    atomR = system->molecules[i].atoms[j].R;
                    for (k = 0; k < 3; k++)
                    {
                        ATOMR(2, k) = 0.5 * h * h * system->molecules[i].atoms[j].force[k]; // assign previously temporarily saved values of accelerations
                    }
                }
            }
        }
        else // if accelerations not given, perform correction as usual
        {
            for (i = 0; i < system->noMolecules; i++)
            {
                Correction(&(system->molecules[i])); // correction after Shake
            }
        }
    }
    else // to be erased in future
    {
        if (accGiven) // if accelerations given, do not perform correction after Shake (not to spoil the values), but multiply the virconstr by sum of corrector coeffs
        {
            system->virconstr = 0.0;
            for (i = 0; i < system->noMolecules; i++)
            {
                for (j = 0; j < system->molecules[i].noAtoms; j++)
                {
                    atomR = system->molecules[i].atoms[j].R;
                    for (k = 0; k < 3; k++)
                    {
                        ATOMR(2, k) = 0.5 * h * h * system->molecules[i].atoms[j].force[k]; // assign previously temporarily saved values of accelerations
                    }
                }
            }
        }
    }

    return 0;
}

// prepare R for .config printing
int AbstractGearIntegrator::PrepareConfig()
{
    int configLevel = size_Shake - 1;

    // actualize potential energy (works only if not dump)
    // if dump, system is incomplete and has no forces (and pairlist == nullptr)
    if (system->pairlist != nullptr)
    {
        system->CalculateForces();
    }

    if (predVel != 1)
    {
        configLevel += 10 * (size - predVel - 1);
    }
    // return (predVel == 1)?r->GetNumberOfRows():(predVel - 1);
    return configLevel;
}

// correction
int AbstractGearIntegrator::Correction(Molecule *mol)
{
    int j;

    for (j = 0; j < mol->noAtoms; j++)
    {
        // correction
        *(mol->atoms[j].R) += CalculateOuterProduct(*r, *(mol->atoms[j].errorG));
        mol->atoms[j].errorG->Clear();
    }
    return 0;
}

// SHAKE hook
double AbstractGearIntegrator::ShakeHook(Molecule *mol)
{
    int i, j, k;
    double local_virial = 0.0;
    double maxRelErr = 0.0;
    double relErr = 0.0;

    for (i = 0; i < mol->noAtoms; i++)
    {
        // prepare errorG by multiplying by sumr
        *(mol->atoms[i].errorG) *= sumr;
        // prepare rPI (not including errorG (actual value to be added when calculating forces))
        (*mol->atoms[i].rPI) = CalculatePartialProduct(*A, *mol->atoms[i].R, 0);
    }

    // iterations
    for (i = 0; i < maxiter_Shake; i++)
    {
        // zero maxRelErr
        maxRelErr = 0.0;
        for (j = 0; j < mol->noConstrBonds; j++)
        {
            // Calculate forces (consider boxscaling) and add them to errorG
            relErr = mol->constrBonds[j].CalculateForces(local_virial, size_Shake, eps_Shake, lambda_Shake, omega_Shake);
            for (k = 0; k < 3; k++)
            {
                mol->atoms[mol->constrBonds[j].atomI].errorG->operator()(k) +=
                    mol->atoms[mol->constrBonds[j].atomI].constraintForces[k];
                mol->atoms[mol->constrBonds[j].atomJ].errorG->operator()(k) +=
                    mol->atoms[mol->constrBonds[j].atomJ].constraintForces[k];
            }
            // Check relErr with current maxRelErr
            maxRelErr = (relErr > maxRelErr) ? (relErr) : (maxRelErr);
        }
        // check if SHAKE convergency reached
        if (maxRelErr < eps_Shake)
        {
            break;
        }
    }

    // update global reached maximum of SHAKE iterations system->maxShakeIter
    (i > system->maxShakeIter) ? (system->maxShakeIter = i) : (true);
    // update system virconstr
    // scaling of the virial according to the comparison of g(r_i), h(r_i) and f(r_i)
#ifdef PARALLEL
#pragma omp atomic
#endif
    system->virconstr += local_virial * 2.0 / sumr / h / h;

    // restore errorG by dividing by sumr
    for (i = 0; i < mol->noAtoms; i++)
    {
        // prepare contraintErrorG to correction (divide by the sum of corrector coefficients not known for SHAKE itself in ConstrainedBond)
        *(mol->atoms[i].errorG) *= (1.0 / sumr);
    }
    // return maximum  relative error
    return maxRelErr;
}
