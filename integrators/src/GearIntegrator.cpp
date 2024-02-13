/*
 * An implementation of Gear integration method for trial simulations of constrained dynamics
 * Based on the general Gear integration method in AbstractGearIntegrator class
 * Author JJ, Date Jan 2021
 */

/*
 * In atoms R following values are stored:
 * R(0, :): x(i) y(i) z(i)             common for all integrators
 * R(1, :): h*vx(i) h*vy(i) h*vz(i)    common for all integrators (h* because of Gear)
 * R(2, :): h^2/2*ax(i)...             Gear specific
 * R(3, :): h^3/6*x'''(i)...           Gear specific
 * R(4, :): h^4/24*x^(iv)(i)...        Gear specific
 *
 * R(0, :) and R(1, :) must be defined before integrator initialization
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
#include "GearIntegrator.hpp"

#define ATOMR(i, j) ((*atomR)(i, j))

// initialization of this integrator
GearIntegrator::GearIntegrator(double step, SimulatedSystem *simsystem, int Rsize, char *integName, Vector trvpB, bool startWithShake)
    : AbstractGearIntegrator{step, simsystem, Rsize, integName, trvpB, startWithShake}
{
    // constructor of AbstractGearIntegrator should be sufficient...
}

// destructor
GearIntegrator::~GearIntegrator()
{
    // should call destructor of AbstractGearIntegrator automatically
}

// calculate right-hand side of the differential equation for physical DoFs
inline int GearIntegrator::CalculateRHS(Atom *at)
{
    int i;

    for(i = 0; i < 3; i++)
    {
        // G0 = h^2/2*f/m - h^2/2*a^P
        at->errorG->operator()(i) = 0.5 * h * h * at->force[i] / at->mass - at->R->operator()(2, i);
    }

    return 0;
}

// perform N  steps of integration
int GearIntegrator::Integrate(int Nsteps)
{
    #ifdef GEARTHERMO
    #undef GEARTHERMO
    #endif
    #ifdef GEARBOXRESC
    #undef GEARBOXRESC
    #endif
    #ifdef GEAREXTDOFS
    #undef GEAREXTDOFS
    #endif
    #include "GearIntegrateMacro.cpp"
}

// other functions are inherited from AbstractGearIntegrator