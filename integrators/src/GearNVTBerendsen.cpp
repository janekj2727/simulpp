/*
 * An implementation of Gear integration method with Berendsen thermostat
 * for trial simulations of constrained dynamics
 * Author JJ, Date Apr 2021
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
 */

/*
 * Berendsen thermostat implemented according to differential equation in MACSIMUS manual p. 162
 * not only velocity rescaling, but also change of higher order derivatives in R...
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
#include "GearNVTBerendsen.hpp"

#define ATOMR(i, j) ((*atomR)(i, j))
// #define DEBUG

// initialization of this integrator (main part is done in AbstractGearIntegrator construction) + TRVP
GearNVTBerendsen::GearNVTBerendsen(double step, SimulatedSystem *simsystem, int Rsize, char*integName, double Tfinal, double tauT, Vector trvpB, bool startWithShake)
    : AbstractGearIntegrator{step, simsystem, Rsize, integName, trvpB, startWithShake}
{
    Tf = Tfinal;
    tau_T = tauT;
    lambda_T = 1.0;
    Tinst = Tf;
}

// empty destructor
GearNVTBerendsen::~GearNVTBerendsen()
{
}

// temperature scaling
inline double GearNVTBerendsen::GetTempScaling()
{
    return 1.0 / (2.0 * tau_T * h) * log(Tinst / Tf); // MACSIMUS
    // lambda = -1/(2*tau_T*h)*(T_f * system->noDegrOfFred/Ekininst - 1) -- original Berendsen article
}

// right-hand side of the differential equation for physical degrees of freedom
inline int GearNVTBerendsen::CalculateRHS(Atom *at)
{
    int i;

    for (i = 0; i < 3; i++)
    {
        // G = h^2/2*(f/m - lambda*v)  - h^2/2*a^P – Berendsen according to MACSIMUS manual p. 162
        at->errorG->operator()(i) = 0.5 * h * h * (at->force[i] / at->mass - lambda_T * at->R->operator()(1, i)) - at->R->operator()(2, i);
    }

    return 0;
}

// perform N steps of integration
int GearNVTBerendsen::Integrate(int Nsteps)
{
    #ifndef GEARTHERMO
    #define GEARTHERMO
    #endif
    #ifdef GEARBOXRESC
    #undef GEARBOXRESC
    #endif
    #ifdef GEAREXTDOFS
    #undef GEAREXTDOFS
    #endif
    #include "GearIntegrateMacro.cpp"
}
