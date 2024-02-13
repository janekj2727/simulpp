/*
 * An implementation of Gear integration method with Berendsen thermostat
 * and volume scaling to final rho
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

/*
 * Requested rho attained in tau.rho ps linearly (the simplest possible solution)
 * Energy increase must be dumped by Berendsen thermostat
 * It is the same as GearNVTBerendsen, the only change is the volume scaling...
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
#include "GearNTVinitBerendsen.hpp"

#define ATOMR(i, j) ((*atomR)(i, j))

// initialization of this integrator + TRVP
GearNTVinitBerendsen::GearNTVinitBerendsen(double step, SimulatedSystem *simsystem, int Rsize, char *integName, double Tfinal, double tauT, double rhofinal, double taurho, Vector trvpB, bool startWithShake)
    : AbstractGearIntegrator{step, simsystem, Rsize, integName, trvpB, startWithShake}
{
    Tf = Tfinal;
    tau_T = tauT;
    lambda_T = 1.0;
    Tinst = Tf;
    Pf = rhofinal;  // rho_final instead of P_f (pressure)
    tau_P = taurho; // tau_rho instead of tau_P
    molecular_based = true;
    lambda_P = pow(system->totalMass / rhofinal / system->CalculateVolume(), h / tau_P / 3.0);
}

// empty destructor
GearNTVinitBerendsen::~GearNTVinitBerendsen()
{
}

// box scaling first rough approximation, but should work
inline double GearNTVinitBerendsen::GetScalingFactor(int order)
{
    double ideal_scaling = pow(system->totalMass / Pf / system->CalculateVolume(), 1.0 / 3.0);

    if (order == 2)
    {
        if (fabs(1.0 - ideal_scaling) < fabs(1.0 - lambda_P))
        {
            return ideal_scaling;
        }
        else
        {
            return lambda_P;
        }
    }
    else
    {
        lambda_P = 1.0;
        return 1.0;
    }

    return 1.0;
}

// temperature scaling
inline double GearNTVinitBerendsen::GetTempScaling()
{
    return 1 / (2 * tau_T * h) * log(Tinst / Tf); // MACSIMUS
    // lambda = -1/(2*tau_T*h)*(T_f * system->noDegrOfFred/Ekininst - 1) -- original Berendsen article
}

// right-hand side of the differential equation for physical degrees of freedom
inline int GearNTVinitBerendsen::CalculateRHS(Atom *at)
{
    int i;

    for (i = 0; i < 3; i++)
    {
        // G = h^2/2*(f/m - lambda*v)  - h^2/2*a^P – Berendsen according to MACSIMUS manual p. 162
        at->errorG->operator()(i) = 0.5 * h * h * (at->force[i] / at->mass - lambda_T * at->R->operator()(1, i)) - at->R->operator()(2, i);
    }

    return 0;
}

// perform N usteps of integration
int GearNTVinitBerendsen::Integrate(int Nsteps)
{
#ifndef GEARTHERMO
#define GEARTHERMO
#endif
#ifndef GEARBOXRESC
#define GEARBOXRESC
#endif
#ifdef GEAREXTDOFS
#undef GEAREXTDOFS
#endif
#include "GearIntegrateMacro.cpp"
}
