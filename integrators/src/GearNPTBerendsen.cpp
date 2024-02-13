/*
 * An implementation of Gear integration method with Berendsen thermostat + friction (Berendsen) barostat
 * for simul++ simulation package
 * Author JJ, Date Feb 2022
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
 * Berendsen (friction) barostat implemented according to DL POLY manual (p. 94)(MACSIMUS manual lack detailed information...)
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
#include "GearNPTBerendsen.hpp"

#define ATOMR(i, j) ((*atomR)(i, j))
// #define DEBUG
#define BETA 0.00625 // isothermal kompresibility (1/bulk modulus) (aproximation for liquid water (very roughly))

// initialization of this integrator (main part is done in AbstractGearIntegrator construction)
GearNPTBerendsen::GearNPTBerendsen(double step, SimulatedSystem *simsystem, int Rsize, char *integName, double Tfinal, double tauT, double Pfinal, double tauP, Vector trvpB, bool startWithShake)
    : AbstractGearIntegrator{step, simsystem, Rsize, integName, trvpB, startWithShake}
{
    Tf = Tfinal;
    tau_T = tauT;
    lambda_T = 1.0;
    Tinst = Tf;
    Pf = Pfinal;
    tau_P = tauP;
    molecular_based = true;
    lambda_P = 1.0;
}

// empty destructor
GearNPTBerendsen::~GearNPTBerendsen()
{
}

// temperature scaling
inline double GearNPTBerendsen::GetTempScaling()
{
    return 1 / (2 * tau_T * h) * log(Tinst / Tf); // MACSIMUS
    // lambda = -1/(2*tau_T*h)*(T_f * system->noDegrOfFred/Ekininst - 1) -- original Berendsen article
}

// box scaling
inline double GearNPTBerendsen::GetScalingFactor(int order)
{
    double Pinst;

    if (order == 2) // after correction
    {
        Pinst = system->CalculatePkin(h) + system->CalculatePconf();
        lambda_P = 1.0 - BETA * h / tau_P * (Pf - Pinst);
        return pow(lambda_P, 1.0 / 3.0);
    }
    else // after prediction
    {
        lambda_P = 1.0;
    }
    return 1.0;
}

inline int GearNPTBerendsen::CalculateRHS(Atom *at)
{
    int i;

    for (i = 0; i < 3; i++)
    {
        // G = h^2/2*(f/m - lambda*v)  - h^2/2*a^P – Berendsen according to MACSIMUS manual p. 162
        at->errorG->operator()(i) = 0.5 * h * h * (at->force[i] / at->mass - lambda_T * at->R->operator()(1, i)) - at->R->operator()(2, i);
    }

    return 0;
}

// perform N  steps of integration
int GearNPTBerendsen::Integrate(int Nsteps)
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
