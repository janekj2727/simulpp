#ifndef VERLETNTVINITBERENDSENHEADER
#define VERLETNTVINITBERENDSENHEADER

#include "AbstractVerletIntegrator.hpp"

class VerletNTVinitBerendsen : public AbstractVerletIntegrator
{
private:
    double tau_T;
    double T_f;            // final (desired) temperature
    double rho_f;          // final (desired) density
    double tau_rho;        // time in which the final density should be reached
    double scaling_factor; // scaling of box in one step...
    int scaling_steps;     // number of steps in which scaling will be done
    int VerletStep(Molecule *mol, double lambda);

public:
    int Integrate(int Nsteps) override;

    VerletNTVinitBerendsen(double step, SimulatedSystem *simsystem, double Tfinal, double tauT, double rhofinal, double taurho, bool startWithShake);
    VerletNTVinitBerendsen(double step, SimulatedSystem *simsystem, Vector trvpB, double Tfinal, double tauT, double rhofinal, double taurho, bool startWithShake);
    ~VerletNTVinitBerendsen();
};

#endif