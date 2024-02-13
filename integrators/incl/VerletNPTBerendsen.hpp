#ifndef VERLETNPTBERENDSENHEADER
#define VERLETNPTBERENDSENHEADER

#include "AbstractVerletIntegrator.hpp"

class VerletNPTBerendsen : public AbstractVerletIntegrator
{
private:
    double tau_T;
    double T_f; // final (desired) temperature
    double tau_P;
    double P_f;             // final (desired) pressure
    double lambda_v; // temperature scaling factor (1 + h/tau[T]*(Tfinal/Tkin - 1))^(1/2)
    double lambda_P; // volume scaling factor to reach desired pressure (1 + BETA*h/tau[P]*(P_f - P_inst)) where BETA is the value of kompresibility
    int VerletStep(Molecule *mol, double lambda_v);

public:
    int Integrate(int Nsteps) override;

    VerletNPTBerendsen(double step, SimulatedSystem *simsystem, double Tfinal, double tauT, double Pfinal, double tauP, bool startWithShake);
    VerletNPTBerendsen(double step, SimulatedSystem *simsystem, Vector trvpB, double Tfinal, double tauT, double Pfinal, double tauP, bool startWithShake);
    ~VerletNPTBerendsen();
};

#endif