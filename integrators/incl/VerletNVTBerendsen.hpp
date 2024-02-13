#ifndef VERLETNVTBERENDSENHEADER
#define VERLETNVTBERENDSENHEADER

#include "AbstractVerletIntegrator.hpp"

class VerletNVTBerendsen: public AbstractVerletIntegrator
{      
    private:
        double tau_T;
        double T_f; // final (desired) temperature 
        double lambda; // temperature scaling factor (1 + h/tau[T]*(Tfinal/Tkin - 1))^(1/2)
        int VerletStep(Molecule *mol, double lambda);

    public:
        int Integrate(int Nsteps) override;

        VerletNVTBerendsen(double step, SimulatedSystem *simsystem, double Tfinal, double tauT, bool startWithShake);
        VerletNVTBerendsen(double step, SimulatedSystem *simsystem, Vector trvpB, double Tfinal, double tauT, bool startWithShake);
        ~VerletNVTBerendsen();

};

#endif