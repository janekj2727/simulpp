#ifndef VERLETINTEGRATORHEADER
#define VERLETINTEGRATORHEADER

#include "AbstractVerletIntegrator.hpp"

class VerletIntegrator: public AbstractVerletIntegrator
{      
    public:
        int Integrate(int Nsteps) override;

        VerletIntegrator(double step, SimulatedSystem *simsystem, bool startWithShake);
        VerletIntegrator(double step, SimulatedSystem *simsystem, Vector trvpB, bool startWithShake);
        ~VerletIntegrator();

    private:
        int VerletStep(Molecule *mol) override;

};

#endif