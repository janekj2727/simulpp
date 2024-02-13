#ifndef GEARINTEGRATORHEADER
#define GEARINTEGRATORHEADER

#include "AbstractIntegrator.hpp"
#include "AbstractGearIntegrator.hpp"
class Matrix;
class Vector;

class GearIntegrator: public AbstractGearIntegrator
{      
    private:
        int Integrate(int Nsteps) override;
        inline int CalculateRHS(Atom *at) override;

    public:
        GearIntegrator(double step, SimulatedSystem *simsystem, int Rsize, char *integName, Vector trvpB, bool startWithShake = true);
        ~GearIntegrator();
        
};

#endif