#ifndef VERLETNVTNOSEITERHEADER
#define VERLETNVTNOSEITERHEADER

#include "AbstractVerletIntegrator.hpp"

class VerletNVTNoseIter: public AbstractVerletIntegrator
{      
    private:
        double tau_T;
        double T_f; // final (desired) temperature 
        Vector *Xi; // extended variable – thermostat
        Vector *Xi_old;
        int noiter; // number of iterations
        int VerletStep(Molecule *mol, double hvXiPredold, bool isfirst);
        int SimplePredict(Molecule *mol);
        int StoreVelocity(Molecule *mol);
        int VelocityToPred(Molecule *mol);

    public:
        int Integrate(int Nsteps) override;
        double EnergyOfExtraDOF() const override;
        int InitXi(int predVel);

        VerletNVTNoseIter(double step, SimulatedSystem *simsystem, double Tfinal, double tauT, bool startWithShake, int noIter);
        VerletNVTNoseIter(double step, SimulatedSystem *simsystem, Vector trvpB, double Tfinal, double tauT, bool startWithShake, int noIter);
        ~VerletNVTNoseIter();
        int PrepareConfig() override;
        int RestoreExtendedDOFs() override;
};

#endif