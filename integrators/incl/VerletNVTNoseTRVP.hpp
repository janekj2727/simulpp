#ifndef VERLETNVTNOSETRVPHEADER
#define VERLETNVTNOSETRVPHEADER

#include "AbstractVerletIntegrator.hpp"

class VerletNVTNoseTRVP: public AbstractVerletIntegrator
{      
    private:
        double tau_T;
        double T_f; // final (desired) temperature 
        Vector *Xi; // extended variable – thermostat
        Vector *Xi_old;
        int VerletStep(Molecule *mol, double hvXiPred);
        VerletNVTNoseTRVP(double step, SimulatedSystem *simsystem, double Tfinal, double tauT, bool startWithShake); // this version is never used (without TRVP)

    public:
        int Integrate(int Nsteps) override;
        double EnergyOfExtraDOF() const override;
        int InitXi(int predVel);


        VerletNVTNoseTRVP(double step, SimulatedSystem *simsystem, Vector trvpB, double Tfinal, double tauT, bool startWithShake);
        ~VerletNVTNoseTRVP();
        int PrepareConfig() override;
        int RestoreExtendedDOFs() override;
};

#endif