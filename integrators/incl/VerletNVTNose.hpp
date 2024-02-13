#ifndef VERLETNVTNOSEHEADER
#define VERLETNVTNOSEHEADER

#include "AbstractIntegrator.hpp"
#include "AbstractVerletIntegrator.hpp"
class Matrix;
class Vector;

class VerletNVTNose: public AbstractVerletIntegrator
{      
    private:
        double T_f;
        double tau_T;
        Vector *Xi;
        Vector *Xi_old;
        int ThermostatHalfStep(double &Ekin);
        int VerletStep(Molecule *p_mol) override;
        int VerletStep2(Molecule *p_mol);
        int TimeShift() override;
        int InitXi();
        double Shake(Molecule *p_mol, double &virial) override;
        double Rattle(Molecule *p_mol, double &virial);
        
    public:
        VerletNVTNose(double step, SimulatedSystem *simsystem, double Tfinal, double tauT, bool startWithShake);
        ~VerletNVTNose();
        int Integrate(int Nsteps) override;
        double EnergyOfExtraDOF() const override;
        int PrepareConfig() override;
        int Initialize(bool startWithShake) override;   
        int RestoreExtendedDOFs() override;     
};

#endif