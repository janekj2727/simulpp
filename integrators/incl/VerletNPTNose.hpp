#ifndef VERLETNPTNOSEHEADER
#define VERLETNPTNOSEHEADER

#include "AbstractIntegrator.hpp"
#include "AbstractVerletIntegrator.hpp"
class Matrix;
class Vector;

class VerletNPTNose: public AbstractVerletIntegrator
{      
    private:
        double T_f;
        double tau_T;
        double P_f;
        double tau_P;
        double M_T; // mass of thermostat: M_T = N_f * k_B * T_f * tau_T^2
        double M_P; // mass of barostat:   M_P = (N_f + 3) * k_B * T_f * tau_P^2
        Vector *Xi;
        Vector *Lambda;
        Vector *Xi_old;
        Vector *Lambda_old;
        double orig_size[3]; // original box_sizes
        int ThermostatQuarterStep(double &Ekin, double &velscale);
        int BarostatHalfStep(double &Ekin, double &velscale);
        int VerletStep(Molecule *p_mol) override;
        int VerletStep2(Molecule *p_mol);
        int TimeShift() override;
        int InitXi();
        int InitLambda();
        double Shake(Molecule *p_mol, double &virial) override;
        double Rattle(Molecule *p_mol, double &virial);
        
    public:
        VerletNPTNose(double step, SimulatedSystem *simsystem, double Tfinal, double tauT, double Pfinal, double tauP, bool startWithShake);
        ~VerletNPTNose();
        int Integrate(int Nsteps) override;
        double EnergyOfExtraDOF() const override;
        int PrepareConfig() override;
        int Initialize(bool startWithShake) override; 
        int RestoreExtendedDOFs() override;        
};

#endif