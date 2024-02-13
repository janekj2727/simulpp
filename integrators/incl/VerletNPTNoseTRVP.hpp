#ifndef VERLETNPTNOSETRVPHEADER
#define VERLETNPTNOSETRVPHEADER

#include "AbstractVerletIntegrator.hpp"

class VerletNPTNoseTRVP : public AbstractVerletIntegrator
{
private:
    double tau_T;
    double T_f; // final (desired) temperature
    double tau_P;
    double P_f;         // final (desired) pressure
    Vector *Xi;         // extended variable – thermostat
    Vector *Lambda;     // extended variable – barostat
    Vector *Xi_old;     // backup of Xi
    Vector *Lambda_old; // backup of Lambda
    double M_T;         // mass of thermostat: M_T = N_f * k_B * T_f * tau_T^2
    double M_P;         // mass of barostat:   M_P = (N_f + 3) * k_B * T_f * tau_P^2
    double recM_T;      // 1/M_T
    double recM_P;      // 1/M_P
    int VerletStep(Molecule *mol, double hvXiPred, double hvLamPred);
    int TimeShiftVec(Vector *vec);
    double TRVPvec(Vector *vec);
    double VelocityVec(Vector *vec);
    int InitXi(int predVel, double Ekininst);
    int InitLambda(int predVel, double Ekininst);
    int TimeShift() override;
    int TRVP() override;
    double VelocityCalculation(Molecule *mol) override;
    VerletNPTNoseTRVP(double step, SimulatedSystem *simsystem, double Tfinal, double tauT, double Pfinal, double tauP, bool startWithShake); // this version is never used (without TRVP)

public:
    int Integrate(int Nsteps) override;
    double EnergyOfExtraDOF() const override;
    int PrepareConfig() override;
    int RestoreExtendedDOFs() override;

    VerletNPTNoseTRVP(double step, SimulatedSystem *simsystem, Vector trvpB, double Tfinal, double tauT, double Pfinal, double tauP, bool startWithShake);
    ~VerletNPTNoseTRVP();
};

#endif