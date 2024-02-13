#ifndef GEARNPTNOSEHEADER
#define GEARNPTNOSEHEADER

#include "AbstractIntegrator.hpp"
#include "AbstractGearIntegrator.hpp"
class Matrix;
class Vector;

class GearNPTNose : public AbstractGearIntegrator
{
private:
    int InitXi(double Ekininst);
    int InitLambda(double Ekininst);
    int Integrate(int Nsteps) override;
    inline int CalculateRHS(Atom *at) override;
    // thermostat variables and functions
    double Tf; // final (desired) temperature
    double tau_T; // thermostat constant
    double lambda_T; // velocity scaling due to thermostat
    double Tinst; // instanteneous temperature
    inline double GetTempScaling();
    // barostat variables and functions
    double Pf; // final (desired) pressure
    double tau_P; // barostat constant
    bool molecular_based; // scaling is molecular/atomic based (NPTBerendsen implemented with molecular scaling -> true)
    double lambda_P; // predicted value of box size
    inline double GetScalingFactor(int order);
    inline double GetShakeScaling();
    // extended degrees of freedom variables and functions
    Matrix *ExtDoFs; // extended degrees of freedom (and their scaled derivatives) (as Atom::R)
    Vector *GForExtDoFs; // prediction error for extra degrees of freedom
    inline double RHSForExtDoF(int extDoFno); // right-hand side of the differential equation for EDoF number @p extDoFno
    inline double EnergyOfExtDoF(); // extra energy (to Hamiltonian) due to EDoFs

public:
    GearNPTNose(double step, SimulatedSystem *simsystem, int Rsize, char *integName, double Tfinal, double tauT, double Pfinal, double tauP, Vector trvpB, bool startWithShake = true);
    ~GearNPTNose();
    int PrepareConfig() override;
};

#endif