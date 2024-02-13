#ifndef GEARNVTNOSEHEADER
#define GEARNVTNOSEHEADER

#include "AbstractIntegrator.hpp"
#include "AbstractGearIntegrator.hpp"
class Matrix;
class Vector;

class GearNVTNose : public AbstractGearIntegrator
{
private:
    int Integrate(int Nsteps) override;
    inline int CalculateRHS(Atom *at) override;
    // thermostat variables and functions
    double Tf;       // final (desired) temperature
    double tau_T;    // thermostat constant
    double lambda_T; // velocity scaling due to thermostat
    double Tinst;    // instanteneous temperature
    inline double GetTempScaling();
    // extended degrees of freedom variables and functions
    Matrix *ExtDoFs;                          // extended degrees of freedom (and their scaled derivatives) (as Atom::R)
    Vector *GForExtDoFs;                      // prediction error for extra degrees of freedom
    inline double RHSForExtDoF(int extDoFno); // right-hand side of the differential equation for EDoF number @p extDoFno
    inline double EnergyOfExtDoF();           // extra energy (to Hamiltonian) due to EDoFs

public:
    GearNVTNose(double step, SimulatedSystem *simsystem, int Rsize, char *integName, double Tfinal, double tauT, Vector trvpB, bool startWithShake = true);
    ~GearNVTNose();
    int InitXi();
    int PrepareConfig() override;
};

#endif