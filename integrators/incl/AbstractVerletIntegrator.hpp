#ifndef ABSTRACTVERLETINTEGRATORHEADER
#define ABSTRACTVERLETINTEGRATORHEADER

#include "AbstractIntegrator.hpp"
class Matrix;
class Vector;
class Molecule;

class AbstractVerletIntegrator: public AbstractIntegrator
{      
    protected:
        // int predVel;
        Vector *trvpCoeff;

        virtual int Integrate(int Nsteps) = 0;
        int Initialize(bool startWithShake = true) override;
        virtual int Initialize(int predVel, bool startWithShake = true);
        int PrepareConfig() override;
        virtual double Shake(Molecule *mol, double &virial);
        virtual double Shake(Molecule *mol, double &virial, double rescaling);
        virtual int TimeShift();
        virtual double VelocityCalculation(Molecule *mol);
        virtual int VerletStep(Molecule *mol);
        virtual int TRVP();

        AbstractVerletIntegrator(double step, SimulatedSystem *simsystem, bool startWithShake);
        AbstractVerletIntegrator(double step, SimulatedSystem *simsystem, bool startWithShake, bool simpleinit);
        AbstractVerletIntegrator(double step, SimulatedSystem *simsystem, Vector trvpB, bool startWithShake);
        ~AbstractVerletIntegrator();

    public:
        virtual double EnergyOfExtraDOF() const;
};

#endif