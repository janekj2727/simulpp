#ifndef CONSTRAINTBONDHEADER
#define CONSTRAINTBONDHEADER

// #include "./mymath/Vector.hpp"

class Molecule;
// class Vector;

class ConstraintBond
{
    private:
        ConstraintBond() {};
        // Vector *rC2, *rP;
        double recMass; // reciprocal mass (1/(m1+m2))
        double redMass; // reduced mass (m1*m2/(m1+m2))
        
    public:
        double length; // constraint bond length
        int atomI; // index of first atom
        int atomJ; // index of second bounded atom
        Molecule *parent; // parent molecule (to know where to search for atoms)
        double relErr; // relative length error
        //bool trueBond; // true if bond, false if used to fix angle

        ConstraintBond(double bondLength, int indexI, int indexJ, Molecule *parentMol); //, bool trueB = true);
        ~ConstraintBond() {};
        double GetCurrentLength() const;
        double CalculateForces(double &virial, int shakeType = 0, double epsc = 1e-6, double scaling = 1.0, double omega = 1.0);
        double GetVelocityAngle() const;

};

#endif