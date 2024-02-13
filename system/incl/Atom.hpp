#ifndef ATOMHEADER
#define ATOMHEADER


#include "Matrix.hpp"
class Vector;

class Atom
{
private:
    Atom(){};

public:
    double mass;                // atom mass
    char name[10];              // atom label/name
    double charge;              // atom partial charge
    int LJid;                   // ID for non-bonded interactions
    Matrix *R;                  // positions velocities...
    double force[3];            // forces on atom (constraintForces not included)
    double constraintForces[3]; // forces caused by constraints multiplied by h^2/m
    Vector *rPI;                // predicted position from intermediate correction (Gear) R^I
    Vector *errorG;             // error of prediction (used for Gear)

    Atom(double atomMass, double atomCharge, char *atomName, int LJid, int Rsize);
    Atom(double atomMass, double atomCharge, char *atomName, int LJid, int Rsize, bool delayedalloc);
    Atom(const Atom &otherAtom);
    ~Atom();
    void ZeroForces();
    inline double GetPosition(int coord) const
    {
        return R->Read(0, coord);
    };
    Atom &operator=(const Atom &otherAtom);
};

#endif