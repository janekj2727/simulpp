/*
 * A class to store atom info for trial simulations of contrained dynamics
 * Author JJ, Date Jan 2021
 */

#include <cstring>
#include <iostream>

#include "Vector.hpp"
#include "Matrix.hpp"
#include "Atom.hpp"

// #define DEBUG

// Default constructor
Atom::Atom(double atomMass, double atomCharge, char *atomName, int LJident, int Rsize)
{
    mass = atomMass;
    charge = atomCharge;
    LJid = LJident;
    strcpy(name, atomName);
    R = new Matrix(Rsize, 3);
    errorG = new Vector(3);
    rPI = new Vector(3);
}

// Constructor with delayed allocation of R
Atom::Atom(double atomMass, double atomCharge, char *atomName, int LJident, int Rsize, bool delayedalloc)
{
    if (delayedalloc)
    {
        mass = atomMass;
        charge = atomCharge;
        LJid = LJident;
        strcpy(name, atomName);
        R = nullptr;
        errorG = new Vector(3);
        rPI = new Vector(3);
    }
    else
    {
        mass = atomMass;
        charge = atomCharge;
        LJid = LJident;
        strcpy(name, atomName);
        R = new Matrix(Rsize, 3);
        errorG = new Vector(3);
        rPI = new Vector(3);
    }
}

// Destructor
Atom::~Atom()
{
    if (R != nullptr)
    {
        delete R;
        R = nullptr;
    }
    if (errorG != nullptr)
    {
        delete errorG;
        errorG = nullptr;
    }
    if (rPI != nullptr)
    {
        delete rPI;
        rPI = nullptr;
    }
}

// Zero forces before next forces evaluation
void Atom::ZeroForces()
{
    int i;

    for (i = 0; i < 3; i++)
    {
        force[i] = 0.0;
    }
}

// Copy constructor
Atom::Atom(const Atom &otherAtom)
{
    mass = otherAtom.mass;
    charge = otherAtom.charge;
    LJid = otherAtom.LJid;
    strcpy(name, otherAtom.name);
    if (otherAtom.R != nullptr)
    {
        R = new Matrix(*otherAtom.R);
    }
    else
    {
        R = nullptr;
    }
    
    errorG = new Vector(*otherAtom.errorG);
    rPI = new Vector(*otherAtom.rPI);
}

// Overloading operator=
Atom &Atom::operator=(const Atom &at)
{
    int k;

    mass = at.mass;
    charge = at.charge;
    LJid = at.LJid;
    strcpy(name, at.name);
    if (at.R != nullptr)
    {
        R = new Matrix(*at.R);
    }
    else
    {
        R = nullptr;
    }

    errorG = new Vector(*at.errorG);
    rPI = new Vector(*at.rPI);

    for (k = 0; k < 3; k++)
    {
        force[k] = at.force[k];
    }

    return *this;
}