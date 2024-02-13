/*
 * A class to store bond info (simul++)
 * Author JJ, Date Jan 2021
 * version 2 (Jun 2022): restructuralization using AbstractIntraMolField
 */

#include <iostream>
#include <cassert>
#include <cmath>

#include "Vector.hpp"
#include "Matrix.hpp"
#include "Atom.hpp"
#include "Molecule.hpp"
#include "LJBond.hpp"

// Default constructor
LJBond::LJBond(Molecule *parentMol, int indexI, int indexJ, double sig, double eps)
    : AbstractIntraMolField{3, parentMol}
{
    sigma = sig;
    epsilon = eps;
    sigma6 = pow(sigma, 6.0);
    atomI = indexI;
    atomJ = indexJ;
}

// Copy constructor
LJBond::LJBond(const LJBond &otherBond)
    : AbstractIntraMolField{otherBond}
{
    sigma = otherBond.sigma;
    epsilon = otherBond.epsilon;
    sigma6 = otherBond.sigma6;
    atomI = otherBond.atomI;
    atomJ = otherBond.atomJ;
}

// Overloading operator =
LJBond &LJBond::operator=(const LJBond &otherBond)
{
    parent = NULL;
    type = otherBond.type;
    sigma = otherBond.sigma;
    sigma6 = otherBond.sigma6;
    epsilon = otherBond.epsilon;
    atomI = otherBond.atomI;
    atomJ = otherBond.atomJ;
    return *this;
}

// cloning
LJBond *LJBond::copy(Molecule *newparent) const
{
    LJBond *newBond = new LJBond(*this);
    newBond->SetParent(newparent);
    return newBond;
}

// Calculate current bond length
double LJBond::GetCurrentLength() const
{
    Vector delta(3);
    int i;

    for (i = 0; i < 3; i++)
    {
        delta[i] = parent->atoms[atomJ].GetPosition(i) - parent->atoms[atomI].GetPosition(i);
    }

    return delta.CalculateNorm();
}

// Calculate forces caused by this bond
double LJBond::CalculateForces()
{
    double l[3];
    double l2, l2i, l6i;
    double Elj, force;
    int i;

    // get current distance vector
    for (i = 0; i < 3; i++)
    {
        l[i] = parent->atoms[atomJ].GetPosition(i) - parent->atoms[atomI].GetPosition(i);
    }

    // calculate square of the atom-atom distance: its norm not needed
    l2 = l[0] * l[0] + l[1] * l[1] + l[2] * l[2];
    // inverse of the square atom distance
    l2i = 1 / l2;
    // 1/(l^6)
    l6i = pow(l2i, 3.0);
    // forces calculation (direction later...)
    force = 48.0 * epsilon * l2i * l6i * sigma6 * (sigma6 * l6i - 0.5);
    // calculate potential energy...
    Elj = 4.0 * epsilon * l6i * sigma6 * (l6i * sigma6 - 1.0);

    // now, calculate forces
    for (i = 0; i < 3; i++)
    {
        parent->atoms[atomJ].force[i] += force * l[i];
        parent->atoms[atomI].force[i] -= force * l[i];
    }

    // update virial
    parent->virial[type] -= force * l2;
    // update potential energy
    parent->Epot[type] += Elj;
    // return potential energy
    return Elj;
}

// Calculate potential energy caused by this bond
double LJBond::CalculateEpot() const
{
    double l[3];
    double l2, l2i, l6i, Elj;
    int i;

    // get current distance vector
    for (i = 0; i < 3; i++)
    {
        l[i] = parent->atoms[atomJ].GetPosition(i) - parent->atoms[atomI].GetPosition(i);
    }

    // calculate square of the atom-atom distance: its norm not needed
    l2 = l[0] * l[0] + l[1] * l[1] + l[2] * l[2];
    // inverse of the square atom distance
    l2i = 1 / l2;
    // 1/(l^6)
    l6i = pow(l2i, 3.0);
    // calculate potential energy...
    Elj = 4.0 * epsilon * l6i * sigma6 * (l6i * sigma6 - 1.0);

    // return potential energy
    return Elj;
}

// Get one of connected atoms
int LJBond::GetAtom(int index) const
{
    if (index == 1)
    {
        return atomI;
    }
    else if (index == 2)
    {
        return atomJ;
    }
    else
    {
        return -1;
    }
}

// Add info to table row for .prt printing
int LJBond::GetParams(std::string &fieldname, std::vector<double> &params) const
{
    params.clear();
    params.push_back(epsilon);
    params.push_back(sigma);
    fieldname = "lj";
    return params.size();    
}