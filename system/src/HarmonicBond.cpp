/*
 * A class to store bond info (simul++)
 * Author JJ, Date Jan 2021
 * version 2 (Jun 2022): restructuralization using AbstractIntraMolField
 */

#include <iostream>
#include <cassert>

#include "Vector.hpp"
#include "Matrix.hpp"
#include "Atom.hpp"
#include "Molecule.hpp"
#include "HarmonicBond.hpp"
#include "MDTable.hpp"

// Default constructor
HarmonicBond::HarmonicBond(Molecule *parentMol, int indexI, int indexJ, double bondLength, double bondStrength)
    : AbstractIntraMolField{0, parentMol}
{
    length = bondLength;
    k = bondStrength;
    atomI = indexI;
    atomJ = indexJ;
}

// Copy constructor
HarmonicBond::HarmonicBond(const HarmonicBond &otherBond)
    : AbstractIntraMolField{otherBond}
{
    length = otherBond.length;
    k = otherBond.k;
    atomI = otherBond.atomI;
    atomJ = otherBond.atomJ;
}

// Overloading operator =
HarmonicBond &HarmonicBond::operator=(const HarmonicBond &otherBond)
{
    parent = NULL;
    type = otherBond.type;
    length = otherBond.length;
    k = otherBond.k;
    atomI = otherBond.atomI;
    atomJ = otherBond.atomJ;
    return *this;
}

// cloning
HarmonicBond* HarmonicBond::copy(Molecule *newparent) const
{
    HarmonicBond* newBond = new HarmonicBond(*this);
    newBond->SetParent(newparent);
    return newBond;
}

// Calculate current bond length
double HarmonicBond::GetCurrentLength() const
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
double HarmonicBond::CalculateForces()
{
    Vector l(3);
    double delta;
    double norm;
    double Epot;
    int i;

    for (i = 0; i < 3; i++)
    {
        l[i] = parent->atoms[atomJ].GetPosition(i) - parent->atoms[atomI].GetPosition(i);
    }

    norm = l.CalculateNorm();
    delta = norm - length;
    l = l * (1 / norm);

    for (i = 0; i < 3; i++)
    {
        parent->atoms[atomI].force[i] += k * delta * l[i];
        parent->atoms[atomJ].force[i] -= k * delta * l[i];
    }

    parent->virial[type] += k * delta * norm;

    Epot = 0.5 * k * delta * delta;
    parent->Epot[type] += Epot;

    return Epot;
}

// Calculate potential energy caused by this bond
double HarmonicBond::CalculateEpot() const
{
    Vector l(3);
    double delta;
    double norm;
    int i;

    for (i = 0; i < 3; i++)
    {
        l[i] = parent->atoms[atomJ].GetPosition(i) - parent->atoms[atomI].GetPosition(i);
    }

    norm = l.CalculateNorm();
    delta = norm - length;

    return 0.5 * k * delta * delta;
}

// Get one of connected atoms
int HarmonicBond::GetAtom(int index) const
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
int HarmonicBond::GetParams(std::string &fieldname, std::vector<double> &params) const
{
    params.clear();
    params.push_back(k);
    params.push_back(length);
    fieldname = "harm";
    return params.size();    
}