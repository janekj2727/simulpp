/*
 * A class to store bond info (simul++)
 * Author JJ, Date Jul 2022
 * version 1
 */

#include <iostream>
#include <cassert>
#include <cmath>

#include "Vector.hpp"
#include "Matrix.hpp"
#include "Atom.hpp"
#include "Molecule.hpp"
#include "ElstatBond.hpp"

// Default constructor
ElstatBond::ElstatBond(Molecule *parentMol, int indexI, int indexJ, double sft, double scaling)
    : AbstractIntraMolField{4, parentMol}
{
    qproduct = scaling * parentMol->atoms[indexI].charge * parentMol->atoms[indexJ].charge;
    shift = sft;
    atomI = indexI;
    atomJ = indexJ;
}

// Copy constructor
ElstatBond::ElstatBond(const ElstatBond &otherBond)
    : AbstractIntraMolField{otherBond}
{
    qproduct = otherBond.qproduct;
    atomI = otherBond.atomI;
    atomJ = otherBond.atomJ;
}

// Overloading operator =
ElstatBond &ElstatBond::operator=(const ElstatBond &otherBond)
{
    parent = NULL;
    type = otherBond.type;
    qproduct = otherBond.qproduct;
    atomI = otherBond.atomI;
    atomJ = otherBond.atomJ;
    return *this;
}

// cloning
ElstatBond *ElstatBond::copy(Molecule *newparent) const
{
    ElstatBond *newBond = new ElstatBond(*this);
    newBond->SetParent(newparent);
    return newBond;
}

// Calculate current bond length
double ElstatBond::GetCurrentLength() const
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
double ElstatBond::CalculateForces()
{
    double l[3];
    double l2, lnorm;
    double Eelstat, force;
    int i;

    // get current distance vector
    for (i = 0; i < 3; i++)
    {
        l[i] = parent->atoms[atomJ].GetPosition(i) - parent->atoms[atomI].GetPosition(i);
    }

    // calculate square of the atom-atom distance: its norm not needed
    l2 = l[0] * l[0] + l[1] * l[1] + l[2] * l[2];
    lnorm = sqrt(l2);
    
    // forces calculation (direction later...)
    force = qproduct / (l2 * lnorm);
    // calculate potential energy...
    Eelstat = qproduct * (1.0 / lnorm - shift);

    // now, calculate forces
    for (i = 0; i < 3; i++)
    {
        parent->atoms[atomJ].force[i] += force * l[i];
        parent->atoms[atomI].force[i] -= force * l[i];
    }

    // update virial
    parent->virial[type] -= force * l2;
    // update potential energy
    parent->Epot[type] += Eelstat;
    // return potential energy
    return Eelstat;
}

// Calculate potential energy caused by this bond
double ElstatBond::CalculateEpot() const
{
    double l[3];
    double lnorm, Eelstat;
    int i;

    // get current distance vector
    for (i = 0; i < 3; i++)
    {
        l[i] = parent->atoms[atomJ].GetPosition(i) - parent->atoms[atomI].GetPosition(i);
    }

    // calculate square of the atom-atom distance: its norm not needed
    lnorm = sqrt(l[0] * l[0] + l[1] * l[1] + l[2] * l[2]);
    
    // calculate potential energy...
    Eelstat = qproduct * (1.0 / lnorm - shift);

    // return potential energy
    return Eelstat;
}

// Get one of connected atoms
int ElstatBond::GetAtom(int index) const
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
int ElstatBond::GetParams(std::string &fieldname, std::vector<double> &params) const
{
    params.clear();
    params.push_back(qproduct);
    params.push_back(shift);
    fieldname = "elstat"; // maybe not needed
    return params.size();    
}