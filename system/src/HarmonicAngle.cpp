/*
 * A class to store harmonic angle (simul++)
 * Author JJ, Date Jun 2022
 */

#include <iostream>
#include <cmath>

#include "Vector.hpp"
#include "Matrix.hpp"
#include "Atom.hpp"
#include "Molecule.hpp"
#include "HarmonicAngle.hpp"

// Default constructor
HarmonicAngle::HarmonicAngle(Molecule *parentMol, int indexI, int indexJ, int indexK, double equilAngle, double angleStrength)
    : AbstractIntraMolField{1, parentMol}
{
    angle = equilAngle;
    k = angleStrength;
    atomI = indexI;
    atomJ = indexJ;
    atomK = indexK;
}

// Copy constructor
HarmonicAngle::HarmonicAngle(const HarmonicAngle &otherAngle)
    : AbstractIntraMolField{otherAngle}
{
    angle = otherAngle.angle;
    k = otherAngle.k;
    atomI = otherAngle.atomI;
    atomJ = otherAngle.atomJ;
    atomK = otherAngle.atomK;
}

// Overloading operator =
HarmonicAngle &HarmonicAngle::operator=(const HarmonicAngle &otherAngle)
{
    parent = NULL;
    type = otherAngle.type;
    angle = otherAngle.angle;
    k = otherAngle.k;
    atomI = otherAngle.atomI;
    atomJ = otherAngle.atomJ;
    atomK = otherAngle.atomK;
    return *this;
}

// cloning
HarmonicAngle *HarmonicAngle::copy(Molecule *newparent) const
{
    HarmonicAngle *newAngle = new HarmonicAngle(*this);
    newAngle->SetParent(newparent);
    return newAngle;
}

// Calculate current bond length
double HarmonicAngle::GetCurrentAngle() const
{
    Vector r21(3), r23(3);
    int i;

    for (i = 0; i < 3; i++)
    {
        r21[i] = parent->atoms[atomI].GetPosition(i) - parent->atoms[atomJ].GetPosition(i);
        r23[i] = parent->atoms[atomK].GetPosition(i) - parent->atoms[atomJ].GetPosition(i);
    }

    return CalculateAngle(r21, r23);
}

// Calculate forces caused by this angle
// DL POLY 4 manual (pp 17–19), source code: angles_forces.f90 
double HarmonicAngle::CalculateForces()
{
    Vector r21(3), r23(3);
    double delta, Epot;
    double r21_norm, r23_norm, rr21_norm, rr23_norm;
    double prefactor, cosinus, sinus, instangle;
    double fI, fK;
    int i;

    for (i = 0; i < 3; i++)
    {
        r21[i] = parent->atoms[atomI].GetPosition(i) - parent->atoms[atomJ].GetPosition(i);
        r23[i] = parent->atoms[atomK].GetPosition(i) - parent->atoms[atomJ].GetPosition(i);
    }

    r21_norm = r21.CalculateNorm();
    r23_norm = r23.CalculateNorm();
    rr21_norm = 1.0/r21_norm;
    rr23_norm = 1.0/r23_norm;
    r21 = r21 * rr21_norm;
    r23 = r23 * rr23_norm;


    cosinus = CalculateScalarProduct(r21, r23);
    sinus = sqrt(1 - cosinus*cosinus);
    instangle = acos(cosinus);
    delta = instangle - angle;
    prefactor = k * delta / sinus;
    
    
    for (i = 0; i < 3; i++)
    {
        fI = prefactor * rr21_norm * (r23(i) - r21(i)*cosinus);
        fK = prefactor * rr23_norm * (r21(i) - r23(i)*cosinus);
        parent->atoms[atomI].force[i] += fI;
        parent->atoms[atomJ].force[i] += -fI - fK;
        parent->atoms[atomK].force[i] += fK;
    }

    parent->virial[type] += 0; // angles does not contribute to virial (pressure)
    Epot = 0.5 * k * delta * delta;
    parent->Epot[type] += Epot;

    return Epot;
}

// Calculate potential energy caused by this angle
double HarmonicAngle::CalculateEpot() const
{
    double diff = GetCurrentAngle() - angle;
    // Epot = k/2 * (theta - theta0)^2
    return 0.5 * k * diff * diff;
}

// Get one of connected atoms
int HarmonicAngle::GetAtom(int index) const
{
    if (index == 1)
    {
        return atomI;
    }
    else if (index == 2)
    {
        return atomJ;
    }
    else if (index == 3)
    {
        return atomK;
    }
    else
    {
        return -1;
    }
}

// Add info to table row for .prt printing
int HarmonicAngle::GetParams(std::string &fieldname, std::vector<double> &params) const
{
    params.clear();
    params.push_back(k);
    params.push_back(angle);
    fieldname = "harm";
    return params.size();    
}