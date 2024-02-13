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
#include "CosineDihedral.hpp"

// Default constructor
CosineDihedral::CosineDihedral(Molecule *parentMol, int indexI, int indexJ, int indexK, int indexL, double Astrength, double d, int mult, bool proper)
    : AbstractIntraMolField{(proper?2:5), parentMol}
{
    A = Astrength;
    delta = d;
    m = mult;
    atomI = indexI;
    atomJ = indexJ;
    atomK = indexK;
    atomL = indexL;
}

// Copy constructor
CosineDihedral::CosineDihedral(const CosineDihedral &otherDihedral)
    : AbstractIntraMolField{otherDihedral}
{
    A = otherDihedral.A;
    delta = otherDihedral.delta;
    m = otherDihedral.m;
    atomI = otherDihedral.atomI;
    atomJ = otherDihedral.atomJ;
    atomK = otherDihedral.atomK;
    atomL = otherDihedral.atomL;
}

// Overloading operator =
CosineDihedral &CosineDihedral::operator=(const CosineDihedral &otherDihedral)
{
    parent = NULL;
    type = otherDihedral.type;
    delta = otherDihedral.delta;
    m = otherDihedral.m;
    atomI = otherDihedral.atomI;
    atomJ = otherDihedral.atomJ;
    atomK = otherDihedral.atomK;
    atomL = otherDihedral.atomL;
    return *this;
}

// cloning
CosineDihedral *CosineDihedral::copy(Molecule *newparent) const
{
    CosineDihedral *newDihedral = new CosineDihedral(*this);
    newDihedral->SetParent(newparent);
    return newDihedral;
}

// Calculate current bond length
double CosineDihedral::GetCurrentAngle() const
{
    Vector rij(3), rjk(3), rkl(3);
    Vector rb(3), rc(3);
    double cosinus, sinus, rrb, rrc, rrjk;
    int i;

    for (i = 0; i < 3; i++)
    {
        rij[i] = parent->atoms[atomJ].GetPosition(i) - parent->atoms[atomI].GetPosition(i);
        rjk[i] = parent->atoms[atomK].GetPosition(i) - parent->atoms[atomJ].GetPosition(i);
        rkl[i] = parent->atoms[atomL].GetPosition(i) - parent->atoms[atomK].GetPosition(i);
    }
    rb = CalculateVectorProduct(rij, rjk);
    rc = CalculateVectorProduct(rjk, rkl);
    rrb = 1.0/rb.CalculateNorm();
    rrc = 1.0/rc.CalculateNorm();
    rrjk = 1.0/rjk.CalculateNorm();

    // cosinus of the dihedral angle
    cosinus = CalculateScalarProduct(rb, rc) * rrb * rrc;
    // sinus of the dihedral angle (with right sign)
    sinus = CalculateScalarProduct(rjk, CalculateVectorProduct(rb, rc)) * rrb * rrc * rrjk;
    // dihedral angle
    return atan2(sinus, cosinus);
}

// Calculate forces caused by this angle
// DL POLY 4 manual (pp 17–19), source code: angles_forces.f90
double CosineDihedral::CalculateForces()
{
    Vector rij(3), rjk(3), rkl(3);
    Vector rb(3), rc(3);
    double Epot;
    double rb_sq, rc_sq, rrb2, rrc2, rrb, rrc, rb_norm, rc_norm, rbc, rrjk;
    double prefactor, cosinus, sinus, phi, term; // prefactor for forces calculation, term = term inside cos in E_pot
    double fI[3], fJ[3], fK[3], fL[3];
    int i;

    // position vectors (NOTE: in DL POLY defined with opposite sign (r_ij = r_i - r_j)!!!)
    for (i = 0; i < 3; i++)
    {
        rij[i] = parent->atoms[atomJ].GetPosition(i) - parent->atoms[atomI].GetPosition(i);
        rjk[i] = parent->atoms[atomK].GetPosition(i) - parent->atoms[atomJ].GetPosition(i);
        rkl[i] = parent->atoms[atomL].GetPosition(i) - parent->atoms[atomK].GetPosition(i);
    }
    // auxiliary vectors b and c (cross-products)
    rb = CalculateVectorProduct(rij, rjk);
    rc = CalculateVectorProduct(rjk, rkl);
    // rbc = rb.rc
    rbc = CalculateScalarProduct(rb, rc);

    // auxiliary numbers (square norm, reciprocal square and reciprocal norm)
    rb_sq = CalculateScalarProduct(rb, rb);
    rc_sq = CalculateScalarProduct(rc, rc);
    rrb2 = 1.0 / rb_sq;
    rrc2 = 1.0 / rc_sq;
    rb_norm = sqrt(rb_sq);
    rc_norm = sqrt(rc_sq);
    rrb = 1.0 / rb_norm;
    rrc = 1.0 / rc_norm;

    // rrjk = 1/norm(rjk)
    rrjk = 1.0 / rjk.CalculateNorm();

    // cosinus of the dihedral angle
    cosinus = rbc * rrb * rrc;
    // sinus of the dihedral angle (with right sign)
    sinus = CalculateScalarProduct(rjk, CalculateVectorProduct(rb, rc)) * rrjk * rrb * rrc;
    // dihedral angle
    phi = atan2(sinus, cosinus);
    // term inside cos in E_pot
    term = phi * m - delta;
    // prefactor = - diff(E_pot, phi)*diff(phi, B)
    prefactor = -A * m / sinus * sin(term) * rrb * rrc;

    // diff(B, r_a^i)
    // B * ||r_c|| / ||r_b|| = (r_b . r_c) * ||r_c|| / (||r_c|| * ||r_b||^2)
    // fI corresponds to fax, fay, faz in DL POLY code
    fI[0] = (rc(1) * rjk(2) - rc(2) * rjk(1)) - rbc * rrb2 * (rb(1) * rjk(2) - rb(2) * rjk(1));
    fI[1] = (rc(2) * rjk(0) - rc(0) * rjk(2)) - rbc * rrb2 * (rb(2) * rjk(0) - rb(0) * rjk(2));
    fI[2] = (rc(0) * rjk(1) - rc(1) * rjk(0)) - rbc * rrb2 * (rb(0) * rjk(1) - rb(1) * rjk(0));

    // fJ corresponds to fb1x, fb1y, fb1z in DL POLY code
    fJ[0] = (rb(1) * rkl(2) - rb(2) * rkl(1)) - rbc * rrc2 * (rc(1) * rkl(2) - rc(2) * rkl(1));
    fJ[1] = (rb(2) * rkl(0) - rb(0) * rkl(2)) - rbc * rrc2 * (rc(2) * rkl(0) - rc(0) * rkl(2));
    fJ[2] = (rb(0) * rkl(1) - rb(1) * rkl(0)) - rbc * rrc2 * (rc(0) * rkl(1) - rc(1) * rkl(0));

    // fK corresponds to fcx, fcy, fcz in DL POLY code
    fK[0] = (rc(1) * rij(2) - rc(2) * rij(1)) - rbc * rrb2 * (rb(1) * rij(2) - rb(2) * rij(1));
    fK[1] = (rc(2) * rij(0) - rc(0) * rij(2)) - rbc * rrb2 * (rb(2) * rij(0) - rb(0) * rij(2));
    fK[2] = (rc(0) * rij(1) - rc(1) * rij(0)) - rbc * rrb2 * (rb(0) * rij(1) - rb(1) * rij(0));

    // fL corresponds to fd1x, fd1y, fd1z in DL POLY code
    fL[0] = (rb(1) * rjk(2) - rb(2) * rjk(1)) - rbc * rrc2 * (rc(1) * rjk(2) - rc(2) * rjk(1));
    fL[1] = (rb(2) * rjk(0) - rb(0) * rjk(2)) - rbc * rrc2 * (rc(2) * rjk(0) - rc(0) * rjk(2));
    fL[2] = (rb(0) * rjk(1) - rb(1) * rjk(0)) - rbc * rrc2 * (rc(0) * rjk(1) - rc(1) * rjk(0));

    for (i = 0; i < 3; i++)
    {
        parent->atoms[atomI].force[i] += prefactor * (fI[i]);
        parent->atoms[atomJ].force[i] += prefactor * (fJ[i] - fI[i] - fK[i]);
        parent->atoms[atomK].force[i] += prefactor * (fK[i] - fJ[i] - fL[i]);
        parent->atoms[atomL].force[i] += prefactor * (fL[i]);
    }

    parent->virial[type] += 0; //dihedrals does not contribute to virial (pressure)
    Epot = A * (1.0 + cos(term));    
    parent->Epot[type] += Epot;

    return Epot;
}

// Calculate potential energy caused by this angle
double CosineDihedral::CalculateEpot() const
{
    double phi = GetCurrentAngle();
    // Epot = A[1 + cos(m*phi - delta)]
    return A * (1.0 + cos(m * phi - delta));
}

// Get one of connected atoms
int CosineDihedral::GetAtom(int index) const
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
    else if (index == 4)
    {
        return atomL;
    }
    else
    {
        return -1;
    }
}

// Add info to table row for .prt printing
int CosineDihedral::GetParams(std::string &fieldname, std::vector<double> &params) const
{
    params.clear();
    params.push_back(A);
    params.push_back(delta);
    params.push_back(m);
    fieldname = "cos";
    return params.size();    
}