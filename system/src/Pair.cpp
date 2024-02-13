/*
 *  Atomic pair for `simul++`
 *
 *  Author: JJ
 *  Date: Dec 2022
 *
 */

#include <cmath>
#include "Pair.hpp"
#include "Atom.hpp"

Pair::Pair()
{
    int i;

    atom1 = nullptr;
    atom2 = nullptr;
    sqdist = 0.0;
    for (i = 0; i < 3; i++)
    {
        distvec[i] = 0.0;
    }
}

Pair::~Pair()
{
}

Pair::Pair(const Pair &origpair)
{
    int i;
    atom1 = origpair.atom1;
    atom2 = origpair.atom2;
    sqdist = origpair.sqdist;

    for (i = 0; i < 3; i++)
    {
        distvec[i] = origpair.distvec[i];
    }
}

Pair::Pair(Atom *at1, Atom *at2)
{
    atom1 = at1;
    atom2 = at2;
    CalculateDistance();
}

Pair::Pair(Atom *at1, Atom *at2, double *box)
{
    atom1 = at1;
    atom2 = at2;
    CalculateDistance(box);
}

Pair &Pair::operator=(const Pair &origpair)
{
    int i;
    atom1 = origpair.atom1;
    atom2 = origpair.atom2;
    sqdist = origpair.sqdist;

    for (i = 0; i < 3; i++)
    {
        distvec[i] = origpair.distvec[i];
    }

    return *this;
}

double Pair::CalculateDistance()
{
    int m;
// calculate distance vector r_kl
#pragma GCC unroll 3
    for (m = 0; m < 3; m++)
    {
        distvec[m] = atom2->GetPosition(m) - atom1->GetPosition(m);
    }
    // calculate square of the atom-atom distance
    sqdist = distvec[0] * distvec[0] + distvec[1] * distvec[1] + distvec[2] * distvec[2];
    return sqdist;
}

double Pair::CalculateDistance(double *box)
{
    int m;
// calculate distance vector r_kl
#pragma GCC unroll 3
    for (m = 0; m < 3; m++)
    {
        distvec[m] = atom2->GetPosition(m) - atom1->GetPosition(m);
        distvec[m] = distvec[m] - box[m] * round(distvec[m] / box[m]);
    }
    // calculate square of the atom-atom distance
    sqdist = distvec[0] * distvec[0] + distvec[1] * distvec[1] + distvec[2] * distvec[2];
    return sqdist;
}

double Pair::CalculateDistance(double sqcutoff)
{
    int m;
    double r[3];
// calculate distance vector r_kl
#pragma GCC unroll 3
    for (m = 0; m < 3; m++)
    {
        r[m] = atom2->GetPosition(m) - atom1->GetPosition(m);
    }
    // calculate square of the atom-atom distance
    sqdist = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
    if (sqdist < sqcutoff)
    {
#pragma GCC unroll 3
        for (m = 0; m < 3; m++)
        {
            distvec[m] = r[m];
        }
    }
    return sqdist;
}

double Pair::CalculateDistance(double sqcutoff, double *box)
{
    int m;
    double r[3];
// calculate distance vector r_kl
#pragma GCC unroll 3
    for (m = 0; m < 3; m++)
    {
        r[m] = atom2->GetPosition(m) - atom1->GetPosition(m);
        r[m] = r[m] - box[m] * round(r[m] / box[m]);
    }
    // calculate square of the atom-atom distance
    sqdist = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
    if (sqdist < sqcutoff)
    {
#pragma GCC unroll 3
        for (m = 0; m < 3; m++)
        {
            distvec[m] = r[m];
        }
    }
    return sqdist;
}

void Pair::SetAtoms(Atom *at1, Atom *at2)
{
    atom1 = at1;
    atom2 = at2;
}

void Pair::SetDistance(double sqdistance)
{
    sqdist = sqdistance;
}

void Pair::SetDistance(double sqdistance, double *distvector)
{
    int i;
    sqdist = sqdistance;

    for (i = 0; i < 3; i++)
    {
        distvec[i] = distvector[i];
    }
}
