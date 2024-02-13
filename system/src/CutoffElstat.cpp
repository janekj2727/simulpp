/*
 * A class to store cutoff electrostatic interactions
 * Author JJ, Date Jul 2022
 */

#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <cstring>
#include <cmath>

#include "Vector.hpp"
#include "Matrix.hpp"
#include "NewtonMethod.hpp"
#include "SimulatedSystem.hpp"
#include "CutoffElstat.hpp"
#include "ElstatBond.hpp"
#include "general_utils.hpp"
#include "AbstractPairList.hpp"

// #define DEBUG

// constructor (allocate matrix of correct sizes)
CutoffElstat::CutoffElstat(SimulatedSystem *parent, double cut, double alp, unsigned int shift)
    : AbstractInterMolField(0, parent)
{
    int i;
    double total_charge = 0.0;

    cutoff = cut;
    cutoff2 = cut * cut;
    alpha = alp;
    Ecutoffcorr = 0.0;
    sft = shift;
    if (sft != 3)
    {
        print_warning(1, "Cutoff electrostatics selected, shifting and/or tail smoothing disabled (sft = " + std::to_string(sft) + ").\n",
        "Continue only if you know what you are doing...\n");
    }
    CalculateABshift();

    // default do not calculate electrostatic forces (until an Atom with charge appears)
    empty = true;

    // check if free boundary conditions or neutral system... + check if any charges
    for (i = 0; i < parent->noMolecules; i++)
    {
        if (parent->molecules[i].HasAnyCharged())
        {
            empty = false;
        }
    }
    if (!empty)
    {
        for (i = 0; i < parent->noMolecules; i++)
        {
            total_charge += parent->molecules[i].TotalCharge();
        }
    }
    if ((parent->boundaryCond != 0) && (fabs(total_charge) > 5.0e-10))
    {
        print_warning(0, "Total charge is not zero (", std::to_string(total_charge), "), simulation in PBC cannot be done!\n");
        parent->error = 49;
    }
}

// copy constructor
CutoffElstat::CutoffElstat(const CutoffElstat &system)
    : AbstractInterMolField(0, system.parent)
{
    // LJmap = system.LJmap;

    cutoff = system.cutoff;
    cutoff2 = cutoff * cutoff;
    empty = system.empty;
    A = system.A;
    B = system.B;
    shift = system.shift;
    alpha = system.alpha;
    Ecutoffcorr = system.Ecutoffcorr;
}

// delete matrices of coefficients
CutoffElstat::~CutoffElstat()
{
    // empty
}

// calculate electrostatic forces
double CutoffElstat::CalculateForces()
{
    double Eel = 0.0;
    double virial = 0.0;

    double alphacut = alpha * cutoff;
    double A3mBcut = 3.0 * A - B * cutoff;
    // double B2 = 2.0 * B;
    double B4 = 4.0 * B;

    if (empty)
    {
        return 0.0;
    }

    /* forces evaluation - electrostatic interactions */
    parent->pairlist->Reset();

#ifdef PARALLEL
#pragma omp parallel num_threads(thread_count) reduction(+                  \
                                                         : Eel) reduction(+ \
                                                                          : virial)
    {
#endif
        int m;
        Vector r(3);
        double r2, rnorm;
        double force;
        double q2;
        Atom *atom1, *atom2;
        Pair pair;
        while ((r2 = parent->pairlist->GetNextPair(pair, cutoff2)) > 0.0)
        {
            atom1 = pair.GetAtom1();
            atom2 = pair.GetAtom2();
            // get charges, if zero, then this couple can be omitted
            q2 = atom1->charge * atom2->charge;
            if (q2 == 0.0)
            {
                continue;
            }

            rnorm = sqrt(r2);
            // if |r| < alpha*cutoff then shifted Coulomb interaction
            if ((rnorm < alphacut) || (parent->boundaryCond == 0))
            {
                force = q2 / r2 / rnorm;
                // calculate potential energy...
                Eel += q2 * (1.0 / rnorm - shift);
            }
            // else MACSIMUS style cutoff intermediate (see MACSIMUS manual) (if periodic bound. cond.)
            else
            {
                force = -q2 * pow(rnorm - cutoff, 2.0) * (A3mBcut / rnorm + B4);
                // calculate potential energy
                Eel += q2 * (pow(rnorm - cutoff, 3.0) * (A + B * rnorm));
            }
            // now, calculate forces
            for (m = 0; m < 3; m++)
            {
#ifdef PARALLEL
#pragma omp atomic
#endif
                atom2->force[m] += force * pair.GetVector(m);
#ifdef PARALLEL
#pragma omp atomic
#endif
                atom1->force[m] -= force * pair.GetVector(m);
            }
            // calculate virial
            virial -= force * r2;
        }
#ifdef PARALLEL
    }
#endif
    parent->virintermol[0] += virial;
    parent->Epotintermol[0] += Eel;
    return Eel;
}

// calculate electrostatic potential energy
double CutoffElstat::CalculateEpot() const
{
    double Eel = 0.0;
    double alphacut = alpha * cutoff;

    if (empty)
    {
        return 0.0;
    }

    /* forces evaluation - intermolecular interactions */
    parent->pairlist->Reset();

#ifdef PARALLEL
#pragma omp parallel num_threads(thread_count) reduction(+ \
                                                         : Eel)
    {
#endif
        Vector r(3);
        double r2, rnorm;
        double q2;
        Atom *atom1, *atom2;
        Pair pair;

        while ((r2 = parent->pairlist->GetNextPair(pair, cutoff2)) > 0.0)
        {
            atom1 = pair.GetAtom1();
            atom2 = pair.GetAtom2();
            // get charges, if zero, then this couple can be omitted
            q2 = atom1->charge * atom2->charge;
            if (q2 == 0.0)
            {
                continue;
            }
            rnorm = sqrt(r2);
            // if |r| < alpha*cutoff then shifted Coulomb interaction
            if ((rnorm < alphacut) || (parent->boundaryCond == 0))
            {
                // calculate potential energy...
                Eel += q2 * (1.0 / rnorm - shift);
            }
            // else MACSIMUS style cutoff intermediate (see MACSIMUS manual) (if periodic bound. cond.)
            else
            {
                // calculate potential energy
                Eel += q2 * (pow(rnorm - cutoff, 3.0) * (A + B * rnorm));
            }
        }
#ifdef PARALLEL
    }
#endif

    return Eel;
}

// calculate C1 and A for each pair
void CutoffElstat::CalculateABshift()
{
    switch (sft)
    {
    case 0: // no shift
        alpha = 1.0;
        shift = 0.0;
        A = 0.0;
        B = 0.0;
        break;
    case 1: // energy not shifted, but truncation done so that forces are continuous (weird)
        shift = 0.0;
        A = (5 * alpha - 2.0) / (cutoff2 * cutoff2 * alpha * pow(alpha * alpha - 2.0 * alpha + 1.0, 2.0));
        B = -(4.0 * alpha - 1) / (alpha * alpha * cutoff2 * cutoff2 * cutoff * (pow(alpha, 4.0) - 4.0 * pow(alpha, 3.0) + 6.0 * alpha * alpha - 4.0 * alpha + 1.0));
        break;
    case 2: // energy shifted but no intermediate function
        alpha = 1.0;
        shift = 1.0 / cutoff;
        A = 0.0;
        B = 0.0;
        break;
    default: // energy shifted and tail smoothed (the original version)
        shift = (10.0 * alpha * alpha - 5.0 * alpha + 1.0) / (6.0 * pow(alpha, 3.0) * cutoff);
        A = -(10.0 * alpha * alpha - 8.0 * alpha + 1.0) / (6.0 * pow(alpha, 3.0) * pow(cutoff, 4.0) * (alpha - 1.0) * (alpha * alpha - 2.0 * alpha + 1.0));
        B = (2.0 * alpha - 1.0) / (2.0 * pow(alpha, 3.0) * pow(cutoff, 5.0) * (pow(alpha, 3.0) - 3.0 * alpha * alpha + 3.0 * alpha - 1.0));
        break;
    }
}

// return energy correction (times volume)
double CutoffElstat::CalculateCutoffCorrection()
{
    Ecutoffcorr = 0.0;
    return 0.0;
}

// set electrostatic interactions cutoff
int CutoffElstat::SetCutoff(double vdwcut, double cut, int boundaryC)
{
    if (boundaryC == 0)
    {
        cutoff = 1000.0; // effectively infinity
        cutoff2 = 1.0e6;
        shift = 0.0; // if free boundary conditions, don't shift the potential to get the original Coulomb
        A = 0.0;
        B = 0.0; // A and B not used if free b.c.
    }
    else
    {
        cutoff = cut;
        cutoff2 = cut * cut;
        CalculateABshift();
    }

    return 0;
}

// cloning
CutoffElstat *CutoffElstat::copy(SimulatedSystem *newparent) const
{
    CutoffElstat *newCutoffElstat = new CutoffElstat(*this);
    newCutoffElstat->SetParent(newparent);
    return newCutoffElstat;
}

// intramol electrostatic interactions
// initialization of intramol bonds
int CutoffElstat::InitializeIntramolBonds(double LJ14factor, double el14factor) const
{
    int noInter = 0;
    int i, j, k;
    Molecule *mol;

    for (i = 0; i < parent->noMolecules; i++)
    {
        mol = &(parent->molecules[i]);
        for (j = 0; j < mol->noAtoms; j++)
        {
            for (k = 0; k < j; k++)
            {
                if (mol->GetDistance(j, k) == 3)
                {
                    mol->intraMolFields.push_back(new ElstatBond(mol, j, k, shift, el14factor));
                    noInter++;
                }
                else if (mol->GetDistance(j, k) > 3)
                {
                    mol->intraMolFields.push_back(new ElstatBond(mol, j, k, shift));
                    noInter++;
                }
            }
        }
    }

    return noInter;
}

// print info about the system to output stream (.prt file)
void CutoffElstat::PrintInfo(std::ofstream &stream, double u_eng, std::string engunit) const
{
    stream << "### Cut & shift electrostatics\n\n";
    stream << "The pair contribution to potential energy due to the cut & shift electrostatics is\n\n";
    stream << "$E_{\\mathrm{cutel}} = \\frac{q_i q_j}{4\\pi \\epsilon_0}\n";
    stream << "    \\begin{cases}\n";
    stream << "        \\left(\\frac{1}{r_{ij}} - \\mathtt{shift}\\right)                 & r_{ij} < \\alpha \\mathtt{cutoff} \\\\\n";
    stream << "        \\left[\\left(r_{ij}-\\mathtt{cutoff}\\right)^3(A+Br_{ij})\\right] & \\alpha \\mathtt{cutoff} \\leq r_{ij} < \\mathtt{cutoff} \\\\\n";
    stream << "        0                                                                  & r_{ij} \\leq \\mathtt{cutoff} \\\\\n";
    stream << "     \\end{cases}$\n\n";
    stream << "where charges $q_i$ and $q_j$ are charges of atoms (can be found above in atom tables in the molecules section),\n";
    stream << "cutoff parameters have this values: $\\alpha = " << alpha << "$, $\\mathtt{cuttof} = " << cutoff << "$ [AA],\n";
    stream << "$A = " << A << "$, $B = " << B << "$ and $\\mathtt{shift} = " << shift << "$ [AA^-1].\n\n";
}
