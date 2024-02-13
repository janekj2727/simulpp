/*
 * A class to store parameters of LJ interaction
 * Author JJ, Date Feb 2021
 */

#include <iostream>
#include <map>
#include <string>
#include <cstring>
#include <cmath>
#include <fstream>

#include "Vector.hpp"
#include "Matrix.hpp"
#include "NewtonMethod.hpp"
#include "SimulatedSystem.hpp"
#include "LJsystem.hpp"
#include "LJBond.hpp"
#include "units.hpp"
#include "AbstractPairList.hpp"

#define DEFAULT_SIGMA 1
// #define DEBUG

// class for C1 and A calculation
class EquationSystem
{
private:
    double sigma;
    double epsilon;
    double sigma6;
    double cut2;

public:
    EquationSystem(){};
    ~EquationSystem(){};
    Vector Evalf(Vector &vec)
    {
        Vector res(2);
        res(0) = epsilon * (6 * sigma6 / pow(vec(0), 7.0) - 12 * sigma6 * sigma6 / pow(vec(0), 13.0)) -
                 vec(1) * (vec(0) * vec(0) - cut2) * vec(0);
        res(1) = 4.0 * epsilon * (sigma6 * sigma6 / pow(vec(0), 12.0) - sigma6 / pow(vec(0), 6.0)) -
                 vec(1) * (vec(0) * vec(0) - cut2) * (vec(0) * vec(0) - cut2);
        return res;
    };
    void SetEpsilonSigma(double eps, double sig) { epsilon = eps, sigma = sig, sigma6 = pow(sigma, 6.0); };
    void SetLJcut(double cutoff2) { cut2 = cutoff2; };
};

EquationSystem *eqsystem;

Vector eqsystemevalf(Vector &vec)
{
    Vector result(2);
    result = eqsystem->Evalf(vec);
    return result;
}

// constructor (allocate matrix of correct sizes)
LJsystem::LJsystem(SimulatedSystem *parent, int siz)
    : AbstractInterMolField(1, parent)
{
    int i;

    size = siz;

    sigmaMatrix = new Matrix(size, size);
    sigma6Matrix = new Matrix(size, size);

    // set default diameter to diagonal terms to avoid atoms of zero radius
    for (i = 0; i < size; i++)
    {
        sigmaMatrix->operator()(i, i) = DEFAULT_SIGMA;
    }

    epsilonMatrix = new Matrix(size, size);
    C1Matrix = new Matrix(size, size);
    AMatrix = new Matrix(size, size);

    // temporarily set cutoff
    cutoff2 = 10000.0;
    cutoff = 100.0;

    // default do not calculate LJ forces (until a pair with non-zero epsilon appears)
    empty = true;

    // zero cutoff correction energy
    Ecutoffcorr = 0.0;
}

// copy constructor
LJsystem::LJsystem(const LJsystem &system)
    : AbstractInterMolField(1, system.parent)
{
    // LJmap = system.LJmap;

    sigmaMatrix = new Matrix(*system.sigmaMatrix);
    sigma6Matrix = new Matrix(*system.sigma6Matrix);
    epsilonMatrix = new Matrix(*system.epsilonMatrix);
    C1Matrix = new Matrix(*system.C1Matrix);
    AMatrix = new Matrix(*system.AMatrix);

    cutoff = system.cutoff;
    cutoff2 = system.cutoff2;
    empty = system.empty;
    Ecutoffcorr = system.Ecutoffcorr;
}

// delete matrices of coefficients
LJsystem::~LJsystem()
{
    delete sigmaMatrix;
    delete sigma6Matrix;
    delete epsilonMatrix;
    delete C1Matrix;
    delete AMatrix;
}

// set sigma of a given pair
int LJsystem::SetSigma(int LJid1, int LJid2, double sigma)
{

    if ((sigma < 0) || (sigma > 1000000))
    {
        return 29;
    }
    sigmaMatrix->operator()(LJid1, LJid2) = sigma;
    sigma6Matrix->operator()(LJid1, LJid2) = pow(sigma, 6.0);
    if (LJid1 != LJid2)
    {
        sigmaMatrix->operator()(LJid2, LJid1) = sigma;
        sigma6Matrix->operator()(LJid2, LJid1) = pow(sigma, 6.0);
    }
    return 0;
}

// set epsilon of a given pair
int LJsystem::SetEpsilon(int LJid1, int LJid2, double epsilon)
{
    if ((epsilon < 0.0) || (epsilon > 1000000.0))
    {
        return 28;
    }

    if (epsilon > 0.0) // forces need to be calculated
    {
        empty = false;
    }

    epsilonMatrix->operator()(LJid1, LJid2) = epsilon;
    if (LJid1 != LJid2)
    {
        epsilonMatrix->operator()(LJid2, LJid1) = epsilon;
    }
    return 0;
}

// set parameters (epsilon and sigma) for the given pair
int LJsystem::SetPairParams(int LJid1, int LJid2, char *params)
{
    double epsilon, sigma;
    sscanf(params, "%lf %lf", &epsilon, &sigma);

    // epsilon
    epsilon /= parent->GetEnergyUnit(); // unit conversion from input to p.u.
    if ((epsilon < 0.0) || (epsilon > 1000000.0))
    {
        return 28;
    }

    if (epsilon > 0.0) // forces need to be calculated
    {
        empty = false;
    }

    epsilonMatrix->operator()(LJid1, LJid2) = epsilon;
    if (LJid1 != LJid2)
    {
        epsilonMatrix->operator()(LJid2, LJid1) = epsilon;
    }

    if ((sigma < 0) || (sigma > 1000000))
    {
        return 29;
    }
    sigmaMatrix->operator()(LJid1, LJid2) = sigma;
    sigma6Matrix->operator()(LJid1, LJid2) = pow(sigma, 6.0);
    if (LJid1 != LJid2)
    {
        sigmaMatrix->operator()(LJid2, LJid1) = sigma;
        sigma6Matrix->operator()(LJid2, LJid1) = pow(sigma, 6.0);
    }
    return 0;
}

// get size of LJsystem
int LJsystem::GetSystemSize() const
{
    return size;
}

// calculate LJ forces
double LJsystem::CalculateForces()
{
    double Elj = 0.0;
    double virial = 0.0;

    if (empty)
    {
        return 0.0;
    }
    parent->pairlist->Reset();

#ifdef PARALLEL
#pragma omp parallel num_threads(thread_count) reduction(+ : Elj) reduction(+ : virial)
    {
#endif
        int m;
        double r2, r2i, r6i;
        double epsilon, C1, A;
        double sigma6, r2mC2;
        double force;
        Atom *atom1, *atom2;
        Pair pair;
        while ((r2 = parent->pairlist->GetNextPair(pair, cutoff2)) > 0.0)
        {
            atom1 = pair.GetAtom1();
            atom2 = pair.GetAtom2();

            // get C1 for the couple
            C1 = GetC1(atom2->LJid, atom1->LJid);

            // if r2 < C1^2 then normal LJ interaction (see Frenkel, Smit)
            if ((r2 < C1 * C1) || (parent->boundaryCond == 0))
            {
                // get interaction parameters
                epsilon = GetEpsilon(atom2->LJid, atom1->LJid);
                sigma6 = GetSigma6(atom2->LJid, atom1->LJid);
                r2i = 1 / r2;
                r6i = pow(r2i, 3.0);
                // sigma6 = pow(sigma, 6.0);
                force = 48.0 * epsilon * r2i * r6i * sigma6 * (sigma6 * r6i - 0.5);
                // calculate potential energy...
                Elj += 4.0 * epsilon * r6i * sigma6 * (r6i * sigma6 - 1.0);
            }
            // else MACSIMUS style cutoff intermediate (see MACSIMUS manual) (if periodic bound. cond.)
            else
            {
                A = GetA(atom2->LJid, atom1->LJid);
                r2mC2 = r2 - cutoff2;
                force = -4 * A * r2mC2;
                // calculate potential energy
                Elj += A * r2mC2 * r2mC2;
            }
            // now, save forces

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
            virial += force * r2;
        }
#ifdef PARALLEL
    }
#endif

    parent->Epotintermol[1] += Elj;
    parent->virintermol[1] -= virial;
    return Elj;
}

// calculate LJ potential energy
double LJsystem::CalculateEpot() const
{
    double Elj = 0.0;

    if (empty)
    {
        return 0.0;
    }

    parent->pairlist->Reset();
#ifdef PARALLEL
#pragma omp parallel num_threads(thread_count) reduction(+ : Elj)
    {
#endif
        Vector r(3);
        double r2, r2i, r6i;
        double epsilon, C1, A;
        double sigma6, r2mC2;
        Atom *atom1, *atom2;
        Pair pair;
        while ((r2 = parent->pairlist->GetNextPair(pair, cutoff2)) > 0.0)
        {
            atom1 = pair.GetAtom1();
            atom2 = pair.GetAtom2();
            // get C1 for the couple
            C1 = GetC1(atom2->LJid, atom1->LJid);

            // if r2 < C1^2 then normal LJ interaction (see Frenkel, Smit)
            if ((r2 < C1 * C1) || (parent->boundaryCond == 0))
            {
                // get interaction parameters
                epsilon = GetEpsilon(atom2->LJid, atom1->LJid);
                sigma6 = GetSigma6(atom2->LJid, atom1->LJid);
                r2i = 1 / r2;
                r6i = pow(r2i, 3.0);
                // calculate potential energy...
                Elj += 4.0 * epsilon * r6i * sigma6 * (r6i * sigma6 - 1.0);
            }
            // else MACSIMUS style cutoff intermediate (see MACSIMUS manual) (if periodic bound. cond.)
            else
            {
                A = GetA(atom2->LJid, atom1->LJid);
                r2mC2 = r2 - cutoff2;
                // calculate potential energy
                Elj += A * r2mC2 * r2mC2;
            }
        }
#ifdef PARALLEL
    }
#endif

    return Elj;
}

// calculate mixed terms
int LJsystem::CalculateMixTerms(int mixingrule)
{
    int nopairscal = 0;
    int i, j;
    for (i = 0; i < size - 1; i++)
    {
        for (j = i + 1; j < size; j++)
        {
            if ((epsilonMatrix->operator()(i, j) == 0.0) && (sigmaMatrix->operator()(i, j) == 0.0))
            {
                switch (mixingrule)
                {
                case 1: // Lorentz–Berthelot (geom energy, aritm sigma)
                    if ((epsilonMatrix->operator()(i, i) > 0.0) && (epsilonMatrix->operator()(j, j) > 0.0))
                    {
                        SetEpsilon(i, j, sqrt(epsilonMatrix->operator()(i, i) * epsilonMatrix->operator()(j, j)));
                    }
                    if ((sigmaMatrix->operator()(i, i) > 0.0) && (sigmaMatrix->operator()(j, j) > 0.0))
                    {
                        // positive values of epsilon are assumed (checked during parameters reading)
                        SetSigma(i, j, 0.5 * (sigmaMatrix->operator()(i, i) + sigmaMatrix->operator()(j, j)));
                    }
                    break;
                default: // no mixing
                    break;
                }
                nopairscal++;
            }
        }
    }

    return nopairscal;
}

// calculate C1 and A for each pair
int LJsystem::CalculateC1A()
{
    int i, j;
    int error = 0;
    double epsilon, sigma, sigma6, cutoff6, x1, x2, C1, C16, A;

    cutoff6 = pow(cutoff, 6.0);
    for (i = 0; i < GetSystemSize(); i++)
    {
        for (j = i; j < GetSystemSize(); j++)
        {
            epsilon = epsilonMatrix->operator()(i, j);
            sigma = sigmaMatrix->operator()(i, j);
            if ((sigma * 1.7818 < cutoff) && (epsilon > 0.0)) // check if cutoff is 2^(5/6) * LJ sigma minimum (else incorrect results)
            {
                sigma6 = pow(sigma, 6.0);
                x1 = cbrt(27.0 * sigma6 * cutoff6 + 160.0 * sigma6 * sigma6 + sqrt(729.0 * cutoff6 * cutoff6 * sigma6 * sigma6 - 24128.0 * cutoff6 * sigma6 * sigma6 * sigma6 + 25600.0 * sigma6 * sigma6 * sigma6 * sigma6));
                x2 = sqrt((640.0 * cutoff2 * sigma6 + 9.0 * cutoff2 * cutoff2 * x1 + 20.0 * x1 * x1) / x1);
                C1 = sqrt((3.0 * cutoff2 + x2 + sqrt(-2.0 * (320.0 * x2 * cutoff2 * sigma6 - 27.0 * x1 * cutoff6 - 1600.0 * x1 * sigma6 - 9.0 * cutoff2 * cutoff2 * x1 * x2 + 10.0 * x2 * x1 * x1) / (x1 * x2))) / 20.0);
                C16 = pow(C1, 6.0);
                A = -4.0 * epsilon * sigma6 * (C16 - sigma6) / (C16 * C16 * (C1 * C1 - cutoff2) * (C1 * C1 - cutoff2));
                // assign results for both (i, j) and (j, i)
                C1Matrix->operator()(i, j) = C1;
                C1Matrix->operator()(j, i) = C1;
                AMatrix->operator()(i, j) = A;
                AMatrix->operator()(j, i) = A;
            }
            else
            {
                // hard cutoff
                C1Matrix->operator()(i, j) = cutoff;
                C1Matrix->operator()(j, i) = cutoff;
                AMatrix->operator()(i, j) = 0.0;
                AMatrix->operator()(j, i) = 0.0;
                if (epsilon > 0.0)
                    error++; // add to error counter
            }
        }
    }
    return error;
}

// return energy correction (times volume)
double LJsystem::CalculateCutoffCorrection()
{
    int multiplicity[GetSystemSize()] = {}; // number of each LJid
    double Ecorr = 0.0;                     // cutoff correction
    double A, C1, epsilon, sigma6;
    double cutoff4 = cutoff2 * cutoff2;
    double cutoff7 = cutoff4 * cutoff2 * cutoff;
    int i, j;

    /*
     *  projit vsechny pary a spocitat, kolik je kazdeho paru
     *  pro kazdy par spocitat korekci za predpokladu ciselne hustoty castic rovnomerne...
     *  melo by vyjit cislo, ktere se pak vydeli objemem a dostanu energii
     */

    for (i = 0; i < parent->noMolecules; i++)
    {
        for (j = 0; j < parent->molecules[i].noAtoms; j++)
        {
            multiplicity[parent->molecules[i].atoms[j].LJid]++;
        }
    }

    for (i = 0; i < GetSystemSize(); i++)
    {
        A = AMatrix->operator()(i, i);
        C1 = C1Matrix->operator()(i, i);
        epsilon = epsilonMatrix->operator()(i, i);
        sigma6 = sigma6Matrix->operator()(i, i);
        Ecorr += 2 * M_PI * multiplicity[i] * (multiplicity[i] - 1) * // multiplicity * (multiplicity - 1)  or multiplicity^2 ??? I think it is not obvious...
                 (A * (C1 * C1 * C1 * (cutoff4 / 3 - 2 * cutoff2 * C1 * C1 / 5 + C1 * C1 * C1 * C1 / 7) - 8 * cutoff7 / 105) +
                  4 * epsilon * sigma6 / (3.0 * C1 * C1 * C1) * (sigma6 / (3.0 * pow(C1, 6.0)) - 1));

        for (j = i + 1; j < GetSystemSize(); j++)
        {
            A = AMatrix->operator()(i, j);
            C1 = C1Matrix->operator()(i, j);
            epsilon = epsilonMatrix->operator()(i, j);
            sigma6 = sigma6Matrix->operator()(i, j);
            Ecorr += 2 * M_PI * multiplicity[i] * multiplicity[j] *
                     (A * (C1 * C1 * C1 * (cutoff4 / 3 - 2 * cutoff2 * C1 * C1 / 5 + C1 * C1 * C1 * C1 / 7) - 8 * cutoff7 / 105) +
                      4 * epsilon * sigma6 / (3.0 * C1 * C1 * C1) * (sigma6 / (3.0 * pow(C1, 6.0)) - 1));
        }
    }

    Ecutoffcorr = Ecorr;
    return Ecorr;
}

// set LJ interactions cutoff
int LJsystem::SetCutoff(double cut, double elcut, int boundaryC)
{
    int error = 0;
    cutoff = cut;
    cutoff2 = cutoff * cutoff;

    if (boundaryC > 0) // if not free b.c.
    {
        error = CalculateC1A();
    }
    return error;
}

// cloning
LJsystem *LJsystem::copy(SimulatedSystem *newparent) const
{
    LJsystem *newLJsystem = new LJsystem(*this);
    newLJsystem->SetParent(newparent);
    return newLJsystem;
}

int LJsystem::GetDiameter(int LJid) const
{
    return (int)round(GetSigma(LJid, LJid) * pow(2.0, -5.0 / 6.0) * 0.7 * 100);
}

// initialization of intramol bonds
int LJsystem::InitializeIntramolBonds(double LJ14factor, double el14factor) const
{
    int noInter = 0;
    int i, j, k;
    Molecule *mol;
    double sigma, epsilon;
    int id1, id2;

    for (i = 0; i < parent->noMolecules; i++)
    {
        mol = &(parent->molecules[i]);
        for (j = 0; j < mol->noAtoms; j++)
        {
            for (k = 0; k < j; k++)
            {
                if (mol->GetDistance(j, k) == 3)
                {
                    id1 = mol->atoms[j].LJid;
                    id2 = mol->atoms[k].LJid;
                    sigma = sigmaMatrix->operator()(id1, id2);
                    epsilon = LJ14factor * epsilonMatrix->operator()(id1, id2); // scaled interaction
                    mol->intraMolFields.push_back(new LJBond(mol, j, k, sigma, epsilon));
                    noInter++;
                }
                else if (mol->GetDistance(j, k) > 3)
                {
                    id1 = mol->atoms[j].LJid;
                    id2 = mol->atoms[k].LJid;
                    sigma = sigmaMatrix->operator()(id1, id2);
                    epsilon = epsilonMatrix->operator()(id1, id2);
                    mol->intraMolFields.push_back(new LJBond(mol, j, k, sigma, epsilon));
                    noInter++;
                }
            }
        }
    }

    return noInter;
}

// print info about the system to output stream (.prt file)
void LJsystem::PrintInfo(std::ofstream &stream, double u_eng, std::string engunit) const
{
    int i, j, k;
    Matrix *auxmatrix; // to save matrices converted to output units

    stream << "### Lennard-Jones interactions\n\n";
    stream << "The contribution to potential energy due to the pair LJ interaction is\n\n";
    if (parent->boundaryCond > 0)
    {
        stream << "$E_{\\mathrm{LJ}} =\n";
        stream << "    \\begin{cases}\n";
        stream << "        4 \\varepsilon_{ij} \\left[\\left(\\frac{\\sigma_{ij}}{r_{ij}}\\right)^{12} - \n";
        stream << "            \\left(\\frac{\\sigma_{ij}}{r_{ij}}\\right)^6\\right]      & \\text{if } r_{ij} \\leq C_{ij} \\\\\n";
        stream << "        A_{ij} (r_{ij}^2 - \\mathtt{cutoff}^2)^2                   & \\text{if } C_{ij} < r_{ij} < \\mathtt{cutoff} \\\\\n";
        stream << "        0                                                         & \\text{otherwise}\n";
        stream << "    \\end{cases}$\n\n";
        stream << "where energy parameters $\\varepsilon_{ij}$ and length parameters $\\sigma_{ij}$ are summarized in the following two tables.\n";
        stream << "The used cutoff value is $\\mathtt{cutoff} = " << cutoff << "$ [AA].\n";
        stream << "The other cutoff parameters ($A_{ij}$ and $C_{ij}$) are listed in the tables below.\n\n";
    }
    else
    {
        stream << "$E_{\\mathrm{LJ}} = 4 \\varepsilon_{ij} \\left[\\left(\\frac{\\sigma_{ij}}{r_{ij}}\\right)^{12} - \\left(\\frac{\\sigma_{ij}}{r_{ij}}\\right)^6\\right]$\n";
        stream << "where energy parameters $\\varepsilon_{ij}$ and length parameters $\\sigma_{ij}$ are summarized in the following two tables.\n\n";
    }

    stream << "Type of atom for this system is deduced from its name as given in the atom tables in molecules section.\n\n";
    // get atom names
    std::string name;
    std::vector<std::string> names;
    for (k = 0; k < size; k++)
    {
        for (i = 0; i < parent->noMolecules; i++)
        {
            for (j = 0; j < parent->molecules[i].noAtoms; j++)
            {
                if (parent->molecules[i].atoms[j].LJid == k)
                    name = (parent->molecules[i].atoms[j].name);
                continue;
            }
            if (name.size() > 0)
                continue;
        }
        names.push_back(name);
        name.clear();
    }

    stream << "**Energy parameters $\\varepsilon$ " << engunit << "**\n\n";
    auxmatrix = new Matrix(*epsilonMatrix);
    *auxmatrix = *auxmatrix * u_eng;
    auxmatrix->PrintMDTable(stream, names, names);
    stream << "\n";

    stream << "**Length parameters $\\sigma$ [AA]**\n\n";
    sigmaMatrix->PrintMDTable(stream, names, names);
    stream << "\n";

    if (parent->boundaryCond > 0)
    {
        stream << "**Cutoff parameters $A$ " << engunit << "**\n\n";
        *auxmatrix = *AMatrix * u_eng;
        auxmatrix->PrintMDTable(stream, names, names);
        stream << "\n";

        stream << "**Length parameters $C$ [AA]**\n\n";
        C1Matrix->PrintMDTable(stream, names, names);
        stream << "\n";

        if (Ecutoffcorr != 0.0)
        {
            stream << "**Cutoff corrections** due to truncation of Lennard-Jones interactions:\n\n";
            stream << "- energy correction:   $E_{\\mathrm{LJ,corr}} = " << Ecutoffcorr * u_eng << " / V$ " << engunit << "\n";
            stream << "- pressure correction: $P_{\\mathrm{LJ,corr}} = " << Ecutoffcorr * U_PRESSURE_PA << " / V^2$ "
                   << "[Pa]\n\n";
            stream << "where volume V is in [AA^3]\n\n";
        }
        else
        {
            stream << "**No cutoff corrections** are added neither to energy, nor to pressure\n\n";
        }
    }

    delete auxmatrix;
}