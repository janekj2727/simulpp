/*
 * A class to compute electrostatic interactions using basic Ewald summation with splines approximation of erfc
 * Author JJ, Date Apr 2023
 *
 * Implementation of CalculateForces was inspired by MACSIMUS code (the order of sums, hyperbolic splines...)
 */

#include <fstream>
#include <map>
#include <string>
#include <cstring>
#include <cmath>

#include "Vector.hpp"
#include "Matrix.hpp"
#include "NewtonMethod.hpp"
#include "SimulatedSystem.hpp"
#include "EwaldSplinesElstat.hpp"
#include "ElstatBond.hpp"
#include "general_utils.hpp"
#include "AbstractPairList.hpp"

#define MIN_QQ_DISTANCE 1.0 // minimum intermolacular charges distance

// constructor (allocate matrix of correct sizes)
EwaldSplinesElstat::EwaldSplinesElstat(SimulatedSystem *parent, double cut, double alp, double kap, unsigned int shift, SplineType type1, unsigned int gridsize1, SplineType type2, unsigned int gridsize2)
    : AbstractInterMolField(0, parent)
{
    int i, j, n;
    double total_charge = 0.0;
    double ksq;
    double box[3];
    double rboxsq[3];
    Ecutoffcorr = 0.0;
    sft = shift;
    cutoff = cut;
    cutoff2 = cut * cut;

    double step1 = 1.0 / gridsize1;
    std::vector<double> erfcvalues;
    std::vector<double> derivatives, derivatives2, derivatives3;
    double min = MIN_QQ_DISTANCE * MIN_QQ_DISTANCE;
    double max = cutoff2 + step1;
    double r2 = min;
    double r = 0.0;
    double step2 = 1.0 / gridsize2;

    alpha = alp;
    alpha42 = 4.0 * alpha * alpha;
    if ((alpha < 0.0) || (alpha > 1e6))
    {
        print_warning(0, "Wrong value of screening parameter alpha for Ewald method (" + std::to_string(alpha) + ")!\n",
                      "    Cannot continue...\n");
        parent->error = 70;
        goto endOfConstructor;
    }
    kappa = kap;
    if ((kappa < 0.0) || (kappa > 1e6))
    {
        print_warning(0, "Wrong value of reciprocal space parameter kappa for Ewald method(" + std::to_string(kappa) + ")!\n",
                      "    Cannot continue...\n");
        parent->error = 70;
        goto endOfConstructor;
    }

    // check if free boundary conditions or neutral system... + check if any charges + count them
    no_ch_atoms = 0;
    for (i = 0; i < parent->noMolecules; i++)
    {
        if (parent->molecules[i].HasAnyCharged())
        {
            for (j = 0; j < parent->molecules[i].noAtoms; j++)
            {
                if (parent->molecules[i].atoms[j].charge != 0.0)
                {
                    no_ch_atoms++;
                }
            }
        }
    }
    if (no_ch_atoms > 1)
    {
        for (i = 0; i < parent->noMolecules; i++)
        {
            total_charge += parent->molecules[i].TotalCharge();
        }
    }
    else
    {
        goto endOfConstructor;
    }
    if (parent->boundaryCond == 0) // cannot use Ewald in free b.c.
    {
        print_warning(0, "Cannot use Ewald summation in free boundary conditions!\n",
                      "    Cannot continue...\n");
        parent->error = 69;
        no_ch_atoms = 0;
        goto endOfConstructor;
    }
    if (fabs(total_charge) > 5.0e-10) // cannot use Ewald for nonneutral system
    {
        print_warning(0, "Total charge is not zero (" + std::to_string(total_charge) + "), simulation in PBC cannot be done!\n");
        parent->error = 49;
        no_ch_atoms = 0;
        goto endOfConstructor;
    }

    // calculate shifts to assure continuity in r-space
    if (sft & 0x02)
        shiftE = erfc(alpha * cutoff) / cutoff;
    else
        shiftE = 0.0;

    if (sft & 0x01)
        shiftf = (M_2_SQRTPI * alpha * exp(-alpha * alpha * cutoff * cutoff) + shiftE) / cutoff;
    else
        shiftf = 0.0;

    // splines initialization
    Erfc = nullptr;
    Derfc = nullptr;
    // erfc calculation
    if (gridsize1 == 0)
    {
        print_warning(1, "Ewald electrostatics with splines requested, but gridsize = 0, using default value gridsize = 256!\n");
        gridsize1 = 256;
    }

    while (r2 <= max)
    {
        r = sqrt(r2);
        erfcvalues.push_back(erfc(alpha * r) / r - shiftE);
        derivatives.push_back(((-M_2_SQRTPI * alpha * exp(-alpha * alpha * r2)) / r - erfc(alpha * r) / r2) / (2 * r));
        r2 += step1;
    }
    // other types of derivatives
    derivatives2.push_back(*(derivatives.begin()));
    derivatives2.push_back(*(derivatives.end() - 1));
    derivatives3.push_back(*(derivatives.end() - 1));

    // interpolator allocations
    switch (type1)
    {
    case hyperbolic:
        Erfc = new MacsimusHyperbSplines<double>(erfcvalues, max, min, derivatives3);
        break;
    case quadratic:
        Erfc = new MacsimusQuadSplines<double>(erfcvalues, max, min, derivatives3);
        break;
    case natural:
        Erfc = new NaturalCubSplines<double>(erfcvalues, max, min, derivatives2);
        break;
    case hermite:
        Erfc = new HermiteCubSplines<double>(erfcvalues, max, min, derivatives);
        break;
    case linear:
        Erfc = new LinearInterpolator<double>(erfcvalues, max, min);
        break;
    case linear3:
        Erfc = new Linear3PInterpolator<double>(erfcvalues, max, min);
        break;
    case linear4:
        Erfc = new Linear4PInterpolator<double>(erfcvalues, max, min);
        break;
    default:
        Erfc = new MacsimusHyperbSplines<double>(erfcvalues, max, min, derivatives3);
        break;
    }
    // extend the interpolator range to zero to speedup calculation
    if (Erfc->ExtendToZero() < 1)
    {
        print_warning(0, "Unable to extend the interpolator Erfc to 0.0!", "Incorrect electrostatic energy and maybe forces expected!");
    }

    // d erfc / d r  * r (forces) calculation
    if (gridsize2 > 0)
    {
        erfcvalues.clear();
        derivatives.clear();
        derivatives2.clear();
        derivatives3.clear();
        max = cutoff2 + step2;
        r2 = min;
        r = 0.0;
        while (r2 <= max)
        {
            r = sqrt(r2);
            erfcvalues.push_back((M_2_SQRTPI * alpha * exp(-alpha * alpha * r2) + erfc(alpha * r) / r) / r2 - shiftf);
            derivatives.push_back(-M_2_SQRTPI * pow(alpha, 3.0) * exp(-alpha * alpha * r2) / r2 + 1.5 * erf(alpha * r) / (r2 * r2 * r) - 1.5 * M_2_SQRTPI * alpha * exp(-alpha * alpha * r2) / (r2 * r2));
            r2 += step2;
        }
        // other types of derivatives
        derivatives2.push_back(*(derivatives.begin()));
        derivatives2.push_back(*(derivatives.end() - 1));
        derivatives3.push_back(*(derivatives.end() - 1));

        // interpolator allocations
        switch (type2)
        {
        case hyperbolic:
            Derfc = new MacsimusHyperbSplines<double>(erfcvalues, max, min, derivatives3);
            break;
        case quadratic:
            Derfc = new MacsimusQuadSplines<double>(erfcvalues, max, min, derivatives3);
            break;
        case natural:
            Derfc = new NaturalCubSplines<double>(erfcvalues, max, min, derivatives2);
            break;
        case hermite:
            Derfc = new HermiteCubSplines<double>(erfcvalues, max, min, derivatives);
            break;
        case linear:
            Derfc = new LinearInterpolator<double>(erfcvalues, max, min);
            break;
        case linear3:
            Derfc = new Linear3PInterpolator<double>(erfcvalues, max, min);
            break;
        case linear4:
            Derfc = new Linear4PInterpolator<double>(erfcvalues, max, min);
            break;
        default:
            Derfc = new MacsimusHyperbSplines<double>(erfcvalues, max, min, derivatives3);
            break;
        }
        // extend the interpolator range to zero to speedup calculation
        if (Derfc->ExtendToZero() < 1)
        {
            print_warning(0, "Unable to extend the interpolator Derfc to 0.0!", "Incorrect electrostatic forces expected!");
        }
    }

    // charged atoms pointers allocation
    charged_atoms = new Atom *[no_ch_atoms];
    n = 0;
    for (i = 0; i < parent->noMolecules; i++)
    {
        if (parent->molecules[i].HasAnyCharged())
        {
            for (j = 0; j < parent->molecules[i].noAtoms; j++)
            {
                if (parent->molecules[i].atoms[j].charge != 0.0)
                {
                    charged_atoms[n] = &(parent->molecules[i].atoms[j]);
                    n++;
                }
            }
        }
    }

    // memory allocation
    lmax = (int)std::min(floor(kappa * parent->GetBox(0)), (double)MAX_K);
    mmax = (int)std::min(floor(kappa * parent->GetBox(1)), (double)MAX_K);
    nmax = (int)std::min(floor(kappa * parent->GetBox(2)), (double)MAX_K);
    zmult = (nmax + 1);
    ymult = zmult * (mmax + 1);
    xmult = (nmax + 1) * (mmax + 1) * (lmax + 1) * 4 + (lmax + 1) * (mmax + 1) * 2 + (lmax + 1) - 4; // -4 for not using 0,0,0 k-vector

    kymax = new int[lmax + 1];
    kzmax = new int[(lmax + 1) * (mmax + 1)];

    sins = new double[xmult * no_ch_atoms];
    coss = new double[xmult * no_ch_atoms];
    expk = new double[ymult * (lmax + 1)];

#ifdef PARALLEL
    real = new double[xmult * thread_count];
    imag = new double[xmult * thread_count];
#else
    real = new double[xmult];
    imag = new double[xmult];
#endif

    // precompute constant values (if box not changing)
    for (i = 0; i < 3; i++)
    {
        box[i] = parent->GetBox(i);
        rboxsq[i] = pow(box[i], -2.0);
    }
    for (i = 0; i <= lmax; i++)
    {
        kx[i] = i * 2.0 * M_PI;
    }
    for (i = 0; i <= mmax; i++)
    {
        ky[i] = i * 2.0 * M_PI;
    }
    for (i = 0; i <= nmax; i++)
    {
        kz[i] = i * 2.0 * M_PI;
    }

    for (i = 0; i < lmax + 1; i++)
    {
        kymax[i] = std::min((int)floor(sqrt(kappa * kappa - i * i * rboxsq[0]) * box[1]), mmax);
        for (j = 0; j < mmax + 1; j++)
        {
            kzmax[i * (lmax + 1) + j] = std::min((int)floor(sqrt(kappa * kappa - i * i * rboxsq[0] - j * j * rboxsq[1]) * box[2]), nmax);
            for (n = 0; n < nmax + 1; n++)
            {
                if ((i == 0) && (j == 0) && (n == 0))
                {
                    expk[(i)*ymult + (j)*zmult + (n)] = 0.0;
                    real[(i)*ymult + (j)*zmult + (n)] = 0.0;
                    imag[(i)*ymult + (j)*zmult + (n)] = 0.0;
                    continue;
                }
                ksq = kx[i] * kx[i] * rboxsq[0] + ky[j] * ky[j] * rboxsq[1] + kz[n] * kz[n] * rboxsq[2];
                expk[(i)*ymult + (j)*zmult + (n)] = exp(-ksq / alpha42) / ksq;
            }
        }
    }

    Eselfcorr = 0.0;
    for (i = 0; i < parent->noMolecules; i++)
    {
        for (j = 0; j < parent->molecules[i].noAtoms; j++)
        {
            Eselfcorr += pow(parent->molecules[i].atoms[j].charge, 2.0);
        }
    }
    Eselfcorr *= -alpha * M_2_SQRTPI * 0.5;

endOfConstructor:;
}

// copy constructor
EwaldSplinesElstat::EwaldSplinesElstat(const EwaldSplinesElstat &system)
    : AbstractInterMolField(0, system.parent)
{
    double ksq;
    int i, j, n;
    double box[3];
    double rboxsq[3];

    for (i = 0; i < 3; i++)
    {
        box[i] = parent->GetBox(i);
        rboxsq[i] = pow(box[i], -2.0);
    }

    cutoff = system.cutoff;
    cutoff2 = cutoff * cutoff;
    alpha = system.alpha;
    alpha42 = 4.0 * alpha * alpha;
    kappa = system.kappa;
    lmax = system.lmax;
    mmax = system.mmax;
    nmax = system.nmax;
    Ecutoffcorr = system.Ecutoffcorr;
    shiftE = system.shiftE;
    shiftf = system.shiftf;
    Erfc = nullptr;
    Derfc = nullptr;

    // allocation
    lmax = (int)floor(kappa * parent->GetBox(0));
    mmax = (int)floor(kappa * parent->GetBox(1));
    nmax = (int)floor(kappa * parent->GetBox(2));
    zmult = (nmax + 1);
    ymult = zmult * (mmax + 1);
    xmult = (nmax + 1) * (mmax + 1) * (lmax + 1) * 4 + (lmax + 1) * (mmax + 1) * 2 + (lmax + 1) - 4; // -4 for not using 0,0,0 k-vector

    kymax = new int[lmax + 1];
    kzmax = new int[(lmax + 1) * (mmax + 1)];

    sins = new double[xmult * no_ch_atoms];
    coss = new double[xmult * no_ch_atoms];
    expk = new double[ymult * (lmax + 1)];

#ifdef PARALLEL
    real = new double[xmult * thread_count];
    imag = new double[xmult * thread_count];
#else
    real = new double[xmult];
    imag = new double[xmult];
#endif

    // charged atoms pointers allocation
    charged_atoms = new Atom *[no_ch_atoms];
    n = 0;
    for (i = 0; i < parent->noMolecules; i++)
    {
        if (parent->molecules[i].HasAnyCharged())
        {
            for (j = 0; j < parent->molecules[i].noAtoms; j++)
            {
                if (parent->molecules[i].atoms[j].charge != 0.0)
                {
                    charged_atoms[n] = &(parent->molecules[i].atoms[j]);
                    n++;
                }
            }
        }
    }

    // precompute constant values
    for (i = 0; i <= lmax; i++)
    {
        kx[i] = i * 2.0 * M_PI;
    }
    for (i = 0; i <= mmax; i++)
    {
        ky[i] = i * 2.0 * M_PI;
    }
    for (i = 0; i <= nmax; i++)
    {
        kz[i] = i * 2.0 * M_PI;
    }

    for (i = 0; i < lmax + 1; i++)
    {
        kymax[i] = std::min((int)floor(sqrt(kappa * kappa - i * i * rboxsq[0]) * box[1]), mmax);
        for (j = 0; j < mmax + 1; j++)
        {
            kzmax[i * (lmax + 1) + j] = std::min((int)floor(sqrt(kappa * kappa - i * i * rboxsq[0] - j * j * rboxsq[1]) * box[2]), nmax);
            for (n = 0; n < nmax + 1; n++)
            {
                if ((i == 0) && (j == 0) && (n == 0))
                {
                    expk[(i)*ymult + (j)*zmult + (n)] = 0.0;
                    real[(i)*ymult + (j)*zmult + (n)] = 0.0;
                    imag[(i)*ymult + (j)*zmult + (n)] = 0.0;
                    continue;
                }
                ksq = kx[i] * kx[i] * rboxsq[0] + ky[j] * ky[j] * rboxsq[1] + kz[n] * kz[n] * rboxsq[2];
                expk[(i)*ymult + (j)*zmult + (n)] = exp(-ksq / alpha42) / ksq;
            }
        }
    }
    Eselfcorr = 0.0;
    for (i = 0; i < parent->noMolecules; i++)
    {
        for (j = 0; j < parent->molecules[i].noAtoms; j++)
        {
            Eselfcorr += pow(parent->molecules[i].atoms[j].charge, 2.0);
        }
    }
    Eselfcorr *= -alpha * M_2_SQRTPI * 0.5;
}

// delete matrices of coefficients
EwaldSplinesElstat::~EwaldSplinesElstat()
{
    // int i, j;
    if (no_ch_atoms > 0)
    {
        delete[] sins;
        sins = nullptr;
        delete[] coss;
        coss = nullptr;
        delete[] imag;
        imag = nullptr;
        delete[] real;
        real = nullptr;
        delete[] expk;
        expk = nullptr;
        delete[] charged_atoms;
        charged_atoms = nullptr;
        delete[] kymax;
        kymax = nullptr;
        delete[] kzmax;
        kzmax = nullptr;
        if (Erfc != nullptr)
        {
            delete Erfc;
            Erfc = nullptr;
        }
        if (Derfc != nullptr)
        {
            delete Derfc;
            Derfc = nullptr;
        }
    }
}

// calculate electrostatic forces
double EwaldSplinesElstat::CalculateForces()
{
    double Eel = 0.0;
    double Erinter = 0.0;
    double Erintra = 0.0;
    double Ek = 0.0;
    double box[3];
    double rboxsq[3];
    double ksq;
    int m, j;
    static double Vol = parent->CalculateVolume();
// int kxmax, kymax, kzmax;
#ifndef PARALLEL
    double rnorm, alphar;
#endif

    if (no_ch_atoms < 1)
    {
        return 0.0;
    }
    for (m = 0; m < 3; m++)
    {
        box[m] = parent->GetBox(m);
        rboxsq[m] = pow(box[m], -2.0);
    }

    /* r-space */
    // all pairs from different molecules
    parent->pairlist->Reset();

#ifdef PARALLEL
#pragma omp parallel num_threads(thread_count) reduction(+ \
                                                         : Erinter)
    {
        int m;
#endif
        double r2;
        double force;
        double q2;
        Atom *atom1, *atom2;
        Pair pair;
        if (Derfc != nullptr)
        {
            while ((r2 = parent->pairlist->GetNextPair(pair, cutoff2)) > 0.0)
            {
                atom1 = pair.GetAtom1();
                atom2 = pair.GetAtom2();
                // get charges (if one of them zero, skip this pair)
                q2 = atom1->charge * atom2->charge;
                if (q2 == 0.0)
                {
                    continue;
                }
                // get energy and force from interpolation
                Erinter += q2 * Erfc->Eval0(r2);
                force = q2 * Derfc->Eval0(r2);

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
            }
        }
        else // Derfc == nullptr, Erfc.Diff() is used
        {
            while ((r2 = parent->pairlist->GetNextPair(pair, cutoff2)) > 0.0)
            {
                atom1 = pair.GetAtom1();
                atom2 = pair.GetAtom2();
                // get charges (if one of them zero, skip this pair)
                q2 = atom1->charge * atom2->charge;
                if (q2 == 0.0)
                {
                    continue;
                }
                // get energy and force
                Erinter += q2 * Erfc->Eval0(r2);
                force = -q2 * (Erfc->Diff0(r2) * 2.0 - shiftf);

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
            }
        }

#ifdef PARALLEL
    }
#endif

// intramolecular interactions (r-space correction)
// all pairs from the same molecule
#ifdef PARALLEL
#pragma omp parallel num_threads(thread_count) reduction(+ \
                                                         : Erintra)
    {
        int m;
        double q2;
        double r2, rnorm;
        double force;
        double alphar; // alpha*r
#endif
        int i, k, l;
        double r[3];
        double erfrr; // erf(alpha*r)/r
#ifdef PARALLEL
#pragma omp for schedule(static, 10)
#endif
        for (i = 0; i < parent->noMolecules; i++)
        {
            for (k = 0; k < parent->molecules[i].noAtoms - 1; k++)
            {
                for (l = k + 1; l < parent->molecules[i].noAtoms; l++)
                {
                    // get charges
                    q2 = parent->molecules[i].atoms[k].charge * parent->molecules[i].atoms[l].charge;
                    if (q2 == 0.0)
                    {
                        continue;
                    }
                    // calculate r_kl (distance vector)
                    for (m = 0; m < 3; m++)
                    {
                        r[m] = parent->molecules[i].atoms[l].GetPosition(m) - parent->molecules[i].atoms[k].GetPosition(m);
                    }
                    // get nearest image... (useless for molecules smaller then half the box size)
                    for (m = 0; m < 3; m++)
                    {
                        r[m] = r[m] - box[m] * round(r[m] / box[m]);
                    }
                    // calculate square of the atom-atom distance
                    r2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
                    rnorm = sqrt(r2);
                    // get energy and force
                    alphar = alpha * rnorm;
                    erfrr = erf(alphar) / rnorm;
                    Erintra -= q2 * erfrr;
                    force = q2 / r2 * (M_2_SQRTPI * alpha * exp(-alphar * alphar) - erfrr);

                    // now, calculate forces
                    // no need for #pragma omp atomic – each thread has its molecule
                    for (m = 0; m < 3; m++)
                    {
                        parent->molecules[i].atoms[l].force[m] += force * r[m];
                        parent->molecules[i].atoms[k].force[m] -= force * r[m];
                    }
                }
            }
        }
#ifdef PARALLEL
    }
    int i, k;
#endif

    /* k-space */
    // update expk if volume changed (uses the symmetry of expk)
    // meanwhile not parallelized
    if (Vol != parent->CalculateVolume())
    {
        Vol = parent->CalculateVolume();
        for (i = 0; i < lmax + 1; i++)
        {
            for (j = 0; j < mmax + 1; j++)
            {
                for (k = 0; k < nmax + 1; k++)
                {
                    if ((i == 0) && (j == 0) && (k == 0))
                    {
                        expk[(i)*ymult + (j)*zmult + (k)] = 0.0;
                        continue;
                    }
                    ksq = kx[i] * kx[i] * rboxsq[0] + ky[j] * ky[j] * rboxsq[1] + kz[k] * kz[k] * rboxsq[2];
                    expk[(i)*ymult + (j)*zmult + (k)] = exp(-ksq / alpha42) / ksq;
                }
            }
        }
    }

// calculate sums (imag and real) - sums are symmetric, kx can be only positive and the result than multiplied by 2
// zero cumulants
#ifdef PARALLEL
    std::memset(imag, 0.0, xmult * thread_count * sizeof(double));
    std::memset(real, 0.0, xmult * thread_count * sizeof(double));
#else
    std::memset(imag, 0.0, xmult * sizeof(double));
    std::memset(real, 0.0, xmult * sizeof(double));
#endif
// kxmax = std::min((int)floor(kappa * box[0]), lmax);

// calculate terms
#ifdef PARALLEL
#pragma omp parallel num_threads(thread_count)
    {
        int i, j, k, m;
        int thread_id = omp_get_thread_num();
#else
    const int thread_id = 0;
#endif
        int w, w0, w00;
        int mxmult;
        double charge;
        double pi2rx_L, pi2ry_L, pi2rz_L; // 2 * pi * x / L_x,...
        double cosx, sinx, cosy, siny, cosz, sinz;
        double A, B, C, D;
#ifdef PARALLEL
#pragma omp for schedule(static, 10)
#endif
        for (m = 0; m < no_ch_atoms; m++)
        {
            charge = charged_atoms[m]->charge;
            mxmult = m * xmult;
            w = mxmult;
            pi2rx_L = 2.0 * M_PI * charged_atoms[m]->GetPosition(0) / box[0];
            pi2ry_L = 2.0 * M_PI * charged_atoms[m]->GetPosition(1) / box[1];
            pi2rz_L = 2.0 * M_PI * charged_atoms[m]->GetPosition(2) / box[2];
            cosx = cos(pi2rx_L);
            sinx = sin(pi2rx_L);
            cosy = cos(pi2ry_L);
            siny = sin(pi2ry_L);
            cosz = cos(pi2rz_L);
            sinz = sin(pi2rz_L);
            for (i = 0; i < lmax + 1; i++)
            {
                // j = 0;
                // k = 0;
                if (i == 0)
                {
                    sins[w] = 0;
                    coss[w] = 1;
                    w0 = w;
                    w++;
                }
                else
                {
                    // b = (1, 0, 0)
                    // a = (kx - 1, 0, 0)
                    // sin(a)sin(b)
                    A = sins[w0] * sinx;
                    // sin(a)cos(b)
                    B = sins[w0] * cosx;
                    // cos(a)sin(b)
                    C = coss[w0] * sinx;
                    // cos(a)cos(b)
                    D = coss[w0] * cosx;
                    w0 = w;
                    // (kx, 0, 0) = (kx - 1, 0, 0) + (1, 0, 0);
                    // sin(a + b) = sin(a)cos(b) + sin(b)cos(a)
                    // cos(a + b) = cos(a)cos(b) - sin(a)sin(b)
                    sins[w] = B + C;
                    coss[w] = D - A;
                    w++;
                }
                for (j = 0; j < kymax[i] + 1; j++)
                {
                    // k = 0;
                    if (j == 0)
                    {
                        w00 = w;
                        sins[w] = sins[w0];
                        coss[w] = coss[w0];
                        w++;
                        sins[w] = sins[w0];
                        coss[w] = coss[w0];
                        w++;
                    }
                    else
                    {
                        // b = (0, 1, 0)
                        // a = (kx, ky - 1, 0)
                        // sin(a)sin(b)
                        A = sins[w00] * siny;
                        // sin(a)cos(b)
                        B = sins[w00] * cosy;
                        // cos(a)sin(b)
                        C = coss[w00] * siny;
                        // cos(a)cos(b)
                        D = coss[w00] * cosy;
                        // (kx, ky, 0) = (kx, ky - 1, 0) + (0, 1, 0);
                        // sin(a + b) = sin(a)cos(b) + sin(b)cos(a)
                        // cos(a + b) = cos(a)cos(b) - sin(a)sin(b)
                        sins[w] = B + C;
                        coss[w] = D - A;
                        // a = (kx, -(ky - 1), 0)
                        // sin(a)sin(b)
                        A = sins[w00 + 1] * siny;
                        // sin(a)cos(b)
                        B = sins[w00 + 1] * cosy;
                        // cos(a)sin(b)
                        C = coss[w00 + 1] * siny;
                        // cos(a)cos(b)
                        D = coss[w00 + 1] * cosy;
                        w00 = w;
                        w++;
                        // (kx, -ky, 0) = (kx, -(ky - 1), 0) - (0, 1, 0);
                        // sin(a - b) = sin(a)cos(b) - sin(b)cos(a)
                        // cos(a - b) = cos(a)cos(b) + sin(a)sin(b)
                        sins[w] = B - C;
                        coss[w] = D + A;
                        w++;
                    }
                    // kzmax = std::min((int)floor(sqrt(kappa * kappa - i * i * rboxsq[0] - j * j * rboxsq[1]) * box[2]), nmax);
                    for (k = (((j == 0) && (i == 0)) ? (1) : 0); k < kzmax[i * (lmax + 1) + j] + 1; k++)
                    {
                        if (k == 0)
                        {
                            sins[w] = sins[w00];
                            coss[w] = coss[w00];
                            w++;
                            sins[w] = sins[w00];
                            coss[w] = coss[w00];
                            w++;
                            sins[w] = sins[w00 + 1];
                            coss[w] = coss[w00 + 1];
                            w++;
                            sins[w] = sins[w00 + 1];
                            coss[w] = coss[w00 + 1];
                            w++;
                        }
                        else if (k == 1)
                        {
                            // b = (0, 0, 1)
                            // a = (kx, ky, 0)
                            // sin(a)sin(b)
                            A = sins[w - 3] * sinz;
                            // sin(a)cos(b)
                            B = sins[w - 3] * cosz;
                            // cos(a)sin(b)
                            C = coss[w - 3] * sinz;
                            // cos(a)cos(b)
                            D = coss[w - 3] * cosz;
                            // (kx, ky, 1) = (kx, ky, 0) + (0, 0, 1);
                            // sin(a + b) = sin(a)cos(b) + sin(b)cos(a)
                            // cos(a + b) = cos(a)cos(b) - sin(a)sin(b)
                            sins[w] = B + C;
                            coss[w] = D - A;
                            w++;
                            // (kx, ky, -1) = (kx, ky, 0) - (0, 0, 1);
                            // sin(a - b) = sin(a)cos(b) - sin(b)cos(a)
                            // cos(a - b) = cos(a)cos(b) + sin(a)sin(b)
                            sins[w] = B - C;
                            coss[w] = D + A;
                            w++;
                            // a = (kx, -ky, 0)
                            // sin(a)sin(b)
                            A = sins[w - 3] * sinz;
                            // sin(a)cos(b)
                            B = sins[w - 3] * cosz;
                            // cos(a)sin(b)
                            C = coss[w - 3] * sinz;
                            // cos(a)cos(b)
                            D = coss[w - 3] * cosz;
                            // (kx, -ky, 1) = (kx, -ky, 0) + (0, 0, 1);
                            // sin(a + b) = sin(a)cos(b) + sin(b)cos(a)
                            // cos(a + b) = cos(a)cos(b) - sin(a)sin(b)
                            sins[w] = B + C;
                            coss[w] = D - A;
                            w++;
                            // (kx, -ky, -1) = (kx, -ky, 0) - (0, 0, 1);
                            // sin(a - b) = sin(a)cos(b) - sin(b)cos(a)
                            // cos(a - b) = cos(a)cos(b) + sin(a)sin(b)
                            sins[w] = B - C;
                            coss[w] = D + A;
                            w++;
                        }
                        else
                        {
                            // b = (0, 0, 1)
                            // a = (kx, ky, kz - 1)
                            // sin(a)sin(b)
                            A = sins[w - 4] * sinz;
                            // sin(a)cos(b)
                            B = sins[w - 4] * cosz;
                            // cos(a)sin(b)
                            C = coss[w - 4] * sinz;
                            // cos(a)cos(b)
                            D = coss[w - 4] * cosz;
                            // (kx, ky, kz) = (kx, ky, kz - 1) + (0, 0, kz);
                            // sin(a + b) = sin(a)cos(b) + sin(b)cos(a)
                            // cos(a + b) = cos(a)cos(b) - sin(a)sin(b)
                            sins[w] = B + C;
                            coss[w] = D - A;
                            w++;
                            // a = (kx, ky, -(kz - 1))
                            // sin(a)sin(b)
                            A = sins[w - 4] * sinz;
                            // sin(a)cos(b)
                            B = sins[w - 4] * cosz;
                            // cos(a)sin(b)
                            C = coss[w - 4] * sinz;
                            // cos(a)cos(b)
                            D = coss[w - 4] * cosz;
                            // (kx, ky, -kz) = (kx, ky, -(kz - 1)) - (0, 0, kz);
                            // sin(a - b) = sin(a)cos(b) - sin(b)cos(a)
                            // cos(a - b) = cos(a)cos(b) + sin(a)sin(b)
                            sins[w] = B - C;
                            coss[w] = D + A;
                            w++;
                            // a = (kx, -ky, kz - 1)
                            // sin(a)sin(b)
                            A = sins[w - 4] * sinz;
                            // sin(a)cos(b)
                            B = sins[w - 4] * cosz;
                            // cos(a)sin(b)
                            C = coss[w - 4] * sinz;
                            // cos(a)cos(b)
                            D = coss[w - 4] * cosz;
                            // (kx, -ky, kz) = (kx, -ky, kz - 1) + (0, 0, 1);
                            // sin(a + b) = sin(a)cos(b) + sin(b)cos(a)
                            // cos(a + b) = cos(a)cos(b) - sin(a)sin(b)
                            sins[w] = B + C;
                            coss[w] = D - A;
                            w++;
                            // a = (kx, -ky, -(kz - 1))
                            // sin(a)sin(b)
                            A = sins[w - 4] * sinz;
                            // sin(a)cos(b)
                            B = sins[w - 4] * cosz;
                            // cos(a)sin(b)
                            C = coss[w - 4] * sinz;
                            // cos(a)cos(b)
                            D = coss[w - 4] * cosz;
                            // (kx, -ky, -kz) = (kx, -ky, -(kz - 1)) - (0, 0, 1);
                            // sin(a - b) = sin(a)cos(b) - sin(b)cos(a)
                            // cos(a - b) = cos(a)cos(b) + sin(a)sin(b)
                            sins[w] = B - C;
                            coss[w] = D + A;
                            w++;
                        }

                        switch ((i > 0) + (j > 0) + (k > 0))
                        {
                        case 3:
                            // (kx, ky, kz), (kx, ky, -kz), (kx, -ky, kz), (kx, -ky, -kz)
                            imag[thread_id * xmult + w - 4 - mxmult] -= charge * sins[w - 4];
                            imag[thread_id * xmult + w - 3 - mxmult] -= charge * sins[w - 3];
                            imag[thread_id * xmult + w - 2 - mxmult] -= charge * sins[w - 2];
                            imag[thread_id * xmult + w - 1 - mxmult] -= charge * sins[w - 1];
                            real[thread_id * xmult + w - 4 - mxmult] += charge * coss[w - 4];
                            real[thread_id * xmult + w - 3 - mxmult] += charge * coss[w - 3];
                            real[thread_id * xmult + w - 2 - mxmult] += charge * coss[w - 2];
                            real[thread_id * xmult + w - 1 - mxmult] += charge * coss[w - 1];
                            break;
                        case 2:
                            if (k == 0)
                            {
                                // (kx, ky, 0), (kx, -ky, 0)
                                imag[thread_id * xmult + w - 4 - mxmult] -= charge * sins[w - 4];
                                imag[thread_id * xmult + w - 2 - mxmult] -= charge * sins[w - 2];
                                real[thread_id * xmult + w - 4 - mxmult] += charge * coss[w - 4];
                                real[thread_id * xmult + w - 2 - mxmult] += charge * coss[w - 2];
                            }
                            else
                            {
                                // (0, ky, kz), (0, ky, -kz) or (kx, 0, kz), (kx, 0, -kz)
                                imag[thread_id * xmult + w - 4 - mxmult] -= charge * sins[w - 4];
                                imag[thread_id * xmult + w - 3 - mxmult] -= charge * sins[w - 3];
                                real[thread_id * xmult + w - 4 - mxmult] += charge * coss[w - 4];
                                real[thread_id * xmult + w - 3 - mxmult] += charge * coss[w - 3];
                            }
                            break;
                        case 1:
                            // (kx, 0, 0) or (0, ky, 0) or (0, 0, kz)
                            imag[thread_id * xmult + w - 4 - mxmult] -= charge * sins[w - 4];
                            real[thread_id * xmult + w - 4 - mxmult] += charge * coss[w - 4];
                            break;
                        default:
                            print_warning(0, "Ewald sum reached point it shouldn't... see code\n");
                            break;
                        }
                    }
                }
            }
        }
#ifdef PARALLEL
#pragma omp barrier

// cummulate imag and real
#pragma omp for private(i, j)
        for (i = 0; i < xmult; i++)
        {
            for (j = 1; j < thread_count; j++)
            {
                imag[i] += imag[j * xmult + i];
                real[i] += real[j * xmult + i];
            }
        }
    }

    int w;
#endif

    // add terms to energy
    // meanwhile not parallelized
    w = 0;
    for (i = 0; i < lmax + 1; i++)
    {
        w++;
        // kymax = std::min((int)floor(sqrt(kappa * kappa - i * i * rboxsq[0]) * box[1]), mmax);
        for (j = 0; j < kymax[i] + 1; j++)
        {
            w += 2;
            // kzmax = std::min((int)floor(sqrt(kappa * kappa - i * i * rboxsq[0] - j * j * rboxsq[1]) * box[2]), nmax);
            for (k = (((j == 0) && (i == 0)) ? (1) : 0); k < kzmax[i * (lmax + 1) + j] + 1; k++)
            {
                w += 4;
                switch ((i > 0) + (j > 0) + (k > 0))
                {
                case 3:
                    // (kx, ky, kz), (kx, ky, -kz), (kx, -ky, kz), (kx, -ky, -kz)
                    Ek += expk[(i)*ymult + (j)*zmult + (k)] * (pow(real[w - 4], 2.0) + pow(imag[w - 4], 2.0) +
                                                               pow(real[w - 3], 2.0) + pow(imag[w - 3], 2.0) +
                                                               pow(real[w - 2], 2.0) + pow(imag[w - 2], 2.0) +
                                                               pow(real[w - 1], 2.0) + pow(imag[w - 1], 2.0));
                    break;
                case 2:
                    if (k == 0)
                    {
                        // (kx, ky, 0), (kx, -ky, 0)
                        Ek += expk[(i)*ymult + (j)*zmult + (k)] * (pow(real[w - 4], 2.0) + pow(imag[w - 4], 2.0) +
                                                                   pow(real[w - 2], 2.0) + pow(imag[w - 2], 2.0));
                    }
                    else
                    {
                        // (0, ky, kz), (0, ky, -kz) or (kx, 0, kz), (kx, 0, -kz)
                        Ek += expk[(i)*ymult + (j)*zmult + (k)] * (pow(real[w - 4], 2.0) + pow(imag[w - 4], 2.0) +
                                                                   pow(real[w - 3], 2.0) + pow(imag[w - 3], 2.0));
                    }
                    break;
                case 1:
                    // (kx, 0, 0) or (0, ky, 0) or (0, 0, kz)
                    Ek += expk[(i)*ymult + (j)*zmult + (k)] * (pow(real[w - 4], 2.0) + pow(imag[w - 4], 2.0));
                    break;
                default:
                    print_warning(0, "Ewald sum reached point it shouldn't... see code\n");
                    break;
                }
            }
        }
    }

    // final reciprocal space energy (2*pi/V*Ek...) (times 2.0 from symmetry (see above))
    Ek *= 4.0 * M_PI / Vol;

// forces from k-space
#ifdef PARALLEL
#pragma omp parallel num_threads(thread_count)
    {
        int i, j, k, l, m, mxmult;
        int w;
        double A, B, C, D;
#endif
        double force3[3];
#ifdef PARALLEL
#pragma omp for schedule(static, 10)
#endif
        for (m = 0; m < no_ch_atoms; m++)
        {
            mxmult = m * xmult;
            w = 0;
            std::fill(&force3[0], &force3[0] + 3, 0.0);
            for (i = 0; i < lmax + 1; i++)
            {
                w++;
                // kymax = std::min((int)floor(sqrt(kappa * kappa - i * i * rboxsq[0]) * box[1]), mmax);
                for (j = 0; j < kymax[i] + 1; j++)
                {
                    w += 2;
                    // kzmax = std::min((int)floor(sqrt(kappa * kappa - i * i * rboxsq[0] - j * j * rboxsq[1]) * box[2]), nmax);
                    for (k = (((j == 0) && (i == 0)) ? (1) : 0); k < kzmax[i * (lmax + 1) + j] + 1; k++)
                    {
                        w += 4;
                        switch ((i > 0) + (j > 0) + (k > 0))
                        {
                        case 3:
                            // (kx, ky, kz), (kx, ky, -kz), (kx, -ky, kz), (kx, -ky, -kz)
                            A = expk[(i)*ymult + (j)*zmult + (k)] * (sins[mxmult + w - 4] * real[w - 4] + coss[mxmult + w - 4] * imag[w - 4]);
                            B = expk[(i)*ymult + (j)*zmult + (k)] * (sins[mxmult + w - 3] * real[w - 3] + coss[mxmult + w - 3] * imag[w - 3]);
                            C = expk[(i)*ymult + (j)*zmult + (k)] * (sins[mxmult + w - 2] * real[w - 2] + coss[mxmult + w - 2] * imag[w - 2]);
                            D = expk[(i)*ymult + (j)*zmult + (k)] * (sins[mxmult + w - 1] * real[w - 1] + coss[mxmult + w - 1] * imag[w - 1]);
                            force3[0] += A * kx[i] / box[0] + B * kx[i] / box[0] + C * kx[i] / box[0] + D * kx[i] / box[0];
                            force3[1] += A * ky[j] / box[1] + B * ky[j] / box[1] - C * ky[j] / box[1] - D * ky[j] / box[1];
                            force3[2] += A * kz[k] / box[2] - B * kz[k] / box[2] + C * kz[k] / box[2] - D * kz[k] / box[2];
                            break;
                        case 2:
                            if (k == 0)
                            {
                                // (kx, ky, 0), (kx, -ky, 0)
                                A = expk[(i)*ymult + (j)*zmult + (k)] * (sins[mxmult + w - 4] * real[w - 4] + coss[mxmult + w - 4] * imag[w - 4]);
                                B = expk[(i)*ymult + (j)*zmult + (k)] * (sins[mxmult + w - 2] * real[w - 2] + coss[mxmult + w - 2] * imag[w - 2]);
                                force3[0] += A * kx[i] / box[0] + B * kx[i] / box[0];
                                force3[1] += A * ky[j] / box[1] - B * ky[j] / box[1];
                                // kz[0] = 0;
                            }
                            else
                            {
                                // (0, ky, kz), (0, ky, -kz) or (kx, 0, kz), (kx, 0, -kz)
                                A = expk[(i)*ymult + (j)*zmult + (k)] * (sins[mxmult + w - 4] * real[w - 4] + coss[mxmult + w - 4] * imag[w - 4]);
                                B = expk[(i)*ymult + (j)*zmult + (k)] * (sins[mxmult + w - 3] * real[w - 3] + coss[mxmult + w - 3] * imag[w - 3]);
                                force3[0] += A * kx[i] / box[0] + B * kx[i] / box[0];
                                force3[1] += A * ky[j] / box[1] + B * ky[j] / box[1];
                                force3[2] += A * kz[k] / box[2] - B * kz[k] / box[2];
                            }
                            break;
                        case 1:
                            // (kx, 0, 0) or (0, ky, 0) or (0, 0, kz)
                            A = expk[(i)*ymult + (j)*zmult + (k)] * (sins[mxmult + w - 4] * real[w - 4] + coss[mxmult + w - 4] * imag[w - 4]);
                            force3[0] += A * kx[i] / box[0];
                            force3[1] += A * ky[j] / box[1];
                            force3[2] += A * kz[k] / box[2];
                            break;
                        default:
                            print_warning(0, "Ewald sum reached point it shouldn't... see code\n");
                            break;
                        }
                    }
                }
            }
            for (l = 0; l < 3; l++)
            {
                // multiplied by 2 for symmetry reasons
                charged_atoms[m]->force[l] += 8.0 * M_PI * charged_atoms[m]->charge * force3[l] / Vol;
            }
        }
#ifdef PARALLEL
    }
#endif

    /* total energy */

#ifdef DEBUG // stored for use in future...
    std::cout << "Forces: Er_inter = " << Erinter << ", Er_intra = " << Erintra << ", Ek = " << Ek << ", Eself = " << Eselfcorr << "\n";
#endif

    Eel = Erinter + Erintra + Ek + Eselfcorr;

    parent->virintermol[type] -= Eel;
    parent->Epotintermol[type] += Eel;
    return Eel;
}

// calculate electrostatic potential energy
double EwaldSplinesElstat::CalculateEpot() const
{
    double Eel = 0.0;
    double Erinter = 0.0;
    double Erintra = 0.0;
    double Ek = 0.0;
    int m;
    double box[3];
    double rboxsq[3];
    double ksq;
    double Vol;
    double localexpk; // local storage for expk (not to change the one used in CalculateForces)
#ifndef PARALLEL
    double rnorm, alphar;
#endif

    if (no_ch_atoms < 1)
    {
        return 0.0;
    }
    for (m = 0; m < 3; m++)
    {
        box[m] = parent->GetBox(m);
        rboxsq[m] = pow(box[m], -2.0);
    }

    /* r-space */
    // all pairs from different molecules
    parent->pairlist->Reset();

#ifdef PARALLEL
#pragma omp parallel num_threads(thread_count) reduction(+ \
                                                         : Erinter)
    {
#endif
        double r2;
        double q2;
        Atom *atom1, *atom2;
        Pair pair;
        while ((r2 = parent->pairlist->GetNextPair(pair, cutoff2)) > 0.0)
        {
            atom1 = pair.GetAtom1();
            atom2 = pair.GetAtom2();
            // get charges
            q2 = atom1->charge * atom2->charge;
            if (q2 == 0.0)
            {
                continue;
            }
            // get energy and force
            Erinter += q2 * Erfc->Eval0(r2);
        }
#ifdef PARALLEL
    }
#endif

// intramolecular interactions (r-space correction)
// all pairs from the same molecule
#ifdef PARALLEL
#pragma omp parallel num_threads(thread_count) reduction(+ \
                                                         : Erintra)
    {
        int m;
        double q2;
        double r2, rnorm;
        double alphar; // alpha*r
#endif
        int i, k, l;
        double r[3];
        double erfrr; // erf(alpha*r)/r
#ifdef PARALLEL
#pragma omp for schedule(static, 10)
#endif
        for (i = 0; i < parent->noMolecules; i++)
        {
            for (k = 0; k < parent->molecules[i].noAtoms - 1; k++)
            {
                for (l = k + 1; l < parent->molecules[i].noAtoms; l++)
                {
                    // get charges
                    q2 = parent->molecules[i].atoms[k].charge * parent->molecules[i].atoms[l].charge;
                    if (q2 == 0.0)
                    {
                        continue;
                    }
                    // calculate distance vectro r_kl
                    for (m = 0; m < 3; m++)
                    {
                        r[m] = parent->molecules[i].atoms[l].GetPosition(m) - parent->molecules[i].atoms[k].GetPosition(m);
                    }
                    // get nearest image... (useless for molecules smaller then half the box size)
                    for (m = 0; m < 3; m++)
                    {
                        r[m] = r[m] - box[m] * round(r[m] / box[m]);
                    }
                    // calculate square of the atom-atom distance
                    r2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
                    rnorm = sqrt(r2);

                    // get energy and force
                    alphar = alpha * rnorm;
                    erfrr = erf(alphar) / rnorm;
                    Erintra -= q2 * erfrr;
                }
            }
        }
#ifdef PARALLEL
    }
#endif
    /* k-space */
    // update volume, expk is not updated in this function (preserve expk to CalculateForces)
    Vol = parent->CalculateVolume();

    // calculate sums (imag and real), sums are symmetric, the first index can span only half the range
    // could be truncated spherically (as forces are) but the resulting PVVC would be quite different
    // thus spherical cutoff in k-space not used here...
    // calculate sums (imag and real) - sums are symmetric, kx can be only positive and the result than multiplied by 2
    // zero cumulants
#ifdef PARALLEL
    std::memset(imag, 0.0, xmult * thread_count * sizeof(double));
    std::memset(real, 0.0, xmult * thread_count * sizeof(double));
#else
    std::memset(imag, 0.0, xmult * sizeof(double));
    std::memset(real, 0.0, xmult * sizeof(double));
#endif
    // kxmax = std::min((int)floor(kappa * box[0]), lmax);
    int kxmax = lmax;
// calculate terms
#ifdef PARALLEL
#pragma omp parallel num_threads(thread_count)
    {
        int i, k, m;
        int thread_id = omp_get_thread_num();
#else
    const int thread_id = 0;
#endif
        int j;
        int w, w0, w00;
        int mxmult;
        double charge;
        double pi2rx_L, pi2ry_L, pi2rz_L; // 2 * pi * x / L_x,...
        double cosx, sinx, cosy, siny, cosz, sinz;
        double A, B, C, D;
        int kymax, kzmax;
#ifdef PARALLEL
#pragma omp for schedule(static, 10)
#endif
        for (m = 0; m < no_ch_atoms; m++)
        {
            charge = charged_atoms[m]->charge;
            mxmult = m * xmult;
            w = mxmult;
            pi2rx_L = 2.0 * M_PI * charged_atoms[m]->GetPosition(0) / box[0];
            pi2ry_L = 2.0 * M_PI * charged_atoms[m]->GetPosition(1) / box[1];
            pi2rz_L = 2.0 * M_PI * charged_atoms[m]->GetPosition(2) / box[2];
            cosx = cos(pi2rx_L);
            sinx = sin(pi2rx_L);
            cosy = cos(pi2ry_L);
            siny = sin(pi2ry_L);
            cosz = cos(pi2rz_L);
            sinz = sin(pi2rz_L);
            for (i = 0; i < kxmax + 1; i++)
            {
                // j = 0;
                // k = 0;
                if (i == 0)
                {
                    sins[w] = 0;
                    coss[w] = 1;
                    w0 = w;
                    w++;
                }
                else
                {
                    // b = (1, 0, 0)
                    // a = (kx - 1, 0, 0)
                    // sin(a)sin(b)
                    A = sins[w0] * sinx;
                    // sin(a)cos(b)
                    B = sins[w0] * cosx;
                    // cos(a)sin(b)
                    C = coss[w0] * sinx;
                    // cos(a)cos(b)
                    D = coss[w0] * cosx;
                    w0 = w;
                    // (kx, 0, 0) = (kx - 1, 0, 0) + (1, 0, 0);
                    // sin(a + b) = sin(a)cos(b) + sin(b)cos(a)
                    // cos(a + b) = cos(a)cos(b) - sin(a)sin(b)
                    sins[w] = B + C;
                    coss[w] = D - A;
                    w++;
                }

                // kymax = std::min((int)floor(sqrt(kappa * kappa - i * i * rboxsq[0]) * box[1]), mmax);
                kymax = mmax;
                for (j = 0; j < kymax + 1; j++)
                {
                    // k = 0;
                    if (j == 0)
                    {
                        w00 = w;
                        sins[w] = sins[w0];
                        coss[w] = coss[w0];
                        w++;
                        sins[w] = sins[w0];
                        coss[w] = coss[w0];
                        w++;
                    }
                    else
                    {
                        // b = (0, 1, 0)
                        // a = (kx, ky - 1, 0)
                        // sin(a)sin(b)
                        A = sins[w00] * siny;
                        // sin(a)cos(b)
                        B = sins[w00] * cosy;
                        // cos(a)sin(b)
                        C = coss[w00] * siny;
                        // cos(a)cos(b)
                        D = coss[w00] * cosy;
                        // (kx, ky, 0) = (kx, ky - 1, 0) + (0, 1, 0);
                        // sin(a + b) = sin(a)cos(b) + sin(b)cos(a)
                        // cos(a + b) = cos(a)cos(b) - sin(a)sin(b)
                        sins[w] = B + C;
                        coss[w] = D - A;
                        // a = (kx, -(ky - 1), 0)
                        // sin(a)sin(b)
                        A = sins[w00 + 1] * siny;
                        // sin(a)cos(b)
                        B = sins[w00 + 1] * cosy;
                        // cos(a)sin(b)
                        C = coss[w00 + 1] * siny;
                        // cos(a)cos(b)
                        D = coss[w00 + 1] * cosy;
                        w00 = w;
                        w++;
                        // (kx, -ky, 0) = (kx, -(ky - 1), 0) - (0, 1, 0);
                        // sin(a - b) = sin(a)cos(b) - sin(b)cos(a)
                        // cos(a - b) = cos(a)cos(b) + sin(a)sin(b)
                        sins[w] = B - C;
                        coss[w] = D + A;
                        w++;
                    }
                    // kzmax = std::min((int)floor(sqrt(kappa * kappa - i * i * rboxsq[0] - j * j * rboxsq[1]) * box[2]), nmax);
                    kzmax = nmax;
                    for (k = (((j == 0) && (i == 0)) ? (1) : 0); k < kzmax + 1; k++)
                    {
                        if (k == 0)
                        {
                            sins[w] = sins[w00];
                            coss[w] = coss[w00];
                            w++;
                            sins[w] = sins[w00];
                            coss[w] = coss[w00];
                            w++;
                            sins[w] = sins[w00 + 1];
                            coss[w] = coss[w00 + 1];
                            w++;
                            sins[w] = sins[w00 + 1];
                            coss[w] = coss[w00 + 1];
                            w++;
                        }
                        else if (k == 1)
                        {
                            // b = (0, 0, 1)
                            // a = (kx, ky, 0)
                            // sin(a)sin(b)
                            A = sins[w - 3] * sinz;
                            // sin(a)cos(b)
                            B = sins[w - 3] * cosz;
                            // cos(a)sin(b)
                            C = coss[w - 3] * sinz;
                            // cos(a)cos(b)
                            D = coss[w - 3] * cosz;
                            // (kx, ky, 1) = (kx, ky, 0) + (0, 0, 1);
                            // sin(a + b) = sin(a)cos(b) + sin(b)cos(a)
                            // cos(a + b) = cos(a)cos(b) - sin(a)sin(b)
                            sins[w] = B + C;
                            coss[w] = D - A;
                            w++;
                            // (kx, ky, -1) = (kx, ky, 0) - (0, 0, 1);
                            // sin(a - b) = sin(a)cos(b) - sin(b)cos(a)
                            // cos(a - b) = cos(a)cos(b) + sin(a)sin(b)
                            sins[w] = B - C;
                            coss[w] = D + A;
                            w++;
                            // a = (kx, -ky, 0)
                            // sin(a)sin(b)
                            A = sins[w - 3] * sinz;
                            // sin(a)cos(b)
                            B = sins[w - 3] * cosz;
                            // cos(a)sin(b)
                            C = coss[w - 3] * sinz;
                            // cos(a)cos(b)
                            D = coss[w - 3] * cosz;
                            // (kx, -ky, 1) = (kx, -ky, 0) + (0, 0, 1);
                            // sin(a + b) = sin(a)cos(b) + sin(b)cos(a)
                            // cos(a + b) = cos(a)cos(b) - sin(a)sin(b)
                            sins[w] = B + C;
                            coss[w] = D - A;
                            w++;
                            // (kx, -ky, -1) = (kx, -ky, 0) - (0, 0, 1);
                            // sin(a - b) = sin(a)cos(b) - sin(b)cos(a)
                            // cos(a - b) = cos(a)cos(b) + sin(a)sin(b)
                            sins[w] = B - C;
                            coss[w] = D + A;
                            w++;
                        }
                        else
                        {
                            // b = (0, 0, 1)
                            // a = (kx, ky, kz - 1)
                            // sin(a)sin(b)
                            A = sins[w - 4] * sinz;
                            // sin(a)cos(b)
                            B = sins[w - 4] * cosz;
                            // cos(a)sin(b)
                            C = coss[w - 4] * sinz;
                            // cos(a)cos(b)
                            D = coss[w - 4] * cosz;
                            // (kx, ky, kz) = (kx, ky, kz - 1) + (0, 0, kz);
                            // sin(a + b) = sin(a)cos(b) + sin(b)cos(a)
                            // cos(a + b) = cos(a)cos(b) - sin(a)sin(b)
                            sins[w] = B + C;
                            coss[w] = D - A;
                            w++;
                            // a = (kx, ky, -(kz - 1))
                            // sin(a)sin(b)
                            A = sins[w - 4] * sinz;
                            // sin(a)cos(b)
                            B = sins[w - 4] * cosz;
                            // cos(a)sin(b)
                            C = coss[w - 4] * sinz;
                            // cos(a)cos(b)
                            D = coss[w - 4] * cosz;
                            // (kx, ky, -kz) = (kx, ky, -(kz - 1)) - (0, 0, kz);
                            // sin(a - b) = sin(a)cos(b) - sin(b)cos(a)
                            // cos(a - b) = cos(a)cos(b) + sin(a)sin(b)
                            sins[w] = B - C;
                            coss[w] = D + A;
                            w++;
                            // a = (kx, -ky, kz - 1)
                            // sin(a)sin(b)
                            A = sins[w - 4] * sinz;
                            // sin(a)cos(b)
                            B = sins[w - 4] * cosz;
                            // cos(a)sin(b)
                            C = coss[w - 4] * sinz;
                            // cos(a)cos(b)
                            D = coss[w - 4] * cosz;
                            // (kx, -ky, kz) = (kx, -ky, kz - 1) + (0, 0, 1);
                            // sin(a + b) = sin(a)cos(b) + sin(b)cos(a)
                            // cos(a + b) = cos(a)cos(b) - sin(a)sin(b)
                            sins[w] = B + C;
                            coss[w] = D - A;
                            w++;
                            // a = (kx, -ky, -(kz - 1))
                            // sin(a)sin(b)
                            A = sins[w - 4] * sinz;
                            // sin(a)cos(b)
                            B = sins[w - 4] * cosz;
                            // cos(a)sin(b)
                            C = coss[w - 4] * sinz;
                            // cos(a)cos(b)
                            D = coss[w - 4] * cosz;
                            // (kx, -ky, -kz) = (kx, -ky, -(kz - 1)) - (0, 0, 1);
                            // sin(a - b) = sin(a)cos(b) - sin(b)cos(a)
                            // cos(a - b) = cos(a)cos(b) + sin(a)sin(b)
                            sins[w] = B - C;
                            coss[w] = D + A;
                            w++;
                        }

                        switch ((i > 0) + (j > 0) + (k > 0))
                        {
                        case 3:
                            // (kx, ky, kz), (kx, ky, -kz), (kx, -ky, kz), (kx, -ky, -kz)
                            imag[thread_id * xmult + w - 4 - mxmult] -= charge * sins[w - 4];
                            imag[thread_id * xmult + w - 3 - mxmult] -= charge * sins[w - 3];
                            imag[thread_id * xmult + w - 2 - mxmult] -= charge * sins[w - 2];
                            imag[thread_id * xmult + w - 1 - mxmult] -= charge * sins[w - 1];
                            real[thread_id * xmult + w - 4 - mxmult] += charge * coss[w - 4];
                            real[thread_id * xmult + w - 3 - mxmult] += charge * coss[w - 3];
                            real[thread_id * xmult + w - 2 - mxmult] += charge * coss[w - 2];
                            real[thread_id * xmult + w - 1 - mxmult] += charge * coss[w - 1];
                            break;
                        case 2:
                            if (k == 0)
                            {
                                // (kx, ky, 0), (kx, -ky, 0)
                                imag[thread_id * xmult + w - 4 - mxmult] -= charge * sins[w - 4];
                                imag[thread_id * xmult + w - 2 - mxmult] -= charge * sins[w - 2];
                                real[thread_id * xmult + w - 4 - mxmult] += charge * coss[w - 4];
                                real[thread_id * xmult + w - 2 - mxmult] += charge * coss[w - 2];
                            }
                            else
                            {
                                // (0, ky, kz), (0, ky, -kz) or (kx, 0, kz), (kx, 0, -kz)
                                imag[thread_id * xmult + w - 4 - mxmult] -= charge * sins[w - 4];
                                imag[thread_id * xmult + w - 3 - mxmult] -= charge * sins[w - 3];
                                real[thread_id * xmult + w - 4 - mxmult] += charge * coss[w - 4];
                                real[thread_id * xmult + w - 3 - mxmult] += charge * coss[w - 3];
                            }
                            break;
                        case 1:
                            // (kx, 0, 0) or (0, ky, 0) or (0, 0, kz)
                            imag[thread_id * xmult + w - 4 - mxmult] -= charge * sins[w - 4];
                            real[thread_id * xmult + w - 4 - mxmult] += charge * coss[w - 4];
                            break;
                        default:
                            print_warning(0, "Ewald sum reached point it shouldn't... see code\n");
                            break;
                        }
                    }
                }
            }
        }
#ifdef PARALLEL
#pragma omp barrier

// cummulate imag and real
#pragma omp for private(i, j)
        for (i = 0; i < xmult; i++)
        {
            for (j = 1; j < thread_count; j++)
            {
                imag[i] += imag[j * xmult + i];
                real[i] += real[j * xmult + i];
            }
        }
    }

    int w;
    int i, j, k, kymax, kzmax;
#endif

    // add terms to energy
    // meanwhile not parallelized
    w = 0;
    for (i = 0; i < kxmax + 1; i++)
    {
        w++;
        // kymax = std::min((int)floor(sqrt(kappa * kappa - i * i * rboxsq[0]) * box[1]), mmax);
        kymax = mmax;
        for (j = 0; j < kymax + 1; j++)
        {
            w += 2;
            // kzmax = std::min((int)floor(sqrt(kappa * kappa - i * i * rboxsq[0] - j * j * rboxsq[1]) * box[2]), nmax);
            kzmax = nmax;
            for (k = (((j == 0) && (i == 0)) ? (1) : 0); k < kzmax + 1; k++)
            {
                w += 4;
                ksq = kx[i] * kx[i] * rboxsq[0] + ky[j] * ky[j] * rboxsq[1] + kz[k] * kz[k] * rboxsq[2];
                localexpk = exp(-ksq / alpha42) / ksq;
                switch ((i > 0) + (j > 0) + (k > 0))
                {
                case 3:
                    // (kx, ky, kz), (kx, ky, -kz), (kx, -ky, kz), (kx, -ky, -kz)
                    Ek += localexpk * (pow(real[w - 4], 2.0) + pow(imag[w - 4], 2.0) +
                                       pow(real[w - 3], 2.0) + pow(imag[w - 3], 2.0) +
                                       pow(real[w - 2], 2.0) + pow(imag[w - 2], 2.0) +
                                       pow(real[w - 1], 2.0) + pow(imag[w - 1], 2.0));
                    break;
                case 2:
                    if (k == 0)
                    {
                        // (kx, ky, 0), (kx, -ky, 0)
                        Ek += localexpk * (pow(real[w - 4], 2.0) + pow(imag[w - 4], 2.0) +
                                           pow(real[w - 2], 2.0) + pow(imag[w - 2], 2.0));
                    }
                    else
                    {
                        // (0, ky, kz), (0, ky, -kz) or (kx, 0, kz), (kx, 0, -kz)
                        Ek += localexpk * (pow(real[w - 4], 2.0) + pow(imag[w - 4], 2.0) +
                                           pow(real[w - 3], 2.0) + pow(imag[w - 3], 2.0));
                    }
                    break;
                case 1:
                    // (kx, 0, 0) or (0, ky, 0) or (0, 0, kz)
                    Ek += localexpk * (pow(real[w - 4], 2.0) + pow(imag[w - 4], 2.0));
                    break;
                default:
                    print_warning(0, "Ewald sum reached point it shouldn't... see code\n");
                    break;
                }
            }
        }
    }

    // final reciprocal space energy (2*pi/V*Ek...) (times 2.0 from symmetry (see above))
    Ek *= 4.0 * M_PI / Vol;

/* total energy */
#ifdef DEBUG
    std::cout << "Epot:   Er_inter = " << Erinter << ", Er_intra = " << Erintra << ", Ek = " << Ek << ", Eself = " << Eselfcorr << "\n";
#endif
    Eel = Erinter + Erintra + Ek + Eselfcorr;

    return Eel;
}

// return energy correction (times volume)
double EwaldSplinesElstat::CalculateCutoffCorrection()
{
    Ecutoffcorr = 0.0;
    return 0.0;
}

// set electrostatic interactions cutoff
int EwaldSplinesElstat::SetCutoff(double vdwcut, double cut, int boundaryC)
{
    if (boundaryC == 0)
    {
        cutoff = 1000.0; // effectively infinity, should never be needed
        cutoff2 = 1e6;
        print_warning(0, "Cannot use Ewald sum in free boundary conditions!!!\n");
    }
    else
    {
        cutoff = cut;
        cutoff2 = cut * cut;
        // calculate shifts to assure continuity in r-space
        if (sft & 0x02)
            shiftE = erfc(alpha * cutoff) / cutoff;
        else
            shiftE = 0.0;

        if (sft & 0x01)
            shiftf = (M_2_SQRTPI * alpha * exp(-alpha * alpha * cutoff * cutoff) + shiftE) / cutoff;
        else
            shiftf = 0.0;
    }

    return 0;
}

// cloning
EwaldSplinesElstat *EwaldSplinesElstat::copy(SimulatedSystem *newparent) const
{
    EwaldSplinesElstat *newEwaldSplinesElstat = new EwaldSplinesElstat(*this);
    newEwaldSplinesElstat->SetParent(newparent);
    return newEwaldSplinesElstat;
}

// intramol electrostatic interactions
// initialization of intramol bonds
int EwaldSplinesElstat::InitializeIntramolBonds(double LJ14factor, double el14factor) const
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
                    mol->intraMolFields.push_back(new ElstatBond(mol, j, k, 0.0, el14factor));
                    noInter++;
                }
                else if (mol->GetDistance(j, k) > 3)
                {
                    mol->intraMolFields.push_back(new ElstatBond(mol, j, k, 0.0));
                    noInter++;
                }
            }
        }
    }

    return noInter;
}

// print info about the system to output stream (.prt file)
void EwaldSplinesElstat::PrintInfo(std::ofstream &stream, double u_eng, std::string engunit) const
{
    stream << "### Electrostatic interactions calculated by Ewald summation\n\n";
    stream << "The pair contribution to potential energy due to the electrostatic interaction is\n\n";
    stream << "$E_{\\mathrm{el}} = \\frac{q_i q_j}{4\\pi \\epsilon_0}\\frac{1}{r_{ij}}$\n\n";
    stream << "where charges $q_i$ and $q_j$ are charges of atoms (can be found above in atom tables in the molecules section).\n\n";
    stream << "Electrostatic interactions are separated into short-range part, which is calculated in $r$-space,\n";
    stream << "and long-range part, which is computed in $k$-space.\n";
    stream << "No other long-range corrections are considered.\n\n";
    stream << "The cutoff distance for short-range part is " << cutoff << " [AA] and screening parameter $\\alpha = " << alpha << "$.\n";
    if (shiftE != 0.0)
    {
        stream << "The real-space energy is shifted by " << shiftE * u_eng << " " << engunit << " multiplied by charges ($q_i q_j$) to avoid discontinuities at cutoff.\n";
    }
    if (shiftf != 0.0)
    {
        stream << "The real-space force is shifted by " << shiftf << " [p.u.] multiplied by charges ($q_i q_j$) to avoid discontinuities at cutoff.\n";
    }
    stream << "The real-space energy and force is calculated by interpolation of the exact formula:\n";
    stream << "- energy: \n";
    Erfc->PrintInfo(stream, 4);
    if (Derfc == nullptr)
    {
        stream << "- force: \n    - function Diff() of the interpolator for energy is used\n";
    }
    else
    {
        stream << "- force:\n";
        Derfc->PrintInfo(stream, 4);
    }
    stream << "\n";
    stream << "Long range parameter $\\kappa = " << kappa << "$ leads to $k$-vector indeces ";
    stream << "$\\left[ l, m, n \\right] = \\left[ -" << lmax << "\\dots +" << lmax << ", -" << mmax << "\\dots +" << mmax << ", -" << nmax << "\\dots +" << nmax << "\\right]$.\n";
    stream << "For the calculation of forces (not for PVVC), the $k$-vector is truncated spherically so that $\\left|\\left|\\mathbf{k}\\right|\\right| < 2 \\pi \\kappa$.\n\n";
}