/*
 * VerletList class is a speed-up implementation of AbstractPairList interface
 * Serves to cycle through all intermolecular atomic pairs
 * Part of the simul++ package by JJ
 * Date: 2022/10/05
 * Modified: 2023/02/15 parallelization
 */

#include <cmath>
#include "VerletList.hpp"
#include "SimulatedSystem.hpp"
#include "Atom.hpp"
#include "Vector.hpp"
#include "Matrix.hpp"

VerletList::VerletList(SimulatedSystem *sys, double cut, double savedist)
    : AbstractPairList{sys, cut}
{
    int i;
    no_atoms = system->noAtomsTotal;
    no_molecules = system->noMolecules;
    total_no_pairs = no_atoms * (no_atoms - 1) / 2;
    origpos = new double[no_atoms][3];
    PBCdisp = new double[no_molecules][3];
    outercut = savedist * cutoff;
    outercut2 = outercut * outercut;
    maxdisp2 = pow((outercut - cutoff) * 0.5, 2.0) * 0.999; // times 0.999 for safety

    for (i = 0; i < no_atoms; i++)
    {
        origpos[i][0] = 0.0;
        origpos[i][1] = 0.0;
        origpos[i][2] = 0.0;
    }
    for (i = 0; i < no_molecules; i++)
    {
        PBCdisp[i][0] = 0.0;
        PBCdisp[i][1] = 0.0;
        PBCdisp[i][2] = 0.0;
    }
    Refresh();
}

VerletList::~VerletList()
{
    delete[] origpos;
    origpos = nullptr;
    delete[] PBCdisp;
    PBCdisp = nullptr;
}

int VerletList::Refresh()
{
    int i, j, k, m;
    double r[3];
    double box[3];
    double displ[3];
    bool newlistneeded = false;
    for (m = 0; m < 3; m++)
    {
        box[m] = system->GetBox(m);
    }

    // check the displacements and decide whether a new list is needed
    k = 0;
    for (i = 0; i < system->noMolecules; i++)
    {
        // if molecule has been moved due to the PBC, move origpos
        for (m = 0; m < 3; m++)
        {
            if (PBCdisp[i][m] != 0.0)
            {
                displ[m] = PBCdisp[i][m];
                PBCdisp[i][m] = 0.0;
            }
            else
            {
                displ[m] = 0.0;
            }
        }
        // for all atoms check the displacement
        for (j = 0; j < system->molecules[i].noAtoms; j++)
        {
            if (newlistneeded)
            {
                break;
            }
            for (m = 0; m < 3; m++)
            {
                origpos[k][m] += displ[m];
                r[m] = origpos[k][m] - system->molecules[i].atoms[j].GetPosition(m);
            }
            if ((r[0] * r[0] + r[1] * r[1] + r[2] * r[2]) > maxdisp2)
            {
                newlistneeded = true;
            }
            k++;
        }
        if (newlistneeded)
        {
            break;
        }
    }

    if (newlistneeded)
    {
        MakeList();
    }
    else
    {
        // cycle through all pairs on the list and recalculate distances
        if (system->boundaryCond > 0)
        {
#ifdef PARALLEL
#pragma omp parallel num_threads(thread_count)
#endif
            for (pair_iterator = pairs.begin(); pair_iterator != pairs.end(); pair_iterator++)
            {
                pair_iterator->CalculateDistance(box);
            }
        }
        else
        {
#ifdef PARALLEL
#pragma omp parallel num_threads(thread_count)
#endif
            for (pair_iterator = pairs.begin(); pair_iterator != pairs.end(); pair_iterator++)
            {
                pair_iterator->CalculateDistance();
            }
        }
    }
    Reset();

    return 0;
}

void VerletList::PrintInfo(std::ofstream &stream) const
{
    stream << "### Verlet list\n\n";
    stream << "List of intermolecular atomic pairs whose distance is shorter than " << outercut << " [AA] is made.\n";
    stream << "Distance between atoms is recalculated every step (before each calcualtion of forces).\n";
    stream << "Displacement of each atom from the position it had when the list was last made is calculated for each atom.\n";
    stream << "If the displacement is more then " << 0.5 * (outercut - cutoff) << " [AA], a new list is made.\n";
    stream << "Cutoff: " << cutoff << " [AA], outer cutoff: " << outercut << " [AA].\n\n";
}

int VerletList::MakeList()
{
    int i, j, k, l, m;
    Pair loc_pair;
    double box[3];
    for (i = 0; i < 3; i++)
    {
        box[i] = system->GetBox(i);
    }

    // save the positions of atoms in time of list making (origpos)
    k = 0;
    for (i = 0; i < system->noMolecules; i++)
    {
        for (m = 0; m < 3; m++)
        {
            PBCdisp[i][m] = 0.0; // zero displacement due to the periodic boundary conditions
        }
        for (j = 0; j < system->molecules[i].noAtoms; j++)
        {
            for (m = 0; m < 3; m++)
            {
                origpos[k][m] = system->molecules[i].atoms[j].GetPosition(m);
            }
            k++;
        }
    }

#ifdef PARALLEL
#pragma omp parallel num_threads(thread_count)
#endif
    pairs.clear();

    if (system->boundaryCond > 0)
    {
#ifdef PARALLEL
#pragma omp parallel for num_threads(thread_count) private(i, j, k, l, loc_pair) schedule(static, 10)
#endif
        for (i = 0; i < system->noMolecules - 1; i++)
        {
            for (j = i + 1; j < system->noMolecules; j++)
            {
                for (k = 0; k < system->molecules[i].noAtoms; k++)
                {
                    for (l = 0; l < system->molecules[j].noAtoms; l++)
                    {
                        loc_pair.SetAtoms(&(system->molecules[i].atoms[k]), &(system->molecules[j].atoms[l]));
                        if (loc_pair.CalculateDistance(box) < outercut2)
                        {
                            pairs.push_back(loc_pair);
                        }
                    }
                }
            }
        }
    }
    else
    {
#ifdef PARALLEL
#pragma omp parallel for num_threads(thread_count) private(i, j, k, l, loc_pair) schedule(static, 10)
#endif
        for (i = 0; i < system->noMolecules - 1; i++)
        {
            for (j = i + 1; j < system->noMolecules; j++)
            {
                for (k = 0; k < system->molecules[i].noAtoms; k++)
                {
                    for (l = 0; l < system->molecules[j].noAtoms; l++)
                    {
                        loc_pair.SetAtoms(&(system->molecules[i].atoms[k]), &(system->molecules[j].atoms[l]));
                        if (loc_pair.CalculateDistance() < outercut2)
                        {
                            pairs.push_back(loc_pair);
                        }
                    }
                }
            }
        }
    }

    Reset();

#ifdef PARALLEL
    int no_pairs = 0;
    std::vector<Pair> feed;
#pragma omp parallel num_threads(thread_count) reduction(+ \
                                                         : no_pairs)
    {
        no_pairs += pairs.size();
    }
    total_no_pairs = no_pairs;
    // load balancing
#pragma omp parallel num_threads(thread_count)
    {
        if (pairs.size() > (size_t) round(total_no_pairs/thread_count))
        {
#pragma omp critical
            while (pairs.size() > (size_t) round(total_no_pairs/thread_count))
            {
                feed.push_back(pairs.back());
                pairs.pop_back();
            }
        }
#pragma omp barrier
        if (pairs.size() < (size_t) round(total_no_pairs/thread_count))
        {
#pragma omp critical
            while (pairs.size() < (size_t) round(total_no_pairs/thread_count))
            {
                if (feed.empty())
                {
                    break;
                }
                pairs.push_back(feed.back());
                feed.pop_back();
            }
        }
#pragma omp barrier
#pragma omp single
        while (!feed.empty())
        {
            pairs.push_back(feed.back());
            feed.pop_back();
        }
    }

#else
    total_no_pairs = pairs.size();
#endif
    system->noPairs = total_no_pairs;
    return 0;
}

int VerletList::MoveMolecule(int molno, int dim, double distance)
{
    PBCdisp[molno][dim] += distance;
    return 0;
}
