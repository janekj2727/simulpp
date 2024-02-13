/*
 * AllPairList class is the simplest implementation of AbstractPairList interface
 * Serves to cycle through all intermolecular atomic pairs
 * Part of the simul++ package by JJ
 * Date: 2022/10/05
 * Modified: 2023/02/15 parallelization
 */

#include <cmath>
#include "AllPairList.hpp"
#include "Atom.hpp"
#include "Vector.hpp"
#include "Matrix.hpp"

std::vector<Pair> AbstractPairList::pairs;
std::vector<Pair>::iterator AbstractPairList::pair_iterator;

AllPairList::AllPairList(SimulatedSystem *sys, double cut)
    : AbstractPairList{sys, cut}
{
    int i, j, k, l;
    Pair loc_pair;
#ifdef PARALLEL
    std::vector<Pair> global_pairs;
#endif
    double box[3];

    for (i = 0; i < 3; i++)
    {
        box[i] = system->GetBox(i);
    }

    if (system->boundaryCond > 0)
    {
        for (i = 0; i < system->noMolecules - 1; i++)
        {
            for (j = i + 1; j < system->noMolecules; j++)
            {
                for (k = 0; k < system->molecules[i].noAtoms; k++)
                {
                    for (l = 0; l < system->molecules[j].noAtoms; l++)
                    {
                        loc_pair.SetAtoms(&(system->molecules[i].atoms[k]), &(system->molecules[j].atoms[l]));
                        loc_pair.CalculateDistance(box);
#ifdef PARALLEL
                        global_pairs.push_back(loc_pair);
#else
                        pairs.push_back(loc_pair);
#endif
                    }
                }
            }
        }
    }
    else
    {
        for (i = 0; i < system->noMolecules - 1; i++)
        {
            for (j = i + 1; j < system->noMolecules; j++)
            {
                for (k = 0; k < system->molecules[i].noAtoms; k++)
                {
                    for (l = 0; l < system->molecules[j].noAtoms; l++)
                    {
                        loc_pair.SetAtoms(&(system->molecules[i].atoms[k]), &(system->molecules[j].atoms[l]));
                        loc_pair.CalculateDistance();
#ifdef PARALLEL
                        global_pairs.push_back(loc_pair);
#else
                        pairs.push_back(loc_pair);
#endif
                    }
                }
            }
        }
    }

#ifdef PARALLEL
    // assign pairs to different threads
    total_no_pairs = global_pairs.size();
    system->noPairs = global_pairs.size();
    int per_thread = (int) round((double) total_no_pairs / thread_count); 
#pragma omp parallel num_threads(thread_count)
    {
        int my_rank = omp_get_thread_num();
        if (my_rank < thread_count - 1)
        {
            pairs.insert(pairs.begin(), global_pairs.begin() + (my_rank * per_thread), global_pairs.begin() + ((my_rank + 1) * per_thread));
        }
        else
        {
            pairs.insert(pairs.begin(), global_pairs.begin() + (my_rank * per_thread), global_pairs.end());
        }      
    }
#else
    total_no_pairs = pairs.size();
    system->noPairs = pairs.size();
#endif
}

AllPairList::~AllPairList()
{
    // nothing needed
}

int AllPairList::Refresh()
{
    int i;
    double box[3];

    for (i = 0; i < 3; i++)
    {
        box[i] = system->GetBox(i);
    }

    // refresh distances
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

    Reset();

    return 0;
}

void AllPairList::PrintInfo(std::ofstream &stream) const
{
    stream << "### All-pairs list\n\n";
    stream << "Distance between atoms is recalculated every step (before each calculation of forces), all interatomic pairs are considered each time.\n";
    stream << "This is the simplest (and default) possibility.\n";
    stream << "The actual number of pairs closer than the largest considered cutoff (" << cutoff << " [AA]) can be monitored using `measure nopairs` directive.\n";
    stream << "Total number of pairs: " << total_no_pairs << ", largest distance considered: " << cutoff << " [AA].\n\n";
}