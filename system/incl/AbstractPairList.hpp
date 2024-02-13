#ifndef ABSTRACTPAIRLISTHEADER
#define ABSTRACTPAIRLISTHEADER

#define MAX_THREAD_COUNT 10

class SimulatedSystem;
class Molecule;
class Atom;
class Vector;

#ifdef PARALLEL
extern int thread_count;
#include <omp.h>
#include <iostream>
#endif

#include <fstream>
#include <list>
#include "Pair.hpp"
#include "SimulatedSystem.hpp"

/**
 * @defgroup PairLists Pair Lists
 * @brief Set of classes providing lists of atom pairs to interact by intermolecular interactions.
 * Intermolecular interactions describes interactions between two atoms closer than certain cutoff distance.
 * Pair list class provide a way to cycle through all relevant atomic pairs.
 * It also enable reuse of once computed interatomic distances.
 * More complex pair lists should reduce the computation complexity to something better than \f$\mathcal{O}\left(N^2\right)\f$.
 *
 * @class AbstractPairList
 * @ingroup PairLists
 * @brief Abstract class (interface) for pair lists.
 *
 * This abstract class provides an interface for pair lists.
 * The simplest pair list should be the AllPairList with complexity \f$\mathcal{0}\left(N^2\right)\f$.
 *
 * Three main functions are prescribed.
 * @ref GetNextPair() returns next pair of atoms for interaction.
 * @ref Reset() is used to set the internal initial state before calling GetNextPairs() in cycle while a negative number is returned, meaning you have reached the end of the list.
 * @ref Refresh() is used to recalculate distances and set up the list.
 *
 * The basic usage may be as follows:
 * 1. Call @ref Refresh() before forces calculation (in each step).
 * 2. Call @ref Reset() before each cycle through atomic pairs to set the iterator to the first pair.
 * 3. Call @ref GetNextPair() to get the next atomic pair until a negative number is returned, meaning you have reached the end.
 * Steps 2 and 3 are called by every single intermolecular field.
 *
 * Serves only for *normal*-space (\f$r\f$-space) calculation.
 * 
 * In the parallel version, each thread has its own vector of @ref pairs.
 * Loops over atomic pairs are then performed simply using the @ref pair_iterator which is also thread-private.
 */
class AbstractPairList
{
private:
    /**
     * @brief Default constructor of a new AbstractPairList object
     *
     * Default constructor made private in order not instantiate empty object.
     */
    AbstractPairList() : total_no_pairs(0), cutoff(0.0), sqcutoff(0.0), system(nullptr){};

protected:
    /**
     * @brief Number of intermolecular atomic pairs currently in the list.
     *
     * Total number of possible intermolecular pairs for @ref AllPairList, current number of pairs on the list otherwise.
     */
    int total_no_pairs;
    /**
     * @brief Vector of intermolecular atomic pairs (treated by one thread).
     */
    static std::vector<Pair> pairs;
    /**
     * @brief Iterator to vector @ref pairs (private for each thread).
     *
     * Stores current iterator for @ref GetNextPair().
     */
    static std::vector<Pair>::iterator pair_iterator;
#ifdef PARALLEL
#pragma omp threadprivate(pairs, pair_iterator)
#endif
    /**
     * @brief Maximum interaction cutoff
     * Pairs with larger distances can be omitted from the pair list.
     */
    double cutoff;
    /**
     * @brief Maximum interaction cutoff squared
     * Square of @ref cutoff
     */
    double sqcutoff;
    /**
     * @brief SimulatedSystem where to get pairs of Atoms.
     * The pointer to the System instance where the involved Molecules (and thus Atoms) are stored.
     */
    SimulatedSystem *system;

public:
    /**
     * @brief Constructor of an AbstractPairList object specifying the involved SimulatedSystem and maximum cutoff distance.
     *
     * @param cut Maximum interactions cutoff (see @ref cutoff)
     * @param sys SimulatedSystem (see @ref system)
     */
    AbstractPairList(SimulatedSystem *sys, double cut) : cutoff(cut), sqcutoff(cut * cut), system(sys){

                                                                                           };
    /**
     * @brief Default destructor for AbstractPairList
     * This destructor is empty, because no dynamic memory is allocated in constructors.
     * Can be overriden in classes that implement this interface (inherit this pure virtual class).
     */
    virtual ~AbstractPairList(){

    };
    /**
     * @brief Get next atomic pair
     *
     * @param pair Next pair.
     * @param sqcutoff Get the pair with squared distance less than sqcutoff.
     * @return Square of the distance between atoms of the pair. Negative number if the last pair has already been reached.
     */
    virtual double GetNextPair(Pair &pair, double sqcutoff)
    {
        while (pair_iterator != pairs.end())
        {
            if (pair_iterator->GetSqDistance() > sqcutoff)
            {
                pair_iterator++;
            }
            else
            {
                pair = *pair_iterator;
                pair_iterator++;
                return pair.GetSqDistance();
            }
        }
        // end of list reached
        return -1.0;
    };
    /**
     * @brief Reset the @ref pair_iterator to the begining of the @ref pairs list
     *
     * @return 0
     */
    int Reset()
    {
#ifdef PARALLEL
#pragma omp parallel num_threads(thread_count)
        {
            pair_iterator = pairs.begin();
        } // this works for AllPairList (tested)
#else
        pair_iterator = pairs.begin();
#endif
        return 0;
    };
    /**
     * @brief Refresh the pair list, evaluate all distances.
     * Reconstruct the entire list (if needed), store distance vectors and distances for further use.
     *
     * @return 0 on success.
     */
    virtual int Refresh() = 0;
    /**
     * @brief Update original positions for atoms in molecule if needed to measure atom displacement
     *
     * @param molno The number of the molecule to be moved
     * @param dim Dimension (0 for \f$x\f$, 1 for \f$y\f$ and 2 for \f$z\f$)
     * @param distance Move molecule by this value (multiple of actual boxsize)
     * @return 0 on success
     */
    virtual int MoveMolecule(int molno, int dim, double distance) { return 0; };
    /**
     * @brief Print informations about the pair list to the @p stream
     *
     * @param stream Output stream.
     */
    virtual void PrintInfo(std::ofstream &stream) const = 0;
};

#endif