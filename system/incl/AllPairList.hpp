#ifndef ALLPAIRLISTHEADER
#define ALLPAIRLISTHEADER

#include "AbstractPairList.hpp"

/**
 * @class AllPairList
 * @ingroup PairLists
 * @brief The simplest implementation of the AbstractPairList interface.
 *
 * The simplest pair list with complexity \f$\mathcal{O}\left(N_\mathrm{atoms}^2\right)\f$.
 * Each @ref Refresh(), all possible pairs are tryied and (squared) distance computed.
 * Number of pairs is calculated using the *old*-style cycle through molecules and atoms.
 * Pointers to atoms are saved to @ref pairs once in the constructor.
 * 3D vector is saved only if the interatomic distance for the given pair is shorter than @ref cutoff.
 * Should be fool-proof and can be used for testing.
 * For computation, @ref VerletList is usually faster and thus recommended.
 *
 * AllPairList is used by default in `simul++`, but can be demanded explicitly in `.control` file by:
 * ```(control)
 * pairlist all
 * ```
 * 
 * In the parallel version, the vector of @ref pairs is distributed evenly over used threads.
 * (Initially a vector `global_pairs` is created from which the the @ref pairs are creating by copying.)
 */
class AllPairList : public AbstractPairList
{
private:
    /**
     * @brief Default constructor of a new AllPairList object
     *
     * Default constructor made private in order not instantiate empty object.
     */
    AllPairList();

public:
    /**
     * @brief Constructor of an AllPairList object specifying the involved SimulatedSystem and maximum cutoff distance.
     *
     * @param cut Maximum interactions cutoff (see @ref cutoff)
     * @param sys SimulatedSystem (see @ref system)
     */
    AllPairList(SimulatedSystem *sys, double cut);
    /**
     * @brief Destroy the AllPairList object
     *
     * Deallocates memory used by this pair list.
     */
    ~AllPairList();
    /**
     * @brief Refresh the pair list, evaluate all distances.
     * Store distances for further use. Distance vectors are saved only for pairs for which
     * distance < cutoff.
     *
     * @return 0 on success.
     */
    int Refresh() override;
    /**
     * @brief Print informations about the pair list to the @p stream
     *
     * @param stream Output stream.
     */
    void PrintInfo(std::ofstream &stream) const;
};

#endif