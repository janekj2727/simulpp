#ifndef VERLETLISTHEADER
#define VERLETLISTHEADER

#include "AbstractPairList.hpp"
using arrvec3 = double (*)[3];

/**
 * @class VerletList
 * @ingroup PairLists
 * @brief The implementation of Verlet list book-keeping method in its simplest form.
 *
 * Verlet list makes use of the fact, that only those intermolecular atom pairs are relevant,
 * whose distance is shorter than the maximum cutoff (see @ref cutoff).
 * Furthermore, in molecular dynamics, atoms move only a little during one step.
 * Hence, only pairs that are not further away than the @ref cutoff plus certain save margin
 * (see @ref outercut) can affect energy in the next few steps.
 * Consequently, the list of pairs within @ref outercut is made and used for energy and forces
 * calculation until the displacement of any Atom is bigger then half the save margin (see @ref maxdisp2), i.e. one half of @ref outercut − @ref cutoff.
 * Then, a new list is made from scratch by @ref MakeList().
 * Atom displacements are calculated during each @ref Refresh() call (method @ref MakeList() is called when needed).
 *
 * Number of pairs is changing during the simulation and can be monitored by adding
 * ```(control)
 * measure nopairs
 * ```
 * to the `.control` file. The Verlet list is activated using
 * ```(control)
 * pairlist verlet SAVE_MARGIN
 * ```
 * `.control` directive, where `SAVE_MARGIN` determines the value of @ref outercut,
 * which is computed by @ref outercut = `SAVE_MARGIN` × @ref cutoff, @ref cutoff being the maximum cutoff distance used by any InterMolecularField.
 * 
 * In the parallel version, each thread constructs its own vector of @ref pairs during @ref MakeList().
 * Then, a load-balancing is done, so that each thread has the same number of @ref pairs.
 */
class VerletList : public AbstractPairList
{
private:
    /**
     * @brief Default constructor of a new VerletList object
     *
     * Default constructor made private in order not instantiate empty object.
     */
    VerletList();
    /**
     * @brief Total number of atoms in the system
     */
    int no_atoms;
    /**
     * @brief Total number of molecules in the system
     */
    int no_molecules;
    /**
     * @brief (Original) positions of atoms when the list was made
     * Saves the original positions to enable calculation of displacement, which determines if the VerletList must be recomputed or not (if @ref MakeList() must be called).
     */
    arrvec3 origpos;
    /**
     * @brief Displacement of @ref origpos due to the periodic boundary conditions
     * When molecule is **moved** due to the periodic boundary conditions, the **original position** saved in @ref origpos should be updated.
     * No boundary checking is then needed when calculating molecular displacements.
     */
    arrvec3 PBCdisp;
    /**
     * @brief Outer (Verlet) cuttof determines if the pair should be on the list
     * Outer cutoff is slightly larger then @ref cutoff to include more pairs (even those who can fall into inner @ref cutoff in the next few steps).
     * The ratio of @ref cutoff and outer cutoff affects the efficiency of the algorithm.
     */
    double outercut;
    /**
     * @brief Outer cutoff squared (see @ref outercut)
     */
    double outercut2;
    /**
     * @brief Squared maximum displacement of one atom before two list recalculations
     * When the displacement of any atom is larger, the whole list is recalculated by @ref MakeList().
     */
    double maxdisp2;
    /**
     * @brief Make a new list from scratch.
     * Every time any atom crosses the maximum displacement (see @ref maxdisp2), a new list is made from scratch.
     * This step has a \f$\mathcal{O}\left(N_{\mathrm{atoms}}^2\right)\f$ complexity.
     * @return 0 on success
     */
    int MakeList();

public:
    /**
     * @brief Constructor of an VerletList object specifying the involved SimulatedSystem, maximum cutoff distance and @ref outercut.
     *
     * @param cut Maximum interactions cutoff (see @ref cutoff)
     * @param sys SimulatedSystem (see @ref system)
     * @param savedist Determines the *save margin*. @ref outercut is calculated by @p savedist × @p cut.
     */
    VerletList(SimulatedSystem *sys, double cut, double savedist);
    /**
     * @brief Destroy the VerletList object
     *
     * Deallocates memory used by this pair list.
     */
    ~VerletList();
    /**
     * @brief Refresh the pair list, evaluate all distances.
     * Reconstruct the entire list (if needed, see @ref MakeList()), store distance vectors and distances for further use. Distance vectors are stored only if distance is shorter then @ref cutoff for the given pair.
     *
     * @return 0 on success.
     */
    int Refresh();
    /**
     * @brief Update @ref origpos for atoms in molecule
     * Useful when applying periodic boundary conditions.
     *
     * @param molno The number of the molecule to be moved
     * @param dim Dimension (0 for \f$x\f$, 1 for \f$y\f$ and 2 for \f$z\f$)
     * @param distance Move molecule by this value (multiple of actual boxsize)
     * @return 0 on success
     */
    int MoveMolecule(int molno, int dim, double distance) override;
    /**
     * @brief Print informations about the pair list to the @p stream
     *
     * @param stream Output stream.
     */
    void PrintInfo(std::ofstream &stream) const;
};

#endif