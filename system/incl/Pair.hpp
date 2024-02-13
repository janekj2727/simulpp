/*
 *  Atomic pair for use in `simul++`
 *
 *  Author: JJ
 *  Date: Dec 2022
 *
 */

#ifndef PAIRHEADER
#define PAIRHEADER

class Atom;

/**
 * @class Pair
 * @ingroup PairLists
 * @brief Intermolecular atomic pair
 *
 * An object representing an intermolecular atomic pair used for pair lists (see @ref AbstractPairList).
 * Pointers to involved atoms, a displacement vector and its squared norm are saved.
 * This structure is saved in @ref AbstractPairList::pairs and handled by pair list routines.
 */
class Pair
{
private:
    /**
     * @brief Pointer to the first Atom in the Pair
     */
    Atom *atom1;
    /**
     * @brief Pointer to the second Atom in the Pair
     */
    Atom *atom2;
    /**
     * @brief The square of the distance between @ref atom1 and @ref atom2
     */
    double sqdist;
    /**
     * @brief The displacement vector between @ref atom1 and @ref atom2
     */
    double distvec[3];

public:
    /**
     * @brief Default constructor
     *
     * No allocation is needed, only set pointers to `nullptr` and numbers to 0.0.
     */
    Pair();
    /**
     * @brief Empty destructor
     *
     * No deallocation needed.
     */
    ~Pair();
    /**
     * @brief Copy constructor
     *
     * Simple copy constructor, no allocation needed.
     * @param origpair An instance of Pair to be copied.
     */
    Pair(const Pair &origpair);
    /**
     * @brief Constructor from pointers to Atom
     * Simpler version for free boundary conditions. Distance is calculated inside.
     */
    Pair(Atom *at1, Atom *at2);
    /**
     * @brief Constructor from pointers to Atom
     * Version for periodic b. c. Distance calculated inside.
     */
    Pair(Atom *at1, Atom *at2, double *box);
    /**
     * @brief Assignment operator
     *
     * Overloading of an assignment operator=.
     * Pointers and numbers are copied to the pair on the left side of assignment.
     *
     * @note This is (obviously) not a constructor!
     *
     * @param origpair An instance of Pair to be copied
     * @return A copy of @p origpair.
     */
    Pair &operator=(const Pair &origpair);
    /**
     * @brief Calculate distance between atoms in free boundary conditions.
     *
     * Recalculate displacement vector between @ref atom1 and @ref atom2 (`atom2 - atom1`)
     * and its square (to @ref sqdist).
     * No boundary conditions applied, so it is suitable for free b. c.
     *
     * @return The squared norm of the displacement vector ( @ref sqdist).
     */
    double CalculateDistance(); // without PBC
    /**
     * @brief Calculate distance between atoms in periodic boundary conditions.
     *
     * Recalculate displacement vector between @ref atom1 and @ref atom2 (`atom2 - atom1`)
     * and its square (to @ref sqdist).
     * Periodic boundary conditions are applied, so the size of the simulation box must be supplemented.
     * The nearest image convention is thus used.
     *
     * @param box Simulation box size (`double[3]`).
     *
     * @return The squared norm of the displacement vector ( @ref sqdist).
     */
    double CalculateDistance(double *box); // with PBC
    /**
     * @brief Calculate distance between atoms in periodic boundary conditions (version 2).
     *
     * Recalculate displacement vector between @ref atom1 and @ref atom2 (`atom2 - atom1`)
     * and its square (to @ref sqdist).
     * Periodic boundary conditions are applied, so the size of the simulation box must be supplemented.
     * The nearest image convention is thus used.
     * The displacement vector is saved only if @ref sqdist < @p sqcutoff.
     *
     * @note Not eventually used, seems to be slower...
     *
     * @param box Simulation box size (`double[3]`).
     * @param sqcutoff Square of the cutoff distance [AA^2].
     *
     * @return The squared norm of the displacement vector ( @ref sqdist).
     */
    double CalculateDistance(double sqcutoff, double *box); // not assign vector if not needed
    /**
     * @brief Calculate distance between atoms in free boundary conditions (version 2).
     *
     * Recalculate displacement vector between @ref atom1 and @ref atom2 (`atom2 - atom1`)
     * and its square (to @ref sqdist).
     * No boundary conditions applied, so it is suitable for free b. c.
     * The displacement vector is saved only if @ref sqdist < @p sqcutoff.
     *
     * @note Not eventually used, seems to be slower...
     *
     * @param sqcutoff Square of the cutoff distance [AA^2].
     *
     * @return The squared norm of the displacement vector ( @ref sqdist).
     */
    double CalculateDistance(double sqcutoff);
    /**
     * @brief Set both pointer to atoms.
     *
     * Write access to the private members @ref atom1 and @ref atom2.
     *
     * @param at1 Pointer to the first atom (saved to @ref atom1)
     * @param at2 Pointer to the second atom (saved to @ref atom2)
     */
    void SetAtoms(Atom *at1, Atom *at2);
    /**
     * @brief Set square distance
     *
     * If the displacement is calculated outside the Pair object, use the result and save it to @ref sqdist.
     * @param sqdistance Squared norm of the displacement vector to be saved into @ref sqdist.
     */
    void SetDistance(double sqdistance);
    /**
     * @brief Set square distance and displacement vector.
     *
     * If the displacement is calculated outside the Pair object, use the result and save it to @ref sqdist and @ref distvec.
     * @param sqdistance Squared norm of the displacement vector to be saved into @ref sqdist.
     * @param distvector Displacement vector to be saved into @ref distvec.
     */
    void SetDistance(double sqdistance, double *distvector);
    /**
     * @brief Get @ref sqdist.
     *
     * Read access to private member @ref sqdist.
     * @return The value of @ref sqdist
     */
    inline double GetSqDistance() const
    {
        return sqdist;
    };
    /**
     * @brief Get @ref distvec entry.
     *
     * Read access to private member @ref distvec.
     * @param coordinate 0 for x, 1 for y and 2 for z coordinate of @ref distvec.
     * @return The value of @ref distvec[coordinate]
     */
    inline double GetVector(int coordinate) const
    {
        return distvec[coordinate];
    };
    /**
     * @brief Get @ref atom1.
     *
     * Read access to private member @ref atom1.
     * @return The pointer saved in @ref atom1
     */
    inline Atom *GetAtom1() const
    {
        return atom1;
    };
    /**
     * @brief Get @ref atom2.
     *
     * Read access to private member @ref atom2.
     * @return The pointer saved in @ref atom2
     */
    inline Atom *GetAtom2() const
    {
        return atom2;
    };
};

#endif
