#ifndef BONDHEADER
#define BONDHEADER

#include "AbstractIntraMolField.hpp"

class Molecule;

/**
 * @class HarmonicBond
 * @ingroup IntraMolField
 * @brief Describes one bond with a harmonic potential.
 *
 * HarmonicBond connects two atoms by a harmonic spring.
 * The interaction potential energy is
 * \f{equation}{
 *     E_{\mathrm{pot}} = \frac{1}{2} k (r_{ij} - r_0)^2
 * \f}
 * where \f$k\f$ is the force constant ( @ref k), \f$r_0\f$ is the equilibrium bond length ( @ref length) and
 * \f$r_{ij}\f$ is the instanteneous distance between atoms connected by the bond ( @ref atomI and @ref atomJ of the @ref parent Molecule).
 *
 * The HarmonicBond is defined in `.field` file in a `bonds` section by line with this form:
 * ```(field)
 * harm atomI atomJ k length
 * ```
 *
 * The class implements (inherits) AbstractIntraMolField interface (class), having @ref type = 0 (belonging to bonds).
 */
class HarmonicBond : public AbstractIntraMolField
{
private:
    /**
     * @brief Default constructor.
     * Made private in order not to be used because an empty bond should not be instantiated.
     */
    HarmonicBond();

public:
    double length; ///< equilibrium length ([AA])
    double k;      ///< force constant (originally in [eng. unit./AA^2], energy units are given in the `.field` file header); inside the program uses program units
    int atomI;     ///< index of first bonded Atom
    int atomJ;     ///< index of second bonded Atom

    /**
     * @brief Constructor of a HarmonicBond object.
     * Constructor uses parameters to set its members to appropriate values.
     *
     * @param parentMol Parent Molecule.
     * @param indexI First bonded Atom index
     * @param indexJ Second bonded Atom index
     * @param bondLength Equilibrium bond length.
     * @param bondStrength HarmonicBond force constant.
     */
    HarmonicBond(Molecule *parentMol, int indexI, int indexJ, double bondLength, double bondStrength);
    /**
     * @brief Copy constructor.
     * Sets all the attributes to be the same as the template @p otherBond, except for @ref parent, which is set to `nullptr`.
     * New HarmonicBond is thus constructed as an *orphan*.
     *
     * @param otherBond Template HarmonicBond.
     */
    HarmonicBond(const HarmonicBond &otherBond);
    /**
     * @brief Overriden assignment operator.
     * Eventually not used in `simul++`.
     *
     * @par See also
     * @ref HarmonicBond(Molecule*, int, int, double, double)
     *
     *  @note This is (obviously) not a constructor. Left-hand side HarmonicBond must be instantiated or declared before.
     *
     * @param otherBond Template HarmonicBond whose attributes (except for @ref parent) are copied to the left-hand-side HarmonicBond.
     * @return Left-hand-side HarmonicBond.
     */
    HarmonicBond &operator=(const HarmonicBond &otherBond);
    /**
     * @brief Default destructor
     *
     * Empty because no memory allocation is done in constructors.
     */
    ~HarmonicBond(){};
    /**
     * @brief Get the current bond length.
     * Returns the current distance between the Atoms connected by this bond.
     *
     * @return Current bond length [AA].
     */
    double GetCurrentLength() const;
    /**
     * @brief Calculate forces caused by this bond.
     * Calculate forces caused by this bond and add them to @ref Atom::force (of @ref atomI and @ref atomJ).
     * Save the potential energy and force virial to the @ref parent Molecule @ref Molecule::Epot and @ref Molecule::virial.
     * Return the potential energy.
     *
     * The potential energy of a HarmonicBond is given by:
     * \f{equation}{
     *     E_{\mathrm{pot}} = \frac{1}{2} k (r_{ij} - r_0)^2
     * \f}
     * where \f$k\f$ is the force constant ( @ref k), \f$r_0\f$ is the equilibrium bond length ( @ref length) and
     * \f$r_{ij}\f$ is the instanteneous distance between atoms connected by the bond ( @ref atomI and @ref atomJ of the @ref parent Molecule).
     * The potential energy is added to the parent Molecule::Epot[0].
     *
     * Forces (3D vector) are given by the derivative of the potential energy:
     * \f{equation}{
     *     \mathbf{f}_i = -\mathbf{f}_j = k (r_{ij} - r_0) \frac{\mathbf{r}_{ij}}{r_{ij}}
     * \f}
     * where \f$\mathbf{r}_{ij}\f$ is the Vector between bonded Atoms defined by:
     * \f{equation}{
     *     \mathbf{r}_{ij} = \mathbf{r}_j - \mathbf{r}_i
     * \f}
     * where \f$\mathbf{r}_i\f$ is the position Vector of Atom @ref atomI and vice versa.
     * The norm of the Vector \f$\mathbf{r}_{ij}\f$ is the instanteneous interatomic distance \f$r_{ij}\f$.
     *
     * The contribution to the viral by this bond is:
     * \f{equation}{
     *     w_{\mathrm{bond}} = k (r_{ij} - r_0) r_{ij}
     * \f}
     * This value is added to the parent Molecule::virial[0].
     *
     * @par See also
     * @ref CalculateEpot() const.
     *
     * @return Potential energy caused by this bond.
     */
    double CalculateForces() override; // must save virial to parent Molecule
    /**
     * @brief Calculate potential energy caused by this bond.
     * Similar to @ref CalculateForces(), but calculates only the potential energy and does not change anything neither in the @ref parent Molecule,
     * nor in the SimulatedSystem.
     * The potential energy is given by:
     * \f{equation}{
     *     E_{\mathrm{pot}} = \frac{1}{2} k (r_{ij} - r_0)^2
     * \f}
     * where \f$k\f$ is the force constant ( @ref k), \f$r_0\f$ is the equilibrium bond length ( @ref length) and
     * \f$r_{ij}\f$ is the instanteneous distance between atoms connected by the bond ( @ref atomI and @ref atomJ of the @ref parent Molecule).
     *
     * @par See also
     * @ref CalculateForces().
     *
     * @return Potential energy caused by this harmonic bond.
     */
    double CalculateEpot() const override; // calculate potential energy of the HarmonicBond
    /**
     * @brief Get the bonded Atom index according to the index value.
     * Returns the label of an Atom valid inside a @ref Molecule.
     * The value of @p index is a label valid inside this HarmonicBond (1 for @ref atomI and 2 for @ref atomJ).
     *
     * @param index Label of an Atom inside this bond class.
     * @return Label of the same Atom inside the @ref parent Molecule.
     */
    int GetAtom(int index) const override; // get one of connected atoms
    /**
     * @brief Clone itself and get adopted by another @ref parent Molecule.
     * Returns a copy of itself and migrate to the @p newparent Molecule.
     * Useful when copying Molecule object.
     * Every HarmonicBond is then cloned by this method and assigned to the new Molecule object.
     *
     * @par See also
     * @ref HarmonicBond(const HarmonicBond&), @ref SetParent(Molecule*)
     *
     * @param newparent Parent Molecule of the newly created HarmonicBond object.
     * @return (Pointer to) the newly created HarmonicBond object.
     */
    HarmonicBond *copy(Molecule *newparent) const override; // clone itself
    /**
     * @brief Get the bond parameters: force constant and length
     * 
     * @param fieldname Human-readable type of this field (`harm` in this case).
     * @param params Force constant and length to be printed in info.
     * @return 2
     */
    int GetParams(std::string &fieldname, std::vector<double> &params) const override;
};

#endif