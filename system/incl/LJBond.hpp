#ifndef LJBONDHEADER
#define LJBONDHEADER

#include "AbstractIntraMolField.hpp"

class Molecule;

/**
 * @class LJBond
 * @ingroup IntraMolField
 * @brief Describes one pseudo-bond with Lennard-Jones (12–6) potential (intramolecular disperse interaction for distances 1–4 and higher).
 *
 * LJBond connects two atoms \f$i\f$ ( @ref atomI) and \f$j\f$ ( @ref atomJ) not connected by a regular bond (minimum distance 1–4).
 * The interaction potential energy is
 * \f{equation}{
 *      E_{\mathrm{LJ}} = 4 \varepsilon \left[\left(\frac{\sigma}{r_{ij}}\right)^{12} - \left(\frac{\sigma}{r_{ij}}\right)^6\right]
 * \f}
 * where \f$r_{ij}\f$ is the instanteneous distance between atoms connected by the bond ( @ref atomI and @ref atomJ of the @ref parent Molecule)
 * and \f$\varepsilon\f$ and \f$\sigma\f$ are Lennard-Jones parameters (see @ref epsilon and @ref sigma).
 *
 * The pseudo-bond is instantiated to realize disperse interactions between atoms of the same parent Molecule, whereas disperse interactions normally belongs to intermolecular interactions.
 * Minimum distance inside the Molecule is calculated and for distances 1–4 and longer, the pseudo-bonds are created.
 *
 * The class implements (inherits) AbstractIntraMolField interface (class), having @ref type = 3 (belonging to disperse interactions).
 */
class LJBond : public AbstractIntraMolField
{
private:
    /**
     * @brief Default constructor.
     * Made private in order not to be used because an empty bond should not be instantiated.
     */
    LJBond();

public:
    double sigma;   ///< LJ \f$\sigma\f$ parameter ([AA])
    double epsilon; ///< LJ \f$\varepsilon\f$ parameter (originally in [eng. unit./AA^2], energy units are given in the `.field` file header); inside the program uses program units; for 1–4 interactions usually scaled down
    double sigma6;  ///< The 6th power of LJ length parameter (\f$\sigma^6\f$).
    int atomI;      ///< index of first bonded Atom
    int atomJ;      ///< index of second bonded Atom

    /**
     * @brief Constructor of a LJBond object.
     * Constructor uses parameters to set its members to appropriate values.
     *
     * @param parentMol Parent Molecule.
     * @param indexI First bonded Atom index
     * @param indexJ Second bonded Atom index
     * @param sig LJ \f$\sigma\f$ parameter.
     * @param eps LJ \f$\varepsilon\f$ parameter.
     */
    LJBond(Molecule *parentMol, int indexI, int indexJ, double sig, double eps);
    /**
     * @brief Copy constructor.
     * Sets all the attributes to be the same as the template @p otherBond, except for @ref parent, which is set to `nullptr`.
     * New LJBond is thus constructed as an *orphan*.
     *
     * @param otherBond Template LJBond.
     */
    LJBond(const LJBond &otherBond);
    /**
     * @brief Overriden assignment operator.
     * Eventually not used in `simul++`.
     *
     * @par See also
     * @ref LJBond(Molecule*, int, int, double, double)
     *
     *  @note This is (obviously) not a constructor. Left-hand side LJBond must be instantiated or declared before.
     *
     * @param otherBond Template LJBond whose attributes (except for @ref parent) are copied to the left-hand-side LJBond.
     * @return Left-hand-side LJBond.
     */
    LJBond &operator=(const LJBond &otherBond);
    /**
     * @brief Default destructor
     *
     * Empty because no memory allocation is done in constructors.
     */
    ~LJBond(){};
    /**
     * @brief Get the current distance between @ref atomI and @ref atomJ.
     * Returns the current distance between the Atoms connected by this pseudo-bond.
     *
     * @return Current distance between the Atoms \f$i\f$ and \f$j\f$ [AA].
     */
    double GetCurrentLength() const;
    /**
     * @brief Calculate forces caused by this bond.
     * Calculate forces caused by this bond and add them to @ref Atom::force (of @ref atomI and @ref atomJ).
     * Save the potential energy and force virial to the @ref parent Molecule @ref Molecule::Epot and @ref Molecule::virial.
     * Return the potential energy.
     *
     * The potential energy of a LJBond is given by:
     * \f{equation}{
     *     E_{\mathrm{LJ}} = 4 \varepsilon \left[\left(\frac{\sigma}{r_{ij}}\right)^{12} - \left(\frac{\sigma}{r_{ij}}\right)^6\right]
     * \f}
     * where \f$r_{ij}\f$ is the instanteneous distance between atoms connected by the bond ( @ref atomI and @ref atomJ of the @ref parent Molecule)
     * and \f$\varepsilon\f$ and \f$\sigma\f$ are LJ parameters (see @ref epsilon, @ref sigma).
     * The potential energy is added to the parent Molecule::Epot[3].
     *
     * Forces (3D vector) are given by the (negative) derivative of the potential energy:
     * \f{equation}{
     *     \mathbf{f}_i = -\mathbf{f}_j = 4 \varepsilon \left(\frac{12\sigma^{12}}{r_{ij}^{13}} -
     *                                      \frac{6\sigma^6}{r_{ij}^7}\right)\frac{\mathbf{r}_{ij}}{r_{ij}}
     * \f}
     * where \f$\mathbf{r}_{ij}\f$ is the Vector between Atoms \f$i\f$ and \f$j\f$ defined by:
     * \f{equation}{
     *     \mathbf{r}_{ij} = \mathbf{r}_j - \mathbf{r}_i
     * \f}
     * where \f$\mathbf{r}_i\f$ is the position Vector of Atom @ref atomI and vice versa.
     * The norm of the Vector \f$\mathbf{r}_{ij}\f$ is the instanteneous interatomic distance \f$r_{ij}\f$.
     *
     * The contribution to the viral by this bond is:
     * \f{equation}{
     *     w_{\mathrm{LJ}} = -4 \varepsilon \left[\left(\frac{12\sigma^{12}}{r_{ij}^{13}}\right) -
     *                       \left(\frac{6\sigma^6}{r_{ij}^7}\right)\right]r_{ij} 
     * \f}
     * This value is added to the parent Molecule::virial[3].
     *
     * @par See also
     * @ref CalculateEpot() const.
     *
     * @return Potential energy caused by this LJ interaction.
     */
    double CalculateForces() override; // must save virial to parent Molecule
    /**
     * @brief Calculate potential energy caused by this LJ interaction (realized as pseudo-bond).
     * Similar to @ref CalculateForces(), but calculates only the potential energy and does not change anything neither in the @ref parent Molecule,
     * nor in the SimulatedSystem.
     * The potential energy is given by:
     * \f{equation}{
     *     E_{\mathrm{pot}} = 4 \varepsilon \left[\left(\frac{\sigma}{r_{ij}}\right)^{12} - \left(\frac{\sigma}{r_{ij}}\right)^6\right]
     * \f}
     * where \f$r_{ij}\f$ is the instanteneous distance between atoms connected by the bond ( @ref atomI and @ref atomJ of the @ref parent Molecule)
     * and \f$\varepsilon\f$ and \f$\sigma\f$ are LJ parameters (see @ref epsilon, @ref sigma).
     *
     * @par See also
     * @ref CalculateForces().
     *
     * @return Potential energy caused by this pair intramolecular LJ interaction.
     */
    double CalculateEpot() const override; // calculate potential energy of the LJBond
    /**
     * @brief Get the involved Atom index according to the index value.
     * Returns the label of an Atom valid inside a @ref Molecule.
     * The value of @p index is a label valid inside this LJBond (1 for @ref atomI and 2 for @ref atomJ).
     *
     * @param index Label of an Atom inside this bond class.
     * @return Label of the same Atom inside the @ref parent Molecule.
     */
    int GetAtom(int index) const override; // get one of connected atoms
    /**
     * @brief Clone itself and get adopted by another @ref parent Molecule.
     * Returns a copy of itself and migrate to the @p newparent Molecule.
     * Useful when copying Molecule object.
     * Every LJBond is then cloned by this method and assigned to the new Molecule object.
     *
     * @par See also
     * @ref LJBond(const LJBond&), @ref SetParent(Molecule*)
     *
     * @param newparent Parent Molecule of the newly created LJBond object.
     * @return (Pointer to) the newly created LJBond object.
     */
    LJBond *copy(Molecule *newparent) const override; // clone itself
    /**
     * @brief Get the bond parameters: @ref epsilon and @ref sigma 
     * 
     * @param fieldname Human-readable type of the field (`lj` in this case)
     * @param params Parameters of LJBond to be printed in info.
     * @return 2
     */
    int GetParams(std::string &fieldname, std::vector<double> &params) const override;
};

#endif