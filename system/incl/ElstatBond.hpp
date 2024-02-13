#ifndef ELSTATBONDHEADER
#define ELSTATBONDHEADER

#include "AbstractIntraMolField.hpp"

class Molecule;

/**
 * @class ElstatBond
 * @ingroup IntraMolField
 * @brief Describes one pseudo-bond with electrostatic potential (intramolecular electrostatic interaction for distances 1–4 and higher).
 *
 * ElstatBond connects two atoms \f$i\f$ ( @ref atomI) and \f$j\f$ ( @ref atomJ) not connected by a regular bond (minimum distance 1–4).
 * The interaction potential energy is
 * \f{equation}{
 *      E_{\mathrm{elstat}} = \frac{\omega q_i q_j}{4\pi\epsilon_0}\left(\frac{1}{r_{ij}}-\mathrm{shift}\right)
 * \f}
 * where \f$r_{ij}\f$ is the instanteneous distance between atoms connected by the bond ( @ref atomI and @ref atomJ of the @ref parent Molecule)
 * and \f$q_i\f$ and \f$q_j\f$ are charges of atoms \f$i\f$ and \f$j\f$.
 * The scaling factor \f$\omega\f$ is usually equal to 1, except for 1–4 interactions.
 * The `shift` parameter is used with cutoff electrostatics (see @ref CutoffElstat).
 *
 * The pseudo-bond is instantiated to realize electrostatic interaction between atoms of the same parent Molecule, whereas electrostatic interactions normally belongs to intermolecular interactions.
 * Minimum distance inside the Molecule is calculated and for distances 1–4 and longer, the pseudo-bonds are created.
 *
 * The class implements (inherits) AbstractIntraMolField interface (class), having @ref type = 4 (belonging to electrostatic interactions).
 */
class ElstatBond : public AbstractIntraMolField
{
private:
    /**
     * @brief Default constructor.
     * Made private in order not to be used because an empty bond should not be instantiated.
     */
    ElstatBond();

public:
    double qproduct; ///< Product of charges \f$q_i q_j\f$, possibly includes scaling factor \f$\omega\f$ (for 1–4 interactions).
    double shift;    ///< Shift in energy (usefull for cutoff electrostatics)
    int atomI;       ///< Index of first Atom \f$i\f$.
    int atomJ;       ///< Index of second Atom \f$j\f$.

    /**
     * @brief Constructor of a ElstatBond object.
     * Constructor creates a new ElstatBond between atoms defined by indeces ( @p indexI and @p indexJ), possibly using scaling factor for 1–4 interactions.
     *
     * @param parentMol Parent Molecule.
     * @param indexI First bonded Atom index.
     * @param indexJ Second bonded Atom index.
     * @param sft Shift for electrostatic energy (usefull for cutoff electrostatics).
     * @param scaling Scaling factor (usefull for 1-4 interactions).
     */
    ElstatBond(Molecule *parentMol, int indexI, int indexJ, double sft = 0.0, double scaling = 1.0);
    /**
     * @brief Copy constructor.
     * Sets all the attributes to be the same as the template @p otherBond, except for @ref parent, which is set to `nullptr`.
     * New ElstatBond is thus constructed as an *orphan*.
     *
     * @param otherBond Template ElstatBond.
     */
    ElstatBond(const ElstatBond &otherBond);
    /**
     * @brief Overriden assignment operator.
     * Eventually not used in `simul++`.
     *
     * @par See also
     * @ref ElstatBond(Molecule*, int, int, double)The scaling factor \f$\omega\f$ is usually equal to 1, except for 1–4 interactions.
     * @param otherBond Template ElstatBond whose attributes (except for @ref parent) are copied to the left-hand-side ElstatBond.
     * @return Left-hand-side ElstatBond.
     */
    ElstatBond &operator=(const ElstatBond &otherBond);
    /**
     * @brief Default destructor
     *
     * Empty because no memory allocation is done in constructors.
     */
    ~ElstatBond(){};
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
     * The potential energy of a ElstatBond is given by:
     * \f{equation}{
     *     E_\mathrm{elstat} = \frac{\omega q_i q_j}{4\pi\epsilon_0}\left(\frac{1}{r_{ij}}-\mathrm{shift}\right)
     * \f}
     * where \f$r_{ij}\f$ is the instanteneous distance between atoms connected by the bond ( @ref atomI and @ref atomJ of the @ref parent Molecule)
     * and \f$q_i\f$ and \f$q_j\f$ are charges of atoms \f$i\f$ and \f$j\f$ (see @ref Atom::charge).
     * The scaling factor \f$\omega\f$ is usually equal to 1, except for 1–4 interactions.
     * The `shift` parameter is used with cutoff electrostatics (see @ref CutoffElstat).
     * The potential energy is added to the parent Molecule::Epot[4].
     *
     * Forces (3D vector) are given by the (negative) derivative of the potential energy:
     * \f{equation}{
     *     \mathbf{f}_j = -\mathbf{f}_i = \frac{\omega q_i q_j}{4\pi \epsilon_0} \frac{1}{r_{ij}^2}\frac{\mathbf{r}_{ij}}{r_{ij}}
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
     *     w_{\mathrm{LJ}} = -\frac{\omega q_i q_j}{4\pi \epsilon_0} \frac{1}{r_{ij}}
     * \f}
     * This value is added to the parent Molecule::virial[4].
     *
     * @par See also
     * @ref CalculateEpot() const.
     *
     * @return Potential energy caused by this electrostatic interaction.
     */
    double CalculateForces() override; // must save virial to parent Molecule
    /**
     * @brief Calculate potential energy caused by this LJ interaction (realized as pseudo-bond).
     * Similar to @ref CalculateForces(), but calculates only the potential energy and does not change anything neither in the @ref parent Molecule,
     * nor in the SimulatedSystem.
     * The potential energy is given by:
     * \f{equation}{
     *     E_{\mathrm{elstat}} = \frac{\omega q_i q_j}{4\pi\epsilon_0}\left(\frac{1}{r_{ij}}-\mathrm{shift}\right)
     * \f}
     * where \f$r_{ij}\f$ is the instanteneous distance between atoms connected by the bond ( @ref atomI and @ref atomJ of the @ref parent Molecule)
     * and \f$q_i\f$ and \f$q_j\f$ are charges of atoms \f$i\f$ and \f$j\f$ (see @ref Atom::charge).
     * The scaling factor \f$\omega\f$ is usually equal to 1, except for 1–4 interactions.
     * The `shift` parameter is used with cutoff electrostatics (see @ref CutoffElstat).
     *
     * @par See also
     * @ref CalculateForces().
     *
     * @return Potential energy caused by this electrostatic interaction.
     */
    double CalculateEpot() const override; // calculate potential energy of the ElstatBond
    /**
     * @brief Get the involved Atom index according to the index value.
     * Returns the label of an Atom valid inside a @ref Molecule.
     * The value of @p index is a label valid inside this ElstatBond (1 for @ref atomI and 2 for @ref atomJ).
     *
     * @param index Label of an Atom inside this bond class.
     * @return Label of the same Atom inside the @ref parent Molecule.
     */
    int GetAtom(int index) const override; // get one of connected atoms
    /**
     * @brief Clone itself and get adopted by another @ref parent Molecule.
     * Returns a copy of itself and migrate to the @p newparent Molecule.
     * Useful when copying Molecule object.
     * Every ElstatBond is then cloned by this method and assigned to the new Molecule object.
     *
     * @par See also
     * @ref ElstatBond(const ElstatBond&), @ref SetParent(Molecule*)
     *
     * @param newparent Parent Molecule of the newly created ElstatBond object.
     * @return (Pointer to) the newly created ElstatBond object.
     */
    ElstatBond *copy(Molecule *newparent) const override; // clone itself
    /**
     * @brief Get the bond parameters: product of charges (including scaling parameter, @ref qproduct) and energy @ref shift
     * 
     * @param fieldname Human-readable type of this field (`elstat` in this case).
     * @param params Parameters of electrostatic bond to be printed in info.
     * @return 2
     */
    int GetParams(std::string &fieldname, std::vector<double> &params) const override;
};

#endif