#ifndef HARMONICANGLEHEADER
#define HARMONICANGLEHEADER

#include "AbstractIntraMolField.hpp"

class Molecule;

/**
 * @class HarmonicAngle
 * @ingroup IntraMolField
 * @brief Describes one angle with a harmonic potential.
 *
 * HarmonicAngle describes a valence angle involving atoms @ref atomI, @ref atomJ and @ref atomK, @ref atomJ being the central Atom.
 * The interaction potential energy is
 * \f{equation}{
 *     E_{\mathrm{pot}} = \frac{1}{2} k (\theta_{ijk} - \theta_0)^2
 * \f}
 * where \f$k\f$ is the force constant ( @ref k), \f$\theta_0\f$ is the equilibrium angle ( @ref angle) and
 * \f$\theta_{ijk}\f$ is the instanteneous angle defined by atoms @ref atomI, @ref atomJ and @ref atomK of the @ref parent Molecule.
 *
 * The HarmonicAngle is defined in `.field` file in an `angles` section by line with this form:
 * ```(field)
 * harm atomI atomJ atomK k angle
 * ```
 *
 * The class implements (inherits) AbstractIntraMolField interface (class), having @ref type = 1 (belonging to angles).
 */
class HarmonicAngle : public AbstractIntraMolField
{
private:
    /**
     * @brief Default constructor.
     * Made private in order not to be used because an empty angle should not be instantiated.
     */
    HarmonicAngle();

public:
    double angle; ///< equilibrium angle ([rad], in `.field` file given in degrees)
    double k;     ///< force constant ([eng. units], inside the program the programming units are used)
    int atomI;    ///< index of the first Atom
    int atomJ;    ///< index of the second (central) Atom
    int atomK;    ///< index of the third atom

    /**
     * @brief Constructor of a HarmonicAngle object.
     * Constructor uses parameters to set its members to appropriate values.
     *
     * @param parentMol Parent Molecule.
     * @param indexI First Atom index.
     * @param indexJ Second (central) Atom index.
     * @param indexK Third Atom index.
     * @param equilAngle Equilibrium angle.
     * @param angleStrength Force constant.
     */
    HarmonicAngle(Molecule *parentMol, int indexI, int indexJ, int indexK, double equilAngle, double angleStrength);
     /**
     * @brief Copy constructor.
     * Sets all the attributes to be the same as the template @p otherAngle, except for @ref parent, which is set to `nullptr`.
     * New HarmonicAngle is thus constructed as an *orphan*.
     *
     * @param otherAngle Template HarmonicAngle.
     */
    HarmonicAngle(const HarmonicAngle &otherAngle);
     /**
     * @brief Overriden assignment operator.
     * Eventually not used in `simul++`.
     *
     * @par See also
     * @ref HarmonicAngle(Molecule*, int, int, int, double, double)
     *
     *  @note This is (obviously) not a constructor. Left-hand side HarmonicAngle must be instantiated or declared before.
     *
     * @param otherAngle Template HarmonicAngle whose attributes (except for @ref parent) are copied to the left-hand-side HarmonicAngle.
     * @return Left-hand-side HarmonicAngle.
     */
    HarmonicAngle &operator=(const HarmonicAngle &otherAngle);
     /**
     * @brief Default destructor
     *
     * Empty because no memory allocation is done in constructors.
     */
    ~HarmonicAngle(){};
    /**
     * @brief Get the instanteneous angle \f$\theta_{ijk}\f$.
     * 
     * @return The instanteneous angle in radians.
     */
    double GetCurrentAngle() const;
    /**
     * @brief Calculate forces caused by this angle.
     * Calculate forces caused by this angle and add them to @ref Atom::force (of @ref atomI, @ref atomJ and @ref atomK in the parent Molecule).
     * Save the potential energy and force virial to the @ref parent Molecule @ref Molecule::Epot and @ref Molecule::virial.
     * Return the potential energy.
     *
     * The potential energy of a HarmonicAngle is given by:
     * \f{equation}{
     *     E_{\mathrm{pot}} = \frac{1}{2} k (\theta_{ijk} - \theta_0)^2
     * \f}
     * where \f$k\f$ is the force constant ( @ref k), \f$\theta_0\f$ is the equilibrium angle ( @ref angle) and
     * \f$\theta_{ijk}\f$ is the instanteneous angle defined by atoms @ref atomI, @ref atomJ and @ref atomK of the @ref parent Molecule.
     * The potential energy is added to the parent Molecule::Epot[1].
     *
     * Forces (3D vector) are given by the derivative of the potential energy. 
     * The calculation is not straightforward and details can be found in @cite dl_poly4-manual, pp 17–19.
     * In this simple case the derivation is briefly:
     * \f{equation}{
     *     \mathbf{f}_a = -\frac{\partial E_{\mathrm{pot}}}{\partial \mathbf{r}_a} = -k (\theta_{ijk} - \theta_0) \frac{\partial \theta_{ijk}}{\partial \mathbf{r}_a}
     * \f}
     * where \f$\mathbf{r}_a\f$ is the position Vector of Atom \f$i, j\f$ or \f$k\f$ (\f$a=i, j\f$ or \f$k\f$). 
     * The instanteneous angle can be calculated by
     * \f{equation}{
     *     \mathbf{\theta}_{ijk} = \arccos \left(\frac{\mathbf{r}_{ji} \cdot \mathbf{r}_{jk}}{r_{ji}r_{jk}}\right)
     * \f}
     * where
     * \f{equation}{
     *     \mathbf{r}_{ij} = \mathbf{r}_j - \mathbf{r}_i
     * \f}
     * is the difference between position vectors, whose norm is the interatomic distance (\f$\left|\left|\mathbf{r}_{ij}\right|\right| = r_{ij}\f$).
     * Since
     * \f{equation}{
     *     \cos \theta_{ijk} = \frac{\mathbf{r}_{ji} \cdot \mathbf{r}_{jk}}{r_{ji}r_{jk}}
     * \f}
     * the derivative can be calculated by the chain rule:
     * \f{equation}{
     *     \frac{\partial \theta_{ijk}}{\partial \mathbf{r}_a} = \frac{\partial \theta_{ijk}}{\partial (\cos \theta_{ijk})} \frac{\partial \left(\frac{\mathbf{r}_{ji} \cdot \mathbf{r}_{jk}}{r_{ji}r_{jk}}\right)}{\partial \mathbf{r}_a} = - \frac{1}{\sin \theta_{ijk}} \frac{\partial \left(\frac{\mathbf{r}_{ji} \cdot \mathbf{r}_{jk}}{r_{ji}r_{jk}}\right)}{\partial \mathbf{r}_a}
     * \f}
     * The most complicated part follows:
     * \f{equation}{
     *     \frac{\partial \left(\frac{\mathbf{r}_{ji} \cdot \mathbf{r}_{jk}}{r_{ji}r_{jk}}\right)}{\partial \mathbf{r}_a} = 
     *     (\delta_{ia} - \delta_{ja}) \frac{\mathbf{r}_{jk}}{r_{jk} r_{ji}} + (\delta_{ka} - \delta_{ja}) \frac{\mathbf{r}_{ji}}{r_{ji} r_{jk}} -
     *     \cos \theta_{ijk} \left[(\delta_{ia} - \delta_{ja}) \frac{\mathbf{r}_{ji}}{r_{ji}^2} + (\delta_{ka} - \delta_{ja}) \frac{\mathbf{r}_{jk}}{r_{jk}^2} \right]
     * \f}
     * which can be used to get the forces.
     *
     * The contribution to the viral by harmonic bond is zero, since the angle value remains the same when the system is shrunk or stretched.
     * Nothing is thus added to the parent Molecule::virial[1].
     *
     * @par See also
     * @ref CalculateEpot() const.
     *
     * @return Potential energy caused by this angle.
     */
    double CalculateForces() override;
    /**
     * @brief Calculate potential energy caused by this angle.
     * Similar to @ref CalculateForces(), but calculates only the potential energy and does not change anything neither in the @ref parent Molecule,
     * nor in the SimulatedSystem.
     * The potential energy is given by:
     * \f{equation}{
     *     E_{\mathrm{pot}} = \frac{1}{2} k (\theta_{ijk} - \theta_0)^2
     * \f}
     * where \f$k\f$ is the force constant ( @ref k), \f$\theta_0\f$ is the equilibrium angle ( @ref angle) and
     * \f$\theta_{ijk}\f$ is the instanteneous angle defined by atoms @ref atomI, @ref atomJ and @ref atomK of the @ref parent Molecule.
     *
     * @par See also
     * @ref CalculateForces().
     *
     * @return Potential energy caused by this harmonic angle.
     */
    double CalculateEpot() const override;                   // calculate potential energy of the HarmonicAngle
    /**
     * @brief Get the involved Atom index according to the index value.
     * Returns the label of an Atom valid inside a @ref Molecule.
     * The value of @p index is a label valid inside this HarmonicAngle (1 for @ref atomI, 2 for @ref atomJ and 3 for @ref atomK).
     *
     * @param index Label of an Atom inside this angle class.
     * @return Label of the same Atom inside the @ref parent Molecule.
     */
    int GetAtom(int index) const override;                   // get one of atoms I, J, K
    /**
     * @brief Clone itself and get adopted by another @ref parent Molecule.
     * Returns a copy of itself and migrate to the @p newparent Molecule.
     * Useful when copying Molecule object.
     * Every HarmonicAngle is then cloned by this method and assigned to the new Molecule object.
     *
     * @par See also
     * @ref HarmonicAngle(const HarmonicAngle&), @ref SetParent(Molecule*)
     *
     * @param newparent Parent Molecule of the newly created HarmonicAngle object.
     * @return (Pointer to) the newly created HarmonicAngle object.
     */
    HarmonicAngle *copy(Molecule *newparent) const override; // clone itself
    /**
     * @brief Get the angle parameters: force constant @ref k and equilibrium angle @ref angle.
     * 
     * @param fieldname Human-readable type of this field (`harm` in this case).
     * @param params Force constant and angle to be printed in info.
     * @return 2
     */
    int GetParams(std::string &fieldname, std::vector<double> &params) const override;
};

#endif