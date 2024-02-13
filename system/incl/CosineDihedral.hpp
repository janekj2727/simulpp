#ifndef COSINEDIHEDRALHEADER
#define COSINEDIHEDRALHEADER

#include "AbstractIntraMolField.hpp"

class Molecule;

/**
 * @class CosineDihedral
 * @ingroup IntraMolField
 * @brief Describes one dihedral with a cosine potential.
 *
 * CosineDihedral describes a dihedral angle involving atoms @ref atomI, @ref atomJ, @ref atomK and @ref atomL.
 * @image html dihedrals.svg
 * The interaction potential energy is
 * \f{equation}{
 *     E_{\mathrm{pot}} = A \left[ 1 + \cos\left(m\varphi_{ijkl} - \delta \right)\right]
 * \f}
 * where \f$A\f$ is the force constant ( @ref A), \f$m\f$ is an integer number ( @ref m) and
 * \f$\delta\f$ is the angle difference ( @ref delta).
 * The instanteneous dihedral angle \f$\varphi_{ijkl}\f$ is given by
 * \f{align}{
 *     B_{ijkl} &= \frac{\left(\mathbf{r}_{ij} \times \mathbf{r}_{jk}\right)\left(\mathbf{r}_{jk} \times \mathbf{r}_{kl}\right)}
 *                 {\left|\left|\mathbf{r}_{ij} \times \mathbf{r}_{jk}\right|\right|\left|\left|\mathbf{r}_{jk} \times \mathbf{r}_{kl}\right|\right|} \\
 *     \varphi_{ijkl} &= \arccos\left(B_{ijkl}\right) \\
 * \f}
 * where \f$\mathbf{r}_{ij} = \mathbf{r}_j - \mathbf{r}_i\f$ (\f$\mathbf{r}_i\f$ is the position Vector of @ref atomI of the @ref parent Molecule), etc.
 * The Vector cross-product is marked by \f$\times\f$ and the Euclidean norm (2-norm) of Vector \f$\mathbf{r}_i\f$ is denoted \f$\left|\left| \mathbf{r}_i \right|\right| = r_i\f$ for short.
 *
 * The CosineDihedral is defined in `.field` file in the `dihedrals` section by line with this form:
 * ```(field)
 * cos atomI atomJ atomK atomL A delta m
 * ```
 *
 * The class implements (inherits) AbstractIntraMolField interface (class), having @ref type = 2 (belonging to dihedrals).
 * 
 * **Improper torsion** can also be described by this dihedral type, the @ref type is then 5 (torsions). 
 * This affects only info printing and connectivity check. 
 * The energy and forces are the same as if it was a normal dihedral.
 * Impropers are marked by an asterisk ('\*') in `.field` file:
 * ```(field)
 * cos* atomI atomJ atomK atomL A delta m
 * ```
 */
class CosineDihedral : public AbstractIntraMolField
{
private:
    /**
     * @brief Default constructor.
     * Made private in order not to be used because an empty dihedral should not be instantiated.
     */
    CosineDihedral();

public:
    double delta; ///< angle difference ([deg], inside the program in [rad])
    double A;     ///< force constant ([eng. units], inside the program the programming units are used)
    int m;        ///< *multiplicity* of the dihedral (e.g. 3 in conformations of ethane)
    int atomI;    ///< index of the first Atom
    int atomJ;    ///< index of the second Atom
    int atomK;    ///< index of the third Atom
    int atomL;    ///< index of the fourth Atom

    /**
     * @brief Constructor of a CosineDihedral object.
     * Constructor uses parameters to set its members to appropriate values.
     *
     * @param parentMol Parent Molecule.
     * @param indexI First Atom index.
     * @param indexJ Second Atom index.
     * @param indexK Third Atom index.
     * @param indexL Fourth Atom index.
     * @param d Angle difference.
     * @param mult Dihedral *multiplicity*.
     * @param Astrength Force constant.
     * @param proper Proper or improper dihedral (improper disables connectivity checks).
     */
    CosineDihedral(Molecule *parentMol, int indexI, int indexJ, int indexK, int indexL, double Astrength, double d, int mult, bool proper);
    /**
     * @brief Copy constructor.
     * Sets all the attributes to be the same as the template @p otherDihedral, except for @ref parent, which is set to `nullptr`.
     * New CosineDihedral is thus constructed as an *orphan*.
     *
     * @param otherDihedral Template CosineDihedral.
     */
    CosineDihedral(const CosineDihedral &otherDihedral);
    /**
     * @brief Overriden assignment operator.
     * Eventually not used in `simul++`.
     *
     * @par See also
     * @ref CosineDihedral(Molecule*, int, int, int, int, double, double, int)
     *
     *  @note This is (obviously) not a constructor. Left-hand side CosineDihedral must be instantiated or declared before.
     *
     * @param otherDihedral Template CosineDihedral whose attributes (except for @ref parent) are copied to the left-hand-side CosineDihedral.
     * @return Left-hand-side CosineDihedral.
     */
    CosineDihedral &operator=(const CosineDihedral &otherDihedral);
    /**
     * @brief Default destructor
     *
     * Empty because no memory allocation is done in constructors.
     */
    ~CosineDihedral(){};
    /**
     * @brief Get the instanteneous dihedral angle \f$\varphi_{ijkl}\f$.
     * The instanteneous dihedral angle \f$\varphi_{ijkl}\f$ is calculated from
     * \f{align}{
     *     \cos\varphi_{ijkl} = B_{ijkl} &= \frac{\left(\mathbf{r}_{ij} \times \mathbf{r}_{jk}\right)\left(\mathbf{r}_{jk} \times \mathbf{r}_{kl}\right)}
     *                 {\left|\left|\mathbf{r}_{ij} \times \mathbf{r}_{jk}\right|\right|\left|\left|\mathbf{r}_{jk} \times \mathbf{r}_{kl}\right|\right|} \\
     *      \sin\varphi_{ijkl} &= \frac{\mathbf{r}_{jk}\cdot \left(\mathbf{r}_b \times \mathbf{r}_c\right)}
     *     {\left|\left|\mathbf{r}_{jk}\right|\right|\left|\left|\mathbf{r}_{b}\right|\right|\left|\left|\mathbf{r}_{c}\right|\right|} \\
     * \f}
     * The `atan2` function (with arguments `sinus` and `cosinus`) is used to determine the angle in the appropriate quadrant.
     *
     * @return The instanteneous dihedral angle in radians \f$[-\pi,\pi]\f$.
     */
    double GetCurrentAngle() const;
    /**
     * @brief Calculate forces caused by this dihedral.
     * Calculate forces caused by this angle and add them to @ref Atom::force (of @ref atomI, @ref atomJ, @ref atomK and @ref atomL in the parent Molecule).
     * Save the potential energy and force virial to the @ref parent Molecule @ref Molecule::Epot and @ref Molecule::virial.
     * Return the potential energy.
     *
     * The potential energy of a CosineDihedral is given by:
     * \f{equation}{
     *      E_{\mathrm{pot}} = A \left[ 1 + \cos\left(m\varphi_{ijkl} - \delta \right)\right]
     * \f}
     * where \f$A\f$ is the force constant ( @ref A), \f$\delta\f$ is the angle difference (shift) ( @ref delta), \f$m\f$ is the dihedral *multiplicity* (number of energy minima, @ref m) and
     * \f$\varphi_{ijkl}\f$ is the instanteneous angle defined by atoms @ref atomI, @ref atomJ, @ref atomK @ref atomL of the @ref parent Molecule.
     * The potential energy is added to the parent Molecule::Epot[1].
     *
     * Forces (3D vectors) are given by the negative derivative of the potential energy.
     * The calculation is not straightforward and details can be found in @cite dl_poly4-manual, pp 20–22.
     * In this case the derivation is briefly:
     * \f{equation}{
     *     f_a^\alpha = -\frac{\partial E_{\mathrm{pot}}}{\partial r_a^\alpha} =
     *                -\frac{\partial E_\mathrm{pot}}{\partial \varphi_{ijkl}}\frac{\partial \varphi_{ijkl}}{\partial B_{ijkl}}\frac{\partial B_{ijkl}}{\partial r_a^\alpha}
     * \f}
     * where \f$r_a^\alpha\f$ is the position coordinate of Atom \f$i, j\f$ or \f$k\f$ (\f$a=i, j\f$ or \f$k\f$; \f$\alpha=x, y\f$ or \f$z\f$).
     * The instanteneous dihedral angle can be calculated by
     * \f{align}{
     *     B_{ijkl} &= \frac{\left(\mathbf{r}_{ij} \times \mathbf{r}_{jk}\right)\left(\mathbf{r}_{jk} \times \mathbf{r}_{kl}\right)}
     *                 {\left|\left|\mathbf{r}_{ij} \times \mathbf{r}_{jk}\right|\right|\left|\left|\mathbf{r}_{jk} \times \mathbf{r}_{kl}\right|\right|} \\
     *     \varphi_{ijkl} &= \arccos\left(B_{ijkl}\right) \\
     * \f}
     * where
     * \f{equation}{
     *     \mathbf{r}_{ij} = \mathbf{r}_j - \mathbf{r}_i
     * \f}
     * is the difference between position vectors, whose norm is the interatomic distance (\f$\left|\left|\mathbf{r}_{ij}\right|\right| = r_{ij}\f$).
     *
     * The first two parts of the chain rule can be readily calculated
     * \f{align}{
     *       \frac{\partial E_{\mathrm{pot}}}{\partial \varphi_{ijkl}} &= -Am\sin\left(m\varphi_{ijkl}-\delta\right) \\
     *       \frac{\partial \varphi_{ijkl}}{\partial B_{ijkl}} &= -\frac{1}{\sin\varphi_{ijkl}} \\
     * \f}
     * The most complicated is the third partial derivative.
     * Let \f$\mathbf{r}_{b} = \mathbf{r}_{ij} \times \mathbf{r}_{jk}\f$ and \f$\mathbf{r}_{c} = \mathbf{r}_{jk} \times \mathbf{r}_{kl}\f$.
     * Then \f$B_{ijkl}\f$ can be written as
     * \f{equation}{
     *     B_{ijkl} = \frac{\mathbf{r}_b \cdot \mathbf{r}_c}{\left|\left|\mathbf{r}_b\right|\right| \left|\left|\mathbf{r}_c\right|\right|}
     * \f}
     * and the derivative is
     * \f{equation}{
     *     \frac{\partial B_{ijkl}}{\partial r_a^\alpha} =
     *         \frac{\left(\frac{\partial (\mathbf{r}_b \cdot \mathbf{r}_c)}{\partial r_a^\alpha}\right) \left|\left|\mathbf{r}_b\right|\right| \left|\left|\mathbf{r}_c\right|\right|
     *         - \mathbf{r}_b \cdot \mathbf{r}_c \left(\frac{\partial (\left|\left|\mathbf{r}_b\right|\right| \left|\left|\mathbf{r}_c\right|\right|)}{\partial r_a^\alpha}\right)}
     *         {\left|\left|\mathbf{r}_b\right|\right|^2 \left|\left|\mathbf{r}_c\right|\right|^2} = 
     *    \frac{\left(\frac{\partial (\mathbf{r}_b \cdot \mathbf{r}_c)}{\partial r_a^\alpha}\right)}
     *         {\left|\left|\mathbf{r}_b\right|\right| \left|\left|\mathbf{r}_c\right|\right|} -
     *    B_{ijkl} \frac{\left(\frac{\partial (\left|\left|\mathbf{r}_b\right|\right| \left|\left|\mathbf{r}_c\right|\right|)}{\partial r_a^\alpha}\right)}
     *         {\left|\left|\mathbf{r}_b\right|\right| \left|\left|\mathbf{r}_c\right|\right|}
     * \label{eq:B_diff}
     * \f}
     * The first derivative in the numerator expands to
     * \f{equation}{
     *     \frac{\partial (\mathbf{r}_b \cdot \mathbf{r}_c)}{\partial r_a^\alpha} =
     *     \frac{\partial \mathbf{r}_b}{\partial r_a^\alpha} \mathbf{r}_c + \mathbf{r}_b \frac{\partial \mathbf{r}_c}{\partial r_a^\alpha}
     * \f}
     * Each derivative is a derivative of a vector product which further expands to
     * \f{equation}{
     *     \frac{\partial \mathbf{r}_b}{\partial r_a^\alpha} = \frac{\partial \mathbf{r}_{ij} \times \mathbf{r}_{jk}}{\partial r_a^\alpha} =
     *     \frac{\partial \mathbf{r}_{ij}}{\partial r_a^\alpha} \times \mathbf{r}_{jk} + \mathbf{r}_{ij}\times\frac{\partial \mathbf{r}_{jk}}{\partial r_a^\alpha}
     * \f}
     * and similarly for \f$\mathbf{r}_c\f$.
     * Finally,
     * \f{equation}{
     *     \frac{\partial \mathbf{r}_{ij}}{\partial r_a^\alpha} = \left(\delta_{aj} - \delta_{ai}\right)
     *     \begin{bmatrix} \delta_{\alpha x} \\ \delta_{\alpha y} \\ \delta_{\alpha z} \end{bmatrix}
     * \f}
     * where \f$\delta_{mn}=1\f$ when \f$m=n\f$ and 0 otherwise.
     * Therefore,
     * \f{multline}{
     *     \frac{\partial \mathbf{r}_b}{\partial r_a^\alpha} = \left(\delta_{aj} - \delta_{ai}\right)
     *     \begin{bmatrix} \delta_{\alpha x} \\ \delta_{\alpha y} \\ \delta_{\alpha z} \end{bmatrix} \times \mathbf{r}_{jk} -
     *     \left(\delta_{ak} - \delta_{aj}\right)
     *     \begin{bmatrix} \delta_{\alpha x} \\ \delta_{\alpha y} \\ \delta_{\alpha z} \end{bmatrix} \times \mathbf{r}_{ij}
     *     = \\
     *     \left(\delta_{aj} - \delta_{ai}\right)
     *     \begin{bmatrix} \delta_{\alpha y} r_{jk}^z - \delta_{\alpha z} r_{jk}^y \\
     *                     \delta_{\alpha z} r_{jk}^x - \delta_{\alpha x} r_{jk}^z \\
     *                     \delta_{\alpha x} r_{jk}^y - \delta_{\alpha y} r_{jk}^x
     *     \end{bmatrix} -
     *     \left(\delta_{ak} - \delta_{aj}\right)
     *     \begin{bmatrix} \delta_{\alpha y} r_{ij}^z - \delta_{\alpha z} r_{ij}^y \\
     *                     \delta_{\alpha z} r_{ij}^x - \delta_{\alpha x} r_{ij}^z \\
     *                     \delta_{\alpha x} r_{ij}^y - \delta_{\alpha y} r_{ij}^x
     *     \end{bmatrix}
     * \f}
     * and eventually,
     * \f{multline}{
     *     \frac{\partial (\mathbf{r}_b \cdot \mathbf{r}_c)}{\partial r_a^\alpha} =
     *     \frac{\partial \mathbf{r}_b}{\partial r_a^\alpha} \mathbf{r}_c + \mathbf{r}_b \frac{\partial \mathbf{r}_c}{\partial r_a^\alpha} = \\
     *     \left(\delta_{aj} - \delta_{ai}\right)
     *       \left[r_c^x (\delta_{\alpha y} r_{jk}^z - \delta_{\alpha z} r_{jk}^y) + 
     *             r_c^y (\delta_{\alpha z} r_{jk}^x - \delta_{\alpha x} r_{jk}^z) +
     *             r_c^z (\delta_{\alpha x} r_{jk}^y - \delta_{\alpha y} r_{jk}^x)
     *       \right] - \\
     *     \left(\delta_{ak} - \delta_{aj}\right)
     *       \left[r_c^x (\delta_{\alpha y} r_{ij}^z - \delta_{\alpha z} r_{ij}^y) +
     *             r_c^y (\delta_{\alpha z} r_{ij}^x - \delta_{\alpha x} r_{ij}^z) +
     *             r_c^z (\delta_{\alpha x} r_{ij}^y - \delta_{\alpha y} r_{ij}^x)
     *       \right] + \\
     *     \left(\delta_{ak} - \delta_{aj}\right)
     *       \left[r_b^x (\delta_{\alpha y} r_{kl}^z - \delta_{\alpha z} r_{kl}^y) +
     *             r_b^y (\delta_{\alpha z} r_{kl}^x - \delta_{\alpha x} r_{kl}^z) +
     *             r_b^z (\delta_{\alpha x} r_{kl}^y - \delta_{\alpha y} r_{kl}^x)
     *       \right] - \\
     *     \left(\delta_{al} - \delta_{ak}\right)
     *       \left[r_b^x (\delta_{\alpha y} r_{jk}^z - \delta_{\alpha z} r_{jk}^y) +
     *             r_b^y (\delta_{\alpha z} r_{jk}^x - \delta_{\alpha x} r_{jk}^z) +
     *             r_b^z (\delta_{\alpha x} r_{jk}^y - \delta_{\alpha y} r_{jk}^x)
     *       \right]
     * \f}
     *
     * The second derivative in \eqref{eq:B_diff} can be evaluated using
     * \f{equation}{
     *     \frac{\partial (\left|\left|\mathbf{r}_b\right|\right| \left|\left|\mathbf{r}_c\right|\right|)}{\partial r_a^\alpha} =
     *     \frac{\partial \left|\left|\mathbf{r}_b\right|\right|}{\partial r_a^\alpha} \left|\left|\mathbf{r}_c\right|\right| +
     *     \left|\left|\mathbf{r}_b\right|\right| \frac{\partial \left|\left|\mathbf{r}_c\right|\right|}{\partial r_a^\alpha}
     * \f}
     * The calculation of the square root connected with the norm can be reduced by realising
     * \f{equation}{
     *     \frac{\partial \left|\left|\mathbf{r}_b\right|\right|}{\partial r_a^\alpha} =
     *     \frac{1}{2 \left|\left|\mathbf{r}_b\right|\right|} \frac{\partial \left|\left|\mathbf{r}_b\right|\right|^2}{\partial r_a^\alpha}
     * \f}
     * Since
     * \f{equation}{
     *     \left|\left|\mathbf{r}_b\right|\right|^2 =
     *     (r_{ij}^y r_{jk}^z - r_{ij}^z r_{jk}^y)^2 + (r_{ij}^z r_{jk}^x - r_{ij}^x r_{jk}^z)^2 + (r_{ij}^x r_{jk}^y - r_{ij}^y r_{jk}^x)^2
     * \f}
     * the derivative is
     * \f{multline}{
     *    \frac{\partial \left|\left|\mathbf{r}_b\right|\right|}{\partial r_a^\alpha}
     *    = \frac{r_b^x}{\left|\left|\mathbf{r}_b\right|\right|}
     *        \left[(\delta_{aj} - \delta_{ai})\left(r_{jk}^z \delta_{\alpha y} - r_{jk}^y \delta_{\alpha z}\right)
     *                     +(\delta_{ak} - \delta_{aj})\left(r_{ij}^y \delta_{\alpha z} - r_{ij}^z \delta_{\alpha y}\right)\right] + \\
     *      \frac{r_b^y}{\left|\left|\mathbf{r}_b\right|\right|}
     *        \left[(\delta_{aj} - \delta_{ai})\left(r_{jk}^x \delta_{\alpha z} - r_{jk}^z \delta_{\alpha x}\right)
     *                     +(\delta_{ak} - \delta_{aj})\left(r_{ij}^z \delta_{\alpha x} - r_{ij}^x \delta_{\alpha z}\right)\right] + \\
     *      \frac{r_b^z}{\left|\left|\mathbf{r}_b\right|\right|}
     *        \left[(\delta_{aj} - \delta_{ai})\left(r_{jk}^y \delta_{\alpha x} - r_{jk}^x \delta_{\alpha y}\right)
     *                     +(\delta_{ak} - \delta_{aj})\left(r_{ij}^x \delta_{\alpha y} - r_{ij}^y \delta_{\alpha x}\right)\right]
     * \f}
     * and the final expression for the second derivative in \eqref{eq:B_diff} is
     * \f{multline}{
     *     \frac{\partial (\left|\left|\mathbf{r}_b\right|\right| \left|\left|\mathbf{r}_c\right|\right|)}{\partial r_a^\alpha} =
     *     \frac{\partial \left|\left|\mathbf{r}_b\right|\right|}{\partial r_a^\alpha} \left|\left|\mathbf{r}_c\right|\right| +
     *     \left|\left|\mathbf{r}_b\right|\right| \frac{\partial \left|\left|\mathbf{r}_c\right|\right|}{\partial r_a^\alpha} = \\
     *       \frac{\left|\left|\mathbf{r}_c\right|\right|}{\left|\left|\mathbf{r}_b\right|\right|}r_b^x
     *        \left[(\delta_{aj} - \delta_{ai})\left(r_{jk}^z \delta_{\alpha y} - r_{jk}^y \delta_{\alpha z}\right)
     *                     +(\delta_{ak} - \delta_{aj})\left(r_{ij}^y \delta_{\alpha z} - r_{ij}^z \delta_{\alpha y}\right)\right] + \\
     *      \frac{\left|\left|\mathbf{r}_c\right|\right|}{\left|\left|\mathbf{r}_b\right|\right|}r_b^y
     *        \left[(\delta_{aj} - \delta_{ai})\left(r_{jk}^x \delta_{\alpha z} - r_{jk}^z \delta_{\alpha x}\right)
     *                     +(\delta_{ak} - \delta_{aj})\left(r_{ij}^z \delta_{\alpha x} - r_{ij}^x \delta_{\alpha z}\right)\right] + \\
     *      \frac{\left|\left|\mathbf{r}_c\right|\right|}{\left|\left|\mathbf{r}_b\right|\right|}r_b^z
     *        \left[(\delta_{aj} - \delta_{ai})\left(r_{jk}^y \delta_{\alpha x} - r_{jk}^x \delta_{\alpha y}\right)
     *                     +(\delta_{ak} - \delta_{aj})\left(r_{ij}^x \delta_{\alpha y} - r_{ij}^y \delta_{\alpha x}\right)\right] + \\
     *      
     *      \frac{\left|\left|\mathbf{r}_b\right|\right|}{\left|\left|\mathbf{r}_c\right|\right|}r_c^x
     *        \left[(\delta_{ak} - \delta_{aj})\left(r_{kl}^z \delta_{\alpha y} - r_{kl}^y \delta_{\alpha z}\right)
     *                     +(\delta_{al} - \delta_{ak})\left(r_{jk}^y \delta_{\alpha z} - r_{jk}^z \delta_{\alpha y}\right)\right] + \\
     *      \frac{\left|\left|\mathbf{r}_b\right|\right|}{\left|\left|\mathbf{r}_c\right|\right|}r_c^y
     *        \left[(\delta_{ak} - \delta_{aj})\left(r_{kl}^x \delta_{\alpha z} - r_{kl}^z \delta_{\alpha x}\right)
     *                     +(\delta_{al} - \delta_{ak})\left(r_{jk}^z \delta_{\alpha x} - r_{jk}^x \delta_{\alpha z}\right)\right] + \\
     *      \frac{\left|\left|\mathbf{r}_b\right|\right|}{\left|\left|\mathbf{r}_c\right|\right|}r_c^z
     *        \left[(\delta_{ak} - \delta_{aj})\left(r_{kl}^y \delta_{\alpha x} - r_{kl}^x \delta_{\alpha y}\right)
     *                     +(\delta_{al} - \delta_{ak})\left(r_{jk}^x \delta_{\alpha y} - r_{jk}^y \delta_{\alpha x}\right)\right]
     * \f}
     * 
     * Expression used for calculations in the program were taken from DL POLY code (`dihedrals_forces.f90`) and checked against the presented formulae.
     * 
     * The final sign of forces depends on the quadrant of the dihedral angle. 
     * Therefore, the final formulae were taken from DL POLY code.
     * Sinus of the dihedral angle must be computed by
     * \f{equation}{
     *     \sin\varphi_{ijkl} = \frac{\mathbf{r}_{jk}\cdot \left(\mathbf{r}_b \times \mathbf{r}_c\right)}
     *     {\left|\left|\mathbf{r}_{jk}\right|\right|\left|\left|\mathbf{r}_{b}\right|\right|\left|\left|\mathbf{r}_{c}\right|\right|}
     * \f}
     * and the dihedral angle by function `atan2`.
     * 
     * The contribution to the viral by the cosine dihedral is zero, since the dihedral angle value remains the same when the system is shrunk or stretched.
     * Nothing is thus added to the parent Molecule::virial[2].
     *
     * @par See also
     * @ref CalculateEpot() const, GetCurrentAngle() const.
     *
     * @return Potential energy caused by this dihedral angle.
     */
    double CalculateForces() override;
    /**
     * @brief Calculate potential energy caused by this dihedral.
     * Similar to @ref CalculateForces(), but calculates only the potential energy and does not change anything neither in the @ref parent Molecule,
     * nor in the SimulatedSystem.
     * The potential energy is given by:
     * \f{equation}{
     *     E_{\mathrm{pot}} = A \left[ 1 + \cos\left(m\varphi_{ijkl} - \delta \right)\right]
     * \f}
     * where \f$A\f$ is the force constant ( @ref A), \f$\delta\f$ is the angle difference ( @ref delta), \f$m\f$ is the dihedral *multiplicity* ( @ref m) and
     * \f$\varphi_{ijkl}\f$ is the instanteneous dihedral angle defined by atoms @ref atomI, @ref atomJ, @ref atomK and @ref atomL of the @ref parent Molecule.
     *
     * @par See also
     * @ref CalculateForces().
     *
     * @return Potential energy caused by this cosine dihedral.
     */
    double CalculateEpot() const override; // calculate potential energy of the CosineDihedral
    /**
     * @brief Get the involved Atom index according to the index value.
     * Returns the label of an Atom valid inside a @ref Molecule.
     * The value of @p index is a label valid inside this CosineDihedral (1 for @ref atomI, 2 for @ref atomJ, 3 for @ref atomK and 4 for @ref atomL).
     *
     * @param index Label of an Atom inside this dihedral class.
     * @return Label of the same Atom inside the @ref parent Molecule.
     */
    int GetAtom(int index) const override; // get one of atoms I, J, K, L
    /**
     * @brief Clone itself and get adopted by another @ref parent Molecule.
     * Returns a copy of itself and migrate to the @p newparent Molecule.
     * Useful when copying Molecule object.
     * Every CosineDihedral is then cloned by this method and assigned to the new Molecule object.
     *
     * @par See also
     * @ref CosineDihedral(const CosineDihedral&), @ref SetParent(Molecule*)
     *
     * @param newparent Parent Molecule of the newly created CosineDihedral object.
     * @return (Pointer to) the newly created CosineDihedral object.
     */
    CosineDihedral *copy(Molecule *newparent) const override; // clone itself
    /**
     * @brief Get the dihedral parameters: force constant @ref A, equilibrium angle @ref delta and multiplicity @ref m
     * 
     * @param fieldname Human-readable type of this field (`cos` in this case).
     * @param params Force constant, angle and multiplicity to be printed in info.
     * @return 3
     */
    int GetParams(std::string &fieldname, std::vector<double> &params) const override;
};

#endif