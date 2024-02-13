#ifndef LJSYSTEMHEADER
#define LJSYSTEMHEADER

#include "AbstractInterMolField.hpp"

#include "Matrix.hpp"
class SimulatedSystem;

#ifdef PARALLEL
extern int thread_count;
#endif

/**
 * @class LJsystem
 * @ingroup InterMolField
 * @brief Description of Lennard-Jones (12–6) interactions between Atoms (from different Molecules).
 *
 * This class implements the interface for intermolecular interactions @ref AbstractInterMolField ( @ref type = 1).
 * Disperse LJ interactions are determined by the Atom unique ID ( @ref Atom::LJid) which is uniquely connected with @ref Atom::name.
 * LJ potential has two terms. The strong Pauli repulsion is approximated by the term proportional to \f$1/r^{12}\f$, the London cohesive potential is proportional to \f$1/r^6\f$, where \f$r\f$ is the interatomic distance.
 * The basic form of the LJ potential is thus:
 * \f{equation}{
 *     E_{\mathrm{LJ}} = \sum_{i,j} 4 \varepsilon_{ij} \left[\left(\frac{\sigma_{ij}}{r_{ij}}\right)^{12} - \left(\frac{\sigma_{ij}}{r_{ij}}\right)^6\right]
 * \f}
 * where \f$r_{ij}\f$ is the 2-norm (see @ref Vector::CalculateNorm(int) const) of the vector \f$\mathbf{r}_{ij} = \mathbf{r}_j - \mathbf{r}_i\f$, where \f$\mathbf{r}_i\f$ is the position of Atom \f$i\f$.
 * The sum runs over all pairs of Atoms in the SimulatedSystem, except for Atoms belonging to the same Molecule.
 * For intramolecular interactions, the corresponding term must be specified as an intramolecular field (see @ref AbstractIntraMolField).
 * The potential has two parameters which are defined for each pair type in the `.field` file.
 * (Therefore, the parameters are the same for all pairs with the same @ref Atom::LJid)
 * The first is the energy parameter \f$\varepsilon_{ij}\f$, which determines the depth of the attractive well.
 * In `.field` file it is given in corresponding energy units (directive `units`), internally it is converted to program units.
 * The second parameter \f$\sigma_{ij}\f$ determines the range of interaction and is given in [AA].
 *
 * The above form is valid in free boundary conditions.
 * If the SimulatedSystem is subjected to periodic boundary conditions, the nearest (minimum) image convention and some form of cutoff must be incorporated.
 * The former is done as usual (see e.g. @cite AllenCompSimLiq, p 29):
 * \f{equation}{
 *     (\mathbf{r}_{ij})_a = \begin{cases} (\mathbf{r}_{ij})_a + L_a/2 & \text{if } (\mathbf{r}_{ij})_a < L_a/2\\
 *                                         (\mathbf{r}_{ij})_a - L_a/2 & \text{if } (\mathbf{r}_{ij})_a > L_a/2\\
 *                                         (\mathbf{r}_{ij})_a       & \text{otherwise}
 *                           \end{cases}
 * \f}
 * where \f$a\f$ runs over \f$x, y\f$ and \f$z\f$ (or 0, 1 and 2, respectively) and \f$(\mathbf{r}_{ij})_a\f$ denotes the \f$a-\f$th component of the Vector \f$\mathbf{r}_{ij}\f$, \f$L_a\f$ is the simulation box size (see @ref SimulatedSystem::boxSize) in the \f$a\f$ direction.
 *
 * The latter can be handled in various ways.
 * The approach used in MACSIMUS (see @cite MACSIMUS_manual p 157) is copied.
 * Beyond @ref cutoff distance, the potential is zero.
 * To reach smoothness in potential and continuity in forces, the intermediate potential is introduced whose parameters (intermediate cutoff distance \f$C\f$ and energy parameter \f$A\f$) are uniquely determined by boundary conditions (satisfying continuity of forces).
 * The overall potential is then:
 * \f{equation}{
 *     E_{\mathrm{LJ}} = \sum_{i,j} \begin{cases} 4 \varepsilon_{ij} \left[\left(\frac{\sigma_{ij}}{r_{ij}}\right)^{12} - \left(\frac{\sigma_{ij}}{r_{ij}}\right)^6\right] & \text{if } r_{ij} \leq C_{ij} \\
 *                                     A_{ij} (r_{ij}^2 - \mathtt{cutoff}^2)^2 & \text{if } C_{ij} < r_{ij} < \mathtt{cutoff} \\
 *                                     0 & \text{otherwise}
 *                                   \end{cases}
 * \f}
 * The intermediate cutoff and parameter \f$A\f$ are calculated for each pair of \f$\varepsilon_{ij}\f$ and \f$\sigma_{ij}\f$ separately.
 * See @ref CalculateC1A() for details.
 *
 * The main purpose of this class, i.e. the calculation of forces and potential energy, is done by implementing @ref CalculateEpot() const which yields only the potential energy of the LJ interactions and
 * @ref CalculateForces() which also modifies the forces for each @ref Atom in the SimulatedSystem.
 */
class LJsystem : public AbstractInterMolField
{
private:
  /**
   * @brief Size of the LJsystem.
   *
   * Number of unique Atoms (with unique name (LJid)).
   */
  int size; // size of the system
  /**
   * @brief Matrix of LJ \f$\sigma_{ij}\f$ for each pair.
   */
  Matrix *sigmaMatrix; // matrix of LJ sigma values (of each pair)
  /**
   * @brief Matrix of LJ \f$\sigma_{ij}^6\f$ for each pair.
   * The sixth power is stored to speed up the evaluation.
   */
  Matrix *sigma6Matrix; // matrix of LJ sigma^6 (for each pair for easier forces calculation)
  /**
   * @brief Matrix of LJ \f$\varepsilon_{ij}\f$ for each pair.
   */
  Matrix *epsilonMatrix; // matrix of LJ epsilon values (of each pair)
  /**
   * @brief Matrix of intermediate cutoffs \f$C_{ij}\f$ for each pair.
   */
  Matrix *C1Matrix; // matrix of MACSIMUS style C1 ("first cutoff")
  /**
   * @brief Matrix of intermediate energy parameters \f$A_{ij}\f$ for each pair.
   */
  Matrix *AMatrix; // matrix of MACSIMUS style A (cutoff energy...)
  /**
   * @brief Default constructor.
   *
   * Made private in order not to be used because an empty system should not be instantiated.
   */
  LJsystem(); // making default constructor private disables its use
  /**
   * @brief If *true*, the system does not contribute to potential energy (all \f$\varepsilon_{ij}\f$ are zero).
   *
   * Setting all \f$\varepsilon_{ij}\f$ to zero can be advantageous if you want to specify the Atom diameter for `.mol` file, but do not want any LJ interactions.
   * The value of @ref empty is then *true* and the force or potential energy calculation is just skipped.
   */
  bool empty; // empty LJ system (no LJ interactions, zero energy)
public:
  /**
   * @brief Primary constructor, sets the parent SimulatedSystem and @ref size.
   *
   * @param parent Parent SimulatedSystem (see @ref AbstractInterMolField::parent)
   * @param siz LJsystem @ref size (number of atoms with unique @ref Atom::LJid)
   */
  LJsystem(SimulatedSystem *parent, int siz); // constructor (must have map form atom names to LJid prepared)
  /**
   * @brief LJsystem copy constructor
   *
   * @param system LJsystem to be copied.
   */
  LJsystem(const LJsystem &system); // copy constructor
  /**
   * @brief LJsystem destructor
   *
   * Destroy all matrices with parameters (see @ref epsilonMatrix, @ref sigmaMatrix, @ref sigma6Matrix, @ref C1Matrix, @ref AMatrix)
   */
  ~LJsystem(); // destructor
  /**
   * @brief Get the \f$\sigma_{ij}\f$ for the given pair (\f$i, j\f$)
   *
   * Provides access to private matrix of \f$\sigma_{ij}\f$, @ref sigmaMatrix.
   *
   * @par See also
   * @ref GetSigma6(int, int) const, @ref GetEpsilon(int, int) const.
   *
   * @param LJid1 @ref Atom::LJid of the first Atom (\f$i\f$).
   * @param LJid2 @ref Atom::LJid of the second Atom (\f$j\f$).
   * @return LJ parameter \f$\sigma_{ij}\f$.
   */
  inline double GetSigma(int LJid1, int LJid2) const
  {
    return sigmaMatrix->operator()(LJid1, LJid2);
  }; // return LJ sigma of given pair
  /**
   * @brief Get the \f$\sigma_{ij}^6\f$ for the given pair (\f$i, j\f$)
   *
   * Provides access to private matrix of \f$\sigma_{ij}^6\f$, @ref sigma6Matrix.
   * The sixth power is stored to speedup the calculation.
   *
   * @par See also
   * @ref GetSigma(int, int) const, @ref GetEpsilon(int, int) const.
   *
   * @param LJid1 @ref Atom::LJid of the first Atom (\f$i\f$).
   * @param LJid2 @ref Atom::LJid of the second Atom (\f$j\f$).
   * @return LJ parameter \f$\sigma_{ij}^6\f$.
   */
  inline double GetSigma6(int LJid1, int LJid2) const
  {
    return sigma6Matrix->operator()(LJid1, LJid2);
  }; // return LJ sigma^6 of given pair
  /**
   * @brief Get the \f$\varepsilon_{ij}\f$ for the given pair (\f$i, j\f$)
   *
   * Provides access to private matrix of \f$\varepsilon_{ij}\f$, @ref epsilonMatrix.
   *
   * @par See also
   * @ref GetSigma(int, int) const, @ref GetSigma6(int, int) const.
   *
   * @param LJid1 @ref Atom::LJid of the first Atom (\f$i\f$).
   * @param LJid2 @ref Atom::LJid of the second Atom (\f$j\f$).
   * @return LJ parameter \f$\varepsilon_{ij}\f$.
   */
  inline double GetEpsilon(int LJid1, int LJid2) const
  {
    return epsilonMatrix->operator()(LJid1, LJid2);
  }; // return LJ epsilon of given pair
  /**
   * @brief Get the cutoff parameter \f$C_{ij}\f$ for the given pair (\f$i, j\f$)
   *
   * Provides access to private matrix of \f$C_{ij}\f$, @ref C1Matrix.
   *
   * @par See also
   * @ref GetA(int, int) const, @ref CalculateC1A().
   *
   * @param LJid1 @ref Atom::LJid of the first Atom (\f$i\f$).
   * @param LJid2 @ref Atom::LJid of the second Atom (\f$j\f$).
   * @return LJ cutoff parameter \f$C_{ij}\f$.
   */
  inline double GetC1(int LJid1, int LJid2) const
  {
    return C1Matrix->operator()(LJid1, LJid2);
  }; // return C1 (MACSIMUS first cutoff) for the given pair
  /**
   * @brief Get the cutoff parameter \f$A_{ij}\f$ for the given pair (\f$i, j\f$)
   *
   * Provides access to private matrix of \f$A_{ij}\f$, @ref AMatrix.
   *
   * @par See also
   * @ref GetC1(int, int) const, @ref CalculateC1A().
   *
   * @param LJid1 @ref Atom::LJid of the first Atom (\f$i\f$).
   * @param LJid2 @ref Atom::LJid of the second Atom (\f$j\f$).
   * @return LJ cutoff parameter \f$A_{ij}\f$.
   */
  inline double GetA(int LJid1, int LJid2) const
  {
    return AMatrix->operator()(LJid1, LJid2);
  }; // retrun A parameter (MACSIMUS style cutoff energy)
  /**
   * @brief Set the LJ parameter \f$\sigma_{ij}\f$ for the given pair (\f$i, j\f$)
   *
   * Provide method to set the member of private @ref sigmaMatrix.
   * Modify also @ref sigma6Matrix.
   *
   * @par See also
   * @ref SetEpsilon(int, int, double).
   *
   * @param LJid1 Atom::LJid of the first Atom (\f$i\f$).
   * @param LJid2 Atom::LJid of the second Atom (\f$i\f$).
   * @param sigma LJ parameter \f$\sigma\f$ for the given pair (\f$\sigma_{ij}\f$)
   * @return 0 if success
   * @return 29 if the desired value of parameter @p sigma is out of expected range.
   */
  int SetSigma(int LJid1, int LJid2, double sigma); // set sigma of given pair (given ID)
  /**
   * @brief Set the LJ parameter \f$\varepsilon_{ij}\f$ for the given pair (\f$i, j\f$)
   *
   * Provide method to set the member of private @ref epsilonMatrix.
   * Set @ref empty to *false* if the value of @p epsilon is greater then zero, because it means that the forces (and energy) are not zero.
   *
   * @par See also
   * @ref SetSigma(int, int, double).
   *
   * @param LJid1 Atom::LJid of the first Atom (\f$i\f$).
   * @param LJid2 Atom::LJid of the second Atom (\f$i\f$).
   * @param epsilon LJ parameter \f$\varepsilon\f$ for the given pair (\f$\varepsilon_{ij}\f$)
   * @return 0 if success
   * @return 28 if the desired value of parameter @p epsilon is out of expected range.
   */
  int SetEpsilon(int LJid1, int LJid2, double epsilon); // set epsilon of given pair (given ID)
  /**
   * @brief Set the LJ parameters \f$\varepsilon_{ij}\f$ and \f$\sigma_{ij}\f$ for the given pair (\f$i, j\f$)
   *
   * Implements the pure virtual method from the parent class @ref AbstractInterMolField.
   * LJ parameters are given as a string in the form used in the `.field` file.
   * The idea is that the concrete class must know how to handle the description in the `.field` file.
   * Both parameters \f$\varepsilon_{ij}\f$ and \f$\sigma_{ij}\f$ are set and appropriate modifications to @ref sigma6Matrix and @ref empty are made.
   *
   * @par See also
   * @ref SetSigma(int, int, double), @ref SetEpsilon(int, int, double).
   *
   * @param LJid1 Atom::LJid of the first Atom (\f$i\f$).
   * @param LJid2 Atom::LJid of the second Atom (\f$i\f$).
   * @param params Line of pair potential parameters from the `.field` file (`(atom1 atom2) epsilon sigma`)
   * @return 0 if success
   * @return 28 if the desired value of parameter epsilon \f$\varepsilon_{ij}\f$ is out of expected range.
   * @return 29 if the desired value of parameter sigma \f$\sigma_{ij}\f$ is out of expected range.
   */
  int SetPairParams(int LJid1, int LJid2, char *params) override; // set parameters (sigma && epsilon for the given pair)
  /**
   * @brief Get the @ref size of this @ref LJsystem (number of Atoms with unique name/LJid).
   *
   * Provides read-only access to private member @ref size.
   *
   * @return The size of the system (see @ref size).
   */
  int GetSystemSize() const; // return number of different LJ atom types
  /**
   * @brief Calculate forces caused by Lennard-Jones interactions (intermolecular).
   * Calculate forces caused by this field (LJsystem) and add them to @ref Atom::force.
   * Save the potential energy and force virial to the @ref parent SimulatedSystem @ref SimulatedSystem::Epotintermol[1] and @ref SimulatedSystem::virintermol[1].
   *
   *  The potential energy caused by Lennard-Jones disperse interactions:
   *  \f{equation}{
   *      E_{\mathrm{LJ}} = \sum_{i,j} \begin{cases} 4 \varepsilon_{ij} \left[\left(\frac{\sigma_{ij}}{r_{ij}}\right)^{12} - \left(\frac{\sigma_{ij}}{r_{ij}}\right)^6\right] & \text{if } r_{ij} \leq C_{ij} \\
   *                                     A_{ij} (r_{ij}^2 - \mathtt{cutoff}^2)^2 & \text{if } C_{ij} < r_{ij} < \mathtt{cutoff} \\
   *                                     0 & \text{otherwise}
   *                                   \end{cases}
   * \f}
   * where \f$r_{ij}\f$ is the 2-norm (see @ref Vector::CalculateNorm(int) const) of the Vector \f$\mathbf{r}_{ij} = \mathbf{r}_j - \mathbf{r}_i\f$, where \f$\mathbf{r}_i\f$ is the position of Atom \f$i\f$.
   * The sum runs over all pairs of Atoms in the SimulatedSystem, except for Atoms belonging to the same Molecule.
   * For intramolecular interactions, the corresponding term must be specified as an intramolecular field (see @ref AbstractIntraMolField).
   * The potential energy is added to the parent SimulatedSystem::Epotintermol[1].
   * In case of free periodic boundary conditions, the energy is the same as if \f$C_{ij}\rightarrow\infty\f$ and \f$\mathtt{cutoff}\rightarrow\infty\f$.
   *
   * Forces (3D vector) are given by the derivative of the potential energy:
   * \f{equation}{
   *     \mathbf{f}_j = -\mathbf{f}_i =
   *          \begin{cases}
   *               4 \varepsilon_{ij} \left[\left(\frac{12\sigma_{ij}^{12}}{r_{ij}^{13}}\right) - \left(\frac{6\sigma_{ij}^6}{r_{ij}^7}\right)\right]\frac{\mathbf{r}_{ij}}{r_{ij}} & \text{if } r_{ij} \leq C_{ij} \\
   *               4 A_{ij} r_{ij} (-r_{ij}^2 + \mathtt{cutoff}^2)\frac{\mathbf{r}_{ij}}{r_{ij}} & \text{if } C_{ij} < r_{ij} < \mathtt{cutoff} \\
   *               0 & \text{otherwise}
   *            \end{cases}
   * \f}
   * Note that only even powers of \f$r_{ij}\f$ are needed throughout the computation.
   *
   * The contribution to the viral by each interaction pair is:
   * \f{equation}{
   *     w_{\mathrm{LJ}} = \mathbf{f}_i \cdot \mathbf{r}_{ij} =
   *         \begin{cases}
   *               -4 \varepsilon_{ij} \left[\left(\frac{12\sigma_{ij}^{12}}{r_{ij}^{13}}\right) - \left(\frac{6\sigma_{ij}^6}{r_{ij}^7}\right)\right]r_{ij} & \text{if } r_{ij} \leq C_{ij} \\
   *               -4 A_{ij} r_{ij}^2 (-r_{ij}^2 + \mathtt{cutoff}^2) & \text{if } C_{ij} < r_{ij} < \mathtt{cutoff} \\
   *               0 & \text{otherwise}
   *         \end{cases}
   * \f}
   * This value is added to the parent @ref SimulatedSystem::virintermol[1].
   *
   * Returns the potential energy.
   *
   * @par See also
   * @ref CalculateEpot() const.
   *
   * @return Potential energy caused by LJ interactions in the @ref parent @ref SimulatedSystem (internally in program units).
   */
  double CalculateForces() override; // calculate forces in given simSystem
  /**
   * @brief Calculate potential energy caused by Lennard-Jones interactions (intermolecular).
   * Similar to @ref CalculateForces(), but calculates only the potential energy and does not change anything neither in the @ref parent SimulatedSystem, nor in the involved Atoms.
   *  The potential energy cause by Lennard-Jones disperse interactions:
   *  \f{equation}{
   *      E_{\mathrm{LJ}} = \sum_{i,j} \begin{cases} 4 \varepsilon_{ij} \left[\left(\frac{\sigma_{ij}}{r_{ij}}\right)^{12} - \left(\frac{\sigma_{ij}}{r_{ij}}\right)^6\right] & \text{if } r_{ij} \leq C_{ij} \\
   *                                     A_{ij} (r_{ij}^2 - \mathtt{cutoff}^2)^2 & \text{if } C_{ij} < r_{ij} < \mathtt{cutoff} \\
   *                                     0 & \text{otherwise}
   *                                   \end{cases}
   * \f}
   * where \f$r_{ij}\f$ is the 2-norm (see @ref Vector::CalculateNorm(int) const) of the Vector \f$\mathbf{r}_{ij} = \mathbf{r}_j - \mathbf{r}_i\f$, where \f$\mathbf{r}_i\f$ is the position of Atom \f$i\f$.
   * The sum runs over all pairs of Atoms in the SimulatedSystem, except for Atoms belonging to the same Molecule.
   * For intramolecular interactions, the corresponding term must be specified as an intramolecular field (see @ref AbstractIntraMolField).
   * In case of free periodic boundary conditions, the energy is the same as if \f$C_{ij}\rightarrow\infty\f$ and \f$\mathtt{cutoff}\rightarrow\infty\f$.
   *
   * @par See also
   * @ref CalculateForces().
   *
   * @return Potential energy caused by LJ interactions in the @ref parent @ref SimulatedSystem (internally in program units).
   */
  double CalculateEpot() const override; // calculate potential energy in given simSystem
  /**
   * @brief Calculate intermediate cutoff parameters \f$C_{ij}\f$ and \f$A_{ij}\f$ for each pair.
   * Intermediate cutoff parameters are calculated in such a way that both the potential and force are continuous.
   * The system of two equations must be solved
   * \f{align}{
   *     4 \varepsilon_{ij} \left[\left(\frac{\sigma_{ij}}{C_{ij}}\right)^{12} - \left(\frac{\sigma_{ij}}{C_{ij}}\right)^6\right] &=
   *               A_{ij} (C_{ij}^2 - \mathtt{cutoff}^2)^2 \\
   *     4 \varepsilon_{ij} \left[\left(\frac{12\sigma_{ij}^{12}}{C_{ij}^{13}}\right) - \left(\frac{6\sigma_{ij}^6}{C_{ij}^7}\right)\right] &=
   *               4 A_{ij} C_{ij} (C_{ij}^2 - \mathtt{cutoff}^2) \\
   * \f}
   * The continuity at `cutoff` is guaranteed by the form of the potential.
   *
   * The system has an analytical solution (yet quite complicated, see Maple `LJforces_cutoff2.mw` in the documentation)
   * if \f$\texttt{cutoff} > 2^{5/6} \sigma_{ij}\f$.
   *
   * The calculated values are saved in matrices @ref C1Matrix and @ref AMatrix.
   * 
   * @return 0 upon success, else number of pairs with too short cutoff
   */
  int CalculateC1A(); // calculate C1 and A for each pair
  /**
   * @brief Set the @ref cutoff
   * Set the distance beyond which the interaction is zero.
   * Applies only in periodic boundary conditions (then @p boundaryC > 0).
   * Modifies also @ref cutoff2 and then uses @ref CalculateC1A() to calculate the parameters of intermediate cutoff modification in MACSIMUS style.
   *
   * @param cut Cutoff distance for disperse interactions (in [AA])
   * @param elcut Cutoff distance for electrostatics (in [AA]) – irrelevant here
   * @param boundaryC Boundary conditions (see @ref SimulatedSystem::boundaryCond)
   * @return 0 if success, else number of pairs with too short cutoff
   */
  int SetCutoff(double cut, double elcut, int boundaryC) override; // set lj interactions spherical cutoff
  /**
   * @brief Clone itself and get adopted by another @ref parent SimulatedSystem.
   * Returns a copy of itself and migrate to the @p newparent SimulatedSystem.
   * Useful when copying SimulatedSystem object.
   * Every intermolecular term is then cloned by this method and assigned to the new SimulatedSystem object.
   *
   * @par See also
   * @ref LJsystem(const LJsystem&), @ref SetParent(SimulatedSystem*), @ref AbstractInterMolField::copy(SimulatedSystem *) const, @ref AbstractIntraMolField::copy(Molecule *) const.
   *
   * @param newparent Parent SimulatedSystem of the newly created LJsystem object.
   * @return (Pointer to) the newly created LJsystem object.
   */
  LJsystem *copy(SimulatedSystem *newparent) const override; // cloning itself
  /**
   * @brief Get the diameter for atoms with given @ref Atom::LJid.
   * Lennard-Jones parameter \f$\sigma_{ii}\f$ can be used to estimate the vdw radius of Atoms of type \f$i\f$.
   * Thus, having defined the self pair interaction, the diameter can be estimated.
   * This estimate is used only for visualisation (for `.mol` file, see @ref FieldFile::PrintMol(char*, SimulatedSystem*))
   *
   * The diameter is calculated by (adopted from MACSIMUS (see @cite MACSIMUS_manual) where it is not documented...)
   * \f{equation}{
   *     70\sigma_{ii} 2^{-5/6}
   * \f}
   *
   * @param LJid ID (@ref Atom::LJid) of the Atom (unique number corresponding to @ref Atom::name)
   * @return Atom diameter estimate rounded to integer number
   */
  int GetDiameter(int LJid) const override; // get atom diameter
  /**
   * @brief Calculate cross-terms for LJ interactions between pairs of different atoms
   * Unspecified mixed pairs parameters are calculated according to specified @p mixingrule.
   * 1. **Lorentz–Berthelot**
   *     - energy parameter \f$\varepsilon_{12}\f$ is calculated as a geometric mean
   *         \f{equation}{
   *             \varepsilon_{ij} = \sqrt{\varepsilon_i \varepsilon_j}
   *         \f}
   *     - aritmetic mean is used for atomic diameters
   *         \f{equation}{
   *             \sigma_{ij} = \frac{i}{j}(\sigma_i + \sigma_j)
   *         \f}
   *
   * @param mixingrule Specification of mixing rule to be used
   * @return Number of pairs whose interaction parameters were determined automatically
   */
  int CalculateMixTerms(int mixingrule) override;
  /**
   * @brief Calculate corrections arising from cutting the potential.
   * Generally we want to avoid consequences of cutting the potential after certain distance.
   * We want to obtain results for the system with the original Lennard-Jones potential without cutoff.
   * However, we measure quantities in a system with different potential (if we use periodic boundary conditions).
   * Therefore, to reach the results valid for the original system, we should try to estimate the discrepancy in the potential energy (and pressure) caused by cutting the potential.
   * Generally, cutoff corrections are given by
   * \f{equation}{
   *      E_\mathrm{LJ,corr} = 2\pi \sum_{i,j} \frac {N_i N_j}{V} \int_{0}^{\infty}\! \left(E_{\mathrm{LJ,id}} - E_{\mathrm{LJ,sim}}\right) g(r) r^2\,\mathrm{d}r
   * \f}
   * where \f$N_i\f$ is the number of Atoms with @ref Atom::LJid = \f$i\f$, \f$E_{\mathrm{LJ,id}}\f$ is the original Lennard-Jones potential, \f$E_{\mathrm{LJ,sim}}\f$ is the truncated potential used in the simulation and \f$g(r)\f$ is the radial distribution function.
   * Essentially, we must integrate the difference between the real and desired potential for all pairs in the system.
   * Having no reasonable estimate for the radial distribution function, \f$g(r)=1\f$ is usually assumed (after the cutoff distance this estimate should be decent).
   * For our version of LJ cutoff, it thus holds
   * \f{multline}{
   *     E_\mathrm{LJ,corr} = 2\pi \sum_{i,j} \frac {N_i N_j}{V}\left(
   * \int_{C_{ij}}^{\mathtt{cutoff}}\! \left\{4 \varepsilon_{ij} \left[\left(\frac{\sigma_{ij}}{r}\right)^{12} - \left(\frac{\sigma_{ij}}{r}\right)^6\right] - A_{ij} \left(r^2 - \mathtt{cutoff}^2\right)^2\right\} r^2\,\mathrm{d}r
   * +  \int_{\mathtt{cutoff}}^\infty \! 4 \varepsilon_{ij} \left[\left(\frac{\sigma_{ij}}{r}\right)^{12} - \left(\frac{\sigma_{ij}}{r}\right)^6\right]r^2\,\mathrm{d}r\right) = \\
   * =  2\pi \sum_{i,j} \frac {N_i N_j}{V}\left\{\left[ 4 \varepsilon_{ij} \sigma_{ij}^6 \left(\frac {\sigma_{ij}^6 r^{-9}}{-9} - \frac{r^{-3}}{-3}\right) -A\left(\frac{r^7}{7} - \frac{2r^5C_{ij}^2}{5}+\frac{r^3C_{ij}^4}{3}\right)\right]_{r=C_{ij}}^{r=\mathtt{cutoff}} + \left[4\varepsilon_{ij}\sigma_{ij}^6 \left(\frac{\sigma_{ij}^6 r^{-9}}{-9} - \frac{r^{-3}}{-3}\right) \right]_{r=\mathtt{cutoff}}^{r=\infty} \right\}= \\
   * =  2\pi \sum_{i,j} \frac {N_i N_j}{V}\left\{ 4\varepsilon_{ij}\sigma_{ij}^6 \left(\frac{\sigma_{ij}^6}{9C_{ij}^9}-\frac{1}{3C_{ij}^3}\right) + A\left[ \mathtt{cutoff}^7\left(\frac{2}{5} - \frac{1}{7} - \frac{1}{3}\right) + \frac{C_{ij}^7}{7} - \frac{2C_{ij}^5\mathtt{cutoff}^2}{5}+\frac{C_{ij}^3 \mathtt{cutoff}^4}{3}\right] \right\} =\\
   * = 2\pi \sum_{i,j} \frac {N_i N_j}{V}\left[ \frac{4\varepsilon_{ij}\sigma_{ij}^6}{3C_{ij}^3} \left(\frac{\sigma_{ij}^6}{3C_{ij}^6}-1\right) + A\left( \frac{C_{ij}^7}{7} -\frac{8\mathtt{cutoff}^7}{105}  - \frac{2C_{ij}^5\mathtt{cutoff}^2}{5}+\frac{C_{ij}^3 \mathtt{cutoff}^4}{3}\right)\right]
   * \f}
   * where \f$\left[f(r)\right]_{r=a}^{r=b}\f$ is a shorthand for \f$f(b)-f(a)\f$ and \f$\frac{N_i N_j}{V}\f$ is the *density of pair* \f$ij\f$.
   * For \f$i=j\f$ we take \f$N_j = N_i - 1\f$, where \f$N_i\f$ is the real number of atoms \f$i\f$ in the system (avoiding self-interactions) but IMHO the choice is not obvious.
   *
   * Correction for pressure is \f$-\frac{\partial E_\mathrm{LJ,corr}}{\partial V} = \frac{E_\mathrm{LJ,corr}}{V}\f$.
   *
   * @return Correction energy times volume (\f$E_\mathrm{LJ,corr}V\f$).
   */
  double CalculateCutoffCorrection() override; // calculate cutoff correction energy (returns energy times volume)
  /**
   * @brief Initialize fake *bonds* used to handle Lennard-Jones interactions between atoms of one Molecule.
   *
   * Lennard-Jones interactions are basically intermolecular, but they occur also in longer molecules.
   * They are realized by bonds of class @ref LJBond.
   * Atoms on distances 1–4 interact by scaled LJ potential,
   * atoms which are furhter from each other interact as if they belonged to different molecules.
   *
   * @param LJ14factor Scaling factor for 1–4 disperse interactions.
   * @param el14factor Scaling factor for 1–4 electrostatic interactions.
   * @return Number of newly initialized bonds.
   */
  int InitializeIntramolBonds(double LJ14factor = 1.0, double el14factor = 1.0) const override; // initialize fake *bonds* used to handle this intermolecular field between atoms in one molecule
  /**
   * @brief Print informations about the LJsystem to the @p stream
   *
   * @param stream Output stream.
   * @param u_eng Energy unit value (see `units.hpp`).
   * @param engunit Energy unit name in brackets (as in `.cpa` file).
   */
  void PrintInfo(std::ofstream &stream, double u_eng, std::string engunit) const override;
  /**
   * @brief System of equations used to determine parameters \f$A_{ij}\f$ and \f$C_{ij}\f$.
   * The use of the class as a friend is needed to enable solving parametrized equations (the same equation with varying \f$\varepsilon_{ij}\f$, \f$\sigma_{ij}\f$ and `cutoff`).
   */
  friend class EquationSystem;
};

#endif