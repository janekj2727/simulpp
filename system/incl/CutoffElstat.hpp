#ifndef CUTOFFELSTATHEADER
#define CUTOFFELSTATHEADER

#include "AbstractInterMolField.hpp"

class Matrix;
class SimulatedSystem;

#ifdef PARALLEL
extern int thread_count;
#endif

/**
 * @class CutoffElstat
 * @ingroup InterMolField
 * @brief The simplest possible description of electrostatic interactions between Atoms (from different Molecules) truncated after certain distance.
 *
 * This class implements the interface for intermolecular interactions @ref AbstractInterMolField ( @ref type = 0).
 * Electrostatic interactions are given by the charge of Atoms (see @ref Atom::charge).
 * The original Coulomb interaction energy is given by
 * \f{equation}{
 *     E_{\mathrm{Coul}} = \sum_{i,j} \frac{q_i q_j}{4\pi \epsilon_0} \frac{1}{r_{ij}}
 * \f}
 * where \f$r_{ij}\f$ is the 2-norm (see @ref Vector::CalculateNorm(int) const) of the vector \f$\mathbf{r}_{ij} = \mathbf{r}_j - \mathbf{r}_i\f$, where \f$\mathbf{r}_i\f$ is the position of Atom \f$i\f$.
 * Charge of Atom \f$i\f$ (see @ref Atom::charge) is denoted as \f$q_i\f$.
 * The sum runs over all pairs of Atoms in the SimulatedSystem, except for Atoms belonging to the same Molecule.
 * For intramolecular interactions, the corresponding term must be specified as an intramolecular field (see @ref AbstractIntraMolField).
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
 * The approach used in MACSIMUS (see @cite MACSIMUS_manual p 206 or @cite Kolafa2008) is copied.
 * In this simplest implementation, electrostatic interactions are truncated after certain @ref cutoff distance.
 * This approach is justified only if the system is neutral and cutoff is much longer than the Debye screening length (see @ref CalculateCutoffCorrection()).
 * To ensure continuous forces (and their derivatives), the original Coulomb interaction is shifted at short separations, neglected beyond @ref cutoff and smoothly interpolated in between.
 * The basic form of the cut-and-shifted electrostatic potential is thus:
 * \f{equation}{
 *     E_{\mathrm{cutel}} = \sum_{i,j} \begin{cases}
 *          \frac{q_i q_j}{4\pi \epsilon_0} \left(\frac{1}{r_{ij}} - \mathtt{shift}\right) & r_{ij} < \alpha \mathtt{cutoff} \\
 *          \frac{q_i q_j}{4\pi \epsilon_0} \left[\left(r_{ij}-\mathtt{cutoff}\right)^3(A+Br_{ij})\right] & \alpha \mathtt{cutoff} \leq r_{ij} < \mathtt{cutoff} \\
 *          0 & r_{ij} \leq \mathtt{cutoff} \\
 *          \end{cases}
 *     \label{eq:elstat_cut_shift}
 * \f}
 * The value of \f$\alpha\f$ must be given by the user (values around 0.7 or 0.666 are recommended by MACSIMUS maual @cite MACSIMUS_manual).
 * Values of parameters \f$A\f$, \f$B\f$ and \f$\mathtt{shift}\f$ are then obtained from the conditions on continuity of the potential and its first two derivatives with respect to \f$r_{ij}\f$.
 * These conditions yield to equations that can be solved analytically to obtain
 * \f{align}{
 *      A = & -\frac{10 \alpha^{2}-8 \alpha +1}{6 \alpha^{3} \mathtt{cutoff}^{4} \left(\alpha -1\right) \left(\alpha^{2}-2 \alpha +1\right)} \label{eq:cutoff_paramsA}\\
 *      B = &  \frac{2 \alpha -1}{2 \alpha^{3} \mathtt{cutoff}^{5} \left(\alpha^{3}-3 \alpha^{2}+3 \alpha -1\right)} \label{eq:cutoff_paramsB}\\
 *      \mathtt{shift} = & \frac{10 \alpha^{2}-5 \alpha +1}{6 \alpha^{3} \mathtt{cutoff}}
 * \label{eq:cutoff_paramsS}
 * \f}
 *
 * The cutoff electrostatics computation is invoked by the following directive in `.control` file:
 * ```(control)
 * elstat cutoff cutoff [alpha=0.7]
 * ```
 * where the first `cutoff` sets the method and the second `cutoff` should be a numerical value. Setting @ref cutoff is obligatory. If no @ref alpha value is given, the default value 0.7 is used.
 * 
 * Since version 1.4, the shift can be turned off in periodic boundary conditions to obtain various version of cutoff behaviour.
 * The following options in `.control` file can be used:
 * 1. `no elstatshift` – any shift and smoothing is turned off, @ref alpha is set to 1 regardless to its initial value
 *     and @ref shift is set to 0. 
 *     Leads to normal Coulomb electrostatics that in most cases does not converge.
 * 2. `no elstatforceshift` – only shift in energy, @ref alpha is set to 1 
 *     and @ref shift is set to \f$1/\mathtt{cutoff}\f$.
 *     Leads to jump in forces.
 * 3. `no elstatenergyshift` – does not make much sense. 
 *     Energy is not shifted, but the tail is made smooth by the transition function.
 *     Continuous potential and forces, but not the second derivative (only two parameters @ref A and @ref B left).
 *     Weird behaviour of forces (non-monotonous). Do not use... 
 *
 * The main purpose of this class, i.e. the calculation of forces and potential energy, is done by implementing @ref CalculateEpot() const which yields only the potential energy of the electrostatic interactions and
 * @ref CalculateForces() which also modifies the forces for each @ref Atom in the SimulatedSystem.
 */
class CutoffElstat : public AbstractInterMolField
{
private:
  /**
   * @brief Cutoff parameter \f$A\f$.
   * Cutoff parameter calculated to ensure smoothness of the potential. See \f$\eqref{eq:cutoff_paramsA}\f$.
   */
  double A;
  /**
   * @brief Cutoff parameter \f$B\f$.
   * Cutoff parameter calculated to ensure smoothness of the potential. See \f$\eqref{eq:cutoff_paramsB}\f$.
   */
  double B;
  /**
   * @brief Cutoff parameter `shift`.
   * Cutoff parameter calculated to ensure smoothness of the potential. See \f$\eqref{eq:cutoff_paramsS}\f$.
   */
  double shift;
  /**
   * @brief User-defined cutoff parameter \f$\alpha\f$.
   * Shifted Coulomb potential is used for separations shorter than \f$\alpha \mathtt{cutoff}\f$, intermediate interpolation is used in between \f$\alpha \mathtt{cutoff}\f$ and `cutoff`.
   * See \f$\eqref{eq:elstat_cut_shift}\f$.
   */
  double alpha;
  /**
   * @brief Default constructor.
   *
   * Made private in order not to be used because an empty electrostatic system should not be instantiated.
   */
  CutoffElstat(); // making default constructor private disables its use
  /**
   * @brief If *true*, the system does not contribute to potential energy (all charges \f$q_i\f$ are zero).
   */
  bool empty;
  /**
   * @brief Shift to ensure energy (if shift & 2) and forces (if shift & 1) continuity (default 3 (both apply)).
   * Saving this value to the class private variable enables to reuse it for @ref SetCutoff(). 
   */
  unsigned int sft;

public:
  /**
   * @brief Primary constructor, sets the parent SimulatedSystem
   *
   * @param parent Parent SimulatedSystem (see @ref AbstractInterMolField::parent).
   * @param cut Electrostatic potential cutoff (see @ref AbstractInterMolField::cutoff).
   * @param alp  Intermediate cutoff parameter (see @ref alpha), default value 0.7 (recommended 0.6–0.7).
   * @param shift Shift to ensure energy (if shift & 2) and forces (if shift & 1) continuity (default 3 (both apply)).
   */
  CutoffElstat(SimulatedSystem *parent, double cut, double alp = 0.7, unsigned int shift = 3); // constructor
  /**
   * @brief CutoffElstat copy constructor
   *
   * @param system CutoffElstat class to be copied.
   */
  CutoffElstat(const CutoffElstat &system); // copy constructor
  /**
   * @brief CutoffElstat destructor
   *
   * Empty, because no allocation is needed for this class.
   */
  ~CutoffElstat(); // destructor
  /**
   * @brief Calculate forces caused by electrostatic interactions (intermolecular).
   * Calculate forces caused by this field (CutoffElstat) and add them to @ref Atom::force.
   * Save the potential energy and force virial to the @ref parent SimulatedSystem @ref SimulatedSystem::Epotintermol[0] and @ref SimulatedSystem::virintermol[0].
   *
   * The potential energy caused by electrostatic interactions (cut, shifted and smoothed):
   * \f{equation}{
   *      E_{\mathrm{cutel}} = \sum_{i,j} \begin{cases}
   *          \frac{q_i q_j}{4\pi \epsilon_0} \left(\frac{1}{r_{ij}} - \mathtt{shift}\right) & r_{ij} < \alpha \mathtt{cutoff} \\
   *          \frac{q_i q_j}{4\pi \epsilon_0} \left[\left(r_{ij}-\mathtt{cutoff}\right)^3(A+Br_{ij})\right] & \alpha \mathtt{cutoff} \leq r_{ij} < \mathtt{cutoff} \\
   *          0 & r_{ij} \leq \mathtt{cutoff} \\
   *          \end{cases}
   * \f}
   * where \f$r_{ij}\f$ is the 2-norm (see @ref Vector::CalculateNorm(int) const) of the vector \f$\mathbf{r}_{ij} = \mathbf{r}_j - \mathbf{r}_i\f$, where \f$\mathbf{r}_i\f$ is the position of Atom \f$i\f$. Intermediate cutoff parameters \f$A\f$, \f$B\f$ and \f$\mathtt{shift}\f$ are determined from the value of `cutoff` and `alpha` (see @ref AbstractInterMolField::cutoff and @ref alpha).
   * The sum runs over all pairs of Atoms in the SimulatedSystem, except for Atoms belonging to the same Molecule.
   * For intramolecular electrostatic interactions, the corresponding term must be specified as an intramolecular field (see @ref AbstractIntraMolField).
   * The potential energy is added to the parent SimulatedSystem::Epotintermol[0] (electrostatic terms).
   * In case of free periodic boundary conditions, the energy is the same as if \f$\mathtt{shift}\rightarrow0\f$ and \f$\mathtt{cutoff}\rightarrow\infty\f$.
   *
   * Forces (3D vector) are given by the derivative of the potential energy:
   * \f{equation}{
   *     \mathbf{f}_j = -\mathbf{f}_i =
   *          \begin{cases}
   *              \frac{q_i q_j}{4\pi \epsilon_0} \frac{1}{r_{ij}^2}\frac{\mathbf{r}_{ij}}{r_{ij}} & \text{if } r_{ij} \leq \alpha \mathtt{cutoff} \\
   *              \frac{q_i q_j}{4\pi \epsilon_0} \left(r_{ij}-\mathtt{cutoff}\right)^2\left(3A - B\mathtt{cutoff} + 4Br_{ij}\right)\frac{\mathbf{r}_{ij}}{r_{ij}} & \text{if } \alpha\mathtt{cutoff} < r_{ij} < \mathtt{cutoff} \\
   *               0 & \text{otherwise}
   *          \end{cases}
   * \f}
   *
   * For precise (not cut) electrostatics, the virial is exactly equal to potential energy.
   * While the agreement is good for methods using Ewald summation, here, the virial theorem does not apply.
   * Better values of pressure are obtained using the virial of forces.
   * The contribution to the viral by each interaction pair is:
   * \f{equation}{
   *     w_{\mathrm{LJ}} = \mathbf{f}_i \cdot \mathbf{r}_{ij} =
   *         \begin{cases}
   *              -\frac{q_i q_j}{4\pi \epsilon_0} \frac{1}{r_{ij}} & \text{if } r_{ij} \leq \alpha \mathtt{cutoff} \\
   *              -\frac{q_i q_j}{4\pi \epsilon_0} \left(r_{ij}-\mathtt{cutoff}\right)^2\left(3A - B\mathtt{cutoff} + 2Br_{ij}\right)r_{ij} & \text{if } \alpha\mathtt{cutoff} < r_{ij} < \mathtt{cutoff} \\
   *               0 & \text{otherwise}
   *         \end{cases}
   * \f}
   * This value is added to the parent @ref SimulatedSystem::virintermol[0].
   *
   * Returns the potential energy.
   *
   * @par See also
   * @ref CalculateEpot() const.
   *
   * @return Potential energy caused by cut-and-shifted electrostatic interactions in the @ref parent @ref SimulatedSystem (internally in program units).
   */
  double CalculateForces() override; // calculate forces in given simSystem
  /**
   * @brief Calculate potential energy caused by electrostatic interactions (intermolecular).
   * Similar to @ref CalculateForces(), but calculates only the potential energy and does not change anything neither in the @ref parent SimulatedSystem, nor in the involved Atoms.
   *
   * The potential energy caused by electrostatic interactions (cut, shifted and smoothed):
   * \f{equation}{
   *      E_{\mathrm{cutel}} = \sum_{i,j} \begin{cases}
   *          \frac{q_i q_j}{4\pi \epsilon_0} \left(\frac{1}{r_{ij}} - \mathtt{shift}\right) & r_{ij} < \alpha \mathtt{cutoff} \\
   *          \frac{q_i q_j}{4\pi \epsilon_0} \left[\left(r_{ij}-\mathtt{cutoff}\right)^3(A+Br_{ij})\right] & \alpha \mathtt{cutoff} \leq r_{ij} < \mathtt{cutoff} \\
   *          0 & r_{ij} \leq \mathtt{cutoff} \\
   *          \end{cases}
   * \f}
   * where \f$r_{ij}\f$ is the 2-norm (see @ref Vector::CalculateNorm(int) const) of the vector \f$\mathbf{r}_{ij} = \mathbf{r}_j - \mathbf{r}_i\f$, where \f$\mathbf{r}_i\f$ is the position of Atom \f$i\f$. Intermediate cutoff parameters \f$A\f$, \f$B\f$ and \f$\mathtt{shift}\f$ are determined from the value of `cutoff` and `alpha` (see @ref AbstractInterMolField::cutoff and @ref alpha).
   * The sum runs over all pairs of Atoms in the SimulatedSystem, except for Atoms belonging to the same Molecule.
   * For intramolecular electrostatic interactions, the corresponding term must be specified as an intramolecular field (see @ref AbstractIntraMolField).
   * In case of free periodic boundary conditions, the energy is the same as if \f$\mathtt{shift}\rightarrow0\f$ and \f$\mathtt{cutoff}\rightarrow\infty\f$.
   *
   * @par See also
   * @ref CalculateForces().
   *
   * @return Potential energy caused by cut-and-shifted electrostatic interactions in the @ref parent @ref SimulatedSystem (internally in program units).
   */
  double CalculateEpot() const override; // calculate potential energy in given simSystem
  /**
   * @brief Calculate intermediate cutoff parameters \f$A\f$, \f$B\f$ and `shift`.
   * Intermediate cutoff parameters are calculated in such a way that the potential, force and its derivative are continuous.
   * The system of three equations is solved
   * \f{align}{
   *    \frac{1}{\alpha  \mathtt{cutoff}}-\mathtt{shift} &=
   *        \mathtt{cutoff}^3 \left(\alpha - 1\right)^{3} \left(B \alpha  \mathtt{cutoff} +A \right) \\
   *  -\frac{1}{\alpha^{2} \mathtt{cutoff}^{2}} &=
   *        3 \mathtt{cutoff}^2 \left(\alpha - 1\right)^{2} \left(B \alpha  \mathtt{cutoff} +A \right)+\mathtt{cutoff}^3\left(\alpha - 1\right)^{3} B \\
   *  \frac{2}{\alpha^{3} \mathtt{cutoff}^{3}} &=
   *        6 \mathtt{cutoff}\left(\alpha-1\right) \left(B \alpha  \mathtt{cutoff} +A \right)+6\mathtt{cutoff}^2 \left(\alpha -1\right)^{2} B \\
   * \f}
   * The continuity at `cutoff` is guaranteed by the form of the potential.
   *
   * Analytical solution can be obtained:
   * \f{align}{
   *      A = & -\frac{10 \alpha^{2}-8 \alpha +1}{6 \alpha^{3} \mathtt{cutoff}^{4} \left(\alpha -1\right) \left(\alpha^{2}-2 \alpha +1\right)} \label{eq:cutoff_paramsA2}\\
   *      B = &  \frac{2 \alpha -1}{2 \alpha^{3} \mathtt{cutoff}^{5} \left(\alpha^{3}-3 \alpha^{2}+3 \alpha -1\right)} \label{eq:cutoff_paramsB2}\\
   *      \mathtt{shift} = & \frac{10 \alpha^{2}-5 \alpha +1}{6 \alpha^{3} \mathtt{cutoff}}
   * \f}
   *
   * The calculated values are saved in @ref A, @ref B and @ref shift.
   */
  void CalculateABshift(); // calculate intermediate cutoff parameters
  /**
   * @brief Set the @ref cutoff
   * Set the distance beyond which the interaction is zero.
   * Provides check if cutoff is less than half the shortest box size.
   *
   * @param vdwcut Cutoff distance for disperse interactions (in [AA]) – irrelevant for electrostatics
   * @param elcut Cutoff distance for electrostatic interactions (in [AA])
   * @param boundaryC Boundary conditions (see @ref SimulatedSystem::boundaryCond)
   * @return 0 if success
   */
  int SetCutoff(double vdwcut, double elcut, int boundaryC) override; // set electrostatic interactions spherical cutoff
  /**
   * @brief Clone itself and get adopted by another @ref parent SimulatedSystem.
   * Returns a copy of itself and migrate to the @p newparent SimulatedSystem.
   * Useful when copying SimulatedSystem object.
   * Every intermolecular term is then cloned by this method and assigned to the new SimulatedSystem object.
   *
   * @par See also
   * @ref CutoffElstat(const CutoffElstat&), @ref SetParent(SimulatedSystem*), @ref AbstractInterMolField::copy(SimulatedSystem *) const, @ref AbstractIntraMolField::copy(Molecule *) const.
   *
   * @param newparent Parent SimulatedSystem of the newly created CutoffElstat object.
   * @return (Pointer to) the newly created CutoffElstat object.
   */
  CutoffElstat *copy(SimulatedSystem *newparent) const override; // cloning itself
  /**
   * @brief Cutoff corrections arising from cutting the electrostatic potential are assumed to be zero.
   * Generally we want to avoid consequences of cutting the potential after certain distance.
   * Cutoff corrections are usually given by
   * \f{equation}{
   *      E_\mathrm{corr} = 2\pi \sum_{i,j} \frac {1}{V} \int_{0}^{\infty}\! \left(E_{\mathrm{id}} - E_{\mathrm{sim}}\right) g(r) r^2\,\mathrm{d}r
   * \f}
   * where \f$E_{\mathrm{id}}\f$ is the original (desired) potential (in our case original Coulomb potential), \f$E_{\mathrm{sim}}\f$ is the truncated potential used in the simulation and \f$g(r)\f$ is the radial distribution function and the sum runs over all charged atoms in the system.
   * Essentially, we must integrate the difference between the real and desired potential for all pairs in the system.
   *
   * In free (vacuum) boundary conditions, the interactions can be calculated explicitely and no cutoff is applied.
   * Therefore the system can be charged and it is still reasonable to simulate it.
   * However, in periodic boundary conditions, the system must be neutral not to have infinite energy.
   * Furthermore, the integral used for correction calculation diverges when the energy decreases as \f$\mathcal{O}(1/r)\f$.
   * The use of cutoff electrostatics is only justified when either each molecule is neutral as a whole (having permanent dipole, quadrupole, etc.) or
   * the cutoff is longer then the Debye screening length of the system (which, nevertheless, cannot be calculated in advance as we do not know the value of dielectric constant) @cite kolafaskriptaCUni pp 81–82.
   * In all cases the cutoff distance should be long enough not to introduce big errors (MACSIMUS manual recommends 15 Å as a minimum, @cite MACSIMUS_manual p 206).
   *
   * Consequently this method returns zero energy correction.
   *
   * @return 0.
   */
  double CalculateCutoffCorrection() override; // calculate cutoff correction energy (returns zero)
  /**
   * @brief Initialize fake *bonds* used to handle electrostatic interactions between atoms of one Molecule.
   *
   * Electrostatic interactions are basically intermolecular, but they occur also in longer molecules.
   * They are realized by bonds of class @ref ElstatBond.
   * Atoms on distances 1–4 interact by scaled Coulombic interaction,
   * atoms which are furhter from each other interact by original Coulombic interaction.
   *
   * @param LJ14factor Scaling factor for 1–4 disperse interactions.
   * @param el14factor Scaling factor for 1–4 electrostatic interactions.
   * @return Number of newly initialized bonds.
   */
  int InitializeIntramolBonds(double LJ14factor = 1.0, double el14factor = 1.0) const override; // initialize fake *bonds* used to handle this intermolecular field between atoms in one molecule
  /**
   * @brief Print informations about the cutoff electrostatics to the @p stream
   *
   * @param stream Output stream.
   * @param u_eng Energy unit value (see `units.hpp`).
   * @param engunit Energy unit name in brackets (as in `.cpa` file).
   */
  void PrintInfo(std::ofstream &stream, double u_eng, std::string engunit) const override;
};

#endif