#ifndef EWALDSPLINESELSTATHEADER
#define EWALDSPLINESELSTATHEADER

#define MAX_K 20

#include "AbstractInterMolField.hpp"
#include "SplineType.hpp"
#include "AbstractInterpolator.hpp"
#include "HermiteCubSplines.hpp"
#include "LinearInterpolator.hpp"
#include "Linear3PInterpolator.hpp"
#include "Linear4PInterpolator.hpp"
#include "MacsimusHyperbSplines.hpp"
#include "MacsimusQuadSplines.hpp"
#include "NaturalCubSplines.hpp"


class Matrix;
class Atom;

#ifdef PARALLEL
extern int thread_count;
#endif

/**
 * @class EwaldSplinesElstat
 * @ingroup InterMolField
 * @brief Standard Ewald summation for charged systems. Includes both the \f$r\f$-space and \f$k\f$-space sums. Uses interpolation for \f$r\f$-space sums.
 *
 * @note This class is equivalent to @ref EwaldElstat, the only change is the use of interpolation instead of erfc function.
 *
 * This class implements the interface for intermolecular interactions @ref AbstractInterMolField ( @ref type = 0).
 * Electrostatic interactions are given by the charge of Atoms (see @ref Atom::charge).
 * The original Coulomb interaction energy is given by
 * \f{equation}{
 *     E_{\mathrm{Coul}} = \sum_{i,j} \frac{q_i q_j}{4\pi \epsilon_0} \frac{1}{r_{ij}}
 * \f}
 * where \f$r_{ij}\f$ is the 2-norm (see @ref Vector::CalculateNorm(int) const) of the vector \f$\mathbf{r}_{ij} = \mathbf{r}_j - \mathbf{r}_i\f$, where \f$\mathbf{r}_i\f$ is the position of Atom \f$i\f$.
 * Charge of Atom \f$i\f$ (see @ref Atom::charge) is denoted as \f$q_i\f$.
 * The sum runs over all pairs of Atoms in the SimulatedSystem and their periodic images, except for Atoms belonging to the same Molecule.
 * For intramolecular interactions, the corresponding term must be specified as an intramolecular field (see @ref AbstractIntraMolField, @ref ElstatBond).
 * In free boundary conditions, the @ref CutoffElstat class must be used to calculate all-pairs Coulombic interactions (cutoff does not apply in free b.c.).
 *
 * In Ewald summation the above described infinite sum is divided into real \f$r\f$-space sum and a sum in a reciprocal space (\f$k\f$-space sum).
 * The main idea can be explained in two ways:
 * 1. **Physical approach** @cite dl_poly4-manual, @cite AllenCompSimLiq
 *     A gaussian distributed screening charge is added with the opposite sign to allow rapid convergence of the real-space sum.
 *     To obtain results for the original distribution, interactions among the screening charge clouds must be subtracted, this is done in reciprocal space,
 *     where the sum converges.
 *     The Fourier transform is involved.
 * 2. **Mathematical approach** @cite kolafaskriptaCUni
 *     The Ewald method can also be explained by purely mathematical arguments.
 *     It uses the trivial (but not obviously useful) identity:
 *     \f{equation}{
 *         \frac{1}{r} = \frac{2}{\sqrt{\pi}}\int_{0}^\infty \! \exp\left(-t^2 r^2\right)\, \mathrm{d}t = \frac{2}{\sqrt{\pi}}\int_{0}^\alpha \! \exp\left(-t^2 r^2\right)\, \mathrm{d}t + \frac{2}{\sqrt{\pi}}\int_{\alpha}^\infty \! \exp\left(-t^2 r^2\right)\, \mathrm{d}t
 *     \f}
 *     The second term leads to complementary error function and rapidly converges in real-space, whereas the first is useful on application of Poisson summation rule and definitely leads to the sum in reciprocal space.
 *
 * Eventually, the potential energy can be computed from @cite dl_poly4-manual
 * \f{equation}{
 *     E_\mathrm{Coul} = \frac{1}{2V\epsilon} \sum_{\mathbf{k}\neq \mathbf{0}}^\mathbf{\infty} \frac{\exp\left(-k^2 / (4\alpha^2)\right)}{k^2}\left|\sum_j^{N_\mathrm{ions}} q_j \exp\left(-\mathrm{i} \mathbf{k}\cdot\mathbf{r_j}\right) \right|^2 + \frac{1}{4\pi\epsilon}\sum_{i,j}^{N_\mathrm{pairs}}\frac{q_i q_j}{r_{ij}}\mathrm{erfc}\left(\alpha r_{ij}\right) - \frac{1}{4\pi\epsilon}\frac{\alpha}{\sqrt{\pi}} \sum_j^{N_\mathrm{ions}} q_j^2 - \frac{1}{4\pi\epsilon} \sum_{\mathrm{molecules}} \sum_{i,j}^{N_\mathrm{inmol}} q_i q_j \frac{\mathrm{erf}\left(\alpha r_{ij}\right)}{r_{ij}}
 * \f}
 * where \f$\alpha\f$ is the screening parameter, \f$N_{\mathrm{ions}}\f$ is the total number of ions (atoms),
 * \f$N_\mathrm{pairs}\f$ means the sum over all pairs of atoms from different molecules, \f$N_\mathrm{inmol}\f$ means summing over all intramolecular pairs.
 * The first term is the reciprocal-space sum, the second is the real-space sum and the last two terms are corrections to exclude self-interactions and intramolecular interactions (which are treated separately in `simul++`).
 * Both the \f$k\f$-space sum and the \f$r\f$-space sum are truncated.
 * Real cutoff must be shorter than one half of the shortest boxsize.
 * Pairs of atoms further than the cutoff distance are not used in the \f$r\f$-space sum.
 * This class uses interpolation by an @ref AbstractInterpolator to speed-up calculation of the erfc() function.
 *
 * Here, we assume the so-called *tin-foil* boundary conditions in infinity.
 * Otherwise the term describing the contribution of the cell dipole moment must have been added (see @cite MACSIMUS_manual).
 *
 * The reciprocal-space vector \f$\mathbf{k}\f$ is given by
 * \f{equation}{
 *     \mathbf{k} = \left[ l \frac{2\pi}{L_0}, m \frac{2\pi}{L_1}, n \frac{2\pi}{L_2}\right]^\mathrm{T}
 * \f}
 * where \f$L_0\f$, \f$L_1\f$ and \f$L_2\f$ are boxsizes and \f$l\f$, \f$m\f$ and \f$n\f$ are integers.
 *
 * The accuracy and convergence of the Ewald method is determined by three parameters.
 * The first one is real-space cutoff ( @ref cutoff), the second is the screening (separation, convergence) factor \f$\alpha\f$ ( @ref alpha) and the third is the reciprocal space cutoff \f$\kappa\f$.
 * The reciprocal-space sum runs over vectors \f$k = \left|\left|\mathbf{k}\right|\right| < 2\pi\kappa \f$.
 * Indices \f$l\f$, \f$m\f$ and \f$n\f$ range from \f$-\kappa L_i\f$ to \f$\kappa L_i\f$ where i is 0 for \f$l\f$ etc. but not all of them are involved in the final sum because the sum is limited to sphere not cube (or cuboid).
 * Nevertheless, this truncation at sphere is used only when calculating forces (see @ref CalculateForces()), for @ref CalculateEpot() (and thus for pressure from virtual volume change (PVVC)) the spherical truncation is not done since it leads to wrong values of PVVC.
 * Therefore, calculation of PVVC is time-consuming and should be used only for initial debugging and tuning, not for productive runs.
 * The recommended values for convergence parameters are @cite MACSIMUS_manual @cite dl_poly4-manual
 * \f{equation}{
 *     \alpha = \kappa = \frac{\pi}{\mathtt{cutoff}}
 * \f}
 * Large values of \f$\alpha\f$ lead to slower convergence in reciprocal-space, lower values of \f$\alpha\f$ lead to slower convergence in real-space.
 * The accuracy is also affected by the used interpolators and their grid density.
 *
 * The Ewald summation for electrostatic interactions with interpolation is invoked by the following directives in `.control` file:
 * ```(control)
 * elstat ewalds cutoff [alpha kappa]
 * splines type gridsize [type gridsize]
 * ```
 * Setting @ref cutoff is obligatory. If no @ref alpha value is given, the default value \f$\pi/\mathtt{cutoff}\f$ is used.
 * The same applies for @ref kappa (default \f$\kappa=\pi/\mathtt{cutoff}\f$).
 * The directive `splines` determines the type and density of the used interpolation scheme,
 * where `type` is one of the following `hyperbolic`, `quadratic`, `natural`, `hermite`, `linear`, `linear3` or `linear4`
 * and `gridsize` is the number of grid points per unity in argument (squared distance) used for interpolation.
 * If `splines` directive is not present, the default value `splines hyperbolic 256` is used.
 *
 * From version 1.4, both energy and forces in real-space are shifted (should be only a small shift) to avoid discontinuities.
 * This is done by subtracting the value of @ref shiftE from energy in \f$r\f$-space and
 * @ref shiftf from forces from \f$r\f$-space.
 * Both can be turned off to mimic the behaviour of the older versions by this `.control` directives:
 * ```(control)
 * no elstatshift       # turns off both (shiftE = 0.0, shiftf = 0.0)
 * no elstatenergyshift # turns off energy shifting (shiftE = 0.0), forces continuous
 * no elstatforceshift  # turns off forces shifting (shiftf = 0.0), energy continuous
 * ```
 *
 * The main purpose of this class, i.e. the calculation of forces and potential energy, is done by implementing @ref CalculateEpot() const which yields only the potential energy of the electrostatic interactions and
 * @ref CalculateForces() which also modifies the forces for each @ref Atom in the SimulatedSystem.
 *
 * Parallelization in r-space is done using by some pairlist (see @ref AbstractPairList),
 * parallelization in k-space is done by splitting the `for` loop over @ref charged_atoms to different threads.
 * To avoid atomic access to @ref imag and @ref real, both are allocated larger to have one copy per thread,
 * after the terms calculation, per-thread values are summed up.
 */
class EwaldSplinesElstat : public AbstractInterMolField
{
private:
  /**
   * @brief Ewald parameter \f$\alpha\f$.
   * Parameter determining the width of the screening Gaussian charges and convergence.
   */
  double alpha;
  /**
   * @brief Ewald parameter \f$\alpha\f$ squared and multiplied by 4.
   * The expression \f$4\alpha^2\f$ is used so frequently, that it would be useless to compute it over again.
   */
  double alpha42;
  /**
   * @brief Ewald parameter \f$\kappa\f$.
   * Parameter determining the maximum size of the reciprocal vector used in the \f$k\f$-space sum.
   */
  double kappa;
  /**
   * @brief Self-correction energy \f$E_{\mathrm{self-corr}}\f$.
   * Correction term cancelling the self-interactions in the reciprocal space.
   */
  double Eselfcorr;
  /**
   * @brief Shift to ensure energy (if shift & 2) and forces (if shift & 1) continuity (default 3 (both apply)).
   * Saving this value to the class private variable enables to reuse it for @ref SetCutoff().
   */
  unsigned int sft;
  /**
   * @brief Energy shift to assure continuity in r-space energy at cutoff
   * Calculated by
   * \f{equation}{
   *     E_{\mathrm{shift}} = \frac{1}{4\pi\epsilon}\frac{\mathrm{erfc}\left(\alpha \mathtt{cutoff}\right)}{\mathtt{cutoff}}
   * \f}
   * The value is multiplied by the factor of charges and subtracted from the value of r-space energy.
   */
  double shiftE;
  /**
   * @brief Shift in forces to assure continuity in r-space forces at cutoff
   * \f{equation}{
   *     f_{\mathrm{shift}} =  \frac{1}{4\pi\epsilon} \frac{2\alpha \mathtt{cutoff}\exp\left(-\alpha^2 \mathtt{cutoff}^2\right) + \sqrt{\pi}\mathrm{erfc}\left(\alpha \mathtt{cutoff}\right)}{\sqrt{\pi}\mathtt{cutoff}^2}
   * \f}
   * The value is multiplied by the factor of charges and subtracted from the value of r-space forces.
   */
  double shiftf;
  /**
   * @brief Vector \f$\mathbf{k}\f$ stores values for given index multiplied by box length \f$L_x\f$.
   * Three arrays ( @ref kx, @ref ky and @ref kz) are needed to enable rectangular box.
   */
  double kx[MAX_K + 1];
  /**
   * @brief Vector \f$\mathbf{k}\f$ stores values for given index multiplied by box length \f$L_y\f$.
   * Three arrays ( @ref kx, @ref ky and @ref kz) are needed to enable rectangular box.
   */
  double ky[MAX_K + 1];
  /**
   * @brief Vector \f$\mathbf{k}\f$ stores values for given index multiplied by box length \f$L_z\f$.
   * Three arrays ( @ref kx, @ref ky and @ref kz) are needed to enable rectangular box.
   */
  double kz[MAX_K + 1];
  /**
   * @brief Precomputed term used in energy calculation.
   * For each vector \f$\mathbf{k}\f$ (with positive indices) the value
   * \f{equation}{
   *     \frac{\exp \left(-k^2/\left(4\alpha^2\right)\right)}{k^2}
   * \f}
   * is stored.
   * This array is used only inside constructors and @ref CalculateForces() method, which is also responsible for any updates due to the change of box size.
   *
   * The (originally and logically 3-dim) array is stored as 1-dim dynamically allocated array to enhance performance.
   * Values of @ref ymult and @ref zmult can be used to index the array.
   * Symbolically: `expk[i][j][k] = expk[i*ymult + j*zmult + k]`.
   */
  double *expk; // [(ymult) * (lmax + 1)];
  /**
   * @brief Maximum value of index \f$l\f$.
   * See the introduction to this class.
   */
  int lmax;
  /**
   * @brief Maximum value of index \f$m\f$.
   * See the introduction to this class.
   */
  int mmax;
  /**
   * @brief Maximum value of index \f$n\f$.
   * See the introduction to this class.
   */
  int nmax;
  ///@{
  /**
   * @brief Auxiliary arrays of maximum indeces of \f$k\f$-vector
   *
   * For each \f$l\f$, the value of maximum \f$m\f$ is stored; the same for \f$m\f$ and \f$n\f$.
   * To avoid evaluating of maximum index in each loop and to conserve the same sphere when volume fluctuates.
   * @note Volume should not fluctuate too much, otherwise the number of considered \f$k\f$-space vectors is not correct.
   */
  int *kymax;
  int *kzmax;
  ///@}
  ///@{
  /**
   * @brief The imaginary ( @ref imag) and real ( @ref real) parts of the sum in the reciprocal space.
   * The two arrays @ref real and @ref imag store the real and the imaginary part of the sum
   * \f{equation}{
   *     \sum_j^{N_\mathrm{ions}} q_j \exp\left(-\mathrm{i} \mathbf{k}\cdot\mathbf{r_j}\right)
   * \f}
   * for each vector \f$\mathbf{k}\f$.
   *
   * The array is stored as 1-dim dynamically allocated array of size @ref xmult to enhance performance.
   * In parallel version, each thread has its own part of array and the size is thus @ref xmult * thread_count.
   * Values are saved in a *weird* order optimizing the cache usage (inspired by MACSIMUS implementation).
   * Schematically the order is: (+, 0, 0), (+, +, 0), (+, −, 0), (+, +, +), (+, +, −), (+, −, +), (+, −, −).
   * For further information, see the implementation.
   */
  double *imag; // [xmult * thread_count];
  double *real; // [xmult * thread_count];
  ///@}
  /**
   * @brief Array of pointers to all charged Atoms in the SimulatedSystem.
   * Defined to enhance computation of \f$k\f$-space sum.
   */
  Atom **charged_atoms;
  ///@{
  /**
   * @brief Arrays storing sine ( @ref sins) and cosine ( @ref coss) terms of the complex exponential \f$\exp\left(-\mathrm{i} \mathbf{k}\cdot\mathbf{r_j}\right)\f$.
   *
   * The (originally and logically 4-dim) array is stored as 1-dim dynamically allocated array to enhance performance.
   * The first index runs over charged atoms, the others over vectors in \f$k\f$-space.
   * The order in \f$k\f$-space is the same as for @ref imag and @ref real.
   * Symbolically: `sins[m][k-space] = sins[m*xmult + index_in_k-space]`.
   */
  double *sins;
  double *coss;
  ///@}
  /**
   * @brief Number of charged Atoms in the SimulatedSystem.
   * The size of array @ref charged_atoms.
   * Replaces also the original `empty` variable (test if electrostatics contribute to the energy and forces).
   */
  int no_ch_atoms;
  AbstractInterpolator<double> *Erfc;
  AbstractInterpolator<double> *Derfc;
  /**
   * @brief Default constructor.
   *
   * Made private in order not to be used because an empty electrostatic system should not be instantiated.
   */
  EwaldSplinesElstat(); // making default constructor private disables its use
  ///@{
  /**
   * @brief Integers used to index @ref expk and find the size of arrays needed, determined by the size of the \f$k\f$-space.
   *
   * It holds: @ref zmult = ( @ref nmax + 1), @ref ymult = @ref zmult * ( @ref mmax + 1)
   * and @ref xmult = ( @ref nmax + 1) * ( @ref mmax + 1) * ( @ref lmax + 1) * 4 + ( @ref lmax + 1) * ( @ref mmax + 1) * 2 + ( @ref lmax + 1) - 4
   */
  int xmult;
  int ymult;
  int zmult;
  ///@}

public:
  /**
   * @brief Primary constructor, sets the parent SimulatedSystem and Ewald parameters.
   *
   * The real space cutoff ( @ref cutoff) is set.
   * If @p alp is not given (default value), the value of screening parameter @ref alpha is set to
   * \f{equation}{
   *     \alpha = \frac{\pi}{\mathtt{cutoff}}
   * \f}
   * Similarly, if @p kap is not given (default value), the value of reciprocal space cutoff @ref kappa is set to
   * \f{equation}{
   *     \kappa = \frac{\pi}{\mathtt{cutoff}}
   * \f}
   *
   * The self-energy correction is calculated by
   * \f{equation}{
   *     E_{\mathrm{self-corr}} = - \frac{1}{4\pi\epsilon}\frac{\alpha}{\sqrt{\pi}} \sum_j^{N_\mathrm{ions}} q_j^2
   * \f}
   * and saved to @ref Eselfcorr.
   *
   * The constructor is also responsible for allocation of arrays @ref kx, @ref ky, @ref kz, @ref expk, @ref real and @ref imag.
   * Values are feeded to the first four of them, the last two must be updated during energy calculation as they depend on ion positions.
   *
   * @param parent Parent SimulatedSystem (see @ref AbstractInterMolField::parent).
   * @param cut Electrostatic potential cutoff (see @ref AbstractInterMolField::cutoff).
   * @param alp  Screening parameter \f$\alpha\f$ (see @ref alpha).
   * @param kap Ewald parameter \f$\kappa\f$ determining the cutoff in reciprocal space.
   * @param shift Shift to ensure energy (if shift & 2) and forces (if shift & 1) continuity (default 3 (both apply)).
   * @param type Type of splines used to calculate values of \f$\mathrm{erfc}(\alpha r)/r\f$ (see @ref SplineType).
   * @param gridsize Grid for interpolation by splines. Points per Angström.
   * @param type2 Type of splines used to calculate derivatives of \f$\mathrm{erfc}(\alpha r)/r\f$ (see @ref SplineType). Used only if @p gridsize2 != 0.
   * @param gridsize2 As @p gridsize, but for derivatives. If 0, then derivatives are calculated from splines for values.
   */
  EwaldSplinesElstat(SimulatedSystem *parent, double cut, double alp = -1.0, double kap = -1.0, unsigned int shift = 3, SplineType type = hyperbolic, unsigned int gridsize = 256, SplineType type2 = hyperbolic, unsigned int gridsize2 = 0); // constructor
  /**
   * @brief EwaldSplinesElstat copy constructor
   *
   * @param system EwaldSplinesElstat class to be copied.
   */
  EwaldSplinesElstat(const EwaldSplinesElstat &system); // copy constructor
  /**
   * @brief EwaldSplinesElstat destructor
   *
   * The destructor is responsible for deallocation of arrays @ref kx, @ref ky, @ref kz, @ref expk, @ref real and @ref imag.
   */
  ~EwaldSplinesElstat(); // destructor
  /**
   * @brief Calculate forces caused by electrostatic interactions (intermolecular).
   * Calculate forces caused by this field (EwaldSplinesElstat) and add them to @ref Atom::force.
   * Save the potential energy and force virial to the @ref parent SimulatedSystem @ref SimulatedSystem::Epotintermol[0] and @ref SimulatedSystem::virintermol[0].
   *
   * The potential energy caused by electrostatic interactions calculated by Ewald method reads
   * Eventually, the potential energy can be computed from @cite dl_poly4-manual
   * \f{equation}{
   *     E_\mathrm{Coul} = \frac{1}{2V\epsilon} \sum_{\mathbf{k}\neq \mathbf{0}}^\mathbf{\infty} \frac{\exp\left(-k^2 / (4\alpha^2)\right)}{k^2}\left|\sum_j^{N_\mathrm{ions}} q_j \exp\left(-\mathrm{i} \mathbf{k}\cdot\mathbf{r_j}\right) \right|^2 + \frac{1}{4\pi\epsilon}\sum_{i,j}^{N_\mathrm{pairs}}\frac{q_i q_j}{r_{ij}}\mathrm{erfc}\left(\alpha r_{ij}\right) - \frac{1}{4\pi\epsilon}\frac{\alpha}{\sqrt{\pi}} \sum_j^{N_\mathrm{ions}} q_j^2 - \frac{1}{4\pi\epsilon} \sum_{\mathrm{molecules}} \sum_{i,j}^{N_\mathrm{inmol}} q_i q_j \frac{\mathrm{erf}\left(\alpha r_{ij}\right)}{r_{ij}}
   * \f}
   * where \f$\alpha\f$ is the screening parameter, \f$N_{\mathrm{ions}}\f$ is the total number of ions (atoms),
   * \f$N_\mathrm{pairs}\f$ means the sum over all pairs of atoms from different molecules, \f$N_\mathrm{inmol}\f$ means summing over all intramolecular pairs.
   * For intramolecular electrostatic interactions, the corresponding term must be specified as an intramolecular field (see @ref AbstractIntraMolField).
   * The potential energy is added to the parent SimulatedSystem::Epotintermol[0] (electrostatic terms).
   *
   * For the sake of computation, it is convenient to split the computation in two parts.
   * In reciprocal space, we are not able to distinguish between atoms.
   * The energy contribution is
   * \f{equation}{
   *     E_{k\text{-}\mathrm{space}} = \frac{1}{2V\epsilon} \sum_{\mathbf{k}\neq \mathbf{0}}^\mathbf{\infty} \frac{\exp\left(-k^2 / (4\alpha^2)\right)}{k^2}\left|\sum_j^{N_\mathrm{ions}} q_j \exp\left(-\mathrm{i} \mathbf{k}\cdot\mathbf{r_j}\right) \right|^2
   * \f}
   * On the other hand, in real space, we can proceed as usually by summing up terms resulting from pair interactions.
   * \f{equation}{
   *     E_{r\text{-}\mathrm{space}} = \frac{1}{4\pi\epsilon}\sum_{i,j}^{N_\mathrm{pairs}}q_i q_j \frac{\mathrm{erfc}\left(\alpha r_{ij}\right)}{r_{ij}} - \frac{1}{4\pi\epsilon} \sum_{\mathrm{molecules}} \sum_{i,j}^{N_\mathrm{inmol}} q_i q_j \frac{\mathrm{erf}\left(\alpha r_{ij}\right)}{r_{ij}}
   * \f}
   * Note, that we must also compute intramolecular corrections yielding from the fact, that we are not able to distinguish between atoms in \f$k\f$-space.
   * Finally, we are missing the constant self-interaction correction
   * \f{equation}{
   *     E_{\mathrm{self-corr}} = - \frac{1}{4\pi\epsilon}\frac{\alpha}{\sqrt{\pi}} \sum_j^{N_\mathrm{ions}} q_j^2
   * \f}
   * which is calculated once in constructor and then only added to the result.
   *
   * Forces (3D vector) are given by the derivative of the potential energy.
   * In real space, the pair-wise approach can be used.
   * For each pair of atoms from different molecules
   * \f{equation}{
   *     \mathbf{f}_j^{r\text{-}\mathrm{space}} = -\mathbf{f}_i^{r\text{-}\mathrm{space}} =
   *     \frac{q_i q_j}{4\pi\epsilon} \frac{2\alpha r_{ij}\exp\left(-\alpha^2 r_{ij}^2\right) + \sqrt{\pi}\mathrm{erfc}\left(\alpha r_{ij}\right)}{\sqrt{\pi}r_{ij}^2} \frac{\mathbf{r}_{ij}}{r_{ij}}
   * \f}
   * and for each intramolecular pair (self-corrections does not contribute to forces)
   * \f{equation}{
   *     \mathbf{f}_j^{r\text{-}\mathrm{space}} = -\mathbf{f}_i^{r\text{-}\mathrm{space}} =
   *     \frac{q_i q_j}{4\pi\epsilon} \frac{2\alpha r_{ij}\exp\left(-\alpha^2 r_{ij}^2\right) - \sqrt{\pi}\mathrm{erf}\left(\alpha r_{ij}\right)}{\sqrt{\pi}r_{ij}^2} \frac{\mathbf{r}_{ij}}{r_{ij}}
   * \f}
   * In reciprocal space, the contribution is @cite MACSIMUS_manual (chacked against @cite Aguado2003-Ewald)
   * \f{equation}{
   *     \mathbf{f}_j^{k\text{-}\mathrm{space}} =
   *     \frac{q_j}{V\epsilon} \sum_{\mathbf{k}\neq\mathbf{0}}^{\mathbf{\infty}} \frac{\exp\left(-k^2/(4\alpha^2)\right)}{k^2} \mathfrak{Im}\left[\exp\left(\mathrm{i}\mathbf{k}\cdot\mathbf{r}_{j}\right) \sum_i^{N_\mathrm{ions}} q_i \exp\left(-\mathrm{i} \mathbf{k}\cdot\mathbf{r_i}\right)\right] \mathbf{k}
   * \f}
   *
   * For precise (not cut) electrostatics, the virial is exactly equal to the negative potential energy.
   * The agreement should be good for Ewald summation @cite Kolafa2008.
   * The expressions for virial can be derived from expressions for stress cited in @cite Aguado2003-Ewald but they are not used here.
   * Instead, the negative value of electrostatic energy itself is added to the parent @ref SimulatedSystem::virintermol[0].
   *
   * Returns the potential energy.
   *
   * @par See also
   * @ref CalculateEpot() const.
   *
   * @return Potential energy caused by electrostatic interactions in the @ref parent @ref SimulatedSystem (internally in program units).
   */
  double CalculateForces() override; // calculate forces in given simSystem
  /**
   * @brief Calculate potential energy caused by electrostatic interactions (intermolecular).
   * Similar to @ref CalculateForces(), but calculates only the potential energy and does not change anything neither in the @ref parent SimulatedSystem, nor in the involved Atoms.
   *
   * The potential energy caused by electrostatic interactions calculated by Ewald method reads
   * Eventually, the potential energy can be computed from @cite dl_poly4-manual
   * \f{equation}{
   *     E_\mathrm{Coul} = \frac{1}{2V\epsilon} \sum_{\mathbf{k}\neq \mathbf{0}}^\mathbf{\infty} \frac{\exp\left(-k^2 / (4\alpha^2)\right)}{k^2}\left|\sum_j^{N_\mathrm{ions}} q_j \exp\left(-\mathrm{i} \mathbf{k}\cdot\mathbf{r_j}\right) \right|^2 + \frac{1}{4\pi\epsilon}\sum_{i,j}^{N_\mathrm{pairs}}\frac{q_i q_j}{r_{ij}}\mathrm{erfc}\left(\alpha r_{ij}\right) - \frac{1}{4\pi\epsilon}\frac{\alpha}{\sqrt{\pi}} \sum_j^{N_\mathrm{ions}} q_j^2 - \frac{1}{4\pi\epsilon} \sum_{\mathrm{molecules}} \sum_{i,j}^{N_\mathrm{inmol}} q_i q_j \frac{\mathrm{erf}\left(\alpha r_{ij}\right)}{r_{ij}}
   * \f}
   * where \f$\alpha\f$ is the screening parameter, \f$N_{\mathrm{ions}}\f$ is the total number of ions (atoms),
   * \f$N_\mathrm{pairs}\f$ means the sum over all pairs of atoms from different molecules, \f$N_\mathrm{inmol}\f$ means summing over all intramolecular pairs.
   * For intramolecular electrostatic interactions, the corresponding term must be specified as an intramolecular field (see @ref AbstractIntraMolField).
   * The potential energy is added to the parent SimulatedSystem::Epotintermol[0] (electrostatic terms).
   *
   * For the sake of computation, it is convenient to split the computation in two parts.
   * In reciprocal space, we are not able to distinguish between atoms.
   * The energy contribution is
   * \f{equation}{
   *     E_{k\text{-}\mathrm{space}} = \frac{1}{2V\epsilon} \sum_{\mathbf{k}\neq \mathbf{0}}^\mathbf{\infty} \frac{\exp\left(-k^2 / (4\alpha^2)\right)}{k^2}\left|\sum_j^{N_\mathrm{ions}} q_j \exp\left(-\mathrm{i} \mathbf{k}\cdot\mathbf{r_j}\right) \right|^2
   * \f}
   * On the other hand, in real space, we can proceed as usually by summing up terms resulting from pair interactions.
   * \f{equation}{
   *     E_{r\text{-}\mathrm{space}} = \frac{1}{4\pi\epsilon}\sum_{i,j}^{N_\mathrm{pairs}}q_i q_j \frac{\mathrm{erfc}\left(\alpha r_{ij}\right)}{r_{ij}} - \frac{1}{4\pi\epsilon} \sum_{\mathrm{molecules}} \sum_{i,j}^{N_\mathrm{inmol}} q_i q_j \frac{\mathrm{erf}\left(\alpha r_{ij}\right)}{r_{ij}}
   * \f}
   * Note, that we must also compute intramolecular corrections yielding from the fact, that we are not able to distinguish between atoms in \f$k\f$-space.
   * Finally, we are missing the constant self-interaction correction
   * \f{equation}{
   *     E_{\mathrm{self-corr}} = - \frac{1}{4\pi\epsilon}\frac{\alpha}{\sqrt{\pi}} \sum_j^{N_\mathrm{ions}} q_j^2
   * \f}
   * which is calculated once in constructor and then only added to the result.
   *
   * For intramolecular electrostatic interactions, the corresponding term must be specified as an intramolecular field (see @ref AbstractIntraMolField).
   * In case of free periodic boundary conditions, the energy is the same as if \f$\mathtt{shift}\rightarrow0\f$ and \f$\mathtt{cutoff}\rightarrow\infty\f$.
   *
   * @par See also
   * @ref CalculateForces().
   *
   * @return Potential energy caused by electrostatic interactions in the @ref parent @ref SimulatedSystem (internally in program units).
   */
  double CalculateEpot() const override;
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
   * @ref EwaldSplinesElstat(const EwaldSplinesElstat&), @ref SetParent(SimulatedSystem*), @ref AbstractInterMolField::copy(SimulatedSystem *) const, @ref AbstractIntraMolField::copy(Molecule *) const.
   *
   * @param newparent Parent SimulatedSystem of the newly created EwaldSplinesElstat object.
   * @return (Pointer to) the newly created EwaldSplinesElstat object.
   */
  EwaldSplinesElstat *copy(SimulatedSystem *newparent) const override; // cloning itself
  /**
   * @brief Ewald corrections arising from cutting the electrostatic potential are assumed to be zero.
   * Generally we want to avoid consequences of cutting the potential after certain distance.
   * Ewald corrections are usually given by
   * \f{equation}{
   *      E_\mathrm{corr} = 2\pi \sum_{i,j} \frac {1}{V} \int_{0}^{\infty}\! \left(E_{\mathrm{id}} - E_{\mathrm{sim}}\right) g(r) r^2\,\mathrm{d}r
   * \f}
   * where \f$E_{\mathrm{id}}\f$ is the original (desired) potential (in our case original Coulomb potential), \f$E_{\mathrm{sim}}\f$ is the truncated potential used in the simulation and \f$g(r)\f$ is the radial distribution function and the sum runs over all charged atoms in the system.
   * Essentially, we must integrate the difference between the real and desired potential for all pairs in the system.
   *
   * In free (vacuum) boundary conditions, the interactions can be calculated explicitely and no cutoff is applied.
   * Ewald summation cannot be used.
   *
   * In periodic boundary conditions the \f$k\f$-space sum can be viewed as an alternative to cutoff corrections.
   * The *normal* cutoff corrections cannot be computed as the integral diverges for slowly decaying electrostatic interactions.
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
   * @brief Print informations about the electrostatics calculated by Ewald summation to the @p stream
   *
   * @param stream Output stream.
   * @param u_eng Energy unit value (see `units.hpp`).
   * @param engunit Energy unit name in brackets (as in `.cpa` file).
   */
  void PrintInfo(std::ofstream &stream, double u_eng, std::string engunit) const override;
};

#endif

/**********************END OF CODE*********************************************
## Notes on parallelization
- `sins, coss`
  - atoms * kspace
- `expk`
  - reduced kspace
- `imag, real`
  - kspace
- `Eel, Erinter, Erintra, Ek`
  - cumulators
- `i, j, k, l, m, n, w, w0, w00`
  - loop indeces
- `r[3]`, `force3[3]`
  - aux. arrays
  - can be made private by omp since their size is known
- `box[3], rbox[3]`
  - *constant* arrays set once and read-only
- `force, q2, alphar, erfcrr, erfrr, A, B, C, D, mxmult, charge, pi2rx_L, ..., cosx, sinx, ...`
  - aux. numbers
- `ksq`
- `atom1, atom2`
  - aux. pointers
- `pair`
  - aux. object
- `no_ch_atoms`
  - class member (not changing)

#### rspace inter
int m; // redeclared (needed prior)
double r2, rnorm;
double force;
double q2;
double alphar; // alpha*r
double erfcrr; // erfc(alpha*r)/r
Atom *atom1, *atom2;
Pair pair;

#### rspace intra (loop through atoms)
int m; // redeclared
int i, k, l;
double q2; // redeclared
double r[3];
double r2, rnorm; // redeclared
double force; // redeclared
double alphar; // alpha*r // redeclared
double erfrr; // erf(alpha*r)/r

#### kspace sums (loop through atoms)
int m; // redeclared
double charge;
double pi2rx_L,...;
double cosx, sinx, ...;
int i, k; // redeclared
int j; // redeclared (needed prior)
int mxmult;
int w, w0, w00;
* double sins, coss; // need not be private (uses different locations for different atoms)
* double imag, real; // **shared!!!** – must be made threadprivate (allocate more space)
double A, B, C, D;

#### kspace forces (loop through atoms)
int i, j, k, l, m, mxmult; // redeclared
int w; // redeclared
double force3[3];
double A, B, C, D; // redeclared

*/