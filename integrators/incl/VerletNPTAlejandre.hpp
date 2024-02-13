#ifndef VERLETNPTALEJANDREHEADER
#define VERLETNPTALEJANDREHEADER

/**
 * @class VerletNPTAlejandre
 * @brief Implementation of Liouville-operator derived measure-preserving integrator for MTK equations of motion generating NPT ensemble @cite Tuckerman2006.
 *
 * Implementation of Alejandre's NPT integrator according to his original article @cite Tuckerman2006.
 * The integrator integrates the MTK equations of motion @cite MTTK1996 (the same as implemented in @ref VerletNPTAlejandre)
 * \f{align}{
 *     \dot{\mathbf{r}}_i(t) &= \mathbf{v}_i(t)+\dot{\lambda}(t)\mathbf{r}_i(t) \\
 *     \dot{\mathbf{v}}_i(t)  &=  \frac{\mathbf{f}_i(t)}{m_i}-\left[\dot{\xi}(t) + \left(1 + \frac{3}{N_\mathrm{f}}\right)\dot{\lambda}(t)\right]\mathbf{v}_i(t) \\
 *     \ddot{\xi} &= \frac{1}{M_T}\left(2E_\mathrm{kin} + M_P\dot{\lambda}(t)^2-(N_\mathrm{r}+1)k_\mathrm{B}T_\mathrm{ext}\right)\\
 *     \ddot{\lambda}(t) &= \frac{3}{M_P}\left[V(P(t)-P_\mathrm{ext})+\frac{2E_\mathrm{kin}}{N_\mathrm{f}}\right]-\dot\xi(t)\dot\lambda(t) \\
 *     V &= L_x L_y L_z \\
 *     L_k &= L_{k,0} \exp⁡(\lambda) \quad \text{where } k=x,y,z \\
 *     E_\mathrm{tot} &= E_\mathrm{kin}+\frac{M_T}{2}\dot{\xi}(t)^2+\frac{M_P}{2}\dot{\lambda}(t)^2+E_\mathrm{pot}+P_\mathrm{ext}V(t)+(N_\mathrm{f}+1)k_\mathrm{B}T_\mathrm{ext}\xi \\
 * \f}
 * where \f$\mathbf{r}_i\f$, \f$\mathbf{v}_i\f$ and \f$\mathbf{f}_i\f$ are vectors of positions and velocities of atom \f$i\f$ and forces on atom \f$i\f$,
 * \f$\lambda\f$ is the extra degree of freedom coupled with barostat (it holds \f$\lambda=\ln(L_k/L_{k,0})\f$),
 * \f$\xi\f$ is the extra degree of freedom connected with thermostat (here we use the same thermostat both for the real degrees of freedom and for \f$\lambda\f$, though originally two thermostats are involved in @cite Tuckerman2006),
 * \f$M_T\f$ is the *mass of the thermostat* and \f$M_P\f$ is the *mass of the piston of the barostat*,
 * \f$T_\mathrm{ext}\f$ and \f$P_\mathrm{ext}\f$ are final (desired) temperature and pressure,
 * \f$L_k\f$ denotes one side of the simulation box (\f$k=x,y,z\f$), \f$V\f$ is the volume of the box,
 * \f$E_\mathrm{kin}\f$ is the kinetic energy of the system (real degrees of freedom only), number of degrees of freedom is called \f$N_\mathrm{f}\f$,
 * \f$E_\mathrm{pot}\f$ is the potential energy of the system and \f$E_\mathrm{tot}\f$ is the total (conserved) energy including the extra degrees of freedom.
 * Currently, only isotopic changes of the simulation box are implemented.
 *
 * The *masses* of the thermostat and barostat are given by
 * \f{align}{
 *     M_T&=N_\mathrm{f}k_\mathrm{B}T_\mathrm{ext}\tau_T^2 \\
 *     M_P&=(N_\mathrm{f}+3)k_\mathrm{B}T_\mathrm{ext}\tau_P^2 \\
 * \f}
 * where \f$\tau_T\f$ and \f$\tau_P\f$ are the relaxation (correlation, dumping) times of the thermostat and barostat, respectively.
 *
 * Integration of these equations of motion is done in a symplectic, time-reversible way, usign the algorithm devised in @cite Tuckerman2006. The integration steps are as follows
 * 1. thermostat half step (\f$\mathrm{i} \mathcal{L}_T(h/2)\f$ according to the article notation)
 * 2. barostat half step (\f$\mathrm{i} \mathcal{L}_{\epsilon,2}(h/2)\f$)
 * 3. Verlet half step (velocities) (\f$\mathrm{i}\mathcal{L}_2(h/2)\f$)
 * 4. barostat full step (\f$\mathrm{i} \mathcal{L}_{\epsilon,1}(h)\f$)
 * 5. Verlet full step (positions) (\f$\mathrm{i} \mathcal{L}_1(h)\f$)
 * 6. SHAKE
 * 7. forces calculation
 * 8. Verlet half step (velocities) (\f$\mathrm{i}\mathcal{L}_2(h/2)\f$)
 * 9. RATTLE
 * 10. barostat half step (\f$\mathrm{i} \mathcal{L}_{\epsilon,2}(h/2)\f$)
 * 11. thermostat half step (\f$\mathrm{i} \mathcal{L}_T(h/2)\f$)
 *
 * All steps are described in appropriate functions.
 * This algorithm is different from the original integration devised in @cite MTTK1996 (and used in @ref VerletNPTNose) and does not comprise quarter steps.
 * It should perform comparably or better.
 *
 */

#include "AbstractIntegrator.hpp"
#include "AbstractVerletIntegrator.hpp"
class Matrix;
class Vector;

class VerletNPTAlejandre : public AbstractVerletIntegrator
{
private:
    /**
     * @brief Final (external, desired) temperature \f$T_\mathrm{ext}\f$
     * Temperature of the thermostat, the *heat bath*.
     */
    double T_f;
    /**
     * @brief Characteristic time of thermostat \f$\tau_T\f$ in ps
     * The shorter the time, the faster the convergence.
     * However, too short times can cause undesired side-effects.
     */
    double tau_T;
    /**
     * @brief Final (external, desired) pressure \f$P_\mathrm{ext}\f$
     * External pressure that should be maintained in the system.
     */
    double P_f;
    /**
     * @brief Characteristic time of barostat \f$\tau_P\f$ in ps
     * Should be several times longer than \f$\tau_T\f$.
     */
    double tau_P;
    /**
     * @brief The mass of the thermostat
     * Thermostat mass is given by \f$M_T = N_\mathrm{f} k_\mathrm{B} T_\mathrm{ext} \tau_T^2\f$.
     */
    double M_T;
    /**
     * @brief The mass of the barostat
     * Barostat mass is given by \f$M_P = (N_\mathrm{f} + 3) k_\mathrm{B} T_\mathrm{ext} \tau_P^2\f$.
     */
    double M_P;
    /**
     * @brief Extra degree of freedom connected with the thermostat
     * Derivatives are saved in the same manner as in @ref Atom::R.
     */
    Vector *Xi;
    /**
     * @brief Extra degree of freedom connected with the barostat
     * Derivatives are saved in the same manner as in @ref Atom::R.
     * This degree of freedom is connected with the box size by \f$\lambda = \ln{L_k/L_{k,0}}\f$.
     */
    Vector *Lambda;
    /**
     * @brief Extra degree of freedom connected with the thermostat – old value
     * Derivatives are saved in the same manner as in @ref Atom::R.
     */
    Vector *Xi_old;
    /**
     * @brief Extra degree of freedom connected with the barostat – old value
     * Derivatives are saved in the same manner as in @ref Atom::R.
     */
    Vector *Lambda_old;
    /**
     * @brief Old box size for SHAKE
     * Stores the value of the box size (\f$x,y,z\f$) for the use in SHAKE.
     */
    double orig_size[3]; // original box_sizes
    /**
     * @brief Half-step of the thermostat connected Liouville operator (\f$\mathrm{i}\mathcal{L}_T(h/2)\f$)
     * The thermostat half-step proceeds as follows.
     * 1. update of \f$\xi\f$
     *     \f{equation}{ \xi(t + h/8) = \xi(t) + \frac{h}{8} \dot\xi(t)\f}
     * 2. *force* on \f$\xi\f$
     *     \f{equation}{\dot\xi(t+h/4) = \dot\xi(t) + \frac{h}{4} \frac{2E_\mathrm{kin}(t) + M_P \dot\lambda(t)^2 - (N_\mathrm{f} + 1) k_\mathrm{B} T_\mathrm{ext}}{M_T}\f}
     * 3. scaling of velocities:
     *     \f{align}{
     *         \mathbf{v}_i(t) &= \mathbf{v}_i(t) \exp \left(-\frac{h}{2} \dot\xi(t+h/4)\right) \\
     *         \lambda(t) &=      \lambda(t) \exp \left(-\frac{h}{2} \dot\xi(t+h/4)\right) \\
     *     \f}
     * 4. update of \f$\xi\f$
     *     \f{equation}{ \xi(t + 3h/8) = \xi(t + h/8) + \frac{h}{4} \dot\xi(t + h/4)\f}
     * 5. *force* on \f$\xi\f$ (kinetic energy has changed)
     *     \f{equation}{\dot\xi(t+h/2) = \dot\xi(t+h/4) + \frac{h}{4} \frac{2E_\mathrm{kin}(t) + M_P \dot\lambda(t)^2 - (N_\mathrm{f} + 1) k_\mathrm{B} T_\mathrm{ext}}{M_T}\f}
     * 6. update of \f$\xi\f$
     *     \f{equation}{\xi(t + h/2) = \xi(t+3h/8) + \frac{h}{8} \dot\xi(t+h/2)\f}
     *
     * @param Ekin System kinetic energy (calculated before and changing inside)
     * @param velscale Velocity scaling factor (to accumulate velocity scaling and done it only when necessary)
     * @return 0 on success
     */
    int ThermostatHalfStep(double &Ekin, double &velscale);
    /**
     * @brief Integration step for barostat connected variable \f$\lambda\f$, \f$\mathrm{i}\mathcal{L}_{\epsilon,1}(h)\f$
     * The barostat integration step is
     * \f{equation}{
     *     \lambda(t+h) = \lambda(t) + h \dot\lambda(t+h/2)
     * \f}
     * Box scaling by \f$\exp (\dot\lambda(t+h/2))\f$ is performed after this step.
     *
     * @return 0 on success
     */
    int BarostatFullStep1();
    /**
     * @brief Half-step integratrion of barostat connected velocity \f$\dot\lambda\f$, \f$\mathrm{i}\mathcal{L}_{\epsilon,2}(h/2)\f$
     * Action of barostat is defined in this step
     * \f{equation}{
     *     \dot\lambda(t+h/2) = \dot\lambda(t) + \frac{h}{2}\frac{3\left[V (P(t) - P_\mathrm{ext}) + 2E_\mathrm{kin}/N_\mathrm{f}\right]}{M_P}
     * \f}
     * Instanteneous pressure \f$P(t)\f$ is calculated inside.
     *
     * @param Ekin Kinetic energy of the system
     * @return 0 on success
     */
    int BarostatHalfStep2(double Ekin);
    /**
     * @brief Verlet step (update of positions), \f$\mathrm{i}\mathcal{L}_1(h)\f$
     * Verlet step is modified as follows
     * \f{equation}{
     *     \mathbf{r}_i(t+h) = \mathbf{r}(t) + h \mathbf{v}_i(t+h/2) \exp(\alpha)\frac{\sinh(\alpha)}{\alpha}
     * \f}
     * where \f$\alpha = h\dot\lambda(t+h/2)/2\f$.
     *
     * @param p_mol Pointer to the molecule whose atom positions should be updated
     * @return 0 on success
     */
    int VerletStep(Molecule *p_mol) override;
    /**
     * @brief Verlet half-step (update of velocities), \f$\mathrm{i}\mathcal{L}_2(h/2)\f$
     * Atom velocities are updated in the following way
     * \f{equation}{
     *     \mathbf{v}_i(t + h/2) = \mathbf{v}_i(t)\exp \left[-\frac{h}{2}(1+3/N_\mathrm{f})\dot\lambda(t+h/2)\right] +
     *     \frac{h}{2}\frac{\mathbf{f}_i(t)}{m_i}\exp(-\alpha)\frac{\sinh(\alpha)}{\alpha}
     * \f}
     * where \f$\alpha = \frac{h}{4}\left(1 + \frac{3}{N_\mathrm{f}}\right)\dot\lambda(t+h/2)\f$ and
     * \f$\dot\lambda(t+h/2)\f$ is still the same (even for the second application of this step).
     *
     * @param p_mol Pointer to the molecule whose atom velocities should be updated
     * @return 0 on success
     */
    int VerletStep2(Molecule *p_mol);
    /**
     * @brief Time shift \f$t+h \rightarrow t\f$.
     * Needed both for the Atom::R and for the extended degrees of freedom.
     *
     * @return 0 on success
     */
    int TimeShift() override;
    /**
     * @brief Initialization of the extended variable \f$\xi\f$
     *
     * @return 0 on success
     */
    int InitXi();
    /**
     * @brief Initialization of the extended variable \f$\lambda\f$
     *
     * @return 0 on success
     */
    int InitLambda();
    /**
     * @brief Shake implementation for this type of integrator.
     * Constrained bonds are corrected to have the right length.
     * The same implementation as for MTTK version (see @ref VerletNPTNose).
     *
     * @param p_mol Pointer to the molecule which should be *shaked*.
     * @param virial Virial of constraint forces to accumulate the contributions from all molecules.
     * @return The highest relative deviation (relative error) of constrained bonds.
     */
    double Shake(Molecule *p_mol, double &virial) override;
    /**
     * @brief Rattle implementation for this kind of integrator.
     * Velocities are made perpendicular to constrained bonds.
     * The same implementation as for MTTK version (see @ref VerletNPTNose).
     *
     * @param p_mol Pointer to the molecule which should be *rattled*.
     * @param virial Virial of constraint forces to accumulate the contributions from all molecules.
     * @return The highest error of the angle between velocity and corresponding constrained bond.
     */
    double Rattle(Molecule *p_mol, double &virial);

public:
    /**
     * @brief Construct a new integrator for NPT Nosé/MTK ensemble.
     * New integrator is created that integrates NPT MTK dynamics by Trotter decomposition introduced by Alejandre (@cite Tuckerman2006).
     *
     * @param step Integration step in ps.
     * @param simsystem Simulated system of atoms.
     * @param Tfinal Temperature of the thermostat (\f$T_\mathrm{ext}\f$)
     * @param tauT Thermostat time constant (\f$\tau_T\f$)
     * @param Pfinal Pressure of the barostat (\f$P_\mathrm{ext}\f$)
     * @param tauP Barostat time constant (\f$\tau_P\f$)
     * @param startWithShake Whether to start by SHAKE or not
     */
    VerletNPTAlejandre(double step, SimulatedSystem *simsystem, double Tfinal, double tauT, double Pfinal, double tauP, bool startWithShake);
    /**
     * @brief Destroy the VerletNPTAlejandre integrator.
     * Delete the extended degrees of freedom variables @ref Xi and @ref Lambda.
     */
    ~VerletNPTAlejandre();
    /**
     * @brief Perform @p Nsteps integration steps.
     * The crucial method of the class. 
     * Performs one integration cycle and save integrator-connected variables to @ref system. 
     * 
     * @param Nsteps Number of integration steps to be performed
     * @return 0 on success
     */
    int Integrate(int Nsteps) override;
    /**
     * @brief Returns the energy of the extra degrees of freedom.
     * The energy of the extra DOFs is
     * \f{equation}{
     *     E_\mathrm{extra} = \frac{M_T}{2}\dot{\xi}(t)^2+\frac{M_P}{2}\dot{\lambda}(t)^2+P_\mathrm{ext}V(t)+(N_\mathrm{f}+1)k_\mathrm{B}T_\mathrm{ext}\xi
     * \f}
     * This value is added to the kinetic and potential energy of the @ref system to obtain the conserved energy.
     * 
     * @return The energy of the extra degrees of freedom
     */
    double EnergyOfExtraDOF() const override;
    /**
     * @brief Prepare values in @ref Atom::R for `.config` printing.
     * Values must be ordered in a unified manner. 
     * Dependence on the integration step is removed.
     * 
     * @see @ref SimulatedSystem::PrintConfig(). 
     * 
     * @return `configlevel` (or `levcfg`) – level of `.config` information
     */
    int PrepareConfig() override;
    /**
     * @brief Initialize values in @ref Atom::R for integration.
     * Values are loaded in a unified manner from `.config` file.
     * This function reorganizes them and scale them by the integration step.
     * 
     * @param startWithShake Indicates whether to perform SHAKE during initialization or not
     * @return 0 on success
     */
    int Initialize(bool startWithShake) override;
};

#endif