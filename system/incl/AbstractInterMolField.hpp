#ifndef ABSTRACTINTERMOLHEADER
#define ABSTRACTINTERMOLHEADER

class SimulatedSystem;

/**
 * @defgroup InterMolField Intermolecular interactions
 * @brief Set of classes describing interactions between molecules.
 * Intermolecular interactions can be devided into two groups.
 *
 * **Electrostatic interactions** are caused by electric **charges** (full or partial) on Atoms.
 * They can be handled in several ways, the simplest being the cutoff approach.
 * More rigorous approaches would involve either Ewald summation or a reaction field.
 * The atomic charges are defined in `.field` file and the computation method is specified in `.control` file.
 *
 * **Disperse (or van der Waals) interactions** originate from quantum interactions between electronic clouds.
 * They usually have two parts.
 * Strong repulsive forces prevail over short distances, long separations are dominated by soft attractions.
 *
 * @class AbstractInterMolField
 * @ingroup InterMolField
 * @brief Abstract class (interface) for potentials (and forces) between Molecules.
 *
 * This abstract class provides an interface for intermolecular interactions.
 * The following constituents are expected to implement this interface:
 * - electrostatic interactions (cutoff, Ewald sum, reaction field,...) – type 0
 * - disperse interactions (Lennard-Jones, Morse,...) – type 1
 *
 * Two main functions are prescribed, that both return potential energy due to the constituent.
 * @ref CalculateForces() also modifies the forces for each @ref Atom in the SimulatedSystem.
 * 
 * Parameters of disperse interactions are specified in `.field` file (see particular interaction type). 
 * Typically, only parameters for homogeneous pairs (both atoms of the same type) are specified 
 * and cross-terms can be calculated automatically using some `vdw mixing` rule in `.control` file (default no mixing).
 * Automatically calculated values can be overriden by specifying the pair manually.
 * This way, the interaction can be supressed for the particular pair by specifying non-zero diameter and zero energy.
 * Currently, the following mixing rules are implemented:
 * - **no mixing** – default, no cross-terms
 * - Lorentz–Berthelot – geometric mean of energy and arithmetic mean of diameters
 * ```(control)
 * vdwmixing lorentz-berthelot
 * ```
 */
class AbstractInterMolField
{
private:
    /**
     * @brief Default constructor of a new AbstractInterMolField object
     *
     * Default constructor made private in order not instantiate empty object.
     */
    AbstractInterMolField() : type(0), parent(nullptr){};

protected:
    /**
     * @brief Type of intermolecular field
     *
     * The following types are allowed now:
     * - electrostatic (type = 0)
     * - disperse (type > 0) – @ref LJsystem(1),..
     */
    int type;
    /**
     * @brief Parent SimulatedSystem
     * The pointer to the System instance this intramolecular field term belongs to.
     * Indicates, where the involved Molecules are stored.
     */
    SimulatedSystem *parent;

public:
    /**
     * @brief Cutoff correction term caused by truncating this field after some cutoff distance.
     * Cutoff correction energy is calculated by dividing this term by volume.
     * Another division by volume yields correction for pressure. 
     */
    double Ecutoffcorr;
    /**
     * @brief Interaction cutoff
     * Maximum distance of interaction. 
     * Beyond the cutoff, there is no interaction between particles.
     */
    double cutoff;
     /**
     * @brief Interaction cutoff squared
     * See @ref cutoff.
     */
    double cutoff2;
    /**
     * @brief Constructor of an AbstractInterMolField object specifying the type and parent SimulatedSystem.
     *
     * @param tp Type of interaction (see @ref type)
     * @param par Parent system (see @ref parent)
     */
    AbstractInterMolField(int tp, SimulatedSystem *par) : type(tp), parent(par){};
    /**
     * @brief Copy constructor of an AbstractInterMolField object
     *
     * @param otherField Other AbstractInterMolField object that should be copied.
     */
    AbstractInterMolField(const AbstractInterMolField &otherField) : type(otherField.type), parent(nullptr){};
    /**
     * @brief Default destructor for AbstractInterMolField
     * This destructor is empty, because no dynamic memory is allocated in constructors.
     * Can be overriden in classes that implement this interface (inherit this pure virtual class).
     */
    virtual ~AbstractInterMolField(){};
    /**
     * @brief Overloading assignment operator=
     * This function assigns values from other AbstractInterMolField object to the left-hand-side object.
     * Eventually not needed in the implementation.
     *
     * @param otherField Object whose properties are copied to the LHS object.
     * @return The LHS object.
     */
    AbstractInterMolField &operator=(const AbstractInterMolField &otherField)
    {
        parent = nullptr;
        type = otherField.type;
        return *this;
    };
    /**
     * @brief Set the @ref parent SimulatedSystem
     * New objects are created as orphans (the pointer @ref parent is `null`). This function is used to set the pointer @ref parent to the new parent SimulatedSystem ( @p par).
     *
     * @param par (Pointer to the) new parent Molecule.
     * @return *true* Setting was successfull (it means that this instance has been an orphan before).
     * @return *false* This instance has already have a valid parent SimulatedSystem and the new parent is not set.
     */
    bool SetParent(SimulatedSystem *par)
    {
        if (parent == nullptr)
        {
            parent = par;
            return true;
        }
        else
            return false;
    };
    /**
     * @brief Calculate forces caused by this intermolecular field.
     * Calculate forces caused by this field and add them to @ref Atom::force.
     * Save the potential energy and force virial to the @ref parent SimulatedSystem @ref SimulatedSystem::Epotintermol and @ref SimulatedSystem::virintermol.
     * Return the potential energy.
     *
     * @par See also
     * @ref CalculateEpot() const.
     *
     * @return Potential energy caused by this intermolecular term (internally in program units).
     */
    virtual double CalculateForces() = 0;
    /**
     * @brief Calculate potential energy caused by this intermolecular field.
     * Similar to @ref CalculateForces(), but calculates only the potential energy and does not change anything neither in the @ref parent SimulatedSystem.
     * Particularly, it does not affect Atom::force.
     *
     * @par See also
     * @ref CalculateForces().
     *
     * @return Potential energy caused by this intermolecular term.
     */
    virtual double CalculateEpot() const = 0;
    /**
     * @brief Clone itself and get adopted by another @ref parent SimulatedSystem.
     * Returns a copy of itself and migrate to the @p newparent SimulatedSystem.
     * Useful when copying SimulatedSystem object.
     * Every intermolecular term is then cloned by this method and assigned to the new SimulatedSystem object.
     *
     * @par See also
     * @ref AbstractInterMolField(const AbstractInterMolField&), @ref SetParent(SimulatedSystem*), @ref AbstractIntraMolField::copy(Molecule *) const.
     *
     * @param newparent Parent SimulatedSystem of the newly created AbstractInterMolField object.
     * @return (Pointer to) the newly created AbstractInterMolField object.
     */
    virtual AbstractInterMolField *copy(SimulatedSystem *newparent) const = 0;
    /**
     * @brief Set the force-field parameters for a given pair
     * All intermolecular terms (except for electrostatic) are set up in `vdw` part of the `.field` file.
     * The common part is the pair specification (by Atom names). 
     * Names are converted to LJ IDs (see @ref Atom::LJid) (integer numbers) inside FieldFile class.
     * Other parameters are to be read and understood by the particular intermolecular field.
     * 
     * @param id1 ID of the first Atom (unique number corresponding to @ref Atom::name)
     * @param id2 ID of the second Atom (unique number corresponding to @ref Atom::name
     * @param params Line of pair potential parameters (starting by some kind of energy term).
     * @return  0 if success
     * @return 28 energy parameter (epsilon,...) incorrect
     * @return 29 other parameter (sigma,...) incorrect
     */
    virtual int SetPairParams(int id1, int id2, char *params) {return 0;};
    /**
     * @brief Get the diameter for atoms with given @ref Atom::LJid.
     * Disperse interactions are used to estimate the vdw radius of Atoms. 
     * Thus, having defined the self pair interaction, the diameter can be estimated. 
     * This estimate is used only for visualisation (for `.mol` file, see @ref FieldFile::PrintMol(char*, SimulatedSystem*))
     * 
     * @param id1 ID of the Atom (unique number corresponding to @ref Atom::name)
     * @return Atom diameter estimate rounded to integer number
     */
    virtual int GetDiameter(int id1) const {return 0;};
    /**
     * @brief Set the cutoff distance.
     * Set the distance beyond which the interaction is zero.
     * Applies only in periodic boundary conditions.
     * Variable @p cutoff corresponds to vdW interactions cutoff, @p elcutoff to electrostatics real cutoff.
     * 
     * @param cutoff Cutoff distance for vdW (disperse) interactions (in [AA])
     * @param elcutoff Cutoff distance for electrostatics (in [AA])
     * @param boundaryC Boundary conditions (see @ref SimulatedSystem::boundaryCond)
     * @return 0 if success, number of pairs leading to error (for LJ system too short cutoff)
     */
    virtual int SetCutoff(double cutoff, double elcutoff, int boundaryC) = 0;
    /**
     * @brief Calculate mixed terms for disperse interactions between atoms of different types
     * Calculate potential parameters for unspecified atomic pairs according to mixing rule @p mixingrule.
     * Possible mixing rules are:
     * - **no mixing** ( @p mixingrule = 0, default) – no implicit mixing for unspecified pairs
     * - Lorentz–Berthelot ( @p mixingrule = 1) – geometric mean of energy and arithmetic mean of diameters
     * 
     * @param mixingrule Value of mixing rule defined in `.control` file
     * @return Number of newly specified mixed terms. 
     */
    virtual int CalculateMixTerms(int mixingrule) {return 0;};
    /**
     * @brief Calculate the correction to energy and pressure resulting from omitting the interactions beyond @ref cutoff
     * Cutoff corrections are calculated under the assumption, that the radial distribution function is one beyond the @ref cutoff distance.
     * The mean *density of pairs* can than be evaluated and the *neglected interaction energy* calculated integrating over the neglected potential from @ref cutoff to infinity.
     * 
     * @par See also
     * @ref Ecutoffcorr, @ref SimulatedSystem::PrintMeasurement(double, double)
     * 
     * @return Correction *energy*, divide by volume to get correction to potential energy and by volume squared to get correction to pressure
     */
    virtual double CalculateCutoffCorrection() = 0; // returns energy (times volume)
    /**
     * @brief Initialize fake *bonds* used to handle this intermolecular field between atoms of one Molecule.
     * 
     * Electrostatic and dispersion interactions are basically intermolecular, but they occur also in longer molecules.
     * Atoms on distances 1–4 interact by scaled versions of the original intermolecular interactions,
     * atoms which are furhter from each other interact as if they belonged to different molecules.
     * Each intermolecular field must know, what kind of bond can replace it as an intramolecular field.
     * 
     * @param LJ14factor Scaling factor for 1–4 disperse interactions.
     * @param el14factor Scaling factor for 1–4 electrostatic interactions.
     * @return Number of newly initialized bonds.
     */
    virtual int InitializeIntramolBonds(double LJ14factor = 1.0, double el14factor = 1.0) const = 0; // initialize fake *bonds* used to handle this intermolecular field between atoms in one molecule
    /**
     * @brief Get the field type.
     * Provides access to the private member @ref type.
     * 
     * @return Type of the intermolecular field.
     */
    int GetType() const { return type; }
    /**
     * @brief Print informations about the force field to the @p stream
     * 
     * @param stream Output stream.
     * @param u_eng Energy unit value (see `units.hpp`).
     * @param engunit Energy unit name in brackets (as in `.cpa` file).
     */
    virtual void PrintInfo(std::ofstream &stream, double u_eng, std::string engunit) const = 0;
};

#endif