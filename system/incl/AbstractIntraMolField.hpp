#ifndef ABSTRACTINTRAMOLHEADER
#define ABSTRACTINTRAMOLHEADER

class Molecule;

/**
 * @defgroup IntraMolField Intramolecular interactions
 * @brief Set of classes forming interactions within a Molecule.
 * Intramolecular interactions comprises bonds ( @ref HarmonicBond), angles ( @ref HarmonicAngle), dihedrals ( @ref CosineDihedral) and other interactions within a Molecule. 
 * Moreover disperse (Lennard-Jones) and electrostatic interactions within a Molecule are realized by classes from this set ( @ref LJBond).
 * They all implement @ref AbstractIntraMolField (as an *interface*).
 * They are defined in `.field` file and they perform forces calculation for the Molecule.
 *
 * @class AbstractIntraMolField
 * @ingroup IntraMolField
 * @brief Abstract class (interface) for potentials (and forces) constituting a Molecule (holding Atoms together).
 *
 * This abstract class provides an interface for all forces within a Molecule (except for constrained bonds ( @ref ConstraintBond),
 * as they are essentially different.
 * The following constituents are expected to implement this interface:
 * - bonds (harmonic ( @ref HarmonicBond), Morse,...) – type 0
 * - angles (harmonic ( @ref HarmonicAngle),...) – type 1
 * - dihedrals (cosine ( @ref CosineDihedral),...) – type 2
 * - improper torsions (cosine ( @ref CosineDihedral),...) – type 5
 * - global pair forces (LJ interactions (type 3, @ref LJBond), electrostatic interactions (type 4),...)
 * Global pair forces are assumed to be similar to intermolecular forces (see @ref AbstractInterMolField), but they can ommit neighbours and 1–3 interactions
 * and modify 1–4 interactions. In practice, they will behave as a (set of) bonds.
 * 
 * Improper torsions share implementation with dihedrals, the only difference is that
 * the involved atoms need not to be directly connected by bonds or constrained bonds.
 *
 * Two main functions are prescribed, that both return potential energy due to the constituent.
 * @ref CalculateForces() also modifies the forces for each @ref Atom in the parent Molecule.
 */
class AbstractIntraMolField
{
private:
    /**
     * @brief Default constructor of a new AbstractIntraMolField object
     *
     * Default constructor made private in order not instantiate empty object.
     */
    AbstractIntraMolField() : type(0), parent(NULL){};

protected:
    /**
     * @brief Type of intramolecular field
     *
     * The following types are allowed now:
     * - bond (type = 0) – @ref HarmonicBond,...
     * - angle (type = 1) – @ref HarmonicAngle,..
     * - dihedral (type = 2) – @ref CosineDihedral,...
     * - improper torsion (type = 5) – @ref CosineDihedral,...
     * - global vdw (type = 3) – 1–4, 1–5,... interactions via vdw forces, @ref LJBond
     * - electrostatic (type = 4) – as previous, but forces caused by charges; not implemented yet
     */
    int type; // bond(0), angle(1), dihedral(2), global (LJ,...) (3), electrostatic(4),...
    /**
     * @brief Parent molecule
     * The pointer to the Molecule instance this intramolecular field term belongs to.
     * Indicates, where the involved atoms are stored.
     */
    Molecule *parent; // parent molecule (to know where to search for atoms)

public:
    /**
     * @brief Constructor of an AbstractIntraMolField object specifying the type and parent Molecule.
     *
     * @param tp Type of interaction (see @ref type)
     * @param par Parent molecule (see @ref parent)
     */
    AbstractIntraMolField(int tp, Molecule *par) : type(tp), parent(par){};
    /**
     * @brief Copy constructor of an AbstractIntraMolField object
     *
     * @param otherField Other AbstractIntraMolField object that should be copied.
     */
    AbstractIntraMolField(const AbstractIntraMolField &otherField) : type(otherField.type), parent(NULL){};
    /**
     * @brief Default destructor for AbstractIntraMolField
     * This destructor is empty, because no dynamic memory is allocated in constructors.
     * Can be overriden in classes that implement this interface (inherit this pure virtual class).
     */
    virtual ~AbstractIntraMolField(){};
    /**
     * @brief Overloading assignment operator=
     * This function assigns values from other AbstractIntraMolField object to the left-hand-side object.
     * Eventually not needed in the implementation.
     *
     * @param otherField Object whose properties are copied to the LHS object.
     * @return The LHS object.
     */
    AbstractIntraMolField &operator=(const AbstractIntraMolField &otherField)
    {
        parent = NULL;
        type = otherField.type;
        return *this;
    };
    /**
     * @brief Set the @ref parent Molecule
     * New objects are created as orphans (the pointer @ref parent is `null`). This function is used to set the pointer @ref parent to the new parent Molecule ( @p par).
     *
     * @param par (Pointer to the) new parent Molecule.
     * @return true Setting was successfull (it means that this instance has been an orphan before).
     * @return false This instance has already have a valid parent Molecule and the new parent is not set.
     */
    bool SetParent(Molecule *par)
    {
        if (parent == NULL)
        {
            parent = par;
            return true;
        }
        else
            return false;
    };
    /**
     * @brief Get the Atom index according to the index value.
     * Returns the label of an Atom valid inside a @ref Molecule.
     * The value of @p index is a label valid inside the itramolecular field class.
     * E.g., in a bond, there are two atoms labeled I and J (or 1 and 2, respectively); this function returns their labels in the parent Molecule.
     *
     * @param index Label of an Atom inside this intramolecular field class.
     * @return int Label of the same Atom inside the @ref parent Molecule.
     */
    virtual int GetAtom(int index) const = 0;
    /**
     * @brief Calculate forces caused by this intramolecular field.
     * Calculate forces caused by this bond/angle/dihedral/... and add them to @ref Atom::force.
     * Save the potential energy and force virial to the @ref parent Molecule @ref Molecule::Epot and @ref Molecule::virial.
     * Return the potential energy.
     *
     * @par See also
     * @ref CalculateEpot() const.
     *
     * @return Potential energy caused by this intramolecular term.
     */
    virtual double CalculateForces() = 0;
    /**
     * @brief Calculate potential energy caused by this intramolecular field.
     * Similar to @ref CalculateForces(), but calculates only the potential energy and does not change anything neither in the @ref parent Molecule,
     * nor in the SimulatedSystem.
     *
     * @par See also
     * @ref CalculateForces().
     *
     * @return Potential energy caused by this intramolecular term.
     */
    virtual double CalculateEpot() const = 0;
    /**
     * @brief Clone itself and get adopted by another @ref parent Molecule.
     * Returns a copy of itself and migrate to the @p newparent Molecule.
     * Useful when copying Molecule object.
     * Every intramolecular term is then cloned by this method and assigned to the new Molecule object.
     *
     * @par See also
     * @ref AbstractIntraMolField(const AbstractIntraMolField&), @ref SetParent(Molecule*)
     *
     * @param newparent Parent Molecule of the newly created AbstractIntraMolField object.
     * @return (Pointer to) the newly created AbstractIntraMolField object.
     */
    virtual AbstractIntraMolField *copy(Molecule *newparent) const = 0;
    /**
     * @brief Get the field type.
     * Provides access to the private member @ref type.
     * 
     * @return Type of the intramolecular field.
     */
    int GetType() const { return type; }
    /**
     * @brief Get the field parameters to print info about.
     * The number of parameters depends on @ref type.
     * 
     * @param fieldname Description of field type (e.g. 'harm' for harmonic).
     * @param params Paramters to be obtained.
     * @return Number of parameters given.
     */
    virtual int GetParams(std::string &fieldname, std::vector<double> &params) const = 0;
};

#endif