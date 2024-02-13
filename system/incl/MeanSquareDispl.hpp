/*
*  Mean square displacement of atoms calculation 
*  Written for simul++ simulation package
*
*  Author: JJ
*  Date: Sep 2023
*
*/

#ifndef MSDHEADER
#define MSDHEADER
#include <cstdio>
#include <vector>
#include <string>
class Vector;
class Molecule;
/**
 * @class MeanSquareDispl
 * @ingroup QuantsMeasurements
 * @brief Mean square displacement online calculation.
 * 
 * This class enables the calculation of mean square displacement (MSD) of individual atoms grouped by their names (for vdW interactions).
 * The calculation is performed online and results are printed to `.cpa` file. 
 */
class MeanSquareDispl
{
private:
    /**
     * @brief Vector of original atom positions.
     * 
     * Size is the number of atoms in the system ( @ref SimulatedSystem::noAtomsTotal) times 3 (coordinates)
     */
    Vector *orig_pos;
    /**
     * @brief Indices of the x-coordinate of the first atom of each molecule in orig_pos.
     * 
     * These indices are needed when moving molecules...
     */
    std::vector<int> mol_index_in_orig;
    /**
     * @brief For each @ref Atom::LJid count the number of atoms. 
     * 
     * Size is thus the size of `LJsystem` or similar vdW system.
     */
    std::vector<int> no_atoms_of_type;
    /**
     * @brief Atom names (LJ names).
     * 
     * Atom names for each @ref Atom::LJid are needed to print header to `.cpa` file.
     * 
     * @see PrintHeader() const; 
     */
    std::vector<std::string> atom_names;
    /**
     * @brief Displacements for each @ref Atom::LJid type.
     * 
     * Size is thus the size of `LJsystem` or similar vdW system.
     */
    Vector *displacements;
    /**
     * @brief Pointer to SimulatedSystem::molecules to be able to get atom positions.
     */
    std::vector<Molecule> *molecules;
    /**
     * @brief Default constructor (disabled).
     * 
     * Default constructor is made private to disable its use. 
     */
    MeanSquareDispl(){};

public:
    /**
     * @brief Construct a new MeanSquareDispl class to measure atom MSD.
     * 
     * Original positions of atoms are stored, number of atom types is stored and used to allocate @ref displacements and @ref no_atoms_of_type.
     * 
     * @param noAtoms Total number of atoms in the system.
     * @param mol Pointer to vector of molecules ( @ref SimulatedSystem::molecules ).
     */
    MeanSquareDispl(int noAtoms, std::vector<Molecule> *mol);
    /**
     * @brief Destroy the MeanSquareDispl object
     * 
     * Deallocate memory used by private members.
     */
    ~MeanSquareDispl();
    /**
     * @brief Copy constructor for this class.
     * 
     * Not needed, but written for clarity.
     * 
     * @param orig Original class instance to be copied. 
     */
    MeanSquareDispl(const MeanSquareDispl& orig);
    /**
     * @brief Assignment operator overloading.
     * 
     * Not needed, but written for clarity.
     * 
     * @param orig Original class to be *copied* to the left-hand side value. 
     * @return MeanSquareDispl& Right-hand side instance is copied to the left-hand side.
     * 
     * @note This is (obviously) **NOT a constructor**. 
     */
    MeanSquareDispl& operator=(const MeanSquareDispl& orig);
    /**
     * @brief Print mean square displacements to `.cpa` file.
     * 
     * @param fcpa Pointer to `.cpa` file (c-style).
     * @param format Format of printed numer(s)
     * @return 0 upon success.
     */
    int PrintMSD(FILE *fcpa, const char *format);
    /**
     * @brief Print header to `.cpa` file.
     * 
     * Print columns header prior to measurements.
     * 
     * @param fcpa Pointer to `.cpa` file (c-style).
     * @param col_no First MSD column number.
     * @param format Format for column name.
     * @return Last MSD column number.
     */
    int PrintHeader(FILE *fcpa, int col_no, const char *format);
    /**
     * @brief Move original position to the opposite way the moleculed is moved upon applying PBC.
     * 
     * When calculating MSD we don't want molecules to return back to the simulation box from the opposite side.
     * One possibility is not to apply periodic boundary conditions and be carefull when calculating distances.
     * The other, here preferred, method is to move the original image of molecule in the opposite direction for MSD measurements.
     * 
     * @param molecule_index Index of molecule which was moved due to PBC.
     * @param direction Direction of move (coordinate): 0 for x, 1 for y and 2 for z.
     * @param distance Distance. Original image will be moved in the opposite direction.
     * @return 0 upon success.
     */
    int MoveMolecule(int molecule_index, int direction, double distance);
    /**
     * @brief Get number of distinct atom types (number of columns with MSD) in `.cpa` file
     * 
     * @return size of @ref atom_names
     */
    int GetNumberOfTypes() const;
};

#endif
