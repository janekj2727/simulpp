#ifndef MOLECULEHEADER
#define MOLECULEHEADER

/*
#include "../mymath/Vector.hpp"
#include "../mymath/Matrix.hpp"
*/

#define INTRA_MOL_TYPES 6

#include <vector>

class ConstraintBond;
class AbstractIntraMolField;

class Molecule
{
private:
    int curNoAtoms;     // current number of atoms in molecule (initialization)
    int curNoCBonds;    // current number of constraints in molecule (initialization)
    int curNoIntraMolF; // current number of flexible bonds (initialization)
    // int **connectivityMatrix; // matrix of shortest paths (minimum number of bonds between atoms)
    Matrix *connectivityMatrix;

public:
    Atom *atoms;                                           // array of atoms in molecule
    int noAtoms;                                           // number of atoms in molecule
    int noConstrBonds;                                     // number of CONSTRAINED bonds
    ConstraintBond *constrBonds;                           // array of constrained bonds
    int noIntraMolFields;                                  // number of bonds, angles, dihedrals,... in the following array
    std::vector<AbstractIntraMolField *> intraMolFields;   // array of bonds, angles, dihedrals,...
    double molMass;                                        // total mass of the molecule
    int noInnerDegF;                                       // number of inner degrees of freedom (except translation of COM)
    std::string name;                                      // name of the molecule in .field file
    int hashname;                                          // hash of the molecule name in .field file
    double virial[INTRA_MOL_TYPES], Epot[INTRA_MOL_TYPES]; // virial and potential of bonds(0), angles(1), dihedrals(2), LJ and other vdW(3) and electrostatic(4) interactions

    double CalculateIntraMolForces(double *Epotintramol, double *virial); // calculate intramolecular forces (LJ, bonds, angles...)
    double CalculateIntraMolEpot() const;                                 // calculate intramolecular potential energy
    // double CalculateConstrBondForces(double &virial, int shakeType = 0, double epsc = 1e-6,  //
    //                                  double rescaling = 1.0, double omega = 1.0);            // calculate SHAKE forces or similar
    int InitShake(int shakeType = 0);                                                        // initialize SHAKE (calculate R^PI...)
    void ZeroForces();                                                                       // zero forces in atoms (before force calculation)
    int AddAtom(double atomMass, double atomCharge, char *atomName, int LJident, int Rsize); // initialization
    int InitConstrBond(int noConstrB);                                                       // initialize constr. bonds array
    int AddConstrBond(double bondLength, int indexI, int indexJ, bool trueBond = true);      // initialization of constrained bonds
    int InitBonds(int noB);                                                                  // initialize flex. bonds array
    int InitAngles(int noA);                                                                 // initialize flex. angles array (reserve space)
    int InitDihedrals(int noD);                                                              // initialize dihedrals array (reserve space)
    int AddBond(char *bondDesciption, double u_eng);                                         // initialization of bond
    int AddAngle(char *angleDescription, double u_eng);                                      // initialization of angle
    int AddDihedral(char *dihedralDescription, double u_eng);                                // initialization of dihedrals
    void CalculateCOM(double *centerOfMass);                                                 // calculate center of mass
    double CalculateEkinTrans(double timeStep);                                              // calculate translational kinetic energy
    double CalculateEkin(double timeStep, int velocityIndex = 1) const;                      // calculate (total) kinetic energy
    void Move(double length, int direction, std::vector<int> &Rrows);                        // move whole molecule in direction by length Angstroms
    double MaxConstrErrSquared() const;                                                      // return maximum relative error of constrained bonds
    double TotalCharge() const;                                                              // calculate total charge of the molecule
    bool HasAnyCharged() const;                                                              // returns true if any atom is charged
    int ConnectivityCheck();                                                                 // calculate shortest distance and check consistency of angles and dihedrals
    int GetDistance(int indexI, int indexJ) const;                                           // get shorters distance between i and j (in number of bonds)
    void PrintInfo(std::ofstream &fprt, int mol_no,
                   int number_of_this, double u_eng, std::string engunit) const; // print info to .prt file
    Molecule(int Natoms, std::string molname);                                   // default constructor
    Molecule(const Molecule &mol);                                               // copy constructor
    ~Molecule();                                                                 // default destructor
};

#endif