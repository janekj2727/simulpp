/*
 * A class to store molecular info for trial simulations of contrained dynamics
 * Author JJ, Date Jan 2021
 */

#include <cstring>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>

#include "Vector.hpp"
#include "Matrix.hpp"
#include "Atom.hpp"
#include "Molecule.hpp"
#include "ConstraintBond.hpp"
#include "HarmonicBond.hpp"
#include "HarmonicAngle.hpp"
#include "CosineDihedral.hpp"
#include "math_utils.hpp"
#include "MDTable.hpp"
#include "units.hpp"
#include "general_utils.hpp"

// #define DEBUG

// #define SUM_SQ3(X) (X[0] * X[0] + X[1] * X[1] + X[2] * X[2]) for Vector not needed

// Default constructor
Molecule::Molecule(int NAtoms, std::string molname)
{
    curNoAtoms = 0;         // current number of atoms in molecule (initialization)
    curNoCBonds = 0;        // current number of constraints in molecule (initialization)
    curNoIntraMolF = 0;     // current number of flexible bonds
    noAtoms = NAtoms;       // number of atoms in molecule
    noConstrBonds = 0;      // number of CONSTRAINED bonds
    constrBonds = NULL;     // array of constrained bonds
    noIntraMolFields = 0;   // number of bonds, angles, dihedrals,...
    intraMolFields.clear(); // array of bonds, angles, dihedrals,...
    noInnerDegF = NAtoms * 3 - 3;
    molMass = 0.0;
    name = molname;
    hashname = 0;
    int i;
    for (i = 0; i < (int)name.length(); i++)
    {
        hashname += i * name[i];
    }
    if ((atoms = (Atom *)calloc(noAtoms, sizeof(Atom))) == NULL)
    {
        print_warning(0, "Allocation of memory failed (atom array)!!!\n");
    }
    connectivityMatrix = new Matrix(noAtoms, noAtoms);
    // connectivityMatrix = new int *[noAtoms];
    // for (int i = 0; i < noAtoms; ++i)
    // {
    //     connectivityMatrix[i] = new int[noAtoms];
    // }
    // for (i = 0; i < noAtoms; i++)
    // {
    //     for (j = 0; j < noAtoms; j++)
    //     {
    //         connectivityMatrix[i][j] = 0;
    //     }
    // }
}

// Copy constructor
Molecule::Molecule(const Molecule &mol)
{
    int i;
    std::vector<AbstractIntraMolField *>::const_iterator it;

    noAtoms = mol.noAtoms;
    curNoAtoms = noAtoms;
    noConstrBonds = mol.noConstrBonds;
    curNoCBonds = noConstrBonds;
    noIntraMolFields = mol.noIntraMolFields;
    curNoIntraMolF = noIntraMolFields;
    name = mol.name;
    hashname = mol.hashname;

    if ((atoms = (Atom *)calloc(noAtoms, sizeof(Atom))) == NULL)
    {
        print_warning(0, "Allocation of memory failed (atom array)!!!\n");
    }

    InitConstrBond(noConstrBonds);
    // InitBonds(noIntraMolFields);

    for (i = 0; i < noAtoms; i++)
    {
        // tempAtom = Atom(mol.atoms[i].mass, mol.atoms[i].charge, mol.atoms[i].name, mol.atoms[i].LJid, mol.atoms[i].R->GetNumberOfRows());
        atoms[i] = mol.atoms[i];
    }

    for (i = 0; i < noConstrBonds; i++)
    {
        constrBonds[i] = ConstraintBond(mol.constrBonds[i].length, mol.constrBonds[i].atomI,
                                        mol.constrBonds[i].atomJ, this);
    }
    intraMolFields.reserve(mol.noIntraMolFields);
    for (it = mol.intraMolFields.begin(); it != mol.intraMolFields.end(); it++)
    {
        intraMolFields.push_back((*it)->copy(this)); // make copy and assign it to the new parent (this Molecule)
    }

    // connectivityMatrix = new int *[noAtoms];
    // for (i = 0; i < noAtoms; i++)
    // {
    //     connectivityMatrix[i] = new int[noAtoms];
    // }
    // for (i = 0; i < noAtoms; i++)
    // {
    //     for (j = 0; j < noAtoms; j++)
    //     {
    //         connectivityMatrix[i][j] = mol.connectivityMatrix[i][j];
    //     }
    // }
    connectivityMatrix = new Matrix(*mol.connectivityMatrix);

    noInnerDegF = mol.noInnerDegF;
    molMass = mol.molMass;
}

// default destructor
Molecule::~Molecule()
{
    int i;
    std::vector<AbstractIntraMolField *>::const_iterator it;
    for (i = 0; i < noAtoms; i++)
    {
        atoms[i].~Atom(); // delete matrices allocated by Atom operator=
    }

    free((void *)atoms); // atoms are allocated by calloc in Molecule constructor, this does not call destructor...!!!
    atoms = NULL;

    if (noConstrBonds != 0)
    {
        for (i = 0; i < noConstrBonds; i++)
        {
            constrBonds[i].~ConstraintBond(); // empty destructor
        }
        free((void *)constrBonds);
        constrBonds = NULL;
    }
    for (it = intraMolFields.begin(); it != intraMolFields.end(); it++)
    {
        delete (*it);
    }
    intraMolFields.clear();

    // for (i = 0; i < noAtoms; i++)
    // {
    //     delete[] connectivityMatrix[i];
    // }
    // delete[] connectivityMatrix;
    // connectivityMatrix = nullptr;
    delete connectivityMatrix;
}

// calculate intramolecular forces – bonds, angles,...
double Molecule::CalculateIntraMolForces(double *Epotintramol, double *virialintramol)
{
    int i;
    double EpotIntra = 0.0;
    std::vector<AbstractIntraMolField *>::const_iterator it;
    for (i = 0; i < INTRA_MOL_TYPES; i++)
    {
        Epot[i] = 0.0;
        virial[i] = 0.0;
    }
    for (it = intraMolFields.begin(); it != intraMolFields.end(); it++)
    {
        EpotIntra += (*it)->CalculateForces();
    }
    for (i = 0; i < INTRA_MOL_TYPES; i++)
    {
        Epotintramol[i] += Epot[i];
        virialintramol[i] += virial[i];
    }

    return EpotIntra;
}

// calculate intramolecular potential energy
double Molecule::CalculateIntraMolEpot() const
{
    double EpotIntra = 0.0;
    std::vector<AbstractIntraMolField *>::const_iterator it;
    for (it = intraMolFields.begin(); it != intraMolFields.end(); it++)
    {
        EpotIntra += (*it)->CalculateEpot();
    }
    return EpotIntra;
}

// initialize SHAKE (calculate R^PI)
int Molecule::InitShake(int shakeType)
{
    // sum of R for Gear R^PI
    int i, k;
    for (i = 0; i < noAtoms; i++)
    {
        for (k = 0; k < 3; k++)
        {
            atoms[i].rPI[k] = atoms[i].R->SumColumn(k, shakeType);
        }
    }
    return 0;
}

// zero forces in atoms (before force calculation)
void Molecule::ZeroForces()
{
    int i;

    for (i = 0; i < noAtoms; i++)
    {
        atoms[i].ZeroForces();
    }
}

// initialization
int Molecule::AddAtom(double atomMass, double atomCharge, char *atomName, int LJident, int Rsize)
{
    // Atom *tempAtom;
    // tempAtom = new Atom(atomMass, atomCharge, atomName, LJident, Rsize);
    Atom tempAtom(atomMass, atomCharge, atomName, LJident, Rsize, true);
    atoms[curNoAtoms] = tempAtom;
    molMass += atomMass;
    return ++curNoAtoms;
}

// initialize constrained bonds array
int Molecule::InitConstrBond(int noConstrB)
{
    noConstrBonds = noConstrB;
    if ((noConstrBonds > 0) && ((constrBonds = (ConstraintBond *)calloc(noConstrBonds, sizeof(ConstraintBond))) == NULL))
    {
        print_warning(0, "Allocation of memory failed (constr bonds array)!!!\n");
        return 1;
    }
    return 0;
}

// initialization of constrained bonds
int Molecule::AddConstrBond(double bondLength, int indexI, int indexJ, bool trueBond)
{
    constrBonds[curNoCBonds] = ConstraintBond(bondLength, indexI, indexJ, this);
    if (trueBond)
    {
        connectivityMatrix->operator()(indexI, indexJ) = 1;
        connectivityMatrix->operator()(indexJ, indexI) = 1;
    }
    return ++curNoCBonds;
}

// initialize array of flexible bonds
int Molecule::InitBonds(int noB)
{
    intraMolFields.reserve(noIntraMolFields + noB);
    return 0;
}

// initialize array of angles
int Molecule::InitAngles(int noA)
{
    intraMolFields.reserve(noIntraMolFields + noA);
    return 0;
}

// initialize array of dihedrals
int Molecule::InitDihedrals(int noD)
{
    intraMolFields.reserve(noIntraMolFields + noD);
    return 0;
}

// add flexible bond
int Molecule::AddBond(char *bondDescription, double u_eng)
{
    char type[10];
    int atom1, atom2;
    double strength, length;
    sscanf(bondDescription, "%s %d %d %lf %lf", type, &atom1, &atom2, &strength, &length);
    // strength p.u. eng/AA^2
    strength /= u_eng;
    if ((length <= 0) || (length > 1000))
    {
        print_warning(0, "Wrong length of (flexible) bond in .field file\n",
                      "    Given '" + std::string(bondDescription) + "', length thus: " + std::to_string(length) + "\n");
        return -20;
    }
    if ((atom1 > noAtoms) || (atom2 > noAtoms) || (atom1 < 1) || (atom2 < 1) || (atom1 == atom2))
    {
        print_warning(0, "Wrong atom number in bond definition in .field file\n",
                      "    Given '" + std::string(bondDescription) + "', atoms: " + std::to_string(atom1) + " and " + std::to_string(atom2) + "\n");
        return -21;
    }
    if (strength < 0.0) // equal to 0.0 enabled to enable fictitious bonds (testing)
    {
        print_warning(0, "Wrong strength of (flexible) bond in .field file\n",
                      "    Given '" + std::string(bondDescription) + "', strength thus: " + std::to_string(strength) + "\n");
        return -22;
    }

    if (strstr(type, "harm") != NULL)
    {
        intraMolFields.push_back(new HarmonicBond(this, atom1 - 1, atom2 - 1, length, strength));
        noIntraMolFields++;
        connectivityMatrix->operator()(atom1 - 1, atom2 - 1) = 1;
        connectivityMatrix->operator()(atom2 - 1, atom1 - 1) = 1;
        return ++curNoIntraMolF;
    }
    else
    {
        print_warning(0, "Unknown type of bond in .field file\n",
                      "    Given '" + std::string(type) + "', expected 'harm'...\n");
        return -50;
    }
    return -50;
}

// add dihedral
int Molecule::AddDihedral(char *dihedralDescription, double u_eng)
{
    char type[10];
    int atom1, atom2, atom3, atom4;
    double A, delta;
    int m;
    bool proper = true; // proper or improper dihedral

    sscanf(dihedralDescription, "%s %d %d %d %d %lf %lf %d", type, &atom1, &atom2, &atom3, &atom4, &A, &delta, &m);

    if ((delta <= -180.0) || (delta > 180))
    {
        print_warning(0, "Wrong delta in dihedral in .field file\n",
                      "    Given '" + std::string(dihedralDescription) + "', delta thus: " + std::to_string(delta) + "\n",
                      "    Expected number between -180 and 180 degrees.\n");
        return -60;
    }
    if ((atom1 > noAtoms) || (atom2 > noAtoms) || (atom3 > noAtoms) || (atom1 < 1) || (atom2 < 1) || (atom3 < 1) || (atom1 == atom2) || (atom2 == atom3) || (atom1 == atom3) || (atom4 > noAtoms) || (atom4 < 1) || (atom3 == atom4) || (atom2 == atom4) || (atom1 == atom4))
    {
        print_warning(0, "Wrong atom number in dihedral definition in .field file\n",
                      "    Given '" + std::string(dihedralDescription) + "', atoms: " + std::to_string(atom1) + ", " + std::to_string(atom2) + ", " + std::to_string(atom3) + " and " + std::to_string(atom4) + "\n");
        return -61;
    }
    if (fabs(A) > 100000.0)
    {
        print_warning(0, "Wrong strength (A) of dihedral in .field file\n",
                      "    Given '" + std::string(dihedralDescription) + "', strength thus: " + std::to_string(A) + "\n");
        return -62;
    }
    if ((m < 0) || (m > 12))
    {
        print_warning(0, "Wrong multiplicity (number of minima, m) of dihedral in .field file\n",
                      "    Given '" + std::string(dihedralDescription) + "', m thus: " + std::to_string(m) + "\n");
        return -63;
    }

    if (strchr(type, '*') != nullptr)
    {
        proper = false;
    }

    if (strstr(type, "cos") != NULL)
    {
        intraMolFields.push_back(new CosineDihedral(this, atom1 - 1, atom2 - 1, atom3 - 1, atom4 - 1, A / u_eng, delta * M_PI / 180.0, m, proper));
        noIntraMolFields++;
        return ++curNoIntraMolF;
    }
    else
    {
        print_warning(0, "Unknown type of dihedral in .field file\n",
                      "    Given '" + std::string(type) + "', expected 'cos'...\n");
        return -64;
    }
    return -64;
}

// add flexible angle
int Molecule::AddAngle(char *angleDescription, double u_eng)
{
    char type[10];
    int atom1, atom2, atom3;
    double strength, angl;
    sscanf(angleDescription, "%s %d %d %d %lf %lf", type, &atom1, &atom2, &atom3, &strength, &angl);

    if ((angl < 0.0) || (angl > 180.0))
    {
        print_warning(0, "Wrong value of equilibrium angle in .field file\n",
                      "    Given '" + std::string(angleDescription) + "', angle thus: " + std::to_string(angl) + "\n");
        return -51;
    }
    angl *= M_PI / 180.0;
    if ((atom1 > noAtoms) || (atom2 > noAtoms) || (atom3 > noAtoms) || (atom1 < 1) || (atom2 < 1) || (atom3 < 1) || (atom1 == atom2) || (atom2 == atom3) || (atom1 == atom3))
    {
        print_warning(0, "Wrong atom number in angle definition in .field file\n",
                      "    Given '" + std::string(angleDescription) + "', atoms: " + std::to_string(atom1) + ", " + std::to_string(atom2) + " and " + std::to_string(atom3) + "\n");
        return -52;
    }
    if (strength <= 0.0)
    {
        print_warning(0, "Wrong strength of (flexible) angle in .field file\n",
                      "    Given '" + std::string(angleDescription) + "', strength thus: " + std::to_string(strength) + "\n");
        return -53;
    }
    // strength p.u. eng/AA^2
    strength /= u_eng;

    if (strstr(type, "harm") != NULL)
    {
        intraMolFields.push_back(new HarmonicAngle(this, atom1 - 1, atom2 - 1, atom3 - 1, angl, strength));
        noIntraMolFields++;
        return ++curNoIntraMolF;
    }
    else
    {
        print_warning(0, "Unknown type of angle in .field file\n",
                      "    Given '" + std::string(type) + "', expected 'harm'...\n");
        return -54;
    }
    return -54;
}

// calculate center of mass
void Molecule::CalculateCOM(double *centerOfMass)
{
    int i, m;

    for (m = 0; m < 3; m++)
    {
        centerOfMass[m] = 0.0;
    }
    for (i = 0; i < noAtoms; i++)
    {
        for (m = 0; m < 3; m++)
        {
            centerOfMass[m] += atoms[i].GetPosition(m) * atoms[i].mass;
        }
    }
    for (m = 0; m < 3; m++)
    {
        centerOfMass[m] /= molMass;
    }
}

// move whole molecule in direction by length Angstroms
void Molecule::Move(double length, int direction, std::vector<int> &Rrows)
{
    int i;
    std::vector<int>::iterator j;

    for (i = 0; i < noAtoms; i++)
    {
        for (j = Rrows.begin(); j != Rrows.end(); j++)
        {
            atoms[i].R->operator()(*j, direction) += length;
        }
    }
}

// calculate translational kinetic energy
double Molecule::CalculateEkinTrans(double timeStep)
{
    double Etrans = 0.0;
    int k, l;
    Vector velocity(3);

    for (k = 0; k < noAtoms; k++)
    {
        for (l = 0; l < 3; l++)
        {
            velocity[l] += atoms[k].mass * atoms[k].R->operator()(1, l) / timeStep;
        }
    }
    Etrans = 0.5 / molMass * CalculateScalarProduct(velocity, velocity);
    return Etrans;
}

// calculate (total) kinetic energy
double Molecule::CalculateEkin(double timeStep, int velocityIndex) const
{
    int j;
    double Ekin = 0.0;

    for (j = 0; j < noAtoms; j++)
    {
        Ekin += atoms[j].mass * atoms[j].R->RowRowDotProduct(velocityIndex, velocityIndex) / (timeStep * timeStep);
        // for (k = 0; k < 3; k++)
        // {
        //     Ekin += atoms[j].mass * atoms[j].R->operator()(1, k) * atoms[j].R->operator()(1, k) / timeStep / timeStep;
        // }
    }

    return 0.5 * Ekin;
}

// return maximum error of constrained bonds
double Molecule::MaxConstrErrSquared() const
{
    int j;
    double maxCErr = 0.0;

    for (j = 0; j < noConstrBonds; j++)
    {
        maxCErr = fmax(maxCErr, pow((constrBonds[j].GetCurrentLength() - constrBonds[j].length) / constrBonds[j].length, 2.0));
    }
    return maxCErr;
}

double Molecule::TotalCharge() const
{
    int i;
    double charge = 0.0;

    for (i = 0; i < noAtoms; i++)
    {
        charge += atoms[i].charge;
    }

    return charge;
}

bool Molecule::HasAnyCharged() const
{
    int i;

    for (i = 0; i < noAtoms; i++)
    {
        if (atoms[i].charge != 0.0)
        {
            return true;
        }
    }

    return false;
}

int Molecule::ConnectivityCheck()
{
    Matrix original(noAtoms, noAtoms);
    Matrix *shortest;
    int i, j, k, l;
    std::vector<AbstractIntraMolField *>::const_iterator it;

    original = *connectivityMatrix;

    shortest = shortestpath(original); // allocation is done inside

    for (i = 0; i < noAtoms; i++)
    {
        for (j = 0; j < noAtoms; j++)
        {
            connectivityMatrix->operator()(i, j) = (int)shortest->operator()(i, j);
            // check if Molecule not divided into separated independent parts
            if ((i != j) && ((connectivityMatrix->operator()(i, j) < 1) || (connectivityMatrix->operator()(i, j) > noAtoms)))
            {
                print_warning(0, "Atoms in one molecule not interconnected by bonds!!!\n",
                              "    Wrong distance between atoms: " + std::to_string(i + 1) + " and " + std::to_string(j + 1) + ".\n");
                return 67;
            }
        }
    }
    delete shortest;
    shortest = nullptr;

    for (it = intraMolFields.begin(); it != intraMolFields.end(); it++)
    {
        if (((*it)->GetType() == 1) && // angles
            ((i = (*it)->GetAtom(1), j = (*it)->GetAtom(2), k = (*it)->GetAtom(3), connectivityMatrix->operator()(i, j) != 1.0) ||
             (connectivityMatrix->operator()(j, k) != 1.0)))
        {
            print_warning(0, "Atoms within one angle are not connected by bonds!!!\n",
                          "Connectivity check fails at angle consisting of atoms: " + std::to_string(i) + ", " + std::to_string(j) + " and " + std::to_string(k) + ".\n");
#ifdef DEBUG
            std::cerr << "Connectivity matrix:\n"
                      << *connectivityMatrix;
#endif
            return 65;
        }
        if (((*it)->GetType() == 2) && // dihedrals
            ((i = (*it)->GetAtom(1), j = (*it)->GetAtom(2), k = (*it)->GetAtom(3), l = (*it)->GetAtom(4), connectivityMatrix->operator()(i, j) != 1.0) ||
             (connectivityMatrix->operator()(j, k) != 1.0) || (connectivityMatrix->operator()(k, l) != 1.0)))
        {
            print_warning(0, "Atoms within one dihedral are not connected by bonds!!!\n",
                          "Connectivity check fails at dihedral consisting of atoms: " + std::to_string(i) + ", " + std::to_string(j) + ", " + std::to_string(k) + " and " + std::to_string(l) + ".\n");
#ifdef DEBUG
            std::cerr << "Connectivity matrix:\n"
                      << *connectivityMatrix;
#endif
            return 66;
        }
    }

    return 0;
}

// get connectivityMatrix
int Molecule::GetDistance(int indexI, int indexJ) const
{
    return connectivityMatrix->operator()(indexI, indexJ);
}

// print molecule info to .prt file
void Molecule::PrintInfo(std::ofstream &fprt, int mol_no, int number_of_this, double u_eng, std::string engunit) const
{
    int i;
    char aux[25];
    std::vector<std::string> atom_names;
    int precision;

    // summary
    fprt << "-----------------------------------------------------------------\n\n";
    fprt << "## Molecule " << mol_no + 1 << ": " << name << "\n\n";
    fprt << "### Summary\n\n";
    int no_intramolf[INTRA_MOL_TYPES] = {0};
    for (auto it = intraMolFields.begin(); it != intraMolFields.end(); it++)
    {
        no_intramolf[(*it)->GetType()]++;
    }

    MDTable<std::string, double> summary({"Variable", "Value"});
    summary.setColumnPrecision({0, 9});

    summary.addRow("Number of atoms", (double)noAtoms);
    summary.addRow("Number of molecules of this type", (double)(number_of_this));
    summary.addRow("Molecular mass [g/mol]", molMass * U_MASS_GMOL);
    summary.addRow("Total charge [e]", TotalCharge() * U_CHARGE_E);
    summary.addRow("Number of constraints (rigid bonds)", noConstrBonds);
    summary.addRow("Number of flexible bonds", no_intramolf[0]);
    summary.addRow("Number of angles", no_intramolf[1]);
    summary.addRow("Number of dihedrals", no_intramolf[2]);
    summary.addRow("Number of improper torsions", no_intramolf[5]);
    summary.addRow("Number of intramolecular dispersive interacions", no_intramolf[3]);
    summary.addRow("Number of intramolecular electrostatic interactions", no_intramolf[4]);

    summary.print(fprt);

    // atom list
    fprt << "\n### Atom list\n\n";

    MDTable<int, std::string, double, double> atomlist({"No.", "Name", "Mass [g/mol]", "Charge [e]"});

    for (i = 0; i < noAtoms; i++)
    {
        atomlist.addRow(i + 1, std::string(atoms[i].name), atoms[i].mass * U_MASS_GMOL, atoms[i].charge * U_CHARGE_E);
    }
    atomlist.print(fprt);

    // connectivity table
    fprt << "\n### Connectivity table\n\n";
    for (i = 0; i < noAtoms; i++)
    {
        sprintf(aux, "%02d: %s", i + 1, atoms[i].name);
        atom_names.push_back(std::string(aux));
    }
    precision = fprt.precision();
    fprt << std::setprecision(3);
    connectivityMatrix->PrintMDTable(fprt, atom_names, atom_names);
    fprt << std::setprecision(precision);

    // constraints if any
    MDTable<std::string, std::string, double> cbonds({"atom 1", "atom 2", "length [AA]"});
    Atom *atI, *atJ;
    if (noConstrBonds > 0)
    {
        fprt << "\n### Constrained (rigid) bonds\n\n";
        for (i = 0; i < noConstrBonds; i++)
        {
            atI = &(atoms[constrBonds[i].atomI]);
            atJ = &(atoms[constrBonds[i].atomJ]);
            cbonds.addRow(std::to_string(constrBonds[i].atomI + 1) + ": " + atI->name, std::to_string(constrBonds[i].atomJ + 1) + ": " + atJ->name, constrBonds[i].length);
        }
        cbonds.print(fprt);
    }

    // intramolecular interactions
    std::vector<double> params;
    std::string fieldname;
    int nopar;
    MDTable<std::string, std::string, std::string, double, double, double, double> bonds({"atom 1", "atom 2", "bond type", "force constant " + engunit, "length [AA]", "par 3", "par 4"});
    MDTable<std::string, std::string, std::string, std::string, double, double, double, double> angles({"atom 1", "atom 2", "atom 3", "bond type", "force constant " + engunit, "angle [rad]", "par 3", "par 4"});
    MDTable<std::string, std::string, std::string, std::string, std::string, double, double, double, double, double> dihedrals({"atom 1", "atom 2", "atom 3", "atom 4", "bond type", "force constant " + engunit, "angle [rad]", "par 3", "par 4", "par 5"});
    MDTable<std::string, std::string, std::string, std::string, std::string, double, double, double, double, double> impropers({"atom 1", "atom 2", "atom 3", "atom 4", "bond type", "force constant " + engunit, "angle [rad]", "par 3", "par 4", "par 5"});
    MDTable<std::string, std::string, std::string, double, double, double, double> disperse({"atom 1", "atom 2", "interaction type", "energy constant " + engunit, "length [AA]", "par 3", "par 4"});
    MDTable<std::string, std::string, double, double> elstat({"atom 1", "atom 2", "product of charges (incl. scaling) ['e^2']", "energy shift " + engunit});
    for (i = 0; i < INTRA_MOL_TYPES; i++)
    {
        if (no_intramolf[i] > 0)
        {
            switch (i)
            {
            case 0:
                fprt << "\n### Flexible bonds\n\n";
                for (auto it = intraMolFields.begin(); it != intraMolFields.end(); it++)
                {
                    if ((*it)->GetType() == 0)
                    {
                        nopar = (*it)->GetParams(fieldname, params);
                        bonds.addRow(std::to_string((*it)->GetAtom(1) + 1) + ": " + atoms[(*it)->GetAtom(1)].name, std::to_string((*it)->GetAtom(2) + 1) + ": " + atoms[(*it)->GetAtom(2)].name, fieldname, params[0] * u_eng, params[1], ((nopar > 2) ? params[2] : 0.0), ((nopar > 3) ? params[3] : 0.0));
                    }
                }
                bonds.print(fprt);
                fprt << "Potential energy formula depends on the bond type:\n";
                fprt << "- `harm`: $E_{\\mathrm{pot}} = \\frac{1}{2} k (r_{ij} - r_0)^2$, where $k$ is the force constant and $r_0$ is the length from the table above\n";
                // next types described here in future...
                break;
            case 1:
                fprt << "\n### Flexible angles\n\n";
                for (auto it = intraMolFields.begin(); it != intraMolFields.end(); it++)
                {
                    if ((*it)->GetType() == 1)
                    {
                        nopar = (*it)->GetParams(fieldname, params);
                        angles.addRow(std::to_string((*it)->GetAtom(1) + 1) + ": " + atoms[(*it)->GetAtom(1)].name, std::to_string((*it)->GetAtom(2) + 1) + ": " + atoms[(*it)->GetAtom(2)].name, std::to_string((*it)->GetAtom(3) + 1) + ": " + atoms[(*it)->GetAtom(3)].name, fieldname, params[0] * u_eng, params[1], ((nopar > 2) ? params[2] : 0.0), ((nopar > 3) ? params[3] : 0.0));
                    }
                }
                angles.print(fprt);
                fprt << "Potential energy formula depends on the angle type:\n";
                fprt << "- `harm`: $E_{\\mathrm{pot}} = \\frac{1}{2} k (\\theta_{ijk} - \\theta_0)^2$, where $k$ is the force constant and $\\theta_0$ is the angle from the table above\n";
                // next types described here in future...
                break;
            case 2:
                fprt << "\n### Dihedrals\n\n";
                for (auto it = intraMolFields.begin(); it != intraMolFields.end(); it++)
                {
                    if ((*it)->GetType() == 2)
                    {
                        nopar = (*it)->GetParams(fieldname, params);
                        dihedrals.addRow(std::to_string((*it)->GetAtom(1) + 1) + ": " + atoms[(*it)->GetAtom(1)].name, std::to_string((*it)->GetAtom(2) + 1) + ": " + atoms[(*it)->GetAtom(2)].name, std::to_string((*it)->GetAtom(3) + 1) + ": " + atoms[(*it)->GetAtom(3)].name, std::to_string((*it)->GetAtom(4) + 1) + ": " + atoms[(*it)->GetAtom(4)].name, fieldname, params[0] * u_eng, params[1], ((nopar > 2) ? params[2] : 0.0), ((nopar > 3) ? params[3] : 0.0), ((nopar > 4) ? params[4] : 0.0));
                    }
                }
                dihedrals.print(fprt);
                fprt << "Potential energy formula depends on the dihedral type:\n";
                fprt << "- `cos`: $E_{\\mathrm{pot}} = A \\left[ 1 + \\cos\\left(m\\varphi_{ijkl} - \\delta \\right)\\right]$, where $A$ is the force constant, $\\delta$ is the angle and $m$ is the 'par 3' from the table above\n";
                // next types described here in future...
                break;
            case 5:
                fprt << "\n### Improper torsions\n\n";
                for (auto it = intraMolFields.begin(); it != intraMolFields.end(); it++)
                {
                    if ((*it)->GetType() == 5)
                    {
                        nopar = (*it)->GetParams(fieldname, params);
                        impropers.addRow(std::to_string((*it)->GetAtom(1) + 1) + ": " + atoms[(*it)->GetAtom(1)].name, std::to_string((*it)->GetAtom(2) + 1) + ": " + atoms[(*it)->GetAtom(2)].name, std::to_string((*it)->GetAtom(3) + 1) + ": " + atoms[(*it)->GetAtom(3)].name, std::to_string((*it)->GetAtom(4) + 1) + ": " + atoms[(*it)->GetAtom(4)].name, fieldname, params[0] * u_eng, params[1], ((nopar > 2) ? params[2] : 0.0), ((nopar > 3) ? params[3] : 0.0), ((nopar > 4) ? params[4] : 0.0));
                    }
                }
                impropers.print(fprt);
                fprt << "Potential energy formula depends on the torsion (dihedral) type:\n";
                fprt << "- `cos`: $E_{\\mathrm{pot}} = A \\left[ 1 + \\cos\\left(m\\varphi_{ijkl} - \\delta \\right)\\right]$, where $A$ is the force constant, $\\delta$ is the angle and $m$ is the 'par 3' from the table above\n";
                // next types described here in future...
                break;
            case 3:
                fprt << "\n### Intramolecular disperse interactions (LJ interactions etc.)\n\n";
                for (auto it = intraMolFields.begin(); it != intraMolFields.end(); it++)
                {
                    if ((*it)->GetType() == 3)
                    {
                        nopar = (*it)->GetParams(fieldname, params);
                        disperse.addRow(std::to_string((*it)->GetAtom(1) + 1) + ": " + atoms[(*it)->GetAtom(1)].name, std::to_string((*it)->GetAtom(2) + 1) + ": " + atoms[(*it)->GetAtom(2)].name, fieldname, params[0] * u_eng, params[1], ((nopar > 2) ? params[2] : 0.0), ((nopar > 3) ? params[3] : 0.0));
                    }
                }
                disperse.print(fprt);
                fprt << "Potential energy formula depends on the interaction type:\n";
                fprt << "- `lj`: $E_{\\mathrm{LJ}} = 4 \\varepsilon \\left[\\left(\\frac{\\sigma}{r_{ij}}\\right)^{12} - \\left(\\frac{\\sigma}{r_{ij}}\\right)^6\\right]$, where $\\varepsilon$ is the energy constant and $\\sigma$ is the length from the table above\n";
                // next types described here in future...
                break;
            case 4:
                fprt << "\n### Intramolecular electrostatic interactions\n\n";
                for (auto it = intraMolFields.begin(); it != intraMolFields.end(); it++)
                {
                    if ((*it)->GetType() == 4)
                    {
                        nopar = (*it)->GetParams(fieldname, params);
                        elstat.addRow(std::to_string((*it)->GetAtom(1) + 1) + ": " + atoms[(*it)->GetAtom(1)].name, std::to_string((*it)->GetAtom(2) + 1) + ": " + atoms[(*it)->GetAtom(2)].name, params[0] * U_CHARGE_E * U_CHARGE_E, params[1] * u_eng);
                    }
                }
                elstat.print(fprt);
                fprt << "Potential energy of electrostatic interactions:\n";
                fprt << "- $E_{\\mathrm{elstat}} = \\frac{\\omega q_i q_j}{4\\pi\\epsilon_0}\\left(\\frac{1}{r_{ij}}-\\mathrm{shift}\\right)$, where $\\omega q_i q_j$ is the charge product ($\\omega$ is the scaling factor) and $\\mathrm{shift}$ is the energy shift from the table above\n";
                // next types described here in future...
                break;
            default:
                break;
            }
        }
    }
}