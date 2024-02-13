/*
 * .plb file utils
 *  read one plb frame, print interatomic distances etc. (for simul++)
 *
 *  Author: JJ
 *  Date:
 *
 */

#include <cstring>
#include <cstdio>
#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <cmath>

#include "file_openclose.hpp"
#include "Vector.hpp"
#include "SimulatedSystem.hpp"
#include "Molecule.hpp"
#include "ConstraintBond.hpp"
#include "AbstractIntraMolField.hpp"
#include "LJsystem.hpp"
#include "FieldFile.hpp"
#include "units.hpp"
#include "general_utils.hpp"

char DetermineColorForMol(char *name);

FieldFile::FieldFile(char *field_name)
{
    errorstatus = 0;
    char rmode[] = "r";
    strcpy(fieldname, field_name);

    if (my_fopen_r(&fieldfile, field_name, rmode) != 0)
    {
        print_warning(0, "Cannot open .field file: " + std::string(field_name) + "\n");
        errorstatus = 7;
        goto endOfConstructor;
    }
    // erase header
    header[0] = '\0';

endOfConstructor:;
}

FieldFile::~FieldFile()
{
    if (my_fclose(&fieldfile, fieldname) != 0)
    {
        print_warning(0, "Error while closing " + std::string(fieldname) + "\n");
        errorstatus = 31;
    };
}

int FieldFile::Read(SimulatedSystem *system)
{
    // initialize simulated systemLJ->GetLJid(auxname) by reading field file
    char auxiliary[MAX_COMMENT + 1];
    char aux2[MAX_COMMENT + 1];
    char engunits[MAX_COMMENT + 1];
    int noMoleculeTypes = 0; // number of molecule types
    int i;
    int noLJ = 0;
    int returncode = 0;

    // reading and printing header
    if (fgets(header, MAX_COMMENT, fieldfile) == nullptr)
    {
        return 46;
    }

    // read units
    while ((fgets(auxiliary, MAX_COMMENT, fieldfile) != NULL) && (feof(fieldfile) == 0))
    {
        if ((auxiliary[0] == '!') || (auxiliary[0] == '#'))
        {
            continue;
        }
        else if (strstr(auxiliary, "units") != NULL)
        {
            sscanf(auxiliary, "%s %s", aux2, engunits);
            if (strstr(engunits, "K") != NULL)
            {
                system->SetEnergyUnit(U_ENERGY_KELVIN);
            }
            else if (strstr(engunits, "kcal") != NULL)
            {
                system->SetEnergyUnit(U_ENERGY_KCALMOL);
            }
            else if (strstr(engunits, "kJ") != NULL)
            {
                system->SetEnergyUnit(U_ENERGY_KJMOL);
            }
            break;
        }
        else
        {
            print_warning(1, "Unexpected directive in .field file.\n", "    Expected 'units', but received " + std::string(auxiliary) + ".\n", "    Units defaulted to K\n");
            system->SetEnergyUnit(U_ENERGY_KELVIN);
            break;
        }
    }

    // read until molecular types
    while ((fgets(auxiliary, MAX_COMMENT, fieldfile) != NULL) && (feof(fieldfile) == 0))
    {
        if (strstr(auxiliary, "molecules") != NULL)
        {
            sscanf(auxiliary, "%s %d", aux2, &noMoleculeTypes);
            break;
        }
    }

    // error if wrong number of molecular types...
    if ((noMoleculeTypes < 1) || (noMoleculeTypes > 1000000))
    {
        print_warning(0, "Wrong number of molecular types ('molecules') in .field file or not given\n",
                      "    Given '" + std::string(auxiliary) + "', converted to number: " + std::to_string(noMoleculeTypes) + "\n");
        return 9;
    }

    // read species
    for (i = 0; i < noMoleculeTypes; i++)
    {
        returncode = ReadMolecules(system);
        if (returncode != 0)
        {
            print_warning(0, "Error while reading molecules from .field file\n");
            return returncode;
        }
    }

    // // vdw - number of LJ terms
    // if ((fgets(auxiliary, MAX_COMMENT, fieldfile) != nullptr) && (strstr(auxiliary, "vdw") != 0))
    // {
    //     sscanf(auxiliary, "%s %d", aux2, &noLJ);
    //     // reading LJ terms
    //     returncode = ReadvdW(noLJ, system);
    //     if (returncode != 0)
    //     {
    //         print_warning(0, "Error during wdW terms reading from field file\n");
    //         return returncode;
    //     }
    // }
    // else if (strstr(auxiliary, "extern") != 0)
    // {
    //     // read extern potential definition...
    // }
    // else
    // {
    //     print_warning(0, "Error in final part of .field file, unknown directive (expected: vdw, extern).\n", "    Obtained: " + std::string(auxiliary) + "\n");
    //     return 24;
    // }

    // now close directive should be there (except for comments)
    while ((fgets(auxiliary, MAX_COMMENT, fieldfile) != nullptr) && (feof(fieldfile) == 0))
    {
        if ((auxiliary[0] == '#') || (auxiliary[0] == '!'))
        {
            continue;
        }
        else if (strstr(auxiliary, "close") != NULL)
        {
            return 0;
        }
        else if ((strstr(auxiliary, "vdw") != 0))
        {
            sscanf(auxiliary, "%s %d", aux2, &noLJ);
            // reading LJ terms
            returncode = ReadvdW(noLJ, system);
            if (returncode != 0)
            {
                print_warning(0, "Error during wdW terms reading from field file\n");
                return returncode;
            }
        }
        else if (strstr(auxiliary, "extern") != 0)
        {
            // read extern potential definition...
            print_warning(1, "Extern potential requested, but this feature has not been implemented yet.\n");
        }
        else
        {
            print_warning(0, "Error in final part of .field file, unknown directive (expected: vdw, extern or close).\n", "    Obtained: " + std::string(auxiliary) + "\n");
            return 24;
        }
    }

    if (strstr(auxiliary, "close") != NULL)
    {
        return 0;
    }
    print_warning(0, "Unexpected line in .field file ('close' statement not reached)!\n",
                  "    Given: " + std::string(auxiliary) + ", expected: close\n");
    return 30;
}

int FieldFile::PrintMol(char *molname, SimulatedSystem *system)
{
    FILE *mol_file;
    int i, j;
    static int k = 0; // to count overall number of atoms
    int LJid;
    int diameter;
    char name[MAX_COMMENT];
    std::string namestr;
    std::vector<std::string> bonded_atoms; // for each atom stores a string with bonded atoms
    std::vector<int> bondNumber;           // for each atom stores number of bonds
    char aux[10];
    int bAtomI, bAtomJ;
    char col;
    int first;
    char writemode[] = "w";

    // open .mol file for reading
    if (my_fopen_w(&mol_file, molname, writemode) != 0)
    {
        print_warning(0, ".mol file " + std::string(molname) + " cannot be opened!\n");
        return 8;
    }

    // print header of species
    fprintf(mol_file, "\nnumber_of_atoms=%d\n\n", system->noAtomsTotal);
    fprintf(mol_file, "atoms\n! i   atom-id  a-type charge chir nbonds bound_atoms\n");

    // for each molecule print .mol entries
    for (i = 0; i < system->noMolecules; i++)
    {
        // number of first atom in the molecule
        first = k;
        // resize bonded_atoms and bondNumber to the size needed for the molecule
        if ((int)bonded_atoms.size() < system->molecules[i].noAtoms)
        {
            bonded_atoms.resize(system->molecules[i].noAtoms);
            bondNumber.resize(system->molecules[i].noAtoms);
        }
        // empty strings before saving new info
        for (j = 0; j < system->molecules[i].noAtoms; j++)
        {
            bonded_atoms[j].clear();
            bondNumber[j] = 0;
        }
        // for each bond (IntraMolField type 0) add the info to both atoms
        for (j = 0; j < system->molecules[i].noIntraMolFields; j++)
        {
            // if not bond, continue...
            if (system->molecules[i].intraMolFields[j]->GetType() != 0)
            {
                continue;
            }
            bAtomI = system->molecules[i].intraMolFields[j]->GetAtom(1);
            bAtomJ = system->molecules[i].intraMolFields[j]->GetAtom(2);
            sprintf(aux, " %d", bAtomJ + first);
            bonded_atoms[bAtomI].append(aux);
            sprintf(aux, " %d", bAtomI + first);
            bonded_atoms[bAtomJ].append(aux);
            bondNumber[bAtomI]++;
            bondNumber[bAtomJ]++;
        }
        // for each constrained bond add the info to both atoms
        for (j = 0; j < system->molecules[i].noConstrBonds; j++)
        {
            bAtomI = system->molecules[i].constrBonds[j].atomI;
            bAtomJ = system->molecules[i].constrBonds[j].atomJ;
            sprintf(aux, " %d", bAtomJ + first);
            bonded_atoms[bAtomI].append(aux);
            sprintf(aux, " %d", bAtomI + first);
            bonded_atoms[bAtomJ].append(aux);
            bondNumber[bAtomI]++;
            bondNumber[bAtomJ]++;
        }
        // write atom lines and update k (first atom number)
        for (j = 0; j < system->molecules[i].noAtoms; j++)
        {
            LJid = system->molecules[i].atoms[j].LJid;
            diameter = diametersMap[LJid];
            GetAtomName(LJid, &namestr);
            strcpy(name, namestr.c_str());
            col = DetermineColorForMol(name);

            fprintf(mol_file, "%d %c%03d/%d-%s%03d %s %10.8f %d %d %s\n", k++, col, diameter,
                    j, name, i, name, system->molecules[i].atoms[j].charge * U_CHARGE_E, 0,
                    bondNumber[j], bonded_atoms[j].c_str());
        }
    }

    // close .mol file
    if (my_fclose(&mol_file, molname) != 0)
    {
        print_warning(0, "Error while closing " + std::string(molname) + "\n");
        return 32;
    }
    // if successful, return 0
    return 0;
}

// Read LJ (VdW) terms from .field and write molecular diameters
int FieldFile::ReadvdW(int novdw, SimulatedSystem *system)
{
    char auxiliary[MAX_COMMENT + 1];
    char aux2[MAX_COMMENT + 1];
    char params[MAX_COMMENT];
    char name1[MAX_COMMENT], name2[MAX_COMMENT];
    int current_type = 0;
    int error = 0;
    int i, j;
    // double diameters[LJmap.size()];

    // check number of LJ terms
    if ((novdw < 0) || (novdw > (int)(LJmap.size() * (LJmap.size() + 1) / 2)))
    {
        print_warning(0, "Wrong number of disperse forces terms ('vdw' directive) (negative or greater then possible combinations)\n", "    Given: '" + std::string(auxiliary) + "', thus number of vdW terms defined: " + std::to_string(novdw) + "\n", "    Number of atoms with different names): " + std::to_string(LJmap.size()) + "\n");
        return 25;
    }

    // read vdw terms line by line
    for (i = 0; i < novdw;)
    {
        if (fgets(auxiliary, MAX_COMMENT, fieldfile) == NULL)
        {
            print_warning(0, "End of .field file reached sooner then expected (from number of 'vdw' terms)\n");
            error = 26;
        }
        if ((auxiliary[0] == '#') || (auxiliary[0] == '!'))
            continue;
        sscanf(auxiliary, "%s %s %s %[^\n]", name1, name2, aux2, params);

        std::string nomen1(name1);
        std::string nomen2(name2);
        int LJid1 = 0, LJid2 = 0;
        try
        {
            LJid1 = LJmap.at(nomen1);
            LJid2 = LJmap.at(nomen2);
        }
        catch (const std::out_of_range &)
        {
            error = 27;
        }

        // check vdw type
        if (strstr(aux2, "lj") != NULL)
        {
            current_type = 0; // 0 means LJ (12–6) interactions
        }
        else
        {
            print_warning(0, "Unknown type of disperse (vdw) interactions in .field file!!!\n",
                          "    Given: " + std::string(aux2) + "; expected: lj\n");
            error = 56;
        }

        try
        {
            j = vdwtypesMap.at(current_type);
        }
        catch (const std::out_of_range &)
        {
            switch (current_type)
            {
            case 0:
                system->interMolFields.push_back(new LJsystem(system, LJmap.size()));
                break;
            }
            j = system->interMolFields.size() - 1;
            vdwtypesMap[current_type] = j;
        }

        error = system->interMolFields.at(j)->SetPairParams(LJid1, LJid2, params);

        if (error != 0)
        {
            switch (error)
            {
            case 28:
                print_warning(0, "Wrong value of energy parameter (first parameter, epsilon) of this vdw pair\n",
                              "    Atom 1: " + std::string(nomen1) + ", atom2: " + std::string(nomen2) + ", parameters: " + std::string(params) + "\n");
                break;
            case 29:
                print_warning(0, "Wrong value of sigma (second parameter) of this vdw pair\n",
                              "    Atom 1: " + std::string(nomen1) + ", atom2: " + std::string(nomen2) + ", parameters: " + std::string(params) + "\n");
                break;
            }
            print_warning(0, "Wrong vdw pair specification\n");
            return error;
        }

        if (LJid1 == LJid2)
        {
            diametersMap[LJid1] = system->interMolFields.at(j)->GetDiameter(LJid1);
        }
        i++;
    }

    return 0;
}

// Read one molecular type from .field file
int FieldFile::ReadMolecules(SimulatedSystem *system)
{
    char auxiliary[MAX_COMMENT + 1];
    char aux2[MAX_COMMENT + 1];
    char name[MAX_COMMENT + 1];
    std::string strName;
    int noThisMolecule = 0; // total number of molecules of this type in simulation
    int noThisAtoms = 0;    // no of atoms in one molecule of this type
    double mass = 0.0, charge = 0.0, length = 0.0;
    int repetition = 0;
    int LJid = 0, atom1, atom2;
    int i, j;
    Molecule *mol; // pointer to molecule prototype
    bool finishreached = false;
    int error = 0;

    // species header (molecule name)
    if (fgets(auxiliary, MAX_COMMENT, fieldfile) == nullptr)
        return 46;

    // save molecule name and its hash
    strName.assign(auxiliary);
    strName.pop_back();

    // nummols - number of molecules
    if ((fgets(auxiliary, MAX_COMMENT, fieldfile) != nullptr) && (strstr(auxiliary, "nummols") != 0))
    {
        sscanf(auxiliary, "%s %d", aux2, &noThisMolecule);
    }
    else
    {
        print_warning(0, "Error nummols (in .field file) not defined or not in right place\n",
                      "    Should obtain nummols and number but got: " + std::string(auxiliary) + "\n");
        return 10;
    }

    if ((noThisMolecule < 1) || (noThisMolecule > 1000000))
    {
        print_warning(0, "Wrong number of molecules of one type ('nummols') in .field file\n",
                      "    Given '" + std::string(auxiliary) + "', converted to number: " + std::to_string(noThisMolecule) + "\n");
        return 11;
    }

    system->noMolecules += noThisMolecule;

    if ((fgets(auxiliary, MAX_COMMENT, fieldfile) != nullptr) && (strstr(auxiliary, "atom") != 0))
    {
        sscanf(auxiliary, "%s %d", aux2, &noThisAtoms);
    }
    else
    {
        print_warning(0, "Error number of atoms ('atoms' in .field file) not defined or not in right place\n", "    Should obtain atoms and number but got: " + std::string(auxiliary) + "\n");
        return 12;
    }

    // return if wrong number of atoms
    if ((noThisAtoms < 1) || (noThisAtoms > 1000000))
    {
        print_warning(0, "Wrong number of atoms in one molecule ('atoms') in .field file\n",
                      "    Given '" + std::string(auxiliary) + "', converted to number: " + std::to_string(noThisAtoms) + "\n");
        return 13;
    }

    // molecule initialization
    mol = new Molecule(noThisAtoms, strName);

    // read atoms part of FIELD
    for (i = 0; i < noThisAtoms;)
    {
        if (fgets(auxiliary, MAX_COMMENT, fieldfile) == nullptr)
            return 46;
        if ((auxiliary[0] == '#') || (auxiliary[1] == '!')) // comment line
            continue;
        sscanf(auxiliary, "%s %lf %lf %d", name, &mass, &charge, &repetition);
        charge /= U_CHARGE_E; // conversion of charge to program units...
        strName.assign(name);
        LJmap.insert(std::make_pair(strName, LJmap.size()));
        LJid = LJmap[strName];
        i += repetition;
        if ((repetition < 0) || (repetition > 1000000) || (i > noThisAtoms))
        {
            print_warning(0, "Wrong number of repetition in 'atoms' part in .field file (or mismatch with total 'atoms' number)\n",
                          "    Given: '" + std::string(auxiliary) + "', this atom repetiton is thus: " + std::to_string(repetition) + "\n");
            return 14;
        }
        if (mass <= 0.0)
        {
            print_warning(0, "Cannot handle massles (or even negative mass) atoms\n",
                          "    In .field file obtained '" + std::string(auxiliary) + "', second word must be positive atom mass\n",
                          "    Translated to mass: " + std::to_string(mass) + "\n");
            return 15;
        }
        for (j = 0; j < repetition; j++)
        {
            mol->AddAtom(mass / U_MASS_GMOL, charge, name, LJid, system->GetRsize());
        }
    }

    // constraints, bonds, angles, dihedrals
    while ((fgets(auxiliary, MAX_COMMENT, fieldfile) != NULL) && (feof(fieldfile) == 0))
    {
        if ((auxiliary[0] == '#') || (auxiliary[0] == '!')) // lines starting with '#' or '!' are treated as comments
        {
            continue;
        }
        // fgets(auxiliary, MAX_COMMENT, field_file);
        if (strstr(auxiliary, "constraints") != 0)
        {
            sscanf(auxiliary, "%s %d", aux2, &repetition);
            if ((repetition < 0) || (repetition > mol->noAtoms * mol->noAtoms))
            {
                print_warning(0, "Wrong number of constraints in .field file\n",
                              "    Given: '" + std::string(auxiliary) + "', converted to number of constraints: " + std::to_string(repetition) + "\n");
                return 16;
            }
            mol->InitConstrBond(repetition);
            mol->noInnerDegF -= repetition;                       // subtract number of constraints from the number of DOF of the molecule
            system->noConstrTotal += repetition * noThisMolecule; // total number of constraints sum...

            for (j = 0; j < repetition;)
            {
                if (fgets(auxiliary, MAX_COMMENT, fieldfile) == nullptr)
                    return 46;
                if ((auxiliary[0] == '#') || (auxiliary[0] == '!'))
                    continue;
                sscanf(auxiliary, "%d %d %lf", &atom1, &atom2, &length);
                if ((length <= 0) || (length > 1000))
                {
                    print_warning(0, "Wrong length of constrained bond in .field file\n",
                                  "    Given '" + std::string(auxiliary) + "', length thus: " + std::to_string(length) + "\n");
                    return 17;
                }
                if ((atom1 > noThisAtoms) || (atom2 > noThisAtoms) || (atom1 < 1) || (atom2 < 1) || (atom1 == atom2))
                {
                    print_warning(0, "Wrong atom number in constained bond definition in .field file\n",
                                  "    Given '" + std::string(auxiliary) + "', atoms: " + std::to_string(atom1) + " and " + std::to_string(atom2) + "\n");
                    return 18;
                }
                mol->AddConstrBond(length, atom1 - 1, atom2 - 1);
                j++;
            }
        }
        else if (strstr(auxiliary, "bonds") != 0)
        {
            sscanf(auxiliary, "%s %d", aux2, &repetition);
            if ((repetition < 0) || (repetition > mol->noAtoms * mol->noAtoms))
            {
                print_warning(0, "Wrong number of bonds in .field file\n",
                              "    Given: '" + std::string(auxiliary) + "', converted to number of bonds: " + std::to_string(repetition) + "\n");
                return 19;
            }

            mol->InitBonds(repetition);

            for (j = 0; j < repetition;)
            {
                if (fgets(auxiliary, MAX_COMMENT, fieldfile) == nullptr)
                    return 46;
                if ((auxiliary[0] == '#') || (auxiliary[0] == '!'))
                    continue;
                if ((error = -mol->AddBond(auxiliary, system->GetEnergyUnit())) > 0)
                    return error;
                j++;
            }
        }
        else if (strstr(auxiliary, "angles") != 0)
        {
            sscanf(auxiliary, "%s %d", aux2, &repetition);
            if ((repetition < 0) || (repetition > mol->noAtoms * mol->noAtoms * mol->noAtoms))
            {
                print_warning(0, "Wrong number of angles in .field file\n",
                              "    Given: '" + std::string(auxiliary) + "', converted to number of angles: " + std::to_string(repetition) + "\n");
                return 55;
            }

            mol->InitAngles(repetition);

            for (j = 0; j < repetition;)
            {
                if (fgets(auxiliary, MAX_COMMENT, fieldfile) == nullptr)
                    return 46;
                if ((auxiliary[0] == '#') || (auxiliary[0] == '!'))
                    continue;
                if ((error = -mol->AddAngle(auxiliary, system->GetEnergyUnit())) > 0)
                    return error;
                j++;
            }
        }
        else if (strstr(auxiliary, "dihedrals") != 0)
        {
            sscanf(auxiliary, "%s %d", aux2, &repetition);
            if ((repetition < 0) || (repetition > mol->noAtoms * mol->noAtoms * mol->noAtoms))
            {
                print_warning(0, "Wrong number of dihedrals in .field file\n",
                              "    Given: '" + std::string(auxiliary) + "', converted to number of dihedrals: " + std::to_string(repetition) + "\n");
                return 59;
            }

            mol->InitDihedrals(repetition);

            for (j = 0; j < repetition;)
            {
                if (fgets(auxiliary, MAX_COMMENT, fieldfile) == nullptr)
                    return 46;
                if ((auxiliary[0] == '#') || (auxiliary[0] == '!'))
                    continue;
                if ((error = -mol->AddDihedral(auxiliary, system->GetEnergyUnit())) > 0)
                    return error;
                j++;
            }
        }
        else if (strstr(auxiliary, "finish") != NULL)
        {
            finishreached = true;
            break;
        }
        // other molecule descriptors (dihedrals, etc.) to be added later here...
    }

    if (!finishreached)
    {
        print_warning(0, "'finish' directive in .field file not reached.\n",
                      "    Molecules cannot be initiated properly (or safely).\n");
        return 23;
    }

    // connectivity checks, distance computing...
    if ((error = mol->ConnectivityCheck()) != 0)
    {
        print_warning(0, "Error occured in molecule: " + std::string(mol->name) + " after connectivity checks\n");
        return error;
    }

    // the most important part - add the prototype molecule to the system
    for (i = 0; i < noThisMolecule; i++)
    {
        system->AddMolecule(mol);
    }

    delete mol;
    mol = nullptr;

    // save boundary between different molecular types for 'nfold' directive
    system->AddMolTypeBound();

    return 0;
}

void FieldFile::GetAtomName(int LJid, std::string *name) const
{
    std::map<std::string, int>::const_iterator it;
    std::string result;

    for (it = LJmap.begin(); it != LJmap.end(); it++)
    {
        if (it->second == LJid)
        {
            result = it->first;
        }
    }

    name->assign(result);
}

char DetermineColorForMol(char *name)
{
    char colour;

    // determine color (for .mol file)
    switch (toupper(name[0])) // colors can be customized here
    {
    case 'C':
        colour = 'C';
        break;
    case 'F':
    case 'X':
        colour = 'G';
        break;
    case 'N':
        colour = 'B';
        break;
    case 'O':
        colour = 'R';
        break;
    case 'H':
        colour = 'W';
        break;
    case 'B':
    case 'S':
        colour = 'Y';
        break;
    case 'M':
        colour = 'M';
        break;
    default:
        colour = 'O';
    }
    return colour;
}
