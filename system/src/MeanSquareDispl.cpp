/*
 *  Mean square displacement class (see header file)
 *
 *  Author: JJ
 *  Date: Sep 2023
 *
 */

#include "MeanSquareDispl.hpp"
#include "Vector.hpp"
#include "Atom.hpp"
#include "Molecule.hpp"
#include <cstring>

MeanSquareDispl::MeanSquareDispl(int noAtoms, std::vector<Molecule> *mol)
{
    int curr_atom_types = -1;
    int i;
    int curr_atom = 0;
    orig_pos = new Vector(noAtoms * 3);
    molecules = mol;
    for (auto it = molecules->begin(); it != molecules->end(); it++)
    {
        mol_index_in_orig.push_back(curr_atom);
        for (i = 0; i < it->noAtoms; i++)
        {
            if (it->atoms[i].LJid > curr_atom_types)
            {
                atom_names.push_back(std::string(it->atoms[i].name));
                no_atoms_of_type.push_back(1);
                curr_atom_types++;
            }
            else
            {
                no_atoms_of_type[it->atoms[i].LJid]++;
            }
            orig_pos->operator()(curr_atom++) = it->atoms[i].GetPosition(0);
            orig_pos->operator()(curr_atom++) = it->atoms[i].GetPosition(1);
            orig_pos->operator()(curr_atom++) = it->atoms[i].GetPosition(2);
        }
    }
    displacements = new Vector(atom_names.size());
}

MeanSquareDispl::~MeanSquareDispl()
{
    if (orig_pos != nullptr)
    {
        delete orig_pos;
    }
    if (displacements != nullptr)
    {
        delete displacements;
    }
    atom_names.clear();
    no_atoms_of_type.clear();
    mol_index_in_orig.clear();
}

MeanSquareDispl::MeanSquareDispl(const MeanSquareDispl &orig)
{
    orig_pos = new Vector(*(orig.orig_pos));
    displacements = new Vector(*(orig.displacements));
    atom_names = orig.atom_names;
    no_atoms_of_type = orig.no_atoms_of_type;
    molecules = orig.molecules;
    mol_index_in_orig = orig.mol_index_in_orig;
}

MeanSquareDispl &MeanSquareDispl::operator=(const MeanSquareDispl &orig)
{
    *orig_pos = *orig.orig_pos;
    *displacements = *displacements;
    atom_names = orig.atom_names;
    no_atoms_of_type = orig.no_atoms_of_type;
    molecules = orig.molecules;
    mol_index_in_orig = orig.mol_index_in_orig;
    return *this;
}

int MeanSquareDispl::PrintMSD(FILE *fcpa, const char *format)
{
    int curr_atom_coord = 0;
    int i, j, id;
    // clear Vector of displacements
    displacements->Clear();
    // for each atom from each molecule...
    for (auto it = molecules->begin(); it != molecules->end(); it++)
    {
        for (i = 0; i < it->noAtoms; i++)
        {
            // ... get LJ id ...
            id = it->atoms[i].LJid;
            // ... and add the square displacement from the original position to displacements vector
            for (j = 0; j < 3; j++)
            {
                displacements->operator()(id) += (it->atoms[i].GetPosition(j) - orig_pos->operator()(curr_atom_coord)) * (it->atoms[i].GetPosition(j) - orig_pos->operator()(curr_atom_coord));
                curr_atom_coord++;
            }
        }
    }
    // for each atom_name (each LJ id) print meand displacement to .cpa file (division by number of atoms of each type is needed)
    for (i = 0; i < (int) atom_names.size(); i++)
    {
        fprintf(fcpa, format, displacements->operator()(i) / (double)no_atoms_of_type[i]);
    }
    // if reached here return 0
    return 0;
}

int MeanSquareDispl::PrintHeader(FILE *fcpa, int col_no, const char *format)
{
    int i;
    int col = col_no - 1; // column number, one must be subtracted to make ++col work few lines below...
    char quantname[15];
    // for each atom_name (each LJ id) print meand displacement to .cpa file (division by number of atoms of each type is needed)
    for (i = 0; i < (int) atom_names.size(); i++)
    {
        // prepare column name (MSD(ATOM_NAME)[AA2]) and print column header in given format to `.cpa` file
        sprintf(quantname, "MSD(%.4s)[AA2]", atom_names[i].c_str());
        fprintf(fcpa, format, ++col, quantname);
    }
    // return last column number
    return col;
}

int MeanSquareDispl::MoveMolecule(int molecule_index, int direction, double distance)
{
    int j;
    // coordinate of first molecule atom + direction -> first index of orig_pos to be changed
    int coord = mol_index_in_orig[molecule_index] + direction;
    for (j = 0; j < molecules->at(molecule_index).noAtoms; j++)
    {
        // original position must be changed
        orig_pos->operator()(coord) += distance;
        coord += 3;                
    }
    return 0;
}

int MeanSquareDispl::GetNumberOfTypes() const
{
    return atom_names.size();
}
