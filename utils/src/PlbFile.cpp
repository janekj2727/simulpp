/*
* .plb file utils 
*  read one plb frame, print interatomic distances etc. (for simul++)
*
*  Author: JJ
*  Date:
*
*/

#include <cstring>
#include <iostream>

#include "file_openclose.hpp"
#include "Vector.hpp"
#include "PlbFile.hpp"
#include "general_utils.hpp"

PlbFile::PlbFile(char *plbfilename)
{
    char rbmode[] = "rb";
    float version_check;
    float N_atoms_f;
    size_t read_size;

    initialized_properly = false;
    frameNo = 0;
    strcpy(filename, plbfilename);

    if ((my_fopen_r(&plbfile, plbfilename, rbmode) != 0) &&
        (strcat(plbfilename, ".plb"), print_warning(2, "Trying to open file with .plb extension added...\n"),
         my_fopen_r(&plbfile, plbfilename, rbmode) != 0))
    {
        print_warning(0, "Cannot open .plb file: " + std::string(plbfilename) + "\n");
        goto endOfConstructor;
    }

    read_size = fread(&N_atoms_f, sizeof(float), 1, plbfile); // read number of sites
    if (read_size > 0)
    {
        N_atoms = (int)N_atoms_f;
    }
    else
    {
        print_warning(0, "Corrupted .plb file: " + std::string(plbfilename) + " cannot read properly\n");
        goto endOfConstructor;
    }

    read_size = fread(&version_check, sizeof(float), 1, plbfile); // read second number (should be -3.0))

    if ((read_size > 0) && (version_check > -2.5))
    {
        print_warning(0, "Old .plb format or corrupted plb file. Cannot read properly...\n");
        N_atoms = 0;
        goto endOfConstructor;
    }

    // allocation of memory for atom coordinates
    positions = new float[N_atoms * 3];

    initialized_properly = true;

endOfConstructor:;
}

PlbFile::~PlbFile()
{
    my_fclose(&plbfile, filename);

    delete[] positions;
    positions = nullptr;
}

/*
PlbFile::PlbFile(const PlbFile& C)
{
    
}

PlbFile& PlbFile::operator=(const PlbFile& C)
{
    return *this;
}
*/

// read one frame
int PlbFile::ReadFrame()
{
    if (N_atoms == 0)
    {
        return 2;
    }
    if (fread(boxsize, sizeof(float), 3, plbfile) != 3)
    {
        if (feof(plbfile))
        {
            return -1; // last frame was the previous
        }
        print_warning(0, "Error during reading of frame no " + std::to_string(frameNo) + " of file: " + std::string(filename) + "!\n");
        return 1;
    }
    if ((int) fread(positions, sizeof(float), 3 * N_atoms, plbfile) != 3 * N_atoms)
    {
        print_warning(0, "Error during reading of frame no " + std::to_string(frameNo) + " of file: " + std::string(filename) + "!\n");
        return 1;
    }

    frameNo++;
    return 0;
}

// return distance between atoms with indices i and j
double PlbFile::AtomDistance(int i, int j)
{
    Vector dist(3);
    int k;

    dist = GetAtomCoordinates(j) - GetAtomCoordinates(i);

    // periodic b.c.
    for (k = 0; k < 3; k++)
    {
        if (dist[k] > (double)boxsize[k] * 0.5)
        {
            dist[k] -= boxsize[k];
        }
        else if (dist[k] < (double)-boxsize[k] * 0.5)
        {
            dist[k] += boxsize[k];
        }
    }

    return dist.CalculateNorm();
}

// get atom coordinates (convert to double)
Vector PlbFile::GetAtomCoordinates(int i)
{
    Vector R(3);
    int j;

    for (j = 0; j < 3; j++)
    {
        R[j] = (double)positions[3 * i + j];
    }

    return R;
}

// get number of atoms
int PlbFile::GetNAtoms() const
{
    return N_atoms;
}