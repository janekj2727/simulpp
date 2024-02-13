/*
 *  plbatomdist – get atom distances from .plb file
 *  part of Simul++
 *
 *  Author: JJ
 *  Date: May 2021
 *
 */

#include <ctype.h>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include "PlbFile.hpp"
#include "general_utils.hpp"

#ifdef PARALLEL
int thread_count = 1;
#endif

int print_help();

int main(int argc, char **argv)
{
    int *atomJ;
    int *atomI;
    char plbfilename[50];
    int Npairs = 0;
    int Natoms = 0;
    int i, k;
    int returnvalue = 0;

    if ((argc < 4) || ((argv[1][0] == '-') && (tolower(argv[1][1]) == 'h')))
    {
        print_help();
        return -1;
    }

    // read commandline arguments
    // 1) filename
    strcpy(plbfilename, argv[1]);
    // open plbfile (use PlbFile class constructor)
    PlbFile plb(plbfilename);
    Natoms = plb.GetNAtoms();

    // number of pairs
    Npairs = (argc - 2) / 2;
    atomI = new int[Npairs];
    atomJ = new int[Npairs];
    k = 2;

    // 2) atom indices of each pair from 0 to Npairs
    for (i = 0; i < Npairs; i++)
    {
        atomI[i] = atoi(argv[k++]);
        if ((atomI[i] >= Natoms) || (atomI[i] < 0))
        {
            print_warning(0, "Requested atom out of range!\n",
                          "Requested: " + std::to_string(atomI[i]) + ", total number of atoms in .plb file: " + std::to_string(Natoms) + "\n");
            return 2;
        }
        atomJ[i] = atoi(argv[k++]);
        if ((atomJ[i] >= Natoms) || (atomJ[i] < 0))
        {
            print_warning(0, "Requested atom out of range!\n",
                          "Requested: " + std::to_string(atomJ[i]) + ", total number of atoms in .plb file: " + std::to_string(Natoms) + "\n");
            return 2;
        }
    }
    std::cout << std::setprecision(12) << std::scientific;

    while ((returnvalue = plb.ReadFrame()) == 0)
    {
        for (i = 0; i < Npairs; i++)
        {
            std::cout << plb.AtomDistance(atomI[i], atomJ[i]) << "  ";
        }
        std::cout << std::endl;
    }

    if (returnvalue != -1)
    {
        print_warning(0, "End of file NOT reached, corrupted .plb file\n");
            return 1;
    }
    else
    {
        std::cout << ".plb file read successfully\n";
    }

    return 0;
}

// print help
int print_help()
{
    std::cout << "plbatomdist - utility to get interatomic distances from .plb file\n";
    std::cout << "originally part of the experimental simul++ package\n";
    std::cout << "Usage:\n";
    std::cout << "    plbatomdist PLBFILE atomI1 atomJ1 [atomI2 atomJ2 ...]\n";
    std::cout << "Where:\n";
    std::cout << "    PLBFILE - .plb file name (only new format with variable box accepted)\n";
    std::cout << "    atomI1 atomJ1 ... - atom indeces\n";
    std::cout << "Returns interatomc distances between atomI1 and atomJ1 (and other specified pairs) in minimum image convention ";
    std::cout << "assuming periodic boundary conditions\n";
    return 0;
}