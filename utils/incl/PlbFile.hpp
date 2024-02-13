/*
* .plb file utils 
* read one plb frame, print interatomic distances etc. (for simul++)
*
*  Author: JJ
*  Date: May 2021
*
*/

#ifndef PLBFILEHEADER
#define PLBFILEHEADER

#include<cstdio>

class Vector;

class PlbFile
{
private:
    PlbFile();
    FILE *plbfile; // pointer to .plb file
    float *positions; // 2-dimensional array of positions saved as 1-dimensional
    float boxsize[3]; // box size
    int N_atoms; // number of atoms
    bool initialized_properly;
    char filename[50];
    int frameNo;

public:
    PlbFile(char *plbfilename);
    ~PlbFile();
    // PlbFile(const PlbFile& C);
    // PlbFile& operator=(const PlbFile& C);
    int ReadFrame();
    double AtomDistance(int i, int j);
    Vector GetAtomCoordinates(int i);
    int GetNAtoms() const;

};

#endif