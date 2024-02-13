/*
* .field file reading 
* read .field file and initialize system force fields
* works as simple factory, replace factory methods from SimulatedSystem
*
*  Author: JJ
*  Date: Jun 2022
*
*/

#ifndef FIELDFILEHEADER
#define FIELDFILEHEADER

#define MAX_COMMENT 100
#define DISPERSE_TYPES 1

#include<cstdio>

class Vector;
class SimulatedSystem;

class FieldFile
{
private:
    FieldFile();
    FILE *fieldfile; // pointer to .field file
    char fieldname[MAX_COMMENT]; // name of .field file
    char header[MAX_COMMENT];
    int ReadMolecules(SimulatedSystem *system); // read molecules part of .field file
    int ReadvdW(int novdw, SimulatedSystem *system); // read vdw terms
    int errorstatus;
    std::map<std::string, int> LJmap;
    std::map<int, int> vdwtypesMap;
    std::map<int, int> diametersMap;

public:
    FieldFile(char *field_name);
    ~FieldFile();
    int Read(SimulatedSystem *system);
    int PrintMol(char *mol_name, SimulatedSystem *system); // open, print and close .mol file
    void GetAtomName(int LJid, std::string *name) const;
    // FieldFile(const FieldFile& C);
    // FieldFile& operator=(const FieldFile& C);
};

#endif