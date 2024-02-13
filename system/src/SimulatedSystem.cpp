/*
 * A class to store simulated system for trial simulations of contrained dynamics
 * Author JJ, Date Jan 2021
 */

#include <cstring>
#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <cmath>
#include <vector>

#include "Vector.hpp"
#include "Matrix.hpp"
#include "math_utils.hpp"
#include "LJsystem.hpp"
#include "Atom.hpp"
#include "Molecule.hpp"
#include "ConstraintBond.hpp"
#include "HarmonicBond.hpp"
#include "SimulatedSystem.hpp"
#include "file_openclose.hpp"
#include "units.hpp"
#include "FieldFile.hpp"
#include "MDTable.hpp"
#include "general_utils.hpp"
#include "AbstractPairList.hpp"

#define ATOMR(i, j) ((*atomR)(i, j))
#define REAL_RAND() (((double)rand() + 1.0) / (RAND_MAX + 2.0))
// #define DEBUG

// local function for writing to .cfg
int VarPut(FILE *file, void *v, int size);

// === LOCAL STRUCTURES FOR MACSIMUS OUTPUT ===

typedef double vector3[3];
typedef struct
{
    int key;
    int intval;
    vector3 vecval;
} REC;

typedef struct
{
    int size;       // whole struct in bytes
    int padd;       // padded to 8 bytes
    double logs;    // log of the Nose variable s
    vector3 lambda; // log(box)
    vector3 shape;  // (reserved for box shape)
    vector3 rp[1];  // contiguous array of all r[ns] (and p[ns])
} A_S;

// default constructor
SimulatedSystem::SimulatedSystem(char *simname, char *config_name, char *field_name,
                                 double LJcutoff, int totalProlif, int sizeOfR, int sizeTRVP)
    : Rsize(sizeOfR), TRVPsize(sizeTRVP)
{
    // prepare .mol file name
    char mol_name[MAX_COMMENT];
    int i;
    FieldFile fieldfile(field_name);

    strcpy(mol_name, simname);
    strcat(mol_name, ".mol");

    // null file *
    fcpa = NULL;
    fplb = NULL;
    fhist = NULL;
    facc = NULL;

    // initialize counters to 0
    noAtomsTotal = 0;
    noConstrTotal = 0;
    noDegrOfFred = 0;
    totalMass = 0.0;
    instEpot = 0.0;
    noMolecules = 0;
    Ecorr = 0.0;
    molecularTypesBound.push_back(0);
    error = 0;
    isExternElliptic = false;
    ellipticK[0] = 0.0;
    ellipticK[1] = 0.0;
    ellipticK[2] = 0.0;
    Eextra = 0.0;
    Xi = 0.0;
    Xivel = 0.0;
    Lambda = 0.0;
    Lambdavel = 0.0;
    cutoffCorrections = true;
    pressureNfCorrection = true;
    noDOFforPressure = 0;
    maxShakeIter = 0;
    vdWcutoff = LJcutoff;
    omegaShake = 1.0;
    Pext = -999999.9;
    extendedDOFs = nullptr;
    pairlist = nullptr;
    noPairs = 0;
    memoryForR = nullptr;
    msd = nullptr;

    for (i = 0; i < INTRA_MOL_TYPES; i++)
    {
        Epotintramol[i] = 0.0;
        virintramol[i] = 0.0;
    }

    for (i = 0; i < INTER_MOL_TYPES; i++)
    {
        Epotintermol[i] = 0.0;
        virintermol[i] = 0.0;
    }

    error = fieldfile.Read(this);
    if (error != 0)
    {
        goto endOfConstructor;
    }

    // allocate Matrices for Atom::R
    error = AllocateR();
    if (error != 0)
    {
        goto endOfConstructor;
    }

    error = ReadConfig(config_name, totalProlif);
    if (error != 0)
    {
        goto endOfConstructor;
    }

    if (fieldfile.PrintMol(mol_name, this) != 0)
    {
        print_warning(0, "Error while printing .mol file, .mol file will be corrupted\n");
        // not a fatal error...
    }

endOfConstructor:;
}

// copy constructor (not needed) – not complete (bonds, constr. bonds missing...)
SimulatedSystem::SimulatedSystem(SimulatedSystem &simSystem, bool withfields)
    : Rsize(simSystem.Rsize), TRVPsize(simSystem.TRVPsize)
{
    int i;
    std::vector<AbstractInterMolField *>::iterator it;

    noMolecules = simSystem.noMolecules;
    for (i = 0; i < noMolecules; i++)
    {
        molecules.push_back(simSystem.molecules.at(i));
    }
    for (i = 0; i < 3; i++)
    {
        boxSize[i] = simSystem.GetBox(i);
    }
    noAtomsTotal = simSystem.noAtomsTotal;
    noConstrTotal = simSystem.noConstrTotal;
    noDegrOfFred = simSystem.noDegrOfFred;
    noDOFforPressure = simSystem.noDOFforPressure;
    totalMass = simSystem.totalMass;
    instEpot = simSystem.instEpot;
    Ecorr = simSystem.Ecorr;
    boundaryCond = simSystem.boundaryCond;
    fplb = NULL;
    fcpa = NULL;
    fhist = NULL;
    facc = NULL;
    molecularTypesBound = simSystem.molecularTypesBound;
    error = simSystem.error;
    isExternElliptic = simSystem.isExternElliptic;
    ellipticK[0] = simSystem.ellipticK[0];
    ellipticK[1] = simSystem.ellipticK[1];
    ellipticK[2] = simSystem.ellipticK[2];
    Eextra = simSystem.Eextra;
    virconstr = simSystem.virconstr;
    Xi = simSystem.Xi;
    Xivel = simSystem.Xivel;
    Lambda = simSystem.Lambda;
    Lambdavel = simSystem.Lambdavel;
    cutoffCorrections = simSystem.cutoffCorrections;
    pressureNfCorrection = simSystem.pressureNfCorrection;
    maxShakeIter = 0;
    vdWcutoff = simSystem.vdWcutoff;
    omegaShake = simSystem.omegaShake;
    Pext = simSystem.Pext;
    extendedDOFs = nullptr;
    pairlist = nullptr;
    noPairs = simSystem.noPairs;
    memoryForR = nullptr;
    msd = nullptr;

    for (i = 0; i < INTRA_MOL_TYPES; i++)
    {
        Epotintramol[i] = simSystem.Epotintramol[i];
        virintramol[i] = simSystem.virintramol[i];
    }

    for (i = 0; i < INTER_MOL_TYPES; i++)
    {
        Epotintermol[i] = simSystem.Epotintermol[i];
        virintermol[i] = simSystem.virintermol[i];
    }

    if (withfields)
    {
        for (it = simSystem.interMolFields.begin(); it != simSystem.interMolFields.end(); it++)
        {
            interMolFields.push_back((*it)->copy(this));
        }
    }
}

// merging constructor (only partial LJsystem and without bonds,... initialization)
SimulatedSystem::SimulatedSystem(SimulatedSystem &simSystem1, SimulatedSystem &simSystem2, int direction)
    : Rsize(simSystem1.Rsize), TRVPsize(simSystem1.TRVPsize)
{
    int i;
    std::vector<int> positionsCoord;
    positionsCoord.push_back(0);
    std::vector<int> molecularTypesVec;
    std::set<int> molecularTypesSet;
    std::vector<int>::iterator it;

    // move systems in order not to overlap
    simSystem1.MoveSystem(-simSystem2.GetBox(direction) * 0.5, direction, positionsCoord);
    simSystem2.MoveSystem(simSystem1.GetBox(direction) * 0.5, direction, positionsCoord);

    noMolecules = simSystem1.noMolecules + simSystem2.noMolecules;
    // molecules must be reordered according to their names
    // create set of unique names (their hashed values) and store them in vector
    for (i = 0; i < simSystem1.noMolecules; i++)
    {
        if (molecularTypesSet.insert(simSystem1.molecules[i].hashname).second)
        {
            molecularTypesVec.push_back(simSystem1.molecules[i].hashname);
        }
    }
    for (i = 0; i < simSystem2.noMolecules; i++)
    {
        if (molecularTypesSet.insert(simSystem2.molecules[i].hashname).second)
        {
            molecularTypesVec.push_back(simSystem2.molecules[i].hashname);
        }
    }

    // for each value in set
    for (it = molecularTypesVec.begin(); it != molecularTypesVec.end(); it++)
    {
        for (i = 0; i < simSystem1.noMolecules; i++)
        {
            if (simSystem1.molecules[i].hashname == *it)
                molecules.push_back(simSystem1.molecules.at(i));
        }
        for (i = 0; i < simSystem2.noMolecules; i++)
        {
            if (simSystem1.molecules[i].hashname == *it)
                molecules.push_back(simSystem2.molecules.at(i));
        }
    }

    // box size of the new system
    for (i = 0; i < 3; i++)
    {
        if (i == direction)
        {
            boxSize[i] = simSystem1.GetBox(i) + simSystem2.GetBox(i);
        }
        else
        {
            boxSize[i] = simSystem1.GetBox(i);
        }
    }

    noAtomsTotal = simSystem1.noAtomsTotal + simSystem2.noAtomsTotal;
    noConstrTotal = simSystem1.noConstrTotal + simSystem2.noConstrTotal;
    noDegrOfFred = simSystem1.noDegrOfFred + simSystem2.noDegrOfFred;
    totalMass = simSystem1.totalMass + simSystem2.totalMass;
    instEpot = 0.0;
    Ecorr = 0.0;
    boundaryCond = simSystem1.boundaryCond;
    isExternElliptic = false;
    ellipticK[0] = 0.0;
    ellipticK[1] = 0.0;
    ellipticK[2] = 0.0;
    Eextra = 0.0;
    fplb = NULL;
    fcpa = NULL;
    fhist = NULL;
    facc = NULL;
    Eextra = 0.0;
    virconstr = 0.0;
    Xi = 0.0;
    Xivel = 0.0;
    Lambda = 0.0;
    Lambdavel = 0.0;
    maxShakeIter = 0;
    vdWcutoff = fmin(simSystem1.vdWcutoff, simSystem2.vdWcutoff);
    omegaShake = 1.0;
    Pext = -999999.9;
    extendedDOFs = nullptr;
    pairlist = nullptr;
    noPairs = 0;
    memoryForR = nullptr;
    msd = nullptr;

    for (i = 0; i < INTRA_MOL_TYPES; i++)
    {
        Epotintramol[i] = 0.0;
        virintramol[i] = 0.0;
    }
    for (i = 0; i < INTER_MOL_TYPES; i++)
    {
        Epotintermol[i] = 0.0;
        virintermol[i] = 0.0;
    }
}

// default destructor
SimulatedSystem::~SimulatedSystem()
{
    char appendix[10];
    strcpy(appendix, ".plb");
    std::vector<AbstractInterMolField *>::iterator it;

    if ((fplb != NULL) && (my_fclose(&fplb, appendix) != 0))
    {
        print_warning(0, "Cannot close .plb file\n");
    }
    strcpy(appendix, ".cpa");
    if ((fcpa != NULL) && (my_fclose(&fcpa, appendix) != 0))
    {
        print_warning(0, "Cannot close .cpa file\n");
    }
    strcpy(appendix, ".history");
    if ((fhist != NULL) && (my_fclose(&fhist, appendix) != 0))
    {
        print_warning(0, "Cannot close .history file\n");
    }
    strcpy(appendix, ".acc");
    if ((facc != NULL) && (my_fclose(&facc, appendix) != 0))
    {
        print_warning(0, "Cannot close .acc file\n");
    }

    for (it = interMolFields.begin(); it != interMolFields.end(); it++)
    {
        delete (*it);
    }
    interMolFields.clear();

    if (extendedDOFs != nullptr)
    {
        delete extendedDOFs;
    }
    if (pairlist != nullptr)
    {
        delete pairlist;
    }
    if (memoryForR != nullptr)
    {
        DeallocateR();
    }
    if (msd != nullptr)
    {
        delete msd;
    }
}

// calculate intermolecular forces (LJ, charges)
// return intermolecular energy
double SimulatedSystem::CalculateInterMolForces()
{
    int i;
    std::vector<AbstractInterMolField *>::const_iterator it;
    double Epot = 0.0;

    // recalculate interatomic distances
    pairlist->Refresh();

    for (i = 0; i < INTER_MOL_TYPES; i++)
    {
        Epotintermol[i] = 0.0;
        virintermol[i] = 0.0;
    }
    for (it = interMolFields.begin(); it != interMolFields.end(); it++)
    {
        Epot += (*it)->CalculateForces();
    }

    if (isExternElliptic)
    {
        CalculateExternalForces(); // if needed virial should be computed...
    }

    return Epot;
}

// calculate intermolecular potential energy
double SimulatedSystem::CalculateInterMolEpot() const
{
    double Eintermol = 0.0;
    std::vector<AbstractInterMolField *>::const_iterator it;

    // recalculate interatomic distances
    pairlist->Refresh();

    for (it = interMolFields.begin(); it != interMolFields.end(); it++)
    {
        Eintermol += (*it)->CalculateEpot();
    }

    if (isExternElliptic)
    {
        Eintermol += CalculateExternalEpot(); // if needed virial should be computed...
    }

    return Eintermol;
}

// initialization
int SimulatedSystem::AddMolecule(Molecule *mol)
{
    molecules.push_back(*mol);
    noAtomsTotal += mol->noAtoms;
    totalMass += mol->molMass;
    noDegrOfFred += mol->noInnerDegF + 3;
    noDOFforPressure += mol->noInnerDegF + 3;
    return molecules.size();
}

// calculate all physical forces (except fictitious constrained forces)
double SimulatedSystem::CalculateForces()
{
    double Epot = 0.0;
    int i;

    for (i = 0; i < INTRA_MOL_TYPES; i++)
    {
        Epotintramol[i] = 0.0;
        virintramol[i] = 0.0;
    }

#ifdef PARALLEL
#pragma omp parallel reduction(+ : Epot) num_threads(thread_count)
    {
        double localEpotintramol[INTRA_MOL_TYPES];
        double localvirintramol[INTRA_MOL_TYPES];
        int j;
        for (j = 0; j < INTRA_MOL_TYPES; j++)
        {
            localEpotintramol[j] = 0.0;
            localvirintramol[j] = 0.0;
        }

#pragma omp for schedule(static, 10)
#endif
        for (i = 0; i < noMolecules; i++)
        {
            molecules[i].ZeroForces();
#ifdef PARALLEL
            Epot += molecules[i].CalculateIntraMolForces(localEpotintramol, localvirintramol);
#else
        Epot += molecules[i].CalculateIntraMolForces(Epotintramol, virintramol);
#endif
        }
#ifdef PARALLEL
#pragma omp critical
        {
            for (j = 0; j < INTRA_MOL_TYPES; j++)
            {
                Epotintramol[j] += localEpotintramol[j];
                virintramol[j] += localvirintramol[j];
            }
        }
    }
#endif

    // now in Epot only EpotIntra...
    Epot += CalculateInterMolForces(); // parallelization inside

    // now EpotInter is added and EpotTotal is returned
    instEpot = Epot;
    virconstr = 0.0;

    return Epot;
}

// print configuration to .plb
int SimulatedSystem::PrintPlayback()
{
    int i, j;
    Matrix *atomR;
    float x, y, z, Lx, Ly, Lz;

    // print configuration to .plb
    Lx = (float)boxSize[0];
    Ly = (float)boxSize[1];
    Lz = (float)boxSize[2];
    fwrite(&Lx, sizeof(float), 1, fplb);
    fwrite(&Ly, sizeof(float), 1, fplb);
    fwrite(&Lz, sizeof(float), 1, fplb);

    for (i = 0; i < noMolecules; i++)
    {
        for (j = 0; j < molecules[i].noAtoms; j++)
        {
            atomR = molecules[i].atoms[j].R;
            x = (float)ATOMR(0, 0) + 0.5 * Lx;
            y = (float)ATOMR(0, 1) + 0.5 * Ly;
            z = (float)ATOMR(0, 2) + 0.5 * Lz;
            fwrite(&x, sizeof(float), 1, fplb);
            fwrite(&y, sizeof(float), 1, fplb);
            fwrite(&z, sizeof(float), 1, fplb);
        }
    }

    return 0;
}

// Read simname.config file to get atom coordinates
int SimulatedSystem::ReadConfig(char *config_name, int totalProlif)
{
    char auxiliary[MAX_COMMENT + 1];
    double auxdouble;
    int auxint;
    char auxname[MAX_COMMENT];
    std::string auxstring;
    FILE *fconfig;                    // .config file
    int levcfg, imcon, noAtomsConfig; // level of info (positions, positions+velocities,...), boundary conditions, total no of atoms
    int normalLevCfg, TRVPlev;        // level of info (standard – DL POLY type), level of TRVP (number of previous velocities stored - 1)
    int extDOFno;                     // number of extended DOFs stored in place of molecule 0
    double x, y, z;
    int i, j, k, l, m;
    int noThisMolecule;
    char read[] = "r";

    // Field file opening...
    if (my_fopen_r(&fconfig, config_name, read) != 0)
    {
        print_warning(0, "Config file ", std::string(config_name), " cannot be opened!\n");
        return 33;
    }
    // protocol part:
    // printf("Reading .config file: %s\n", config_name);

    // reading and printing header
    if (fgets(auxiliary, MAX_COMMENT, fconfig) == nullptr)
        return 47;
    // protocol part:
    // printf("Header .config: %s\n", auxiliary);

    // reading .config level, imcon and noAtomsConfig
    if (fgets(auxiliary, MAX_COMMENT, fconfig) == nullptr)
        return 47;
    sscanf(auxiliary, "%d %d %d", &levcfg, &imcon, &noAtomsConfig);
    // levcfg = xy
    // y - 0 positions only, 1 - positions + velocities, ...
    // x - trvp size
    // imcon - 0 free, 1 - cubic, 2 - rectangular box
    boundaryCond = imcon;
    if (boundaryCond == 2)
    {
        noDegrOfFred -= 3; // subtract 3 from total number of degrees of freedom for periodic b. c.
    }
    // these values determines the content of the .config file that is being read...
    // number of extended degrees of freedom in .config file
    extDOFno = (int)floor(levcfg / 100.0);
    levcfg -= 100 * extDOFno;
    // TRVP level in .config file
    TRVPlev = (int)floor(levcfg / 10.0);
    // normal cfg level in .config file (0 .. positions only, 1 .. positions + velocities, ...)
    normalLevCfg = levcfg - 10 * TRVPlev;
    // check total no Atoms
    if (noAtomsConfig * totalProlif != noAtomsTotal)
    {
        print_warning(0, "Bad number of atoms in " + std::string(config_name) + ": " + std::to_string(noAtomsConfig) + "(total proliferation " + std::to_string(totalProlif) + ")\n", "    mismatch with .field file: " + std::to_string(noAtomsTotal) + "\n");
        return 34;
    }

    // read box sizes (if not free b.c.) else boxSize = 0 (only for .plb file)
    if (imcon != 0)
    {
        if (fgets(auxiliary, MAX_COMMENT, fconfig) == nullptr)
            return 47;
        sscanf(auxiliary, "%lf %lf %lf", &boxSize[0], &auxdouble, &auxdouble);
        if (fgets(auxiliary, MAX_COMMENT, fconfig) == nullptr)
            return 47;
        sscanf(auxiliary, "%lf %lf %lf", &auxdouble, &boxSize[1], &auxdouble);
        if (fgets(auxiliary, MAX_COMMENT, fconfig) == nullptr)
            return 47;
        sscanf(auxiliary, "%lf %lf %lf", &auxdouble, &auxdouble, &boxSize[2]);
    }
    else // imcon = 0 .. free boundary condition
    {
        for (i = 0; i < 3; i++)
        {
            boxSize[i] = 0.0;
        }
    }

    // if extended DOFs stored, read them first
    if (extDOFno > 0)
    {
        // allocate system->extendedDOFs (deallocation in integrators' InitXi or InitLambda)
        // if not deallocated in integrators, deallocation occurs in ~SimulatedSystem...
        extendedDOFs = new Matrix(Rsize, extDOFno);
        // read header and check DOF name
        if ((fgets(auxiliary, MAX_COMMENT, fconfig) != NULL) && (feof(fconfig) == 0))
        {
            if ((strstr(auxiliary, "Xi") == NULL) && (strstr(auxiliary, "Lambda") == NULL))
            {
                print_warning(0, "Error unexpected name of extra degree of freedom in .config file.\n",
                              "    Expected: Xi or Lambda, but recieved " + std::string(auxname) + "\n");
                return 36;
            }
            else
            {
                strcpy(extDOFnames, auxiliary);
            }
        }
        else
        {
            print_warning(0, "End of .config file reached sooner then expected.\n");
            return 35;
        }
        // read positions, velocities, forces, ... (according to levcfg (normalLevCfg),
        // but save only those needed (max to Rsize - TRVPsize - 1)
        for (l = 0; l < normalLevCfg + 1; l++)
        {
            if ((fgets(auxiliary, MAX_COMMENT, fconfig) != NULL) && (feof(fconfig) == 0))
            {
                if (l < Rsize - ((TRVPsize == 0) ? 0 : (TRVPsize + 1))) // if there is a place in Rsize to store these values (in DOFs there is always a place, but not so in the final destination in integrator)
                {
                    m = sscanf(auxiliary, "%lf %lf %lf", &x, &y, &z); // now max 3 DOFs can be read
                    if (m > 0)
                    {
                        extendedDOFs->operator()(l, 0) = x;
                    }
                    if ((extDOFno > 1) && (m > 1))
                    {
                        extendedDOFs->operator()(l, 1) = y;
                    }
                    if ((extDOFno > 2) && (m > 2))
                    {
                        extendedDOFs->operator()(l, 2) = z;
                    }
                }
            }
            else
            {
                print_warning(0, "End of .config file reached sooner then expected.\n");
                return 35;
            }
        }
        // read stored values of past velocities for TRVP
        if (TRVPlev != 0)
        {
            for (l = 0; l < TRVPlev + 1; l++)
            {

                if ((fgets(auxiliary, MAX_COMMENT, fconfig) != NULL) && (feof(fconfig) == 0))
                {
                    if (l < ((TRVPsize == 0) ? 0 : (TRVPsize + 1))) // if there is a place in Rsize to store these values
                    {
                        m = sscanf(auxiliary, "%lf %lf %lf", &x, &y, &z); // now max 3 DOFs can be read
                        if (m > 0)
                        {
                            extendedDOFs->operator()(l + Rsize - TRVPsize - 1, 0) = x;
                        }
                        if ((extDOFno > 1) && (m > 1))
                        {
                            extendedDOFs->operator()(l + Rsize - TRVPsize - 1, 1) = y;
                        }
                        if ((extDOFno > 2) && (m > 2))
                        {
                            extendedDOFs->operator()(l + Rsize - TRVPsize - 1, 2) = z;
                        }
                    }
                }
                else
                {
                    print_warning(0, "End of .config file reached sooner then expected.\n");
                    return 35;
                }
            }
        }
    }

    // read atoms position (and velocities and ...)
    // little complicated due to correct nfold handling
    // first cycle through molecular types
    for (k = 1; k < (int)molecularTypesBound.size(); k++)
    {
        noThisMolecule = (int)round((molecularTypesBound.at(k) - molecularTypesBound.at(k - 1)) / totalProlif);
        // read positions of this type (cycle through molecules of this type)
        for (i = molecularTypesBound.at(k - 1); i < molecularTypesBound.at(k - 1) + noThisMolecule; i++)
        {
            // cycle through atoms in molecule
            for (j = 0; j < molecules[i].noAtoms; j++)
            {
                // read atom header and check atom name
                if ((fgets(auxiliary, MAX_COMMENT, fconfig) != NULL) && (feof(fconfig) == 0))
                {
                    sscanf(auxiliary, "%s %d", auxname, &auxint);
                    if (strstr(auxname, molecules[i].atoms[j].name) == NULL)
                    {
                        print_warning(0, "Error wrong atom type in .config file (molecule " + std::to_string(i + 1) + ", atom " + std::to_string(j + 1) + ")\n", "    Expected: " + std::string(molecules[i].atoms[j].name) + ", but recieved " + std::string(auxname) + "\n");
                        return 36;
                    }
                }
                else
                {
                    print_warning(0, "End of .config file reached sooner then expected.\n");
                    return 35;
                }
                // read positions, velocities, forces, ... (according to levcfg (normalLevCfg),
                // but save only those needed (max to Rsize - TRVPsize - 1)
                for (l = 0; l < normalLevCfg + 1; l++)
                {
                    if ((fgets(auxiliary, MAX_COMMENT, fconfig) != NULL) && (feof(fconfig) == 0))
                    {
                        if (l < Rsize - ((TRVPsize == 0) ? 0 : (TRVPsize + 1))) // if there is a place in Rsize to store these values
                        {
                            if (sscanf(auxiliary, "%lf %lf %lf", &x, &y, &z) != 3)
                            {
                                print_warning(0, "Error in .config file. Corrupted atom coordinates\n.");
                                return 35;
                            }
                            molecules[i].atoms[j].R->operator()(l, 0) = x;
                            molecules[i].atoms[j].R->operator()(l, 1) = y;
                            molecules[i].atoms[j].R->operator()(l, 2) = z;
                        }
                    }
                    else
                    {
                        print_warning(0, "End of .config file reached sooner then expected.\n");
                        return 35;
                    }
                }
                // read stored values of past velocities for TRVP
                if (TRVPlev != 0)
                {
                    for (l = 0; l < TRVPlev + 1; l++)
                    {

                        if ((fgets(auxiliary, MAX_COMMENT, fconfig) != NULL) && (feof(fconfig) == 0))
                        {
                            if (l < ((TRVPsize == 0) ? 0 : (TRVPsize + 1))) // if there is a place in Rsize to store these values
                            {
                                sscanf(auxiliary, "%lf %lf %lf", &x, &y, &z);
                                molecules[i].atoms[j].R->operator()(l + Rsize - TRVPsize - 1, 0) = x;
                                molecules[i].atoms[j].R->operator()(l + Rsize - TRVPsize - 1, 1) = y;
                                molecules[i].atoms[j].R->operator()(l + Rsize - TRVPsize - 1, 2) = z;
                            }
                        }
                        else
                        {
                            print_warning(0, "End of .config file reached sooner then expected.\n");
                            return 35;
                        }
                    }
                }
            }
        }
    }
    // .config file closing...
    if (my_fclose(&fconfig, config_name) != 0)
    {
        print_warning(0, "Failed closing " + std::string(config_name) + "...\n");
        return 37;
    }

    return 0;
}

// get Box size element
double SimulatedSystem::GetBox(int coord) const
{
    return boxSize[coord];
}

// set Box size directly
int SimulatedSystem::SetBox(double size, int coord)
{
    if ((size < 0.0) || (size > 100000))
    {
        print_warning(1, "Trying to set box size to unphysical number: " + std::to_string(size) + "\n",
                      "    Not allowing this operation!\n");
    }
    else
    {
        boxSize[coord] = size;
    }

    return 0;
}

// initialization of .plb file
int SimulatedSystem::InitPlbFile(char *plb_name, bool append)
{
    float N_atom_f;
    float var_box;
    char writebinary[] = "wb";
    char read[] = "rb";

    // if append then append binary and check existence
    if (append)
    {
        writebinary[0] = 'a';
        if (my_fopen_r(&fplb, plb_name, read) != 0)
        {
            print_warning(0, "Cannot open " + std::string(plb_name) + " file to append playback\n", "    File does not exist?\n");
            return 38;
        }
        else
        {
            my_fclose(&fplb, plb_name);
        }
    }

    // open .plb file
    if (my_fopen_w(&fplb, plb_name, writebinary) != 0)
    {
        print_warning(0, "Cannot open " + std::string(plb_name) + " file to write\n");
        return 38;
    }

    // write header to .plb file in not append
    if (!append)
    {
        /* Write first line to .plb file (binary file, first line contains
         *  two float numbers - first being the number of atoms and second the box size
         *  (-3) in case of variable box size, which is our case, because this is the de-
         *  fault value for the newer versions of MACSIMUS) */
        N_atom_f = (float)noAtomsTotal;
        var_box = -3.0;
        fwrite(&N_atom_f, sizeof(float), 1, fplb);
        fwrite(&var_box, sizeof(float), 1, fplb);
    }

    return 0;
}

// initialization of .history file
int SimulatedSystem::InitHistoryFile(char *hist_name, int no_frames, int hist_level, bool append)
{
    char write[] = "w";
    char read[] = "r";

    // if append then append instead of write but first check existence
    if (append)
    {
        write[0] = 'a';
        if (my_fopen_r(&fhist, hist_name, read) != 0)
        {
            print_warning(0, "Cannot open " + std::string(hist_name) + " file to append history\n", "    File does not exist?\n");
            return 41;
        }
        else
        {
            my_fclose(&fhist, hist_name);
        }
    }

    // open .history file
    if (my_fopen_w(&fhist, hist_name, write) != 0)
    {
        print_warning(0, "Cannot open " + std::string(hist_name) + " file to write\n");
        return 41;
    }

    if (!append)
    {
        // write header to .history file
        /* Write first line to .history file – comment */
        fprintf(fhist, "Trajectory of simulation. View to check velocities or visualize by VMD.\n");
        /* Write second line to .history file – keytraj, imcon, megatm, frame, records
         *  information level (0 = positions only, 1 = positions + velocities,...
         *  periodic image convention (0 (free), 1 (cubic), 2 (periodic)),
         *  megatm – number of atoms in system – totalNoAtoms,
         *  frame – number of frames in .history file – no_frames
         *  records – number of lines in .history file (2 for header, each frame: 4 for header, (2+keytraj) for each atom */

        fprintf(fhist, "%d %d %d %d %d\n", hist_level, boundaryCond, noAtomsTotal, no_frames,
                2 + no_frames * (4 + noAtomsTotal * (2 + hist_level)));
    }
    return 0;
}

// apply periodic boundary conditions
int SimulatedSystem::ApplyPeriodicBC(std::vector<int> &Rrows)
{
    int moved = 0;

#ifdef PARALLEL
#pragma omp parallel num_threads(thread_count) reduction(+ : moved)
    {
#endif
        double COM[3];
        int i, m;
#ifdef PARALLEL
#pragma omp for
#endif
        for (i = 0; i < noMolecules; i++)
        {
            molecules[i].CalculateCOM(COM);

            for (m = 0; m < 3; m++)
            {
                if (COM[m] > 0.53 * boxSize[m])
                {
                    molecules[i].Move(-boxSize[m], m, Rrows);
                    pairlist->MoveMolecule(i, m, -boxSize[m]);
                    if (msd != nullptr)
                        msd->MoveMolecule(i, m, -boxSize[m]);
                    moved++;
                }
                else if (COM[m] < -0.53 * boxSize[m])
                {
                    molecules[i].Move(boxSize[m], m, Rrows);
                    pairlist->MoveMolecule(i, m, boxSize[m]);
                    if (msd != nullptr)
                        msd->MoveMolecule(i, m, boxSize[m]);
                    moved++;
                }
            }
        }
#ifdef PARALLEL
    }
#endif
    return moved;
}

// print .config file to enable restart
int SimulatedSystem::PrintConfig(int configLevel, char *simname, double h)
{
    // implement print of .config file with configLevel 0--n (for gear...)
    // if configLevel > 10, then trvp size = floor(configLevel/10)
    FILE *fconfig;
    char config_name[MAX_COMMENT];
    int i, j, k;
    int curratom = 1;
    char write[] = "w";
    Matrix *atomR;
    int normalConfigLevel, trvpLevel, trvpStart, noExtDOFs;

    // .config file name
    if (strstr(simname, ".dump") == nullptr) // if called from dump, it has name already set
    {
        strcpy(config_name, simname);
        strcat(config_name, ".configfinal");
    }
    else
    {
        strcpy(config_name, simname);
    }

    // open .config file
    if (my_fopen_w(&fconfig, config_name, write) != 0)
    {
        print_warning(1, "Cannot open " + std::string(config_name) + " file to write\n", "    .config file will not be made!!!\n");
        // this is not a fatal error, thus continue simulation (or its closing) without it
        return -1;
    }

    // parse configLevel
    noExtDOFs = (int)floor(configLevel / 100.0);
    trvpLevel = (int)floor((configLevel - 100.0 * noExtDOFs) / 10.0);
    normalConfigLevel = configLevel - 10 * trvpLevel - 100 * noExtDOFs;
    trvpStart = Rsize - trvpLevel - 1;

    // header of .config file, hopefully compatible with DL_POLY, potential energy instead of total
    fprintf(fconfig, "Final configuration of %s\n", simname);
    fprintf(fconfig, "%9d %9d %9d %19f\n", configLevel, boundaryCond, noAtomsTotal, instEpot);

    // cell geometry to CONFIG if not free b.c.
    if (boundaryCond != 0)
    {
        fprintf(fconfig, "%19.12E %19.12E %19.12E \n", boxSize[0], 0.0, 0.0);
        fprintf(fconfig, "%19.12E %19.12E %19.12E \n", 0.0, boxSize[1], 0.0);
        fprintf(fconfig, "%19.12E %19.12E %19.12E \n", 0.0, 0.0, boxSize[2]);
    }

    // if any extended DOFs should be printed, print them now
    if (noExtDOFs > 0)
    {
        fprintf(fconfig, "%-17s%1d\n", extDOFnames, 0);
        for (k = 0; k <= normalConfigLevel; k++)
        {
            for (i = 0; i < noExtDOFs; i++)
            {
                fprintf(fconfig, "%19.12E ", extendedDOFs->operator()(k, i) * fact(k) / pow(h, (double)k));
            }
            fprintf(fconfig, "\n");
        }
        if (trvpLevel > 0)
        {
            for (k = 0; k <= trvpLevel; k++)
            {
                for (i = 0; i < noExtDOFs; i++)
                {
                    fprintf(fconfig, "%19.12E ", extendedDOFs->operator()(k + trvpStart, i) / h);
                }
                fprintf(fconfig, "\n");
            }
        }
    }

    // print atom positions, velocities etc.
    for (i = 0; i < noMolecules; i++)
    {
        for (j = 0; j < molecules[i].noAtoms; j++)
        {
            fprintf(fconfig, "%-8s%10d\n", molecules[i].atoms[j].name, curratom++);
            atomR = molecules[i].atoms[j].R;
            for (k = 0; k <= normalConfigLevel; k++)
            {
                fprintf(fconfig, "%19.12E %19.12E %19.12E \n",
                        ATOMR(k, 0) * fact(k) / pow(h, (double)k),
                        ATOMR(k, 1) * fact(k) / pow(h, (double)k),
                        ATOMR(k, 2) * fact(k) / pow(h, (double)k));
            }
            if (trvpLevel > 0)
            {
                for (k = 0; k <= trvpLevel; k++)
                {
                    fprintf(fconfig, "%19.12E %19.12E %19.12E \n",
                            ATOMR(k + trvpStart, 0) / h,
                            ATOMR(k + trvpStart, 1) / h,
                            ATOMR(k + trvpStart, 2) / h);
                }
            }
        }
    }

    // close .config file
    if (my_fclose(&fconfig, config_name) != 0)
    {
        print_warning(1, "Error during .config file closing\n", "    file can be corrupted!!!\n");
        // not a fatal error
        return -2;
    }

    return 0;
}

// print configuration to .history file
int SimulatedSystem::PrintHistory(double h, int hist_level, double curr_time)
{
    // implement print of 1 frame to .history file with configLevel = hist_level
    // designed for Gear-family integrators, not showing actual acceleration etc. for Verlet
    int i, j, k;
    int curratom = 1;
    Matrix *atomR;
    int step_no = (int)round(curr_time / h);

    // timestep and frame info timestep
    fprintf(fhist, "timestep %d %d %d %d %f %f\n", step_no, noAtomsTotal,
            hist_level, boundaryCond, h, curr_time);

    // cell geometry
    fprintf(fhist, "%19f %19f %19f \n", boxSize[0], 0.0, 0.0);
    fprintf(fhist, "%19f %19f %19f \n", 0.0, boxSize[1], 0.0);
    fprintf(fhist, "%19f %19f %19f \n", 0.0, 0.0, boxSize[2]);

    // print atom positions, velocities etc.
    for (i = 0; i < noMolecules; i++)
    {
        for (j = 0; j < molecules[i].noAtoms; j++)
        {
            fprintf(fhist, "%-8s%10d %f %f %f\n", molecules[i].atoms[j].name, curratom++, molecules[i].atoms[j].mass * U_MASS_GMOL, molecules[i].atoms[j].charge * U_CHARGE_E, 0.0);
            // last number should be displacement from position in t=0, but here we don't provide this information...
            // last number provided due to potential VMD support (file is recognized as DL_POLY 4 HISTORY...)
            atomR = molecules[i].atoms[j].R;
            for (k = 0; k <= hist_level; k++) // optimized for Gear-type integrators, for Verlet-type integrators scaling factors are wierd
            {
                fprintf(fhist, "%19.12E %19.12E %19.12E \n",
                        ATOMR(k, 0) * fact(k) / pow(h, (double)k),
                        ATOMR(k, 1) * fact(k) / pow(h, (double)k),
                        ATOMR(k, 2) * fact(k) / pow(h, (double)k));
            }
        }
    }

    return 0;
}

// proliferate simulated system in x y z directions
int SimulatedSystem::Proliferate(int nfoldx, int nfoldy, int nfoldz)
{
    std::vector<int> Rcoord;
    int i, j, k, l, m;
    int noThisMolecule;
    int currentMolecule;

    // which coordinates should be moved (only R(0))... we are before integrator initialization
    Rcoord.push_back(0);

    // cykle through each molecular type
    for (k = 1; k < (int)molecularTypesBound.size(); k++)
    {
        // original number of the molecular type in config
        noThisMolecule = (int)round((molecularTypesBound.at(k) - molecularTypesBound.at(k - 1)) /
                                    (nfoldx * nfoldy * nfoldz));
        // move original box to the leftmost corner
        for (i = molecularTypesBound.at(k - 1); i < noThisMolecule; i++)
        {
            molecules[i].Move(-0.5 * (nfoldx - 1) * boxSize[0], 0, Rcoord);
            molecules[i].Move(-0.5 * (nfoldy - 1) * boxSize[1], 1, Rcoord);
            molecules[i].Move(-0.5 * (nfoldz - 1) * boxSize[2], 2, Rcoord);
        }
        // assign (leftmost) coordinates to images of the original box
        for (j = molecularTypesBound.at(k - 1) + noThisMolecule;
             j < molecularTypesBound.at(k); j++)
        {
            for (l = 0; l < molecules[j].noAtoms; l++)
            {
                molecules[j].atoms[l].R->operator()(0, 0) =
                    molecules[j - noThisMolecule].atoms[l].R->operator()(0, 0);
                molecules[j].atoms[l].R->operator()(0, 1) =
                    molecules[j - noThisMolecule].atoms[l].R->operator()(0, 1);
                molecules[j].atoms[l].R->operator()(0, 2) =
                    molecules[j - noThisMolecule].atoms[l].R->operator()(0, 2);
            }
        }
        // move images to final positions
        currentMolecule = molecularTypesBound.at(k - 1);
        for (i = 0; i < nfoldx; i++)
        {
            for (j = 0; j < nfoldy; j++)
            {
                for (l = 0; l < nfoldz; l++)
                {
                    for (m = 0; m < noThisMolecule; m++)
                    {
                        molecules[currentMolecule].Move(i * boxSize[0], 0, Rcoord);
                        molecules[currentMolecule].Move(j * boxSize[1], 1, Rcoord);
                        molecules[currentMolecule].Move(l * boxSize[2], 2, Rcoord);
                        currentMolecule++;
                    }
                }
            }
        }
    }

    // new box size
    boxSize[0] *= nfoldx;
    boxSize[1] *= nfoldy;
    boxSize[2] *= nfoldz;

    return 0;
}

// assign random velocities (-0.5, 0.5) to atoms
int SimulatedSystem::RandomVelocities()
{
    int i, j, k, m;
    Vector sum_mv(3); // total momentum is cumulated here
    Matrix *atomR;

    /* Initializing atom coordinates and velocities and writing coordinates to .plb file */
    srand((unsigned int)time(NULL));

    for (i = 0; i < noMolecules; i++) // for all molecules
    {
        for (j = 0; j < molecules[i].noAtoms; j++) // for each atom
        {
            atomR = molecules[i].atoms[j].R;
            for (k = 0; k < 3; k++)
            {
                ATOMR(1, k) = REAL_RAND() - 0.5;                       // assign random velocity between -0.5 and 0.5
                sum_mv(k) += molecules[i].atoms[j].mass * ATOMR(1, k); // sum momenta to subtract box COM movement
            }
            // sum_mv2 += molecules[i].atoms[j].mass * atomR->RowRowDotProduct(1, 1); // sum squares of velocities to get correct rescaling
        }
    }

    /* Total momentum to zero and higher derivatives in R set to 0.0 */
    sum_mv = sum_mv * (1 / totalMass); // velocity of COM to be subtracted from all velocities
    for (i = 0; i < noMolecules; i++)  // for all molecules
    {
        for (j = 0; j < molecules[i].noAtoms; j++) // for each atom
        {
            atomR = molecules[i].atoms[j].R;
            for (k = 0; k < 3; k++)
            {
                ATOMR(1, k) -= sum_mv(k);
                for (m = 2; m < Rsize; m++)
                {
                    ATOMR(m, k) = 0.0;
                }
            }
        }
    }

    return 0;
}

// rescale velocities to desired temperature (work with R(1,:) = v, not hv)
double SimulatedSystem::SetTemperature(double temp)
{
    int i, j;
    double sum_mv2 = 0; // 2Ekin is cumulated here
    double v_scf = 1;   // velocity scaling factor

    // calculate Ekin (2Ekin stored in sum_mv2)
    for (i = 0; i < noMolecules; i++)
    {
        for (j = 0; j < molecules[i].noAtoms; j++)
        {
            sum_mv2 += molecules[i].atoms[j].mass * molecules[i].atoms[j].R->RowRowDotProduct(1, 1);
        }
    }

    /* Set desired temperature */
    v_scf = sqrt(noDegrOfFred * temp / (sum_mv2)); // velocity scaling factor
    RescaleVelocities(v_scf);

    return v_scf;
}

// simply rescale velocities by a constant factor 'scaling'
int SimulatedSystem::RescaleVelocities(double scaling)
{
    int i, j, k;

    for (i = 0; i < noMolecules; i++)
    {
        for (j = 0; j < molecules[i].noAtoms; j++)
        {
            for (k = 0; k < 3; k++)
            {
                molecules[i].atoms[j].R->operator()(1, k) *= scaling;
            }
        }
    }
    return 0;
}

// rescale boxsize and positions (based on center of mass by default)
int SimulatedSystem::RescaleBox(double boxscaling, std::vector<int> positionsAt, bool molecularbased)
{
    int i, j, k;
    std::vector<int>::iterator it;

    if (molecularbased) // molecular based rescaling (protects bond lengths), complicated for large molecules
    {
        double center_of_mass[3];

        for (i = 0; i < noMolecules; i++)
        {
            molecules[i].CalculateCOM(center_of_mass);
            for (j = 0; j < 3; j++)
            {
                center_of_mass[j] = (-1 + boxscaling) * center_of_mass[j];
                molecules[i].Move(center_of_mass[j], j, positionsAt);
                // pairlist->MoveMolecule(i, j, center_of_mass[j]); // should not be done (??)
            }
        }
    }
    else // atom-based rescaling needed in different barostats ...
    {
        for (i = 0; i < noMolecules; i++)
        {
            for (j = 0; j < molecules[i].noAtoms; j++)
            {
                for (k = 0; k < 3; k++)
                {
                    for (it = positionsAt.begin(); it != positionsAt.end(); it++)
                    {
                        molecules[i].atoms[j].R->operator()(*it, k) *= boxscaling;
                    }
                }
            }
        }
        // pairlist->RescaleOrigPos(boxscaling); // should not be done (??)
    }

    for (j = 0; j < 3; j++)
    {
        boxSize[j] *= boxscaling;
    }
    return 0;
}

// rescale boxsize and positions (based on center of mass)
int SimulatedSystem::RescaleBox3D(Vector boxscaling, std::vector<int> positionsAt)
{
    int i, j;
    double center_of_mass[3];

    for (i = 0; i < noMolecules; i++)
    {
        molecules[i].CalculateCOM(center_of_mass);
        for (j = 0; j < 3; j++)
        {
            center_of_mass[j] = (-1 + boxscaling[j]) * center_of_mass[j];
            molecules[i].Move(center_of_mass[j], j, positionsAt);
        }
    }
    for (j = 0; j < 3; j++)
    {
        boxSize[j] *= boxscaling[j];
    }
    return 0;
}

// get Rsize
int SimulatedSystem::GetRsize() const
{
    return Rsize;
}

// calculate external forces (elliptic potential, gravity,...)
double SimulatedSystem::CalculateExternalForces()
{
    Eext = 0.0;
    int i, j, k;

    // external elliptic potential
    if (isExternElliptic)
    {
        for (i = 0; i < noMolecules; i++)
        {
            for (j = 0; j < molecules[i].noAtoms; j++)
            {
                for (k = 0; k < 3; k++)
                {
                    molecules[i].atoms[j].force[k] += -molecules[i].atoms[j].R->operator()(0, k) * ellipticK[k];
                    Eext += 0.5 * ellipticK[k] * (molecules[i].atoms[j].R->operator()(0, k) * molecules[i].atoms[j].R->operator()(0, k));
                }
            }
        }
    }

    return Eext;
}

// calculate potential energy caused by external forces (elliptic potential, gravity,...)
double SimulatedSystem::CalculateExternalEpot() const
{
    double Eextern = 0.0;
    int i, j, k;

    // external elliptic potential
    if (isExternElliptic)
    {
        for (i = 0; i < noMolecules; i++)
        {
            for (j = 0; j < molecules[i].noAtoms; j++)
            {
                for (k = 0; k < 3; k++)
                {
                    Eextern += 0.5 * ellipticK[k] * (molecules[i].atoms[j].R->operator()(0, k) * molecules[i].atoms[j].R->operator()(0, k));
                }
            }
        }
    }

    return Eextern;
}

// invert whole configuration (inversion center [0,0,0])
int SimulatedSystem::InvertBox(std::vector<int> positionsAt)
{
    int i, j, k;
    std::vector<int>::iterator l;
    for (i = 0; i < noMolecules; i++)
    {
        for (j = 0; j < molecules[i].noAtoms; j++)
        {
            for (k = 0; k < 3; k++)
            {
                for (l = positionsAt.begin(); l != positionsAt.end(); l++)
                {
                    molecules[i].atoms[j].R->operator()(*l, k) *= -1.0;
                }
            }
        }
    }

    return 0;
}

// reflect simulation along mirror in the middle of the "direction" axis
int SimulatedSystem::ReflectBox(std::vector<int> positionsAt, int direction)
{
    int i, j;
    std::vector<int>::iterator l;
    for (i = 0; i < noMolecules; i++)
    {
        for (j = 0; j < molecules[i].noAtoms; j++)
        {
            for (l = positionsAt.begin(); l != positionsAt.end(); l++)
            {
                molecules[i].atoms[j].R->operator()(*l, direction) *= -1.0;
            }
        }
    }

    return 0;
}

// add velocity to the whole configuration
int SimulatedSystem::AddVelocity(double velocity, int direction)
{
    int i, j;
    for (i = 0; i < noMolecules; i++)
    {
        for (j = 0; j < molecules[i].noAtoms; j++)
        {
            molecules[i].atoms[j].R->operator()(1, direction) += velocity;
        }
    }

    return 0;
}

// move whole configuration
int SimulatedSystem::MoveSystem(double displacement, int direction, std::vector<int> positionsAt)
{
    int i;

    for (i = 0; i < noMolecules; i++)
    {
        molecules[i].Move(displacement, direction, positionsAt);
    }

    return 0;
}

// print MACSIMUS style .cfg file (final configuration...)
int SimulatedSystem::PrintCfg(int configLevel, char *simname, double h) const
{
    // inspired by KlimaM output.cu (which was taken from MACSIMUS utilities)
    // writes binary file with system configuration and simulation properties...
    // not all properties are relevant and some may be used in other way then in the original...

    FILE *fcfg;
    char cfg_name[MAX_COMMENT];
    int i, j, k;
    char write[] = "wb";
    REC rec;
    A_S *a[3];
    double Ekin = 0.0;

    // .config file name
    strcpy(cfg_name, simname);
    strcat(cfg_name, ".cfg");

    // open .config file
    if (my_fopen_w(&fcfg, cfg_name, write) != 0)
    {
        print_warning(1, "Cannot open " + std::string(cfg_name) + " file to write (binary)\n", "    .cfg file will not be made!!!\n");
        // this is not a fatal error, thus continue simulation (or its closing) without it
        return -1;
    }

    //  atm_position - data_p pointer to atom positions [Angstrom]
    //  atm_velocity - data_p pointer to atom velocities [macsimus velocity units]
    //  Etot - double total energy [K]
    //  h - double time step [ps]
    //  sim_time - double running simulation time [ps]
    //  box_edge_length - double length of simulation box edge [Angstrom] (cubic box expected)
    //  h= time step [ps]
    //  sim_time= running simulation time [ps]
    //  box_edge_length= box sizes [Angstrom]
    //  Etot= total energy [K]

    int size = 0;

    // parse configLevel
    if ((int)floor(configLevel / 10.0) > 0)
    {
        print_warning(1, "writing .cfg with trvplevel>0 requested\n", "    This should never happen...\n");
        return -2;
    }

    ////First part of header
    int cfgkey = 1; // 1 for nonpolar version, 2 for polar version
    // MACSIMUS options -@, -a,...
    int options[32] = {-9, 100, 1, 9,                // irrelevant, polarizabilities scaling (%), check for .stp file each cycle, measuring of constraints errors
                       0, -1, 0, -32333,             // check site-site dist., items in .cp file (!?!), playback mode force calculation, -g (REMOVED)
                       0, 2, 100, -9,                // center atoms of given valence, interrupt option (2 = finish step, save and exit), LJ scaling (%), force constant for keeping sites in place (-k)
                       0, configLevel, 0, 1,         // dump of .asc, .vel, .for, etc.; integrator type (in MACSIMUS), write .vlb, ask for .loc file
                       229, 100, 1, 0,               // polar version option/trvp order (check again for NH?!?), charges rescaling, config level for dt.cfg, store fully restartable configuration
                       0, 0, 99999, 1,               // interactive mode, runtime measurements, fix bonds (-u option), verbosity level
                       6, 0, 2147483647, 0,          // write configuration check, water models options, write playback of first no molecules, seed for generator
                       -32333, 0, 0, -32333};        // POLAR only, NSLOTS (removed), water-protein, debugging option
    int nspec = (int)molecularTypesBound.size() - 1; // number of molecule types

    int *species_specification; // for each molecule type
    if ((species_specification = (int *)malloc(2 * nspec * sizeof(int))) == NULL)
    {
        print_warning(1, "Error during memory allocation in .cfg file print\n", "    .cfg file corrupted\n");
        fclose(fcfg);
        return -3;
    }

    j = 0;
    for (i = 0; i < nspec; i++)
    {
        species_specification[j++] = molecularTypesBound[i + 1] - molecularTypesBound[i]; // number of molecules
        species_specification[j++] = molecules[molecularTypesBound[i]].noAtoms;           // number of sites in molecule
    }

    VarPut(fcfg, &cfgkey, sizeof(cfgkey));
    VarPut(fcfg, options, sizeof(options));
    VarPut(fcfg, &nspec, sizeof(nspec));
    VarPut(fcfg, species_specification, nspec * 2 * sizeof(int));

    free(species_specification);
    species_specification = nullptr;

    ////First type record specifications
    // energy calculation
    std::vector<Molecule>::const_iterator itmol;
    for (itmol = molecules.begin(); itmol != molecules.end(); itmol++)
    {
        Ekin += itmol->CalculateEkin(h, 1);
    }
    // first record
    rec.key = 1;                              // key for record type specification
    rec.intval = noMolecules;                 // Number of molecules
    rec.vecval[0] = 0.0;                      // running simulation time [ps] (here zero, because the value is stored elsewhere)
    rec.vecval[1] = h;                        // time step [ps]
    rec.vecval[2] = instEpot + Eextra + Ekin; // total energy [K]

    VarPut(fcfg, &rec, sizeof(rec));

    ////Second type record specifications
    rec.key = 2;               // key for record type specification
    rec.intval = noAtomsTotal; // Number of atoms
    rec.vecval[0] = GetBox(0); // box sizes [Angstrom]
    rec.vecval[1] = GetBox(1); // box sizes [Angstrom]
    rec.vecval[2] = GetBox(2); // box sizes [Angstrom]
    VarPut(fcfg, &rec, sizeof(rec));

    ////Thermostat specifications (stored elsewhere, add later...)
    rec.key = 4;       // key for record type specification
    rec.intval = 0;    // Thermostat type
    rec.vecval[0] = 0; // RvdW (for van der Waals radius setup)
    rec.vecval[1] = 0; // tau.T
    rec.vecval[2] = 0; // tau.P

    VarPut(fcfg, &rec, sizeof(rec));

    ////alocation of structures a[i]
    size = sizeof(*a[0]) + sizeof(vector3) * (noAtomsTotal - 1);
    for (i = 0; i < 3; i++)
    {
        if (((a[i]) = (A_S *)malloc(size)) == nullptr)
        {
            print_warning(1, "Error during memory allocation in .cfg file print\n", "    .cfg file corrupted\n");
            fclose(fcfg);
            return -3;
        }
        memset(a[i], 0, size);
        a[i]->size = size;
    }

    ////logs // Nosé–Hoover Xi (later)
    a[0]->logs = 0;
    a[1]->logs = 0;
    a[2]->logs = 0;
    ////ln(L)=lambda
    a[0]->lambda[0] = log(boxSize[0]);
    a[0]->lambda[1] = log(boxSize[1]);
    a[0]->lambda[2] = log(boxSize[2]);
    // lambda derivatives * h
    a[1]->lambda[0] = 0;
    a[1]->lambda[1] = 0;
    a[1]->lambda[2] = 0;
    // lambda second derivatives * h^2/2
    a[2]->lambda[0] = 0;
    a[2]->lambda[1] = 0;
    a[2]->lambda[2] = 0;

    ////setting of positions, velocities and accelerations
    for (j = 0, i = 0; j < noMolecules; j++)
    {
        for (k = 0; k < molecules[j].noAtoms; k++)
        {
            //// positions x,y,z
            a[0]->rp[i][0] = molecules[j].atoms[k].GetPosition(0) + 0.5 * boxSize[0];
            a[0]->rp[i][1] = molecules[j].atoms[k].GetPosition(1) + 0.5 * boxSize[1];
            a[0]->rp[i][2] = molecules[j].atoms[k].GetPosition(2) + 0.5 * boxSize[2];
            //// h*velocities x,y,z
            a[1]->rp[i][0] = molecules[j].atoms[k].R->operator()(1, 0);
            a[1]->rp[i][1] = molecules[j].atoms[k].R->operator()(1, 1);
            a[1]->rp[i][2] = molecules[j].atoms[k].R->operator()(1, 2);
            //// h*h*accelerations x,y,z
            a[2]->rp[i][0] = molecules[j].atoms[k].R->operator()(2, 0);
            a[2]->rp[i][1] = molecules[j].atoms[k].R->operator()(2, 1);
            a[2]->rp[i][2] = molecules[j].atoms[k].R->operator()(2, 2);
            i++;
        }
    }

    //// write into cfg file
    VarPut(fcfg, a[0], size);
    VarPut(fcfg, a[1], size);
    VarPut(fcfg, a[2], size);
    free(a[2]);
    free(a[1]);
    free(a[0]);

    // EOF marked by zero and time in ms (current time)
    i = 0;
    if (fwrite(&i, sizeof(i), 1, fcfg) != 1)
    {
        print_warning(1, "Error during .cfg printing...\n");
        return -3;
    }

    // current time

    // close .config file
    if (my_fclose(&fcfg, cfg_name) != 0)
    {
        print_warning(1, "Error during .cfg file closing\n", "    file can be corrupted!!!\n");
        // not a fatal error
        return -2;
    }

    return 0;
}

// auxiliary function for PrintCfg, inspired by the MACSIMUS function of the same name
int VarPut(FILE *file, void *v, int size)
{
    // file is opened in wb (or ab) mode
    // void *v is anything you want to write byte by byte
    // size is int size of v
    int i, s;
    unsigned char a;
    int written = 0;

    if (!size)
        return -1; // nothing to write

    written += (int)fwrite(&size, sizeof(size), 1, file);
    s = 0x4A4B; // magic constant (no idea why 0x4a4b)

    for (i = 0; i < size; i++) // originally loop(I,FROM,TO) for ((I)=(FROM); (I)<(TO); (I)++)
    {
        a = *((unsigned char *)v + i); // by explicit type conversion (casting) of pointer v (which can point to any date type) to *unsigned char we are reading v byte by byte
        s += a;
        if (putc(a, file) != EOF)
            written++; // putc (or fputc) casts automatically the argument to unsigned char (formally, its 1st argument is 'int')
    }
    written += (int)fwrite(&s, sizeof(s), 1, file);

    return written;
}

// set corrections that should be done (see above)
int SimulatedSystem::CalculateMixTerms(int mixrule)
{
    std::vector<AbstractInterMolField *>::iterator it;

    for (it = interMolFields.begin(); it != interMolFields.end(); it++)
    {
        if ((*it)->GetType() == 0)
        {
            continue;
        }
        (*it)->CalculateMixTerms(mixrule);
    }

    return 0;
}

// set corrections that should be done (see above)
int SimulatedSystem::SetCorrections(bool cutoff, bool pressure)
{
    std::vector<AbstractInterMolField *>::iterator it;

    pressureNfCorrection = pressure;
    cutoffCorrections = cutoff;

    double aux;

    // if periodic boundary conditions then calculate cutoff related issues
    // C1 and A - parameters of cutoff correction
    if ((boundaryCond > 0) && ((vdWcutoff > 0.5 * boxSize[0]) || (vdWcutoff > 0.5 * boxSize[1]) || (vdWcutoff > 0.5 * boxSize[2])))
    {
        aux = vdWcutoff;
        vdWcutoff = 0.5 * fmin(fmin(boxSize[0], boxSize[1]), boxSize[2]);
        print_warning(1, "Requested vdW (LJ) cutoff bigger than half the smallest box size\n",
                      "    Requested: " + std::to_string(aux) + ", box size in .config file: [",
                      std::to_string(boxSize[0]) + ", " + std::to_string(boxSize[1]) + ", " + std::to_string(boxSize[2]) + "]\n",
                      "    Setting vdW cutoff to: " + std::to_string(vdWcutoff) + " [AA]\n");
    }

    if ((boundaryCond > 0) && (elstcutoff != -999.0) && ((elstcutoff > 0.5 * boxSize[0]) || (elstcutoff > 0.5 * boxSize[1]) || (elstcutoff > 0.5 * boxSize[2])))
    {
        aux = elstcutoff;
        elstcutoff = 0.5 * fmin(fmin(boxSize[0], boxSize[1]), boxSize[2]);
        print_warning(1, "Requested electrostatics cutoff bigger than half the smallest box size\n",
                      "    Requested: " + std::to_string(aux) + ", box size in .config file: [",
                      std::to_string(boxSize[0]) + ", " + std::to_string(boxSize[1]) + ", " + std::to_string(boxSize[2]) + "]\n",
                      "    Setting elstat cutoff to: " + std::to_string(elstcutoff) + " [AA]\n");
    }

    if (vdWcutoff < 1.0)
    {
        print_warning(1, "Requested vdW (LJ) cutoff smaller than 1.0 [AA] (probably a typo?).\n",
                      "    Requested: " + std::to_string(vdWcutoff) + ", box size in .config file: [",
                      std::to_string(boxSize[0]) + ", " + std::to_string(boxSize[1]) + ", " + std::to_string(boxSize[2]) + "]\n");
    }

    if ((elstcutoff != -999.0) && (elstcutoff < 1.0))
    {
        print_warning(1, "Requested electrostatics cutoff smaller than 1.0 [AA] (probably a typo?).\n",
                      "    Requested: " + std::to_string(elstcutoff) + ", box size in .config file: [",
                      std::to_string(boxSize[0]) + ", " + std::to_string(boxSize[1]) + ", " + std::to_string(boxSize[2]) + "]\n");
    }

    for (it = interMolFields.begin(); it != interMolFields.end(); it++)
    {
        if (error = (*it)->SetCutoff(vdWcutoff, elstcutoff, boundaryCond), error != 0)
        {
            print_warning(1, "Problems when assigning cutoff (in SimulatedSystem.cpp::SetCorrections()).\n",
                          "Maybe too short cutoff for some pairs. Error number: " + std::to_string(error) + ".\n");
        }
    }

    if (!pressureNfCorrection)
    {
        noDOFforPressure = noDegrOfFred;
    }

    if (!cutoffCorrections)
    {
        Ecorr = 0.0;
    }
    else if (boundaryCond > 0)
    {
        Ecorr = 0.0;
        for (it = interMolFields.begin(); it != interMolFields.end(); it++)
        {
            Ecorr += (*it)->CalculateCutoffCorrection();
        }
    }

    return 0;
}

// set energy unit u_eng (which is private)
int SimulatedSystem::SetEnergyUnit(const double ueng)
{
    u_eng = ueng;
    return 0;
}

// add boundary to molboundTypes
int SimulatedSystem::AddMolTypeBound()
{
    molecularTypesBound.push_back(molecules.size());
    return 0;
}

// add fake bonds to handle disperse and electrostatic interactions inside Molecule
// set corrections that should be done (see above)
int SimulatedSystem::AddIntramolInteractions(double lj14, double elstat14)
{
    std::vector<AbstractInterMolField *>::iterator it;

    for (it = interMolFields.begin(); it != interMolFields.end(); it++)
    {
        (*it)->InitializeIntramolBonds(lj14, elstat14);
    }
    return 0;
}

int SimulatedSystem::SetPressForH(double constP)
{
    Pext = constP;
    return 0;
}

int SimulatedSystem::PrintDump(char *cpaname, char *plbname, char *histname)
{
    char auxname[100];

    if (fcpa != nullptr)
    {
        if (fclose(fcpa) != 0)
        {
            return 1;
        }
        strcpy(auxname, cpaname);
        strcat(auxname, ".dump");
        my_copy(cpaname, auxname);
        if ((fcpa = fopen(cpaname, "a")) == nullptr)
        {
            return 2;
        }
    }

    if (fplb != nullptr)
    {
        if (fclose(fplb) != 0)
        {
            return 1;
        }
        strcpy(auxname, plbname);
        strcat(auxname, ".dump");
        my_copy(plbname, auxname);
        if ((fplb = fopen(plbname, "ab")) == nullptr)
        {
            return 2;
        }
    }

    if (fhist != nullptr)
    {
        if (fclose(fhist) != 0)
        {
            return 1;
        }
        strcpy(auxname, histname);
        strcat(auxname, ".dump");
        my_copy(histname, auxname);
        if ((fhist = fopen(histname, "a")) == nullptr)
        {
            return 2;
        }
    }

    return 0;
}

// print system info to .prt file
void SimulatedSystem::PrintInfo(std::ofstream &fprt, std::string &engunit) const
{
    int i;
    int no_inter[INTER_MOL_TYPES] = {0};
    // energy unit
    if (fabs(u_eng - U_ENERGY_KCALMOL) < 1e-12)
    {
        engunit = "[kcal]";
    }
    else if (fabs(u_eng - U_ENERGY_KJMOL) < 1e-12)
    {
        engunit = "[kJ]";
    }
    else
    {
        engunit = "[K]";
    }
    fprt << "----------------------------------------------------------------------\n\n";
    fprt << "# System\n\n";
    fprt << "## Summary\n\n";

    MDTable<std::string, double> summary({"Variable", "Value"});
    summary.setColumnPrecision({0, 9});

    summary.addRow("Number of molecules", (double)noMolecules);
    summary.addRow("Types of molecules", (double)molecularTypesBound.size() - 1);
    summary.addRow("Total number of atoms", (double)noAtomsTotal);
    summary.addRow("Total mass [g/mol]", totalMass * U_MASS_GMOL);
    summary.addRow("Total number of constraints", (double)noConstrTotal);
    summary.addRow("Boundaries (0 = free, 2 = periodic)", boundaryCond);
    summary.addRow("Number of degrees of freedom (DOF)*", (double)noDegrOfFred);
    summary.addRow("Number of DOF used to pressure calc.", noDOFforPressure);
    summary.addRow("Initial system volume [AA^3]", CalculateVolume());
    summary.addRow("Initial density [kg/m3]", totalMass / CalculateVolume() * U_DENSITY_KGM3);

    summary.print(fprt);

    if (boundaryCond > 0)
    {
        fprt << "\\*) 3 degrees of freedom are subtracted for the total momentum conservation\n\n";
    }
    else
    {
        fprt << "* no degrees of freedom subtracted in free boundary conditions\n\n";
    }

    fprt << "To calculate density [kg/m^3] from volume [AA^3] use: $\\rho = " << std::fixed << std::setprecision(12) << totalMass * U_DENSITY_KGM3 << " / V$\n\n";

    fprt << "Initial positions of the very first atom (provided to check which `.config` file was used):\n\n";
    fprt << "```(config)\n";
    for (i = 0; i < 3; i++)
    {
        fprt << std::scientific << std::setprecision(12) << molecules[0].atoms[0].GetPosition(i) << "  ";
    }
    fprt << std::defaultfloat << "\n```\n\n";

    // Molecules info

    for (i = 0; i < (int)(molecularTypesBound.size() - 1); i++)
    {
        molecules[molecularTypesBound[i]].PrintInfo(fprt, i, molecularTypesBound[i + 1] - molecularTypesBound[i], u_eng, engunit);
    }

    // Intermol. fields info
    if (interMolFields.size() > 0)
    {
        fprt << "------------------------------------------------------------\n\n";
        fprt << "\n## Intermolecular interactions\n\n";
    }
    // Pair list info
    fprt << "### Pair list\n\n";
    pairlist->PrintInfo(fprt);
    fprt << "------------------------------------------------------------\n\n";
    // get number of intermol fields according to their type
    for (auto it = interMolFields.begin(); it != interMolFields.end(); it++)
    {
        no_inter[(*it)->GetType()]++;
    }
    // disperse first
    if (no_inter[1] > 0)
    {
        fprt << "### Disperse (van der Waals) interactions\n\n";
        for (auto it = interMolFields.begin(); it != interMolFields.end(); it++)
        {
            if ((*it)->GetType() == 1)
            {
                (*it)->PrintInfo(fprt, u_eng, engunit);
            }
        }
    }

    // electrostatics now
    if (no_inter[0] > 0)
    {
        if (no_inter[1] > 0)
        {
            fprt << "------------------------------------------------------------\n\n";
        }
        fprt << "### Electrostatic interactions\n\n";
        for (auto it = interMolFields.begin(); it != interMolFields.end(); it++)
        {
            if ((*it)->GetType() == 0)
            {
                (*it)->PrintInfo(fprt, u_eng, engunit);
            }
        }
    }
}

// close opened files
void SimulatedSystem::CloseFiles()
{
    char appendix[10];
    strcpy(appendix, ".plb");

    if ((fplb != NULL) && (my_fclose(&fplb, appendix) != 0))
    {
        print_warning(0, "Cannot close .plb file\n");
    }
    else
    {
        fplb = NULL;
    }
    strcpy(appendix, ".cpa");
    if ((fcpa != NULL) && (my_fclose(&fcpa, appendix) != 0))
    {
        print_warning(0, "Cannot close .cpa file\n");
    }
    else
    {
        fcpa = NULL;
    }
    strcpy(appendix, ".history");
    if ((fhist != NULL) && (my_fclose(&fhist, appendix) != 0))
    {
        print_warning(0, "Cannot close .history file\n");
    }
    else
    {
        fhist = NULL;
    }
    strcpy(appendix, ".acc");
    if ((facc != NULL) && (my_fclose(&facc, appendix) != 0))
    {
        print_warning(0, "Cannot close .acc file\n");
    }
    else
    {
        facc = NULL;
    }
}

// get number of measured quantities
int SimulatedSystem::GetNoQuant() const
{
    return measuredQuant.size();
}

// allocate Atom::R
int SimulatedSystem::AllocateR()
{
    double *p_memory;
    int i, j;

    p_memory = new double[Rsize * 3 * noAtomsTotal];
    memoryForR = p_memory;
    for (i = 0; i < noMolecules; i++)
    {
        for (j = 0; j < molecules[i].noAtoms; j++)
        {
            molecules[i].atoms[j].R = new Matrix(Rsize, 3, p_memory);
            p_memory += 3 * Rsize;
        }
    }

    return 0;
}

// deallocate Atom::R
int SimulatedSystem::DeallocateR()
{
    int i, j;

    if (memoryForR == nullptr)
    {
        return -1;
    }

    delete[] memoryForR;
    for (i = 0; i < noMolecules; i++)
    {
        for (j = 0; j < molecules[i].noAtoms; j++)
        {
            molecules[i].atoms[j].R->DeleteMData();
        }
    }
    memoryForR = nullptr;

    return 0;
}