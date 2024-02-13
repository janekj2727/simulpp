#ifndef SIMULATEDSYSTEMHEADER
#define SIMULATEDSYSTEMHEADER

#define MAX_COMMENT 100
#define INTER_MOL_TYPES 2

#include <vector>
#include <set>

/*
#include "../mymath/Vector.hpp"
#include "../mymath/Matrix.hpp"
*/
#include "Atom.hpp"
#include "Molecule.hpp"
#include "MeanSquareDispl.hpp"
class AbstractInterMolField;
class AbstractPairList;

#ifdef PARALLEL
extern int thread_count;
#include <omp.h>
#endif

// functions related to measurement are now in Measurement.cpp (not in SimulatedSystem.cpp as one would expect)
class SimulatedSystem
{
private:
    SimulatedSystem() : Rsize(5), TRVPsize(0){}; // disable constructor without arguments
    FILE *fplb;
    FILE *fcpa;
    FILE *fhist;
    FILE *facc;
    int ReadConfig(char *config_name, int totalProlif);
    double CalculateMaxConstrErr() const;    // calculate maximum relative error of constrained bonds
    double CalculateMaxConstrVAngle() const; // calculate maximum deviation from pi/2 angle between velocity and constraint bond
    double ForceFromEpot(Atom *atom, int direction,
                         double displacement = 0.0001) const; // calculate force on atom in direction (x,y,z) by diferentiation of Epot
    int PrintAccelsFromEpot(double time) const;               // print acceleration of atoms to file ffor
    int AllocateR();                                          // allocate Atom::R
    int DeallocateR();                                        // deallocate Atom::R

    double boxSize[3];   // system box size
    double ellipticK[3]; // force constants of elliptic external potential
    bool isExternElliptic;
    const int Rsize;                      // size of Atom's R (position, velocity, ... matrix)
    const int TRVPsize;                   // size of TRVP (determines the content of Atom::R)
    double Ecorr;                         // long-range interaction cutoff correction
    std::vector<int> molecularTypesBound; // boundaries between different molecular types in molecules
    double Epotintramol[INTRA_MOL_TYPES]; // energy of intramolecular interactions bonds(0), angles(1), dihedrals(2), lj(3), elstat(4)
    double Eext;                          // potential energy of external forces
    std::vector<int> measuredQuant;       // measured quantities
    double u_eng;                         // energy unit conversion
    bool pressureNfCorrection;            // calculate pressure correction according to MACSIMUS manual eq. 15.15?
    bool cutoffCorrections;               // calculate cutoff corrections to energy and pressure?
    double Pext;                          // constant pressure −> if true, H = U + P_ext <V>, else H = U + <P> V; where U = Epot + Ekin + Ecorr
    double *memoryForR;                   // memory block for Atom::R

public:
    int noMolecules;                                     // number of molecules in the system
    std::vector<Molecule> molecules;                     // array of molecules in the system
    int noAtomsTotal;                                    // total number of atoms in the system
    int noConstrTotal;                                   // total number of CONSTRAINTS
    int noDegrOfFred;                                    // total number of degrees of freedom
    double totalMass;                                    // total mass of the simulated system
    double instEpot;                                     // instantenous potential energy
    int boundaryCond;                                    // boundary condition (imcon key), 0 (free), 1 (cubic), 2 (periodic)
    int error;                                           // stores error indication for program-flow control
    double epsShake;                                     // maximum relative error of constr. bonds
    double Eextra;                                       // energy of extra degrees of freedom (calculated by integrator)
    double virintermol[INTER_MOL_TYPES];                 // intermolecular virial: electrostatic(0), LJ forces(1)
    double Epotintermol[INTER_MOL_TYPES];                // intermolecular interactions energy (electrostatic(0), Lennard-Jones(1))
    double virintramol[INTRA_MOL_TYPES];                 // virial: bonds(0), angles(1),... (see Epotintramol)
    double virconstr;                                    // contrained bond virial
    double Xi;                                           // variable connected with thermostat (NVT Nosé & NPT Nosé)
    double Xivel;                                        // velocity of Xi (see above)
    double Lambda;                                       // variable connected with barostat (NPT Nosé)
    double Lambdavel;                                    // velocity of Lambda (see above)
    int noPairs;                                         // current number of considered atomic pairs (intermolecular)
    int noDOFforPressure;                                // number of degrees of freedom used to calculate pressure (see MACSIMUS equation 15.15)
    int maxShakeIter;                                    // maximum number of SHAKE iteration from the last printing (during the whole 'cycle' (MACSIMUS meaning))
    double vdWcutoff;                                    // disperse (vdw) interactions cutoff
    double elstcutoff;                                   // electrostatic interactions cutoff
    double omegaShake;                                   // overrelaxation for SHAKE
    Matrix *extendedDOFs;                                // matrix of extended degrees of freedom (for .config reading and printing)
    char extDOFnames[MAX_COMMENT];                       // names of extended degrees of freedom (to distinguish between xi and lambda)
    std::vector<AbstractInterMolField *> interMolFields; // intermolecular interactions
    AbstractPairList *pairlist;                          // list of intermolecular pairs
    MeanSquareDispl *msd;                                // mean square displacement calculation

    int InitPlbFile(char *plb_name, bool append);                                        // initialize .plb file
    int InitCpaFile(char *cpa_name, std::set<int> measuredQ, bool append, double &time); // initialize .cpa file
    int InitHistoryFile(char *hist_name, int no_frames, int hist_level, bool append);    // initialize .history file
    double CalculateInterMolForces();                                                    // calculate intermolecular forces (LJ, charges) - to be implemented later; returns Epotinter
    double CalculateInterMolEpot() const;                                                // calculate intermolecular potential energy
    double CalculateForces();                                                            // calculate all forces in the SimulatedSystem, returns total potential energy
    double CalculateEpot() const;                                                        // calculate potential energy (without forces calculation)
    double CalculateEkin(double h, int velocityIndex = 1) const;                         // calculate kinetic energy of the system
    double CalculateVolume() const;                                                      // calculate system volume
    double CalculateExternalForces();                                                    // calculate external forces (elliptic potential, gravity...)
    double CalculateExternalEpot() const;                                                // calculate potential energy caused by external field
    int AddMolecule(Molecule *mol);                                                      // initialization
    int PrintMeasurement(double time, double h);                                         // print measured quantities
    int PrintPlayback();                                                                 // print configuration to plb
    int PrintHistory(double h, int hist_level, double curr_time);                        // print configuration to .history
    double GetBox(int coord) const;                                                      // get box size element
    int SetBox(double size, int coord);                                                  // set box size directly
    int ApplyPeriodicBC(std::vector<int> &Rrows);                                        // move molecules to remain in box vicinity...
    int PrintConfig(int configLevel, char *simname, double h);                           // print final configuration to .config file
    int PrintCfg(int configLevel, char *simname, double h) const;                        // print MACSIMUS-style .cfg file (alternative to .config file)
    int Proliferate(int nfoldx, int nfoldy, int nfoldz);                                 // proliferate system
    int RandomVelocities();                                                              // assign random velocities (-0.5, 0.5)
    double SetTemperature(double temp);                                                  // rescale velocities according to given temperature (work with v not h*v, returns scaling factor)
    int RescaleVelocities(double scaling);                                               // multiply all velocities by factor 'scaling'
    int RescaleBox(double boxscaling, std::vector<int> positionsAt,                      //
                   bool molecularbased = true);                                          // rescale box size and positions (based on COM)
    int RescaleBox3D(Vector boxscaling, std::vector<int> positionsAt);                   // rescale box separately in each direction
    int InvertBox(std::vector<int> positionsAt);                                         // invert configuration (inversion center [0, 0, 0])
    int ReflectBox(std::vector<int> positionsAt, int direction);                         // reflect box along mirror in the middle of the "direction" axis
    int AddVelocity(double velocity, int direction);                                     // add velocity in direction given by 2nd argument
    int MoveSystem(double displacement, int direction, std::vector<int> positionsAt);    // move whole configuration
    double EnergyCutoffCorr() const;                                                     // get energy cutoff correction
    double GetTotalVirial() const;                                                       // get the sum of all virial terms
    double CalculatePconf() const;                                                       // calculate configurational (nonideal) part of pressure from virial
    double CalculatePconfVVC();                                                          // calculate configurational (nonideal) pressure using virtual volume change
    double CalculatePkin(double h, double Ekininst = 0.0) const;                         // calculate kinetic (ideal) part of pressure
    int SetCorrections(bool cutoff, bool pressure);                                      // set corrections that should be done (see above)
    int CalculateMixTerms(int mixrule);                                                  // calculate mixing terms for vdw interactions (between different types of atoms)
    int SetEnergyUnit(const double ueng);                                                // set unit of energy to handle input and output conversions
    double GetEnergyUnit() const;                                                        // get energy unit ...
    int GetRsize() const;                                                                // get size of Matrix R storing Atoms positions, velocities,...
    int AddMolTypeBound();                                                               // add index dividing different types of molecules to molecularTypesBound
    int AddIntramolInteractions(double lj14, double elstat14);                           // add fake bonds to handle disperse and elstat interaction inside Molecule
    int SetPressForH(double constP);                                                     // set private member Pext for enthalpy measurement
    double CalculateEnthalpy(double Ekin) const;                                         // calculate enthalpy of the system
    int PrintDump(char *cpaname, char *plbname, char *histname);                         // close .cpa, .plb and .history, copy them to .dump files and open them for append
    void PrintInfo(std::ofstream &fprt, std::string &engunit) const;                     // print system info to .prt file
    void CloseFiles();                                                                   // close opened files
    int GetNoQuant() const;                                                              // get number of measured quantities
    double MeanVelocity(double h, int dim = 3) const;                                    // get mean velocity (should be 0)
    double VarianceOfVelocity(double h, int dim = 3) const;                              // variance of velocity distribution
    double KurtosisOfVelocity(double h, int dim = 3) const;                              // (excess) kurtosis of velocity distribution (0 for canonical)
    double LeapFrogTtr(double h, double Ekintot) const;                                  // LF (VERLET3) translational temperature
    double LeapFrogTin(double h, double Ekintot) const;                                  // LF (VERLET3) internal temperature

    SimulatedSystem(char *simname, char *config_name,
                    char *field_name, double LJcutoff,
                    int totalProlif, int sizeOfR = 5,
                    int sizeTRVP = 0);                                                        // default constructor
    SimulatedSystem(SimulatedSystem &simSystem, bool withfields = true);                      // copy constructor
    SimulatedSystem(SimulatedSystem &simSystem1, SimulatedSystem &simSystem2, int direction); // merging constructor
    ~SimulatedSystem();                                                                       // default destructor
};

#endif