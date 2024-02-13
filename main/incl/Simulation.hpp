#ifndef SIMULATIONHEADER
#define SIMULATIONHEADER

/*
#include "Vector.hpp"
#include "Matrix.hpp"
*/

#include "SplineType.hpp"
class AbstractIntegrator;

class Simulation
{
private:
    int Rsize;                // size of atom R needed
    double Tfinal;            // temperature set
    double tauT;              // thermostat constant
    double rhofinal;          // final density
    double taurho;            // time to reach final density
    double Pfinal;            // pressure set
    double tauP;              // barostat constant
    int ensembleLabel;        // ensemble label (nve (0), nvt berendsen (1), ntvinit berendsen (2), nvt nose (3), npt nose (4)...)
    int integratorLabel;      // integrator used label (verlet(0), gear(1), verlet+TRVP(2), gear+TRVP(3), verlet+iter(4), verlet+TRVP+iter(6), alejandre(-1, NPT nosé/mtk/hoover only) (,...?))
    int iterLabel;            // number of iterations
    char integratorName[30];  // integrator name ("\x1b[%dm"human-readable)
    int initLabel;            // initialization label (continue(0), restart(1 – default), scalevel(2), randomvel(3), random(4))
    int elstatlabel;          // type of elstat interactions (no(0), cutoff(1), ewald(2), ewald+splines(3))
    SplineType splineslabel1; // type of splines for ewald+splines see SplineType
    unsigned int gridsize1;   // gridsize for interpolators (energy)
    SplineType splineslabel2; // type of splines for ewald+splines forces
    unsigned int gridsize2;   // gridsize for interpolators (forces)
    int pairlistLabel;        // type of pair list (allpairs(0))
    double shakeTol;          // relative tolerance of SHAKE
    double shakeOmega;        // overrelaxation for SHAKE
    Vector *trvpCoeff;           // TRVP coefficients (for both Verlet and Gear formalism)
    char simName[MAX_COMMENT];   // simulation (and .control) name
    int nFold[3];                // proliferate system in x y z
    int historyLevel;            // information level in .history file
    std::set<int> measuredQuant; // set of measured quantities
    bool printcfg;               // print final configuration also in MACSIMUS .cfg format (beside the DL_POLY .config...)
    bool pressureNfCorrection;   // calculate pressure correction according to MACSIMUS manual eq. 15.15?
    bool cutoffCorrections;      // calculate cutoff corrections to energy and pressure?
    unsigned int elstatshift;    // electrostatics shift (&1 force shift, &2 energy shift)
    double elalpha;              // electrostatics parameter alpha
    double elkappa;              // electrostatics parameter kappa
    double lj14;                 // scaling factor for 1–4 disperse interactions
    double elstat14;             // scaling factor for 1–4 electrostatic interactions
    bool startWithShake;         // perform Shake during initialization (after first forces calculation), default true
    double vel_scaling;          // velocity scaling factor (in Initialize, reported in PrintInfo)
    double list_save;            // parameter for pair lists (save border)
    int mixingrule;              // mixing rule for vdw interactions

public:
    double h;    // integration step
    int noSteps; // number of integration steps
    double simulationTime;
    int noMeasurementsDone; // number of measurements already done
    double measurementFreq; // frequency of measurement (time between two consecutive measurements)/
    double plbFreq;         // frequency of .plb writings
    double histFreq;        // frequency of .history writings
    double dumpFreq;
    int freqMeasPlbHistDump[4];     // number of steps between two measurements, .plb write, .history write and .dump write
    int freqSteps;                  // greatest common divisor of the three frequencies above
    char fieldName[MAX_COMMENT];    // .field file name
    char configName[MAX_COMMENT];   // .config file name
    double LJcutoff;                // vdw interactions spherical cutoff
    double elcutoff;                // electrostatics interaction cutoff
    int error;                      // stores error information to handle
    SimulatedSystem *system;        // simulated system of molecules...
    AbstractIntegrator *integrator; // specific integrator used

    Simulation(char *simName, std::ofstream &fprt); // default constructor
    ~Simulation();                                  // default destructor
    int Execute();                                  // perform simulation and measurement
    void SetSystem(SimulatedSystem *simSystem);     // assign system
    int Initialize();                               // initialize integrator, ...
    int ReadControl(std::ofstream &fprt);           // read .control file
    int SelectIntegrator(std::string integrator);   // select integrator (parsing .control keyword "integrator")
    int SelectEnsemble(std::string ensemble);       // select ensemble (parsing .control keyword "ensemble")
    int SelectHistory(std::string historyline);
    int SelectInit(std::string initialization);                           // parse 'init' command for ReadControl()
    int AddMeasuredQuantities(std::string quant);                         // parse 'measure' command for ReadControl()
    int GetRsize() const;                                                 // get Rsize for atom initialization
    int GetTRVPsize() const;                                              // get TRVP size for system initialization
    int GetTotalProliferation() const;                                    // get total proliferation number (nFold[0]*nFold[1]*nFold[2])
    int SaveNfold(std::string nfold);                                     // save values of nfold from .control
    int PrintInfo(std::ofstream &fprt, std::string engunit);              // print simulation parameters to .prt file
    int SelectElectrostatics(std::string elstat);                         // electrostatics handling class selection
    int SelectPairList(std::string pairlist);                             // select pair list
    int SetScaling14(std::string scaleline);                              // set scaling for 1–4 interactions
    int SelectMixingRule(std::string mixline);                            // select mixing rule for vdw interactions
    int SelectSplines(std::string splinesline);                       // select type of interpolation for ewald+splines
    std::vector<std::string> GetMeasuredQuant(std::string engunit) const; // get measured quantities names
};

#endif