/*
 * A class to store simulated system for trial simulations of contrained dynamics
 * Author JJ, Date Jan 2021
 */

#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <string>
#include <set>
#include <utility>
// #include <cstdio>

#include "Vector.hpp"
#include "math_utils.hpp"
#include "general_utils.hpp"
#include "SimulatedSystem.hpp"
#include "AbstractIntegrator.hpp"
#include "VerletIntegrator.hpp"
#include "GearIntegrator.hpp"
#include "VerletNVTBerendsen.hpp"
#include "GearNVTBerendsen.hpp"
#include "VerletNTVinitBerendsen.hpp"
#include "GearNTVinitBerendsen.hpp"
#include "GearNVTNose.hpp"
#include "VerletNVTNose.hpp"
#include "VerletNVTNoseTRVP.hpp"
#include "VerletNPTBerendsen.hpp"
#include "GearNPTBerendsen.hpp"
#include "VerletNPTNose.hpp"
#include "VerletNPTNoseTRVP.hpp"
#include "GearNPTNose.hpp"
#include "VerletNVTNoseIter.hpp"
#include "VerletNPTNoseIter.hpp"
#include "VerletNPTAlejandre.hpp"
#include "units.hpp"
#include "Simulation.hpp"
#include "CutoffElstat.hpp"
#include "EwaldElstat.hpp"
#include "EwaldSplinesElstat.hpp"
#include "Quantities.hpp"
#include "Quantities.cpp"
#include "MDTable.hpp"
#include "AllPairList.hpp"
#include "VerletList.hpp"

// #define DEBUG

// default constructor
Simulation::Simulation(char *simname, std::ofstream &fprt)
{
    int i;

    strcpy(simName, simname);
    h = 0.0;
    noSteps = 0;
    simulationTime = 0.0;
    noMeasurementsDone = 0;
    measurementFreq = 0.0;
    plbFreq = 0.0;
    histFreq = 0.0;
    dumpFreq = 0.0;
    integratorLabel = 0; // 0 Verlet
    iterLabel = 0;
    ensembleLabel = 0; // 0 NVE
    pairlistLabel = 0; // all pairs
    historyLevel = 0;  // keytraj for history
    Rsize = 4;         // Verlet by default
    nFold[0] = 1;
    nFold[1] = 1;
    nFold[2] = 1;
    LJcutoff = 10.0;
    integrator = nullptr;
    initLabel = 1;     // restart (positions and velocities from CONFIG)
    shakeTol = 1.0e-6; // relative tolerance for SHAKE
    taurho = 0.0;
    Tfinal = 300.0;
    tauT = 0.0;
    Pfinal = 101325 / U_PRESSURE_PA;
    tauP = 0.0;
    rhofinal = 1000.0;
    trvpCoeff = nullptr;
    printcfg = false;            // do not print .cfg file by default
    cutoffCorrections = true;    // calculate cutoff corrections to energy and pressure by default
    pressureNfCorrection = true; // calculate correction to pressure (number of degrees of freedom)
    elstatshift = 3;             // default to shift both energy and forces to avoid discontinuities
    lj14 = 1.0;
    elstat14 = 1.0;
    elstatlabel = 0;       // no electrostatics by default
    startWithShake = true; // start with SHAKE by default
    mixingrule = 0;        // no implicit mixing unless specified
    shakeOmega = 1.0;      // default - no overrelaxation

    splineslabel1 = hyperbolic;
    splineslabel2 = hyperbolic;
    gridsize1 = 256;
    gridsize2 = 0;

    // measured quantities - add Etot, Epot, Ekin, Tkin, Ttr, Tin, P (default measured)
    for (i = 1; i < 7; i++)
    {
        measuredQuant.insert(i);
    }
    measuredQuant.insert(15); // pressure is measured by default (if not free b.c.)

    error = ReadControl(fprt);
    // in case no error occured
    if (error == 0)
    {
        print_warning(2, ".control file read successfully\n");
    }
}

// destructor
Simulation::~Simulation()
{
    if (integrator != nullptr)
    {
        delete integrator;
    }
    if (trvpCoeff != nullptr)
    {
        delete trvpCoeff;
    }
}

// perform simulation
int Simulation::Execute()
{
    double curTime = 0.0; // current time
    int curStep = 0;      // current step
    int configLevel;      // .config file level
    int this_error = 0;
    int no_frames = 0;
    int next_measurement = 0;
    int next_plb = 0;
    int next_history = 0;
    int next_dump = freqMeasPlbHistDump[3];
    bool append = false;
    SimulatedSystem *system2;
    std::ifstream stopfile;

    char plb_name[MAX_COMMENT];
    char cpa_name[MAX_COMMENT];
    char hist_name[MAX_COMMENT];
    char dump_name[MAX_COMMENT];
    char stop_name[MAX_COMMENT];

    strcpy(plb_name, simName);
    strcat(plb_name, ".plb");
    strcpy(cpa_name, simName);
    strcat(cpa_name, ".cpa");
    strcpy(hist_name, simName);
    strcat(hist_name, ".history");
    strcpy(dump_name, simName);
    strcat(dump_name, ".config.dump");
    strcpy(stop_name, simName);
    strcat(stop_name, ".stp");

    // if free b.c. it does not make sense to measure pressure or volume or density
    if (system->boundaryCond == 0)
    {
        measuredQuant.erase(15); // pressure
        measuredQuant.erase(14); // pressure from VVC
        measuredQuant.erase(7);  // volume
        measuredQuant.erase(8);  // density
    }

    // the initial configuration is always measured except for 'init continue'
    // if 'init' == 'continue' then append .cpa file and not print initial measurement
    if (initLabel == 0)
    {
        append = true;
    }
    if ((error = system->InitCpaFile(cpa_name, measuredQuant, append, curTime)) != 0)
    {
        return error;
    }
    if (!append)
    {
        system->PrintMeasurement(curTime, h);
    }

    next_measurement += freqMeasPlbHistDump[0];

    if ((freqMeasPlbHistDump[1] != 0))
    {
        if ((error = system->InitPlbFile(plb_name, append)) == 0)
        {
            if (!append)
            {
                system->PrintPlayback();
            }
            next_plb += freqMeasPlbHistDump[1];
        }
        else
        {
            print_warning(0, "ERROR: Cannot open .plb file: ", std::string(plb_name));
            return error;
        }
    }
    if ((freqMeasPlbHistDump[2] != 0))
    {
        no_frames = (int)(floor(noSteps / freqMeasPlbHistDump[2]) + 1);
        if ((error = system->InitHistoryFile(hist_name, no_frames, historyLevel, append)) == 0)
        {
            if (!append)
            {
                system->PrintHistory(h, historyLevel, curTime);
            }
            next_history += freqMeasPlbHistDump[2];
        }
        else
        {
            print_warning(1, "Cannot open .history file: ", std::string(hist_name));
            return error;
        }
    }
#ifndef NOPROGRESSBAR
    print_progress_bar(0);
#endif

    while (curStep < noSteps)
    {
        integrator->Integrate(freqSteps);
        curTime += freqSteps * h;
        curStep += freqSteps;
#ifdef PARALELL
#pragma omp parallel num_threads(thread_count) sections
        {
#pragma omp section
#endif
            if (curStep == next_measurement)
            {
                error = system->PrintMeasurement(curTime, h);
                next_measurement += freqMeasPlbHistDump[0];
            }
#ifdef PARALELL
#pragma omp section
#endif
            if (curStep == next_plb)
            {
                system->PrintPlayback();
                next_plb += freqMeasPlbHistDump[1];
            }
#ifdef PARALELL
#pragma omp section
#endif
            if (curStep == next_history)
            {
                system->PrintHistory(h, historyLevel, curTime);
                next_history += freqMeasPlbHistDump[2];
            }
#ifdef PARALELL
        }
#endif
        if ((curStep == next_dump) && (curStep < noSteps))
        {
            system->PrintDump(cpa_name, plb_name, hist_name);

            // copy system to enable preparation of config (new system is created withoud any force field)
            system2 = new SimulatedSystem(*system, false);
            integrator->SetSystem(system2);

            // prepare Atom::R to Print Config
            // returns configLevel of the .configfinal file
            configLevel = integrator->PrepareConfig();

            // final configuration to .config file (with config level 2)
            error = system2->PrintConfig(configLevel, dump_name, h); // error is 0 if success., negative if problems
            integrator->SetSystem(system);
            delete system2;
            system2 = nullptr;
            integrator->RestoreExtendedDOFs();
            next_dump += freqMeasPlbHistDump[3];
        }
#ifndef NOPROGRESSBAR
        print_progress_bar((int)floor(100 * curStep / noSteps));
#endif
        if (error != 0)
        {
            return error;
        }

        if (stopfile.open(stop_name, std::ifstream::in), stopfile.good())
        {
            print_warning(1, "Stop file " + std::string(stop_name) + " found.\n",
                          "Configuration will be printed to .config (and .cfg) file and simulation then exited!\n");
            error = -999;
            this_error = error;
            break;
        }
    }

    // prepare Atom::R to Print Config
    // returns configLevel of the .configfinal file
    configLevel = integrator->PrepareConfig();

    // final configuration to .config file (with config level 2)
    error = system->PrintConfig(configLevel, simName, h); // error is 0 if success., negative if problems
    // if MACSIMUS .cfg file is requested
    if (printcfg)
    {
        error = system->PrintCfg(configLevel, simName, h); // print final configuration in MACSIMUS .cfg format
    }
    return ((error == 0) ? this_error : error);
}

// assign system
void Simulation::SetSystem(SimulatedSystem *simSystem)
{
    system = simSystem;
}

// read .control file
int Simulation::ReadControl(std::ofstream &fprt)
{
    bool finishReached = false;
    int error = 0;
    std::pair<std::set<std::string>::const_iterator, bool> insertion;
    std::set<std::string> userDefinedOptions;

    std::fstream controlFile;
    std::string controlName, line, keyword;
    std::size_t position;

    // name of .control file
    controlName.append(simName);
    controlName.append(".control");

    // default names of .field and .config
    strcpy(fieldName, simName);
    strcat(fieldName, ".field");
    strcpy(configName, simName);
    strcat(configName, ".config");

    // open .control file
    controlFile.open(controlName, std::ios::in);
    if (!controlFile)
    {
        print_warning(0, "Cannot open .control file (" + controlName + "). No such file!!!\n");
        return 1;
    }

    // print header to .control file replication in .prt file
    fprt << "`" << controlName << "` file:\n```(control)\n";

    // read .control file
    while (!controlFile.eof())
    {
        std::getline(controlFile, line);
        fprt << line << "\n";
        if ((line.length() == 0) || (line.at(0) == '#') || (line.at(0) == '!'))
        {
            // lines beginning with # or ! are treated as comments
            // empty lines are also ignored...
            continue;
        }
        // first word is the main keyword which determines later handling
        position = line.find(" ");
        keyword = line.substr(0, position);
        line.erase(0, position + 1);
        // comments on line ends must be stripped
        if ((position = line.find("#")) != std::string::npos)
        {
            line.erase(position, std::string::npos);
        }
        if ((position = line.find("!")) != std::string::npos)
        {
            line.erase(position, std::string::npos);
        }
        // check if not already defined, if yes print warning and ignore
        insertion = userDefinedOptions.insert(keyword);
        if (!insertion.second && !(keyword.compare("no") == 0))
        {
            print_warning(1, "Trying to set previously set option '" + keyword + "'\n",
                          "    The new value: " + line + " will be ignored.\n");
            continue;
        }
        // parse keyword (one line in .config file)
        if (keyword.compare("temperature") == 0)
        {
            Tfinal = stod(line);
            if ((Tfinal <= 0) || (Tfinal > 10000))
            {
                print_warning(1, "Invalid value of 'temperature' in .control file.\n",
                              "    Given: '" + line + "', converted to: " + std::to_string(Tfinal) + "\n",
                              "    Default value 300 K used\n");
                Tfinal = 300.0;
            }
        }
        else if (keyword.compare("pressure") == 0)
        {
            Pfinal = stod(line);
            if ((Pfinal <= -1e9) || (Pfinal > 1e10))
            {
                print_warning(1, "Invalid value of 'pressure' in .control file.\n",
                              "    Given: '" + line + "', converted to: " + std::to_string(Pfinal) + "\n",
                              "    Default value 101 325 Pa used\n");
                Pfinal = 101325.0;
            }
            Pfinal /= U_PRESSURE_PA;
        }
        else if ((keyword.compare("density") == 0) || (keyword.compare("rho") == 0))
        {
            rhofinal = stod(line) / U_DENSITY_KGM3;
            if ((rhofinal <= 0) || (rhofinal > 20))
            {
                print_warning(1, "Invalid value of 'density' in .control file.\n",
                              "    Given: '" + line + "', converted to: " + std::to_string(rhofinal * U_DENSITY_KGM3) + "kg/m^3\n",
                              "    Default value 1000 kg/m3 used\n");
                rhofinal = 1000.0 / U_DENSITY_KGM3;
            }
        }
        else if (keyword.compare("timestep") == 0)
        {
            h = stod(line);
        }
        else if (keyword.compare("steps") == 0)
        {
            noSteps = stoi(line);
            if (noSteps < 0)
            {
                print_warning(0, "Cannot perform negative number of steps.\n",
                              "    Wrong value of steps directive in .control file\n",
                              "    Given: '" + line + "', converted to: " + std::to_string(noSteps) + "\n");
                return 4;
            } // calculate greatest common divisor of three integers in vector abc (zeros ignored)
        }
        else if (keyword.compare("stats") == 0)
        {
            measurementFreq = stod(line);
            if (measurementFreq <= 0.0)
            {
                print_warning(0, "Negative or zero measurement frequency ('stats' directive in .control)\n",
                              "    Given: '" + line + "', converted to: " + std::to_string(measurementFreq) + "\n");
                return 6;
            }
        }
        else if (keyword.compare("playback") == 0)
        {
            plbFreq = stod(line);
            if (plbFreq < 0.0)
            {
                print_warning(0, "Negative plb frequency ('playback' directive in .control)\n",
                              "    Given: '" + line + "', converted to: " + std::to_string(plbFreq) + "\n");
                return 6;
            }
        }
        else if (keyword.compare("history") == 0)
        {
            if (SelectHistory(line) != 0)
            {
                print_warning(1, "Unknown setting in history directive\n",
                              "    Given: '" + line + "', no history will be produced...\n");
                histFreq = 0.0;
            }
        }
        else if (keyword.compare("dump") == 0)
        {
            dumpFreq = stod(line);
            if (dumpFreq < 0.0)
            {
                print_warning(0, "Negative dump frequency ('dump' directive in .control)\n",
                              "    Given: '" + line + "', converted to: " + std::to_string(dumpFreq) + "\n");
                return 6;
            }
        }
        else if (keyword.compare("ljcutoff") == 0)
        {
            LJcutoff = stod(line);
            if ((LJcutoff < 0.0) || (LJcutoff > 1000000))
            {
                print_warning(0, "Invalid 'ljcutoff' directive in .control file\n",
                              "    Given: '" + line + "', converted to: " + std::to_string(LJcutoff) + "\n");
                return 5;
            }
        }
        else if (keyword.compare("shake") == 0)
        {
            if (sscanf(line.c_str(), "%lf %lf", &shakeTol, &shakeOmega) < 2)
            {
                shakeOmega = 1.0;
            }
            if ((shakeTol <= 0.0) || (shakeTol > 0.5))
            {
                print_warning(1, "Shake tolerance (in 'shake' directive) in .control file out of 'normal values'\n",
                              "    Given: '" + line + "', converted to: " + std::to_string(shakeTol) + "\n",
                              "    Defaulted back to 1e-6\n");
                shakeTol = 1.0e-6;
            }
            if ((shakeOmega < 1.0) || (shakeOmega > 2.5))
            {
                print_warning(1, "Shake overrelaxation parameter (in 'shake' directive) in .control file out of 'normal values'\n",
                              "    Given: '" + line + "', converted to: " + std::to_string(shakeOmega) + "\n",
                              "    Defaulted back to 1.0\n");
                shakeOmega = 1.0;
            }
        }
        else if (keyword.compare("config") == 0)
        {
            strcpy(configName, line.c_str());
        }
        else if (keyword.compare("field") == 0)
        {
            strcpy(fieldName, line.c_str());
        }
        else if (keyword.compare("integrator") == 0)
        {
            if ((error = SelectIntegrator(line)) != 0)
            {
                print_warning(1, "Unknown integrator in .control file\n",
                              "    Given: '" + line + "', defaulted to Verlet...\n");
            }
        }
        else if ((keyword.compare("vdwmixing") == 0) || (keyword.compare("ljmixing") == 0))
        {
            if ((error = SelectMixingRule(line)) != 0)
            {
                print_warning(1, "Unknown mixing rule for vdw interactions in .control file\n",
                              "    Given: '" + line + "'\n",
                              "No cross-terms will be calculated!\n");
                mixingrule = 0;
            }
        }
        else if (keyword.compare("ensemble") == 0)
        {
            if ((error = SelectEnsemble(line)) > 0)
            {
                return error;
            }
        }
        else if (keyword.compare("nfold") == 0)
        {
            if (SaveNfold(line) != 0)
            {
                return 3;
            }
        }
        else if (keyword.compare("init") == 0)
        {
            initLabel = SelectInit(line);
            if ((initLabel < 0) || (initLabel > 4))
            {
                print_warning(1, "Unknown 'init' selection\n",
                              "     Given: '" + line + "', defaulted to 'restart'\n");
                initLabel = 1;
            }
        }
        else if (keyword.compare("measure") == 0)
        {
            if (AddMeasuredQuantities(line) != 0)
            {
                print_warning(1, "Unknown quantity in 'measure' command\n",
                              "    Given: " + line + "\n");
            }
        }
        else if (keyword.compare("cfg") == 0)
        {
            printcfg = true; // print final configuration in MACSIMUS .cfg format
        }
        else if (keyword.compare("no") == 0)
        {
            if (line.compare("cutoffcorr") == 0)
            {
                cutoffCorrections = false;
            }
            else if (line.compare("pressurecorr") == 0)
            {
                pressureNfCorrection = false;
            }
            else if ((line.compare("elstat") == 0) || (line.compare("elec") == 0))
            {
                elstatlabel = 0; // already default, but directive enabled for explicit statement
            }
            else if ((line.compare("shakeinit") == 0) || (line.compare("initshake") == 0))
            {
                startWithShake = false; // do not perform SHAKE during initialzation (after first forces calculation)
            }
            else if ((line.compare("elstatenergyshift") == 0) || (line.compare("elstatshiftenergy") == 0) || (line.compare("ewaldenergyshift") == 0))
            {
                // do not shift elstat energy (lead to jump in r-space electrostatics)
                elstatshift &= 0x01;
            }
            else if ((line.compare("elstatforceshift") == 0) || (line.compare("elstatshiftforce") == 0) || (line.compare("ewaldforceshift") == 0))
            {
                // do not shift elstat force (lead to jump in r-space electrostatics)
                elstatshift &= 0x02;
            }
            else if ((line.compare("elstatshift") == 0) || (line.compare("elstatshift") == 0) || (line.compare("ewaldshift") == 0))
            {
                // do not shift elstat energy (lead to jump in r-space electrostatics)
                elstatshift = 0;
            }
            else
            {
                print_warning(1, "Unknown 'no' directive in .control file\n",
                              "    Given: '" + keyword + " " + line + "'\n",
                              "    This directive will be ignored...\n");
            }
        }
        else if (keyword.compare("elstat") == 0)
        {
            if ((error = SelectElectrostatics(line)) != 0)
            {
                print_warning(0, "Unknown electrostatics in .control file\n",
                              "    Given: '" + line + "', cannot continue...\n");
                return error;
            }
        }
        else if (keyword.compare("pairlist") == 0)
        {
            if ((error = SelectPairList(line)) != 0)
            {
                print_warning(1, "Unknown pair list in .control file\n",
                              "    Given: '" + line + "'\n",
                              "Defaulted to all pairs...\n");
                pairlistLabel = 0;
            }
        }
        else if (keyword.compare("scale14") == 0)
        {
            if ((error = SetScaling14(line)) != 0)
            {
                print_warning(0, "Corrupted line 'scale14' in .control file\n",
                              "    Given: '" + line + "', cannot continue...\n");
            }
        }
        else if (keyword.compare("splines") == 0)
        {
            SelectSplines(line);
        }
        else if (keyword.compare("finish") == 0) // if "finish" then last line reached
        {
            finishReached = true;
            break;
        }
        else
        {
            print_warning(1, "Unknown directive in .control file\n",
                          "    Given: '" + keyword + "'\n",
                          "    This option will be ignored!\n");
        }
    }

    if (!finishReached)
    {
        print_warning(1, "'finish' keyword not reached during reading .control file\n");
    }

    // print footer to .control file replication in .prt file
    fprt << "```\n\n";

    if (h > 0)
    {
        freqMeasPlbHistDump[0] = (int)floor(measurementFreq / h);
        freqMeasPlbHistDump[1] = (int)floor(plbFreq / h);
        freqMeasPlbHistDump[2] = (int)floor(histFreq / h);
        freqMeasPlbHistDump[3] = (int)floor(dumpFreq / h);
        freqSteps = gcd4(freqMeasPlbHistDump); // greatest common divisor
        if (dumpFreq == 0.0)
        {
            freqMeasPlbHistDump[3] = 0;
        }
        simulationTime = noSteps * h;
    }
    else
    {
        print_warning(0, "Negative value of 'timestep' in .control file or not given at all.\n",
                      "    Cannot perform simulation with negative or zero integration step.\n");
        return 2;
    }

    if ((noSteps > 0) && (freqSteps == 0))
    {
        freqSteps = 10;
        print_warning(1, "Zero frequency (cycle length) while non-zero number of integration steps.\n",
                      "    cycle length defaulted to 10.\n");
        print_warning(1, "No measurements nor configuration playback for the whole simulation...\n");
    }

    if ((histFreq > 0.0) && (historyLevel >= Rsize))
    {
        print_warning(1, "History information level (keytraj) greater than size of atoms coordinates (Rsize) (" + std::to_string(historyLevel),
                      " >= " + std::to_string(Rsize) + ")\n",
                      "    historyLevel changed to Rsize - 1 = " + std::to_string(Rsize - 1) + "\n");
        historyLevel = Rsize - 1;
    }

    if ((taurho > simulationTime) && (ensembleLabel == 2))
    {
        print_warning(1, "Density reaching time exceeds simulation time\n",
                      "    Final density will not be reached\n");
    }

    if (ensembleLabel == 5)
    {
        if (pressureNfCorrection == true)
        {
            print_warning(1, "Pressure correction (number of degrees of freedom) set to true (requested)\n",
                          "    but ensemble NPT Nosé (MTK) selected which does not work with it.\n",
                          "    Pressure correction will NOT be done!\n");
        }
        pressureNfCorrection = false;
    }

    return 0;
}

// proliferate system, assign random velocities (and/or positions) if needed, initialize integrator, ...
int Simulation::Initialize()
{
    std::vector<AbstractIntraMolField *>::iterator it;

    // proliferate system if needed
    if (GetTotalProliferation() > 1)
    {
        system->Proliferate(nFold[0], nFold[1], nFold[2]);
    }

    system->vdWcutoff = LJcutoff;
    system->elstcutoff = elcutoff;
    system->omegaShake = shakeOmega;

    // electrostatics initialization
    if (elstatlabel == 0)
    {
        // no electrostatic interactions, do nothing...
        // signal for SimulatedSystem::SetCorrection...
        system->elstcutoff = -999.0;
        elcutoff = -999.0;
    }
    else if (elstatlabel == 1)
    {
        // cutoff elstat
        system->interMolFields.push_back(new CutoffElstat(system, elcutoff, elalpha, elstatshift));
        if (system->error != 0)
        {
            return system->error;
        }
    }
    else if (elstatlabel == 2)
    {
        // Ewald sumation elstat
        system->interMolFields.push_back(new EwaldElstat(system, elcutoff, elalpha, elkappa, elstatshift));
        if (system->error != 0)
        {
            return system->error;
        }
    }
    else if (elstatlabel == 3)
    {
        // Ewald sumation elstat
        system->interMolFields.push_back(new EwaldSplinesElstat(system, elcutoff, elalpha, elkappa, elstatshift, splineslabel1, gridsize1, splineslabel2, gridsize2));
        if (system->error != 0)
        {
            return system->error;
        }
    }

    // calculate mixing vdw interactions
    system->CalculateMixTerms(mixingrule);

    // set corrections to System
    system->SetCorrections(cutoffCorrections, pressureNfCorrection);

    // add intramol elstat and disperse bonds
    system->AddIntramolInteractions(lj14, elstat14);

    // constraint maximum relative error for SHAKE
    system->epsShake = shakeTol;

    // assign random positions or velocities if needed (according to 'init' directive in .control)
    if ((initLabel == 0) || (initLabel == 1)) // continue and restart (noscale)
    {
        // do nothing...
    }
    else if (initLabel == 2)
    {
        // scale velocities
        vel_scaling = system->SetTemperature(Tfinal);
    }
    else if (initLabel == 3)
    {
        // initialize random velocities
        system->RandomVelocities();
        vel_scaling = system->SetTemperature(Tfinal);
    }
    else if (initLabel == 4)
    {
        // initialize random positions
        // to be implemented later ...

        // initialize random velocities
        system->RandomVelocities();
        vel_scaling = system->SetTemperature(Tfinal);
    }

    // initialize pairlist
    if (pairlistLabel == 1) // Verlet list
    {
        system->pairlist = new VerletList(system, std::max(LJcutoff, elcutoff), list_save);
    }
    else // all pairs list as default
    {
        system->pairlist = new AllPairList(system, std::max(LJcutoff, elcutoff));
    }

    // choose and initialize correct integrator
    if ((integratorLabel == -1) && (ensembleLabel != 5))
    {
        integratorLabel = 0;
        strcpy(integratorName, "Verlet");
        print_warning(1, "Alejandre integrator chosen, but not in NPT MTK ensemble, defaulted to Verlet\n");
    }
    if (ensembleLabel == 0) // NVE ensemble
    {
        if (integratorLabel > 3) // if iter demanded (no effect)
        {
            integratorLabel -= 4;
        }
        if (integratorLabel == 0) // Verlet NVE
        {
            integrator = new VerletIntegrator(h, system, startWithShake);
        }
        else if (integratorLabel == 1) // Gear NVE
        {
            integrator = new GearIntegrator(h, system, Rsize, integratorName, *trvpCoeff, startWithShake);
        }
        else if (integratorLabel == 2) // Verlet NVE + TRVP (unnecessary, but possible)
        {
            integrator = new VerletIntegrator(h, system, *trvpCoeff, startWithShake);
        }
        else // Gear + TRVP (dtto)
        {
            integrator = new GearIntegrator(h, system, Rsize, integratorName, *trvpCoeff, startWithShake);
        }
    }
    else if (ensembleLabel == 1) // NVT Berendsen
    {
        if (integratorLabel > 3) // if iter demanded (no effect)
        {
            integratorLabel -= 4;
        }
        if (integratorLabel == 0) // Verlet NVT Berendsen
        {
            integrator = new VerletNVTBerendsen(h, system, Tfinal, tauT, startWithShake);
        }
        else if (integratorLabel == 1) // Gear NVT Berendsen
        {
            integrator = new GearNVTBerendsen(h, system, Rsize, integratorName, Tfinal, tauT, *trvpCoeff, startWithShake);
        }
        else if (integratorLabel == 2) // Verlet NVT Berendsen + TRVP (unnecessary but possible)
        {
            integrator = new VerletNVTBerendsen(h, system, *trvpCoeff, Tfinal, tauT, startWithShake);
        }
        else // Gear NVT Berendsen + TRVP (unnecessary but possible)
        {
            integrator = new GearNVTBerendsen(h, system, Rsize, integratorName, Tfinal, tauT, *trvpCoeff, startWithShake);
        }
    }
    else if (ensembleLabel == 2) // NTVinit Berendsen (Berendsen thermostat + volume scaling to final density)
    {
        if (integratorLabel > 3) // if iter demanded (no effect)
        {
            integratorLabel -= 4;
        }
        if (integratorLabel == 0) // Verlet NTVinit Berendsen
        {
            integrator = new VerletNTVinitBerendsen(h, system, Tfinal, tauT, rhofinal, taurho, startWithShake);
        }
        else if (integratorLabel == 1) // Gear NTVinit Berendsen
        {
            integrator = new GearNTVinitBerendsen(h, system, Rsize, integratorName, Tfinal, tauT, rhofinal, taurho, *trvpCoeff, startWithShake);
        }
        else if (integratorLabel == 2) // Verlet NTVinit Berendsen + TRVP (unnecessary but possible)
        {
            integrator = new VerletNTVinitBerendsen(h, system, *trvpCoeff, Tfinal, tauT, rhofinal, taurho, startWithShake);
        }
        else // Gear NTVinit Berendsen + TRVP (unnecessary but possible)
        {
            integrator = new GearNTVinitBerendsen(h, system, Rsize, integratorName, Tfinal, tauT, rhofinal, taurho, *trvpCoeff, startWithShake);
        }
    }
    else if (ensembleLabel == 3) // NVT Nose–Hoover
    {
        if (integratorLabel == 0) // Verlet NVT Nose–Hoover (MTTK algorithm)
        {
            integrator = new VerletNVTNose(h, system, Tfinal, tauT, startWithShake);
        }
        else if (integratorLabel == 1) // Gear NVT Nose–Hoover
        {
            integrator = new GearNVTNose(h, system, Rsize, integratorName, Tfinal, tauT, *trvpCoeff, startWithShake);
        }
        else if (integratorLabel == 2) // Verlet NVT Nose–Hoover + TRVP
        {
            integrator = new VerletNVTNoseTRVP(h, system, *trvpCoeff, Tfinal, tauT, startWithShake); // totally different from VerletNVTNose (not MTTK)
        }
        else if (integratorLabel == 3) // Gear NVT Nose–Hoover + TRVP (neccessary – depends on the specific integrator)
        {
            integrator = new GearNVTNose(h, system, Rsize, integratorName, Tfinal, tauT, *trvpCoeff, startWithShake);
        }
        else if (integratorLabel == 4) // Verlet NVT Nose–Hoover (iterations)
        {
            integrator = new VerletNVTNoseIter(h, system, Tfinal, tauT, startWithShake, iterLabel);
        }
        else // Verlet NVT Nose–Hoover + TRVP (iterations)
        {
            integrator = new VerletNVTNoseIter(h, system, *trvpCoeff, Tfinal, tauT, startWithShake, iterLabel);
        }
    }
    else if (ensembleLabel == 4) // NPT Berendsen
    {
        system->SetPressForH(Pfinal);
        if (integratorLabel > 3) // if iter demanded (no effect)
        {
            integratorLabel -= 4;
        }
        if (integratorLabel == 0) // Verlet NPT Berendsen
        {
            integrator = new VerletNPTBerendsen(h, system, Tfinal, tauT, Pfinal, tauP, startWithShake);
        }
        else if (integratorLabel == 1) // Gear NPT Berendsen
        {
            integrator = new GearNPTBerendsen(h, system, Rsize, integratorName, Tfinal, tauT, Pfinal, tauP, *trvpCoeff, startWithShake);
        }
        else if (integratorLabel == 2) // Verlet NPT Berendsen + TRVP (unnecessary and currently not used)
        {
            integrator = new VerletNPTBerendsen(h, system, *trvpCoeff, Tfinal, tauT, Pfinal, tauP, startWithShake);
        }
        else // Gear NPT Berendsen + TRVP
        {
            integrator = new GearNPTBerendsen(h, system, Rsize, integratorName, Tfinal, tauT, Pfinal, tauP, *trvpCoeff, startWithShake);
        }
    }
    else if (ensembleLabel == 5) // NPT Nose-Hoover
    {
        system->SetPressForH(Pfinal);
        if (integratorLabel == -1)
        {
            integrator = new VerletNPTAlejandre(h, system, Tfinal, tauT, Pfinal, tauP, startWithShake);
        }
        else if (integratorLabel == 0) // Verlet NPT Nose-Hoover (MTTK formulation Trotter decomposition)
        {
            integrator = new VerletNPTNose(h, system, Tfinal, tauT, Pfinal, tauP, startWithShake);
        }
        else if (integratorLabel == 1) // Gear NPT Nose-Hoover
        {
            integrator = new GearNPTNose(h, system, Rsize, integratorName, Tfinal, tauT, Pfinal, tauP, *trvpCoeff, startWithShake);
        }
        else if (integratorLabel == 2) // Verlet NPT Nose-Hoover + TRVP (MACSIMUS style)
        {
            integrator = new VerletNPTNoseTRVP(h, system, *trvpCoeff, Tfinal, tauT, Pfinal, tauP, startWithShake);
        }
        else if (integratorLabel == 3) // Gear NPT Nose-Hoover + TRVP
        {
            integrator = new GearNPTNose(h, system, Rsize, integratorName, Tfinal, tauT, Pfinal, tauP, *trvpCoeff, startWithShake);
        }
        else if (integratorLabel == 4) // Verlet NVT Nose–Hoover (iterations)
        {
            integrator = new VerletNPTNoseIter(h, system, Tfinal, tauT, Pfinal, tauP, startWithShake, iterLabel);
        }
        else // Verlet NVT Nose–Hoover + TRVP (iterations)
        {
            integrator = new VerletNPTNoseIter(h, system, *trvpCoeff, Tfinal, tauT, Pfinal, tauP, startWithShake, iterLabel);
        }
    }
    else
    {
        return 40;
    }

    return 0;
}

// select integrator (parsing .control keyword "integrator")
int Simulation::SelectIntegrator(std::string integrator)
{
    size_t position;
    int size, trvp, iter; // size (k), order (m) (and trvp order(k in MACSIMUS)), number of iterations
    char auxc;
    int i;
    std::string trvpstring;
    std::string iterstring;
    Rsize = 0; // reset Rsize (to be updated in this function)

    // first word determines integrator, second TRVP (if any), third iter (if any)
    position = integrator.find(" ");
    if (position != std::string::npos)
    {
        trvpstring = integrator.substr(position + 1, std::string::npos);
        integrator.erase(position);
    }
    position = trvpstring.find(" ");
    if (position != std::string::npos)
    {
        iterstring = trvpstring.substr(position + 1, std::string::npos);
        trvpstring.erase(position);
    }

    // check if the second string is not iterstring (no trvpstring)
    if (iterstring.empty() && (trvpstring.find("iter") != std::string::npos))
    {
        iterstring = trvpstring;
        trvpstring.clear();
    }

    // integrator type init 0 (verlet)
    integratorLabel = 0;

    // if iter, parse iter
    iter = 0;
    if ((!iterstring.empty()) && (iterstring.find("iter") != std::string::npos))
    {
        if (iterstring.substr(4, std::string::npos).length() == 0)
        {
            print_warning(1, "Iteration demanded but no order given\n",
                          "    Defaulted to 3\n");
            iter = 3;
        }
        else
        {
            iter = std::stoi(iterstring.substr(4, std::string::npos));
        }
        if ((iter < 1) || (iter > 20))
        {
            print_warning(1, "Wrong number of iterations requested: " + std::to_string(iter) + "\n",
                          "    Defaulted to 3\n");
            iter = 3;
        }
        integratorLabel += 4; // integrator with iter has label +4 from default values
        iterLabel = iter;
    }

    // if trvp, parse trvp
    trvp = 0;
    if ((!trvpstring.empty()) && (trvpstring.find("trvp") != std::string::npos))
    {
        if (trvpstring.substr(4, std::string::npos).length() == 0)
        {
            print_warning(1, "TRVP demanded but no order given\n",
                          "    Defaulted to 2\n");
            trvp = 2;
        }
        else
        {
            trvp = std::stoi(trvpstring.substr(4, std::string::npos));
        }
        if ((trvp < 1) || (trvp > 4))
        {
            print_warning(1, "Wrong value of TRVP order: " + std::to_string(trvp) + "\n",
                          "    Defaulted to 2\n");
            trvp = 2;
        }
        Rsize += 1 + trvp;                // store history (trvp values) and predicted velocity (1 extra value)
        trvpCoeff = new Vector(trvp + 1); // TRVP coefficients for velocity prediction recursively (B in Kolafa+Lisal:JCTC(2011), trvp = k)
        trvpCoeff->operator()(0) = (2.0 * trvp + 1.0) / (trvp + 1.0);
        for (i = 0; i < trvp; i++)
        {
            trvpCoeff->operator()(i + 1) = -trvpCoeff->operator()(i) * (trvp - i) / (trvp + i + 2.0);
        }
        integratorLabel += 2; // integrator with trvp has label +2 from default values
    }
    else // empty trvpCoeff...
    {
        trvpCoeff = new Vector(0);
    }

    // verlet or leapfrog (meanwhile the same)
    if ((integrator.compare("leapfrog") == 0) || (integrator.compare("verlet") == 0))
    {
        integratorLabel += 0; // verlet integrator (even labels)
        strcpy(integratorName, "Verlet");
        Rsize += 4;
    }
    // alejandre (verlet for NPT nosé/mtk/hoover)
    else if (integrator.find("alej") == 0)
    {
        integratorLabel = -1;
        strcpy(integratorName, "Alejandre for NPT MTK");
        Rsize += 4;
    }
    // if not verlet and does not contain 'k' -> unknown integrator
    else if ((position = integrator.find("k")) != 0)
    {
        print_warning(1, "Unknown integrator selected: " + integrator + "\n",
                      "    Defaulted to Verlet\n");
        // default: verlet
        integratorLabel += 0;
        Rsize += 4;
        return -1;
    }
    // parse Gear integrator
    else
    {
        // check if Gear integrator name given (kXmXX)
        if (integrator.find("m") != position + 2)
        {
            print_warning(1, "Unknown integrator selected: " + integrator + "\n");
            // default: verlet
            integratorLabel += 0;
            Rsize += 4;
            return -1;
        }
        integratorLabel += 1; // gear integrator (odd labels)
        strcpy(integratorName, integrator.c_str());
        // method size (number after k)
        auxc = integrator.at(1);
        // convert number in auxc (ASCII) to integer (method size (not including TRVP))
        size = (int)auxc - '0'; // convert char to int
        // this number must be between 3 and 9
        if ((size < 3) || (size > 9))
        {
            size = 3; // defaulted to Verlet
            strcpy(integratorName, "k3m2e");
            print_warning(1, "Unknown integrator selected (wrong number after 'k'): " + integrator + "\n", "Deafaulted to Verlet.\n");
            // default: verlet
            integratorLabel += 0;
            Rsize += 4;
            return -1;
        }
        Rsize += size;
        // predictor and corrector are constructed in AbstractGearIntegrator constructor from integratorName
    }

    return 0;
}

// select ensemble (parsing .control keyword "ensemble")
int Simulation::SelectEnsemble(std::string ensemble)
{
    std::size_t position;
    std::string auxstring;

    if (ensemble.find("nve") != std::string::npos)
    {
        // NVE ensemble
        ensembleLabel = 0;
    }
    else if (ensemble.find("nvt berendsen") != std::string::npos)
    {
        // NVT Berendsen ensemble
        ensembleLabel = 1;
        position = ensemble.find("sen ");
        auxstring = ensemble.substr(position + 4, std::string::npos);
        if (auxstring.empty() || (position == std::string::npos))
        {
            print_warning(0, "Berendsen thermostat requested, but tau.T not specified\n",
                          "    Given: '" + ensemble + "', expected tau.T in [ps]\n");
            return 48;
        }
        tauT = atof(auxstring.c_str());
        if ((tauT < 0.0) || (tauT > 1000.0))
        {
            print_warning(1, "NVT Berendsen ensemble – thermostat tau.T outside limits (0, 1000)\n",
                          "    Obtained: '" + ensemble + "', converted to tau.T = " + std::to_string(tauT) + "\n",
                          "    Defaulted to 1 ps...\n");
            tauT = 1.0;
        }
    }
    else if (ensemble.find("ntvinit berendsen") != std::string::npos)
    {
        // NVT Berendsen + rescaling to final density ensemble
        ensembleLabel = 2;
        position = ensemble.find("sen ");
        auxstring = ensemble.substr(position + 4, std::string::npos);
        if (auxstring.empty() || (position == std::string::npos))
        {
            print_warning(0, "NTVinit ensemble requested, but tau.T and tau.rho not specified\n",
                          "    Given: '" + ensemble + "', expected tau.T and tau.rho in [ps]\n");
            return 48;
        }
        if (sscanf(auxstring.c_str(), "%lf %lf", &tauT, &taurho) < 2)
        {
            print_warning(0, "NTVinit ensemble requested, tau.rho not specified\n",
                          "    Given: '" + ensemble + "', expected tau.T and tau.rho in [ps]\n");
            return 48;
        }
        if ((tauT < 0.0) || (tauT > 1000.0))
        {
            print_warning(1, "NTVinit Berendsen ensemble – thermostat tau.T outside limits (0, 1000)\n",
                          "    Obtained: '" + ensemble + "', converted to tau.T = " + std::to_string(tauT) + "\n",
                          "    Defaulted to 1 ps...\n");
            tauT = 1.0;
        }
        if ((taurho < 0.0) || (taurho > 1.0e6))
        {
            print_warning(1, "NTVinit Berendsen ensemble – time to reach final density (tau[rho]) outside limits (0, 1e6)\n",
                          "    Obtained: '" + ensemble + "', converted to tau.rho = " + std::to_string(taurho) + "\n",
                          "    Defaulted to 1 ps...\n");
            taurho = 1.0;
        }
        // measure volume and density
        measuredQuant.insert(7);
        measuredQuant.insert(8);
    }
    else if (ensemble.find("nvt nose") != std::string::npos)
    {
        // NVT Nose–Hoover ensemble
        ensembleLabel = 3;
        position = ensemble.find("ose ");
        auxstring = ensemble.substr(position + 4, std::string::npos);
        if (auxstring.empty() || (position == std::string::npos))
        {
            print_warning(0, "Nose–Hoover thermostat requested, but tau.T not specified\n",
                          "    Given: '" + ensemble + "', expected tau.T in [ps]\n");
            return 48;
        }
        tauT = atof(auxstring.c_str());
        if ((tauT < 0.0) || (tauT > 1000.0))
        {
            print_warning(1, "NVT Nose–Hoover ensemble – thermostat tau.T outside limits (0, 1000)\n",
                          "    Obtained: '" + ensemble + "', converted to tau.T = " + std::to_string(tauT) + "\n",
                          "    Defaulted to 1 ps...\n");
            tauT = 1.0;
        }
    }
    else if (ensemble.find("npt berendsen") != std::string::npos)
    {
        // NPT Berendsen – Berendsen (friction) thermostat + friction barostat
        ensembleLabel = 4;
        position = ensemble.find("sen ");
        auxstring = ensemble.substr(position + 4, std::string::npos);
        if (auxstring.empty() || (position == std::string::npos))
        {
            print_warning(0, "NPT Berendsen ensemble requested, but tau.T and tau.P not specified\n",
                          "    Given: '" + ensemble + "', expected tau.T and tau.P in [ps]\n");
            return 48;
        }
        if (sscanf(auxstring.c_str(), "%lf %lf", &tauT, &tauP) < 2)
        {
            print_warning(0, "NPT Berendsen ensemble requested, tau.P not specified\n",
                          "    Given: '" + ensemble + "', expected tau.T and tau.P in [ps]\n");
            return 48;
        }
        if ((tauT < 0.0) || (tauT > 1000.0))
        {
            print_warning(1, "NPT Berendsen ensemble – thermostat tau.T outside limits (0, 1000)\n",
                          "    Obtained: '" + ensemble + "', converted to tau.T = " + std::to_string(tauT) + "\n",
                          "    Defaulted to 1 ps...\n");
            tauT = 1.0;
        }
        if ((tauP < 0.0) || (tauP > 1.0e6))
        {
            print_warning(1, "NPT Berendsen ensemble – barostat tau.P outside limits (0, 1e6)\n",
                          "    Obtained: '" + ensemble + "', converted to tau.P = " + std::to_string(tauP) + "\n",
                          "    Defaulted to 1 ps...\n");
            tauP = 1.0;
        }
        // measure volume and density and pressure
        measuredQuant.insert(7);
        measuredQuant.insert(8);
        measuredQuant.insert(15); // pressure form virial
    }
    else if ((ensemble.find("npt nose") != std::string::npos) || (ensemble.find("npt hoover") != std::string::npos) || (ensemble.find("npt mtk") != std::string::npos))
    {
        // NPT Nosé–Hoover in MTK fomulation
        ensembleLabel = 5;
        position = ensemble.find(" ", 4); // search from the first letter after the first space
        auxstring = ensemble.substr(position + 1, std::string::npos);
        if (auxstring.empty() || (position == std::string::npos))
        {
            print_warning(0, "NPT Nose-Hoover ensemble requested, but tau.T and tau.P not specified\n",
                          "    Given: '" + ensemble + "', expected tau.T and tau.P in [ps]\n");
            return 48;
        }
        if (sscanf(auxstring.c_str(), "%lf %lf", &tauT, &tauP) < 2)
        {
            print_warning(0, "NPT Nose-Hoover ensemble requested, tau.P not specified\n",
                          "    Given: '" + ensemble + "', expected tau.T and tau.P in [ps]\n");
            return 48;
        }
        if ((tauT < 0.0) || (tauT > 1000.0))
        {
            print_warning(1, "NPT Nose-Hoover ensemble – thermostat tau.T outside limits (0, 1000)\n",
                          "    Obtained: '" + ensemble + "', converted to tau.T = " + std::to_string(tauT) + "\n",
                          "    Defaulted to 1 ps...\n");
            tauT = 1.0;
        }
        if ((tauP < 0.0) || (tauP > 1.0e3))
        {
            print_warning(1, "NPT Nose-Hoover ensemble – barostat tau.P outside limits (0, 1000)\n",
                          "    Obtained: '" + ensemble + "', converted to tau.P = " + std::to_string(tauP) + "\n",
                          "    Defaulted to 1 ps...\n");
            tauP = 1.0;
        }
        // measure volume and density and pressure
        measuredQuant.insert(7);
        measuredQuant.insert(8);
        measuredQuant.insert(15); // pressure form virial
    }
    else
    {
        print_warning(1, "Unknown ensemble in .control file\n",
                      "    Given: '" + ensemble + "', defaulted to nve...\n");
        ensembleLabel = 0;
        return -1;
    }

    return 0;
}

// get Rsize for atom initialization
int Simulation::GetRsize() const
{
    return Rsize;
}

// get total proliferation number (nFold[0]*nFold[1]*nFold[2])
int Simulation::GetTotalProliferation() const
{
    return nFold[0] * nFold[1] * nFold[2];
}

// parse nfold (system multiplication) from .config
int Simulation::SaveNfold(std::string nfold)
{
    int i;
    sscanf(nfold.c_str(), "%d %d %d", &nFold[0], &nFold[1], &nFold[2]);
    for (i = 0; i < 3; i++)
    {
        if ((nFold[i] < 1) || (nFold[i] > 100000))
        {
            print_warning(0, "Error during .control reading\n",
                          "    wrong value of nfold[" + std::to_string(i) + "] = " + std::to_string(nFold[i]) + "\n");
            return 3;
        }
    }
    return 0;
}

// parse history directive in .control file
int Simulation::SelectHistory(std::string line)
{
    sscanf(line.c_str(), "%lf %d", &histFreq, &historyLevel);
    if ((histFreq < 0) || (histFreq > 9999999))
    {
        print_warning(0, "Wrong frequency of writing to .history file\n",
                      "    Option history got: " + line + "\n",
                      "    Parsed to histFreq = " + std::to_string(histFreq) + "\n");
        return 42;
    }
    if ((historyLevel < 0) || (historyLevel > 10))
    {
        print_warning(0, "Wrong history information level (keytraj) to .history file\n",
                      "    Option history got: " + line + "\n",
                      "    Parsed to historyLevel(keytraj) = " + std::to_string(historyLevel) + "\n");
        return 43;
    }
    return 0;
}

// parse init directive in .control file
int Simulation::SelectInit(std::string initialization)
{
    if (initialization.compare("continue") == 0)
    {
        return 0;
    }
    else if (initialization.compare("restart") == 0)
    {
        return 1;
    }
    else if (initialization.compare("scalevel") == 0)
    {
        return 2;
    }
    else if (initialization.compare("randomvel") == 0)
    {
        return 3;
    }
    else if (initialization.compare("random") == 0)
    {
        return 4;
    }
    else
    {
        return -1;
    }
}

// parse pairlist directive in .control file
int Simulation::SelectPairList(std::string pairlist)
{
    size_t position;
    std::string auxstring;
    std::string label;
    char auxch[30];

    position = pairlist.find(" ");
    label = pairlist.substr(0, position);
    auxstring = pairlist.substr(position + 1, std::string::npos);
    if (label.compare("all") == 0)
    {
        pairlistLabel = 0; // all pairs
        return 0;
    }
    else if (label.compare("verlet") == 0)
    {
        pairlistLabel = 1;
        if (sscanf(auxstring.c_str(), "%lf %s", &list_save, auxch) < 1)
        {
            list_save = 1.2;
            print_warning(1, "Verlet list requested, but no save margin given.\n",
                          "Defaulted to 1.2!\n");
        }
        if ((list_save < 1.0) || (list_save > 2.0))
        {
            list_save = 1.2;
            print_warning(1, "Verlet list requested, but incorrect save margin given.\n",
                          "Given: " + std::to_string(list_save) + ", expected value between 1.0 and 2.0\n",
                          "Defaulted to 1.2!\n");
        }
        return 0;
    }
    return -1;
}

// parse 'measure' directive in .control file
int Simulation::AddMeasuredQuantities(std::string quant)
{
    int error = 0;
    // Used to split string around spaces.
    std::istringstream ss(quant);

    std::string word; // for storing each word

    bool del = false;
    int quantindex;

    // Traverse through all words
    // while loop till we get
    // strings to store in string word
    while (ss >> word)
    {
        if (word[0] == '-')
        {
            del = true;       // quantity will be deleted instead of added
            word.erase(0, 1); // to recognize quantity
        }
        quantindex = quant::idx_from_code(word);
        if (quantindex < 0)
        {
            print_warning(1, "Unknown quantity in 'measure' directive: " + word + "\n");
            error++;
        }
        else if (del)
        {
            measuredQuant.erase(quantindex + 1);
            del = false;
        }
        else
        {
            measuredQuant.insert(quantindex + 1);
        }
    }
    return error;
}

// get TRVP size for System initialization (needed in System::ReadConfig)
int Simulation::GetTRVPsize() const
{
    return (trvpCoeff == NULL) ? 0 : (trvpCoeff->GetSize() - 1);
}

// select class to handle electrostatics
int Simulation::SelectElectrostatics(std::string line)
{
    int numbersread = 0;
    double param1, param2, param3;
    char aux[MAX_COMMENT];
    int localerror = 0;

    numbersread = sscanf(line.c_str(), "%s %lf %lf %lf", aux, &param1, &param2, &param3) - 1;
    if (strstr(aux, "no") != nullptr)
    {
        elstatlabel = 0; // no electrostatics
    }
    else if (strstr(aux, "cut") != nullptr)
    {
        elstatlabel = 1; // cut-and-shifted electrostatics
        if (numbersread < 1)
        {
            localerror = 57;
            error = 57;
            print_warning(0, "Cutoff for cutoff electrostatics must be given!!!\n",
                          "    Given: " + line + ", expected cutoff CUTOFF_VALUE [ALPHA_VALUE]\n");
        }
        else if (numbersread < 2)
        {
            elcutoff = param1;
            elalpha = 0.7; // default value of alpha
        }
        else
        {
            elcutoff = param1;
            elalpha = param2;
        }
    }
    else if (strstr(aux, "ewalds") != nullptr)
    {
        elstatlabel = 3; // basic Ewald summation using splines for erfc interpolation
        if (numbersread < 1)
        {
            localerror = 57;
            error = 57;
            print_warning(0, "Cutoff for Ewald electrostatics must be given!!!\n",
                          "    Given: " + line + ", expected ewalds CUTOFF_VALUE [ALPHA_VALUE KAPPA_VALUE]\n");
        }
        else if (numbersread < 2)
        {
            elcutoff = param1;
            elalpha = M_PI / elcutoff; // default value for alpha
            elkappa = M_PI / elcutoff; // default value for kappa
        }
        else if (numbersread < 3)
        {
            elcutoff = param1;
            elalpha = param2;
            elkappa = M_PI / elcutoff; // default value for kappa
        }
        else
        {
            elcutoff = param1;
            elalpha = param2;
            elkappa = param3;
        }
    }
    else if (strstr(aux, "ewald") != nullptr)
    {
        elstatlabel = 2; // basic Ewald summation
        if (numbersread < 1)
        {
            localerror = 57;
            error = 57;
            print_warning(0, "Cutoff for Ewald electrostatics must be given!!!\n",
                          "    Given: " + line + ", expected ewald CUTOFF_VALUE [ALPHA_VALUE KAPPA_VALUE]\n");
        }
        else if (numbersread < 2)
        {
            elcutoff = param1;
            elalpha = M_PI / elcutoff; // default value for alpha
            elkappa = M_PI / elcutoff; // default value for kappa
        }
        else if (numbersread < 3)
        {
            elcutoff = param1;
            elalpha = param2;
            elkappa = M_PI / elcutoff; // default value for kappa
        }
        else
        {
            elcutoff = param1;
            elalpha = param2;
            elkappa = param3;
        }
    }
    else
    {
        elstatlabel = 0; // unknown electrostatics -> no electrostatics
        localerror = 58;
        error = 58;
        print_warning(0, "Unknown version of electrostatics in .control file!!!\n",
                      "    Given: " + line + ", expected 'elstat no' or 'elstat cutoff ...'\n");
    }

    return localerror;
}

// set scaling for 1–4 interactions
int Simulation::SetScaling14(std::string scalingline)
{
    int numbersread;
    double param1, param2;
    int err = 0;
    numbersread = sscanf(scalingline.c_str(), "%lf %lf", &param1, &param2);

    if (numbersread < 1)
    {
        err = 68;
        print_warning(0, "At least one number must be given in scale14!!!\n",
                      "    Given: '" + scalingline + "', expected 'scale DISPERSE14SCALING [ELSTAT14SCALING]\n");
    }
    else if (numbersread < 2)
    {
        lj14 = param1;
        elstat14 = lj14;
    }
    else
    {
        lj14 = param1;
        elstat14 = param2;
    }

    return err;
}

// prit resume to .prt file
int Simulation::PrintInfo(std::ofstream &fprt, std::string engunit)
{
    std::set<int>::const_iterator it;
    std::string qname;
    fprt << "------------------------------------------------------------------\n\n";
    fprt << "\n# Simulation\n\n";
    fprt << "## Summary\n\n";
    // timestep, number of steps
    MDTable<std::string, double> summary({"Simulation property", "value"});
    summary.addRow("Simulaton timestep [ps]", h);
    summary.addRow("Number of steps", noSteps);
    summary.addRow("Total simulation time [ps]", h * noSteps);
    if (freqMeasPlbHistDump[0] > 0)
    {
        summary.addRow("Number of steps between two consecutive measuremet (`.cpa`)", freqMeasPlbHistDump[0]);
    }
    if (freqMeasPlbHistDump[1] > 0)
    {
        summary.addRow("Number of steps between two playback frames (`.plb`)", freqMeasPlbHistDump[1]);
    }
    if (freqMeasPlbHistDump[2] > 0)
    {
        summary.addRow("Number of steps between two trajectory frames (`.history`)", freqMeasPlbHistDump[2]);
        summary.addRow("Trajectory level in `.history`\\*", historyLevel);
    }
    if (freqMeasPlbHistDump[3] > 0)
    {
        summary.addRow("Number of steps between two dumps (backups) (`.dump`)", freqMeasPlbHistDump[3]);
    }
    summary.print(fprt);
    if (freqMeasPlbHistDump[2] > 0)
    {
        fprt << "\\*) 0 = positions, 1 = positions + velocities,...\n\n";
    }
    else
    {
        fprt << "\n";
    }

    // // filenames
    // fprintf(prtfile, ".config file name: %s\n", configName);
    // fprintf(prtfile, ".field file name: %s\n", fieldName);

    // measured quantities
    if (measurementFreq != 0)
    {
        fprt << "## Measured quantities (`.cpa`)\n\n";
        MDTable<std::string, std::string> measurements({"Quantity [unit]", "Description"});

        for (it = measuredQuant.begin(); it != measuredQuant.end(); it++)
        {
            qname = quant::name[(*it) - 1];
            if (quant::iseng[(*it) - 1])
                qname.append(engunit);
            measurements.addRow(qname, quant::description[(*it) - 1]);
        }
        measurements.print(fprt);
        fprt << "\n";
    }

    // ensemble and related values (temperature, thermostat constant...)
    fprt << "## Statistical ensemble\n\n";
    switch (ensembleLabel)
    {
    case 0:
        fprt << "- **Microcanonical ensemble (NVE)**\n\n";
        break;
    case 1:
        fprt << "- **Berendsen (friction) thermostat (NVT, not precisely)**\n";
        fprt << "- Target **temperature**: " << Tfinal << " K\n";
        fprt << "- Thermostat constant: " << tauT << " ps\n\n";
        break;
    case 2:
        fprt << "- **Berendsen (friction) thermostat with volume changing to desired density** (*NTVinit*)\n";
        fprt << "- Target **temperature**: " << Tfinal << " K\n";
        fprt << "- Thermostat constant: " << tauT << " ps\n";
        fprt << "- Target **density**: " << rhofinal * U_DENSITY_KGM3 << " kg m^(-3)\n";
        fprt << "- Time to reach target density: " << taurho << " ps\n\n";
        break;
    case 3:
        fprt << "- **Nosé–Hoover thermostat (NVT, extra degree of freedom)**\n";
        fprt << "- Target **temperature**: " << Tfinal << " K\n";
        fprt << "- Thermostat constant: " << tauT << " ps\n\n";
        break;
    case 4:
        fprt << "- **Berendsen (friction) thermostat & friction barostat (NPT, not precisely)**\n";
        fprt << "- Target **temperature**: " << Tfinal << " K\n";
        fprt << "- Thermostat constant: " << tauT << " ps\n\n";
        fprt << "- Target pressure: " << Pfinal * U_PRESSURE_PA << " Pa\n";
        fprt << "- Barostat constant: " << tauP << " ps\n\n";
        break;
    case 5:
        fprt << "- **Nosé–Hoover themostat and barostat based on extended Lagrangian (MTK formulation) (NPT)**\n";
        fprt << "- Target **temperature**: " << Tfinal << " K\n";
        fprt << "- Thermostat constant: " << tauT << " ps\n\n";
        fprt << "- Target pressure: " << Pfinal * U_PRESSURE_PA << " Pa\n";
        fprt << "- Barostat constant: " << tauP << " ps\n\n";
        break;
    }

    // integrator
    fprt << "## Integrator\n\n";
    switch (integratorLabel)
    {
    case -1:
        fprt << "- **Verlet – Alejandre's version of Trotter decomposition**\n\n";
        break;
    case 0:
        fprt << "- **Verlet/Leapfrog**\n\n";
        break;
    case 1:
        fprt << "- **" << integratorName << "**\n\n";
        break;
    case 2:
        fprt << "- **Verlet + TimeReversibleVelocityPredictor(TRVP)**\n\n";
        break;
    case 3:
        fprt << "- **" << integratorName << "**\n\n";
        break;
    case 4:
        fprt << "- **Verlet + self-consistent iterations**\n\n";
        break;
    case 6:
        fprt << "- **Verlet + TimeReversibleVelocityPredictor(TRVP) + self-consistent iterations**\n\n";
        break;
    }

    // SHAKE info
    if (system->noConstrTotal > 0)
    {
        fprt << "## SHAKE parameters\n\n";
        fprt << "- **relative tolerance**: " << shakeTol << "\n";
        fprt << "- overrelaxation parameter: " << shakeOmega << "\n\n";
    }

    // initialization
    fprt << "## Initialization\n\n";
    switch (initLabel)
    {
    case 0:
        fprt << "`init continue`\n";
        fprt << " -`.config` file read and **velocities NOT rescaled**\n";
        fprt << "- **measurements** in `.cpa` file, playback in `.plb` file (if present) and history in `.history file` (if present)\n";
        fprt << "**will be appended** to the files from previous simulation\n\n";
        fprt << "**WARNING**: No consistency checks done (number of molecules, measured quantities etc. not compared)\n\n";
        break;
    case 1:
        fprt << "`init restart`\n";
        fprt << "- `.config` file read and **velocities NOT rescaled**\n";
        fprt << "- **measurements**, playback and history **will start anew** (default option)\n\n";
        break;
    case 2:
        fprt << "`init scalevel`\n";
        fprt << "- `.config` file read and **velocities RESCALED (by factor " << vel_scaling << ") to desired temperature " << Tfinal << " K**\n";
        fprt << "- **measurements**, playback and history **will start anew**\n\n";
        break;
    case 3:
        fprt << "`init randomvel`\n";
        fprt << "- `.config` file read for positions\n";
        fprt << "- **velocities chosen randomly** from uniform distribution **and scaled (by factor " << vel_scaling << ") to desired temperature " << Tfinal << " K\n";
        fprt << "- **measurements**, playback and history **will start anew**\n\n";
        break;
    default:
        fprt << "`init random`\n\n";
        fprt << "**WARNING: NOT IMPLEMENTED YET...**\n\n";
        break;
    }

    return 0;
}

// get names of measured quantities
std::vector<std::string> Simulation::GetMeasuredQuant(std::string engunit) const
{
    std::string qname;
    std::vector<std::string> vec;
    int i;
    for (auto it = measuredQuant.begin(); it != measuredQuant.end(); it++)
    {
        qname = quant::name[(*it) - 1];
        if (quant::iseng[(*it) - 1])
            qname.append(engunit);
        if ((*it) == 38) // MSD
        {
            for (i = 0; i < system->msd->GetNumberOfTypes(); i++)
            {
                vec.push_back(qname);
            }
        }
        else
        {
            vec.push_back(qname);
        }
    }
    return vec;
}

int Simulation::SelectMixingRule(std::string line)
{
    if (line.compare("lorentz-berthelot") == 0)
    {
        mixingrule = 1;
        return 0;
    }
    else
    {
        return -1;
    }

    return 0;
}

int Simulation::SelectSplines(std::string splinesline)
{
    char type1[MAX_COMMENT], type2[MAX_COMMENT];
    int grid1, grid2;
    int numbersread = sscanf(splinesline.c_str(), "%s %d %s %d", type1, &grid1, type2, &grid2);

    if (numbersread < 2)
    {
        print_warning(1, "Too few parameters for keyword `splines`.", "Given: " + splinesline + "Expected at least two parameters (type and gridsize).");
    }
    else
    {
        if (strstr(type1, "hyperbolic,") != nullptr)
        {
            splineslabel1 = hyperbolic;
        }
        else if (strstr(type1, "quadratic") != nullptr)
        {
            splineslabel1 = quadratic;
        }
        else if (strstr(type1, "natural") != nullptr)
        {
            splineslabel1 = natural;
        }
        else if (strstr(type1, "hermite") != nullptr)
        {
            splineslabel1 = hermite;
        }
        else if (strstr(type1, "linear") != nullptr)
        {
            splineslabel1 = linear;
        }
        else if (strstr(type1, "linear3") != nullptr)
        {
            splineslabel1 = linear3;
        }
        else if (strstr(type1, "linear4") != nullptr)
        {
            splineslabel1 = linear4;
        }
        gridsize1 = grid1;
    }
    if (numbersread == 4)
    {
        if (strstr(type2, "hyperbolic,") != nullptr)
        {
            splineslabel2 = hyperbolic;
        }
        else if (strstr(type2, "quadratic") != nullptr)
        {
            splineslabel2 = quadratic;
        }
        else if (strstr(type2, "natural") != nullptr)
        {
            splineslabel2 = natural;
        }
        else if (strstr(type2, "hermite") != nullptr)
        {
            splineslabel2 = hermite;
        }
        else if (strstr(type2, "linear") != nullptr)
        {
            splineslabel2 = linear;
        }
        else if (strstr(type2, "linear3") != nullptr)
        {
            splineslabel2 = linear3;
        }
        else if (strstr(type2, "linear4") != nullptr)
        {
            splineslabel2 = linear4;
        }
        gridsize2 = grid2;
    }
    return 0;
}
