/*
 * Main program for testing simulations with constrained bonds
 *
 * Author: JJ
 * Date: Jan 2021
 *
 */

#include <iostream>
#include <fstream>
// #include <string>
#include <cstring>
#include <sys/time.h>

#include "Vector.hpp"
#include "Matrix.hpp"
#include "ConstraintBond.hpp"
#include "Atom.hpp"
#include "Molecule.hpp"
#include "SimulatedSystem.hpp"
#include "file_openclose.hpp"
#include "Simulation.hpp"
#include "version.hpp"
#include "general_utils.hpp"

#ifndef VERLET
#define VERLET 1
#endif

#ifndef ERRORSPATH
#define ERRORSPATH "/home/janekj/Programs/simul++/documentation/errors.txt"
#endif

#ifdef PARALLEL
#include <omp.h>
int thread_count = 1;
#define MAX_THREAD_COUNT 10
#endif

void printHelp();
void printHeader(timeval start, char *simname, std::ofstream &fprt);
void printError(int errnumber, std::ofstream &fprt);
void printConclusion(timeval end, double duration, std::ofstream &fprt);

int main(int argc, char **argv)
{
    char simname[MAX_COMMENT];
    char prtname[MAX_COMMENT];
    char cpaname[MAX_COMMENT];
    char prtback[MAX_COMMENT];
    int error = 0;
    int i;
    bool simnameset = false;         // simname given (if not print help and exit)
    timeval start_time, final_time;  // start and final time
    gettimeofday(&start_time, NULL); // save start time
    double computing_time = 0.0;     // total computing (CPU) time

    std::ofstream fprt;

    // if no argument or -h as an argument, then print help
    if ((argc < 2) || ((argv[1][0] == '-') && (argv[1][1]) == 'h'))
    {
        printHelp();
        return -1;
    }

    // read commandline arguments
    for (i = 1; i < argc; i++)
    {
        if (argv[i][0] == '-') // options
        {
            switch (argv[i][1])
            {
#ifdef PARALLEL
            case 'p': // parallel: number of threads
                thread_count = std::stoi(argv[++i]);
                if ((thread_count > MAX_THREAD_COUNT) || (thread_count < 1))
                {
                    print_warning(1, "Wrong number of threads for OpenMP.\n",
                                  "    Given: " + std::to_string(argv[i][0]) + ", converted to: " + std::to_string(thread_count) + "\n",
                                  "    Defaulted to thread_count = 1\n");
                    thread_count = 1;
                }
                break;
#endif
            case 'h': // help
                printHelp();
                return -1;
            case 'v': // version
                print_warning(2, "This is `simul++` v" + std::string(VERSION) + "\n");
                break;
            default: // unknown option
                print_warning(1, "Skipping unknown option: " + std::to_string(argv[i][1]) + "\n");
                break;
            }
        }
        else // simname
        {
            // get simname (name of .control and default name of .field and .config files)
            strcpy(simname, argv[i]);
            simnameset = true;
        }
    }

    if (!simnameset)
    {
        print_warning(0, "Simulation name not given, cannot choose .control file...\n");
        return 1;
    }

    // open protocol file
    strcpy(prtname, simname);
    strcat(prtname, ".prt");
    fprt.open(prtname, std::ios::in);
    if (fprt)
    {
        print_warning(2, "File " + std::string(prtname) + " exists, backup " + std::string(prtname) + "~ will be made\n");
        strcpy(prtback, prtname);
        strcat(prtback, "~");
        std::rename(prtname, prtback);
        fprt.close();
    }
    fprt.open(prtname, std::ios::out);
    if (!fprt)
    {
        print_warning(0, " Protocol file cannot be opened!\n",
                      "    Simulation ended unsuccessfully...\n");
        error = 44;
        return error;
    }

    printHeader(start_time, simname, fprt); // program version, simulation name

    Simulation simulation(simname, fprt); // read .control, provide Rsize and total proliferation number; print .control to .prt
    error = simulation.error;
    if (error != 0)
    {
        printError(error, fprt);
        fprt.close();
        return error;
    }
    SimulatedSystem system(simname, simulation.configName,
                           simulation.fieldName, simulation.LJcutoff,
                           simulation.GetTotalProliferation(),
                           simulation.GetRsize(), simulation.GetTRVPsize()); // read .field and .config
    error = system.error;
    if (error != 0)
    {
        printError(error, fprt);
        fprt.close();
        return error;
    }
    simulation.SetSystem(&system);   // set system to simulation
    error = simulation.Initialize(); // choose correct integrator and prepare simulation
    if (error != 0)
    {
        printError(error, fprt);
        fprt.close();
        return error;
    }
    // print system info to .prt
    std::string engunit;
    simulation.system->PrintInfo(fprt, engunit);

    // print simulation info to .prt
    simulation.PrintInfo(fprt, engunit);

    error = simulation.Execute(); // execute simulation and print measurements
    if (error > 0)
    {
        printError(error, fprt);
        fprt.close();
        return error;
    }

    // print statistics to .prt (from .cpa) (.cpa would have to be closed before by SimulatedSystem destructor...)
    simulation.system->CloseFiles();
    strcpy(cpaname, simname);
    strcat(cpaname, ".cpa");
    basic_statistics(fprt, cpaname, simulation.GetMeasuredQuant(engunit));

    if (error == -999) // simulation stopped by .stp file
    {
        fprt << "**Simulation stopped prematurely by `.stp` file!**\n";
        fprt << "To continue it use `init continue` in `.control` file.\n";
    }

    gettimeofday(&final_time, NULL); // save final time
    computing_time = (final_time.tv_sec - start_time.tv_sec) + ((double)final_time.tv_usec - (double)start_time.tv_usec + 500) / 1000000;
    printConclusion(final_time, computing_time, fprt);

    fprt.close();

    return 0;
}

// print help - to be extended
void printHelp()
{
    std::cout << TEXTBF << COL_ESCAPE_IN << fg_blue << COL_ESCAPE_OUT << "simul++" << COL_RESET << " by JJ version " << VERSION << " compiled on " << DATE << "\n";
    std::cout << "    velocity version for Verlet-type integrators: VERLET=" << VERLET << "\n";
#ifdef PARALLEL
    std::cout << "    parallelization using OpenMP supported, use -p NUMBER_OF_THREADS to run simul++ in parallel\n";
#endif
    std::cout << TEXTBF << "    usage:" << COL_RESET << " simul++ CONTROLFILENAME\n";
}

// print protocol header
void printHeader(timeval start, char *simname, std::ofstream &fprt)
{
    std::string currtime(ctime((time_t *)&(start.tv_sec)));
    currtime.replace(currtime.find("\n"), std::string::npos, "\0");

    fprt << "# " << simname << "\n";
    fprt << "## `simul++` - experimental package for MD simulations by JJ\n\n";
    fprt << "- `simul++` v" << VERSION << "\n- compiled on " << DATE << "\n";
#ifdef PARALLEL
    fprt << "- parallel version (OpenMP); number of threads: " << thread_count << "\n\n";
#endif
    fprt << "-------------------------------------------------------------\n\n";
    fprt << "Simulation **started at " << currtime << "**\n\n";
}

// print final info to prt file
void printConclusion(timeval end, double duration, std::ofstream &fprt)
{
    std::string endtime(ctime((time_t *)&(end.tv_sec)));
    endtime.replace(endtime.find("\n"), std::string::npos, "\0");
    fprt << "\n--------------------------------------------------------------\n\n";
    fprt << "Simulation **ended** successfully **at: " << endtime << "**\n\n";
    fprt << "Total computing time: " << duration << " s\n";
}

// print error description to protocol file
void printError(int errnumber, std::ofstream &fprt)
{
    FILE *errfile;
    char errfilename[] = ERRORSPATH;
    char readmode[] = "r";
    int err;
    char message[200];
    char line[210];

    if (my_fopen_r(&errfile, errfilename, readmode) == 0)
    {
        while (fgets(line, 210, errfile) != NULL)
        {
            sscanf(line, "%d %s", &err, message);
            if (err == errnumber)
            {
                fprt << "**ERROR:** " << line << "Cannot continue...\n";
                break;
            }
        }
        my_fclose(&errfile, errfilename);
    }
    else
    {
        fprt << "**ERROR:** " << errnumber << " occured, but cannot find the reference in `documentation/errors.txt`";
    }
}
