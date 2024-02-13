/*
 * Merge simul++ (or DL_POLY) .config files
 * originally written for simul++ simulation package
 * Author: JJ
 * Date: Jan 2022
 */

/*
 * Features and options:
 *     load each .config file to SimulatedSystem class (.field needed – need not be a problem)
 *     move configurations (already implemented in Simulated System)
 *     print merged config
 *     1D, 2D and 3D mixing (linear, square, cubic) – half or mixed.
 *     if non-equal box-sizes, the smaller can be expanded or not (only centered)
 *     velocity can be added to all atoms in one .config file, impact momentum (one number)
 *     larger box-size can be specified
 *     reduce the number of configlevel to the least of the original files (max 2, positions and velocities, because it does not make sense to continue directly the simulation and forces are different after merge)
 *     in case of velocities missing, they can be assigned according to the selected temperature
 * For options list see function print_help()...
 */

#include <cstring>
#include <iostream>
#include <cmath>
#include "file_openclose.hpp"
#include "SimulatedSystem.hpp"
#include "Vector.hpp"
#include "general_utils.hpp"

#define MAXFILENAME 100

#ifdef PARALLEL
int thread_count = 1;
#endif

void print_help();

int main(int argc, char **argv)
{
    int i, j;
    double momentum = 0.0;                    // impact momentum, default 0
    bool mix = false;                         // default -H HALF
    bool expand = false;                      // default -C CENTER
    int direction[3] = {-1, -1, -1};          // dir: z, zy, zyx
    int dimension = 1;                        // default -D: 1 (linear merge)
    int no_config = 0;                        // number of .config files
    double temperature = 0.0;                 // if -T given, random velocities are assigned
    double new_box_size[3] = {0.0, 0.0, 0.0}; // if -L given, bigger box
    double finalbox[3] = {0.0, 0.0, 0.0};
    char config1[MAXFILENAME];
    config1[0] = '\0';
    char config2[MAXFILENAME];
    config2[0] = '\0';
    char field1[MAXFILENAME];
    field1[0] = '\0';
    char field2[MAXFILENAME];
    field2[0] = '\0';
    char output[MAXFILENAME];
    SimulatedSystem *system1, *system2, *system3;
    Vector scaling1(3), scaling2(3);

    system3 = nullptr; // just to avoid compiler warnings

    strcpy(output, "output");

    // if too few arguments given, print help...
    if ((argc < 2) || ((argv[1][0] == '-') && (argv[1][1]) == 'h'))
    {
        print_help();
        return -1;
    }

    // read commandline arguments
    for (i = 1; i < argc; i++)
    {
        if (argv[i][0] == '-') // options
        {
            switch (argv[i][1])
            {
            case 'o':
                if (argc <= (i + 1))
                {
                    print_warning(0, "Unexpected end of input (option -o value expected)\n");
                    return 5;
                }
                strcpy(output, argv[i + 1]);
                i++;
                break;
            case 'f':
                if (argc <= (i + no_config))
                {
                    print_warning(0, "Unexpected end of input (option -f values expected)\n");
                    return 5;
                }
                if (no_config > 0)
                {
                    strcpy(field1, argv[i + 1]);
                    i++;
                }
                if (no_config > 1)
                {
                    strcpy(field2, argv[i + 1]);
                    i++;
                }
                break;
            case 'T':
                if (argc <= (i + 1))
                {
                    print_warning(0, "Unexpected end of input (option -T value expected)\n");
                    return 5;
                }
                temperature = atof(argv[i + 1]);
                i++;
                break;
            case 'D':
                if (argc <= (i + 1))
                {
                    print_warning(0, "Unexpected end of input (option -D value expected)\n");
                    return 5;
                }
                dimension = atoi(argv[i + 1]);
                if ((dimension < 1) || (dimension > 3))
                {
                    print_warning(0, "Wrong dimension of merge given (option -D)!!!\n",
                                  "    Given " + std::string(argv[i + 1]) + ", expected 1, 2 or 3...\n",
                                  "    Cannot continue...\n");
                    return 1;
                }
                i++;
                break;
            case 'd':
                if (argc <= (i + 1))
                {
                    print_warning(0, "Unexpected end of input (option -d value expected)\n");
                    return 5;
                }
                for (j = 0; j < (int)strlen(argv[i + 1]); j++)
                {
                    direction[j] = (int)(argv[i + 1][j] - 'x');
                    if ((direction[j] < 0) || (direction[j] > 2))
                    {
                        print_warning(0, "Wrong direction of merge given (option -d)!!!\n",
                                      "    Given " + std::string(argv[i + 1]) + ", expected x, y or z\n",
                                      "    Cannot continue...\n");
                        return 2;
                    }
                }
                i++;
                break;
            case 'L':
                if (argc <= (i + 1))
                {
                    print_warning(0, "Unexpected end of input (option -L value expected)\n");
                    return 5;
                }
                new_box_size[0] = atof(argv[i + 1]);
                i++;
                j = 1;
                while (((i + 1) < argc) && (isdigit(argv[i + 1][0]) != 0))
                {
                    new_box_size[j] = atof(argv[i + 1]);
                    j++;
                    i++;
                }
                break;
            case 'p':
                if (argc <= (i + 1))
                {
                    print_warning(0, "Unexpected end of input (option -p value expected)\n");
                    return 5;
                }
                momentum = atof(argv[i + 1]);
                i++;
                break;
            case 'H':
                mix = false;
                break;
            case 'M':
                mix = true;
                break;
            case 'E':
                expand = true;
                break;
            case 'C':
                expand = false;
                break;
            case 'h':
                print_help();
                return -1;
            default:
                print_warning(1, "Skipping unknown option: " + std::to_string(argv[i][1]) + "\n");
                break;
            }
        }
        else // .config files
        {
            /* assign .config file names */
            if (no_config < 1)
            {
                strcpy(config1, argv[i]);
                no_config++;
            }
            else if (no_config < 2)
            {
                strcpy(config2, argv[i]);
                no_config++;
            }
            else
            {
                print_warning(1, "Too many .config files given accepting only the first two...\n",
                              "    Unused: " + std::string(argv[i]) + "\n");
            }
        }
    }

    // check if at least one .config file was specified
    if (no_config < 2)
    {
        print_warning(0, "No or only one .config file given, nothing to do...\n");
        return 3;
    }

    // if .field files not given they must be derived from .config files
    if (strstr(config1, ".config") == NULL)
    {
        if (strcmp(field1, "") == 0)
        {
            strcpy(field1, config1);
            strcat(field1, ".field");
        }
        strcat(config1, ".config");
    }
    else
    {
        if (strcmp(field1, "") == 0)
        {
            print_warning(0, ".config file specified with extension, but .field file not given\n",
                          "    Cannot continue without .field file...\n");
            return 4;
        }
    }

    if (strstr(config2, ".config") == NULL)
    {
        if (strlen(field2) == 0)
        {
            strcpy(field2, config2);
            strcat(field2, ".field");
        }
        strcat(config2, ".config");
    }
    else
    {
        {
            print_warning(0, ".config file specified with extension, but .field file not given\n",
                          "    Cannot continue without .field file...\n");
            return 4;
        }
    }

    // directions must be specified up to dimension
    for (i = 0; i < dimension; i++)
    {
        if (direction[i] == -1)
        {
            for (j = 2; j >= 0; j--)
            {
                if ((direction[0] != j) && (direction[1] != j) && (direction[2] != j))
                {
                    direction[i] = j;
                    break;
                }
            }
        }
    }

    // temporal check settings: in future summary what will be done...
    std::cout << "dimension = " << dimension << "\n";
    std::cout << "config1 = " << config1 << "\n";
    std::cout << "config2 = " << config2 << "\n";
    std::cout << "field1 = " << field1 << "\n";
    std::cout << "field2 = " << field2 << "\n";

    // std::cout << "directions: " << direction[0] << ", " << direction[1] << ", " << direction[2] << "\n";

    // now load configurations from files to class SimulatedSystem
    // specify simulation names and load configurations
    char simname1[MAXFILENAME];
    char simname2[MAXFILENAME];

    strcpy(simname1, output);
    strcpy(simname2, output);
    strcat(simname1, "1");
    strcat(simname2, "2");

    system1 = new SimulatedSystem(simname1, config1, field1, 0.0, 1, 3);
    system2 = new SimulatedSystem(simname2, config2, field2, 0.0, 1, 3);

    // memory allocation check (existence check)
    if ((system1 == nullptr) || (system2 == nullptr))
    {
        print_warning(0, "Error during input systems initialization ( memory allocation??? )\n",
                      "    This should never happen. Exiting...\n");
        return 6;
    }

    // assign random velocities if needed (option -T given) and set temperature to desired value
    if (temperature > 0.0)
    {
        system1->RandomVelocities();
        system2->RandomVelocities();
        system1->SetTemperature(temperature); // work with v not h*v in R
        system2->SetTemperature(temperature);
    }

    // assign impact velocities if needed (option -p) and not mixing configurations
    if ((momentum > 0.0) && (!mix))
    {
        system1->AddVelocity(momentum * 0.5 / system1->totalMass, direction[0]);
        system2->AddVelocity(-momentum * 0.5 / system2->totalMass, direction[0]);
    }

    // variables for scaling
    std::vector<int> positionsCoords;
    positionsCoords.push_back(0);

    // each direction
    for (j = 0; j < dimension; j++)
    {
        // box-size balancing
        for (i = 0; i < 3; i++)
        {
            if (direction[j] == i)
            {
                finalbox[i] = fmax(system1->GetBox(i) + system2->GetBox(i), new_box_size[j]);
                scaling1[i] = finalbox[i] / (system1->GetBox(i) + system2->GetBox(i));
                scaling2[i] = finalbox[i] / (system1->GetBox(i) + system2->GetBox(i));
            }
            else
            {
                finalbox[i] = fmax(system1->GetBox(i), system2->GetBox(i));
                scaling1[i] = finalbox[i] / system1->GetBox(i);
                scaling2[i] = finalbox[i] / system2->GetBox(i);
            }
        }

        // expand box or only change box size (E/C options)
        if (expand) // option -E
        {
            system1->RescaleBox3D(scaling1, positionsCoords);
            system2->RescaleBox3D(scaling2, positionsCoords);
        }
        else // default, option -C
        {
            for (i = 0; i < 3; i++)
            {
                system1->SetBox(scaling1[i] * system1->GetBox(i), i);
                system2->SetBox(scaling2[i] * system2->GetBox(i), i);
            }
        }

        // invert 2nd system (if needed; option -M MIX)
        if (mix && (j > 1))
        {
            system2->ReflectBox(positionsCoords, direction[1]);
        }
        else if (mix)
        {
            system2->InvertBox(positionsCoords);
        }

        // merge two systems
        if (system3 != nullptr)
        {
            delete system3;
            system3 = nullptr;
        }
        system3 = new SimulatedSystem(*system1, *system2, direction[j]);

        delete system1;
        delete system2;
        system1 = new SimulatedSystem(*system3);
        system2 = new SimulatedSystem(*system3);
    }

    delete system1;
    delete system2;

    // print final configuration (timestep is 1.0 not to have impact on the configuration velocity)
    if (system3 != nullptr) // just to be on a save side
    {
        system3->ApplyPeriodicBC(positionsCoords);
        system3->PrintConfig(1, output, 1.0);
        delete system3;
        system3 = nullptr;
    }
    else
    {
        print_warning(0, "Something weird happened, no final SimulatedSystem initiated!\n",
                      "    Exiting...\n");
        return 5;
    }

    return 0;
}

// print help
void print_help()
{
    std::cout << "mergeconfig utility from simul++ package by JJ\n";
    std::cout << "==============================================\n";
    std::cout << "    Merge two .config files to one\n";
    std::cout << "    Add impact momentum to the two merged systems\n\n";
    std::cout << "USAGE:\n";
    std::cout << "    mergeconfig [OPTIONS(-f not possible)] CONFIG1 CONFIG2 [OPTIONS(-f possible)]\n";
    std::cout << "          when CONFIG1 and CONFIG2 are given without extension,\n";
    std::cout << "          the existence of the .field files with the same base name is expected\n";
    std::cout << "OPTIONS:\n";
    std::cout << "  Options without value:\n";
    std::cout << "    -h    show this help\n";
    std::cout << "    -H    'Half'   - keep the configurations separated (default)\n";
    std::cout << "    -M    'Mix'    - mix the configurations as much as possible\n";
    std::cout << "    -C    'Center' - keep the configuration in the middle of the box, if enlarging it (default)\n";
    std::cout << "    -E    'Expand' - expand the configuration if enlarging the box\n";
    std::cout << "  Options in pairs 'option value':\n";
    std::cout << "    -L    change (enlarge) final box size in directions specified by -d\n";
    std::cout << "    -T    set temperature and force random velocities assigning\n";
    std::cout << "    -d    directions of merge (in their respective order), e.g. 'x', 'yz' or 'xzy' (default 'z', 'zy' or 'zyx')\n";
    std::cout << "    -D    dimensions of merge: 1, 2 or 3 (linear, in 'square' or 'cube'), (default 1)\n";
    std::cout << "          each system is once, twice or four times in the final config (must be specified in advance in .field files)\n";
    std::cout << "    -o    output file name (without extension) (default output)\n";
    std::cout << "    -p    impact momentum; assign velocities to the both configuration (mass momentum conserved)\n";
    std::cout << "          applies only with '-H', not '-M'\n";
    std::cout << "    -f    specify field files (one for each .config files)\n";
}