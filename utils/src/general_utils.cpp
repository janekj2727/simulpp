/*
 * General utils originally for Simul++
 * Author: JJ
 * Date: Apr 2021
 */

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <cstring>

#include "general_utils.hpp"
#include "file_openclose.hpp"
#include "Matrix.hpp"
#include "Vector.hpp"
#include "MDTable.hpp"

// print progres bar [|||||    ] 54 % with \r (work on console output only - not a bug, it is a feature... :-)
void print_progress_bar(int percentage)
{
    int i;
    int j = (int)floor(percentage / 2);

    std::cout << "  [" << COL_ESCAPE_IN << 7 << COL_ESCAPE_OUT;

    for (i = 0; i < j; i++)
    {
        std::cout << " ";
    }
    std::cout << COL_RESET;
    for (; i < 50; i++)
    {
        std::cout << " ";
    }
    std::cout << "]" << std::setw(3) << percentage << " %\r" << std::flush;
}

// print warning to prt or to std::cerr (defined on compile time)
void print_warning(int error_level, std::string str1, std::string str2, std::string str3, std::string str4)
{
    fg_cols fgcol;
    bg_cols bgcol;
    std::string intro;

    switch (error_level)
    {
    case 0:
        intro = " ERR: ";
        fgcol = fg_red;
        bgcol = bg_red;
        break;
    case 1:
        intro = " WRN: ";
        fgcol = fg_yellow;
        bgcol = bg_yellow;
        break;
    default:
        intro = " INF: ";
        fgcol = fg_white;
        bgcol = bg_blue;
        break;
    }

    std::cerr << COL_ESCAPE_IN << bgcol << COL_ESCAPE_OUT << TEXTBF << intro << COL_RESET << " ";
    std::cerr << COL_ESCAPE_IN << fgcol << COL_ESCAPE_OUT << str1 << str2 << str3 << str4 << COL_RESET;
}

int basic_statistics(std::ofstream &fprt, char *cpa_name, std::vector<std::string> qnames)
{
    FILE *fcpa;
    char read[2] = "r";
    size_t siz = 1000;
    char *line;
    char *p_line;
    int i, j, currline = 0;
    double aux, aux1, aux2, mean, var;
    std::vector<std::vector<double>> stdmat;
    Vector *auxvec;
    // number of columns based on qnames size
    int no_cols = (int)qnames.size();

    MDTable<std::string, double, double, double, double, double, double> tab({"Quantity [unit]", "Mean value", "Stderr estim.\\*", "Min. value", "Max. value", "Variance", "Error of var\\**"});

    if ((line = (char *)calloc(siz, sizeof(char))) == NULL)
    {
        print_warning(1, "Cannot allocate memory for final statistics!\n");
        return -2;
    }

    if (my_fopen_r(&fcpa, cpa_name, read) != 0)
    {
        print_warning(1, "File " + std::string(cpa_name) + " cannot be opened for reading!\n", "    No statistics will be printed to protocol.\n");
        return -1;
    }

    for (i = 0; i < no_cols; i++)
    {
        stdmat.push_back(std::vector<double>());
    }

    while (getline(&line, &siz, fcpa) > 0)
    {
        if ((line[0] == '#') || (line[0] == '!')) // comment line
        {
            continue;
        }
        currline++;
        p_line = &(line[0]);
        while (*(++p_line) == ' ') // to avoid whitespaces at the beginning of the line
            ;
        for (i = 0; i < no_cols; i++)
        {
            if (sscanf(p_line, "%lf", &aux) < 0)
            {
                print_warning(1, "Line around " + std::to_string(currline + 2) + " in " + std::string(cpa_name) + " is too short!");
                return -3;
            }
            p_line = strchr(p_line, ' ');
            while (*(++p_line) == ' ')
                ;
            stdmat[i].push_back(aux);
        }
    }

    if (my_fclose(&fcpa, cpa_name) != 0)
    {
        print_warning(-1, "File " + std::string(cpa_name) + " cannot be closed after final statistics!\n");
        return -3;
    }

    auxvec = new Vector(currline);

    for (i = 0; i < no_cols; i++)
    {
        for (j = 0; j < currline; j++)
        {
            auxvec->operator()(j) = stdmat[i][j];
        }
        mean = auxvec->MeanWithErrorEstimate(aux1);
        var = auxvec->VarianceWithError(aux2, std::min(10, (int)ceil(currline / 10.0)));
        tab.addRow(qnames[i], mean, aux1, auxvec->Min(), auxvec->Max(), var, aux2);
    }
    tab.setColumnPrecision(std::vector<int>({0, 9, 4, 9, 9, 9, 4}));

    fprt << "\n-----------------------------------------------------------------------------------------------\n";
    fprt << "\n# Simulation results\n\n";
    fprt << "Total number of " << std::to_string(currline) << " measurements was done. \n";
    fprt << "Short statistics follows.\n\n";
    tab.print(fprt);
    fprt << "\\*) Standard error of mean is estimated using the 0.6*(max−min) formula form cumulative average\n\n";
    fprt << "\\**) Error of variance is calculated from 10 blocks and is relevant only for large data sets\n\n";

    delete auxvec;
    free(line);
    line = nullptr;

    return currline;
}
