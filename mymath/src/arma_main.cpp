/*
*  Time series analysis – example main
*
*  Author: JJ
*  Date:
*
*/

#include "ARMA.hpp"
#include "Vector.hpp"
#include "file_openclose.hpp"
#include <iostream>
#include <random>
#include <chrono>
#include <cstring>

#define DEFAULT_OUTPUT "vector.txt"

int main(int argc, char **argv)
{
    ARMA *aa;
    int p = 1, q = 0;
    Vector theta(q + 1), phi(p + 1);
    theta(0) = -1.0;
    // theta(1) = 0.3;
    // theta(2) = 0.1;
    phi(0) = -1.0;
    phi(1) = 0.98;
    // phi(2) = 0.2;
    char filename[100];
    char writemode[] = "w";
    char correltxtname[] = "correl.txt";
    

    if (argc < 2)
    {
        strcpy(filename, DEFAULT_OUTPUT);
    }
    else
    {
        strcpy(filename, argv[1]);
    }

    aa = new ARMA(p, q, &phi, &theta);

    Vector realisation(10000);
    Vector corr(10000);
    realisation = aa->GetSeries(10000, 1000);
    corr = realisation.CorrelationCoefficients();
    aa->PrintSummary(stdout);

    FILE *fout;
    my_fopen_r(&fout, filename, writemode);
    realisation.PrintToFile(fout);
    my_fclose(&fout, filename);

    my_fopen_r(&fout, correltxtname, writemode);
    corr.PrintToFile(fout);
    my_fclose(&fout, correltxtname);
    
    delete aa;

    return 0;
}

// Cryer, Chan pp. 77

