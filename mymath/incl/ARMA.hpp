/*
*  Mixed Autoregressive Moving Average Model
*  Time series analysis
*  According to Time Series Analysis with Applications in R by J.D.Cryer and K.-S.Chan
*
*  Author: JJ
*  Date: Sep 2021
*
*/

#ifndef ARMAHEADER
#define ARMAHEADER

#include "Vector.hpp"


class ARMA
{
private:
    int m_p;
    int m_q;
    Vector *m_phi;
    Vector *m_theta;
    Vector *m_psi;
    Vector *m_gamma;
    ARMA();
public:
    ARMA(int p, int q, Vector *phi, Vector *theta);
    ~ARMA();
    ARMA(const ARMA& C);
    ARMA& operator=(const ARMA& C);
    Vector GetSeries(int length, int preplength);
    int PrintSummary(FILE *fout);
    Vector GetRho(int length) const;
    double GetTheoreticalVariance() const;
};

#endif
