/*
*  Mixed Autoregressive Moving Average Model
*  Time series analysis
*  According to Time Series Analysis with Applications in R by J.D.Cryer and K.-S.Chan
*
*  Author: JJ
*  Date: Sep 2021
*
*/

#include "ARMA.hpp"
#include "Matrix.hpp"
#include "LinearSystem.hpp"
#include "math_utils.hpp"
#include <iostream>
#include <cmath>

#define MAX_PSI 1000
#define SIGMA 1.0

ARMA::ARMA()
{
}

ARMA::ARMA(int p, int q, Vector *phi, Vector *theta)
{
    int i, j;
    double psi_i;

    m_p = p;
    m_q = q;
    m_phi = phi;
    m_theta = theta;
    m_psi = new Vector(MAX_PSI);
    m_gamma = new Vector(MAX_PSI);

    m_psi->operator()(0) = 1;
    // m_psi->operator()(1) = -m_theta->operator()(1) + m_phi->operator()(1);
    for (j = 1; j <= m_q; j++)
    {
        m_psi->operator()(j) = -m_theta->operator()(j);
        for (i = m_p; i > 0; i--)
        {
            psi_i = ((j - i) >= 0) ? m_psi->operator()(j - i) : 0;
            m_psi->operator()(j) += m_phi->operator()(i) * psi_i;
        }
    }
    for (j = m_q + 1; j < MAX_PSI; j++)
    {
        m_psi->operator()(j) = 0.0;
        for (i = m_p; i > 0; i--)
        {
            psi_i = ((j - i) >= 0) ? m_psi->operator()(j - i) : 0;
            m_psi->operator()(j) += m_phi->operator()(i) * psi_i;
        }
    }

    // std::cout << "theta = " << *m_theta << "\n";
    // std::cout << "phi = " << *m_phi << "\n";
    // std::cout << "psi = " << *m_psi << "\n";

    // gamma evaluation (Appendix C p. 85)
    Matrix A(m_p + 1, m_p + 1);
    for (i = 0; i <= m_p; i++)
    {
        A(i, i) = -1.0;
        for (j = 1; j <= m_p; j++)
        {
            A(i, abs(j - i)) += m_phi->operator()(j);
        }
    }
    Vector b(m_p + 1);
    for (i = 0; i <= m_p; i++)
    {
        for (j = i; j <= m_q; j++)
        {
            b(i) += m_theta->operator()(j) * m_psi->operator()(j - i);
        }
        b(i) *= SIGMA * SIGMA;
    }
    LinearSystem linsysgamma(A, b);
    Vector gamma(m_p + 1);
    gamma = linsysgamma.SolveGJ();

    // std::cout << "gamma = " << gamma << "\n";
    for (j = 0; j < m_p + 1; j++)
    {
        m_gamma->operator()(j) = gamma(j);
    }
    for (j = m_p + 1; j < MAX_PSI; j++)
    {
        for (i = 1; i <= m_p; i++)
        {
            m_gamma->operator()(j) += m_phi->operator()(i) * m_gamma->operator()(j - i);
        }
    }
    // std::cout << "gamma = " << *m_gamma << "\n";
}

ARMA::~ARMA()
{
    delete m_psi;
    delete m_gamma;
}

ARMA::ARMA(const ARMA &C)
{
    m_p = C.m_p;
    m_q = C.m_q;
    m_phi = C.m_phi;
    m_theta = C.m_theta;
    m_psi = C.m_psi;
    m_gamma = C.m_gamma;
}

ARMA &ARMA::operator=(const ARMA &C)
{
    m_p = C.m_p;
    m_q = C.m_q;
    m_phi = C.m_phi;
    m_theta = C.m_theta;
    m_psi = C.m_psi;
    m_gamma = C.m_gamma;

    return *this;
}

// Calculate one realization of the ARMA series
// Length of the Vector returned = length
// preplength is the length of the 'cooking' before collecting results
Vector ARMA::GetSeries(int length, int preplength)
{
    Vector result(length);
    Vector tempor(preplength);
    Vector e(m_q + 1);
    int i, j;

    for (i = 0; i < m_p; i++)
    {
        tempor(i) = randgauss();
    }
    for (i = 0; i <= m_q; i++)
    {
        e(i) = randgauss();
    }
    for (i = m_p; i < preplength; i++)
    {
        for (j = 1; j <= m_p; j++)
        {
            tempor(i) += m_phi->operator()(j) * tempor(i - j);
        }
        for (j = 0; j <= m_q; j++)
        {
            tempor(i) -= m_theta->operator()(j) * e(j);
        }
        for (j = m_q; j > 0; j--)
        {
            e(j) = e(j - 1);
        }
        e(0) = randgauss();
    }
    for (i = 0; i < m_p; i++)
    {
        result(i) = tempor(preplength - m_p + i);
    }
    for (i = m_p; i < length; i++)
    {
        for (j = 1; j <= m_p; j++)
        {
            result(i) += m_phi->operator()(j) * result(i - j);
        }
        for (j = 0; j <= m_q; j++)
        {
            result(i) -= m_theta->operator()(j) * e(j);
        }
        for (j = m_q; j > 0; j--)
        {
            e(j) = e(j - 1);
        }
        e(0) = randgauss();
    }

    return result;
}

// Print ARMA info and psi, gamma and rho
int ARMA::PrintSummary(FILE *fout)
{
    int i;

    // general parameters (user-defined)
    fprintf(fout, "# ARMA model\n");
    fprintf(fout, "# p = %d, phi = [ ", m_p);
    for (i = 0; i <= m_p; i++)
    {
        fprintf(fout, "%f  ", m_phi->operator()(i));
    }
    fprintf(fout, "]\n");
    fprintf(fout, "# q = %d, theta = [ ", m_q);
    for (i = 0; i <= m_q; i++)
    {
        fprintf(fout, "%f  ", m_theta->operator()(i));
    }
    fprintf(fout, "]\n");
    // calculated parameters
    Vector sumRho(m_gamma->GetSize());
    sumRho = m_gamma->CumulativeSum() * (1/m_gamma->operator()(0));
    fprintf(fout, "# 1: psi        2: gamma       3: rho         4: sum(rho, 0..n)\n");
    for (i = 0; i < MAX_PSI; i++)
    {
        fprintf(fout, " %14.10f %14.10f %14.10f %14.10f\n", m_psi->operator()(i), m_gamma->operator()(i), m_gamma->operator()(i) / m_gamma->operator()(0), sumRho(i));
    }
    return 0;
}

// return autocorrelation coeffs
Vector ARMA::GetRho(int length) const
{
    Vector rho(length);
    int i;
    for (i = 0; i < length; i++)
    {
        rho(i) = m_gamma->operator()(i) / m_gamma->operator()(0);
    }
    return rho;
}

// return theoretical value of variance
double ARMA::GetTheoreticalVariance() const
{
    return m_gamma->operator()(0);
}