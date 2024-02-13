/*
* A class to store constraint bond info for trial simulations of contrained dynamics
* Author JJ, Date Jan 2021
*/

#include <iostream>
#include <cmath>

#include "Vector.hpp"
#include "Matrix.hpp"
// #include "NewtonMethod.hpp"
#include "Atom.hpp"
#include "Molecule.hpp"
#include "ConstraintBond.hpp"

// Default (and only) constructor
ConstraintBond::ConstraintBond(double bondLength, int indexI, int indexJ, Molecule *parentMol) //, bool trueB)
{
    length = bondLength;
    atomI = indexI;
    atomJ = indexJ;
    parent = parentMol;
    recMass = 1 / (parent->atoms[atomI].mass + parent->atoms[atomJ].mass);
    redMass = (parent->atoms[atomI].mass * parent->atoms[atomJ].mass) * recMass;
    //trueBond = trueB;
}

// Calculate current bond length
double ConstraintBond::GetCurrentLength() const
{
    Vector delta(3);
    int i;

    for (i = 0; i < 3; i++)
    {
        delta[i] = parent->atoms[atomI].GetPosition(i) - parent->atoms[atomJ].GetPosition(i);
    }

    return delta.CalculateNorm();
}

double ConstraintBond::CalculateForces(double &virial, int shakeType, double epsc, double rescaling, double omega)
{
    /*
    * shakeType = 0 for Verlet, type SHAKE:
    *     implemented according to Nezbeda, Kolafa, Kotrla - Úvod do počítčových simulací (p. 122)
    * shakeType = -1 for Verlet MTTK type – RATTLE1:
    *     implemented according to standard SHAKE (see above) and DL POLY manual (p. 75) and Andersens article (1983)
    * shakeType = -2 for Verlet MTTK type – RATTLE2:
    *     implemented according Andersen (1983) and to DL POLY manual (p. 75)
    * shakeType = -3 for Verlet, type SHAKE, with predicted scaling:
    *     as shakeType = 0, but using predicted box rescaling...
    * shakeType = k > 0 for Gear:
    *     new method documented in simul++/documentation/SHAKEforGear/SHAKEforGear.pdf
    *     first k values in R are summed during prediction (omitting trvp values (past velocities...) in R)
    * epsc: maximum allowed relative error of constraint
    * returns relative error of constr. bond length (or 0 if unknown shakeType...)
    */
    int k;
    double lambda;
    double scalProd = 0.0;
    // double length2 = 0.0; // for RATTLE2 real length, but not better
    Vector rt(3);
    Vector rtph(3);
    Vector hvtph(3);

    switch (shakeType)
    {
    case 0: // standard SHAKE for Verlet
        for (k = 0; k < 3; k++)
        {
            rt(k) = parent->atoms[atomJ].R->operator()(0, k) - parent->atoms[atomI].R->operator()(0, k);
            rtph(k) = parent->atoms[atomJ].R->operator()(3, k) - parent->atoms[atomI].R->operator()(3, k);
        }
        // lambda = (CalculateScalarProduct(rtph, rtph) - CalculateScalarProduct(rt, rt))/(2 * CalculateScalarProduct(rtph, rt));
        // lambda = (CalculateScalarProduct(rtph, rtph) - length * length)/(2 * CalculateScalarProduct(rtph, rt));
        scalProd = CalculateScalarProduct(rtph, rtph);
        relErr = fabs(sqrt(scalProd) - length) / length; // rel. bond length error to return
        if (relErr < epsc)                               // if error smaller than the allowed maximum value, do not continue
        {
            for (k = 0; k < 3; k++)
            {
                parent->atoms[atomJ].constraintForces[k] = 0.0;
                parent->atoms[atomI].constraintForces[k] = 0.0;
            }
            return relErr;
        }
        lambda = omega * (scalProd - length * length) / (2 * length * length);
        // recMass = 1 / (parent->atoms[atomI].mass + parent->atoms[atomJ].mass);
        for (k = 0; k < 3; k++)
        {
            parent->atoms[atomJ].constraintForces[k] = -lambda * recMass * parent->atoms[atomI].mass * rt(k);
            parent->atoms[atomI].constraintForces[k] = +lambda * recMass * parent->atoms[atomJ].mass * rt(k);
        }
        virial += lambda * recMass * parent->atoms[atomJ].mass * parent->atoms[atomI].mass * length * length;
        return relErr;
    case -1: // RATTLE 1 for VerletNVTNose (and similar) (similar as SHAKE but time-shifted)
        for (k = 0; k < 3; k++)
        {
            rt(k) = parent->atoms[atomJ].R->operator()(2, k) - parent->atoms[atomI].R->operator()(2, k);
            rtph(k) = parent->atoms[atomJ].R->operator()(0, k) - parent->atoms[atomI].R->operator()(0, k);
        }
        // lambda = (CalculateScalarProduct(rtph, rtph) - CalculateScalarProduct(rt, rt))/(2 * CalculateScalarProduct(rtph, rt));
        // lambda = (CalculateScalarProduct(rtph, rtph) - length * length)/(2 * CalculateScalarProduct(rtph, rt));
        scalProd = CalculateScalarProduct(rtph, rtph);
        relErr = fabs(sqrt(scalProd) - length) / length; // rel. bond length error to return
        if (relErr < epsc)                               // if error smaller than the allowed maximum value, do not continue
        {
            for (k = 0; k < 3; k++)
            {
                parent->atoms[atomJ].constraintForces[k] = 0.0;
                parent->atoms[atomI].constraintForces[k] = 0.0;
            }
            return relErr;
        }
        lambda = omega * (scalProd - length * length) / (2 * length * length);
        // recMass = 1 / (parent->atoms[atomI].mass + parent->atoms[atomJ].mass);
        for (k = 0; k < 3; k++)
        {
            parent->atoms[atomJ].constraintForces[k] = -lambda * recMass * parent->atoms[atomI].mass * rt(k);
            parent->atoms[atomI].constraintForces[k] = +lambda * recMass * parent->atoms[atomJ].mass * rt(k);
        }
        virial += lambda * recMass * parent->atoms[atomJ].mass * parent->atoms[atomI].mass * length * length;
        return relErr;
    case -2: // RATTLE 2 for VerletNVTNose (and similar) (according to Andersen 1983)
        // hv_i(t+h)^C = hv_i(t+h) − h*k*r_ij(t+h)/m_i
        // hv_j(t+h)^C = hv_j(t+h) + h*k*r_ij(t+h)/m_j
        // k = r_ij . (v_i(t+h) - v_j(t+h))/d_ij^2 * (1/m_i + 1/m_j)
        // instead of k (in Andersen, here we use lambda...)
        for (k = 0; k < 3; k++)
        {
            rtph(k) = parent->atoms[atomJ].R->operator()(0, k) - parent->atoms[atomI].R->operator()(0, k);
            hvtph(k) = parent->atoms[atomJ].R->operator()(1, k) - parent->atoms[atomI].R->operator()(1, k);
        }
        scalProd = CalculateScalarProduct(rtph, hvtph); // h * r_ij . v_ij
        // length2 = CalculateScalarProduct(rtph, rtph); // not work better
        // what should be the value of relErr??? (26/01/2023 consulted with JK and average velocity added)
        relErr = fabs(scalProd) / length / hvtph.CalculateNorm(); 
        if (relErr < epsc)                // if error smaller than the allowed maximum value, do not continue
        {
            for (k = 0; k < 3; k++)
            {
                parent->atoms[atomJ].constraintForces[k] = 0.0;
                parent->atoms[atomI].constraintForces[k] = 0.0;
            }
            return relErr;
        }
        lambda = omega * scalProd / (length * length); // having tried real length – not work better
        // recMass = 1 / (parent->atoms[atomI].mass + parent->atoms[atomJ].mass);
        for (k = 0; k < 3; k++)
        {
            parent->atoms[atomJ].constraintForces[k] = -lambda * recMass * parent->atoms[atomI].mass * rtph(k);
            parent->atoms[atomI].constraintForces[k] = +lambda * recMass * parent->atoms[atomJ].mass * rtph(k);
        }
        virial += lambda * recMass * parent->atoms[atomJ].mass * parent->atoms[atomI].mass * length * length;
        return relErr;
    case -3: // standard SHAKE for Verlet with box rescaling
        for (k = 0; k < 3; k++)
        {
            rt(k) = parent->atoms[atomJ].R->operator()(0, k) - parent->atoms[atomI].R->operator()(0, k);
            rtph(k) = parent->atoms[atomJ].R->operator()(3, k) - parent->atoms[atomI].R->operator()(3, k);
        }
        // lambda = (CalculateScalarProduct(rtph, rtph) - CalculateScalarProduct(rt, rt))/(2 * CalculateScalarProduct(rtph, rt));
        // lambda = (CalculateScalarProduct(rtph, rtph) - length * length)/(2 * CalculateScalarProduct(rtph, rt));
        scalProd = CalculateScalarProduct(rtph, rtph);
        relErr = fabs(sqrt(scalProd) - length / rescaling) / length * rescaling; // rel. bond length error to return
        if (relErr < epsc)                                                       // if error smaller than the allowed maximum value, do not continue
        {
            for (k = 0; k < 3; k++)
            {
                parent->atoms[atomJ].constraintForces[k] = 0.0;
                parent->atoms[atomI].constraintForces[k] = 0.0;
            }
            return relErr;
        }
        lambda = omega * (scalProd - length * length / rescaling / rescaling) / (2 * length * length / rescaling);
        // recMass = 1 / (parent->atoms[atomI].mass + parent->atoms[atomJ].mass);
        for (k = 0; k < 3; k++)
        {
            parent->atoms[atomJ].constraintForces[k] = -lambda * recMass * parent->atoms[atomI].mass * rt(k) * rescaling;
            parent->atoms[atomI].constraintForces[k] = +lambda * recMass * parent->atoms[atomJ].mass * rt(k) * rescaling;
        }
        virial += lambda * recMass * parent->atoms[atomJ].mass * parent->atoms[atomI].mass * length * length;
        return relErr;
    default: // for Gear
        // (rtph - lambda*rt)^2 - d^2 = 0
        for (k = 0; k < 3; k++)
        {
            // interatomic distances (vectors)
            rt(k) = parent->atoms[atomJ].R->operator()(0, k) - parent->atoms[atomI].R->operator()(0, k); // calculated from predicted values in t
            rtph(k) = (parent->atoms[atomJ].rPI->operator()(k) + parent->atoms[atomJ].errorG->operator()(k)) * rescaling -
                   (parent->atoms[atomI].rPI->operator()(k) + parent->atoms[atomI].errorG->operator()(k)) * rescaling; // calculated from values predicted from t to t+h (with all known forces)
        }
        scalProd = CalculateScalarProduct(rtph, rtph);
        relErr = fabs(sqrt(scalProd) - length) / length; // rel. bond length error to return
        if (relErr < epsc)                               // if error smaller than the allowed maximum value, do not continue
        {
            for (k = 0; k < 3; k++)
            {
                parent->atoms[atomJ].constraintForces[k] = 0.0;
                parent->atoms[atomI].constraintForces[k] = 0.0;
            }
            return relErr;
        }
        lambda = omega * (scalProd - length * length) / (2 * length * length); // ala SHAKE version
        // ConstraintBond::recMass = 1 / (parent->atoms[atomI].mass + parent->atoms[atomJ].mass);
        for (k = 0; k < 3; k++)
        {
            parent->atoms[atomJ].constraintForces[k] = -lambda * recMass * parent->atoms[atomI].mass * rt(k);
            parent->atoms[atomI].constraintForces[k] = +lambda * recMass * parent->atoms[atomJ].mass * rt(k);
        }
        virial += lambda * recMass * parent->atoms[atomJ].mass * parent->atoms[atomI].mass * length * length;
        return relErr;
    }
}

// calculate angle between the constraint bond and velocity
double ConstraintBond::GetVelocityAngle() const
{
    int k;
    Vector r(3), hv(3);

    for (k = 0; k < 3; k++)
    {
        r(k) = parent->atoms[atomJ].R->operator()(0, k) - parent->atoms[atomI].R->operator()(0, k);
        hv(k) = parent->atoms[atomJ].R->operator()(1, k) - parent->atoms[atomI].R->operator()(1, k);
    }
    return CalculateAngle(r, hv);
}
