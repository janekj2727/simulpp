/*
 * A class to store simulated system for trial simulations of contrained dynamics
 * Author JJ, Date Jan 2021
 */

#include <cstring>  //
#include <iostream> //
#include <map>      //
#include <cmath>    //
#include <vector>   //
#include <cassert>  //

#include "Vector.hpp"
#include "Matrix.hpp"
#include "math_utils.hpp"
#include "LJsystem.hpp"
#include "Atom.hpp"            //
#include "Molecule.hpp"        //
#include "ConstraintBond.hpp"  //
#include "SimulatedSystem.hpp" //
#include "file_openclose.hpp"  //
#include "units.hpp"           //
#include "Quantities.hpp"      //
#include "general_utils.hpp"

#define ATOMR(i, j) ((*atomR)(i, j))
#define REAL_RAND() (((double)rand() + 1.0) / (RAND_MAX + 2.0))

// calculate potential energy (without forces)
double SimulatedSystem::CalculateEpot() const
{
    double Epot = 0.0;
    int i;

    for (i = 0; i < noMolecules; i++)
    {
        Epot += molecules[i].CalculateIntraMolEpot();
    }
    // now in Epot only EpotIntra...
    Epot += CalculateInterMolEpot();

    return Epot;
}

// calculate kinetic energy
double SimulatedSystem::CalculateEkin(double h, int velocityIndex) const
{
    double Ekininst = 0.0;
    int i;
    const Molecule *p_mol;

    for (i = 0; i < noMolecules; i++)
    {
        p_mol = &(molecules[i]);
        Ekininst += p_mol->CalculateEkin(h, velocityIndex);
    }
    return Ekininst;
}

// calculate system volume
double SimulatedSystem::CalculateVolume() const
{
    return boxSize[0] * boxSize[1] * boxSize[2];
}

// return maximum error of constrained bonds
double SimulatedSystem::CalculateMaxConstrErr() const
{
    double maxCErr = 0.0;
    int i;

    for (i = 0; i < noMolecules; i++)
    {
        maxCErr = fmax(maxCErr, molecules[i].MaxConstrErrSquared());
    }

    return sqrt(maxCErr);
}

// return maximum error of angle between constrained bond and velocities of connected atoms from 90° (pi/2)
double SimulatedSystem::CalculateMaxConstrVAngle() const
{
    double maxCErr = 0.0;
    int i, j;

    for (i = 0; i < noMolecules; i++)
    {
        for (j = 0; j < molecules[i].noConstrBonds; j++)
        {
            maxCErr = fmax(maxCErr, fabs(molecules[i].constrBonds[j].GetVelocityAngle() - M_PI / 2.0));
        }
    }

    return maxCErr;
}

// get energy cutoff correction
double SimulatedSystem::EnergyCutoffCorr() const
{
    return Ecorr / CalculateVolume();
}

// calculate configurational part of pressure using virtual volume change
double SimulatedSystem::CalculatePconfVVC()
{
    double P = 0.0;            // resulting pressure
    double E1 = 0.0, E2 = 0.0; // energy of smaller and bigger system
    double V1 = 1.0, V2 = 1.0; // volumes of smaller and bigger system
    std::vector<int> posAt;
    posAt.push_back(0); // positions in R to be scaled (only positions stored at index 0 should be sufficient)
    // static const double dL = 0.00333;                                                  // delta L (box size) (in % of current size)
    static const double dL = 0.001;
    static const double scaleUp = 1.0 + dL;                                            // scaling parameter
    static const double scaleDown = pow(2.0 - scaleUp * scaleUp * scaleUp, 1.0 / 3.0); // to maintain the 'center' at the current value

    // configurational pressure is minus the derivative of potential energy with volume
    // the derivative is calculated numerically using RescaleBox()

    RescaleBox(scaleUp, posAt); // rescale box up
    E2 = CalculateEpot();       // energy of the bigger system
    V2 = CalculateVolume();     // volume of the bigger system

    RescaleBox(scaleDown / scaleUp, posAt); // rescale box down
    E1 = CalculateEpot();                   // energy of the smaller system
    V1 = CalculateVolume();                 // volume of the smaller system

    RescaleBox(1 / scaleDown, posAt);

    P = ((E1 + Ecorr / V1) - (E2 + Ecorr / V2)) / (V2 - V1); // minus numerical derivative, central difference

    return P; // return pressure in program units, conversion must be done elsewhere (when printing)
}

// calculate configurational part of pressure using force virial
double SimulatedSystem::CalculatePconf() const
{
    double Pconf = 0.0; // resulting pressure
    double V = CalculateVolume();

    Pconf = -GetTotalVirial() / (3.0 * V) + (Ecorr / V / V);

    return Pconf; // return pressure in program units, conversion must be done elsewhere (when printing)
}

// calculate kinetic (ideal) part of pressure
double SimulatedSystem::CalculatePkin(double h, double Ekininst) const
{
    std::vector<Molecule>::const_iterator itmol;
    double Ekin = Ekininst;
    if (Ekininst == 0.0)
    {
        for (itmol = molecules.begin(); itmol != molecules.end(); itmol++)
        {
            Ekin += itmol->CalculateEkin(h, 1);
        }
    }
    return 2.0 * noDOFforPressure / noDegrOfFred * Ekin / 3.0 / CalculateVolume();
}

// sum virial terms and return total virial
double SimulatedSystem::GetTotalVirial() const
{
    double virtotal = 0.0;
    int i;
    for (i = 0; i < INTRA_MOL_TYPES; i++)
    {
        virtotal += virintramol[i];
    }
    for (i = 0; i < INTER_MOL_TYPES; i++)
    {
        virtotal += virintermol[i];
    }
    return virtotal + virconstr;
}

// get energy unit u_eng
double SimulatedSystem::GetEnergyUnit() const
{
    return u_eng;
}

// get force on atom by direct differentiation of Epot
double SimulatedSystem::ForceFromEpot(Atom *atom, int direction, double displacement) const
{
    double force = 0.0;
    double Ep, Em;

    // displace atom in positive direction
    atom->R->operator()(0, direction) += displacement;
    // measure system energy
    Ep = CalculateEpot();

    // displace atom in negative direction
    atom->R->operator()(0, direction) -= 2.0 * displacement;
    // measure system energy
    Em = CalculateEpot();

    // move atom back to its original position
    atom->R->operator()(0, direction) += displacement;

    // calculate force by central difference
    force = -(Ep - Em) / (2.0 * displacement);

    return force;
}

// calculate acceleration (a) of atoms by differences of Epot and print to .for file
int SimulatedSystem::PrintAccelsFromEpot(double time) const
{
    int i, j, k;
    double force;
    double mass;

    if (facc == nullptr)
    {
        return 71;
    }
    // print time
    fprintf(facc, "# time = %16.9E\n", time);
    // for each atom calculate force form Epot and print h^2/2 * acceleration (f/m)
    for (i = 0; i < noMolecules; i++)
    {
        for (j = 0; j < molecules[i].noAtoms; j++)
        {
            mass = molecules[i].atoms[j].mass;
            for (k = 0; k < 3; k++)
            {
                force = ForceFromEpot(&(molecules[i].atoms[j]), k);
                fprintf(facc, "%16.9E  ", force / mass);
            }
            fprintf(facc, "\n");
        }
    }
    // return 0 on success
    return 0;
}

// enthalpy calculation
double SimulatedSystem::CalculateEnthalpy(double Ekin) const
{
    double V = CalculateVolume();
    double enth = instEpot + Ekin + Ecorr / V; // internal energy
    double Penth = 0.0;
    if (Pext == -999999.9)
    {
        Penth = CalculatePkin(1.0, Ekin) + CalculatePconf();
    }
    else
    {
        Penth = Pext;
    }
    enth += Penth * V;
    return enth;
}

// mean velocity
double SimulatedSystem::MeanVelocity(double h, int dim) const
{
    assert(dim > -1);
    Vector vel(noAtomsTotal);
    Vector vel3(3 * noAtomsTotal);
    int i, j, k, l;
    l = 0;
    for (i = 0; i < noMolecules; i++)
    {
        for (j = 0; j < molecules[i].noAtoms; j++)
        {
            if (dim < 2)
            {
                vel(l++) = molecules[i].atoms[j].R->Read(1, dim) / h;
            }
            else
            {
                for (k = 0; k < 3; k++)
                {
                    vel3(l++) = molecules[i].atoms[j].R->Read(1, k) / h;
                }
            }
        }
    }

    if (dim < 2)
    {
        return vel.Mean();
    }
    else
    {
        return vel3.Mean();
    }

    return 9999999.99999; // error (should not reach this...)
}

// variance of velocity distribution
double SimulatedSystem::VarianceOfVelocity(double h, int dim) const
{
    assert(dim > -1);
    Vector vel(noAtomsTotal);
    Vector vel3(3 * noAtomsTotal);
    int i, j, k, l;
    l = 0;
    for (i = 0; i < noMolecules; i++)
    {
        for (j = 0; j < molecules[i].noAtoms; j++)
        {
            if (dim < 2)
            {
                vel(l++) = molecules[i].atoms[j].R->Read(1, dim) / h;
            }
            else
            {
                for (k = 0; k < 3; k++)
                {
                    vel3(l++) = molecules[i].atoms[j].R->Read(1, k) / h;
                }
            }
        }
    }

    if (dim < 2)
    {
        return vel.Variance();
    }
    else
    {
        return vel3.Variance();
    }

    return 9999999.99999; // error (should not reach this...)
}

// (excess) kurtosis of velocity distribution
double SimulatedSystem::KurtosisOfVelocity(double h, int dim) const
{
    assert(dim > -1);
    Vector vel(noAtomsTotal);
    Vector vel3(3 * noAtomsTotal);
    int i, j, k, l;
    l = 0;
    for (i = 0; i < noMolecules; i++)
    {
        for (j = 0; j < molecules[i].noAtoms; j++)
        {
            if (dim < 2)
            {
                vel(l++) = molecules[i].atoms[j].R->Read(1, dim) / h;
            }
            else
            {
                for (k = 0; k < 3; k++)
                {
                    vel3(l++) = molecules[i].atoms[j].R->Read(1, k) / h;
                }
            }
        }
    }

    if (dim < 2)
    {
        return vel.ExcessKurtosis();
    }
    else
    {
        return vel3.ExcessKurtosis();
    }

    return 9999999.99999; // error (should not reach this...)
}

double SimulatedSystem::LeapFrogTtr(double h, double Ekintot) const
{
    double Etr2 = 0.0; // 2 * translational energy
    std::vector<Molecule>::const_iterator itmol;
    int i, k;
    Vector ptotmh(3), ptotph(3); // translational momentum of molecule (i−h/2, i+h/2)
    Matrix *atomR;
    double mass;
    for (itmol = molecules.begin(); itmol != molecules.end(); itmol++)
    {
        ptotmh.Clear();
        ptotph.Clear();
        for (i = 0; i < itmol->noAtoms; i++)
        {
            atomR = itmol->atoms[i].R;
            mass = itmol->atoms[i].mass;
            for (k = 0; k < 3; k++)
            {
                if (Pext == -999999.9) // not barostat
                {
                    ptotmh(k) += mass * (atomR->operator()(0, k) - atomR->operator()(2, k));
                }
                else // with barostat (different values in ATOMR(2, :))
                {
                    ptotmh(k) += mass * atomR->operator()(2, k);
                }
                ptotph(k) += mass * (atomR->operator()(3, k) - atomR->operator()(0, k));
            }
        }
        Etr2 += (CalculateScalarProduct(ptotmh, ptotmh) + CalculateScalarProduct(ptotph, ptotph)) / itmol->molMass;
    }
    Etr2 /= 2.0 * h * h;

    // translation 3*noMolecules - systemConserved (systemConserved = 3*N - noDegrOfFred - noTotalConstr)
    return Etr2 / (3 * noMolecules - 3 * noAtomsTotal + noDegrOfFred + noConstrTotal);
}

double SimulatedSystem::LeapFrogTin(double h, double Ekintot) const
{
    double Etr2 = 0.0; // 2 * translational energy
    std::vector<Molecule>::const_iterator itmol;
    int i, k;
    Vector ptotmh(3), ptotph(3); // translational momentum of molecule (i−h/2, i+h/2)
    Matrix *atomR;
    double mass;
    for (itmol = molecules.begin(); itmol != molecules.end(); itmol++)
    {
        ptotmh.Clear();
        ptotph.Clear();
        for (i = 0; i < itmol->noAtoms; i++)
        {
            atomR = itmol->atoms[i].R;
            mass = itmol->atoms[i].mass;
            for (k = 0; k < 3; k++)
            {
                if (Pext == -999999.9) // not barostat
                {
                    ptotmh(k) += mass * (atomR->operator()(0, k) - atomR->operator()(2, k));
                }
                else // with barostat (different values in ATOMR(2, :))
                {
                    ptotmh(k) += mass * atomR->operator()(2, k);
                }
                ptotph(k) += mass * (atomR->operator()(3, k) - atomR->operator()(0, k));
            }
        }
        Etr2 += (CalculateScalarProduct(ptotmh, ptotmh) + CalculateScalarProduct(ptotph, ptotph)) / itmol->molMass;
    }
    Etr2 /= 2.0 * h * h;
    // internal temperature
    return (2 * Ekintot - Etr2) / (3 * noAtomsTotal - 3 * noMolecules - noConstrTotal);
}

// initialization of .cpa file to write statistics
int SimulatedSystem::InitCpaFile(char *cpa_name, std::set<int> measuredQ, bool append, double &time)
{
    char write[] = "w";
    char read[] = "r";
    char line[MAX_COMMENT * 10]; // should be sufficient for line length in .cpa file
    char buff[MAX_COMMENT * 10];
    char *last_no;
    char engunit[7];
    char *acc_name;
    char *acc_ext;
    double last_time = 0.0; // the time of the last measurement

    if (append)
    {
        write[0] = 'a';
        // get current time from cpa file
        if (my_fopen_r(&fcpa, cpa_name, read) != 0)
        {
            print_warning(0, "Cannot open %s file to append measurement\n" + std::string(cpa_name) + "\n", "    File does not exist?\n");
            return 39;
        }
        else
        {
            while (fgets(line, MAX_COMMENT * 10, fcpa) != nullptr)
            {
                if (strlen(line) > 10)
                    strcpy(buff, line);
            }
            last_no = strrchr(buff, ' ');
            if (last_no == nullptr)
            {
                print_warning(1, "Cannot get last measurement time from " + std::string(cpa_name) + " file\n", "    Setting current time to 0.0\n");
            }
            if (sscanf(last_no, "%lf", &last_time) < 1)
            {
                print_warning(1, "Cannot get last measurement time from " + std::string(cpa_name) + " file\n", "    Setting current time to 0.0\n");
                last_time = 0.0;
            }
            else
            {
                if (last_time < 0.0)
                {
                    print_warning(1, "Cannot get last measurement time from " + std::string(cpa_name) + " file, maybe corrupted\n", "    Setting current time to 0.0\n");
                    last_time = 0.0;
                }
            }
            my_fclose(&fcpa, cpa_name);
        }
    }
    time = last_time;

    // copy measured quantities set from simulation (provided as argument) to system (this)
    measuredQuant.assign(measuredQ.begin(), measuredQ.end());

    // open .cpa file
    if (my_fopen_w(&fcpa, cpa_name, write) != 0)
    {
        if (write[0] == 'a')
        {
            print_warning(0, "Cannot open " + std::string(cpa_name) + " file to append measurement\n");
        }
        else
        {
            print_warning(0, "Cannot open " + std::string(cpa_name) + " file to write\n");
        }
        return 39;
    }

    if (fabs(u_eng - U_ENERGY_KCALMOL) < 1e-12)
    {
        strcpy(engunit, "[kcal]");
    }
    else if (fabs(u_eng - U_ENERGY_KJMOL) < 1e-12)
    {
        strcpy(engunit, "[kJ]");
    }
    else
    {
        strcpy(engunit, "[K]");
    }

    if (!append)
    {
        // write header to .cpa file
        fprintf(fcpa, "# block = 1; Statistics of simulation.\n#");
    }
    // for each quantity in measuredQuant add title (name of quantity)
    int i = 0;
    std::vector<int>::const_iterator it;
    std::string qname;
    for (it = measuredQuant.begin(); it != measuredQuant.end(); it++)
    {
        if (((*it) == 31) || ((*it) == 38))
        {
            continue; // accelerations and MSD later
        }
        qname = quant::name.at((*it) - 1);
        if (quant::iseng.at((*it) - 1))
            qname.append(engunit);

        if (!append)
        {
            fprintf(fcpa, "%2d:%-14s ", ++i, qname.c_str());
        }
    }
    // mean square displacements if requested
    if (measuredQ.count(38) > 0)
    {
        msd = new MeanSquareDispl(noAtomsTotal, &molecules);
        i = msd->PrintHeader(fcpa, ++i, "%2d:%-14s ");
    }

    if (!append)
    {
        fprintf(fcpa, "%2d:time[ps]\n", ++i);
    }

    // if force from Epot should be measured, init facc
    if (measuredQ.count(31) > 0)
    {
        // set .for file name
        acc_name = (char *)malloc((strlen(cpa_name) + 1) * sizeof(char));
        strcpy(acc_name, cpa_name);
        acc_ext = strstr(acc_name, ".cpa");
        strcpy(acc_ext, ".acc");

        // open .for file
        if (my_fopen_w(&facc, acc_name, write) != 0)
        {
            if (write[0] == 'a')
            {
                print_warning(0, "Cannot open " + std::string(acc_name) + " file to append accelerations\n");
            }
            else
            {
                print_warning(0, "Cannot open " + std::string(acc_name) + " file to write\n");
            }
            return 72;
        }
    }

    return 0;
}

// print measured quantities
int SimulatedSystem::PrintMeasurement(double time, double h)
{
    // when changing measured quantities, you must change also header in InitCpaFile... (see control.md)
    // int i;
    std::vector<Molecule>::iterator itmol;
    double Ekin = 0.0;
    double Tkin = 0.0;
    double Ekintr = 0.0;
    double Ttr = 0.0;
    double Tin = 0.0;
    double V = 0.0;
    double Eljintra = 0.0;
    double Eelstintra = 0.0;

    // calculate Ekin
    for (itmol = molecules.begin(); itmol != molecules.end(); itmol++)
    {
        Ekintr += itmol->CalculateEkinTrans(h);
        Ekin += itmol->CalculateEkin(h, 1);
        Eljintra += itmol->Epot[3];
        Eelstintra += itmol->Epot[4];
    }

    Tkin = 2 * Ekin / noDegrOfFred;
    // translation 3*noMolecules - systemConserved (systemConserved = 3*N - noDegrOfFred - noTotalConstr)
    Ttr = 2 * Ekintr / (3 * noMolecules - 3 * noAtomsTotal + noDegrOfFred + noConstrTotal);
    Tin = 2 * (Ekin - Ekintr) / (3 * noAtomsTotal - 3 * noMolecules - noConstrTotal);
    V = CalculateVolume();

    // for each measured quantity print the value
    std::vector<int>::const_iterator it;
    double quantity;
    for (it = measuredQuant.begin(); it != measuredQuant.end(); it++)
    {
        switch (*it)
        {
        case 1:
            if (boundaryCond != 0)
            {
                quantity = (instEpot + Ekin + Eextra + Ecorr / V) * u_eng;
            }
            else
            {
                quantity = (instEpot + Ekin + Eextra) * u_eng;
            }
            break;
        case 2:
            quantity = Ekin * u_eng;
            break;
        case 3:
            if (boundaryCond != 0)
            {
                quantity = (instEpot + Ecorr / V) * u_eng;
            }
            else
            {
                quantity = instEpot * u_eng;
            }
            break;
        case 4:
            quantity = Tkin;
            break;
        case 5:
            quantity = Ttr;
            break;
        case 6:
            quantity = Tin;
            break;
        case 7: // volume
            quantity = V;
            break;
        case 8: // density
            quantity = totalMass / V * U_DENSITY_KGM3;
            break;
        case 9: // energy of LJ interactions without correction
            quantity = (Epotintermol[1] + Eljintra) * u_eng;
            break;
        case 10: // energy of flexible bonds
            quantity = Epotintramol[0] * u_eng;
            break;
        case 11: // maximum error of constrained bond
            quantity = CalculateMaxConstrErr();
            break;
        case 12: // energy of the extra degrees of freedom (extended DOF)
            quantity = Eextra * u_eng;
            break;
        case 13: // maximum angle difference of constraint bond and velocity from perpendicular
            quantity = CalculateMaxConstrVAngle();
            break;
        case 14: // pressure from virtual volume change (molecular based) (number of degrees of freedom ??? )
            if (pressureNfCorrection)
            {
                quantity = (Ttr * noMolecules / V + CalculatePconfVVC()) * U_PRESSURE_PA; // best agreement with MACSIMUS for atoms
            }
            else
            {
                quantity = (2.0 * Ekintr / 3.0 / V + CalculatePconfVVC()) * U_PRESSURE_PA; // good agreement with MACSIMUS
            }
            break;
        case 15:
            // pressure from virial (number of degrees of freedom according to MACSIMUS)
            quantity = (CalculatePkin(h, Ekin) + CalculatePconf()) * U_PRESSURE_PA; // best agreement with DL POLY if noDOfaccPressure = noDegrOfFred (pressureNfCorrection = false)
            break;
        case 16: // virial itself (total)
            quantity = GetTotalVirial();
            break;
        case 17: // virial of constr.
            quantity = virconstr;
            break;
        case 18: // virial of bonds (flexible)
            quantity = virintramol[0];
            break;
        case 19: // virial of LJ interactions
            quantity = virintermol[1] + virintramol[3];
            break;
        case 20: // extended variable Xi connected with thermostat
            quantity = Xi;
            break;
        case 21: // velocity of Xi
            quantity = Xivel;
            break;
        case 22: // extended variable Lambda connected with barostat
            quantity = Lambda;
            break;
        case 23: // velocity of Lambda
            quantity = Lambdavel;
            break;
        case 24: // maximum number of SHAKE iterations
            quantity = maxShakeIter;
            maxShakeIter = 0; // zero before the next cycle
            break;
        case 25: // angle potential energy
            quantity = Epotintramol[1] * u_eng;
            break;
        case 26: // virial of angles (flexible)
            quantity = virintramol[1];
            break;
        case 27: // electrostatic energy
            quantity = (Epotintermol[0] + Eelstintra) * u_eng;
            break;
        case 28: // electrostatic virial
            quantity = virintermol[0] + virintramol[4];
            break;
        case 29: // dihedral potential energy
            quantity = Epotintramol[2] * u_eng;
            break;
        case 30: // virial of dihedralsPrintCpa
            quantity = virintramol[2];
            break;
        case 31: // forces (accelerations) from Epot to .acc file
            PrintAccelsFromEpot(time);
            break;
        case 32: // intramolecular electrostatic energy
            quantity = Eelstintra * u_eng;
            break;
        case 33: // enthalpy
            quantity = CalculateEnthalpy(Ekin) * u_eng;
            break;
        case 34: // mean velocity (should be 0)
            quantity = MeanVelocity(h);
            break;
        case 35: // variance of velocity distribution
            quantity = VarianceOfVelocity(h);
            break;
        case 36: // (excess) kurtosis of velocity distribution
            quantity = KurtosisOfVelocity(h);
            break;
        case 37: // number of intermolecular terms currently considered
            quantity = noPairs;
            break;
        case 39: // translational temperature according to leap-frog
            quantity = LeapFrogTtr(h, Ekin);
            break;
        case 40: // internal temperature according to leap-frog
            quantity = LeapFrogTin(h, Ekin);
            break;
        }
        if ((*it) != 38) // if not MSD print (if yes print nothing)
        {
            fprintf(fcpa, "%16.9E  ", quantity);
        }
    }

    // if MSD requested print MSD of each atom type
    if (msd != nullptr)
    {
        msd->PrintMSD(fcpa, "%16.9E  ");
    }

    fprintf(fcpa, "%16.9E\n", time);

    if (std::isnan(Ekin) || std::isnan(instEpot))
    {
        print_warning(0, "Energy is NaN = simulation crashed\n");
        return 45;
    }

    return 0;
}
