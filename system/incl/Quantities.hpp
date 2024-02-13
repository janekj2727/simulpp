#include <vector>
#include <string>

#ifndef QUANT_NAMESPACE
#define QUANT_NAMESPACE
// measured quantities code, name and description
/**
 * @defgroup QuantsMeasurements Quantities and Measurement
 * @brief Group of classes, namespaces and other utilities for measurement.
 * 
 * Increasing number of measured quantities lead to more systematic approach to measurements in `simul++`.
 * Firstly, a set of possible quantities was described in a dedicated namespace.
 * Then, some quantities needed separate classes to store necessary data.
 * In future, maybe all quantities should be objects (classes) of itself...
 */
/** @file Quantities.hpp
 *  @ingroup QuantsMeasurements
*/
/**
 * @addtogroup QuantsMeasurements
 * @{
 * @namespace quant
 * @brief Namespace defining codes for quantities in `.control` file, their name in `.cpa` file and their description in `.prt` file
 * @ingroup QuantsMeasurements
 * 
 * Useful for the command `measure` in `.control` file.
 * In future, each quantity should be a class of itself... * 
 */
namespace quant
{
    /**
     * @brief Code of quantity for `.control` directive `measure`.
     */
    static const std::vector<std::string> code = {"Etot", "Ekin", "Epot", "Tkin", "Ttr", "Tin",
                                                  "V", "rho", "Elj", "Ebond", "MaxErrC", "Eextra",
                                                  "MaxAngleC", "PVVC", "P", "virial", "virconstr", "virbond",
                                                  "virLJ", "Xi", "dXidt", "Lambda", "dLambdadt", "MaxShakeIter",
                                                  "virangle", "Eangle", "Eelst", "virelst", "Edih", "virdih",
                                                  "accels", "Eelin", "H", "meanvel", "varvel", "kurtvel", "nopairs", "msd", "TtrLF", "TinLF"};
    /**
     * @brief Name of quantity for `.cpa` file.
     */
    static const std::vector<std::string> name = {"Etot", "Ekin", "Epot", "Tkin[K]", "Ttr[K]", "Tin[K]",
                                                  "V[AA3]", "rho[kg/m3]", "Elj", "Ebond", "MaxErrC", "Eextra",
                                                  "MaxAngleC", "PVVC[Pa]", "P[Pa]", "virial[p.u.]", "virconstr[p.u.]", "virbond[p.u.]",
                                                  "virLJ[p.u.]", "Xi[1]", "dXidt[ps-1]", "Lambda[1]", "dLambdadt[ps-1]", "MaxShakeIter",
                                                  "virangle[p.u.]", "Eangle", "Eelst", "virelst[p.u.]", "Edih", "virdih[p.u.]",
                                                  "accels", "Eelin", "H", "<v>[AA/ps]", "Var(v)[AA2/ps2]", "Kurt(v)[1]", "NoAtomPairs", "MSD[AA2]", "TtrLF[K]", "TinLF[K]"};
    /**
     * @brief Indication of whether the quantity is "energy" (has units of energy). 
     */
    static const std::vector<bool> iseng = {true, true, true, false, false, false,
                                            false, false, true, true, false, true,
                                            false, false, false, false, false, false,
                                            false, false, false, false, false, false,
                                            false, true, true, false, true, false,
                                            false, true, true, false, false, false, false, false, false, false};
    /**
     * @brief Description of quantity for `.prt` file. 
     */
    static const std::vector<std::string> description = {
        "total energy including corrections and energy of extra degrees of freedom",
        "kinetic energy",
        "potential energy including corrections",
        "kinetic temperature",
        "translational temperature (center-of-mass translation)",
        "internal temperature (rotational, vibrational, ...)",
        "volume",
        "density",
        "energy of disperse (vdw) interactions (intramolecular included) without (cutoff) corrections",
        "energy of flexible bonds",
        "maximum relative error of constrained bonds",
        "energy of extra degrees of freedom",
        "maximum error of angle between velocity and constrained bond from perpendicular",
        "pressure from virtual volume change (molecular based)",
        "pressure from virials",
        "total virial of forces",
        "virial of constraining forces",
        "virial of flexible bonds",
        "virial of disperse (vdw) interactions (intramolecular included)",
        "extended degree of freedom connected with thermostat (variable $\\xi$)",
        "time derivative of variable $\\xi$",
        "extended degree of freedom connected with barostat (variable $\\lambda$)",
        "time derivative of variable $\\lambda$",
        "maximum number of SHAKE iterations",
        "potential energy of angles",
        "virial of angles",
        "electrostatic energy",
        "virial of electrostatic forces",
        "enegy of dihedrals",
        "virial of dihedrals",
        "accelerations from potential energy differentiation (to `.acc` file)",
        "intramolecular electrostatic energy",
        "enthalpy",
        "mean velocity (average over $v_x$, $v_y$ and $v_z$ over all atoms)",
        "variance of the velocity distribution",
        "kurtosis of the velocity distribution",
        "number of intermolecular atom pairs currently used for force calculation",
        "mean square displacement of atoms",
        "translational temperature according to leap-frog (VERLET3) energy",
        "internal temperature according to leap-frog (VERLET3) energy"};
    /**
     * @brief Get index of quantity from its @ref code .
     * 
     * @param word Code of the quantity in `.control` directive `measure`.
     * @return int Index of the quantity in vectors of @ref quant namespace.
     */
    int idx_from_code(std::string word);

}
/**
 * @}
 */

#endif