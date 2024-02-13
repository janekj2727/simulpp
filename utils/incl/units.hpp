/* Units for simul++
* Unit conversion constants for simul++
* From input to program: x_program = x_in / U_quantity_unit
* To output from program: x_out = x_program * U_quantity_unit
* 
* Version 1.0 (June 2021)
* Author JJ
*/

#ifndef UNITSHEADER
#define UNITSHEADER

#include <cmath>

#define U_AVOGADRO 6.02214076e23  // SI 2019
#define U_BOLTZMANN 1.380649e-23  // SI 2019 (kB = 1 in program units)
#define U_EPSILON0 8.8541878128e-12 // CODATA value
#define U_ELEMENTARY_CHARGE 1.602176634e-19 // SI 2019
#define U_LENGTH_M 1.0e-10        // program units: Angstroms
#define U_VOLUME_M3 1.0e-30       // program units: Angstroms^3
#define U_ENERGY_JMOL (U_BOLTZMANN*U_AVOGADRO)   // program units: Kelvins (kB*K)
#define U_ENERGY_KJMOL (U_BOLTZMANN*U_AVOGADRO/1000.0)   // program units: Kelvins (kB*K)
#define U_ENERGY_KCALMOL (U_BOLTZMANN*U_AVOGADRO/4184.0)   // program units: Kelvins (kB*K)
#define U_ENERGY_KELVIN 1.0   // program units: Kelvins (kB*K)
#define U_MASS_GMOL (U_BOLTZMANN*U_AVOGADRO*0.1)       // program units: kB*K/(AA/ps)^2
#define U_PRESSURE_PA (U_BOLTZMANN/U_VOLUME_M3)  // program units 1 p.u. = 13.8 MPa
#define U_PRESSURE_KPA (U_BOLTZMANN/U_VOLUME_M3/1.0e3)  // program units 1 p.u. = 13.8 MPa
#define U_PRESSURE_MPA (U_BOLTZMANN/U_VOLUME_M3/1.0e6)  // program units 1 p.u. = 13.8 MPa
#define U_DENSITY_KGM3 (0.001*U_MASS_GMOL/(U_AVOGADRO*U_VOLUME_M3))     // program units g/mol/AA3      
#define U_CHARGE_C (sqrt(4.0*M_PI*U_EPSILON0*U_BOLTZMANN*U_LENGTH_M*U_ENERGY_KELVIN)) // program units sqrt(4*pi*epsilon0*k*K*AA)  
#define U_CHARGE_E (U_CHARGE_C/U_ELEMENTARY_CHARGE) // program units are complicated (sqrt(4*pi*epsilon0*k*K*AA))


#endif