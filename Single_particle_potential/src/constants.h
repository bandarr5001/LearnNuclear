#ifndef constants_h
#define constants_h

#include <cmath>



//////////////////////////////////////////////////////////////////////////////////////////
// Conversion constants

// from PDG: https://pdg.lbl.gov/2019/reviews/rpp2018-rev-phys-constants.pdf
const double hBarC = 197.326980; //in MeV fm



//////////////////////////////////////////////////////////////////////////////////////////
// Particle properties

// particle masses, in MeV
const double mNucleon = 938.918; // mean proton + neutron mass, in MeV



//////////////////////////////////////////////////////////////////////////////////////////
// Nuclear matter properties

// saturation density
const double nSat_in_1overfm3 = 0.160; // in fm^-3
const double nSat_in_MeV3 = nSat_in_1overfm3 * hBarC * hBarC * hBarC; // in MeV^3

// binding energy
const double eBinding = -16.3; // in MeV



#endif
