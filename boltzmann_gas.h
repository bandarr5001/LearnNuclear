#ifndef boltzmann_gas_h
#define boltzmann_gas_h

#include <stdexcept> // for using std::runtime_error
#include <iostream>
#include <functional>
#include <cmath>
#include <gsl/gsl_vector.h>
#include "./core_statmech_functions.h"

//HIGH TEMPERATURE LIMIT ******************************************************

double get_density_integrand_boltzmann(double degeneracy_g, double t, 
    double temperature, double energy_epsilon, double chem_potent_mu);

double get_density_integral_boltzmann(double chem_potent_mu, double mass,
    double degeneracy_g, double temperature);

#endif