#ifndef bose_gas_h
#define bose_gas_h

#include <stdexcept> // for using std::runtime_error
#include <iostream>
#include <functional>
#include <cmath>
#include <gsl/gsl_vector.h>
#include "./core_statmech_functions.h"

//Bose-Einsteins Distribution *************************************************

//Calculates the integrand for the density n_B of a Bose gas
double get_density_integrand_bose(double degeneracy_g, double t,
    double temperature, double energy_epsilon, double chem_potent_mu);

//Calculates the integral for the density n_B of a Bose gas
double get_density_integral_bose(double chem_potent_mu, double mass,
    double degeneracy_g, double temperature);

//Calculates the integrand of the pressure of
double get_P_spherical_integrand_bose(double degeneracy_g, double t, 
    double energy_epsilon, double temperature, double chem_potent_mu);

//calculates the total pressure for an ideal fermi gas using the 
double bose_gas_pressure_integral_spherical(double temperature, double chem_potent_mu, double degeneracy_g,
    double mass);

#endif