#ifndef fermi_gas_h
#define fermi_gas_h

#include <stdexcept> // for using std::runtime_error
#include <iostream>
#include <functional>
#include <cmath>
#include <gsl/gsl_vector.h>
#include "./core_statmech_functions.h"

//T=0 *************************************************************************

double fermi_gas_press_temp0_pFermi(double degeneracy_g, double chem_potent_mu,
    double p_momentum_fermi, double mass);

//calculates the fermi gas pressure at T=0 as a function of chemical potential
double fermi_gas_press_temp0_chem(double chem_potent_mu, double degeneracy_g,
    double mass);

//Calculates the fermi gas pressure at T=0 as a function of the density
double fermi_gas_press_temp0_density(double chem_potent_mu,
    double degeneracy_g, double mass);

//FINITE TEMPERATURE **********************************************************

/*All of these formulas are with respect to t rather than using p = (1-t)/t
due to the 1/t^2 term. While there is an equivalent form for implementing p
in this way, it does not work if I actually input a momentum. The mapping takes
the form of:
integral(0 -> infinity){f(x)} = integral(0 -> 1){f((1-t)/t) * (1/t^2)}
There is no way that I can make the 1/t^2 term disappear if I were to use the
momentum p as an input due to the nature of the mapping
(Unless I'm misunderstanding what is being asked of me and it's just a notation
issue rather than an actual implementation issue) 
*/

//Calculates the particle energy epsilon
double get_energy_epsilon(double mass, double t);

//Fermi-Dirac Distribution ****************************************************
//Calculates the Fermi-Dirac distribution
double get_FermiDirac_dist(double temperature, double energy_epsilon,
    double chem_potent_mu);

//Calculates the integrand for the density n_B of a Fermi gas
double get_density_integrand_fermi(double t, double energy_epsilon, double chem_potent_mu, Parameters& params);

//Calculates the integral for the density n_B of a Fermi gas
double get_density_integral_fermi(Parameters& params, double chem_potent_mu);

double get_density_integral_fermi_gsl(double chem_potent, void *p);

double get_density_derivative_integrand_fermi(double t, double energy_epsilon, double chem_potent_mu, Parameters& params);

double get_density_derivative_integral_fermi_numerical(double chem_potent_mu, Parameters& params);

double get_density_derivative_integral_fermi(double chem_potent_mu, Parameters& params);

double get_density_derivative_integral_fermi_gsl_numerical(double chem_potent, void *p);

double get_density_derivative_integral_fermi_gsl(double chem_potent, void *p);

//Calculates the integrand for the pressure of a spherical Fermi gas
double get_P_spherical_integrand_fermi(double t, 
    double energy_epsilon, double chem_potent_mu, Parameters& params);

//calculates the total pressure for an ideal fermi gas
double fermi_gas_pressure_integral_spherical(Parameters& params, double chem_potent_mu);

#endif