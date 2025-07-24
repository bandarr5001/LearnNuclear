#ifndef core_statmech_functions_h
#define core_statmech_functions_h

#include <cmath>
#include <cstdlib>
#include <vector>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_vector.h>
#include <iomanip>

struct Parameters {
     double degeneracy_g_;
     double mass_;
     double density_;
     double temperature_;
};
extern Parameters params;

double get_energy_epsilon(double mass, double p_momentum);

//Calculates the Fermi-Dirac distribution
double get_FermiDirac_dist(double temperature, double energy_epsilon,
     double chem_potent_mu);

//Calculates the Bose-Einstein distribution
double get_BoseEinstein_dist(double temperature, double energy_epsilon,
     double chem_potent_mu);

//Calculates the Boltzmann distribution
double get_Boltzmann_dist(double temperature, double energy_epsilon, double chem_potent_mu);

#endif