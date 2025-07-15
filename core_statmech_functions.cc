#include "./core_statmech_functions.h"
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>

double get_energy_epsilon(double mass, double p_momentum) {
    //particle energy at the lower bound and mapped to a range of (0,1)
    return std::sqrt(std::pow(mass,2.0) + std::pow(p_momentum,2.0));
}

//Fermi-Dirac Distribution
double get_FermiDirac_dist(double temperature, double energy_epsilon,
     double chem_potent_mu) {
    return 1.0 / (std::exp((energy_epsilon - chem_potent_mu)/temperature)
    + 1.0);
}

//Bose-Einstein Distribution
double get_BoseEinstein_dist(double temperature, double energy_epsilon, 
    double chem_potent_mu) {
    return 1.0
    / (std::exp((energy_epsilon - chem_potent_mu)/temperature) - 1.0);
}

//Boltzmann Distribution
double get_Boltzmann_dist(double temperature, double energy_epsilon,
     double chem_potent_mu) {
    return 1.0 / (std::exp((energy_epsilon - chem_potent_mu)/temperature));
}