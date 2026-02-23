#ifndef single_particle_potential_h
#define single_particle_potential_h

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_vector.h>
#include <iomanip>
#include "./core_statmech_functions.h"
#include "./constants.h"

struct potential_parameters {
    const double degeneracy_g_;
    const double mass_;
    const double density_;
    const double binding_energy_;
    double incompress_at_satdense_;
};
extern potential_parameters pparams;

struct potential_results {
    double A;
    double B;
    double tau;
};
extern potential_results parameter_results;

double fermi_momentum(void *p);

double fermi_energy(void *p);

//Calculates the energy density of a spherical, ideal, noninteracting fermi gas
//at T = 0
double energy_density_Fermigas(void *p);

double fermi_momentum_satdense(void *p);

double fermi_energy_satdense(void *p);

//Calculates the energy density of a spherical, ideal, noninteracting fermi gas
//at T = 0
double energy_density_Fermigas_satdense(void *p);

//Contains equations 25, 30, and 33 from Interactions in nuclear matter
//by Agnieszka Sorensen
int conditions(const gsl_vector* x, void* p,
					     gsl_vector* fvec);

//Uses GSL multiroot to calculate the parameters which all for the calculation 
//of the single particle potential 
potential_results get_parameters(potential_results parameter_results, double tolerance, potential_parameters pparams);

double get_single_particle_potential(potential_results parameter_results,
     potential_parameters pparams);

#endif
