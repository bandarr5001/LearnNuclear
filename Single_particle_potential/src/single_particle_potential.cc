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
#include "./single_particle_potential.h"


double fermi_momentum(Parameters params) {
    return std::pow(6*std::pow(M_PI,2.0)
    *params.density_/params.degeneracy_g_,1.0/3.0);
}

double fermi_energy(Parameters params) {
    double p_momentum_fermi = fermi_momentum(params);
    return std::sqrt(params.mass_*params.mass_ + 
        p_momentum_fermi*p_momentum_fermi);
}

//Calculates the energy density of a spherical, ideal, noninteracting fermi gas
double energy_density_Fermigas(Parameters params){

    //calculates the fermi momentum using the density, n
    double p_momentum_fermi = fermi_momentum(params);

    //calculates the chemical potential using the fermi momentum
    double chem_potent_mu = std::sqrt(std::pow(params.mass_,2.0) 
    + std::pow(p_momentum_fermi,2.0));

    double fermi_energy_epsilon = fermi_energy(params);

    //Calculates fermi pressure at T = 0
    return (params.degeneracy_g_ / (16.0 * std::pow(M_PI,2.0)))
    * (2.0 * fermi_energy_epsilon*fermi_energy_epsilon*fermi_energy_epsilon
    * p_momentum_fermi) - (fermi_energy_epsilon*p_momentum_fermi*params.mass_)
    * (params.mass_*params.mass_*params.mass_*params.mass_
    * std::log((fermi_energy_epsilon + p_momentum_fermi)/params.mass_));
}

//Contains equations 25, 30, and 33 from Interactions in nuclear matter
//by Agnieszka Sorensen
int conditions(const gsl_vector* x, void* p,
					     gsl_vector* f) {
    auto* par = static_cast<potential_parameters*>(p);
    double A = gsl_vector_get(x,0);
    double B = gsl_vector_get(x,1);
    double tau = gsl_vector_get(x,2);
    const double energy_density_FG = energy_density_Fermigas(params);
    const double energy_epsilon_Fermi = fermi_energy(params);
    const double p_momentum_Fermi = fermi_momentum(params); 


    const double condition1 = (energy_density_FG/nSat_in_MeV3) 
        + ((A/2) + (B/tau))*nSat_in_MeV3 - par->binding_energy_;

    gsl_vector_set(f, 0, condition1);


    const double condition2 = energy_epsilon_Fermi 
        - (energy_density_FG/nSat_in_MeV3) + (A/2) + B*((tau - 1)/tau);

    gsl_vector_set(f, 1, condition2);
//Need to solve equations 31 and 32 explicitly
    const double K = 9 * ;
    const double K_0 = p_momentum_Fermi*p_momentum_Fermi
        * ;

    const double condition3 = K - K_0;

    gsl_vector_set(f, 2, condition3);

    return GSL_SUCCESS;
}

potential_results get_parameters(double guess_A, double guess_B, double guess_tau, double tolerance, potential_parameters pparams) {

  const size_t dim = 3;

  gsl_multiroot_function F;

  F.f = &conditions;
  F.n = dim;
  F.params = &params;
    
  gsl_vector* f = gsl_vector_alloc(dim);
  gsl_vector* x = gsl_vector_alloc(dim);
  gsl_vector_set(x, 0, guess_A);
  gsl_vector_set(x, 1, guess_B);
  gsl_vector_set(x, 2, guess_tau);

  const gsl_multiroot_fsolver_type* T = gsl_multiroot_fsolver_hybrids;
  gsl_multiroot_fsolver* s = gsl_multiroot_fsolver_alloc(T, dim);
  gsl_multiroot_fsolver_set(s, &F, x);

  int test_status;

  size_t i = 0, iMAX = 100;

  do {
    i++;
    test_status = gsl_multiroot_fsolver_iterate(s);
    test_status = gsl_multiroot_test_residual(s->f, tolerance);
  } while (test_status == GSL_CONTINUE && i < iMAX);

  parameter_results.A = gsl_vector_get(s->x,0);
  parameter_results.B = gsl_vector_get(s->x,1);
  parameter_results.tau = gsl_vector_get(s->x,2);
    
  if (test_status == GSL_SUCCESS) {
    std::cout << "Single particle potential parameters for V(" << params.density_ <<") at T = 0: A = " 
	      << parameter_results.A << " B = " << parameter_results.B << " tau = " << parameter_results.tau << " MeV" << std::endl;
  }

  else {
    std::cerr << "Failed to converge:" << gsl_strerror(test_status) 
	      << std::endl;
  }
  
  return ;
}
