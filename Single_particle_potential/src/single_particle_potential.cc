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


double fermi_momentum(void *p) {
    auto* par = static_cast<potential_parameters*>(p);
    return std::pow(6.0*std::pow(M_PI,2.0)
    *par->density_/par->degeneracy_g_,1.0/3.0);
}

double fermi_energy(void *p) {
    auto* par = static_cast<potential_parameters*>(p);
    double p_momentum_fermi = fermi_momentum(p);
    return std::sqrt(par->mass_*par->mass_ + 
        p_momentum_fermi*p_momentum_fermi);
}

//Calculates the energy density of a spherical, ideal, noninteracting fermi gas
double energy_density_Fermigas(void *p){
    auto* par = static_cast<potential_parameters*>(p);

    //calculates the fermi momentum using the density, n
    double p_momentum_fermi = fermi_momentum(p);

    double fermi_energy_epsilon = fermi_energy(p);

    double term1 = 2.0*fermi_energy_epsilon*fermi_energy_epsilon*fermi_energy_epsilon*p_momentum_fermi;

    double term2 = fermi_energy_epsilon*p_momentum_fermi*par->mass_*par->mass_;

    double term3 = par->mass_*par->mass_*par->mass_*par->mass_*std::log((fermi_energy_epsilon + p_momentum_fermi)/par->mass_);

    //Calculates energy density at T = 0
    return (par->degeneracy_g_ / (16.0 * M_PI*M_PI))
    * (term1 - term2 - term3);
}

double fermi_momentum_satdense(void *p) {
    auto* par = static_cast<potential_parameters*>(p);
    return std::pow(6.0*std::pow(M_PI,2.0)
    *nSat_in_MeV3/par->degeneracy_g_,1.0/3.0);
}

double fermi_energy_satdense(void *p) {
    auto* par = static_cast<potential_parameters*>(p);
    double p_momentum_fermi = fermi_momentum_satdense(p);
    return std::sqrt(par->mass_*par->mass_ + 
        p_momentum_fermi*p_momentum_fermi);
}

double energy_density_Fermigas_satdense(void *p){
    auto* par = static_cast<potential_parameters*>(p);

    //calculates the fermi momentum using the density, n
    double p_momentum_fermi = fermi_momentum_satdense(p);

    double fermi_energy_epsilon = fermi_energy_satdense(p);

    double term1 = 2.0*fermi_energy_epsilon*fermi_energy_epsilon*fermi_energy_epsilon*p_momentum_fermi;

    double term2 = fermi_energy_epsilon*p_momentum_fermi*par->mass_*par->mass_;

    double term3 = par->mass_*par->mass_*par->mass_*par->mass_*std::log((fermi_energy_epsilon + p_momentum_fermi)/par->mass_);

    //Calculates energy density at T = 0
    return (par->degeneracy_g_ / (16.0 * M_PI*M_PI))
    * (term1 - term2 - term3);
}

//Contains equations 25, 30, and 33 from Interactions in nuclear matter
//by Agnieszka Sorensen
int conditions(const gsl_vector* x, void* p,
					     gsl_vector* f) {
    auto* par = static_cast<potential_parameters*>(p);
    double A = gsl_vector_get(x,0);
    double B = gsl_vector_get(x,1);
    double tau = gsl_vector_get(x,2);
    const double energy_density_FG = energy_density_Fermigas(p);
    const double energy_epsilon_Fermi = fermi_energy(p);
    const double p_momentum_Fermi = fermi_momentum(p); 
    const double energy_density_FG_satdense = energy_density_Fermigas_satdense(p);
    const double energy_epsilon_Fermi_satdense = fermi_energy_satdense(p);
    const double p_momentum_Fermi_satdense = fermi_momentum_satdense(p); 


    double condition1 = (energy_density_FG/nSat_in_MeV3) 
        + ((A/2.0) + (B/tau))*nSat_in_MeV3 - par->mass_ - par->binding_energy_;

    gsl_vector_set(f, 0, condition1);


    double condition2 = energy_epsilon_Fermi_satdense 
        - (energy_density_FG_satdense/nSat_in_MeV3) + (A/2) + B*((tau - 1)/tau);

    gsl_vector_set(f, 1, condition2);

    double term1 = (par->degeneracy_g_*energy_epsilon_Fermi_satdense
        *p_momentum_Fermi_satdense*p_momentum_Fermi_satdense*p_momentum_Fermi_satdense)/(9*nSat_in_MeV3
        *M_PI*M_PI);

    double term2 = A;//*(par->density_/nSat_in_MeV3);

    double term3 = B*(tau-1.0);//*std::pow((par->density_/nSat_in_MeV3),tau - 1.0);

    double K = 9.0 * (term1 + term2 + term3);

    double condition3 = K - par->incompress_at_satdense_;

    gsl_vector_set(f, 2, condition3);

    std::cout << std::setprecision(10);
    std::cout << "[DEBUG] Inputs: A = " << A << ", B = " << B << ", tau = " << tau << std::endl;
    std::cout << "[DEBUG] Nuclear Matter Properties: g = " << par->degeneracy_g_ << ", mass = " << par->mass_ << ", density = " << par->density_ << ", B_0 = " << par->binding_energy_ << ", K_0 = " << par->incompress_at_satdense_ << std::endl;
    std::cout << "[DEBUG] fermi_momentum = " << p_momentum_Fermi << std::endl;
    std::cout << "[DEBUG] fermi_energy = " << energy_epsilon_Fermi << std::endl;
    std::cout << "[DEBUG] energy_density_FG = " << energy_density_FG << std::endl;
    std::cout << "[DEBUG] condition1 = " << condition1 << std::endl;
    std::cout << "[DEBUG] condition2 = " << condition2 << std::endl;
    std::cout << "[DEBUG] condition3 = " << condition3 << std::endl;

    if (!std::isfinite(condition1) || !std::isfinite(condition2) || !std::isfinite(condition3)) {
        std::cerr << "NaN or Inf in conditions!" << std::endl;
    }

    return GSL_SUCCESS;
}

potential_results get_parameters(potential_results parameter_results, double tolerance, potential_parameters pparams) {

  const size_t dim = 3;

  gsl_multiroot_function F;

  F.f = &conditions;
  F.n = dim;
  F.params = &pparams;
    
  //gsl_vector* f = gsl_vector_alloc(dim);
  gsl_vector* x = gsl_vector_alloc(dim);
  gsl_vector_set(x, 0, parameter_results.A);
  gsl_vector_set(x, 1, parameter_results.B);
  gsl_vector_set(x, 2, parameter_results.tau);

  const gsl_multiroot_fsolver_type* T = gsl_multiroot_fsolver_hybrids;
  gsl_multiroot_fsolver* s = gsl_multiroot_fsolver_alloc(T, dim);
  gsl_multiroot_fsolver_set(s, &F, x);

  int test_status;

  size_t i = 0, iMAX = 100;

  do {
    i++;
    test_status = gsl_multiroot_fsolver_iterate(s);
    if (test_status != GSL_SUCCESS && test_status != GSL_CONTINUE) {
        std::cerr << "Iteration failed: " << gsl_strerror(test_status) << std::endl;
    break;
    }
    test_status = gsl_multiroot_test_residual(s->f, tolerance);
    
    std::cout << "i: " << i << std::endl;
  } while (test_status == GSL_CONTINUE && i < iMAX);

  parameter_results.A = gsl_vector_get(s->x,0);
  parameter_results.B = gsl_vector_get(s->x,1);
  parameter_results.tau = gsl_vector_get(s->x,2);
    
  if (test_status == GSL_SUCCESS) {
    std::cout << "Single particle potential parameters for V(" << pparams.density_ <<") at T = 0: A = " 
	      << parameter_results.A << " B = " << parameter_results.B << " tau = " << parameter_results.tau << " MeV" << std::endl;
  }

  else {
    std::cerr << "Failed to converge:" << gsl_strerror(test_status) 
	      << std::endl;
  }
  
  return parameter_results;
}

double get_single_particle_potential(potential_results parameter_results,
     potential_parameters pparams) {
    double potential = parameter_results.A*(pparams.density_/nSat_in_MeV3) + parameter_results.B*std::pow(pparams.density_/nSat_in_MeV3,parameter_results.tau - 1.0);
    std::cout << "Single particle potential parameters for V(" << pparams.density_ << ") = " << potential << " MeV" << std::endl;
    return potential; 
}
