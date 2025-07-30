#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_vector.h>
#include <iomanip>
#include "./single_particle_potential.h"

int main() {
    std::cout << "\n\n*********************************************************************"
	    << "\n* This code evaluates the potential parameters A, B, and tau for single particle potentials *"
	    << "\n*********************************************************************\n"
	    << std::endl;

    std::cout << std::fixed;
    std::cout << std::setprecision(6);

    try{

      potential_parameters pparams = {
        .degeneracy_g_ = 4.0,     //spin 1/2 for a fermion
        .mass_ = mNucleon,         // mean proton + neutron mass in MeV
        .saturation_density_ = nSat_in_MeV3, //density in MeV^3
        .binding_energy_ = -16.3, //in MeV
        .incompress_at_satdense_ = 240.0 //in MeV
      };

      potential_results parameter_results = {
        .A = -200.0,
        .B = 150.0,
        .tau = 3.0
      };

      const double tolerance = 1e-8;

      parameter_results = get_parameters(parameter_results, tolerance, pparams);

      get_single_particle_potential(parameter_results, pparams);

      std::cout << "\n\n*********************************************************************"
	    << "\n* Finished evaluating the potential parameters A, B, and tau for single particle potentials *"
	    << "\n*********************************************************************\n"
	    << std::endl;
    } catch (const std::runtime_error& e) {
    std::cerr << "Caught exception: " << e.what() << std::endl;
  }
}