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

    try{

      potential_parameters pparams = {
        .degeneracy_g_ = 2.0,     //spin 1/2 for a fermion
        .mass_ = mNucleon,         // mean proton + neutron mass in MeV
        .density_ = nSat_in_MeV3, //density in MeV^3
        .binding_energy_ = -16.3, //in MeV
        .incompress_at_satdense_ = 240.0 //in MeV
      };

      potential_results parameter_results = {
        .A = -218.699,
        .B = 166.215,
        .tau = 2.33376
      };

      double tolerance = 1e-2;

      get_parameters(parameter_results, tolerance, pparams);

      std::cout << "\n\n*********************************************************************"
	    << "\n* Finished evaluating the potential parameters A, B, and tau for single particle potentials *"
	    << "\n*********************************************************************\n"
	    << std::endl;
    } catch (const std::runtime_error& e) {
    std::cerr << "Caught exception: " << e.what() << std::endl;
  }
}