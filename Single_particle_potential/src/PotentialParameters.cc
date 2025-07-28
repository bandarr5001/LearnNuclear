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


      // Thermodynamic parameters
      Parameters params = {
        .degeneracy_g_ = 2.0,     //spin 1/2 for a fermion
        .mass_ = 938.918,         // mean proton + neutron mass in MeV
        .density_ = nSat_in_MeV3, //density in MeV^3
        .temperature_ = 200.0     //temperature T in MeV  
      };

      potential_parameters pparams = {
        .density_ = nSat_in_MeV3,
        .binding_energy_ = -16.3, //in MeV
        .incompress_at_satdense_ = 240.0 //in MeV
      };

      double guess_A =  -200.0;
      double guess_B = 150.0;
      double guess_tau = 2.0;
      double tolerance = 1e-6;

      get_parameters(guess_A, guess_B, guess_tau, tolerance, pparams);

      std::cout << "\n\n*********************************************************************"
	    << "\n* Finished evaluating the potential parameters A, B, and tau for single particle potentials *"
	    << "\n*********************************************************************\n"
	    << std::endl;
    } catch (const std::runtime_error& e) {
    std::cerr << "Caught exception: " << e.what() << std::endl;
  }
}