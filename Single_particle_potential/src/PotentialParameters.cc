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
  std::cout << "\n*********************************************************************"
	    << "\n* This code finds nuclear matter potential parameters A, B, & tau   *"
	    << "\n*********************************************************************"
	    << "\n" << std::endl;

  std::cout << std::fixed;
  std::cout << std::setprecision(6);

  try{

    PotentialParameters parameter_set(4.0,     //spin 1/2 for a fermion
				      mNucleon,   // mean proton + neutron mass in MeV
				      nSat_in_MeV3, //density in MeV^3
				      -16.3, //in MeV
				      240.0 //in MeV
				      );
    const double tolerance = 1e-8;

    parameter_set.get_parameters(tolerance);

    parameter_set.get_single_particle_potential(nSat_in_MeV3);

    std::cout << "\n*********************************************************************"
	      << "\n* Finished evaluating the potential parameters A, B, & tau          *"
	      << "\n*********************************************************************"
	      << "\n"
	      << std::endl;
  } catch (const std::runtime_error& e) {
    std::cerr << "Caught exception: " << e.what() << std::endl;
  }
}
