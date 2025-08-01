#include <iomanip>   // to control the number of significant digits in std::cout
#include <iostream>  // for std::cout

#include "./constants.h"
#include "./potential_parameters.h"



int main() {
  std::cout << "\n*********************************************************************"
	    << "\n* This code finds nuclear matter potential parameters A, B, & tau   *"
	    << "\n*********************************************************************"
	    << "\n" << std::endl;

  std::cout << std::fixed;
  std::cout << std::setprecision(6);

  try{
    // Choose nuclear matter properties and initialize PotentialParameters object
    const double degeneracy_g = 4.0; // 2.0 = protons/neutrons only; 4.0 = nuclear matter
    const double incompressibility = 240.0;
    PotentialParameters parameter_set(degeneracy_g, mNucleon, nSat_in_MeV3, eBinding,
				      incompressibility);
    const double tolerance = 1e-8;

    // Find potential parameters reproducing nuclear matter properties
    parameter_set.get_parameters(tolerance);

    // Print out the single-particle potential at saturation density
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
