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


        std::cout << "\n\n*********************************************************************"
	      << "\n* Finished evaluating the potential parameters A, B, and tau for single particle potentials *"
	      << "\n*********************************************************************\n"
	      << std::endl;
    } catch (const std::runtime_error& e) {
    std::cerr << "Caught exception: " << e.what() << std::endl;
  }
}