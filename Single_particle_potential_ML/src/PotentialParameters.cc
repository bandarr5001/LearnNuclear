#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_vector.h>
#include <iomanip>
#include <fstream>
#include "./read_and_write_4c.h"
#include "./single_particle_potential.h"

int main() {
    std::cout << "\n\n*********************************************************************"
	    << "\n* This code evaluates the potential parameters A, B, and tau for single particle potentials *"
	    << "\n*********************************************************************\n"
	    << std::endl;

    std::cout << std::fixed;
    std::cout << std::setprecision(6);

    double incompress_array[100];
    double A_array[100];
    double B_array[100];
    double tau_array[100];

    try{
      potential_parameters pparams = {
        .degeneracy_g_ = 4.0,     //spin 1/2 for a fermion
        .mass_ = mNucleon,         // mean proton + neutron mass in MeV
        .density_ = nSat_in_MeV3, //density in MeV^3
        .binding_energy_ = -16.3, //in MeV
        .incompress_at_satdense_ = 180.0 //in MeV
      };

      potential_results parameter_results = {
        .A = -200.0,
        .B = 150.0,
        .tau = 3.0
      };

      double tolerance = 1e-8;
      double step = (420.0-180.0)/99.0;
      for(int i=0; i<100; i++) {
      	parameter_results = get_parameters(parameter_results, tolerance, pparams);
	std::cout << pparams.incompress_at_satdense_ << std::endl;
        std::cout << parameter_results.A << std::endl;        
        std::cout << parameter_results.B << std::endl;          
        std::cout << parameter_results.tau << std::endl;          

	incompress_array[i] = pparams.incompress_at_satdense_;
        A_array[i] = parameter_results.A;
        B_array[i] = parameter_results.B;
        tau_array[i] = parameter_results.tau;

      	get_single_particle_potential(parameter_results, pparams);

      	std::cout << "\n\n*********************************************************************"
	   	 << "\n* Finished evaluating the potential parameters A, B, and tau for single particle potentials *"
	   	 << "\n*********************************************************************\n"
		    << std::endl;
        pparams.incompress_at_satdense_ = pparams.incompress_at_satdense_ + step;
      
     }} catch (const std::runtime_error& e) {
    	    std::cerr << "Caught exception: " << e.what() << std::endl;
        }
     write_array_4c(100, incompress_array, A_array, B_array, tau_array);
}
