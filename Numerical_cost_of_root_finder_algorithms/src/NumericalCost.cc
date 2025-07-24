#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>

#include "./constants.h"
#include "./chemical_potential.h"
#include "./fermi_gas.h"



int main() {
  ////////////////////////////////////////////////////////////////////////////////////////
  // Let the user know what you will do
  std::cout << "\n\n*********************************************************************"
	    << "\n* This code evaluates the numerical cost of root-finding algorithms *"
	    << "\n*********************************************************************\n"
	    << std::endl;


  
  ////////////////////////////////////////////////////////////////////////////////////////
  // Define constants

  // Themordynamic parameters
  Parameters params = {
    .degeneracy_g_ = 2.0,     //spin 1/2 for a fermion
    .mass_ = 938.918,         // mean proton + neutron mass in MeV
    .density_ = nSat_in_MeV3, //density in MeV^3
    .temperature_ = 200.0     //temperature T in MeV  
  };

  std::cout << "\nLooking for the chemical potential muB for:"
	    << "\n     degeneracy = " << params.degeneracy_g_
	    << "\n           mass = " << params.mass_ << " [MeV]"
	    << "\n baryon density = " << params.density_ << " [MeV^3]"
	    << " = " << params.density_ / std::pow(hBarC, 3.0) << " [fm^-3]"
	    << " = " << params.density_ / nSat_in_MeV3 << " [n_0]"
	    << "\n    temperature = " << params.temperature_ << " [MeV]"
	    << std::endl;



  // Root-finding algorithm parameters
  const double tolerance = 1e-6;
  const double guess = 900.0;

  std::cout << "\nRoot-finding algorithms will use:"
	    << "\n      tolerance = " << tolerance
	    << "\n          guess = " << guess << " [MeV]" << "\n\n" << std::endl;


  
  ////////////////////////////////////////////////////////////////////////////////////////
  // Solve the root equation for the chemical potential using a number of methods

  ///////////////////////////////////////////
  // My bisection method
  auto start_bisect = std::chrono::high_resolution_clock::now();

  get_chem_potent_mu_fermi_bisection(tolerance,500.0,1500.0,params);

  auto stop_bisect = std::chrono::high_resolution_clock::now();
  auto duration_bisect = 
    std::chrono::duration_cast<std::chrono::milliseconds>(stop_bisect 
							  - start_bisect);
  std::cout << "My bisection method took " << duration_bisect.count()
	    << " milliseconds\n" << std::endl;


  
  ///////////////////////////////////////////
  // My Newton-Raphson method
  auto start_newt = std::chrono::high_resolution_clock::now();

  get_chem_potent_mu_fermi_newton(tolerance, guess, params);

  auto stop_newt = std::chrono::high_resolution_clock::now();
  auto duration_newt = 
    std::chrono::duration_cast<std::chrono::milliseconds>(stop_newt 
							  - start_newt);
  std::cout << "My Newton-Raphson method took " << duration_newt.count()
	    << " milliseconds\n" << std::endl;


  
  ///////////////////////////////////////////
  // My Newton-Raphson method with finite-difference derivative
  auto start_newt_num = std::chrono::high_resolution_clock::now();

  get_chem_potent_mu_fermi_newton_numerical(tolerance, guess, params);

  auto stop_newt_num = std::chrono::high_resolution_clock::now();
  auto duration_newt_num = 
    std::chrono::duration_cast<std::chrono::milliseconds>(stop_newt_num 
							  - start_newt_num);
  std::cout << "My numerical derivative Newton-Raphson method took " << duration_newt_num.count() << 
    " milliseconds\n" << std::endl;


  
  ///////////////////////////////////////////
  // GSL bisection
  auto start_gsl_bisect = std::chrono::high_resolution_clock::now();

  gsl_chem_potent_fermi_bisection(tolerance, params);

  auto stop_gsl_bisect = std::chrono::high_resolution_clock::now();
  auto duration_gsl_bisect = 
    std::chrono::duration_cast<std::chrono::milliseconds>(stop_gsl_bisect 
							  - start_gsl_bisect);
  std::cout << "GSL's bisection method took " << duration_gsl_bisect.count()
	    << " milliseconds\n" << std::endl;


  
  ///////////////////////////////////////////
  // GSL Newton-Rapshon method
  auto start_gsl_newt = std::chrono::high_resolution_clock::now();

  gsl_chem_potent_fermi_newton(tolerance, guess, params);

  auto stop_gsl_newt = std::chrono::high_resolution_clock::now();
  auto duration_gsl_newt = 
    std::chrono::duration_cast<std::chrono::milliseconds>(stop_gsl_newt - 
							  start_gsl_newt);
  std::cout << "GSL's Newton-Raphson method took " << 
    duration_gsl_newt.count() << " milliseconds\n" << std::endl;


  
  ///////////////////////////////////////////
  // GSL Newton-Rapshon method with finite-difference derivative
  auto start_gsl_newt_num = std::chrono::high_resolution_clock::now();

  gsl_chem_potent_fermi_newton_numerical(tolerance, guess, params);

  auto stop_gsl_newt_num = std::chrono::high_resolution_clock::now();
  auto duration_gsl_newt_num = 
    std::chrono::duration_cast<std::chrono::milliseconds>(stop_gsl_newt_num - 
							  start_gsl_newt_num);
  std::cout << "GSL's Newton-Raphson method with numerical derivative took " << 
    duration_gsl_newt_num.count() << " milliseconds\n" << std::endl;


  
  ///////////////////////////////////////////
  // GSL multiroot hybrid solver

  // Fermi gas
  auto start_gsl_multi_fermi = std::chrono::high_resolution_clock::now();

  gsl_chem_potent_fermi_multiroot(guess, tolerance, params);

  auto stop_gsl_multi_fermi = std::chrono::high_resolution_clock::now();
  auto duration_gsl_multi_fermi = 
    std::chrono::duration_cast<std::chrono::milliseconds>
    (stop_gsl_multi_fermi - start_gsl_multi_fermi);
  std::cout << "GSL's multiroot method for Fermi gas took "
	    << duration_gsl_multi_fermi.count() << " milliseconds\n" << std::endl;

  // Bose gas
  auto start_gsl_multi_bose = std::chrono::high_resolution_clock::now();

  get_chem_potent_mu_bose_multiroot(guess, tolerance, params);

  auto stop_gsl_multi_bose = std::chrono::high_resolution_clock::now();
  auto duration_gsl_multi_bose = 
    std::chrono::duration_cast<std::chrono::milliseconds>
    (stop_gsl_multi_bose - start_gsl_multi_bose);
  std::cout << "GSL's multiroot method for Bose gas took "
	    << duration_gsl_multi_bose.count() << " milliseconds\n" << std::endl;

  // Boltzmann gas
  auto start_gsl_multi_boltzmann = std::chrono::high_resolution_clock::now();

  get_chem_potent_mu_boltzmann_multiroot(guess, tolerance, params);

  auto stop_gsl_multi_boltzmann = std::chrono::high_resolution_clock::now();
  auto duration_gsl_multi_boltzmann = 
    std::chrono::duration_cast<std::chrono::milliseconds>
    (stop_gsl_multi_boltzmann - start_gsl_multi_boltzmann);
  std::cout << "GSL's multiroot method for Boltzmann gas took "
	    << duration_gsl_multi_boltzmann.count() << " milliseconds\n" << std::endl;


  
  ////////////////////////////////////////////////////////////////////////////////////////
  // // Let the user know everything went well
  std::cout << "\n\n*********************************************************************"
	    << "\n* Finished evaluating the numerical cost of root-finding algorithms *"
	    << "\n*********************************************************************\n"
	    << std::endl;
  return 0;
}
