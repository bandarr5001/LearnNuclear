#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>

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

  // Use try{} catch{} to handle errors
  try {

  
    ////////////////////////////////////////////////////////////////////////////////////////
    // Define constants

    // Thermodynamic parameters
    Parameters params = {
      .degeneracy_g_ = 2.0,     //spin 1/2 for a fermion
      .mass_ = 938.918,         // mean proton + neutron mass in MeV
      .density_ = nSat_in_MeV3, //density in MeV^3
      .temperature_ = 200.0     //temperature T in MeV  
    };

    std::cout << "\nLooking for the chemical potential muB for:"
	      << "\n              degeneracy = " << params.degeneracy_g_
	      << "\n                    mass = " << params.mass_ << " [MeV]"
	      << "\n          baryon density = " << params.density_ << " [MeV^3]"
	      << " = " << params.density_ / std::pow(hBarC, 3.0) << " [fm^-3]"
	      << " = " << params.density_ / nSat_in_MeV3 << " [n_0]"
	      << "\n             temperature = " << params.temperature_ << " [MeV]"
	      << std::endl;



    // Root-finding algorithm parameters
    const double tolerance = 1e-6;
    const double bisection_interval_low = 100.0;
    const double bisection_interval_high = 1500.0;
    const double muB_guess = 900.0;

    std::cout << "\nRoot-finding algorithms will use:"
	      << "\n  bisection_interval_low = " << bisection_interval_low << " [MeV]"
	      << "\n bisection_interval_high = " << bisection_interval_high << " [MeV]"
	      << "\n               muB_guess = " << muB_guess << " [MeV]" << "\n\n"
	      << "\n               tolerance = " << tolerance << std::endl;


  
    ////////////////////////////////////////////////////////////////////////////////////////
    // Solve the root equation for the chemical potential using a number of methods

    ///////////////////////////////////////////
    // My bisection method
    auto start_bisect = std::chrono::high_resolution_clock::now();

    get_chem_potent_mu_fermi_bisection
      (tolerance, bisection_interval_low, bisection_interval_high, params);

    auto stop_bisect = std::chrono::high_resolution_clock::now();
    auto duration_bisect = std::chrono::duration_cast<std::chrono::milliseconds>
      (stop_bisect - start_bisect);
    std::cout << "My bisection method took "
	      << duration_bisect.count()
	      << " milliseconds\n" << std::endl;


  
    ///////////////////////////////////////////
    // My Newton-Raphson method
    auto start_newt = std::chrono::high_resolution_clock::now();

    get_chem_potent_mu_fermi_newton(tolerance, muB_guess, params);

    auto stop_newt = std::chrono::high_resolution_clock::now();
    auto duration_newt = std::chrono::duration_cast<std::chrono::milliseconds>
      (stop_newt - start_newt);
    std::cout << "My Newton-Raphson method took "
	      << duration_newt.count() << " milliseconds\n" << std::endl;


  
    ///////////////////////////////////////////
    // My Newton-Raphson method with finite-difference derivative
    auto start_newt_num = std::chrono::high_resolution_clock::now();

    get_chem_potent_mu_fermi_newton_numerical(tolerance, muB_guess, params);

    auto stop_newt_num = std::chrono::high_resolution_clock::now();
    auto duration_newt_num = std::chrono::duration_cast<std::chrono::milliseconds>
      (stop_newt_num - start_newt_num);
    std::cout << "My numerical derivative Newton-Raphson method took "
	      << duration_newt_num.count() << " milliseconds\n" << std::endl;


  
    ///////////////////////////////////////////
    // GSL bisection
    auto start_gsl_bisect = std::chrono::high_resolution_clock::now();

    gsl_chem_potent_fermi_bisection(tolerance, params);

    auto stop_gsl_bisect = std::chrono::high_resolution_clock::now();
    auto duration_gsl_bisect = std::chrono::duration_cast<std::chrono::milliseconds>
      (stop_gsl_bisect - start_gsl_bisect);
    std::cout << "GSL's bisection method took "
	      << duration_gsl_bisect.count() << " milliseconds\n" << std::endl;


  
    ///////////////////////////////////////////
    // GSL Newton-Rapshon method
    auto start_gsl_newt = std::chrono::high_resolution_clock::now();

    gsl_chem_potent_fermi_newton(tolerance, muB_guess, params);

    auto stop_gsl_newt = std::chrono::high_resolution_clock::now();
    auto duration_gsl_newt = std::chrono::duration_cast<std::chrono::milliseconds>
      (stop_gsl_newt - start_gsl_newt);
    std::cout << "GSL's Newton-Raphson method took " << 
      duration_gsl_newt.count() << " milliseconds\n" << std::endl;


  
    ///////////////////////////////////////////
    // GSL Newton-Rapshon method with finite-difference derivative
    auto start_gsl_newt_num = std::chrono::high_resolution_clock::now();

    gsl_chem_potent_fermi_newton_numerical(tolerance, muB_guess, params);

    auto stop_gsl_newt_num = std::chrono::high_resolution_clock::now();
    auto duration_gsl_newt_num = std::chrono::duration_cast<std::chrono::milliseconds>
      (stop_gsl_newt_num - start_gsl_newt_num);
    std::cout << "GSL's Newton-Raphson method with numerical derivative took " << 
      duration_gsl_newt_num.count() << " milliseconds\n" << std::endl;


  
    ///////////////////////////////////////////
    // GSL multiroot hybrid solver

    // Fermi gas
    auto start_gsl_multi_fermi = std::chrono::high_resolution_clock::now();

    gsl_chem_potent_fermi_multiroot(muB_guess, tolerance, params);

    auto stop_gsl_multi_fermi = std::chrono::high_resolution_clock::now();
    auto duration_gsl_multi_fermi = std::chrono::duration_cast<std::chrono::milliseconds>
      (stop_gsl_multi_fermi - start_gsl_multi_fermi);
    std::cout << "GSL's multiroot method for Fermi gas took "
	      << duration_gsl_multi_fermi.count() << " milliseconds\n" << std::endl;

    // Bose gas
    auto start_gsl_multi_bose = std::chrono::high_resolution_clock::now();

    get_chem_potent_mu_bose_multiroot(muB_guess, tolerance, params);

    auto stop_gsl_multi_bose = std::chrono::high_resolution_clock::now();
    auto duration_gsl_multi_bose = std::chrono::duration_cast<std::chrono::milliseconds>
      (stop_gsl_multi_bose - start_gsl_multi_bose);
    std::cout << "GSL's multiroot method for Bose gas took "
	      << duration_gsl_multi_bose.count() << " milliseconds\n" << std::endl;

    // Boltzmann gas
    auto start_gsl_multi_boltz = std::chrono::high_resolution_clock::now();

    get_chem_potent_mu_boltzmann_multiroot(muB_guess, tolerance, params);

    auto stop_gsl_multi_boltz = std::chrono::high_resolution_clock::now();
    auto duration_gsl_multi_boltz = std::chrono::duration_cast<std::chrono::milliseconds>
      (stop_gsl_multi_boltz - start_gsl_multi_boltz);
    std::cout << "GSL's multiroot method for Boltzmann gas took "
	      << duration_gsl_multi_boltz.count() << " milliseconds\n" << std::endl;


  
    ////////////////////////////////////////////////////////////////////////////////////////
    // // Let the user know everything went well
    std::cout << "\n\n*********************************************************************"
	      << "\n* Finished evaluating the numerical cost of root-finding algorithms *"
	      << "\n*********************************************************************\n"
	      << std::endl;
  } catch (const std::runtime_error& e) {
    std::cerr << "Caught exception: " << e.what() << std::endl;
  }
 
  return 0;
}
