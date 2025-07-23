#include <iostream>
#include <cmath>
#include <cstdlib>
#include <chrono>
#include "./fermi_gas.h"
#include "./chemical_potential.h"

int main() {
    Parameters params;
    params.degeneracy_g_ = 2.0; //spin 1/2 for a fermion
    params.mass_ = 938.918; // mean proton + neutron mass in MeV
    params.temperature_ = 200.0; //temperature T in MeV
    params.density_ = 122856.691; //density n_0
    double tolerance = 1e-6;
    double guess = 900.0;

    std::cout << "Tolerance is: " << tolerance << std::endl;

    std::cout << "Guess is: " << guess << std::endl;

    auto start_bisect = std::chrono::high_resolution_clock::now();

    get_chem_potent_mu_fermi_bisection(tolerance,500.0,1500.0,params);

    auto stop_bisect = std::chrono::high_resolution_clock::now();
    auto duration_bisect = 
    std::chrono::duration_cast<std::chrono::milliseconds>(stop_bisect 
        - start_bisect);
    std::cout << "My bisection method took " << duration_bisect.count()
     << " milliseconds\n" << std::endl;

    auto start_newt = std::chrono::high_resolution_clock::now();

    get_chem_potent_mu_fermi_newton(tolerance, guess, params);

    auto stop_newt = std::chrono::high_resolution_clock::now();
    auto duration_newt = 
    std::chrono::duration_cast<std::chrono::milliseconds>(stop_newt 
        - start_newt);
    std::cout << "My Newton-Raphson method took " << duration_newt.count() << 
    " milliseconds\n" << std::endl;


    auto start_newt_num = std::chrono::high_resolution_clock::now();

    get_chem_potent_mu_fermi_newton_numerical(tolerance, guess, params);

    auto stop_newt_num = std::chrono::high_resolution_clock::now();
    auto duration_newt_num = 
    std::chrono::duration_cast<std::chrono::milliseconds>(stop_newt_num 
        - start_newt_num);
    std::cout << "My numerical derivative Newton-Raphson method took " << duration_newt_num.count() << 
    " milliseconds\n" << std::endl;


    auto start_gsl_bisect = std::chrono::high_resolution_clock::now();

    gsl_chem_potent_fermi_bisection(tolerance, params);

    auto stop_gsl_bisect = std::chrono::high_resolution_clock::now();
    auto duration_gsl_bisect = 
    std::chrono::duration_cast<std::chrono::milliseconds>(stop_gsl_bisect 
        - start_gsl_bisect);
    std::cout << "GSL's bisection method took " << duration_gsl_bisect.count()
     << " milliseconds\n" << std::endl;



    auto start_gsl_newt = std::chrono::high_resolution_clock::now();

    gsl_chem_potent_fermi_newton(tolerance, guess, params);

    auto stop_gsl_newt = std::chrono::high_resolution_clock::now();
    auto duration_gsl_newt = 
    std::chrono::duration_cast<std::chrono::milliseconds>(stop_gsl_newt - 
        start_gsl_newt);
    std::cout << "GSL's Newton-Raphson method took " << 
    duration_gsl_newt.count() << " milliseconds\n" << std::endl;


    auto start_gsl_newt_num = std::chrono::high_resolution_clock::now();

    gsl_chem_potent_fermi_newton_numerical(tolerance, guess, params);

    auto stop_gsl_newt_num = std::chrono::high_resolution_clock::now();
    auto duration_gsl_newt_num = 
    std::chrono::duration_cast<std::chrono::milliseconds>(stop_gsl_newt_num - 
        start_gsl_newt_num);
    std::cout << "GSL's Newton-Raphson method with numerical derivative took " << 
    duration_gsl_newt_num.count() << " milliseconds\n" << std::endl;


    auto start_gsl_multi = std::chrono::high_resolution_clock::now();

    gsl_chem_potent_fermi_multiroot(guess, tolerance, params);

    auto stop_gsl_multi = std::chrono::high_resolution_clock::now();
    auto duration_gsl_multi = 
    std::chrono::duration_cast<std::chrono::milliseconds>(stop_gsl_multi - 
        start_gsl_multi);
    std::cout << "GSL's multiroot method took " << duration_gsl_multi.count() 
    << " milliseconds\n" << std::endl;

    auto start_bose_multi = std::chrono::high_resolution_clock::now();

    get_chem_potent_mu_bose_multiroot(guess, tolerance, params);

    auto stop_bose_multi = std::chrono::high_resolution_clock::now();
    auto duration_bose_multi = 
    std::chrono::duration_cast<std::chrono::milliseconds>(stop_bose_multi - 
        start_bose_multi);
    std::cout << "GSL's multiroot method took " << duration_bose_multi.count() 
    << " milliseconds\n" << std::endl;

    auto start_boltzmann_multi = std::chrono::high_resolution_clock::now();

    get_chem_potent_mu_boltzmann_multiroot(guess, tolerance, params);

    auto stop_boltzmann_multi = std::chrono::high_resolution_clock::now();
    auto duration_boltzmann_multi = 
    std::chrono::duration_cast<std::chrono::milliseconds>(stop_boltzmann_multi - 
        start_boltzmann_multi);
    std::cout << "GSL's multiroot method took " << duration_boltzmann_multi.count() 
    << " milliseconds\n" << std::endl;
}