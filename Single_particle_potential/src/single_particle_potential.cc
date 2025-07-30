#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_vector.h>
#include <iomanip>
#include "./core_statmech_functions.h"
#include "./constants.h"
#include "./single_particle_potential.h"

//Calculates the fermi momentum at density n_B
double fermi_momentum(void *p) {
    auto* par = static_cast<potential_parameters*>(p);
    return std::pow(6.0*std::pow(M_PI,2.0)
    *par->density_/par->degeneracy_g_,1.0/3.0);
}

//Calculates the fermi energy at density n_B
double fermi_energy(void *p) {
    auto* par = static_cast<potential_parameters*>(p);
    double p_momentum_fermi = fermi_momentum(p);
    return std::sqrt(par->mass_*par->mass_ + 
        p_momentum_fermi*p_momentum_fermi);
}

//Calculates the energy density of a spherical, ideal, 
//noninteracting fermi gas at density n_B at T=0
double energy_density_Fermigas(void *p){
    auto* par = static_cast<potential_parameters*>(p);

    //calculates the fermi momentum using the density, n
    double p_momentum_fermi = fermi_momentum(p);

    double fermi_energy_epsilon = fermi_energy(p);

    double term1 = 2.0*fermi_energy_epsilon*fermi_energy_epsilon
    *fermi_energy_epsilon*p_momentum_fermi;

    double term2 = fermi_energy_epsilon*p_momentum_fermi*par->mass_*par->mass_;

    double term3 = par->mass_*par->mass_*par->mass_*par->mass_
    *std::log((fermi_energy_epsilon + p_momentum_fermi)/par->mass_);

    //Calculates energy density at T = 0
    return (par->degeneracy_g_ / (16.0 * M_PI*M_PI))
    * (term1 - term2 - term3);
}

//Calculates the fermi momentum at saturation density in MeV^3
double fermi_momentum_satdense(void *p) {
    auto* par = static_cast<potential_parameters*>(p);
    return std::pow(6.0*std::pow(M_PI,2.0)
    *nSat_in_MeV3/par->degeneracy_g_,1.0/3.0);
}

//Calculates the fermi energy at saturation density in MeV^3
double fermi_energy_satdense(void *p) {
    auto* par = static_cast<potential_parameters*>(p);
    double p_momentum_fermi = fermi_momentum_satdense(p);
    return std::sqrt(par->mass_*par->mass_ + 
        p_momentum_fermi*p_momentum_fermi);
}

//Calculates the energy density of a spherical, ideal, 
//noninteracting fermi gas at saturation density in MeV^3 at T=0
double energy_density_Fermigas_satdense(void *p){
    auto* par = static_cast<potential_parameters*>(p);

    double p_momentum_fermi = fermi_momentum_satdense(p);

    double fermi_energy_epsilon = fermi_energy_satdense(p);

    //3 terms to the equation, separated out for readability

    double term1 = 2.0*fermi_energy_epsilon*fermi_energy_epsilon
    *fermi_energy_epsilon*p_momentum_fermi;

    double term2 = fermi_energy_epsilon*p_momentum_fermi*par->mass_*par->mass_;

    double term3 = par->mass_*par->mass_*par->mass_*par->mass_
    *std::log((fermi_energy_epsilon + p_momentum_fermi)/par->mass_);

    return (par->degeneracy_g_ / (16.0 * M_PI*M_PI))
    * (term1 - term2 - term3);
}

//Contains equations 25, 30, and 33 from Interactions in nuclear matter
//Calculates the conditions that all three parameters must satisfy
int conditions(const gsl_vector* x, void* p,
					     gsl_vector* f) {
    auto* par = static_cast<potential_parameters*>(p);
    
    //pulling in the current guess for each parameter
    double A = gsl_vector_get(x,0);
    double B = gsl_vector_get(x,1);
    double tau = gsl_vector_get(x,2);

    //Calculating the components of each condition
    const double energy_density_FG = energy_density_Fermigas(p);
    const double energy_density_FG_satdense = 
        energy_density_Fermigas_satdense(p);
    const double energy_epsilon_Fermi_satdense = fermi_energy_satdense(p);
    const double p_momentum_Fermi_satdense = fermi_momentum_satdense(p); 

    //Equation 25, written out explicitly and equal to zero
    //Condition 1 is at density n_B
    double condition1 = (energy_density_FG/nSat_in_MeV3) + ((A/2.0) + (B/tau))
     - par->mass_ - par->binding_energy_;

    gsl_vector_set(f, 0, condition1);

    //Equation 30, written out explicitly and equal to zero
    //Condition 2 is at saturation density
    double condition2 = (energy_epsilon_Fermi_satdense/nSat_in_MeV3)
        - (energy_density_FG_satdense/(nSat_in_MeV3*nSat_in_MeV3)) 
        + ((A/2.0)*(1/nSat_in_MeV3)) + B*((tau - 1.0)/tau)*(1/nSat_in_MeV3);

    gsl_vector_set(f, 1, condition2);

    //Equation 33, written out explicitly and equal to zero
    //Condition 3 is at saturation density

    //Term 1 is the derivative of P_FG 
    //Calculated using Mathematica for time, numerically identical to given solution
    double term1 = (par->degeneracy_g_*p_momentum_Fermi_satdense*
        p_momentum_Fermi_satdense*p_momentum_Fermi_satdense*
        p_momentum_Fermi_satdense*p_momentum_Fermi_satdense)
        /(18.0*energy_epsilon_Fermi_satdense*M_PI*M_PI*nSat_in_MeV3);

    double term2 = A;

    double term3 = B*(tau-1.0);

    double K = 9.0 * (term1 + term2 + term3);

    double condition3 = K - par->incompress_at_satdense_;

    gsl_vector_set(f, 2, condition3);

    //Makes sure none of the conditions break
    if (!std::isfinite(condition1)) {
        std::cerr << "NaN or Inf in condition 1" << std::endl;
    }
    else if (!std::isfinite(condition2)) {
        std::cerr << "NaN or Inf in condition 2" << std::endl;
    }
    else if (!std::isfinite(condition3)) {
        std::cerr << "NaN or Inf in condition 3" << std::endl;
    }
    
    return GSL_SUCCESS;
}

//Uses GSL multiroot to calculate the parameters A, B, and tau for the 
//single particle potential 
potential_results get_parameters(potential_results parameter_results,
    double tolerance, potential_parameters pparams) {
    //Number of parameters being solved for
    const size_t dim = 3;

    gsl_multiroot_function F;

    F.f = &conditions;
    F.n = dim;
    F.params = &pparams;
        
    //Sets the parameters vector to the initial guesses for A, B, and tau
    gsl_vector* x = gsl_vector_alloc(dim);
    gsl_vector_set(x, 0, parameter_results.A);
    gsl_vector_set(x, 1, parameter_results.B);
    gsl_vector_set(x, 2, parameter_results.tau);

    //Sets the multiroot solver type
    const gsl_multiroot_fsolver_type* T = gsl_multiroot_fsolver_hybrids;
    gsl_multiroot_fsolver* s = gsl_multiroot_fsolver_alloc(T, dim);
    gsl_multiroot_fsolver_set(s, &F, x);

    int test_status;

    size_t i = 0, iMAX = 1000; //starting iteration, max number of iterations


    do {
        //iterates the root solver
        i++;
        test_status = gsl_multiroot_fsolver_iterate(s);
        //Checks for a failed iteration
        if (test_status != GSL_SUCCESS && test_status != GSL_CONTINUE) {
            std::cerr << "Iteration " << i << " failed: " 
            << gsl_strerror(test_status) << std::endl;
            break;
        }
        //Checks if the root solver has reached the roots
        test_status = gsl_multiroot_test_residual(s->f, tolerance);
    } while (test_status == GSL_CONTINUE && i < iMAX);

    //Sets the results struct to the corresponding roots
    parameter_results.A = gsl_vector_get(s->x,0);
    parameter_results.B = gsl_vector_get(s->x,1);
    parameter_results.tau = gsl_vector_get(s->x,2);
    
    //Tells you what the roots are when they are solved for
    if (test_status == GSL_SUCCESS) {
        std::cout << "Single particle potential parameters for V("
        << pparams.density_ <<") at T = 0: A = " 
        << parameter_results.A << " B = " << parameter_results.B << " tau = " 
        << parameter_results.tau << " MeV" << std::endl;
    }
    //Tells you if the root solver fails to converge
    else {
        std::cerr << "Failed to converge:" << gsl_strerror(test_status) 
            << std::endl;
    }
    //Returns the solved roots
    return parameter_results;
}

//Calculates the single particle potential with the calculated parameters
double get_single_particle_potential(potential_results parameter_results,
    potential_parameters pparams) {
    //Equation 15 in Interactions in Nuclear matter
    double potential = parameter_results.A*(pparams.density_/nSat_in_MeV3) +
     parameter_results.B*
     std::pow(pparams.density_/nSat_in_MeV3,parameter_results.tau - 1.0);

    std::cout << "Single particle potential for V(n_B = " << pparams.density_ 
    << ") = " << potential << " MeV" << std::endl;
    //Returns the potential for possible future use
    return potential; 
}
