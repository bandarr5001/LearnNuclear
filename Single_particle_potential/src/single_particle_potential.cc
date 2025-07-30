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

// Constructor
FermiWrapper::FermiWrapper(const potential_parameters& p) : params(p) {}

//Calculates the fermi momentum at density n_B
double FermiWrapper::fermi_momentum() const {
    return std::pow(6.0*std::pow(M_PI,2.0)*
    params.density_/params.degeneracy_g_,1.0/3.0);
}

//Calculates the fermi energy at density n_B
double FermiWrapper::fermi_energy() const {
    double p_momentum_fermi = fermi_momentum();
    return std::sqrt(params.mass_*params.mass_ + 
        p_momentum_fermi*p_momentum_fermi);
}

//Calculates the energy density of a spherical, ideal, 
//noninteracting fermi gas at density n_B at T=0
double FermiWrapper::energy_density_Fermigas() const {
    //calculates the fermi momentum using the density, n
    double p_F = fermi_momentum();
    double epsilon_F = fermi_energy();

    double term1 = 2.0*epsilon_F*epsilon_F*epsilon_F*p_F;

    double term2 =epsilon_F*p_F*params.mass_*params.mass_;

    double term3 = params.mass_*params.mass_*params.mass_*params.mass_
    *std::log((epsilon_F + p_F)/params.mass_);

    //Calculates energy density at T = 0
    return (params.degeneracy_g_ / (16.0 * M_PI*M_PI))
    * (term1 - term2 - term3);
}

//Calculates the fermi momentum at saturation density in MeV^3
double FermiWrapper::fermi_momentum_satdense() const {
    return std::pow(6.0*std::pow(M_PI,2.0)
    *nSat_in_MeV3/params.degeneracy_g_,1.0/3.0);
}

//Calculates the fermi energy at saturation density in MeV^3
double FermiWrapper::fermi_energy_satdense() const {
    double p_F = fermi_momentum_satdense();
    return std::sqrt(params.mass_*params.mass_ + p_F*p_F);
}

//Calculates the energy density of a spherical, ideal, 
//noninteracting fermi gas at saturation density in MeV^3 at T=0
double FermiWrapper::energy_density_Fermigas_satdense() const {
    double p_F = fermi_momentum_satdense();
    double epsilon_F = fermi_energy_satdense();

    //3 terms to the equation, separated out for readability

    double term1 = 2.0*epsilon_F*epsilon_F*epsilon_F*p_F;

    double term2 = epsilon_F*p_F*params.mass_*params.mass_;

    double term3 = params.mass_*params.mass_*params.mass_*params.mass_
    *std::log((epsilon_F + p_F)/params.mass_);

    return (params.degeneracy_g_ / (16.0 * M_PI*M_PI))
    * (term1 - term2 - term3);
}

//static callbacks
double FermiWrapper::fermi_momentum_callback(void* p) {
    auto* self = static_cast<FermiWrapper*>(p);
    return self->fermi_momentum();
}

double FermiWrapper::energy_density_callback(void* p) {
    auto* self = static_cast<FermiWrapper*>(p);
    return self->energy_density_Fermigas();
}

double FermiWrapper::fermi_momentum_callback_satdense(void* p) {
    auto* self = static_cast<FermiWrapper*>(p);
    return self->fermi_momentum_satdense();
}

double FermiWrapper::fermi_energy_callback_satdense(void* p) {
    auto* self = static_cast<FermiWrapper*>(p);
    return self->fermi_energy_satdense();
}

//Contains equations 25, 30, and 33 from Interactions in nuclear matter
//Calculates the conditions that all three parameters must satisfy
int conditions(const gsl_vector* x, void* p,
					     gsl_vector* f) {
    auto* wrapper = static_cast<FermiWrapper*>(p);
    auto* par = static_cast<potential_parameters*>(p);
    
    //pulling in the current guess for each parameter
    double A = gsl_vector_get(x,0);
    double B = gsl_vector_get(x,1);
    double tau = gsl_vector_get(x,2);

    //Calculating the components of each condition
    const double energy_density_FG = wrapper->energy_density_Fermigas();
    const double energy_density_FG_satdense = 
        wrapper->energy_density_Fermigas_satdense();
    const double energy_epsilon_Fermi_satdense = wrapper->fermi_energy_satdense();
    const double p_momentum_Fermi_satdense = wrapper->fermi_momentum_satdense(); 

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
