#include "./potential_parameters.h"

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

#include "./Fermi_gas_functions_Tzero.h"



// Condition for the value of the binding energy
// (equation 25 in "Interactions in nuclear matter" notes)
double PotentialParameters::binding_energy_condition(double A, double B, double tau) {
  const double energy_density_FG =
      Fermi_gas_energy_density(degeneracy_, mass_, saturation_density_);
  const double energy_per_particle_at_saturation =
    energy_density_FG/saturation_density_ + A/2.0 + B/tau - mass_;

  return energy_per_particle_at_saturation - binding_energy_;
}



// Condition for the location of the minimum in the energy per particle, i.e., value of
// the saturation density (equation 30 in "Interactions in nuclear matter" notes)
double PotentialParameters::energy_minimum_condition(double A, double B, double tau) {
  const double energy_epsilon_Fermi =
    fermi_energy(degeneracy_, mass_, saturation_density_);
  const double energy_density_FG =
      Fermi_gas_energy_density(degeneracy_, mass_, saturation_density_);

  return (energy_epsilon_Fermi/saturation_density_)
    - (energy_density_FG / std::pow(saturation_density_, 2.0)) 
    + ( A/2.0 + B * (tau - 1.0)/tau ) / saturation_density_;
}



// Condition for the value of the incompressibility at saturation density
// (equation 33 in "Interactions in nuclear matter" notes)
double PotentialParameters::incompressibility_condition(double A, double B, double tau) {
  const double d_PFG_dnB = d_PFG_d_nB(degeneracy_, mass_, saturation_density_) ;
  // derivative of the interaction contribution to pressure at saturation density
  const double d_Pint_dnB = A + B * (tau - 1.0);

  return 9.0 * ( d_PFG_dnB + d_Pint_dnB ) - incompress_at_satdense_;
}



// Evaluates equations 25, 30, and 33 from "Interactions in nuclear matter" notes, which
// are conditions that all three potential parameters A, B, tau must satisfy
int PotentialParameters::conditions(const gsl_vector* x, void* p, gsl_vector* f) {
  // This member function is declared static, which means that "this" (the pointer to the
  // current class instance) is not available in this scope. As a result, we cannot
  // directly access member variables like mass_ or call member functions like
  // Fermi_momentum().
  // However, the GSL interface allows us to pass a user-defined pointer via the void* p
  // argument. In this case, we pass a pointer to the current class instance ("this"),
  // cast to void*. To recover access to class members and functions, we static_cast the
  // void* p back to a pointer to our class, which gives us the same members and methods
  // we would have via "this".
  // Since "this" is a reserved keyword and not usable here, we conventionally call the
  // recovered pointer "self" to imply that it plays the same role.
    auto* self = static_cast<PotentialParameters*>(p);

    
    // Pull in the current guess for each parameter
    double A = gsl_vector_get(x,0);
    double B = gsl_vector_get(x,1);
    double tau = gsl_vector_get(x,2);

    /////////////////////////////////////////
    // Evaluate the conditions for current parameter guesses   
    const double condition1 = self->binding_energy_condition(A, B, tau);
    gsl_vector_set(f, 0, condition1);
    
    const double condition2 = self->energy_minimum_condition(A, B, tau);
    gsl_vector_set(f, 1, condition2);
    
    const double condition3 = self->incompressibility_condition(A, B, tau);
    gsl_vector_set(f, 2, condition3);

    /////////////////////////////////////////
    // Make sure none of the conditions break
    if (!std::isfinite(condition1)) {
      throw std::runtime_error("NaN or Inf in condition 1");
    }
    else if (!std::isfinite(condition2)) {
      throw std::runtime_error("NaN or Inf in condition 2");
    }
    else if (!std::isfinite(condition3)) {
      throw std::runtime_error("NaN or Inf in condition 3");
    }
    
    return GSL_SUCCESS;
}



// Uses GSL multiroot to find the parameters A, B, and tau of the single-particle
// potential for given nuclear matter properties (saturation_density_, binding_energy_,
// incompress_at_satdense_)
void PotentialParameters::get_parameters(double tolerance, double A_guess, double B_guess,
					 double tau_guess) {
    // Number of parameters being solved for
    const size_t dim = 3;

    /////////////////////////////////////////
    // Initialize the struct needed for multidimensional root finding
    gsl_multiroot_function Root_function_struct;

    Root_function_struct.f = &PotentialParameters::conditions;
    Root_function_struct.n = dim;
    // static_cast "this" (the pointer to the current instance of the class) onto void*,
    // required by GSL:
    Root_function_struct.params = static_cast<void*>(this);

    /////////////////////////////////////////
    // Initialize the roots vector to the initial guesses for A, B, and tau
    gsl_vector* roots = gsl_vector_alloc(dim);
    // hard-code initial guesses for now
    gsl_vector_set(roots, 0, A_guess);
    gsl_vector_set(roots, 1, B_guess);
    gsl_vector_set(roots, 2, tau_guess);

    /////////////////////////////////////////
    // Set the multiroot solver type
    const gsl_multiroot_fsolver_type* type = gsl_multiroot_fsolver_hybrids;
    gsl_multiroot_fsolver* solver = gsl_multiroot_fsolver_alloc(type, dim);
    gsl_multiroot_fsolver_set(solver, &Root_function_struct, roots);

    /////////////////////////////////////////
    // Find roots
    int test_status{};
    int i{};   // to count iterations
    const int iMAX = 1000;   // max number of iterations

    do {
        // IteratE the root solver
        i++;
        test_status = gsl_multiroot_fsolver_iterate(solver);
        // Checs for a failed iteration
        if (test_status != GSL_SUCCESS && test_status != GSL_CONTINUE) {
	  std::cerr << "Iteration " << i << " failed: " 
		    << gsl_strerror(test_status) << std::endl;
	  break;
        }
        //Checks if the root solver has reached the roots
        test_status = gsl_multiroot_test_residual(solver->f, tolerance);
    } while (test_status == GSL_CONTINUE && i < iMAX);

    
    // If roots are found, set the potential parameters and print out the values
    if (test_status == GSL_SUCCESS) {
      // Set the potential parameter members of the class
      A_ = gsl_vector_get(solver->x, 0);
      B_ = gsl_vector_get(solver->x, 1);
      tau_ = gsl_vector_get(solver->x, 2);
      
      std::cout << "Single particle potential parameters: A = " << A_
		<< " [MeV], B = " << B_ << " [MeV], tau = " << tau_ << " [1]"
		<< std::endl;
    }
    // Notify the user if the root solver fails to converge
    else {
      std::cerr << "Failed to converge:" << gsl_strerror(test_status) 
		<< std::endl;
    }    
    
    // FREE POINTERS!!!!!!
    gsl_multiroot_fsolver_free(solver);
    gsl_vector_free(roots);
}



// Calculates the single particle potential with the calculated parameters
double PotentialParameters::get_single_particle_potential(double density) {
  // Equation 15 in "Interactions in nuclear matter" notes
  double potential = A_ * (density/saturation_density_) +
    B_ * std::pow(density/saturation_density_, tau_ - 1.0);

  std::cout << "Single particle potential for V(n_B = " << density 
	    << ") = " << potential << " MeV" << std::endl;
    
  //Returns the potential for possible future use
  return potential; 
}
