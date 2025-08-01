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
// FermiWrapper::FermiWrapper(nuclear_properties* p) : props(p) {}

// Calculates the Fermi momentum at density n_B
double PotentialParameters::fermi_momentum() {
  return std::pow((6.0 * std::pow(M_PI,2.0) * this->saturation_density() /
		   this->degeneracy() ), 1.0/3.0);
}

//Calculates the fermi energy at density n_B
double PotentialParameters::fermi_energy() {
    double p_momentum_fermi = fermi_momentum();
    return std::sqrt(this->mass() * this->mass() + p_momentum_fermi * p_momentum_fermi);
}

//Calculates the energy density of a spherical, ideal, 
//noninteracting fermi gas at density n_B at T=0
double PotentialParameters::energy_density_Fermigas() {
    //calculates the fermi momentum using the density, n
    double p_F = fermi_momentum();
    double epsilon_F = fermi_energy();
    const double mass = this->mass();

    double term1 = 2.0*epsilon_F*epsilon_F*epsilon_F*p_F;

    double term2 =epsilon_F * p_F * std::pow(mass, 2.0);

    double term3 = std::pow(mass, 4.0) * std::log((epsilon_F + p_F)/mass);

    //Calculates energy density at T = 0
    return (this->degeneracy() / (16.0 * M_PI*M_PI))
    * (term1 - term2 - term3);
}

//static callbacks
double PotentialParameters::fermi_momentum_callback(void* p) {
    auto* self = static_cast<PotentialParameters*>(p);
    return self->fermi_momentum();
}

double PotentialParameters::energy_density_callback(void* p) {
    auto* self = static_cast<PotentialParameters*>(p);
    return self->energy_density_Fermigas();
}

//Contains equations 25, 30, and 33 from Interactions in nuclear matter
//Calculates the conditions that all three parameters must satisfy
int PotentialParameters::conditions(const gsl_vector* x, void* p,
				    gsl_vector* f) {
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

    
    //pulling in the current guess for each parameter
    double A = gsl_vector_get(x,0);
    double B = gsl_vector_get(x,1);
    double tau = gsl_vector_get(x,2);

    //Calculating the components of each condition
    const double energy_density_FG = self->energy_density_Fermigas();
    const double energy_epsilon_Fermi = self->fermi_energy();
    const double p_momentum_Fermi = self->fermi_momentum(); 

    //Equation 25, written out explicitly and equal to zero
    //Condition 1 is at density n_B
    double condition1 = (energy_density_FG/nSat_in_MeV3) + ((A/2.0) + (B/tau))
      - self->mass() - self->binding_energy();

    gsl_vector_set(f, 0, condition1);

    //Equation 30, written out explicitly and equal to zero
    //Condition 2 is at saturation density
    double condition2 = (energy_epsilon_Fermi/nSat_in_MeV3)
        - (energy_density_FG/(nSat_in_MeV3*nSat_in_MeV3)) 
        + ((A/2.0)*(1/nSat_in_MeV3)) + B*((tau - 1.0)/tau)*(1/nSat_in_MeV3);

    gsl_vector_set(f, 1, condition2);

    //Equation 33, written out explicitly and equal to zero
    //Condition 3 is at saturation density

    //Term 1 is the derivative of P_FG 
    //Calculated using Mathematica for time, numerically identical to given solution
    double term1 = (self->degeneracy()*p_momentum_Fermi*
        p_momentum_Fermi*p_momentum_Fermi*
        p_momentum_Fermi*p_momentum_Fermi)
        /(18.0*energy_epsilon_Fermi*M_PI*M_PI*nSat_in_MeV3);

    double term2 = A;

    double term3 = B*(tau-1.0);

    double K = 9.0 * (term1 + term2 + term3);

    double condition3 = K - self->incompress_at_satdense();

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
void PotentialParameters::get_parameters(double tolerance) {
    //Number of parameters being solved for
    const size_t dim = 3;

    gsl_multiroot_function F;

    F.f = &PotentialParameters::conditions;
    F.n = dim;
    // static_cast "this" (the pointer to the current instance of the class) onto void*,
    // required by GSL:
    F.params = static_cast<void*>(this);
        
    //Sets the parameters vector to the initial guesses for A, B, and tau
    gsl_vector* x = gsl_vector_alloc(dim);
    // hard-code initial guesses for now
    gsl_vector_set(x, 0, 200.0);
    gsl_vector_set(x, 1, 100.0);
    gsl_vector_set(x, 2, 3.0);

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
    A_ = gsl_vector_get(s->x,0);
    B_ = gsl_vector_get(s->x,1);
    tau_ = gsl_vector_get(s->x,2);
    
    //Tells you what the roots are when they are solved for
    if (test_status == GSL_SUCCESS) {
        std::cout << "Single particle potential parameters for V(n_0) at T = 0: A = " 
        << A_ << ", B = " << B_ << ", tau = " 
        << tau_ << " MeV" << std::endl;
    }
    //Tells you if the root solver fails to converge
    else {
        std::cerr << "Failed to converge:" << gsl_strerror(test_status) 
            << std::endl;
    }


    // FREE POINTERS!!!!!!
    gsl_multiroot_fsolver_free(s);
    gsl_vector_free(x);
}

//Calculates the single particle potential with the calculated parameters
double PotentialParameters::get_single_particle_potential(double density) {
    //Equation 15 in Interactions in Nuclear matter
    double potential = A_ * (density/nSat_in_MeV3) +
     B_ * std::pow(density/nSat_in_MeV3, tau_ - 1.0);

    std::cout << "Single particle potential for V(n_B = " << density 
    << ") = " << potential << " MeV" << std::endl;
    //Returns the potential for possible future use
    return potential; 
}
