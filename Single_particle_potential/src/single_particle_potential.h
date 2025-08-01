#ifndef single_particle_potential_h
#define single_particle_potential_h

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_vector.h>
#include <iomanip>
#include <memory>   // for std::unique_ptr
#include "./core_statmech_functions.h"
#include "./constants.h"



class PotentialParameters {
public:
  ///////////////////////////////////////////
  // Constructor
  PotentialParameters(double degeneracy,
		      double mass,
		      double saturation_density,
		      double binding_energy,
		      double incompress_at_satdense)
    : degeneracy_(degeneracy),
      mass_(mass),
      saturation_density_(saturation_density),
      binding_energy_(binding_energy),
      incompress_at_satdense_(incompress_at_satdense) {
    // (nothing to do here)
  }

  // Default destructor
  virtual ~PotentialParameters() {}


  
  ///////////////////////////////////////////
  // Member functions

  // Condition for the value of the binding energy
  // (equation 25 in "Interactions in nuclear matter" notes)
  double binding_energy_condition(double A, double B, double tau);

  // Condition for the location of the minimum in the energy per particle, i.e., value of
  // the saturation density (equation 30 in "Interactions in nuclear matter" notes)
  double energy_minimum_condition(double A, double B, double tau);

  // Condition for the value of the incompressibility at saturation density
  // (equation 33 in "Interactions in nuclear matter" notes)
  double incompressibility_condition(double A, double B, double tau);

  // Evaluates equations 25, 30, and 33 from "Interactions in nuclear matter" notes, which
  // are conditions that all three potential parameters A, B, tau must satisfy.
  // Note: this function needs to be declared static to satisfy the required GSL
  // signature of a function pointer passed to the gsl_multiroot_funtion F: F.f needs a
  // fuction pointer of the type int (*)(const gsl_vector*, void*, gsl_vector*), a plain
  // C-style function pointer without "this". However, a non-static class member function
  // implicitly would have the signature
  // int (PotentialParameters::*)(const gsl_vector*, void*, gsl_vector*).
  // Hence, declaring conditions() as static is necessary.
  static int conditions(const gsl_vector* x, void* p,
		 gsl_vector* fvec);

  // Uses GSL multiroot to find the parameters A, B, and tau of the single-particle
  // potential for given nuclear matter properties (saturation_density_, binding_energy_,
  // incompress_at_satdense_); default values of initial parameter guesses are provided
  void get_parameters(double tolerance, double A_guess = 200.0, double B_guess = 100.0,
		      double tau_guess = 3.0);

  // Calculates the single particle potential with the calculated parameters
  double get_single_particle_potential(double density);


  
  ///////////////////////////////////////////
  // Return functions for class members
  double degeneracy() { return degeneracy_; }
  double mass() { return mass_; }
  double saturation_density() { return saturation_density_; }
  double binding_energy() { return binding_energy_; }
  double incompress_at_satdense() { return incompress_at_satdense_; } 

  

private:
  ///////////////////////////////////////////
  // Nuclear matter properties:
  const double degeneracy_;
  const double mass_;
  const double saturation_density_;
  const double binding_energy_;
  const double incompress_at_satdense_;
  ///////////////////////////////////////////
  // Potential parameters, fit to reproduce nuclear matter properties:
  double A_{};
  double B_{};
  double tau_{};
};



#endif
