#include "./Fermi_gas_functions_Tzero.h"

#include <iostream>
#include <cmath>



// Kinetic energy of a relativistic particle
// (not a Fermi gas only quantity but it'd be a waste to have a separate file for it)
double E_kinetic(double mass, double momentum) {
  return std::sqrt(mass * mass + momentum * momentum);
}



// Fermi momentum at baryon density n_B (T=0)
double fermi_momentum(double degeneracy, double density) {
  return std::pow((6.0 * std::pow(M_PI, 2.0) * density / degeneracy), 1.0/3.0);
}



// Fermi energy at baryon density n_B (T=0)
double fermi_energy(double degeneracy, double mass, double density) {
  const double p_F = fermi_momentum(degeneracy, density);
  return E_kinetic(mass, p_F);
}



// Energy density of an ideal (non-interacting) Fermi gas at baryon density n_B (T=0)
double Fermi_gas_energy_density(double degeneracy, double mass, double density) {
  // Calculate the Fermi momentum and energy:
  const double p_F = fermi_momentum(degeneracy, density);
  const double e_F = fermi_energy(degeneracy, mass, density);
  
  const double term1 = 2.0 * std::pow(e_F, 3.0) * p_F;
  const double term2 = e_F * p_F * std::pow(mass, 2.0);
  const double term3 = std::pow(mass, 4.0) * std::log((e_F + p_F)/mass);

  return (degeneracy / (16.0 * M_PI* M_PI)) * (term1 - term2 - term3);
}



// Derivative of the Fermi pressure w.r.t. baryon density nB (T=0)
double d_PFG_d_nB(double degeneracy, double mass, double density) {
  // Calculate the Fermi momentum and energy:
  const double p_F = fermi_momentum(degeneracy, density);
  const double e_F = fermi_energy(degeneracy, mass, density);

  return std::pow(p_F, 2.0) / (3.0 * e_F);
}


