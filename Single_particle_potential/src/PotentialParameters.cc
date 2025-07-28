#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_vector.h>
#include <iomanip>
#include "./single_particle_potential.h"

int main() {
    std::cout << "\n\n*********************************************************************"
	    << "\n* This code evaluates the potential parameters A, B, and tau for single particle potentials *"
	    << "\n*********************************************************************\n"
	    << std::endl;

    try{

      potential_parameters pparams = {
        .degeneracy_g_ = 4.0,     //spin 1/2 for a fermion
        .mass_ = mNucleon,         // mean proton + neutron mass in MeV
        .density_ = nSat_in_MeV3, //density in MeV^3
        .binding_energy_ = -16.3, //in MeV
        .incompress_at_satdense_ = 240.0 //in MeV
      };

      potential_results parameter_results = {
        .A = -218.699,
        .B = 166.215,
        .tau = 2.33376
      };

    double p_momentum_Fermi = std::pow(6.0*std::pow(M_PI,2.0)*pparams.density_/pparams.degeneracy_g_,1.0/3.0); 
    double energy_epsilon_Fermi = std::sqrt(pparams.mass_*pparams.mass_ + p_momentum_Fermi*p_momentum_Fermi);

    double term1_1 = 2.0*energy_epsilon_Fermi*energy_epsilon_Fermi*energy_epsilon_Fermi*p_momentum_Fermi;

    double term2_1 = energy_epsilon_Fermi*p_momentum_Fermi*pparams.mass_*pparams.mass_;

    double term3_1 = pparams.mass_*pparams.mass_*pparams.mass_*pparams.mass_*std::log((energy_epsilon_Fermi + p_momentum_Fermi)/pparams.mass_);

    double energy_density_FG = (pparams.degeneracy_g_ / (16.0 * M_PI*M_PI))
    * (term1_1 - term2_1 - term3_1);

    double p_momentum_Fermi_satdense = std::pow(6.0*std::pow(M_PI,2.0)
    *nSat_in_MeV3/pparams.degeneracy_g_,1.0/3.0);

    double energy_epsilon_Fermi_satdense = std::sqrt(pparams.mass_*pparams.mass_ + 
        p_momentum_Fermi_satdense*p_momentum_Fermi_satdense);

    double term1_2 = 2.0*energy_epsilon_Fermi_satdense*energy_epsilon_Fermi_satdense*energy_epsilon_Fermi_satdense*p_momentum_Fermi_satdense;

    double term2_2 = energy_epsilon_Fermi_satdense*p_momentum_Fermi_satdense*pparams.mass_*pparams.mass_;

    double term3_2 = pparams.mass_*pparams.mass_*pparams.mass_*pparams.mass_*std::log((energy_epsilon_Fermi_satdense + p_momentum_Fermi_satdense)/pparams.mass_);

    double energy_density_FG_satdense = (pparams.degeneracy_g_ / (16.0 * M_PI*M_PI))
    * (term1_2 - term2_2 - term3_2);

    double condition1 = (energy_density_FG/nSat_in_MeV3) 
        + ((parameter_results.A/2.0) + (parameter_results.B/parameter_results.tau))*nSat_in_MeV3 - pparams.mass_ - pparams.binding_energy_;

    double condition2 = energy_epsilon_Fermi_satdense 
        - (energy_density_FG_satdense/nSat_in_MeV3) + (parameter_results.A/2) + parameter_results.B*((parameter_results.tau - 1)/parameter_results.tau);

      double term1_3 = (pparams.degeneracy_g_*energy_epsilon_Fermi_satdense
        *p_momentum_Fermi_satdense*p_momentum_Fermi_satdense*p_momentum_Fermi_satdense)/(9*nSat_in_MeV3
        *M_PI*M_PI);

    double term2_3 = parameter_results.A;//*(par->density_/nSat_in_MeV3);

    double term3_3 = parameter_results.B*(parameter_results.tau-1.0);//*std::pow((par->density_/nSat_in_MeV3),tau - 1.0);

    double K = 9.0 * (term1_3 + term2_3 + term3_3);

    double condition3 = K - pparams.incompress_at_satdense_;

    std::cout << "TEST OF CONDITIONS" << std::endl;
    std::cout << "condition1 = " << condition1 << std::endl;
    std::cout << "condition2 = " << condition2 << std::endl;
    std::cout << "condition3 = " << condition3 << std::endl;

      double tolerance = 1e-6;

      get_parameters(parameter_results, tolerance, pparams);

      std::cout << "\n\n*********************************************************************"
	    << "\n* Finished evaluating the potential parameters A, B, and tau for single particle potentials *"
	    << "\n*********************************************************************\n"
	    << std::endl;
    } catch (const std::runtime_error& e) {
    std::cerr << "Caught exception: " << e.what() << std::endl;
  }
}