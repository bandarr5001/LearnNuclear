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
#include "./core_statmech_functions.h"
#include "./constants.h"

//struct for the parameters and nuclear properties for the calculation
struct potential_parameters {
    const double degeneracy_g_;
    const double mass_;
    const double density_;
    const double binding_energy_;
    const double incompress_at_satdense_;
};
extern potential_parameters pparams;

//struct for the results of the root solver
struct potential_results {
    double A;
    double B;
    double tau;
};
extern potential_results parameter_results;

struct FermiWrapper {
    potential_parameters params;

    //Constructor
    FermiWrapper(const potential_parameters& p);

    //Member Functions

    //Calculates the fermi momentum at density n_B
    double fermi_momentum() const;

    //Calculates the fermi energy at density n_B
    double fermi_energy() const;

    //Calculates the energy density of a spherical, ideal, 
    //noninteracting fermi gas at density n_B at T=0
    double energy_density_Fermigas() const;

    //Calculates the fermi momentum at saturation density in MeV^3
    double fermi_momentum_satdense() const;

    //Calculates the fermi energy at saturation density in MeV^3
    double fermi_energy_satdense() const;

    //Calculates the energy density of a spherical, ideal, 
    //noninteracting fermi gas at saturation density in MeV^3 at T=0
    double energy_density_Fermigas_satdense() const;

    //Static callbacks for GSL
    static double fermi_momentum_callback(void *p);
    static double energy_density_callback(void *p);

    static double fermi_momentum_callback_satdense(void *p);
    static double fermi_energy_callback_satdense(void *p);
};

//Contains equations 25, 30, and 33 from Interactions in nuclear matter
//Calculates the conditions that all three parameters must satisfy
int conditions(const gsl_vector* x, void* p,
					     gsl_vector* fvec);

//Uses GSL multiroot to calculate the parameters A, B, and tau for the 
//single particle potential 
potential_results get_parameters(potential_results parameter_results, 
    double tolerance, potential_parameters pparams);

//Calculates the single particle potential with the calculated parameters
double get_single_particle_potential(potential_results parameter_results,
     potential_parameters pparams);

#endif