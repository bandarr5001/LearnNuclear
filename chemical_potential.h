#ifndef chemical_potential_h
#define chemical_potential_h

#include "./fermi_gas.h"
#include "./bose_gas.h"
#include "./boltzmann_gas.h"
#include "./core_statmech_functions.h"
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_vector.h>
#include <iomanip>


//Fermi Gas *******************************************************************
double root_function_for_mu_fermi(Parameters& params, double chem_potent_mu);

double root_function_for_mu_fermi_gsl(double chem_potent, void *p);

double get_chem_potent_mu_fermi_bisection(double tolerance, double low, 
    double high, Parameters& params);

double get_chem_potent_mu_fermi_newton(double tolerance, double guess,
     Parameters& params);

double get_chem_potent_mu_fermi_newton_numerical(double tolerance, double guess,
     Parameters& params);

double gsl_chem_potent_fermi_bisection(double tolerance, Parameters& params);

void gsl_fdf_fermi_newton(double mu, void *p, double *f, double *df);

double gsl_chem_potent_fermi_newton(double tolerance, double guess, 
    Parameters& params);

double gsl_chem_potent_fermi_newton_numerical(double tolerance, double guess, 
    Parameters& params);

int root_function_for_mu_fermi_gsl_multiroot(const gsl_vector* x, void* p,
     gsl_vector* fvec);

double gsl_chem_potent_fermi_multiroot(double guess, double tolerance, 
    Parameters& params);

//Bose Gas ********************************************************************

//Calculates the root function to find the chemical potential of a Bose gas
//Uses the density integral
double root_function_for_mu_bose(double density, double chem_potent_mu,
    double mass, double degeneracy_g, double temperature);

//Calculates the chemical potential of a Bose gas using the bisection method
double get_chem_potent_mu_bose_bisection(double tolerance, double low, double high, double temperature, double density,
    double degeneracy_g, double mass);

//High Temperature Limit ******************************************************
double root_function_for_mu_boltzmann(double density, double chem_potent_mu,
    double mass, double degeneracy_g, double temperature);

//Energy needed to add a particle to the system
double get_chem_potent_mu_boltzmann_bisection(double tolerance, double low,
     double high, double temperature, double density, double degeneracy_g,
     double mass);

#endif