#include "./bose_gas.h"
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

//Bose-Einstein Distribution **************************************************

double get_density_integrand_bose(double degeneracy_g, double t, 
    double temperature, double energy_epsilon, double chem_potent_mu) {

//density of a fermi gas and mapped to a bound of (0,1)
    return ((degeneracy_g/(2.0*std::pow(M_PI,2.0))) * 
    (std::pow((1 - t)/t,2.0) * 
    get_BoseEinstein_dist(temperature, energy_epsilon, chem_potent_mu))) 
    / std::pow(t,2.0);
}

double get_density_integral_bose(double chem_potent_mu, double mass,
    double degeneracy_g, double temperature) {
    int intervals = 10000;
    double h = 1.0/(1.0*intervals);
    double sum = 0.0;
    //Loop over the intervals
    for(int i = 0; i < intervals; i++) {
    //edges of current interval
        double t_lower = i * h;
        double t_upper = (i*1.0 + 1.0) * h;
            
        double energy_epsilon_lower = 0.0;
        double density_lower = 0.0;

//At t = 0, the 1/t^2 term diverges to infinity, but the argument of the
//exponential diverges as well and more quickly than the two previous terms.
//Since the exponential is in the denominator, the 
//integrand goes to 0 as t goes to 0.

        if(t_lower != 0.0) {
            //particle energy at the lower bound
            energy_epsilon_lower = 
            get_energy_epsilon(mass, (t_lower-1)/t_lower);
            //density of a fermi gas at lower bound
            density_lower = 
            get_density_integrand_bose(degeneracy_g, t_lower, temperature,
                energy_epsilon_lower, chem_potent_mu);
        }

        double energy_epsilon_upper = 0.0;
        double density_upper = 0.0;

        if(t_upper != 1.0) {
            //particle energy at the upper bound
            energy_epsilon_upper = 
            get_energy_epsilon(mass, (t_upper-1)/t_upper);
            //density of a fermi gas at upper bound
            density_upper = 
            get_density_integrand_bose(degeneracy_g, t_upper, temperature,
             energy_epsilon_upper, chem_potent_mu);
        }
            
        double area = 0.5 * h * (density_lower + density_upper);
        //std::cout << "Area: " << area << std::endl;
        sum += area;
    }
    return sum;
}

double get_P_spherical_integrand_bose(double degeneracy_g, double t, 
    double energy_epsilon, double temperature, double chem_potent_mu) {

//integrand for the fermi pressure of an ideal Fermi gas at the upper bound and
//mapped to a range of (0,1)
    return ((degeneracy_g / (6.0 * std::pow(M_PI,2.0))) * 
    (std::pow((1 - t)/t,4.0) / energy_epsilon) 
    * get_BoseEinstein_dist(temperature, energy_epsilon, chem_potent_mu))
     / std::pow(t,2.0);
}

//fermi pressure of an ideal Bose gas at the upper bound and
//mapped to a range of (0,1)
double bose_gas_pressure_integral_spherical(double temperature,
    double chem_potent_mu, double degeneracy_g, double mass) {

    int intervals = 10000;

    //Using a trapezoidal Riemann sum
    const double h = 1.0/(1.0*intervals); //h is the height of each trapezoid

    double sum = 0.0;

    for(int i = 0; i < intervals; i++) {

        double t_lower = i * h; 
        double t_upper = (i*1.0 + 1.0) * h; 

        double energy_epsilon_lower = 0.0;
        double P_spherical_lower = 0.0;

/*At t = 0, the 1/t^2 term diverges to infinity as well as the [(1-t)/t]^4
term, but the argument of the exponential diverges as well and more quickly
than the two previous terms. Since the exponential is in the denominator, the 
integrand goes to 0 as t goes to 0.*/

        if(t_lower != 0.0) {
//particle energy at the lower bound and mapped to a range of (0,1)
            energy_epsilon_lower = 
            get_energy_epsilon(mass, (t_lower-1)/t_lower);

//fermi pressure of a spherical, isotropic fermi gas at the lower bound
// and mapped to a range of (0,1)
            P_spherical_lower =
            get_P_spherical_integrand_bose(degeneracy_g, t_lower,
                energy_epsilon_lower, temperature, chem_potent_mu);
        }
        
        double energy_epsilon_upper = 0.0;
        double P_spherical_upper = 0.0;
        
/*At t = 1, the integrand is 0 due to the 1 - t in the numerator
 (see the function get_P_spherical_integrand)*/

        if(t_lower != 1.0) {
//particle energy at the upper bound and mapped to a range of (0,1)
            energy_epsilon_upper =
                get_energy_epsilon(mass, (t_upper-1)/t_upper);
        
//fermi pressure of a spherical, isotropic fermi gas at the upper bound
// and mapped to a range of (0,1)
            P_spherical_upper = 
                get_P_spherical_integrand_bose(degeneracy_g, t_upper,
                energy_epsilon_upper, temperature, chem_potent_mu);
        }
    
        //Area of the trapezoid
        const double area = 0.5 * h * (P_spherical_lower + P_spherical_upper);

        sum += area; //sums all of the areas step by step

        //std::cout << "area: " << area << std::endl;
        //std::cout << "t_upper" << t_upper << std::endl;
    }

    //hbar*c = 197.3 MeV fm
    //Must divide sum by (hbar*c)^3
    double Pressure = sum/std::pow(197.3,3);

    return Pressure; //Returns the area under the curve in MeV/fm^3
}