#include "./fermi_gas.h"
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

//T=0 *************************************************************************
//At T=0, the chemical potential is the fermi energy

double fermi_gas_press_temp0_pFermi(double degeneracy_g, double chem_potent_mu,
    double p_momentum_fermi, double mass){

    return (degeneracy_g / (16.0 * std::pow(M_PI,2.0)))
    * ((2.0/3.0)*chem_potent_mu*std::pow(p_momentum_fermi,3.0)
    - chem_potent_mu*p_momentum_fermi*std::pow(mass,2.0)
    + std::pow(mass,4.0)*std::log((chem_potent_mu + p_momentum_fermi)/mass));
}

double fermi_gas_press_temp0_chem(double chem_potent_mu, double degeneracy_g,
    double mass) { 
    
    //Calculates the ferm momentum using the chemical potential
    double p_momentum_fermi = std::sqrt(std::pow(chem_potent_mu,2.0)
     - std::pow(mass,2.0));

    //calculates the fermi pressure at T = 0

    return fermi_gas_press_temp0_pFermi(degeneracy_g, chem_potent_mu, 
        p_momentum_fermi, mass);
}

double fermi_gas_press_temp0_dense(double density_n, double degeneracy_g,
    double mass){

    //calculates the fermi momentum using the density, n
    double p_momentum_fermi = std::pow(6*std::pow(M_PI,2.0)
    *density_n/degeneracy_g,1.0/3.0);

    //calculates the chemical potential using the fermi momentum
    double chem_potent_mu = std::sqrt(std::pow(mass,2.0) 
    + std::pow(p_momentum_fermi,2.0));

    //Calculates fermi pressure at T = 0
    return fermi_gas_press_temp0_pFermi(degeneracy_g, chem_potent_mu, 
        p_momentum_fermi, mass);

}

//FINITE TEMPERATURE **********************************************************

double get_density_integrand_fermi(double t, double energy_epsilon, 
    double chem_potent_mu, Parameters& params) {

//density of a fermi gas and mapped to a bound of (0,1)
    return ((params.degeneracy_g_/(2.0*std::pow(M_PI,2.0))) * 
    (std::pow((1.0 - t)/t,2.0) * 
    get_FermiDirac_dist(params.temperature_, energy_epsilon, chem_potent_mu))) 
    / std::pow(t,2.0);
}

double get_P_spherical_integrand_fermi(double t, 
    double energy_epsilon, double chem_potent_mu, Parameters& params) {

//integrand for the fermi pressure of an ideal Fermi gas at the upper bound and
//mapped to a range of (0,1)
    return ((params.degeneracy_g_ / (6.0 * std::pow(M_PI,2.0))) * 
    (std::pow((1.0 - t)/t,4.0) / energy_epsilon) 
    * get_FermiDirac_dist(params.temperature_, energy_epsilon, chem_potent_mu))
     / std::pow(t,2.0);
}

double get_density_integral_fermi_gsl(double chem_potent, void *p) {
    auto *par = static_cast<Parameters*>(p);

    if (std::isnan(chem_potent)) {
        std::cerr << "ERROR: mu is NaN in get_density_integral_fermi_gsl"
         << std::endl;
        //return 0.0;
    }
    
    int intervals = 100000;
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
            get_energy_epsilon(par->mass_, (t_lower-1)/t_lower);
            //density of a fermi gas at lower bound
            density_lower = 
            get_density_integrand_fermi(t_lower,
                energy_epsilon_lower, chem_potent, *par);
            //NANs are caused by rounding errors due to the code thinking that
            //it is dividing by zero even though it is not. The values around
            //these go to 0 so these values are set to zero.
            if (std::isnan(density_lower)) {
                density_lower = 0.0;
            }
        }

        double energy_epsilon_upper = 0.0;
        double density_upper = 0.0;

        if(t_upper != 1.0) {
            //particle energy at the upper bound
            energy_epsilon_upper = 
            get_energy_epsilon(par->mass_, (t_upper-1)/t_upper);
            //density of a fermi gas at upper bound
            density_upper = 
            get_density_integrand_fermi(t_upper, energy_epsilon_upper, 
                chem_potent, *par);
            
            //NANs are caused by rounding errors due to the code thinking that
            //it is dividing by zero even though it is not. The values around
            //these go to 0 so these values are set to zero.
            if (std::isnan(density_upper)) {
                density_upper = 0.0;
            }
        }
            
        double area = 0.5 * h * (density_lower + density_upper);
        //std::cout << "Area: " << area << std::endl;
        sum += area;
    }
    //std::cout << "sum: " << sum << std::endl;

    if (std::isnan(sum) || !std::isfinite(sum)) {
        std::cerr << "ERROR: Integral sum is invalid (NaN or Inf), returning 0"
         << std::endl;
        return 0.0;  // or return a penalty value, or throw
    }
    return sum;
}

double get_density_integral_fermi(Parameters& params, double chem_potent_mu) {
    int intervals = 100000;
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
            get_energy_epsilon(params.mass_, (t_lower-1)/t_lower);
            //density of a fermi gas at lower bound
            density_lower = 
            get_density_integrand_fermi(t_lower,
                energy_epsilon_lower, chem_potent_mu, params);
            //NANs are caused by rounding errors due to the code thinking that
            //it is dividing by zero even though it is not. The values around
            //these go to 0 so these values are set to zero.
            if (std::isnan(density_lower)) {
                density_lower = 0.0;
            }
        }

        double energy_epsilon_upper = 0.0;
        double density_upper = 0.0;

        if(t_upper != 1.0) {
            //particle energy at the upper bound
            energy_epsilon_upper = 
            get_energy_epsilon(params.mass_, (t_upper-1)/t_upper);
            //density of a fermi gas at upper bound
            density_upper = 
            get_density_integrand_fermi(t_upper, energy_epsilon_upper, 
                chem_potent_mu, params);
            
            //NANs are caused by rounding errors due to the code thinking that
            //it is dividing by zero even though it is not. The values around
            //these go to 0 so these values are set to zero.
            if (std::isnan(density_upper)) {
                density_upper = 0.0;
            }
        }
            
        double area = 0.5 * h * (density_lower + density_upper);
        //std::cout << "Area: " << area << std::endl;
        sum += area;
    }
    //std::cout << "sum: " << sum << std::endl;
    return sum;
}

double get_density_derivative_integrand_fermi(double t, 
double energy_epsilon, double chem_potent_mu, Parameters& params) {

    double p_momentum = (1.0 - t)/t;
    double f_FD_dist = get_FermiDirac_dist(params.temperature_,energy_epsilon,chem_potent_mu);
    double exponential_components = f_FD_dist - f_FD_dist*f_FD_dist;

    //std::cout << "right: " << 1.0/(std::exp(term) + 1.0) << " wrong: " << 1.0/(std::exp(term) - 1.0) << " f_FD^2: " << f_FD_dist*f_FD_dist << std::endl;

    if(isinf(exponential_components)) {
        std::cerr << "Exponential is infinite" << std::endl;
    }
    else if(isnan(exponential_components)) {
        std::cerr << "Exponential is nan!" << std::endl;
    }

    double numerator = params.degeneracy_g_ * 
    exponential_components * p_momentum * p_momentum;

    double denominator = 2.0 * M_PI * params.temperature_ * t * t;

    double integrand = -numerator / denominator;

    if (isnan(integrand)) {
        std::cout << "The integrand is not a number! " << std::endl;
        std::cout << "Numerator: " << numerator << std::endl;
        std::cout << "exponential: " << exponential_components << std::endl;
        std::cout << "energy_epsilon: " << energy_epsilon << std::endl;
        std::cout << "exponent: " << (energy_epsilon - chem_potent_mu)
        /params.temperature_ << std::endl;
        std::cout << "temperature: " << params.temperature_ << std::endl;
        std::cout << "Denominator: " << denominator << std::endl;
        return 1e-8;
    }
    return integrand;
}

double get_density_derivative_integral_fermi_gsl_numerical(double chem_potent_mu, void *p) {

    double h = 1e-8;

    double fxplush = get_density_integral_fermi_gsl(chem_potent_mu + h, p);

    double fxminush = get_density_integral_fermi_gsl(chem_potent_mu - h, p);

    double fcalcprime = -(fxplush - fxminush)/(2*h);

    return fcalcprime;
}

double get_density_derivative_integral_fermi_gsl(double chem_potent, void *p) {

    auto *par = static_cast<Parameters*>(p);

    int intervals = 100000;
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
            get_energy_epsilon(par->mass_, (t_lower-1)/t_lower);
            //density of a fermi gas at lower bound
            density_lower = 
            get_density_derivative_integrand_fermi(t_lower,
                energy_epsilon_lower, chem_potent, *par);
            
            //NANs are caused by rounding errors due to the code thinking that
            //it is dividing by zero even though it is not. The values around
            //these go to 0 so these values are set to zero.
            if (std::isnan(density_lower)) {
                density_lower = 0.0;
            }
        }
        //std::cout << "density_lower (" << t_lower << "): " << density_lower
        // << std::endl;

        double energy_epsilon_upper = 0.0;
        double density_upper = 0.0;

        if(t_upper != 1.0) {
            //particle energy at the upper bound
            energy_epsilon_upper = 
            get_energy_epsilon(par->mass_, (t_upper-1)/t_upper);
            //density of a fermi gas at upper bound
            density_upper = 
            get_density_derivative_integrand_fermi(t_upper,
                 energy_epsilon_upper, chem_potent, *par);

            //NANs are caused by rounding errors due to the code thinking that
            //it is dividing by zero even though it is not. The values around
            //these go to 0 so these values are set to zero.
            if (std::isnan(density_upper)) {
                density_upper = 0.0;
            }
        }
        //std::cout << "density_upper (" << t_upper << "): " << density_upper 
        //<< std::endl;
            
        double area = 0.5 * h * (density_lower + density_upper);
        //std::cout << "Area: " << area << std::endl;
        sum += area;
    }
    //std::cout << "sum derivative: " << sum << std::endl;
    return sum;
}

//Taking the numerical derivative
double get_density_derivative_integral_fermi_numerical(double chem_potent_mu, Parameters& params) {
    double h = 1e-8;

    double fxplush = get_density_integral_fermi(params,chem_potent_mu + h);

    double fxminush = get_density_integral_fermi(params,chem_potent_mu - h);

    double fcalcprime = -(fxplush - fxminush)/(2*h);

    return fcalcprime;
}

double get_density_derivative_integral_fermi(double chem_potent_mu, 
    Parameters& params) {

    int intervals = 100000;
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
            get_energy_epsilon(params.mass_, (t_lower-1)/t_lower);
            //density of a fermi gas at lower bound
            density_lower = 
            get_density_derivative_integrand_fermi(t_lower,
                energy_epsilon_lower, chem_potent_mu, params);
            
            //NANs are caused by rounding errors due to the code thinking that
            //it is dividing by zero even though it is not. The values around
            //these go to 0 so these values are set to zero.
            if (std::isnan(density_lower)) {
                density_lower = 0.0;
            }
        }
        //std::cout << "density_lower (" << t_lower << "): " << density_lower
        // << std::endl;

        double energy_epsilon_upper = 0.0;
        double density_upper = 0.0;

        if(t_upper != 1.0) {
            //particle energy at the upper bound
            energy_epsilon_upper = 
            get_energy_epsilon(params.mass_, (t_upper-1)/t_upper);
            //density of a fermi gas at upper bound
            density_upper = 
            get_density_derivative_integrand_fermi(t_upper,
                 energy_epsilon_upper, chem_potent_mu, params);

            //NANs are caused by rounding errors due to the code thinking that
            //it is dividing by zero even though it is not. The values around
            //these go to 0 so these values are set to zero.
            if (std::isnan(density_upper)) {
                density_upper = 0.0;
            }
        }
        //std::cout << "density_upper (" << t_upper << "): " << density_upper 
        //<< std::endl;
            
        double area = 0.5 * h * (density_lower + density_upper);
        //std::cout << "Area: " << area << std::endl;
        sum += area;
    }
    //std::cout << "sum derivative: " << sum << std::endl;
    return sum;
}

//fermi pressure of an ideal Fermi gas at the upper bound and
//mapped to a range of (0,1)
double fermi_gas_pressure_integral_spherical(Parameters& params,
     double chem_potent_mu) {
    
    int intervals = 10000;
    //Using a trapezoidal Riemann sum
    const double h = 1.0/(1.0*intervals); //h is the height of each trapezoid

    double sum = 0.0;

    for(int i = 0; i < intervals; i ++) {

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
            get_energy_epsilon(params.mass_, (t_lower-1)/t_lower);

//fermi pressure of a spherical, isotropic fermi gas at the lower bound
// and mapped to a range of (0,1)
            P_spherical_lower =
            get_P_spherical_integrand_fermi(t_lower,
                energy_epsilon_lower, chem_potent_mu, params);
        }
        
        double energy_epsilon_upper = 0.0;
        double P_spherical_upper = 0.0;
        
/*At t = 1, the integrand is 0 due to the 1 - t in the numerator
 (see the function get_P_spherical_integrand)*/

        if(t_lower != 1.0) {
//particle energy at the upper bound and mapped to a range of (0,1)
            energy_epsilon_upper =
                get_energy_epsilon(params.mass_, (t_upper-1)/t_upper);
        
//fermi pressure of a spherical, isotropic fermi gas at the upper bound
// and mapped to a range of (0,1)
            P_spherical_upper = 
                get_P_spherical_integrand_fermi(t_upper,
                energy_epsilon_upper, chem_potent_mu, params);
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