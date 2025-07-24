#include "./chemical_potential.h"
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


double root_function_for_mu_fermi(Parameters& params, double chem_potent_mu) {
    return params.density_ - 
    get_density_integral_fermi(params, chem_potent_mu);
}

double root_function_for_mu_fermi_gsl(double chem_potent, void *p) {
    auto *par = static_cast<Parameters*>(p);

    //std::cout << "density in root_function_for_mu_fermi_gsl: "
    // << par->density_ << std::endl;

    double residual = par->density_ - 
    get_density_integral_fermi(*par, chem_potent);

    //std::cout << "residual in root_function_for_mu_fermi_gsl: "
    // << residual << std::endl;

    return residual;
}

//Energy needed to add a particle to the system
double get_chem_potent_mu_fermi_bisection(double tolerance, double low, 
    double high, Parameters& params) {
    //Using the bisection method
    double f_low = 
    root_function_for_mu_fermi(params, low);
    double f_high = 
    root_function_for_mu_fermi(params, high);

    // Check if there's a root in the interval
    if (f_low * f_high > 0) {
        std::cerr << "Error: f_low and f_high have the same sign" << std::endl;
        std::cout << "f_low: " << f_low << " f_high: " << f_high << std::endl;
        std::cout << "temperature: " << params.temperature_ << " density: " 
        << params.density_ 
        << " degeneracy_g: " << params.degeneracy_g_ << " mass: " 
        << params.mass_ << std::endl;
        return NAN;
    }

    double mid = 0.5 * (low + high);


    while(std::abs(high - low)/2 > tolerance) {
        mid = 0.5 * (low + high);
        double f_mid = root_function_for_mu_fermi(params, mid);

        if (std::abs(f_mid) < tolerance) return mid;

        if (f_low * f_mid < 0.0) {
            high = mid;
            f_high = f_mid;
        }
        
        else {
            low = mid;
            f_low = f_mid;
        }
    }

    double root = 0.5 * (low + high);

    std::cout << "Fermi gas chemical potential bisection: mu = " << 
    std::setprecision (15) << root << " MeV" << std::endl;

    return root;
}

double get_chem_potent_mu_fermi_newton_numerical(double tolerance, double guess,
     Parameters& params) {

    double root_condition = 1.0;

    double chem_potent_mu = guess;

    double delta_mu = 0.0;

    double density_deriv_integral = 0.0;

    int i = 0, iMAX = 100;

    while (std::abs(root_condition) > tolerance && i++ < iMAX) {
        root_condition = root_function_for_mu_fermi(params, chem_potent_mu);
        //std::cout << "root_condition: " << root_condition << std::endl;

        density_deriv_integral = 
        get_density_derivative_integral_fermi_numerical(chem_potent_mu, params);

        //std::cout << "density_integral: " << density_integral << 
        //" density_deriv_integral: " << density_deriv_integral << std::endl;

        delta_mu = root_condition/density_deriv_integral;

        if(std::isnan(delta_mu)) {
            std::cerr << "delta_mu is nan in get_chem_potent_mu_fermi_newton"
             << std::endl;
            return chem_potent_mu;
        }

        //std::cout << "delta_mu: " << delta_mu << std::endl;

        chem_potent_mu -= delta_mu;
        //std::cout << "chem_potent_mu: " << chem_potent_mu << std::endl;

    }
    std::cout << "Fermi gas chemical potential numerical Newton-Raphson: mu = " << 
    std::setprecision (15) << chem_potent_mu
     << " MeV" << std::endl;
    return chem_potent_mu;
}

double get_chem_potent_mu_fermi_newton(double tolerance, double guess,
     Parameters& params) {

    double root_condition = 1.0;

    double chem_potent_mu = guess;

    double delta_mu = 0.0;

    double density_deriv_integral = 0.0;

    int i = 0, iMAX = 100;

    while (std::abs(root_condition) > tolerance && i++ < iMAX) {
        root_condition = root_function_for_mu_fermi(params, chem_potent_mu);
        //std::cout << "root_condition: " << root_condition << std::endl;

        density_deriv_integral = 
        get_density_derivative_integral_fermi(chem_potent_mu, params);

        //std::cout << "density_integral: " << density_integral << 
        //" density_deriv_integral: " << density_deriv_integral << std::endl;

        delta_mu = root_condition/density_deriv_integral;

        if(std::isnan(delta_mu)) {
            std::cerr << "delta_mu is nan in get_chem_potent_mu_fermi_newton"
             << std::endl;
            return chem_potent_mu;
        }

        //std::cout << "delta_mu: " << delta_mu << std::endl;

        chem_potent_mu -= delta_mu;
        //std::cout << "chem_potent_mu: " << chem_potent_mu << std::endl;

    }
    std::cout << "Fermi gas chemical potential Newton-Raphson: mu = " << chem_potent_mu
    << std::setprecision (15) << " MeV" << std::endl;
    return chem_potent_mu;
}

double gsl_chem_potent_fermi_bisection(double tolerance, Parameters& params) {

    gsl_set_error_handler_off();

    gsl_function F;
    F.function = &root_function_for_mu_fermi_gsl;
    F.params = &params;

    const gsl_root_fsolver_type *T = gsl_root_fsolver_bisection;
    gsl_root_fsolver *solver = gsl_root_fsolver_alloc(T);

    double lower = 500;
    double upper = 1500;

    gsl_root_fsolver_set(solver, &F, lower, upper);

    //std::cout << "HERE" << std::endl;

    int test_status;
    int i = 0;
    int iMAX = 100;
    double root;
    do {
        i++;
        gsl_root_fsolver_iterate(solver);             // Do one step
        root = gsl_root_fsolver_root(solver);         // Get current guess
        lower = gsl_root_fsolver_x_lower(solver);   // Update bracket
        upper = gsl_root_fsolver_x_upper(solver);

        // Stop when close enough
        test_status = gsl_root_test_interval(lower, upper, 0.0, tolerance);

        //std::cout << "Iter " << i << ": x = " << root << "\n";
    } while (test_status == GSL_CONTINUE && i < iMAX);

    if (test_status == GSL_SUCCESS) {
        std::cout << "Fermi gas chemical potential bisection method GSL: mu = "
        << std::setprecision (15) << root << " MeV" << std::endl;
    }
    else {
        std::cerr << "Failed to converge" << std::endl;
    }

    return root;
}

void gsl_fdf_fermi_newton(double mu, void *p, double *f, double *df) {
    auto *par = static_cast<Parameters*>(p);

    *f  = root_function_for_mu_fermi_gsl(mu,par);
    *df = get_density_derivative_integral_fermi_gsl_numerical(mu,par);

    //std::cout << "f(mu) = " << *f << " df = " << *df << std::endl;
}

double gsl_chem_potent_fermi_newton(double tolerance, double guess, 
    Parameters& params) {
    //gsl_set_error_handler_off();
    gsl_function_fdf F;

    F.f = &root_function_for_mu_fermi_gsl;
    F.df = &get_density_derivative_integral_fermi_gsl;
    F.fdf = &gsl_fdf_fermi_newton;
    F.params = &params;


    const gsl_root_fdfsolver_type* T = gsl_root_fdfsolver_newton;
    gsl_root_fdfsolver* s = gsl_root_fdfsolver_alloc(T);

    gsl_root_fdfsolver_set(s,&F,guess);

    int test_status;
    int i = 0;
    int iMAX = 100;
    double root = guess;
    do {   
        i++;
        double root_prev = root; 
        test_status = gsl_root_fdfsolver_iterate(s);
        root = gsl_root_fdfsolver_root(s);
        test_status =  gsl_root_test_delta(root, root_prev, 0, tolerance);

        /*std::cout << "iter " << i 
          << ": mu = " << root 
          << ", f(mu) = " << root_function_for_mu_fermi_gsl(root, &params)
          << ", df(mu) = " << 
          get_density_derivative_integral_fermi_gsl(root, &params)
          << std::endl;*/
    } while(test_status == GSL_CONTINUE && i < iMAX);
    
    if (test_status == GSL_SUCCESS) {
        std::cout << "Fermi gas chemical potential Newton-Raphson GSL: mu = "
         << root << " MeV" << std::endl;
    }
    else {
        std::cerr << "Failed to converge:" << gsl_strerror(test_status) 
        << std::endl;
    }
    gsl_root_fdfsolver_free(s);

    return root;
}

double gsl_chem_potent_fermi_newton_numerical(double tolerance, double guess, 
    Parameters& params) {
    //gsl_set_error_handler_off();
    gsl_function_fdf F;

    F.f = &root_function_for_mu_fermi_gsl;
    F.df = &get_density_derivative_integral_fermi_gsl_numerical;
    F.fdf = &gsl_fdf_fermi_newton;
    F.params = &params;


    const gsl_root_fdfsolver_type* T = gsl_root_fdfsolver_newton;
    gsl_root_fdfsolver* s = gsl_root_fdfsolver_alloc(T);

    gsl_root_fdfsolver_set(s,&F,guess);

    int test_status;
    int i = 0;
    int iMAX = 100;
    double root = guess;
    do {   
        i++;
        double root_prev = root; 
        test_status = gsl_root_fdfsolver_iterate(s);
        root = gsl_root_fdfsolver_root(s);
        test_status =  gsl_root_test_delta(root, root_prev, 0, tolerance);

        /*std::cout << "iter " << i 
          << ": mu = " << root 
          << ", f(mu) = " << root_function_for_mu_fermi_gsl(root, &params)
          << ", df(mu) = " << 
          get_density_derivative_integral_fermi_gsl(root, &params)
          << std::endl;*/
    } while(test_status == GSL_CONTINUE && i < iMAX);
    
    if (test_status == GSL_SUCCESS) {
        std::cout << "Fermi gas chemical potential numerical Newton-Raphson GSL: mu = "
         << root << " MeV" << std::endl;
    }
    else {
        std::cerr << "Failed to converge:" << gsl_strerror(test_status) 
        << std::endl;
    }
    gsl_root_fdfsolver_free(s);

    return root;
}

int root_function_for_mu_fermi_gsl_multiroot(const gsl_vector* x, void* p,
     gsl_vector* fvec) {
    auto* par = static_cast<Parameters*>(p);
    
    double chem_potent_mu = gsl_vector_get(x, 0);

    if (std::isnan(chem_potent_mu)) {
        std::cerr << "FATAL: μ is NaN inside root_function: " << chem_potent_mu
         << std::endl;
        gsl_vector_set(fvec, 0, 1e10); // Large residual penalty
        return GSL_SUCCESS;
    }

    double rho = get_density_integral_fermi_gsl(chem_potent_mu, par);

    if (std::isnan(rho) || !std::isfinite(rho)) {
        std::cerr << "FATAL: ρ(μ) returned NaN or Inf!" << std::endl;
        gsl_vector_set(fvec, 0, 1e10);
        return GSL_SUCCESS;
    }

    double residual = par->density_ - rho;

    /*std::cout << "root function call: μ = " << chem_potent_mu
              << ", ρ(μ) = " << rho
              << ", residual = " << residual << std::endl;*/
    gsl_vector_set(fvec, 0, residual);
    return GSL_SUCCESS;
}

double gsl_chem_potent_fermi_multiroot(double guess, double tolerance, 
    Parameters& params) {

    const size_t dim = 1;

    gsl_multiroot_function F;

    F.f = &root_function_for_mu_fermi_gsl_multiroot;
    F.n = dim;
    F.params = &params;
    
    gsl_vector* x = gsl_vector_alloc(dim);
    gsl_vector_set(x, 0, guess);

    const gsl_multiroot_fsolver_type* T = gsl_multiroot_fsolver_hybrids;
    gsl_multiroot_fsolver* s = gsl_multiroot_fsolver_alloc(T, dim);
    gsl_multiroot_fsolver_set(s, &F, x);

    int test_status;

    size_t i = 0, iMAX = 100;

    //std::cout << "iter μ          f(μ)            |Δμ|" << std::endl;

    do {
        i++;
        //double root_prev = gsl_vector_get(s->x,0);
        //std::cout << "previous root:" << root_prev << std::endl;

        test_status = gsl_multiroot_fsolver_iterate(s);
        
        //double root = gsl_vector_get(s->x,0);

        //std::cout << "current root:" << root << std::endl;
        //double residual = gsl_vector_get(s->f,0);

        /*std::cout << i << "  " << std::setw(12) << root
                            << "  " << std::setw(12) << residual
                            << "  " << std::setw(12) 
                            << std::fabs(root-root_prev)
                            << std::endl;*/

        test_status = gsl_multiroot_test_residual(s->f, tolerance);


    } while (test_status == GSL_CONTINUE && i < iMAX);

    double result = gsl_vector_get(s->x,0);


    
    if (test_status == GSL_SUCCESS) {
        std::cout << "Fermi gas chemical potential multiroot GSL: mu = " 
        << result << " MeV" << std::endl;
    }

    else {
        std::cerr << "Failed to converge:" << gsl_strerror(test_status) 
        << std::endl;
    }

    //gsl_multiroot_fsolver_free(s);
    //gsl_vector_free(x);

    return result;
}

//Bose Gas ********************************************************************

int root_function_for_mu_bose_multiroot(const gsl_vector* x, void* p,
     gsl_vector* fvec) {
    auto* par = static_cast<Parameters*>(p);
    
    double chem_potent_mu = gsl_vector_get(x, 0);

    if (std::isnan(chem_potent_mu)) {
        std::cerr << "FATAL: μ is NaN inside root_function: " << chem_potent_mu
         << std::endl;
        gsl_vector_set(fvec, 0, 1e10); // Large residual penalty
        return GSL_SUCCESS;
    }

    double rho = get_density_integral_bose(chem_potent_mu, par);

    if (std::isnan(rho) || !std::isfinite(rho)) {
        std::cerr << "FATAL: ρ(μ) returned NaN or Inf!" << std::endl;
        gsl_vector_set(fvec, 0, 1e10);
        return GSL_SUCCESS;
    }

    double residual = par->density_ - rho;

    /*std::cout << "root function call: μ = " << chem_potent_mu
              << ", ρ(μ) = " << rho
              << ", residual = " << residual << std::endl;*/
    gsl_vector_set(fvec, 0, residual);
    return GSL_SUCCESS;
}

//Energy needed to add a particle to the system
double get_chem_potent_mu_bose_multiroot(double guess, double tolerance, 
    Parameters& params) {
    const size_t dim = 1;

    gsl_multiroot_function F;

    F.f = &root_function_for_mu_bose_multiroot;
    F.n = dim;
    F.params = &params;
    
    gsl_vector* x = gsl_vector_alloc(dim);
    gsl_vector_set(x, 0, guess);

    const gsl_multiroot_fsolver_type* T = gsl_multiroot_fsolver_hybrids;
    gsl_multiroot_fsolver* s = gsl_multiroot_fsolver_alloc(T, dim);
    gsl_multiroot_fsolver_set(s, &F, x);

    int test_status;

    size_t i = 0, iMAX = 100;

    //std::cout << "iter μ          f(μ)            |Δμ|" << std::endl;

    do {
        i++;
        //double root_prev = gsl_vector_get(s->x,0);
        //std::cout << "previous root:" << root_prev << std::endl;

        test_status = gsl_multiroot_fsolver_iterate(s);
        
        //double root = gsl_vector_get(s->x,0);

        //std::cout << "current root:" << root << std::endl;
        //double residual = gsl_vector_get(s->f,0);

        /*std::cout << i << "  " << std::setw(12) << root
                            << "  " << std::setw(12) << residual
                            << "  " << std::setw(12) 
                            << std::fabs(root-root_prev)
                            << std::endl;*/

        test_status = gsl_multiroot_test_residual(s->f, tolerance);


    } while (test_status == GSL_CONTINUE && i < iMAX);

    double result = gsl_vector_get(s->x,0);


    
    if (test_status == GSL_SUCCESS) {
        std::cout << "Bose gas chemical potential multiroot GSL: mu = " 
        << result << " MeV" << std::endl;
    }

    else {
        std::cerr << "Failed to converge:" << gsl_strerror(test_status) 
        << std::endl;
    }

    //gsl_multiroot_fsolver_free(s);
    //gsl_vector_free(x);

    return result;
}

//High Temperature Limit ******************************************************

int root_function_for_mu_boltzmann_multiroot(const gsl_vector* x, void* p,
     gsl_vector* fvec) {
    auto* par = static_cast<Parameters*>(p);
    
    double chem_potent_mu = gsl_vector_get(x, 0);

    if (std::isnan(chem_potent_mu)) {
        std::cerr << "FATAL: μ is NaN inside root_function: " << chem_potent_mu
         << std::endl;
        gsl_vector_set(fvec, 0, 1e10); // Large residual penalty
        return GSL_SUCCESS;
    }

    double rho = get_density_integral_boltzmann(chem_potent_mu, par);

    if (std::isnan(rho) || !std::isfinite(rho)) {
        std::cerr << "FATAL: ρ(μ) returned NaN or Inf!" << std::endl;
        gsl_vector_set(fvec, 0, 1e10);
        return GSL_SUCCESS;
    }

    double residual = par->density_ - rho;

    /*std::cout << "root function call: μ = " << chem_potent_mu
              << ", ρ(μ) = " << rho
              << ", residual = " << residual << std::endl;*/
    gsl_vector_set(fvec, 0, residual);
    return GSL_SUCCESS;
}

//Energy needed to add a particle to the system
double get_chem_potent_mu_boltzmann_multiroot(double guess, double tolerance, 
    Parameters& params) {

    const size_t dim = 1;

    gsl_multiroot_function F;

    F.f = &root_function_for_mu_boltzmann_multiroot;
    F.n = dim;
    F.params = &params;
    
    gsl_vector* x = gsl_vector_alloc(dim);
    gsl_vector_set(x, 0, guess);

    const gsl_multiroot_fsolver_type* T = gsl_multiroot_fsolver_hybrids;
    gsl_multiroot_fsolver* s = gsl_multiroot_fsolver_alloc(T, dim);
    gsl_multiroot_fsolver_set(s, &F, x);

    int test_status;

    size_t i = 0, iMAX = 100;

    //std::cout << "iter μ          f(μ)            |Δμ|" << std::endl;

    do {
        i++;
        //double root_prev = gsl_vector_get(s->x,0);
        //std::cout << "previous root:" << root_prev << std::endl;

        test_status = gsl_multiroot_fsolver_iterate(s);
        
        //double root = gsl_vector_get(s->x,0);

        //std::cout << "current root:" << root << std::endl;
        //double residual = gsl_vector_get(s->f,0);

        /*std::cout << i << "  " << std::setw(12) << root
                            << "  " << std::setw(12) << residual
                            << "  " << std::setw(12) 
                            << std::fabs(root-root_prev)
                            << std::endl;*/

        test_status = gsl_multiroot_test_residual(s->f, tolerance);


    } while (test_status == GSL_CONTINUE && i < iMAX);

    double result = gsl_vector_get(s->x,0);


    
    if (test_status == GSL_SUCCESS) {
        std::cout << "Boltzmann gas chemical potential multiroot GSL: mu = " 
        << result << " MeV" << std::endl;
    }

    else {
        std::cerr << "Failed to converge:" << gsl_strerror(test_status) 
        << std::endl;
    }

    //gsl_multiroot_fsolver_free(s);
    //gsl_vector_free(x);

    return result;
}