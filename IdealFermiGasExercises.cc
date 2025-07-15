#include <iostream>
#include <cmath>
#include <cstdlib>
#include "./fermi_gas.h"
#include "./bose_gas.h"
#include "./boltzmann_gas.h"
#include "./chemical_potential.h"
#include "./core_statmech_functions.h"
#include "./read_and_write.h"

int main() {
    int intervals = 1000; //number of intervals for the trapezoidal rule
    double degen = 2.0; //spin 1/2 for a fermion
    double mass = 938.918; // mean proton + neutron mass in MeV
    double temperature = 20.0; //temperature T in MeV
    double density = 122856.691; //density n_0
    const double tolerance = 1e-8;

    params.degeneracy_g_ = degen;
    params.density_ = density;
    params.mass_ = mass;
    params.temperature_ = temperature;

    double mu = -1000.0;

    for(int i = 0; i < 200; i++){
        std::cout << root_function_for_mu_bose(density,mu,mass,degen,temperature) << std::endl;
        mu -= i*100;
    }

    double chem_potent_mu_fermi = get_chem_potent_mu_fermi_bisection(tolerance,500.0,1500.0,params);
    std::cout << "Fermi Gas Chemical potential: " << chem_potent_mu_fermi << " MeV" << 
    std::endl;

    //Calculates the Fermi gas pressure using equation (68) and (72) of 
    //"Thermodynamics of the ideal Fermi gas" by Agnieszka Sorensen
    std::cout << "Fermi Gas Pressure = " << fermi_gas_pressure_integral_spherical(params, chem_potent_mu_fermi) << " MeV/fm^3" << std::endl;

   double chem_potent_mu_bose =
    get_chem_potent_mu_bose_bisection(tolerance,-1000,1500,temperature,density,degen,mass);
    std::cout << "Bose Gas Chemical potential: " << chem_potent_mu_bose << " MeV" << 
    std::endl;

    //Calculates the Fermi gas pressure using equations (68) and (78) of 
    //"Thermodynamics of the ideal Fermi gas" by Agnieszka Sorensen
    std::cout << "Bose Gas Pressure = " << 
    bose_gas_pressure_integral_spherical(temperature,
        chem_potent_mu_bose, degen, mass) << " MeV/fm^3" << std::endl;


    double read_col1[intervals];
    double read_col2[intervals];

    //Problem 4
    /*std::cout << "Problem 4" << std::endl;
    //Calculates
    //and writes it into a file
    double temps[15] = {1.0,10.0,25.0,50.0,75.0,100.0,125.0,150.0,175.0,200.0,225.0,250.0,275.0,300.0,325.0};
    double mu_T[15] = {0.0};

    for(int j = 0; j < 15; j++) {
        mu_T[j] = get_chem_potent_mu_fermi_bisection(tolerance,-100000.0,1500.0,temps[j],density,degen,mass);
        
    }

    write_array(15,temps,mu_T);
    read_array(15,read_col1,read_col2);
    */
    //Problem 6
    std::cout << "Problem 6" << std::endl;

    std::cout << "High Temperature Limit Chemical Potential:" << get_chem_potent_mu_boltzmann_bisection(tolerance,-100000.0,1500.0,temperature,density,degen,mass) << std::endl;


    //Problem 7
    std::cout << "Problem 7" << std::endl;
    //Calculating the chemical potential for various densities at
    // T = 1, 10, 50, 100, 200 MeV
    
    //Calculates the chemical potential at a temperature of T = 1 MeV
    //and writes it into a file
    double temp[5] = {1.0,10.0,50.0,100.0,200.0};
    double dense[7] = {0.0, 122856.691, 2457713.382, 3686570.073, 4915426.764, 6144283.455, 7373140.146}; //density n_B
    double mu_T_fermi[5] = {0.0};
    double mu_T_bose[5] = {0.0};
    double mu_highT[5] = {0.0};
    /*
    std::cout << "Ideal Fermi Gas Chemical Potentials" << std::endl;
    for(int j = 0; j < 5; j++) {
        for(int i = 0; i < 7; i++) {
            mu_T_fermi[i] = get_chem_potent_mu_fermi_bisection(tolerance, -10000.0,1500.0,temp[j],dense[i],degen,mass);
            //std::cout << "T = " << temp[j] << ", n = " << dense[i] << ", mu = " << mu_T_fermi[i] << std::endl;
        }
        //write_array(7,dense,mu_T_fermi);
        //read_array(7,read_col1,read_col2);
    }*/

    double dense_bose[7] = {0.0, 122856.691, 2457713.382, 3686570.073, 4915426.764, 6144283.455, 7373140.146}; //density n_B
    std::cout << "Ideal Bose Gas Chemical Potentials" << std::endl;
    for(int j = 0; j < 5; j++) {
        for(int i = 0; i < 7; i++) {
            mu_T_bose[i] = get_chem_potent_mu_bose_bisection(tolerance,-1500.0,1000.0,temp[j],dense_bose[i],degen,mass);
            std::cout << "T = " << temp[j] << ", n = " << dense_bose[i] << ", mu = " << mu_T_bose[i] << std::endl;
        }
        //write_array(7,dense_bose,mu_T_bose);
        //read_array(7,read_col1,read_col2);
    }
}