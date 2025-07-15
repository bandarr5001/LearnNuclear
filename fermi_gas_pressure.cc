#include <iostream>
#include <cmath>
#include <cstdlib>
#include "./fermi_gas.h"
#include "./read_and_write.h"

int main() {
    int intervals = 1000; //number of intervals for the trapezoidal rule
    double temp = 1;//temperature T in MeV
    double degen = 2.0; //spin 1/2 for a fermion
    double mass = 938.918; // mean proton + neutron mass in MeV
    double density = 1.0; //density n_B
    
   

    double chem_potent_mu_fermi = get_chem_potent_mu_fermi(temp,density,degen);
    std::cout << "chemical potential: " << chem_potent_mu_fermi << " MeV" << 
    std::endl;

    std::cout << "Fermi Gas Pressure = " << 
    fermi_gas_pressure_integral_spherical(intervals, temp, chem_potent_mu_fermi,
        degen, mass) << " MeV/fm^3" << std::endl;

    //Use equation 46 to calculate the fermi gas pressure as a function of the 
    //chemical potential mu
    
    double p_fermi = 200.0;

    double h = p_fermi/intervals;
    //std::cout << "h = " << h << std::endl;
    double P_T0[intervals];
    double chem_potential_[intervals];

    double mu = 0.0;

    for(int i = 0; i < intervals; i += 1) {
        chem_potential_[i] = mu;
        P_T0[i] = fermi_gas_press_temp0_chem(mu, degen, mass);
        //std::cout << P_T0[i] << std::endl;
        mu += h;
    }

    double read_col1_[intervals];
    double read_col2_[intervals];

    write_array(intervals,chem_potential_,P_T0);
    read_array(intervals,read_col1_,read_col2_);
}