#ifndef Fermi_gas_functions_Tzero_h
#define Fermi_gas_functions_Tzero_h



double E_kinetic(double mass, double momentum);

// Fermi momentum at baryon density n_B (T=0)
double fermi_momentum(double degeneracy, double density);

// Fermi energy at baryon density n_B (T=0)
double fermi_energy(double degeneracy, double mass, double density);

// Energy density of an ideal (non-interacting) Fermi gas at baryon density n_B (T=0)
double Fermi_gas_energy_density(double degeneracy, double mass, double density);

// Derivative of the Fermi pressure w.r.t. baryon density nB (T=0)
double d_PFG_d_nB(double degeneracy, double mass, double density);


#endif
