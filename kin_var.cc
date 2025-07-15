#include <iostream>    // standard input-output stream
#include <fstream>     // handles file input-output: file creation, write to files, etc.
#include <iomanip>     // manipulates the input-output; here used for setprecision 

#include "./kinematic_variables.h"



int main () {

  std::cout << "\nHello there! This is my kinematic variables code\n" << std::endl;

  ////////////////////////////////////////////////////////////////////////////////////////
  // Define your variables
  double mass = 0.938918; // mean proton + neutron mass

  // Simply hard-code what you provide; simple solution: provide one or the other
  // (but you can surely come up with something better)
  double sqrts = 14000.0;
  double Ekin = 0.0;

  ////////////////////////////////////////////////////////////////////////////////////////
  // Calculate kinematic variables

  // Construct a Kinematic_Variables object
  KinematicVariables HADES_experiment(mass, sqrts, Ekin);
  HADES_experiment.print_variables();

  return 0;
}

