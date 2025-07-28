//////////////////////////////////////////////////////////////////////////////////////////
// This class is used to calculate kinematic variables in heavy-ion collisions.
//////////////////////////////////////////////////////////////////////////////////////////

#ifndef kinematic_variables_h
#define kinematic_variables_h

#include <stdexcept> // for using std::runtime_error
#include <iostream>



class KinematicVariables {
public:
  ////////////////////////////////////////////////////////////////////////////////////////
  // Default constructor; the constrctor uses provided variables to construct the object.
  KinematicVariables(double mass,
		     double sqrts,
		     double Ekin) {
    // Assign mass
    mass_ = mass;
    // Check which of the beam energy variables is provided
    if ( sqrts > 0.0 ) {
      // Check that Ekin is not provided
      if ( Ekin > 0.0 ) {
	throw std::runtime_error("Provide EITHER sqrts OR Ekin.");
      } else {
	sqrts_ = sqrts;
	//pbeam_CM = here is the routine calculating p_beam in the CM frame
	get_pbeamCM_from_sqrts();
    //Ebeam_CM = here is the routine calculating E_beam in the CM frame
    get_EbeamCM_from_sqrts();
    //vbeam_CM = here is the routine calculating v_beam in the CM frame
    get_vbeamCM_from_sqrts();
    //ybeam_CM = here is the routine calculating y_beam in the CM frame
	get_ybeamCM_from_sqrts();
    //ytarg_CM = here is the routine calculating y_targ in the CM frame
    get_y_coverage_CM();
    //gamma_CM = here is the routine calculating gamma in the CM frame
    get_gammaCM_from_sqrts();


	//Ekin_ = here is the routine calculating Ekin
	get_Ekin_from_sqrts();
    //Ebeam_FXT = here is the routine calculating E_beam in the FXT frame
    get_EbeamFXT_from_Ekin();
	//pbeam_FXT = here is the routine calculating p_beam in the FXT frame
	get_pbeamFXT_from_Ekin();
    //vbeam_FXT = here is the routine calculating v_beam in the FXT frame
    get_vbeamFXT_from_Ekin();
    //ybeam_FXT = here is the routine calculating v_beam in the FXT frame
    get_ybeamFXT_from_Ekin();
    //gamma_FXT = here is the routine calculating gamma in the FXT frame
    get_gammaFXT_from_Ekin();
      }    
    } else if ( Ekin > 0.0 ) {
      // Check that sqrts is not provided
      if ( sqrts > 0.0 ) {
	throw std::runtime_error("Provide EITHER sqrts OR Ekin.");
      } else {
	Ekin_ = Ekin;
	//Ebeam_FXT = here is the routine calculating E_beam in the FXT frame
    get_EbeamFXT_from_Ekin();
	//pbeam_FXT = here is the routine calculating p_beam in the FXT frame
	get_pbeamFXT_from_Ekin();
    //vbeam_FXT = here is the routine calculating v_beam in the FXT frame
    get_vbeamFXT_from_Ekin();
    //ybeam_FXT = here is the routine calculating v_beam in the FXT frame
    get_ybeamFXT_from_Ekin();
    //gamma_FXT = here is the routine calculating gamma in the FXT frame
    get_gammaFXT_from_Ekin();


	//sqrts_ = here is the routine calculating sqrts
    get_sqrts_from_Ekin();
	//pbeam_CM = here is the routine calculating p_beam in the CM frame
	get_pbeamCM_from_sqrts();
    //Ebeam_CM = here is the routine calculating E_beam in the CM frame
    get_EbeamCM_from_sqrts();
    //vbeam_CM = here is the routine calculating v_beam in the CM frame
    get_vbeamCM_from_sqrts();
    //ybeam_CM = here is the routine calculating y_beam in the CM frame
	get_ybeamCM_from_sqrts();
    //ytarg_CM = here is the routine calculating y_targ in the CM frame
    get_y_coverage_CM();
    //gamma_CM = here is the routine calculating gamma in the CM frame
    get_gammaCM_from_sqrts();    
      }      
    }
  }

  // Default destructor (releases the memory associated with the object when it is,
  // destroyed, i.e., when it goes out of scope)
  virtual ~KinematicVariables() {}


  
  ////////////////////////////////////////////////////////////////////////////////////////
  // Member functions of the class  
  ////////////////////////////////////////////////////////////////////////////////////////

  // Function that calculates Ekin_ given sqrts_
  void get_Ekin_from_sqrts();

  // Function that calculates p_beam_CM_ given sqrts_
  void get_pbeamCM_from_sqrts();
  
  // Function that calculates E_beam_CM_ given sqrts_
  void get_EbeamCM_from_sqrts();

  // Function that calculates v_beam_CM_ given sqrts_
  void get_vbeamCM_from_sqrts();

  // Function that calculates y_beam_CM_ given sqrts_
  void get_ybeamCM_from_sqrts();

  // Function that calculates y_beam - y_targ in CM frame
  void get_y_coverage_CM();

  // Function that calculates gamma_CM_ given sqrts_
  void get_gammaCM_from_sqrts();
  
  //Function that calculates E_beam_FXT_ given Ekin_
  void get_EbeamFXT_from_Ekin();

  // Function that calculates sqrts_ given Ekin_
  void get_sqrts_from_Ekin();

  // Function that calculates p_beam_FXT_ given Ekin_
  void get_pbeamFXT_from_Ekin();

  // Function that calculates v_beam_FXT_ given Ekin_
  void get_vbeamFXT_from_Ekin();

  // Function that calculates y_beam_FXT_ given Ekin_
  void get_ybeamFXT_from_Ekin();
  
  // Function that calculates y_beam - y_targ in FXT frame
  void get_y_coverage_FXT();

  // Function that calculates gamma_FXT_ given Ekin_
  void get_gammaFXT_from_Ekin();

  void print_variables();

  
  
  ////////////////////////////////////////////////////////////////////////////////////////
  // Return functions for private class members
  ////////////////////////////////////////////////////////////////////////////////////////

  double sqrts() { return sqrts_; }
  double y_beam_CM() { return y_beam_CM_; }
  double p_beam_CM() { return p_beam_CM_; }
  // etc., no need to define all of them -- only those you use outside of the class

  
 private:
  // Provided in the constructor
  double mass_{}; // brace initialization to a default value (0.0 for a double)
  double sqrts_{};
  double Ekin_{};
  // Calculated through the class methods:
  double p_beam_CM_{};
  double E_beam_CM_{};
  double v_beam_CM_{};
  double y_beam_CM_{};
  double y_coverage_CM_{};
  double gamma_CM_{};
  double p_beam_FXT_{};
  double E_beam_FXT_{};
  double v_beam_FXT_{};
  double y_beam_FXT_{};
  double y_coverage_FXT_{};
  double gamma_FXT_{};
};

#endif
