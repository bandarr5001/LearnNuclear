//////////////////////////////////////////////////////////////////////////////////////////
// This class is used to calculate kinematic variables in heavy-ion collisions.
//////////////////////////////////////////////////////////////////////////////////////////

#include "./kinematic_variables.h"

#include <iostream> // input, output (std::cout etc.)
#include <cstdlib>
#include <cmath>



void KinematicVariables::get_Ekin_from_sqrts() {
    Ekin_ = (sqrts_ * sqrts_)/(2.0 * mass_) - 2.0 * mass_;
}

void KinematicVariables::get_pbeamCM_from_sqrts() {
    p_beam_CM_ = sqrt((pow(sqrts_,2)/4) - pow(mass_,2));
}

void KinematicVariables::get_EbeamCM_from_sqrts() {
    E_beam_CM_ = sqrts_/2;
}

void KinematicVariables::get_vbeamCM_from_sqrts() {
    v_beam_CM_ = p_beam_CM_/E_beam_CM_;
}

void KinematicVariables::get_ybeamCM_from_sqrts() {
    y_beam_CM_ = 0.5*log((1 + v_beam_CM_) / (1 - v_beam_CM_));
}

void KinematicVariables::get_y_coverage_CM() {
    double y_targ_CM_ = 0.5*log((1 - v_beam_CM_) / (1 + v_beam_CM_));
    y_coverage_CM_ = y_beam_CM_ - y_targ_CM_;
}

void KinematicVariables::get_gammaCM_from_sqrts() {
    gamma_CM_ = 1 / sqrt(1 - pow(v_beam_CM_,2));
}

void KinematicVariables::get_EbeamFXT_from_Ekin() {
    E_beam_FXT_ = Ekin_ + mass_;
}

void KinematicVariables::get_sqrts_from_Ekin() {
    sqrts_ = sqrt(2*pow(mass_,2) + 2*E_beam_FXT_*mass_);
}

void KinematicVariables::get_pbeamFXT_from_Ekin() {
    p_beam_FXT_ = sqrt(pow(E_beam_FXT_,2) - pow(mass_,2));
}

void KinematicVariables::get_vbeamFXT_from_Ekin() {
    v_beam_FXT_ = p_beam_FXT_/E_beam_FXT_;
}

void KinematicVariables::get_ybeamFXT_from_Ekin() {
    y_beam_FXT_ = 0.5*log((1 + v_beam_FXT_) / (1 - v_beam_FXT_));
}

void KinematicVariables::get_gammaFXT_from_Ekin() {
    gamma_FXT_ = 1 / sqrt(1 - pow(v_beam_FXT_,2));
}

//Also need to get, Lab/Total Rapidity Coverage and y_proj - y_targ (for both FXT and CM),

void KinematicVariables::print_variables() {
  std::cout << "***********************************************************" << std::endl;
  std::cout << "\n          sqrts = " << sqrts_
        << "\n      E_beam_CM = " << E_beam_CM_
	    << "\n      p_beam_CM = " << p_beam_CM_
	    << "\n      v_beam_CM = " << v_beam_CM_
	    << "\n      y_beam_CM = " << y_beam_CM_
        << "\n  y_coverage_CM = " << y_coverage_CM_
	    << "\n       gamma_CM = " << gamma_CM_
	    << "\n"
	    << "\n           Ekin = " << Ekin_
        << "\n     E_beam_FXT = " << E_beam_FXT_
	    << "\n     p_beam_FXT = " << p_beam_FXT_
	    << "\n     v_beam_FXT = " << v_beam_FXT_
	    << "\n     y_beam_FXT = " << y_beam_FXT_
        << "\n y_coverage_FXT = " << y_beam_FXT_
	    << "\n      gamma_FXT = " << gamma_FXT_
	    << "\n" << std::endl;
}

