#include <iostream> //gives standard C++ I/O: cin, cout
#include <cstdlib> //gives access to the C standard library
#include <cmath>
using namespace std; //don't have to put std in front of everything

struct results
{
    float E_kin; //Beam kinetic energy in fixed target frame
    float sqrts; //sqrt{s}, relates the energy in the collider frame to the energy in the fixed target frame
    float E_cm; //Energy of the particle in collider frame
    float E_fixed; //Energy of the beam in fixed target frame
    float p_beamCM; //Beam Momentum in center-of-mass frame
    float p_fixed;
    float beta_beamCM; //beam velocity in center-of-mass frame
    float beta_fixed;
    float y_beamCM; //beam rapidity in center-of-mass frame
    float lorentzFactorCM; //Lorentz contraction factor gamma
};

int main() {
    
    float SfromKinetic(float E_kin); //Calculates sqrt(S) given beam kinetic energy E_kin (not including mass)
    float KineticfromS(float sqrts); //Calculates E_kin given sqrt(S)
    results KinematicVariablesCollider(); //Calculates beam momentum, beam velocity, beam rapidity, and Lorentz contraction factor gamma when given sqrt(S) and E_kin in the Collider mode
    results KinematicVariablesFXT(); //Calculates beam momentum, beam velocity, beam rapidity, and Lorentz contraction factor gamma when given sqrt(S) and E_kin in the FXT mode

    results collide = KinematicVariablesCollider();
    results fixed = KinematicVariablesFXT();

    //cout << "sqrt(s) calculated from given E_kin: " << SfromKinetic(1.24) << endl;
    //cout << "E_kin calculated from given sqrt(s): " << KineticfromS(2.42) << endl;
}

float SfromKinetic(float E_kin){ //Calculates sqrt(S) given beam kinetic energy E_kin (not including mass)
    //Used for target mode
    //Assuming particles are protons
    float sqrts = 0.0;
    float E = 0.0;
    float m = 0.938272; //in GeV
    E = E_kin + m;
    cout << "E_fixed = " <<  E << endl;
    sqrts = sqrt(pow(m,2) + pow(m,2) + (2*E*m)); //Gives the square root of the total energy in the center-of-mass frame in a 2 particle system
    return sqrts;
}

float KineticfromS(float sqrts){ //Calculates E_kin given sqrt(S)
    //Used for collider mode
    //Assuming particles are protons
    float E_kin = 0.0;
    float m = 0.938272; //in GeV
    //Knowing that sqrt(S) = 2E (assuming m1 = m2 and p1 = p2)
    cout << "E_cm = " << sqrts/2 << endl;
    E_kin = (pow(sqrts,2)/(2*m)) - 2*m; //Gives the kinetic energy for one particle, for the total kinetic energy multiply by 2
    return E_kin;
}

results KinematicVariablesCollider(){ //Calculates beam momentum, beam velocity, beam rapidity, and Lorentz contraction factor gamma when given sqrt(S) and E_kin
    //Used for collider mode
    //Assuming particles are protons
    float m = 0.938272; //in GeV
    float sqrts = 0.0;
    cout << "Enter sqrts: ";
    cin >> sqrts;
    results collide;

    //Need to add: Total rapidity Coverage, y_proj - y_target

    collide.sqrts = sqrts;
    //Knowing that sqrt(S) = 2E (assuming m1 = m2 and p1 = p2)
    collide.E_cm = sqrts/2;
    collide.E_kin = (pow(sqrts,2)/(2*m)) - 2*m; //Gives the kinetic energy for one particle, for the total kinetic energy multiply by 2
    collide.E_fixed = collide.E_kin + m;
    collide.p_beamCM = sqrt((pow(sqrts,2)/4) - pow(m,2)); //Beam Momentum in center-of-mass frame
    collide.beta_beamCM = collide.p_beamCM/collide.E_cm; //beam velocity in center-of-mass frame
    collide.y_beamCM = 0.5*log((1 + collide.beta_beamCM)/(1 - collide.beta_beamCM));//beam rapidity in center-of-mass frame
    collide.lorentzFactorCM = 1/sqrt(1 - pow(collide.beta_beamCM,2)); //Lorentz factor gamma

    cout << "For fixed collder experiments, given sqrts = " << collide.sqrts << ", find: " << endl;
    cout << "E_kin: " << collide.E_kin << endl;
    cout << "E_fixed: " << collide.E_fixed << endl;
    cout << "E_cm: " << collide.E_cm << endl;
    cout << "p_beamCM: " << collide.p_beamCM << endl;
    cout << "beta_beamCM: " << collide.beta_beamCM << endl;
    cout << "y_beamCM: " << collide.y_beamCM << endl;
    cout << "lorentzFactorCM: " <<collide.lorentzFactorCM << endl;

    return collide;
}

results KinematicVariablesFXT(){ //Calculates beam momentum, beam velocity, beam rapidity, and Lorentz contraction factor gamma when given sqrt(S) and E_kin
    float m = 0.938272; //in GeV
    float E_kin = 0.0;
    cout << "Enter E_kin: ";
    cin >> E_kin;
    results fixed;
    
    //Need to add: Total rapidity Coverage, y_proj - y_target, beta_fixed

    fixed.E_kin = E_kin;
    fixed.E_fixed = fixed.E_kin + m;
    fixed.sqrts = sqrt(pow(m,2) + pow(m,2) + (2*fixed.E_fixed*m)); //Gives the square root of the total energy in the center-of-mass frame in a 2 particle system;
    fixed.E_cm = fixed.sqrts/2; //Energy of 

    fixed.p_beamCM = sqrt(pow(E_kin,2) + 2*E_kin*m); //Beam Momentum in center-of-mass frame
    
    fixed.p_fixed = (fixed.p_beamCM*m)/fixed.sqrts; // Beam Momentum in fixed target frame
    std::cout << fixed.p_fixed << std::endl;
    std::cin.get();
    std::cin.get();
    std::cin.get();
    fixed.beta_beamCM = fixed.p_beamCM/(fixed.E_cm); //beam velocity in center-of-mass frame
    
    fixed.beta_fixed = fixed.p_fixed/fixed.E_fixed; //Why does the velocity of the beam use the fixed frame energy
    fixed.y_beamCM = 0.5*log((1 + fixed.beta_beamCM)/(1 - fixed.beta_beamCM));//beam rapidity in center-of-mass frame
    fixed.lorentzFactorCM = 1/sqrt(1 - pow(fixed.beta_beamCM,2)); //Lorentz factor gamma

    cout << "For fixed target experiments, given E_kin = " << E_kin << ", find: " << endl;
    cout << "E_fixed: " << fixed.E_fixed << endl;
    cout << "sqrt(s): " << fixed.sqrts << endl;
    cout << "E_cm: " << fixed.E_cm << endl;
    cout << "p_beamCM: " << fixed.p_beamCM << endl;
    cout << "p_fixed:  " << fixed.p_fixed << endl;
    cout << "beta_beamCM: " << fixed.beta_beamCM << endl;
    cout << "beta_fixed: " << fixed.beta_fixed << endl;
    cout << "y_beamCM: " << fixed.y_beamCM << endl;
    cout << "lorentzFactorCM: " << fixed.lorentzFactorCM << endl;

    return fixed;
}