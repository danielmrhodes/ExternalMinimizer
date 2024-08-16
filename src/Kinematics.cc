#include "Kinematics.h"
#include "TMath.h"

Kinematics::Kinematics() {;}
Kinematics::~Kinematics() {;}

double Kinematics::Theta_CM_FP(double ThetaLAB, bool sol2) {
  
  double tau = (beam_mass/targ_mass)/std::sqrt(1 - (Ex/Ep)*(1 + beam_mass/targ_mass));
  if(std::sin(ThetaLAB) > 1.0/tau) {

    ThetaLAB = std::asin(1.0/tau);
    if(ThetaLAB < 0)
      ThetaLAB += TMath::Pi();  

    return std::asin(tau*std::sin(ThetaLAB)) + ThetaLAB;
  }

  if(sol2)
    return std::asin(tau*std::sin(-ThetaLAB)) + ThetaLAB + TMath::Pi();
  
  return std::asin(tau*std::sin(ThetaLAB)) + ThetaLAB;
  
}

double Kinematics::Recoil_Theta_LAB(double thetaCM) {

  double tau = 1.0/std::sqrt(1 - (Ex/Ep)*(1 + beam_mass/targ_mass));
  double tanTheta = std::sin(TMath::Pi() - thetaCM)/(std::cos(TMath::Pi() - thetaCM) + tau);
  
  return std::atan(tanTheta);
  
}

double Kinematics::Dsig(double thetaP, bool sol2, bool targ_det) {
  
  double thetaCM = Theta_CM_FP(thetaP,sol2);
  double r3; //jacobian
  if(targ_det) {
    double zcm = TMath::Pi() - thetaCM;
    double zlb = Recoil_Theta_LAB(thetaCM);

    double term1 = TMath::Power(TMath::Sin(zlb)/TMath::Sin(zcm),2.0);
    double term2 = TMath::Abs(TMath::Cos(zcm - zlb));
    r3 = 1.0/(term1*term2);
  }
  else {
    double num = TMath::Power(TMath::Sin(thetaCM)/TMath::Sin(thetaP),2.0);
    double denom = TMath::Abs(TMath::Cos(thetaCM-thetaP));
    r3 = num/denom; 
  }
  
  double ared = 1.0 + beam_mass/targ_mass;
  double dista = 0.0719949*ared*beam_Z*targ_Z/Ep;
  double eps = 1.0/TMath::Sin(thetaCM/2.0);

  return 250.0*r3*sqrt(Ep/(Ep - ared*Ex))*dista*dista*pow(eps,4);
  
  //Need sin(theta_det)*Dsig(thetaP) for scaling parameters
  
}
