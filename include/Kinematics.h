#ifndef Kinematics_h
#define Kinematics_h 1

class Kinematics {
  
 public:
  
  Kinematics();
  ~Kinematics();
  
  double Theta_CM_FP(double ThetaLAB, bool sol2 = false);
  double Recoil_Theta_LAB(double thetaCM);
  double Dsig(double thetaP, bool sol2 = false, bool targ_det = false);
  
  /*
  double Theta_CM_FR(double ThetaLAB, bool sol2 = false);
  double Theta_LAB(double thetaCM);
  double Theta_LAB_Max();
  
  double Recoil_Theta_LAB_Max();
  double KE_LAB(double thetaCM);
  double Recoil_KE_LAB(double thetaCM);
  double RuthCM(double thetaCM);
  double RuthLAB(double thetaLAB, bool sol2=false);
  double Beta_LAB(double thetaCM);
  double Recoil_Beta_LAB(double thetaCM);
  */
  
  void SetBeamZ(double z) {beam_Z = z;}
  void SetBeamMass(double ms) {beam_mass = ms;}
  void SetTargetZ(double z) {targ_Z = z;}
  void SetTargetMass(double ms) {targ_mass = ms;}
  void SetBeamEnergy(double en) {Ep = en;}
  void SetEx(double ex) {Ex = ex;}

 private:
  
  double beam_Z, beam_mass, targ_Z, targ_mass, Ep, Ex;
  
};

#endif
