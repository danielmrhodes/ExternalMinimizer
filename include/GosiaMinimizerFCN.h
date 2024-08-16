#ifndef GosiaMinimizerFCN_h
#define GosiaMinimizerFCN_h 1

#include "Nucleus.h"
#include "ExperimentalData.h"
#include "TFile.h"
#include "TNtuple.h"

class GosiaMinimizerFCN {
  
 public:
  
  GosiaMinimizerFCN();
  ~GosiaMinimizerFCN();

  double ErrorDef() const {return chi_cut;}
  double GetErrorDef() const {return chi_cut;}
  
  double up() const {return chi_cut;}
  double operator() (const double* pars);
  void SetErrorDef(double err) {chi_cut = err;}
  
  friend class GosiaMinimizer;

 private:
  
  void CreateNTuple();
  void FillFactors();
  
  double chi_cut;

  std::string beam_name;
  Nucleus* beam;
  Literature* beam_lit;
  ExperimentalData* beam_data;

  TFile* outFile;
  TNtuple* ntuple;
  
  std::vector<std::string> target_names;
  std::vector<Nucleus*> targets;
  std::vector<Literature*> target_lits;
  std::vector<ExperimentalData*> target_data;

  std::vector<int> indices; //list of which free scaling parameter goes with which experiments
  std::vector<double> scalings; //the fixed relative scaling between experiments
  std::vector<double> factors; //sin(theta_det)*Dsig(thetaP) factor needed for coupling experiments

};

#endif
