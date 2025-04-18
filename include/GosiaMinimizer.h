#ifndef GosiaMinimizer_h
#define GosiaMinimizer_h 1

#include "Nucleus.h"
#include "GosiaMinimizerFCN.h"
#include "Math/MinimizerOptions.h"
#include "Math/Minimizer.h"

class GosiaMinimizer {

 public:
  
  GosiaMinimizer(std::string meth, std::string alg);
  GosiaMinimizer(std::string meth);
  GosiaMinimizer();
  ~GosiaMinimizer();
  
  void Minimize();
  void ValidateErrors() {mini->SetValidError(true); validate = true;}
  void CalculateMINOSErrors() {calc_unc = true;}
  void WriteParameterSpace() {write = true;}
  void LimitMatrixElements() {limited = true;}
  void IncludeRelativeCS() {relCS = true;}
  void FixScalingParameters() {fixed = true;}
  
  void Scan(int index, double llim, double ulim, int nstep, bool inti = false);
  void Scan2D(int indexX, double llimX, double ulimX, int nstepX, int indexY, double llimY, double ulimY, int nstepY, bool inti = false);
  
  void SetMaximumIterations(int max) {maxIter = max;}
  void SetMaximumCalls(int max) {maxCalls = max;}
  void SetNumTrys(int num) {numTrys = num;}
  void SetFitTolerance(double tol) {fitTol = tol;}
  void SetChi2Cut(double cut) {chi_cut = cut; theFCN->chi_cut = cut;}
  
  void SetMethod(std::string str) {method = str;}
  void SetAlgorithm(std::string str) {algorithm = str;}
  void LinkBeamExperiments(int exp1, int exp2, double rel = 1.0);
  void LinkBeamTargetExperiments(int exp1B, int numE, int target_num);
  
  void SetBeam(Nucleus* nuc);
  void AddTarget(Nucleus* nuc);

  void SetBeamLiterature(Literature* lit) {theFCN->beam_lit = lit;}
  void AddTargetLiterature(Literature* lit) {theFCN->target_lits.push_back(lit);}

  void SetBeamData(ExperimentalData* dat);
  void AddTargetData(ExperimentalData* dat);
  
  void SetPrintLevel(int lvl) {mini->SetPrintLevel(lvl);}
  void Print() const {mini->PrintResults();}

 private:
  
  void InitialSetup();
  void LinkExperiments(int num1, int num2, double rel);
  
  void Resize(int size);
  void SetupParameters(std::vector<double> scales);
  std::vector<double> FindInitialScalings();
  void UpdateScalings(const double* min);

  bool limited; //flag to apply upper and lower limits to free MEs
  bool validate; //flag to validate symmetric errors (perform Hessian calculation)
  bool calc_unc; //flag to calculate asymmetric MINOS erros
  bool write; //flag to write parameter space to a ROOT file
  bool relCS; //flag to include RuthCS factor in relative couplings as GOSIA does
  bool fixed; //flag to fix the scaling parameters at their initial values

  int maxIter;
  int maxCalls;
  int numTrys;
  double fitTol;
  double chi_cut;

  GosiaMinimizerFCN* theFCN;

  std::vector<std::array<int,3>> targ_info;

  std::string method;
  std::string algorithm;
  ROOT::Math::Minimizer* mini;

};

#endif
