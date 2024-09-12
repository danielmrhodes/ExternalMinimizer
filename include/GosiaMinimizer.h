#ifndef GosiaMinimizer_h
#define GosiaMinimizer_h 1

#include "Nucleus.h"
#include "GosiaMinimizerFCN.h"
#include "Math/MinimizerOptions.h"
#include "Math/Minimizer.h"

class GosiaMinimizer {

 public:
  
  GosiaMinimizer();
  GosiaMinimizer(std::string meth, std::string alg);
  ~GosiaMinimizer();
  
  void Minimize();
  void CalculateFullUncertainty() {calc_unc = true;}
  void WriteParameterSpace() {write = true;}
  void LimitMatrixElements() {limited = true;}
  
  void Print() const {mini->PrintResults();}
  
  void SetMaximumIterations(int max) {maxIter = max;}
  void SetMaximumCalls(int max) {maxCalls = max;}
  void SetNumTrys(int num) {numTrys = num;}
  void SetFitTolerance(double tol) {fitTol = tol;}
  void SetChi2Cut(double cut) {chi_cut = cut; theFCN->chi_cut = cut;}
  
  void SetMethod(std::string str) {method = str;}
  void SetAlgorithm(std::string str) {algorithm = str;}
  void LinkBeamExperiments(int exp1, int exp2, double rel);
  void LinkBeamTargetExperiments(int exp1B, int numE, int target_num);
  
  void SetBeam(Nucleus* nuc);
  void AddTarget(Nucleus* nuc);

  void SetBeamLiterature(Literature* lit) {theFCN->beam_lit = lit;}
  void AddTargetLiterature(Literature* lit) {theFCN->target_lits.push_back(lit);}

  void SetBeamData(ExperimentalData* dat);
  void AddTargetData(ExperimentalData* dat);
  
 private:
  
  void LinkExperiments(int num1, int num2, double rel);
  
  void Resize(int size);
  void SetupParameters(std::vector<double> scales);
  std::vector<double> FindInitialScalings();
  void UpdateScalings();

  bool limited; //flag for applying upper and lower limits to free MEs
  bool calc_unc;
  bool write;

  int maxIter;
  int maxCalls;
  int numTrys;
  double fitTol;
  double chi_cut;

  GosiaMinimizerFCN* theFCN;
  
  std::string method;
  std::string algorithm;
  ROOT::Math::Minimizer* mini;

};

#endif
