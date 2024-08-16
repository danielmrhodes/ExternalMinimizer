#ifndef Experiment_h
#define Experiment_h 1

#include "Yield.h"
#include <vector>
#include <string>

class Experiment {
  
 public:

  Experiment();
  ~Experiment();
  
  int GetTotalSize() const {return inti_yields_all.size();}
  int GetRawDataSize() const {return data_raw.size();}
  void PrintAllYields(std::string name) const;
  void PrintRawData(std::string name) const;
  void PrintRawAndCorrData(std::string name) const;
  void PrintComparison(std::string name, double scale, double factor) const;

  void SetNumber(int nm) {num = nm;}
  void SetEnergyMin(double en) {emin = en;}
  void SetEnergyMax(double en) {emax = en;}
  void SetThetaMin(double th) {tmin = th;}
  void SetThetaMax(double th) {tmax = th;}
  void SetRutherfordCS(double xs) {cs = xs;}
  
  int GetNumber() const {return num;}
  double GetEnergyMin() const {return emin;}
  double GetEnergyMax() const {return emax;}
  double GetThetaMin() const {return tmin;}
  double GetThetaMax() const {return tmax;}
  double GetRutherfordCS() const {return cs;}

  Yield* GetPointYield(int ni, int nf) const;
  Yield* GetIntiYield(int ni, int nf) const;

  std::vector<YieldError*> GetRawData() const {return data_raw;}
  std::vector<YieldError*> GetCorrectedData() const {return data_corr;}
  std::vector<Yield*> GetPointYields() const {return point_yields;}
  
  std::vector<Yield*> GetAllPointYields() const {return point_yields_all;}
  std::vector<Yield*> GetAllIntiYields() const {return inti_yields_all;}
  std::vector<double> GetCorrectionFactors() const {return corr_facs;}
  
  friend class ExperimentalData;
  friend class GosiaMinimizerFCN;
  
 private:

  void AddRawData(YieldError* yld) {data_raw.push_back(yld);}
  void AddCorrData(YieldError* yld) {data_corr.push_back(yld);}
  
  void AddIntiYield_All(Yield* yld) {inti_yields_all.push_back(yld);}
  void UpdateIntiYield_All(int index, double val) {inti_yields_all[index]->SetValue(val);}

  void AddPointYield_All(Yield* yld) {point_yields_all.push_back(yld);} 
  void UpdatePointYield_All(int index, double val) {point_yields_all[index]->SetValue(val);}

  void AddPointYield(Yield* yld) {point_yields.push_back(yld);} 
  void UpdatePointYield(int index, double val) {point_yields[index]->SetValue(val);}  

  void ClearCorrectionFactors() {corr_facs.clear();}
  void ClearAllIntiYields();
  void ClearAllPointYields();
  void ClearPointYields();
  void ClearRawData();
  void ClearCorrectedData();

  void ZeroAllIntiYields();
  void ZeroAllPointYields();
  void ZeroPointYields();
  
  //bool targ_det;
  int num, beam_Z, targ_Z;
  double emin, emax, tmin, tmax, cs, Ep, thP, beam_mass, targ_mass;
  
  std::vector<YieldError*> data_raw, data_corr;
  std::vector<Yield*> point_yields, point_yields_all, inti_yields_all;
  std::vector<double> corr_facs;
  
};

#endif
