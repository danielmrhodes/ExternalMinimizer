#ifndef ExperimentalData_h
#define ExperimentalData_h 1

#include "Experiment.h"
#include <vector>
#include <string>

class ExperimentalData {

 public:
  
  ExperimentalData(std::string nm);
  ~ExperimentalData();

  int Size() const {return data.size();}
  std::string GetName() const {return name;}
  std::vector<double> GetScalings() const {return scalings;}
  Experiment* GetExperiment(int num) const {return data.at(num-1);}

  void SetScalings(std::vector<double> scls) {scalings = scls;}
  void RunGosiaCommands();
  void ReadGosiaFiles();
  void ReadDataFile();
  void Correct();
  void RecorrectExp(int exp, double scale);

  void PrintComparison() const;
  void PrintExpRawData(int index) const;
  void PrintRawData() const;
  void PrintRawAndCorrData() const;
  
  void PrintAllExpYields(int index) const;
  void PrintAllYields() const;

  friend class GosiaMinimizerFCN;
  //friend class GosiaMinimizer;
  
 private:

  void CreateFromGosiaInput();
  void FillFromIntiOutput(bool create = false);
  void FillFromPoinOutput(bool create = false);
  
  void DeriveCorrectionFactors();
  void SetFactors(std::vector<double> facs) {factors = facs;}
  
  std::string name;
  std::vector<Experiment*> data;
  std::vector<double> scalings; //fixed coupling x fitted scaling parameters
  std::vector<double> factors; //sin(theta_det)*Dsig(thetaP) factor needed for coupling experiments
  
};

#endif
