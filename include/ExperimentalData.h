#ifndef ExperimentalData_h
#define ExperimentalData_h 1

#include "Experiment.h"
#include "Nucleus.h"
#include "TF1.h"
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

  void SetNormalizationIndex(int index) {nrm_index = index;}
  int GetNormalizationIndex() const {return nrm_index;}

  void SetScalings(std::vector<double> scls) {scalings = scls;}
  void RunGosiaCommands();

  void ReadGosiaFiles(std::string poin_inp, std::string poin_out, std::string inti_out);
  void ReadGosiaFiles();

  void ReadDataFile(std::string file_name);
  void ReadDataFile();
  
  void GenerateData(std::string file_name, Nucleus* nuc, TF1* eff, double scale);
  void GenerateData(Nucleus* nuc, TF1* eff, double scale);
  void GenerateAllData(Nucleus* nuc, TF1* eff, double scale);

  void WriteDataFile(std::string file_name,int A, int Z);  
  void WriteDataFile(int A, int Z);

  void Correct();
  void RecorrectExp(int exp, double scale);

  void PrintComparison() const;
  void PrintExpRawData(int index) const;
  void PrintRawData() const;
  void PrintRawAndCorrData() const;
  
  void PrintAllExpYields(int index) const;
  void PrintAllYields() const;

  friend class GosiaMinimizerFCN;
  
 private:

  void CreateFromGosiaInput(std::string file_name);
  void FillFromIntiOutput(std::string file_name, bool create = false);
  void FillFromPoinOutput(std::string file_name, bool create = false);

  void CreateFromGosiaInput();
  void FillFromIntiOutput(bool create = false);
  void FillFromPoinOutput(bool create = false);

  void DeriveCorrectionFactors();
  void SetFactors(std::vector<double> facs) {factors = facs;}
  
  int nrm_index;
  std::string name;
  std::vector<Experiment*> data;
  std::vector<double> scalings; //fixed coupling x fitted scaling parameters
  std::vector<double> factors; //sin(theta_det)*Dsig(thetaP) factor needed for coupling experiments
  
};

#endif
