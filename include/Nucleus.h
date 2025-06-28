#ifndef Nucleus_h
#define Nucleus_h 1

#include "Literature.h"
#include "MatrixElement.h"
#include "Level.h"
#include "TSpline.h"
#include <vector>
#include <map>

class Nucleus {
  
 public:
  
  Nucleus(std::string nm);
  ~Nucleus();
  
  void CreateFromGosiaOutputFile(std::string file_name);
  void CreateFromGosiaInputFile(std::string file_name);
  void FillFromBSTFile(std::string file_name);

  void CreateFromGosiaOutputFile();
  void CreateFromGosiaInputFile();
  void FillFromBSTFile();
  
  void ReleaseMatrixElement(int n1, int n2, int mt);
  void FixMatrixElement(int n1, int n2, int mt);
  void ReleaseAll();
  void FixAll();

  void SetConverionCoefficients(int mult, std::vector<double> ens, std::vector<double> ccs);

  void PrintLevels() const;
  void PrintLifetimes() const;
  void PrintBranchProbabilities() const;
  void PrintBranchingRatios() const;
  void PrintMixingRatios() const;
  void PrintMatrixElements() const;
  void PrintGammaTransitions() const;
  void PrintAll() const; 
  void PrintComparison(const Literature* lit) const;

  void Write(std::string file_name, int A, int Z) const;
  void WriteLevelScheme(std::string file_name) const;
  void Write(int A, int Z) const;
  void WriteLevelScheme() const;

  bool CheckMatrixElements() const;

  static double Moment(double dme, double spin, int mult);
  double DecayLifetime(int mult, double me, double spin, double egam) const;

  int GetNumMatrixElements() const {return matrix_elements.size();}
  int GetNumLevels() const {return levels.size();}
  int GetNumFree() const;

  double GetLevelEnergy(int index) const {return levels.at(index-1)->GetEnergy();}
  
  double CompareWithLiterature(const Literature* lit) const;
  double CompareLifetime(const LitVal* lv) const;
  double CompareBranchingRatio(const LitVal* lv) const;
  double CompareMixingRatio(const LitVal* lv) const;
  double CompareMatrixElement(const LitVal* lv) const;

  double CalculateLifetime(int index) const;
  double CalculateBranchProbability(int ni, int nf) const;
  double CalculateBranchingRatio(int ni, int nf1, int nf2) const;
  double CalculateMixingRatio(int ni, int nf) const;
  
  double ConvCoef(int mult, double egam) const {return convCoeffs[mult-1]->Eval(egam);}
  std::vector<int> GetTransitionMultipolarities(int ni, int nf) const;
  std::vector<Level*> GetLevels() const {return levels;}
  
  MatrixElement* GetMatrixElement(int index) const {return matrix_elements[index];}
  MatrixElement* GetMatrixElement2(int n1, int n2, int mt) const {return GetMatrixElement(n1-1,n2-1,mt);}
  std::vector<MatrixElement*> GetMatrixElements() const {return matrix_elements;}
  
  float GetMatrixElementValue(int index) const {return matrix_elements[index]->GetValue();}
  std::vector<float> GetMatrixElementValues() const;
  
  std::string GetName() const {return name;}
  void SetName(std::string nm) {name = nm;}

  friend class NucleusPlotter;

  private:
  
  MatrixElement* GetMatrixElement(int n1, int n2, int mt) const;

  std::array<TGraph*,8> convCoeffs;
  
  std::string name;
  std::vector<MatrixElement*> matrix_elements;
  std::vector<Level*> levels;
  
};

#endif
