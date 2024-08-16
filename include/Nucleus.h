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
  void PrintAll() const; 
  void PrintComparison(const Literature* lit) const;

  inline double DecayLifetime(int mult, double me, double spin, double egam) const;

  int GetNumMatrixElements() const {return matrix_elements.size();}
  int GetNumLevels() const {return levels.size();}
  int GetNumFree() const;
  
  double CompareWithLiterature(const Literature* lit) const;
  double CompareLifetime(const LitVal* lv) const;
  double CompareBranchingRatio(const LitVal* lv) const;
  double CompareMixingRatio(const LitVal* lv) const;
  double CompareMatrixElement(const LitVal* lv) const;

  double CalculateLifetime(int index) const;
  double CalculateBranchProbability(int ni, int nf) const;
  double CalculateBranchingRatio(int ni, int nf1, int nf2) const;
  double CalculateMixingRatio(int ni, int nf) const;

  std::vector<Level*> GetLevels() const {return levels;}
  
  MatrixElement* GetMatrixElement(int index) const {return matrix_elements[index];}
  //std::vector<MatrixElement*> GetMatrixElements() const {return matrix_elements;}
  float GetMatrixElementValue(int index) const {return matrix_elements[index]->GetValue();}
  std::vector<float> GetMatrixElementValues() const;
  
  std::string GetName() const {return name;}
  void SetName(std::string nm) {name = nm;}

  friend class NucleusPlotter;

  private:
  
  MatrixElement* GetMatrixElement(int n1, int n2, int mt) const;
  //void Sort();

  std::array<TSpline3*,8> convCoeffs;

  std::string name;
  //std::map<std::array<int,3>,int> the_map; //{mult,n1,n2} <-> vector index
  std::vector<MatrixElement*> matrix_elements;
  std::vector<Level*> levels;
  
};

#endif
