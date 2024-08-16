#ifndef Literature_h
#define Literature_h 1

#include <vector>
#include <string>

class LitVal {

 public:

  LitVal();
  LitVal(int indI, double val, double erU, double erD);
  LitVal(int indI, int indF, double val, double erU, double erD);
  LitVal(int indI, int indF1, int indF2, double val, double erU, double erD);
  ~LitVal();
  
 int GetInitialIndex() const {return indices[0];}
 int GetFinalIndex() const {return indices[1];}
 int GetFinalIndexTwo() const {return indices[2];}
 
 double GetValue() const {return value;}
 double GetErrorUp() const {return errUp;}
 double GetErrorDown() const {return errDn;}
 double Compare(double calc) const;
 double NSigma(double calc) const;

 /*
 void SetInitialIndex(int ind) {indices[0] = ind;}
 void SetFinalIndex(int ind) {indices[1] = ind;}
 void SetFinalIndexTwo(int ind) {indices[2] = ind;}
 void SetValue(double val) {value = val;}
 void SetErrorUp(double err) {errUp = err;}
 void SetErrorEn(double err) {errDn = err;}
 */
 
 private:
 
 std::vector<int> indices;
 double value, errUp, errDn; 
  
};

class Literature {

 public:
  
  Literature(std::string nm);
  ~Literature();

  void CreateFromFile();
  
  void AddLifetime(int index, double val, double erU, double erD);
  void AddBranchingRatio(int init, int fin1, int fin2, double val, double erU, double erD);
  void AddMixingRatio(int init, int fin, double val, double erU, double erD);
  void AddMatrixElement(int ind1, int ind2, int mult, double val, double erU, double erD);
  
  int Size() const;
  std::vector<double> GetWeights() const {return weights;}
  std::vector<LitVal*> GetLifetimes() const {return lifetimes;}
  std::vector<LitVal*> GetBranchingRatios() const {return branching_ratios;}
  std::vector<LitVal*> GetMixingRatios() const {return mixing_ratios;}
  std::vector<LitVal*> GetMatrixElements() const {return matrix_elements;}

 private:
  
  std::string name;
  std::vector<double> weights;
  std::vector<LitVal*> lifetimes;
  std::vector<LitVal*> branching_ratios;
  std::vector<LitVal*> mixing_ratios;
  std::vector<LitVal*> matrix_elements;
  
};

#endif
