#ifndef MatrixElement_h
#define MatrixElement_h 1

#include <string>

class MatrixElement {

 public:
  
  MatrixElement();
  ~MatrixElement();
  
  bool IsFixed() const {return fixed;}
  bool IsLimited() const {return limited;}
  int GetIndex1() const {return index1;}
  int GetIndex2() const {return index2;}
  int GetMultipolarity() const {return mult;}
  std::string GetMultS() const;
  double GetValue() const {return value;}
  double GetUpperLimit() const {return ulim;}
  double GetLowerLimit() const {return llim;}

  void Limit() {limited = true;}
  void Fix() {fixed = true;}
  void Release() {fixed = false;}
  void SetIndex1(int ind) {index1 = ind;}
  void SetIndex2(int ind) {index2 = ind;}
  void SetMultipolarity(int mlt) {mult = mlt;} 
  void SetValue(double v) {value = v;}

  void ApplyLimits(double ll, double ul) {limited = true; llim = ll; ulim = ul;}
  void SetLimits(double ll, double ul) {llim = ll; ulim = ul;}
  void SetUpperLimit(double lim) {ulim = lim;}
  void SetLowerLimit(double lim) {llim = lim;}

  static std::string GetMultS(int l);
  static bool Compare(const MatrixElement* me1, const MatrixElement* me2);
  
 private:

  bool fixed, limited;
  int index1, index2, mult;
  double value, ulim, llim;
  
};

#endif
