#ifndef Yield_h
#define Yield_h 1

#include <vector>

class Yield {

 public:
  
  Yield();
  ~Yield();

  void AddInitialIndex(int ind) {nis.push_back(ind);}
  void AddFinalIndex(int ind) {nfs.push_back(ind);}
  void SetValue(double yl) {yld = yl;}
  void SetThreshold(double th) {thresh = th;}
  void SetObserved() {observed = true;}

  bool IsObserved() const {return observed;}
  int Size() const {return nis.size();}
  double GetValue() const {return yld;}
  double GetThreshold() const {return thresh;}

  std::vector<int> GetInitialIndices() const {return nis;}
  std::vector<int> GetFinalIndices() const {return nfs;}
  
    
 private:

  bool observed;
  double yld, thresh;
  std::vector<int> nis, nfs;

};

class YieldError : public Yield {

 public:
  
  YieldError();
  ~YieldError();

  void SetErrorUp(double er) {erUp = er;}
  void SetErrorDown(double er) {erDn = er;}
  void SetWeight(double wt) {weight = wt;}
  
  double GetErrorUp() const {return erUp;}
  double GetErrorDown() const {return erDn;}
  double GetWeight() const {return weight;}
  
 private:

  double erUp, erDn, weight;

};

#endif
