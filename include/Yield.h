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

  int Size() const {return nis.size();}
  std::vector<int> GetInitialIndices() const {return nis;}
  std::vector<int> GetFinalIndices() const {return nfs;}
  double GetValue() const {return yld;}
    
 private:

  std::vector<int> nis, nfs;
  double yld;

};

class YieldError : public Yield {

 public:
  
  YieldError() : Yield() {;}
  ~YieldError() {;}

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
