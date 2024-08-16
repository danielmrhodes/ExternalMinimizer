#ifndef Level_h
#define Level_h 1

#include <string>

class Level {
  
 public:
  
  Level();
  ~Level();

  int GetIndex() const {return index;}
  int GetParity() const {return parity;}
  double GetSpin() const {return spin;}
  double GetEnergy() const {return energy;}

  void SetIndex(int ind) {index = ind;}
  void SetParity(int par) {parity =par;}
  void SetSpin(int sp) {spin = sp;}
  void SetEnergy(double en) {energy = en;}

  std::string GetParityS() const {
    if(parity == 1)
      return "+";
    return "-";
  }
  
 private:
  
  int index, parity;
  double spin, energy;
  
};

#endif
