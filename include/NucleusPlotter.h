#ifndef NucleusPlotter_h
#define NucleusPlotter_h 1

#include "Nucleus.h"
//#include "TH1D.h"

class NucleusPlotter {
  
 public:
  
  NucleusPlotter();
  NucleusPlotter(Nucleus* nuc);  
  ~NucleusPlotter();
  
  void PlotLifetime(int ni, int nbins);
  void PlotBranchingRatio(int ni, int nf1, int nf2, int nbins);
  void PlotMixingRatio(int ni, int nf, int nbins);
  void PlotMatrixElement(int n1, int n2, int mult, int nbins);

  std::string GetMinimum() const;
  Nucleus* GetNucleus() const {return nucleus;}

  void SetNucleus(Nucleus* nuc) {nucleus = nuc;}
  void SetFileName(std::string nm) {file_name = nm;}
  void SetFlag(bool flg) {flag = flg;}

 private:
  
  std::string MakeName(int ni, int nf1, int nf2, int type) const;
  std::string MakeTitle(int ni, int nf1, int nf2, int type) const;
  
  double GetValue(int ni, int nf1, int nf2, int type) const;
  void Plot(int ni, int nf1, int nf2, int nbins, int type);

  Nucleus* nucleus;
  std::string file_name;
  bool flag; //Flag for file_name to refer to a Parameter_Space.root file

};

#endif
