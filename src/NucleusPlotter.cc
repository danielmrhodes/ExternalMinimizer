#include "NucleusPlotter.h"
#include "TH1D.h"
#include "TString.h"

#include <iostream>
#include <fstream>
#include <sstream>

NucleusPlotter::NucleusPlotter() {
  
  nucleus = NULL;
  file_name = "";
  flag = false;
  
  return;
}

NucleusPlotter::NucleusPlotter(Nucleus* nuc) : NucleusPlotter() {
  
  nucleus = nuc;
  return;
}

NucleusPlotter::~NucleusPlotter() {;}

void NucleusPlotter::Plot(int ni, int nf1, int nf2, int nbins, int type) {

  if(!nucleus) {
    std::cout << "Nucleus not set!" << std::endl;
    return;
  }

  std::string minimum = GetMinimum();
  if(minimum == "")
    return;
  
  double min_ch2;
  std::stringstream ss(minimum);
  ss >> min_ch2;
  
  int num = nucleus->GetNumMatrixElements();
  for(int i=0;i<num;++i) {
    
    double val;
    ss >> val;

    nucleus->GetMatrixElement(i)->SetValue(val);

  }

  double best = GetValue(ni,nf1,nf2,type);
  std::vector<double> vals;
  std::vector<double> vals_1s;

  std::ifstream file(file_name.c_str(),std::ios::in);
  
  std::string line, word;
  while(std::getline(file,line)) {
    
    double ch2;
    std::stringstream ss1(line);
    ss1 >> ch2;
    
    for(int i=0;i<num;++i) {
      
      double val;
      ss1 >> val;
      
      nucleus->GetMatrixElement(i)->SetValue(val);
      
    }

    double val = GetValue(ni,nf1,nf2,type);
    vals.push_back(val);

    if(ch2 < min_ch2 + 1.0)
      vals_1s.push_back(val);

  }

  double min = *std::min_element(vals.begin(),vals.end());
  double max = *std::max_element(vals.begin(),vals.end());
  double eU = max - best;
  double eD = best - min;
  
  std::cout << best << " +" << eU << " -" << eD << std::endl;;

  std::string name = MakeName(ni,nf1,nf2,type);
  std::string title = MakeTitle(ni,nf1,nf2,type);

  double step = (max-min)/10.0;
  TH1D* h = new TH1D(name.c_str(),title.c_str(),nbins,min-step,max+step);
  TH1D* h_1s = new TH1D((name+"_1s").c_str(),title.c_str(),nbins,min-step,max+step);
  h_1s->SetLineColor(kRed);
  
  for(double vl : vals)
    h->Fill(vl);

  for(double vl : vals_1s)
    h_1s->Fill(vl);
  
  h->Draw("hist");
  h_1s->Draw("hist same");

  return;

}

std::string NucleusPlotter::MakeName(int ni, int nf1, int nf2, int type) const {

  switch(type) {
    
  case 0:
    return Form("hTau_%d",ni);
    
  case 1:
    return Form("hBR_%d_%d_%d",ni,nf1,nf2);

  case 2:
    return Form("hDel_%d_%d",ni,nf1);

  case 3:
    return Form("hMR_%d_%d_%d",ni,nf1,nf2);
    
  }
  
  return "";
}

std::string NucleusPlotter::MakeTitle(int ni, int nf1, int nf2, int type) const {
  
  switch(type) {
    
  case 0:
    return Form("State %d Lifetime Distribution",ni);
    
  case 1:
    return Form("%d->%d / %d->%d Branching Ratio Distribution",ni,nf1,ni,nf2);

  case 2:
    return Form("%d->%d Transition Mixing Ratio Distribution",ni,nf1);

  case 3:
    return Form("<%d||%d||%d> Matrix Element Distribution",ni,nf2,nf1);
    
  }

  return "";
}

double NucleusPlotter::GetValue(int ni, int nf1, int nf2, int type) const {

  switch(type) {
    
  case 0:
    return nucleus->CalculateLifetime(ni);
    
  case 1:
    return nucleus->CalculateBranchingRatio(ni,nf1,nf2);

  case 2:
    return nucleus->CalculateMixingRatio(ni,nf1);

  case 3:
    return nucleus->GetMatrixElement(ni-1,nf1-2,nf2)->GetValue();
    
  }

  return 0;
  
}

std::string NucleusPlotter::GetMinimum() const {

  if(file_name == "") {
    std::cout << "Chi2 file name not set!" << std::endl;
    return "";
  }
  
  std::ifstream file(file_name.c_str(),std::ios::in);
  if(!(file.is_open())) {
    std::cout << "FAILED TO OPEN FILE " << file_name << std::endl;
    return "";
  }
  
  double min = 1E10;
  std::string line, word, best_line;
  while(std::getline(file,line)) {
    
    double ch2;
    std::stringstream ss(line);
    ss >> ch2;

    if(ch2 < min) {
      min = ch2;
      best_line = line;
    }

  }

  return best_line;

}
