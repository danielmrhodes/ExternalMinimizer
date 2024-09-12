#include "Literature.h"
#include "TMath.h"
#include <iostream>
#include <fstream>
#include <sstream>

LitVal::LitVal() {;}
LitVal::~LitVal() {;}
LitVal::LitVal(int ind, double val, double erU, double erD) {

  indices.push_back(ind);
  
  value = val;
  errUp = erU;
  errDn = erD;
    
    return;
}

LitVal::LitVal(int indI, int indF, double val, double erU, double erD) : LitVal(indI,val,erU,erD) {

  indices.push_back(indF);

}

LitVal::LitVal(int indI, int indF1, int indF2, double val, double erU, double erD) 
  : LitVal(indI,indF1,val,erU,erD) {

  indices.push_back(indF2);

}

double LitVal::Compare(double calc) const{
  
  double diff = value - calc;
  if(diff > 0.0)
    return TMath::Power(diff/GetErrorDown(),2.0);

  return TMath::Power(diff/GetErrorUp(),2.0);

}

double LitVal::NSigma(double calc) const{
  
  double diff = value - calc;
  if(diff > 0.0)
    return diff/GetErrorDown();

  return diff/GetErrorUp();

}

Literature::Literature(std::string nm) {
  
  name = nm;
  for(int i=0;i<4;++i)
    weights.push_back(1.0);
  
}

Literature::~Literature() {

  for(LitVal* lv : lifetimes)
    delete lv;
  lifetimes.clear();

  for(LitVal* lv : branching_ratios)
    delete lv;
  branching_ratios.clear();

  for(LitVal* lv : mixing_ratios)
    delete lv;
  mixing_ratios.clear();

  for(LitVal* lv : matrix_elements)
    delete lv;
  matrix_elements.clear();

  return;
}

void Literature::CreateFromFile() {
  
  std::ifstream inFile((name+".lit").c_str());
  if(!inFile.is_open()) {
    std::cout << "Could not open literature file " << name+".lit" << "!" << std::endl;
    return;
  }

  weights.clear();
  lifetimes.clear();
  branching_ratios.clear();
  mixing_ratios.clear();
  matrix_elements.clear();
  
  std::string line;
  int num;
  double weight;

  std::getline(inFile,line);
  std::stringstream ss1(line);
  ss1 >> num >> weight;
  weights.push_back(weight);
  
  for(int i=0;i<num;++i) {
    
    std::getline(inFile,line);
    int ni, nf1, nf2;
    double val, err;

    std::stringstream ss(line);
    ss >> ni >> nf1 >> nf2 >> nf2 >> val >> err;
    AddBranchingRatio(ni,nf1,nf2,val,err,err);

  }

  std::getline(inFile,line);
  std::stringstream ss2(line);
  ss2 >> num >> weight;
  weights.push_back(weight);

  for(int i=0;i<num;++i) {
    
    std::getline(inFile,line);
    int ni;
    double val, err;

    std::stringstream ss(line);
    ss >> ni >> val >> err;
    AddLifetime(ni,val,err,err);

  }

  std::getline(inFile,line);
  std::stringstream ss3(line);
  ss3 >> num >> weight;
  weights.push_back(weight);
  
  for(int i=0;i<num;++i) {
    
    std::getline(inFile,line);
    int ni, nf;
    double val, err;
    
    std::stringstream ss(line);
    ss >> ni >> nf >> val >> err;
    AddMixingRatio(ni,nf,val,err,err);
    
  }

  std::getline(inFile,line);
  std::stringstream ss4(line);
  ss4 >> num >> weight;
  weights.push_back(weight);
  
  for(int i=0;i<num;++i) {
    
    std::getline(inFile,line);
    int mult, ni, nf;
    double val, err;
    
    std::stringstream ss(line);
    ss >> mult >> ni >> nf >> val >> err;
    AddMatrixElement(ni,nf,mult,val,err,err);
    
  }

  return;
}

void Literature::AddLifetime(int index, double val, double erU, double erD) {

  lifetimes.push_back(new LitVal(index-1,val,erU,erD));

  return;
}

void Literature::AddBranchingRatio(int init, int fin1, int fin2, double val, double erU, double erD) {
  
  branching_ratios.push_back(new LitVal(init-1,fin1-1,fin2-1,val,erU,erD));
  
  return;
}


void Literature::AddMixingRatio(int init, int fin, double val, double erU, double erD) {
  
  mixing_ratios.push_back(new LitVal(init-1,fin-1,val,erU,erD));

  return;
}

void Literature::AddMatrixElement(int ind1, int ind2, int mult, double val, double erU, double erD) {

  matrix_elements.push_back(new LitVal(ind1-1,ind2-1,mult,val,erU,erD));
  return;
}

int Literature::Size() const {
  
  return lifetimes.size() + branching_ratios.size() + mixing_ratios.size() + matrix_elements.size();

}
