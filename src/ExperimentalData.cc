#include "ExperimentalData.h"
#include "TMath.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

ExperimentalData::ExperimentalData(std::string nm) : name(nm) {;}
ExperimentalData::~ExperimentalData() {

  for(Experiment* exp : data)
    delete exp;

  return;
}

void ExperimentalData::ReadGosiaFiles() {

  CreateFromGosiaInput();
  if(data.size() == 0) {
    std::cout << "Did not create any experiments for " << name << "!" << std::endl;
    return;
  }
  
  FillFromIntiOutput(true);
  FillFromPoinOutput(true);

  int size0 = data.at(0)->GetTotalSize();
  for(Experiment* exp : data) {
    
    int size = exp->GetTotalSize();
    if(size != size0)
      std::cout << "Not all " << name << " experimets have the same size!" << std::endl;
    
    exp->ClearCorrectionFactors();
    exp->corr_facs.resize(size);
  }

  DeriveCorrectionFactors();

  return;
}

void ExperimentalData::ReadDataFile() {

  std::ifstream file((name+".yld").c_str());
  if(!file.is_open()) {
    std::cout << "Could not open data file " << name+".yld" << "!" << std::endl;
    return;
  }

  for(Experiment* exp : data)
    exp->ClearRawData();

  std::string line, word;
  int index = 0;
  while(std::getline(file,line)) {
    if(line.empty())
      break;

    Experiment* exp = data.at(index);

    int num_y;
    double weight;
    std::stringstream ss(line);
    ss >> num_y >> num_y >> num_y >> num_y >> weight >> num_y >> weight;

    for(int i=0;i<num_y;++i) {
      std::getline(file,line);
      
      int indsI, indsF;
      double yl, er;
      
      std::stringstream ss1(line);
      ss1 >> indsI >> indsF >> yl >> er;
      
      YieldError* yld = new YieldError();
      if(indsI < 100) {
	yld->AddInitialIndex(indsI-1);
	yld->AddFinalIndex(indsF-1);
      }
      else {
	
	int ni1 = indsI%100;
	int ni2 = indsI/100;
	int nf1 = indsF%100;
	int nf2 = indsF/100;

	yld->AddInitialIndex(ni1-1);
	yld->AddFinalIndex(nf1-1);
	yld->AddInitialIndex(ni2-1);
	yld->AddFinalIndex(nf2-1);
      
      }
      
      yld->SetValue(yl);
      yld->SetErrorUp(er);
      yld->SetErrorDown(er);
      yld->SetWeight(weight);
      
      exp->AddRawData(yld);

    }
    ++index;
  }

  return;
}

void ExperimentalData::CreateFromGosiaInput() {

  std::ifstream file((name+".POIN.inp").c_str());
  if(!file.is_open()) {
    std::cout << "Could not open Gosia file " << name+".POIN.inp" << "!" << std::endl;
    return;
  }

  for(Experiment* exp : data)
    delete exp;
  
  scalings.clear();
  factors.clear();
  
  std::string qry = "EXPT";
  std::string line;
  while(std::getline(file,line)) {
    if(line.compare(0,qry.size(),qry))
      continue;
    
    std::getline(file,line);
    break;
  }
  
  int num_exp, z1, a1;
  std::stringstream ss(line);
  ss >> num_exp >> z1 >> a1;
  
  for(int i=0;i<num_exp;++i) {
    std::getline(file,line);
    
    Experiment* exp = new Experiment();
    exp->SetNumber(i+1);

    int z2, a2;
    double en, th;

    std::stringstream ss1(line);
    ss1 >> z2 >> a2 >> en >> th;

    if(z2 < 0) {
      exp->beam_Z = z1;
      exp->beam_mass = a1;
      exp->targ_Z = z2;
      exp->targ_mass = a2;
    }
    else {
      exp->beam_Z = z2;
      exp->beam_mass = a2;
      exp->targ_Z = z1;
      exp->targ_mass = a1;
    }

    exp->Ep = en;
    exp->thP = th*TMath::DegToRad();
    
    data.push_back(exp);
    scalings.push_back(1.0);
    factors.push_back(1.0);
  
  }
  
  return;
}

void ExperimentalData::FillFromIntiOutput(bool create) {
  
  std::ifstream file((name+".INTI.out").c_str());
  if(!file.is_open()) {
    std::cout << "Could not open Gosia output file " << name+".INTI.out" << "!" 
	      << std::endl;
    return;
  }

   if(create) {
    for(Experiment* exp : data) {
      exp->ClearAllIntiYields();
    }
  }
  else {
    for(Experiment* exp : data) {
      exp->ZeroAllIntiYields();
    }
  }

  std::string qry = "     INTEGRATED RUTHERFORD CROSS SECTION=";

  std::string line, word;
  while(std::getline(file,line)) {
    if(line.compare(0,qry.size(),qry))
      continue;
   
    word = line.substr(qry.size(),9);
    
    double cs;
    std::stringstream ss(word);
    ss >> cs;

    word = line.substr(line.size()-2,2);
    
    int num;
    std::stringstream ss1(word);
    ss1 >> num;
    
    for(int i=0;i<9;i++)
      std::getline(file,line);

    word = line.substr(19,7);

    double emin, emax, tmin, tmax;
    
    std::stringstream ss2(word);
    ss2 >> emin;

    word = line.substr(30,7);
    std::stringstream ss3(word);
    ss3 >> emax;

    word = line.substr(67,7);
    std::stringstream ss4(word);
    ss4 >> tmin;

    word = line.substr(77,7);
    std::stringstream ss5(word);
    ss5 >> tmax;
    
    Experiment* exp = data.at(num-1);
    exp->SetRutherfordCS(cs);
    exp->SetEnergyMin(emin);
    exp->SetEnergyMax(emax);
    exp->SetThetaMin(tmin*TMath::DegToRad());
    exp->SetThetaMax(tmax*TMath::DegToRad());

    //Jump to yield printout for this experiment
    for(int i=0;i<3;i++)
      std::getline(file,line);

    int indexY = 0;
    while(std::getline(file,line)) { //Add all yields

      if(line.empty())
	break;
      
      if(line.find("********* END OF EXECUTION **********") != std::string::npos)
        break;

      int ni, nf;
      double ji, jf, yl;

      std::stringstream ss6(line);
      ss6 >> ni >> nf >> ji >> jf >> yl;
      
      if(create) {
	
	Yield* yld = new Yield();
	yld->AddInitialIndex(ni-1);
	yld->AddFinalIndex(nf-1);
	yld->SetValue(yl);
	
	exp->AddIntiYield_All(yld);
	
      }
      else {
	exp->UpdateIntiYield_All(indexY,yl);
      }
      ++indexY;
    }
  }
  
  return;
}

void ExperimentalData::FillFromPoinOutput(bool create) {

  std::ifstream file((name+".POIN.out").c_str());
  if(!file.is_open()) {
    std::cout << "Could not open Gosia output file " << name+".POIN.out" << "!" << std::endl;
    return;
  }

  if(create) {
    for(Experiment* exp : data) {
      exp->ClearAllPointYields();
      exp->ClearPointYields();
    }
  }
  else {
    for(Experiment* exp : data) {
      exp->ZeroAllPointYields();
      exp->ZeroPointYields();
    }
  }
  
  std::string qry = "                                                  CALCULATED YIELDS";
  std::string line, word;
  int indexE = 0;
  while(std::getline(file,line)) {
    if(line.compare(0,qry.size(),qry))
      continue;
    
    //Jump to yield printout for this experiment
    for(int i=0;i<6;++i)
      std::getline(file,line);

    Experiment* exp = data[indexE];
    int indexY = 0;
    while(std::getline(file,line)) { //Add all yields

      if(line.empty())
	break;
      
      if(line.find("********* END OF EXECUTION **********") != std::string::npos)
        break;

      int ni, nf;
      double ji, jf, yl;

      std::stringstream ss6(line);
      ss6 >> ni >> nf >> ji >> jf >> yl;
      
      if(create) {
	
	Yield* yld = new Yield();
	yld->AddInitialIndex(ni-1);
	yld->AddFinalIndex(nf-1);
	yld->SetValue(yl);
	
	exp->AddPointYield_All(yld);
      }
      else {
	
	exp->UpdatePointYield_All(indexY,yl);
	
	for(int i=0;i<exp->point_yields.size();++i) {

	  Yield* yldP = exp->point_yields[i];
	  double val = yldP->GetValue();
	  
	  std::vector<int> nis = yldP->GetInitialIndices();
	  std::vector<int> nfs = yldP->GetFinalIndices();
	  for(int j=0;j<nis.size();++j) {
	    if(nis[j] == ni-1 && nfs[j] == nf-1)
	      yldP->SetValue(val + yl);
	  }
	}
      }
      ++indexY;
    }
    ++indexE;
  }

  return;
}

void ExperimentalData::DeriveCorrectionFactors() {

  const int size = data.at(0)->GetTotalSize();
  const double scale = data.at(0)->GetAllIntiYields().at(size-1)->GetValue() / data.at(0)->GetAllPointYields().at(size-1)->GetValue();
  
  for(Experiment* exp : data) {
    
    std::vector<Yield*> yldsI = exp->GetAllIntiYields();
    std::vector<Yield*> yldsP = exp->GetAllPointYields();
    for(int i=0;i<size;++i) {
      
      double cc = scale*yldsP.at(i)->GetValue()/yldsI.at(i)->GetValue();
      exp->corr_facs.at(i) = cc;
      
    }
  }

  return;
}

void ExperimentalData::Correct() {

  const int sizeT = data.at(0)->GetTotalSize();
  double scale = data.at(0)->GetAllIntiYields().at(sizeT-1)->GetValue() / data.at(0)->GetAllPointYields().at(sizeT-1)->GetValue();

  for(Experiment* exp : data) {

    exp->ClearCorrectedData();
    exp->ClearPointYields();
    for(YieldError* yldR : exp->GetRawData()) {
      
      YieldError* yldC = new YieldError();
      Yield* yldP = new Yield();
      double valP = 0.0;
      
      std::vector<int> nis = yldR->GetInitialIndices();
      std::vector<int> nfs = yldR->GetFinalIndices();
      
      double corF = 0.0;
      double sumW = 0.0;

      int num = nis.size();
      for(int j=0;j<num;++j) {
	int ni = nis.at(j);
	int nf = nfs.at(j);

	yldC->AddInitialIndex(ni);
	yldC->AddFinalIndex(nf);
	yldP->AddInitialIndex(ni);
	yldP->AddFinalIndex(nf);
	
	Yield* ylp = exp->GetPointYield(ni,nf);
	Yield* yli = exp->GetIntiYield(ni,nf);
	
	double yp = ylp->GetValue();
	double yi = yli->GetValue();
	
	valP += yp;
	double cc = scale*yp/yi;

	corF += yi*cc;
	sumW += yi;
      }
      corF /= sumW;      
      
      yldC->SetValue(yldR->GetValue() * corF);
      yldC->SetErrorUp(yldR->GetErrorUp() * corF);
      yldC->SetErrorDown(yldR->GetErrorDown() * corF);
      yldC->SetWeight(yldR->GetWeight());
      exp->AddCorrData(yldC);
      
      yldP->SetValue(valP);
      exp->AddPointYield(yldP);
      
    }
  }

  return;
}

void ExperimentalData::RecorrectExp(int exp_num, double scale) {

  Experiment* exp = data.at(exp_num-1);
  for(int i=0;i<exp->GetCorrectedData().size();++i) {
    
    YieldError* yldC = exp->GetCorrectedData().at(i);
    std::vector<int> nis = yldC->GetInitialIndices();
    std::vector<int> nfs = yldC->GetFinalIndices();
      
    double corF = 0.0;
    double sumW = 0.0;

    int num = nis.size();
    for(int j=0;j<num;++j) {
      
      int ni = nis.at(j);
      int nf = nfs.at(j);
	
      Yield* ylp = exp->GetPointYield(ni,nf);
      Yield* yli = exp->GetIntiYield(ni,nf);
	
      double yp = ylp->GetValue();
      double yi = yli->GetValue();
      double cc = scale*yp/yi;
      
      corF += yi*cc;
      sumW += yi;
      
    }
    corF /= sumW;      
      
    YieldError* yldR = exp->GetRawData().at(i);  
    yldC->SetValue(yldR->GetValue() * corF);
    yldC->SetErrorUp(yldR->GetErrorUp() * corF);
    yldC->SetErrorDown(yldR->GetErrorDown() * corF);
      
  }

  const int size = exp->GetTotalSize();
  std::vector<Yield*> yldsI = exp->GetAllIntiYields();
  std::vector<Yield*> yldsP = exp->GetAllPointYields();
  for(int i=0;i<size;++i) {
    double cc = scale*yldsP.at(i)->GetValue()/yldsI.at(i)->GetValue();
    exp->corr_facs.at(i) = cc;
  }
  
  return;
}

void ExperimentalData::PrintComparison() const {

  for(int i=0;i<data.size();++i) {
    Experiment* exp = data.at(i);
    exp->PrintComparison(name,scalings.at(i),factors.at(i));
  }
  
  return;
}

void ExperimentalData::PrintAllYields() const {

  for(Experiment* exp : data)
    exp->PrintAllYields(name);

  return;
}

void ExperimentalData::PrintAllExpYields(int num) const {

  int size = data.size();
  if(num < 1 || num > size) {
    std::cout << "Experiment number must be between 1 and " << size << " for " << name << " data" 
	      << std::endl;
    return;
  }
  
  Experiment* exp = data.at(num-1);
  exp->PrintAllYields(name);

  return;
}

void ExperimentalData::PrintRawData() const {

  for(Experiment* exp : data)
    exp->PrintRawData(name);

  return;
}

void ExperimentalData::PrintRawAndCorrData() const {

  for(Experiment* exp : data)
    exp->PrintRawAndCorrData(name);

  return;
}

void ExperimentalData::PrintExpRawData(int num) const {

  int size = data.size();
  if(num < 1 || num > size) {
    std::cout << "Experiment number must be between 1 and " << size << " for " << name << " data" 
	      << std::endl;
    return;
  }
  
  Experiment* exp = data.at(num-1);
  exp->PrintRawData(name);

  return;
}

