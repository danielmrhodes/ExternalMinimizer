#include "ExperimentalData.h"
#include "TMath.h"
#include "TRandom3.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

ExperimentalData::ExperimentalData(std::string nm) : name(nm), nrm_index(-1) {;}
ExperimentalData::~ExperimentalData() {

  for(Experiment* exp : data)
    delete exp;
  data.clear();
  
  return;
}

void ExperimentalData::ReadGosiaFiles() {
  ReadGosiaFiles(name+".POIN.inp",name+".POIN.out",name+".INTI.out");
  return;
}

void ExperimentalData::ReadGosiaFiles(std::string poin_inp, std::string poin_out, std::string inti_out) {

  CreateFromGosiaInput(poin_inp);
  if(data.size() == 0) {
    std::cout << "Did not create any experiments for " << name << "!" << std::endl;
    return;
  }
  
  FillFromIntiOutput(inti_out,true);
  FillFromPoinOutput(poin_out,true);

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
  ReadDataFile(name+".yld");
  return;
}

void ExperimentalData::ReadDataFile(std::string file_name) {

  std::ifstream file(file_name.c_str());
  if(!file.is_open()) {
    std::cout << "Could not open data file " << file_name << "!" << std::endl;
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
	
	int ni1 = indsI/100;
	int ni2 = indsI%100;
	int nf1 = indsF/100;
	int nf2 = indsF%100;

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

void ExperimentalData::GenerateData(Nucleus* nuc, TF1* f, double scale) {
  GenerateData(name+".sim",nuc,f,scale);
  return;
}

void ExperimentalData::GenerateData(std::string file_name, Nucleus* nuc, TF1* eff_curve, double scale) {
  
  std::ifstream file(file_name.c_str());
  if(!file.is_open()) {
    std::cout << "Could not open file " << file_name << "!" << std::endl;
    return;
  }

  for(Experiment* exp : data)
    exp->ClearRawData();

  TRandom3 rand;

  int index = 0;
  std::string line, word;
  while(std::getline(file,line)) {
    if(line.empty())
      break;

    Experiment* exp = data.at(index);

    int num_y;
    double weight;
    std::stringstream ss(line);
    ss >> num_y >> weight;

    for(int i=0;i<num_y;++i) {
      std::getline(file,line);
     
      int indsI, indsF;
      double yl, er;

      std::stringstream ss1(line);
      ss1 >> indsI >> indsF;

      YieldError* yldD = new YieldError();
      yldD->SetWeight(weight);      

      if(indsI < 100) {
       
	double egam = 1000.0*(nuc->GetLevelEnergy(indsI) - nuc->GetLevelEnergy(indsF));
	double eff = eff_curve->Eval(egam); 

	Yield* yldI = exp->GetIntiYield(indsI-1,indsF-1);
	double yl = scale*yldI->GetValue()*eff;
	
	double ydat = rand.PoissonD(yl);
	//double ydat = yl;
	double err = TMath::Sqrt(ydat);

	ydat /= eff;
	err /= eff;	

	yldD->AddInitialIndex(indsI-1);
        yldD->AddFinalIndex(indsF-1);
	yldD->SetValue(ydat);
	yldD->SetErrorUp(err);
	yldD->SetErrorDown(err);

      }
      else {

        int ni1 = indsI/100;
        int ni2 = indsI%100;
        int nf1 = indsF/100;
        int nf2 = indsF%100;

	double egam1 = 1000.0*(nuc->GetLevelEnergy(ni1) - nuc->GetLevelEnergy(nf1));
	double eff1 = eff_curve->Eval(egam1);
	
	Yield* yldI1 = exp->GetIntiYield(ni1-1,nf1-1);
        double y1 = scale*yldI1->GetValue()*eff1;        

	double ydat1 = rand.PoissonD(y1);
	//double ydat1 = y1;
	double err1 = TMath::Sqrt(ydat1);

	ydat1 /= eff1;
        err1 /= eff1;
	
	double egam2 = 1000.0*(nuc->GetLevelEnergy(ni2) - nuc->GetLevelEnergy(nf2));
	double eff2 = eff_curve->Eval(egam2);

	Yield* yldI2 = exp->GetIntiYield(ni2-1,nf2-1);
        double y2 = scale*yldI2->GetValue()*eff2;
	
	double ydat2 = rand.PoissonD(y2);
	//double ydat2 = y2;
	double err2 = TMath::Sqrt(ydat2);
	
	ydat2 /= eff2;
        err2 /= eff2;

	double val = ydat1 + ydat2;
	double err = TMath::Sqrt(err1*err1 + err2*err2);
	
        yldD->AddInitialIndex(ni1-1);
        yldD->AddFinalIndex(nf1-1);
        yldD->AddInitialIndex(ni2-1);
        yldD->AddFinalIndex(nf2-1);

	yldD->SetValue(val);
        yldD->SetErrorUp(err);
        yldD->SetErrorDown(err);

      }

      exp->AddRawData(yldD);
       
    } //End loop over one experiment's yields
    	
    ++index;
  } //End loop over all experiments

  return;
}

void ExperimentalData::GenerateAllData(Nucleus* nuc, TF1* eff_curve, double scale) {
 
  TRandom3 rand;
  for(Experiment* exp : data) {
    exp->ClearRawData();
  
    for(Yield* yldI : exp->GetAllIntiYields()) {

      std::vector<int> indsI = yldI->GetInitialIndices();
      std::vector<int> indsF = yldI->GetFinalIndices();
      int ni = indsI.at(0);
      int nf = indsF.at(0);

      double egam = 1000.0*(nuc->GetLevelEnergy(ni) - nuc->GetLevelEnergy(nf));
      double eff = eff_curve->Eval(egam);
      
      double yl = scale*yldI->GetValue()*eff;
      double ydat = rand.PoissonD(yl);
      double err = TMath::Sqrt(ydat);

      ydat /= eff;
      err /= eff;

      YieldError* yldD = new YieldError();
      yldD->AddInitialIndex(ni);
      yldD->AddFinalIndex(nf);
      yldD->SetValue(ydat);
      yldD->SetErrorUp(err);
      yldD->SetErrorDown(err);
      yldD->SetWeight(1.0);
    
      exp->AddRawData(yldD);
    } //End loop over inti yields  
  } //End loop over experiments
  
  return;
}

void ExperimentalData::WriteDataFile(int A, int Z) {
  WriteDataFile(name+".yld",A,Z);
  return;
}

void ExperimentalData::WriteDataFile(std::string file_name, int A, int Z) {
  
  std::ofstream file(file_name.c_str());
  if(!file.is_open()) {
    std::cout << "Could not open file " << file_name << "!" << std::endl;
    return;
  } 
  
  for(Experiment* exp : data) {
    
    int num = exp->GetNumber();
    int en = 0.5*(exp->GetEnergyMax() + exp->GetEnergyMin());
    int num_y = exp->GetRawData().size();

    file << num << " " << 1 << " " << Z << " " << A << " " << en << " " << num_y << " " << 1.0 << "\n";
    for(YieldError* yldR : exp->GetRawData()) {

      std::vector<int> indsI = yldR->GetInitialIndices();
      std::vector<int> indsF = yldR->GetFinalIndices();

      int ni = indsI.at(0);
      int nf = indsF.at(0);
      if(indsI.size() > 1) {

	ni = 100*(ni+1);
	nf = 100*(nf+1);
	ni += indsI.at(1);
	nf += indsF.at(1);
	
      }

      double val = yldR->GetValue();
      double err = 0.5*(yldR->GetErrorUp() + yldR->GetErrorDown());

      file << ni+1 << " " << nf+1 << " " << val << " " << err << "\n";

    }

  }

  return;
}

void ExperimentalData::CreateFromGosiaInput() {
  CreateFromGosiaInput(name+".POIN.inp");
  return;
}

void ExperimentalData::CreateFromGosiaInput(std::string file_name) {

  std::ifstream file(file_name.c_str());
  if(!file.is_open()) {
    std::cout << "Could not open Gosia file " << file_name << "!" << std::endl;
    return;
  }

  for(Experiment* exp : data)
    delete exp;
  data.clear();
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
  FillFromIntiOutput(name+".INTI.out",create);
  return;
}

void ExperimentalData::FillFromIntiOutput(std::string file_name, bool create) {
  
  std::ifstream file(file_name.c_str());
  if(!file.is_open()) {
    std::cout << "Could not open Gosia output file " << file_name << "!" 
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
  FillFromPoinOutput(name+".POIN.out",create);
  return;
}

void ExperimentalData::FillFromPoinOutput(std::string file_name, bool create) {

  std::ifstream file(file_name.c_str());
  if(!file.is_open()) {
    std::cout << "Could not open Gosia output file " << file_name << "!" << std::endl;
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
  if(nrm_index < 0)
    nrm_index = size-1;

  const double scale = data.at(0)->GetAllIntiYields().at(nrm_index)->GetValue() / data.at(0)->GetAllPointYields().at(nrm_index)->GetValue();
  
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
  if(nrm_index < 0)
    nrm_index = sizeT-1;

  double scale = data.at(0)->GetAllIntiYields().at(nrm_index)->GetValue() / data.at(0)->GetAllPointYields().at(nrm_index)->GetValue();

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

int ExperimentalData::GetNumRawData() {

  int size = 0;
  for(Experiment* exp : data)
    size += exp->GetRawDataSize();
  
  return size;

}
