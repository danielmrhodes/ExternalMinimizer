#include "Nucleus.h"
#include "TMath.h"
#include "TString.h"
#include "Math/SpecFuncMathMore.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

Nucleus::Nucleus(std::string nm) {
  
  name = nm;
  
  std::vector<double> ens = {0.01,0.015,1.0,5.5,6.0};
  std::vector<double> ccs = {0.0,0.0,0.0,0.0,0.0};

  for(int i=0;i<8;++i) {
    convCoeffs[i] = new TGraph(5,&ens[0],&ccs[0]);
    convCoeffs[i]->SetName(Form("CC%d",i+1));
    convCoeffs[i]->SetTitle(Form("CC%d",i+1));
  }

  return;

}

Nucleus::~Nucleus() {

  for(MatrixElement* me : matrix_elements)
    delete me;
  matrix_elements.clear();

  for(Level* lvl : levels)
    delete lvl;
  levels.clear();

  for(TGraph* g : convCoeffs)
    delete g;
  
  return;

}

double Nucleus::DecayLifetime(int mult, double me, double spin, double egam) const {
  
  double bsl = me*me/(2.0*spin + 1.0);
  switch(mult) {

  case 1:
    
    //return 0.62757/(TMath::Power(egam,3.0)*bsl*100000.0);
    return 0.629/(TMath::Power(egam,3.0)*bsl*100000.0);

  case 2:

    //return 813.68/(TMath::Power(egam,5.0)*bsl*10000.0);
    return 816.0/(TMath::Power(egam,5.0)*bsl*10000.0);
    
  case 3:

    //return 1748.55/(TMath::Power(egam,7.0)*bsl);
    return 1760.0/(TMath::Power(egam,7.0)*bsl);
    
  case 7:

    //return 56.842/(TMath::Power(egam,3.0)*bsl*1000.0);
    return 56.8/(TMath::Power(egam,3.0)*bsl*1000.0);
    
  case 8:
    
    //return 738.66/(TMath::Power(egam,5.0)*bsl);
    return 738.1/(TMath::Power(egam,5.0)*bsl);
    
  }
  
  std::cout << "Unsupported Multipolarity " << mult << std::endl;
  return 0.0;
  
}

double Nucleus::CalculateLifetime(int index) const {
  
  double sum = 0.0;
  for(MatrixElement* me : matrix_elements) {
    
    int n1 = me->GetIndex1();
    int n2 = me->GetIndex2();
    if(n1 == n2 || (n1 != index && n2 != index))
      continue;
    
    double egam = levels[n2]->GetEnergy() - levels[n1]->GetEnergy();
    double val = me->GetValue();
    int mult = me->GetMultipolarity();
    double cc = convCoeffs[mult-1]->Eval(std::abs(egam));

    if(n2 == index && egam > 0.0)
      sum += (cc + 1.0)/DecayLifetime(mult,val,levels[n2]->GetSpin(),egam);
    else if(n1 == index && egam < 0.0)
      sum += (cc + 1.0)/DecayLifetime(mult,val,levels[n1]->GetSpin(),-egam);

  }
  
  return 1.0/sum;
}


double Nucleus::CalculateBranchProbability(int ni, int nf) const {

  double egam = levels[ni]->GetEnergy() - levels[nf]->GetEnergy();
  if(egam < 0.0)
    return 0.0;
  
  double sum = 0.0;
  for(MatrixElement* me : matrix_elements) {
    
    int n1 = me->GetIndex1();
    int n2 = me->GetIndex2();
    if((ni == n1 && nf == n2) || (ni == n2 && nf == n1)) {

      int mult = me->GetMultipolarity();
      double cc = convCoeffs[mult-1]->Eval(egam);
      
      sum += (cc + 1.0)/DecayLifetime(mult,me->GetValue(),levels[ni]->GetSpin(),egam);

    }
  
  }

  return sum*CalculateLifetime(ni);

}

double Nucleus::CalculateBranchingRatio(int ni, int nf1, int nf2) const {

  double egam1 = levels[ni]->GetEnergy() - levels[nf1]->GetEnergy();
  if(egam1 < 0.0)
    return 0.0;

  double egam2 = levels[ni]->GetEnergy() - levels[nf2]->GetEnergy();
  if(egam2 < 0.0)
    return 0.0;
  
  double spin = levels[ni]->GetSpin();
  double sum1 = 0.0;
  double sum2 = 0.0;
  
  for(MatrixElement* me : matrix_elements) {

    int n1 = me->GetIndex2();
    int n2 = me->GetIndex1();
    int mult = me->GetMultipolarity();
    
    if((ni == n1 && nf1 == n2) || (ni == n2 && nf1 == n1)) {

      
      double cc = convCoeffs[mult-1]->Eval(egam1);
      sum1 += (cc + 1.0)/DecayLifetime(mult,me->GetValue(),spin,egam1);
      
    }
    else if((ni == n1 && nf2 == n2) || (ni == n2 && nf2 == n1)) {
      
      double cc = convCoeffs[mult-1]->Eval(egam2);
      sum2 += (cc + 1.0)/DecayLifetime(mult,me->GetValue(),spin,egam2);
      
    }
  }

  return sum1/sum2;

}

double Nucleus::CalculateMixingRatio(int ni, int nf) const {

  //Only E2/M1 supported

  double egam = levels[ni]->GetEnergy() - levels[nf]->GetEnergy();
  if(egam < 0.0)
    return 0.0;
  
  double valE2 = 0.0;
  double valM1 = 0.0;
  for(MatrixElement* me : matrix_elements) {

    int n1 = me->GetIndex1();
    int n2 = me->GetIndex2();
    if((ni == n1 && nf == n2) || (ni == n2 && nf == n1)) {

      int mult = me->GetMultipolarity();
      if(mult == 2)
	valE2 = me->GetValue();
      else if(mult == 7)
	valM1 = me->GetValue();
    }
    
  }

  if(std::abs(valE2) < 1E-18 || std::abs(valM1) < 1E-18)
    return 0.0;

  return 0.835*egam*valE2/valM1;

}

std::vector<int> Nucleus::GetTransitionMultipolarities(int ni, int nf) const {

  std::vector<int> lvals;
  if(ni == nf)
    return lvals;
  
  if(levels[ni]->GetEnergy() < levels[nf]->GetEnergy())
    return lvals;

  for(MatrixElement* me : matrix_elements) {

    int n1 = me->GetIndex1();
    int n2 = me->GetIndex2();
    if((ni == n1 && nf == n2) || (ni == n2 && nf == n1))
      lvals.push_back(me->GetMultipolarity());
    
  }

  return lvals;

}

double Nucleus::CompareLifetime(const LitVal* lv) const {
  
  double tau = CalculateLifetime(lv->GetInitialIndex());
  return lv->Compare(tau);

}

double Nucleus::CompareBranchingRatio(const LitVal* lv) const {
  
  double br = CalculateBranchingRatio(lv->GetInitialIndex(),lv->GetFinalIndex(),
				       lv->GetFinalIndexTwo()); 
  return lv->Compare(br);

}

double Nucleus::CompareMixingRatio(const LitVal* lv) const {
  
  double mx = CalculateMixingRatio(lv->GetInitialIndex(),lv->GetFinalIndex()); 
  return lv->Compare(mx);

}

double Nucleus::CompareMatrixElement(const LitVal* lv) const {

  double me = GetMatrixElement(lv->GetInitialIndex(),lv->GetFinalIndex(),
			       lv->GetFinalIndexTwo())->GetValue();
  return lv->Compare(me);
}

double Nucleus::CompareWithLiterature(const Literature* lit) const {
  
  std::vector<double> ws = lit->GetWeights();
  double chi2 = 0.0;
  
  for(LitVal* lv : lit->GetLifetimes())
    chi2 += ws[0]*CompareLifetime(lv);
  
  for(LitVal* lv : lit->GetBranchingRatios())
    chi2 += ws[1]*CompareBranchingRatio(lv);

  for(LitVal* lv : lit->GetMixingRatios())
    chi2 += ws[2]*CompareMixingRatio(lv);

  for(LitVal* lv : lit->GetMatrixElements())
    chi2 += ws[3]*CompareMatrixElement(lv);
  
  return chi2;
}

/*
void Nucleus::CreateTheMap() {

  the_map.clear();
  for(int i=0;i<matrix_elements.size();++i) {
    
    MatrixElement* me = matrix_elements.at(i);
    int n1 = me->GetIndex1();
    int n2 = me->GetIndex2();
    int mult = me->GetMultipolarity();
    
    the_map[{mult,n1,n2}] = i;
    the_map[{mult,n2,n1}] = i;
  }

}
*/

MatrixElement* Nucleus::GetMatrixElement(int n1, int n2, int mult) const {
  
  //int index = the_map[{mult,n1,n2}];
  //return matrix_elements[index];

  for(MatrixElement* me : matrix_elements) {
    if(me->GetMultipolarity() != mult)
      continue;
    
    int ni = me->GetIndex1();
    int nf = me->GetIndex2();
    
    if((ni == n1 && nf == n2) || (ni == n2 && nf == n1))
      return me;
    
  }
  
  std::cout << name << " does not have <" << n1 << "||" << mult << "||" << n2 << "> matrix element!" << std::endl;  

  return NULL;
  
}

void Nucleus::ReleaseMatrixElement(int n1, int n2, int mult) {

  MatrixElement* me = GetMatrixElement(n1-1,n2-1,mult);
  me->Release();
  
  if(me->GetUpperLimit() - me->GetLowerLimit() < 1E-10)
    me->SetLowerLimit(-1.0*me->GetUpperLimit());

  return;
}

void Nucleus::FixMatrixElement(int n1, int n2, int mult) {

  MatrixElement* me = GetMatrixElement(n1-1,n2-1,mult);
  me->Fix();

  return;
}

void Nucleus::ReleaseAll() {

  for(MatrixElement* me : matrix_elements) {
    
    me->Release();
    if(me->GetUpperLimit() - me->GetLowerLimit() < 1E-10)
      me->SetLowerLimit(-1.0*me->GetUpperLimit());
  
  }
  
  return;
}

void Nucleus::FixAll() {
  
  for(MatrixElement* me : matrix_elements)  
    me->Fix();
  
  return;
}

/*
void Nucleus::Sort() {

  std::sort(matrix_elements.begin(),matrix_elements.end(),MatrixElement::Compare);
  //CreateTheMap();
  //for(int i=0;i<matrix_elements.size();++i)
  //matrix_elements[i]->SetVectorIndex(i);

  return;

}
*/

int Nucleus::GetNumFree() const {

  int num = 0;
  for(MatrixElement* me : matrix_elements)
    if(!(me->IsFixed()))
      num++;
  
  return num;
}

void Nucleus::CreateFromGosiaInputFile(std::string file_name) {
  
  levels.clear();
  matrix_elements.clear();
  
  //std::ifstream inFile((name+".POIN.inp").c_str());
  std::ifstream inFile(file_name.c_str());
  if(!inFile.is_open()) {
    std::cout << "Could not open Gosia input file " << name+".POIN.inp" << "!" << std::endl;
    return;
  }

  bool gosi = true;
  std::string line;
  
  while(std::getline(inFile,line)) {
    if(line.find("OP,COUL") != std::string::npos)
      gosi = false;
    if(line.find("LEVE") != std::string::npos)
      break;
  }
  
  while(std::getline(inFile,line)) {
    if(line.find("0 0 0 0") != std::string::npos)
      break;

    int index, par;
    double energy, spin;

    std::stringstream ss(line);
    ss >> index >> par >> spin >> energy;
    
    Level* lvl = new Level();
    lvl->SetIndex(index-1);
    lvl->SetParity(par);
    lvl->SetSpin(spin);
    lvl->SetEnergy(energy);
    
    levels.push_back(lvl);
    
  }
  
  std::string flag, end;
  if(gosi) {
    flag = " 0 0 0 0";
    end = "0 0 0 0 0";
  }
  else {
    flag = " 0 0";
    end = "0 0 0";
  }

  //int index = 0;
  int mult = 0;
  std::getline(inFile,line); //ME
  while(std::getline(inFile,line)) {
    if(line.find(end) != std::string::npos)
      break;
    
    if(line.find(flag) != std::string::npos) {
      std::stringstream ss(line);
      ss >> mult;
      continue;
    }
    
    int n1, n2;
    double val, llim, ulim;
    
    std::stringstream ss(line);
    ss >> n1 >> n2 >> val;
    if(gosi)
      ss >> llim >> ulim;
    
    MatrixElement* me = new MatrixElement();
    //me->SetVectorIndex(index);
    me->SetIndex1(n1-1);
    me->SetIndex2(n2-1);
    me->SetMultipolarity(mult);
    me->SetValue(val);

    if(gosi) {
      me->SetLowerLimit(llim);
      me->SetUpperLimit(ulim);
    }
    else {
      me->SetLowerLimit(val - 3.0);
      me->SetUpperLimit(val + 3.0);
    }
    
    matrix_elements.push_back(me);
    //++index;
  }

  //CreateTheMap();

  return;
}

void Nucleus::CreateFromGosiaOutputFile(std::string file_name) {

  levels.clear();
  matrix_elements.clear();
  
  std::ifstream inFile(file_name);
  if(!inFile.is_open()) {
    std::cout << "Could not open Gosia output file " << name + ".out" << "!" << std::endl;
    return;
  }

  std::string qry = "                                        LEVELS";

  std::string line;
  while(std::getline(inFile,line))
    if(line.compare(qry) == 0)
      break;

  std::getline(inFile,line);
  std::getline(inFile,line);
  while(std::getline(inFile,line)) {
    if(line.empty())
      break;
    
    int index;
    std::string parity;
    double energy, spin;

    std::stringstream ss(line);
    ss >> index >> parity >> spin >> energy;

    int par = 1;
    if(parity.compare("-")==0)
      par=-1;
    
    Level* lvl = new Level();
    lvl->SetIndex(index-1);
    lvl->SetParity(par);
    lvl->SetSpin(spin);
    lvl->SetEnergy(energy);
    
    levels.push_back(lvl);
    
  }

  int count = 0;
  const int nBlock = 1;
  qry = "                                        MATRIX ELEMENTS";
  while(std::getline(inFile,line)) { 
    if(line.compare(qry) != 0)
      continue;
    
    if(count < nBlock) {
      count++;
      continue;
    }
    
    std::getline(inFile,line);
    std::getline(inFile,line);

    int mult;
    while(std::getline(inFile,line)) {
      if(line.find("MULTIPOLARITY") != std::string::npos) {
	std::string multS(1,line.back());
	std::stringstream ss(multS);
	ss >> mult;
	std::getline(inFile,line);
	continue;
      }

      if(line.empty() ||
	 line.find("END OF EXECUTION") != std::string::npos ||
	 line.find("CONVERGENCE ACHIEVED") != std::string::npos)
	return;

      int index, n1, n2;
      double val;

      std::stringstream ss1(line);
      ss1 >> index >> n1 >> n2 >> val;

      MatrixElement* me = new MatrixElement();
      //me->SetVectorIndex(index-1);
      me->SetIndex1(n1-1);
      me->SetIndex2(n2-1);
      me->SetMultipolarity(mult);
      me->SetValue(val);
        
      matrix_elements.push_back(me);
    }
    
    //CreateTheMap();
    return;
  }
  
  //CreateTheMap();
  return;
}

void Nucleus::FillFromBSTFile(std::string file_name) {
  
  std::ifstream inFile(file_name);
  if(!inFile.is_open()) {
    std::cout << "Could not open Gosia file " << name+".bst" << "!" << std::endl;
    return;
  }
  
  std::string line;
  for(int i=0;i<matrix_elements.size();++i) {
    std::getline(inFile,line);

    double val;
    std::stringstream ss(line);
    ss >> val;

    MatrixElement* me = matrix_elements[i];
    me->SetValue(val);

    double ul = me->GetUpperLimit();
    double ll = me->GetLowerLimit();
    double relU = std::abs((val - ul)/val);
    double relL = std::abs((val - ll)/val);

    if(val > ul || val < ll || relU < 0.01 || relL < 0.01)  {
      me->SetLowerLimit(val - 3.0);
      me->SetUpperLimit(val + 3.0);
    }
    
  }

  return;
}

void Nucleus::FillFromBSTFile() {
  FillFromBSTFile(name + ".bst");
  return;
}

void Nucleus::CreateFromGosiaInputFile() {
  CreateFromGosiaInputFile(name + ".POIN.inp");
  return;
}

void Nucleus::CreateFromGosiaOutputFile() {
  CreateFromGosiaOutputFile(name + ".POIN.out");
  return;
}

void Nucleus::SetConverionCoefficients(int mult, std::vector<double> ens, std::vector<double> ccs) {

  int size = ens.size();
  if(size < 5) {
    std::cout << "Give at least 5 points for conversion coefficients (Mult=" << mult << ")" << std::endl;
    return;
  }

  if(ccs.size() != size) {
    std::cout << "Must give same number of energies (" << size << ") as values (" 
	      << ccs.size() << ") when setting conversion coefficients (Mult=" << mult << ")" 
	      << std::endl;
    return;
  }

  TGraph* g = convCoeffs[mult-1];
  delete g;
  g = NULL;

  convCoeffs[mult-1] = new TGraph(size,&ens[0],&ccs[0]);
  convCoeffs[mult-1]->SetName(Form("CC%d",mult));
  convCoeffs[mult-1]->SetTitle(Form("CC%d",mult));

  return;
}

std::vector<float> Nucleus::GetMatrixElementValues() const {

  std::vector<float> vals;
  
  int size = matrix_elements.size();
  vals.resize(size);
  
  for(int i=0;i<size;++i)
    vals[i] = float(matrix_elements[i]->GetValue());

  return vals;
}

void Nucleus::PrintLevels() const {
  
  int intw = 7;
  int valw = 12;
  
  std::cout << std::left;
  std::cout << "\nLevels:\n";
  std::cout << std::setw(intw) << "Index" 
	    << std::setw(intw) << "Jpi" 
	    << "Energy (MeV)\n";
    
  for(int i=0;i<levels.size();i++) {

    Level* lvl = levels.at(i);

    int index = lvl->GetIndex();
    double en = lvl->GetEnergy();
    double sp = lvl->GetSpin();
    std::string pr = lvl->GetParityS();

    std::string jpi = Form("%.1f",sp) + pr;
    if(std::floor(sp) == sp)
      jpi = Form("%.0f",sp) + pr;
    
    std::cout << std::setw(intw) << index+1
	      << std::setw(intw) << jpi
	      << std::setw(valw) << en << "\n"; 
  }
  
  return;
  
}

void Nucleus::PrintLifetimes() const {

  std::cout << "\nLifetimes:\n";
  for(int i=1;i<levels.size();i++) {
    
    int index = levels.at(i)->GetIndex();
    std::cout << "State " << index+1 << ": " << CalculateLifetime(index) << " ps\n";

  }

  return;
}

void Nucleus::PrintBranchProbabilities() const {
  
  int intw = 7;
  std::cout << std::left;
  std::cout << "\nBranch Probabilities:\n";
  std::cout << std::setw(intw) << "ni" 
	    << std::setw(intw) << "nf" 
	    << "prob\n";

  int size = levels.size();
  for(int i=0;i<size;++i) {
    
    int ni = levels[size-i-1]->GetIndex();
    for(int j=0;j<size;++j) {
      if(i==j)
	continue;
      
      int nf = levels[size-j-1]->GetIndex();
      double prob = CalculateBranchProbability(ni,nf);
      if(prob < 1E-18 || std::isnan(prob))
	continue;
      
      std::cout << std::setw(intw) << ni+1  
		<< std::setw(intw) << nf+1 
		<< prob << "\n";
    }
  }
  
  return;

}

void Nucleus::PrintBranchingRatios() const {
  
  int intw = 7;
  int valw = 12;
  
  std::cout << std::left;
  std::cout << "\nBranching Ratios:\n";
  std::cout << std::setw(intw) << "ni" 
	    << std::setw(intw) << "nf1" 
	    << std::setw(intw) << "nf2" 
	    << "ratio\n";
  
  int size = levels.size();
  for(int i=0;i<size;++i) {
    
    std::vector<int> nfs;
    std::vector<double> probs;
    double max = 0.0;
    int nfm = 0;
    
    int ni = levels.at(size-i-1)->GetIndex();
    for(int j=0;j<size;++j) {
      if(i==j)
	continue;
      
      int nf = levels.at(size-j-1)->GetIndex();
      double prob = CalculateBranchProbability(ni,nf);
      if(prob < 1E-18)
	continue;

      nfs.push_back(nf);
      probs.push_back(prob);

      if(prob > max) {
	max = prob;
	nfm = nf;
      }
      
    }
    
    for(int j=0;j<probs.size();++j)
      if(nfm != nfs[j])
	std::cout << std::setw(intw) << ni+1 
		  << std::setw(intw) << nfs[j]+1 
		  << std::setw(intw) << nfm+1
		  << probs[j]/max << "\n";
  }
  
  return;

}

void Nucleus::PrintMixingRatios() const {

  int intw = 7;
  std::cout << std::left;
  std::cout << "\nMixing Ratios:\n";
  std::cout << std::setw(intw) << "ni" 
	    << std::setw(intw) << "nf" 
	    << "value\n";

  int size = levels.size();
  for(int i=0;i<size;++i) {
    
    int ni = levels[size-i-1]->GetIndex();
    for(int j=i+1;j<size;++j) {
      
      int nf = levels[size-j-1]->GetIndex();
      double del = CalculateMixingRatio(ni,nf);
      if(std::abs(del) < 1E-18)
	continue;
      
      std::cout << std::setw(intw) << ni+1
		<< std::setw(intw) << nf+1
		<< del << "\n";
    }
  }
  
  return;

}

double Nucleus::Moment(double dme, double spin, int mult) {
  
  int isp2 = int(2.0*spin + 0.001);
  double w3j = ROOT::Math::wigner_3j(isp2,2*mult,isp2,-isp2,0,isp2);

  return TMath::Sqrt(16.0*TMath::Pi()/5.0)*w3j*dme;
}

void Nucleus::PrintMatrixElements() const {

  int intw = 7;
  int valw = 12;

  std::cout << std::left;
  std::cout << "\nMatrix Elements:\n";
  std::cout << std::setw(intw) << "n1" 
	    << std::setw(intw) << "n2" 
	    << std::setw(intw) << "mult" 
	    << std::setw(valw) << "value"
	    << "BSL/Q\n";
  
  for(int i=0;i<matrix_elements.size();i++) {
    
    MatrixElement* me = matrix_elements.at(i);
    
    int n1 = me->GetIndex1();
    int n2 = me->GetIndex2();
    int l = me->GetMultipolarity();
    if(l == 7)
      l=1;
    if(l == 8)
      l=2;

    double val = me->GetValue();
    std::string mult = me->GetMultS();

    double en1 = levels[n1]->GetEnergy();
    double en2 = levels[n2]->GetEnergy();
    
    double spin;
    if(en1 - en2 > 0 || n1 == n2)
      spin = levels[n1]->GetSpin();
    else
      spin = levels[n2]->GetSpin();
    
    double bsl_q;
    if(n1 != n2)
      bsl_q = val*val/(2.0*spin + 1.0);
    else
      bsl_q = Moment(val,spin,l);
      
    std::cout << std::setw(intw) << n1+1
	      << std::setw(intw) << n2+1
	      << std::setw(intw) << mult
	      << std::setw(valw) << val 
	      << bsl_q << "\n";

  }

  return;
}

void Nucleus::PrintGammaTransitions() const {

  int intw = 7;
  int valw = 12;

  std::cout << std::left;
  std::cout << "\nGamma-Ray Transitions:\n";
  std::cout << std::setw(intw) << "ni" 
	    << std::setw(intw) << "nf" 
	    << std::setw(valw) << "energy" 
	    << std::setw(valw) << "mult" 
	    << std::setw(valw) << "delta"
	    << "cc\n";
  
  int size = levels.size();
  for(int i=0;i<size;++i) {
    
    int ni = levels[size-i-1]->GetIndex();
    double tau = CalculateLifetime(ni);
    double spin = levels[size-i-1]->GetSpin();
    double eni = levels[size-i-1]->GetEnergy();
    
    for(int j=0;j<size;++j) {
      if(i==j)
	continue;

      int nf = levels[size-j-1]->GetIndex();
      double enf = levels[size-j-1]->GetEnergy();
      
      double egam = eni - enf;
      if(egam < 0.0)
	continue;
      
      double prob = CalculateBranchProbability(ni,nf);
      if(prob < 1E-18)
	continue;
      
      std::vector<int> lvals = GetTransitionMultipolarities(ni,nf);
      if(lvals.size() == 0)
	continue;

      std::vector<double> taups;
      std::string mults;
      int sign = 1;
      
      for(int l : lvals) {
	
	MatrixElement* me = GetMatrixElement(ni,nf,l);
	mults += me->GetMultS() + " ";

	double val = me->GetValue();
	if(val < 0.)
	  sign *= -1;

	double taup = DecayLifetime(l,val,spin,egam);
	taups.push_back(taup);

      }
      
      double cc = convCoeffs[lvals.at(0)-1]->Eval(egam);
      double del = 0.0;
      
      if(lvals.size() == 2) {
	
	int mult0 = lvals.at(0);
	int mult1 = lvals.at(1);

	for(int& l : lvals) {
	  if(l == 7)
	    l=1;
	  if(l == 8)
	    l=2;
	}
	
	if(lvals.at(0) > lvals.at(1))
	  del = sign*TMath::Sqrt(taups.at(1)/taups.at(0));
	else
	  del = sign*TMath::Sqrt(taups.at(0)/taups.at(1));

	double frac0 = del*del/(1. + del*del);
	double frac1 = 1./(1. + del*del);

	cc = frac0*convCoeffs[mult0-1]->Eval(egam) + frac1*convCoeffs[mult1-1]->Eval(egam);
      }
      
      std::cout << std::setw(intw) << ni+1
		<< std::setw(intw) << nf+1
		<< std::setw(valw) << 1000.*egam 
		<< std::setw(valw) << mults
		<< std::setw(valw) << del
		<< cc << "\n";
    }
  }

  return;
}

void Nucleus::PrintAll() const {

  PrintLevels();
  PrintLifetimes();
  PrintBranchProbabilities();
  PrintBranchingRatios();
  //PrintMixingRatios();
  PrintMatrixElements();
  PrintGammaTransitions();
  
  return;
}

void Nucleus::PrintComparison(const Literature* lit) const {

  std::vector<double> ws = lit->GetWeights();
  
  int intw = 7;
  int valw = 12;
  int errw = 18;
  
  std::cout << std::left;
  std::cout << "\nBranching Ratios (w=" << ws[0] << ")\n";
  std::cout << std::setw(intw) << "ni" 
	    << std::setw(intw) << "nf1" 
	    << std::setw(intw) << "nf2" 
	    << std::setw(valw) << "Calc" 
	    << std::setw(valw) << "Lit" 
	    << std::setw(errw) << "Err"  
	    << std::setw(errw) << "nSig"  
	    << "Chi2\n";
  
  for(LitVal* lv : lit->GetBranchingRatios()) {
    
    int ni = lv->GetInitialIndex();
    int nf1 = lv->GetFinalIndex();
    int nf2 = lv->GetFinalIndexTwo();

    double br =  CalculateBranchingRatio(ni,nf1,nf2);
    double nSig = lv->NSigma(br);
    double chi2 = ws[0]*lv->Compare(br);

    double val = lv->GetValue();
    double eU = lv->GetErrorUp();
    double eD = lv->GetErrorDown();
    
    std::string err = Form("+%.3f-%.3f",eU,eD);
    std::cout << std::setw(intw) << ni+1
	      << std::setw(intw) << nf1+1
	      << std::setw(intw) << nf2+1
	      << std::setw(valw) << br
	      << std::setw(valw) << val
	      << std::setw(errw) << err
	      << std::setw(errw) << nSig
	      << chi2 << "\n";
  }
 
  std::cout << "\nLifetimes (w=" << ws[1] << ")\n";
  std::cout << std::setw(intw) << "index" 
	    << std::setw(valw) << "Calc" 
	    << std::setw(valw) << "Lit" 
	    << std::setw(errw) << "Err"  
	    << std::setw(errw) << "nSig"  
	    << "Chi2\n";
  
  for(LitVal* lv : lit->GetLifetimes()) {
    
    int index = lv->GetInitialIndex();
    
    double tau = CalculateLifetime(index);
    double nSig = lv->NSigma(tau);
    double chi2 = ws[1]*lv->Compare(tau);
    
    double val = lv->GetValue();
    double eU = lv->GetErrorUp();
    double eD = lv->GetErrorDown();
    
    std::string err = Form("+%.3f-%.3f",eU,eD);
    std::cout << std::setw(intw) << index+1
	      << std::setw(valw) << tau
	      << std::setw(valw) << val
	      << std::setw(errw) << err
	      << std::setw(errw) << nSig
	      << chi2 << "\n";
  }
  
  std::cout << "\nMixing Ratios (w=" << ws[2] << ")\n";
  std::cout << std::setw(intw) << "ni" 
	    << std::setw(intw) << "nf" 
	    << std::setw(valw) << "Calc" 
	    << std::setw(valw) << "Lit" 
	    << std::setw(errw) << "Err"  
	    << std::setw(errw) << "nSig"  
	    << "Chi2\n";
  
  for(LitVal* lv : lit->GetMixingRatios()) {
    
    int ni = lv->GetInitialIndex();
    int nf = lv->GetFinalIndex();

    double mx =  CalculateMixingRatio(ni,nf);
    double nSig = lv->NSigma(mx);
    double chi2 = ws[2]*lv->Compare(mx);

    double val = lv->GetValue();
    double eU = lv->GetErrorUp();
    double eD = lv->GetErrorDown();

    std::string err = Form("+%.3f-%.3f",eU,eD);
    std::cout << std::setw(intw) << ni+1
	      << std::setw(intw) << nf+1
	      << std::setw(valw) << mx
	      << std::setw(valw) << val
	      << std::setw(errw) << err
	      << std::setw(errw) << nSig
	      << chi2 << "\n";
  }

  std::cout << "\nMatrix Elements (w=" << ws[3] << ")\n";
  std::cout << std::setw(intw) << "mult" 
	    << std::setw(intw) << "n1" 
	    << std::setw(intw) << "n2" 
	    << std::setw(valw) << "Calc" 
	    << std::setw(valw) << "Lit" 
	    << std::setw(errw) << "Err"  
	    << std::setw(errw) << "nSig"  
	    << "Chi2\n";
  
  for(LitVal* lv : lit->GetMatrixElements()) {

    int n1 = lv->GetInitialIndex();
    int n2 = lv->GetFinalIndex();
    int mult = lv->GetFinalIndexTwo();

    double me = GetMatrixElement(n1,n2,mult)->GetValue();
    double nSig = lv->NSigma(me);
    double chi2 = ws[3]*lv->Compare(me);

    double val = lv->GetValue();
    double eU = lv->GetErrorUp();
    double eD = lv->GetErrorDown();
    
    std::string err = Form("+%.3f-%.3f",eU,eD);
    std::cout << std::setw(intw) << mult
	      << std::setw(intw) << n1+1
	      << std::setw(intw) << n2+1
	      << std::setw(valw) << me
	      << std::setw(valw) << val
	      << std::setw(errw) << err
	      << std::setw(errw) << nSig
	      << chi2 << "\n";
  }

  return;
}

void Nucleus::Write(int A, int Z) const {
  Write(name+".txt",A,Z);
  return;
}

void Nucleus::Write(std::string file_name, int A, int Z) const {
  
  //This function writes a Cynus nucleus file
  std::ofstream outFile(file_name.c_str());
  
  outFile << Z << "\t" << A << "\t" << levels.size() << "\t7\n";
  for(Level* lvl : levels)
    outFile << lvl->GetIndex() << "\t" << lvl->GetEnergy() << "\t" << lvl->GetSpin() << "\t" << lvl->GetParity() << "\n";
  
  for(MatrixElement* me : matrix_elements)
    outFile << me->GetIndex1() << "\t" << me->GetIndex2() << "\t" << me->GetValue() << "\t\t" << me->GetMultS() << "\n";
  
  return;
}

void Nucleus::WriteLevelScheme() const {
  WriteLevelScheme(name+".lvl");
  return;
}

void Nucleus::WriteLevelScheme(std::string file_name) const {

  //This function writes a level scheme file for G4CLX
  std::ofstream outFile(file_name.c_str());
  
  for(Level* lvl : levels) {
    
    int index = lvl->GetIndex();
    if(index == 0) //skip ground state
      continue;

    double en = lvl->GetEnergy();
    double sp = lvl->GetSpin();
    double lt = CalculateLifetime(index);

    std::vector<int> indices;
    for(Level* lvl2 : levels) {

      int index2 = lvl2->GetIndex();
      double prob = CalculateBranchProbability(index,index2);
      if(prob > 1E-18)
	indices.push_back(index2);
      
    }

    int num = indices.size();
    outFile << index << " " << 1000*en << " " << sp << " " << lt << " " << num << "\n";
    
    for(int i=0;i<num;++i) {
     
      int index2 = indices.at(i);
      double egam = en - levels.at(index2)->GetEnergy();
      double prob = CalculateBranchProbability(index,index2);
      std::vector<int> mults = GetTransitionMultipolarities(index,index2);

      int l1, l2;
      double mx = 0.0;
      double cc = 0.0;
      if(mults.size() == 1) {
	l1 = mults.at(0);
	l2 = 0;
	
	cc = convCoeffs[l1-1]->Eval(egam);
      }
      else {
	
	l1 = mults.at(0);
	l2 = mults.at(1);
	
	//Only E2/M1 mixing ratio implemented right now
	if((l1 == 2 && l2 == 7) || (l1 == 7 && l2 == 2)) {
	  mx = CalculateMixingRatio(index,index2);

	  double fe2 = mx*mx/(mx*mx + 1.0);
	  double fm1 = 1.0/(mx*mx + 1.0);
	  cc = fe2*convCoeffs[1]->Eval(egam) + fm1*convCoeffs[6]->Eval(egam);;
	}
      }

      if(l1 == 7)
	l1 = 1;
      if(l1 == 8)
	l1 = 2;
      if(l2 == 7)
	l2 = 1;
      if(l2 == 8)
	l2 = 2;

      if(l2 > 0 && l2 < l1)
	std::swap(l1,l2);

      outFile << " " << index2 << " " << prob << " " << l1 << " " << l2 << " " << mx << " " << cc << "\n";
    
    }
    
  }

  return;
}
