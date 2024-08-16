#include "Nucleus.h"
#include "TMath.h"
#include "TString.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

Nucleus::Nucleus(std::string nm) {
  
  name = nm;
  
  std::vector<double> ens = {0.01,0.02,1.0,3.5,4.0};
  std::vector<double> ccs = {0.0,0.0,0.0,0.0,0.0};

  for(int i=0;i<8;++i)
    convCoeffs[i] = new TSpline3(Form("CC%d",i+1),&ens[0],&ccs[0],5);

  return;

}

Nucleus::~Nucleus() {

  for(MatrixElement* me : matrix_elements)
    delete me;

  for(Level* lvl : levels)
    delete lvl;

  for(int i=0;i<8;++i)
    delete convCoeffs[i];

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
    int mult = me->GetMultipolarity();
    
    if(n2 == index && egam > 0.0) {  
      
      double cc = convCoeffs[mult-1]->Eval(egam);
      sum += (cc + 1.0)/DecayLifetime(mult,me->GetValue(),levels[n2]->GetSpin(),egam);
    }
    else if(n1 == index && egam < 0.0) {
      
      double cc = convCoeffs[mult-1]->Eval(-egam);
      sum += (cc + 1.0)/DecayLifetime(me->GetMultipolarity(),me->GetValue(),levels[n1]->GetSpin(),-egam);
    }

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
    if((ni == n1 && nf == n2) || (ni == n2 && nf == n1))
      sum += 1.0/DecayLifetime(me->GetMultipolarity(),me->GetValue(),levels[ni]->GetSpin(),egam);
  
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
  
  double sum1 = 0.0;
  double sum2 = 0.0;
  for(MatrixElement* me : matrix_elements) {

    int n1 = me->GetIndex2();
    int n2 = me->GetIndex1();
   
    if((ni == n1 && nf1 == n2) || (ni == n2 && nf1 == n1))
      sum1 += 1.0/DecayLifetime(me->GetMultipolarity(),me->GetValue(),levels[ni]->GetSpin(),egam1);
    else if((ni == n1 && nf2 == n2) || (ni == n2 && nf2 == n1))
      sum2 += 1.0/DecayLifetime(me->GetMultipolarity(),me->GetValue(),levels[ni]->GetSpin(),egam2);
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

  if(std::abs(valE2) < 1E-10 || std::abs(valM1) < 1E-10)
    return 0.0;

  return 0.835*egam*valE2/valM1;

}

double Nucleus::CompareLifetime(const LitVal* lv) const {
  
  double tau =  CalculateLifetime(lv->GetInitialIndex());
  return lv->Compare(tau);

}

double Nucleus::CompareBranchingRatio(const LitVal* lv) const {
  
  double br =  CalculateBranchingRatio(lv->GetInitialIndex(),lv->GetFinalIndex(),
				       lv->GetFinalIndexTwo()); 
  return lv->Compare(br);

}

double Nucleus::CompareMixingRatio(const LitVal* lv) const {
  
  double mx =  CalculateMixingRatio(lv->GetInitialIndex(),lv->GetFinalIndex()); 
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

void Nucleus::CreateFromGosiaInputFile() {
  
  levels.clear();
  matrix_elements.clear();
  
  std::ifstream inFile((name+".POIN.inp").c_str());
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
  
  /*
  std::string flag = " 0 0 0 0";
  if(!gosi)
    flag = " 0 0";
  std::string end = flag + " 0";
  */
  
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

void Nucleus::CreateFromGosiaOutputFile() {

  levels.clear();
  matrix_elements.clear();
  
  std::ifstream inFile((name + ".out").c_str());
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

void Nucleus::FillFromBSTFile() {
  
  std::ifstream inFile((name + ".bst").c_str());
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

void Nucleus::SetConverionCoefficients(int mult, std::vector<double> ens, std::vector<double> ccs) {

  int size = ens.size();
  if(size < 5) {
    std::cout << "Give at least 5 spline points for conversion coefficients (Mult=" << mult << ")" << std::endl;
    return;
  }

  if(ccs.size() != size) {
    std::cout << "Must give same number of energies (" << size << ") as values (" 
	      << ccs.size() << ") when setting conversion coefficients (Mult=" << mult << ")" 
	      << std::endl;
    return;
  }

  TSpline3* spline = convCoeffs[mult-1];
  delete spline;
  spline = NULL;
  
  convCoeffs[mult-1] = new TSpline3(Form("CC%d",mult),&ens[0],&ccs[0],size);

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
	    << "Prob\n";

  int size = levels.size();
  for(int i=0;i<size;++i) {
    
    int ni = levels[size-i-1]->GetIndex();
    for(int j=0;j<size;++j) {
      if(i==j)
	continue;
      
      int nf = levels[size-j-1]->GetIndex();
      double prob = CalculateBranchProbability(ni,nf);
      if(prob < 1E-10)
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
	    << "Ratio\n";
  
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
      if(prob < 1E-10)
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
	    << "Value\n";

  int size = levels.size();
  for(int i=0;i<size;++i) {
    
    int ni = levels[size-i-1]->GetIndex();
    for(int j=i+1;j<size;++j) {
      
      int nf = levels[size-j-1]->GetIndex();
      double del = CalculateMixingRatio(ni,nf);
      if(std::abs(del) < 1E-10)
	continue;
      
      std::cout << std::setw(intw) << ni+1
		<< std::setw(intw) << nf+1
		<< del << "\n";
    }
  }
  
  return;

}

void Nucleus::PrintMatrixElements() const {

  int intw = 7;
  int valw = 12;

  std::cout << std::left;
  std::cout << "\nMatrix Elements:\n";
  std::cout << std::setw(intw) << "mult" 
	    << std::setw(intw) << "n1" 
	    << std::setw(intw) << "n2" 
	    << "Value\n";
  
  for(int i=0;i<matrix_elements.size();i++) {
    
    MatrixElement* me = matrix_elements.at(i);
    
    int n1 = me->GetIndex1();
    int n2 = me->GetIndex2();
    int mult = me->GetMultipolarity();
    double val = me->GetValue();

    std::cout << std::setw(intw) << mult
	      << std::setw(intw) << n1+1
	      << std::setw(intw) << n2+1
	      << std::setw(valw) << val<< "\n";

  }

  return;
}

void Nucleus::PrintAll() const {

  PrintLevels();
  PrintLifetimes();
  PrintBranchProbabilities();
  PrintBranchingRatios();
  PrintMixingRatios();
  PrintMatrixElements();

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
