#include "Experiment.h"
#include <iostream>
#include <iomanip>
#include "TMath.h"
#include "TString.h"

Experiment::Experiment() {;}
Experiment::~Experiment() {
  
  ClearAllIntiYields();
  ClearAllPointYields();
  ClearPointYields();
  ClearRawData();
  ClearCorrectedData();
  
  return;
}

void Experiment::ClearAllIntiYields() {

  for(Yield* yld : inti_yields_all)
    delete yld;
  inti_yields_all.clear();
  
  return;
}

void Experiment::ClearAllPointYields() {

  for(Yield* yld : point_yields_all)
    delete yld;
  point_yields_all.clear();
  
  return;
}

void Experiment::ClearPointYields() {

  for(Yield* yld : point_yields)
    delete yld;
  point_yields.clear();
  
  return;
}

void Experiment::ClearRawData() {

  for(YieldError* yld : data_raw)
    delete yld;
  data_raw.clear();
  
  return;
}

void Experiment::ClearCorrectedData() {

  for(YieldError* yld : data_corr)
    delete yld;
  data_corr.clear();
  
  return;
}

void Experiment::ZeroAllIntiYields() {

  for(Yield* yld : inti_yields_all)
    yld->SetValue(0.0);

  return;
}

void Experiment::ZeroAllPointYields() {

  for(Yield* yld : point_yields_all)
    yld->SetValue(0.0);
}

void Experiment::ZeroPointYields() {

  for(Yield* yld : point_yields)
    yld->SetValue(0.0);
  
  return;
}

Yield* Experiment::GetPointYield(int ni, int nf) const {

  for(Yield* y : point_yields_all)
    if(y->GetInitialIndices()[0] == ni && y->GetFinalIndices()[0] == nf)
      return y;
  
  std::cout << "No yield between states " << ni+1 << " and " << nf+1 << std::endl;
  return NULL;
  
}

Yield* Experiment::GetIntiYield(int ni, int nf) const {

  for(Yield* y : inti_yields_all)
    if(y->GetInitialIndices()[0] == ni && y->GetFinalIndices()[0] == nf)
      return y;
  
  std::cout << "No yield between states " << ni+1 << " and " << nf+1 << std::endl;
  return NULL;
  
}

void Experiment::PrintAllYields(std::string name) const {
  
  std::vector<Yield*> yldsP = GetAllPointYields();
  std::vector<Yield*> yldsI = GetAllIntiYields();
  std::vector<double> corrF = GetCorrectionFactors();
  
  int size1 = yldsI.size();
  if(yldsP.size() != size1) {
    std::cout << "Bad sizes for all yields in " << name << " Experiment " << num << "!" << std::endl;
    return;
  }

  int intw = 7;
  int valw = 20;
  double r2d = TMath::RadToDeg(); 
 
  std::cout << std::left;
  std::cout << "\nAll Yields for " << name << " Experiment " << num << "\n"
	    << "Energy Range: " << emin << " to " << emax << " MeV\n"
	    << "Angle Range: " << tmin*r2d << " to " << tmax*r2d << " deg\n"
	    << "Total Rutherford Cross Section: " << cs << " mb\n" 
	    << std::endl;

  std::cout << std::setw(intw) << "ni" 
	    << std::setw(intw) << "nf" 
	    << std::setw(valw) << "Point Yield" 
	    << std::setw(valw) << "Integrated Yield"
	    << "Correction Factor\n";
  
  for(int i=0;i<size1;++i) {
    
    Yield* ylI = yldsI.at(i);
    Yield* ylP = yldsP.at(i);
    double cc = corrF.at(i);

    //No multiplets for these
    int ni = ylI->GetInitialIndices().at(0);
    int nf = ylI->GetFinalIndices().at(0);
    double yi = ylI->GetValue();
    double yp = ylP->GetValue();

    std::cout << std::setw(intw) << ni+1
	      << std::setw(intw) << nf+1
	      << std::setw(valw) << yp
	      << std::setw(valw) << yi 
	      << cc << "\n";
  }

  return;

}

void Experiment::PrintRawData(std::string name) const {

  int size = data_raw.size();

  int intw = 8;
  int valw = 20;
  double r2d = TMath::RadToDeg();
   
  std::cout << std::left;
  std::cout << "\nRaw Data for " << name << " Experiment " << num << "\n"
	    << "Energy Range: " << emin << " to " << emax << " MeV\n"
	    << "Angle Range: " << tmin*r2d << " to " << tmax*r2d << " deg\n"
	    << "Total Rutherford Cross Section: " << cs << " mb\n" 
	    << std::endl;
  
  std::cout << std::setw(intw) << "ni" 
	    << std::setw(intw) << "nf" 
	    << std::setw(valw) << "Raw Yield" 
	    << std::setw(valw) << "Error"
	    << "Weight\n";

  for(int i=0;i<size;++i) {
    YieldError* yld = data_raw.at(i);

    std::vector<int> nis = yld->GetInitialIndices();
    std::vector<int> nfs = yld->GetFinalIndices();
    double yl = yld->GetValue();
    double eU = yld->GetErrorUp();
    double eD = yld->GetErrorDown();
    double wt = yld->GetWeight();

    std::string err = Form("+%.0f-%.0f",eU,eD);
    std::string ni = Form("%d",nis[0]+1);
    for(int j=1;j<nis.size();j++)
      ni += Form("+%d",nis[j]+1);

    std::string nf = Form("%d",nfs[0]+1);
    for(int j=1;j<nfs.size();j++)
      nf += Form("+%d",nfs[j]+1);

    std::cout << std::setw(intw) << ni
	      << std::setw(intw) << nf
	      << std::setw(valw) << yl
	      << std::setw(valw) << err
	      << wt << "\n";
  }
  
  return;
}

void Experiment::PrintRawAndCorrData(std::string name) const {

  int size = data_raw.size();

  int intw = 8;
  int valw = 20;
  double r2d = TMath::RadToDeg();
  
  std::cout << std::left;
  std::cout << "\nRaw and Corrected Data for " << name << " Experiment " << num << "\n"
	    << "Energy Range: " << emin << " to " << emax << " MeV\n"
	    << "Angle Range: " << tmin*r2d << " to " << tmax*r2d << " deg\n"
	    << "Total Rutherford Cross Section: " << cs << " mb\n" 
	    << std::endl;
  
  std::cout << std::setw(intw) << "ni" 
	    << std::setw(intw) << "nf" 
	    << std::setw(valw) << "Raw Yield"
	    << std::setw(valw) << "Error"
	    << std::setw(valw) << "Corr Yield"
	    << std::setw(valw) << "Error"
	    << "Weight\n";

  for(int i=0;i<size;++i) {
    
    YieldError* yldR = data_raw.at(i);
    std::vector<int> nis = yldR->GetInitialIndices();
    std::vector<int> nfs = yldR->GetFinalIndices();
    
    double ylR = yldR->GetValue();
    double eUR = yldR->GetErrorUp();
    double eDR = yldR->GetErrorDown();
    double wtR = yldR->GetWeight();

    YieldError* yldC = data_corr.at(i);
    double ylC = yldC->GetValue();
    double eUC = yldC->GetErrorUp();
    double eDC = yldC->GetErrorDown();
    double wtC = yldC->GetWeight();
    
    std::string errR = Form("+%.0f-%.0f",eUR,eDR);
    std::string errC = Form("+%.0f-%.0f",eUC,eDC);
    
    std::string ni = Form("%d",nis[0]+1);
    for(int j=1;j<nis.size();j++)
      ni += Form("+%d",nis[j]+1);

    std::string nf = Form("%d",nfs[0]+1);
    for(int j=1;j<nfs.size();j++)
      nf += Form("+%d",nfs[j]+1);

    std::cout << std::setw(intw) << ni
	      << std::setw(intw) << nf
	      << std::setw(valw) << ylR
	      << std::setw(valw) << errR
	      << std::setw(valw) << ylC
	      << std::setw(valw) << errC
	      << wtR << "\n";
  }
  
  return;
}

void Experiment::PrintComparison(std::string name, double scale, double factor) const {

  int size = data_raw.size();

  int intw = 8;
  int valw = 20;
  double r2d = TMath::RadToDeg();
   
  std::cout << std::left;
  std::cout << "\nYield Comparison for " << name << " Experiment " << num << "\n"
	    << "Energy Range: " << emin << " to " << emax << " MeV\n"
	    << "Angle Range: " << tmin*r2d << " to " << tmax*r2d << " deg\n"
	    << "Total Rutherford Cross Section: " << cs << " mb\n" 
	    << "Coupling Factor: " << factor << "\n"
	    << "Scaling Parameter: " << scale << "\n" << std::endl;
  
  std::cout << std::setw(intw) << "ni" 
	    << std::setw(intw) << "nf" 
	    << std::setw(valw) << "Corr Yield"
	    << std::setw(valw) << "Error"
	    << std::setw(intw) << "Weight"
	    << std::setw(valw) << "Calc Yield"
	    << std::setw(valw) << "nSig"
	    << "Chi2\n";

  for(int i=0;i<size;++i) {

    YieldError* yldC = data_corr.at(i);
    std::vector<int> nis = yldC->GetInitialIndices();
    std::vector<int> nfs = yldC->GetFinalIndices();
    
    double ylC = yldC->GetValue();
    double eUC = yldC->GetErrorUp();
    double eDC = yldC->GetErrorDown();
    double wtC = yldC->GetWeight();

    Yield* yldP = point_yields.at(i);
    double ylP = yldP->GetValue();

    double diff = ylC - scale*ylP/factor;
    double nSig;
    if(diff > 0.0)
      nSig = diff/eDC;
    else
      nSig = diff/eUC;
    
    double chi2 = wtC*nSig*nSig;

    std::string errC = Form("+%.0f-%.0f",eUC,eDC);
    std::string ni = Form("%d",nis[0]+1);
    for(int j=1;j<nis.size();j++)
      ni += Form("+%d",nis[j]+1);

    std::string nf = Form("%d",nfs[0]+1);
    for(int j=1;j<nfs.size();j++)
      nf += Form("+%d",nfs[j]+1);

    std::cout << std::setw(intw) << ni
	      << std::setw(intw) << nf
	      << std::setw(valw) << ylC
	      << std::setw(valw) << errC
	      << std::setw(intw) << wtC
	      << std::setw(valw) << scale*ylP/factor
	      << std::setw(valw) << nSig
	      << chi2 << "\n";

  }
  
  return;
}
