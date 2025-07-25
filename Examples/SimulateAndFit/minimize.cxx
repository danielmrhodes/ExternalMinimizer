#include "Nucleus.h"
#include "ExperimentalData.h"
#include "GosiaMinimizer.h"

int main(int argc, char** argv) {

  std::string beam_name = "sm148";
  
  Nucleus* beam = new Nucleus(beam_name);
  beam->CreateFromGosiaInputFile(); //beam_name.POIN.inp
  beam->FillFromBSTFile(); //beam_name.bst
  beam->CheckMatrixElements();
  
  std::vector<double> ensB = {0.05,0.075,0.10,0.15,0.20,0.30,0.45,0.60,0.80,1.0,1.25,1.50,1.75,2.0,2.5,3.0};
  std::vector<double> ccB1 = {1.732,0.599,0.276,0.0916,0.0423,0.01471,0.00549,0.00288,0.001588,0.001032,0.000736,0.000709,0.000783,
                              0.000892,0.001130,0.001354};
  std::vector<double> ccB2 = {34.1,6.80,2.34,0.556,0.209,0.0562,0.01707,0.00801,0.00405,0.00249,0.001586,0.001177,0.000997,0.000932,
                              0.000956,0.001064};
  std::vector<double> ccB3 = {1760,170.5,35.5,4.56,1.209,0.224,0.0516,0.0209,0.00934,0.00531,0.00315,0.00214,0.001605,0.001315,0.001078,
                              0.001044};

  beam->SetConverionCoefficients(1,ensB,ccB1); //mult, energies, coefficients
  beam->SetConverionCoefficients(2,ensB,ccB2);
  beam->SetConverionCoefficients(3,ensB,ccB3);
  
  beam->PrintAll();
  //beam->Write(148,62);
  //beam->WriteLevelScheme();
  
  Literature* beam_lit = new Literature(beam_name);
  beam_lit->CreateFromFile(); //beam_name.lit
  
  ExperimentalData* beam_data = new ExperimentalData(beam_name);
  beam_data->ReadGosiaFiles(); //beam_name.POIN.inp, beam_name.POIN.out, beam_name.INTI.out
	
  /* 
  /////////////////
  //Uncomment this block to generate data from the calculations and write it to file
  //Gamma-ray detection efficiency curve. Energy (x) in keV.
  TF1* eff = new TF1("eff","6.798*TMath::Power(x + 100.,-0.621)",10.,5000.);

  //Experimental parameters
  double Nav = 6.02214076*TMath::Power(10,23); //Avagadro's number
  double targ_mm = 207.98; // target molar mass (g/mol)
  double peff = 0.8; //particle detection efficiency
  double coincLT = 1.0; //Coincidence livetime
  double mtm = 1.0; //measurement time (days)
  double rate = 6.0*TMath::Power(10,8); //beam rate (pps)
  double sa = 4.0*TMath::Pi(); //germanium solid angle (sr)
  
  double scale = TMath::Power(10,-30)*24*60*60*mtm*rate*sa*peff*coincLT*Nav/targ_mm;
  beam_data->GenerateAllData(beam,eff,scale,2000);
  beam_data->WriteDataFile(148,62); //writes to beam_name.yld
  
  delete beam;
  delete beam_lit;
  delete beam_data;
  delete eff;
  
  return 1;
  /////////////////
  */
 
  //Read in data
  beam_data->ReadDataFile(); //beam_name.yld
  beam_data->Correct();

  int numE = beam_data->Size();
  int numY = beam_data->GetNumRawData();
  int numC = beam_lit->Size();
 
  std::cout << "\n\n######Begin Fitting Routine######\n\n";
  std::cout << "Fitting " << numY << " experimental yields (" << numE << " data sets) and " << numC << " literature constraints\n";

  GosiaMinimizer* mini = new GosiaMinimizer(); 
  mini->SetFitTolerance(0.00001);
  //mini->SetNumTrys(10);
  //mini->ValidateErrors(); //Flag to perform Hessian calculation
  //mini->LimitMatrixElements(); //Apply upper and lower limits to matrix elements
  //mini->WriteParameterSpace(); //Flag to write the parameter space to a ROOT file
  //mini->CalculateMINOSErrors(); //Flag for MINOS uncertainties. Only for free beam matrix elements. Very time consuming.
  //mini->SetPrintLevel(1000);

  mini->SetBeam(beam);
  mini->SetBeamLiterature(beam_lit);
  mini->SetBeamData(beam_data);
  
  //Leaving out Rutherford CS term, which is different from how GOSIA does it
  mini->LinkBeamExperiments(1,2);
  mini->LinkBeamExperiments(1,3);
  mini->LinkBeamExperiments(1,4);
  mini->LinkBeamExperiments(1,5);
  mini->LinkBeamExperiments(1,6);
  mini->LinkBeamExperiments(1,7);
  mini->LinkBeamExperiments(1,8);
  mini->LinkBeamExperiments(1,9);
  mini->LinkBeamExperiments(1,10);

  /*  
  //gsb
  beam->ReleaseMatrixElement(1,2,2); //index1, index2, mult
  beam->ReleaseMatrixElement(2,2,2);
  beam->ReleaseMatrixElement(2,3,2);
  beam->ReleaseMatrixElement(3,3,2);
  
  //negative parity band
  beam->ReleaseMatrixElement(1,10,3);
  beam->ReleaseMatrixElement(2,12,3);
  beam->ReleaseMatrixElement(10,12,2);
  beam->ReleaseMatrixElement(10,10,2);
  */
  mini->Minimize();
  
  std::cout << "\n######End Fitting Routine######\n\n"; 
  
  //Print final yield and literature constraint comparisons.
  beam_data->PrintComparison();
  beam->PrintComparison(beam_lit); 
   
  delete beam;
  delete beam_lit;
  delete beam_data;
  delete mini;
  
  return 0;
}
