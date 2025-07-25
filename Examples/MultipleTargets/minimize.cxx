#include "Nucleus.h"
#include "ExperimentalData.h"
#include "GosiaMinimizer.h"

int main(int argc, char** argv) {

  std::string beam_name = "kr78";
  Nucleus* beam = new Nucleus(beam_name);
  beam->CreateFromGosiaInputFile(); //beam_name.POIN.inp
  beam->FillFromBSTFile(); //beam_name.bst
  
  //For conversion coefficients
  std::vector<double> ensB = {0.30,0.35,0.40,0.45,0.50,0.60,0.70,0.80,0.90,1.00,1.25,1.50,1.75,2.00,2.25,2.50,2.75,3.00};
  std::vector<double> ccB1 = {0.00370,0.00243,0.001708,0.001263,0.000971,0.000629,0.000444,0.000334,0.000263,0.000214,0.000219,
			      0.000352,0.000525,0.000693,0.000852,0.001003,0.001140,0.001261};
  std::vector<double> ccB2 = {0.01700,0.01005,0.00646,0.00443,0.00319,0.00186,0.001214,0.000855,0.000637,0.000495,0.000318,0.000291,
			      0.000344,0.000427,0.000527,0.000631,0.000735,0.000837};
  std::vector<double> ccB7 = {0.00728,0.00499,0.00362,0.00274,0.00214,0.001410,0.000999,0.000746,0.000579,0.000463,0.000306,0.000272,
			      0.000311,0.000384,0.000473,0.000570,0.000669,0.000766};
  
  beam->SetConverionCoefficients(1,ensB,ccB1); //mult, energies, coefficients
  beam->SetConverionCoefficients(2,ensB,ccB2);
  beam->SetConverionCoefficients(7,ensB,ccB7); 
  beam->PrintAll();
  
  Literature* beam_lit = new Literature(beam_name);
  beam_lit->CreateFromFile(); //beam_name.lit
  
  ExperimentalData* beam_data = new ExperimentalData(beam_name);
  beam_data->ReadGosiaFiles(); //beam_name.POIN.inp, beam_name.POIN.out, beam_name.INTI.out
  beam_data->ReadDataFile(); //beam_name.yld
  beam_data->Correct();
  
  std::string targ1_name = "pt196";
  Nucleus* targ1 = new Nucleus(targ1_name);
  targ1->CreateFromGosiaInputFile();
  targ1->FillFromBSTFile();

  std::vector<double> ensT1 = {0.25,0.30,0.325,0.35,0.375,0.40,0.45,0.50,0.55,0.60,0.70,0.8,0.9,1.0};
  std::vector<double> ccT12 = {0.1724,0.0985,0.0779,0.0631,0.0521,0.0438,0.0322,0.0248,0.0197,0.01610,0.01143,0.00860,0.00676,0.00547};
  std::vector<double> ccT17 = {0.503,0.305,0.246,0.201,0.1672,0.1407,0.1029,0.0780,0.0607,0.0484,0.0325,0.0231,0.0171,0.01309};
  
  targ1->SetConverionCoefficients(2,ensT1,ccT12);
  targ1->SetConverionCoefficients(7,ensT1,ccT17);
  targ1->PrintAll();

  Literature* targ1_lit = new Literature(targ1_name);
  targ1_lit->CreateFromFile();
  
  ExperimentalData* targ1_data = new ExperimentalData(targ1_name);
  targ1_data->ReadGosiaFiles();
  targ1_data->ReadDataFile();
  targ1_data->Correct();

  std::string targ2_name = "pt194";
  Nucleus* targ2 = new Nucleus(targ2_name);
  targ2->CreateFromGosiaInputFile();
  targ2->FillFromBSTFile();
  
  targ2->SetConverionCoefficients(2,ensT1,ccT12);
  targ2->SetConverionCoefficients(7,ensT1,ccT17);
  targ2->PrintAll();
  
  Literature* targ2_lit = new Literature(targ2_name);
  targ2_lit->CreateFromFile();

  ExperimentalData* targ2_data = new ExperimentalData(targ2_name);
  targ2_data->ReadGosiaFiles();
  targ2_data->ReadDataFile();
  targ2_data->Correct();

  std::string targ3_name = "pd110";
  Nucleus* targ3 = new Nucleus(targ3_name);
  targ3->CreateFromGosiaInputFile();
  targ3->FillFromBSTFile();

  std::vector<double> ensT3 = {0.30,0.35,0.40,0.45,0.50,0.60,0.70,0.80,0.90,1.00,1.25,1.50,1.75,2.00};
  std::vector<double> ccT32 = {0.0296,0.0179,0.01171,0.00816,0.00598,0.00358,0.00238,0.001699,0.001282,0.001006,0.000635,0.000511,0.000504,0.000551};
  std::vector<double> ccT37 = {0.0198,0.01339,0.00961,0.00720,0.00558,0.00361,0.00252,0.00186,0.001425,0.001128,0.000708,0.000546,0.000514,0.000545};

  targ3->SetConverionCoefficients(2,ensT3,ccT32);
  targ3->SetConverionCoefficients(7,ensT3,ccT37);
  targ3->PrintAll();
  
  Literature* targ3_lit = new Literature(targ3_name);
  targ3_lit->CreateFromFile();

  ExperimentalData* targ3_data = new ExperimentalData(targ3_name);
  targ3_data->ReadGosiaFiles();
  targ3_data->ReadDataFile();
  targ3_data->Correct();

  GosiaMinimizer* mini = new GosiaMinimizer();
  mini->SetFitTolerance(0.00001);
  mini->SetNumTrys(5);
  //mini->IncludeRelativeCS(); //Includes relative Rutherford cross section factor (Dsig) for coupling experiments
  //mini->ValidateErrors(); //Flag to perform Hessian calculation
  //mini->LimitMatrixElements(); //Apply upper and lower limits to matrix elements
  //mini->WriteParameterSpace(); //Flag to write the parameter space to a ROOT file
  //mini->CalculateMINOSErrors(); //Flag for MINOS uncertainties. Only for free beam matrix elements. Very time consuming.
  //mini->FixScalingParameters(); //Fixes the the scaling of each experiment at its initial value

  mini->SetBeam(beam);
  mini->SetBeamLiterature(beam_lit);
  mini->SetBeamData(beam_data);
  
  mini->AddTarget(targ1);
  mini->AddTargetLiterature(targ1_lit);
  mini->AddTargetData(targ1_data);
  
  mini->AddTarget(targ2);
  mini->AddTargetLiterature(targ2_lit);
  mini->AddTargetData(targ2_data);

  mini->AddTarget(targ3);
  mini->AddTargetLiterature(targ3_lit);
  mini->AddTargetData(targ3_data);
  
  
  /*
  //With relative CS
  mini->LinkBeamTargetExperiments(1,5,1); //first beam expiment, number of experiments, target number
  mini->LinkBeamExperiments(6,7,0.39178260281956); //experiment1, experiment2, relative scaling (scales exp2 to exp1). exp2 > exp1
  mini->LinkBeamExperiments(6,8,0.21274922666436);
  mini->LinkBeamExperiments(6,9,0.13805731574278);
  mini->LinkBeamTargetExperiments(13,5,2);
  mini->LinkBeamTargetExperiments(24,3,3);
  */
  
  //Make sure you have added all the data (beam and targets) to the minimizer before you begin linking experiments
  //Without relative CS
  mini->LinkBeamTargetExperiments(1,5,1); //first beam expiment, number of experiments, target number
  mini->LinkBeamExperiments(6,7); //experiment 1, experiment 2, relative scaling (scales exp2 to exp1). exp2 > exp1
  mini->LinkBeamExperiments(6,8);
  mini->LinkBeamExperiments(6,9);
  mini->LinkBeamExperiments(6,10);
  mini->LinkBeamExperiments(6,11);
  //mini->LinkBeamExperiments(6,12);
  mini->LinkBeamTargetExperiments(13,5,2);
  mini->LinkBeamExperiments(18,19);
  mini->LinkBeamExperiments(18,20);
  mini->LinkBeamExperiments(21,22);
  mini->LinkBeamExperiments(21,23);
  mini->LinkBeamTargetExperiments(24,3,3); 
  
  /*  
  //gsb
  beam->ReleaseMatrixElement(1,2,2); //index1, index2, mult
  beam->ReleaseMatrixElement(2,2,2);
  beam->ReleaseMatrixElement(2,3,2);
  beam->ReleaseMatrixElement(3,3,2);
  beam->ReleaseMatrixElement(3,4,2);
  beam->ReleaseMatrixElement(4,4,2);
  beam->ReleaseMatrixElement(4,5,2);
  beam->ReleaseMatrixElement(5,5,2);
  
  //2_2+
  beam->ReleaseMatrixElement(1,8,2);
  beam->ReleaseMatrixElement(2,8,2);
  beam->ReleaseMatrixElement(8,8,2);
  beam->ReleaseMatrixElement(2,8,7);
  
  //3_1+
  beam->ReleaseMatrixElement(2,9,2);
  beam->ReleaseMatrixElement(3,9,2);
  beam->ReleaseMatrixElement(8,9,2);
  //beam->ReleaseMatrixElement(9,9,2);
  beam->ReleaseMatrixElement(2,9,7);
  beam->ReleaseMatrixElement(3,9,7);
  beam->ReleaseMatrixElement(8,9,7);
  
  //4_2+
  beam->ReleaseMatrixElement(2,10,2);
  beam->ReleaseMatrixElement(3,10,2);
  beam->ReleaseMatrixElement(8,10,2);
  beam->ReleaseMatrixElement(10,10,2);
  beam->ReleaseMatrixElement(3,10,7);
  
  //0_2+
  beam->ReleaseMatrixElement(2,13,2);
  //beam->ReleaseMatrixElement(13,25,7);
  
  //2_3+
  beam->ReleaseMatrixElement(1,14,2);
  beam->ReleaseMatrixElement(2,14,2);
  beam->ReleaseMatrixElement(3,14,2);
  beam->ReleaseMatrixElement(8,14,2);
  beam->ReleaseMatrixElement(13,14,2);
  beam->ReleaseMatrixElement(14,14,2);
  beam->ReleaseMatrixElement(2,14,7);
  beam->ReleaseMatrixElement(8,14,7);
  
  //Weird interference terms
  beam->ReleaseMatrixElement(3,8,2);
  beam->ReleaseMatrixElement(8,13,2);
  
  //Negative parity band
  //beam->ReleaseMatrixElement(2,15,1);
  //beam->ReleaseMatrixElement(3,16,1);
  //beam->ReleaseMatrixElement(4,16,1);
  beam->ReleaseMatrixElement(15,15,2);
  beam->ReleaseMatrixElement(15,16,2);
  //beam->ReleaseMatrixElement(16,16,2);
  beam->ReleaseMatrixElement(1,15,3);
  beam->ReleaseMatrixElement(2,16,3);
  //beam->ReleaseMatrixElement(3,17,3); 
  */
  
  /*
  //Weakly populated states
  beam->ReleaseMatrixElement(1,18,2);
  beam->ReleaseMatrixElement(2,18,2);
  beam->ReleaseMatrixElement(8,18,2);
  beam->ReleaseMatrixElement(14,18,2);
  beam->ReleaseMatrixElement(18,18,2);
  beam->ReleaseMatrixElement(2,18,7);
  beam->ReleaseMatrixElement(8,18,7);
  beam->ReleaseMatrixElement(14,18,7);

  beam->ReleaseMatrixElement(1,19,2);
  beam->ReleaseMatrixElement(2,19,2);
  beam->ReleaseMatrixElement(9,19,2);
  beam->ReleaseMatrixElement(19,19,2);
  beam->ReleaseMatrixElement(2,19,7);
  beam->ReleaseMatrixElement(9,19,7);
  
  beam->ReleaseMatrixElement(1,20,2);
  beam->ReleaseMatrixElement(3,20,2);
  beam->ReleaseMatrixElement(8,20,2);
  beam->ReleaseMatrixElement(9,20,2);
  beam->ReleaseMatrixElement(20,20,2);
  beam->ReleaseMatrixElement(8,20,7);
  beam->ReleaseMatrixElement(9,20,7);

  beam->ReleaseMatrixElement(2,21,2);
  beam->ReleaseMatrixElement(3,21,2);
  beam->ReleaseMatrixElement(14,21,2);
  beam->ReleaseMatrixElement(21,21,2);
  beam->ReleaseMatrixElement(3,21,2);

  beam->ReleaseMatrixElement(2,22,2);
  */
  
  /*  
  targ1->ReleaseMatrixElement(1,2,2);
  targ1->ReleaseMatrixElement(2,2,2);
  targ2->ReleaseMatrixElement(1,2,2);
  targ2->ReleaseMatrixElement(2,2,2);
  targ3->ReleaseMatrixElement(1,2,2);
  targ3->ReleaseMatrixElement(2,2,2);
  */
  
  mini->Minimize();
  //mini->Scan(11,-0.8,-0.7,2,true);
  //mini->Scan(11,-0.85,-0.73,5,true);
  //mini->Scan2D(11,-0.86,-0.72,29,5,0.782,0.794,25,false);
  //mini->Scan2D(80,-0.25,0.25,26,13,0.50,0.65,16,false);
  //mini->Scan2D(80,-0.22,0.23,31,13,0.51,0.61,21,false);
 
  //Print final yield and literature constraint comparisons
  targ1->PrintComparison(targ1_lit);
  targ2->PrintComparison(targ2_lit);
  targ3->PrintComparison(targ3_lit);
  
  targ1_data->PrintComparison();
  targ2_data->PrintComparison();
  targ3_data->PrintComparison();
  
  beam->PrintComparison(beam_lit);
  beam_data->PrintComparison();
  
  delete beam;
  delete beam_lit;
  delete beam_data;
  delete targ1;
  delete targ1_lit;
  delete targ1_data;
  delete targ2;
  delete targ2_lit;
  delete targ2_data;
  delete targ3;
  delete targ3_lit;
  delete targ3_data;
  delete mini;
  
  return 0;
}
