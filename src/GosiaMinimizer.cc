#include "GosiaMinimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TString.h"
#include <fstream>

GosiaMinimizer::GosiaMinimizer() : GosiaMinimizer("Minuit2","Migrad") {;}
GosiaMinimizer::GosiaMinimizer(std::string meth) : GosiaMinimizer(meth,"") {;}

GosiaMinimizer::GosiaMinimizer(std::string meth, std::string alg) {
  
  chi_cut = 1.0;
  maxIter = 10000;
  maxCalls = 50000;
  numTrys = 1;
  fitTol = 0.00001;

  validate = false;
  write = false;
  calc_unc = false;
  limited = false;
  relCS = false;

  SetMethod(meth);
  SetAlgorithm(alg);

  theFCN = new GosiaMinimizerFCN();

  mini = ROOT::Math::Factory::CreateMinimizer(meth,alg); 
  if(!mini) {
    std::cout << "Either method or algorithm is unknown!\nmethod: " << method 
	      << "\nalgorithm: " << algorithm << std::endl;
  }
  
  return;
}

GosiaMinimizer::~GosiaMinimizer() {

  delete theFCN;

  return;
}

void GosiaMinimizer::SetBeam(Nucleus* nuc) {

  theFCN->beam = nuc;
  theFCN->beam_name = nuc->GetName();
  
  return;
}

void GosiaMinimizer::AddTarget(Nucleus* nuc) {

  theFCN->targets.push_back(nuc);
  theFCN->target_names.push_back(nuc->GetName());
  
  return;
}

void GosiaMinimizer::SetBeamData(ExperimentalData* data) {

  theFCN->beam_data = data;
  
  int num_exp = data->Size();
  Resize(num_exp);
  
  return;

}

void GosiaMinimizer::AddTargetData(ExperimentalData* data) {

  theFCN->target_data.push_back(data);
  
  int num_exp = data->Size();
  Resize(num_exp);
  
  return;

}

void GosiaMinimizer::Resize(int size) {

  int cur_size = theFCN->scalings.size();
  for(int i=0;i<size;++i) {
    theFCN->indices.push_back(i + cur_size);
    theFCN->scalings.push_back(1.0);
    theFCN->factors.push_back(1.0);
  }
  
  return;
}

void GosiaMinimizer::SetupParameters(std::vector<double> scales) {

  Nucleus* beam = theFCN->beam;
  
  int par_num = 0;
  for(int i=0;i<beam->GetNumMatrixElements();++i) {
    
    MatrixElement* me = beam->GetMatrixElement(i);
    if(me->IsFixed())
      continue;
    
    double val = me->GetValue();
    double ulim = me->GetUpperLimit();
    double llim = me->GetLowerLimit();    
    double st = TMath::Abs(val/1000.0);
    
    if(limited)
      mini->SetLimitedVariable(par_num,Form("beamME%02d",i+1),val,st,llim,ulim);
    else
      mini->SetVariable(par_num,Form("beamME%02d",i+1),val,st);
    ++par_num;
    
  }

  for(int i=0;i<theFCN->targets.size();++i) {
    
    Nucleus* targ = theFCN->targets[i];
    for(int j=0;j<targ->GetNumMatrixElements();++j) {

      MatrixElement* me = targ->GetMatrixElement(j);
      if(me->IsFixed())
	continue;

      double val = me->GetValue();
      double ulim = me->GetUpperLimit();
      double llim = me->GetLowerLimit();    
      double st = TMath::Abs(val/1000.0);
      
      if(limited)
	mini->SetLimitedVariable(par_num,Form("targ%dME%02d",i+1,j+1),val,st,llim,ulim);
      else
	mini->SetVariable(par_num,Form("targ%dME%02d",i+1,j+1),val,st);
      ++par_num;
      
    }
  }

  if(!theFCN->beam_data)
    return;
  
  //Minimum one free scaling parameter for experimental data
  int scale_num = 0;
  mini->SetVariable(par_num,"Scaling00",scales.at(scale_num),scales.at(scale_num)/100.0);
  ++par_num;
  ++scale_num;
  
  int size = theFCN->indices.size();
  for(int i=1;i<size;++i) {
    
    int num = theFCN->indices.at(i);
    if(num == scale_num) {

      double scl = scales.at(scale_num);
      mini->SetVariable(par_num,Form("Scaling%02d",num),scl,scl/100.0);
      ++par_num;
      ++scale_num;
    }
  }

  return;
}

void GosiaMinimizer::LinkBeamExperiments(int exp1, int exp2, double rel) {

  LinkExperiments(exp1,exp2,rel);

  return;
}

void GosiaMinimizer::LinkBeamTargetExperiments(int exp1B, int numE, int target_num) {

  ExperimentalData* dataB = theFCN->beam_data;

  int first_targ_exp = dataB->Size() + 1;
  for(int i=0;i<target_num-1;++i)
    first_targ_exp += theFCN->target_data[i]->Size();
  
  int nrm_index = dataB->GetNormalizationIndex();
  double scale0 = dataB->GetExperiment(exp1B)->GetAllIntiYields().at(nrm_index)->GetValue() / 
    dataB->GetExperiment(exp1B)->GetAllPointYields().at(nrm_index)->GetValue();
  
  for(int i=0;i<numE;++i) {
    
    dataB->RecorrectExp(exp1B + i,scale0);
    LinkExperiments(exp1B + i,first_targ_exp + i,1.0);
    
  }
  
  return;
  
}

void GosiaMinimizer::LinkExperiments(int exp1, int exp2, double rel) {
  
  std::cout << "Linking Experiments " << exp1 << " and " << exp2 << std::endl;
  if(exp2 < exp1) {
    std::cout << "Must link larger experiment number to smaller experiment number" << std::endl;
    return;
  }

  int old_ind = theFCN->indices.at(exp2-1);

  theFCN->indices.at(exp2-1)= theFCN->indices.at(exp1-1);
  theFCN->scalings.at(exp2-1) = rel;

  for(int i=0;i<theFCN->indices.size();++i)
    if(theFCN->indices.at(i) >= old_ind)
      theFCN->indices.at(i) -= 1;

  return;
}

std::vector<double> GosiaMinimizer::FindInitialScalings() {

  //Fix beam matrix elements at initial values
  std::vector<int> beam_inds;
  for(int i=0;i<theFCN->beam->GetNumMatrixElements();++i) {
    MatrixElement* me = theFCN->beam->GetMatrixElement(i);
    if(!(me->IsFixed())) {
      beam_inds.push_back(i);
      me->Fix();
    }
    
  }
  
  //Fix target matrix elements at initial values
  std::vector<std::vector<int>> targ_inds;
  for(Nucleus* targ : theFCN->targets) {
    targ_inds.emplace_back();
    for(int i=0;i<targ->GetNumMatrixElements();++i) {
      MatrixElement* me = targ->GetMatrixElement(i);
      if(!(me->IsFixed())) {
	targ_inds.back().push_back(i);
	me->Fix();
      }
    }
  }
  
  //Make a temporary minimizer to find initial scaling parameters for the GOSIA experiments
  ROOT::Math::Minimizer* tmp_mini = ROOT::Math::Factory::CreateMinimizer(method,algorithm); 
  if(!tmp_mini)
    std::cout << "Either method or algorithm is unknown!\nmethod: " << method 
	      << "\nalgorithm: " << algorithm << std::endl;
  
  tmp_mini->SetErrorDef(chi_cut);
  tmp_mini->SetMaxFunctionCalls(maxCalls);
  tmp_mini->SetMaxIterations(maxIter);
  tmp_mini->SetTolerance(fitTol);

  //Minimum one free scaling parameter
  int scale_num = 0;
  tmp_mini->SetVariable(scale_num,"Scaling00",1000,100);
  ++scale_num;
  
  int size = theFCN->indices.size();
  for(int i=1;i<size;++i) {
    
    int num = theFCN->indices.at(i);
    
    if(num == scale_num) {
      tmp_mini->SetVariable(scale_num,Form("Scaling%02d",num),1000,100);
      scale_num++;
    }
  }
  
  //Make the fuction and do the minimization
  ROOT::Math::Functor func(*theFCN,scale_num);
  tmp_mini->SetFunction(func);

  tmp_mini->Minimize();
  const double* x = tmp_mini->X();
  std::cout << "\nInitial Chi2: " << tmp_mini->MinValue();

  //Record scaling parameters for future use
  std::vector<double> pars;
  for(int i=0;i<scale_num;++i)
    pars.push_back(x[i]);

  //Release the previously fixed matrix elements
  for(int i : beam_inds)
    theFCN->beam->GetMatrixElement(i)->Release();

  for(int i=0;i<theFCN->targets.size();++i) {
    
    std::vector<int> inds = targ_inds.at(i);
    Nucleus* targ = theFCN->targets.at(i);
    
    for(int j : inds)
      targ->GetMatrixElement(j)->Release();
  }
  
  delete tmp_mini;
  return pars;
  
  /*
  int num = theFCN->beam->GetNumFree();
  for(Nucleus* targ : theFCN->targets)
    num += targ->GetNumFree();

  for(int i=0;i<num;++i)
    mini->FixVariable(i);
    
  mini->Minimize();

  for(int i=0;i<num;++i)
    mini->ReleaseVariable(i);

  return;
  
  */

}

void GosiaMinimizer::Minimize() {

  if(!theFCN->beam) {
    std::cout << "No beam nucleus set!" << std::endl;
    return;
  }

  if(write)
    theFCN->CreateNTuple();
  
  std::vector<double> scales;
  if(theFCN->beam_data) {
    
    if(relCS)
      theFCN->FillFactors();
    
    scales = FindInitialScalings();
    
  }
  SetupParameters(scales);

  int size = mini->NFree();
  ROOT::Math::Functor func(*theFCN,size);
  
  mini->SetFunction(func);
  mini->SetErrorDef(chi_cut);
  mini->SetMaxFunctionCalls(maxCalls);
  mini->SetMaxIterations(maxIter);
  mini->SetTolerance(fitTol);
  
  std::cout << "\n" << size << " free parameters" << std::endl;;
  for(int i=0;i<size;++i)
    std::cout << " Par" << i << " " << mini->VariableName(i) << ": " << mini->X()[i] << std::endl;
  std::cout << "Minimizing..." << std::flush;

  int status = 3;
  int count = 0;
  while(count < numTrys && status != 0) {
    mini->Minimize();
    status = mini->Status();
    ++count;
  }
  
  if(status == 0) {
    std::cout << " Done! (" << count << " attempt"; 
    if(count > 1)
      std::cout << "s";
    std::cout << ".)" << std::endl;
  }
  else {
    std::cout << " Failed. (" << count << " attempt";
    if(count > 1)
      std::cout << "s";
    std::cout << ".)" << std::endl;
  }

  Nucleus* beam = theFCN->beam;
  int numM = beam->GetNumMatrixElements();

  //Write all matrix elements to file (a GOSIA matrix element file)
  std::ofstream meFileB((theFCN->beam_name + ".bst-minimized").c_str());
  for(int i=0;i<numM;++i) {
    double val = beam->GetMatrixElementValue(i);
    meFileB << std::setprecision(19) << val << "\n";
  }
  meFileB.close();

  //Perform error validation if requested
  if(validate && status == 0) {
    
    std::cout << "Validating Errors..." << std::flush;
    mini->Hesse();
    
    if(mini->Status() == 0)
      std::cout << " Done!" << std::endl;
    else
      std::cout << " Failed." << std::endl;
  }
  
  std::cout << "\n";
  Print();
  
  //Write Migrad (symmetric) uncertainties to file
  //All free parameters
  const double* min = mini->X();
  const double* errors = mini->Errors();

  std::ofstream uncFile("Symmetric_Errors.txt");
  for(int i=0;i<size;++i) {
    
     double val = min[i];
     double err = errors[i];
     std::string name = mini->VariableName(i);
     
     uncFile << name << ": " << val << " +- " << err << "\n";

  }
  uncFile.close();
  
  //Calculate MINOS (asymmetric) uncertainties if requested and write them to file
  //Only for free beam matrix elements
  if(calc_unc) {
    
    std::cout << "\nCalculating MINOS uncertainties for free beam matrix elements..." << std::endl;
    std::ofstream minosFile("MINOS_Errors.txt");

    int numF = beam->GetNumFree();
    for(int i=0;i<numF;++i) {

      double eU, eD;
      mini->GetMinosError(i,eD,eU);
      int statusM = mini->Status();
      
      double val = mini->X()[i];
      std::string name = mini->VariableName(i);

      std::cout << name << ": " << val << " +" << eU << " " << eD << " (status = " << statusM << ")" << std::endl;
      minosFile << name << ": " << val << " +" << eU << " " << eD << " (status = " << statusM << ")" << "\n";
    }
    
    minosFile.close();
    beam->FillFromBSTFile(theFCN->beam_name + ".bst-minimized"); //Reset MEs to minimum   
  
  }
  
  if(write) {
    
    TNtuple* tpl = theFCN->ntuple;
    if(tpl) {
      TFile* outFile = new TFile("ParameterSpace.root","RECREATE"); 
      theFCN->ntuple->Write();
      outFile->Close();
    }
  }
 
  if(theFCN->beam_data)
    UpdateScalings(min);
  
  return;
}

void GosiaMinimizer::UpdateScalings(const double* x) {
  
  int numM = theFCN->beam->GetNumFree();
  for(Nucleus* targ : theFCN->targets)
    numM += targ->GetNumFree();

  int numE_beam = 0;
  if(theFCN->beam_data) {
    numE_beam += theFCN->beam_data->Size();

    std::vector<double> scales;
    for(int i=0;i<numE_beam;++i)
      scales.push_back(theFCN->scalings.at(i)*x[numM + theFCN->indices.at(i)]);
    
    theFCN->beam_data->SetScalings(scales);

  }

  int numE_targs = 0;
  for(int i=0;i<theFCN->target_data.size();++i) {
    
    ExperimentalData* dataT = theFCN->target_data.at(i);
    int size = dataT->Size();
    
    std::vector<double> scales;
    for(int j=0;j<size;++j)
      scales.push_back(theFCN->scalings.at(numE_beam + numE_targs + j)*x[numM + theFCN->indices.at(numE_beam + numE_targs + j)]);
    
    dataT->SetScalings(scales);
    
    numE_targs += size;
  }

  return;
}
