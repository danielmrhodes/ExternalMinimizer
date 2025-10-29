#include "GosiaMinimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TString.h"
#include "TMatrixD.h"
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
  fixed = false;
  calculated = false;
  
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
  
  return;
}

void GosiaMinimizer::AddTarget(Nucleus* nuc) {

  theFCN->targets.push_back(nuc);
  
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

  int cur_size = theFCN->rel_scalings.size();
  for(int i=0;i<size;++i) {
    theFCN->indices.push_back(i + cur_size);
    theFCN->rel_scalings.push_back(1.0);
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

  if(!theFCN->beam_data || fixed || calculated)
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

  targ_info.push_back({exp1B,numE,target_num});
  
  ExperimentalData* dataB = theFCN->beam_data;
  ExperimentalData* dataT = theFCN->target_data[target_num-1];

  int first_targ_exp = dataB->Size() + 1;
  for(int i=0;i<target_num-1;++i)
    first_targ_exp += theFCN->target_data[i]->Size();
  
  int nrm_index = dataB->GetNormalizationIndex();
  double scale0 = dataB->GetExperiment(exp1B)->GetAllIntiYields().at(nrm_index)->GetValue() / 
    dataB->GetExperiment(exp1B)->GetAllPointYields().at(nrm_index)->GetValue();
  
  for(int i=0;i<numE;++i) {
    
    dataB->RecorrectExp(exp1B + i,scale0);
    dataT->RecorrectExp(i+1,scale0);
    
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

  theFCN->indices.at(exp2-1) = theFCN->indices.at(exp1-1);
  theFCN->rel_scalings.at(exp2-1) = rel;

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
  //UpdateScalings(x);
    
  std::cout << "\nInitial Chi2: " << tmp_mini->MinValue();
  
  //Record scaling parameters for future use
  std::vector<double> pars;
  for(int i=0;i<scale_num;++i) {
    pars.push_back(x[i]);
    if(fixed || calculated)
      std::cout << Form("\n Scaling%02d: ",i) << pars.back();
  }
  std::cout << "\n" << std::endl;

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

}

void GosiaMinimizer::InitialSetup() {

  if(write)
    theFCN->CreateNTuple();
  
  std::vector<double> scales;
  if(theFCN->beam_data) {
    
    if(relCS)
      theFCN->FillFactors();
    
    scales = FindInitialScalings();
    if(fixed) {
      theFCN->fixed = true;
      theFCN->fixed_scales = scales;
    }
    if(calculated) {
      theFCN->calculated = true;
      theFCN->fixed_scales = scales;
    }
  }
  SetupParameters(scales);
  
  mini->SetErrorDef(chi_cut);
  mini->SetMaxFunctionCalls(maxCalls);
  mini->SetMaxIterations(maxIter);
  mini->SetTolerance(fitTol);

  int size = mini->NFree();
  if(size == 0) {
    mini->SetValidError(false); 
    validate = false;
  }

  return;
}

void GosiaMinimizer::Scan(int index, double llim, double ulim, int nstep, bool inti) {

  Nucleus* beam = theFCN->beam;
  if(!beam) {
    std::cout << "No beam nucleus set!" << std::endl;
    return;
  }
  std::string name = beam->GetName();

  MatrixElement* me = beam->GetMatrixElement(index - 1);
  if(!me->IsFixed()) {
    std::cout << "Sanning matrix element is not fixed." << std::endl;
    return;
  }
  
  InitialSetup();
  
  int size = mini->NFree();
  ROOT::Math::Functor func(*theFCN,size);
  mini->SetFunction(func);
  
  ExperimentalData* data = theFCN->beam_data;
  std::ofstream outFile("scan.txt");
  outFile << index << " " << llim << " " << ulim << " " << nstep << " " << int(inti) << "\n";
  
  double step = (ulim - llim)/double(nstep - 1);
  for(int i=0;i<nstep;++i) {

    double val = llim + i*step;
    me->SetValue(val);

    std::ofstream meFile((name + ".bst").c_str());
    for(MatrixElement* me1 : beam->GetMatrixElements())
      meFile << std::setprecision(19) << me1->GetValue() << "\n";
    meFile.close();

    if(data) {
      
      std::string cmd = "gosia < " + name + ".POIN.inp > /dev/null 2>&1";
      system(cmd.c_str());
      
      data->FillFromPoinOutput(); 

      if(inti) {
	
	cmd = "gosia < " + name + ".INTI.inp > /dev/null 2>&1";
	system(cmd.c_str());
	
	data->FillFromIntiOutput();
	data->Correct();

	//Recorrect target-linked experiments
	for(std::array<int,3> info : targ_info) {

	  int exp1B = info[0];
	  int numE = info[1];
	  int targ_num = info[2];
	  ExperimentalData* dataT = theFCN->target_data[targ_num-1];
	  
	  int nrm_index = data->GetNormalizationIndex();
	  double scale0 = data->GetExperiment(exp1B)->GetAllIntiYields().at(nrm_index)->GetValue() / 
	    data->GetExperiment(exp1B)->GetAllPointYields().at(nrm_index)->GetValue();

	  for(int j=0;j<numE;++j)  {
	    data->RecorrectExp(exp1B + j,scale0);
	    dataT->RecorrectExp(j+1,scale0);
	  }
	  
	}
      }
    }
    
    mini->Minimize();
    
    double chi2 = mini->MinValue();
    std::cout << val << " " << chi2 << std::endl;
    //Print();

    outFile << chi2;
    for(MatrixElement* me1 : beam->GetMatrixElements())
      outFile << " " << me1->GetValue();
    outFile << "\n";
    
  }
  
  outFile.close();
  
  return;
}

void GosiaMinimizer::Scan2D(int indexX, double llimX, double ulimX, int nstepX, int indexY, double llimY, double ulimY, int nstepY, bool inti) {

  Nucleus* beam = theFCN->beam;
  if(!beam) {
    std::cout << "No beam nucleus set!" << std::endl;
    return;
  }
  std::string name = beam->GetName();

  MatrixElement* meX = beam->GetMatrixElement(indexX - 1);
  MatrixElement* meY = beam->GetMatrixElement(indexY - 1);
  if(!meX->IsFixed() || !meY->IsFixed()) {
    std::cout << "Sanning matrix elements are not fixed." << std::endl;
    return;
  }

  InitialSetup();

  int size = mini->NFree();
  ROOT::Math::Functor func(*theFCN,size);
  mini->SetFunction(func);
  
  ExperimentalData* data = theFCN->beam_data;
  std::ofstream outFile("scan2D.txt");
  outFile << indexX << " " << llimX << " " << ulimX << " " << nstepX << " " <<  indexY
	  << " " << llimY << " " << ulimY << " " <<  nstepY << " " << int(inti) << "\n";

  double stepX = (ulimX - llimX)/double(nstepX - 1);
  double stepY = (ulimY - llimY)/double(nstepY - 1);  
  for(int i=0;i<nstepX;++i) {
    
    double valX = llimX + i*stepX;
    meX->SetValue(valX);
    
    for(int j=0;j<nstepY;++j) {
      
      double valY = llimY + j*stepY;
      if(i%2)
	valY = ulimY - j*stepY;
      meY->SetValue(valY);
      
      std::ofstream meFile((name + ".bst").c_str());
      for(MatrixElement* me1 : beam->GetMatrixElements())
	meFile << std::setprecision(19) << me1->GetValue() << "\n";
      meFile.close();
      
      if(data) {
	
	std::string cmd = "gosia < " + name + ".POIN.inp > /dev/null 2>&1";
	system(cmd.c_str());
	
	data->FillFromPoinOutput(); 
	
	if(inti) {
	  
	  cmd = "gosia < " + name + ".INTI.inp > /dev/null 2>&1";
	  system(cmd.c_str());
	  
	  data->FillFromIntiOutput();
	  data->Correct();
	  
	  //Recorrect target-linked experiments
	  for(std::array<int,3> info : targ_info) {
	    
	    int exp1B = info[0];
	    int numE = info[1];
	    int targ_num = info[2];
	    ExperimentalData* dataT = theFCN->target_data[targ_num-1];

	    int nrm_index = data->GetNormalizationIndex();
	    double scale0 = data->GetExperiment(exp1B)->GetAllIntiYields().at(nrm_index)->GetValue() / 
	      data->GetExperiment(exp1B)->GetAllPointYields().at(nrm_index)->GetValue();
	    
	    for(int k=0;k<numE;++k) { 
	      data->RecorrectExp(exp1B + k,scale0);
	      dataT->RecorrectExp(k+1,scale0);
	    }
	    
	  }
	}
      }
      
      mini->Minimize();

      double chi2 = mini->MinValue();
      std::cout << valX << " " << valY << " " << chi2 << std::endl;
      //Print();
      
      outFile << chi2;
      for(MatrixElement* me1 : beam->GetMatrixElements())
	outFile << " " << me1->GetValue();
      outFile << "\n";
      
    }
  }
  
  outFile.close();

  return;
}

void GosiaMinimizer::Minimize() {

  if(!theFCN->beam) {
    std::cout << "No beam nucleus set!" << std::endl;
    return;
  }
  if(!(theFCN->beam->CheckMatrixElements()))
    return;

  InitialSetup();
  
  int size = mini->NFree();
  ROOT::Math::Functor func(*theFCN,size);
  mini->SetFunction(func);
  
  std::cout << size << " free parameters" << std::endl;;
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

  std::vector<double> calc_scales; 
  if(calculated)
     calc_scales = theFCN->CalculateScalingParameters();
  
  Nucleus* beam = theFCN->beam;
  std::string beam_name = beam->GetName();
  
  //Write all matrix elements in the beam to file (a GOSIA matrix element file)
  std::ofstream meFileB((beam_name + ".bst-minimized").c_str());
  for(MatrixElement* me : beam->GetMatrixElements()) {
    double val = me->GetValue();
    meFileB << std::setprecision(19) << val << "\n";
  }
  meFileB.close();

  //Write all matrix elements in all the targets to file (a GOSIA matrix element file)
  for(int i=0;i<theFCN->targets.size();i++) {
    
    Nucleus* targ = theFCN->targets[i];
    
    std::ofstream meFileT((targ->GetName() + ".bst-minimized").c_str());
    for(MatrixElement* me : targ->GetMatrixElements()) {
      double val = me->GetValue();
      meFileT << std::setprecision(19) << val << "\n";
    }
    meFileT.close();

  }

  //Perform error validation if requested
  if(validate && status == 0) {
    
    std::cout << "Validating Errors..." << std::flush;
    mini->Hesse();
    
    status = mini->Status();
    if(status == 0)
      std::cout << " Done!" << std::endl;
    else
      std::cout << " Failed." << std::endl;
  }
  
  std::cout << "\n";
  if(calculated) { 

    int sizeBE = theFCN->beam_data->Size();
    int scale_num = 0;

    std::cout << "Final Scalings\n";
    for(int i=0;i<sizeBE;++i) {
      int num = theFCN->indices[i];
      if(num == scale_num) {
	std::cout << Form(" Scaling%02d: ",num) << calc_scales[num] << "\n";
	scale_num++;
      }
    }
    std::cout << std::endl;
  }
  
  Print();

  const double* min = mini->X();
  const double* errors = mini->Errors();
  if(validate && status == 0) {
    
    TMatrixD res;
    res.ResizeTo(2,size);
    for(int i=0;i<size;++i) {
      res[0][i] = min[i];
      res[1][i] = errors[i];
    }
    
    TMatrixD covM;
    TMatrixD corM;
    covM.ResizeTo(size,size);
    corM.ResizeTo(size,size);
    
    for(int i=0;i<size;++i) {
      for(int j=i;j<size;++j) {
	covM[i][j] = mini->CovMatrix(i,j);
	covM[j][i] = mini->CovMatrix(j,i);
	corM[i][j] = mini->Correlation(i,j);
	corM[j][i] = mini->Correlation(j,i);
      }
    }

    //Save minimum, uncertainties, and correlations to a ROOT file
    TFile* fileM = new TFile("ResultFile.root","RECREATE");
    fileM->cd();
    
    res.Write("Minimum");
    covM.Write("Covariance_Matrix");
    corM.Write("Correlation_Matrix");
    
    fileM->Close();
    delete fileM;

    //std::cout << "\n******** Covariance matrix **********";
    //covM.Print();
    
    //std::cout << "\n******** Correlation matrix **********";
    //corM.Print();
  }

  //Write Migrad (symmetric) uncertainties to a text file
  std::ofstream uncFile("Symmetric_Errors.txt");
  for(int i=0;i<size;++i) {
    
     double val = min[i];
     double err = errors[i];
     std::string name = mini->VariableName(i);
     
     uncFile << "Par " << i << " " << name << ": " << val << " +- " << err << "\n";

  }
  uncFile.close();
  
  //Calculate MINOS (asymmetric) uncertainties if requested and write them to a text file
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

      std::cout << "Par " << i << " " <<  name << ": " << val << " +" << eU << " " << eD << " (status = " << statusM << ")" << std::endl;
      minosFile << "Par " << i << " " <<  name << ": " << val << " +" << eU << " " << eD << " (status = " << statusM << ")" << "\n";
    }
    minosFile.close();
  }
  
  if(write)
    theFCN->Write();
 
  beam->FillFromBSTFile(beam_name + ".bst-minimized"); //Reset MEs to minimum  
  //if(theFCN->beam_data && !fixed && !calculated)
  //UpdateScalings(min);

  return;
}

/*
void GosiaMinimizer::UpdateScalings(const double* x) {
  
  int numM = theFCN->beam->GetNumFree();
  for(Nucleus* targ : theFCN->targets)
    numM += targ->GetNumFree();
  
  int numE_beam = theFCN->beam_data->Size();
  
  std::vector<double> scales;
  for(int i=0;i<numE_beam;++i)
    scales.push_back(theFCN->rel_scalings.at(i)*x[numM + theFCN->indices.at(i)]);
  
  theFCN->beam_data->SetScalings(scales);

  int numE_targs = 0;
  for(int i=0;i<theFCN->target_data.size();++i) {
    
    ExperimentalData* dataT = theFCN->target_data.at(i);
    int size = dataT->Size();
    
    std::vector<double> scales;
    for(int j=0;j<size;++j)
      scales.push_back(theFCN->rel_scalings.at(numE_beam + numE_targs + j)*x[numM + theFCN->indices.at(numE_beam + numE_targs + j)]);
    
    dataT->SetScalings(scales);
    
    numE_targs += size;
  }

  return;
}
*/
