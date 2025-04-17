#include "GosiaMinimizerFCN.h"
#include "Kinematics.h"
#include "TMath.h"
#include <iostream>
#include <fstream>

GosiaMinimizerFCN::GosiaMinimizerFCN() {
  
  chi_cut = 1.0;
  fixed = false;

  beam_name = "";
  beam = NULL;
  beam_lit = NULL;
  beam_data = NULL;
  ntuple = NULL;
  
  target_names.clear();
  targets.clear();
  target_lits.clear();
  target_data.clear();

  indices.clear();
  scalings.clear();
  factors.clear();

  return;

}
GosiaMinimizerFCN::~GosiaMinimizerFCN() {;}

double GosiaMinimizerFCN::operator() (const double* pars) {

  int par_num = 0;
  
  if(beam->GetNumFree() > 0) {
      
    //Update free beam matrix elements, write all of them to file
    const int sizeB = beam->GetNumMatrixElements();
    std::ofstream meFileB((beam_name + ".bst").c_str());
    
    for(int i=0;i<sizeB;++i) {

      MatrixElement* me = beam->GetMatrixElement(i);
      if(me->IsFixed()) {
	
	meFileB << std::setprecision(19) << me->GetValue() << "\n";
	continue;
      }
    
      double val = pars[par_num];
      
      meFileB << std::setprecision(19) << val << "\n";
      me->SetValue(val);
      ++par_num;

    }
    meFileB.close();

    if(beam_data) { 

      //Calculate beam point yields with new matrix elements
      std::string cmdB = "gosia < " + beam_name + ".POIN.inp > /dev/null 2>&1";
      system(cmdB.c_str());
      
      //Read in the new point yields
      beam_data->FillFromPoinOutput(); 
    }
  }
  
  //Update free target matrix elements, write all of them to file
  for(int i=0;i<targets.size();++i) {

    Nucleus* targ = targets[i];
    if(targ->GetNumFree() == 0)
      continue;
    
    const int sizeT = targ->GetNumMatrixElements();
    std::ofstream meFileT((target_names[i] + ".bst").c_str());
    
    for(int j=0;j<sizeT;++j) {
      
      MatrixElement* me = targ->GetMatrixElement(j);
      if(me->IsFixed()) {
	meFileT << std::setprecision(19) << me->GetValue() << "\n";
	continue;
      }
      
      double val = pars[par_num];
      meFileT << std::setprecision(19) << val << "\n";
      me->SetValue(val);
      ++par_num;
      
    }
    meFileT.close();

    if(i < target_data.size()) { 
      
      //Calculate target point yields with new matrix elements
      std::string cmdT = "gosia < " + target_names[i] + ".POIN.inp > /dev/null 2>&1";
      system(cmdT.c_str());
      
      //Read in the new point yields
      target_data[i]->FillFromPoinOutput();
    }
  }

  //Calculate yield chi2 for beam
  double yld_chi2_beam = 0.0;
  int numE_beam = 0;
  if(beam_data) {
    numE_beam += beam_data->Size();
    
    for(int i=0;i<numE_beam;++i) {
      
      Experiment* exp = beam_data->data[i];
    
      double rel = scalings[i];
      int num = indices[i];
      double fac = factors[i];
      
      double scale;
      if(fixed)
	scale = rel*fixed_scales[num];
      else
	scale = rel*pars[par_num + num];    

      for(int j=0;j<exp->data_corr.size();++j) {
      
	YieldError* yldC = exp->data_corr[j];
	Yield* yldP = exp->point_yields[j];  
      
	double diff = yldC->GetValue() - scale*yldP->GetValue()/fac;
      
	double nSig = 0.0;
	if(diff > 0.0)
	  nSig = diff/yldC->GetErrorDown();
	else
	  nSig = diff/yldC->GetErrorUp();
      
	yld_chi2_beam += nSig*nSig*yldC->GetWeight();

      }
    }
  }

  //Calculate yield chi2 for targets
  double yld_chi2_targs = 0.0;
  int numE_targs = 0;
  for(int i=0;i<target_data.size();++i) {
    
    ExperimentalData* dataT = target_data[i];
    int size = dataT->Size();

    for(int j=0;j<size;++j) {
      
      Experiment* exp = dataT->data[j];
    
      double rel = scalings[numE_beam + numE_targs + j];
      int num = indices[numE_beam + numE_targs + j];
      double fac = factors[numE_beam + numE_targs + j];
      
      double scale;
      if(fixed)
	scale = rel*fixed_scales[num];
      else
	scale = rel*pars[par_num + num]; 

      for(int k=0;k<exp->data_corr.size();++k) {
      
	YieldError* yldC = exp->data_corr[k];
	Yield* yldP = exp->point_yields[k];  
      
	double diff = yldC->GetValue() - scale*yldP->GetValue()/fac;
      
	double nSig = 0.0;
	if(diff > 0.0)
	  nSig = diff/yldC->GetErrorDown();
	else
	  nSig = diff/yldC->GetErrorUp();
      
	yld_chi2_targs += nSig*nSig*yldC->GetWeight();

      }
    }
    numE_targs += size;
  }
  
  double lit_chi2_beam = 0.0;
  if(beam_lit)
    lit_chi2_beam += beam->CompareWithLiterature(beam_lit);
  
  double lit_chi2_targs = 0.0;
  for(int i=0;i<targets.size();++i)
    if(i < target_lits.size())
      lit_chi2_targs += targets[i]->CompareWithLiterature(target_lits[i]);
  
  if(ntuple) {
    
    std::vector<float> output = beam->GetMatrixElementValues();
    output.insert(output.begin(),lit_chi2_targs);
    output.insert(output.begin(),yld_chi2_targs);
    output.insert(output.begin(),lit_chi2_beam);
    output.insert(output.begin(),yld_chi2_beam);
    
    ntuple->Fill(&output[0]);
  }
  
  //std::cout << yld_chi2_beam << std::endl;
  return yld_chi2_beam + lit_chi2_beam + yld_chi2_targs + lit_chi2_targs;

}

void GosiaMinimizerFCN::CreateNTuple() {

  if(beam->GetNumFree() == 0) {
    std::cout << "No free matrix elements in the beam. Will not write parameters space file." 
	      << std::endl;

      return;
  }

  std::string varlist = "yldChi2B:litChi2B:yldChi2T:litChi2T";
  
  int size = beam->GetNumMatrixElements();
  for(int i=1;i<size+1;++i)
    varlist += Form(":ME%d",i);
  
  ntuple = new TNtuple("ps","Parameter Space",varlist.c_str());

  return;
}

void GosiaMinimizerFCN::Write() {

  TFile* outFile = new TFile("ParameterSpace.root","RECREATE");
  
  ntuple->Write();
  outFile->Close();
  
  delete outFile;

  return;
}

void GosiaMinimizerFCN::FillFactors() {  
  
  if(!beam_data) 
    return; 
  
  factors.clear();
  Kinematics kin;
    
  double beam_ex = beam->GetLevels().at(1)->GetEnergy();
  for(Experiment* exp : beam_data->data) {
    
    kin.SetBeamZ(exp->beam_Z);
    kin.SetBeamMass(exp->beam_mass);
    kin.SetTargetZ(exp->targ_Z);
    kin.SetTargetMass(exp->targ_mass);
    kin.SetBeamEnergy(exp->Ep);
    kin.SetEx(beam_ex);

    bool target_det = false;
    double thc = exp->thP;
    if(thc < 0.0) {
      target_det = true;
      thc *= -1.0;
    }

    double thd = thc;
    if(target_det)
      thd = kin.Recoil_Theta_LAB(kin.Theta_CM_FP(thc));

    double fac = TMath::Sin(thd)*kin.Dsig(thc,false,target_det);
    factors.push_back(fac);
  }
  beam_data->SetFactors(factors);
  
  for(int i=0;i<target_data.size();++i) {
    
    ExperimentalData* dataT = target_data.at(i);
    double targ_ex = targets.at(i)->GetLevels()[1]->GetEnergy();
  
    std::vector<double> tmp_facs;
    for(Experiment* exp : dataT->data) {

      kin.SetBeamZ(exp->beam_Z);
      kin.SetBeamMass(exp->beam_mass);
      kin.SetTargetZ(exp->targ_Z);
      kin.SetTargetMass(exp->targ_mass);
      kin.SetBeamEnergy(exp->Ep);
      kin.SetEx(targ_ex);

      bool target_det = false;
      double thc = exp->thP;
      if(thc < 0.0) {
	target_det = true;
	thc *= -1.0;
      }

      double thd = thc;
      if(target_det)
	thd = kin.Recoil_Theta_LAB(kin.Theta_CM_FP(thc));

      double fac = TMath::Sin(thd)*kin.Dsig(thc,false,target_det);
      factors.push_back(fac);
      tmp_facs.push_back(fac);

    }
    dataT->SetFactors(tmp_facs);
  }

  return;
}
