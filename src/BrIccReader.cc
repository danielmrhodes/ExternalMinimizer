#include <iostream>
#include <cmath>
#include <cstdio>
#include <vector>

#include "BrIccReader.h"

BrIccReader::BrIccReader(std::string idxpath, std::string iccpath) {
  idxfile = fopen(idxpath.c_str(), "r");
  iccfile = fopen(iccpath.c_str(), "r");
  //std::cout << "Conversion electron coefficients from BrIcc database. \nPlease cite T. Kibedi et al. NIM A589 (2009) 202-229 for the conversion coefficients!" << std::endl;
}

void BrIccReader::Open(std::string idxpath, std::string iccpath) {
  idxfile = fopen(idxpath.c_str(), "r");
  iccfile = fopen(iccpath.c_str(), "r");
  //std::cout << "Conversion electron coefficients from BrIcc database. \nPlease cite T. Kibedi et al. NIM A589 (2009) 202-229 for the conversion coefficients!" << std::endl;
}

void BrIccReader::Close() {
  if (idxfile) {
    fclose(idxfile);
  }
  if (iccfile) {
    fclose(iccfile);
  }
  idxfile=NULL;
  iccfile=NULL;
}

double BrIccReader::GetTotalCC(int Z, double egamma, int mult) { //note that this mult is defined differently from GOSIAReader mult
  if (!idxfile || !iccfile) { std::cout << "Open BrIcc data base before getting conversion coefficients!" << std::endl; return -1; }
  fseek(idxfile, 2048*(Z-1), SEEK_SET);
  elem e;
  fread(&e, 12, 1, idxfile);

  if (egamma > 6000) { return 0; }

  std::vector<BrIccReader::eshell> shell;
  for (int i=0; i<37; ++i) {
    eshell rec;
    fread(&rec, 32, 1, idxfile);
    if ((int)rec.exist > 0) { rec.exist = true; }
    else { rec.exist = false; }
    shell.push_back(rec);
  }

  if (egamma < 1100) { 
    shell[36].exist = false; //no IPC below 1100 keV
  }
  if ((mult > 2 && mult < 5) || (mult > 7 && mult < 10)) {
    shell[36].exist = false; //no IPC for E4,E5,M4,M5
  }

  float energy[37][200];
  float cc[37][200][10];
  for (int i = 0; i < 37; ++i) {
    if (egamma < shell[i].be) { shell[i].exist = false; }
    if (!shell[i].exist) { continue; }
    int iEn = 0;
    for (int rec = shell[i].nrec - 1; rec < shell[i].nrec + shell[i].nmesh - 1; ++rec) {
      fseek(iccfile, rec*44, SEEK_SET);
      icc coef;
      fread(&coef, 44, 1, iccfile);
      energy[i][iEn] = coef.energy;
      if (iEn > 200) { std::cout << "SEVERE ERROR! " << iEn << std::endl; }
      for (int imult=0; imult<10; ++imult) { cc[i][iEn][imult] = coef.icc[imult]; }
      ++iEn;
    }        
  }

  double total_cc = 0;
  for (int i=0; i<37; ++i) {
    double this_cc = 0;
    if (!shell[i].exist) { continue; }
    if (egamma < energy[i][0]) {
      this_cc = cc[i][0][mult];
      std::cout << "Warning! Egamma = " << egamma << " is in the regime where solid state effects dominate, conversion coefficients for shell ";
      for (int j=0; j<8; ++j) { std::cout << shell[i].shell[j]; }
      std::cout << std::endl;      
    }
    else if (egamma > energy[i][shell[i].nmesh-1]) {
      this_cc = cc[i][shell[i].nmesh-1][mult];
      //std::cout << "Warning! Egamma = " << egamma << " exceeds the range of conversion coefficients table for shell ";
      //for (int j=0; j<8; ++j) { std::cout << shell[i].shell[j]; }
      //std::cout << std::endl;      
    }
    else {
      int n = shell[i].nmesh;
      double *x = new double[n]();
      double *y = new double[n]();

      for (int iEn = 0; iEn < n; ++iEn) {
        x[iEn] = std::log(energy[i][iEn]);
        y[iEn] = cc[i][iEn][mult];
      }
      BrIccReader::Spline(x,y,n,std::log(egamma),this_cc);

      delete[] x;
      delete[] y;
    }    
    total_cc += this_cc;
  }

  return total_cc;
 
}

double BrIccReader::GetTotalOmg(int Z, double egamma) { //Gets Omega for E0 transitions
  if (!idxfile || !iccfile) { std::cout << "Open BrIcc data base before getting conversion coefficients!" << std::endl; return -1; }
  fseek(idxfile, 2048*(Z-1), SEEK_SET);
  elem e;
  fread(&e, 12, 1, idxfile);

  if (egamma > 6000) { return 0; }

  std::vector<BrIccReader::eshell> shell;
  for (int i=0; i<41; ++i) {
    eshell rec;
    fread(&rec, 32, 1, idxfile);
    if ((int)rec.exist > 0) { rec.exist = true; }
    else { rec.exist = false; }
    shell.push_back(rec);
  }

  if (egamma < 1100) { 
    shell[40].exist = false; //no IPC below 1100 keV
  }

  float energy[4][200];
  float omg[4][200];
  for (int i = 37; i < 41; ++i) {
    if (egamma < shell[i].be) { shell[i].exist = false; }
    if (!shell[i].exist) { continue; }
    int iEn = 0;
    for (int rec = shell[i].nrec - 1; rec < shell[i].nrec + shell[i].nmesh - 1; ++rec) {
      fseek(iccfile, rec*44, SEEK_SET);
      icc coef;
      fread(&coef, 44, 1, iccfile);
      energy[i-37][iEn] = coef.energy;
      if (iEn > 200) { std::cout << "SEVERE ERROR! " << iEn << std::endl; }
      omg[i-37][iEn] = coef.icc[0];
      ++iEn;
    }        
  }

  double total_omg = 0;
  for (int i=37; i<41; ++i) {
    double this_omg = 0;
    if (!shell[i].exist) { continue; }
    if (egamma < energy[i-37][0]) {
      this_omg = omg[i-37][0];
      std::cout << "Warning! Egamma = " << egamma << " is in the regime where solid state effects dominate, conversion coefficients for shell ";
      for (int j=0; j<8; ++j) { std::cout << shell[i].shell[j]; }
      std::cout << std::endl;      
    }
    else if (egamma > energy[i-37][shell[i].nmesh-1]) {
      this_omg = omg[i-37][shell[i].nmesh-1];
      std::cout << "Warning! Egamma = " << egamma << " exceeds the range of conversion coefficients table for shell ";
      for (int j=0; j<8; ++j) { std::cout << shell[i].shell[j]; }
      std::cout << std::endl;      
    }
    else {
      int n = shell[i].nmesh;
      double *x = new double[n]();
      double *y = new double[n]();

      for (int iEn = 0; iEn < n; ++iEn) {
        x[iEn] = std::log(energy[i-37][iEn]);
        y[iEn] = omg[i-37][iEn];
      }
      BrIccReader::Spline(x,y,n,std::log(egamma),this_omg);

      delete[] x;
      delete[] y;
    }    
    total_omg += this_omg;
  }

  return total_omg;
 
}

void BrIccReader::SplineFit(double *x, double *y, int n, double dy1, double dyn, double* ddy) {
  double y2[n];
  double u[1500];
  double qn;
  double un;

  if (dy1 > 1e30) {
    y2[0] = 0;
    u[0] = 1500;
  }
  else {
    y2[0] = -0.5;
    u[0] = (3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-dy1);
  }
  for (int i=1; i<n-1; ++i) {
    double sig = (x[i] - x[i-1])/(x[i+1]-x[i-1]);
    double p = sig * ddy[i-1] + 2;
    y2[i] = (sig-1)/p;
    u[i] = (6 * ((y[i+1]-y[i])/(x[i+1] - x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]))/(x[i+1]-x[i-1]) - sig*u[i-1])/p;    
  }
  if (dyn > 1e30) {
    qn = 0;
    un = 0;
  }
  else {
    qn = 0.5;
    un = (3/(x[n-1]-x[n-2])) * (dyn - (y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
  }
  ddy[n-1] = (un-qn*u[n-2])/(qn*ddy[n-2] + 1.);
  for (int k=n-2; k>=0; --k) {
    ddy[k] = ddy[k]*ddy[k+1] + u[k];
  }
}

void BrIccReader::SplineEval(double *x, double *y, double *ddy, int n, double &xval, double &yval) {
  int klo = 0;
  int khi = n-1;
  int k;
  while (khi-klo > 1) {
    k = (khi + klo)/2;
    if (x[k] > xval) { khi = k; }
    else { klo = k; }
  }
  double h = x[khi] - x[klo];

  double a;
  double b;
  if (std::abs(h) < 1e-9) {
    if ( xval < x[0] ) { h = 1; }
    if ( xval > x[n-1] ) { k = n-2; }
    a = (y[k] - y[k+1])/(x[k] - x[k+1]);
    b = (y[k] + y[k+1] - a*(x[k] + x[k+1]))*0.5;
    yval = a*xval + b;
    std::cout << "SplineEval " << xval << " " << yval << " extrapolation" << std::endl;
    return;
  }
  a = (x[khi] - xval)/h;
  b = (xval - x[klo])/h;
  yval = a*y[klo] + b*y[khi] + ((std::pow(a,3) - a)*ddy[klo] + (std::pow(b,3) - b)*ddy[khi])*(h*h/6.);
                                                     
  return;
}

void BrIccReader::Spline(double *x, double *y, int n, double xval, double &yval) {
  double *y_loc = new double[n]();
  for (int i=0; i<n; ++i) {
    y_loc[i] = std::log(y[i]);
  }
  double dy1 = (y_loc[1]-y_loc[0])/(x[1]-x[0]);
  double dyn = (y_loc[n-1]-y_loc[n-2])/(x[n-1]-x[n-2]);

  double *ddy = new double[n]();
  SplineFit(x, y_loc, n, dy1, dyn, ddy);

  SplineEval(x, y_loc, ddy, n, xval, yval);

  yval = std::exp(yval);

  delete[] y_loc;
  delete[] ddy;
  return;
}
