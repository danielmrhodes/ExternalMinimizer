#ifndef BrIccReader_h
#define BrIccReader_h 1

#include <iostream>

class BrIccReader {
  
  //these are structures for reading the BrIcc database format
  struct eshell {
  public:
    bool exist;
    float be;
    char tab[8];
    char shell[8];
    int nmesh;
    int nrec;

    void print() {
      std::cout << "Exists = " << exist << std::endl;
      std::cout << "Binding energy: " << be << std::endl;
      std::cout << "Tabulation    : ";
      for (int i=0; i<8; ++i) { std::cout << tab[i]; }
      std::cout << std::endl;
      std::cout << "Electron shell: ";
      for (int i=0; i<8; ++i) { std::cout << shell[i]; }
      std::cout << std::endl;
      std::cout << "Number of mesh: " << nmesh << std::endl;
      std::cout << "Record offset : " << nrec << std::endl;
    }    
  };

  struct elem {
  public:
    int Z;
    char symb[4];
    int A;

    void print() {
      std::cout << "Z : " << Z << std::endl;
      std::cout << "Symbol : ";
      for (int i=0; i<4; ++i) { std::cout << symb[i]; }
      std::cout << std::endl;
      std::cout << "A : " << A << std::endl;
    }
  };

  struct icc {
  public:
    float energy;
    float icc[10];  
  };
  
 public:
 BrIccReader() : idxfile(NULL), iccfile(NULL) { };
 BrIccReader(std::string idxpath, std::string iccpath);
 ~BrIccReader() { Close(); }  

  void Open(std::string idxpath, std::string iccpath);
  void Close();

  double GetTotalCC(int Z, double egamma, int mult); //note that this mult is defined differently from GOSIAReader mult
  double GetTotalOmg(int Z, double egamma); //Gets Omega for E0 transitions

  static void Spline(double *x, double *y, int n, double xval, double &yval);
  static void SplineEval(double *x, double *y, double *ddy, int n, double &xval, double &yval);
  static void SplineFit(double *x, double *y, int n, double dy1, double dyn, double* ddy);

 private:
  FILE *idxfile;
  FILE *iccfile;
   
};

#endif
