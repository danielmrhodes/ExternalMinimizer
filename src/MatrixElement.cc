#include "MatrixElement.h"
#include <iostream>

MatrixElement::MatrixElement() {

  fixed = true;
  limited = false;
  
  index1 = -1; 
  index2 = -1; 
  mult = -1;

  value = 0.0;
  ulim = 4.0; 
  llim = -4.0;

}
MatrixElement::~MatrixElement() {;}

std::string MatrixElement::GetMultS(int lval) {

  switch(lval) {

  case 1:
    return "E1";

  case 2:
    return "E2";

  case 3:
    return "E3";
  
  case 4:
    return "E4";

  case 5:
    return "E5";

  case 6:
    return "E6";

  case 7:
    return "M1";

  case 8:
    return "M2";
    
  }

  std::cout << "Unknown multipolarity: " << lval << std::endl; 

  return "??";
}

std::string MatrixElement::GetMultS() const {

  switch(mult) {

  case 1:
    return "E1";

  case 2:
    return "E2";

  case 3:
    return "E3";
  
  case 4:
    return "E4";

  case 5:
    return "E5";

  case 6:
    return "E6";

  case 7:
    return "M1";

  case 8:
    return "M2";
    
  }

  std::cout << "Unknown multipolarity: " << mult << std::endl; 

  return "??";
}

bool MatrixElement::Compare(const MatrixElement* me1, const MatrixElement* me2) {

  if(me1->GetMultipolarity() != me2->GetMultipolarity())
   return me1->GetMultipolarity() < me2->GetMultipolarity();
  
  if(me1->GetIndex1() != me2->GetIndex1())
    return me1->GetIndex1() < me2->GetIndex1();
  
  return me1->GetIndex2() < me2->GetIndex2();
  
}
