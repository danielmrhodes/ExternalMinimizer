#include "MatrixElement.h"

MatrixElement::MatrixElement() {

  fixed = true;
  ulim = 4.0; 
  llim = -4.0;

}
MatrixElement::~MatrixElement() {;}

bool MatrixElement::Compare(const MatrixElement* me1, const MatrixElement* me2) {

  if(me1->GetMultipolarity() != me2->GetMultipolarity())
   return me1->GetMultipolarity() < me2->GetMultipolarity();
  
  if(me1->GetIndex1() != me2->GetIndex1())
    return me1->GetIndex1() < me2->GetIndex1();
  
  return me1->GetIndex2() < me2->GetIndex2();
  
}
