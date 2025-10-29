#include "Yield.h"

Yield::Yield() {
  yld = 0.0;
  thresh = 0.5;
}
Yield::~Yield() {;}

YieldError::YieldError() {
  Yield();
  erUp = 0.0;
  erDn = 0.0; 
  weight = 1.0;
}
YieldError::~YieldError() {;}
