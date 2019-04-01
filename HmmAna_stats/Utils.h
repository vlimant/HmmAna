#ifndef _Utils_h_
#define _Utils_h_
#include "CommonHead.h"
#include "RooFitHead.h"

using namespace std;
using namespace RooFit;

namespace Utils{

vector<TString> SplitString(const TString& theOpt, const char separator );

void setDefaultMinimize();

float minimizeTest(RooNLLVar* nll,float neps);

float minimizeMinosTest(RooNLLVar* nll,RooRealVar* poi,float neps);

};

#endif
