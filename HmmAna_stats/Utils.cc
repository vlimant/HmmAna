#include "Utils.h"

using namespace std;
using namespace RooFit;

vector<TString> Utils::SplitString(const TString& theOpt, const char separator )
{
   vector<TString> splitV;
   TString splitOpt(theOpt);
   splitOpt.ReplaceAll("\n"," ");
   splitOpt = splitOpt.Strip(TString::kBoth,separator);
   while (splitOpt.Length()>0) {
      if ( !splitOpt.Contains(separator) ) {
         splitV.push_back(splitOpt);
         break;
      }
      else {
         TString toSave = splitOpt(0,splitOpt.First(separator));
         splitV.push_back(toSave);
         splitOpt = splitOpt(splitOpt.First(separator),splitOpt.Length());
      }
      splitOpt = splitOpt.Strip(TString::kLeading,separator);
   }
   return splitV;
}


void Utils::setDefaultMinimize(){
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
  ROOT::Math::MinimizerOptions::SetDefaultStrategy(1);
  ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(-1);
}

float Utils::minimizeTest(RooNLLVar* nll,float neps)
{
    nll->enableOffsetting(true);
    int printLevel = ROOT::Math::MinimizerOptions::DefaultPrintLevel();
    RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
    int strat = ROOT::Math::MinimizerOptions::DefaultStrategy();

    RooMinimizer minim(*nll);
    minim.setStrategy(strat);
    minim.setPrintLevel(1);
    minim.setProfile();
    minim.setEps(neps);
    int status = minim.minimize(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());

    return status;
}

float Utils::minimizeMinosTest(RooNLLVar* nll,RooRealVar* poi, float neps)
{
    nll->enableOffsetting(true);
    int printLevel = ROOT::Math::MinimizerOptions::DefaultPrintLevel();
    RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
    int strat = ROOT::Math::MinimizerOptions::DefaultStrategy();

    RooMinimizer minim(*nll);
    minim.setStrategy(strat);
    minim.setPrintLevel(1);
    minim.setProfile();
    minim.setEps(neps);
    minim.migrad(); 
    minim.minos(*poi);
    //pos_mean = ml->getVal();
    int status = minim.minimize(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());
    
    return status;
}

