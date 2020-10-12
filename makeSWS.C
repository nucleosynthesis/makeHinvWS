void makeSWS(TString year, TString region){
  gROOT->ProcessLine(".L  ../makeSignalAndMCBackgroundWS.C++");
  gROOT->ProcessLine("makeSignalAndMCBackgroundWS(\"" + year + "\",\"" + region + "\")");
  gROOT->ProcessLine(".q");
}

