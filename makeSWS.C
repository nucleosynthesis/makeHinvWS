void makeSWS(TString year, TString region){
  gROOT->ProcessLine(".L  ../makeSignalWS.C++");
  gROOT->ProcessLine("makeSignalWS(\"" + year + "\",\"" + region + "\")");
  gROOT->ProcessLine(".q");
}

