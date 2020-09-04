void makeallWS(TString year, TString region){
  gROOT->ProcessLine(".L  ../makeWS_percategory.C++");
  gROOT->ProcessLine("makeWS_percategory(\"" + year + "\",\"" + region + "\")");
  gROOT->ProcessLine(".q");
}
