#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "TSystem.h"
#include "TFile.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooArgList.h"
#include "TH1F.h"
#include "RooDataHist.h"
#include "RooFormulaVar.h"
#include "/afs/cern.ch/user/v/vmilosev/CMSSW_10_2_13/src/HiggsAnalysis/CombinedLimit/interface/RooParametricHist.h"
#include "RooAddition.h"

enum PROCESS{
  data = 0,
  VBFH = 1,
  ggH = 2,
  QCDZnunu = 3,
  EWKZnunu = 4,
  QCDW = 5,
  EWKW = 6,
  QCDDYll = 7,
  EWKZll = 8
};


int makeWS(){
    // As usual, load the combine library to get access to the RooParametricHist
    gSystem->Load("libHiggsAnalysisCombinedLimit.so");

    ///////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////
    //Hardcoded input parameters
    //To be adjusted to channel !

    std::string lChannel = "VBF";
    std::string lOutFileName = "param_ws_"+lChannel+".root";

    //define the variable that is fitted
    RooRealVar lVarFit("mjj","M_{jj} (GeV)",200,5000);
    std::string lVarLabel = "Mjj";

    // Regions
    // Use same order! SR=region 0.
    const unsigned nR = 5;
    std::string lRegions[5] = {"SR","Wenu","Wmunu","Zee","Zmumu"};
    //processes
    //use same order: data= process 0, signal = process 1, QCD Z+Jets in SR = 2, etc....
    const unsigned nP = 12;
    std::string lProcs[nP] = {"data_obs","VBFHtoInv","GluGluHtoInv",
			      "ZJETS","EWKZNUNU",
			      "WJETS","EWKW",
			      "DY","EWKZll",
			      "TOP","VV","QCD"};

    const unsigned nN = 12;
    std::string lNuis[12] = {"bjet_veto","pileup","tau_veto",
			     "eventVetoVEleIdIso","eventVetoLMuId","eventVetoLMuIso",
			     "eventSelTEleIdIso","eventSelTMuId","eventSelTMuIso",
      			     "eventSelVEleIdIso","eventSelLMuId","eventSelLMuIso"
    };

    const unsigned nS = 2*nN+1;
    std::string lSysts[25];
    for (unsigned iS(0); iS<nS; ++iS){
      if (iS==0) lSysts[iS] = "";
      else lSysts[iS] = (iS-1)%2==0? lNuis[(iS-1)/2]+"Up" : lNuis[(iS-1)/2]+"Down";
      //std::cout << lSysts[iS] << std::endl;
    }
    const bool isSRsyst[25] = {1,1,1,1,1,
			       1,1,1,1,1,
			       1,1,1,0,0,
			       0,0,0,0,0,
			       0,0,0,0,0};
    const bool isCRWsyst[25] = {1,1,1,1,1,
				1,1,0,0,0,
				0,0,0,1,1,
				1,1,1,1,0,
				0,0,0,0,0};
    const bool isCRZsyst[25] = {1,1,1,1,1,
				1,1,0,0,0,
				0,0,0,1,1,
				1,1,1,1,1,
				1,1,1,1,1};
			       


    //input file path and name
    //input file is expected to contain one directory per region with names as in lRegions,
    //and one histogram per process with name as in lProcs with shape of the variable that is fitted.
    std::string lInFileName[nR];
    for (unsigned iR(0); iR<nR; ++iR){     
      lInFileName[iR] = "out_VBF_ana_"+lRegions[iR]+"_2017_vtest_jpt2_40/VBF_shapes.root";
    }
    
    //indices of QCD or EWK Z/W processes in SR/CR in lProcs array. 
    const unsigned nT = 2;
    std::string lType[nT] = {"QCD","EWK"};
    unsigned vproc[nT][5] = {
      {PROCESS::QCDZnunu,PROCESS::QCDW,PROCESS::QCDW,PROCESS::QCDDYll,PROCESS::QCDDYll},
      {PROCESS::EWKZnunu,PROCESS::EWKW,PROCESS::EWKW,PROCESS::EWKZll,PROCESS::EWKZll}
    };


    //values of nuisances
    double WZratioSyst = 1.12;

    double jesWWSyst = 1.02;
    double jesZZSyst = 1.01;

    /////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////
    // Output file and workspace

    TFile *fOut = new TFile(lOutFileName.c_str(),"RECREATE");
    RooWorkspace wspace("wspace","wspace");

    RooArgList vars(lVarFit);

    //-- Define the shapes and binning -- read from input files

    TFile *finput[nR];
    TH1F *histos[nR][nP][nS];
    for (unsigned iR(0); iR<nR; ++iR){
      for (unsigned iP(0); iP<nP; ++iP){
	for (unsigned iS(0); iS < nS; ++iS){
	  histos[iR][iP][iS] = 0;
	}
      }
    }

    for (unsigned iR(0); iR<nR; ++iR){
      
      finput[iR] = TFile::Open(lInFileName[iR].c_str());
      finput[iR]->cd(lRegions[iR].c_str());
      
      for (unsigned iP(0); iP<nP; ++iP){
	for (unsigned iS(0); iS < nS; ++iS){
	  if (iS==0 || 
	      (iR==0 && !isSRsyst[iS]) ||
	      (iR>0 && iR<3 && !isCRWsyst[iS]) ||
	      (iR>2 && !isCRZsyst[iS])
	      ) histos[iR][iP][iS] = (TH1F*)gDirectory->Get(lProcs[iP].c_str());
	  else if (iP>0) {
	    histos[iR][iP][iS] = (TH1F*)gDirectory->Get((lProcs[iP]+"_"+lSysts[iS]).c_str());
	  }
	  if (!histos[iR][iP][iS]) {
	    if (iP>0){
	      std::cout<< " Histo not found for region " << lRegions[iR] << " process " << lProcs[iP] << " syst " << lSysts[iS] << std::endl;
	      histos[iR][iP][iS] = (TH1F*)gDirectory->Get(lProcs[iP].c_str());
	    }
	    if (iS==0) return 1;
	    continue;
	  }
	  std::cout << " --- histos " << histos[iR][iP][iS]->GetName() << " " << histos[iR][iP][iS]->GetEntries() << " " << histos[iR][iP][iS]->Integral() << std::endl;
	}
      }
    }

    const unsigned nB = histos[0][PROCESS::data][0]->GetNbinsX();//-4;
    double bins[nB+1];
    for (unsigned iB(0); iB<nB; ++iB){
      bins[iB] = histos[0][PROCESS::data][0]->GetXaxis()->GetBinLowEdge(iB+1);//+2);
      //std::cout << " bin " << iB << " Low edge " <<  bins[iB] << std::endl;
    }
    bins[nB] = histos[0][PROCESS::data][0]->GetXaxis()->GetBinLowEdge(nB+1);
    //std::cout << " bin " << nB << " Low edge " <<  bins[nB] << std::endl;

    TH1F dummyHist("dummyHist","Dummy hist for binning",nB,bins);

    std::cout << nB << " - Check Binning: ";
    for (unsigned iB(0); iB<nB+1; ++iB){
      std::cout << bins[iB] << " ";
    }
    std::cout << std::endl;

    RooDataHist *data_hist[nR];
    RooDataHist *bkg_hist[nR][nP][nS];
    for (unsigned iR(0); iR<nR; ++iR){
      
      data_hist[iR] = new RooDataHist(("data_obs_"+lRegions[iR]).c_str(),"Data observed",vars,histos[iR][PROCESS::data][0]);
      wspace.import(*data_hist[iR]);
      
      for (unsigned iP(1); iP<nP; ++iP){
	for (unsigned iS(0); iS < nS; ++iS){
	  std::ostringstream label;
	  label << lProcs[iP] << "_hist_" << lRegions[iR];
	  if (iS>0) {
	    label << "_" << lSysts[iS];
	  }
	  if (histos[iR][iP][iS] && histos[iR][iP][iS]->GetEntries()>0){
	    bkg_hist[iR][iP][iS] = new RooDataHist(label.str().c_str(),"Background",vars,histos[iR][iP][iS]);
	    wspace.import(*bkg_hist[iR][iP][iS]);
	  }
	}
      }
    }

    // In the signal region, W and Z are tied together by their ratio from MC + theory uncertainty.
    // SR/CR remaining non-cancellations from JES uncertainties.
    RooRealVar *wzratio[nT];
    RooRealVar *jesWW[nT];
    RooRealVar *jesZZ[nT];
    for (unsigned iT(0); iT<nT; ++iT){     
      wzratio[iT] = new RooRealVar(("wzratio"+lType[iT]).c_str(), (lType[iT]+" W/Z ratio nuisance parameter").c_str(),0);
      jesWW[iT] = new RooRealVar(("jesWW"+lType[iT]).c_str(), (lType[iT]+" W/W JES nuisance parameter").c_str(),0);
      jesZZ[iT] = new RooRealVar(("jesZZ"+lType[iT]).c_str(), (lType[iT]+" Z/Z JES nuisance parameter").c_str(),0);
    }


    RooRealVar *ewkqcdratiostat[nB];
    RooRealVar *wzratiostat[nT][nB];
    RooRealVar *TFstat[nT][nB][nR-1];
    RooRealVar *TFsysts[nN];
    //same nuisance name for different CR, correlated accross CR or given different systs names already....??
    for (unsigned iN(0); iN < nN; ++iN){
      std::ostringstream lname;
      lname.str("");
      lname << "TF_syst_" << lNuis[iN];
      TFsysts[iN] = new RooRealVar(lname.str().c_str(),"CR/SR ratio syst nuisance parameter",0);
    }


    // Create one parameter per bin representing the yield. (note of course we can have multiple processes like this)
    for (unsigned iB(1); iB<nB+1; ++iB){
      std::cout << " -- Processing bin " << iB << std::endl;
      std::ostringstream lname;

      RooFormulaVar *EWKQCDbin = 0;
      
      for (unsigned iT(0); iT<nT; ++iT){     

	std::cout << " --- Processing type " << lType[iT] << std::endl;

	lname.str("");
	lname << lType[iT] << "Z_SR_bin" << iB;
	RooRealVar binParZ(lname.str().c_str(),(lType[iT]+" Z+jets yield in signal region, per bin").c_str(),histos[0][iT==0?PROCESS::QCDZnunu:PROCESS::EWKZnunu][0]->GetBinContent(iB),0,10*histos[0][iT==0?PROCESS::QCDZnunu:PROCESS::EWKZnunu][0]->GetBinContent(iB));


	if (iT==0) {
	  wspace.import(binParZ,RooFit::RecycleConflictNodes());
	  
	  //program link between QCD and EWK yields in SR
	  lname.str("");
	  lname << "ewkqcdratio_stat_bin" << iB;
	  ewkqcdratiostat[iB] = new RooRealVar(lname.str().c_str()," EWK/QCD ratio stat nuisance parameter",0);
	  lname.str("");
	  lname << "TF_EWKQCDSR_bin" << iB;
	  std::ostringstream lFormula;
	  double ratio = histos[0][PROCESS::EWKZnunu][0]->GetBinContent(iB) / histos[0][PROCESS::QCDZnunu][0]->GetBinContent(iB); 
	  double ratiostat = 1+sqrt(pow(histos[0][PROCESS::EWKZnunu][0]->GetBinError(iB)/histos[0][PROCESS::EWKZnunu][0]->GetBinContent(iB),2)+pow(histos[0][PROCESS::QCDZnunu][0]->GetBinError(iB)/histos[0][PROCESS::QCDZnunu][0]->GetBinContent(iB),2));
	  lFormula << ratio;
	  lFormula << "*TMath::Power(" << ratiostat << ",@0)";
	  
	  RooFormulaVar TFEWKQCD(lname.str().c_str(),"Transfer factor EWK/QCD Z",lFormula.str().c_str(),RooArgList(*(ewkqcdratiostat[iB])) );
	  wspace.import(TFEWKQCD,RooFit::RecycleConflictNodes());
	  lname.str("");
	  lname << "EWKQCD_SR_bin" << iB;
	  EWKQCDbin = new RooFormulaVar(lname.str().c_str(),"EWK Z+jets yield in signal regions from QCD Z yield, per bin","@0*@1",RooArgList(TFEWKQCD,binParZ));
	  wspace.import((*EWKQCDbin),RooFit::RecycleConflictNodes());
	}

	//lname.str("");
	//lname << "EWKQCD_SR_bin" << iB;
	//RooFormulaVar EWKQCDbin = *wspace.function(lname.str().c_str());
	lname.str("");
	lname << lType[iT] << "wzratio_stat_bin" << iB;
	wzratiostat[iT][iB] = new RooRealVar(lname.str().c_str(),"W/Z ratio stat nuisance parameter",0);
	lname.str("");
	lname << lType[iT] << "TF_WZSR_bin" << iB;
	std::ostringstream lFormula;
	lFormula.str("");
	double ratio = 0;
	double ratiostat = 0;
	if (iT==0) {
	  ratio = histos[0][PROCESS::QCDW][0]->GetBinContent(iB) / histos[0][PROCESS::QCDZnunu][0]->GetBinContent(iB);
	  ratiostat = 1+sqrt(pow(histos[0][PROCESS::QCDW][0]->GetBinError(iB)/histos[0][PROCESS::QCDW][0]->GetBinContent(iB),2)+pow(histos[0][PROCESS::QCDZnunu][0]->GetBinError(iB)/histos[0][PROCESS::QCDZnunu][0]->GetBinContent(iB),2));
	} else{
	  ratio = histos[0][PROCESS::EWKW][0]->GetBinContent(iB) / histos[0][PROCESS::EWKZnunu][0]->GetBinContent(iB);
	  ratiostat = 1+sqrt(pow(histos[0][PROCESS::EWKW][0]->GetBinError(iB)/histos[0][PROCESS::EWKW][0]->GetBinContent(iB),2)+pow(histos[0][PROCESS::EWKZnunu][0]->GetBinError(iB)/histos[0][PROCESS::EWKZnunu][0]->GetBinContent(iB),2));
	}
	lFormula << ratio;
	lFormula << "*TMath::Power(" << WZratioSyst << ",@0)*TMath::Power(" << ratiostat << ",@1)";

	std::cout << " ---- Check stat error " << lType[iT] << " WZratio: " << ratiostat << std::endl;
	
	RooFormulaVar TFWZ(lname.str().c_str(),"Transfer factor W/Z",lFormula.str().c_str(),RooArgList(*(wzratio[iT]),*(wzratiostat[iT][iB])) );
	wspace.import(TFWZ,RooFit::RecycleConflictNodes());

	lname.str("");
	lname << lType[iT] << "WZ_SR_bin" << iB;
	RooFormulaVar WZbin(lname.str().c_str(),(lType[iT]+" W+jets yield in signal regions from Z yield, per bin").c_str(),"@0*@1",iT==0?RooArgList(TFWZ,binParZ):RooArgList(TFWZ,*(EWKQCDbin)));
	wspace.import(WZbin,RooFit::RecycleConflictNodes());
	
	if (iT==0) std::cout << " --- SR QCD Z yield = " << histos[0][PROCESS::QCDZnunu][0]->GetBinContent(iB) << " SR QCD W yield = " << histos[0][PROCESS::QCDW][0]->GetBinContent(iB) << " SR signal yield = " << histos[0][PROCESS::VBFH][0]->GetBinContent(iB) << std::endl;
	else std::cout << " --- SR EWK Z yield = " << histos[0][PROCESS::EWKZnunu][0]->GetBinContent(iB) << " SR EWK W yield = " << histos[0][PROCESS::EWKW][0]->GetBinContent(iB) << " SR signal yield = " << histos[0][PROCESS::VBFH][0]->GetBinContent(iB) << std::endl;
	
	for (unsigned iR(1); iR<nR; ++iR){
	  
	  std::cout << " ---- Doing control regions: " << lRegions[iR] << std::endl;
	  
	  lname.str("");
	  lname << lType[iT] << "TF_" << lRegions[iR] << "_stat_bin" << iB;
	  TFstat[iT][iB][iR-1] = new RooRealVar(lname.str().c_str(),"CR/SR ratio stat nuisance parameter",0);
	  ratio = 0;
	  ratiostat = 0;
	  std::cout << " ----- Check process names for CR/SR ratio: " << lProcs[vproc[iT][iR]] << " with " ;
	  if (iR<3) std::cout << lProcs[vproc[iT][iR]];
	  else std::cout << lProcs[vproc[iT][0]];
	  std::cout << std::endl;
	  if (iR<3) {
	    //WCR / WSR
	    ratio = histos[iR][vproc[iT][iR]][0]->GetBinContent(iB) / histos[0][vproc[iT][iR]][0]->GetBinContent(iB);
	    ratiostat = 1+sqrt(pow(histos[iR][vproc[iT][iR]][0]->GetBinError(iB)/histos[iR][vproc[iT][iR]][0]->GetBinContent(iB),2)+pow(histos[0][vproc[iT][iR]][0]->GetBinError(iB)/histos[0][vproc[iT][iR]][0]->GetBinContent(iB),2));
	  }
	  else {
	    //ZCR / ZSR
	    ratio = histos[iR][vproc[iT][iR]][0]->GetBinContent(iB) / histos[0][vproc[iT][0]][0]->GetBinContent(iB);
	    ratiostat = 1+sqrt(pow(histos[iR][vproc[iT][iR]][0]->GetBinError(iB)/histos[iR][vproc[iT][iR]][0]->GetBinContent(iB),2)+pow(histos[0][vproc[iT][0]][0]->GetBinError(iB)/histos[0][vproc[iT][0]][0]->GetBinContent(iB),2));
	  }
	  lFormula.str("");
	  lFormula << ratio;
	  lFormula << "*TMath::Power(";
	  if (iR<3) lFormula << jesWWSyst;
	  else lFormula << jesZZSyst;
	  lFormula << ",@0)*TMath::Power(" << ratiostat << ",@1)";

	  RooArgList nuisances;
	  nuisances.add(iR<3? *(jesWW[iT]) : *(jesZZ[iT]));
	  nuisances.add(*(TFstat[iT][iB][iR-1]));
	  for (unsigned iN(0); iN < nN; ++iN){
	    unsigned iSyst = 2+iN;
	    double ratiovar[2];
	    double ratiosyst[2];
	    //get up and down variations
	    for (unsigned iV(0); iV<2; ++iV){
	      unsigned iS = 2*iN+1+iV;
	      ratiovar[iV] = iR<3 ? histos[iR][vproc[iT][iR]][iS]->GetBinContent(iB) / histos[0][vproc[iT][iR]][iS]->GetBinContent(iB):
		histos[iR][vproc[iT][iR]][iS]->GetBinContent(iB) / histos[0][vproc[iT][0]][iS]->GetBinContent(iB);
	      ratiosyst[iV] = 1+(ratiovar[iV]-ratio)/ratio;
	      if (ratiosyst[iV] < 0){
		std::cout << " -- ERROR in systematics variations! For process " << lProcs[vproc[iT][iR]] << " region " << lRegions[iR] << " syst " << lSysts[iS] << " bin " << iB << ": ratio = " << ratio << " ratiovar = " << ratiovar[iV] << std::endl;
		return 1;
	      }
	      std::cout << std::setprecision(10) << " ------ bin " << iB << " type " << lType[iT] << " region " << lRegions[iR] << " Check syst " << lSysts[iS] << " " << iV << " " << ratiosyst[iV] << std::endl;
	    }
	    if (ratiovar[0]<ratio){
	      std::cout << " -- INFO: up variations is actually giving smaller ratio...." << lProcs[vproc[iT][iR]] << " region " << lRegions[iR] << " syst " << lNuis[iN] << " bin " << iB << ": ratio = " << ratio << " ratiovar = " << ratiovar[0] << std::endl;
	    }
	    if (ratiovar[1]>ratio){
	      std::cout << " -- INFO: down variations is actually giving larger ratio...." << lProcs[vproc[iT][iR]] << " region " << lRegions[iR] << " syst " << lNuis[iN] << " bin " << iB << ": ratio = " << ratio << " ratiovar = " << ratiovar[1] << std::endl;
	    }
	    //take sign of @i to decide up / down...
	    lFormula << "*( (@" << iSyst << ">=0)*TMath::Power(" << ratiosyst[0] << ",@" << iSyst << ")+(@" << iSyst << "<0)*TMath::Power(" << 1./ratiosyst[1] << ",@" << iSyst << "))";
	    nuisances.add(*(TFsysts[iN]));
	  }


	  lname.str("");
	  lname << lType[iT] << "TF_" << lRegions[iR] << "_bin" << iB;
	  RooFormulaVar TF(lname.str().c_str(),"Transfer factor CR/SR",lFormula.str().c_str(),nuisances);
	  wspace.import(TF,RooFit::RecycleConflictNodes());
	  lname.str("");
	  lname << lType[iT] << "V_" << lRegions[iR] << "_bin" << iB;
	  RooFormulaVar CRbin(lname.str().c_str(),(lType[iT]+" V+jets yield in control regions, per bin").c_str(),"@0*@1",iR<3?RooArgList(TF,WZbin):iT==0?RooArgList(TF,binParZ):RooArgList(TF,*(EWKQCDbin)));
	  wspace.import(CRbin,RooFit::RecycleConflictNodes());
	  
	  if (iR==3) std::cout << " ---- CR " << lRegions[iR] << " V yield = " << histos[iR][vproc[iT][iR]][0]->GetBinContent(iB) << " data yield " << histos[iR][0][0]->GetBinContent(iB) << std::endl;
	  
	  
	}//loop on regions
      }//loop on types


    }//loop on bins

    std::cout << " - Adding to lists: " << std::endl;
    RooArgList Z_SR_bins[nT];
    RooArgList W_SR_bins[nT];
    RooArgList V_CR_bins[nT][nR-1];
    for (unsigned iT(0); iT<nT; ++iT){     
      for (unsigned iB(0); iB<nB; ++iB){
	std::ostringstream lname;
	if (iT==0){
	  lname << lType[iT] << "Z_SR_bin" << iB+1;
	} else {
	  lname << "EWKQCD_SR_bin" << iB+1;
	}
	if ( (iT==0 && !wspace.var(lname.str().c_str())) ||
	     (iT==1 && !wspace.function(lname.str().c_str()))
	     ) {
	  std::cout << "Error for " << lType[iT] << " Z bin " << iB << " " << lname.str() << std::endl;
	  return 1;
	}
	if (iT==0) Z_SR_bins[iT].add(*wspace.var(lname.str().c_str()));
	else Z_SR_bins[iT].add(*wspace.function(lname.str().c_str()));
	lname.str("");
	lname << lType[iT] << "WZ_SR_bin" << iB+1;
	if (!wspace.function(lname.str().c_str())) {
	  std::cout << "Error for " << lType[iT] << " W bin " << iB << " " << lname.str() << std::endl;
	  return 1;
	}
	W_SR_bins[iT].add(*wspace.function(lname.str().c_str()));
	for (unsigned iR(1); iR<nR; ++iR){
	  lname.str("");
	  lname << lType[iT] << "V_" << lRegions[iR] << "_bin" << iB+1;
	  if (!wspace.function(lname.str().c_str())) {
	    std::cout << "Error for  " << lType[iT] << " Z bin " << iB << " region " << lRegions[iR] << " " << lname.str() << std::endl;
	    return 1;
	  }
	  V_CR_bins[iT][iR-1].add(*wspace.function(lname.str().c_str()));
	}
      }//loop on bins
    }//loop on types
    
    std::cout << " - Creating the SR parametric hists" << std::endl;
    
    // Create a RooParametericHist which contains those yields, last argument is just for the binning,
    // can use the data TH1 for that

    for (unsigned iT(0); iT<nT; ++iT){     
      std::cout << " -- Processing type " << lType[iT] << std::endl;

      RooParametricHist p_Z((lType[iT]+"Z_SR").c_str(), (lType[iT]+"Z+jets PDF in signal region").c_str(),lVarFit,Z_SR_bins[iT],dummyHist);
      // Always include a _norm term which should be the sum of the yields (thats how combine likes to play with pdfs)
      RooAddition p_Z_norm((lType[iT]+"Z_SR_norm").c_str(),("Total Number of events from "+lType[iT]+" Z+jets in signal region").c_str(),Z_SR_bins[iT]);
      
      RooParametricHist p_W((lType[iT]+"W_SR").c_str(), (lType[iT]+"W+jets PDF in signal region").c_str(),lVarFit,W_SR_bins[iT],dummyHist);
      RooAddition p_W_norm((lType[iT]+"W_SR_norm").c_str(),("Total Number of events from "+lType[iT]+" W+jets in signal region").c_str(),W_SR_bins[iT]);
    
      std::cout << " -- Importing the parametric hists" << std::endl;
    
      // import the pdfs
      wspace.import(p_Z);
      wspace.import(p_Z_norm,RooFit::RecycleConflictNodes());
      wspace.import(p_W);
      wspace.import(p_W_norm,RooFit::RecycleConflictNodes());
      
      std::cout << " -- Creating and importing the CR parametric hists" << std::endl;
    
      for (unsigned iR(1); iR<nR; ++iR){
	std::ostringstream lname;
	lname << lType[iT] << "V_" << lRegions[iR];
	RooParametricHist p_CRV(lname.str().c_str(), "Background PDF in control region",lVarFit,V_CR_bins[iT][iR-1],dummyHist);
	lname << "_norm";
	RooAddition p_CRV_norm(lname.str().c_str(),("Total Number of events from "+lType[iT]+" V+jets background in control region").c_str(),V_CR_bins[iT][iR-1]);
	wspace.import(p_CRV);
	wspace.import(p_CRV_norm,RooFit::RecycleConflictNodes());
      }
    }//loop on types
    
    std::cout << " - Printing workspace." << std::endl;
    wspace.Print();
    
    fOut->cd();
    wspace.Write();
    
    std::cout << " - Workspace " << wspace.GetName() << " written." << std::endl;
    
    // Clean up
    fOut->Close();
    fOut->Delete();
    
    return 0;
    
}
