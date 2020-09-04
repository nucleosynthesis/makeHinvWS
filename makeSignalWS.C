#include <iostream>
#include <iomanip>
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
//#include "../../interface/RooParametricHist.h"
#include "../HiggsAnalysis/CombinedLimit/interface/RooParametricHist.h"
//#include "/afs/cern.ch/user/v/vmilosev/CMSSW_10_2_13/src/HiggsAnalysis/CombinedLimit/interface/RooParametricHist.h"
#include "RooAddition.h"

enum PROCESS{
  VBFH = 0,
  ggH = 1,
};
void makeSignalWS(std::string year="2017", std::string cat="MTR"){


    const bool is2017 = year=="2017";
    std::string lChannel = "VBF";
    std::string lCategory = cat+"_";
    std::string lYear = year+"_";
    std::string lOutFileName = "signal_ws_"+lCategory+lYear+lChannel+".root";

    //define the variable that is fitted
    RooRealVar lVarFit(("mjj_"+cat+"_"+year).c_str(),"M_{jj} (GeV)",200,5000);
    std::string lVarLabel = "Mjj";
    
    TFile *fOut = new TFile(lOutFileName.c_str(),"RECREATE");
    RooWorkspace wspace("wspace_signal","wspace_signal");
    RooArgList vars(lVarFit);
    
    std::string lRegions = "SR";
    std::string lInFileName = "out_VBF_ana_"+lRegions+"_"+year+"_v"+cat+"_"+year+"_200109/VBF_shapes.root";

    TFile *finput = TFile::Open(lInFileName.c_str());
    //finput->cd(lRegions.c_str());
    
    TFile *finputJES = TFile::Open("../vbf_shape_jes_uncs.root");

    const unsigned nP = 2;
    std::string lProcs[nP] = {"VBFHtoInv","GluGluHtoInv"};


    const unsigned nN = 7;
    std::string lSysts[nN] = {"bjet_veto","pileup","tau_veto",
			     "eventVetoVEleIdIso","eventVetoLMuId","eventVetoLMuIso",
                             "eventVetoVEleReco",			      
    };

    const unsigned nJ = 11;
    std::string lJes[nJ] = {
		  "jesAbsolute"
		 , Form("jesAbsolute_%s",year.c_str())
		 , "jesBBEC1"
		 , Form("jesBBEC1_%s",year.c_str())
		 , "jesEC2"
		 , Form("jesEC2_%s",year.c_str())
		 , "jesFlavorQCD"
		 , "jesHF"
		 , Form("jesHF_%s",year.c_str())
		 , "jesRelativeBal"
		 , Form("jesRelativeSample_%s",year.c_str())
    };


   for (unsigned iP(0); iP<nP; ++iP){
    
       TH1F* Thist = (TH1F*)finput->Get(Form("%s/%s",lRegions.c_str(),lProcs[iP].c_str()));
       RooDataHist *hist = new RooDataHist((lProcs[iP]+"_hist_"+lRegions).c_str(),"Signal proces",vars,Thist);
       wspace.import(*hist);

       for (unsigned iS(0); iS < nN; ++iS){
         
	 TH1F* ThistSU = (TH1F*)finput->Get(Form("%s/%s_%sUp",lRegions.c_str(),lProcs[iP].c_str(),lSysts[iS].c_str()));
         RooDataHist *histSU = new RooDataHist((lProcs[iP]+"_hist_"+lRegions+"_"+lSysts[iS]+"Up").c_str(),"Signal proces",vars,ThistSU);
         wspace.import(*histSU);	    

         TH1F* ThistSD = (TH1F*)finput->Get(Form("%s/%s_%sDown",lRegions.c_str(),lProcs[iP].c_str(),lSysts[iS].c_str()));
         RooDataHist *histSD = new RooDataHist((lProcs[iP]+"_hist_"+lRegions+"_"+lSysts[iS]+"Down").c_str(),"Signal proces",vars,ThistSD);
         wspace.import(*histSD);	     

       }

       // now create the variations due to the different sources of JES
       std::cout << finputJES->GetName() <<std::endl;
       for (unsigned iJ(0); iJ < nJ; ++iJ){

	
	  TH1F *hSUp   = (TH1F*)finputJES->Get(Form("VBF%s_%sUp_smoothed",year.c_str(),lJes[iJ].c_str()));
	  TH1F *hSDown = (TH1F*)finputJES->Get(Form("VBF%s_%sDown_smoothed",year.c_str(),lJes[iJ].c_str()));

	  TH1F *hSUpnew = (TH1F*)Thist->Clone(); hSUpnew->SetName(Form("%s/%s_%sUp",lRegions.c_str(),lProcs[iP].c_str(),lJes[iJ].c_str()));
	  TH1F *hSDownnew = (TH1F*)Thist->Clone(); hSDownnew->SetName(Form("%s/%s_%sDown",lRegions.c_str(),lProcs[iP].c_str(),lJes[iJ].c_str()));

	  for (int b=1; b <= Thist->GetNbinsX() ; b++){
	  	double xv = hSUpnew->GetBinCenter(b);
		double yv = hSUpnew->GetBinContent(b);
		if ( xv > hSUp->GetBinLowEdge(hSUp->GetNbinsX()) ) xv = hSUp->GetBinLowEdge(hSUp->GetNbinsX())+0.5*hSUp->GetBinWidth(hSUp->GetNbinsX());

		hSUpnew->SetBinContent(b,yv*hSUp->GetBinContent(hSUp->FindBin(xv)));
		hSDownnew->SetBinContent(b,yv*hSDown->GetBinContent(hSDown->FindBin(xv)));
	  }
         
	 RooDataHist *histSU = new RooDataHist((lProcs[iP]+"_hist_"+lRegions+"_CMS_scale_j_"+lJes[iJ]+"Up").c_str(),"Signal proces",vars,hSUpnew);
         RooDataHist *histSD = new RooDataHist((lProcs[iP]+"_hist_"+lRegions+"_CMS_scale_j_"+lJes[iJ]+"Down").c_str(),"Signal proces",vars,hSDownnew);
	 wspace.import(*histSU);
	 wspace.import(*histSD);

	}

   }


   fOut->WriteTObject(&wspace);

};

