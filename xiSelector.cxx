#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TSystem.h"
#include "TStopwatch.h"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RResultPtr.hxx"

#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TString.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TFile.h" 
#include "HarryPlotter.hh"

std::vector<float> ptbins = HarryPlotter::Getptbins(); 
std::vector<float> layerPos = HarryPlotter::GetposLayers();

/* 
   workflow: 
   1) take the xis for kind of granted -> QA cuts
   -> additional hits + pT cut 
   -> transverse Radius (difference between truth and reco, causality) + Daughter DCA (?)
   -> invariant Mass cut (?)
	  
   2) Ultimately: Select the Xi_c from Xi_cc 
   -> QA Cuts: transverse radius, dca among daughters, dca to PV ... invariant mass
   -> Pion DCA distribution... 
	  
   3) Select the Xi_cc 
   -> Xi_c from Xi_cc Cuts: DCA to PV 
   What this needs to be able: 
   1) process signal + background files 

   Some random documentation: 

   Things to think about: 
   1) Do we compare Topo & Strangeness Tracking? 
   2) Switch between signal and background
*/ 


int main(int argc, char **argv) {
  TStopwatch timer;
  timer.Start();
  
  const char* fileName = argv[1]; 
  const char* outAddon = (argv[2])?argv[2]:""; 
  
  int xiDec = argv[3]?atoi(argv[3]):0; 
  bool ExclusiveSignal = (xiDec == 0)?false:true;
  
  int pTbin = argv[4]?atoi(argv[4]):-1; 
  if (pTbin >= 0 && pTbin > ptbins.size()-1) { 
    std::cout << "Crashing cause the requested pT bin is out of range ( requested = " << pTbin << " , max. available = " << ptbins.size()-1 << " )\n";
    return -999; 
  }
  
  HarryPlotter::StyleBox(); 
 
  ROOT::EnableImplicitMT(); // Tell ROOT you want to go parallel          
  TString filePath = TString::Format("%s", fileName); 
  
  auto h_cand_counter = new TH1D("df_xi_c_candCounter", "candCounter", 1, 0, 1); 

  auto h_gen_xi_pt_eta_counter = new TH2D("ptetaXiGen", "candCounter", 100, 0, 10, 30, -1.5, 1.5); 
  auto h_gen_xi_pt_y_counter = new TH2D("ptyXiGen", "candCounter", 100, 0, 10, 30, -1.5, 1.5); 
   
  TChain input("fTreeCandidates"); 
  int inputFiles = 0; 
  int inputFailures = 0; 
  if (filePath.Contains(".root")) { 
    TFile *inFile = TFile::Open(filePath);
    TH1D* evtCounter = (TH1D*)inFile->Get("hEventCounter"); 
    TH2D* ptXipteta = (TH2D*)inFile->Get("hPtEtaGeneratedXi"); 
    TH2D* ptXipty = (TH2D*)inFile->Get("hPtYGeneratedXi"); 

    if (!evtCounter||!ptXipteta||!ptXipty) { 
      inputFailures++; 
      std::cout << "Zis is vehry bad, ze generation histograms are missing, Guenther! No histogram, no chain! \n"; 
    } else { 
      h_cand_counter->Add(evtCounter); 
	
      h_gen_xi_pt_eta_counter->Add(ptXipteta); 
      h_gen_xi_pt_y_counter->Add(ptXipty); 
	
      input.Add(filePath); 
      inputFiles++;
    }
    
  } else { 
    std::cout << "Input seems not to be a file, trying to build a chain from sudirs...\n"; 
    // if file does not exist, loop over subdirs to find TFiles 
    TSystemDirectory dir("ZeDirectory", filePath);
    auto files = dir.GetListOfFiles();
    for (auto fileObj : *files)  {
      auto file = (TSystemFile*) fileObj;
      TString inSubDirFile = TString::Format("%s/%s/treeoutput.root", filePath.Data(), file->GetName()).Data(); 
      if (!gSystem->AccessPathName(inSubDirFile)) { 
	TFile *inFile = TFile::Open(inSubDirFile);
	if (!inFile) { 
	  inputFailures++;
	  continue;
	}
	if (inFile->IsZombie()) { 
	  inFile->Close();
	  inputFailures++;
	  continue; 
	} 
	TH1D* evtCounter = (TH1D*)inFile->Get("hEventCounter"); 
	TH2D* ptXipteta = (TH2D*)inFile->Get("hPtEtaGeneratedXi"); 
	TH2D* ptXipty = (TH2D*)inFile->Get("hPtYGeneratedXi"); 

	if (!evtCounter||!ptXipteta||!ptXipty) { 
	  inputFailures++;
	  continue; 
	}
	h_cand_counter->Add(evtCounter); 

	h_gen_xi_pt_eta_counter->Add(ptXipteta); 
	h_gen_xi_pt_y_counter->Add(ptXipty); 
	
	input.Add(inSubDirFile);
	inputFiles++;
	inFile->Close();
      }
    }
  }
  
  std::cout << "Added " << inputFiles << " files to the chain, failed in " << inputFailures << " cases \n"; 

  //Lmb cuts
  float radDiffMax = 1.; //cm 

  float lmbMass = 1.116; 
  float invMassDiffLmb = 0.012; //8 MeV/c2 mass window 
  auto invMassLmbCut = [&invMassDiffLmb, &lmbMass](float invMass) { return (TMath::Abs(invMass-lmbMass) < invMassDiffLmb); }; 
  auto radCut = [&radDiffMax](float radDiff) { return (radDiff < radDiffMax);};
  auto decLengthLmb = [&lmbMass](float len, float mom){ return TMath::Abs(mom)>1e-4?len*lmbMass/mom:-999; }; 
  
  //Xi cuts
  float xiMass = 1.322; 
  float invMassDiffXi = 0.012; 
  int addedHitsMin = !ExclusiveSignal?0:1; 
  
  auto invMassXiCut = [&invMassDiffXi, &xiMass](float invMass) { return (TMath::Abs(invMass-xiMass) < invMassDiffXi); }; 
  auto hitsCut = [&addedHitsMin](int AddedHits) { return (AddedHits >= addedHitsMin); }; 
  
  auto decLengthXi = [&xiMass](float len, float mom){ return TMath::Abs(mom)>1e-4?len*xiMass/mom:-999; }; 

  float xipTmin = pTbin >=0?ptbins[pTbin]:0; 
  float xipTmax = pTbin >=0?ptbins[pTbin+1]:20.0; 
  
  std::cout << "Xi pT range set to " << xipTmin << " to " << xipTmax << std::endl; 
  
  auto pTCut = [&xipTmin, &xipTmax](float pT) { return (xipTmin < pT)&&(pT < xipTmax); }; 

  ROOT::RDataFrame df(input);
  
  if (ExclusiveSignal) {
    std::cout << "Looking Exclusively for Signal!\n"; 
  } else { 
    std::cout << "Rejecting Signal!\n"; 
  }


  if (!ExclusiveSignal) { 
    std::cout << "Rejecting Xis, make sure you know what you doing!\n"; 
  }
  
  auto df_ForceXi = (!ExclusiveSignal)?df.Filter("!fTrueXi","FakeXis"):df.Filter("fTrueXi","TrueXis");

  auto df_in = df_ForceXi
    .Filter("TMath::Abs(fXiCCEta)<0.5")
    .Filter(pTCut, {"fXiPtStraTrack"}, "pTXi")
    ;
  
  auto h_df_in_im_xi_mass_stra = df_in.Histo1D({"h_df_in_im_xi_mass_stra", "xi_cc inv mass", 700, 2.6, 4.6}, "fXiMass"); 
  auto in_counter = df_in.Count(); 
  
  auto df_in_qa = df_in; 
  
  //towards the Lambda for the Xi 
  auto h_df_in_qa_pos_pt = df_in_qa.Histo1D({"df_in_qa_pos_pt", "pos pt", 100, 0, 10}, "fPositivePt"); 
  auto h_df_in_qa_pos_dca_xy = df_in_qa.Histo1D({"df_in_qa_pos_dca_xy", "pos dca xy pv", 1000, -1e4, 1e4}, "fPositiveDCAxy"); 
  auto h_df_in_qa_pos_dca_z = df_in_qa.Histo1D({"df_in_qa_pos_dca_z", "pos dca z pv", 1000, -1e4, 1e4}, "fPositiveDCAz"); 
  auto h_df_in_qa_pos_dca_xy_wide = df_in_qa.Histo1D({"df_in_qa_pos_dca_xy_wide", "pos dca xy pv", 1000, -1e5, 1e5}, "fPositiveDCAxy"); 
  auto h_df_in_qa_pos_dca_z_wide = df_in_qa.Histo1D({"df_in_qa_pos_dca_z_wide", "pos dca z pv", 1000, -1e5, 1e5}, "fPositiveDCAz"); 
  auto h_df_in_qa_pos_hits = df_in_qa.Histo1D({"df_in_qa_pos_hits", "pos hits", 15, 0, 15}, "fPositiveClusters"); 
  auto h_df_in_qa_pos_chisq = df_in_qa.Histo1D({"df_in_qa_pos_chisq", "pos chisq", 100, 0, 100}, "fPositiveChisquare"); 
  auto h_df_in_qa_pos_chisqhits = df_in_qa.Define("fPositiveChisquareOverHits", "fPositiveChisquare/fPositiveClusters").Histo1D({"df_in_qa_pos_chisqhits", "pos chisq over hits", 100, 0, 50}, "fPositiveChisquareOverHits"); 
  
  auto h_df_in_qa_pos_dca_xy_vs_hits = df_in_qa.Histo2D<float,int>({"df_in_qa_pos_dca_xy_vs_hits", "pos dca xy pv vs hits", 1000, -10000, 10000, 15, -0.5, 14.5}, "fPositiveDCAxy", "fPositiveClusters"); 
  auto h_df_in_qa_pos_dca_xy_vs_lmb_dl_pv = df_in_qa.Define("fLmbInvDecayLengthToPV", decLengthLmb, {"fV0DecayLength", "fV0TotalMomentum"}).Histo2D<float,float>({"df_in_qa_pos_dca_xy_vs_lmb_dl_pv", "pos dca xy pv vs lmb decay length", 1000, -10000, 10000, 100, 0, 40}, "fPositiveDCAxy", "fLmbInvDecayLengthToPV"); 
  auto h_df_in_qa_pos_dca_xy_vs_pos_pt = df_in_qa.Histo2D<float,float>({"h_df_in_qa_pos_dca_xy_vs_pos_pt", "pos dca xy pv vs pos pt", 1000, -10000, 10000, 100, 0, 10}, "fPositiveDCAxy", "fPositivePt"); 
    
  auto h_df_in_qa_neg_pt = df_in_qa.Histo1D({"df_in_qa_neg_pt", "neg pt", 100, 0, 10}, "fNegativePt"); 
  auto h_df_in_qa_neg_dca_xy = df_in_qa.Histo1D({"df_in_qa_neg_dca_xy", "neg dca xy pv", 1000, -1e4, 1e4}, "fNegativeDCAxy"); 
  auto h_df_in_qa_neg_dca_z = df_in_qa.Histo1D({"df_in_qa_neg_dca_z", "neg dca z pv", 1000, -1e4, 1e4}, "fNegativeDCAz"); 
  auto h_df_in_qa_neg_dca_xy_wide = df_in_qa.Histo1D({"df_in_qa_neg_dca_xy_wide", "neg dca xy pv", 1000, -1e5, 1e5}, "fNegativeDCAxy"); 
  auto h_df_in_qa_neg_dca_z_wide = df_in_qa.Histo1D({"df_in_qa_neg_dca_z_wide", "neg dca z pv", 1000, -1e5, 1e5}, "fNegativeDCAz"); 
  auto h_df_in_qa_neg_hits = df_in_qa.Histo1D({"df_in_qa_neg_hits", "neg hits", 15, 0, 15}, "fNegativeClusters"); 
  auto h_df_in_qa_neg_chisq = df_in_qa.Histo1D({"df_in_qa_neg_chisq", "neg chisq", 100, 0, 100}, "fNegativeChisquare"); 
  auto h_df_in_qa_neg_chisqhits = df_in_qa.Define("fNegativeChisquareOverHits", "fNegativeChisquare/fNegativeClusters").Histo1D({"df_in_qa_neg_chisqhits", "neg chisq over hits", 100, 0, 50}, "fNegativeChisquareOverHits"); 
  
  auto h_df_in_qa_neg_dca_xy_vs_hits = df_in_qa.Histo2D<float,int>({"df_in_qa_neg_dca_xy_vs_hits", "neg dca xy pv vs hits", 1000, -10000, 10000, 15, -0.5, 14.5}, "fNegativeDCAxy", "fNegativeClusters"); 
  auto h_df_in_qa_neg_dca_xy_vs_lmb_dl_pv = df_in_qa.Define("fLmbInvDecayLengthToPV", decLengthLmb, {"fV0DecayLength", "fV0TotalMomentum"}).Histo2D<float,float>({"df_in_qa_neg_dca_xy_vs_lmb_dl_pv", "neg dca xy pv vs lmb decay length", 1000, -10000, 10000, 100, 0, 40}, "fNegativeDCAxy", "fLmbInvDecayLengthToPV"); 
  auto h_df_in_qa_neg_dca_xy_vs_neg_pt = df_in_qa.Histo2D<float,float>({"h_df_in_qa_neg_dca_xy_vs_neg_pt", "neg dca xy pv vs neg pt", 1000, -10000, 10000, 100, 0, 10}, "fNegativeDCAxy", "fNegativePt"); 

  auto h_df_in_qa_neg_hits_vs_pos_hits = df_in_qa.Histo2D<int,int>({"h_df_in_qa_neg_hits_vs_pos_hits", "neg hits vs pos hits", 15, -0.5, 14.5, 15, -0.5, 14.5}, "fNegativeClusters", "fPositiveClusters"); 
  auto h_df_in_qa_neg_dca_xy_vs_pos_dca_xy = df_in_qa.Histo2D<float,float>({"df_in_qa_neg_dca_xy_vs_pos_dca_xy", "neg dca xy pv vs pos dca xy pv", 1000, -10000, 10000, 1000, -10000, 10000}, "fNegativeDCAxy", "fPositiveDCAxy"); 

  auto h_df_in_qa_lmb_dca_xy = df_in_qa.Histo1D({"df_in_qa_lmb_dca_xy", "lmb dca xy pv", 1000, -1e4, 1e4}, "fV0DCAxyToPV");
  auto h_df_in_qa_lmb_dca_z = df_in_qa.Histo1D({"df_in_qa_lmb_dca_z", "lmb dca z pv", 1000, -1e4, 1e4}, "fV0DCAzToPV");
  auto h_df_in_qa_lmb_dca_xy_wide = df_in_qa.Histo1D({"df_in_qa_lmb_dca_xy_wide", "lmb dca xy pv", 1000, -1e5, 1e5}, "fV0DCAxyToPV");
  auto h_df_in_qa_lmb_dca_z_wide = df_in_qa.Histo1D({"df_in_qa_lmb_dca_z_wide", "lmb dca z pv", 1000, -1e5, 1e5}, "fV0DCAzToPV");
  
  auto h_df_in_qa_lmb_totp = df_in_qa.Histo1D({"df_in_qa_lmb_totp", "lmb tot p", 100, 0, 10}, "fV0TotalMomentum");   
  auto h_df_in_qa_lmb_ddist_pv = df_in_qa.Define("fLmbInvDecayLengthToPV", decLengthLmb, {"fV0DecayLength", "fV0TotalMomentum"}).Histo1D({"df_in_qa_lmb_ddist_pv", "lmb ddist to pv", 1000, 0, 40}, "fLmbInvDecayLengthToPV");
  auto h_df_in_qa_lmb_ddca = df_in_qa.Histo1D({"df_in_qa_lmb_ddca", "lmb prong dca", 500, 0, 100}, "fXiV0DauDCA"); 
  auto h_df_in_qa_lmb_ddca_wide = df_in_qa.Histo1D({"df_in_qa_lmb_ddca_wide", "lmb prong dca", 500, 0, 5000}, "fXiV0DauDCA"); 
  auto h_df_in_qa_lmb_trad = df_in_qa.Histo1D({"df_in_qa_lmb_trad", "lmb trad", 600, 0, 30}, "fV0DecayRadius"); 
  auto h_df_in_qa_lmb_trad_mc = df_in_qa.Histo1D({"df_in_qa_lmb_trad_mc", "lmb trad", 600, 0, 30}, "fV0DecayRadiusMC"); 
  auto h_df_in_qa_lmb_trad_diff = df_in_qa.Define("XiV0DecayRadDiff", "TMath::Abs(fV0DecayRadiusMC-fV0DecayRadius)").Histo1D({"df_in_qa_lmb_trad_diff", "lmb trad", 500, 0, 2}, "XiV0DecayRadDiff"); 
  
  auto h_df_in_qa_lmb_mass = df_in_qa.Histo1D({"df_in_qa_lmb_mass", "lmb inv mass", 750, 1., 1.8}, "fLambdaMass"); 
  auto h_df_in_qa_xi_mass = df_in_qa.Histo1D({"df_in_qa_xi_mass", "xi inv mass", 750, 1.2, 2}, "fXiMass"); 
  
  //Define future variables and select Lambdas 

  auto df_lmb_im = df_in
    //.Define("XiV0DecayRadDiff", "TMath::Abs(fV0DecayRadiusMC-fV0DecayRadius)")
    .Define("XiDecayRadDiff", "TMath::Abs(fXiDecayRadiusMC-fXiDecayRadius)")
    .Define("XiLmbDecayRadDiff", "fV0DecayRadius-fXiDecayRadius")
    .Define("XiXicDecayRadDiffTopo", "fXiDecayRadius-fXicDecayRadiusTopo")
    .Define("XiXicDecayRadDiffStra", "fXiDecayRadius-fXicDecayRadiusStraTrack")
    .Define("fLmbInvDecayLengthToPV", decLengthLmb, {"fV0DecayLength", "fV0TotalMomentum"})
    .Define("fXiInvDecayLengthToPV", decLengthXi, {"fXiDecayLength", "fXiTotalMomentum"})
    //.Filter(radCut, {"XiV0DecayRadDiff"})
    // .Filter("TMath::Abs(fV0DCAxyToPV) < 5000", "fV0DCAxyToPV")
    // .Filter("TMath::Abs(fV0DCAzToPV) < 7000", "fV0DCAzToPV")
    // .Filter("fXiV0DauDCA < 2000","fXiV0DauDCA")
    // .Filter("fV0DecayRadius > 0.5","fV0DecayRadius")
    // .Filter("fLmbInvDecayLengthToPV > 0.04","fLmbInvDecayLengthToPV")
    // .Filter("TMath::Abs(fPositiveDCAxy) > 50","fPositiveDCAxy")
    // .Filter("TMath::Abs(fPositiveDCAz) > 40","fPositiveDCAz")
    // .Filter("TMath::Abs(fNegativeDCAxy) > 100","fNegativeDCAxy")
    // .Filter("TMath::Abs(fNegativeDCAz) > 50","fNegativeDCAz")
    ;

  auto h_df_lmb_im_lmb_mass = df_lmb_im.Histo1D({"df_lmb_im_lmb_mass", "lmb inv mass", 750, 1., 1.8}, "fLambdaMass"); 

  auto df_lmb =  df_lmb_im
    .Filter(invMassLmbCut, {"fLambdaMass"}, "fLambdaMass")
    ;

  //Towards the Xi for Xic->Xi+2pi 
  //Fill some Histograms
  
  auto df_lmb_qa = df_lmb;

  auto h_df_lmb_qa_bach_pt = df_lmb_qa.Histo1D({"df_lmb_qa_bach_pt", "bach pt", 100, 0, 10}, "fBachelorPt"); 
  auto h_df_lmb_qa_bach_dca_xy = df_lmb_qa.Histo1D({"df_lmb_qa_bach_dca_xy", "bach dca xy pv", 1000, -1e4, 1e4}, "fBachelorDCAxy"); 
  auto h_df_lmb_qa_bach_dca_z = df_lmb_qa.Histo1D({"df_lmb_qa_bach_dca_z", "bach dca z pv", 1000, -1e4, 1e4}, "fBachelorDCAz"); 
  auto h_df_lmb_qa_bach_dca_xy_wide = df_in_qa.Histo1D({"df_in_qa_bach_dca_xy_wide", "bach dca xy pv", 1000, -1e5, 1e5}, "fBachelorDCAxy"); 
  auto h_df_lmb_qa_bach_dca_z_wide = df_in_qa.Histo1D({"df_in_qa_bach_dca_z_wide", "bach dca z pv", 1000, -1e5, 1e5}, "fBachelorDCAz"); 
  auto h_df_lmb_qa_bach_hits = df_lmb_qa.Histo1D({"df_lmb_qa_bach_hits", "bach hits", 15, 0, 15}, "fBachelorClusters"); 
  auto h_df_lmb_qa_bach_chisq = df_lmb_qa.Histo1D({"df_lmb_qa_bach_chisq", "bach chisq", 100, 0, 100}, "fBachelorChisquare"); 
  auto h_df_lmb_qa_bach_chisqhits = df_lmb_qa.Define("fBachelorChisquareOverHits", "fBachelorChisquare/fBachelorClusters").Histo1D({"df_lmb_qa_bach_chisqhits", "bach chisq over hits", 100, 0, 50}, "fBachelorChisquareOverHits"); 

  auto h_df_in_qa_bach_dca_xy_vs_hits = df_in_qa.Histo2D<float,int>({"df_in_qa_bach_dca_xy_vs_hits", "bach dca xy pv vs hits", 1000, -10000, 10000, 15, -0.5, 14.5}, "fBachelorDCAxy", "fBachelorClusters"); 
  auto h_df_in_qa_bach_dca_xy_vs_lmb_dl_pv = df_in_qa.Define("fLmbInvDecayLengthToPV", decLengthLmb, {"fV0DecayLength", "fV0TotalMomentum"}).Histo2D<float,float>({"df_in_qa_bach_dca_xy_vs_lmb_dl_pv", "bach dca xy pv vs lmb decay length", 1000, -10000, 10000, 100, 0, 40}, "fBachelorDCAxy", "fLmbInvDecayLengthToPV"); 
  auto h_df_in_qa_bach_dca_xy_vs_bach_pt = df_in_qa.Histo2D<float,float>({"h_df_in_qa_bach_dca_xy_vs_bach_pt", "bach dca xy pv vs bach pt", 1000, -10000, 10000, 100, 0, 10}, "fBachelorDCAxy", "fBachelorPt"); 

  auto h_df_lmb_qa_xi_pt = df_lmb_qa.Histo1D({"df_lmb_qa_xi_pt", "xi pt", 100, 0, 10}, "fXiPtStraTrack"); 
  auto h_df_lmb_qa_xi_ddist_pv = df_lmb_qa.Histo1D({"df_lmb_qa_xi_ddist_pv", "xi ddist to pv", 1000, 0, 20}, "fXiInvDecayLengthToPV");
  auto h_df_lmb_qa_xi_ddca = df_lmb_qa.Histo1D({"df_lmb_qa_xi_ddca", "xi prong dca", 500, 0, 2000}, "fXiCascDauDCA"); 
  auto h_df_lmb_qa_xi_trad = df_lmb_qa.Histo1D({"df_lmb_qa_xi_trad", "xi trad", 250, 0, 10}, "fXiDecayRadius"); 
  auto h_df_lmb_qa_xi_trad_mc = df_lmb_qa.Histo1D({"df_lmb_qa_xi_trad_mc", "xi trad", 250, 0, 10}, "fXiDecayRadiusMC"); 
  auto h_df_lmb_qa_xi_trad_diff = df_lmb_qa.Histo1D({"df_lmb_qa_xi_trad_diff", "xi trad", 250, 0, 2}, "XiDecayRadDiff"); 

  auto h_df_lmb_qa_trad_diff_lmb_xi = df_lmb_qa.Histo1D({"df_lmb_qa_trad_diff_lmb_xi", "lmb-xi trad", 500, -100, 150}, "XiLmbDecayRadDiff"); 

  auto h_df_lmb_qa_xi_mass = df_lmb_qa.Histo1D({"df_lmb_qa_xi_mass", "xi inv mass", 750, 1.2, 2}, "fXiMass"); 
  
  //Select Xis excluding hits to avoid cheating 
  auto df_xi_sel = df_lmb
    .Filter("fXiDecayRadius > 0.5","fXiDecayRadius")
    .Filter("XiLmbDecayRadDiff > 0","XiLmbDecayRadDiff")
    // .Filter("fXiCascDauDCA> 4 && fXiCascDauDCA < 1400")
    //.Filter("fXiDecayLength > 0.02","fXiDecayLength")
    //.Filter("TMath::Abs(fBachelorDCAxy) > 40","fBachelorDCAxy")
    //.Filter("TMath::Abs(fBachelorDCAz) > 40","fBachelorDCAz")
    ;

  auto h_df_xi_sel_xi_mass = df_xi_sel.Histo1D({"df_xi_sel_xi_mass", "xi inv mass", 750, 1.2, 2}, "fXiMass"); 

  //Select Xis including Hits
  auto df_xi_im = df_xi_sel.
    Filter(hitsCut, {"fXiHitsAdded"},"fXiHitsAdded")
    ; 
  
  auto h_df_xi_im_xi_mass = df_xi_im.Histo1D({"df_xi_im_xi_mass", "xi inv mass", 750, 1.2, 2}, "fXiMass"); 
  
  auto h_df_xi_im_xi_ddist_dv_topo = df_xi_im.Histo1D({"df_xi_im_xi_dist_dv_topo", "xi decay dist", 1500, 0, 30}, "fXiInvDecayLengthToDVTopo"); 
  auto h_df_xi_im_xi_ddist_dv_stra = df_xi_im.Histo1D({"df_xi_im_xi_dist_dv_stra", "xi decay dist", 1500, 0, 30}, "fXiInvDecayLengthToDVStra"); 

  auto df_xi = df_xi_im
    .Filter(invMassXiCut, {"fXiMass"},"fXiMass")
    ;
  auto h_df_xi_qa_xi_dca_xy_topo = df_xi_im.Histo1D({"df_xi_qa_xi_dca_xy_topo", "xi dca xy topo", 1000, -500, 500}, "fXiDCAxyToPVTopo");  
  auto h_df_xi_qa_xi_dca_z_topo = df_xi_im.Histo1D({"df_xi_qa_xi_dca_z_topo", "xi dca z topo", 1000, -500, 500}, "fXiDCAzToPVTopo");  
  
  auto h_df_xi_qa_xi_dca_xy_stra = df_xi_im.Histo1D({"df_xi_qa_xi_dca_xy_stra", "xi dca xy stra", 1000, -500, 500}, "fXiDCAxyToPVStraTrack");  
  auto h_df_xi_qa_xi_dca_z_stra = df_xi_im.Histo1D({"df_xi_qa_xi_dca_z_stra", "xi dca z stra", 1000, -500, 500}, "fXiDCAzToPVStraTrack");  
  

  auto out_counter_c1 = df_xi_im.Count(); 

  std::cout <<"==========================================================================================================\n";
  std::cout <<"You might want to close your eyes in case you are not really interesterd .... Very not nice cut overview: \n";
  std::cout <<"==========================================================================================================\n";
  
  auto cutsRepo = df.Report(); 
  cutsRepo->Print(); 
 
   
  std::cout <<"==========================================================================================================\n";
  std::cout <<"You can open your eyes again! \n";
  std::cout <<"==========================================================================================================\n";
  
 
  auto cutCounter  = new TH1D("cutCounter", "cutCounter", 20, 0, 20); 
  double inputCounter = *in_counter; 

  double reduction = *out_counter_c1/inputCounter;  
  std::cout << reduction << std::endl; 
  cutCounter->SetBinContent(1, reduction);
  
  TString outName = TString::Format("outxiccSelector_%s.root",outAddon )  ; 
  TFile* out = TFile::Open(outName.Data(), "recreate");
  
  out->cd(); 
  //this is my start
  HarryPlotter::CheckAndStore(out, h_df_in_im_xi_mass_stra);
  //to the lmb
  HarryPlotter::CheckAndStore(out, h_df_in_qa_pos_pt);
  HarryPlotter::CheckAndStore(out, h_df_in_qa_pos_dca_xy);
  HarryPlotter::CheckAndStore(out, h_df_in_qa_pos_dca_z);
  HarryPlotter::CheckAndStore(out, h_df_in_qa_pos_dca_xy_wide);
  HarryPlotter::CheckAndStore(out, h_df_in_qa_pos_dca_z_wide);
  HarryPlotter::CheckAndStore(out, h_df_in_qa_pos_hits);
  HarryPlotter::CheckAndStore(out, h_df_in_qa_pos_chisq);
  HarryPlotter::CheckAndStore(out, h_df_in_qa_pos_chisqhits);

  HarryPlotter::CheckAndStore(out, h_df_in_qa_neg_pt);
  HarryPlotter::CheckAndStore(out, h_df_in_qa_neg_dca_xy);
  HarryPlotter::CheckAndStore(out, h_df_in_qa_neg_dca_z);
  HarryPlotter::CheckAndStore(out, h_df_in_qa_neg_dca_xy_wide);
  HarryPlotter::CheckAndStore(out, h_df_in_qa_neg_dca_z_wide);
  HarryPlotter::CheckAndStore(out, h_df_in_qa_neg_hits);
  HarryPlotter::CheckAndStore(out, h_df_in_qa_neg_chisq);
  HarryPlotter::CheckAndStore(out, h_df_in_qa_neg_chisqhits);
  
  HarryPlotter::CheckAndStore(out, h_df_in_qa_lmb_dca_xy);
  HarryPlotter::CheckAndStore(out, h_df_in_qa_lmb_dca_z);
  HarryPlotter::CheckAndStore(out, h_df_in_qa_lmb_dca_xy_wide);
  HarryPlotter::CheckAndStore(out, h_df_in_qa_lmb_dca_z_wide);
  
  HarryPlotter::CheckAndStore(out, h_df_in_qa_lmb_totp);
  HarryPlotter::CheckAndStore(out, h_df_in_qa_lmb_ddist_pv);
  HarryPlotter::CheckAndStore(out, h_df_in_qa_lmb_ddca);
  HarryPlotter::CheckAndStore(out, h_df_in_qa_lmb_ddca_wide);
  HarryPlotter::CheckAndStore(out, h_df_in_qa_lmb_trad); 
  HarryPlotter::CheckAndStore(out, h_df_in_qa_lmb_trad_mc);
  HarryPlotter::CheckAndStore(out, h_df_in_qa_lmb_trad_diff);
  HarryPlotter::CheckAndStore(out, h_df_in_qa_lmb_mass);
  HarryPlotter::CheckAndStore(out, h_df_in_qa_xi_mass);
  
  //to xi 
  HarryPlotter::CheckAndStore(out, h_df_lmb_im_lmb_mass);
  
  HarryPlotter::CheckAndStore(out, h_df_lmb_qa_bach_pt);
  HarryPlotter::CheckAndStore(out, h_df_lmb_qa_bach_dca_xy);
  HarryPlotter::CheckAndStore(out, h_df_lmb_qa_bach_dca_z);
  HarryPlotter::CheckAndStore(out, h_df_lmb_qa_bach_dca_xy_wide);
  HarryPlotter::CheckAndStore(out, h_df_lmb_qa_bach_dca_z_wide);
  HarryPlotter::CheckAndStore(out, h_df_lmb_qa_bach_hits);
  HarryPlotter::CheckAndStore(out, h_df_lmb_qa_bach_chisq);
  HarryPlotter::CheckAndStore(out, h_df_lmb_qa_bach_chisqhits);
  
  HarryPlotter::CheckAndStore(out, h_df_lmb_qa_xi_pt);
  HarryPlotter::CheckAndStore(out, h_df_lmb_qa_xi_ddist_pv);
  HarryPlotter::CheckAndStore(out, h_df_lmb_qa_xi_ddca);
  HarryPlotter::CheckAndStore(out, h_df_lmb_qa_xi_trad); 
  HarryPlotter::CheckAndStore(out, h_df_lmb_qa_xi_trad_mc);
  HarryPlotter::CheckAndStore(out, h_df_lmb_qa_xi_trad_diff);
  
  HarryPlotter::CheckAndStore(out, h_df_lmb_qa_xi_mass);

  HarryPlotter::CheckAndStore(out, h_df_xi_im_xi_ddist_dv_topo); 
  HarryPlotter::CheckAndStore(out, h_df_xi_im_xi_ddist_dv_stra); 
  
  HarryPlotter::CheckAndStore(out, h_df_xi_sel_xi_mass);
  
  HarryPlotter::CheckAndStore(out, h_df_xi_qa_xi_dca_xy_topo);
  HarryPlotter::CheckAndStore(out, h_df_xi_qa_xi_dca_z_topo);
   
  HarryPlotter::CheckAndStore(out, h_df_xi_qa_xi_dca_xy_stra);
  HarryPlotter::CheckAndStore(out, h_df_xi_qa_xi_dca_z_stra);

  HarryPlotter::CheckAndStore(out, h_cand_counter); 
  HarryPlotter::CheckAndStore(out, h_gen_xi_pt_eta_counter); 
  HarryPlotter::CheckAndStore(out, h_gen_xi_pt_y_counter); 
  
  HarryPlotter::CheckAndStore(out, cutCounter);
  
  HarryPlotter::CheckAndStore(out, h_df_in_qa_pos_dca_xy_vs_hits);   
  HarryPlotter::CheckAndStore(out, h_df_in_qa_pos_dca_xy_vs_lmb_dl_pv);   
  HarryPlotter::CheckAndStore(out, h_df_in_qa_pos_dca_xy_vs_pos_pt);   
  
  HarryPlotter::CheckAndStore(out, h_df_in_qa_neg_dca_xy_vs_hits);   
  HarryPlotter::CheckAndStore(out, h_df_in_qa_neg_dca_xy_vs_lmb_dl_pv);   
  HarryPlotter::CheckAndStore(out, h_df_in_qa_neg_dca_xy_vs_neg_pt);   
  
  HarryPlotter::CheckAndStore(out, h_df_in_qa_neg_hits_vs_pos_hits);   
  HarryPlotter::CheckAndStore(out, h_df_in_qa_neg_dca_xy_vs_pos_dca_xy);   

  HarryPlotter::CheckAndStore(out, h_df_in_qa_bach_dca_xy_vs_hits);   
  HarryPlotter::CheckAndStore(out, h_df_in_qa_bach_dca_xy_vs_lmb_dl_pv);   
  HarryPlotter::CheckAndStore(out, h_df_in_qa_bach_dca_xy_vs_bach_pt);   
  
  out->Close(); 
  timer.Stop();
  timer.Print();
  return 0; 
}
