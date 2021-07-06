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
  
  int xiccDec = argv[3]?atoi(argv[3]):0; 
  bool ExclusiveSignal = (xiccDec == 0)?false:true;
  
  int pTbin = argv[4]?atoi(argv[4]):-1; 
  if (pTbin >= 0 && pTbin > ptbins.size()-1) { 
    std::cout << "Crashing cause the requested pT bin is out of range ( requested = " << pTbin << " , max. available = " << ptbins.size()-1 << " )\n";
    return -999; 
  }
  int wrongAssociationMode = argv[5]?atoi(argv[5]):0; 
  
  int noXi = argv[6]?atoi(argv[6]):0; 
  bool ForceNoXi = (noXi==0)?false:true; 
  
  HarryPlotter::StyleBox(); 
 
  ROOT::EnableImplicitMT(); // Tell ROOT you want to go parallel          
  TString filePath = TString::Format("%s", fileName); 
  
  auto h_cand_counter = new TH1D("df_xi_c_candCounter", "candCounter", 1, 0, 1); 
  auto h_gen_xi_c_counter = new TH1D("ptXicGen", "candCounter", 200, 0, 20); 
  auto h_gen_xi_c_pt_eta_counter = new TH2D("ptetaXicGen", "candCounter", 200, 0, 20, 30, -1.5, 1.5); 
  auto h_gen_xi_c_pt_y_counter = new TH2D("ptyXicGen", "candCounter", 200, 0, 20, 30, -1.5, 1.5); 
    
  auto h_gen_xi_cc_counter = new TH1D("ptXiccGen", "candCounter", 200, 0, 20); 
  auto h_gen_xi_cc_pt_eta_counter = new TH2D("ptetaXiccGen", "candCounter", 200, 0, 20, 30, -1.5, 1.5); 
  auto h_gen_xi_cc_pt_y_counter = new TH2D("ptyXiccGen", "candCounter", 200, 0, 20, 30, -1.5, 1.5); 
   
  TChain input("fTreeCandidates"); 
  int inputFiles = 0; 
  int inputFailures = 0; 
  if (filePath.Contains(".root")) { 
    TFile *inFile = TFile::Open(filePath);
    TH1D* evtCounter = (TH1D*)inFile->Get("hEventCounter"); 
    TH1D* ptXicGen = (TH1D*)inFile->Get("hXiCGeneratedPt"); 
    TH1D* ptXiccGen = (TH1D*)inFile->Get("hXiCCGeneratedPt"); 
    TH2D* ptXicpteta = (TH2D*)inFile->Get("hPtEtaGeneratedXiC"); 
    TH2D* ptXicpty = (TH2D*)inFile->Get("hPtYGeneratedXiC"); 
    TH2D* ptXiccpteta = (TH2D*)inFile->Get("hPtEtaGeneratedXiCC"); 
    TH2D* ptXiccpty = (TH2D*)inFile->Get("hPtYGeneratedXiCC"); 

    if (!evtCounter||!ptXicGen||!ptXicpteta||!ptXicpty||!ptXiccGen||!ptXiccpteta||!ptXiccpty) { 
      inputFailures++; 
      std::cout << "Zis is vehry bad, ze generation histograms are missing, Guenther! No histogram, no chain! \n"; 
    } else { 
      h_cand_counter->Add(evtCounter); 
      h_gen_xi_c_counter->Add(ptXicGen); 
      h_gen_xi_cc_counter->Add(ptXiccGen); 
	
      h_gen_xi_c_pt_eta_counter->Add(ptXicpteta); 
      h_gen_xi_c_pt_y_counter->Add(ptXicpty); 
	
      h_gen_xi_cc_pt_eta_counter->Add(ptXiccpteta); 
      h_gen_xi_cc_pt_y_counter->Add(ptXiccpty); 
	
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
	TH1D* ptXicGen = (TH1D*)inFile->Get("hXiCGeneratedPt"); 
	TH1D* ptXiccGen = (TH1D*)inFile->Get("hXiCCGeneratedPt"); 
	TH2D* ptXicpteta = (TH2D*)inFile->Get("hPtEtaGeneratedXiC"); 
	TH2D* ptXicpty = (TH2D*)inFile->Get("hPtYGeneratedXiC"); 
	TH2D* ptXiccpteta = (TH2D*)inFile->Get("hPtEtaGeneratedXiCC"); 
	TH2D* ptXiccpty = (TH2D*)inFile->Get("hPtYGeneratedXiCC"); 

	if (!evtCounter||!ptXicGen||!ptXicpteta||!ptXicpty||!ptXiccGen||!ptXiccpteta||!ptXiccpty) { 
	  inputFailures++;
	  continue; 
	}
	h_cand_counter->Add(evtCounter); 
	h_gen_xi_c_counter->Add(ptXicGen); 
	h_gen_xi_cc_counter->Add(ptXiccGen); 
	
	h_gen_xi_c_pt_eta_counter->Add(ptXicpteta); 
	h_gen_xi_c_pt_y_counter->Add(ptXicpty); 
	
	h_gen_xi_cc_pt_eta_counter->Add(ptXiccpteta); 
	h_gen_xi_cc_pt_y_counter->Add(ptXiccpty); 
	
	input.Add(inSubDirFile);
	inputFiles++;
	inFile->Close();
      }
    }
  }
  
  std::cout << "Added " << inputFiles << " files to the chain, failed in " << inputFailures << " cases \n"; 

  //TOF things
  float clight = 0.0299; //cm/ps
  auto beta = [&clight](float length, float time) { return length/(clight*time);};

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
  int addedHitsMin = ForceNoXi?0:1; 
  
  auto invMassXiCut = [&invMassDiffXi, &xiMass](float invMass) { return (TMath::Abs(invMass-xiMass) < invMassDiffXi); }; 
  auto hitsCut = [&addedHitsMin](int AddedHits) { return (AddedHits >= addedHitsMin); }; 
  
  auto decLengthXi = [&xiMass](float len, float mom){ return TMath::Abs(mom)>1e-4?len*xiMass/mom:-999; }; 

  //Xic cuts
  float xicMass = 2.468; 
  float invMassDiffXic = 0.12; //8 MeV/c2 mass window 
  auto decLengthXic = [&xicMass](float len, float mom){ return TMath::Abs(mom)>1e-4?len*xicMass/mom:-999; }; 
  
  auto invMassXicCut = [&invMassDiffXic, &xicMass](float invMass) { return (TMath::Abs(invMass-xicMass) < invMassDiffXic); }; 
  

  //Xicc cuts
  float xiccMass = 3.621;   
  float invMassDiffXicc = 0.120; //8 MeV/c2 mass window 
  auto decLengthXicc = [&xiccMass](float len, float mom){ return TMath::Abs(mom)>1e-4?len*xiccMass/mom:-999; }; 
  auto invMassXiccCut = [&invMassDiffXicc, &xiccMass](float invMass) { return (TMath::Abs(invMass-xiccMass) < invMassDiffXicc); }; 

  float xiccpTmin = pTbin >=0?ptbins[pTbin]:0; 
  float xiccpTmax = pTbin >=0?ptbins[pTbin+1]:20.0; 
  
  std::cout << "Xicc pT range set to " << xiccpTmin << " to " << xiccpTmax << std::endl; 
  
  auto pTCut = [&xiccpTmin, &xiccpTmax](float pT) { return (xiccpTmin < pT)&&(pT < xiccpTmax); }; 
  if (ExclusiveSignal) { 
    wrongAssociationMode = 0; 
  }
  auto associations = [&wrongAssociationMode](bool isXicc, bool isXic, bool piccUsed, int picUsed) { 
    bool out = true; 
    if (wrongAssociationMode == 1) { 
      //pions come from the correct decay chain but are not correctly associated
      out = (!isXicc && !isXic && piccUsed & (picUsed == 2)); 
    } else if (wrongAssociationMode == 2) { 
      //Xic and its pions are good but the pion from the xicc is replaced by a pythia pion
      out = (!isXicc && isXic && !piccUsed & (picUsed == 2)); 
    } else if (wrongAssociationMode == 3) { 
      //Xic and Xicc are wrong but the two pions from the xic are used only the xicc pion is replaced by a pythia pion
      out = (!isXicc && !isXic && !piccUsed & (picUsed == 2)); 
    } else if (wrongAssociationMode == 4) { 
      //Xic and Xicc are wrong but one xic pion and the xicc pion is replaced by a pythia pion
      out = (!isXicc && !isXic && !piccUsed & (picUsed == 1)); 
    } else if (wrongAssociationMode == 5) { 
      //Xic and Xicc are wrong but both xic pions are replaced by pythia pions
      out = (!isXicc && !isXic && piccUsed & (picUsed == 0)); 
    } else if (wrongAssociationMode == 6) { 
      //All three pions are replaced by pythia pions 
      out = (!isXicc && !isXic && !piccUsed & (picUsed == 0)); 
    } 
    return out;
  }; 

  ROOT::RDataFrame df(input);
  
  if (ExclusiveSignal) {
    std::cout << "Looking Exclusively for Signal!\n"; 
  } else { 
    std::cout << "Rejecting Signal!\n"; 
  }

  std::cout << "Wrong association mode set to " << wrongAssociationMode << std::endl;

  if (ForceNoXi) { 
    std::cout << "Rejecting Xis, make sure you know what you doing!\n"; 
  }
  
  auto df_ForceXi = ForceNoXi?df.Filter("!fTrueXi","noTrueXis"):df.Filter("fTrueXi||!fTrueXi","TrueAndFalseXis");

  auto df_in = (
		ExclusiveSignal?
		df_ForceXi
		.Filter("fTrueXicc")
		:df_ForceXi
		.Filter("!fTrueXicc")
		)
    .Filter(associations,{"fTrueXicc", "fTrueXic", "fPiccUsed", "fPicUsed"})
    .Define("fXiccPDGMass", [&xiccMass]() {return xiccMass;})
    .Define("fXiccY", HarryPlotter::YFromMomentum, {"lPXiCCStraTrack", "lPtXiCCStraTrack", "fXiccPDGMass", "fXiCCEta"})
    .Define("fNegBeta", beta, {"fNegativeLength", "fNegativeTOFSignal"})
    .Define("fPosBeta", beta, {"fPositiveLength", "fPositiveTOFSignal"})
    .Define("fBachBeta", beta, {"fBachelorLength", "fBachelorTOFSignal"}) 
    .Define("fNegExpectedBeta", beta, {"fNegativeLength", "fNegativeExpectedSignal"})
    .Define("fPosExpectedBeta", beta, {"fPositiveLength", "fPositiveExpectedSignal"})
    .Define("fBachExpectedBeta", beta, {"fBachelorLength", "fBachelorExpectedSignal"})     // .Define("fPic1Beta", beta, {"fNegativeLength", "fNegativeTOFSignal"})
    .Filter("TMath::Abs(fXiCCEta)<0.5")
    .Filter(pTCut, {"lPtXiCCStraTrack"}, "pTXicc")
    ;

  auto h_df_in_im_xi_cc_mass_stra = df_in.Histo1D({"h_df_in_im_xi_cc_mass_stra", "xi_cc inv mass", 700, 2.6, 4.6}, "fXiccMassStraTrack"); 
  auto in_counter = df_in.Count(); 

  auto df_in_qa = df_in.Filter("fFirstCandidateXiCC","df_in_h_bool"); 

  //towards the Lambda for the Xi 
  auto h_df_in_qa_pos_pt = df_in_qa.Histo1D({"df_in_qa_pos_pt", "pos pt", 200, 0, 20}, "fPositivePt"); 
  auto h_df_in_qa_pos_dca_xy = df_in_qa.Histo1D({"df_in_qa_pos_dca_xy", "pos dca xy pv", 1000, -1e4, 1e4}, "fPositiveDCAxy"); 
  auto h_df_in_qa_pos_dca_z = df_in_qa.Histo1D({"df_in_qa_pos_dca_z", "pos dca z pv", 1000, -1e4, 1e4}, "fPositiveDCAz"); 
  auto h_df_in_qa_pos_dca_xy_wide = df_in_qa.Histo1D({"df_in_qa_pos_dca_xy_wide", "pos dca xy pv", 1000, -1e5, 1e5}, "fPositiveDCAxy"); 
  auto h_df_in_qa_pos_dca_z_wide = df_in_qa.Histo1D({"df_in_qa_pos_dca_z_wide", "pos dca z pv", 1000, -1e5, 1e5}, "fPositiveDCAz"); 
  auto h_df_in_qa_pos_hits = df_in_qa.Histo1D({"df_in_qa_pos_hits", "pos hits", 15, 0, 15}, "fPositiveClusters"); 
  auto h_df_in_qa_pos_chisq = df_in_qa.Histo1D({"df_in_qa_pos_chisq", "pos chisq", 200, 0, 200}, "fPositiveChisquare"); 
  auto h_df_in_qa_pos_chisqhits = df_in_qa.Define("fPositiveChisquareOverHits", "fPositiveChisquare/fPositiveClusters").Histo1D({"df_in_qa_pos_chisqhits", "pos chisq over hits", 100, 0, 50}, "fPositiveChisquareOverHits"); 

  auto h_df_in_qa_pos_dca_xy_vs_hits = df_in_qa.Histo2D<float,int>({"df_in_qa_pos_dca_xy_vs_hits", "pos dca xy pv vs hits", 1000, -10000, 10000, 15, -0.5, 14.5}, "fPositiveDCAxy", "fPositiveClusters"); 
  auto h_df_in_qa_pos_dca_xy_vs_lmb_dl_pv = df_in_qa.Define("fLmbInvDecayLengthToPV", decLengthLmb, {"fV0DecayLength", "fV0TotalMomentum"}).Histo2D<float,float>({"df_in_qa_pos_dca_xy_vs_lmb_dl_pv", "pos dca xy pv vs lmb decay length", 1000, -10000, 10000, 100, 0, 40}, "fPositiveDCAxy", "fLmbInvDecayLengthToPV"); 
  auto h_df_in_qa_pos_dca_xy_vs_pos_pt = df_in_qa.Histo2D<float,float>({"h_df_in_qa_pos_dca_xy_vs_pos_pt", "pos dca xy pv vs pos pt", 1000, -10000, 10000, 200, 0, 20}, "fPositiveDCAxy", "fPositivePt"); 

  auto h_df_in_qa_pos_beta_vs_pos_pt = df_in_qa.Histo2D({"h_df_in_qa_pos_beta_vs_pos_pt", "pos beta vs pos pt", 200, 0, 20, 1200, 0, 1.2}, "fPositivePt", "fPosBeta"); 

  auto h_df_in_qa_neg_pt = df_in_qa.Histo1D({"df_in_qa_neg_pt", "neg pt", 200, 0, 20}, "fNegativePt"); 
  auto h_df_in_qa_neg_dca_xy = df_in_qa.Histo1D({"df_in_qa_neg_dca_xy", "neg dca xy pv", 1000, -1e4, 1e4}, "fNegativeDCAxy"); 
  auto h_df_in_qa_neg_dca_z = df_in_qa.Histo1D({"df_in_qa_neg_dca_z", "neg dca z pv", 1000, -1e4, 1e4}, "fNegativeDCAz"); 
  auto h_df_in_qa_neg_dca_xy_wide = df_in_qa.Histo1D({"df_in_qa_neg_dca_xy_wide", "neg dca xy pv", 1000, -1e5, 1e5}, "fNegativeDCAxy"); 
  auto h_df_in_qa_neg_dca_z_wide = df_in_qa.Histo1D({"df_in_qa_neg_dca_z_wide", "neg dca z pv", 1000, -1e5, 1e5}, "fNegativeDCAz"); 
  auto h_df_in_qa_neg_hits = df_in_qa.Histo1D({"df_in_qa_neg_hits", "neg hits", 15, 0, 15}, "fNegativeClusters"); 
  auto h_df_in_qa_neg_chisq = df_in_qa.Histo1D({"df_in_qa_neg_chisq", "neg chisq", 200, 0, 200}, "fNegativeChisquare"); 
  auto h_df_in_qa_neg_chisqhits = df_in_qa.Define("fNegativeChisquareOverHits", "fNegativeChisquare/fNegativeClusters").Histo1D({"df_in_qa_neg_chisqhits", "neg chisq over hits", 100, 0, 50}, "fNegativeChisquareOverHits"); 
  
  auto h_df_in_qa_neg_dca_xy_vs_hits = df_in_qa.Histo2D<float,int>({"df_in_qa_neg_dca_xy_vs_hits", "neg dca xy pv vs hits", 1000, -10000, 10000, 15, -0.5, 14.5}, "fNegativeDCAxy", "fNegativeClusters"); 
  auto h_df_in_qa_neg_dca_xy_vs_lmb_dl_pv = df_in_qa.Define("fLmbInvDecayLengthToPV", decLengthLmb, {"fV0DecayLength", "fV0TotalMomentum"}).Histo2D<float,float>({"df_in_qa_neg_dca_xy_vs_lmb_dl_pv", "neg dca xy pv vs lmb decay length", 1000, -10000, 10000, 100, 0, 40}, "fNegativeDCAxy", "fLmbInvDecayLengthToPV"); 
  auto h_df_in_qa_neg_dca_xy_vs_neg_pt = df_in_qa.Histo2D<float,float>({"h_df_in_qa_neg_dca_xy_vs_neg_pt", "neg dca xy pv vs neg pt", 1000, -10000, 10000, 200, 0, 20}, "fNegativeDCAxy", "fNegativePt"); 
  auto h_df_in_qa_neg_hits_vs_pos_hits = df_in_qa.Histo2D<int,int>({"h_df_in_qa_neg_hits_vs_pos_hits", "neg hits vs pos hits", 15, -0.5, 14.5, 15, -0.5, 14.5}, "fNegativeClusters", "fPositiveClusters"); 
  auto h_df_in_qa_neg_dca_xy_vs_pos_dca_xy = df_in_qa.Histo2D<float,float>({"df_in_qa_neg_dca_xy_vs_pos_dca_xy", "neg dca xy pv vs pos dca xy pv", 1000, -10000, 10000, 1000, -10000, 10000}, "fNegativeDCAxy", "fPositiveDCAxy"); 
  auto h_df_in_qa_neg_beta_vs_neg_pt = df_in_qa.Histo2D({"h_df_in_qa_neg_beta_vs_neg_pt", "neg beta vs neg pt", 200, 0, 20, 1200, 0, 1.2}, "fNegativePt", "fNegBeta"); 
  auto h_df_in_qa_lmb_dca_xy = df_in_qa.Histo1D({"df_in_qa_lmb_dca_xy", "lmb dca xy pv", 1000, -1e4, 1e4}, "fV0DCAxyToPV");
  auto h_df_in_qa_lmb_dca_z = df_in_qa.Histo1D({"df_in_qa_lmb_dca_z", "lmb dca z pv", 1000, -1e4, 1e4}, "fV0DCAzToPV");
  auto h_df_in_qa_lmb_dca_xy_wide = df_in_qa.Histo1D({"df_in_qa_lmb_dca_xy_wide", "lmb dca xy pv", 1000, -1e5, 1e5}, "fV0DCAxyToPV");
  auto h_df_in_qa_lmb_dca_z_wide = df_in_qa.Histo1D({"df_in_qa_lmb_dca_z_wide", "lmb dca z pv", 1000, -1e5, 1e5}, "fV0DCAzToPV");
  
  auto h_df_in_qa_lmb_totp = df_in_qa.Histo1D({"df_in_qa_lmb_totp", "lmb tot p", 200, 0, 20}, "fV0TotalMomentum");   
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
    .Define("DCAxyProd"  , "fXiDCAxyToPVStraTrack*fXicPionDCAxyToPV1*fXicPionDCAxyToPV2*fXicDCAxyToPVStraTrack*fPicDCAxyToPVStraTrack") //O(1e7)
    .Define("DCAxySqSum"  , "fXiDCAxyToPVStraTrack*fXiDCAxyToPVStraTrack+fXicPionDCAxyToPV1*fXicPionDCAxyToPV1+fXicPionDCAxyToPV2*fXicPionDCAxyToPV2+fXicDCAxyToPVStraTrack*fXicDCAxyToPVStraTrack+fPicDCAxyToPVStraTrack*fPicDCAxyToPVStraTrack") //O(1e7)
    .Define("DCAzProd"   , "fXiDCAzToPVStraTrack*fXicPionDCAzToPV1*fXicPionDCAzToPV2*fXicDCAzToPVStraTrack*fPicDCAzToPVStraTrack")
    .Define("DCAzSqSum"  , "fXiDCAzToPVStraTrack*fXiDCAzToPVStraTrack+fXicPionDCAzToPV1*fXicPionDCAzToPV1+fXicPionDCAzToPV2*fXicPionDCAzToPV2+fXicDCAzToPVStraTrack*fXicDCAzToPVStraTrack+fPicDCAzToPVStraTrack*fPicDCAzToPVStraTrack") //O(1e7)
    .Define("DCAProdAdd" , "DCAxyProd*DCAxyProd+DCAzProd*DCAzProd")
    .Define("DCAProdMult", "DCAxyProd*DCAzProd")
    .Define("DCAxyProdXicPi", "fPicDCAxyToPVStraTrack*fXicDCAxyToPVStraTrack")
    .Define("DCAzProdXicPi", "fPicDCAzToPVStraTrack*fXicDCAzToPVStraTrack")
    .Define("fLmbInvDecayLengthToPV", decLengthLmb, {"fV0DecayLength", "fV0TotalMomentum"})
    .Define("fXicInvDecayLengthToPVTopo", decLengthXic, {"fXicDecayDistanceFromPVTopo", "lPXiCTopo"}) //Xi_c decay length from the PV topological
    .Define("fXiccInvDecayLengthToPVTopo", decLengthXicc, {"fXiccDecayDistanceFromPVTopo", "lPXiCCTopo"}) // Xi_cc decay length from the PV topological
    .Define("fXicInvDecayLengthToPVStra", decLengthXic, {"fXicDecayDistanceFromPVStraTrack", "lPXiCStraTrack"})   // Xi_c decay length from the PV stra tracking 
    .Define("fXiccInvDecayLengthToPVStra", decLengthXicc, {"fXiccDecayDistanceFromPVStraTrack", "lPXiCCStraTrack"}) // Xi_cc decay length form the PV stra tracking = this is the Xi_cc decay length! 
    .Define("fXiInvDecayLengthToPV", decLengthXi, {"fXiDecayLength", "fXiTotalMomentum"})
    .Define("fXiInvDecayLengthToDVTopo", decLengthXi, {"fXiCtoXiLengthTopo", "fXiTotalMomentum"}) //this is the Xi decay length 
    .Define("fXicInvDecayLengthToDVTopo", decLengthXic, {"fXiCCtoXiCLengthTopo", "lPXiCTopo"}) //this is the Xi_c decay length 
    .Define("fXiInvDecayLengthToDVStra", decLengthXi, {"fXiCtoXiLengthStraTrack", "fXiTotalMomentum"}) //this is the Xi decay length 
    .Define("fXicInvDecayLengthToDVStra", decLengthXic, {"fXiCCtoXiCLengthStraTrack", "lPXiCStraTrack"})    //this is the Xi_c decay length 
    .Define("fPosExpPVDV", "fPositiveExpectedSignal - fPositiveExpectedSignalFromPV")
    .Define("fNegExpPVDV", "fNegativeExpectedSignal - fNegativeExpectedSignalFromPV")
    .Define("fBachExpPVDV", "fBachelorExpectedSignal - fBachelorExpectedSignalFromPV")
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

  auto h_df_lmb_im_lmb_mass = df_lmb_im.Filter("fFirstCandidateXiCC","df_lmb_im_h_bool").Histo1D({"df_lmb_im_lmb_mass", "lmb inv mass", 750, 1., 1.8}, "fLambdaMass"); 

  auto df_lmb =  df_lmb_im
    .Filter(invMassLmbCut, {"fLambdaMass"}, "fLambdaMass")
    ;

  //Towards the Xi for Xic->Xi+2pi 
  //Fill some Histograms
  auto h_df_lmb_xi_cc_mass_stra = df_lmb.Histo1D({"h_df_lmb_xi_cc_mass_stra", "xi_cc inv mass", 700, 2.6, 4.6}, "fXiccMassStraTrack"); 

  auto df_lmb_qa = df_lmb.Filter("fFirstCandidateXiCC", "df_lmb_h_bool");
  
  auto h_df_lmb_qa_pos_beta_expected_vs_pos_pt = df_lmb_qa.Histo2D({"h_df_lmb_qa_pos_beta_expected_vs_pos_pt", "pos beta exp vs pos pt", 200, 0, 20, 1200, 0, 1.2}, "fPositivePt", "fPosExpectedBeta"); 
  auto h_df_lmb_qa_pos_tof_dv_expected_vs_measured = df_lmb_qa.Histo2D({"h_df_lmb_qa_pos_tof_dv_expected_vs_measured", "beta expected vs measured", 2000, 2500, 6500, 2000, 2500, 6500}, "fPositiveTOFSignal", "fPositiveExpectedSignal"); 
  auto h_df_lmb_qa_pos_tof_pv_expected_vs_measured = df_lmb_qa.Histo2D({"h_df_lmb_qa_pos_tof_pv_expected_vs_measured", "beta expected vs measured", 2000, 2500, 6500, 2000, 2500, 6500}, "fPositiveTOFSignal", "fPositiveExpectedSignalFromPV"); 
  auto h_df_lmb_qa_pos_tof_dv_expected_vs_pv_expected = df_lmb_qa.Histo2D({"h_df_lmb_qa_pos_tof_dv_expected_vs_pv_expected", "beta expected vs measured", 2000, 2500, 6500, 2000, 2500, 6500}, "fPositiveExpectedSignal","fPosExpPVDV");

  auto h_df_lmb_qa_neg_beta_expected_vs_neg_pt = df_lmb_qa.Histo2D({"h_df_lmb_qa_neg_beta_expected_vs_neg_pt", "neg beta exp vs neg pt", 200, 0, 20, 1200, 0, 1.2}, "fNegativePt", "fNegExpectedBeta"); 
  auto h_df_lmb_qa_neg_tof_dv_expected_vs_measured = df_lmb_qa.Histo2D({"h_df_lmb_qa_neg_tof_dv_expected_vs_measured", "beta expected vs measured", 2000, 2500, 6500, 2000, 2500, 6500}, "fNegativeTOFSignal", "fNegativeExpectedSignal"); 
  auto h_df_lmb_qa_neg_tof_pv_expected_vs_measured = df_lmb_qa.Histo2D({"h_df_lmb_qa_neg_tof_pv_expected_vs_measured", "beta expected vs measured", 2000, 2500, 6500, 2000, 2500, 6500}, "fNegativeTOFSignal", "fNegativeExpectedSignalFromPV"); 
  auto h_df_lmb_qa_neg_tof_dv_expected_vs_pv_expected = df_lmb_qa.Histo2D({"h_df_lmb_qa_neg_tof_dv_expected_vs_pv_expected", "beta expected vs measured", 2000, 2500, 6500, 2000, 2500, 6500}, "fNegativeExpectedSignal","fNegExpPVDV");

  auto h_df_lmb_qa_bach_beta_vs_bach_pt = df_lmb_qa.Histo2D({"h_df_lmb_qa_bach_beta_vs_bach_pt", "bach beta vs bach pt", 200, 0, 20, 1200, 0, 1.2}, "fBachelorPt", "fBachBeta"); 
  auto h_df_lmb_qa_bach_beta_expected_vs_bach_pt = df_lmb_qa.Histo2D({"h_df_lmb_qa_bach_beta_expected_vs_bach_pt", "bach beta exp vs bach pt", 200, 0, 20, 1200, 0, 1.2}, "fBachelorPt", "fBachExpectedBeta"); 
  auto h_df_lmb_qa_bach_tof_dv_expected_vs_measured = df_lmb_qa.Histo2D({"h_df_lmb_qa_bach_tof_dv_expected_vs_measured", "beta expected vs measured", 2000, 2500, 6500, 2000, 2500, 6500}, "fBachelorTOFSignal", "fBachelorExpectedSignal"); 
  auto h_df_lmb_qa_bach_tof_pv_expected_vs_measured = df_lmb_qa.Histo2D({"h_df_lmb_qa_bach_tof_pv_expected_vs_measured", "beta expected vs measured", 2000, 2500, 6500, 2000, 2500, 6500}, "fBachelorTOFSignal", "fBachelorExpectedSignalFromPV"); 
  auto h_df_lmb_qa_bach_tof_dv_expected_vs_pv_expected = df_lmb_qa.Histo2D({"h_df_lmb_qa_bach_tof_dv_expected_vs_pv_expected", "beta expected vs measured", 2000, 2500, 6500, 2000, 2500, 6500}, "fBachelorExpectedSignal", "fBachExpPVDV"); 


  auto h_df_lmb_qa_bach_pt = df_lmb_qa.Histo1D({"df_lmb_qa_bach_pt", "bach pt", 200, 0, 20}, "fBachelorPt"); 
  auto h_df_lmb_qa_bach_dca_xy = df_lmb_qa.Histo1D({"df_lmb_qa_bach_dca_xy", "bach dca xy pv", 1000, -1e4, 1e4}, "fBachelorDCAxy"); 
  auto h_df_lmb_qa_bach_dca_z = df_lmb_qa.Histo1D({"df_lmb_qa_bach_dca_z", "bach dca z pv", 1000, -1e4, 1e4}, "fBachelorDCAz"); 
  auto h_df_lmb_qa_bach_dca_xy_wide = df_lmb_qa.Histo1D({"df_lmb_qa_bach_dca_xy_wide", "bach dca xy pv", 1000, -1e5, 1e5}, "fBachelorDCAxy"); 
  auto h_df_lmb_qa_bach_dca_z_wide = df_lmb_qa.Histo1D({"df_lmb_qa_bach_dca_z_wide", "bach dca z pv", 1000, -1e5, 1e5}, "fBachelorDCAz"); 
  auto h_df_lmb_qa_bach_hits = df_lmb_qa.Histo1D({"df_lmb_qa_bach_hits", "bach hits", 15, 0, 15}, "fBachelorClusters"); 
  auto h_df_lmb_qa_bach_chisq = df_lmb_qa.Histo1D({"df_lmb_qa_bach_chisq", "bach chisq", 200, 0, 200}, "fBachelorChisquare"); 
  auto h_df_lmb_qa_bach_chisqhits = df_lmb_qa.Define("fBachelorChisquareOverHits", "fBachelorChisquare/fBachelorClusters").Histo1D({"df_lmb_qa_bach_chisqhits", "bach chisq over hits", 100, 0, 50}, "fBachelorChisquareOverHits"); 

  auto h_df_lmb_qa_bach_dca_xy_vs_hits = df_lmb_qa.Histo2D<float,int>({"df_lmb_qa_bach_dca_xy_vs_hits", "bach dca xy pv vs hits", 1000, -10000, 10000, 15, -0.5, 14.5}, "fBachelorDCAxy", "fBachelorClusters"); 
  auto h_df_lmb_qa_bach_dca_xy_vs_lmb_dl_pv = df_lmb_qa.Histo2D<float,float>({"df_lmb_qa_bach_dca_xy_vs_lmb_dl_pv", "bach dca xy pv vs lmb decay length", 1000, -10000, 10000, 100, 0, 40}, "fBachelorDCAxy", "fLmbInvDecayLengthToPV"); 
  auto h_df_lmb_qa_bach_dca_xy_vs_bach_pt = df_lmb_qa.Histo2D<float,float>({"h_df_lmb_qa_bach_dca_xy_vs_bach_pt", "bach dca xy pv vs bach pt", 1000, -10000, 10000, 200, 0, 20}, "fBachelorDCAxy", "fBachelorPt"); 
  
    
  auto h_df_lmb_qa_xi_pt = df_lmb_qa.Histo1D({"df_lmb_qa_xi_pt", "xi pt", 200, 0, 20}, "fXiPtStraTrack"); 
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

  auto h_df_xi_sel_xi_mass = df_xi_sel.Filter("fFirstCandidateXiCC","df_xi_sel_h_bool").Histo1D({"df_xi_sel_xi_mass", "xi inv mass", 750, 1.2, 2}, "fXiMass"); 

  //Select Xis including Hits
  auto df_xi_im = df_xi_sel.
    Filter(hitsCut, {"fXiHitsAdded"},"fXiHitsAdded")
    ; 
  
  auto h_df_xi_im_xi_mass = df_xi_im.Filter("fFirstCandidateXiCC","df_xi_im_h_bool").Histo1D({"df_xi_im_xi_mass", "xi inv mass", 750, 1.2, 2}, "fXiMass"); 
  
  auto h_df_xi_im_xi_ddist_dv_topo = df_xi_im.Filter("fFirstCandidateXiCC","df_xi_im_h_bool").Histo1D({"df_xi_im_xi_dist_dv_topo", "xi decay dist", 1500, 0, 30}, "fXiInvDecayLengthToDVTopo"); 
  auto h_df_xi_im_xi_ddist_dv_stra = df_xi_im.Filter("fFirstCandidateXiCC","df_xi_im_h_bool").Histo1D({"df_xi_im_xi_dist_dv_stra", "xi decay dist", 1500, 0, 30}, "fXiInvDecayLengthToDVStra"); 

  auto df_xi = df_xi_im
    .Filter(invMassXiCut, {"fXiMass"},"fXiMass")
    ;
  auto h_df_xi_xi_cc_mass_stra = df_xi.Histo1D({"h_df_xi_xi_cc_mass_stra", "xi_cc inv mass", 700, 2.6, 4.6}, "fXiccMassStraTrack"); 

  //Towards the Xi_c for Xi_cc->Xi_c+pi: 
  //Fill some histograms 
  //Topo
  
  auto df_xi_qa = df_xi.Filter("fFirstCandidateXiCC","df_xi_h_bool");

  auto h_df_xi_qa_trad_diff_xi_xi_c_topo = df_xi_qa.Histo1D({"df_xi_qa_trad_diff_xi_xi_c_topo", "xi-xi_c trad", 500, -50, 100}, "XiXicDecayRadDiffTopo") ;

  auto h_df_xi_qa_xi_c_mass_topo = df_xi_qa.Histo1D({"df_xi_qa_xi_c_mass_topo", "xi_c inv mass", 700, 1.6, 3.2}, "fXicMassTopo"); 
  
  auto h_df_xi_qa_xi_c_ddca_topo = df_xi_qa.Histo1D({"df_xi_qa_xi_c_ddca_topo", "xi_c prong dca", 500, 0, 100}, "fXicDaughterDCATopo"); 
  auto h_df_xi_qa_xi_c_ddist_pv_topo = df_xi_qa.Histo1D({"df_xi_qa_xi_c_dist_pv_topo", "xi_c decay dist", 1500, 0, 0.30}, "fXicInvDecayLengthToPVTopo"); 
  auto h_df_xi_qa_xi_c_ddist_dv_topo = df_xi_qa.Histo1D({"df_xi_qa_xi_c_dist_dv_topo", "xi_c decay dist", 1500, 0, 0.30}, "fXicInvDecayLengthToDVTopo");   
  auto h_df_xi_qa_xi_c_trad_topo = df_xi_qa.Histo1D({"df_xi_qa_xi_c_trad_topo", "xi_c trad", 2000, 0, 0.4}, "fXicDecayRadiusTopo"); 
  
  auto h_df_xi_qa_xi_dca_xy_topo = df_xi_qa.Histo1D({"df_xi_qa_xi_dca_xy_topo", "xi dca xy topo", 1000, -500, 500}, "fXiDCAxyToPVTopo");  
  auto h_df_xi_qa_xi_dca_z_topo = df_xi_qa.Histo1D({"df_xi_qa_xi_dca_z_topo", "xi dca z topo", 1000, -500, 500}, "fXiDCAzToPVTopo");  

  //Stra

  auto h_df_xi_qa_trad_diff_xi_xi_c_stra = df_xi_qa.Histo1D({"df_xi_qa_trad_diff_xi_xi_c_stra", "xi-xi_c trad", 500, -50, 200}, "XiXicDecayRadDiffStra") ;
  
  auto h_df_xi_qa_xi_c_mass_stra = df_xi_qa.Histo1D({"df_xi_qa_xi_c_mass_stra", "xi_c inv mass", 700, 1.6, 3.2}, "fXicMassStraTrack"); 
  auto h_df_xi_qa_xi_c_pt = df_xi_qa.Histo1D({"df_xi_qa_xi_c_pt", "xi_c pt", 200, 0, 20}, "fXicPtStraTrack");  
  
  auto h_df_xi_qa_xi_c_ddca_stra = df_xi_qa.Histo1D({"df_xi_qa_xi_c_ddca_stra", "xi_c prong dca", 500, 0, 100}, "fXicDaughterDCAStraTrack"); 
  auto h_df_xi_qa_xi_c_ddist_pv_stra = df_xi_qa.Histo1D({"df_xi_qa_xi_c_dist_pv_stra", "xi_c decay dist", 1500, 0, 0.30}, "fXicInvDecayLengthToPVStra"); //redefine as mL/p 
  auto h_df_xi_qa_xi_c_ddist_dv_stra = df_xi_qa.Histo1D({"df_xi_qa_xi_c_dist_dv_stra", "xi_c decay dist", 1500, 0, 0.30}, "fXicInvDecayLengthToDVStra"); 
  auto h_df_xi_qa_xi_c_trad_stra = df_xi_qa.Histo1D({"df_xi_qa_xi_c_trad_stra", "xi_c trad", 2000, 0, 0.4}, "fXicDecayRadiusStraTrack"); 
  
  auto h_df_xi_qa_xi_dca_xy_stra = df_xi_qa.Histo1D({"df_xi_qa_xi_dca_xy_stra", "xi dca xy stra", 1000, -500, 500}, "fXiDCAxyToPVStraTrack");  
  auto h_df_xi_qa_xi_dca_z_stra = df_xi_qa.Histo1D({"df_xi_qa_xi_dca_z_stra", "xi dca z stra", 1000, -500, 500}, "fXiDCAzToPVStraTrack");  
  
  //Pions (from the xic) 
  
  auto h_df_xi_qa_pi_one_dca_xy = df_xi_qa.Histo1D({"df_xi_qa_pi_one_dca_xy", "pi1 dca xy stra", 1000, -500, 500}, "fXicPionDCAxyToPV1");  
  auto h_df_xi_qa_pi_one_dca_z = df_xi_qa.Histo1D({"df_xi_qa_pi_one_dca_z", "pi1 dca z stra", 1000, -500, 500}, "fXicPionDCAzToPV1");  
  
  auto h_df_xi_qa_pi_two_dca_xy = df_xi_qa.Histo1D({"df_xi_qa_pi_two_dca_xy", "pi2 dca xy stra", 1000, -500, 500}, "fXicPionDCAxyToPV2");  
  auto h_df_xi_qa_pi_two_dca_z = df_xi_qa.Histo1D({"df_xi_qa_pi_two_dca_z", "pi2 dca z stra", 1000, -500, 500}, "fXicPionDCAzToPV2");  
  
  auto h_df_xi_qa_pi_one_pt = df_xi_qa.Histo1D({"df_xi_qa_pi_one_pt", "pi c pt", 200, 0, 20}, "fPiC1Pt");  
  auto h_df_xi_qa_pi_two_pt = df_xi_qa.Histo1D({"df_xi_qa_pi_two_pt", "pi c pt", 200, 0, 20}, "fPiC2Pt");  

  //study a bit decay length + trad 
  //cut on Trad and check decay lengths to pv and to dv 
  auto h_df_xi_qa_trad_xi_c_ddist_pv_stra = df_xi_qa.Filter("fXicDecayRadiusStraTrack > 0.004","corrStudy_fXicDecayRadiusStraTrack").Histo1D({"df_xi_qa_trad_xi_c_dist_pv_stra", "xi_c decay dist", 1500, 0, 0.30}, "fXicInvDecayLengthToPVStra"); 
  auto h_df_xi_qa_trad_xi_c_ddist_dv_stra = df_xi_qa.Filter("fXicDecayRadiusStraTrack > 0.004","corrStudy_fXicDecayRadiusStraTrack").Histo1D({"df_xi_qa_trad_xi_c_dist_dv_stra", "xi_c decay dist", 1500, 0, 0.30}, "fXicInvDecayLengthToDVStra"); 
  
  //cut on dl to pv 
  auto h_df_xi_qa_dl_pv_xi_c_ddist_dv_stra =  df_xi_qa.Filter("fXicInvDecayLengthToPVStra > 0.002","corrStudy_fXicInvDecayLengthToPVStra").Histo1D({"df_xi_qa_dl_pv_xi_c_dist_dv_stra", "xi_c decay dist", 1500, 0, 0.30}, "fXicInvDecayLengthToDVStra");
  auto h_df_xi_qa_dl_pv_xi_c_trad_stra = df_xi_qa.Filter("fXicInvDecayLengthToPVStra > 0.002","corrStudy_fXicInvDecayLengthToPVStra").Histo1D({"df_xi_qa_dl_pv_xi_c_trad_stra", "xi_c trad", 2000, 0, 0.4}, "fXicDecayRadiusStraTrack");
  
  //cut on dl to dv 
  auto h_df_xi_qa_dl_dv_xi_c_ddist_dv_stra =  df_xi_qa.Filter("fXicInvDecayLengthToDVStra < 0.1","corrStudy_fXicInvDecayLengthToDVStra").Histo1D({"df_xi_qa_dl_dv_xi_c_dist_pv_stra", "xi_c decay dist", 1500, 0, 0.30}, "fXicInvDecayLengthToPVStra");
  auto h_df_xi_qa_dl_dv_xi_c_trad_stra = df_xi_qa.Filter("fXicInvDecayLengthToDVStra < 0.1","corrStudy_fXicInvDecayLengthToDVStra").Histo1D({"df_xi_qa_dl_dv_xi_c_trad_stra", "xi_c trad", 2000, 0, 0.4}, "fXicDecayRadiusStraTrack");

  //Select the Xi_c 

  auto df_xi_c_im = df_xi
    .Define("XicXiccDecayRadDiffTopo", "fXicDecayRadiusTopo-fXiccDecayRadiusTopo")
    .Define("XicXiccDecayRadDiffStra", "fXicDecayRadiusStraTrack-fXiccDecayRadiusStraTrack")
    .Filter("XiXicDecayRadDiffStra > 0","XiXicDecayRadDiffStra")
    ;
  
  auto h_df_xi_c_im_xi_c_mass_stra = df_xi_c_im.Filter("fFirstCandidateXiCC","df_xi_c_im_h_bool").Histo1D({"df_xi_c_im_xi_c_mass_stra", "xi_c inv mass", 700, 1.6, 3.2}, "fXicMassStraTrack"); 

  auto df_xi_c = df_xi_c_im.Filter(invMassXicCut, {"fXicMassStraTrack"},"fXicMassStraTrack");
  
  auto df_xi_c_qa = df_xi_c.Filter("fFirstCandidateXiCC","df_xi_c_h_bool"); 
  
  //Towards the actual Xi_cc
  //Fill some histograms 
  
  //Topo
  
  auto h_df_xi_c_qa_trad_diff_xi_xi_c_topo = df_xi_c_qa.Histo1D({"df_xi_c_qa_trad_diff_xi_c_xi_cc_topo", "xi_c-xi_cc trad", 500, -100, 150}, "XicXiccDecayRadDiffTopo") ;

  auto h_df_xi_c_qa_xi_cc_ddca_topo = df_xi_c_qa.Histo1D({"df_xi_c_qa_xi_cc_ddca_topo", "xi_cc prong dca", 500, 0, 500}, "fXiccDaughterDCATopo"); 
  auto h_df_xi_c_qa_xi_cc_ddist_pv_topo = df_xi_c_qa.Histo1D({"df_xi_c_qa_xi_cc_dist_pv_topo", "xi_cc decay dist", 3000, 0, 0.5}, "fXiccInvDecayLengthToPVTopo"); 
  auto h_df_xi_c_qa_xi_cc_trad_topo = df_xi_c_qa.Histo1D({"df_xi_c_qa_xi_cc_trad_topo", "xi_cc trad", 2000, 0, 0.5}, "fXiccDecayRadiusTopo"); 

  auto h_df_xi_c_qa_xi_c_dca_xy_topo = df_xi_c_qa.Histo1D({"df_xi_c_qa_xi_c_dca_xy_topo", "xi_c dca xy topo", 1000, -500, 500}, "fXicDCAxyToPVTopo");  
  auto h_df_xi_c_qa_xi_c_dca_z_topo  = df_xi_c_qa.Histo1D({"df_xi_c_qa_xi_c_dca_z_topo", "xi_c dca z topo", 1000, -500, 500}, "fXicDCAzToPVTopo");  

  auto h_df_xi_c_qa_xi_cc_dca_xy_topo = df_xi_c_qa.Histo1D({"df_xi_c_qa_xi_cc_dca_xy_topo", "xi_cc dca xy topo", 1000, -500, 500}, "fXiccDCAxyToPVTopo");  
  auto h_df_xi_c_qa_xi_cc_dca_z_topo  = df_xi_c_qa.Histo1D({"df_xi_c_qa_xi_cc_dca_z_topo", "xi_cc dca z topo", 1000, -500, 500}, "fXiccDCAzToPVTopo");  
  
  //Stra

  auto h_df_xi_c_qa_trad_diff_xi_xi_c_stra = df_xi_c_qa.Histo1D({"df_xi_c_qa_trad_diff_xi_c_xi_cc_stra", "xi_c-xi_cc trad", 500, -100, 150}, "XicXiccDecayRadDiffStra") ;

  auto h_df_xi_c_qa_xi_cc_ddca_stra = df_xi_c_qa.Histo1D({"df_xi_c_qa_xi_cc_ddca_stra", "xi_cc prong dca", 500, 0, 500}, "fXiccDaughterDCAStraTrack"); 
  auto h_df_xi_c_qa_xi_cc_ddist_pv_stra = df_xi_c_qa.Histo1D({"df_xi_c_qa_xi_cc_dist_pv_stra", "xi_cc decay dist", 3000, 0, 0.50}, "fXiccInvDecayLengthToPVStra"); 
  auto h_df_xi_c_qa_xi_cc_trad_stra = df_xi_c_qa.Histo1D({"df_xi_c_qa_xi_cc_trad_stra", "xi_cc trad", 2000, 0, 0.5}, "fXiccDecayRadiusStraTrack"); 
  
  auto h_df_xi_c_qa_xi_c_dca_xy_stra = df_xi_c_qa.Histo1D({"df_xi_c_qa_xi_c_dca_xy_stra", "xi_c dca xy stra", 1000, -500, 500}, "fXicDCAxyToPVStraTrack");  
  auto h_df_xi_c_qa_xi_c_dca_z_stra  = df_xi_c_qa.Histo1D({"df_xi_c_qa_xi_c_dca_z_stra", "xi_c dca z stra", 1000, -500, 500}, "fXicDCAzToPVStraTrack");  
  
  auto h_df_xi_c_qa_xi_cc_dca_xy_stra = df_xi_c_qa.Histo1D({"df_xi_c_qa_xi_cc_dca_xy_stra", "xi_cc dca xy stra", 1000, -500, 500}, "fXiccDCAxyToPVStraTrack");  
  auto h_df_xi_c_qa_xi_cc_dca_z_stra  = df_xi_c_qa.Histo1D({"df_xi_c_qa_xi_cc_dca_z_stra", "xi_cc dca z stra", 1000, -500, 500}, "fXiccDCAzToPVStraTrack");  

  //Pions (from the xicc)
 
  auto h_df_xi_c_qa_pi_dca_xy = df_xi_c_qa.Histo1D({"df_xi_c_qa_pi_dca_xy", "xi_c dca xy stra", 1000, -500, 500}, "fPicDCAxyToPVStraTrack");  
  auto h_df_xi_c_qa_pi_dca_z  = df_xi_c_qa.Histo1D({"df_xi_c_qa_pi_dca_z", "xi_c dca z stra", 1000, -500, 500}, "fPicDCAzToPVStraTrack");  
  
  auto h_df_xi_c_qa_xi_cc_pt = df_xi_c_qa.Histo1D({"df_xi_c_qa_xi_cc_pt", "xi_cc pt", 200, 0, 20}, "lPtXiCCStraTrack");  
  auto h_df_xi_c_qa_pi_pt = df_xi_c_qa.Histo1D({"df_xi_c_qa_pi_pt", "pi cc pt", 200, 0, 20}, "fPiCCPt");  

  //cut on Trad and check decay lengths 
  auto h_df_xi_c_qa_trad_xi_cc_ddist_pv_stra = df_xi_c_qa.Filter("fXiccDecayRadiusStraTrack > 0.004","corrStudy_fXiccDecayRadiusStraTrack").Histo1D({"df_xi_c_qa_trad_xi_cc_dist_pv_stra", "xi_cc decay dist", 3000, 0, 0.50}, "fXiccInvDecayLengthToPVStra"); 
  
  //cut on dl to pv 
  auto h_df_xi_c_qa_dl_pv_xi_cc_trad_stra = df_xi_c_qa.Filter("fXiccInvDecayLengthToPVStra > 0.004","corrStudy_fXiccInvDecayLengthToPVStra").Histo1D({"df_xi_c_qa_dl_pv_xi_cc_trad_stra", "xi_cc trad", 2000, 0, 0.4}, "fXiccDecayRadiusStraTrack");
  
  //some product magic .. 

  auto h_df_xi_c_qa_dcaxyProd   = df_xi_c_qa.Histo1D({"df_xi_c_qa_dca_xy_prod", "dca xy prod", 1000, -1e8, 1e8},"DCAxyProd"  ); 
  auto h_df_xi_c_qa_dcazProd    = df_xi_c_qa.Histo1D({"df_xi_c_qa_dca_z_prod", "dca z prod", 1000, -1e8, 1e8},"DCAzProd"   ); 
  auto h_df_xi_c_qa_dcaProdAdd  = df_xi_c_qa.Histo1D({"df_xi_c_qa_dca_prod_add", "dca prod add", 1000, 0, 1e15},"DCAProdAdd" ); 
  auto h_df_xi_c_qa_dcaProdMult = df_xi_c_qa.Histo1D({"df_xi_c_qa_dca_prod_mult", "dca prod mult", 1000, -1e9, 1e9},"DCAProdMult"); 
  auto h_df_xi_c_qa_dcaxySqSum   = df_xi_c_qa.Histo1D({"df_xi_c_qa_dca_xy_sqsum", "dca xy sqsum", 1000, 0, 1000000},"DCAxySqSum"  ); 
  auto h_df_xi_c_qa_dcazSqSum   = df_xi_c_qa.Histo1D({"df_xi_c_qa_dca_z_sqsum", "dca z sqsum", 1000, 0, 1000000},"DCAzSqSum"  ); 
  
  auto h_df_xi_c_qa_dcaxy_prod_xic_pi   = df_xi_c_qa.Histo1D({"df_xi_c_qa_dca_xy_prod_xic_pi", "dca xy prod", 1000, -1000, 1000},"DCAxyProdXicPi"  ); 
  auto h_df_xi_c_qa_dcaz_prod_xic_pi   = df_xi_c_qa.Histo1D({"df_xi_c_qa_dca_z_prod_xic_pi", "dca z prod", 1000, -1000, 1000},"DCAzProdXicPi"  ); 
  
  auto h_df_xi_c_xi_cc_mass_stra = df_xi_c.Filter("XicXiccDecayRadDiffStra > 0","XicXiccDecayRadDiffStra").Histo1D({"df_xi_c_xi_cc_mass_stra", "xi_cc inv mass", 700, 2.6, 4.6}, "fXiccMassStraTrack"); 
  
  //Select the Xi_cc
  auto df_xi_cc_im_c1 = df_xi_c
    .Filter("XicXiccDecayRadDiffStra > 0","c1_XicXiccDecayRadDiffStra")
    .Filter("fXicDaughterDCAStraTrack < 20","c1_fXicDaughterDCAStraTrack")
    .Filter("fXicDecayRadiusStraTrack > 0.004","c1_fXicDecayRadiusStraTrack")
    .Filter("fXicInvDecayLengthToPVStra > 0.002","c1_fXicInvDecayLengthToPVStra")
    .Filter("fXicInvDecayLengthToDVStra < 0.1","c1_fXicInvDecayLengthToDVStra")
    .Filter("fXiccDaughterDCAStraTrack < 20","c1_fXiccDaughterDCAStraTrack")
    .Filter("fXiccDecayRadiusStraTrack > 0.005","c1_fXiccDecayRadiusStraTrack")
    .Filter("fXiccInvDecayLengthToPVStra > 0.004","c1_fXiccInvDecayLengthToPVStra")
    .Filter("TMath::Abs(fXiDCAxyToPVStraTrack) > 5","c1_fXiDCAxyToPVStraTrack")
    .Filter("TMath::Abs(fXiDCAzToPVStraTrack) > 10","c1_fXiDCAzToPVStraTrack")
    .Filter("TMath::Abs(fXicDCAxyToPVStraTrack) > 15","c1_fXicDCAxyToPVStraTrack")
    .Filter("TMath::Abs(fXicDCAzToPVStraTrack) > 10","c1_fXicDCAzToPVStraTrack")
    .Filter("TMath::Abs(fXicPionDCAxyToPV1) > 10","c1_fXicPionDCAxyToPV1")
    .Filter("TMath::Abs(fXicPionDCAzToPV1) > 15","c1_fXicPionDCAzToPV1")
    .Filter("TMath::Abs(fXicPionDCAxyToPV2) > 10","c1_fXicPionDCAxyToPV2")
    .Filter("TMath::Abs(fXicPionDCAzToPV2) > 20","c1_fXicPionDCAzToPV2")
    .Filter("TMath::Abs(fPicDCAxyToPVTopo) > 20","c1_fPicDCAxyToPVTopo")
    ;
  
  //Fill some final histograms  
  auto h_df_xi_cc_im_xi_cc_mass_stra_c1 = df_xi_cc_im_c1.Histo1D({"df_xi_cc_im_xi_cc_mass_stra_c1", "xi_cc inv mass", 700, 2.6, 4.6}, "fXiccMassStraTrack"); 
  auto h_df_xi_cc_im_xi_cc_pt_c1 = df_xi_cc_im_c1.Filter(invMassXiccCut, {"fXiccMassStraTrack"}).Histo2D({"df_xi_cc_im_xi_cc_pt_vs_y_c1", "pt selected", 200, 0, 20, 30, -1.5, 1.5}, "lPtMCXiCC", "fXiccY"); 
  auto out_counter_c1 = df_xi_cc_im_c1.Filter(invMassXiccCut, {"fXiccMassStraTrack"}).Count(); 
  
  //Select the Xi_cc
  auto df_xi_cc_im_c2 = df_xi_c
    .Filter("XicXiccDecayRadDiffStra > 0","c2_XicXiccDecayRadDiffStra")
    .Filter("fXicDaughterDCAStraTrack < 20","c2_fXicDaughterDCAStraTrack")
    //.Filter("fXicDecayRadiusStraTrack > 0.004","c2_fXicDecayRadiusStraTrack")
    //.Filter("fXicInvDecayLengthToPVStra > 0.002","c2_fXicInvDecayLengthToPVStra")
    .Filter("fXicInvDecayLengthToDVStra < 0.1","c2_fXicInvDecayLengthToDVStra")
    .Filter("fXiccDaughterDCAStraTrack < 20","c2_fXiccDaughterDCAStraTrack")
    //.Filter("fXiccDecayRadiusStraTrack > 0.005","c2_fXiccDecayRadiusStraTrack")
    //.Filter("fXiccInvDecayLengthToPVStra > 0.004","c2_fXiccInvDecayLengthToPVStra")
    .Filter("TMath::Abs(fXiDCAxyToPVStraTrack) > 5","c2_fXiDCAxyToPVStraTrack")
    .Filter("TMath::Abs(fXiDCAzToPVStraTrack) > 10","c2_fXiDCAzToPVStraTrack")
    .Filter("TMath::Abs(fXicDCAxyToPVStraTrack) > 15","c2_fXicDCAxyToPVStraTrack")
    .Filter("TMath::Abs(fXicDCAzToPVStraTrack) > 10","c2_fXicDCAzToPVStraTrack")
    .Filter("TMath::Abs(fXicPionDCAxyToPV1) > 10","c2_fXicPionDCAxyToPV1")
    .Filter("TMath::Abs(fXicPionDCAzToPV1) > 15","c2_fXicPionDCAzToPV1")
    .Filter("TMath::Abs(fXicPionDCAxyToPV2) > 10","c2_fXicPionDCAxyToPV2")
    .Filter("TMath::Abs(fXicPionDCAzToPV2) > 20","c2_fXicPionDCAzToPV2")
    .Filter("TMath::Abs(fPicDCAxyToPVTopo) > 20","c2_fPicDCAxyToPVTopo")
    ;
  
  //Fill some final histograms  
  auto h_df_xi_cc_im_xi_cc_mass_stra_c2 = df_xi_cc_im_c2.Histo1D({"df_xi_cc_im_xi_cc_mass_stra_c2", "xi_cc inv mass", 700, 2.6, 4.6}, "fXiccMassStraTrack"); 
  auto h_df_xi_cc_im_xi_cc_pt_c2 = df_xi_cc_im_c2.Filter(invMassXiccCut, {"fXiccMassStraTrack"}).Histo2D({"df_xi_cc_im_xi_cc_pt_vs_y_c2", "pt selected", 200, 0, 20, 30, -1.5, 1.5}, "lPtMCXiCC", "fXiccY"); 
  
  auto out_counter_c2 = df_xi_cc_im_c2.Filter(invMassXiccCut, {"fXiccMassStraTrack"}).Count(); 
  
  //Select the Xi_cc
  auto df_xi_cc_im_c3 = df_xi_c
    .Filter("XicXiccDecayRadDiffStra > 0","c3_XicXiccDecayRadDiffStra")
    .Filter("fXicDaughterDCAStraTrack < 20","c3_fXicDaughterDCAStraTrack")
    .Filter("fXicDecayRadiusStraTrack > 0.004","c3_fXicDecayRadiusStraTrack")
    .Filter("fXicInvDecayLengthToPVStra > 0.002","c3_fXicInvDecayLengthToPVStra")
    .Filter("fXicInvDecayLengthToDVStra < 0.1","c3_fXicInvDecayLengthToDVStra")
    .Filter("fXiccDaughterDCAStraTrack < 20","c3_fXiccDaughterDCAStraTrack")
    //.Filter("fXiccDecayRadiusStraTrack > 0.005","c3_fXiccDecayRadiusStraTrack")
    //.Filter("fXiccInvDecayLengthToPVStra > 0.004","c3_fXiccInvDecayLengthToPVStra")
    .Filter("TMath::Abs(fXiDCAxyToPVStraTrack) > 5","c3_fXiDCAxyToPVStraTrack")
    .Filter("TMath::Abs(fXiDCAzToPVStraTrack) > 10","c3_fXiDCAzToPVStraTrack")
    .Filter("TMath::Abs(fXicDCAxyToPVStraTrack) > 15","c3_fXicDCAxyToPVStraTrack")
    .Filter("TMath::Abs(fXicDCAzToPVStraTrack) > 10","c3_fXicDCAzToPVStraTrack")
    .Filter("TMath::Abs(fXicPionDCAxyToPV1) > 10","c3_fXicPionDCAxyToPV1")
    .Filter("TMath::Abs(fXicPionDCAzToPV1) > 15","c3_fXicPionDCAzToPV1")
    .Filter("TMath::Abs(fXicPionDCAxyToPV2) > 10","c3_fXicPionDCAxyToPV2")
    .Filter("TMath::Abs(fXicPionDCAzToPV2) > 20","c3_fXicPionDCAzToPV2")
    .Filter("TMath::Abs(fPicDCAxyToPVTopo) > 20","c3_fPicDCAxyToPVTopo")
    ;
  
  //Fill some final histograms  
  auto h_df_xi_cc_im_xi_cc_mass_stra_c3 = df_xi_cc_im_c3.Histo1D({"df_xi_cc_im_xi_cc_mass_stra_c3", "xi_cc inv mass", 700, 2.6, 4.6}, "fXiccMassStraTrack"); 
  auto h_df_xi_cc_im_xi_cc_pt_c3 = df_xi_cc_im_c3.Filter(invMassXiccCut, {"fXiccMassStraTrack"}).Histo2D({"df_xi_cc_im_xi_cc_pt_vs_y_c3", "pt selected", 200, 0, 20, 30, -1.5, 1.5}, "lPtMCXiCC", "fXiccY"); 
  
  auto out_counter_c3 = df_xi_cc_im_c3.Filter(invMassXiccCut, {"fXiccMassStraTrack"}).Count(); 
  
  auto df_xi_cc_im_c4 = df_xi_c
    .Filter("XicXiccDecayRadDiffStra > 0","c4_XicXiccDecayRadDiffStra")
    .Filter("XicXiccDecayRadDiffStra > 0","c4_XicXiccDecayRadDiffStra")
    .Filter("XicXiccDecayRadDiffStra > 0","c4_XicXiccDecayRadDiffStra")
    .Filter("fXicDaughterDCAStraTrack < 20","c4_fXicDaughterDCAStraTrack")
    .Filter("fXicDecayRadiusStraTrack > 0.004","c4_fXicDecayRadiusStraTrack")
    .Filter("fXicInvDecayLengthToPVStra > 0.002","c4_fXicInvDecayLengthToPVStra")
    .Filter("fXicInvDecayLengthToDVStra < 0.1","c4_fXicInvDecayLengthToDVStra")
    .Filter("fXiccDaughterDCAStraTrack < 20","c4_fXiccDaughterDCAStraTrack")
    .Filter("fXiccDecayRadiusStraTrack > 0.005","c4_fXiccDecayRadiusStraTrack")
    .Filter("fXiccInvDecayLengthToPVStra > 0.004","c4_fXiccInvDecayLengthToPVStra")
    .Filter("TMath::Abs(fXiDCAxyToPVStraTrack) > 5","c4_fXiDCAxyToPVStraTrack")
    .Filter("TMath::Abs(fXiDCAzToPVStraTrack) > 10","c4_fXiDCAzToPVStraTrack")
    .Filter("TMath::Abs(fXicDCAxyToPVStraTrack) > 20","c4_fXicDCAxyToPVStraTrack")
    .Filter("TMath::Abs(fXicDCAzToPVStraTrack) > 20","c4_fXicDCAzToPVStraTrack")
    .Filter("TMath::Abs(fXicPionDCAxyToPV1) > 10","c4_fXicPionDCAxyToPV1")
    .Filter("TMath::Abs(fXicPionDCAzToPV1) > 15","c4_fXicPionDCAzToPV1")
    .Filter("TMath::Abs(fXicPionDCAxyToPV2) > 10","c4_fXicPionDCAxyToPV2")
    .Filter("TMath::Abs(fXicPionDCAzToPV2) > 20","c4_fXicPionDCAzToPV2")
    .Filter("TMath::Abs(fPicDCAxyToPVTopo) > 20","c4_fPicDCAxyToPVTopo")
    .Filter("TMath::Abs(fPicDCAzToPVTopo) > 20","c4_fPicDCAzToPVTopo")
    .Filter("TMath::Abs(fXiccDCAxyToPVStraTrack) < 50","c4_fXiccDCAxyToPVStraTrack")
    .Filter("TMath::Abs(fXiccDCAzToPVStraTrack) < 50","c4_fXiccDCAzToPVStraTrack")    
    ;
  
  //Fill some final histograms  
  auto h_df_xi_cc_im_xi_cc_mass_stra_c4 = df_xi_cc_im_c4.Histo1D({"df_xi_cc_im_xi_cc_mass_stra_c4", "xi_cc inv mass", 700, 2.6, 4.6}, "fXiccMassStraTrack"); 
  auto h_df_xi_cc_im_xi_cc_pt_c4 = df_xi_cc_im_c4.Filter(invMassXiccCut, {"fXiccMassStraTrack"}).Histo2D({"df_xi_cc_im_xi_cc_pt_vs_y_c4", "pt selected", 200, 0, 20, 30, -1.5, 1.5}, "lPtMCXiCC", "fXiccY"); 
  
  auto out_counter_c4 = df_xi_cc_im_c4.Filter(invMassXiccCut, {"fXiccMassStraTrack"}).Count();
 
  auto df_xi_cc_im_c4_2Hit = df_xi_c
    .Filter("fXiHitsAdded > 1")
    .Filter("XicXiccDecayRadDiffStra > 0","c4_2Hit_XicXiccDecayRadDiffStra")
    .Filter("XicXiccDecayRadDiffStra > 0","c4_2Hit_XicXiccDecayRadDiffStra")
    .Filter("XicXiccDecayRadDiffStra > 0","c4_2Hit_XicXiccDecayRadDiffStra")
    .Filter("fXicDaughterDCAStraTrack < 20","c4_2Hit_fXicDaughterDCAStraTrack")
    .Filter("fXicDecayRadiusStraTrack > 0.004","c4_2Hit_fXicDecayRadiusStraTrack")
    .Filter("fXicInvDecayLengthToPVStra > 0.002","c4_2Hit_fXicInvDecayLengthToPVStra")
    .Filter("fXicInvDecayLengthToDVStra < 0.1","c4_2Hit_fXicInvDecayLengthToDVStra")
    .Filter("fXiccDaughterDCAStraTrack < 20","c4_2Hit_fXiccDaughterDCAStraTrack")
    .Filter("fXiccDecayRadiusStraTrack > 0.005","c4_2Hit_fXiccDecayRadiusStraTrack")
    .Filter("fXiccInvDecayLengthToPVStra > 0.004","c4_2Hit_fXiccInvDecayLengthToPVStra")
    .Filter("TMath::Abs(fXiDCAxyToPVStraTrack) > 5","c4_2Hit_fXiDCAxyToPVStraTrack")
    .Filter("TMath::Abs(fXiDCAzToPVStraTrack) > 10","c4_2Hit_fXiDCAzToPVStraTrack")
    .Filter("TMath::Abs(fXicDCAxyToPVStraTrack) > 20","c4_2Hit_fXicDCAxyToPVStraTrack")
    .Filter("TMath::Abs(fXicDCAzToPVStraTrack) > 20","c4_2Hit_fXicDCAzToPVStraTrack")
    .Filter("TMath::Abs(fXicPionDCAxyToPV1) > 10","c4_2Hit_fXicPionDCAxyToPV1")
    .Filter("TMath::Abs(fXicPionDCAzToPV1) > 15","c4_2Hit_fXicPionDCAzToPV1")
    .Filter("TMath::Abs(fXicPionDCAxyToPV2) > 10","c4_2Hit_fXicPionDCAxyToPV2")
    .Filter("TMath::Abs(fXicPionDCAzToPV2) > 20","c4_2Hit_fXicPionDCAzToPV2")
    .Filter("TMath::Abs(fPicDCAxyToPVTopo) > 20","c4_2Hit_fPicDCAxyToPVTopo")
    .Filter("TMath::Abs(fPicDCAzToPVTopo) > 20","c4_2Hit_fPicDCAzToPVTopo")
    .Filter("TMath::Abs(fXiccDCAxyToPVStraTrack) < 50","c4_2Hit_fXiccDCAxyToPVStraTrack")
    .Filter("TMath::Abs(fXiccDCAzToPVStraTrack) < 50","c4_2Hit_fXiccDCAzToPVStraTrack")    
    ;
  
  //Fill some final histograms  
  auto h_df_xi_cc_im_xi_cc_mass_stra_c4_2Hit = df_xi_cc_im_c4_2Hit.Histo1D({"df_xi_cc_im_xi_cc_mass_stra_c4_2Hit", "xi_cc inv mass", 700, 2.6, 4.6}, "fXiccMassStraTrack"); 
  auto h_df_xi_cc_im_xi_cc_pt_c4_2Hit = df_xi_cc_im_c4_2Hit.Filter(invMassXiccCut, {"fXiccMassStraTrack"}).Histo2D({"df_xi_cc_im_xi_cc_pt_vs_y_c4_2Hit", "pt selected", 200, 0, 20, 30, -1.5, 1.5}, "lPtMCXiCC", "fXiccY"); 
  auto out_counter_c4_2Hit = df_xi_cc_im_c4_2Hit.Filter(invMassXiccCut, {"fXiccMassStraTrack"}).Count();
  
  auto df_xi_cc_im_c4_3Hit = df_xi_c
    .Filter("fXiHitsAdded > 2")
    .Filter("XicXiccDecayRadDiffStra > 0","c4_3Hit_XicXiccDecayRadDiffStra")
    .Filter("XicXiccDecayRadDiffStra > 0","c4_3Hit_XicXiccDecayRadDiffStra")
    .Filter("XicXiccDecayRadDiffStra > 0","c4_3Hit_XicXiccDecayRadDiffStra")
    .Filter("fXicDaughterDCAStraTrack < 20","c4_3Hit_fXicDaughterDCAStraTrack")
    .Filter("fXicDecayRadiusStraTrack > 0.004","c4_3Hit_fXicDecayRadiusStraTrack")
    .Filter("fXicInvDecayLengthToPVStra > 0.002","c4_3Hit_fXicInvDecayLengthToPVStra")
    .Filter("fXicInvDecayLengthToDVStra < 0.1","c4_3Hit_fXicInvDecayLengthToDVStra")
    .Filter("fXiccDaughterDCAStraTrack < 20","c4_3Hit_fXiccDaughterDCAStraTrack")
    .Filter("fXiccDecayRadiusStraTrack > 0.005","c4_3Hit_fXiccDecayRadiusStraTrack")
    .Filter("fXiccInvDecayLengthToPVStra > 0.004","c4_3Hit_fXiccInvDecayLengthToPVStra")
    .Filter("TMath::Abs(fXiDCAxyToPVStraTrack) > 5","c4_3Hit_fXiDCAxyToPVStraTrack")
    .Filter("TMath::Abs(fXiDCAzToPVStraTrack) > 10","c4_3Hit_fXiDCAzToPVStraTrack")
    .Filter("TMath::Abs(fXicDCAxyToPVStraTrack) > 20","c4_3Hit_fXicDCAxyToPVStraTrack")
    .Filter("TMath::Abs(fXicDCAzToPVStraTrack) > 20","c4_3Hit_fXicDCAzToPVStraTrack")
    .Filter("TMath::Abs(fXicPionDCAxyToPV1) > 10","c4_3Hit_fXicPionDCAxyToPV1")
    .Filter("TMath::Abs(fXicPionDCAzToPV1) > 15","c4_3Hit_fXicPionDCAzToPV1")
    .Filter("TMath::Abs(fXicPionDCAxyToPV2) > 10","c4_3Hit_fXicPionDCAxyToPV2")
    .Filter("TMath::Abs(fXicPionDCAzToPV2) > 20","c4_3Hit_fXicPionDCAzToPV2")
    .Filter("TMath::Abs(fPicDCAxyToPVTopo) > 20","c4_3Hit_fPicDCAxyToPVTopo")
    .Filter("TMath::Abs(fPicDCAzToPVTopo) > 20","c4_3Hit_fPicDCAzToPVTopo")
    .Filter("TMath::Abs(fXiccDCAxyToPVStraTrack) < 50","c4_3Hit_fXiccDCAxyToPVStraTrack")
    .Filter("TMath::Abs(fXiccDCAzToPVStraTrack) < 50","c4_3Hit_fXiccDCAzToPVStraTrack")    
    ;
  //Fill some final histograms  
  auto h_df_xi_cc_im_xi_cc_mass_stra_c4_3Hit = df_xi_cc_im_c4_3Hit.Histo1D({"df_xi_cc_im_xi_cc_mass_stra_c4_3Hit", "xi_cc inv mass", 700, 2.6, 4.6}, "fXiccMassStraTrack"); 
  auto h_df_xi_cc_im_xi_cc_pt_c4_3Hit = df_xi_cc_im_c4_3Hit.Filter(invMassXiccCut, {"fXiccMassStraTrack"}).Histo2D({"df_xi_cc_im_xi_cc_pt_vs_y_c4_3Hit", "pt selected", 200, 0, 20, 30, -1.5, 1.5}, "lPtMCXiCC", "fXiccY"); 
  auto out_counter_c4_3Hit = df_xi_cc_im_c4_3Hit.Filter(invMassXiccCut, {"fXiccMassStraTrack"}).Count();

  //Select the Xi_cc
  auto df_xi_cc_im_c5 = df_xi_c
    .Filter("XicXiccDecayRadDiffStra > 0","c5_XicXiccDecayRadDiffStra")
    .Filter("fXicDaughterDCAStraTrack < 20","c5_fXicDaughterDCAStraTrack")
    //.Filter("fXicDecayRadiusStraTrack > 0.004","c5_fXicDecayRadiusStraTrack")
    //.Filter("fXicInvDecayLengthToPVStra > 0.002","c5_fXicInvDecayLengthToPVStra")
    .Filter("fXicInvDecayLengthToDVStra < 0.1","c5_fXicInvDecayLengthToDVStra")
    .Filter("fXiccDaughterDCAStraTrack < 20","c5_fXiccDaughterDCAStraTrack")
    //.Filter("fXiccDecayRadiusStraTrack > 0.005","c5_fXiccDecayRadiusStraTrack")
    //.Filter("fXiccInvDecayLengthToPVStra > 0.004","c5_fXiccInvDecayLengthToPVStra")
    .Filter("TMath::Abs(fXiDCAxyToPVStraTrack) > 5","c5_fXiDCAxyToPVStraTrack")
    .Filter("TMath::Abs(fXiDCAzToPVStraTrack) > 10","c5_fXiDCAzToPVStraTrack")
    //.Filter("TMath::Abs(fXicDCAxyToPVStraTrack) > 15","c5_fXicDCAxyToPVStraTrack")
    //.Filter("TMath::Abs(fXicDCAzToPVStraTrack) > 10","c5_fXicDCAzToPVStraTrack")
    .Filter("TMath::Abs(fXicPionDCAxyToPV1) > 10","c5_fXicPionDCAxyToPV1")
    .Filter("TMath::Abs(fXicPionDCAzToPV1) > 15","c5_fXicPionDCAzToPV1")
    .Filter("TMath::Abs(fXicPionDCAxyToPV2) > 10","c5_fXicPionDCAxyToPV2")
    .Filter("TMath::Abs(fXicPionDCAzToPV2) > 20","c5_fXicPionDCAzToPV2")
    .Filter("TMath::Abs(fPicDCAxyToPVTopo) > 20","c5_fPicDCAxyToPVTopo")
    ;
  
  //Fill some final histograms  
  auto h_df_xi_cc_im_xi_cc_mass_stra_c5 = df_xi_cc_im_c5.Histo1D({"df_xi_cc_im_xi_cc_mass_stra_c5", "xi_cc inv mass", 700, 2.6, 4.6}, "fXiccMassStraTrack"); 
  auto h_df_xi_cc_im_xi_cc_pt_c5 = df_xi_cc_im_c5.Filter(invMassXiccCut, {"fXiccMassStraTrack"}).Histo2D({"df_xi_cc_im_xi_cc_pt_vs_y_c5", "pt selected", 200, 0, 20, 30, -1.5, 1.5}, "lPtMCXiCC", "fXiccY"); 
  
  auto out_counter_c5 = df_xi_cc_im_c5.Filter(invMassXiccCut, {"fXiccMassStraTrack"}).Count(); 

  //Select the Xi_cc
  auto df_xi_cc_im_c6 = df_xi_c
    .Filter("XicXiccDecayRadDiffStra > 0","c6_XicXiccDecayRadDiffStra")
    .Filter("fXicDaughterDCAStraTrack < 20","c6_fXicDaughterDCAStraTrack")
    .Filter("fXicDecayRadiusStraTrack > 0.004","c6_fXicDecayRadiusStraTrack")
    .Filter("fXicInvDecayLengthToPVStra > 0.002","c6_fXicInvDecayLengthToPVStra")
    .Filter("fXicInvDecayLengthToDVStra < 0.1","c6_fXicInvDecayLengthToDVStra")
    .Filter("fXiccDaughterDCAStraTrack < 20","c6_fXiccDaughterDCAStraTrack")
    //.Filter("fXiccDecayRadiusStraTrack > 0.005","c6_fXiccDecayRadiusStraTrack")
    //.Filter("fXiccInvDecayLengthToPVStra > 0.004","c6_fXiccInvDecayLengthToPVStra")
    .Filter("TMath::Abs(fXiDCAxyToPVStraTrack) > 5","c6_fXiDCAxyToPVStraTrack")
    .Filter("TMath::Abs(fXiDCAzToPVStraTrack) > 10","c6_fXiDCAzToPVStraTrack")
    //.Filter("TMath::Abs(fXicDCAxyToPVStraTrack) > 15","c6_fXicDCAxyToPVStraTrack")
    //.Filter("TMath::Abs(fXicDCAzToPVStraTrack) > 10","c6_fXicDCAzToPVStraTrack")
    .Filter("TMath::Abs(fXicPionDCAxyToPV1) > 10","c6_fXicPionDCAxyToPV1")
    .Filter("TMath::Abs(fXicPionDCAzToPV1) > 15","c6_fXicPionDCAzToPV1")
    .Filter("TMath::Abs(fXicPionDCAxyToPV2) > 10","c6_fXicPionDCAxyToPV2")
    .Filter("TMath::Abs(fXicPionDCAzToPV2) > 20","c6_fXicPionDCAzToPV2")
    .Filter("TMath::Abs(fPicDCAxyToPVTopo) > 20","c6_fPicDCAxyToPVTopo")
    ;
  
  //Fill some final histograms  
  auto h_df_xi_cc_im_xi_cc_mass_stra_c6 = df_xi_cc_im_c6.Histo1D({"df_xi_cc_im_xi_cc_mass_stra_c6", "xi_cc inv mass", 700, 2.6, 4.6}, "fXiccMassStraTrack"); 
  auto h_df_xi_cc_im_xi_cc_pt_c6 = df_xi_cc_im_c6.Filter(invMassXiccCut, {"fXiccMassStraTrack"}).Histo2D({"df_xi_cc_im_xi_cc_pt_vs_y_c6", "pt selected", 200, 0, 20, 30, -1.5, 1.5}, "lPtMCXiCC", "fXiccY"); 

  auto out_counter_c6 = df_xi_cc_im_c6.Filter(invMassXiccCut, {"fXiccMassStraTrack"}).Count(); 
  
  
  //Select the Xi_cc
  auto df_xi_cc_im_c7 = df_xi_c
    .Filter("XicXiccDecayRadDiffStra > 0","c7_XicXiccDecayRadDiffStra")
    .Filter("fXicDaughterDCAStraTrack < 20","c7_fXicDaughterDCAStraTrack")
    .Filter("fXicDecayRadiusStraTrack > 0.004","c7_fXicDecayRadiusStraTrack")
    .Filter("fXicInvDecayLengthToPVStra > 0.002","c7_fXicInvDecayLengthToPVStra")
    .Filter("fXicInvDecayLengthToDVStra < 0.1","c7_fXicInvDecayLengthToDVStra")
    .Filter("fXiccDaughterDCAStraTrack < 20","c7_fXiccDaughterDCAStraTrack")
    .Filter("fXiccDecayRadiusStraTrack > 0.005","c7_fXiccDecayRadiusStraTrack")
    .Filter("fXiccInvDecayLengthToPVStra > 0.004","c7_fXiccInvDecayLengthToPVStra")
    .Filter("TMath::Abs(fXiDCAxyToPVStraTrack) > 5","c7_fXiDCAxyToPVStraTrack")
    .Filter("TMath::Abs(fXiDCAzToPVStraTrack) > 10","c7_fXiDCAzToPVStraTrack")
    //.Filter("TMath::Abs(fXicDCAxyToPVStraTrack) > 15","c7_fXicDCAxyToPVStraTrack")
    //.Filter("TMath::Abs(fXicDCAzToPVStraTrack) > 10","c7_fXicDCAzToPVStraTrack")
    .Filter("TMath::Abs(fXicPionDCAxyToPV1) > 10","c7_fXicPionDCAxyToPV1")
    .Filter("TMath::Abs(fXicPionDCAzToPV1) > 15","c7_fXicPionDCAzToPV1")
    .Filter("TMath::Abs(fXicPionDCAxyToPV2) > 10","c7_fXicPionDCAxyToPV2")
    .Filter("TMath::Abs(fXicPionDCAzToPV2) > 20","c7_fXicPionDCAzToPV2")
    .Filter("TMath::Abs(fPicDCAxyToPVTopo) > 20","c7_fPicDCAxyToPVTopo")
    ;
  
  //Fill some final histograms  
  auto h_df_xi_cc_im_xi_cc_mass_stra_c7 = df_xi_cc_im_c7.Histo1D({"df_xi_cc_im_xi_cc_mass_stra_c7", "xi_cc inv mass", 700, 2.6, 4.6}, "fXiccMassStraTrack"); 
  auto h_df_xi_cc_im_xi_cc_pt_c7 = df_xi_cc_im_c7.Filter(invMassXiccCut, {"fXiccMassStraTrack"}).Histo2D({"df_xi_cc_im_xi_cc_pt_vs_y_c7", "pt selected", 200, 0, 20, 30, -1.5, 1.5}, "lPtMCXiCC", "fXiccY"); 
  auto out_counter_c7 = df_xi_cc_im_c7.Filter(invMassXiccCut, {"fXiccMassStraTrack"}).Count(); 
  
  //Select the Xi_cc
  auto df_xi_cc_im_c8 = df_xi_c
    .Filter("XicXiccDecayRadDiffStra > 0","c8_XicXiccDecayRadDiffStra")
    .Filter("fXicDaughterDCAStraTrack < 20","c8_fXicDaughterDCAStraTrack")
    .Filter("fXicDecayRadiusStraTrack > 0.004","c8_fXicDecayRadiusStraTrack")
    .Filter("fXicInvDecayLengthToPVStra > 0.002","c8_fXicInvDecayLengthToPVStra")
    .Filter("fXicInvDecayLengthToDVStra < 0.1","c8_fXicInvDecayLengthToDVStra")
    .Filter("fXiccDaughterDCAStraTrack < 20","c8_fXiccDaughterDCAStraTrack")
    .Filter("fXiccDecayRadiusStraTrack > 0.005","c8_fXiccDecayRadiusStraTrack")
    .Filter("fXiccInvDecayLengthToPVStra > 0.004","c8_fXiccInvDecayLengthToPVStra")
    .Filter("TMath::Abs(fXiDCAxyToPVStraTrack) > 5","c8_fXiDCAxyToPVStraTrack")
    .Filter("TMath::Abs(fXiDCAzToPVStraTrack) > 10","c8_fXiDCAzToPVStraTrack")
    .Filter("TMath::Abs(fXicDCAxyToPVStraTrack) > 15","c8_fXicDCAxyToPVStraTrack")
    .Filter("TMath::Abs(fXicDCAzToPVStraTrack) > 10","c8_fXicDCAzToPVStraTrack")
    .Filter("TMath::Abs(fXicPionDCAxyToPV1) > 10","c8_fXicPionDCAxyToPV1")
    .Filter("TMath::Abs(fXicPionDCAzToPV1) > 15","c8_fXicPionDCAzToPV1")
    .Filter("TMath::Abs(fXicPionDCAxyToPV2) > 10","c8_fXicPionDCAxyToPV2")
    .Filter("TMath::Abs(fXicPionDCAzToPV2) > 20","c8_fXicPionDCAzToPV2")
    .Filter("TMath::Abs(fPicDCAxyToPVTopo) > 20","c8_fPicDCAxyToPVTopo")
    .Filter("fPiC1Pt > 0.35", "c8_fPiC1Pt")
    .Filter("fPiC2Pt > 0.35", "c8_fPiC2Pt")
    .Filter("fPiCCPt > 0.35", "c8_fPiCCPt")
    ;
  
  //Fill some final histograms  
  auto h_df_xi_cc_im_xi_cc_mass_stra_c8 = df_xi_cc_im_c8.Histo1D({"df_xi_cc_im_xi_cc_mass_stra_c8", "xi_cc inv mass", 700, 2.6, 4.6}, "fXiccMassStraTrack"); 
  auto h_df_xi_cc_im_xi_cc_pt_c8 = df_xi_cc_im_c8.Filter(invMassXiccCut, {"fXiccMassStraTrack"}).Histo2D({"df_xi_cc_im_xi_cc_pt_vs_y_c8", "pt selected", 200, 0, 20, 30, -1.5, 1.5}, "lPtMCXiCC", "fXiccY"); 
  auto out_counter_c8 = df_xi_cc_im_c8.Filter(invMassXiccCut, {"fXiccMassStraTrack"}).Count(); 
  
  auto df_xi_cc_im_c9 = df_xi_c
    .Filter("XicXiccDecayRadDiffStra > 0","c9_XicXiccDecayRadDiffStra")
    .Filter("XicXiccDecayRadDiffStra > 0","c9_XicXiccDecayRadDiffStra")
    .Filter("XicXiccDecayRadDiffStra > 0","c9_XicXiccDecayRadDiffStra")
    .Filter("fXicDaughterDCAStraTrack < 20","c9_fXicDaughterDCAStraTrack")
    .Filter("fXicDecayRadiusStraTrack > 0.004","c9_fXicDecayRadiusStraTrack")
    .Filter("fXicInvDecayLengthToPVStra > 0.002","c9_fXicInvDecayLengthToPVStra")
    .Filter("fXicInvDecayLengthToDVStra < 0.1","c9_fXicInvDecayLengthToDVStra")
    .Filter("fXiccDaughterDCAStraTrack < 20","c9_fXiccDaughterDCAStraTrack")
    .Filter("fXiccDecayRadiusStraTrack > 0.005","c9_fXiccDecayRadiusStraTrack")
    .Filter("fXiccInvDecayLengthToPVStra > 0.004","c9_fXiccInvDecayLengthToPVStra")
    .Filter("TMath::Abs(fXiDCAxyToPVStraTrack) > 5","c9_fXiDCAxyToPVStraTrack")
    .Filter("TMath::Abs(fXiDCAzToPVStraTrack) > 10","c9_fXiDCAzToPVStraTrack")
    .Filter("TMath::Abs(fXicDCAxyToPVStraTrack) > 20","c9_fXicDCAxyToPVStraTrack")
    .Filter("TMath::Abs(fXicDCAzToPVStraTrack) > 20","c9_fXicDCAzToPVStraTrack")
    .Filter("TMath::Abs(fXicPionDCAxyToPV1) > 10","c9_fXicPionDCAxyToPV1")
    .Filter("TMath::Abs(fXicPionDCAzToPV1) > 15","c9_fXicPionDCAzToPV1")
    .Filter("TMath::Abs(fXicPionDCAxyToPV2) > 10","c9_fXicPionDCAxyToPV2")
    .Filter("TMath::Abs(fXicPionDCAzToPV2) > 20","c9_fXicPionDCAzToPV2")
    .Filter("TMath::Abs(fPicDCAxyToPVTopo) > 20","c9_fPicDCAxyToPVTopo")
    .Filter("TMath::Abs(fPicDCAzToPVTopo) > 20","c9_fPicDCAzToPVTopo")
    .Filter("TMath::Abs(fXiccDCAxyToPVStraTrack) < 50","c9_fXiccDCAxyToPVStraTrack")
    .Filter("TMath::Abs(fXiccDCAzToPVStraTrack) < 50","c9_fXiccDCAzToPVStraTrack")    
    .Filter("fPiC1Pt > 0.35", "c9_fPiC1Pt")
    .Filter("fPiC2Pt > 0.35", "c9_fPiC2Pt")
    .Filter("fPiCCPt > 0.35", "c9_fPiCCPt")
    ;
  
  //Fill some final histograms  
  auto h_df_xi_cc_im_xi_cc_mass_stra_c9 = df_xi_cc_im_c9.Histo1D({"df_xi_cc_im_xi_cc_mass_stra_c9", "xi_cc inv mass", 700, 2.6, 4.6}, "fXiccMassStraTrack"); 
  auto h_df_xi_cc_im_xi_cc_pt_c9 = df_xi_cc_im_c9.Filter(invMassXiccCut, {"fXiccMassStraTrack"}).Histo2D({"df_xi_cc_im_xi_cc_pt_vs_y_c9", "pt selected", 200, 0, 20, 30, -1.5, 1.5}, "lPtMCXiCC", "fXiccY"); 
  auto out_counter_c9 = df_xi_cc_im_c9.Filter(invMassXiccCut, {"fXiccMassStraTrack"}).Count();


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
  
  reduction = *out_counter_c2/inputCounter;
  std::cout << reduction << std::endl; 
  cutCounter->SetBinContent(2, reduction);
  
  reduction = *out_counter_c3/inputCounter;
  std::cout << reduction << std::endl; 
  cutCounter->SetBinContent(3, reduction);
  
  reduction = *out_counter_c4/inputCounter;
  std::cout << reduction << std::endl; 
  cutCounter->SetBinContent(4, reduction);
  
  reduction = *out_counter_c4_2Hit/inputCounter;
  std::cout << reduction << std::endl; 
  cutCounter->SetBinContent(5, reduction);
  
  reduction = *out_counter_c4_3Hit/inputCounter;
  std::cout << reduction << std::endl; 
  cutCounter->SetBinContent(6, reduction);
  
  reduction = *out_counter_c5/inputCounter;
  std::cout << reduction << std::endl; 
  cutCounter->SetBinContent(7, reduction);
  
  reduction = *out_counter_c6/inputCounter;
  std::cout << reduction << std::endl; 
  cutCounter->SetBinContent(8, reduction);
  
  reduction = *out_counter_c7/inputCounter;
  std::cout << reduction << std::endl; 
  cutCounter->SetBinContent(9, reduction);
  
  reduction = *out_counter_c8/inputCounter;
  std::cout << reduction << std::endl; 
  cutCounter->SetBinContent(10, reduction);
  
  reduction = *out_counter_c9/inputCounter;
  std::cout << reduction << std::endl; 
  cutCounter->SetBinContent(10, reduction);
  
  
  TString outName = TString::Format("outxiccSelector_%s.root",outAddon )  ; 
  TFile* out = TFile::Open(outName.Data(), "recreate");
  
  out->cd(); 
  //this is my start
  HarryPlotter::CheckAndStore(out, h_df_in_im_xi_cc_mass_stra);
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
  HarryPlotter::CheckAndStore(out, h_df_lmb_xi_cc_mass_stra); 

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
  
  HarryPlotter::CheckAndStore(out, h_df_xi_xi_cc_mass_stra);

  HarryPlotter::CheckAndStore(out, h_df_lmb_qa_trad_diff_lmb_xi); 
  
  HarryPlotter::CheckAndStore(out, h_df_xi_sel_xi_mass);
  
  //to xi_c
  HarryPlotter::CheckAndStore(out, h_df_xi_im_xi_mass);

  HarryPlotter::CheckAndStore(out, h_df_xi_qa_xi_c_mass_topo); 
  HarryPlotter::CheckAndStore(out, h_df_xi_qa_trad_diff_xi_xi_c_topo);
 
  HarryPlotter::CheckAndStore(out, h_df_xi_qa_xi_c_ddca_topo); 
  HarryPlotter::CheckAndStore(out, h_df_xi_qa_xi_c_ddist_pv_topo); 
  HarryPlotter::CheckAndStore(out, h_df_xi_qa_xi_c_ddist_dv_topo);
  HarryPlotter::CheckAndStore(out, h_df_xi_qa_xi_c_trad_topo); 
    
  HarryPlotter::CheckAndStore(out, h_df_xi_qa_xi_dca_xy_topo);
  HarryPlotter::CheckAndStore(out, h_df_xi_qa_xi_dca_z_topo);
   
  HarryPlotter::CheckAndStore(out, h_df_xi_qa_xi_c_mass_stra); 
  HarryPlotter::CheckAndStore(out, h_df_xi_qa_xi_c_pt); 
  
  HarryPlotter::CheckAndStore(out, h_df_xi_qa_trad_diff_xi_xi_c_stra);
  HarryPlotter::CheckAndStore(out, h_df_xi_qa_xi_c_ddca_stra); 
  HarryPlotter::CheckAndStore(out, h_df_xi_qa_xi_c_ddist_pv_stra); 
  HarryPlotter::CheckAndStore(out, h_df_xi_qa_xi_c_ddist_dv_stra);
  HarryPlotter::CheckAndStore(out, h_df_xi_qa_xi_c_trad_stra); 

  HarryPlotter::CheckAndStore(out, h_df_xi_qa_xi_dca_xy_stra);
  HarryPlotter::CheckAndStore(out, h_df_xi_qa_xi_dca_z_stra);

  HarryPlotter::CheckAndStore(out, h_df_xi_qa_pi_one_pt);
  HarryPlotter::CheckAndStore(out, h_df_xi_qa_pi_two_pt);
    
  HarryPlotter::CheckAndStore(out, h_df_xi_qa_pi_one_dca_xy);
  HarryPlotter::CheckAndStore(out, h_df_xi_qa_pi_one_dca_z);
  HarryPlotter::CheckAndStore(out, h_df_xi_qa_pi_two_dca_xy);
  HarryPlotter::CheckAndStore(out, h_df_xi_qa_pi_two_dca_z);

  HarryPlotter::CheckAndStore(out, h_df_xi_qa_trad_xi_c_ddist_pv_stra);
  HarryPlotter::CheckAndStore(out, h_df_xi_qa_trad_xi_c_ddist_dv_stra);
  HarryPlotter::CheckAndStore(out, h_df_xi_qa_dl_pv_xi_c_ddist_dv_stra);
  HarryPlotter::CheckAndStore(out, h_df_xi_qa_dl_pv_xi_c_trad_stra);
  HarryPlotter::CheckAndStore(out, h_df_xi_qa_dl_dv_xi_c_ddist_dv_stra);
  HarryPlotter::CheckAndStore(out, h_df_xi_qa_dl_dv_xi_c_trad_stra);
  
  HarryPlotter::CheckAndStore(out, h_df_xi_c_im_xi_c_mass_stra); 
  
  //to xi_cc
  HarryPlotter::CheckAndStore(out, h_df_xi_c_qa_trad_diff_xi_xi_c_topo);
  HarryPlotter::CheckAndStore(out, h_df_xi_c_qa_xi_cc_ddca_topo);
  HarryPlotter::CheckAndStore(out, h_df_xi_c_qa_xi_cc_ddist_pv_topo);
  HarryPlotter::CheckAndStore(out, h_df_xi_c_qa_xi_cc_trad_topo);
  HarryPlotter::CheckAndStore(out, h_df_xi_c_qa_xi_c_dca_xy_topo);
  HarryPlotter::CheckAndStore(out, h_df_xi_c_qa_xi_c_dca_z_topo);
  HarryPlotter::CheckAndStore(out, h_df_xi_c_qa_xi_cc_dca_xy_topo);
  HarryPlotter::CheckAndStore(out, h_df_xi_c_qa_xi_cc_dca_z_topo);


  HarryPlotter::CheckAndStore(out, h_df_xi_c_qa_trad_diff_xi_xi_c_stra);
  HarryPlotter::CheckAndStore(out, h_df_xi_c_qa_xi_cc_ddca_stra);
  HarryPlotter::CheckAndStore(out, h_df_xi_c_qa_xi_cc_ddist_pv_stra);
  HarryPlotter::CheckAndStore(out, h_df_xi_c_qa_xi_cc_trad_stra);
  HarryPlotter::CheckAndStore(out, h_df_xi_c_qa_xi_c_dca_xy_stra);
  HarryPlotter::CheckAndStore(out, h_df_xi_c_qa_xi_c_dca_z_stra);
  HarryPlotter::CheckAndStore(out, h_df_xi_c_qa_xi_cc_dca_xy_stra);
  HarryPlotter::CheckAndStore(out, h_df_xi_c_qa_xi_cc_dca_z_stra);
    
  HarryPlotter::CheckAndStore(out,h_df_xi_c_qa_xi_cc_pt); 
  HarryPlotter::CheckAndStore(out,h_df_xi_c_qa_pi_pt); 

  HarryPlotter::CheckAndStore(out,h_df_xi_c_qa_pi_dca_xy);
  HarryPlotter::CheckAndStore(out,h_df_xi_c_qa_pi_dca_z);
  
  HarryPlotter::CheckAndStore(out,h_df_xi_c_qa_trad_xi_cc_ddist_pv_stra);
  HarryPlotter::CheckAndStore(out,h_df_xi_c_qa_dl_pv_xi_cc_trad_stra);
  
  HarryPlotter::CheckAndStore(out,h_df_xi_c_qa_dcaxyProd  );
  HarryPlotter::CheckAndStore(out,h_df_xi_c_qa_dcazProd   );
  HarryPlotter::CheckAndStore(out,h_df_xi_c_qa_dcaProdAdd );
  HarryPlotter::CheckAndStore(out,h_df_xi_c_qa_dcaProdMult);
  
  HarryPlotter::CheckAndStore(out,h_df_xi_c_qa_dcaxySqSum);
  HarryPlotter::CheckAndStore(out,h_df_xi_c_qa_dcazSqSum);

  HarryPlotter::CheckAndStore(out,h_df_xi_c_qa_dcaxy_prod_xic_pi);
  HarryPlotter::CheckAndStore(out,h_df_xi_c_qa_dcaz_prod_xic_pi);
  
  HarryPlotter::CheckAndStore(out,h_df_xi_c_xi_cc_mass_stra); 
  
  //xi_cc selected
  HarryPlotter::CheckAndStore(out, h_df_xi_cc_im_xi_cc_mass_stra_c1); 
  HarryPlotter::CheckAndStore(out, h_df_xi_cc_im_xi_cc_pt_c1); 
  HarryPlotter::CheckAndStore(out, h_df_xi_cc_im_xi_cc_mass_stra_c2); 
  HarryPlotter::CheckAndStore(out, h_df_xi_cc_im_xi_cc_pt_c2); 
  HarryPlotter::CheckAndStore(out, h_df_xi_cc_im_xi_cc_mass_stra_c3); 
  HarryPlotter::CheckAndStore(out, h_df_xi_cc_im_xi_cc_pt_c3); 
  HarryPlotter::CheckAndStore(out, h_df_xi_cc_im_xi_cc_mass_stra_c4); 
  HarryPlotter::CheckAndStore(out, h_df_xi_cc_im_xi_cc_pt_c4); 
  HarryPlotter::CheckAndStore(out, h_df_xi_cc_im_xi_cc_mass_stra_c4_2Hit); 
  HarryPlotter::CheckAndStore(out, h_df_xi_cc_im_xi_cc_pt_c4_2Hit); 
  HarryPlotter::CheckAndStore(out, h_df_xi_cc_im_xi_cc_mass_stra_c4_3Hit); 
  HarryPlotter::CheckAndStore(out, h_df_xi_cc_im_xi_cc_pt_c4_3Hit); 
  HarryPlotter::CheckAndStore(out, h_df_xi_cc_im_xi_cc_mass_stra_c5); 
  HarryPlotter::CheckAndStore(out, h_df_xi_cc_im_xi_cc_pt_c5); 
  HarryPlotter::CheckAndStore(out, h_df_xi_cc_im_xi_cc_mass_stra_c6); 
  HarryPlotter::CheckAndStore(out, h_df_xi_cc_im_xi_cc_pt_c6); 
  HarryPlotter::CheckAndStore(out, h_df_xi_cc_im_xi_cc_mass_stra_c7); 
  HarryPlotter::CheckAndStore(out, h_df_xi_cc_im_xi_cc_pt_c7); 
  HarryPlotter::CheckAndStore(out, h_df_xi_cc_im_xi_cc_mass_stra_c8); 
  HarryPlotter::CheckAndStore(out, h_df_xi_cc_im_xi_cc_pt_c8); 
  HarryPlotter::CheckAndStore(out, h_df_xi_cc_im_xi_cc_mass_stra_c9); 
  HarryPlotter::CheckAndStore(out, h_df_xi_cc_im_xi_cc_pt_c9); 


  HarryPlotter::CheckAndStore(out, h_cand_counter); 
  HarryPlotter::CheckAndStore(out, h_gen_xi_c_counter); 
  HarryPlotter::CheckAndStore(out, h_gen_xi_c_pt_eta_counter); 
  HarryPlotter::CheckAndStore(out, h_gen_xi_c_pt_y_counter); 
  
  HarryPlotter::CheckAndStore(out, h_gen_xi_cc_counter); 
  HarryPlotter::CheckAndStore(out, h_gen_xi_cc_pt_eta_counter); 
  HarryPlotter::CheckAndStore(out, h_gen_xi_cc_pt_y_counter); 
 
  HarryPlotter::CheckAndStore(out, cutCounter);

			    
  HarryPlotter::CheckAndStore(out, h_df_in_qa_pos_dca_xy_vs_hits);   
  HarryPlotter::CheckAndStore(out, h_df_in_qa_pos_dca_xy_vs_lmb_dl_pv);   
  HarryPlotter::CheckAndStore(out, h_df_in_qa_pos_dca_xy_vs_pos_pt);   

  HarryPlotter::CheckAndStore(out, h_df_in_qa_pos_beta_vs_pos_pt);   
  HarryPlotter::CheckAndStore(out, h_df_lmb_qa_pos_beta_expected_vs_pos_pt);   
  HarryPlotter::CheckAndStore(out, h_df_lmb_qa_pos_tof_dv_expected_vs_measured);   
  HarryPlotter::CheckAndStore(out, h_df_lmb_qa_pos_tof_pv_expected_vs_measured);   
  HarryPlotter::CheckAndStore(out, h_df_lmb_qa_pos_tof_dv_expected_vs_pv_expected);   
    
  HarryPlotter::CheckAndStore(out, h_df_in_qa_neg_dca_xy_vs_hits);   
  HarryPlotter::CheckAndStore(out, h_df_in_qa_neg_dca_xy_vs_lmb_dl_pv);   
  HarryPlotter::CheckAndStore(out, h_df_in_qa_neg_dca_xy_vs_neg_pt);   

  HarryPlotter::CheckAndStore(out, h_df_in_qa_neg_beta_vs_neg_pt);   
  HarryPlotter::CheckAndStore(out, h_df_lmb_qa_neg_beta_expected_vs_neg_pt);   
  HarryPlotter::CheckAndStore(out, h_df_lmb_qa_neg_tof_dv_expected_vs_measured);   
  HarryPlotter::CheckAndStore(out, h_df_lmb_qa_neg_tof_pv_expected_vs_measured);   
  HarryPlotter::CheckAndStore(out, h_df_lmb_qa_neg_tof_dv_expected_vs_pv_expected);   
  
  HarryPlotter::CheckAndStore(out, h_df_in_qa_neg_hits_vs_pos_hits);   
  HarryPlotter::CheckAndStore(out, h_df_in_qa_neg_dca_xy_vs_pos_dca_xy);   

  HarryPlotter::CheckAndStore(out, h_df_lmb_qa_bach_dca_xy_vs_hits);   
  HarryPlotter::CheckAndStore(out, h_df_lmb_qa_bach_dca_xy_vs_lmb_dl_pv);   
  HarryPlotter::CheckAndStore(out, h_df_lmb_qa_bach_dca_xy_vs_bach_pt);   

  HarryPlotter::CheckAndStore(out, h_df_lmb_qa_bach_beta_vs_bach_pt);   
  HarryPlotter::CheckAndStore(out, h_df_lmb_qa_bach_beta_expected_vs_bach_pt);   
  HarryPlotter::CheckAndStore(out, h_df_lmb_qa_bach_tof_dv_expected_vs_measured);   
  HarryPlotter::CheckAndStore(out, h_df_lmb_qa_bach_tof_pv_expected_vs_measured);   
  HarryPlotter::CheckAndStore(out, h_df_lmb_qa_bach_tof_dv_expected_vs_pv_expected);   
  
  
  out->Close(); 
  timer.Stop();
  timer.Print();
  return 0; 
}






  /*
  auto df_xi_cc_im_c1 = df_xi_c
    .Filter("XicXiccDecayRadDiffStra > 0")
    .Filter("TMath::Abs(fXicPionDCAxyToPV1) > 30")
    .Filter("TMath::Abs(fXicPionDCAzToPV1) > 35")
    .Filter("TMath::Abs(fXicPionDCAxyToPV2) > 30")
    .Filter("TMath::Abs(fXicPionDCAzToPV2) > 35")
    .Filter("TMath::Abs(fPicDCAxyToPVTopo) > 30")
    .Filter("TMath::Abs(fPicDCAzToPVTopo) > 35")
    ;
  
  //Fill some final histograms  
  auto h_df_xi_cc_im_xi_cc_mass_stra_c1 = df_xi_cc_im_c1.Histo1D({"df_xi_cc_im_xi_cc_mass_stra_c1", "xi_cc inv mass", 700, 2.6, 4.6}, "fXiccMassStraTrack"); 
  auto h_df_xi_cc_im_xi_cc_pt_c1 = df_xi_cc_im_c1.Filter(invMassXiccCut, {"fXiccMassStraTrack"}).Histo1D({"df_xi_cc_im_xi_cc_pt_c1", "pt selected", 200, 0, 20}, "lPtMCXiCC"); 
  auto out_counter_c1 = df_xi_cc_im_c1.Filter(invMassXiccCut, {"fXiccMassStraTrack"}).Count(); 

  //soft selection all cuts 
  auto df_xi_cc_im_c2 = df_xi_c
    .Filter("XicXiccDecayRadDiffStra > 0")
    .Filter("fXicDaughterDCAStraTrack < 20")
    .Filter("fXicDecayRadiusStraTrack > 0.004")
    .Filter("fXicInvDecayLengthToPVStra > 0.002")
    .Filter("fXiccDaughterDCAStraTrack < 20")
    .Filter("fXiccDecayRadiusStraTrack > 0.003")
    .Filter("fXiccInvDecayLengthToPVStra > 0.003")
    .Filter("TMath::Abs(fXiDCAxyToPVStraTrack) > 5")
    .Filter("TMath::Abs(fXiDCAzToPVStraTrack) > 10")
    .Filter("TMath::Abs(fXicDCAxyToPVStraTrack) > 10")
    .Filter("TMath::Abs(fXicDCAzToPVStraTrack) > 5")
    .Filter("TMath::Abs(fXicPionDCAxyToPV1) > 10")
    .Filter("TMath::Abs(fXicPionDCAzToPV1) > 15")
    .Filter("TMath::Abs(fXicPionDCAxyToPV2) > 10")
    .Filter("TMath::Abs(fXicPionDCAzToPV2) > 15")
    .Filter("TMath::Abs(fPicDCAxyToPVTopo) > 10")
    ;
  
  //Fill some final histograms  
  auto h_df_xi_cc_im_xi_cc_mass_stra_c2 = df_xi_cc_im_c2.Histo1D({"df_xi_cc_im_xi_cc_mass_stra_c2", "xi_cc inv mass", 700, 2.6, 4.6}, "fXiccMassStraTrack"); 
  auto h_df_xi_cc_im_xi_cc_pt_c2 = df_xi_cc_im_c2.Filter(invMassXiccCut, {"fXiccMassStraTrack"}).Histo1D({"df_xi_cc_im_xi_cc_pt_c2", "pt selected", 200, 0, 20}, "lPtMCXiCC"); 
  auto out_counter_c2 = df_xi_cc_im_c2.Filter(invMassXiccCut, {"fXiccMassStraTrack"}).Count(); 

  //soft selection all cuts expereimental xicc dca cuts
  auto df_xi_cc_im_c3 = df_xi_c
    .Filter("XicXiccDecayRadDiffStra > 0")
    .Filter("fXicDaughterDCAStraTrack < 20")
    .Filter("fXicDecayRadiusStraTrack > 0.004")
    .Filter("fXicInvDecayLengthToPVStra > 0.002")
    .Filter("fXiccDaughterDCAStraTrack < 20")
    .Filter("fXiccDecayRadiusStraTrack > 0.003")
    .Filter("fXiccInvDecayLengthToPVStra > 0.003")
    .Filter("TMath::Abs(fXiDCAxyToPVStraTrack) > 5")
    .Filter("TMath::Abs(fXiDCAzToPVStraTrack) > 10")
    .Filter("TMath::Abs(fXicDCAxyToPVStraTrack) > 10")
    .Filter("TMath::Abs(fXicDCAzToPVStraTrack) > 5")
    .Filter("TMath::Abs(fXicPionDCAxyToPV1) > 10")
    .Filter("TMath::Abs(fXicPionDCAzToPV1) > 15")
    .Filter("TMath::Abs(fXicPionDCAxyToPV2) > 10")
    .Filter("TMath::Abs(fXicPionDCAzToPV2) > 15")
    .Filter("TMath::Abs(fPicDCAxyToPVTopo) > 10")
    .Filter("TMath::Abs(fXiccDCAzToPVStraTrack) < 200")
    .Filter("TMath::Abs(fXiccDCAzToPVStraTrack) < 200")
    ;
  
  //Fill some final histograms  
  auto h_df_xi_cc_im_xi_cc_mass_stra_c3 = df_xi_cc_im_c3.Histo1D({"df_xi_cc_im_xi_cc_mass_stra_c3", "xi_cc inv mass", 700, 2.6, 4.6}, "fXiccMassStraTrack"); 
  auto h_df_xi_cc_im_xi_cc_pt_c3 = df_xi_cc_im_c3.Filter(invMassXiccCut, {"fXiccMassStraTrack"}).Histo1D({"df_xi_cc_im_xi_cc_pt_c3", "pt selected", 200, 0, 20}, "lPtMCXiCC"); 
  auto out_counter_c3 = df_xi_cc_im_c3.Filter(invMassXiccCut, {"fXiccMassStraTrack"}).Count(); 

  //harsh selection all cuts
  auto df_xi_cc_im_c4 = df_xi_c
    .Filter("XicXiccDecayRadDiffStra > 0")
    .Filter("fXicDaughterDCAStraTrack < 20")
    .Filter("fXicDecayRadiusStraTrack > 0.01")
    .Filter("fXicInvDecayLengthToPVStra > 0.01")
    .Filter("fXiccDaughterDCAStraTrack < 20")
    .Filter("fXiccDecayRadiusStraTrack > 0.003")
    .Filter("fXiccInvDecayLengthToPVStra > 0.003")
    .Filter("TMath::Abs(fXiDCAxyToPVStraTrack) > 5")
    .Filter("TMath::Abs(fXiDCAzToPVStraTrack) > 10")
    .Filter("TMath::Abs(fXicDCAxyToPVStraTrack) > 15")
    .Filter("TMath::Abs(fXicDCAzToPVStraTrack) > 10")
    .Filter("TMath::Abs(fXicPionDCAxyToPV1) > 10")
    .Filter("TMath::Abs(fXicPionDCAzToPV1) > 15")
    .Filter("TMath::Abs(fXicPionDCAxyToPV2) > 10")
    .Filter("TMath::Abs(fXicPionDCAzToPV2) > 15")
    .Filter("TMath::Abs(fPicDCAxyToPVTopo) > 10")
    ;
  
  //Fill some final histograms  
  auto h_df_xi_cc_im_xi_cc_mass_stra_c4 = df_xi_cc_im_c4.Histo1D({"df_xi_cc_im_xi_cc_mass_stra_c4", "xi_cc inv mass", 700, 2.6, 4.6}, "fXiccMassStraTrack"); 
  auto h_df_xi_cc_im_xi_cc_pt_c4 = df_xi_cc_im_c4.Filter(invMassXiccCut, {"fXiccMassStraTrack"}).Histo1D({"df_xi_cc_im_xi_cc_pt_c4", "pt selected", 200, 0, 20}, "lPtMCXiCC"); 
  auto out_counter_c4 = df_xi_cc_im_c4.Filter(invMassXiccCut, {"fXiccMassStraTrack"}).Count(); 

  
  auto df_xi_cc_im_c6 = df_xi_c
    .Filter("XicXiccDecayRadDiffStra > 0")
    .Filter("fXicDaughterDCAStraTrack < 20")
    .Filter("fXicDecayRadiusStraTrack > 0.004")
    .Filter("fXicInvDecayLengthToPVStra > 0.002")
    .Filter("fXicInvDecayLengthToDVStra < 0.1")
    .Filter("fXiccDaughterDCAStraTrack < 20")
    .Filter("fXiccDecayRadiusStraTrack > 0.005")
    .Filter("fXiccInvDecayLengthToPVStra > 0.004")
    .Filter("TMath::Abs(fXiDCAxyToPVStraTrack) > 5")
    .Filter("TMath::Abs(fXiDCAzToPVStraTrack) > 10")
    .Filter("TMath::Abs(fXicDCAxyToPVStraTrack) > 15")
    .Filter("TMath::Abs(fXicDCAzToPVStraTrack) > 10")
    .Filter("TMath::Abs(fXicPionDCAxyToPV1) > 10")
    .Filter("TMath::Abs(fXicPionDCAzToPV1) > 15")
    .Filter("TMath::Abs(fXicPionDCAxyToPV2) > 10")
    .Filter("TMath::Abs(fXicPionDCAzToPV2) > 20")
    .Filter("TMath::Abs(fPicDCAxyToPVTopo) > 20")
    .Filter("TMath::Abs(fXiccDCAzToPVStraTrack) < 200")
    .Filter("TMath::Abs(fXiccDCAzToPVStraTrack) < 200")    
    ;
  
  //Fill some final histograms  
  auto h_df_xi_cc_im_xi_cc_mass_stra_c6 = df_xi_cc_im_c6.Histo1D({"df_xi_cc_im_xi_cc_mass_stra_c6", "xi_cc inv mass", 700, 2.6, 4.6}, "fXiccMassStraTrack"); 
  auto h_df_xi_cc_im_xi_cc_pt_c6 = df_xi_cc_im_c6.Filter(invMassXiccCut, {"fXiccMassStraTrack"}).Histo1D({"df_xi_cc_im_xi_cc_pt_c6", "pt selected", 200, 0, 20}, "lPtMCXiCC"); 
  auto out_counter_c6 = df_xi_cc_im_c6.Filter(invMassXiccCut, {"fXiccMassStraTrack"}).Count(); 

  auto df_xi_cc_im_c7 = df_xi_c
    .Filter("XicXiccDecayRadDiffStra > 0")
    .Filter("fXicDaughterDCAStraTrack < 20")
    .Filter("fXicDecayRadiusStraTrack > 0.004")
    .Filter("fXicInvDecayLengthToPVStra > 0.002")
    .Filter("fXicInvDecayLengthToDVStra < 0.1")
    .Filter("fXiccDaughterDCAStraTrack < 20")
    .Filter("fXiccDecayRadiusStraTrack > 0.005")
    .Filter("fXiccInvDecayLengthToPVStra > 0.004")
    .Filter("TMath::Abs(fXiDCAxyToPVStraTrack) > 5")
    .Filter("TMath::Abs(fXiDCAzToPVStraTrack) > 10")
    .Filter("TMath::Abs(fXicDCAxyToPVStraTrack) > 15")
    .Filter("TMath::Abs(fXicDCAzToPVStraTrack) > 10")
    .Filter("TMath::Abs(fXicPionDCAxyToPV1) > 10")
    .Filter("TMath::Abs(fXicPionDCAzToPV1) > 15")
    .Filter("TMath::Abs(fXicPionDCAxyToPV2) > 10")
    .Filter("TMath::Abs(fXicPionDCAzToPV2) > 20")
    .Filter("TMath::Abs(fPicDCAxyToPVTopo) > 20")
    .Filter("TMath::Abs(fXiccDCAxyToPVStraTrack) < 200")
    .Filter("TMath::Abs(fXiccDCAzToPVStraTrack) < 200")    
    ;
  
  //Fill some final histograms  
  auto h_df_xi_cc_im_xi_cc_mass_stra_c7 = df_xi_cc_im_c7.Histo1D({"df_xi_cc_im_xi_cc_mass_stra_c7", "xi_cc inv mass", 700, 2.6, 4.6}, "fXiccMassStraTrack"); 
  auto h_df_xi_cc_im_xi_cc_pt_c7 = df_xi_cc_im_c7.Filter(invMassXiccCut, {"fXiccMassStraTrack"}).Histo1D({"df_xi_cc_im_xi_cc_pt_c7", "pt selected", 200, 0, 20}, "lPtMCXiCC"); 
  auto out_counter_c7 = df_xi_cc_im_c7.Filter(invMassXiccCut, {"fXiccMassStraTrack"}).Count(); 

  auto df_xi_cc_im_c8 = df_xi_c
    .Filter("XicXiccDecayRadDiffStra > 0")
    .Filter("XicXiccDecayRadDiffStra > 0")
    .Filter("XicXiccDecayRadDiffStra > 0")
    .Filter("fXicDaughterDCAStraTrack < 20")
    .Filter("fXicDecayRadiusStraTrack > 0.004")
    .Filter("fXicInvDecayLengthToPVStra > 0.002")
    .Filter("fXicInvDecayLengthToDVStra < 0.1")
    .Filter("fXiccDaughterDCAStraTrack < 20")
    .Filter("fXiccDecayRadiusStraTrack > 0.005")
    .Filter("fXiccInvDecayLengthToPVStra > 0.004")
    .Filter("TMath::Abs(fXiDCAxyToPVStraTrack) > 5")
    .Filter("TMath::Abs(fXiDCAzToPVStraTrack) > 10")
    .Filter("TMath::Abs(fXicDCAxyToPVStraTrack) > 15")
    .Filter("TMath::Abs(fXicDCAzToPVStraTrack) > 10")
    .Filter("TMath::Abs(fXicPionDCAxyToPV1) > 10")
    .Filter("TMath::Abs(fXicPionDCAzToPV1) > 15")
    .Filter("TMath::Abs(fXicPionDCAxyToPV2) > 10")
    .Filter("TMath::Abs(fXicPionDCAzToPV2) > 20")
    .Filter("TMath::Abs(fPicDCAxyToPVTopo) > 20")
    .Filter("TMath::Abs(fPicDCAzToPVTopo) > 20")
    .Filter("TMath::Abs(fXiccDCAxyToPVStraTrack) < 200")
    .Filter("TMath::Abs(fXiccDCAzToPVStraTrack) < 200")    
    ;
  
  //Fill some final histograms  
  auto h_df_xi_cc_im_xi_cc_mass_stra_c8 = df_xi_cc_im_c8.Histo1D({"df_xi_cc_im_xi_cc_mass_stra_c8", "xi_cc inv mass", 700, 2.6, 4.6}, "fXiccMassStraTrack"); 
  auto h_df_xi_cc_im_xi_cc_pt_c8 = df_xi_cc_im_c8.Filter(invMassXiccCut, {"fXiccMassStraTrack"}).Histo1D({"df_xi_cc_im_xi_cc_pt_c8", "pt selected", 200, 0, 20}, "lPtMCXiCC"); 
  auto out_counter_c8 = df_xi_cc_im_c8.Filter(invMassXiccCut, {"fXiccMassStraTrack"}).Count(); 

  auto df_xi_cc_im_c9 = df_xi_c
    .Filter("XicXiccDecayRadDiffStra > 0")
    .Filter("XicXiccDecayRadDiffStra > 0")
    .Filter("XicXiccDecayRadDiffStra > 0")
    .Filter("fXicDaughterDCAStraTrack < 20")
    .Filter("fXicDecayRadiusStraTrack > 0.004")
    .Filter("fXicInvDecayLengthToPVStra > 0.002")
    .Filter("fXicInvDecayLengthToDVStra < 0.1")
    .Filter("fXiccDaughterDCAStraTrack < 20")
    .Filter("fXiccDecayRadiusStraTrack > 0.005")
    .Filter("fXiccInvDecayLengthToPVStra > 0.004")
    .Filter("TMath::Abs(fXiDCAxyToPVStraTrack) > 5")
    .Filter("TMath::Abs(fXiDCAzToPVStraTrack) > 10")
    .Filter("TMath::Abs(fXicDCAxyToPVStraTrack) > 15")
    .Filter("TMath::Abs(fXicDCAzToPVStraTrack) > 10")
    .Filter("TMath::Abs(fXicPionDCAxyToPV1) > 10")
    .Filter("TMath::Abs(fXicPionDCAzToPV1) > 15")
    .Filter("TMath::Abs(fXicPionDCAxyToPV2) > 10")
    .Filter("TMath::Abs(fXicPionDCAzToPV2) > 20")
    .Filter("TMath::Abs(fPicDCAxyToPVTopo) > 20")
    .Filter("TMath::Abs(fPicDCAzToPVTopo) > 20")
    .Filter("TMath::Abs(fXiccDCAxyToPVStraTrack) < 50")
    .Filter("TMath::Abs(fXiccDCAzToPVStraTrack) < 50")    
    ;
  
  //Fill some final histograms  
  auto h_df_xi_cc_im_xi_cc_mass_stra_c9 = df_xi_cc_im_c9.Histo1D({"df_xi_cc_im_xi_cc_mass_stra_c9", "xi_cc inv mass", 700, 2.6, 4.6}, "fXiccMassStraTrack"); 
  auto h_df_xi_cc_im_xi_cc_pt_c9 = df_xi_cc_im_c9.Filter(invMassXiccCut, {"fXiccMassStraTrack"}).Histo1D({"df_xi_cc_im_xi_cc_pt_c9", "pt selected", 200, 0, 20}, "lPtMCXiCC"); 
  auto out_counter_c9 = df_xi_cc_im_c9.Filter(invMassXiccCut, {"fXiccMassStraTrack"}).Count();
    
  auto df_xi_cc_im_c9_2Hit = df_xi_c
    .Filter("fXiHitsAdded > 1")
    .Filter("XicXiccDecayRadDiffStra > 0")
    .Filter("XicXiccDecayRadDiffStra > 0")
    .Filter("XicXiccDecayRadDiffStra > 0")
    .Filter("fXicDaughterDCAStraTrack < 20")
    .Filter("fXicDecayRadiusStraTrack > 0.004")
    .Filter("fXicInvDecayLengthToPVStra > 0.002")
    .Filter("fXicInvDecayLengthToDVStra < 0.1")
    .Filter("fXiccDaughterDCAStraTrack < 20")
    .Filter("fXiccDecayRadiusStraTrack > 0.005")
    .Filter("fXiccInvDecayLengthToPVStra > 0.004")
    .Filter("TMath::Abs(fXiDCAxyToPVStraTrack) > 5")
    .Filter("TMath::Abs(fXiDCAzToPVStraTrack) > 10")
    .Filter("TMath::Abs(fXicDCAxyToPVStraTrack) > 15")
    .Filter("TMath::Abs(fXicDCAzToPVStraTrack) > 10")
    .Filter("TMath::Abs(fXicPionDCAxyToPV1) > 10")
    .Filter("TMath::Abs(fXicPionDCAzToPV1) > 15")
    .Filter("TMath::Abs(fXicPionDCAxyToPV2) > 10")
    .Filter("TMath::Abs(fXicPionDCAzToPV2) > 20")
    .Filter("TMath::Abs(fPicDCAxyToPVTopo) > 20")
    .Filter("TMath::Abs(fPicDCAzToPVTopo) > 20")
    .Filter("TMath::Abs(fXiccDCAxyToPVStraTrack) < 50")
    .Filter("TMath::Abs(fXiccDCAzToPVStraTrack) < 50")    
    ;
  
  //Fill some final histograms  
  auto h_df_xi_cc_im_xi_cc_mass_stra_c9_2Hit = df_xi_cc_im_c9_2Hit.Histo1D({"df_xi_cc_im_xi_cc_mass_stra_c9_2Hit", "xi_cc inv mass", 700, 2.6, 4.6}, "fXiccMassStraTrack"); 
  auto h_df_xi_cc_im_xi_cc_pt_c9_2Hit = df_xi_cc_im_c9_2Hit.Filter(invMassXiccCut, {"fXiccMassStraTrack"}).Histo1D({"df_xi_cc_im_xi_cc_pt_c9_2Hit", "pt selected", 200, 0, 20}, "lPtMCXiCC"); 
  auto out_counter_c9_2Hit = df_xi_cc_im_c9_2Hit.Filter(invMassXiccCut, {"fXiccMassStraTrack"}).Count();
  
  auto df_xi_cc_im_c9_3Hit = df_xi_c
    .Filter("fXiHitsAdded > 2")
    .Filter("XicXiccDecayRadDiffStra > 0")
    .Filter("XicXiccDecayRadDiffStra > 0")
    .Filter("fXicDaughterDCAStraTrack < 20")
    .Filter("fXicDecayRadiusStraTrack > 0.004")
    .Filter("fXicInvDecayLengthToPVStra > 0.002")
    .Filter("fXicInvDecayLengthToDVStra < 0.1")
    .Filter("fXiccDaughterDCAStraTrack < 20")
    .Filter("fXiccDecayRadiusStraTrack > 0.005")
    .Filter("fXiccInvDecayLengthToPVStra > 0.004")
    .Filter("TMath::Abs(fXiDCAxyToPVStraTrack) > 5")
    .Filter("TMath::Abs(fXiDCAzToPVStraTrack) > 10")
    .Filter("TMath::Abs(fXicDCAxyToPVStraTrack) > 15")
    .Filter("TMath::Abs(fXicDCAzToPVStraTrack) > 10")
    .Filter("TMath::Abs(fXicPionDCAxyToPV1) > 10")
    .Filter("TMath::Abs(fXicPionDCAzToPV1) > 15")
    .Filter("TMath::Abs(fXicPionDCAxyToPV2) > 10")
    .Filter("TMath::Abs(fXicPionDCAzToPV2) > 20")
    .Filter("TMath::Abs(fPicDCAxyToPVTopo) > 20")
    .Filter("TMath::Abs(fPicDCAzToPVTopo) > 20")
    .Filter("TMath::Abs(fXiccDCAxyToPVStraTrack) < 50")
    .Filter("TMath::Abs(fXiccDCAzToPVStraTrack) < 50")    
    ;
  
  //Fill some final histograms  
  auto h_df_xi_cc_im_xi_cc_mass_stra_c9_3Hit = df_xi_cc_im_c9_3Hit.Histo1D({"df_xi_cc_im_xi_cc_mass_stra_c9_3Hit", "xi_cc inv mass", 700, 2.6, 4.6}, "fXiccMassStraTrack"); 
  auto h_df_xi_cc_im_xi_cc_pt_c9_3Hit = df_xi_cc_im_c9_3Hit.Filter(invMassXiccCut, {"fXiccMassStraTrack"}).Histo1D({"df_xi_cc_im_xi_cc_pt_c9_3Hit", "pt selected", 200, 0, 20}, "lPtMCXiCC"); 
  auto out_counter_c9_3Hit = df_xi_cc_im_c9_3Hit.Filter(invMassXiccCut, {"fXiccMassStraTrack"}).Count();

  
  reduction = *out_counter_c3/inputCounter;
  std::cout << reduction << std::endl; 
  cutCounter->SetBinContent(3, reduction);
  reduction = *out_counter_c4/inputCounter;
  std::cout << reduction << std::endl; 
  cutCounter->SetBinContent(4, reduction);
  reduction = *out_counter_c5/inputCounter;
  std::cout << reduction << std::endl; 
  cutCounter->SetBinContent(5, reduction);
  reduction = *out_counter_c6/inputCounter;
  std::cout << reduction << std::endl; 
  cutCounter->SetBinContent(6, reduction);
  reduction = *out_counter_c7/inputCounter;
  std::cout << reduction << std::endl; 
  cutCounter->SetBinContent(7, reduction);
  reduction = *out_counter_c8/inputCounter ;
  std::cout << reduction << std::endl; 
  cutCounter->SetBinContent(8, reduction);
  reduction = *out_counter_c9/inputCounter ;
  std::cout << reduction << std::endl; 
  cutCounter->SetBinContent(9, reduction);
  reduction = *out_counter_c9_2Hit/inputCounter ;
  std::cout << reduction << std::endl; 
  cutCounter->SetBinContent(10, reduction);
  reduction = *out_counter_c9_3Hit/inputCounter ;
  std::cout << reduction << std::endl; 
  cutCounter->SetBinContent(11, reduction);
  reduction = *out_counter_c10/inputCounter;
  std::cout << reduction << std::endl; 
  cutCounter->SetBinContent(12, reduction);

    HarryPlotter::CheckAndStore(out, h_df_xi_cc_im_xi_cc_mass_stra_c1); 
  HarryPlotter::CheckAndStore(out, h_df_xi_cc_im_xi_cc_pt_c1); 
  HarryPlotter::CheckAndStore(out, h_df_xi_cc_im_xi_cc_mass_stra_c2); 
  HarryPlotter::CheckAndStore(out, h_df_xi_cc_im_xi_cc_pt_c2); 
  HarryPlotter::CheckAndStore(out, h_df_xi_cc_im_xi_cc_mass_stra_c3); 
  HarryPlotter::CheckAndStore(out, h_df_xi_cc_im_xi_cc_pt_c3); 
  HarryPlotter::CheckAndStore(out, h_df_xi_cc_im_xi_cc_mass_stra_c4); 
  HarryPlotter::CheckAndStore(out, h_df_xi_cc_im_xi_cc_pt_c4); 
  HarryPlotter::CheckAndStore(out, h_df_xi_cc_im_xi_cc_mass_stra_c6); 
  HarryPlotter::CheckAndStore(out, h_df_xi_cc_im_xi_cc_pt_c6); 
  HarryPlotter::CheckAndStore(out, h_df_xi_cc_im_xi_cc_mass_stra_c7); 
  HarryPlotter::CheckAndStore(out, h_df_xi_cc_im_xi_cc_pt_c7); 
  HarryPlotter::CheckAndStore(out, h_df_xi_cc_im_xi_cc_mass_stra_c8); 
  HarryPlotter::CheckAndStore(out, h_df_xi_cc_im_xi_cc_pt_c8); 
  HarryPlotter::CheckAndStore(out, h_df_xi_cc_im_xi_cc_mass_stra_c9); 
  HarryPlotter::CheckAndStore(out, h_df_xi_cc_im_xi_cc_pt_c9); 
  HarryPlotter::CheckAndStore(out, h_df_xi_cc_im_xi_cc_mass_stra_c9_2Hit); 
  HarryPlotter::CheckAndStore(out, h_df_xi_cc_im_xi_cc_pt_c9_2Hit); 
  HarryPlotter::CheckAndStore(out, h_df_xi_cc_im_xi_cc_mass_stra_c9_3Hit); 
  HarryPlotter::CheckAndStore(out, h_df_xi_cc_im_xi_cc_pt_c9_3Hit); 

 */
