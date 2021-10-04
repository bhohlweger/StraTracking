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
  
  int addedHitsMin = argv[5]?atoi(argv[5]):1; 
  
  int noXi = argv[6]?atoi(argv[6]):0; 
  bool ForceNoXi = (noXi==0)?false:true; 
      
  int wrongAssociationMode = argv[7]?atoi(argv[7]):0; 
    
  HarryPlotter::StyleBox(); 
 
  ROOT::EnableImplicitMT(5); // Tell ROOT you want to go parallel          
  TString filePath = TString::Format("%s", fileName); 
  
  auto h_LongTracks = new TH1D("NLongTrack", "NLongTrack", 10000, 0, 10000); 
  
  auto h_cand_counter = new TH1D("df_xi_c_candCounter", "candCounter", 1, 0, 1); 

  auto h_gen_xi_pt_eta_counter = new TH2D("ptetaXiGen", "candCounter", 200, 0, 20, 30, -1.5, 1.5); 
  auto h_gen_xi_pt_y_counter = new TH2D("ptyXiGen", "candCounter", 200, 0, 20, 30, -1.5, 1.5); 

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
    TH1D* hLongTracks = (TH1D*)inFile->Get("hNLongTracks");
    TH1D* evtCounter = (TH1D*)inFile->Get("hEventCounter"); 

    TH2D* ptXipteta = (TH2D*)inFile->Get("hPtEtaGeneratedXi"); 
    TH2D* ptXipty = (TH2D*)inFile->Get("hPtYGeneratedXi"); 
    
    TH1D* ptXicGen = (TH1D*)inFile->Get("hXiCGeneratedPt");
    TH2D* ptXicpteta = (TH2D*)inFile->Get("hPtEtaGeneratedXiC"); 
    TH2D* ptXicpty = (TH2D*)inFile->Get("hPtYGeneratedXiC"); 
    
    TH1D* ptXiccGen = (TH1D*)inFile->Get("hXiCCGeneratedPt");
    TH2D* ptXiccpteta = (TH2D*)inFile->Get("hPtEtaGeneratedXiCC"); 
    TH2D* ptXiccpty = (TH2D*)inFile->Get("hPtYGeneratedXiCC"); 

    if (!hLongTracks||!evtCounter||
	!ptXipteta||!ptXipty||
	!ptXicGen||!ptXicpteta||!ptXicpty||
	!ptXiccGen||!ptXiccpteta||!ptXiccpty) { 
      inputFailures++; 
      std::cout << "Zis is vehry bad, ze generation histograms are missing, Guenther! No histogram, no chain! \n"; 
      
      if(!hLongTracks )std::cout << " Missing hLongTracks \n"; 
      if(!evtCounter  )std::cout << " Missing evtCounter  \n";
      if(!ptXipteta   )std::cout << " Missing ptXipteta   \n";
      if(!ptXipty     )std::cout << " Missing ptXipty     \n";
      if(!ptXicGen    )std::cout << " Missing ptXicGen    \n";
      if(!ptXicpteta  )std::cout << " Missing ptXicpteta  \n";
      if(!ptXicpty    )std::cout << " Missing ptXicpty    \n";
      if(!ptXiccGen   )std::cout << " Missing ptXiccGen   \n";
      if(!ptXiccpteta )std::cout << " Missing ptXiccpteta \n";
      if(!ptXiccpty   )std::cout << " Missing ptXiccpty   \n";
    } else { 
      h_LongTracks->Add(hLongTracks);
      
      h_cand_counter->Add(evtCounter); 
      
      h_gen_xi_pt_eta_counter->Add(ptXipteta); 
      h_gen_xi_pt_y_counter->Add(ptXipty); 
      
      h_gen_xi_c_counter->Add(ptXicGen); 	
      h_gen_xi_c_pt_eta_counter->Add(ptXicpteta); 
      h_gen_xi_c_pt_y_counter->Add(ptXicpty); 
      
      h_gen_xi_cc_counter->Add(ptXiccGen); 	
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
    bool oneTimeError = true; 
    for (auto fileObj : *files)  {
      auto file = (TSystemFile*) fileObj;
      TString inSubDirFile = TString::Format("%s/%s/xicc.treeoutput.root", filePath.Data(), file->GetName()).Data(); 
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

	TH1D* hLongTracks = (TH1D*)inFile->Get("hNLongTracks");
	TH1D* evtCounter = (TH1D*)inFile->Get("hEventCounter"); 

	TH2D* ptXipteta = (TH2D*)inFile->Get("hPtEtaGeneratedXi"); 
	TH2D* ptXipty = (TH2D*)inFile->Get("hPtYGeneratedXi"); 
    
	TH1D* ptXicGen = (TH1D*)inFile->Get("hXiCGeneratedPt");
	TH2D* ptXicpteta = (TH2D*)inFile->Get("hPtEtaGeneratedXiC"); 
	TH2D* ptXicpty = (TH2D*)inFile->Get("hPtYGeneratedXiC"); 
    
	TH1D* ptXiccGen = (TH1D*)inFile->Get("hXiCCGeneratedPt");
	TH2D* ptXiccpteta = (TH2D*)inFile->Get("hPtEtaGeneratedXiCC"); 
	TH2D* ptXiccpty = (TH2D*)inFile->Get("hPtYGeneratedXiCC"); 
	
	if (!hLongTracks||!evtCounter||
	    !ptXipteta||!ptXipty||
	    !ptXicGen||!ptXicpteta||!ptXicpty||
	    !ptXiccGen||!ptXiccpteta||!ptXiccpty) { 
	  if (oneTimeError) { 
	    if(!hLongTracks )std::cout << " Missing hLongTracks \n"; 
	    if(!evtCounter  )std::cout << " Missing evtCounter  \n";
	    if(!ptXipteta   )std::cout << " Missing ptXipteta   \n";
	    if(!ptXipty     )std::cout << " Missing ptXipty     \n";
	    if(!ptXicGen    )std::cout << " Missing ptXicGen    \n";
	    if(!ptXicpteta  )std::cout << " Missing ptXicpteta  \n";
	    if(!ptXicpty    )std::cout << " Missing ptXicpty    \n";
	    if(!ptXiccGen   )std::cout << " Missing ptXiccGen   \n";
	    if(!ptXiccpteta )std::cout << " Missing ptXiccpteta \n";
	    if(!ptXiccpty   )std::cout << " Missing ptXiccpty   \n";
	    oneTimeError = false; 
	  }
	  inFile->Close(); 
	  inputFailures++; 
	  continue; 
	} else { 
	  h_LongTracks->Add(hLongTracks);
      
	  h_cand_counter->Add(evtCounter); 
      
	  h_gen_xi_pt_eta_counter->Add(ptXipteta); 
	  h_gen_xi_pt_y_counter->Add(ptXipty); 
      
	  h_gen_xi_c_counter->Add(ptXicGen); 	
	  h_gen_xi_c_pt_eta_counter->Add(ptXicpteta); 
	  h_gen_xi_c_pt_y_counter->Add(ptXicpty); 
      
	  h_gen_xi_cc_counter->Add(ptXiccGen); 	
	  h_gen_xi_cc_pt_eta_counter->Add(ptXiccpteta); 
	  h_gen_xi_cc_pt_y_counter->Add(ptXiccpty); 

	  input.Add(inSubDirFile);
	  inputFiles++;
	  inFile->Close();
	}
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
  float invMassDiffLmb = 0.005; //8 MeV/c2 mass window 
  auto invMassLmbCut = [&invMassDiffLmb, &lmbMass](float invMass) { return (TMath::Abs(invMass-lmbMass) < invMassDiffLmb); }; 
  auto radCut = [&radDiffMax](float radDiff) { return (radDiff < radDiffMax);};
  auto decLengthLmb = [&lmbMass](float len, float mom){ return TMath::Abs(mom)>1e-4?len*lmbMass/mom:-999; }; 
  
  //Xi cuts
  float xiMass = 1.322; 
  float invMassDiffXi = 0.005; 
  
  auto invMassXiCut = [&invMassDiffXi, &xiMass](float invMass) { return (TMath::Abs(invMass-xiMass) < invMassDiffXi); }; 
  auto hitsCut = [&addedHitsMin](int AddedHits) { return (AddedHits >= addedHitsMin); }; 
  
  auto decLengthXi = [&xiMass](float len, float mom){ return TMath::Abs(mom)>1e-4?len*xiMass/mom:-999; }; 

  //Xic cuts
  float xicMass = 2.46793; 
  float invMassDiffXic = 0.021; //8 MeV/c2 mass window 
  auto decLengthXic = [&xicMass](float len, float mom){ return TMath::Abs(mom)>1e-4?len*xicMass/mom:-999; }; 
  
  auto invMassXicCut = [&invMassDiffXic, &xicMass](float invMass) { return (TMath::Abs(invMass-xicMass) < invMassDiffXic); }; 
  

  //Xicc cuts
  float xiccMass = 3.621;   
  float invMassDiffXicc = 0.030; //8 MeV/c2 mass window 
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
  
  std::cout << "Utilizing the full beauty of strangeness hits, requested nHits = "<<  addedHitsMin << "\n";
  
  auto h_df_identified = df
    .Histo2D({"df_pt_vs_eta_ident", "pt selected", 200, 0, 20, 30, -1.5, 1.5}, "fPtMCXiCC", "fXiCCEta"); 

  auto df_precut = df
    .Filter("fV0DecayRadius > 0.45", "safeteyPrecuts_1")
    .Filter("fXiDecayRadius > 0.45", "safeteyPrecuts_2")
    .Filter("TMath::Abs(fPositiveDCAxy)> 50", "safeteyPrecuts_3")
    .Filter("TMath::Abs(fPositiveDCAz )> 40", "safeteyPrecuts_5")
    .Filter("TMath::Abs(fNegativeDCAxy)> 100", "safeteyPrecuts_6")
    .Filter("TMath::Abs(fNegativeDCAz )> 50", "safeteyPrecuts_7")
    .Filter("TMath::Abs(fBachelorDCAxy)> 40", "safeteyPrecuts_8")
    .Filter("TMath::Abs(fBachelorDCAz )> 40", "safeteyPrecuts_9")
    .Filter("TMath::Abs(fPic1DCAxyToPV         )>10", "safeteyPrecuts_10")
    .Filter("TMath::Abs(fPic1DCAzToPV          )>10", "safeteyPrecuts_11")
    .Filter("TMath::Abs(fPic2DCAxyToPV         )>10", "safeteyPrecuts_12")
    .Filter("TMath::Abs(fPic2DCAzToPV          )>10", "safeteyPrecuts_13")
    .Filter("TMath::Abs(fPiccDCAxyToPV         )>10", "safeteyPrecuts_14")
    .Filter("TMath::Abs(fPiccDCAzToPV          )>10", "safeteyPrecuts_15")
    .Filter("fPiC1Pt > 0.15", "safeteyPrecuts_16")
    .Filter("fPiC2Pt > 0.15", "safeteyPrecuts_17")
    .Filter("fPiCCPt > 0.3", "safeteyPrecuts_18")
    .Filter("TMath::Abs(fLambdaMass-1.116) < 0.005", "safeteyPrecuts_19")
    .Filter("fPtXi > 0.", "safeteyPrecuts_20")
    .Filter("fV0DauDCA   < 400", "safeteyPrecuts_21")
    .Filter("TMath::Abs(fXiMass-1.322) < 0.005", "safeteyPrecuts_22")
    .Filter("fXiDauDCA < 400", "safeteyPrecuts_23")
    .Filter("TMath::Abs(fXicMass-2.468) < 0.030", "safeteyPrecuts_24")
    .Define("fPosTOFDiffInner",  "fPositiveInnerTOF20Signal-fPositiveInnerExpectedSignal")
    .Filter("fPosTOFDiffInner < 75", "safeteyPrecuts_25")
    .Define("fNegTOFDiffInner",  "fNegativeInnerTOF20Signal-fNegativeInnerExpectedSignal")
    .Filter("fNegTOFDiffInner < 75", "safeteyPrecuts_26")
    .Define("fBachTOFDiffInner", "fBachelorInnerTOF20Signal-fBachelorInnerExpectedSignal")
    .Filter("fBachTOFDiffInner < 75", "safeteyPrecuts_27")
    .Define("fPic1TOFDiffInner", "fPic1InnerTOF20Signal-fPic1InnerExpectedSignal")
    .Filter("fPic1TOFDiffInner < 75", "safeteyPrecuts_28")
    .Define("fPic2TOFDiffInner", "fPic2InnerTOF20Signal-fPic2InnerExpectedSignal")
    .Filter("fPic2TOFDiffInner < 75", "safeteyPrecuts_29")
    .Define("fPiccTOFDiffInner", "fPiccInnerTOF20Signal-fPiccInnerExpectedSignal")
    .Filter("fPiccTOFDiffInner < 75", "safeteyPrecuts_30")
    .Filter("fXicDaughterDCA < 50", "safeteyPrecuts_31")
    .Filter("fXiccDaughterDCA < 30", "safeteyPrecuts_32")
    .Filter("TMath::Abs(fXiccMass-3.621) < 0.65", "safeteyPrecutsOut")
    ; 

  auto df_ForceXi = ForceNoXi?df_precut.Filter("!fTrueXi","noTrueXis"):df_precut.Filter("fTrueXi||!fTrueXi","TrueAndFalseXis");

  auto df_in = (
		ExclusiveSignal?
		df_ForceXi
		.Filter("fTrueXicc","trueXiccs")
		:df_ForceXi
		.Filter("!fTrueXicc","fakeXiccs")
		)
    .Filter(associations,{"fTrueXicc", "fTrueXic", "fPiccUsed", "fPicUsed"}, "Associations")
    .Define("fXiccPDGMass", [&xiccMass]() {return xiccMass;})
    .Define("fXiccY", HarryPlotter::YFromMomentum, {"fPXiCC", "fPtXiCC", "fXiccPDGMass", "fXiCCEta"})
    .Filter("TMath::Abs(fXiCCEta)<1.5", "XiccEta")
    .Filter(pTCut, {"fPtXiCC"}, "pTXicc")
    ;

  auto h_df_in_im_xi_cc_mass = df_in.Histo1D({"h_df_in_im_xi_cc_mass", "xi_cc inv mass", 700, 2.6, 4.6}, "fXiccMass"); 
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

  auto h_df_in_qa_lmb_dca_xy = df_in_qa.Histo1D({"df_in_qa_lmb_dca_xy", "lmb dca xy pv", 1000, -1e4, 1e4}, "fV0DCAxyToPV");
  auto h_df_in_qa_lmb_dca_z = df_in_qa.Histo1D({"df_in_qa_lmb_dca_z", "lmb dca z pv", 1000, -1e4, 1e4}, "fV0DCAzToPV");
  auto h_df_in_qa_lmb_dca_xy_wide = df_in_qa.Histo1D({"df_in_qa_lmb_dca_xy_wide", "lmb dca xy pv", 1000, -1e5, 1e5}, "fV0DCAxyToPV");
  auto h_df_in_qa_lmb_dca_z_wide = df_in_qa.Histo1D({"df_in_qa_lmb_dca_z_wide", "lmb dca z pv", 1000, -1e5, 1e5}, "fV0DCAzToPV");
  
  auto h_df_in_qa_lmb_totp = df_in_qa.Histo1D({"df_in_qa_lmb_totp", "lmb tot p", 200, 0, 20}, "fV0TotalMomentum");   
  auto h_df_in_qa_lmb_ddist_pv = df_in_qa.Define("fLmbInvDecayLengthToPV", decLengthLmb, {"fV0DecayLength", "fV0TotalMomentum"}).Histo1D({"df_in_qa_lmb_ddist_pv", "lmb ddist to pv", 1000, 0, 40}, "fLmbInvDecayLengthToPV");
  auto h_df_in_qa_lmb_ddca = df_in_qa.Histo1D({"df_in_qa_lmb_ddca", "lmb prong dca", 500, 0, 100}, "fV0DauDCA"); 
  auto h_df_in_qa_lmb_ddca_wide = df_in_qa.Histo1D({"df_in_qa_lmb_ddca_wide", "lmb prong dca", 500, 0, 5000}, "fV0DauDCA"); 
  auto h_df_in_qa_lmb_trad = df_in_qa.Histo1D({"df_in_qa_lmb_trad", "lmb trad", 600, 0, 30}, "fV0DecayRadius"); 
  auto h_df_in_qa_lmb_trad_mc = df_in_qa.Histo1D({"df_in_qa_lmb_trad_mc", "lmb trad", 600, 0, 30}, "fV0DecayRadiusMC"); 
  auto h_df_in_qa_lmb_trad_diff = df_in_qa.Define("XiV0DecayRadDiff", "TMath::Abs(fV0DecayRadiusMC-fV0DecayRadius)").Histo1D({"df_in_qa_lmb_trad_diff", "lmb trad", 500, 0, 2}, "XiV0DecayRadDiff"); 
  
  auto h_df_in_qa_lmb_mass = df_in_qa.Histo1D({"df_in_qa_lmb_mass", "lmb inv mass", 750, 1., 1.8}, "fLambdaMass"); 
  auto h_df_in_qa_xi_mass = df_in_qa.Histo1D({"df_in_qa_xi_mass", "xi inv mass", 750, 1.2, 2}, "fXiMass"); 
  
  //Define future variables and select Lambdas 

  auto df_lmb_im = df_in
    .Define("XiDecayRadDiff", "TMath::Abs(fXiDecayRadiusMC-fXiDecayRadius)")
    .Define("XiLmbDecayRadDiff", "fV0DecayRadius-fXiDecayRadius")
    .Define("XiXicDecayRadDiff", "fXiDecayRadius-fXicDecayRadius")
    .Define("fLmbInvDecayLengthToPV", decLengthLmb, {"fV0DecayLength", "fV0TotalMomentum"})
    .Define("fXiInvDecayLengthToPV", decLengthXi, {"fXiDecayLength", "fXiTotalMomentum"})
    .Define("fXicInvDecayLengthToPV", decLengthXic, {"fXicDecayDistanceFromPV", "fPXiC"})   // Xi_c decay length from the PV stra tracking 
    .Define("fXiccInvDecayLengthToPV", decLengthXicc, {"fXiccDecayDistanceFromPV", "fPXiCC"}) // Xi_cc decay length form the PV stra tracking = this is the Xi_cc decay length! 

    .Define("fXiInvDecayLengthToDV", decLengthXi, {"fXiCtoXiLength", "fXiTotalMomentum"}) //this is the Xi decay length 
    .Define("fXicInvDecayLengthToDV", decLengthXic, {"fXiCCtoXiCLength", "fPXiC"})    //this is the Xi_c decay length 
    
    .Filter("TMath::Abs(fV0DCAxyToPV) < 5000", "fV0DCAxyToPV")
    .Filter("TMath::Abs(fV0DCAzToPV) < 7000", "fV0DCAzToPV")
    .Filter("fV0DauDCA < 300","fV0DauDCA")
    .Filter("fV0DecayRadius > 0.5","fV0DecayRadius")
    .Filter("fLmbInvDecayLengthToPV <  15","fLmbInvDecayLengthToPV")
    .Filter("TMath::Abs(fPositiveDCAxy) > 50","fPositiveDCAxy")
    .Filter("TMath::Abs(fPositiveDCAz) > 40","fPositiveDCAz")
    .Filter("TMath::Abs(fNegativeDCAxy) > 100","fNegativeDCAxy")
    .Filter("TMath::Abs(fNegativeDCAz) > 50","fNegativeDCAz")
    .Filter("TMath::Abs(fPosTOFDiffInner) < 50", "PosTOF")
    .Filter("TMath::Abs(fNegTOFDiffInner) < 50", "NegTOF")
    ;

  auto h_df_lmb_im_lmb_mass = df_lmb_im.Filter("fFirstCandidateXiCC","df_lmb_im_h_bool").Histo1D({"df_lmb_im_lmb_mass", "lmb inv mass", 750, 1., 1.8}, "fLambdaMass"); 

  auto df_lmb =  df_lmb_im
    .Filter(invMassLmbCut, {"fLambdaMass"}, "fLambdaMass")
    ;

  //Towards the Xi for Xic->Xi+2pi 
  //Fill some Histograms
  auto h_df_lmb_xi_cc_mass = df_lmb.Histo1D({"h_df_lmb_xi_cc_mass", "xi_cc inv mass", 700, 2.6, 4.6}, "fXiccMass"); 

  auto df_lmb_qa = df_lmb.Filter("fFirstCandidateXiCC", "df_lmb_h_bool");
  
  auto h_df_lmb_qa_pos_tof_diff_inner = df_lmb_qa.Histo1D({"h_df_lmb_qa_pos_tof_diff_inner", "beta expected vs measured", 2000, -2000, 2000}, "fPosTOFDiffInner");
  auto h_df_lmb_qa_neg_tof_diff_inner = df_lmb_qa.Histo1D({"h_df_lmb_qa_neg_tof_diff_inner", "beta expected vs measured", 2000, -2000, 2000}, "fNegTOFDiffInner");
  auto h_df_lmb_qa_bach_tof_diff_inner = df_lmb_qa.Histo1D({"h_df_lmb_qa_bach_tof_diff_inner", "beta expected vs measured", 2000, -2000, 2000}, "fBachTOFDiffInner");

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
  
  auto h_df_lmb_qa_xi_pt = df_lmb_qa.Histo1D({"df_lmb_qa_xi_pt", "xi pt", 200, 0, 20}, "fPtXi"); 
  auto h_df_lmb_qa_xi_ddist_pv = df_lmb_qa.Histo1D({"df_lmb_qa_xi_ddist_pv", "xi ddist to pv", 1000, 0, 20}, "fXiInvDecayLengthToPV");
  auto h_df_lmb_qa_xi_ddca = df_lmb_qa.Histo1D({"df_lmb_qa_xi_ddca", "xi prong dca", 500, 0, 2000}, "fXiDauDCA"); 
  auto h_df_lmb_qa_xi_trad = df_lmb_qa.Histo1D({"df_lmb_qa_xi_trad", "xi trad", 250, 0, 10}, "fXiDecayRadius"); 
  auto h_df_lmb_qa_xi_trad_mc = df_lmb_qa.Histo1D({"df_lmb_qa_xi_trad_mc", "xi trad", 250, 0, 10}, "fXiDecayRadiusMC"); 
  auto h_df_lmb_qa_xi_trad_diff = df_lmb_qa.Histo1D({"df_lmb_qa_xi_trad_diff", "xi trad", 250, 0, 2}, "XiDecayRadDiff"); 

  auto h_df_lmb_qa_trad_diff_lmb_xi = df_lmb_qa.Histo1D({"df_lmb_qa_trad_diff_lmb_xi", "lmb-xi trad", 500, -100, 150}, "XiLmbDecayRadDiff"); 

  auto h_df_lmb_qa_xi_mass = df_lmb_qa.Histo1D({"df_lmb_qa_xi_mass", "xi inv mass", 750, 1.2, 2}, "fXiMass"); 
  
  //Select Xis excluding hits to avoid cheating 
  auto df_xi_sel = df_lmb
    .Filter("fXiDecayRadius > 0.5","fXiDecayRadius")
    .Filter("XiLmbDecayRadDiff > 0","XiLmbDecayRadDiff")
    .Filter("fXiDauDCA < 300","fXiDauDCA")
    .Filter("fXiInvDecayLengthToPV < 12","fXiInvDecayLengthToPV")
    .Filter("TMath::Abs(fBachelorDCAxy) > 40","fBachelorDCAxy")
    .Filter("TMath::Abs(fBachelorDCAz) > 40","fBachelorDCAz")
    .Filter("TMath::Abs(fBachTOFDiffInner) < 50", "BachelorTOF")
    ;

  auto h_df_xi_sel_xi_mass = df_xi_sel.Filter("fFirstCandidateXiCC","df_xi_sel_h_bool").Histo1D({"df_xi_sel_xi_mass", "xi inv mass", 750, 1.2, 2}, "fXiMass"); 

  //Select Xis including Hits
  auto df_xi_im = df_xi_sel.
    Filter(hitsCut, {"fXiHitsAdded"},"fXiHitsAdded")
    ; 
  
  auto h_df_xi_im_xi_mass = df_xi_im.Filter("fFirstCandidateXiCC","df_xi_im_h_bool").Histo1D({"df_xi_im_xi_mass", "xi inv mass", 750, 1.2, 2}, "fXiMass"); 
  
  auto h_df_xi_im_xi_ddist_dv = df_xi_im.Filter("fFirstCandidateXiCC","df_xi_im_h_bool").Histo1D({"df_xi_im_xi_dist_dv", "xi decay dist", 1500, 0, 30}, "fXiInvDecayLengthToDV"); 

  auto df_xi = df_xi_im
    .Filter(invMassXiCut, {"fXiMass"},"fXiMass")
    ;
  auto h_df_xi_xi_cc_mass = df_xi.Histo1D({"h_df_xi_xi_cc_mass", "xi_cc inv mass", 700, 2.6, 4.6}, "fXiccMass"); 

  //Towards the Xi_c for Xi_cc->Xi_c+pi: 
  //Fill some histograms 
  
  auto df_xi_qa = df_xi.Filter("fFirstCandidateXiCC","df_xi_h_bool");

  auto h_df_xi_qa_trad_diff_xi_xi_c = df_xi_qa.Histo1D({"df_xi_qa_trad_diff_xi_xi_c", "xi-xi_c trad", 500, -50, 200}, "XiXicDecayRadDiff") ;
  
  auto h_df_xi_qa_xi_c_mass = df_xi_qa.Histo1D({"df_xi_qa_xi_c_mass", "xi_c inv mass", 700, 1.6, 3.2}, "fXicMass"); 
  auto h_df_xi_qa_xi_c_pt = df_xi_qa.Histo1D({"df_xi_qa_xi_c_pt", "xi_c pt", 200, 0, 20}, "fPXiC");  
  
  auto h_df_xi_qa_xi_c_ddca = df_xi_qa.Histo1D({"df_xi_qa_xi_c_ddca", "xi_c prong dca", 500, 0, 100}, "fXicDaughterDCA"); 
  auto h_df_xi_qa_xi_c_ddist_pv = df_xi_qa.Histo1D({"df_xi_qa_xi_c_dist_pv", "xi_c decay dist", 1500, 0, 0.30}, "fXicInvDecayLengthToPV"); //redefine as mL/p 
  auto h_df_xi_qa_xi_c_ddist_dv = df_xi_qa.Histo1D({"df_xi_qa_xi_c_dist_dv", "xi_c decay dist", 1500, 0, 0.30}, "fXicInvDecayLengthToDV"); 
  auto h_df_xi_qa_xi_c_trad = df_xi_qa.Histo1D({"df_xi_qa_xi_c_trad", "xi_c trad", 2000, 0, 0.4}, "fXicDecayRadius"); 
  
  auto h_df_xi_qa_xi_dca_xy = df_xi_qa.Histo1D({"df_xi_qa_xi_dca_xy", "xi dca xy", 1000, -500, 500}, "fXiDCAxyToPV");  
  auto h_df_xi_qa_xi_dca_z = df_xi_qa.Histo1D({"df_xi_qa_xi_dca_z", "xi dca z", 1000, -500, 500}, "fXiDCAzToPV");  
  
  //Pions (from the xic) 
  
  auto h_df_xi_qa_pi_one_dca_xy = df_xi_qa.Histo1D({"df_xi_qa_pi_one_dca_xy", "pi1 dca xy stra", 1000, -500, 500}, "fPic1DCAxyToPV");  
  auto h_df_xi_qa_pi_one_dca_z = df_xi_qa.Histo1D({"df_xi_qa_pi_one_dca_z", "pi1 dca z stra", 1000, -500, 500}, "fPic1DCAzToPV");  
  
  auto h_df_xi_qa_pi_two_dca_xy = df_xi_qa.Histo1D({"df_xi_qa_pi_two_dca_xy", "pi2 dca xy stra", 1000, -500, 500}, "fPic2DCAxyToPV");  
  auto h_df_xi_qa_pi_two_dca_z = df_xi_qa.Histo1D({"df_xi_qa_pi_two_dca_z", "pi2 dca z stra", 1000, -500, 500}, "fPic2DCAzToPV");  
  
  auto h_df_xi_qa_pi_one_pt = df_xi_qa.Histo1D({"df_xi_qa_pi_one_pt", "pi c pt", 200, 0, 20}, "fPiC1Pt");  
  auto h_df_xi_qa_pi_two_pt = df_xi_qa.Histo1D({"df_xi_qa_pi_two_pt", "pi c pt", 200, 0, 20}, "fPiC2Pt");  

  auto h_df_xi_qa_pic1_tof_diff_inner = df_xi_qa.Histo1D({"h_df_xi_qa_pic1_tof_diff_inner", "beta expected vs measured", 2000, -2000, 2000}, "fPic1TOFDiffInner");
  auto h_df_xi_qa_pic2_tof_diff_inner = df_xi_qa.Histo1D({"h_df_xi_qa_pic2_tof_diff_inner", "beta expected vs measured", 2000, -2000, 2000}, "fPic2TOFDiffInner");

  //study a bit decay length + trad 
  //cut on Trad and check decay lengths to pv and to dv 
  auto h_df_xi_qa_trad_xi_c_ddist_pv = df_xi_qa.Filter("fXicDecayRadius > 0.01","corrStudy_fXicDecayRadius").Histo1D({"df_xi_qa_trad_xi_c_dist_pv", "xi_c decay dist", 1500, 0, 0.30}, "fXicInvDecayLengthToPV"); 
  auto h_df_xi_qa_trad_xi_c_ddist_dv = df_xi_qa.Filter("fXicDecayRadius > 0.01","corrStudy_fXicDecayRadius").Histo1D({"df_xi_qa_trad_xi_c_dist_dv", "xi_c decay dist", 1500, 0, 0.30}, "fXicInvDecayLengthToDV"); 
  
  //cut on dl to pv 
  auto h_df_xi_qa_dl_pv_xi_c_ddist_dv =  df_xi_qa.Filter("fXicInvDecayLengthToPV > 0.003","corrStudy_fXicInvDecayLengthToPV").Histo1D({"df_xi_qa_dl_pv_xi_c_dist_dv", "xi_c decay dist", 1500, 0, 0.30}, "fXicInvDecayLengthToDV");
  auto h_df_xi_qa_dl_pv_xi_c_trad = df_xi_qa.Filter("fXicInvDecayLengthToPV > 0.003","corrStudy_fXicInvDecayLengthToPV").Histo1D({"df_xi_qa_dl_pv_xi_c_trad", "xi_c trad", 2000, 0, 0.4}, "fXicDecayRadius");
  
  //cut on dl to dv 
  auto h_df_xi_qa_dl_dv_xi_c_ddist_dv =  df_xi_qa.Filter("fXicInvDecayLengthToDV > 0.002","corrStudy_fXicInvDecayLengthToDV").Histo1D({"df_xi_qa_dl_dv_xi_c_dist_pv", "xi_c decay dist", 1500, 0, 0.30}, "fXicInvDecayLengthToPV");
  auto h_df_xi_qa_dl_dv_xi_c_trad = df_xi_qa.Filter("fXicInvDecayLengthToDV > 0.002","corrStudy_fXicInvDecayLengthToDV").Histo1D({"df_xi_qa_dl_dv_xi_c_trad", "xi_c trad", 2000, 0, 0.4}, "fXicDecayRadius");

  //Select the Xi_c 

  auto df_xi_c_im = df_xi
    .Define("XicXiccDecayRadDiff", "fXicDecayRadius-fXiccDecayRadius")
    .Filter("XiXicDecayRadDiff > 0","XiXicDecayRadDiff")
    ;
  
  auto h_df_xi_c_im_xi_c_mass = df_xi_c_im.Filter("fFirstCandidateXiCC","df_xi_c_im_h_bool").Histo1D({"df_xi_c_im_xi_c_mass", "xi_c inv mass", 700, 1.6, 3.2}, "fXicMass"); 
  auto df_xi_c = df_xi_c_im
    .Filter(invMassXicCut, {"fXicMass"},"fXicMass")
    .Filter("TMath::Abs(fPic1TOFDiffInner) < 50", "Pic1TOF")
    .Filter("TMath::Abs(fPic2TOFDiffInner) < 50", "Pic2TOF")
    .Filter("TMath::Abs(fPiccTOFDiffInner) < 50", "PiccTOF")
    ;
  
  auto df_xi_c_qa = df_xi_c.Filter("fFirstCandidateXiCC","df_xi_c_h_bool"); 
  
  //Towards the actual Xi_cc
  //Fill some histograms 
  
  auto h_df_xi_c_qa_trad_diff_xi_xi_c = df_xi_c_qa.Histo1D({"df_xi_c_qa_trad_diff_xi_c_xi_cc", "xi_c-xi_cc trad", 500, -100, 150}, "XicXiccDecayRadDiff") ;

  auto h_df_xi_c_qa_xi_cc_ddca = df_xi_c_qa.Histo1D({"df_xi_c_qa_xi_cc_ddca", "xi_cc prong dca", 500, 0, 100}, "fXiccDaughterDCA"); 
  auto h_df_xi_c_qa_xi_cc_ddist_pv = df_xi_c_qa.Histo1D({"df_xi_c_qa_xi_cc_dist_pv", "xi_cc decay dist", 3000, 0, 0.50}, "fXiccInvDecayLengthToPV"); 
  auto h_df_xi_c_qa_xi_cc_trad = df_xi_c_qa.Histo1D({"df_xi_c_qa_xi_cc_trad", "xi_cc trad", 2000, 0, 0.5}, "fXiccDecayRadius"); 
  
  auto h_df_xi_c_qa_xi_c_dca_xy = df_xi_c_qa.Histo1D({"df_xi_c_qa_xi_c_dca_xy", "xi_c dca xy stra", 1000, -500, 500}, "fXicDCAxyToPV");  
  auto h_df_xi_c_qa_xi_c_dca_z  = df_xi_c_qa.Histo1D({"df_xi_c_qa_xi_c_dca_z", "xi_c dca z stra", 1000, -500, 500}, "fXicDCAzToPV");  
  
  auto h_df_xi_c_qa_xi_cc_dca_xy = df_xi_c_qa.Histo1D({"df_xi_c_qa_xi_cc_dca_xy", "xi_cc dca xy stra", 1000, -500, 500}, "fXiccDCAxyToPV");  
  auto h_df_xi_c_qa_xi_cc_dca_z  = df_xi_c_qa.Histo1D({"df_xi_c_qa_xi_cc_dca_z", "xi_cc dca z stra", 1000, -500, 500}, "fXiccDCAzToPV");  

  //Pions (from the xicc)
 
  auto h_df_xi_c_qa_pi_dca_xy = df_xi_c_qa.Histo1D({"df_xi_c_qa_pi_dca_xy", "xi_c dca xy stra", 1000, -500, 500}, "fPiccDCAxyToPV");  
  auto h_df_xi_c_qa_pi_dca_z  = df_xi_c_qa.Histo1D({"df_xi_c_qa_pi_dca_z", "xi_c dca z stra", 1000, -500, 500}, "fPiccDCAzToPV");  

  auto h_df_xi_c_qa_picc_tof_diff_inner = df_xi_c_qa.Histo1D({"h_df_xi_c_qa_picc_tof_diff_inner", "beta expected vs measured", 2000, -2000, 2000}, "fPiccTOFDiffInner");
  
  auto df_xi_cc_qa = df_xi_c_qa
    .Filter("TMath::Abs(fPic1DCAxyToPV) > 10")
    .Filter("TMath::Abs(fPic1DCAzToPV) > 10")
    .Filter("TMath::Abs(fPic2DCAxyToPV) > 10")
    .Filter("TMath::Abs(fPic2DCAzToPV) > 10")
    .Filter("TMath::Abs(fPiccDCAxyToPV) > 10")
    .Filter("TMath::Abs(fPiccDCAzToPV) > 10"); 


  auto h_df_xi_cc_qa_xi_cc_pt = df_xi_cc_qa.Histo1D({"df_xi_cc_qa_xi_cc_pt", "xi_cc pt", 200, 0, 20}, "fPtXiCC");  
  
  auto h_df_xi_cc_qa_pi_c1_pt = df_xi_cc_qa.Histo2D({"df_xi_cc_qa_pi_c1_pt", "pi c1 pt", 800, 0, 20, 200, 0, 20}, "fPiC1Pt", "fPtXiCC");  
  auto h_df_xi_cc_qa_pi_c2_pt = df_xi_cc_qa.Histo2D({"df_xi_cc_qa_pi_c2_pt", "pi c2 pt", 800, 0, 20, 200, 0, 20}, "fPiC2Pt", "fPtXiCC");  
  auto h_df_xi_cc_qa_pi_cc_pt = df_xi_cc_qa.Histo2D({"df_xi_cc_qa_pi_cc_pt", "pi cc pt", 800, 0, 20, 200, 0, 20}, "fPiCCPt", "fPtXiCC");  

  //cut on Trad and check decay lengths 
  auto h_df_xi_c_qa_trad_xi_cc_ddist_pv = df_xi_c_qa.Filter("fXiccDecayRadius > 0.005","corrStudy_fXiccDecayRadius").Histo1D({"df_xi_c_qa_trad_xi_cc_dist_pv", "xi_cc decay dist", 3000, 0, 0.50}, "fXiccInvDecayLengthToPV"); 
  //cut on dl to pv 
  auto h_df_xi_c_qa_dl_pv_xi_cc_trad = df_xi_c_qa.Filter("(fXiccInvDecayLengthToPV > 0.003) && (fXiccInvDecayLengthToPV <  0.06)","corrStudy_fXiccInvDecayLengthToPV").Histo1D({"df_xi_c_qa_dl_pv_xi_cc_trad", "xi_cc trad", 2000, 0, 0.4}, "fXiccDecayRadius");
    
  auto h_df_xi_c_xi_cc_mass = df_xi_c.Filter("XicXiccDecayRadDiff > 0","XicXiccDecayRadDiff").Histo1D({"df_xi_c_xi_cc_mass", "xi_cc inv mass", 700, 2.6, 4.6}, "fXiccMass"); 
  
  //Select the Xi_cc
  auto df_xi_cc_im_c1 = df_xi_c
    .Filter("XicXiccDecayRadDiff > 0","c1_XicXiccDecayRadDiffStra")
    .Filter("fXicDaughterDCA < 20","c1_fXicDaughterDCA")
    .Filter("fXicDecayRadius > 0.004","c1_fXicDecayRadius")
    .Filter("fXicInvDecayLengthToPV > 0.002","c1_fXicInvDecayLengthToPVStra")
    .Filter("fXicInvDecayLengthToDV < 0.1","c1_fXicInvDecayLengthToDVStra")
    .Filter("fXiccDaughterDCA < 20","c1_fXiccDaughterDCA")
    .Filter("fXiccDecayRadius > 0.005","c1_fXiccDecayRadius")
    .Filter("fXiccInvDecayLengthToPV > 0.004","c1_fXiccInvDecayLengthToPVStra")
    .Filter("TMath::Abs(fXiDCAxyToPV) > 5","c1_fXiDCAxyToPV")
    .Filter("TMath::Abs(fXiDCAzToPV) > 10","c1_fXiDCAzToPV")
    .Filter("TMath::Abs(fXicDCAxyToPV) > 15","c1_fXicDCAxyToPV")
    .Filter("TMath::Abs(fXicDCAzToPV) > 10","c1_fXicDCAzToPV")
    .Filter("TMath::Abs(fPic1DCAxyToPV) > 10","c1_fXicPionDCAxyToPV1")
    .Filter("TMath::Abs(fPic1DCAzToPV) > 15","c1_fXicPionDCAzToPV1")
    .Filter("TMath::Abs(fPic2DCAxyToPV) > 10","c1_fXicPionDCAxyToPV2")
    .Filter("TMath::Abs(fPic2DCAzToPV) > 20","c1_fXicPionDCAzToPV2")
    .Filter("TMath::Abs(fPiccDCAxyToPV) > 20","c1_fPicDCAxyToPVTopo")
    .Filter("TMath::Abs(fPiccDCAzToPV) > 10","c1_fPicDCAzToPV")
    ;
  
  //Fill some final histograms  
  auto h_df_xi_cc_im_xi_cc_mass_c1 = df_xi_cc_im_c1.Histo1D({"df_xi_cc_im_xi_cc_mass_c1", "xi_cc inv mass", 700, 2.6, 4.6}, "fXiccMass"); 
  auto h_df_xi_cc_im_xi_cc_pt_c1 = df_xi_cc_im_c1.Filter(invMassXiccCut, {"fXiccMass"}, "c1_Ivm").Histo2D({"df_xi_cc_im_xi_cc_pt_vs_eta_c1", "pt selected", 200, 0, 20, 30, -1.5, 1.5}, "fPtMCXiCC", "fXiCCEta"); 
  auto h_df_xi_cc_im_xi_cc_bkg_pt_c1 = df_xi_cc_im_c1.Filter(invMassXiccCut, {"fXiccMass"}, "c1_Ivm").Histo1D({"df_xi_cc_im_xi_cc_pt_bkg_c1", "pt selected", 200, 0, 20}, "fPtXiCC"); 
  auto out_counter_c1 = df_xi_cc_im_c1.Filter(invMassXiccCut, {"fXiccMass"}).Count(); 
  
  //Select the Xi_cc
  auto df_xi_cc_im_c2 = df_xi_c
    .Filter("XicXiccDecayRadDiff > 0","c2_XicXiccDecayRadDiff")
    
    .Filter("fXicDaughterDCA < 12","c2_fXicDaughterDCA")
    .Filter("fXicDecayRadius > 0.004","c2_fXicDecayRadius")
    
    .Filter("fXicInvDecayLengthToPV > 0.002","c2_fXicInvDecayLengthToPVStra")
    .Filter("fXicInvDecayLengthToDV < 0.06","c2_fXicInvDecayLengthToDVStra")
    
    .Filter("fXiccDaughterDCA < 5","c2_fXiccDaughterDCA")
    .Filter("fXiccDecayRadius > 0.005","c2_fXiccDecayRadius")
    .Filter("(fXiccInvDecayLengthToPV > 0.004) && (fXiccInvDecayLengthToPV < 0.06)","c2_fXiccInvDecayLengthToPVStra")
    
    .Filter("TMath::Abs(fXiDCAxyToPV) > 5","c2_fXiDCAxyToPV")
    .Filter("TMath::Abs(fXiDCAzToPV) > 10","c2_fXiDCAzToPV")
    
    .Filter("TMath::Abs(fXicDCAxyToPV) > 15","c2_fXicDCAxyToPV")
    .Filter("TMath::Abs(fXicDCAzToPV) > 10","c2_fXicDCAzToPV")
    
    .Filter("TMath::Abs(fXiccDCAxyToPV) < 20","c2_fXiccDCAxyToPV")
    .Filter("TMath::Abs(fXiccDCAzToPV) < 20","c2_fXiccDCAzToPV")    
    
    .Filter("TMath::Abs(fPic1DCAxyToPV) > 10","c2_fXicPionDCAxyToPV1")
    .Filter("TMath::Abs(fPic1DCAzToPV) > 20","c2_fXicPionDCAzToPV1")
    .Filter("TMath::Abs(fPic2DCAxyToPV) > 10","c2_fXicPionDCAxyToPV2")
    .Filter("TMath::Abs(fPic2DCAzToPV) > 20","c2_fXicPionDCAzToPV2")
    .Filter("TMath::Abs(fPiccDCAxyToPV) > 20","c2_fPicDCAxyToPVTopo")
    .Filter("TMath::Abs(fPiccDCAzToPV) > 10","c2_fPicDCAzToPV")
    ;
  
  //Fill some final histograms  
  auto h_df_xi_cc_im_xi_cc_mass_c2 = df_xi_cc_im_c2.Histo1D({"df_xi_cc_im_xi_cc_mass_c2", "xi_cc inv mass", 700, 2.6, 4.6}, "fXiccMass"); 
  auto h_df_xi_cc_im_xi_cc_pt_c2 = df_xi_cc_im_c2.Filter(invMassXiccCut, {"fXiccMass"}, "c2_Ivm").Histo2D({"df_xi_cc_im_xi_cc_pt_vs_eta_c2", "pt selected", 200, 0, 20, 30, -1.5, 1.5}, "fPtMCXiCC", "fXiCCEta"); 
  auto h_df_xi_cc_im_xi_cc_bkg_pt_c2 = df_xi_cc_im_c2.Filter(invMassXiccCut, {"fXiccMass"}, "c2_Ivm").Histo1D({"df_xi_cc_im_xi_cc_pt_bkg_c2", "pt selected", 200, 0, 20}, "fPtXiCC"); 
  auto out_counter_c2 = df_xi_cc_im_c2.Filter(invMassXiccCut, {"fXiccMass"}).Count(); 
  
  //Select the Xi_cc
  auto df_xi_cc_im_c3 = df_xi_c
    .Filter("XicXiccDecayRadDiff > 0","c3_XicXiccDecayRadDiffStra")
    .Filter("fXicDaughterDCA < 12","c3_fXicDaughterDCA")
    .Filter("fXicDecayRadius > 0.01","c3_fXicDecayRadius")
    .Filter("fXicInvDecayLengthToPV > 0.003","c3_fXicInvDecayLengthToPVStra")
    .Filter("fXicInvDecayLengthToDV > 0.002","c3_fXicInvDecayLengthToDVStra")
    .Filter("fXiccDaughterDCA < 8","c3_fXiccDaughterDCA")
    .Filter("fXiccDecayRadius > 0.005","c3_fXiccDecayRadius")
    .Filter("fXiccInvDecayLengthToPV > 0.003","c3_fXiccInvDecayLengthToPVStra")
    .Filter("TMath::Abs(fXiDCAxyToPV) > 10","c3_fXiDCAxyToPV")
    .Filter("TMath::Abs(fXiDCAzToPV) > 10","c3_fXiDCAzToPV")
    .Filter("TMath::Abs(fXicDCAxyToPV) > 10","c3_fXicDCAxyToPV")
    .Filter("TMath::Abs(fXicDCAzToPV) > 10","c3_fXicDCAzToPV")
    .Filter("TMath::Abs(fXiccDCAxyToPV) < 15","c3_fXiccDCAxyToPV")
    .Filter("TMath::Abs(fXiccDCAzToPV) < 15","c3_fXiccDCAzToPV")    
    .Filter("TMath::Abs(fPic1DCAxyToPV) > 10","c3_fXicPionDCAxyToPV1")
    .Filter("TMath::Abs(fPic1DCAzToPV) > 15","c3_fXicPionDCAzToPV1")
    .Filter("TMath::Abs(fPic2DCAxyToPV) > 10","c3_fXicPionDCAxyToPV2")
    .Filter("TMath::Abs(fPic2DCAzToPV) > 20","c3_fXicPionDCAzToPV2")
    .Filter("TMath::Abs(fPiccDCAxyToPV) > 20","c3_fPicDCAxyToPVTopo")
    .Filter("TMath::Abs(fPiccDCAzToPV) > 10","c3_fPicDCAzToPV")
    ;
  
  //Fill some final histograms  
  auto h_df_xi_cc_im_xi_cc_mass_c3 = df_xi_cc_im_c3.Histo1D({"df_xi_cc_im_xi_cc_mass_c3", "xi_cc inv mass", 700, 2.6, 4.6}, "fXiccMass"); 
  auto h_df_xi_cc_im_xi_cc_pt_c3 = df_xi_cc_im_c3.Filter(invMassXiccCut, {"fXiccMass"}, "c3_Ivm").Histo2D({"df_xi_cc_im_xi_cc_pt_vs_eta_c3", "pt selected", 200, 0, 20, 30, -1.5, 1.5}, "fPtMCXiCC", "fXiCCEta"); 
  auto h_df_xi_cc_im_xi_cc_bkg_pt_c3 = df_xi_cc_im_c3.Filter(invMassXiccCut, {"fXiccMass"}, "c3_Ivm").Histo1D({"df_xi_cc_im_xi_cc_pt_bkg_c3", "pt selected", 200, 0, 20}, "fPtXiCC"); 
  auto out_counter_c3 = df_xi_cc_im_c3.Filter(invMassXiccCut, {"fXiccMass"}).Count(); 
  
  auto df_xi_cc_im_c4 = df_xi_c
    .Filter("XicXiccDecayRadDiff > 0","c4_XicXiccDecayRadDiff")
    .Filter("fXicDaughterDCA < 20","c4_fXicDaughterDCA")
    .Filter("fXicDecayRadius > 0.004","c4_fXicDecayRadius")
    .Filter("fXicInvDecayLengthToPV > 0.002","c4_fXicInvDecayLengthToPVStra")
    .Filter("fXicInvDecayLengthToDV < 0.1","c4_fXicInvDecayLengthToDVStra")
    .Filter("fXiccDaughterDCA < 20","c4_fXiccDaughterDCA")
    .Filter("fXiccDecayRadius > 0.005","c4_fXiccDecayRadius")
    .Filter("fXiccInvDecayLengthToPV > 0.004","c4_fXiccInvDecayLengthToPVStra")
    .Filter("TMath::Abs(fXiDCAxyToPV) > 5","c4_fXiDCAxyToPV")
    .Filter("TMath::Abs(fXiDCAzToPV) > 10","c4_fXiDCAzToPV")
    .Filter("TMath::Abs(fXicDCAxyToPV) > 20","c4_fXicDCAxyToPV")
    .Filter("TMath::Abs(fXicDCAzToPV) > 20","c4_fXicDCAzToPV")
    .Filter("TMath::Abs(fPic1DCAxyToPV) > 10","c4_fXicPionDCAxyToPV1")
    .Filter("TMath::Abs(fPic1DCAzToPV) > 15","c4_fXicPionDCAzToPV1")
    .Filter("TMath::Abs(fPic2DCAxyToPV) > 10","c4_fXicPionDCAxyToPV2")
    .Filter("TMath::Abs(fPic2DCAzToPV) > 20","c4_fXicPionDCAzToPV2")
    .Filter("TMath::Abs(fPiccDCAxyToPV) > 20","c4_fPicDCAxyToPV")
    .Filter("TMath::Abs(fPiccDCAzToPV) > 20","c4_fPicDCAzToPV")
    .Filter("TMath::Abs(fXiccDCAxyToPV) < 50","c4_fXiccDCAxyToPV")
    .Filter("TMath::Abs(fXiccDCAzToPV) < 50","c4_fXiccDCAzToPV")    
    ;
  
  //Fill some final histograms  
  auto h_df_xi_cc_im_xi_cc_mass_c4 = df_xi_cc_im_c4.Histo1D({"df_xi_cc_im_xi_cc_mass_c4", "xi_cc inv mass", 700, 2.6, 4.6}, "fXiccMass"); 
  auto h_df_xi_cc_im_xi_cc_pt_c4 = df_xi_cc_im_c4.Filter(invMassXiccCut, {"fXiccMass"}, "c4_Ivm").Histo2D({"df_xi_cc_im_xi_cc_pt_vs_eta_c4", "pt selected", 200, 0, 20, 30, -1.5, 1.5}, "fPtMCXiCC", "fXiCCEta"); 
  auto h_df_xi_cc_im_xi_cc_bkg_pt_c4 = df_xi_cc_im_c4.Filter(invMassXiccCut, {"fXiccMass"}, "c4_Ivm").Histo1D({"df_xi_cc_im_xi_cc_pt_bkg_c4", "pt selected", 200, 0, 20}, "fPtXiCC"); 
  auto out_counter_c4 = df_xi_cc_im_c4.Filter(invMassXiccCut, {"fXiccMass"}).Count();
 
  auto df_xi_cc_im_c4_2Hit = df_xi_c
    .Filter("fXiHitsAdded > 1")
    .Filter("XicXiccDecayRadDiff > 0","c4_2Hit_XicXiccDecayRadDiff")
    .Filter("fXicDaughterDCA < 20","c4_2Hit_fXicDaughterDCA")
    .Filter("fXicDecayRadius > 0.004","c4_2Hit_fXicDecayRadius")
    .Filter("fXicInvDecayLengthToPV > 0.002","c4_2Hit_fXicInvDecayLengthToPVStra")
    .Filter("fXicInvDecayLengthToDV < 0.1","c4_2Hit_fXicInvDecayLengthToDVStra")
    .Filter("fXiccDaughterDCA < 20","c4_2Hit_fXiccDaughterDCA")
    .Filter("fXiccDecayRadius > 0.005","c4_2Hit_fXiccDecayRadius")
    .Filter("fXiccInvDecayLengthToPV > 0.004","c4_2Hit_fXiccInvDecayLengthToPVStra")
    .Filter("TMath::Abs(fXiDCAxyToPV) > 5","c4_2Hit_fXiDCAxyToPV")
    .Filter("TMath::Abs(fXiDCAzToPV) > 10","c4_2Hit_fXiDCAzToPV")
    .Filter("TMath::Abs(fXicDCAxyToPV) > 20","c4_2Hit_fXicDCAxyToPV")
    .Filter("TMath::Abs(fXicDCAzToPV) > 20","c4_2Hit_fXicDCAzToPV")
    .Filter("TMath::Abs(fPic1DCAxyToPV) > 10","c4_2Hit_fXicPionDCAxyToPV1")
    .Filter("TMath::Abs(fPic1DCAzToPV) > 15","c4_2Hit_fXicPionDCAzToPV1")
    .Filter("TMath::Abs(fPic2DCAxyToPV) > 10","c4_2Hit_fXicPionDCAxyToPV2")
    .Filter("TMath::Abs(fPic2DCAzToPV) > 20","c4_2Hit_fXicPionDCAzToPV2")
    .Filter("TMath::Abs(fPiccDCAxyToPV) > 20","c4_2Hit_fPicDCAxyToPVTopo")
    .Filter("TMath::Abs(fPiccDCAzToPV) > 20","c4_2Hit_fPicDCAzToPV")
    .Filter("TMath::Abs(fXiccDCAxyToPV) < 50","c4_2Hit_fXiccDCAxyToPV")
    .Filter("TMath::Abs(fXiccDCAzToPV) < 50","c4_2Hit_fXiccDCAzToPV")    
    ;
  
  //Fill some final histograms  
  auto h_df_xi_cc_im_xi_cc_mass_c4_2Hit = df_xi_cc_im_c4_2Hit.Histo1D({"df_xi_cc_im_xi_cc_mass_c4_2Hit", "xi_cc inv mass", 700, 2.6, 4.6}, "fXiccMass"); 
  auto h_df_xi_cc_im_xi_cc_pt_c4_2Hit = df_xi_cc_im_c4_2Hit.Filter(invMassXiccCut, {"fXiccMass"}, "c4_2Hit_Ivm").Histo2D({"df_xi_cc_im_xi_cc_pt_vs_eta_c4_2Hit", "pt selected", 200, 0, 20, 30, -1.5, 1.5}, "fPtMCXiCC", "fXiCCEta"); 
  auto h_df_xi_cc_im_xi_cc_bkg_pt_c4_2Hit = df_xi_cc_im_c4_2Hit.Filter(invMassXiccCut, {"fXiccMass"}, "c4_2Hit_Ivm").Histo1D({"df_xi_cc_im_xi_cc_pt_bkg_c4_2Hit", "pt selected", 200, 0, 20}, "fPtXiCC"); 
  auto out_counter_c4_2Hit = df_xi_cc_im_c4_2Hit.Filter(invMassXiccCut, {"fXiccMass"}).Count();
  
  auto df_xi_cc_im_c4_3Hit = df_xi_c
    .Filter("fXiHitsAdded > 2")
    .Filter("XicXiccDecayRadDiff > 0","c4_3Hit_XicXiccDecayRadDiff")
    .Filter("fXicDaughterDCA < 20","c4_3Hit_fXicDaughterDCA")
    .Filter("fXicDecayRadius > 0.004","c4_3Hit_fXicDecayRadius")
    .Filter("fXicInvDecayLengthToPV > 0.002","c4_3Hit_fXicInvDecayLengthToPVStra")
    .Filter("fXicInvDecayLengthToDV < 0.1","c4_3Hit_fXicInvDecayLengthToDVStra")
    .Filter("fXiccDaughterDCA < 20","c4_3Hit_fXiccDaughterDCA")
    .Filter("fXiccDecayRadius > 0.005","c4_3Hit_fXiccDecayRadius")
    .Filter("fXiccInvDecayLengthToPV > 0.004","c4_3Hit_fXiccInvDecayLengthToPVStra")
    .Filter("TMath::Abs(fXiDCAxyToPV) > 5","c4_3Hit_fXiDCAxyToPV")
    .Filter("TMath::Abs(fXiDCAzToPV) > 10","c4_3Hit_fXiDCAzToPV")
    .Filter("TMath::Abs(fXicDCAxyToPV) > 20","c4_3Hit_fXicDCAxyToPV")
    .Filter("TMath::Abs(fXicDCAzToPV) > 20","c4_3Hit_fXicDCAzToPV")
    .Filter("TMath::Abs(fPic1DCAxyToPV) > 10","c4_3Hit_fXicPionDCAxyToPV1")
    .Filter("TMath::Abs(fPic1DCAzToPV) > 15","c4_3Hit_fXicPionDCAzToPV1")
    .Filter("TMath::Abs(fPic2DCAxyToPV) > 10","c4_3Hit_fXicPionDCAxyToPV2")
    .Filter("TMath::Abs(fPic2DCAzToPV) > 20","c4_3Hit_fXicPionDCAzToPV2")
    .Filter("TMath::Abs(fPiccDCAxyToPV) > 20","c4_3Hit_fPicDCAxyToPV")
    .Filter("TMath::Abs(fPiccDCAzToPV) > 20","c4_3Hit_fPicDCAzToPV")
    .Filter("TMath::Abs(fXiccDCAxyToPV) < 50","c4_3Hit_fXiccDCAxyToPV")
    .Filter("TMath::Abs(fXiccDCAzToPV) < 50","c4_3Hit_fXiccDCAzToPV")    
    ;
  //Fill some final histograms  
  auto h_df_xi_cc_im_xi_cc_mass_c4_3Hit = df_xi_cc_im_c4_3Hit.Histo1D({"df_xi_cc_im_xi_cc_mass_c4_3Hit", "xi_cc inv mass", 700, 2.6, 4.6}, "fXiccMass"); 
  auto h_df_xi_cc_im_xi_cc_pt_c4_3Hit = df_xi_cc_im_c4_3Hit.Filter(invMassXiccCut, {"fXiccMass"}, "c4_3Hit_Ivm").Histo2D({"df_xi_cc_im_xi_cc_pt_vs_eta_c4_3Hit", "pt selected", 200, 0, 20, 30, -1.5, 1.5}, "fPtMCXiCC", "fXiCCEta"); 
  auto h_df_xi_cc_im_xi_cc_bkg_pt_c4_3Hit = df_xi_cc_im_c4_3Hit.Filter(invMassXiccCut, {"fXiccMass"}, "c4_3Hit_Ivm").Histo1D({"df_xi_cc_im_xi_cc_pt_bkg_c4_3Hit", "pt selected", 200, 0, 20}, "fPtXiCC"); 
  auto out_counter_c4_3Hit = df_xi_cc_im_c4_3Hit.Filter(invMassXiccCut, {"fXiccMass"}).Count();

  //Select the Xi_cc
  auto df_xi_cc_im_c5 = df_xi_c
    .Filter("XicXiccDecayRadDiff > 0","c5_XicXiccDecayRadDiff")
    .Filter("fXicDaughterDCA < 10","c5_fXicDaughterDCA")
    .Filter("fXicDecayRadius > 0.01","c5_fXicDecayRadius")
    .Filter("fXicInvDecayLengthToPV > 0.003","c5_fXicInvDecayLengthToPVStra")
    .Filter("fXicInvDecayLengthToDV > 0.002","c5_fXicInvDecayLengthToDVStra")
    .Filter("fXiccDaughterDCA < 5","c5_fXiccDaughterDCA")
    .Filter("fXiccDecayRadius > 0.005","c5_fXiccDecayRadius")
    .Filter("(fXiccInvDecayLengthToPV > 0.003) && (fXiccInvDecayLengthToPV < 0.06)","c5_fXiccInvDecayLengthToPVStra")
    .Filter("TMath::Abs(fXiDCAxyToPV) > 10","c5_fXiDCAxyToPV")
    .Filter("TMath::Abs(fXiDCAzToPV) > 10","c5_fXiDCAzToPV")
    .Filter("TMath::Abs(fXicDCAxyToPV) > 15","c5_fXicDCAxyToPV")
    .Filter("TMath::Abs(fXicDCAzToPV) > 15","c5_fXicDCAzToPV")
    .Filter("TMath::Abs(fPic1DCAxyToPV) > 10","c5_fXicPionDCAxyToPV1")
    .Filter("TMath::Abs(fPic1DCAzToPV) > 15","c5_fXicPionDCAzToPV1")
    .Filter("TMath::Abs(fPic2DCAxyToPV) > 10","c5_fXicPionDCAxyToPV2")
    .Filter("TMath::Abs(fPic2DCAzToPV) > 15","c5_fXicPionDCAzToPV2")
    .Filter("TMath::Abs(fPiccDCAxyToPV) > 20","c5_fPicDCAxyToPV")
    .Filter("TMath::Abs(fPiccDCAzToPV) > 20","c5_fPicDCAzToPV")
    .Filter("TMath::Abs(fXiccDCAxyToPV) < 10","c5_fXiccDCAxyToPV")
    .Filter("TMath::Abs(fXiccDCAzToPV) < 10","c5_fXiccDCAzToPV")    
    ;
  
  //Fill some final histograms  
  auto h_df_xi_cc_im_xi_cc_mass_c5 = df_xi_cc_im_c5.Histo1D({"df_xi_cc_im_xi_cc_mass_c5", "xi_cc inv mass", 700, 2.6, 4.6}, "fXiccMass"); 
  auto h_df_xi_cc_im_xi_cc_pt_c5 = df_xi_cc_im_c5.Filter(invMassXiccCut, {"fXiccMass"}, "c5_Ivm").Histo2D({"df_xi_cc_im_xi_cc_pt_vs_eta_c5", "pt selected", 200, 0, 20, 30, -1.5, 1.5}, "fPtMCXiCC", "fXiCCEta"); 
  auto h_df_xi_cc_im_xi_cc_bkg_pt_c5 = df_xi_cc_im_c5.Filter(invMassXiccCut, {"fXiccMass"}, "c5_Ivm").Histo1D({"df_xi_cc_im_xi_cc_pt_bkg_c5", "pt selected", 200, 0, 20}, "fPtXiCC");   
  auto out_counter_c5 = df_xi_cc_im_c5.Filter(invMassXiccCut, {"fXiccMass"}).Count(); 


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
  
  TString outName = TString::Format("outxiccSelector_%s.root",outAddon )  ; 
  TFile* out = TFile::Open(outName.Data(), "recreate");
  
  out->cd(); 
  //this is my start
  HarryPlotter::CheckAndStore(out, h_df_in_im_xi_cc_mass);
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
  HarryPlotter::CheckAndStore(out, h_df_lmb_xi_cc_mass); 

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

  HarryPlotter::CheckAndStore(out, h_df_xi_im_xi_ddist_dv); 
  
  HarryPlotter::CheckAndStore(out, h_df_xi_xi_cc_mass);

  HarryPlotter::CheckAndStore(out, h_df_lmb_qa_trad_diff_lmb_xi); 
  
  HarryPlotter::CheckAndStore(out, h_df_xi_sel_xi_mass);
  
  //to xi_c
  HarryPlotter::CheckAndStore(out, h_df_xi_im_xi_mass);
   HarryPlotter::CheckAndStore(out, h_df_xi_qa_xi_c_mass); 
  HarryPlotter::CheckAndStore(out, h_df_xi_qa_xi_c_pt); 
  
  HarryPlotter::CheckAndStore(out, h_df_xi_qa_trad_diff_xi_xi_c);
  HarryPlotter::CheckAndStore(out, h_df_xi_qa_xi_c_ddca); 
  HarryPlotter::CheckAndStore(out, h_df_xi_qa_xi_c_ddist_pv); 
  HarryPlotter::CheckAndStore(out, h_df_xi_qa_xi_c_ddist_dv);
  HarryPlotter::CheckAndStore(out, h_df_xi_qa_xi_c_trad); 

  HarryPlotter::CheckAndStore(out, h_df_xi_qa_xi_dca_xy);
  HarryPlotter::CheckAndStore(out, h_df_xi_qa_xi_dca_z);

  HarryPlotter::CheckAndStore(out, h_df_xi_qa_pi_one_pt);
  HarryPlotter::CheckAndStore(out, h_df_xi_qa_pi_two_pt);
    
  HarryPlotter::CheckAndStore(out, h_df_xi_qa_pi_one_dca_xy);
  HarryPlotter::CheckAndStore(out, h_df_xi_qa_pi_one_dca_z);
  HarryPlotter::CheckAndStore(out, h_df_xi_qa_pi_two_dca_xy);
  HarryPlotter::CheckAndStore(out, h_df_xi_qa_pi_two_dca_z);

  HarryPlotter::CheckAndStore(out, h_df_xi_qa_trad_xi_c_ddist_pv);
  HarryPlotter::CheckAndStore(out, h_df_xi_qa_trad_xi_c_ddist_dv);
  HarryPlotter::CheckAndStore(out, h_df_xi_qa_dl_pv_xi_c_ddist_dv);
  HarryPlotter::CheckAndStore(out, h_df_xi_qa_dl_pv_xi_c_trad);
  HarryPlotter::CheckAndStore(out, h_df_xi_qa_dl_dv_xi_c_ddist_dv);
  HarryPlotter::CheckAndStore(out, h_df_xi_qa_dl_dv_xi_c_trad);
  
  HarryPlotter::CheckAndStore(out, h_df_xi_c_im_xi_c_mass); 
  
  //to xi_cc
 
  HarryPlotter::CheckAndStore(out, h_df_xi_c_qa_trad_diff_xi_xi_c);
  HarryPlotter::CheckAndStore(out, h_df_xi_c_qa_xi_cc_ddca);
  HarryPlotter::CheckAndStore(out, h_df_xi_c_qa_xi_cc_ddist_pv);
  HarryPlotter::CheckAndStore(out, h_df_xi_c_qa_xi_cc_trad);
  HarryPlotter::CheckAndStore(out, h_df_xi_c_qa_xi_c_dca_xy);
  HarryPlotter::CheckAndStore(out, h_df_xi_c_qa_xi_c_dca_z);
  HarryPlotter::CheckAndStore(out, h_df_xi_c_qa_xi_cc_dca_xy);
  HarryPlotter::CheckAndStore(out, h_df_xi_c_qa_xi_cc_dca_z);    

  HarryPlotter::CheckAndStore(out,h_df_xi_c_qa_pi_dca_xy);
  HarryPlotter::CheckAndStore(out,h_df_xi_c_qa_pi_dca_z);
  
  HarryPlotter::CheckAndStore(out,h_df_xi_c_qa_trad_xi_cc_ddist_pv);
  HarryPlotter::CheckAndStore(out,h_df_xi_c_qa_dl_pv_xi_cc_trad);
  
  HarryPlotter::CheckAndStore(out,h_df_xi_c_xi_cc_mass); 
  
  HarryPlotter::CheckAndStore(out,h_df_identified); 
  auto h_df_xi_c_efficiency = h_df_identified->ProjectionX(TString::Format("EfficiencyNoCutting"), h_df_identified->GetYaxis()->FindBin(-1.5),h_df_identified->GetYaxis()->FindBin(+1.5));
  auto h_df_pT_Generated = h_gen_xi_cc_pt_eta_counter->ProjectionX("pTXiccGenerados",h_gen_xi_cc_pt_eta_counter->GetYaxis()->FindBin(-1.5),h_gen_xi_cc_pt_eta_counter->GetYaxis()->FindBin(+1.5));
  h_df_xi_c_efficiency->Sumw2(); 
  h_df_pT_Generated->Sumw2(); 
  h_df_pT_Generated->Rebin(4);
  h_df_xi_c_efficiency->Rebin(4);
  h_df_xi_c_efficiency->Divide(h_df_pT_Generated); 
  
  HarryPlotter::CheckAndStore(out, h_df_pT_Generated); 
  HarryPlotter::CheckAndStore(out, h_df_xi_c_efficiency); 

  //xi_cc selected
  HarryPlotter::CheckAndStore(out,h_df_xi_cc_qa_xi_cc_pt); 
  
  HarryPlotter::CheckAndStore(out,h_df_xi_cc_qa_pi_c1_pt); 
  auto h_df_xi_cc_qa_pi_c1_pt_1  = h_df_xi_cc_qa_pi_c1_pt->ProjectionX("h_df_xicc_qa_pi_c1_pt_1", h_df_xi_cc_qa_pi_c1_pt->GetYaxis()->FindBin(0.) ,h_df_xi_cc_qa_pi_c1_pt->GetYaxis()->FindBin(4.5));
  auto h_df_xi_cc_qa_pi_c1_pt_2  = h_df_xi_cc_qa_pi_c1_pt->ProjectionX("h_df_xicc_qa_pi_c1_pt_2", h_df_xi_cc_qa_pi_c1_pt->GetYaxis()->FindBin(4.5) ,h_df_xi_cc_qa_pi_c1_pt->GetYaxis()->FindBin(20));
  HarryPlotter::CheckAndStore(out,h_df_xi_cc_qa_pi_c1_pt_1); 
  HarryPlotter::CheckAndStore(out,h_df_xi_cc_qa_pi_c1_pt_2); 

  HarryPlotter::CheckAndStore(out,h_df_xi_cc_qa_pi_c2_pt); 
  auto h_df_xi_cc_qa_pi_c2_pt_1  = h_df_xi_cc_qa_pi_c2_pt->ProjectionX("h_df_xicc_qa_pi_c2_pt_1", h_df_xi_cc_qa_pi_c2_pt->GetYaxis()->FindBin(0.) ,h_df_xi_cc_qa_pi_c2_pt->GetYaxis()->FindBin(4.5));
  auto h_df_xi_cc_qa_pi_c2_pt_2  = h_df_xi_cc_qa_pi_c2_pt->ProjectionX("h_df_xicc_qa_pi_c2_pt_2", h_df_xi_cc_qa_pi_c2_pt->GetYaxis()->FindBin(4.5) ,h_df_xi_cc_qa_pi_c2_pt->GetYaxis()->FindBin(20));
  HarryPlotter::CheckAndStore(out,h_df_xi_cc_qa_pi_c2_pt_1); 
  HarryPlotter::CheckAndStore(out,h_df_xi_cc_qa_pi_c2_pt_2); 
  
  HarryPlotter::CheckAndStore(out,h_df_xi_cc_qa_pi_cc_pt); 
  auto h_df_xi_cc_qa_pi_cc_pt_1  = h_df_xi_cc_qa_pi_cc_pt->ProjectionX("h_df_xicc_qa_pi_cc_pt_1", h_df_xi_cc_qa_pi_cc_pt->GetYaxis()->FindBin(0.) ,h_df_xi_cc_qa_pi_cc_pt->GetYaxis()->FindBin(4.5));
  auto h_df_xi_cc_qa_pi_cc_pt_2  = h_df_xi_cc_qa_pi_cc_pt->ProjectionX("h_df_xicc_qa_pi_cc_pt_2", h_df_xi_cc_qa_pi_cc_pt->GetYaxis()->FindBin(4.5) ,h_df_xi_cc_qa_pi_cc_pt->GetYaxis()->FindBin(20));
  HarryPlotter::CheckAndStore(out,h_df_xi_cc_qa_pi_cc_pt_1); 
  HarryPlotter::CheckAndStore(out,h_df_xi_cc_qa_pi_cc_pt_2);
  
  HarryPlotter::CheckAndStore(out, h_df_xi_cc_im_xi_cc_mass_c1); 
  HarryPlotter::CheckAndStore(out, h_df_xi_cc_im_xi_cc_pt_c1); 
  HarryPlotter::CheckAndStore(out, h_df_xi_cc_im_xi_cc_bkg_pt_c1); 
  HarryPlotter::CheckAndStore(out, h_df_xi_cc_im_xi_cc_mass_c2); 
  HarryPlotter::CheckAndStore(out, h_df_xi_cc_im_xi_cc_pt_c2); 
  HarryPlotter::CheckAndStore(out, h_df_xi_cc_im_xi_cc_bkg_pt_c2); 
  HarryPlotter::CheckAndStore(out, h_df_xi_cc_im_xi_cc_mass_c3); 
  HarryPlotter::CheckAndStore(out, h_df_xi_cc_im_xi_cc_pt_c3); 
  HarryPlotter::CheckAndStore(out, h_df_xi_cc_im_xi_cc_bkg_pt_c3); 
  HarryPlotter::CheckAndStore(out, h_df_xi_cc_im_xi_cc_mass_c4); 
  HarryPlotter::CheckAndStore(out, h_df_xi_cc_im_xi_cc_pt_c4); 
  HarryPlotter::CheckAndStore(out, h_df_xi_cc_im_xi_cc_bkg_pt_c4); 
  HarryPlotter::CheckAndStore(out, h_df_xi_cc_im_xi_cc_mass_c4_2Hit); 
  HarryPlotter::CheckAndStore(out, h_df_xi_cc_im_xi_cc_pt_c4_2Hit); 
  HarryPlotter::CheckAndStore(out, h_df_xi_cc_im_xi_cc_bkg_pt_c4_2Hit); 
  HarryPlotter::CheckAndStore(out, h_df_xi_cc_im_xi_cc_mass_c4_3Hit); 
  HarryPlotter::CheckAndStore(out, h_df_xi_cc_im_xi_cc_pt_c4_3Hit); 
  HarryPlotter::CheckAndStore(out, h_df_xi_cc_im_xi_cc_bkg_pt_c4_3Hit); 
  HarryPlotter::CheckAndStore(out, h_df_xi_cc_im_xi_cc_mass_c5); 
  HarryPlotter::CheckAndStore(out, h_df_xi_cc_im_xi_cc_pt_c5); 
  HarryPlotter::CheckAndStore(out, h_df_xi_cc_im_xi_cc_bkg_pt_c5); 

  HarryPlotter::CheckAndStore(out, h_LongTracks); 
  
  HarryPlotter::CheckAndStore(out, h_cand_counter); 
  
  HarryPlotter::CheckAndStore(out, h_gen_xi_pt_eta_counter); 
  HarryPlotter::CheckAndStore(out, h_gen_xi_pt_y_counter); 
  
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

  HarryPlotter::CheckAndStore(out, h_df_lmb_qa_pos_tof_diff_inner);   
      
  HarryPlotter::CheckAndStore(out, h_df_in_qa_neg_dca_xy_vs_hits);   
  HarryPlotter::CheckAndStore(out, h_df_in_qa_neg_dca_xy_vs_lmb_dl_pv);   
  HarryPlotter::CheckAndStore(out, h_df_in_qa_neg_dca_xy_vs_neg_pt);   

  HarryPlotter::CheckAndStore(out, h_df_lmb_qa_neg_tof_diff_inner);   
  
  HarryPlotter::CheckAndStore(out, h_df_in_qa_neg_hits_vs_pos_hits);   
  HarryPlotter::CheckAndStore(out, h_df_in_qa_neg_dca_xy_vs_pos_dca_xy);   

  HarryPlotter::CheckAndStore(out, h_df_lmb_qa_bach_dca_xy_vs_hits);   
  HarryPlotter::CheckAndStore(out, h_df_lmb_qa_bach_dca_xy_vs_lmb_dl_pv);   
  HarryPlotter::CheckAndStore(out, h_df_lmb_qa_bach_dca_xy_vs_bach_pt);   

  HarryPlotter::CheckAndStore(out, h_df_lmb_qa_bach_tof_diff_inner);   
  
  HarryPlotter::CheckAndStore(out, h_df_xi_qa_pic1_tof_diff_inner);
  HarryPlotter::CheckAndStore(out, h_df_xi_qa_pic2_tof_diff_inner);
  
  HarryPlotter::CheckAndStore(out, h_df_xi_c_qa_picc_tof_diff_inner);
  
  out->Close(); 
  timer.Stop();
  timer.Print();
  return 0; 
}
