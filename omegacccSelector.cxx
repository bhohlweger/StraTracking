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
   1) take the omegas for kind of granted -> QA cuts
   -> additional hits + pT cut 
   -> transverse Radius (difference between truth and reco, causality) + Daughter DCA (?)
   -> invariant Mass cut (?)
	  
   2) Ultimately: Select the Omega_c from Omega_cc 
   -> QA Cuts: transverse radius, dca among daughters, dca to PV ... invariant mass
   -> Pion DCA distribution... 
	  
   3) Select the Omega_cc 
   -> Omega_c from Omega_cc Cuts: DCA to PV 
   What this needs to be able: 
   1) process signal + background files 

   Some random documentation: 

   Things to think about: 
   1) Do we compare Topo & Strangeness Tracking? 
   2) Switch between signal and background
*/ 

//TODO: Changes ranges of invariant mass plots!
int main(int argc, char **argv) {
  TStopwatch timer;
  timer.Start();
  
  const char* fileName = argv[1]; 
  const char* outAddon = (argv[2])?argv[2]:""; 
  
  int omegaccDec = argv[3]?atoi(argv[3]):0; 
  bool ExclusiveSignal = (omegaccDec == 0)?false:true;
  
  int pTbin = argv[4]?atoi(argv[4]):-1; 
  if (pTbin >= 0 && pTbin > ptbins.size()-1) { 
    std::cout << "Crashing cause the requested pT bin is out of range ( requested = " << pTbin << " , max. available = " << ptbins.size()-1 << " )\n";
    return -999; 
  }
  int wrongAssociationMode = argv[5]?atoi(argv[5]):0; 
  
  int noOmega = argv[6]?atoi(argv[6]):0; 
  bool ForceNoOmega = (noOmega==0)?false:true; 
  
  HarryPlotter::StyleBox(); 
 
  ROOT::EnableImplicitMT(); // Tell ROOT you want to go parallel          
  TString filePath = TString::Format("%s", fileName); 
  
  auto h_LongTracks = new TH1D("NLongTrack", "NLongTrack", 10000, 0, 10000); 
  
  auto h_cand_counter = new TH1D("df_omega_c_candCounter", "candCounter", 1, 0, 1); 

  auto h_gen_omega_pt_eta_counter = new TH2D("ptetaOmegaGen", "candCounter", 200, 0, 20, 30, -1.5, 1.5); 
  auto h_gen_omega_pt_y_counter = new TH2D("ptyOmegaGen", "candCounter", 200, 0, 20, 30, -1.5, 1.5); 

  auto h_gen_omega_c_counter = new TH1D("ptOmegacGen", "candCounter", 200, 0, 20); 
  auto h_gen_omega_c_pt_eta_counter = new TH2D("ptetaOmegacGen", "candCounter", 200, 0, 20, 30, -1.5, 1.5); 
  auto h_gen_omega_c_pt_y_counter = new TH2D("ptyOmegacGen", "candCounter", 200, 0, 20, 30, -1.5, 1.5); 
    
  auto h_gen_omega_cc_counter = new TH1D("ptOmegaccGen", "candCounter", 200, 0, 20); 
  auto h_gen_omega_cc_pt_eta_counter = new TH2D("ptetaOmegaccGen", "candCounter", 200, 0, 20, 30, -1.5, 1.5); 
  auto h_gen_omega_cc_pt_y_counter = new TH2D("ptyOmegaccGen", "candCounter", 200, 0, 20, 30, -1.5, 1.5); 
   
  auto h_gen_omega_ccc_counter = new TH1D("ptOmegacccGen", "candCounter", 200, 0, 20); 
  auto h_gen_omega_ccc_pt_eta_counter = new TH2D("ptetaOmegacccGen", "candCounter", 200, 0, 20, 30, -1.5, 1.5); 
  auto h_gen_omega_ccc_pt_y_counter = new TH2D("ptyOmegacccGen", "candCounter", 200, 0, 20, 30, -1.5, 1.5); 
   
  TChain input("fTreeCandidates"); 
  int inputFiles = 0; 
  int inputFailures = 0; 
  if (filePath.Contains(".root")) { 
    TFile *inFile = TFile::Open(filePath);
    TH1D* hLongTracks = (TH1D*)inFile->Get("hNLongTracks");
    TH1D* evtCounter = (TH1D*)inFile->Get("hEventCounter"); 

    TH2D* ptOmegapteta = (TH2D*)inFile->Get("hPtEtaGeneratedOmega"); 
    TH2D* ptOmegapty = (TH2D*)inFile->Get("hPtYGeneratedOmega"); 
    
    TH1D* ptOmegacGen = (TH1D*)inFile->Get("hOmegaCGeneratedPt");
    TH2D* ptOmegacpteta = (TH2D*)inFile->Get("hPtEtaGeneratedOmegaC"); 
    TH2D* ptOmegacpty = (TH2D*)inFile->Get("hPtYGeneratedOmegaC"); 
    
    TH1D* ptOmegaccGen = (TH1D*)inFile->Get("hOmegaCCGeneratedPt");
    TH2D* ptOmegaccpteta = (TH2D*)inFile->Get("hPtEtaGeneratedOmegaCC"); 
    TH2D* ptOmegaccpty = (TH2D*)inFile->Get("hPtYGeneratedOmegaCC"); 

    TH1D* ptOmegacccGen = (TH1D*)inFile->Get("hOmegaCCCGeneratedPt");
    TH2D* ptOmegacccpteta = (TH2D*)inFile->Get("hPtEtaGeneratedOmegaCCC"); 
    TH2D* ptOmegacccpty = (TH2D*)inFile->Get("hPtYGeneratedOmegaCCC"); 
    
    if (!hLongTracks||!evtCounter||
	!ptOmegapteta||!ptOmegapty||
	!ptOmegacGen||!ptOmegacpteta||!ptOmegacpty||
	!ptOmegaccGen||!ptOmegaccpteta||!ptOmegaccpty||
	!ptOmegacccGen||!ptOmegacccpteta||!ptOmegacccpty) { 
      inputFailures++; 
      std::cout << "Zis is vehry bad, ze generation histograms are missing, Guenther! No histogram, no chain! \n"; 
      
      if(!hLongTracks )std::cout << " Missing hLongTracks \n"; 
      if(!evtCounter  )std::cout << " Missing evtCounter  \n";
      if(!ptOmegapteta   )std::cout << " Missing ptOmegapteta   \n";
      if(!ptOmegapty     )std::cout << " Missing ptOmegapty     \n";
      if(!ptOmegacGen    )std::cout << " Missing ptOmegacGen    \n";
      if(!ptOmegacpteta  )std::cout << " Missing ptOmegacpteta  \n";
      if(!ptOmegacpty    )std::cout << " Missing ptOmegacpty    \n";
      if(!ptOmegaccGen   )std::cout << " Missing ptOmegaccGen   \n";
      if(!ptOmegaccpteta )std::cout << " Missing ptOmegaccpteta \n";
      if(!ptOmegaccpty   )std::cout << " Missing ptOmegaccpty   \n";
      if(!ptOmegacccGen   )std::cout << " Missing ptOmegacccGen   \n";
      if(!ptOmegacccpteta )std::cout << " Missing ptOmegacccpteta \n";
      if(!ptOmegacccpty   )std::cout << " Missing ptOmegacccpty   \n";
    } else { 
      h_LongTracks->Add(hLongTracks);
      
      h_cand_counter->Add(evtCounter); 
      
      h_gen_omega_pt_eta_counter->Add(ptOmegapteta); 
      h_gen_omega_pt_y_counter->Add(ptOmegapty); 
      
      h_gen_omega_c_counter->Add(ptOmegacGen); 	
      h_gen_omega_c_pt_eta_counter->Add(ptOmegacpteta); 
      h_gen_omega_c_pt_y_counter->Add(ptOmegacpty); 
      
      h_gen_omega_cc_counter->Add(ptOmegaccGen); 	
      h_gen_omega_cc_pt_eta_counter->Add(ptOmegaccpteta); 
      h_gen_omega_cc_pt_y_counter->Add(ptOmegaccpty); 
	
      h_gen_omega_ccc_counter->Add(ptOmegacccGen); 	
      h_gen_omega_ccc_pt_eta_counter->Add(ptOmegacccpteta); 
      h_gen_omega_ccc_pt_y_counter->Add(ptOmegacccpty); 

      input.Add(filePath); 
      inputFiles++;
    }
    
  } else { 
    std::cout << "Input seems not to be a file, trying to build a chain from sudirs...\n"; 
    // if file does not eomegast, loop over subdirs to find TFiles 
    TSystemDirectory dir("ZeDirectory", filePath);
    auto files = dir.GetListOfFiles();
    bool oneTimeError = true; 
    for (auto fileObj : *files)  {
      auto file = (TSystemFile*) fileObj;
      TString inSubDirFile = TString::Format("%s/%s/omegaccc.treeoutput.root", filePath.Data(), file->GetName()).Data(); 
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

	TH2D* ptOmegapteta = (TH2D*)inFile->Get("hPtEtaGeneratedOmega"); 
	TH2D* ptOmegapty = (TH2D*)inFile->Get("hPtYGeneratedOmega"); 
    
	TH1D* ptOmegacGen = (TH1D*)inFile->Get("hOmegaCGeneratedPt");
	TH2D* ptOmegacpteta = (TH2D*)inFile->Get("hPtEtaGeneratedOmegaC"); 
	TH2D* ptOmegacpty = (TH2D*)inFile->Get("hPtYGeneratedOmegaC"); 
    
	TH1D* ptOmegaccGen = (TH1D*)inFile->Get("hOmegaCCGeneratedPt");
	TH2D* ptOmegaccpteta = (TH2D*)inFile->Get("hPtEtaGeneratedOmegaCC"); 
	TH2D* ptOmegaccpty = (TH2D*)inFile->Get("hPtYGeneratedOmegaCC"); 
	
	TH1D* ptOmegacccGen = (TH1D*)inFile->Get("hOmegaCCCGeneratedPt");
	TH2D* ptOmegacccpteta = (TH2D*)inFile->Get("hPtEtaGeneratedOmegaCCC"); 
	TH2D* ptOmegacccpty = (TH2D*)inFile->Get("hPtYGeneratedOmegaCCC"); 
	
	if (!hLongTracks||!evtCounter||
	    !ptOmegapteta||!ptOmegapty||
	    !ptOmegacGen||!ptOmegacpteta||!ptOmegacpty||
	    !ptOmegaccGen||!ptOmegaccpteta||!ptOmegaccpty||
	    !ptOmegacccGen||!ptOmegacccpteta||!ptOmegacccpty) { 
	  if (oneTimeError) { 
	    if(!hLongTracks )std::cout << " Missing hLongTracks \n"; 
	    if(!evtCounter  )std::cout << " Missing evtCounter  \n";
	    if(!ptOmegapteta   )std::cout << " Missing ptOmegapteta   \n";
	    if(!ptOmegapty     )std::cout << " Missing ptOmegapty     \n";
	    if(!ptOmegacGen    )std::cout << " Missing ptOmegacGen    \n";
	    if(!ptOmegacpteta  )std::cout << " Missing ptOmegacpteta  \n";
	    if(!ptOmegacpty    )std::cout << " Missing ptOmegacpty    \n";
	    if(!ptOmegaccGen   )std::cout << " Missing ptOmegaccGen   \n";
	    if(!ptOmegaccpteta )std::cout << " Missing ptOmegaccpteta \n";
	    if(!ptOmegaccpty   )std::cout << " Missing ptOmegaccpty   \n";
	    if(!ptOmegacccGen   )std::cout << " Missing ptOmegacccGen   \n";
	    if(!ptOmegacccpteta )std::cout << " Missing ptOmegacccpteta \n";
	    if(!ptOmegacccpty   )std::cout << " Missing ptOmegacccpty   \n";
	    oneTimeError = false; 
	  }
	  inFile->Close(); 
	  inputFailures++; 
	  continue; 
	} else { 
	  h_LongTracks->Add(hLongTracks);
      
	  h_cand_counter->Add(evtCounter); 
      
	  h_gen_omega_pt_eta_counter->Add(ptOmegapteta); 
	  h_gen_omega_pt_y_counter->Add(ptOmegapty); 
      
	  h_gen_omega_c_counter->Add(ptOmegacGen); 	
	  h_gen_omega_c_pt_eta_counter->Add(ptOmegacpteta); 
	  h_gen_omega_c_pt_y_counter->Add(ptOmegacpty); 
      
	  h_gen_omega_cc_counter->Add(ptOmegaccGen); 	
	  h_gen_omega_cc_pt_eta_counter->Add(ptOmegaccpteta); 
	  h_gen_omega_cc_pt_y_counter->Add(ptOmegaccpty); 

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
  
  //Omega cuts
  float omegaMass = 1.672; 
  float invMassDiffOmega = 0.005; 
  int addedHitsMin = ForceNoOmega?0:1; 
  
  auto invMassOmegaCut = [&invMassDiffOmega, &omegaMass](float invMass) { return (TMath::Abs(invMass-omegaMass) < invMassDiffOmega); }; 
  auto hitsCut = [&addedHitsMin](int AddedHits) { return (AddedHits >= addedHitsMin); }; 
  
  auto decLengthOmega = [&omegaMass](float len, float mom){ return TMath::Abs(mom)>1e-4?len*omegaMass/mom:-999; }; 

  //Omegac cuts
  float omegacMass = 2.695; 
  float invMassDiffOmegac = 0.024; //8 MeV/c2 mass window 
  auto decLengthOmegac = [&omegacMass](float len, float mom){ return TMath::Abs(mom)>1e-4?len*omegacMass/mom:-999; }; 
  
  auto invMassOmegacCut = [&invMassDiffOmegac, &omegacMass](float invMass) { return (TMath::Abs(invMass-omegacMass) < invMassDiffOmegac); }; 
  
  //omegacc cuts
  float omegaccMass = 3.746;   
  float invMassDiffOmegacc = 0.024; //8 MeV/c2 mass window 
  auto decLengthOmegacc = [&omegaccMass](float len, float mom){ return TMath::Abs(mom)>1e-4?len*omegaccMass/mom:-999; }; 
  auto invMassOmegaccCut = [&invMassDiffOmegacc, &omegaccMass](float invMass) { return (TMath::Abs(invMass-omegaccMass) < invMassDiffOmegacc); }; 

  //omegaccc cuts
  float omegacccMass = 4.797;   
  float invMassDiffOmegaccc = 0.024; //8 MeV/c2 mass window 
  auto decLengthOmegaccc = [&omegacccMass](float len, float mom){ return TMath::Abs(mom)>1e-4?len*omegacccMass/mom:-999; }; 
  auto invMassOmegacccCut = [&invMassDiffOmegaccc, &omegacccMass](float invMass) { return (TMath::Abs(invMass-omegacccMass) < invMassDiffOmegaccc); }; 

  float omegacccpTmin = pTbin >=0?ptbins[pTbin]:0; 
  float omegacccpTmax = pTbin >=0?ptbins[pTbin+1]:20.0; 
  
  std::cout << "Omegaccc pT range set to " << omegacccpTmin << " to " << omegacccpTmax << std::endl; 
  
  auto pTCut = [&omegacccpTmin, &omegacccpTmax](float pT) { return (omegacccpTmin < pT)&&(pT < omegacccpTmax); }; 

  if (ExclusiveSignal) { 
    wrongAssociationMode = 0; 
  }

  auto associations = [&wrongAssociationMode](bool isOmegacc, bool isOmegac, bool piccUsed, int picUsed) { 
    bool out = true; 
    if (wrongAssociationMode == 1) { 
      //pions come from the correct decay chain but are not correctly associated
      out = (!isOmegacc && !isOmegac && piccUsed & (picUsed == 2)); 
    } else if (wrongAssociationMode == 2) { 
      //Omegac and its pions are good but the pion from the omegacc is replaced by a pythia pion
      out = (!isOmegacc && isOmegac && !piccUsed & (picUsed == 2)); 
    } else if (wrongAssociationMode == 3) { 
      //Omegac and Omegacc are wrong but the two pions from the omegac are used only the omegacc pion is replaced by a pythia pion
      out = (!isOmegacc && !isOmegac && !piccUsed & (picUsed == 2)); 
    } else if (wrongAssociationMode == 4) { 
      //Omegac and Omegacc are wrong but one omegac pion and the omegacc pion is replaced by a pythia pion
      out = (!isOmegacc && !isOmegac && !piccUsed & (picUsed == 1)); 
    } else if (wrongAssociationMode == 5) { 
      //Omegac and Omegacc are wrong but both omegac pions are replaced by pythia pions
      out = (!isOmegacc && !isOmegac && piccUsed & (picUsed == 0)); 
    } else if (wrongAssociationMode == 6) { 
      //All three pions are replaced by pythia pions 
      out = (!isOmegacc && !isOmegac && !piccUsed & (picUsed == 0)); 
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

  if (ForceNoOmega) { 
    std::cout << "Rejecting Omegas, make sure you know what you doing!\n"; 
  } else { 
    std::cout << "Utilizing the full beauty of strangeness hits, requested nHits = "<<  addedHitsMin << "\n";
  }
  
  auto h_df_identified = df
    .Histo2D({"df_pt_vs_eta_ident", "pt selected", 200, 0, 20, 30, -1.5, 1.5}, "fOmegacccPtMC", "fOmegacccEta"); 
  
  //TODO: REDO 
  auto df_precut = df
    .Define("fPosTOFDiffInner", "fPositiveInnerTOF20Signal-fPositiveInnerExpectedSignal")
    .Define("fNegTOFDiffInner", "fNegativeInnerTOF20Signal-fNegativeInnerExpectedSignal")
    .Define("fBachTOFDiffInner", "fBachelorInnerTOF20Signal-fBachelorInnerExpectedSignal")
    .Define("fPicTOFDiffInner", "fPicInnerTOF20Signal-fPicInnerExpectedSignal")
    .Define("fPiccTOFDiffInner", "fPiccInnerTOF20Signal-fPiccInnerExpectedSignal")
    .Define("fPicccTOFDiffInner", "fPicccInnerTOF20Signal-fPicccInnerExpectedSignal")

    .Filter("fV0DecayRadius > 0.475", "safeteyPrecuts_1")
    .Filter("fOmegaDecayRadius > 0.475", "safeteyPrecuts_2")
    .Filter("TMath::Abs(fPositiveDCAxy)> 60", "safeteyPrecuts_3")
    .Filter("TMath::Abs(fPositiveDCAxy)< 12000", "safeteyPrecuts_4")
    .Filter("TMath::Abs(fPositiveDCAz )> 60", "safeteyPrecuts_5")
    .Filter("TMath::Abs(fNegativeDCAxy)> 60", "safeteyPrecuts_6")
    .Filter("TMath::Abs(fNegativeDCAz )> 60", "safeteyPrecuts_7")
    .Filter("TMath::Abs(fBachelorDCAxy)> 40", "safeteyPrecuts_8")
    .Filter("TMath::Abs(fBachelorDCAz )> 40", "safeteyPrecuts_9")
    .Filter("TMath::Abs(fPicDCAxyToPV         )> 10", "safeteyPrecuts_10")
    .Filter("TMath::Abs(fPicDCAzToPV          )> 10", "safeteyPrecuts_11")
    .Filter("TMath::Abs(fPiccDCAxyToPV        )> 5", "safeteyPrecuts_12")
    .Filter("TMath::Abs(fPiccDCAzToPV         )> 5", "safeteyPrecuts_13")
    .Filter("TMath::Abs(fPicccDCAxyToPV       )> 5", "safeteyPrecuts_14")
    .Filter("TMath::Abs(fPicccDCAzToPV        )> 5", "safeteyPrecuts_15")
    .Filter("TMath::Abs(fLambdaMass-1.116) < 0.012", "safeteyPrecuts_16")
    .Filter("fOmegaPt > 0.4", "safeteyPrecuts_17")
    .Filter("fV0DauDCA   < 1500", "safeteyPrecuts_18")
    .Filter("TMath::Abs(fOmegaMass-1.672) < 0.012", "safeteyPrecuts_19")
    .Filter("fOmegaDauDCA < 1200", "safeteyPrecuts_20")
    .Filter("TMath::Abs(fPosTOFDiffInner) < 100", "safeteyPrecuts_21")
    .Filter("TMath::Abs(fNegTOFDiffInner) < 100", "safeteyPrecuts_22")
    .Filter("TMath::Abs(fBachTOFDiffInner) < 100", "safeteyPrecuts_23")
    .Filter("TMath::Abs(fPicTOFDiffInner) < 100", "safeteyPrecuts_24")
    .Filter("TMath::Abs(fPiccTOFDiffInner) < 100", "safeteyPrecuts_25")
    .Filter("TMath::Abs(fPicccTOFDiffInner) < 100", "safeteyPrecuts_26")
    .Filter("TMath::Abs(fOmegacccDCAxyToPV) < 200", "safeteyPrecuts_27")
    .Filter("TMath::Abs(fOmegacccDCAzToPV) < 200", "safeteyPrecuts_28")
    .Filter("fOmegacDauDCA < 200", "safeteyPrecuts_29") 
    .Filter("fOmegaccDauDCA < 150", "safeteyPrecuts_30") 
    .Filter("fOmegacccDauDCA < 125", "safeteyPrecuts_31") 
    .Filter("fOmegacDecayDistanceFromPV > 0.0050", "safeteyPrecuts_32") 
    .Filter("TMath::Abs(fOmegacMass-2.695) < 0.06", "safeteyPrecuts_33")
    .Filter("TMath::Abs(fOmegaccMass-3.746) < 0.06", "safeteyPrecuts_34")
    .Filter("TMath::Abs(fOmegacccMass-4.797) < 0.6", "safeteyPrecutsOut")
    ; 

  //auto df_ForceOmega = ForceNoOmega?df_precut.Filter("!fTrueOmega","noTrueOmegas"):df_precut.Filter("fTrueOmega||!fTrueOmega","TrueAndFalseOmegas");
  auto df_ForceOmega = df_precut.Filter("fTrueOmega||!fTrueOmega","TrueAndFalseOmegas");

  auto df_in = (
		ExclusiveSignal?
		df_ForceOmega
		.Filter("fTrueOmegaccc","trueOmegacccs")
		:df_ForceOmega
		.Filter("!fTrueOmegaccc","fakeOmegacccs")
		)
    .Define("fOmegacccPDGMass", [&omegacccMass]() {return omegacccMass;})
    .Define("fOmegacccY", HarryPlotter::YFromMomentum, {"fOmegacccP", "fOmegacccPt", "fOmegacccPDGMass", "fOmegacccEta"})
    .Filter("TMath::Abs(fOmegacccEta)<1.5", "OmegacccEta")
    .Filter(pTCut, {"fOmegacccPt"}, "pTOmegaccc")
    ;

  auto h_df_in_im_omega_ccc_mass = df_in.Histo1D({"h_df_in_im_omega_ccc_mass", "omega_ccc inv mass", 700, 3.6, 5.8}, "fOmegacccMass"); 
  auto in_counter = df_in.Count(); 

  auto df_in_qa = df_in.Filter("fFirstCandidateOmegaCCC","df_in_h_bool"); 

  //towards the Lambda for the Omega 
  auto h_df_in_qa_pos_pt = df_in_qa.Histo1D({"df_in_qa_pos_pt", "pos pt", 200, 0, 20}, "fPositivePt"); 
  auto h_df_in_qa_pos_dca_xy = df_in_qa.Histo1D({"df_in_qa_pos_dca_xy", "pos dca xy pv", 1000, -1e4, 1e4}, "fPositiveDCAxy"); 
  auto h_df_in_qa_pos_dca_z = df_in_qa.Histo1D({"df_in_qa_pos_dca_z", "pos dca z pv", 1000, -1e4, 1e4}, "fPositiveDCAz"); 
  auto h_df_in_qa_pos_dca_xy_wide = df_in_qa.Histo1D({"df_in_qa_pos_dca_xy_wide", "pos dca xy pv", 1000, -1e5, 1e5}, "fPositiveDCAxy"); 
  auto h_df_in_qa_pos_dca_z_wide = df_in_qa.Histo1D({"df_in_qa_pos_dca_z_wide", "pos dca z pv", 1000, -1e5, 1e5}, "fPositiveDCAz"); 
  auto h_df_in_qa_pos_hits = df_in_qa.Histo1D({"df_in_qa_pos_hits", "pos hits", 15, 0, 15}, "fPositiveClusters"); 
  auto h_df_in_qa_pos_chisq = df_in_qa.Histo1D({"df_in_qa_pos_chisq", "pos chisq", 200, 0, 200}, "fPositiveChisquare"); 
  auto h_df_in_qa_pos_chisqhits = df_in_qa.Define("fPositiveChisquareOverHits", "fPositiveChisquare/fPositiveClusters").Histo1D({"df_in_qa_pos_chisqhits", "pos chisq over hits", 100, 0, 50}, "fPositiveChisquareOverHits"); 

  auto h_df_in_qa_pos_dca_xy_vs_hits = df_in_qa.Histo2D<float,int>({"df_in_qa_pos_dca_xy_vs_hits", "pos dca xy pv vs hits", 1000, -10000, 10000, 15, -0.5, 14.5}, "fPositiveDCAxy", "fPositiveClusters"); 
  auto h_df_in_qa_pos_dca_xy_vs_lmb_dl_pv = df_in_qa.Define("fLmbInvDecayLengthToPV", decLengthLmb, {"fV0DecayLength", "fLambdaP"}).Histo2D<float,float>({"df_in_qa_pos_dca_xy_vs_lmb_dl_pv", "pos dca xy pv vs lmb decay length", 1000, -10000, 10000, 100, 0, 40}, "fPositiveDCAxy", "fLmbInvDecayLengthToPV"); 
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
  auto h_df_in_qa_neg_dca_xy_vs_lmb_dl_pv = df_in_qa.Define("fLmbInvDecayLengthToPV", decLengthLmb, {"fV0DecayLength", "fLambdaP"}).Histo2D<float,float>({"df_in_qa_neg_dca_xy_vs_lmb_dl_pv", "neg dca xy pv vs lmb decay length", 1000, -10000, 10000, 100, 0, 40}, "fNegativeDCAxy", "fLmbInvDecayLengthToPV"); 
  auto h_df_in_qa_neg_dca_xy_vs_neg_pt = df_in_qa.Histo2D<float,float>({"h_df_in_qa_neg_dca_xy_vs_neg_pt", "neg dca xy pv vs neg pt", 1000, -10000, 10000, 200, 0, 20}, "fNegativeDCAxy", "fNegativePt"); 
  auto h_df_in_qa_neg_hits_vs_pos_hits = df_in_qa.Histo2D<int,int>({"h_df_in_qa_neg_hits_vs_pos_hits", "neg hits vs pos hits", 15, -0.5, 14.5, 15, -0.5, 14.5}, "fNegativeClusters", "fPositiveClusters"); 
  auto h_df_in_qa_neg_dca_xy_vs_pos_dca_xy = df_in_qa.Histo2D<float,float>({"df_in_qa_neg_dca_xy_vs_pos_dca_xy", "neg dca xy pv vs pos dca xy pv", 1000, -10000, 10000, 1000, -10000, 10000}, "fNegativeDCAxy", "fPositiveDCAxy"); 

  auto h_df_in_qa_lmb_dca_xy = df_in_qa.Histo1D({"df_in_qa_lmb_dca_xy", "lmb dca xy pv", 1000, -1e4, 1e4}, "fV0DCAxyToPV");
  auto h_df_in_qa_lmb_dca_z = df_in_qa.Histo1D({"df_in_qa_lmb_dca_z", "lmb dca z pv", 1000, -1e4, 1e4}, "fV0DCAzToPV");
  auto h_df_in_qa_lmb_dca_xy_wide = df_in_qa.Histo1D({"df_in_qa_lmb_dca_xy_wide", "lmb dca xy pv", 1000, -1e5, 1e5}, "fV0DCAxyToPV");
  auto h_df_in_qa_lmb_dca_z_wide = df_in_qa.Histo1D({"df_in_qa_lmb_dca_z_wide", "lmb dca z pv", 1000, -1e5, 1e5}, "fV0DCAzToPV");
  
  auto h_df_in_qa_lmb_totp = df_in_qa.Histo1D({"df_in_qa_lmb_totp", "lmb tot p", 200, 0, 20}, "fLambdaP");   
  auto h_df_in_qa_lmb_ddist_pv = df_in_qa.Define("fLmbInvDecayLengthToPV", decLengthLmb, {"fV0DecayLength", "fLambdaP"}).Histo1D({"df_in_qa_lmb_ddist_pv", "lmb ddist to pv", 1000, 0, 40}, "fLmbInvDecayLengthToPV");
  auto h_df_in_qa_lmb_ddca = df_in_qa.Histo1D({"df_in_qa_lmb_ddca", "lmb prong dca", 500, 0, 100}, "fV0DauDCA"); 
  auto h_df_in_qa_lmb_ddca_wide = df_in_qa.Histo1D({"df_in_qa_lmb_ddca_wide", "lmb prong dca", 500, 0, 5000}, "fV0DauDCA"); 
  auto h_df_in_qa_lmb_trad = df_in_qa.Histo1D({"df_in_qa_lmb_trad", "lmb trad", 600, 0, 30}, "fV0DecayRadius"); 
  auto h_df_in_qa_lmb_trad_mc = df_in_qa.Histo1D({"df_in_qa_lmb_trad_mc", "lmb trad", 600, 0, 30}, "fV0DecayRadiusMC"); 
  auto h_df_in_qa_lmb_trad_diff = df_in_qa.Define("OmegaV0DecayRadDiff", "TMath::Abs(fV0DecayRadiusMC-fV0DecayRadius)").Histo1D({"df_in_qa_lmb_trad_diff", "lmb trad", 500, 0, 2}, "OmegaV0DecayRadDiff"); 
  
  auto h_df_in_qa_lmb_mass = df_in_qa.Histo1D({"df_in_qa_lmb_mass", "lmb inv mass", 750, 1., 1.4}, "fLambdaMass"); 
  auto h_df_in_qa_omega_mass = df_in_qa.Histo1D({"df_in_qa_omega_mass", "omega inv mass", 750, 1.4, 2.}, "fOmegaMass"); 
  
  //Define future variables and select Lambdas 

  auto df_lmb_im = df_in
    .Define("OmegaDecayRadDiff", "TMath::Abs(fOmegaDecayRadiusMC-fOmegaDecayRadius)")
    .Define("OmegaLmbDecayRadDiff", "fV0DecayRadius-fOmegaDecayRadius")
    .Define("OmegaOmegacDecayRadDiff", "fOmegaDecayRadius-fOmegacDecayRadius")
    .Define("fLmbInvDecayLengthToPV", decLengthLmb, {"fV0DecayLength", "fLambdaP"})
    .Define("fOmegaInvDecayLengthToPV", decLengthOmega, {"fOmegaDecayLength", "fOmegaP"})
    .Define("fOmegacInvDecayLengthToPV", decLengthOmegac, {"fOmegacDecayDistanceFromPV", "fOmegacP"})   // Omega_c decay length from the PV stra tracking 
    .Define("fOmegaccInvDecayLengthToPV", decLengthOmegacc, {"fOmegaccDecayDistanceFromPV", "fOmegaccP"})   // Omega_cc decay length from the PV stra tracking 
    .Define("fOmegacccInvDecayLengthToPV", decLengthOmegaccc, {"fOmegacccDecayDistanceFromPV", "fOmegacccP"}) // Omega_ccc decay length form the PV stra tracking = this is the Omega_ccc decay length! 

    .Define("fOmegaInvDecayLengthToDV", decLengthOmega, {"fOmegaCtoOmegaLength", "fOmegaP"}) //this is the Omega decay length 
    .Define("fOmegacInvDecayLengthToDV", decLengthOmegac, {"fOmegaCCtoOmegaCLength", "fOmegacP"})    //this is the Omega_c decay length 
    .Define("fOmegaccInvDecayLengthToDV", decLengthOmegacc, {"fOmegaCCCtoOmegaCCLength", "fOmegaccP"})    //this is the Omega_cc decay length 
    
    
    .Filter("TMath::Abs(fV0DCAxyToPV) < 5000", "fV0DCAxyToPV")
    .Filter("TMath::Abs(fV0DCAzToPV) < 7000", "fV0DCAzToPV")
    .Filter("fV0DauDCA < 200","fV0DauDCA")
    .Filter("fV0DecayRadius > 0.5","fV0DecayRadius")
    .Filter("fLmbInvDecayLengthToPV > 0.04","fLmbInvDecayLengthToPV")
    .Filter("TMath::Abs(fPositiveDCAxy) > 50","fPositiveDCAxy")
    .Filter("TMath::Abs(fPositiveDCAz) > 40","fPositiveDCAz")
    .Filter("TMath::Abs(fNegativeDCAxy) > 100","fNegativeDCAxy")
    .Filter("TMath::Abs(fNegativeDCAz) > 50","fNegativeDCAz")
    .Filter("TMath::Abs(fPosTOFDiffInner) < 50", "PosTOF")
    .Filter("TMath::Abs(fNegTOFDiffInner) < 50", "NegTOF")
    ;

  auto h_df_lmb_im_lmb_mass = df_lmb_im.Filter("fFirstCandidateOmegaCCC","df_lmb_im_h_bool").Histo1D({"df_lmb_im_lmb_mass", "lmb inv mass", 750, 1., 1.4}, "fLambdaMass"); 

  auto df_lmb =  df_lmb_im
    .Filter(invMassLmbCut, {"fLambdaMass"}, "fLambdaMass")
    ;

  //Towards the Omega for Omegac->Omega+2pi 
  //Fill some Histograms
  auto h_df_lmb_omega_ccc_mass = df_lmb.Histo1D({"h_df_lmb_omega_ccc_mass", "omega_cc inv mass", 700, 3.6, 5.8}, "fOmegacccMass"); 

  auto df_lmb_qa = df_lmb.Filter("fFirstCandidateOmegaCCC", "df_lmb_h_bool");
  
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
  
  auto h_df_lmb_qa_omega_pt = df_lmb_qa.Histo1D({"df_lmb_qa_omega_pt", "omega pt", 200, 0, 20}, "fOmegaPt"); 
  auto h_df_lmb_qa_omega_ddist_pv = df_lmb_qa.Histo1D({"df_lmb_qa_omega_ddist_pv", "omega ddist to pv", 1000, 0, 20}, "fOmegaInvDecayLengthToPV");
  auto h_df_lmb_qa_omega_ddca = df_lmb_qa.Histo1D({"df_lmb_qa_omega_ddca", "omega prong dca", 500, 0, 2000}, "fOmegaDauDCA"); 
  auto h_df_lmb_qa_omega_trad = df_lmb_qa.Histo1D({"df_lmb_qa_omega_trad", "omega trad", 250, 0, 10}, "fOmegaDecayRadius"); 
  auto h_df_lmb_qa_omega_trad_mc = df_lmb_qa.Histo1D({"df_lmb_qa_omega_trad_mc", "omega trad", 250, 0, 10}, "fOmegaDecayRadiusMC"); 
  auto h_df_lmb_qa_omega_trad_diff = df_lmb_qa.Histo1D({"df_lmb_qa_omega_trad_diff", "omega trad", 250, 0, 2}, "OmegaDecayRadDiff"); 

  auto h_df_lmb_qa_trad_diff_lmb_omega = df_lmb_qa.Histo1D({"df_lmb_qa_trad_diff_lmb_omega", "lmb-omega trad", 500, -100, 150}, "OmegaLmbDecayRadDiff"); 

  auto h_df_lmb_qa_omega_mass = df_lmb_qa.Histo1D({"df_lmb_qa_omega_mass", "omega inv mass", 750, 1.4, 2.}, "fOmegaMass"); 
  
  //Select Omegas excluding hits to avoid cheating 
  auto df_omega_sel = df_lmb
    .Filter("fOmegaDecayRadius > 0.5","fOmegaDecayRadius")
    .Filter("OmegaLmbDecayRadDiff > 0","OmegaLmbDecayRadDiff")
    .Filter("fOmegaDauDCA < 200","fOmegaDauDCA")
    .Filter("fOmegaInvDecayLengthToPV > 0.3","")
    .Filter("fOmegaDecayRadius > 0.5","")
    .Filter("TMath::Abs(fBachelorDCAxy) > 40","fBachelorDCAxy")
    .Filter("TMath::Abs(fBachelorDCAz) > 40","fBachelorDCAz")
    .Filter("TMath::Abs(fBachTOFDiffInner) < 50", "BachelorTOF")
    ;

  auto h_df_omega_sel_omega_mass = df_omega_sel.Filter("fFirstCandidateOmegaCCC","df_omega_sel_h_bool").Histo1D({"df_omega_sel_omega_mass", "omega inv mass", 750, 1.4, 2.}, "fOmegaMass"); 

  //Select Omegas including Hits
  auto df_omega_im = df_omega_sel.
    Filter(hitsCut, {"fOmegaHitsAdded"},"fOmegaHitsAdded")
    ; 
  
  auto h_df_omega_im_omega_mass = df_omega_im.Filter("fFirstCandidateOmegaCCC","df_omega_im_h_bool").Histo1D({"df_omega_im_omega_mass", "omega inv mass", 750, 1.4, 2.}, "fOmegaMass"); 
  
  auto h_df_omega_im_omega_ddist_dv = df_omega_im.Filter("fFirstCandidateOmegaCCC","df_omega_im_h_bool").Histo1D({"df_omega_im_omega_dist_dv", "omega decay dist", 1500, 0, 30}, "fOmegaInvDecayLengthToDV"); 

  auto df_omega = df_omega_im
    .Filter(invMassOmegaCut, {"fOmegaMass"},"fOmegaMass")
    ;
  auto h_df_omega_omega_ccc_mass = df_omega.Histo1D({"h_df_omega_omega_ccc_mass", "omega_cc inv mass", 700, 3.6, 5.8}, "fOmegacccMass"); 

  //Towards the Omega_c for Omega_cc->Omega_c+pi: 
  //Fill some histograms 
  
  auto df_omega_qa = df_omega.Filter("fFirstCandidateOmegaCCC","df_omega_h_bool");

  auto h_df_omega_qa_trad_diff_omega_omega_c = df_omega_qa.Histo1D({"df_omega_qa_trad_diff_omega_omega_c", "omega-omega_c trad", 500, -50, 200}, "OmegaOmegacDecayRadDiff") ;
  
  auto h_df_omega_qa_omega_c_mass = df_omega_qa.Histo1D({"df_omega_qa_omega_c_mass", "omega_c inv mass", 700, 1.9, 3.4}, "fOmegacMass"); 
  auto h_df_omega_qa_omega_c_pt = df_omega_qa.Histo1D({"df_omega_qa_omega_c_pt", "omega_c pt", 200, 0, 20}, "fOmegacP");  
  
  auto h_df_omega_qa_omega_c_ddca = df_omega_qa.Histo1D({"df_omega_qa_omega_c_ddca", "omega_c prong dca", 500, 0, 100}, "fOmegacDauDCA"); 
  auto h_df_omega_qa_omega_c_ddist_pv = df_omega_qa.Histo1D({"df_omega_qa_omega_c_dist_pv", "omega_c decay dist", 1500, 0, 0.30}, "fOmegacInvDecayLengthToPV"); //redefine as mL/p 
  auto h_df_omega_qa_omega_c_ddist_dv = df_omega_qa.Histo1D({"df_omega_qa_omega_c_dist_dv", "omega_c decay dist", 1500, 0, 0.30}, "fOmegacInvDecayLengthToDV"); 
  auto h_df_omega_qa_omega_c_trad = df_omega_qa.Histo1D({"df_omega_qa_omega_c_trad", "omega_c trad", 2000, 0, 0.4}, "fOmegacDecayRadius"); 
  
  auto h_df_omega_qa_omega_dca_xy = df_omega_qa.Histo1D({"df_omega_qa_omega_dca_xy", "omega dca xy", 1000, -500, 500}, "fOmegaDCAxyToPV");  
  auto h_df_omega_qa_omega_dca_z = df_omega_qa.Histo1D({"df_omega_qa_omega_dca_z", "omega dca z", 1000, -500, 500}, "fOmegaDCAzToPV");  
  
  //Pion (from the omegac) 
  
  auto h_df_omega_qa_pi_dca_xy = df_omega_qa.Histo1D({"df_omega_qa_pi_dca_xy", "pi dca xy stra", 1000, -500, 500}, "fPicDCAxyToPV");  
  auto h_df_omega_qa_pi_dca_z = df_omega_qa.Histo1D({"df_omega_qa_pi_dca_z", "pi dca z stra", 1000, -500, 500}, "fPicDCAzToPV");  
  
  auto h_df_omega_qa_pi_pt = df_omega_qa.Histo1D({"df_omega_qa_pi_pt", "pi c pt", 200, 0, 20}, "fPicPt");  

  auto h_df_omega_qa_pic_tof_diff_inner = df_omega_qa.Histo1D({"h_df_omega_qa_pic_tof_diff_inner", "beta expected vs measured", 2000, -2000, 2000}, "fPicTOFDiffInner");

  //study a bit decay length + trad 
  //cut on Trad and check decay lengths to pv and to dv 
  auto h_df_omega_qa_trad_omega_c_ddist_pv = df_omega_qa.Filter("fOmegacDecayRadius > 0.003","corrStudy_fOmegacDecayRadius").Histo1D({"df_omega_qa_trad_omega_c_dist_pv", "omega_c decay dist", 1500, 0, 0.30}, "fOmegacInvDecayLengthToPV"); 
  auto h_df_omega_qa_trad_omega_c_ddist_dv = df_omega_qa.Filter("fOmegacDecayRadius > 0.003","corrStudy_fOmegacDecayRadius").Histo1D({"df_omega_qa_trad_omega_c_dist_dv", "omega_c decay dist", 1500, 0, 0.30}, "fOmegacInvDecayLengthToDV"); 
  
  //cut on dl to pv 
  auto h_df_omega_qa_dl_pv_omega_c_ddist_dv =  df_omega_qa.Filter("fOmegacInvDecayLengthToPV > 0.002","corrStudy_fOmegacInvDecayLengthToPV").Histo1D({"df_omega_qa_dl_pv_omega_c_dist_dv", "omega_c decay dist", 1500, 0, 0.30}, "fOmegacInvDecayLengthToDV");
  auto h_df_omega_qa_dl_pv_omega_c_trad = df_omega_qa.Filter("fOmegacInvDecayLengthToPV > 0.002","corrStudy_fOmegacInvDecayLengthToPV").Histo1D({"df_omega_qa_dl_pv_omega_c_trad", "omega_c trad", 2000, 0, 0.4}, "fOmegacDecayRadius");
  
  //cut on dl to dv 
  auto h_df_omega_qa_dl_dv_omega_c_ddist_dv =  df_omega_qa.Filter("fOmegacInvDecayLengthToDV < 0.06","corrStudy_fOmegacInvDecayLengthToDV").Histo1D({"df_omega_qa_dl_dv_omega_c_dist_pv", "omega_c decay dist", 1500, 0, 0.30}, "fOmegacInvDecayLengthToPV");
  auto h_df_omega_qa_dl_dv_omega_c_trad = df_omega_qa.Filter("fOmegacInvDecayLengthToDV < 0.06","corrStudy_fOmegacInvDecayLengthToDV").Histo1D({"df_omega_qa_dl_dv_omega_c_trad", "omega_c trad", 2000, 0, 0.4}, "fOmegacDecayRadius");

  //Select the Omega_c 

  auto df_omega_c_im = df_omega
    .Define("OmegacOmegaccDecayRadDiff", "fOmegacDecayRadius-fOmegaccDecayRadius")
    .Filter("TMath::Abs(fPicTOFDiffInner) < 50", "PicTOF")
    .Filter("OmegaOmegacDecayRadDiff > 0","OmegaOmegacDecayRadDiff")
    ;
  
  auto h_df_omega_c_im_omega_c_mass = df_omega_c_im.Filter("fFirstCandidateOmegaCCC","df_omega_c_im_h_bool").Histo1D({"df_omega_c_im_omega_c_mass", "omega_c inv mass", 700, 1.9, 3.4}, "fOmegacMass"); 
  
  auto df_omega_c = df_omega_c_im.
    Filter(invMassOmegacCut, {"fOmegacMass"},"fOmegacMass");
  
  auto df_omega_c_qa = df_omega_c.Filter("fFirstCandidateOmegaCCC","df_omega_c_h_bool"); 

  //Towards the Omega_c for Omega_cc->Omega_c+pi: 
  //Fill some histograms 
  
  auto h_df_omega_c_qa_trad_diff_omega_omega_cc = df_omega_c_qa.Histo1D({"df_omega_c_qa_trad_diff_omega_c_omega_cc", "omega_c-omega_cc trad", 500, -50, 200}, "OmegacOmegaccDecayRadDiff") ;
  
  auto h_df_omega_c_qa_omega_cc_mass = df_omega_c_qa.Histo1D({"df_omega_c_qa_omega_cc_mass", "omega_cc inv mass", 700, 1.6, 3.2}, "fOmegaccMass"); 
  auto h_df_omega_c_qa_omega_cc_pt = df_omega_c_qa.Histo1D({"df_omega_c_qa_omega_cc_pt", "omega_cc pt", 200, 0, 20}, "fOmegaccP");  
  
  auto h_df_omega_c_qa_omega_cc_ddca = df_omega_c_qa.Histo1D({"df_omega_c_qa_omega_cc_ddca", "omega_cc prong dca", 500, 0, 100}, "fOmegaccDauDCA"); 
  auto h_df_omega_c_qa_omega_cc_ddist_pv = df_omega_c_qa.Histo1D({"df_omega_c_qa_omega_cc_dist_pv", "omega_cc decay dist", 1500, 0, 0.30}, "fOmegaccInvDecayLengthToPV"); //redefine as mL/p 
  auto h_df_omega_c_qa_omega_cc_ddist_dv = df_omega_c_qa.Histo1D({"df_omega_c_qa_omega_cc_dist_dv", "omega_cc decay dist", 1500, 0, 0.30}, "fOmegaccInvDecayLengthToDV"); 
  auto h_df_omega_c_qa_omega_cc_trad = df_omega_c_qa.Histo1D({"df_omega_c_qa_omega_cc_trad", "omega_cc trad", 2000, 0, 0.4}, "fOmegaccDecayRadius"); 
  
  auto h_df_omega_c_qa_omega_c_dca_xy = df_omega_c_qa.Histo1D({"df_omega_c_qa_omega_c_dca_xy", "omega dca xy", 1000, -500, 500}, "fOmegacDCAxyToPV");  
  auto h_df_omega_c_qa_omega_c_dca_z = df_omega_c_qa.Histo1D({"df_omega_c_qa_omega_c_dca_z", "omega dca z", 1000, -500, 500}, "fOmegacDCAzToPV");  
  
  //Pion (from the omegacc) 
  
  auto h_df_omega_c_qa_pi_dca_xy = df_omega_c_qa.Histo1D({"df_omega_c_qa_pi_dca_xy", "pi dca xy stra", 1000, -500, 500}, "fPiccDCAxyToPV");  
  auto h_df_omega_c_qa_pi_dca_z = df_omega_c_qa.Histo1D({"df_omega_c_qa_pi_dca_z", "pi dca z stra", 1000, -500, 500}, "fPiccDCAzToPV");  
  
  auto h_df_omega_c_qa_pi_pt = df_omega_c_qa.Histo1D({"df_omega_c_qa_pi_pt", "pi c pt", 200, 0, 20}, "fPiccPt");  

  auto h_df_omega_c_qa_picc_tof_diff_inner = df_omega_c_qa.Histo1D({"h_df_omega_c_qa_picc_tof_diff_inner", "beta expected vs measured", 2000, -2000, 2000}, "fPiccTOFDiffInner");

  //study a bit decay length + trad 
  //cut on Trad and check decay lengths to pv and to dv 
  auto h_df_omega_c_qa_trad_omega_cc_ddist_pv = df_omega_c_qa.Filter("fOmegaccDecayRadius > 0.003","corrStudy_fOmegaccDecayRadius").Histo1D({"df_omega_c_qa_trad_omega_cc_dist_pv", "omega_cc decay dist", 1500, 0, 0.30}, "fOmegaccInvDecayLengthToPV"); 
  auto h_df_omega_c_qa_trad_omega_cc_ddist_dv = df_omega_c_qa.Filter("fOmegaccDecayRadius > 0.003","corrStudy_fOmegaccDecayRadius").Histo1D({"df_omega_c_qa_trad_omega_cc_dist_dv", "omega_cc decay dist", 1500, 0, 0.30}, "fOmegaccInvDecayLengthToDV"); 
  
  //cut on dl to pv 
  auto h_df_omega_c_qa_dl_pv_omega_cc_ddist_dv =  df_omega_c_qa.Filter("fOmegaccInvDecayLengthToPV > 0.002","corrStudy_fOmegaccInvDecayLengthToPV").Histo1D({"df_omega_c_qa_dl_pv_omega_cc_dist_dv", "omega_cc decay dist", 1500, 0, 0.30}, "fOmegaccInvDecayLengthToDV");
  auto h_df_omega_c_qa_dl_pv_omega_cc_trad = df_omega_c_qa.Filter("fOmegaccInvDecayLengthToPV > 0.002","corrStudy_fOmegaccInvDecayLengthToPV").Histo1D({"df_omega_c_qa_dl_pv_omega_cc_trad", "omega_cc trad", 2000, 0, 0.4}, "fOmegaccDecayRadius");
  
  //cut on dl to dv 
  auto h_df_omega_c_qa_dl_dv_omega_cc_ddist_dv =  df_omega_c_qa.Filter("fOmegaccInvDecayLengthToDV < 0.06","corrStudy_fOmegaccInvDecayLengthToDV").Histo1D({"df_omega_c_qa_dl_dv_omega_cc_dist_pv", "omega_cc decay dist", 1500, 0, 0.30}, "fOmegaccInvDecayLengthToPV");
  auto h_df_omega_c_qa_dl_dv_omega_cc_trad = df_omega_c_qa.Filter("fOmegaccInvDecayLengthToDV < 0.06","corrStudy_fOmegaccInvDecayLengthToDV").Histo1D({"df_omega_c_qa_dl_dv_omega_cc_trad", "omega_cc trad", 2000, 0, 0.4}, "fOmegaccDecayRadius");

  //Select the omega_cc 

  auto df_omega_cc_im = df_omega_c
    .Define("OmegaccOmegacccDecayRadDiff", "fOmegaccDecayRadius-fOmegacccDecayRadius")
    .Filter("TMath::Abs(fPiccTOFDiffInner) < 50", "PiccTOF")
    .Filter("OmegacOmegaccDecayRadDiff > 0","OmegacOmegaccDecayRadDiff")
    ;
  
  auto h_df_omega_cc_im_omega_cc_mass = df_omega_cc_im.Filter("fFirstCandidateOmegaCCC","df_omega_cc_im_h_bool").Histo1D({"df_omega_cc_im_omega_cc_mass", "omega_cc inv mass", 700, 2.9, 4.5}, "fOmegaccMass"); 
  
  auto df_omega_cc = df_omega_cc_im.
    Filter(invMassOmegaccCut, {"fOmegaccMass"},"fOmegaccMass");
  
  auto df_omega_cc_qa = df_omega_cc.Filter("fFirstCandidateOmegaCCC","df_omega_cc_h_bool"); 

  /////////////

  //Towards the actual Omega_ccc
  //Fill some histograms 
  
  auto h_df_omega_cc_qa_trad_diff_omega_cc_omega_ccc = df_omega_cc_qa.Histo1D({"df_omega_cc_qa_trad_diff_omega_cc_omega_ccc", "omega_cc-omega_ccc trad", 500, -100, 150}, "OmegaccOmegacccDecayRadDiff") ;

  auto h_df_omega_cc_qa_omega_ccc_ddca = df_omega_cc_qa.Histo1D({"df_omega_cc_qa_omega_ccc_ddca", "omega_ccc prong dca", 500, 0, 100}, "fOmegacccDauDCA"); 
  auto h_df_omega_cc_qa_omega_ccc_ddist_pv = df_omega_cc_qa.Histo1D({"df_omega_cc_qa_omega_ccc_dist_pv", "omega_ccc decay dist", 3000, 0, 0.50}, "fOmegacccInvDecayLengthToPV"); 
  auto h_df_omega_cc_qa_omega_ccc_trad = df_omega_cc_qa.Histo1D({"df_omega_cc_qa_omega_ccc_trad", "omega_ccc trad", 2000, 0, 0.5}, "fOmegacccDecayRadius"); 
  
  auto h_df_omega_cc_qa_omega_cc_dca_xy = df_omega_cc_qa.Histo1D({"df_omega_cc_qa_omega_cc_dca_xy", "omega_cc dca xy stra", 1000, -500, 500}, "fOmegaccDCAxyToPV");  
  auto h_df_omega_cc_qa_omega_cc_dca_z  = df_omega_cc_qa.Histo1D({"df_omega_cc_qa_omega_cc_dca_z", "omega_cc dca z stra", 1000, -500, 500}, "fOmegaccDCAzToPV");  
  
  auto h_df_omega_cc_qa_omega_ccc_dca_xy = df_omega_cc_qa.Histo1D({"df_omega_cc_qa_omega_ccc_dca_xy", "omega_ccc dca xy stra", 1000, -500, 500}, "fOmegacccDCAxyToPV");  
  auto h_df_omega_cc_qa_omega_ccc_dca_z  = df_omega_cc_qa.Histo1D({"df_omega_cc_qa_omega_ccc_dca_z", "omega_ccc dca z stra", 1000, -500, 500}, "fOmegacccDCAzToPV");  

  //Pions (from the omegaccc)
 
  auto h_df_omega_cc_qa_pi_dca_xy = df_omega_cc_qa.Histo1D({"df_omega_cc_qa_pi_dca_xy", "omega_ccc dca xy stra", 1000, -500, 500}, "fPicccDCAxyToPV");  
  auto h_df_omega_cc_qa_pi_dca_z  = df_omega_cc_qa.Histo1D({"df_omega_cc_qa_pi_dca_z", "omega_ccc dca z stra", 1000, -500, 500}, "fPicccDCAzToPV");  

  auto h_df_omega_cc_qa_piccc_tof_diff_inner = df_omega_cc_qa.Histo1D({"h_df_omega_cc_qa_piccc_tof_diff_inner", "beta expected vs measured", 2000, -2000, 2000}, "fPicccTOFDiffInner");
  
  auto h_df_omega_cc_qa_omega_ccc_pt = df_omega_cc_qa
    .Filter("TMath::Abs(fPicDCAxyToPV) > 10")
    .Filter("TMath::Abs(fPicDCAzToPV) > 10")
    .Filter("TMath::Abs(fPiccDCAxyToPV) > 10")
    .Filter("TMath::Abs(fPiccDCAzToPV) > 10")
    .Filter("TMath::Abs(fPicccDCAxyToPV) > 10")
    .Filter("TMath::Abs(fPicccDCAzToPV) > 10")
    .Filter("TMath::Abs(fPicTOFDiffInner) < 50", "PicTOF")
    .Filter("TMath::Abs(fPiccTOFDiffInner) < 50", "PiccTOF")
    .Filter("TMath::Abs(fPicccTOFDiffInner) < 50", "PicccTOF")
    .Histo1D({"df_omega_cc_qa_omega_ccc_pt", "omega_ccc pt", 200, 0, 20}, "fOmegacccPt");  
	    
  auto h_df_omega_cc_qa_pi_pt = df_omega_cc_qa.Histo1D({"df_omega_cc_qa_pi_pt", "pi ccc pt", 200, 0, 20}, "fPicccPt");  

  //cut on Trad and check decay lengths 
  auto h_df_omega_cc_qa_trad_omega_ccc_ddist_pv = df_omega_cc_qa.Filter("fOmegacccDecayRadius > 0.003","corrStudy_fOmegacccDecayRadius").Histo1D({"df_omega_cc_qa_trad_omega_ccc_dist_pv", "omega_ccc decay dist", 3000, 0, 0.50}, "fOmegacccInvDecayLengthToPV"); 
  //cut on dl to pv 
  auto h_df_omega_cc_qa_dl_pv_omega_ccc_trad = df_omega_cc_qa.Filter("(fOmegacccInvDecayLengthToPV > 0.002) && (fOmegacccInvDecayLengthToPV <  0.06)","corrStudy_fOmegacccInvDecayLengthToPV").Histo1D({"df_omega_cc_qa_dl_pv_omega_ccc_trad", "omega_ccc trad", 2000, 0, 0.4}, "fOmegacccDecayRadius");
    
  auto h_df_omega_cc_omega_ccc_mass = df_omega_cc.Filter("OmegaccOmegacccDecayRadDiff > 0","OmegaccOmegacccDecayRadDiff").Histo1D({"df_omega_cc_omega_ccc_mass", "omega_ccc inv mass", 700, 3.6, 5.8}, "fOmegacccMass"); 
  
  //Select the Omega_cc
  auto df_omega_ccc_im_c1 = df_omega_cc
    .Filter("OmegaccOmegacccDecayRadDiff > 0","c1_OmegacOmegaccDecayRadDiffStra")
    .Filter("TMath::Abs(fOmegaDCAxyToPV) > 10","c1_fOmegaDCAxyToPV")
    .Filter("TMath::Abs(fOmegaDCAzToPV) > 10","c1_fOmegaDCAzToPV")
    
    .Filter("fOmegacDauDCA < 20","c1_fOmegacDaughterDCA")
    .Filter("fOmegacDecayRadius > 0.006","c1_fOmegacDecayRadius")
    .Filter("fOmegacInvDecayLengthToPV > 0.006","c1_fOmegacInvDecayLengthToPVStra")
    .Filter("TMath::Abs(fOmegacDCAxyToPV) > 10","c1_fOmegacDCAxyToPV")
    .Filter("TMath::Abs(fOmegacDCAzToPV) > 10","c1_fOmegacDCAzToPV")
    
    .Filter("fOmegaccDauDCA < 10","c1_fOmegaccDaughterDCA")
    .Filter("fOmegaccDecayRadius > 0.002","c1_fOmegaccDecayRadius")
    .Filter("fOmegaccInvDecayLengthToPV > 0.002","c1_fOmegaccInvDecayLengthToPVStra")

    .Filter("fOmegacccDauDCA < 10","c1_fOmegacccDaughterDCA")
    .Filter("fOmegacccDecayRadius > 0.002","c1_fOmegacccDecayRadius")
    .Filter("fOmegacccInvDecayLengthToPV > 0.002","c1_fOmegacccInvDecayLengthToPVStra")
    .Filter("TMath::Abs(fOmegacccDCAxyToPV) < 20","c1_fOmegaDCAxyToPV")
    .Filter("TMath::Abs(fOmegacccDCAzToPV) < 20","c1_fOmegaDCAzToPV")

    .Filter("TMath::Abs(fPicDCAxyToPV) > 10","c1_fPicDCAxyToPV")
    .Filter("TMath::Abs(fPicDCAzToPV) > 10","c1_fPicDCAzToPV")
    .Filter("TMath::Abs(fPiccDCAxyToPV) > 10","c1_fPiccDCAxyToPV")
    .Filter("TMath::Abs(fPiccDCAzToPV) > 10","c1_fPiccDCAzToPV")
    .Filter("TMath::Abs(fPicccDCAxyToPV) > 10","c1_fPicccDCAxyToPV")
    .Filter("TMath::Abs(fPicccDCAzToPV) > 10","c1_fPicccDCAzToPV")
    ;
  
  //Fill some final histograms  
  auto h_df_omega_ccc_im_omega_ccc_mass_c1 = df_omega_ccc_im_c1.Histo1D({"df_omega_cc_im_omega_cc_mass_c1", "omega_cc inv mass", 700, 3.6, 5.8}, "fOmegacccMass"); 
  auto h_df_omega_ccc_im_omega_ccc_pt_c1 = df_omega_ccc_im_c1
    .Filter(invMassOmegacccCut, {"fOmegacccMass"}, "c1_Ivm")
    .Histo2D({"df_omega_cc_im_omega_ccc_pt_vs_eta_c1", "pt selected", 200, 0, 20, 30, -1.5, 1.5}, "fOmegacccPtMC", "fOmegacccEta"); 
  
  auto out_counter_c1 = df_omega_ccc_im_c1.Filter(invMassOmegacccCut, {"fOmegacccMass"}).Count(); 


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
  
  TString outName = TString::Format("outomegacccSelector_%s.root",outAddon )  ; 
  TFile* out = TFile::Open(outName.Data(), "recreate");
  
  out->cd(); 
  //this is my start
  HarryPlotter::CheckAndStore(out, h_df_in_im_omega_ccc_mass);
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
  HarryPlotter::CheckAndStore(out, h_df_in_qa_omega_mass);
  
  //to omega 
  HarryPlotter::CheckAndStore(out, h_df_lmb_im_lmb_mass);
  HarryPlotter::CheckAndStore(out, h_df_lmb_omega_ccc_mass); 

  HarryPlotter::CheckAndStore(out, h_df_lmb_qa_bach_pt);
  HarryPlotter::CheckAndStore(out, h_df_lmb_qa_bach_dca_xy);
  HarryPlotter::CheckAndStore(out, h_df_lmb_qa_bach_dca_z);
  HarryPlotter::CheckAndStore(out, h_df_lmb_qa_bach_dca_xy_wide);
  HarryPlotter::CheckAndStore(out, h_df_lmb_qa_bach_dca_z_wide);
  HarryPlotter::CheckAndStore(out, h_df_lmb_qa_bach_hits);
  HarryPlotter::CheckAndStore(out, h_df_lmb_qa_bach_chisq);
  HarryPlotter::CheckAndStore(out, h_df_lmb_qa_bach_chisqhits);
  
  HarryPlotter::CheckAndStore(out, h_df_lmb_qa_omega_pt);
  HarryPlotter::CheckAndStore(out, h_df_lmb_qa_omega_ddist_pv);
  HarryPlotter::CheckAndStore(out, h_df_lmb_qa_omega_ddca);
  HarryPlotter::CheckAndStore(out, h_df_lmb_qa_omega_trad); 
  HarryPlotter::CheckAndStore(out, h_df_lmb_qa_omega_trad_mc);
  HarryPlotter::CheckAndStore(out, h_df_lmb_qa_omega_trad_diff);
  
  HarryPlotter::CheckAndStore(out, h_df_lmb_qa_omega_mass);

  HarryPlotter::CheckAndStore(out, h_df_omega_im_omega_ddist_dv); 
  
  HarryPlotter::CheckAndStore(out, h_df_omega_omega_ccc_mass);

  HarryPlotter::CheckAndStore(out, h_df_lmb_qa_trad_diff_lmb_omega); 
  
  HarryPlotter::CheckAndStore(out, h_df_omega_sel_omega_mass);
	    
  //to omega_c
  HarryPlotter::CheckAndStore(out, h_df_omega_im_omega_mass);
  HarryPlotter::CheckAndStore(out, h_df_omega_qa_omega_c_mass); 
  HarryPlotter::CheckAndStore(out, h_df_omega_qa_omega_c_pt); 
  
  HarryPlotter::CheckAndStore(out, h_df_omega_qa_trad_diff_omega_omega_c);
  HarryPlotter::CheckAndStore(out, h_df_omega_qa_omega_c_ddca); 
  HarryPlotter::CheckAndStore(out, h_df_omega_qa_omega_c_ddist_pv); 
  HarryPlotter::CheckAndStore(out, h_df_omega_qa_omega_c_ddist_dv);
  HarryPlotter::CheckAndStore(out, h_df_omega_qa_omega_c_trad); 
	    
  HarryPlotter::CheckAndStore(out, h_df_omega_qa_omega_dca_xy);
  HarryPlotter::CheckAndStore(out, h_df_omega_qa_omega_dca_z);

  HarryPlotter::CheckAndStore(out, h_df_omega_qa_pi_pt);
	    
  HarryPlotter::CheckAndStore(out, h_df_omega_qa_pi_dca_xy);
  HarryPlotter::CheckAndStore(out, h_df_omega_qa_pi_dca_z);

  HarryPlotter::CheckAndStore(out, h_df_omega_qa_trad_omega_c_ddist_pv);
  HarryPlotter::CheckAndStore(out, h_df_omega_qa_trad_omega_c_ddist_dv);
  HarryPlotter::CheckAndStore(out, h_df_omega_qa_dl_pv_omega_c_ddist_dv);
  HarryPlotter::CheckAndStore(out, h_df_omega_qa_dl_pv_omega_c_trad);
  HarryPlotter::CheckAndStore(out, h_df_omega_qa_dl_dv_omega_c_ddist_dv);
  HarryPlotter::CheckAndStore(out, h_df_omega_qa_dl_dv_omega_c_trad);
  
  HarryPlotter::CheckAndStore(out, h_df_omega_c_im_omega_c_mass); 
  

  //to omega_cc
  HarryPlotter::CheckAndStore(out, h_df_omega_c_qa_omega_cc_mass); 
  HarryPlotter::CheckAndStore(out, h_df_omega_c_qa_omega_cc_pt); 

  HarryPlotter::CheckAndStore(out, h_df_omega_c_qa_trad_diff_omega_omega_cc); 
  HarryPlotter::CheckAndStore(out, h_df_omega_c_qa_omega_cc_ddca); 
  HarryPlotter::CheckAndStore(out, h_df_omega_c_qa_omega_cc_ddist_pv); 
  HarryPlotter::CheckAndStore(out, h_df_omega_c_qa_omega_cc_ddist_dv);
  HarryPlotter::CheckAndStore(out, h_df_omega_c_qa_omega_cc_trad); 

  HarryPlotter::CheckAndStore(out, h_df_omega_c_qa_omega_c_dca_xy);
  HarryPlotter::CheckAndStore(out, h_df_omega_c_qa_omega_c_dca_z);

  HarryPlotter::CheckAndStore(out, h_df_omega_c_qa_pi_pt);

  HarryPlotter::CheckAndStore(out, h_df_omega_c_qa_pi_dca_xy);
  HarryPlotter::CheckAndStore(out, h_df_omega_c_qa_pi_dca_z);

  HarryPlotter::CheckAndStore(out, h_df_omega_c_qa_trad_omega_cc_ddist_pv);
  HarryPlotter::CheckAndStore(out, h_df_omega_c_qa_trad_omega_cc_ddist_dv);
  HarryPlotter::CheckAndStore(out, h_df_omega_c_qa_dl_pv_omega_cc_ddist_dv);
  HarryPlotter::CheckAndStore(out, h_df_omega_c_qa_dl_pv_omega_cc_trad);
  HarryPlotter::CheckAndStore(out, h_df_omega_c_qa_dl_dv_omega_cc_ddist_dv);
  HarryPlotter::CheckAndStore(out, h_df_omega_c_qa_dl_dv_omega_cc_trad);
  
  HarryPlotter::CheckAndStore(out, h_df_omega_cc_im_omega_cc_mass); 
  
  //to omega_cc
 
  HarryPlotter::CheckAndStore(out, h_df_omega_cc_qa_trad_diff_omega_cc_omega_ccc);
  HarryPlotter::CheckAndStore(out, h_df_omega_cc_qa_omega_ccc_ddca);
  HarryPlotter::CheckAndStore(out, h_df_omega_cc_qa_omega_ccc_ddist_pv);
  HarryPlotter::CheckAndStore(out, h_df_omega_cc_qa_omega_ccc_trad);
  HarryPlotter::CheckAndStore(out, h_df_omega_cc_qa_omega_cc_dca_xy);
  HarryPlotter::CheckAndStore(out, h_df_omega_cc_qa_omega_cc_dca_z);
  HarryPlotter::CheckAndStore(out, h_df_omega_cc_qa_omega_ccc_dca_xy);
  HarryPlotter::CheckAndStore(out, h_df_omega_cc_qa_omega_ccc_dca_z);
    
  HarryPlotter::CheckAndStore(out,h_df_omega_cc_qa_omega_ccc_pt); 
  HarryPlotter::CheckAndStore(out,h_df_omega_cc_qa_pi_pt); 

  HarryPlotter::CheckAndStore(out,h_df_omega_cc_qa_pi_dca_xy);
  HarryPlotter::CheckAndStore(out,h_df_omega_cc_qa_pi_dca_z);
  
  HarryPlotter::CheckAndStore(out,h_df_omega_cc_qa_trad_omega_ccc_ddist_pv);
  HarryPlotter::CheckAndStore(out,h_df_omega_cc_qa_dl_pv_omega_ccc_trad);
  
  HarryPlotter::CheckAndStore(out,h_df_omega_cc_omega_ccc_mass); 

  HarryPlotter::CheckAndStore(out,h_df_identified); 
  auto h_df_omega_c_efficiency = h_df_identified->ProjectionX(TString::Format("EfficiencyNoCutting"), 0, 30);
  auto h_df_pT_Generated = h_gen_omega_cc_pt_eta_counter->ProjectionX("pTOmegaccGenerados",h_gen_omega_cc_pt_eta_counter->GetYaxis()->FindBin(-1.5),h_gen_omega_cc_pt_eta_counter->GetYaxis()->FindBin(+1.5));

  h_df_pT_Generated->Sumw2(); 
  h_df_pT_Generated->Rebin(4);
  h_df_omega_c_efficiency->Sumw2(); 
  h_df_omega_c_efficiency->Rebin(4);
  h_df_omega_c_efficiency->Divide(h_df_pT_Generated); 
  HarryPlotter::CheckAndStore(out, h_df_pT_Generated); 
  HarryPlotter::CheckAndStore(out, h_df_omega_c_efficiency); 

  //omega_ccc selected
  HarryPlotter::CheckAndStore(out, h_df_omega_ccc_im_omega_ccc_mass_c1); 
  HarryPlotter::CheckAndStore(out, h_df_omega_ccc_im_omega_ccc_pt_c1); 
  HarryPlotter::CheckAndStore(out, h_LongTracks); 
  
  HarryPlotter::CheckAndStore(out, h_cand_counter); 
  
  HarryPlotter::CheckAndStore(out, h_gen_omega_pt_eta_counter); 
  HarryPlotter::CheckAndStore(out, h_gen_omega_pt_y_counter); 
  
  HarryPlotter::CheckAndStore(out, h_gen_omega_c_counter); 
  HarryPlotter::CheckAndStore(out, h_gen_omega_c_pt_eta_counter); 
  HarryPlotter::CheckAndStore(out, h_gen_omega_c_pt_y_counter); 
  
  HarryPlotter::CheckAndStore(out, h_gen_omega_cc_counter); 
  HarryPlotter::CheckAndStore(out, h_gen_omega_cc_pt_eta_counter); 
  HarryPlotter::CheckAndStore(out, h_gen_omega_cc_pt_y_counter); 
 
  HarryPlotter::CheckAndStore(out, h_gen_omega_ccc_counter); 
  HarryPlotter::CheckAndStore(out, h_gen_omega_ccc_pt_eta_counter); 
  HarryPlotter::CheckAndStore(out, h_gen_omega_ccc_pt_y_counter); 

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
  
  HarryPlotter::CheckAndStore(out, h_df_omega_qa_pic_tof_diff_inner);
  HarryPlotter::CheckAndStore(out, h_df_omega_c_qa_picc_tof_diff_inner);
  HarryPlotter::CheckAndStore(out, h_df_omega_cc_qa_piccc_tof_diff_inner);
  
  out->Close(); 
  timer.Stop();
  timer.Print();
  return 0; 
}
