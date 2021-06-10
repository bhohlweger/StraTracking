#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TSystem.h"
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

std::vector<float> ptbins = {1.,6.};//HarryPlotter::Getptbins(); 
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
  
  const char* fileName = argv[1]; 
  const char* outAddon = (argv[2])?argv[2]:""; 
  
  HarryPlotter::StyleBox(); 
 
  ROOT::EnableImplicitMT(); // Tell ROOT you want to go parallel          
  TString filePath = TString::Format("%s", fileName); 
  
  auto h_cand_counter = new TH1D("df_xi_c_candCounter", "candCounter", 1, 0, 1); 
  auto h_gen_xi_c_counter = new TH1D("ptXicGen", "candCounter", 100, 0, 10); 
  auto h_gen_xi_cc_counter = new TH1D("ptXiccGen", "candCounter", 100, 0, 10); 
   
  TChain input("fTreeCandidates"); 
  int inputFiles = 0; 
  if (filePath.Contains(".root")) { 
    input.Add(filePath); 
    inputFiles++;
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
	  continue;
	}
	if (inFile->IsZombie()) { 
	  inFile->Close();
	  continue; 
	} 
	TH1D* evtCounter = (TH1D*)inFile->Get("hEventCounter"); 
	TH1D* ptXicGen = (TH1D*)inFile->Get("hXiCGeneratedPt"); 
	TH1D* ptXiccGen = (TH1D*)inFile->Get("hXiCCGeneratedPt"); 
	if (!evtCounter||!ptXicGen||!ptXiccGen) { 
	  continue; 
	}
	h_cand_counter->Add(evtCounter); 
	h_gen_xi_c_counter->Add(ptXicGen); 
	h_gen_xi_cc_counter->Add(ptXiccGen); 
	input.Add(inSubDirFile);
	inputFiles++;
	inFile->Close();
      }
    }
  }
  
  std::cout << "Added " << inputFiles << " files to the chain \n"; 
  /*
  auto Baryons = [](float PDGCode) { 

  }; 
  */
  ROOT::RDataFrame df(input);
  
  auto PDGLead = [](int motherPDGs) {return motherPDGs;}; 

  auto df_in_qa = df
    .Filter("fFirstCandidateXiCC")
    ; 
  auto in_counter = df_in_qa.Count(); 
  
  auto h_PDGCodeDCAxy = df_in_qa.Histo2D({"PDGCodeDCAxy", "PDGCodeDCAxy", 10000, 0, 10000, 600, -300, 300}, "fPiccMotherPDG", "fPicDCAxyToPVTopo"); 

  TString outName = TString::Format("outpionDCA_%s.root",outAddon )  ; 
  TFile* out = TFile::Open(outName.Data(), "recreate");
  
  out->cd(); 
  //this is my start

  HarryPlotter::CheckAndStore(out, h_cand_counter); 
  HarryPlotter::CheckAndStore(out, h_gen_xi_c_counter); 
  HarryPlotter::CheckAndStore(out, h_gen_xi_cc_counter); 
   
  HarryPlotter::CheckAndStore(out, h_PDGCodeDCAxy); 

  out->Close(); 
  return 0; 
}

