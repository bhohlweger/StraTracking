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
  
  auto Baryons_u_d = [](int PDGCode) { 
    if((PDGCode > 1000) && (PDGCode < 3000)){
      return true; 
    } else { 
      return false; 
    }
  };
  
  auto Baryons_s = [](int PDGCode) { 
    if(PDGCode == 3122){ //Lambda
      return true; 
    } else if(PDGCode == 3222){ //Sigma+ 
      return true; 
    } else if(PDGCode == 3112){ //Sigma-
      return true; 
    } else if(PDGCode == 3322){ //Xi0
      return true; 
    } else if(PDGCode == 3312){ //Xi-
      return true; 
    } else if(PDGCode == 3334){ //Omega
      return true; 
    } else { 
      return false; 
    }
  };

  auto Baryons_c = [](int PDGCode) { 
    if(PDGCode == 4122){ //Lambda c
      return true; 
    } else if(PDGCode == 4132){ //Xi0 c
      return true; 
    } else if(PDGCode == 4232){ //Xi- c
      return true; 
    } else if(PDGCode == 4332){ //Omega c
      return true; 
    } else { 
      return false; 
    }
  };
  
  auto Baryons_b = [](int PDGCode) { 
    if(PDGCode == 5122){ //Lambda b 
      return true; 
    } else if(PDGCode == 5132){ //Xi- b
      return true; 
    } else if(PDGCode == 5232){ //Xi0 b
      return true; 
    } else { 
      return false; 
    }
  };
  
  auto Mesons_u_d = [](int PDGCode) { 
    if((PDGCode > 100) && (PDGCode < 300) && !PDGCode == 130){
      return true; 
    } else if (PDGCode == 333) { //Phi
      return true; 
    } else if (PDGCode == 130){ //K0L
      return false; 
    } else {
      return true; 
    }
  };
  
  auto Mesons_s = [](int PDGCode) { 
    if(PDGCode == 130){ //K0L
      return true; 
    } else if(PDGCode == 310){ //K0S
      return true; 
    } else if(PDGCode == 321){ //K+
      return true; 
    } else { 
      return false; 
    }
  };

  auto Mesons_c = [](int PDGCode) { 
    if(PDGCode == 411){ //D+
      return true; 
    } else if(PDGCode == 421){ //D0
      return true; 
    } else if(PDGCode == 431){ //Ds+
      return true; 
    } else { 
      return false; 
    }
  };

  auto Mesons_b = [](int PDGCode) { 
    if(PDGCode == 511){ //B0
      return true; 
    } else if(PDGCode == 521){ //B+
      return true; 
    } else if(PDGCode == 531){ //B0s
      return true; 
    } else if(PDGCode == 541){ //B+c
      return true; 
    } else { 
      return false; 
    }
  };
  
  auto nonOfTheAbove = [&Mesons_u_d, &Mesons_s, &Mesons_b, &Baryons_u_d, &Baryons_s, &Baryons_c, &Baryons_b] (int PDGCode) { 
    return !(Mesons_u_d(PDGCode)||&Mesons_s||Mesons_b||Baryons_u_d||Baryons_s||Baryons_c||Baryons_b); 
  }; 


  auto makeMeAbsolute = [] (int value) { return (int)TMath::Abs(value);}; 

  ROOT::RDataFrame df(input);
  
  auto PDGLead = [](int motherPDGs) {return motherPDGs;}; 

  auto df_in_qa = df
    .Filter("fFirstCandidateXiCC")
    .Define("absfPiccMotherPDG", makeMeAbsolute, {"fPiccMotherPDG"})
    ; 
  auto in_counter = df_in_qa.Count(); 
  
  auto h_PDGCode = df_in_qa.Histo1D({"PDGCode", "PDGCode", 20000, -10000, 10000}, "fPiccMotherPDG"); 
  auto h_PDGCodeWide = df_in_qa.Histo1D({"PDGCodeWide", "PDGCodeWide", 20000, -1e8, 1e8}, "fPiccMotherPDG"); 

  auto df_meson_ud = df_in_qa.Filter(Mesons_u_d, {"absfPiccMotherPDG"}); 
  auto dca_xy_meson_ud = df_meson_ud.Histo1D({"dca_xy_meson_ud", "dca_xy_meson_ud", 1000, -500, 500}, {"fPicDCAxyToPVTopo"}); 
  auto dca_z_meson_ud = df_meson_ud.Histo1D({"dca_z_meson_ud", "dca_z_meson_ud", 1000, -500, 500}, {"fPicDCAzToPVTopo"}); 

  auto df_meson_s = df_in_qa.Filter(Mesons_s, {"absfPiccMotherPDG"}); 
  auto dca_xy_meson_s = df_meson_s.Histo1D({"dca_xy_meson_s", "dca_xy_meson_s", 1000, -500, 500}, {"fPicDCAxyToPVTopo"}); 
  auto dca_z_meson_s = df_meson_s.Histo1D({"dca_z_meson_s", "dca_z_meson_s", 1000, -500, 500}, {"fPicDCAzToPVTopo"}); 

  auto df_meson_c = df_in_qa.Filter(Mesons_c, {"absfPiccMotherPDG"}); 
  auto dca_xy_meson_c = df_meson_c.Histo1D({"dca_xy_meson_c", "dca_xy_meson_c", 1000, -500, 500}, {"fPicDCAxyToPVTopo"}); 
  auto dca_z_meson_c = df_meson_c.Histo1D({"dca_z_meson_c", "dca_z_meson_c", 1000, -500, 500}, {"fPicDCAzToPVTopo"}); 

  auto df_meson_b = df_in_qa.Filter(Mesons_b, {"absfPiccMotherPDG"}); 
  auto dca_xy_meson_b = df_meson_b.Histo1D({"dca_xy_meson_b", "dca_xy_meson_b", 1000, -500, 500}, {"fPicDCAxyToPVTopo"}); 
  auto dca_z_meson_b = df_meson_b.Histo1D({"dca_z_meson_b", "dca_z_meson_b", 1000, -500, 500}, {"fPicDCAzToPVTopo"}); 
  
  auto df_baryon_ud = df_in_qa.Filter(Baryons_u_d, {"absfPiccMotherPDG"});
  auto dca_xy_baryon_ud = df_baryon_ud.Histo1D({"dca_xy_baryon_ud", "dca_xy_baryon_ud", 1000, -500, 500}, {"fPicDCAxyToPVTopo"}); 
  auto dca_z_baryon_ud = df_baryon_ud.Histo1D({"dca_z_baryon_ud", "dca_z_baryon_ud", 1000, -500, 500}, {"fPicDCAzToPVTopo"}); 

  auto df_baryon_s = df_in_qa.Filter(Baryons_s, {"absfPiccMotherPDG"}); 
  auto dca_xy_baryon_s = df_baryon_s.Histo1D({"dca_xy_baryon_s", "dca_xy_baryon_s", 1000, -500, 500}, {"fPicDCAxyToPVTopo"}); 
  auto dca_z_baryon_s = df_baryon_s.Histo1D({"dca_z_baryon_s", "dca_z_baryon_s", 1000, -500, 500}, {"fPicDCAzToPVTopo"}); 

  auto df_baryon_c = df_in_qa.Filter(Baryons_c, {"absfPiccMotherPDG"}); 
  auto dca_xy_baryon_c = df_baryon_c.Histo1D({"dca_xy_baryon_c", "dca_xy_baryon_c", 1000, -500, 500}, {"fPicDCAxyToPVTopo"}); 
  auto dca_z_baryon_c = df_baryon_c.Histo1D({"dca_z_baryon_c", "dca_z_baryon_c", 1000, -500, 500}, {"fPicDCAzToPVTopo"}); 

  auto df_baryon_b = df_in_qa.Filter(Baryons_b, {"absfPiccMotherPDG"}); 
  auto dca_xy_baryon_b = df_baryon_b.Histo1D({"dca_xy_baryon_b", "dca_xy_baryon_b", 1000, -500, 500}, {"fPicDCAxyToPVTopo"}); 
  auto dca_z_baryon_b = df_baryon_b.Histo1D({"dca_z_baryon_b", "dca_z_baryon_b", 1000, -500, 500}, {"fPicDCAzToPVTopo"}); 

  auto df_other = df_in_qa.Filter(nonOfTheAbove, {"absfPiccMotherPDG"}); 
  auto dca_xy_other = df_other.Histo1D({"dca_xy_other", "dca_xy_other", 1000, -500, 500}, {"fPicDCAxyToPVTopo"}); 
  auto dca_z_other = df_other.Histo1D({"dca_z_other", "dca_z_other", 1000, -500, 500}, {"fPicDCAzToPVTopo"}); 
  
  TString outName = TString::Format("outpionDCA_%s.root",outAddon )  ; 
  TFile* out = TFile::Open(outName.Data(), "recreate");
  
  out->cd(); 
  //this is my start

  HarryPlotter::CheckAndStore(out, h_cand_counter); 
  HarryPlotter::CheckAndStore(out, h_gen_xi_c_counter); 
  HarryPlotter::CheckAndStore(out, h_gen_xi_cc_counter); 
   
  HarryPlotter::CheckAndStore(out, h_PDGCode); 
  HarryPlotter::CheckAndStore(out, h_PDGCodeWide); 

  HarryPlotter::CheckAndStore(out, dca_xy_meson_ud);
  HarryPlotter::CheckAndStore(out, dca_z_meson_ud);

  HarryPlotter::CheckAndStore(out, dca_xy_meson_s);
  HarryPlotter::CheckAndStore(out, dca_z_meson_s);

  HarryPlotter::CheckAndStore(out, dca_xy_meson_c);
  HarryPlotter::CheckAndStore(out, dca_z_meson_c);

  HarryPlotter::CheckAndStore(out, dca_xy_meson_b);
  HarryPlotter::CheckAndStore(out, dca_z_meson_b);

  HarryPlotter::CheckAndStore(out, dca_xy_baryon_ud);
  HarryPlotter::CheckAndStore(out, dca_z_baryon_ud);

  HarryPlotter::CheckAndStore(out, dca_xy_baryon_s);
  HarryPlotter::CheckAndStore(out, dca_z_baryon_s);

  HarryPlotter::CheckAndStore(out, dca_xy_baryon_c);
  HarryPlotter::CheckAndStore(out, dca_z_baryon_c);

  HarryPlotter::CheckAndStore(out, dca_xy_baryon_b);
  HarryPlotter::CheckAndStore(out, dca_z_baryon_b);

  HarryPlotter::CheckAndStore(out, dca_xy_other);
  HarryPlotter::CheckAndStore(out, dca_z_other);

  out->Close(); 
  return 0; 
}

