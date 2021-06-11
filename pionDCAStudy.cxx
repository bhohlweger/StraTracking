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



template <typename T> bool Baryons_u_d (T PDGCode) { 
  if((PDGCode > 1000) && (PDGCode < 3000)){
    return true; 
  } else { 
    return false; 
  }
};


template <typename T>  bool Baryons_s(T PDGCode) { 
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

template <typename T>  bool Baryons_c(T PDGCode) { 
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
  
template <typename T>  bool Baryons_b(T PDGCode) { 
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
  
template <typename T>  bool Mesons_u_d(T PDGCode) { 
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
  
template <typename T>  bool Mesons_s(T PDGCode) { 
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

template <typename T>  bool Mesons_c(T PDGCode) { 
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

template <typename T>  bool Mesons_b(T PDGCode) { 
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
	/*
	  if (!evtCounter||!ptXicGen||!ptXiccGen) { 
	  continue; 
	  }
	  h_cand_counter->Add(evtCounter); 
	  h_gen_xi_c_counter->Add(ptXicGen); 
	  h_gen_xi_cc_counter->Add(ptXiccGen); 
	*/
	input.Add(inSubDirFile);
	inputFiles++;
	inFile->Close();
      }
    }
  }
  
  std::cout << "Added " << inputFiles << " files to the chain \n"; 
  
  auto anyWeakDecay = [] (int PDGCode) { 
    return (Mesons_s(PDGCode)||Mesons_c(PDGCode)||Mesons_b(PDGCode)||Baryons_s(PDGCode)||Baryons_c(PDGCode)||Baryons_b(PDGCode)); 
  };   
  auto nonWeakDecay = [] (int PDGCode) { 
    return (Mesons_u_d(PDGCode)||Baryons_u_d(PDGCode));
  }; 
  auto otherDecays = [] (int PDGCode) { 
    return !(Mesons_u_d(PDGCode)||Mesons_s(PDGCode)||Mesons_c(PDGCode)||Mesons_b(PDGCode)||Baryons_u_d(PDGCode)||Baryons_s(PDGCode)||Baryons_c(PDGCode)||Baryons_b(PDGCode)); 
  };
  auto findLeadingWeakInBaryonChain = [](ROOT::RVec<int> &PDGs) {
    auto weakPDGs = ROOT::VecOps::Filter(PDGs, [](int pdg){return Baryons_s(pdg)||Baryons_c(pdg)||Baryons_b(pdg);}); 
    auto sorted = ROOT::VecOps::Reverse(ROOT::VecOps::Sort(weakPDGs));
    return sorted.at(0); 
  }; 
  
  auto findLeadingWeakInMesonChain = [](ROOT::RVec<int> &PDGs) { 
    auto weakBaryonPDGs = ROOT::VecOps::Filter(PDGs, [](int pdg){return Baryons_s(pdg)||Baryons_c(pdg)||Baryons_b(pdg);}); 
    auto sortedBaryons = ROOT::VecOps::Reverse(ROOT::VecOps::Sort(weakBaryonPDGs));
    auto weakMesonPDGs = ROOT::VecOps::Filter(PDGs, [](int pdg){return Mesons_s(pdg)||Mesons_c(pdg)||Mesons_b(pdg);}); 
    auto sortedMesons = ROOT::VecOps::Reverse(ROOT::VecOps::Sort(weakMesonPDGs));
    if (sortedBaryons.size() < 1) { 
      return sortedMesons.at(0); 
    } else { 
      auto checkBaryons = sortedBaryons/1000; 
      auto checkMesons = sortedMesons/100; 
      if (checkBaryons.at(0) >= checkMesons.at(0)) { 
	return sortedBaryons.at(0); 
      } else { 
	return sortedMesons.at(0); 
      }
    }
  };
  
  auto bar_s = [](int pdg) { return Baryons_s(pdg);}; 
  auto bar_c = [](int pdg) { return Baryons_c(pdg);}; 
  auto bar_b = [](int pdg) { return Baryons_b(pdg);}; 
  
  auto mes_s = [](int pdg) { return Mesons_s(pdg);}; 
  auto mes_c = [](int pdg) { return Mesons_c(pdg);}; 
  auto mes_b = [](int pdg) { return Mesons_b(pdg);}; 

  

  //auto findWeakMesonInChain; 

  auto makeMeAbsolute = [] (int value) { return (int)TMath::Abs(value);}; 
  auto givemyintback = [] (int value) {return (int)value;}; 
  auto givemyfloatback = [] (float value) {return (float)value;}; 
  auto givemyRVecBack = [] (ROOT::RVec<int> array) { return array;};
  ROOT::RDataFrame df_in(input);
  
  auto df = df_in
    .Define("fPiccMotherPDG", givemyintback, {"fPionMotherPDG"})
    .Define("fPiccMotherNChain", givemyintback, {"fPionMotherNChain"})
    .Define("fPiccMotherChain", givemyRVecBack, {"fPionMotherChain"})
    .Define("fPicDCAxyToPVTopo", givemyfloatback, {"fPionDCAxy"})
    .Define("fPicDCAzToPVTopo", givemyfloatback, {"fPionDCAz"})
    ;

  auto FromBaryon = [](int pdgCodes) { 
    bool out = false; 
    auto checkDigits = pdgCodes/(int)1000; 
    if ((checkDigits > 0) && (checkDigits < 10)) {
      out = true;
    }
    return out; 
  };
 
  auto FromMeson = [](int pdgCodes) { 
    bool out = false; 
    auto checkDigits = pdgCodes/(int)100; 
    if ((checkDigits > 0) && (checkDigits < 10)) {
      out = true;
    }
    return out; 
  };
  
  //let's tackle this someway different ... 
  
  //Filter primary 
  auto df_prim = df.Filter("fPiccMotherNChain < 1"); 
  auto h_PDGCode_prim = df_prim.Histo1D({"PDGCodePrim", "PDGCode", 20000, -10000, 10000}, "fPiccMotherPDG"); 
  auto h_DCA_prim = df_prim.Histo1D({"DCAPrim", "DCA", 300, -1500, 1500}, "fPicDCAxyToPVTopo"); 

  auto df_strong = df.Filter(nonWeakDecay, {"fPiccMotherPDG"}); 
  auto h_PDGCode_strong = df_strong.Histo1D({"PDGCodeStrong", "PDGCode", 20000, -10000, 10000}, "fPiccMotherPDG"); 
  auto h_DCA_strong = df_strong.Histo1D({"DCAStrong", "DCA", 300, -1500, 1500}, "fPicDCAxyToPVTopo"); 

  auto df_other = df.Filter(otherDecays, {"fPiccMotherPDG"}); 
  auto h_PDGCode_other = df_other.Histo1D({"PDGCodeOther", "PDGCode", 20000, -10000, 10000}, "fPiccMotherPDG"); 
  auto h_DCA_other = df_other.Histo1D({"DCAOther", "DCA", 300, -1500, 1500}, "fPicDCAxyToPVTopo"); 

  //Filter non-primary 
  auto df_dec = df.Filter("fPiccMotherNChain > 0").Filter(anyWeakDecay, {"fPiccMotherPDG"});   
  //From Baryons 
  auto df_baryon = df_dec.Filter(FromBaryon, {"fPiccMotherPDG"}).Define("LeadingBaryon", findLeadingWeakInBaryonChain, {"fPiccMotherChain"});
  auto h_PDGCode_baryon = df_baryon.Histo1D({"PDGCodeBaryon", "PDGCode", 20000, -10000, 10000}, "fPiccMotherPDG"); 
  auto h_PDGCode_lead_baryon = df_baryon.Histo1D({"PDGCodeLeadBaryon", "PDGCode", 20000, -10000, 10000}, "LeadingBaryon"); 
  auto h_DCA_baryon = df_baryon.Histo1D({"DCABaryon", "DCA", 300, -1500, 1500}, "fPicDCAxyToPVTopo"); 
  auto h_DCA_baryon_s = df_baryon.Filter(bar_s, {"fPiccMotherPDG"}).Histo1D({"DCABaryon_s", "DCABaryon_s", 300, -1500, 1500}, "fPicDCAxyToPVTopo"); 
  auto h_DCA_baryon_c = df_baryon.Filter(bar_c, {"fPiccMotherPDG"}).Histo1D({"DCABaryon_c", "DCABaryon_c", 300, -1500, 1500}, "fPicDCAxyToPVTopo"); 
  auto h_DCA_baryon_b = df_baryon.Filter(bar_b, {"fPiccMotherPDG"}).Histo1D({"DCABaryon_b", "DCABaryon_d", 300, -1500, 1500}, "fPicDCAxyToPVTopo"); 

  auto df_baryon_prompt = df_baryon.Filter("fPiccMotherPDG==LeadingBaryon"); 
  auto h_PDGCode_baryon_prompt = df_baryon_prompt.Histo1D({"PDGCodeBaryonPrompt", "PDGCode", 20000, -10000, 10000}, "fPiccMotherPDG"); 
  auto h_DCA_baryon_prompt = df_baryon_prompt.Histo1D({"DCABaryon_Prompt", "DCABaryon_Prompt", 300, -1500, 1500}, "fPicDCAxyToPVTopo"); 
  auto h_DCA_baryon_prompt_s = df_baryon_prompt.Filter(bar_s, {"fPiccMotherPDG"}).Histo1D({"DCABaryonPrompt_s", "DCABaryonPrompt_s", 300, -1500, 1500}, "fPicDCAxyToPVTopo"); 
  auto h_DCA_baryon_prompt_c = df_baryon_prompt.Filter(bar_c, {"fPiccMotherPDG"}).Histo1D({"DCABaryonPrompt_c", "DCABaryonPrompt_c", 300, -1500, 1500}, "fPicDCAxyToPVTopo"); 
  auto h_DCA_baryon_prompt_b = df_baryon_prompt.Filter(bar_b, {"fPiccMotherPDG"}).Histo1D({"DCABaryonPrompt_b", "DCABaryonPrompt_b", 300, -1500, 1500}, "fPicDCAxyToPVTopo"); 

  auto df_baryon_non_prompt = df_baryon.Filter("fPiccMotherPDG!=LeadingBaryon"); 
  auto h_PDGCode_baryon_non_prompt = df_baryon_non_prompt.Histo1D({"PDGCodeBaryonNonPrompt", "PDGCode", 20000, -10000, 10000}, "LeadingBaryon"); 
  auto h_DCA_baryon_non_prompt = df_baryon_non_prompt.Histo1D({"DCABaryonNonPrompt", "DCABaryonNonPrompt", 300, -1500, 1500}, "fPicDCAxyToPVTopo"); 
  auto h_DCA_baryon_non_prompt_s = df_baryon_prompt.Filter(bar_s, {"fPiccMotherPDG"}).Histo1D({"DCABaryonNonPrompt_s", "DCABaryonNonPrompt_s", 300, -1500, 1500}, "fPicDCAxyToPVTopo"); 
  auto h_DCA_baryon_non_prompt_c = df_baryon_prompt.Filter(bar_c, {"fPiccMotherPDG"}).Histo1D({"DCABaryonNonPrompt_c", "DCABaryonNonPrompt_c", 300, -1500, 1500}, "fPicDCAxyToPVTopo"); 
  auto h_DCA_baryon_non_prompt_b = df_baryon_prompt.Filter(bar_b, {"fPiccMotherPDG"}).Histo1D({"DCABaryonNonPrompt_b", "DCABaryonNonPrompt_b", 300, -1500, 1500}, "fPicDCAxyToPVTopo"); 

  //Mesons 
  auto df_meson = df_dec.Filter(FromMeson, {"fPiccMotherPDG"}).Define("LeadingHadron", findLeadingWeakInMesonChain, {"fPiccMotherChain"});
  auto h_PDGCode = df_meson.Histo1D({"PDGCode", "PDGCode", 20000, -10000, 10000}, "fPiccMotherPDG"); 
  auto h_PDGCode_meson = df_meson.Histo1D({"PDGCodeMeson", "PDGCode", 20000, -10000, 10000}, "fPiccMotherPDG"); 
  auto h_DCA_meson = df_meson.Histo1D({"DCAMeson", "DCAMeson", 300, -1500, 1500}, "fPicDCAxyToPVTopo"); 
  auto h_DCA_meson_s = df_meson.Filter(mes_s, {"fPiccMotherPDG"}).Histo1D({"DCAMeson_s", "DCAMeson_s", 300, -1500, 1500}, "fPicDCAxyToPVTopo"); 
  auto h_DCA_meson_c = df_meson.Filter(mes_c, {"fPiccMotherPDG"}).Histo1D({"DCAMeson_c", "DCAMeson_c", 300, -1500, 1500}, "fPicDCAxyToPVTopo"); 
  auto h_DCA_meson_b = df_meson.Filter(mes_b, {"fPiccMotherPDG"}).Histo1D({"DCAMeson_b", "DCAMeson_b", 300, -1500, 1500}, "fPicDCAxyToPVTopo"); 

  auto df_meson_prompt = df_meson.Filter("fPiccMotherPDG==LeadingHadron"); 
  auto h_PDGCode_meson_prompt = df_meson_prompt.Histo1D({"PDGCodeMesonPrompt", "PDGCode", 20000, -10000, 10000}, "fPiccMotherPDG"); 
  auto h_DCA_meson_prompt = df_meson_prompt.Histo1D({"DCAMesonPrompt", "DCAMesonPrompt", 300, -1500, 1500}, "fPicDCAxyToPVTopo"); 
  auto h_DCA_meson_prompt_s = df_meson_prompt.Filter(bar_s, {"fPiccMotherPDG"}).Histo1D({"DCAMesonPrompt_s", "DCAMesonPrompt_s", 300, -1500, 1500}, "fPicDCAxyToPVTopo"); 
  auto h_DCA_meson_prompt_c = df_meson_prompt.Filter(bar_c, {"fPiccMotherPDG"}).Histo1D({"DCAMesonPrompt_c", "DCAMesonPrompt_c", 300, -1500, 1500}, "fPicDCAxyToPVTopo"); 
  auto h_DCA_meson_prompt_b = df_meson_prompt.Filter(bar_b, {"fPiccMotherPDG"}).Histo1D({"DCAMesonPrompt_b", "DCAMesonPrompt_b", 300, -1500, 1500}, "fPicDCAxyToPVTopo"); 

  
  auto df_meson_non_prompt = df_meson.Filter("fPiccMotherPDG!=LeadingHadron"); 
  auto h_PDGCode_meson_non_prompt = df_meson_non_prompt.Histo1D({"PDGCodeMesonNonPrompt", "PDGCode", 20000, -10000, 10000}, "LeadingHadron"); 
  auto h_DCA_meson_non_prompt = df_meson_non_prompt.Histo1D({"DCAMesonNonPrompt", "DCAMesonNonPrompt", 300, -1500, 1500}, "fPicDCAxyToPVTopo"); 
  auto h_DCA_meson_non_prompt_s = df_meson_prompt.Filter(bar_s, {"fPiccMotherPDG"}).Histo1D({"DCAMesonNonPrompt_s", "DCAMesonNonPrompt_s", 300, -1500, 1500}, "fPicDCAxyToPVTopo"); 
  auto h_DCA_meson_non_prompt_c = df_meson_prompt.Filter(bar_c, {"fPiccMotherPDG"}).Histo1D({"DCAMesonNonPrompt_c", "DCAMesonNonPrompt_c", 300, -1500, 1500}, "fPicDCAxyToPVTopo"); 
  auto h_DCA_meson_non_prompt_b = df_meson_prompt.Filter(bar_b, {"fPiccMotherPDG"}).Histo1D({"DCAMesonNonPrompt_b", "DCAMesonNonPrompt_b", 300, -1500, 1500}, "fPicDCAxyToPVTopo"); 
  
  auto df_meson_non_prompt_baryon = df_meson_non_prompt.Filter(FromBaryon, {"LeadingHadron"}); 
  auto h_PDGCode_meson_non_prompt_baryon = df_meson_non_prompt_baryon.Histo1D({"PDGCodeMesonNonPromptBaryon", "PDGCode", 20000, -10000, 10000}, "LeadingHadron"); 
  auto h_DCA_meson_non_prompt_baryon = df_meson_non_prompt_baryon.Histo1D({"DCAMesonNonPromptBaryon", "DCAMesonNonPromptBaryon", 300, -1500, 1500}, "fPicDCAxyToPVTopo"); 

  
  auto df_meson_non_prompt_meson = df_meson_non_prompt.Filter(FromMeson, {"LeadingHadron"}); 
  auto h_PDGCode_meson_non_prompt_meson = df_meson_non_prompt_meson.Histo1D({"PDGCodeMesonNonPromptMeson", "PDGCode", 20000, -10000, 10000}, "LeadingHadron"); 
  auto h_DCA_meson_non_prompt_meson = df_meson_non_prompt_meson.Histo1D({"DCAMesonNonPromptMeson", "DCAMesonNonPromptMeson", 300, -1500, 1500}, "fPicDCAxyToPVTopo"); 
  auto h_DCA_meson_non_prompt_meson_2D = df_meson_non_prompt_meson.Histo2D({"NonPromptMesonOrg", "bla", 2000, -1000, 1000, 2000, -1000, 1000}, "fPiccMotherPDG", "LeadingHadron"); 
  
  TString outName = TString::Format("outpionDCA_%s.root",outAddon )  ; 
  TFile* out = TFile::Open(outName.Data(), "recreate");
  
  out->cd(); 
  //this is my start

  HarryPlotter::CheckAndStore(out, h_cand_counter); 
  HarryPlotter::CheckAndStore(out, h_gen_xi_c_counter); 
  HarryPlotter::CheckAndStore(out, h_gen_xi_cc_counter); 
  
  HarryPlotter::CheckAndStore(out, h_PDGCode_prim); 
  HarryPlotter::CheckAndStore(out, h_DCA_prim); 
  HarryPlotter::CheckAndStore(out, h_PDGCode_strong); 
  HarryPlotter::CheckAndStore(out, h_DCA_strong); 
  HarryPlotter::CheckAndStore(out, h_PDGCode_other); 
  HarryPlotter::CheckAndStore(out, h_DCA_other); 
  HarryPlotter::CheckAndStore(out, h_PDGCode_baryon); 
  HarryPlotter::CheckAndStore(out, h_DCA_baryon); 
  HarryPlotter::CheckAndStore(out, h_DCA_baryon_s); 
  HarryPlotter::CheckAndStore(out, h_DCA_baryon_c); 
  HarryPlotter::CheckAndStore(out, h_DCA_baryon_b); 

  HarryPlotter::CheckAndStore(out, h_PDGCode_lead_baryon); 
  HarryPlotter::CheckAndStore(out, h_DCA_baryon_prompt); 
  HarryPlotter::CheckAndStore(out, h_DCA_baryon_prompt_s); 
  HarryPlotter::CheckAndStore(out, h_DCA_baryon_prompt_b); 
  HarryPlotter::CheckAndStore(out, h_DCA_baryon_prompt_c); 

  HarryPlotter::CheckAndStore(out, h_PDGCode_baryon_prompt); 
  HarryPlotter::CheckAndStore(out, h_DCA_baryon_non_prompt); 
  HarryPlotter::CheckAndStore(out, h_DCA_baryon_non_prompt_s); 
  HarryPlotter::CheckAndStore(out, h_DCA_baryon_non_prompt_c); 
  HarryPlotter::CheckAndStore(out, h_DCA_baryon_non_prompt_b); 

  HarryPlotter::CheckAndStore(out, h_PDGCode_baryon_non_prompt); 
  HarryPlotter::CheckAndStore(out, h_DCA_baryon_non_prompt); 
  HarryPlotter::CheckAndStore(out, h_DCA_baryon_non_prompt_s); 
  HarryPlotter::CheckAndStore(out, h_DCA_baryon_non_prompt_c); 
  HarryPlotter::CheckAndStore(out, h_DCA_baryon_non_prompt_b); 

  HarryPlotter::CheckAndStore(out, h_PDGCode_meson); 
  HarryPlotter::CheckAndStore(out, h_DCA_meson); 
  HarryPlotter::CheckAndStore(out, h_DCA_meson_s); 
  HarryPlotter::CheckAndStore(out, h_DCA_meson_c); 
  HarryPlotter::CheckAndStore(out, h_DCA_meson_b); 

  HarryPlotter::CheckAndStore(out, h_PDGCode_meson_prompt); 
  HarryPlotter::CheckAndStore(out, h_DCA_meson_prompt); 
  HarryPlotter::CheckAndStore(out, h_DCA_meson_prompt_s); 
  HarryPlotter::CheckAndStore(out, h_DCA_meson_prompt_c); 
  HarryPlotter::CheckAndStore(out, h_DCA_meson_prompt_b); 

  HarryPlotter::CheckAndStore(out, h_PDGCode_meson_non_prompt); 
  HarryPlotter::CheckAndStore(out, h_DCA_meson_non_prompt); 
  HarryPlotter::CheckAndStore(out, h_DCA_meson_non_prompt_s); 
  HarryPlotter::CheckAndStore(out, h_DCA_meson_non_prompt_c); 
  HarryPlotter::CheckAndStore(out, h_DCA_meson_non_prompt_b); 

  HarryPlotter::CheckAndStore(out, h_PDGCode_meson_non_prompt_baryon); 
  HarryPlotter::CheckAndStore(out, h_DCA_meson_non_prompt_baryon); 

  HarryPlotter::CheckAndStore(out, h_PDGCode_meson_non_prompt_meson); 
  HarryPlotter::CheckAndStore(out, h_DCA_meson_non_prompt_meson); 
  HarryPlotter::CheckAndStore(out, h_DCA_meson_non_prompt_meson_2D); 
  out->Close(); 
  return 0; 
}


/*
  auto df_in_qa = df
  .Define("fPiCCPt", givemyfloatback, {"fPionPt"})
  .Define("absfPiccMotherPDG", makeMeAbsolute, {"fPiccMotherPDG"})
  ; 
  auto in_counter = df_in_qa.Count(); 
  
  auto h_PDGCode = df_in_qa.Histo1D({"PDGCode", "PDGCode", 20000, -10000, 10000}, "fPiccMotherPDG"); 
  auto h_PDGCodeWide = df_in_qa.Histo1D({"PDGCodeWide", "PDGCodeWide", 20000, -1e8, 1e8}, "fPiccMotherPDG"); 
  
  auto h_dca_xy = df_in_qa.Histo1D({"dca_xy", "dca_xy", 3000, -1500, 1500}, {"fPicDCAxyToPVTopo"}); 
  auto h_dca_z = df_in_qa.Histo1D({"dca_z", "dca_z", 3000, -1500, 1500}, {"fPicDCAzToPVTopo"}); 
  
  auto h_dca_xy_pion_130 = df_in_qa.Filter("TMath::Abs(fPiccMotherPDG) == 130").Histo1D({"dca_xy_130", "dca_xy", 3000, -1500, 1500}, {"fPicDCAxyToPVTopo"}); 
  auto h_dca_z_pion_130 = df_in_qa.Filter("TMath::Abs(fPiccMotherPDG) == 130").Histo1D({"dca_z_130", "dca_z", 3000, -1500, 1500}, {"fPicDCAzToPVTopo"}); 
  
  auto h_dca_xy_pion_310 = df_in_qa.Filter("TMath::Abs(fPiccMotherPDG) == 310").Histo1D({"dca_xy_310", "dca_xy", 3000, -1500, 1500}, {"fPicDCAxyToPVTopo"}); 
  auto h_dca_z_pion_310 = df_in_qa.Filter("TMath::Abs(fPiccMotherPDG) == 310").Histo1D({"dca_z_310", "dca_z", 3000, -1500, 1500}, {"fPicDCAzToPVTopo"}); 

  auto h_dca_xy_pion_321 = df_in_qa.Filter("TMath::Abs(fPiccMotherPDG) == 321").Histo1D({"dca_xy_321", "dca_xy", 3000, -1500, 1500}, {"fPicDCAxyToPVTopo"}); 
  auto h_dca_z_pion_321 = df_in_qa.Filter("TMath::Abs(fPiccMotherPDG) == 321").Histo1D({"dca_z_321", "dca_z", 3000, -1500, 1500}, {"fPicDCAzToPVTopo"}); 
  
  auto h_dca_xy_pion_411 = df_in_qa.Filter("TMath::Abs(fPiccMotherPDG) == 411").Histo1D({"dca_xy_411", "dca_xy", 3000, -1500, 1500}, {"fPicDCAxyToPVTopo"}); 
  auto h_dca_z_pion_411 = df_in_qa.Filter("TMath::Abs(fPiccMotherPDG) == 411").Histo1D({"dca_z_411", "dca_z", 3000, -1500, 1500}, {"fPicDCAzToPVTopo"});   

  auto h_dca_xy_pion_421 = df_in_qa.Filter("TMath::Abs(fPiccMotherPDG) == 421").Histo1D({"dca_xy_421", "dca_xy", 3000, -1500, 1500}, {"fPicDCAxyToPVTopo"}); 
  auto h_dca_z_pion_421 = df_in_qa.Filter("TMath::Abs(fPiccMotherPDG) == 421").Histo1D({"dca_z_421", "dca_z", 3000, -1500, 1500}, {"fPicDCAzToPVTopo"}); 

  auto h_dca_xy_pion_431 = df_in_qa.Filter("TMath::Abs(fPiccMotherPDG) == 431").Histo1D({"dca_xy_431", "dca_xy", 3000, -1500, 1500}, {"fPicDCAxyToPVTopo"}); 
  auto h_dca_z_pion_431 = df_in_qa.Filter("TMath::Abs(fPiccMotherPDG) == 431").Histo1D({"dca_z_431", "dca_z", 3000, -1500, 1500}, {"fPicDCAzToPVTopo"}); 
  
  auto h_dca_xy_pion_511 = df_in_qa.Filter("TMath::Abs(fPiccMotherPDG) == 511").Histo1D({"dca_xy_511", "dca_xy", 3000, -1500, 1500}, {"fPicDCAxyToPVTopo"}); 
  auto h_dca_z_pion_511 = df_in_qa.Filter("TMath::Abs(fPiccMotherPDG) == 511").Histo1D({"dca_z_511", "dca_z", 3000, -1500, 1500}, {"fPicDCAzToPVTopo"}); 
  
  auto h_dca_xy_pion_521 = df_in_qa.Filter("TMath::Abs(fPiccMotherPDG) == 521").Histo1D({"dca_xy_521", "dca_xy", 3000, -1500, 1500}, {"fPicDCAxyToPVTopo"}); 
  auto h_dca_z_pion_521 = df_in_qa.Filter("TMath::Abs(fPiccMotherPDG) == 521").Histo1D({"dca_z_521", "dca_z", 3000, -1500, 1500}, {"fPicDCAzToPVTopo"}); 
  
  auto h_dca_xy_pion_531 = df_in_qa.Filter("TMath::Abs(fPiccMotherPDG) == 531").Histo1D({"dca_xy_531", "dca_xy", 3000, -1500, 1500}, {"fPicDCAxyToPVTopo"}); 
  auto h_dca_z_pion_531 = df_in_qa.Filter("TMath::Abs(fPiccMotherPDG) == 531").Histo1D({"dca_z_531", "dca_z", 3000, -1500, 1500}, {"fPicDCAzToPVTopo"}); 
  
  auto h_dca_xy_pion_541 = df_in_qa.Filter("TMath::Abs(fPiccMotherPDG) == 541").Histo1D({"dca_xy_541", "dca_xy", 3000, -1500, 1500}, {"fPicDCAxyToPVTopo"}); 
  auto h_dca_z_pion_541 = df_in_qa.Filter("TMath::Abs(fPiccMotherPDG) == 541").Histo1D({"dca_z_541", "dca_z", 3000, -1500, 1500}, {"fPicDCAzToPVTopo"}); 
  
  auto h_dca_xy_pion_3122 = df_in_qa.Filter("TMath::Abs(fPiccMotherPDG) == 3122").Histo1D({"dca_xy_3122", "dca_xy", 3000, -1500, 1500}, {"fPicDCAxyToPVTopo"}); 
  auto h_dca_z_pion_3122 = df_in_qa.Filter("TMath::Abs(fPiccMotherPDG) == 3122").Histo1D({"dca_z_3122", "dca_z", 3000, -1500, 1500}, {"fPicDCAzToPVTopo"}); 

  auto h_dca_xy_pion_3112 = df_in_qa.Filter("TMath::Abs(fPiccMotherPDG) == 3112").Histo1D({"dca_xy_3112", "dca_xy", 3000, -1500, 1500}, {"fPicDCAxyToPVTopo"}); 
  auto h_dca_z_pion_3112 = df_in_qa.Filter("TMath::Abs(fPiccMotherPDG) == 3112").Histo1D({"dca_z_3112", "dca_z", 3000, -1500, 1500}, {"fPicDCAzToPVTopo"}); 

  auto h_dca_xy_pion_3222 = df_in_qa.Filter("TMath::Abs(fPiccMotherPDG) == 3222").Histo1D({"dca_xy_3222", "dca_xy", 3000, -1500, 1500}, {"fPicDCAxyToPVTopo"}); 
  auto h_dca_z_pion_3222 = df_in_qa.Filter("TMath::Abs(fPiccMotherPDG) == 3222").Histo1D({"dca_z_3222", "dca_z", 3000, -1500, 1500}, {"fPicDCAzToPVTopo"}); 

  auto h_dca_xy_pion_3322 = df_in_qa.Filter("TMath::Abs(fPiccMotherPDG) == 3322").Histo1D({"dca_xy_3322", "dca_xy", 3000, -1500, 1500}, {"fPicDCAxyToPVTopo"}); 
  auto h_dca_z_pion_3322 = df_in_qa.Filter("TMath::Abs(fPiccMotherPDG) == 3322").Histo1D({"dca_z_3322", "dca_z", 3000, -1500, 1500}, {"fPicDCAzToPVTopo"}); 
  
  auto h_dca_xy_pion_3312 = df_in_qa.Filter("TMath::Abs(fPiccMotherPDG) == 3312").Histo1D({"dca_xy_3312", "dca_xy", 3000, -1500, 1500}, {"fPicDCAxyToPVTopo"}); 
  auto h_dca_z_pion_3312 = df_in_qa.Filter("TMath::Abs(fPiccMotherPDG) == 3312").Histo1D({"dca_z_3312", "dca_z", 3000, -1500, 1500}, {"fPicDCAzToPVTopo"}); 
  
  auto h_dca_xy_pion_3334 = df_in_qa.Filter("TMath::Abs(fPiccMotherPDG) == 3334").Histo1D({"dca_xy_3334", "dca_xy", 3000, -1500, 1500}, {"fPicDCAxyToPVTopo"}); 
  auto h_dca_z_pion_3334 = df_in_qa.Filter("TMath::Abs(fPiccMotherPDG) == 3334").Histo1D({"dca_z_3334", "dca_z", 3000, -1500, 1500}, {"fPicDCAzToPVTopo"}); 
  
  auto h_dca_xy_pion_4122 = df_in_qa.Filter("TMath::Abs(fPiccMotherPDG) == 4122").Histo1D({"dca_xy_4122", "dca_xy", 3000, -1500, 1500}, {"fPicDCAxyToPVTopo"}); 
  auto h_dca_z_pion_4122 = df_in_qa.Filter("TMath::Abs(fPiccMotherPDG) == 4122").Histo1D({"dca_z_4122", "dca_z", 3000, -1500, 1500}, {"fPicDCAzToPVTopo"}); 
  
  auto h_dca_xy_pion_4132 = df_in_qa.Filter("TMath::Abs(fPiccMotherPDG) == 4132").Histo1D({"dca_xy_4132", "dca_xy", 3000, -1500, 1500}, {"fPicDCAxyToPVTopo"}); 
  auto h_dca_z_pion_4132 = df_in_qa.Filter("TMath::Abs(fPiccMotherPDG) == 4132").Histo1D({"dca_z_4132", "dca_z", 3000, -1500, 1500}, {"fPicDCAzToPVTopo"}); 
  
  auto h_dca_xy_pion_4232 = df_in_qa.Filter("TMath::Abs(fPiccMotherPDG) == 4232").Histo1D({"dca_xy_4232", "dca_xy", 3000, -1500, 1500}, {"fPicDCAxyToPVTopo"}); 
  auto h_dca_z_pion_4232 = df_in_qa.Filter("TMath::Abs(fPiccMotherPDG) == 4232").Histo1D({"dca_z_4232", "dca_z", 3000, -1500, 1500}, {"fPicDCAzToPVTopo"}); 
  
  auto h_dca_xy_pion_4332 = df_in_qa.Filter("TMath::Abs(fPiccMotherPDG) == 4332").Histo1D({"dca_xy_4332", "dca_xy", 3000, -1500, 1500}, {"fPicDCAxyToPVTopo"}); 
  auto h_dca_z_pion_4332 = df_in_qa.Filter("TMath::Abs(fPiccMotherPDG) == 4332").Histo1D({"dca_z_4332", "dca_z", 3000, -1500, 1500}, {"fPicDCAzToPVTopo"}); 

  auto h_dca_xy_pion_5122 = df_in_qa.Filter("TMath::Abs(fPiccMotherPDG) == 5122").Histo1D({"dca_xy_5122", "dca_xy", 3000, -1500, 1500}, {"fPicDCAxyToPVTopo"}); 
  auto h_dca_z_pion_5122 = df_in_qa.Filter("TMath::Abs(fPiccMotherPDG) == 5122").Histo1D({"dca_z_5122", "dca_z", 3000, -1500, 1500}, {"fPicDCAzToPVTopo"}); 

  auto h_dca_xy_pion_5132 = df_in_qa.Filter("TMath::Abs(fPiccMotherPDG) == 5132").Histo1D({"dca_xy_5132", "dca_xy", 3000, -1500, 1500}, {"fPicDCAxyToPVTopo"}); 
  auto h_dca_z_pion_5132 = df_in_qa.Filter("TMath::Abs(fPiccMotherPDG) == 5132").Histo1D({"dca_z_5132", "dca_z", 3000, -1500, 1500}, {"fPicDCAzToPVTopo"}); 
  
  auto h_dca_xy_pion_5232 = df_in_qa.Filter("TMath::Abs(fPiccMotherPDG) == 5232").Histo1D({"dca_xy_5232", "dca_xy", 3000, -1500, 1500}, {"fPicDCAxyToPVTopo"}); 
  auto h_dca_z_pion_5232 = df_in_qa.Filter("TMath::Abs(fPiccMotherPDG) == 5232").Histo1D({"dca_z_5232", "dca_z", 100, 0, 10}, {"fPicDCAzToPVTopo"}); 
 

  auto h_pT_pion_130 = df_in_qa.Filter("TMath::Abs(fPiccMotherPDG) == 130").Histo1D({"pT_130", "pT", 100, 0, 10}, {"fPiCCPt"}); 
  
  auto h_pT_pion_310 = df_in_qa.Filter("TMath::Abs(fPiccMotherPDG) == 310").Histo1D({"pT_310", "pT", 100, 0, 10}, {"fPiCCPt"}); 

  auto h_pT_pion_321 = df_in_qa.Filter("TMath::Abs(fPiccMotherPDG) == 321").Histo1D({"pT_321", "pT", 100, 0, 10}, {"fPiCCPt"}); 

  auto h_pT_pion_411 = df_in_qa.Filter("TMath::Abs(fPiccMotherPDG) == 411").Histo1D({"pT_411", "pT", 100, 0, 10}, {"fPiCCPt"}); 

  auto h_pT_pion_421 = df_in_qa.Filter("TMath::Abs(fPiccMotherPDG) == 421").Histo1D({"pT_421", "pT", 100, 0, 10}, {"fPiCCPt"}); 

  auto h_pT_pion_431 = df_in_qa.Filter("TMath::Abs(fPiccMotherPDG) == 431").Histo1D({"pT_431", "pT", 100, 0, 10}, {"fPiCCPt"}); 
  
  auto h_pT_pion_511 = df_in_qa.Filter("TMath::Abs(fPiccMotherPDG) == 511").Histo1D({"pT_511", "pT", 100, 0, 10}, {"fPiCCPt"}); 
  
  auto h_pT_pion_521 = df_in_qa.Filter("TMath::Abs(fPiccMotherPDG) == 521").Histo1D({"pT_521", "pT", 100, 0, 10}, {"fPiCCPt"}); 
  
  auto h_pT_pion_531 = df_in_qa.Filter("TMath::Abs(fPiccMotherPDG) == 531").Histo1D({"pT_531", "pT", 100, 0, 10}, {"fPiCCPt"}); 
  
  auto h_pT_pion_541 = df_in_qa.Filter("TMath::Abs(fPiccMotherPDG) == 541").Histo1D({"pT_541", "pT", 100, 0, 10}, {"fPiCCPt"}); 
  
  auto h_pT_pion_3122 = df_in_qa.Filter("TMath::Abs(fPiccMotherPDG) == 3122").Histo1D({"pT_3122", "pT", 100, 0, 10}, {"fPiCCPt"}); 

  auto h_pT_pion_3112 = df_in_qa.Filter("TMath::Abs(fPiccMotherPDG) == 3112").Histo1D({"pT_3112", "pT", 100, 0, 10}, {"fPiCCPt"}); 

  auto h_pT_pion_3222 = df_in_qa.Filter("TMath::Abs(fPiccMotherPDG) == 3222").Histo1D({"pT_3222", "pT", 100, 0, 10}, {"fPiCCPt"}); 

  auto h_pT_pion_3322 = df_in_qa.Filter("TMath::Abs(fPiccMotherPDG) == 3322").Histo1D({"pT_3322", "pT", 100, 0, 10}, {"fPiCCPt"}); 
  
  auto h_pT_pion_3312 = df_in_qa.Filter("TMath::Abs(fPiccMotherPDG) == 3312").Histo1D({"pT_3312", "pT", 100, 0, 10}, {"fPiCCPt"}); 
  
  auto h_pT_pion_3334 = df_in_qa.Filter("TMath::Abs(fPiccMotherPDG) == 3334").Histo1D({"pT_3334", "pT", 100, 0, 10}, {"fPiCCPt"}); 
  
  auto h_pT_pion_4122 = df_in_qa.Filter("TMath::Abs(fPiccMotherPDG) == 4122").Histo1D({"pT_4122", "pT", 100, 0, 10}, {"fPiCCPt"}); 
  
  auto h_pT_pion_4132 = df_in_qa.Filter("TMath::Abs(fPiccMotherPDG) == 4132").Histo1D({"pT_4132", "pT", 100, 0, 10}, {"fPiCCPt"}); 
  
  auto h_pT_pion_4232 = df_in_qa.Filter("TMath::Abs(fPiccMotherPDG) == 4232").Histo1D({"pT_4232", "pT", 100, 0, 10}, {"fPiCCPt"}); 
  
  auto h_pT_pion_4332 = df_in_qa.Filter("TMath::Abs(fPiccMotherPDG) == 4332").Histo1D({"pT_4332", "pT", 100, 0, 10}, {"fPiCCPt"}); 

  auto h_pT_pion_5122 = df_in_qa.Filter("TMath::Abs(fPiccMotherPDG) == 5122").Histo1D({"pT_5122", "pT", 100, 0, 10}, {"fPiCCPt"}); 

  auto h_pT_pion_5132 = df_in_qa.Filter("TMath::Abs(fPiccMotherPDG) == 5132").Histo1D({"pT_5132", "pT", 100, 0, 10}, {"fPiCCPt"}); 
  
  auto h_pT_pion_5232 = df_in_qa.Filter("TMath::Abs(fPiccMotherPDG) == 5232").Histo1D({"pT_5232", "pT", 100, 0, 10}, {"fPiCCPt"}); 
 

  auto df_meson_ud = df_in_qa.Filter(Mesons_u_d, {"absfPiccMotherPDG"}); 
  auto dca_xy_meson_ud = df_meson_ud.Histo1D({"dca_xy_meson_ud", "dca_xy_meson_ud", 3000, -1500, 1500}, {"fPicDCAxyToPVTopo"}); 
  auto dca_z_meson_ud = df_meson_ud.Histo1D({"dca_z_meson_ud", "dca_z_meson_ud", 3000, -1500, 1500}, {"fPicDCAzToPVTopo"}); 

  auto df_meson_s = df_in_qa.Filter(Mesons_s, {"absfPiccMotherPDG"}); 
  auto dca_xy_meson_s = df_meson_s.Histo1D({"dca_xy_meson_s", "dca_xy_meson_s", 3000, -1500, 1500}, {"fPicDCAxyToPVTopo"}); 
  auto dca_z_meson_s = df_meson_s.Histo1D({"dca_z_meson_s", "dca_z_meson_s", 3000, -1500, 1500}, {"fPicDCAzToPVTopo"}); 

  auto df_meson_c = df_in_qa.Filter(Mesons_c, {"absfPiccMotherPDG"}); 
  auto dca_xy_meson_c = df_meson_c.Histo1D({"dca_xy_meson_c", "dca_xy_meson_c", 3000, -1500, 1500}, {"fPicDCAxyToPVTopo"}); 
  auto dca_z_meson_c = df_meson_c.Histo1D({"dca_z_meson_c", "dca_z_meson_c", 3000, -1500, 1500}, {"fPicDCAzToPVTopo"}); 

  auto df_meson_b = df_in_qa.Filter(Mesons_b, {"absfPiccMotherPDG"}); 
  auto dca_xy_meson_b = df_meson_b.Histo1D({"dca_xy_meson_b", "dca_xy_meson_b", 3000, -1500, 1500}, {"fPicDCAxyToPVTopo"}); 
  auto dca_z_meson_b = df_meson_b.Histo1D({"dca_z_meson_b", "dca_z_meson_b", 3000, -1500, 1500}, {"fPicDCAzToPVTopo"}); 
  
  auto df_baryon_ud = df_in_qa.Filter(Baryons_u_d, {"absfPiccMotherPDG"});
  auto dca_xy_baryon_ud = df_baryon_ud.Histo1D({"dca_xy_baryon_ud", "dca_xy_baryon_ud", 3000, -1500, 1500}, {"fPicDCAxyToPVTopo"}); 
  auto dca_z_baryon_ud = df_baryon_ud.Histo1D({"dca_z_baryon_ud", "dca_z_baryon_ud", 3000, -1500, 1500}, {"fPicDCAzToPVTopo"}); 

  auto df_baryon_s = df_in_qa.Filter(Baryons_s, {"absfPiccMotherPDG"}); 
  auto dca_xy_baryon_s = df_baryon_s.Histo1D({"dca_xy_baryon_s", "dca_xy_baryon_s", 3000, -1500, 1500}, {"fPicDCAxyToPVTopo"}); 
  auto dca_z_baryon_s = df_baryon_s.Histo1D({"dca_z_baryon_s", "dca_z_baryon_s", 3000, -1500, 1500}, {"fPicDCAzToPVTopo"}); 

  auto df_baryon_c = df_in_qa.Filter(Baryons_c, {"absfPiccMotherPDG"}); 
  auto dca_xy_baryon_c = df_baryon_c.Histo1D({"dca_xy_baryon_c", "dca_xy_baryon_c", 3000, -1500, 1500}, {"fPicDCAxyToPVTopo"}); 
  auto dca_z_baryon_c = df_baryon_c.Histo1D({"dca_z_baryon_c", "dca_z_baryon_c", 3000, -1500, 1500}, {"fPicDCAzToPVTopo"}); 

  auto df_baryon_b = df_in_qa.Filter(Baryons_b, {"absfPiccMotherPDG"}); 
  auto dca_xy_baryon_b = df_baryon_b.Histo1D({"dca_xy_baryon_b", "dca_xy_baryon_b", 3000, -1500, 1500}, {"fPicDCAxyToPVTopo"}); 
  auto dca_z_baryon_b = df_baryon_b.Histo1D({"dca_z_baryon_b", "dca_z_baryon_b", 3000, -1500, 1500}, {"fPicDCAzToPVTopo"}); 

  auto df_other = df_in_qa.Filter(nonOfTheAbove, {"absfPiccMotherPDG"}); 
  auto dca_xy_other = df_other.Histo1D({"dca_xy_other", "dca_xy_other", 3000, -1500, 1500}, {"fPicDCAxyToPVTopo"}); 
  auto dca_z_other = df_other.Histo1D({"dca_z_other", "dca_z_other", 3000, -1500, 1500}, {"fPicDCAzToPVTopo"}); 
*/



/*
  HarryPlotter::CheckAndStore(out, h_PDGCode); 
  HarryPlotter::CheckAndStore(out, h_PDGCodeWide); 
  
  HarryPlotter::CheckAndStore(out, h_dca_xy);
  HarryPlotter::CheckAndStore(out, h_dca_z);

  HarryPlotter::CheckAndStore(out, h_dca_xy_pion_130);
  HarryPlotter::CheckAndStore(out, h_dca_z_pion_130);

  HarryPlotter::CheckAndStore(out, h_dca_xy_pion_310);
  HarryPlotter::CheckAndStore(out, h_dca_z_pion_310);

  HarryPlotter::CheckAndStore(out, h_dca_xy_pion_321);
  HarryPlotter::CheckAndStore(out, h_dca_z_pion_321);

  HarryPlotter::CheckAndStore(out, h_dca_xy_pion_411);
  HarryPlotter::CheckAndStore(out, h_dca_z_pion_411);

  HarryPlotter::CheckAndStore(out, h_dca_xy_pion_421);
  HarryPlotter::CheckAndStore(out, h_dca_z_pion_421);

  HarryPlotter::CheckAndStore(out, h_dca_xy_pion_431);
  HarryPlotter::CheckAndStore(out, h_dca_z_pion_431);

  HarryPlotter::CheckAndStore(out, h_dca_xy_pion_411);
  HarryPlotter::CheckAndStore(out, h_dca_z_pion_411);

  HarryPlotter::CheckAndStore(out, h_dca_xy_pion_511);
  HarryPlotter::CheckAndStore(out, h_dca_z_pion_511);

  HarryPlotter::CheckAndStore(out, h_dca_xy_pion_521);
  HarryPlotter::CheckAndStore(out, h_dca_z_pion_521);

  HarryPlotter::CheckAndStore(out, h_dca_xy_pion_531);
  HarryPlotter::CheckAndStore(out, h_dca_z_pion_531);

  HarryPlotter::CheckAndStore(out, h_dca_xy_pion_541);
  HarryPlotter::CheckAndStore(out, h_dca_z_pion_541);

  HarryPlotter::CheckAndStore(out, h_dca_xy_pion_3122);
  HarryPlotter::CheckAndStore(out, h_dca_z_pion_3122);

  HarryPlotter::CheckAndStore(out, h_dca_xy_pion_3112);
  HarryPlotter::CheckAndStore(out, h_dca_z_pion_3112);

  HarryPlotter::CheckAndStore(out, h_dca_xy_pion_3222);
  HarryPlotter::CheckAndStore(out, h_dca_z_pion_3222);

  HarryPlotter::CheckAndStore(out, h_dca_xy_pion_3322);
  HarryPlotter::CheckAndStore(out, h_dca_z_pion_3322);

  HarryPlotter::CheckAndStore(out, h_dca_xy_pion_3312);
  HarryPlotter::CheckAndStore(out, h_dca_z_pion_3312);

  HarryPlotter::CheckAndStore(out, h_dca_xy_pion_3334);
  HarryPlotter::CheckAndStore(out, h_dca_z_pion_3334);

  HarryPlotter::CheckAndStore(out, h_dca_xy_pion_4122);
  HarryPlotter::CheckAndStore(out, h_dca_z_pion_4122);

  HarryPlotter::CheckAndStore(out, h_dca_xy_pion_4132);
  HarryPlotter::CheckAndStore(out, h_dca_z_pion_4132);

  HarryPlotter::CheckAndStore(out, h_dca_xy_pion_4232);
  HarryPlotter::CheckAndStore(out, h_dca_z_pion_4232);

  HarryPlotter::CheckAndStore(out, h_dca_xy_pion_4332);
  HarryPlotter::CheckAndStore(out, h_dca_z_pion_4332);

  HarryPlotter::CheckAndStore(out, h_dca_xy_pion_5122);
  HarryPlotter::CheckAndStore(out, h_dca_z_pion_5122);

  HarryPlotter::CheckAndStore(out, h_dca_xy_pion_5132);
  HarryPlotter::CheckAndStore(out, h_dca_z_pion_5132);

  HarryPlotter::CheckAndStore(out, h_dca_xy_pion_5232);
  HarryPlotter::CheckAndStore(out, h_dca_z_pion_5232);
  
  
  HarryPlotter::CheckAndStore(out, h_pT_pion_130);

  HarryPlotter::CheckAndStore(out, h_pT_pion_310);

  HarryPlotter::CheckAndStore(out, h_pT_pion_321);

  HarryPlotter::CheckAndStore(out, h_pT_pion_411);

  HarryPlotter::CheckAndStore(out, h_pT_pion_421);

  HarryPlotter::CheckAndStore(out, h_pT_pion_431);

  HarryPlotter::CheckAndStore(out, h_pT_pion_411);

  HarryPlotter::CheckAndStore(out, h_pT_pion_511);

  HarryPlotter::CheckAndStore(out, h_pT_pion_521);

  HarryPlotter::CheckAndStore(out, h_pT_pion_531);

  HarryPlotter::CheckAndStore(out, h_pT_pion_541);

  HarryPlotter::CheckAndStore(out, h_pT_pion_3122);

  HarryPlotter::CheckAndStore(out, h_pT_pion_3112);

  HarryPlotter::CheckAndStore(out, h_pT_pion_3222);

  HarryPlotter::CheckAndStore(out, h_pT_pion_3322);

  HarryPlotter::CheckAndStore(out, h_pT_pion_3312);

  HarryPlotter::CheckAndStore(out, h_pT_pion_3334);

  HarryPlotter::CheckAndStore(out, h_pT_pion_4122);

  HarryPlotter::CheckAndStore(out, h_pT_pion_4132);

  HarryPlotter::CheckAndStore(out, h_pT_pion_4232);

  HarryPlotter::CheckAndStore(out, h_pT_pion_4332);

  HarryPlotter::CheckAndStore(out, h_pT_pion_5122);

  HarryPlotter::CheckAndStore(out, h_pT_pion_5132);

  HarryPlotter::CheckAndStore(out, h_pT_pion_5232);


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
*/
