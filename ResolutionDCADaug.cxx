
#include <iostream>

#include "TFile.h"
#include "TTree.h"
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

int main(int argc, char **argv) {
  HarryPlotter::StyleBox(); 

  ROOT::EnableImplicitMT(); // Tell ROOT you want to go parallel          
  TString filePath = TString::Format("/localstore/alice/hohlweger/analysis/StrangnessTracking/210322"); 
  TFile *file_c = new TFile(filePath+"/omegaccc.root", "READ");

  ROOT::RDataFrame df_ca_c("fTreeTripleC", file_c);

  //add classic quality cuts ... 
  //e.g. fOmegaCPtMC>1.0&&fOmegaCPtMC<3.0 
  //Strangeness tracking: fOmegaHitsAdded >= 3
  std::vector<TString> out_names = {"_integrated"};  
  
  ROOT::RDF::RResultPtr<TH1D> dummy;  
  std::vector<ROOT::RDF::RResultPtr<::TH1D>*> ca_c_dca_daug_topo;  
  std::vector<ROOT::RDF::RResultPtr<::TH1D>*> om_c_dca_daug_topo;  
  std::vector<ROOT::RDF::RResultPtr<::TH1D>*> ca_c_dca_daug_stra;  
  std::vector<ROOT::RDF::RResultPtr<::TH1D>*> om_c_dca_daug_stra;  
  
  //define filter 
  auto dfil_om_c_topo = df_ca_c.
    Filter(HarryPlotter::TopoCuts_om, {"fOmegaCPtMC"}).
    Filter("fCorrectPionFromOmegaC&&fCorrectPionFromOmegaCC&&fCorrectPionFromOmegaCCC");
  
  auto dfil_ca_c_topo = df_ca_c.
    Filter(HarryPlotter::TopoCuts_om, {"fOmegaCPtMC"}).
    Filter("!(fCorrectPionFromOmegaC&&fCorrectPionFromOmegaCC&&fCorrectPionFromOmegaCCC)"); 

  auto dfil_om_c_stra = df_ca_c.
    Filter(HarryPlotter::StraCuts_om, {"fOmegaCPtMC", "fOmegaHitsAdded"}).
    Filter("fCorrectPionFromOmegaC&&fCorrectPionFromOmegaCC&&fCorrectPionFromOmegaCCC");     
  auto dfil_ca_c_stra = df_ca_c.
    Filter(HarryPlotter::StraCuts_om, {"fOmegaCPtMC", "fOmegaHitsAdded"}).
    Filter("!(fCorrectPionFromOmegaC&&fCorrectPionFromOmegaCC&&fCorrectPionFromOmegaCCC)");  
  
  ca_c_dca_daug_topo.emplace_back(new ROOT::RDF::RResultPtr<::TH1D>(dfil_ca_c_topo.Histo1D({"ca_c_dca_daug_topo","#Omega_{c}#rightarrow#Omega+#pi", 500,0, 500}, "fOmegaCDecayDCATopo"))); 
  om_c_dca_daug_topo.emplace_back(new ROOT::RDF::RResultPtr<::TH1D>(dfil_om_c_topo.Histo1D({"om_c_dca_daug_topo","#Omega_{c}#rightarrow#Omega+#pi", 500,0, 500}, "fOmegaCDecayDCATopo"))); 
  
  ca_c_dca_daug_stra.emplace_back(new ROOT::RDF::RResultPtr<::TH1D>(dfil_ca_c_stra.Histo1D({"ca_c_dca_daug_stra","#Omega_{c}#rightarrow#Omega+#pi", 500,0, 500}, "fOmegaCDecayDCAStraTrack"))); 
  om_c_dca_daug_stra.emplace_back(new ROOT::RDF::RResultPtr<::TH1D>(dfil_om_c_stra.Histo1D({"om_c_dca_daug_stra","#Omega_{c}#rightarrow#Omega+#pi", 500,0, 500}, "fOmegaCDecayDCAStraTrack"))); 
  
  
  for (int ipt = 0; ipt < ptbins.size()-1; ++ipt) { 
    TString binstring = TString::Format("_pT_%.1f_%.1f", ptbins[ipt], ptbins[ipt+1]); 
    TString cutstring = TString::Format("(%.1f < fOmegaCPtMC) && (fOmegaCPtMC < %.1f)", ptbins[ipt], ptbins[ipt+1]); 
    out_names.push_back(binstring); 
    
    auto dfil_om_c_topo_pT = dfil_om_c_topo.Filter(cutstring.Data()); 
    auto dfil_ca_c_topo_pT = dfil_ca_c_topo.Filter(cutstring.Data()); 
    
    ca_c_dca_daug_topo.emplace_back(new ROOT::RDF::RResultPtr<::TH1D>(dfil_ca_c_topo_pT.Histo1D({"ca_c_dca_daug_topo"+binstring,"#Omega_{c}#rightarrow#Omega+#pi", 500,0, 500}, "fOmegaCDecayDCATopo"))); 
    om_c_dca_daug_topo.emplace_back(new ROOT::RDF::RResultPtr<::TH1D>(dfil_om_c_topo_pT.Histo1D({"om_c_dca_daug_topo"+binstring,"#Omega_{c}#rightarrow#Omega+#pi", 500,0, 500}, "fOmegaCDecayDCATopo"))); 
    
    auto dfil_om_c_stra_pT = dfil_om_c_stra.Filter(cutstring.Data()); 
    auto dfil_ca_c_stra_pT = dfil_ca_c_stra.Filter(cutstring.Data()); 
    
    ca_c_dca_daug_stra.emplace_back(new ROOT::RDF::RResultPtr<::TH1D>(dfil_ca_c_stra_pT.Histo1D({"ca_c_dca_daug_stra"+binstring,"#Omega_{c}#rightarrow#Omega+#pi", 500,0, 500}, "fOmegaCDecayDCAStraTrack"))); 
    om_c_dca_daug_stra.emplace_back(new ROOT::RDF::RResultPtr<::TH1D>(dfil_om_c_stra_pT.Histo1D({"om_c_dca_daug_stra"+binstring,"#Omega_{c}#rightarrow#Omega+#pi", 500,0, 500}, "fOmegaCDecayDCAStraTrack"))); 
    
  }
  //Write histos 
  TFile *out = TFile::Open("outResolutionsDCADaug.root", "recreate"); 
  for (int ipt = 0; ipt < ptbins.size(); ++ipt) { 
      
    ca_c_dca_daug_topo[ipt]->GetPtr()->Scale(1./ca_c_dca_daug_topo[ipt]->GetPtr()->Integral()); 
    om_c_dca_daug_topo[ipt]->GetPtr()->Scale(1./om_c_dca_daug_topo[ipt]->GetPtr()->Integral()); 
   
    HarryPlotter::CheckAndStore(out, ca_c_dca_daug_topo[ipt]->GetPtr()); 
    HarryPlotter::CheckAndStore(out, om_c_dca_daug_topo[ipt]->GetPtr()); 
    
    ca_c_dca_daug_stra[ipt]->GetPtr()->Scale(1./ca_c_dca_daug_stra[ipt]->GetPtr()->Integral()); 
    om_c_dca_daug_stra[ipt]->GetPtr()->Scale(1./om_c_dca_daug_stra[ipt]->GetPtr()->Integral()); 
    
    HarryPlotter::CheckAndStore(out, ca_c_dca_daug_stra[ipt]->GetPtr()); 
    HarryPlotter::CheckAndStore(out, om_c_dca_daug_stra[ipt]->GetPtr()); 
    
    out->cd(); 

  } 
  out->Close(); 
} 

