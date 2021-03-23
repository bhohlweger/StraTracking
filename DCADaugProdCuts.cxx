
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
 
  TString filePath = TString::Format("/localstore/alice/hohlweger/analysis/StrangnessTracking/210323"); 
  TFile *file_c = new TFile(filePath+"/omegaccc.root", "READ");

  //TH1D* om_evt_counter = (TH1D*)file->Get("hEventCounter"); 
  //TH1D* om_c_evt_counter = (TH1D*)file_c->Get("hEventCounter"); 
  
  // double sampleRatio = om_evt_counter->GetEntries()/om_c_evt_counter->GetEntries(); 
  // double abundance = 500*20; // 500 for the production ratio between omega and omega_c, 20 for the branching ratio 
  //double scaleFactor = abundance/sampleRatio; 
    
  ROOT::RDataFrame df_ca_c("fTreeTripleC", file_c);
  
  //add classic quality cuts ... 
  ROOT::RDF::RResultPtr<TH1D> dummy;  
  
  //Define filter 
  //Input filter
 
  auto dfil_om_c_topo = df_ca_c.
    Filter(HarryPlotter::TopoCuts_om, {"fOmegaCPtMC"}).
    Filter("fCorrectPionFromOmegaC&&fCorrectPionFromOmegaCC&&fCorrectPionFromOmegaCCC&&fFirstCombinationCCC").
    Define("dcaDaugProd", "fOmegaCDecayDCATopo*fOmegaCCDecayDCATopo*fOmegaCCCDecayDCATopo");
  
  auto dfil_ca_c_topo = df_ca_c.
    Filter(HarryPlotter::TopoCuts_om, {"fOmegaCPtMC"}).
    Filter("!(fCorrectPionFromOmegaC||fCorrectPionFromOmegaCC||fCorrectPionFromOmegaCCC)").
    Define("dcaDaugProd", "fOmegaCDecayDCATopo*fOmegaCCDecayDCATopo*fOmegaCCCDecayDCATopo"); 

  auto dfil_om_c_stra = df_ca_c.
    Filter(HarryPlotter::StraCuts_om, {"fOmegaCPtMC", "fOmegaHitsAdded"}).
    Filter("fCorrectPionFromOmegaC&&fCorrectPionFromOmegaCC&&fCorrectPionFromOmegaCCC&&fFirstCombinationCCC").
    Define("dcaDaugProd", "fOmegaCDecayDCAStraTrack*fOmegaCCDecayDCAStraTrack*fOmegaCCCDecayDCAStraTrack");     
  auto dfil_ca_c_stra = df_ca_c.
    Filter(HarryPlotter::StraCuts_om, {"fOmegaCPtMC", "fOmegaHitsAdded"}).
    Filter("!(fCorrectPionFromOmegaC||fCorrectPionFromOmegaCC||fCorrectPionFromOmegaCCC)").
    Define("dcaDaugProd", "fOmegaCDecayDCAStraTrack*fOmegaCCDecayDCAStraTrack*fOmegaCCCDecayDCAStraTrack");  

  //Plots 
  //DCA Product 
  auto ca_dca_daug_prod_topo = dfil_ca_c_topo.Histo1D({"ca_dca_daug_prod_topo","#Omega_{c}#rightarrow#Omega+#pi", 10000,0, 0.5e6}, "dcaDaugProd");
  auto ca_dca_daug_prod_stra = dfil_ca_c_stra.Histo1D({"ca_dca_daug_prod_stra","#Omega_{c}#rightarrow#Omega+#pi", 10000,0, 0.5e6}, "dcaDaugProd");
  auto om_c_dca_daug_prod_topo = dfil_om_c_topo.Histo1D({"om_c_dca_daug_prod_topo","#Omega_{c}#rightarrow#Omega+#pi", 10000,0, 0.5e6}, "dcaDaugProd");
  auto om_c_dca_daug_prod_stra = dfil_om_c_stra.Histo1D({"om_c_dca_daug_prod_stra","#Omega_{c}#rightarrow#Omega+#pi", 10000,0, 0.5e6}, "dcaDaugProd");
  
  //MASS
  auto ca_mass_topo = dfil_ca_c_topo.Histo1D({"ca_mass_topo","#Omega_{c}#rightarrow#Omega+#pi", 200, 4.2, 5.2}, "fOmegaCCCMassTopo"); 
  auto ca_mass_stra = dfil_ca_c_stra.Histo1D({"ca_mass_stra","#Omega_{c}#rightarrow#Omega+#pi", 200, 4.2, 5.2}, "fOmegaCCCMassTopo"); 
  auto om_c_mass_topo = dfil_om_c_topo.Histo1D({"om_c_mass_topo","#Omega_{c}#rightarrow#Omega+#pi", 200, 4.2, 5.2}, "fOmegaCCCMassStraTrack"); 
  auto om_c_mass_stra = dfil_om_c_stra.Histo1D({"om_c_mass_stra","#Omega_{c}#rightarrow#Omega+#pi", 200, 4.2, 5.2}, "fOmegaCCCMassStraTrack");   
  /*
  auto ca_mass_topo_cut = dfil_ca_topo_cut.Histo1D({"ca_mass_topo_cut","#Omega_{c}#rightarrow#Omega+#pi", 200, 2.4, 3.0}, "fOmegacMassTopo"); 
  auto ca_mass_stra_cut = dfil_ca_stra_cut.Histo1D({"ca_mass_stra_cut","#Omega_{c}#rightarrow#Omega+#pi", 200, 2.4, 3.0}, "fOmegacMassStraTrack"); 
  auto om_c_mass_topo_cut = dfil_om_c_topo_cut.Histo1D({"om_c_mass_topo_cut","#Omega_{c}#rightarrow#Omega+#pi", 200, 2.4, 3.0}, "fOmegacMassTopo"); 
  auto om_c_mass_stra_cut = dfil_om_c_stra_cut.Histo1D({"om_c_mass_stra_cut","#Omega_{c}#rightarrow#Omega+#pi", 200, 2.4, 3.0}, "fOmegacMassStraTrack");   
  */
  //Write histos 
  TFile *out = TFile::Open("outDCADaugProductCuts.root", "recreate"); 
  
  ca_dca_daug_prod_topo->Scale(1./ca_dca_daug_prod_topo->Integral()); 
  ca_dca_daug_prod_stra->Scale(1./ca_dca_daug_prod_stra->Integral()); 
  om_c_dca_daug_prod_topo->Scale(1./om_c_dca_daug_prod_topo->Integral()); 
  om_c_dca_daug_prod_stra->Scale(1./om_c_dca_daug_prod_stra->Integral()); 
  
  HarryPlotter::CheckAndStore(out, ca_dca_daug_prod_topo); 
  HarryPlotter::CheckAndStore(out, ca_dca_daug_prod_stra); 
  HarryPlotter::CheckAndStore(out, om_c_dca_daug_prod_topo); 
  HarryPlotter::CheckAndStore(out, om_c_dca_daug_prod_stra); 
  
  HarryPlotter::CumulateAndStore(out, ca_dca_daug_prod_topo); 
  HarryPlotter::CumulateAndStore(out, ca_dca_daug_prod_stra); 
  HarryPlotter::CumulateAndStore(out, om_c_dca_daug_prod_topo); 
  HarryPlotter::CumulateAndStore(out, om_c_dca_daug_prod_stra); 
    
  // ca_mass_topo->Scale(scaleFactor); 
  // ca_mass_stra->Scale(scaleFactor); 
  HarryPlotter::CheckAndStore(out, ca_mass_topo); 
  HarryPlotter::CheckAndStore(out, ca_mass_stra); 
  HarryPlotter::CheckAndStore(out, om_c_mass_topo); 
  HarryPlotter::CheckAndStore(out, om_c_mass_stra); 
  /*
  ca_mass_topo_cut->Scale(scaleFactor); 
  ca_mass_stra_cut->Scale(scaleFactor); 
  HarryPlotter::CheckAndStore(out, ca_mass_topo_cut); 
  HarryPlotter::CheckAndStore(out, ca_mass_stra_cut); 
  HarryPlotter::CheckAndStore(out, om_c_mass_topo_cut); 
  HarryPlotter::CheckAndStore(out, om_c_mass_stra_cut); 
  */
  HarryPlotter::AverageBackground(out, om_c_mass_topo, ca_mass_topo); 
  HarryPlotter::AverageBackground(out, om_c_mass_stra, ca_mass_stra); 
  /*
  HarryPlotter::AverageBackground(out, om_c_mass_topo_cut, ca_mass_topo_cut); 
  HarryPlotter::AverageBackground(out, om_c_mass_stra_cut, ca_mass_stra_cut); 
  */
  out->Close(); 
  return 0; 
} 

