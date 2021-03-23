
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
  
  //TString filePath = TString::Format("/localstore/alice/hohlweger/analysis/StrangnessTracking/210323/omegaccc.root"); 
  TString filePath = TString::Format("/data/alice/bhohlweger/test/gimmezetree.root"); 
  //  TString filePath = TString::Format("/localstore/alice/hohlweger/analysis/StrangnessTracking/210323/gimmezetree.root"); 
  TFile *file_c = new TFile(filePath, "READ");

  ROOT::RDataFrame df_ca_c("fTreeTripleC", file_c);

  //add classic quality cuts ... 
  //e.g. fOmegaCPtMC>1.0&&fOmegaCPtMC<3.0 
  //Strangeness tracking: fOmegaHitsAdded >= 3

  TString selection_signal = 
    TString::Format("fCorrectPionFromOmegaC&&fCorrectPionFromOmegaCC&&fCorrectPionFromOmegaCCC"); 
  
  TString selection_background = 
    TString::Format("fFirstCombinationC&&fFirstCombinationCC&&fFirstCombinationCCC&&!fCorrectPionFromOmegaCCC"); 

  //TH1D* om_evt_counter = (TH1D*)file->Get("hEventCounter"); 
  //TH1D* om_c_evt_counter = (TH1D*)file_c->Get("hEventCounter"); 
  
  double sampleRatio = 1.;//om_evt_counter->GetEntries()/om_c_evt_counter->GetEntries(); 
  double abundance = 1; // 500 for the production ratio between omega and omega_c, 20 for the branching ratio 
  double scaleFactor = abundance/sampleRatio; 
    
  //Define filter 
  //Input filter
  auto dfil_om_c_topo = df_ca_c
    .Filter(HarryPlotter::TopoCuts_om, {"fOmegaCPtMC"})
    .Filter(selection_signal.Data())
    .Define("prod_xy_om_c_topo", "fOmegaDCAxyToPVTopo*fOmegaCPionDCAxyToPV")
    .Define("prod_z_om_c_topo", "fOmegaDCAzToPVTopo*fOmegaCPionDCAzToPV");
  
  auto dfil_ca_c_topo = df_ca_c
    .Filter(HarryPlotter::TopoCuts_om, {"fOmegaCPtMC"})
    .Filter(selection_background.Data())
    .Define("prod_xy_ca_c_topo", "fOmegaDCAxyToPVTopo*fOmegaCPionDCAxyToPV")
    .Define("prod_z_ca_c_topo", "fOmegaDCAzToPVTopo*fOmegaCPionDCAzToPV");
    
  auto dfil_om_c_stra = df_ca_c
    .Filter(HarryPlotter::StraCuts_om, {"fOmegaCPtMC", "fOmegaHitsAdded"})
    .Filter(selection_signal.Data())
    .Define("prod_xy_om_c_stra", "fOmegaDCAxyToPVStraTrack*fOmegaCPionDCAxyToPV")
    .Define("prod_z_om_c_stra", "fOmegaDCAzToPVStraTrack*fOmegaCPionDCAzToPV");

  auto dfil_ca_c_stra = df_ca_c
    .Filter(HarryPlotter::StraCuts_om, {"fOmegaCPtMC", "fOmegaHitsAdded"})
    .Filter(selection_background.Data())
    .Define("prod_xy_ca_c_stra", "fOmegaDCAxyToPVStraTrack*fOmegaCPionDCAxyToPV")
    .Define("prod_z_ca_c_stra", "fOmegaDCAzToPVStraTrack*fOmegaCPionDCAzToPV");

  //Cut filter   
  auto dfil_ca_c_topo_cut = dfil_ca_c_topo.Filter(HarryPlotter::TopoCut_ccc_dca, {"prod_xy_ca_c_topo"});
  auto dfil_ca_c_stra_cut = dfil_ca_c_stra.Filter(HarryPlotter::StraCut_ccc_dca, {"prod_xy_ca_c_stra"});

  auto dfil_om_c_topo_cut = dfil_om_c_topo.Filter(HarryPlotter::TopoCut_ccc_dca, {"prod_xy_om_c_topo"}); 
  auto dfil_om_c_stra_cut = dfil_om_c_stra.Filter(HarryPlotter::StraCut_ccc_dca, {"prod_xy_om_c_stra"});
  
  //Plots 
  //DCA Product 
  auto ca_c_dca_om_dca_pi_topo = dfil_ca_c_topo.Histo1D({"ca_dca_om_dca_pi_topo","#Omega_{c}#rightarrow#Omega+#pi", 5000,-5000, 5000}, "prod_xy_ca_c_topo");
  auto ca_c_dca_om_dca_pi_stra = dfil_ca_c_stra.Histo1D({"ca_dca_om_dca_pi_stra","#Omega_{c}#rightarrow#Omega+#pi", 5000,-5000, 5000}, "prod_xy_ca_c_stra");
  auto om_c_dca_om_dca_pi_topo = dfil_om_c_topo.Histo1D({"om_c_dca_om_dca_pi_topo","#Omega_{c}#rightarrow#Omega+#pi", 5000,-5000, 5000}, "prod_xy_om_c_topo");
  auto om_c_dca_om_dca_pi_stra = dfil_om_c_stra.Histo1D({"om_c_dca_om_dca_pi_stra","#Omega_{c}#rightarrow#Omega+#pi", 5000,-5000, 5000}, "prod_xy_om_c_stra");
  
  //MASS
  auto ca_c_mass_topo = dfil_ca_c_topo.Histo1D({"ca_mass_topo","#Omega_{c}#rightarrow#Omega+#pi", 200, 4.2, 5.2}, "fOmegaCCCMassTopo"); 
  auto ca_c_mass_stra = dfil_ca_c_stra.Histo1D({"ca_mass_stra","#Omega_{c}#rightarrow#Omega+#pi", 200, 4.2, 5.2}, "fOmegaCCCMassStraTrack"); 
  auto om_c_mass_topo = dfil_om_c_topo.Histo1D({"om_c_mass_topo","#Omega_{c}#rightarrow#Omega+#pi", 200, 4.2, 5.2}, "fOmegaCCCMassTopo"); 
  auto om_c_mass_stra = dfil_om_c_stra.Histo1D({"om_c_mass_stra","#Omega_{c}#rightarrow#Omega+#pi", 200, 4.2, 5.2}, "fOmegaCCCMassStraTrack");   
  
  auto ca_c_mass_topo_cut = dfil_ca_c_topo_cut.Histo1D({"ca_mass_topo_cut","#Omega_{c}#rightarrow#Omega+#pi", 200, 4.2, 5.2}, "fOmegaCCCMassTopo"); 
  auto ca_c_mass_stra_cut = dfil_ca_c_stra_cut.Histo1D({"ca_mass_stra_cut","#Omega_{c}#rightarrow#Omega+#pi", 200, 4.2, 5.2}, "fOmegaCCCMassStraTrack"); 
  auto om_c_mass_topo_cut = dfil_om_c_topo_cut.Histo1D({"om_c_mass_topo_cut","#Omega_{c}#rightarrow#Omega+#pi", 200, 4.2, 5.2}, "fOmegaCCCMassTopo"); 
  auto om_c_mass_stra_cut = dfil_om_c_stra_cut.Histo1D({"om_c_mass_stra_cut","#Omega_{c}#rightarrow#Omega+#pi", 200, 4.2, 5.2}, "fOmegaCCCMassStraTrack");   
  
  //Write histos 
  TFile *out = TFile::Open("outDCAProductCuts_ccc.root", "recreate"); 
  

  ca_c_dca_om_dca_pi_topo->Scale(1./ca_c_dca_om_dca_pi_topo->Integral()); 
  ca_c_dca_om_dca_pi_stra->Scale(1./ca_c_dca_om_dca_pi_stra->Integral()); 
  om_c_dca_om_dca_pi_topo->Scale(1./om_c_dca_om_dca_pi_topo->Integral()); 
  om_c_dca_om_dca_pi_stra->Scale(1./om_c_dca_om_dca_pi_stra->Integral()); 
  
  HarryPlotter::CheckAndStore(out, ca_c_dca_om_dca_pi_topo); 
  HarryPlotter::CheckAndStore(out, ca_c_dca_om_dca_pi_stra); 
  HarryPlotter::CheckAndStore(out, om_c_dca_om_dca_pi_topo); 
  HarryPlotter::CheckAndStore(out, om_c_dca_om_dca_pi_stra); 
  
  HarryPlotter::CumulateAndStore(out, ca_c_dca_om_dca_pi_topo); 
  HarryPlotter::CumulateAndStore(out, ca_c_dca_om_dca_pi_stra); 
  HarryPlotter::CumulateAndStore(out, om_c_dca_om_dca_pi_topo); 
  HarryPlotter::CumulateAndStore(out, om_c_dca_om_dca_pi_stra); 
  
  ca_c_mass_topo->Scale(scaleFactor); 
  ca_c_mass_stra->Scale(scaleFactor); 
  HarryPlotter::CheckAndStore(out, ca_c_mass_topo); 
  HarryPlotter::CheckAndStore(out, ca_c_mass_stra); 
  HarryPlotter::CheckAndStore(out, om_c_mass_topo); 
  HarryPlotter::CheckAndStore(out, om_c_mass_stra); 
  
  ca_c_mass_topo_cut->Scale(scaleFactor); 
  ca_c_mass_stra_cut->Scale(scaleFactor); 
  HarryPlotter::CheckAndStore(out, ca_c_mass_topo_cut); 
  HarryPlotter::CheckAndStore(out, ca_c_mass_stra_cut); 
  HarryPlotter::CheckAndStore(out, om_c_mass_topo_cut); 
  HarryPlotter::CheckAndStore(out, om_c_mass_stra_cut); 
  
  HarryPlotter::AverageBackground(out, om_c_mass_topo, ca_c_mass_topo); 
  HarryPlotter::AverageBackground(out, om_c_mass_stra, ca_c_mass_stra); 
  
  HarryPlotter::AverageBackground(out, om_c_mass_topo_cut, ca_c_mass_topo_cut); 
  HarryPlotter::AverageBackground(out, om_c_mass_stra_cut, ca_c_mass_stra_cut); 
  
  out->Close(); 
  return 0; 
} 

