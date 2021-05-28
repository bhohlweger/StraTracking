
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
   
  TString filePathSignal = TString::Format("/localstore/alice/hohlweger/analysis/StrangnessTracking/210323/gimmezetree.root"); 
  TString filePathBkg = TString::Format("%s", HarryPlotter::FilePathPythia()); 

  TFile *file_sgn = new TFile(filePathSignal, "READ");
  TFile *file_bkg = new TFile(filePathBkg, "READ");

  ROOT::RDataFrame df_sgn("fTreeTripleC", file_sgn);
  ROOT::RDataFrame df_bkg("fTreeTripleC", file_bkg);

  //add classic quality cuts ... 
  //e.g. fOmegaCPtMC>1.0&&fOmegaCPtMC<3.0 
  //Strangeness tracking: fOmegaHitsAdded >= 3

  TString selection_signal = 
    TString::Format("fCorrectPionFromOmegaC&&fCorrectPionFromOmegaCC&&fCorrectPionFromOmegaCCC"); 
  
  TString selection_background = 
    TString::Format("fFirstCombinationC&&fFirstCombinationCC&&fFirstCombinationCCC&&!fCorrectPionFromOmegaC"); 
  
  double sampleRatio = 1.;//om_evt_counter->GetEntries()/om_c_evt_counter->GetEntries(); 
  double abundance = 1; // 500 for the production ratio between omega and omega_c, 20 for the branching ratio 
  double scaleFactor = abundance/sampleRatio; 
    
  //Define filter 
  //Topological 
  auto dfil_om_topo = df_sgn
    .Filter(HarryPlotter::TopoCuts_om, {"fOmegaCCCPtTopo"})
    .Filter(selection_signal.Data())
    .Define("prod_xy_om_c_topo", "fOmegaDCAxyToPVTopo*fOmegaCPionDCAxyToPV")
    .Define("prod_z_om_c_topo", "fOmegaDCAzToPVTopo*fOmegaCPionDCAzToPV")
    
    .Define("prod_xy_om_cc_topo", "fOmegaCDCAxyToPVTopo*fOmegaCCPionDCAxyToPV")
    .Define("prod_z_om_cc_topo", "fOmegaCDCAzToPVTopo*fOmegaCCPionDCAzToPV")
    
    .Define("prod_xy_om_ccc_topo", "fOmegaCCDCAxyToPVTopo*fOmegaCCCPionDCAxyToPV")
    .Define("prod_z_om_ccc_topo", "fOmegaCCDCAzToPVTopo*fOmegaCCCPionDCAzToPV");

  auto dfil_ca_topo = df_bkg
    .Filter(HarryPlotter::TopoCuts_om, {"fOmegaCCCPtTopo"})
    //.Filter(selection_background.Data())
    .Define("prod_xy_ca_c_topo", "fOmegaDCAxyToPVTopo*fOmegaCPionDCAxyToPV")
    .Define("prod_z_ca_c_topo", "fOmegaDCAzToPVTopo*fOmegaCPionDCAzToPV")
    
    .Define("prod_xy_ca_cc_topo", "fOmegaCDCAxyToPVTopo*fOmegaCCPionDCAxyToPV")
    .Define("prod_z_ca_cc_topo", "fOmegaCDCAzToPVTopo*fOmegaCCPionDCAzToPV")
    
    .Define("prod_xy_ca_ccc_topo", "fOmegaCCDCAxyToPVTopo*fOmegaCCCPionDCAxyToPV")
    .Define("prod_z_ca_ccc_topo", "fOmegaCCDCAzToPVTopo*fOmegaCCCPionDCAzToPV");

  //Strangeness Tracking 
  auto dfil_om_stra = df_sgn
    .Filter(HarryPlotter::StraCuts_om, {"fOmegaCCCPtTopo", "fOmegaHitsAdded"})
    .Filter(selection_signal.Data())
    .Define("prod_xy_om_c_stra", "fOmegaDCAxyToPVStraTrack*fOmegaCPionDCAxyToPV")
    .Define("prod_z_om_c_stra", "fOmegaDCAzToPVStraTrack*fOmegaCPionDCAzToPV")
    
    .Define("prod_xy_om_cc_stra", "fOmegaCDCAxyToPVStraTrack*fOmegaCCPionDCAxyToPV")
    .Define("prod_z_om_cc_stra", "fOmegaCDCAzToPVStraTrack*fOmegaCCPionDCAzToPV")
    
    .Define("prod_xy_om_ccc_stra", "fOmegaCCDCAxyToPVStraTrack*fOmegaCCCPionDCAxyToPV")
    .Define("prod_z_om_ccc_stra", "fOmegaCCDCAzToPVStraTrack*fOmegaCCCPionDCAzToPV");

  auto dfil_ca_stra = df_bkg
    .Filter(HarryPlotter::TopoCuts_om, {"fOmegaCCCPtTopo"})
    //.Filter(HarryPlotter::StraCuts_om, {"fOmegaCCCPtTopo", "fOmegaHitsAdded"})
    //.Filter(selection_background.Data())
    .Define("prod_xy_ca_c_stra", "fOmegaDCAxyToPVStraTrack*fOmegaCPionDCAxyToPV")
    .Define("prod_z_ca_c_stra", "fOmegaDCAzToPVStraTrack*fOmegaCPionDCAzToPV")
    
    .Define("prod_xy_ca_cc_stra", "fOmegaCDCAxyToPVStraTrack*fOmegaCCPionDCAxyToPV")
    .Define("prod_z_ca_cc_stra", "fOmegaCDCAzToPVStraTrack*fOmegaCCPionDCAzToPV")
    
    .Define("prod_xy_ca_ccc_stra", "fOmegaCCDCAxyToPVStraTrack*fOmegaCCCPionDCAxyToPV")
    .Define("prod_z_ca_ccc_stra", "fOmegaCCDCAzToPVStraTrack*fOmegaCCCPionDCAzToPV");

  //Cut filter   
  auto dfil_om_topo_cut = dfil_om_topo
    .Filter(HarryPlotter::TopoCut_c_dca, {"prod_xy_om_c_topo"}) 
    .Filter(HarryPlotter::TopoCut_cc_dca, {"prod_xy_om_cc_topo"}) 
    .Filter(HarryPlotter::TopoCut_ccc_dca, {"prod_xy_om_ccc_topo"}); 

  auto dfil_ca_topo_cut = dfil_ca_topo
    .Filter(HarryPlotter::TopoCut_c_dca, {"prod_xy_ca_c_topo"}) 
    .Filter(HarryPlotter::TopoCut_cc_dca, {"prod_xy_ca_cc_topo"}) 
    .Filter(HarryPlotter::TopoCut_ccc_dca, {"prod_xy_ca_ccc_topo"}); 

  auto dfil_om_stra_cut = dfil_om_stra
    .Filter(HarryPlotter::StraCut_c_dca, {"prod_xy_om_c_stra"}) 
    .Filter(HarryPlotter::StraCut_cc_dca, {"prod_xy_om_cc_stra"}) 
    .Filter(HarryPlotter::StraCut_ccc_dca, {"prod_xy_om_ccc_stra"}); 

  auto dfil_ca_stra_cut = dfil_ca_stra
    .Filter(HarryPlotter::StraCut_c_dca, {"prod_xy_ca_c_stra"}) 
    .Filter(HarryPlotter::StraCut_cc_dca, {"prod_xy_ca_cc_stra"}) 
    .Filter(HarryPlotter::StraCut_ccc_dca, {"prod_xy_ca_ccc_stra"}); 

  //Plots 
  //DCA Product 
  //omega_c
  auto om_c_dca_om_dca_pi_topo = dfil_om_topo.Histo1D({"om_c_dca_om_dca_pi_topo","#Omega_{c}#rightarrow#Omega+#pi", 750,-5000, 5000}, "prod_xy_om_c_topo");
  auto ca_c_dca_om_dca_pi_topo = dfil_ca_topo.Histo1D({"ca_c_dca_om_dca_pi_topo","#Omega_{c}#rightarrow#Omega+#pi", 750,-5000, 5000}, "prod_xy_ca_c_topo");
  
  auto om_c_dca_om_dca_pi_stra = dfil_om_stra.Histo1D({"om_c_dca_om_dca_pi_stra","#Omega_{c}#rightarrow#Omega+#pi", 750,-5000, 5000}, "prod_xy_om_c_stra");
  auto ca_c_dca_om_dca_pi_stra = dfil_ca_stra.Histo1D({"ca_c_dca_om_dca_pi_stra","#Omega_{c}#rightarrow#Omega+#pi", 750,-5000, 5000}, "prod_xy_ca_c_stra");

  //omega_cc
  auto om_cc_dca_om_dca_pi_topo = dfil_om_topo.Histo1D({"om_cc_dca_om_dca_pi_topo","#Omega_{c}#rightarrow#Omega+#pi", 750,-5000, 5000}, "prod_xy_om_cc_topo");
  auto ca_cc_dca_om_dca_pi_topo = dfil_ca_topo.Histo1D({"ca_cc_dca_om_dca_pi_topo","#Omega_{c}#rightarrow#Omega+#pi", 750,-5000, 5000}, "prod_xy_ca_cc_topo");
  
  auto om_cc_dca_om_dca_pi_stra = dfil_om_stra.Histo1D({"om_cc_dca_om_dca_pi_stra","#Omega_{c}#rightarrow#Omega+#pi", 750,-5000, 5000}, "prod_xy_om_cc_stra");
  auto ca_cc_dca_om_dca_pi_stra = dfil_ca_stra.Histo1D({"ca_cc_dca_om_dca_pi_stra","#Omega_{c}#rightarrow#Omega+#pi", 750,-5000, 5000}, "prod_xy_ca_cc_stra");
  
  //omega_ccc
  auto om_ccc_dca_om_dca_pi_topo = dfil_om_topo.Histo1D({"om_ccc_dca_om_dca_pi_topo","#Omega_{c}#rightarrow#Omega+#pi", 750,-5000, 5000}, "prod_xy_om_ccc_topo");
  auto ca_ccc_dca_om_dca_pi_topo = dfil_ca_topo.Histo1D({"ca_ccc_dca_om_dca_pi_topo","#Omega_{c}#rightarrow#Omega+#pi", 750,-5000, 5000}, "prod_xy_ca_ccc_topo");
  
  auto om_ccc_dca_om_dca_pi_stra = dfil_om_stra.Histo1D({"om_ccc_dca_om_dca_pi_stra","#Omega_{c}#rightarrow#Omega+#pi", 750,-5000, 5000}, "prod_xy_om_ccc_stra");
  auto ca_ccc_dca_om_dca_pi_stra = dfil_ca_stra.Histo1D({"ca_ccc_dca_om_dca_pi_stra","#Omega_{c}#rightarrow#Omega+#pi", 750,-5000, 5000}, "prod_xy_ca_ccc_stra");
  
  //MASS
  auto om_c_mass_topo = dfil_om_topo.Histo1D({"om_c_mass_topo","#Omega_{c}#rightarrow#Omega+#pi", 400, 3.2, 6.2}, "fOmegaCCCMassTopo"); 
  auto ca_c_mass_topo = dfil_ca_topo.Histo1D({"ca_mass_topo","#Omega_{c}#rightarrow#Omega+#pi", 400, 3.2, 6.2}, "fOmegaCCCMassTopo"); 
  
  auto ca_c_mass_stra = dfil_ca_stra.Histo1D({"ca_mass_stra","#Omega_{c}#rightarrow#Omega+#pi", 400, 3.2, 6.2}, "fOmegaCCCMassStraTrack"); 
  auto om_c_mass_stra = dfil_om_stra.Histo1D({"om_c_mass_stra","#Omega_{c}#rightarrow#Omega+#pi", 400, 3.2, 6.2}, "fOmegaCCCMassStraTrack");   
  
  auto om_c_mass_topo_cut = dfil_om_topo_cut.Histo1D({"om_c_mass_topo_cut","#Omega_{c}#rightarrow#Omega+#pi", 400, 3.2, 6.2}, "fOmegaCCCMassTopo"); 
  auto ca_c_mass_topo_cut = dfil_ca_topo_cut.Histo1D({"ca_mass_topo_cut","#Omega_{c}#rightarrow#Omega+#pi", 400, 3.2, 6.2}, "fOmegaCCCMassTopo"); 
  
  auto ca_c_mass_stra_cut = dfil_ca_stra_cut.Histo1D({"ca_mass_stra_cut","#Omega_{c}#rightarrow#Omega+#pi", 400, 3.2, 6.2}, "fOmegaCCCMassStraTrack"); 
  auto om_c_mass_stra_cut = dfil_om_stra_cut.Histo1D({"om_c_mass_stra_cut","#Omega_{c}#rightarrow#Omega+#pi", 400, 3.2, 6.2}, "fOmegaCCCMassStraTrack");   
  
  //Write histos 
  TFile *out = TFile::Open("outDCAProductCuts_ccc.root", "recreate"); 
  

  ca_c_dca_om_dca_pi_topo->Scale(1./ca_c_dca_om_dca_pi_topo->Integral()); 
  om_c_dca_om_dca_pi_topo->Scale(1./om_c_dca_om_dca_pi_topo->Integral()); 
  ca_c_dca_om_dca_pi_stra->Scale(1./ca_c_dca_om_dca_pi_stra->Integral()); 
  om_c_dca_om_dca_pi_stra->Scale(1./om_c_dca_om_dca_pi_stra->Integral()); 
  
  ca_cc_dca_om_dca_pi_topo->Scale(1./ca_cc_dca_om_dca_pi_topo->Integral()); 
  om_cc_dca_om_dca_pi_topo->Scale(1./om_cc_dca_om_dca_pi_topo->Integral()); 
  ca_cc_dca_om_dca_pi_stra->Scale(1./ca_cc_dca_om_dca_pi_stra->Integral()); 
  om_cc_dca_om_dca_pi_stra->Scale(1./om_cc_dca_om_dca_pi_stra->Integral()); 
  
  ca_ccc_dca_om_dca_pi_topo->Scale(1./ca_ccc_dca_om_dca_pi_topo->Integral()); 
  om_ccc_dca_om_dca_pi_topo->Scale(1./om_ccc_dca_om_dca_pi_topo->Integral()); 
  ca_ccc_dca_om_dca_pi_stra->Scale(1./ca_ccc_dca_om_dca_pi_stra->Integral()); 
  om_ccc_dca_om_dca_pi_stra->Scale(1./om_ccc_dca_om_dca_pi_stra->Integral()); 
  
  HarryPlotter::CheckAndStore(out,    ca_c_dca_om_dca_pi_topo); 
  HarryPlotter::CheckAndStore(out,    om_c_dca_om_dca_pi_topo); 
  
  HarryPlotter::CumulateAndStore(out, ca_c_dca_om_dca_pi_topo); 
  HarryPlotter::CumulateAndStore(out, om_c_dca_om_dca_pi_topo); 
  
  HarryPlotter::CheckAndStore(out,    ca_cc_dca_om_dca_pi_topo); 
  HarryPlotter::CheckAndStore(out,    om_cc_dca_om_dca_pi_topo); 
  
  HarryPlotter::CumulateAndStore(out, ca_cc_dca_om_dca_pi_topo); 
  HarryPlotter::CumulateAndStore(out, om_cc_dca_om_dca_pi_topo); 
  
  HarryPlotter::CheckAndStore(out,    ca_ccc_dca_om_dca_pi_topo); 
  HarryPlotter::CheckAndStore(out,    om_ccc_dca_om_dca_pi_topo); 
  
  HarryPlotter::CumulateAndStore(out, ca_ccc_dca_om_dca_pi_topo);   
  HarryPlotter::CumulateAndStore(out, om_ccc_dca_om_dca_pi_topo); 
  
  
  HarryPlotter::CheckAndStore(out,    ca_c_dca_om_dca_pi_stra); 
  HarryPlotter::CheckAndStore(out,    om_c_dca_om_dca_pi_stra); 
  
  HarryPlotter::CumulateAndStore(out, ca_c_dca_om_dca_pi_stra); 
  HarryPlotter::CumulateAndStore(out, om_c_dca_om_dca_pi_stra); 
  
  HarryPlotter::CheckAndStore(out,    ca_cc_dca_om_dca_pi_stra); 
  HarryPlotter::CheckAndStore(out,    om_cc_dca_om_dca_pi_stra); 
  
  HarryPlotter::CumulateAndStore(out, ca_cc_dca_om_dca_pi_stra); 
  HarryPlotter::CumulateAndStore(out, om_cc_dca_om_dca_pi_stra); 

  HarryPlotter::CheckAndStore(out,    ca_ccc_dca_om_dca_pi_stra); 
  HarryPlotter::CheckAndStore(out,    om_ccc_dca_om_dca_pi_stra); 

  HarryPlotter::CumulateAndStore(out, ca_ccc_dca_om_dca_pi_stra); 
  HarryPlotter::CumulateAndStore(out, om_ccc_dca_om_dca_pi_stra); 

  ca_c_mass_topo->Scale(scaleFactor); 
  ca_c_mass_stra->Scale(scaleFactor); 
  HarryPlotter::CheckAndStore(out, ca_c_mass_topo); 
  HarryPlotter::CheckAndStore(out, om_c_mass_topo); 
  
  HarryPlotter::CheckAndStore(out, ca_c_mass_stra); 
  HarryPlotter::CheckAndStore(out, om_c_mass_stra); 
  
  ca_c_mass_topo_cut->Scale(scaleFactor); 
  ca_c_mass_stra_cut->Scale(scaleFactor); 
  HarryPlotter::CheckAndStore(out, ca_c_mass_topo_cut); 
  HarryPlotter::CheckAndStore(out, om_c_mass_topo_cut); 
  
  HarryPlotter::CheckAndStore(out, ca_c_mass_stra_cut); 
  HarryPlotter::CheckAndStore(out, om_c_mass_stra_cut); 
  
  HarryPlotter::AverageBackground(out, om_c_mass_topo, ca_c_mass_topo); 
  HarryPlotter::AverageBackground(out, om_c_mass_stra, ca_c_mass_stra); 
  
  HarryPlotter::AverageBackground(out, om_c_mass_topo_cut, ca_c_mass_topo_cut); 
  HarryPlotter::AverageBackground(out, om_c_mass_stra_cut, ca_c_mass_stra_cut); 
  
  out->Close(); 
  return 0; 
} 

