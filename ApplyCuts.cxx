
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

template <typename T> 
void CheckAndStore(TFile *out, T key) { 
  TObject* obj = out->FindKey(key->GetName());
  if (!obj) { 
    out->cd(); 
    key->Write(); 
  }
}

template <typename T, typename F> 
void NormalizeAndStore(TFile *out, T one, F &normLmbd) { 
  normLmbd(one); 
  CheckAndStore(out, one); 
}

template <typename T, typename F> 
void PlotAndStore(TString name, TFile* out, T one, T two, T three, T four, F &normLmbd) { 
  auto c = new TCanvas(name, ""); 
  c->cd(); 
  
  one->SetLineColor(kPink+7); 
  one->SetMarkerColor(kPink+7); 
  normLmbd(one); 
  one->SetTitle(TString::Format("%s;%s;", one->GetTitle(), name.Data()).Data()); 
  one->Draw(); 
  CheckAndStore(out, one); 
  
  if (two) {
    two->SetLineColor(kAzure-4);
    two->SetMarkerColor(kAzure-4);   
    two->SetLineStyle(2);
    normLmbd(two); 
    two->Draw("same"); 
    CheckAndStore(out, two); 
  }
  if (three) { 
    three->SetLineColor(kOrange+7); 
    three->SetMarkerColor(kOrange+7); 
    normLmbd(three); 
    three->Draw("same"); 
    CheckAndStore(out, three); 
  }
  
  if (four) { 
    four->SetLineColor(kBlue); 
    four->SetMarkerColor(kBlue); 
    four->SetLineStyle(2);
    normLmbd(four);
    four->Draw("same"); 
    CheckAndStore(out, four); 
  }
  c->BuildLegend(); 
  out->cd(); 
  c->Write(); 
}

template <typename T,typename F>  
void CumulatePlotAndStore(TString name, TFile* out, T one, T two, T three, T four, F &normLmbd) { 
  auto c = new TCanvas(name, ""); 
  c->cd(); 
  
  normLmbd(one); 
  TH1D* cu_one = (TH1D*)one->GetCumulative(); 
  cu_one->SetLineColor(kPink+7); 
  cu_one->SetMarkerColor(kPink+7); 
  cu_one->SetTitle(TString::Format("%s;%s;", cu_one->GetTitle(), name.Data()).Data()); 
  cu_one->Draw(); 
  CheckAndStore(out, cu_one); 
  
  if (two) {
    normLmbd(two); 
    TH1D* cu_two = (TH1D*)two->GetCumulative(); 
    cu_two->SetName(TString::Format("%s_cumulative", two->GetName()).Data()); 
    cu_two->SetLineColor(kAzure-4);
    cu_two->SetMarkerColor(kAzure-4);   
    cu_two->SetLineStyle(2);
    cu_two->SetTitle(two->GetTitle()); 
    cu_two->Draw("same"); 
    CheckAndStore(out, cu_two); 
  }
  if (three) { 
    normLmbd(three); 
    TH1D* cu_three = (TH1D*)three->GetCumulative(); 
    cu_three->SetLineColor(kOrange+7); 
    cu_three->SetMarkerColor(kOrange+7); 
    cu_three->SetTitle(three->GetTitle()); 
    cu_three->Draw("same"); 
    CheckAndStore(out, cu_three); 
  }
  
  if (four) { 
    normLmbd(four);
    TH1D* cu_four = (TH1D*)four->GetCumulative(); 
    cu_four->SetLineColor(kBlue); 
    cu_four->SetMarkerColor(kBlue); 
    cu_four->SetLineStyle(2);
    cu_four->SetTitle(four->GetTitle());
    cu_four->Draw("same"); 
    CheckAndStore(out, cu_four); 
  }
  c->BuildLegend();
  out->cd(); 
  c->Write(); 
}


int main(int argc, char **argv) {
  StyleBox(); 

  ROOT::EnableImplicitMT(); // Tell ROOT you want to go parallel          
  TFile *file = new TFile("/localstore/alice/hohlweger/analysis/StrangnessTracking/210315/omeganew.root", "READ");
  TFile *file_c = new TFile("/localstore/alice/hohlweger/analysis/StrangnessTracking/210315/omegacnew.root", "READ");
  
  ROOT::RDataFrame df_om("fTreeOmega", file);
  ROOT::RDataFrame df_ca("fTreeCandidates", file);
  
  ROOT::RDataFrame df_om_c("fTreeOmega", file_c);
  ROOT::RDataFrame df_ca_c("fTreeCandidates", file_c);
  
  //add classic quality cuts ... 
  //e.g. fOmegaPtMC>1.0&&fOmegaPtMC<3.0 

  //Strangeness tracking: fOmegaHitsAdded >= 3
  auto TopoCuts = [] (float pT, float decayRad) {return ((1.0 < pT)&&(pT < 3.0)&&(decayRad>1e-12));}; 
  auto StraCuts = [] (float pT, int iAdded, float decayRad) {return ((1.0 < pT)&&(pT < 3.0)&&(iAdded >= 3)&&(decayRad>1e-12));}; 
  auto Topo_dcaxy_z_Cuts = [] (float dcaxy, float dcaz) { double prod = dcaxy*dcaz; return ((prod) <  -900);}; 
  auto Stra_dcaxy_z_Cuts = [] (float dcaxy, float dcaz) { double prod = dcaxy*dcaz; return ((prod) <  -600);}; 
  
  ROOT::RDF::RResultPtr<TH1D> dummy; //hopefully points to nowhere xD
  
  //define filter 

  //omega 
  auto dfil_om_topo = df_om.Filter(TopoCuts, {"fOmegaPtMC", "fOmegacDecayRadiusTopo"}).Define("om_dca_omXdca_pi_topo","fOmegaDCAxyToPVTopo*fOmegacPionDCAxyToPV"); 
  auto dfil_ca_topo = df_ca.Filter(TopoCuts, {"fOmegaPtMC", "fOmegacDecayRadiusTopo"}).Filter("!fTrueOmegac").Define("om_dca_omXdca_pi_topo","fOmegaDCAxyToPVTopo*fOmegacPionDCAxyToPV"); 
  
  auto dfil_om_stra = df_om.Filter(StraCuts, {"fOmegaPtMC", "fOmegaHitsAdded", "fOmegacDecayRadiusStraTrack"}).Define("om_dca_omXdca_pi_stra","fOmegaDCAxyToPVStraTrack*fOmegacPionDCAxyToPV");
  auto dfil_ca_stra = df_ca.Filter(StraCuts, {"fOmegaPtMC", "fOmegaHitsAdded", "fOmegacDecayRadiusStraTrack"}).Filter("!fTrueOmegac").Define("om_dca_omXdca_pi_stra","fOmegaDCAxyToPVStraTrack*fOmegacPionDCAxyToPV");
 
  //omega_c
  auto dfil_om_c_topo = df_om_c.Filter(TopoCuts, {"fOmegaPtMC", "fOmegacDecayRadiusTopo"}).Define("om_c_dca_omXdca_pi_topo","fOmegaDCAxyToPVTopo*fOmegacPionDCAxyToPV"); 
  auto dfil_ca_c_topo = df_ca_c.Filter(TopoCuts, {"fOmegaPtMC", "fOmegacDecayRadiusTopo"}).Filter("!fTrueOmegac").Define("om_c_dca_omXdca_pi_topo","fOmegaDCAxyToPVTopo*fOmegacPionDCAxyToPV"); 
  
  auto dfil_om_c_topo_dca = dfil_om_c_topo.Filter(Topo_dcaxy_z_Cuts, {"fOmegaDCAxyToPVTopo", "fOmegacPionDCAxyToPV"});
  auto dfil_ca_c_topo_dca = dfil_ca_c_topo.Filter(Topo_dcaxy_z_Cuts, {"fOmegaDCAxyToPVTopo", "fOmegacPionDCAxyToPV"});
  
  auto dfil_om_c_stra = df_om_c.Filter(StraCuts, {"fOmegaPtMC", "fOmegaHitsAdded", "fOmegacDecayRadiusStraTrack"}).Define("om_c_dca_omXdca_pi_stra","fOmegaDCAxyToPVStraTrack*fOmegacPionDCAxyToPV"); 
  auto dfil_ca_c_stra = df_ca_c.Filter(StraCuts, {"fOmegaPtMC", "fOmegaHitsAdded", "fOmegacDecayRadiusStraTrack"}).Filter("!fTrueOmegac").Define("om_c_dca_omXdca_pi_stra","fOmegaDCAxyToPVStraTrack*fOmegacPionDCAxyToPV"); 

  auto dfil_om_c_stra_dca = dfil_om_c_stra.Filter(Stra_dcaxy_z_Cuts, {"fOmegaDCAxyToPVStraTrack", "fOmegacPionDCAxyToPV"});
  auto dfil_ca_c_stra_dca = dfil_ca_c_stra.Filter(Stra_dcaxy_z_Cuts, {"fOmegaDCAxyToPVStraTrack", "fOmegacPionDCAxyToPV"});

  //define histos 
  //Omega 
  //topo 
  auto om_dca_xy_topo = dfil_om_topo.Histo1D({"om_dca_xy_topo","#Omega", 200, -150, 150}, "fOmegaDCAxyToPVTopo");
  auto om_dca_z_topo = dfil_om_topo.Histo1D({"om_dca_z_topo","#Omega", 200, -150, 150}, "fOmegaDCAzToPVTopo");
  auto ca_dca_om_dca_pi_topo = dfil_ca_topo.Histo1D({"ca_dca_om_dca_pi_topo","#Omega+#pi Candidate", 1000,-50000, 50000}, "om_dca_omXdca_pi_topo");
  
  //stra 
  auto om_dca_xy_stra = dfil_om_stra.Histo1D({"om_dca_xy_stra","#Omega", 200, -150, 150}, "fOmegaDCAxyToPVStraTrack");
  auto om_dca_z_stra = dfil_om_stra.Histo1D({"om_dca_z_stra","#Omega", 200, -150, 150}, "fOmegaDCAzToPVStraTrack");
  auto ca_dca_om_dca_pi_stra = dfil_ca_stra.Histo1D({"ca_dca_om_dca_pi_stra","#Omega+#pi Candidate", 1000,-50000, 50000}, "om_dca_omXdca_pi_stra");
  
  //Omega_c 
  //topo 
  auto om_c_dca_xy_topo = dfil_om_c_topo.Histo1D({"om_c_dca_xy_topo","#Omega_{c}#rightarrow#Omega+#pi", 200, -150, 150}, "fOmegaDCAxyToPVTopo");
  auto ca_c_dca_xy_topo = dfil_ca_c_topo.Histo1D({"ca_c_dca_xy_topo","#Omega_{c}#rightarrow#Omega+#pi Candiates", 200, -150, 150}, "fOmegaDCAxyToPVTopo");

  auto om_c_dca_z_topo = dfil_om_c_topo.Histo1D({"om_c_dca_z_topo","#Omega_{c}#rightarrow#Omega+#pi", 200, -150, 150}, "fOmegaDCAzToPVTopo");
  auto ca_c_dca_z_topo = dfil_ca_c_topo.Histo1D({"ca_c_dca_z_topo","#Omega_{c}#rightarrow#Omega+#pi Candiates", 200, -150, 150}, "fOmegaDCAzToPVTopo");

  auto om_c_dca_om_dca_pi_topo = dfil_om_c_topo.Histo1D({"om_c_dca_om_dca_pi_topo","#Omega_{c}#rightarrow#Omega+#pi", 1000,-50000, 50000}, "om_c_dca_omXdca_pi_topo");
  
  auto ca_c_dca_om_dca_pi_topo = dfil_ca_c_topo.Histo1D({"ca_c_dca_om_dca_pi_topo","#Omega_{c}#rightarrow#Omega+#pi Candiates", 1000,-50000, 50000}, "om_c_dca_omXdca_pi_topo");

  auto om_c_mass_topo = dfil_om_c_topo_dca.Histo1D({"om_c_mass_topo","#Omega_{c}#rightarrow#Omega+#pi", 200, 2.4, 3.0}, "fOmegacMassTopo"); 
  auto ca_c_mass_topo = dfil_ca_c_topo_dca.Histo1D({"ca_c_mass_topo","#Omega_{c}#rightarrow#Omega+#pi Candidates", 200, 2.4, 3.0}, "fOmegacMassTopo"); 
  
  auto om_c_rad_topo = dfil_om_c_topo.Histo1D({"om_c_rad_topo", "#Omega_{c}#rightarrow#Omega+#pi", 2000, 0, 0.25}, "fOmegacDecayRadiusTopo"); 
  auto ca_c_rad_topo = dfil_ca_c_topo.Histo1D({"ca_c_rad_topo", "#Omega_{c}#rightarrow#Omega+#pi Candidates", 2000, 0, 0.25}, "fOmegacDecayRadiusTopo"); 
  
  auto om_c_dist_topo = dfil_om_c_topo.Histo1D({"om_c_dist_topo", "#Omega_{c}#rightarrow#Omega+#pi", 2000, 0, 0.25}, "fOmegacDecayDistanceFromPVTopo"); 
  auto ca_c_dist_topo = dfil_ca_c_topo.Histo1D({"ca_c_dist_topo", "#Omega_{c}#rightarrow#Omega+#pi Candidates", 2000, 0, 0.25}, "fOmegacDecayDistanceFromPVTopo");
  
  auto om_c_daug_dca_topo = dfil_om_c_topo.Histo1D({"om_c_daug_dca_topo", "#Omega_{c}#rightarrow#Omega+#pi", 1000, 0, 400}, "fOmegacDaughterDCATopo"); 
  auto ca_c_daug_dca_topo = dfil_ca_c_topo.Histo1D({"ca_c_daug_dca_topo", "#Omega_{c}#rightarrow#Omega+#pi Candidates", 1000, 0, 400}, "fOmegacDaughterDCATopo");
  
  auto om_c_rad_dist_topo = dfil_om_c_topo.Histo2D({"om_c_rad_dist_topo", "#Omega_{c}#rightarrow#Omega+#pi", 500, 0, 0.25, 500, 0, 0.25}, "fOmegacDecayRadiusTopo", "fOmegacDecayDistanceFromPVTopo"); 
  auto ca_c_rad_dist_topo = dfil_ca_c_topo.Histo2D({"ca_c_rad_dist_topo", "#Omega_{c}#rightarrow#Omega+#pi", 500, 0, 0.25, 500, 0, 0.25}, "fOmegacDecayRadiusTopo", "fOmegacDecayDistanceFromPVTopo"); 
  
  auto om_c_rad_daug_dca_topo = dfil_om_c_topo.Histo2D({"om_c_rad_daug_dca_topo", "#Omega_{c}#rightarrow#Omega+#pi", 500, 0, 0.25, 500, 0, 400}, "fOmegacDecayRadiusTopo", "fOmegacDaughterDCATopo"); 
  auto ca_c_rad_daug_dca_topo = dfil_ca_c_topo.Histo2D({"ca_c_rad_daug_dca_topo", "#Omega_{c}#rightarrow#Omega+#pi", 500, 0, 0.25, 500, 0, 400}, "fOmegacDecayRadiusTopo", "fOmegacDaughterDCATopo"); 
  
  auto om_c_dist_daug_dca_topo = dfil_om_c_topo.Histo2D({"om_c_dist_daug_dca_topo", "#Omega_{c}#rightarrow#Omega+#pi", 500, 0, 0.25, 500, 0, 400}, "fOmegacDecayDistanceFromPVTopo", "fOmegacDaughterDCATopo"); 
  auto ca_c_dist_daug_dca_topo = dfil_ca_c_topo.Histo2D({"ca_c_dist_daug_dca_topo", "#Omega_{c}#rightarrow#Omega+#pi", 500, 0, 0.25, 500, 0, 400}, "fOmegacDecayDistanceFromPVTopo", "fOmegacDaughterDCATopo"); 
  
  //stra
  auto om_c_dca_xy_stra = dfil_om_c_stra.Histo1D({"om_c_dca_xy_stra","#Omega_{c}#rightarrow#Omega+#pi", 200, -150, 150}, "fOmegaDCAxyToPVStraTrack");
  auto ca_c_dca_xy_stra = dfil_ca_c_stra.Histo1D({"ca_c_dca_xy_stra","#Omega_{c}#rightarrow#Omega+#pi Candiates", 200, -150, 150}, "fOmegaDCAxyToPVStraTrack");
  
  auto om_c_dca_z_stra = dfil_om_c_stra.Histo1D({"om_c_dca_z_stra","#Omega_{c}#rightarrow#Omega+#pi", 200, -150, 150}, "fOmegaDCAzToPVStraTrack");
  auto ca_c_dca_z_stra = dfil_ca_c_stra.Histo1D({"ca_c_dca_z_stra","#Omega_{c}#rightarrow#Omega+#pi Candiates", 200, -150, 150}, "fOmegaDCAzToPVStraTrack");
  
  auto om_c_dca_om_dca_pi_stra = dfil_om_c_stra.Histo1D({"om_c_dca_om_dca_pi_stra","#Omega_{c}#rightarrow#Omega+#pi", 1000,-50000, 50000}, "om_c_dca_omXdca_pi_stra");
  auto ca_c_dca_om_dca_pi_stra = dfil_ca_c_stra.Histo1D({"ca_c_dca_om_dca_pi_stra","#Omega_{c}#rightarrow#Omega+#pi Candiates", 1000,-50000, 50000}, "om_c_dca_omXdca_pi_stra");

  auto om_c_mass_stra = dfil_om_c_stra_dca.Histo1D({"om_c_mass_stra","#Omega_{c}#rightarrow#Omega+#pi", 200, 2.4, 3.0}, "fOmegacMassStraTrack"); 
  auto ca_c_mass_stra = dfil_ca_c_stra_dca.Histo1D({"ca_c_mass_stra","#Omega_{c}#rightarrow#Omega+#pi Candidates", 200, 2.4, 3.0}, "fOmegacMassStraTrack"); 
  
  auto om_c_rad_stra = dfil_om_c_stra.Histo1D({"om_c_rad_stra", "#Omega_{c}#rightarrow#Omega+#pi", 2000, 0, 2.5}, "fOmegacDecayRadiusStraTrack"); 
  auto ca_c_rad_stra = dfil_ca_c_stra.Histo1D({"ca_c_rad_stra", "#Omega_{c}#rightarrow#Omega+#pi Candidates", 2000, 0, 2.5}, "fOmegacDecayRadiusStraTrack"); 
  
  auto om_c_dist_stra = dfil_om_c_stra.Histo1D({"om_c_dist_stra", "#Omega_{c}#rightarrow#Omega+#pi", 2000, 0, 2.5}, "fOmegacDecayDistanceFromPVStraTrack"); 
  auto ca_c_dist_stra = dfil_ca_c_stra.Histo1D({"ca_c_dist_stra", "#Omega_{c}#rightarrow#Omega+#pi Candidates", 2000, 0, 2.5}, "fOmegacDecayDistanceFromPVStraTrack");
  
  auto om_c_daug_dca_stra = dfil_om_c_stra.Histo1D({"om_c_daug_dca_stra", "#Omega_{c}#rightarrow#Omega+#pi", 1000, 0, 400}, "fOmegacDaughterDCAStraTrack"); 
  auto ca_c_daug_dca_stra = dfil_ca_c_stra.Histo1D({"ca_c_daug_dca_stra", "#Omega_{c}#rightarrow#Omega+#pi Candidates", 1000, 0, 400}, "fOmegacDaughterDCAStraTrack");
    
  auto om_c_rad_dist_stra = dfil_om_c_stra.Histo2D({"om_c_rad_dist_stra", "#Omega_{c}#rightarrow#Omega+#pi", 500, 0, 0.25, 500, 0, 0.25}, "fOmegacDecayRadiusStraTrack", "fOmegacDecayDistanceFromPVStraTrack");  
  auto ca_c_rad_dist_stra = dfil_ca_c_stra.Histo2D({"ca_c_rad_dist_stra", "#Omega_{c}#rightarrow#Omega+#pi", 500, 0, 0.25, 500, 0, 0.25}, "fOmegacDecayRadiusStraTrack", "fOmegacDecayDistanceFromPVStraTrack"); 
  
  auto om_c_rad_daug_dca_stra = dfil_om_c_stra.Histo2D({"om_c_rad_daug_dca_stra", "#Omega_{c}#rightarrow#Omega+#pi", 500, 0, 0.25, 500, 0, 400}, "fOmegacDecayRadiusStraTrack", "fOmegacDaughterDCAStraTrack"); 
  auto ca_c_rad_daug_dca_stra = dfil_ca_c_stra.Histo2D({"ca_c_rad_daug_dca_stra", "#Omega_{c}#rightarrow#Omega+#pi", 500, 0, 0.25, 500, 0, 400}, "fOmegacDecayRadiusStraTrack", "fOmegacDaughterDCAStraTrack"); 
  
  auto om_c_dist_daug_dca_stra = dfil_om_c_stra.Histo2D({"om_c_dist_daug_dca_stra", "#Omega_{c}#rightarrow#Omega+#pi", 500, 0, 0.25, 500, 0, 400}, "fOmegacDecayDistanceFromPVStraTrack", "fOmegacDaughterDCAStraTrack"); 
  auto ca_c_dist_daug_dca_stra = dfil_ca_c_stra.Histo2D({"ca_c_dist_daug_dca_stra", "#Omega_{c}#rightarrow#Omega+#pi", 500, 0, 0.25, 500, 0, 400}, "fOmegacDecayDistanceFromPVStraTrack", "fOmegacDaughterDCAStraTrack"); 
    
  TFile *out = TFile::Open("out.root", "recreate"); 
  auto NormalizeToMaximum = [] (ROOT::RDF::RResultPtr<TH1D> hist) { hist->Scale(1./hist->GetMaximum());};
  auto NormalizeToEntries = [] (ROOT::RDF::RResultPtr<TH1D> hist) { hist->Scale(1./hist->GetEntries());};
  auto Normalize2DToEntries = [] (ROOT::RDF::RResultPtr<TH2D> hist) { hist->Scale(1./hist->GetEntries());};
  auto DoNotNormalize = [] (ROOT::RDF::RResultPtr<TH1D> hist) { };
  
  //Just some dca distributions 
  
  PlotAndStore("dca_{om,xy}", out, om_dca_xy_topo, dummy, om_dca_xy_stra, dummy, NormalizeToEntries); 
  PlotAndStore("dca_{om,z}", out, om_dca_z_topo, dummy, om_dca_z_stra, dummy, NormalizeToMaximum); 
  
  PlotAndStore("dca_{omc,xy}", out, om_c_dca_xy_topo, ca_c_dca_xy_topo, om_c_dca_xy_stra, ca_c_dca_xy_stra, NormalizeToMaximum); 
  PlotAndStore("dca_{omc,z}", out, om_c_dca_z_topo, ca_c_dca_z_topo, om_c_dca_z_stra, ca_c_dca_z_stra, NormalizeToMaximum); 
  
  PlotAndStore("dca_{om vs omc,xy}", out, om_dca_xy_stra, dummy, om_c_dca_xy_stra, ca_c_dca_xy_stra, DoNotNormalize); 
  PlotAndStore("dca_{om vs omc,z}", out, om_dca_z_stra, dummy, om_c_dca_z_stra, ca_c_dca_z_stra, DoNotNormalize);
  
  //Decay related distributions

  PlotAndStore("om_vs_ca_c_radius_topo", out, om_c_rad_topo, ca_c_rad_topo, dummy, dummy, NormalizeToEntries); 
  PlotAndStore("om_vs_ca_c_dist_topo", out, om_c_dist_topo, ca_c_dist_topo, dummy, dummy, NormalizeToEntries); 
  PlotAndStore("om_vs_ca_c_daugh_dca_topo", out, om_c_daug_dca_topo, ca_c_daug_dca_topo, dummy, dummy, NormalizeToEntries); 
  
  PlotAndStore("om_vs_ca_c_radius_stra", out, om_c_rad_stra, ca_c_rad_stra, dummy, dummy, NormalizeToEntries); 
  PlotAndStore("om_vs_ca_c_dist_stra", out, om_c_dist_stra, ca_c_dist_stra, dummy, dummy, NormalizeToEntries); 
  PlotAndStore("om_vs_ca_c_daugh_dca_stra", out, om_c_daug_dca_stra, ca_c_daug_dca_stra, dummy, dummy, NormalizeToEntries); 
  
  CumulatePlotAndStore("cu_om_vs_ca_c_radius_topo", out, om_c_rad_topo, ca_c_rad_topo, dummy, dummy, DoNotNormalize); 
  CumulatePlotAndStore("cu_om_vs_ca_c_dist_topo", out, om_c_dist_topo, ca_c_dist_topo, dummy, dummy, DoNotNormalize); 
  CumulatePlotAndStore("cu_om_vs_ca_c_daugh_dca_topo", out, om_c_daug_dca_topo, ca_c_daug_dca_topo, dummy, dummy, DoNotNormalize); 
  
  CumulatePlotAndStore("cu_om_vs_ca_c_radius_stra", out, om_c_rad_stra, ca_c_rad_stra, dummy, dummy, DoNotNormalize); 
  CumulatePlotAndStore("cu_om_vs_ca_c_dist_stra", out, om_c_dist_stra, ca_c_dist_stra, dummy, dummy, DoNotNormalize); 
  CumulatePlotAndStore("cu_om_vs_ca_c_daugh_dca_stra", out, om_c_daug_dca_stra, ca_c_daug_dca_stra, dummy, dummy, DoNotNormalize); 
  
  NormalizeAndStore(out,om_c_rad_dist_topo, Normalize2DToEntries); 
  NormalizeAndStore(out,ca_c_rad_dist_topo, Normalize2DToEntries); 
  NormalizeAndStore(out,om_c_rad_daug_dca_topo , Normalize2DToEntries); 
  NormalizeAndStore(out,ca_c_rad_daug_dca_topo , Normalize2DToEntries); 
  NormalizeAndStore(out,om_c_dist_daug_dca_topo, Normalize2DToEntries); 
  NormalizeAndStore(out,ca_c_dist_daug_dca_topo, Normalize2DToEntries); 
  
  NormalizeAndStore(out,om_c_rad_dist_stra, Normalize2DToEntries); 
  NormalizeAndStore(out,ca_c_rad_dist_stra, Normalize2DToEntries); 
  NormalizeAndStore(out,om_c_rad_daug_dca_stra , Normalize2DToEntries); 
  NormalizeAndStore(out,ca_c_rad_daug_dca_stra , Normalize2DToEntries); 
  NormalizeAndStore(out,om_c_dist_daug_dca_stra, Normalize2DToEntries); 
  NormalizeAndStore(out,ca_c_dist_daug_dca_stra, Normalize2DToEntries); 
  
  //D Meson like cut
  
  PlotAndStore("om_dca_omXdca_pi", out, ca_dca_om_dca_pi_stra, dummy, ca_dca_om_dca_pi_topo, dummy, NormalizeToEntries); 
  PlotAndStore("om_cdca_omXdca_pi", out, om_c_dca_om_dca_pi_stra, ca_c_dca_om_dca_pi_stra, om_c_dca_om_dca_pi_topo,ca_c_dca_om_dca_pi_topo, NormalizeToEntries); 
  
  PlotAndStore("om_vs_om_c_dca_omXdca_pi_topo", out, om_c_dca_om_dca_pi_topo, ca_c_dca_om_dca_pi_topo, dummy, ca_dca_om_dca_pi_topo, DoNotNormalize); 
  PlotAndStore("om_vs_om_c_dca_omXdca_pi_stra", out, om_c_dca_om_dca_pi_stra, ca_c_dca_om_dca_pi_stra, dummy, ca_dca_om_dca_pi_stra, DoNotNormalize); 
  
  CumulatePlotAndStore("cu_om_vs_om_c_dca_omXdca_pi_topo", out, om_c_dca_om_dca_pi_topo, ca_c_dca_om_dca_pi_topo, dummy, ca_dca_om_dca_pi_topo, DoNotNormalize); 
  CumulatePlotAndStore("cu_om_vs_om_c_dca_omXdca_pi_stra", out, om_c_dca_om_dca_pi_stra, ca_c_dca_om_dca_pi_stra, dummy, ca_dca_om_dca_pi_stra, DoNotNormalize); 
  

  //Invariant mass 
  PlotAndStore("mass_{omc stra}", out, om_c_mass_stra, ca_c_mass_stra, dummy, dummy, NormalizeToEntries); 
  PlotAndStore("mass_{omc topo}", out, om_c_mass_topo, ca_c_mass_topo, dummy, dummy, NormalizeToEntries); 

  out->Close(); 
} 

