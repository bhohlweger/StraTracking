
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
  TString filePath = TString::Format("%s",HarryPlotter::FilePath()); 
  TFile *file = new TFile(filePath+"/omega.root", "READ");
  TFile *file_c = new TFile(filePath+"/omegac.root", "READ");

  ROOT::RDataFrame df_om("fTreeOmega", file);
  ROOT::RDataFrame df_ca("fTreeCandidates", file);

  ROOT::RDataFrame df_om_c("fTreeOmega", file_c);
  ROOT::RDataFrame df_ca_c("fTreeCandidates", file_c);

  //add classic quality cuts ... 
  //e.g. fOmegaPtMC>1.0&&fOmegaPtMC<3.0 
  //Strangeness tracking: fOmegaHitsAdded >= 3
  std::vector<TString> out_names = {"_integrated"};  
  
  ROOT::RDF::RResultPtr<TH1D> dummy;  
  std::vector<ROOT::RDF::RResultPtr<::TH2D>*> res_dca_xy_stra;  
  std::vector<ROOT::RDF::RResultPtr<::TH2D>*> res_dca_z_stra;  
  std::vector<ROOT::RDF::RResultPtr<::TH1D>*> res_dca_xy_topo;  
  std::vector<ROOT::RDF::RResultPtr<::TH1D>*> res_dca_z_topo;  
  
  //define filter 
  auto dfil_om_topo = df_om_c.Filter(HarryPlotter::TopoCuts_om_c, {"fOmegaPtMC", "fOmegacDecayRadiusTopo"});
  auto dfil_ca_topo = df_ca_c.Filter(HarryPlotter::TopoCuts_ca_c, {"fOmegaPtMC", "fOmegacDecayRadiusTopo"}).Filter("!fTrueOmegac"); 

  auto dfil_om_stra = df_om_c.Filter(HarryPlotter::StraCuts_om_c, {"fOmegaPtMC", "fOmegaHitsAdded", "fOmegacDecayRadiusStraTrack"}); 
  auto dfil_ca_stra = df_ca_c.Filter(HarryPlotter::StraCuts_ca_c, {"fOmegaPtMC", "fOmegaHitsAdded", "fOmegacDecayRadiusStraTrack"}).Filter("!fTrueOmegac"); 
  
  auto pv_res = dfil_om_stra.Define("delta_pv", HarryPlotter::Distance, {"fPrimaryVertexX", "fPrimaryVertexY", "fPrimaryVertexZ"});  
  auto res_pv = pv_res.Histo1D({"res_pv", "res_pv", 400, 0, 200}, "delta_pv");  

  res_dca_xy_stra.emplace_back(new ROOT::RDF::RResultPtr<::TH2D>(dfil_om_stra.Histo2D({"dca_xy_vs_AddedHits_stra", "dca_xy_vs_AddedHits_stra", 10, 0, 10, 400, -200, 200 }, "fOmegaHitsAdded", "fOmegaDCAxyToPVStraTrack"))); 
  res_dca_z_stra.emplace_back(new ROOT::RDF::RResultPtr<::TH2D>(dfil_om_stra.Histo2D({"dca_z_vs_AddedHits_stra", "dca_z_vs_AddedHits_stra", 10, 0 , 10, 400, -200, 200 }, "fOmegaHitsAdded", "fOmegaDCAzToPVStraTrack"))); 
  
  res_dca_xy_topo.emplace_back(new ROOT::RDF::RResultPtr<::TH1D>(dfil_om_topo.Histo1D({"dcay_xy_topo", "dca_xy_topo", 400, -200, 200} , "fOmegaDCAxyToPVTopo"))); 
  res_dca_z_topo.emplace_back(new ROOT::RDF::RResultPtr<::TH1D>(dfil_om_topo.Histo1D({"dcay_z_topo", "dca_z_topo", 400, -200, 200} , "fOmegaDCAzToPVTopo"))); 
  
  
  for (int ipt = 0; ipt < ptbins.size()-1; ++ipt) { 
    TString binstring = TString::Format("_pT_%.1f_%.1f", ptbins[ipt], ptbins[ipt+1]); 
    TString cutstring = TString::Format("(%.1f < fOmegaPtMC) && (fOmegaPtMC < %.1f)", ptbins[ipt], ptbins[ipt+1]); 
    out_names.push_back(binstring); 
    
    auto dv_res_stra_pt_layer = dfil_om_stra.Filter(cutstring.Data()); 
    res_dca_xy_stra.emplace_back(new ROOT::RDF::RResultPtr<::TH2D>(dv_res_stra_pt_layer.Histo2D({"dca_xy_vs_AddedHits_stra"+binstring, "dca_xy_vs_AddedHits_stra"+binstring, 10, 0, 10, 400, -200, 200 }, "fOmegaHitsAdded", "fOmegaDCAxyToPVStraTrack"))); 
    res_dca_z_stra.emplace_back(new ROOT::RDF::RResultPtr<::TH2D>(dv_res_stra_pt_layer.Histo2D({"dca_z_vs_AddedHits_stra"+binstring, "dca_z_vs_AddedHits_stra"+binstring, 10, 0, 10, 400, -200, 200 }, "fOmegaHitsAdded", "fOmegaDCAzToPVStraTrack"))); 
    
    auto dv_res_topo_pt_layer = dfil_om_topo.Filter(cutstring.Data()); 
    res_dca_xy_topo.emplace_back(new ROOT::RDF::RResultPtr<::TH1D>(dv_res_topo_pt_layer.Histo1D({"dca_xy_topo"+binstring, "dca_xy_topo"+binstring, 400, -200, 200 }, "fOmegaDCAxyToPVTopo"))); 
    res_dca_z_topo.emplace_back(new ROOT::RDF::RResultPtr<::TH1D>(dv_res_topo_pt_layer.Histo1D({"dca_z_topo"+binstring, "dca_z_topo"+binstring, 400, -200, 200 }, "fOmegaDCAzToPVTopo"))); 
    
  }

  //Write histos 
  TFile *out = TFile::Open("outResolutionsDCA.root", "recreate"); 
  for (int ipt = 0; ipt < ptbins.size(); ++ipt) { 
    
    TH1F* width_dca_xy = new TH1F("width_dca_xy"+out_names[ipt], "Width DCA_{xy}", 11, -1, 10); 
    width_dca_xy->GetYaxis()->SetTitle("#sigma_{DCAxy} (#mum)"); 
    TH1F* width_dca_z = new TH1F("width_dca_z" + out_names[ipt], "Width DCA_{z}", 11, -1, 10);  
    width_dca_z->GetYaxis()->SetTitle("#sigma_{DCAz} (#mum)"); 
    
    HarryPlotter::Normalize2DBinByBin(res_dca_xy_stra[ipt]->GetPtr()); 
    HarryPlotter::Normalize2DBinByBin(res_dca_z_stra[ipt]->GetPtr()); 

    HarryPlotter::CheckAndStore(out, res_dca_xy_stra[ipt]->GetPtr()); 
    HarryPlotter::CheckAndStore(out, res_dca_z_stra[ipt]->GetPtr()); 
    
    res_dca_xy_topo[ipt]->GetPtr()->Scale(1./res_dca_xy_topo[ipt]->GetPtr()->Integral()); 
    res_dca_z_topo[ipt]->GetPtr()->Scale(1./res_dca_z_topo[ipt]->GetPtr()->Integral()); 
    
    HarryPlotter::CheckAndStore(out, res_dca_xy_topo[ipt]->GetPtr()); 
    HarryPlotter::CheckAndStore(out, res_dca_z_topo[ipt]->GetPtr()); 

    for (int iBinX = 1; iBinX <= res_dca_xy_stra[ipt]->GetPtr()->GetNbinsX(); ++iBinX) { 
      res_dca_xy_stra[ipt]->GetPtr()->GetXaxis()->SetRange(iBinX, iBinX); 
      width_dca_xy->SetBinContent(iBinX+1, res_dca_xy_stra[ipt]->GetPtr()->GetStdDev(2)); 
      width_dca_xy->GetXaxis()->SetBinLabel(iBinX+1, TString::Format("%d added Hits", iBinX).Data()); 
    }      
    width_dca_xy->SetBinContent(1, res_dca_xy_topo[ipt]->GetPtr()->GetStdDev(1)); 
    width_dca_xy->GetXaxis()->SetBinLabel(1, "Topological Tracking"); 
    
    for (int iBinX = 1; iBinX <= res_dca_z_stra[ipt]->GetPtr()->GetNbinsX(); ++iBinX) { 
      res_dca_z_stra[ipt]->GetPtr()->GetXaxis()->SetRange(iBinX, iBinX); 
      width_dca_z->SetBinContent(iBinX+1, res_dca_z_stra[ipt]->GetPtr()->GetStdDev(2)); 
      width_dca_z->GetXaxis()->SetBinLabel(iBinX+1, TString::Format("%d added Hits", iBinX).Data()); 
    }
    width_dca_z->SetBinContent(1, res_dca_z_topo[ipt]->GetPtr()->GetStdDev(1)); 
    width_dca_z->GetXaxis()->SetBinLabel(1, "Topological Tracking"); 
    
    out->cd(); 
    width_dca_xy->Write(); 
    width_dca_z->Write(); 
  } 
  out->Close(); 
} 

