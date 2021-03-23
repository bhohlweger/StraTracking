
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
  
  std::vector<TString> out_names = {"_integrated"};  
  
  ROOT::RDF::RResultPtr<TH1D> dummy;  
  std::vector<ROOT::RDF::RResultPtr<::TH2D>*> res_dca_x_stra;  
  std::vector<ROOT::RDF::RResultPtr<::TH2D>*> res_dca_y_stra;  
  std::vector<ROOT::RDF::RResultPtr<::TH2D>*> res_dca_z_stra;  

  std::vector<ROOT::RDF::RResultPtr<::TH1D>*> res_dca_x_topo;  
  std::vector<ROOT::RDF::RResultPtr<::TH1D>*> res_dca_y_topo;  
  std::vector<ROOT::RDF::RResultPtr<::TH1D>*> res_dca_z_topo;  
  
  //define filter 
  auto dfil_om_topo = df_om_c.Filter(HarryPlotter::TopoCuts_om_c, {"fOmegacPtMC", "fOmegacDecayRadiusTopo"})
    .Define("dca_xy_topo", HarryPlotter::Distance_xy, {"fOmegacVertexXTopo","fOmegacVertexYTopo", "fOmegacVertexXMC", "fOmegacVertexYMC"})
    .Define("dca_x_topo", HarryPlotter::Distance_z, {"fOmegacVertexXTopo", "fOmegacVertexXMC"})
    .Define("dca_y_topo", HarryPlotter::Distance_z, {"fOmegacVertexYTopo", "fOmegacVertexYMC"})
    .Define("dca_z_topo", HarryPlotter::Distance_z, {"fOmegacVertexZTopo", "fOmegacVertexZMC"});
			  
  //auto dfil_ca_topo = df_ca.Filter(TopoCuts, {"fOmegaPtMC", "fOmegacDecayRadiusTopo"}).Filter("!fTrueOmegac"); 

  auto dfil_om_stra = df_om_c.Filter(HarryPlotter::StraCuts_om_c, {"fOmegacPtMC", "fOmegaHitsAdded", "fOmegacDecayRadiusStraTrack"})
    .Define("dca_x_stra", HarryPlotter::Distance_z, {"fOmegacVertexXStraTrack", "fOmegacVertexXMC"})
    .Define("dca_y_stra", HarryPlotter::Distance_z, {"fOmegacVertexYStraTrack", "fOmegacVertexYMC"})
    .Define("dca_z_stra", HarryPlotter::Distance_z, {"fOmegacVertexZStraTrack", "fOmegacVertexZMC"});   
  //auto dfil_ca_stra = df_ca.Filter(StraCuts, {"fOmegaPtMC", "fOmegaHitsAdded", "fOmegacDecayRadiusStraTrack"}).Filter("!fTrueOmegac"); 
  
  auto res_pv_x = dfil_om_stra.Histo1D({"res_pv_x", "res_pv_x", 400, -0.002, 0.002}, "fPrimaryVertexX");  
  auto res_pv_y = dfil_om_stra.Histo1D({"res_pv_y", "res_pv_y", 400, -0.002, 0.002}, "fPrimaryVertexY");  
  auto res_pv_z = dfil_om_stra.Histo1D({"res_pv_z", "res_pv_z", 400, -0.002, 0.002}, "fPrimaryVertexZ");  

  res_dca_x_stra.emplace_back(new ROOT::RDF::RResultPtr<::TH2D>(dfil_om_stra.Histo2D({"dca_x_vs_AddedHits_stra", "dca_x_vs_AddedHits_stra", 10, 0, 10, 400, -200e-4, 200e-4 }, "fOmegaHitsAdded", "dca_x_stra"))); 
  res_dca_y_stra.emplace_back(new ROOT::RDF::RResultPtr<::TH2D>(dfil_om_stra.Histo2D({"dca_y_vs_AddedHits_stra", "dca_y_vs_AddedHits_stra", 10, 0, 10, 400, -200e-4, 200e-4 }, "fOmegaHitsAdded", "dca_y_stra"))); 
  res_dca_z_stra.emplace_back(new ROOT::RDF::RResultPtr<::TH2D>(dfil_om_stra.Histo2D({"dca_z_vs_AddedHits_stra", "dca_z_vs_AddedHits_stra", 10, 0 , 10, 400, -200e-4, 200e-4 }, "fOmegaHitsAdded", "dca_z_stra"))); 
  
  res_dca_x_topo.emplace_back(new ROOT::RDF::RResultPtr<::TH1D>(dfil_om_topo.Histo1D({"dcay_x_topo", "dca_x_topo", 400, -200e-4, 200e-4} , "dca_x_topo"))); 
  res_dca_y_topo.emplace_back(new ROOT::RDF::RResultPtr<::TH1D>(dfil_om_topo.Histo1D({"dcay_y_topo", "dca_y_topo", 400, -200e-4, 200e-4} , "dca_y_topo"))); 
  res_dca_z_topo.emplace_back(new ROOT::RDF::RResultPtr<::TH1D>(dfil_om_topo.Histo1D({"dcay_z_topo", "dca_z_topo", 400, -200e-4, 200e-4} , "dca_z_topo"))); 
  
  
  for (int ipt = 0; ipt < ptbins.size()-1; ++ipt) { 
    TString binstring = TString::Format("_pT_%.1f_%.1f", ptbins[ipt], ptbins[ipt+1]); 
    TString cutstring = TString::Format("(%.1f < fOmegacPtMC) && (fOmegacPtMC < %.1f)", ptbins[ipt], ptbins[ipt+1]); 
    out_names.push_back(binstring); 
    
    auto dv_res_stra_pt_layer = dfil_om_stra.Filter(cutstring.Data()); 
    res_dca_x_stra.emplace_back(new ROOT::RDF::RResultPtr<::TH2D>(dv_res_stra_pt_layer.Histo2D({"dca_x_vs_AddedHits_stra"+binstring, "dca_x_vs_AddedHits_stra"+binstring, 10, 0, 10, 400, -200e-4, 200e-4 }, "fOmegaHitsAdded", "dca_x_stra"))); 
    res_dca_y_stra.emplace_back(new ROOT::RDF::RResultPtr<::TH2D>(dv_res_stra_pt_layer.Histo2D({"dca_y_vs_AddedHits_stra"+binstring, "dca_y_vs_AddedHits_stra"+binstring, 10, 0, 10, 400, -200e-4, 200e-4 }, "fOmegaHitsAdded", "dca_y_stra"))); 
    res_dca_z_stra.emplace_back(new ROOT::RDF::RResultPtr<::TH2D>(dv_res_stra_pt_layer.Histo2D({"dca_z_vs_AddedHits_stra"+binstring, "dca_z_vs_AddedHits_stra"+binstring, 10, 0, 10, 400, -200e-4, 200e-4 }, "fOmegaHitsAdded", "dca_z_stra"))); 
    
    auto dv_res_topo_pt_layer = dfil_om_topo.Filter(cutstring.Data()); 
    res_dca_x_topo.emplace_back(new ROOT::RDF::RResultPtr<::TH1D>(dv_res_topo_pt_layer.Histo1D({"dca_x_topo"+binstring, "dca_x_topo"+binstring, 400, -200e-4, 200e-4 }, "dca_x_topo"))); 
    res_dca_y_topo.emplace_back(new ROOT::RDF::RResultPtr<::TH1D>(dv_res_topo_pt_layer.Histo1D({"dca_y_topo"+binstring, "dca_y_topo"+binstring, 400, -200e-4, 200e-4 }, "dca_y_topo"))); 
    res_dca_z_topo.emplace_back(new ROOT::RDF::RResultPtr<::TH1D>(dv_res_topo_pt_layer.Histo1D({"dca_z_topo"+binstring, "dca_z_topo"+binstring, 400, -200e-4, 200e-4 }, "dca_z_topo"))); 
    
  }

  //Write histos 
  TFile *out = TFile::Open("outResolutionsDV.root", "recreate"); 
  for (int ipt = 0; ipt < ptbins.size(); ++ipt) { 
    TH1F* width_dca_x = new TH1F("width_dca_x"+out_names[ipt], "x", 11, -1, 10); 
    width_dca_x->GetYaxis()->SetTitle("#sigma_{#Omega_{c} DV} (#mum)"); 
    TH1F* width_dca_y = new TH1F("width_dca_y"+out_names[ipt], "y", 11, -1, 10); 
    width_dca_y->GetYaxis()->SetTitle("#sigma_{#Omega_{c} DV} (#mum)"); 
    TH1F* width_dca_z = new TH1F("width_dca_z" + out_names[ipt], "z", 11, -1, 10);  
    width_dca_z->GetYaxis()->SetTitle("#sigma_{#Omega_{c} DV} (#mum)"); 
    
       
    HarryPlotter::Normalize2DBinByBin(res_dca_x_stra[ipt]->GetPtr()); 
    HarryPlotter::Normalize2DBinByBin(res_dca_y_stra[ipt]->GetPtr()); 
    HarryPlotter::Normalize2DBinByBin(res_dca_z_stra[ipt]->GetPtr()); 

    HarryPlotter::CheckAndStore(out, res_dca_x_stra[ipt]->GetPtr()); 
    HarryPlotter::CheckAndStore(out, res_dca_y_stra[ipt]->GetPtr()); 
    HarryPlotter::CheckAndStore(out, res_dca_z_stra[ipt]->GetPtr()); 
    
    res_dca_x_topo[ipt]->GetPtr()->Scale(1./res_dca_x_topo[ipt]->GetPtr()->Integral()); 
    res_dca_y_topo[ipt]->GetPtr()->Scale(1./res_dca_y_topo[ipt]->GetPtr()->Integral()); 
    res_dca_z_topo[ipt]->GetPtr()->Scale(1./res_dca_z_topo[ipt]->GetPtr()->Integral()); 
    
    HarryPlotter::CheckAndStore(out, res_dca_x_topo[ipt]->GetPtr()); 
    HarryPlotter::CheckAndStore(out, res_dca_y_topo[ipt]->GetPtr()); 
    HarryPlotter::CheckAndStore(out, res_dca_z_topo[ipt]->GetPtr()); 
    
    for (int iBinX = 1; iBinX <= res_dca_x_stra[ipt]->GetPtr()->GetNbinsX(); ++iBinX) { 
      res_dca_x_stra[ipt]->GetPtr()->GetXaxis()->SetRange(iBinX, iBinX); 
      width_dca_x->SetBinContent(iBinX+1, res_dca_x_stra[ipt]->GetPtr()->GetStdDev(2)); 
      width_dca_x->GetXaxis()->SetBinLabel(iBinX+1, TString::Format("#Omega %d Hits", iBinX).Data()); 
      width_dca_x->SetBinError(iBinX+1, 0); 
    }      
    width_dca_x->SetBinContent(1, res_dca_x_topo[ipt]->GetPtr()->GetStdDev(1)); 
    width_dca_x->SetBinError(1, 0); 
    width_dca_x->GetXaxis()->SetBinLabel(1, "Topological Tracking"); 

    for (int iBinX = 1; iBinX <= res_dca_y_stra[ipt]->GetPtr()->GetNbinsX(); ++iBinX) { 
      res_dca_y_stra[ipt]->GetPtr()->GetXaxis()->SetRange(iBinX, iBinX); 
      width_dca_y->SetBinContent(iBinX+1, res_dca_y_stra[ipt]->GetPtr()->GetStdDev(2)); 
      width_dca_y->GetXaxis()->SetBinLabel(iBinX+1, TString::Format("#Omega %d Hits", iBinX).Data()); 
      width_dca_y->SetBinError(iBinX+1, 0); 
    }      
    width_dca_y->SetBinContent(1, res_dca_y_topo[ipt]->GetPtr()->GetStdDev(1)); 
    width_dca_y->SetBinError(1, 0);
    width_dca_y->GetXaxis()->SetBinLabel(1, "#Omega Topological"); 
    
    for (int iBinX = 1; iBinX <= res_dca_z_stra[ipt]->GetPtr()->GetNbinsX(); ++iBinX) { 
      res_dca_z_stra[ipt]->GetPtr()->GetXaxis()->SetRange(iBinX, iBinX); 
      width_dca_z->SetBinContent(iBinX+1, res_dca_z_stra[ipt]->GetPtr()->GetStdDev(2)); 
      width_dca_z->GetXaxis()->SetBinLabel(iBinX+1, TString::Format("#Omega %d Hits", iBinX).Data()); 
      width_dca_z->SetBinError(iBinX+1, 0); 
    }
    width_dca_z->SetBinContent(1, res_dca_z_topo[ipt]->GetPtr()->GetStdDev(1)); 
    width_dca_z->SetBinError(1, 0); 
    width_dca_z->GetXaxis()->SetBinLabel(1, "#Omega Topological"); 
    
    out->cd(); 
    width_dca_x->Scale(1./1e-4); 
    width_dca_y->Scale(1./1e-4); 
    width_dca_z->Scale(1./1e-4); 

    width_dca_x->Write(); 
    width_dca_y->Write(); 
    width_dca_z->Write(); 
  } 
  res_pv_x->Scale(1./res_pv_x->Integral());
  res_pv_y->Scale(1./res_pv_y->Integral());
  res_pv_z->Scale(1./res_pv_z->Integral());
  HarryPlotter::CheckAndStore(out, res_pv_x); 
  HarryPlotter::CheckAndStore(out, res_pv_y); 
  HarryPlotter::CheckAndStore(out, res_pv_z); 

  out->Close(); 
} 

