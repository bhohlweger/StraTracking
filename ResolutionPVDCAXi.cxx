
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
  const char* fileName = argv[1]; 
  const char* outAddon = (argv[2])?argv[2]:""; 
  HarryPlotter::StyleBox(); 

  ROOT::EnableImplicitMT(); // Tell ROOT you want to go parallel          
  TString filePath = TString::Format("%s", fileName); 
  TFile *file = new TFile(filePath, "READ");
  //TFile *file_c = new TFile(filePath+"/xiegac.root", "READ");

  ROOT::RDataFrame df_xi_c("fTreeXi", file);
  //ROOT::RDataFrame df_ca_c("fTreeCandidates", file);

  // ROOT::RDataFrame df_xi_c("fTreeXi", file_c);
  // ROOT::RDataFrame df_ca_c("fTreeCandidates", file_c);

  //add classic quality cuts ... 
  //e.g. fXiPtMC>1.0&&fXiPtMC<3.0 
  //Strangeness tracking: fXiHitsAdded >= 3
  std::vector<TString> out_names = {"_integrated"};  
  
  ROOT::RDF::RResultPtr<TH1D> dummy;  
  std::vector<ROOT::RDF::RResultPtr<::TH2D>*> res_dca_xy_stra;  
  std::vector<ROOT::RDF::RResultPtr<::TH2D>*> res_dca_z_stra;  
  std::vector<ROOT::RDF::RResultPtr<::TH1D>*> res_dca_xy_topo;  
  std::vector<ROOT::RDF::RResultPtr<::TH1D>*> res_dca_z_topo;  
  
  //define filter 
  auto dfil_xi_topo = df_xi_c.Filter(HarryPlotter::TopoCuts_om, {"fXiPtMC"}).Filter("fTrueXi"); 
  //auto dfil_ca_topo = df_ca_c.Filter(HarryPlotter::TopoCuts_ca, {"fXiPtMC", "fXicDecayRadiusTopo"}).Filter("!fTrueXic"); 
  auto dfil_xi_topo_bench = dfil_xi_topo.Filter("fXiDecayRadiusMC < 0.5"); 

  auto dfil_xi_stra = df_xi_c.Filter(HarryPlotter::TopoCuts_om, {"fXiPtMC"}).Filter("fTrueXi"); 
  //auto dfil_ca_stra = df_ca_c.Filter(HarryPlotter::StraCuts_ca_c, {"fXiPtMC", "fXiHitsAdded", "fXicDecayRadiusStraTrack"}).Filter("!fTrueXic"); 
  
  auto pv_res = dfil_xi_stra.Define("delta_pv", HarryPlotter::Distance, {"fPrimaryVertexX", "fPrimaryVertexY", "fPrimaryVertexZ"});  
  auto res_pv = pv_res.Histo1D({"res_pv", "res_pv", 400, 0, 200}, "delta_pv");  

  res_dca_xy_stra.emplace_back(new ROOT::RDF::RResultPtr<::TH2D>(dfil_xi_stra.Histo2D({"dca_xy_vs_AddedHits_stra", "dca_xy_vs_AddedHits_stra", 10, 0, 10, 400, -200, 200 }, "fXiHitsAdded", "fXiDCAxyToPVStraTrack"))); 
  res_dca_z_stra.emplace_back(new ROOT::RDF::RResultPtr<::TH2D>(dfil_xi_stra.Histo2D({"dca_z_vs_AddedHits_stra", "dca_z_vs_AddedHits_stra", 10, 0 , 10, 400, -200, 200 }, "fXiHitsAdded", "fXiDCAzToPVStraTrack"))); 
  
  res_dca_xy_topo.emplace_back(new ROOT::RDF::RResultPtr<::TH1D>(dfil_xi_topo.Histo1D({"dca_xy_topo", "dca_xy_topo", 400, -200, 200} , "fXiDCAxyToPVTopo"))); 
  res_dca_z_topo.emplace_back(new ROOT::RDF::RResultPtr<::TH1D>(dfil_xi_topo.Histo1D({"dca_z_topo", "dca_z_topo", 400, -200, 200} , "fXiDCAzToPVTopo"))); 
  
  res_dca_xy_topo.emplace_back(new ROOT::RDF::RResultPtr<::TH1D>(dfil_xi_topo_bench.Histo1D({"dca_xy_topo_bench", "dca_xy_topo", 400, -200, 200} , "fXiDCAxyToPVTopo"))); 
  res_dca_z_topo.emplace_back(new ROOT::RDF::RResultPtr<::TH1D>(dfil_xi_topo_bench.Histo1D({"dca_z_topo_bench", "dca_z_topo", 400, -200, 200} , "fXiDCAzToPVTopo"))); 
  
  for (int ipt = 0; ipt < ptbins.size()-1; ++ipt) { 
    TString binstring = TString::Format("_pT_%.1f_%.1f", ptbins[ipt], ptbins[ipt+1]); 
    TString cutstring = TString::Format("(%.1f < fXiPtMC) && (fXiPtMC < %.1f)", ptbins[ipt], ptbins[ipt+1]); 
    out_names.push_back(binstring); 
    
    auto dv_res_stra_pt_layer = dfil_xi_stra.Filter(cutstring.Data()); 
    res_dca_xy_stra.emplace_back(new ROOT::RDF::RResultPtr<::TH2D>(dv_res_stra_pt_layer.Histo2D({"dca_xy_vs_AddedHits_stra"+binstring, "dca_xy_vs_AddedHits_stra"+binstring, 10, 0, 10, 400, -200, 200 }, "fXiHitsAdded", "fXiDCAxyToPVStraTrack"))); 
    res_dca_z_stra.emplace_back(new ROOT::RDF::RResultPtr<::TH2D>(dv_res_stra_pt_layer.Histo2D({"dca_z_vs_AddedHits_stra"+binstring, "dca_z_vs_AddedHits_stra"+binstring, 10, 0, 10, 400, -200, 200 }, "fXiHitsAdded", "fXiDCAzToPVStraTrack"))); 
    
    auto dv_res_topo_pt_layer = dfil_xi_topo.Filter(cutstring.Data()); 
    res_dca_xy_topo.emplace_back(new ROOT::RDF::RResultPtr<::TH1D>(dv_res_topo_pt_layer.Histo1D({"dca_xy_topo"+binstring, "dca_xy_topo"+binstring, 400, -200, 200 }, "fXiDCAxyToPVTopo"))); 
    res_dca_z_topo.emplace_back(new ROOT::RDF::RResultPtr<::TH1D>(dv_res_topo_pt_layer.Histo1D({"dca_z_topo"+binstring, "dca_z_topo"+binstring, 400, -200, 200 }, "fXiDCAzToPVTopo"))); 
    
    auto dv_res_topo_pt_layer_bench = dfil_xi_topo_bench.Filter(cutstring.Data()); 
    res_dca_xy_topo.emplace_back(new ROOT::RDF::RResultPtr<::TH1D>(dv_res_topo_pt_layer_bench.Histo1D({"dca_xy_topo_bench"+binstring, "dca_xy_topo_bench"+binstring, 400, -200, 200 }, "fXiDCAxyToPVTopo"))); 
    res_dca_z_topo.emplace_back(new ROOT::RDF::RResultPtr<::TH1D>(dv_res_topo_pt_layer_bench.Histo1D({"dca_z_topo_bench"+binstring, "dca_z_topo_bench"+binstring, 400, -200, 200 }, "fXiDCAzToPVTopo"))); 

  }

  //Write histos 
  TFile *out = TFile::Open(TString::Format("outResolutionsDCA_xi%s.root", outAddon).Data(), "recreate"); 
  for (int ipt = 0; ipt < ptbins.size(); ++ipt) { 
  //for (int ipt = 0; ipt < 1; ++ipt) { 
    
    TH1F* width_dca_xy = new TH1F("width_dca_xy"+out_names[ipt], "Width DCA_{xy}", 12, -2, 10); 
    width_dca_xy->GetYaxis()->SetTitle("#sigma_{DCAxy} (#mum)"); 
    TH1F* width_dca_z = new TH1F("width_dca_z" + out_names[ipt], "Width DCA_{z}", 12, -2, 10);  
    width_dca_z->GetYaxis()->SetTitle("#sigma_{DCAz} (#mum)"); 
    
    // HarryPlotter::Normalize2DBinByBin(res_dca_xy_stra[ipt]->GetPtr()); 
    // HarryPlotter::Normalize2DBinByBin(res_dca_z_stra[ipt]->GetPtr()); 
    
    for (int iBinX = 1; iBinX <= res_dca_xy_stra[ipt]->GetPtr()->GetNbinsX(); ++iBinX) { 
      auto dca_xy_stra_proj = res_dca_xy_stra[ipt]->GetPtr()->ProjectionY(TString::Format("%s_%d", res_dca_xy_stra[ipt]->GetPtr()->GetName(), iBinX).Data(), iBinX, iBinX); 
      if(dca_xy_stra_proj->Integral() > 2000) { 
	width_dca_xy->SetBinContent(iBinX+2, HarryPlotter::FitDCA(dca_xy_stra_proj)); 
	HarryPlotter::CheckAndStore(out, dca_xy_stra_proj); 
      }
      width_dca_xy->GetXaxis()->SetBinLabel(iBinX+2, TString::Format("%d added Hits", iBinX-1).Data()); 
    }      
    
    width_dca_xy->SetBinContent(1, HarryPlotter::FitDCA(res_dca_xy_topo[2*ipt]->GetPtr())); 
    width_dca_xy->GetXaxis()->SetBinLabel(1, "Topological Tracking"); 
    
    width_dca_xy->SetBinContent(2, HarryPlotter::FitDCA(res_dca_xy_topo[2*ipt+1]->GetPtr())); 
    width_dca_xy->GetXaxis()->SetBinLabel(2, "Topo. Tra. (r_{Decay}<0.5cm)"); 

    for (int iBinX = 1; iBinX <= res_dca_z_stra[ipt]->GetPtr()->GetNbinsX(); ++iBinX) { 
      auto dca_z_stra_proj = res_dca_z_stra[ipt]->GetPtr()->ProjectionY(TString::Format("%s_%d", res_dca_z_stra[ipt]->GetPtr()->GetName(), iBinX).Data(), iBinX, iBinX); 
      if(dca_z_stra_proj->Integral() > 2000) { 
	width_dca_z->SetBinContent(iBinX+2, HarryPlotter::FitDCA(dca_z_stra_proj)); 
	HarryPlotter::CheckAndStore(out, dca_z_stra_proj); 
      }
      width_dca_z->GetXaxis()->SetBinLabel(iBinX+2, TString::Format("%d added Hits", iBinX-1).Data()); 
    }      
    
    width_dca_z->SetBinContent(1, HarryPlotter::FitDCA(res_dca_z_topo[2*ipt]->GetPtr()));     
    width_dca_z->GetXaxis()->SetBinLabel(1, "Topological Tracking"); 
    
    width_dca_z->SetBinContent(2, HarryPlotter::FitDCA(res_dca_z_topo[2*ipt+1]->GetPtr())); 
    width_dca_z->GetXaxis()->SetBinLabel(2, "Topo. Tra. (r_{Decay}<0.5cm)"); 
    
    out->cd(); 
    
    HarryPlotter::CheckAndStore(out, res_dca_xy_stra[ipt]->GetPtr()); 
    HarryPlotter::CheckAndStore(out, res_dca_z_stra[ipt]->GetPtr()); 
        
    HarryPlotter::CheckAndStore(out, res_dca_xy_topo[2*ipt]->GetPtr()); 
    HarryPlotter::CheckAndStore(out, res_dca_z_topo[2*ipt]->GetPtr());

    HarryPlotter::CheckAndStore(out, res_dca_xy_topo[2*ipt+1]->GetPtr()); 
    HarryPlotter::CheckAndStore(out, res_dca_z_topo[2*ipt+1]->GetPtr());

    width_dca_xy->Write(); 
    width_dca_z->Write(); 
  } 
  out->Close(); 
} 

