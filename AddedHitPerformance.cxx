
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

  ROOT::RDF::RResultPtr<TH1D> dummy;  
  std::vector<ROOT::RDF::RResultPtr<::TH1D>*> layers_Stra;  

  auto dfil_om_c_topo = df_om_c.Filter(HarryPlotter::TopoCuts_om_c, {"fOmegaPtMC", "fOmegacDecayRadiusTopo"});
  auto dfil_ca_c_topo = df_ca_c.Filter(HarryPlotter::TopoCuts_ca_c, {"fOmegaPtMC", "fOmegacDecayRadiusTopo"}); 

  auto dfv_res_topo = dfil_om_c_topo.Define("dlength_true", HarryPlotter::Distance, {"fOmegacVertexXMC", "fOmegacVertexYMC", "fOmegacVertexZMC"}).Define("delta_dlength", HarryPlotter::RelErr, {"dlength_true", "fOmegacDecayDistanceFromPVTopo"});  

  auto dfil_om_c_stra = df_om_c.Filter(HarryPlotter::StraCuts_om_c, {"fOmegaPtMC", "fOmegaHitsAdded", "fOmegacDecayRadiusStraTrack"}); 
  auto dfil_ca_c_stra = df_ca_c.Filter(HarryPlotter::StraCuts_ca_c, {"fOmegaPtMC", "fOmegaHitsAdded", "fOmegacDecayRadiusStraTrack"}); 

  auto pv_res = dfil_om_c_stra.Define("delta_pv", HarryPlotter::Distance, {"fPrimaryVertexX", "fPrimaryVertexY", "fPrimaryVertexZ"}); 
  auto h_pv_res = pv_res.Histo1D({"pv_resolution", "PV resolution", 150, 0, 3}, "delta_pv"); 

  auto dfv_res_stra = dfil_om_c_stra.Define("dlength_true", HarryPlotter::Distance, {"fOmegacVertexXMC", "fOmegacVertexYMC", "fOmegacVertexZMC"}).Define("delta_dlength", HarryPlotter::RelErr, {"dlength_true", "fOmegacDecayDistanceFromPVStraTrack"});  

  //begin loop over pT bins and Radius bins to define histograms 
  for (int ipt = 0; ipt < ptbins.size()-1; ++ipt) { 
    for (int ilayer = 0; ilayer < layerPos.size()-1; ++ilayer) { 
      TString binstring = TString::Format("_pT_%.1f_%.1f_decay_before_layer_%d", ptbins[ipt], ptbins[ipt+1], ilayer+1); 
      TString cutstring = TString::Format("(%.1f < fOmegaPtMC) && (fOmegaPtMC < %.1f) && (%.1f < fOmegaDecayRadiusMC) && (fOmegaDecayRadiusMC < %.1f)", ptbins[ipt], ptbins[ipt+1], layerPos[ilayer], layerPos[ilayer+1]); 

      auto dfv_res_stra_pt_layer = dfv_res_stra.Filter(cutstring.Data()); 

      auto hit_perf = dfv_res_stra_pt_layer.Define("hit_perf", TString::Format("%d - fOmegaHitsAdded", ilayer).Data()); 
      layers_Stra.emplace_back(new ROOT::RDF::RResultPtr<::TH1D>(hit_perf.Histo1D({"hit_perf"+binstring, "hit_perf"+binstring, 11, -1 , 10 }, "hit_perf"))); 
    }
  } 

  //Write histos 
  TFile *out = TFile::Open("outAddedHitPerformance.root", "recreate"); 

  TH2D* AddedHitPerformance = nullptr; 
  TH2D* AddedHitPerformance_integrated = new TH2D("AddedHitPerformance_integrated", "Added Hit Performance", 12, -0.5, 11.5, 10, -0.5, 9.5); 
  AddedHitPerformance_integrated->GetYaxis()->SetTitle("Crossed Layers#minusAdded Hits"); 
  AddedHitPerformance_integrated->GetXaxis()->SetTitle("Crossed Layers"); 
  int ilayer =layerPos.size() -1; 
  int ipt = 0; 
  for ( auto it : layers_Stra) { 
    if (ilayer == layerPos.size() -1 ) { 
      if (AddedHitPerformance) { 
	out->cd(); 
	AddedHitPerformance->Write();
      }
      TString binstring = TString::Format("_pT_%.1f_%.1f", ptbins[ipt], ptbins[ipt+1]); 
      AddedHitPerformance = new TH2D("AddedHitPerformance"+binstring, "Added Hit Performance", 12, -0.5, 11.5, 10, -0.5, 9.5); 
      AddedHitPerformance->GetYaxis()->SetTitle("Crossed Layers#minusAdded Hits"); 
      AddedHitPerformance->GetXaxis()->SetTitle("Crossed Layers"); 
      ilayer = 0; 
      ipt++; 
    }
    it->GetPtr()->Scale(1./(it->GetPtr()->Integral())); 
    HarryPlotter::CheckAndStore(out, *it); 
    for (int ibin=2; ibin<=it->GetPtr()->GetNbinsX(); ++ibin) { 
      AddedHitPerformance->SetBinContent(ilayer+1, ibin-1, it->GetPtr()->GetBinContent(ibin)); 
      AddedHitPerformance_integrated->SetBinContent(ilayer+1, ibin-1, AddedHitPerformance_integrated->GetBinContent(ilayer+1, ibin-1) + it->GetPtr()->GetBinContent(ibin)); 
    }
    ilayer++; 
  }
  HarryPlotter::Normalize2DBinByBin(AddedHitPerformance_integrated);
  AddedHitPerformance->Write();
  AddedHitPerformance_integrated->Write();
  out->Close(); 
} 

