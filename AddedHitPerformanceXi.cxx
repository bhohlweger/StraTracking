
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

std::vector<float> ptbins = {1.,6.};//HarryPlotter::Getptbins(); 
std::vector<float> layerPos = HarryPlotter::GetposLayers();

int main(int argc, char **argv) {
  const char* fileName = argv[1]; 
  const char* outAddon = (argv[2])?argv[2]:""; 
  HarryPlotter::StyleBox(); 

  ROOT::EnableImplicitMT(); // Tell ROOT you want to go parallel          
  TString filePath = TString::Format("%s", fileName); 
  TFile *file = new TFile(filePath, "READ");
  
  ROOT::RDataFrame df_xi("fTreeXi", file);
  
  ROOT::RDF::RResultPtr<TH1D> dummy;  
  std::vector<ROOT::RDF::RResultPtr<::TH1D>*> layers_Stra;  

  auto dfil_xi_stra = df_xi.Filter(HarryPlotter::TopoCuts_om, {"fXiPtMC"}).Filter("fTrueXi");

  //begin loop over pT bins and Radius bins to define histograms 
  for (int ipt = 0; ipt < ptbins.size()-1; ++ipt) { 
    for (int ilayer = 0; ilayer < layerPos.size()-1; ++ilayer) { 
      TString binstring = TString::Format("_pT_%.1f_%.1f_decay_before_layer_%d", ptbins[ipt], ptbins[ipt+1], ilayer+1); 
      TString cutstring = TString::Format("(%.1f < fXiPtMC) && (fXiPtMC < %.1f) && (%.1f < fXiDecayRadiusMC) && (fXiDecayRadiusMC < %.1f)", ptbins[ipt], ptbins[ipt+1], layerPos[ilayer], layerPos[ilayer+1]); 

      auto dfv_res_stra_pt_layer = dfil_xi_stra.Filter(cutstring.Data()); 

      //auto hit_perf = dfv_res_stra_pt_layer.Define("hit_perf", TString::Format("%d - fXiHitsAdded", ilayer).Data()); 
      layers_Stra.emplace_back(new ROOT::RDF::RResultPtr<::TH1D>(dfv_res_stra_pt_layer.Histo1D({"hit_perf"+binstring, "hit_perf"+binstring, 11, -1 , 10 }, "fXiHitsAdded"))); 
    }
  } 

  //Write histos 
  TFile *out = TFile::Open("outAddedHitPerformanceXi.root", "recreate"); 

  TH2D* AddedHitPerformance = nullptr; 
  TH2D* AddedHitPerformance_integrated = new TH2D("AddedHitPerformance_integrated", "Added Hit Performance", 12, -0.5, 11.5, 10, -0.5, 9.5); 
  AddedHitPerformance_integrated->GetYaxis()->SetTitle("Added Hits"); 
  AddedHitPerformance_integrated->GetXaxis()->SetTitle("Crossed Layers"); 
  int ilayer =layerPos.size() -1; 
  int ipt = 0; 
  double norm = 0; 
  for ( auto it : layers_Stra) { 
    if (ilayer == layerPos.size() -1 ) { 
      if (AddedHitPerformance) { 
	out->cd(); 
	AddedHitPerformance->Write();
      }
      TString binstring = TString::Format("_pT_%.1f_%.1f", ptbins[ipt], ptbins[ipt+1]); 
      AddedHitPerformance = new TH2D("AddedHitPerformance"+binstring, "Added Hit Performance", 12, -0.5, 11.5, 10, -0.5, 9.5); 
      AddedHitPerformance->GetYaxis()->SetTitle("Added Hits"); 
      AddedHitPerformance->GetXaxis()->SetTitle("Crossed Layers"); 
      ilayer = 0; 
      ipt++; 
    }
    norm+=it->GetPtr()->Integral(); 
    HarryPlotter::CheckAndStore(out, *it); 
    for (int ibin=2; ibin<=it->GetPtr()->GetNbinsX(); ++ibin) { 
      AddedHitPerformance->SetBinContent(ilayer+1, ibin-1, it->GetPtr()->GetBinContent(ibin)); 
      AddedHitPerformance_integrated->SetBinContent(ilayer+1, ibin-1, AddedHitPerformance_integrated->GetBinContent(ilayer+1, ibin-1) + it->GetPtr()->GetBinContent(ibin)); 
    }
    ilayer++; 
  }
  //HarryPlotter::Normalize2DBinByBin(AddedHitPerformance_integrated);
  AddedHitPerformance_integrated->Scale(1./norm); 
  AddedHitPerformance->Write();
  AddedHitPerformance_integrated->Write();
  out->Close(); 
} 

