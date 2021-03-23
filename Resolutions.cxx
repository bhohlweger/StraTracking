
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


std::vector<float> ptbins = {1.0, 3.0}; //{1.0, 1.5, 2.0, 2.5, 3.0}; 
std::vector<float> layerPos = {0, 0.5, 1.2, 2.5, 3.75, 7.0, 12, 20, 30, 45, 60, 80, 100}; 

void StyleBox() { 

  const int NCont = 255;
  //gROOT->ForceStyle();                                                                                                                                                                                                   
  gStyle->Reset("Plain");
  gStyle->SetNumberContours(NCont);

  gStyle->SetTitleBorderSize(0);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(8, 0);
  gStyle->SetCanvasBorderMode(0);

  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.1);
  gStyle->SetPadTopMargin(0.15);
  gStyle->SetPadRightMargin(0.15);

  gStyle->SetFrameLineWidth(1);
  gStyle->SetHistLineWidth(2);
  gStyle->SetFuncWidth(2);
  gStyle->SetLineWidth(2);
  gStyle->SetHatchesLineWidth(1);
  gStyle->SetHatchesSpacing(0.5);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadColor(0);
  gStyle->SetHistLineColor(kBlack);
  gStyle->SetFuncColor(kGreen);
  gStyle->SetLabelColor(kBlack, "xyz");
  gStyle->SetPalette(kCividis);

  gStyle->SetTextFont(43);
  gStyle->SetLabelFont(43, "xyz");
  gStyle->SetTitleFont(43, "xyz");
  gStyle->SetLegendFont(43);

  gStyle->SetTextSizePixels(28);
  gStyle->SetLabelSize(28, "xyz");
  gStyle->SetTitleSize(28, "xyz");
  gStyle->SetLegendTextSize(28);

  gStyle->SetLabelOffset(0.01, "xy");
  gStyle->SetTitleOffset(1.25, "y");
  gStyle->SetTitleOffset(1.25, "x");

  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  gStyle->SetLegendBorderSize(0);

}

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
  TFile *file = new TFile("/localstore/alice/hohlweger/analysis/StrangnessTracking/210316/omega.root", "READ");
  TFile *file_c = new TFile("/localstore/alice/hohlweger/analysis/StrangnessTracking/210316/omegac.root", "READ");
  
  ROOT::RDataFrame df_om("fTreeOmega", file);
  ROOT::RDataFrame df_ca("fTreeCandidates", file);
  
  ROOT::RDataFrame df_om_c("fTreeOmega", file_c);
  ROOT::RDataFrame df_ca_c("fTreeCandidates", file_c);
  
  //add classic quality cuts ... 
  //e.g. fOmegaPtMC>1.0&&fOmegaPtMC<3.0 
  //Strangeness tracking: fOmegaHitsAdded >= 3
  auto TopoCuts = [] (float pT, float decayRad) {return ((1.0 < pT)&&(pT < 3.0)&&(decayRad>1e-12));}; 
  auto StraCuts = [] (float pT, int iAdded, float decayRad) {return ((1.0 < pT)&&(pT < 3.0)&&(decayRad>1e-12));}; 

  auto Length = [] (float x1, float y1, float z1, float x2, float y2, float z2) { 
    return (float)TMath::Sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
  }; 
    
  auto Distance = [] (float x1, float y1, float z1) {
    return (float)TMath::Sqrt(x1*x1+y1*y1+z1*z1);
  };  
	 
  auto RelErr = [] (float generated, float measured) {
    return (float)(generated - measured)/generated;
  }; 
  
  ROOT::RDF::RResultPtr<TH1D> dummy;  
  std::vector<ROOT::RDF::RResultPtr<::TH1D>*> resolutions_Stra;  
  std::vector<ROOT::RDF::RResultPtr<::TH1D>*> layers_Stra;  
  std::vector<ROOT::RDF::RResultPtr<::TH1D>*> resolutions_Topo;  
  

  //define filter 
			  
  //omega_c
  auto dfil_om_c_topo = df_om_c.Filter(TopoCuts, {"fOmegaPtMC", "fOmegacDecayRadiusTopo"});
  auto dfil_ca_c_topo = df_ca_c.Filter(TopoCuts, {"fOmegaPtMC", "fOmegacDecayRadiusTopo"}).Filter("!fTrueOmegac"); 
  auto dfv_res_topo = dfil_om_c_topo.Define("dlength_true", Distance, {"fOmegacVertexXMC", "fOmegacVertexYMC", "fOmegacVertexZMC"}).Define("delta_dlength", RelErr, {"dlength_true", "fOmegacDecayDistanceFromPVTopo"});  
  
  auto dfil_om_c_stra = df_om_c.Filter(StraCuts, {"fOmegaPtMC", "fOmegaHitsAdded", "fOmegacDecayRadiusStraTrack"}); 
  auto dfil_ca_c_stra = df_ca_c.Filter(StraCuts, {"fOmegaPtMC", "fOmegaHitsAdded", "fOmegacDecayRadiusStraTrack"}).Filter("!fTrueOmegac"); 
  
  auto pv_res = dfil_om_c_stra.Define("delta_pv", Distance, {"fPrimaryVertexX", "fPrimaryVertexY", "fPrimaryVertexZ"}); 
  auto h_pv_res = pv_res.Histo1D({"pv_resolution", "PV resolution", 150, 0, 3}, "delta_pv"); 
  
  auto dfv_res_stra = dfil_om_c_stra.Define("dlength_true", Distance, {"fOmegacVertexXMC", "fOmegacVertexYMC", "fOmegacVertexZMC"}).Define("delta_dlength", RelErr, {"dlength_true", "fOmegacDecayDistanceFromPVStraTrack"});  
  
  //begin loop over pT bins and Radius bins to define histograms 
  for (int ipt = 0; ipt < ptbins.size()-1; ++ipt) { 
    for (int ilayer = 0; ilayer < layerPos.size()-1; ++ilayer) { 
      
      TString binstring = TString::Format("_pT_%.1f_%.1f_decay_before_layer_%d", ptbins[ipt], ptbins[ipt+1], ilayer+1); 
      TString cutstring = TString::Format("(%.1f < fOmegaPtMC) && (fOmegaPtMC < %.1f) && (%.1f < fOmegaDecayRadiusMC) && (fOmegaDecayRadiusMC < %.1f)", ptbins[ipt], ptbins[ipt+1], layerPos[ilayer], layerPos[ilayer+1]); 
      
      auto dfv_res_stra_pt_layer = dfv_res_stra.Filter(cutstring.Data()); 
      resolutions_Stra.emplace_back(new ROOT::RDF::RResultPtr<::TH1D>(dfv_res_stra_pt_layer.Histo1D({"dv_res_stra"+binstring, "dv_res_stra"+binstring, 200, -1.5 , 1.5 }, "delta_dlength"))); 
      
      auto hit_perf = dfv_res_stra_pt_layer.Define("hit_perf", TString::Format("%d - fOmegaHitsAdded", ilayer).Data()); 
      layers_Stra.emplace_back(new ROOT::RDF::RResultPtr<::TH1D>(hit_perf.Histo1D({"hit_perf"+binstring, "hit_perf"+binstring, 11, -1 , 10 }, "hit_perf"))); 
      
      auto dfv_res_topo_pt_layer = dfv_res_topo.Filter(cutstring.Data()); 
      resolutions_Topo.emplace_back(new ROOT::RDF::RResultPtr<::TH1D>(dfv_res_topo_pt_layer.Histo1D({"dv_res_topo"+binstring, "dv_res_topo"+binstring, 200, -1.5 , 1.5 }, "delta_dlength"))); 
    }
  } 

  //Write histos 
  TFile *out = TFile::Open("outResolutions.root", "recreate"); 
  auto NormalizeToMaximum = [] (ROOT::RDF::RResultPtr<TH1D> hist) { hist->Scale(1./hist->GetMaximum());};
  auto NormalizeToEntries = [] (ROOT::RDF::RResultPtr<TH1D> hist) { hist->Scale(1./hist->GetEntries());};
  auto NormalizeToIntegral = [] (ROOT::RDF::RResultPtr<TH1D> hist) { if (hist->Integral() > 0) hist->Scale(1./hist->Integral());};
  auto Normalize2DToEntries = [] (ROOT::RDF::RResultPtr<TH2D> hist) { hist->Scale(1./hist->GetEntries());};
  auto DoNotNormalize = [] (ROOT::RDF::RResultPtr<TH1D> hist) { };
  for ( auto it : resolutions_Stra) { 
    NormalizeAndStore(out, *it, NormalizeToIntegral); 
  }
  TH2D* AddedHitPerformance = new TH2D("AddedHitPerformance", "Added Hit Performance", 12, -0.5, 11.5, 10, -0.5, 9.5); 
  AddedHitPerformance->GetYaxis()->SetTitle("Crossed Layers#minusAdded Hits"); 
  AddedHitPerformance->GetXaxis()->SetTitle("Crossed Layers"); 
  int ilayer =1; 
  for ( auto it : layers_Stra) { 
    NormalizeAndStore(out, *it, NormalizeToIntegral); 
    for (int ibin=2; ibin<=it->GetPtr()->GetNbinsX(); ++ibin) { 
      AddedHitPerformance->SetBinContent(ilayer, ibin-1, it->GetPtr()->GetBinContent(ibin)); 
    }
    ilayer++; 
  }
  out->cd(); 
  AddedHitPerformance->Write();
  h_pv_res->Write();
  for ( auto it : resolutions_Topo ) { 
    NormalizeAndStore(out, *it, NormalizeToIntegral); 
  }
  out->Close(); 
} 

