
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
  std::vector<ROOT::RDF::RResultPtr<::TH2D>*> layers_dca_xy;  
  std::vector<ROOT::RDF::RResultPtr<::TH2D>*> layers_dca_z;  

  std::vector<ROOT::RDF::RResultPtr<::TH2D>*> layers_dca_xy_refit;  
  std::vector<ROOT::RDF::RResultPtr<::TH2D>*> layers_dca_z_refit;  
  
  auto dfil_xi_stra = df_xi.Filter(HarryPlotter::TopoCuts_om, {"fXiPtMC"}).Filter("fTrueXi");
  
  //begin loop over pT bins and Radius bins to define histograms 
  for (int ilayer = 0; ilayer < layerPos.size()-1; ++ilayer) { 
    TString binstring = TString::Format("_decay_before_layer_%d", ilayer+1); 
    TString cutstring = TString::Format(" (%.1f < fXiDecayRadiusMC) && (fXiDecayRadiusMC < %.1f)", layerPos[ilayer], layerPos[ilayer+1]); 
    auto dfv_res_stra_pt_layer = dfil_xi_stra.Filter(cutstring.Data()); 

    layers_dca_xy.emplace_back(new ROOT::RDF::RResultPtr<::TH2D>(dfv_res_stra_pt_layer.Histo2D({"dca_xy"+binstring, "dca_xy"+binstring, 14, -1 , 13, 400, -200, 200 }, "lHitLayer", "lDCAxyToHit"))); 
    
    layers_dca_z.emplace_back(new ROOT::RDF::RResultPtr<::TH2D>(dfv_res_stra_pt_layer.Histo2D({"dca_z"+binstring, "dca_z"+binstring, 14, -1 , 13, 400, -200, 200 }, "lHitLayer", "lDCAzToHit"))); 

    layers_dca_xy_refit.emplace_back(new ROOT::RDF::RResultPtr<::TH2D>(dfv_res_stra_pt_layer.Histo2D({"dca_xy_refit"+binstring, "dca_xy_refit"+binstring, 14, -1 , 13, 400, -200, 200 }, "lHitLayer", "lDCAxyToHitInward"))); 
    
    layers_dca_z_refit.emplace_back(new ROOT::RDF::RResultPtr<::TH2D>(dfv_res_stra_pt_layer.Histo2D({"dca_z_refit"+binstring, "dca_z_refit"+binstring, 14, -1 , 13, 400, -200, 200 }, "lHitLayer", "lDCAzToHitInward"))); 
  }
  //Write histos 
  TFile *out = TFile::Open("hitstudiesxi.root", "recreate"); 
  
  TList* out_dca_xy = new TList(); 
  out_dca_xy->SetName("out_dca_xy"); 
  out_dca_xy->SetOwner(); 
  for (auto itHist : layers_dca_xy)  {
    out_dca_xy->Add(itHist->GetPtr()); 
    
    TList* outList = new TList(); 
    outList->SetOwner(); 
    outList->SetName(TString::Format("%s_proj", itHist->GetPtr()->GetName()).Data()); 
    out_dca_xy->Add(outList); 
    TString widthName = TString::Format("%s_width", itHist->GetPtr()->GetName()); 
    TH1D* width = new TH1D(widthName.Data(), widthName.Data(), 12, 0, 12);
    width->GetXaxis()->SetTitle("lHitLayer"); 
    width->GetYaxis()->SetTitle("#sigma_{DCA}"); 
    outList->Add(width); 
    for (int iLyr = 2; iLyr < 14; ++iLyr) { 
      auto proj = itHist->GetPtr()->ProjectionY(TString::Format("%s_lHitLayer_%d", itHist->GetPtr()->GetName(), iLyr-1).Data(), iLyr, iLyr); 
      if (proj->Integral() > 100) width->SetBinContent(iLyr-1, HarryPlotter::FitDCA(proj)); 
      outList->Add(proj); 
    }     
  }
  out->cd(); 
  out_dca_xy->Write(out_dca_xy->GetName(), 1); 
  
  TList* out_dca_z = new TList(); 
  out_dca_z->SetName("out_dca_z"); 
  out_dca_z->SetOwner(); 
  for (auto itHist : layers_dca_z)  {
    out_dca_z->Add(itHist->GetPtr()); 
    
    TList* outList = new TList(); 
    outList->SetOwner(); 
    outList->SetName(TString::Format("%s_proj", itHist->GetPtr()->GetName()).Data()); 
    out_dca_z->Add(outList); 
    TString widthName = TString::Format("%s_width", itHist->GetPtr()->GetName()); 
    TH1D* width = new TH1D(widthName.Data(), widthName.Data(), 12, 0, 12);
    width->GetXaxis()->SetTitle("lHitLayer"); 
    width->GetYaxis()->SetTitle("#sigma_{DCA}"); 
    outList->Add(width); 
    for (int iLyr = 2; iLyr < 14; ++iLyr) { 
      auto proj = itHist->GetPtr()->ProjectionY(TString::Format("%s_lHitLayer_%d", itHist->GetPtr()->GetName(), iLyr-1).Data(), iLyr, iLyr); 
      if (proj->Integral() > 100) width->SetBinContent(iLyr-1, HarryPlotter::FitDCA(proj)); 
      outList->Add(proj); 
    }     
  }
  out_dca_z->Write(out_dca_z->GetName(), 1); 

  TList* out_dca_xy_refit = new TList(); 
  out_dca_xy_refit->SetName("out_dca_xy_refit"); 
  out_dca_xy_refit->SetOwner(); 
  for (auto itHist : layers_dca_xy_refit)  {
    out_dca_xy_refit->Add(itHist->GetPtr()); 

    TList* outList = new TList(); 
    outList->SetOwner(); 
    outList->SetName(TString::Format("%s_proj", itHist->GetPtr()->GetName()).Data()); 
    out_dca_xy_refit->Add(outList); 
    TString widthName = TString::Format("%s_width", itHist->GetPtr()->GetName()); 
    TH1D* width = new TH1D(widthName.Data(), widthName.Data(), 12, 0, 12);
    width->GetXaxis()->SetTitle("lHitLayer"); 
    width->GetYaxis()->SetTitle("#sigma_{DCA}"); 
    outList->Add(width); 
    for (int iLyr = 2; iLyr < 14; ++iLyr) { 
      auto proj = itHist->GetPtr()->ProjectionY(TString::Format("%s_lHitLayer_%d", itHist->GetPtr()->GetName(), iLyr-1).Data(), iLyr, iLyr); 
      if (proj->Integral() > 100) width->SetBinContent(iLyr-1, HarryPlotter::FitDCA(proj)); 
      outList->Add(proj); 
    }     
  }
  out_dca_xy_refit->Write(out_dca_xy_refit->GetName(),1); 
  
  TList* out_dca_z_refit = new TList(); 
  out_dca_z_refit->SetName("out_dca_z_refit"); 
  out_dca_z_refit->SetOwner(); 
  for (auto itHist : layers_dca_z_refit)  {
    out_dca_z_refit->Add(itHist->GetPtr()); 
    
    TList* outList = new TList(); 
    outList->SetOwner(); 
    outList->SetName(TString::Format("%s_proj", itHist->GetPtr()->GetName()).Data()); 
    out_dca_z_refit->Add(outList); 
    TString widthName = TString::Format("%s_width", itHist->GetPtr()->GetName()); 
    TH1D* width = new TH1D(widthName.Data(), widthName.Data(), 12, 0, 12);
    width->GetXaxis()->SetTitle("lHitLayer"); 
    width->GetYaxis()->SetTitle("#sigma_{DCA}"); 
    outList->Add(width); 
    for (int iLyr = 2; iLyr < 14; ++iLyr) { 
      auto proj = itHist->GetPtr()->ProjectionY(TString::Format("%s_lHitLayer_%d", itHist->GetPtr()->GetName(), iLyr-1).Data(), iLyr, iLyr); 
      if (proj->Integral() > 100) width->SetBinContent(iLyr-1, HarryPlotter::FitDCA(proj)); 
      outList->Add(proj); 
    }     
  }
  out_dca_z_refit->Write(out_dca_z_refit->GetName(),1); 
  
  out->Close(); 
} 

