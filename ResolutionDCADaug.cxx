
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

void DCADaugCuts(bool signal, const char* howmanycs) { 
  HarryPlotter::StyleBox(); 

  ROOT::EnableImplicitMT(); // Tell ROOT you want to go parallel          

  TString filePath; 

  if (signal) { 
    filePath = TString::Format("/localstore/alice/hohlweger/analysis/StrangnessTracking/210323/gimmezetree.root"); 
  } else { 
    filePath = TString::Format("%s",HarryPlotter::FilePathPythia()); 
  }
  TFile *file_c = new TFile(filePath, "READ");

  ROOT::RDataFrame df_ca_c("fTreeTripleC", file_c);

  //add classic quality cuts ... 
  //e.g. fOmegaCPtMC>1.0&&fOmegaCPtMC<3.0 
  //Strangeness tracking: fOmegaHitsAdded >= 3
  std::vector<TString> out_names = {"_integrated"};  

  TString selection_signal = 
    TString::Format("fCorrectPionFromOmegaC&&fCorrectPionFromOmegaCC&&fCorrectPionFromOmegaCCC"); 
  
  TString selection_background = 
    TString::Format("fFirstCombinationC&&fFirstCombinationCC&&fFirstCombinationCCC&&!fCorrectPionFromOmega%s", howmanycs); 
  
  TString variableTopo = TString::Format("fOmega%sDecayDCATopo", howmanycs); 
  TString variableStra = TString::Format("fOmega%sDecayDCAStraTrack", howmanycs); 

  std::cout << "selection_signal: " << selection_signal.Data() << std::endl;
  std::cout << "selection_background: " << selection_background.Data() << std::endl;  
  std::cout << "variableTopo: " << variableTopo.Data() << std::endl; 
  std::cout << "variableStra: " << variableStra.Data() << std::endl; 
  
  ROOT::RDF::RResultPtr<TH1D> dummy;  
  std::vector<ROOT::RDF::RResultPtr<::TH1D>*> ca_c_dca_daug_topo;  
  std::vector<ROOT::RDF::RResultPtr<::TH1D>*> om_c_dca_daug_topo;  
  
  std::vector<ROOT::RDF::RResultPtr<::TH1D>*> ca_c_dca_daug_stra;  
  std::vector<ROOT::RDF::RResultPtr<::TH1D>*> om_c_dca_daug_stra;  
  
  std::vector<ROOT::RDF::RResultPtr<::TH2D>*> ca_c_dca_daug_stra_addedHits;  
  std::vector<ROOT::RDF::RResultPtr<::TH2D>*> om_c_dca_daug_stra_addedHits;  
  
  //define filter 
  auto dfil_om_c_topo = df_ca_c.
    Filter(HarryPlotter::TopoCuts_om, {"fOmegaCPtMC"}).
    Filter(selection_signal.Data());
  
  auto dfil_ca_c_topo = df_ca_c.
    Filter(HarryPlotter::TopoCuts_om, {"fOmegaCPtMC"}).
    Filter(selection_background.Data()); 
  
  auto dfil_om_c_stra = df_ca_c.
    Filter(HarryPlotter::StraCuts_om, {"fOmegaCPtMC", "fOmegaHitsAdded"}).
    Filter(selection_signal.Data()); 
  
  auto dfil_ca_c_stra = df_ca_c.
    Filter(HarryPlotter::StraCuts_om, {"fOmegaCPtMC", "fOmegaHitsAdded"}).
    Filter(selection_background.Data()); 
  
  auto dfil_om_c_stra_addedHits = df_ca_c.//use the topo cuts here to keep also the ones with 1 added hit etc. 
    Filter(HarryPlotter::TopoCuts_om, {"fOmegaCPtMC"}).
    Filter(selection_signal.Data()); 

  auto dfil_ca_c_stra_addedHits = df_ca_c.
    Filter(HarryPlotter::TopoCuts_om, {"fOmegaCPtMC"}).
    Filter(selection_background.Data()); 
  ca_c_dca_daug_topo.emplace_back(new ROOT::RDF::RResultPtr<::TH1D>(dfil_ca_c_topo.Histo1D({"ca_c_dca_daug_topo","#Omega_{c}#rightarrow#Omega+#pi", 300 , 0, 200}, variableTopo.Data()))); 
  om_c_dca_daug_topo.emplace_back(new ROOT::RDF::RResultPtr<::TH1D>(dfil_om_c_topo.Histo1D({"om_c_dca_daug_topo","#Omega_{c}#rightarrow#Omega+#pi",300, 0, 200}, variableTopo.Data()))); 
  
  ca_c_dca_daug_stra.emplace_back(new ROOT::RDF::RResultPtr<::TH1D>(dfil_ca_c_stra.Histo1D({"ca_c_dca_daug_stra","#Omega_{c}#rightarrow#Omega+#pi",300, 0, 200}, variableStra.Data()))); 
  om_c_dca_daug_stra.emplace_back(new ROOT::RDF::RResultPtr<::TH1D>(dfil_om_c_stra.Histo1D({"om_c_dca_daug_stra","#Omega_{c}#rightarrow#Omega+#pi",300, 0, 200}, variableStra.Data()))); 
  
  ca_c_dca_daug_stra_addedHits.emplace_back(new ROOT::RDF::RResultPtr<::TH2D>(dfil_ca_c_stra_addedHits.Histo2D({"ca_c_dca_daug_stra_addedHits","#Omega_{c}#rightarrow#Omega+#pi", 10, 0, 10,300, 0, 200}, "fOmegaHitsAdded", variableStra.Data()))); 
  om_c_dca_daug_stra_addedHits.emplace_back(new ROOT::RDF::RResultPtr<::TH2D>(dfil_om_c_stra_addedHits.Histo2D({"om_c_dca_daug_stra_addedHits","#Omega_{c}#rightarrow#Omega+#pi", 10, 0, 10,300, 0, 200}, "fOmegaHitsAdded", variableStra.Data()))); 
    
  for (int ipt = 0; ipt < ptbins.size()-1; ++ipt) { 
    TString binstring = TString::Format("_pT_%.1f_%.1f", ptbins[ipt], ptbins[ipt+1]); 
    TString cutstring = TString::Format("(%.1f < fOmegaCPtMC) && (fOmegaCPtMC < %.1f)", ptbins[ipt], ptbins[ipt+1]); 
    out_names.push_back(binstring); 
    
    auto dfil_om_c_topo_pT = dfil_om_c_topo.Filter(cutstring.Data()); 
    auto dfil_ca_c_topo_pT = dfil_ca_c_topo.Filter(cutstring.Data()); 
    
    ca_c_dca_daug_topo.emplace_back(new ROOT::RDF::RResultPtr<::TH1D>(dfil_ca_c_topo_pT.Histo1D({"ca_c_dca_daug_topo"+binstring,"#Omega_{c}#rightarrow#Omega+#pi", 300, 0, 200}, variableTopo.Data()))); 
    om_c_dca_daug_topo.emplace_back(new ROOT::RDF::RResultPtr<::TH1D>(dfil_om_c_topo_pT.Histo1D({"om_c_dca_daug_topo"+binstring,"#Omega_{c}#rightarrow#Omega+#pi", 300, 0, 200}, variableTopo.Data()))); 
    
    auto dfil_om_c_stra_pT = dfil_om_c_stra.Filter(cutstring.Data()); 
    auto dfil_ca_c_stra_pT = dfil_ca_c_stra.Filter(cutstring.Data()); 
    
    ca_c_dca_daug_stra.emplace_back(new ROOT::RDF::RResultPtr<::TH1D>(dfil_ca_c_stra_pT.Histo1D({"ca_c_dca_daug_stra"+binstring,"#Omega_{c}#rightarrow#Omega+#pi", 300, 0, 200}, variableStra.Data()))); 
    om_c_dca_daug_stra.emplace_back(new ROOT::RDF::RResultPtr<::TH1D>(dfil_om_c_stra_pT.Histo1D({"om_c_dca_daug_stra"+binstring,"#Omega_{c}#rightarrow#Omega+#pi", 300, 0, 200}, variableStra.Data()))); 
    
    auto dfil_om_c_stra_pT_addedHits = dfil_om_c_stra_addedHits.Filter(cutstring.Data()); 
    auto dfil_ca_c_stra_pT_addedHits = dfil_ca_c_stra_addedHits.Filter(cutstring.Data()); 
    
    ca_c_dca_daug_stra_addedHits.emplace_back(new ROOT::RDF::RResultPtr<::TH2D>(dfil_ca_c_stra_pT_addedHits.Histo2D({"ca_c_dca_daug_stra_addedHits"+binstring,"#Omega_{c}#rightarrow#Omega+#pi", 10, 0, 10, 300, 0, 200}, "fOmegaHitsAdded",variableStra.Data()))); 
    om_c_dca_daug_stra_addedHits.emplace_back(new ROOT::RDF::RResultPtr<::TH2D>(dfil_om_c_stra_pT_addedHits.Histo2D({"om_c_dca_daug_stra_addedHits"+binstring,"#Omega_{c}#rightarrow#Omega+#pi", 10, 0, 10, 300, 0, 200}, "fOmegaHitsAdded", variableStra.Data()))); 
    
  }
  //Write histos 
  TFile *out = TFile::Open(TString::Format("outResolutionsDCADaug_%s_%s.root", howmanycs, signal?"signal":"bkg"), "recreate"); 			   


  for (int ipt = 0; ipt < ptbins.size(); ++ipt) { 
    
    ca_c_dca_daug_topo[ipt]->GetPtr()->Scale(1./ca_c_dca_daug_topo[ipt]->GetPtr()->Integral()); 
    if (signal) om_c_dca_daug_topo[ipt]->GetPtr()->Scale(1./om_c_dca_daug_topo[ipt]->GetPtr()->Integral()); 
   
    HarryPlotter::CheckAndStore(out, ca_c_dca_daug_topo[ipt]->GetPtr()); 
    if (signal) HarryPlotter::CheckAndStore(out, om_c_dca_daug_topo[ipt]->GetPtr()); 
    
    ca_c_dca_daug_stra[ipt]->GetPtr()->Scale(1./ca_c_dca_daug_stra[ipt]->GetPtr()->Integral()); 
    if (signal) om_c_dca_daug_stra[ipt]->GetPtr()->Scale(1./om_c_dca_daug_stra[ipt]->GetPtr()->Integral()); 
    
    HarryPlotter::CheckAndStore(out, ca_c_dca_daug_stra[ipt]->GetPtr()); 
    if (signal) HarryPlotter::CheckAndStore(out, om_c_dca_daug_stra[ipt]->GetPtr()); 

    TH1F* width_dca_daug_om_c = new TH1F("width_dca_daug_om_c"+out_names[ipt], "Width DCA_{om_c}", 11, -1, 10); 
    width_dca_daug_om_c->GetYaxis()->SetTitle("#sigma_{DCA} (#mum)"); 
    TH1F* width_dca_daug_ca_c = new TH1F("width_dca_daug_ca_c" + out_names[ipt], "Width DCA_{ca_c}", 11, -1, 10);  
    width_dca_daug_ca_c->GetYaxis()->SetTitle("#sigma_{DCA} (#mum)"); 
    
    if (signal) HarryPlotter::Normalize2DBinByBin(om_c_dca_daug_stra_addedHits[ipt]->GetPtr()); 
    HarryPlotter::Normalize2DBinByBin(ca_c_dca_daug_stra_addedHits[ipt]->GetPtr()); 

    if (signal) HarryPlotter::CheckAndStore(out, om_c_dca_daug_stra_addedHits[ipt]->GetPtr()); 
    HarryPlotter::CheckAndStore(out, ca_c_dca_daug_stra_addedHits[ipt]->GetPtr()); 
    
    for (int iBinX = 1; iBinX <= om_c_dca_daug_stra_addedHits[ipt]->GetPtr()->GetNbinsX(); ++iBinX) { 
      om_c_dca_daug_stra_addedHits[ipt]->GetPtr()->GetXaxis()->SetRange(iBinX, iBinX); 
      width_dca_daug_om_c->SetBinContent(iBinX+1, om_c_dca_daug_stra_addedHits[ipt]->GetPtr()->GetStdDev(2)); 
      width_dca_daug_om_c->GetXaxis()->SetBinLabel(iBinX+1, TString::Format("%d added Hits", iBinX).Data()); 
    }      
    width_dca_daug_om_c->SetBinContent(1, om_c_dca_daug_topo[ipt]->GetPtr()->GetStdDev(1)); 
    width_dca_daug_om_c->GetXaxis()->SetBinLabel(1, "Topological Tracking"); 
    
    for (int iBinX = 1; iBinX <= ca_c_dca_daug_stra_addedHits[ipt]->GetPtr()->GetNbinsX(); ++iBinX) { 
      ca_c_dca_daug_stra_addedHits[ipt]->GetPtr()->GetXaxis()->SetRange(iBinX, iBinX); 
      width_dca_daug_ca_c->SetBinContent(iBinX+1, ca_c_dca_daug_stra_addedHits[ipt]->GetPtr()->GetStdDev(2)); 
      width_dca_daug_ca_c->GetXaxis()->SetBinLabel(iBinX+1, TString::Format("%d added Hits", iBinX).Data()); 
    }
    width_dca_daug_ca_c->SetBinContent(1, ca_c_dca_daug_topo[ipt]->GetPtr()->GetStdDev(1)); 
    width_dca_daug_ca_c->GetXaxis()->SetBinLabel(1, "Topological Tracking"); 
    
    out->cd(); 
    if (signal) width_dca_daug_om_c->Write(); 
    width_dca_daug_ca_c->Write(); 
  } 
  out->Close(); 
} 


int main(int argc, char **argv) {
  DCADaugCuts(true, "C"); 
  DCADaugCuts(true, "CC"); 
  DCADaugCuts(true, "CCC"); 

  DCADaugCuts(false, "C"); 
  DCADaugCuts(false, "CC"); 
  DCADaugCuts(false, "CCC"); 
}
