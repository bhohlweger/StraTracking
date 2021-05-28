
template<typename T> 
void HarryPlotter::CheckAndStore(TFile *out, T key) { 
  TObject* obj = out->FindKey(key->GetName());
  if (!obj) { 
    out->cd(); 
    key->Write(); 
  }
}

template<typename T> 
void HarryPlotter::NormToMaxCheckAndStore(TFile *out, T key) { 
  key->Scale(1./key->GetMaximum()); 
  HarryPlotter::CheckAndStore(out, key);
}

template <typename T> void HarryPlotter::CumulateAndStore(TFile *out, T one) { 
  TH1D* cu_one = (TH1D*)one->GetCumulative(); 
  HarryPlotter::CheckAndStore(out, cu_one); 
};

template <typename T, typename F> 
void HarryPlotter::PlotAndStore(TString name, TFile* out, T one, T two, T three, T four) { 
  auto c = new TCanvas(name, ""); 
  c->cd(); 
  one->SetLineColor(kPink+7); 
  one->SetMarkerColor(kPink+7);											  
  one->SetTitle(TString::Format("%s;%s;", one->GetTitle(), name.Data()).Data()); 
  one->Draw(); 
  HarryPlotter::CheckAndStore(out, one); 

  if (two) {
    two->SetLineColor(kAzure-4);
    two->SetMarkerColor(kAzure-4);   
    two->SetLineStyle(2);
    two->Draw("same"); 
    HarryPlotter::CheckAndStore(out, two); 
  }
  if (three) { 
    three->SetLineColor(kOrange+7); 
    three->SetMarkerColor(kOrange+7); 
    three->Draw("same"); 
    HarryPlotter::CheckAndStore(out, three); 
  }

  if (four) { 
    four->SetLineColor(kBlue); 
    four->SetMarkerColor(kBlue); 
    four->SetLineStyle(2);
    four->Draw("same"); 
    HarryPlotter::CheckAndStore(out, four); 
  }
  c->BuildLegend(); 
  out->cd(); 
  c->Write(); 
}

template <typename T,typename F> void HarryPlotter::CumulatePlotAndStore(TString name, TFile* out, T one, T two, T three, T four) { 
  auto c = new TCanvas(name, ""); 
  c->cd(); 
  TH1D* cu_one = (TH1D*)one->GetCumulative(); 
  cu_one->SetLineColor(kPink+7); 
  cu_one->SetMarkerColor(kPink+7); 
  cu_one->SetTitle(TString::Format("%s;%s;", cu_one->GetTitle(), name.Data()).Data()); 
  cu_one->Draw(); 
  HarryPlotter::CheckAndStore(out, cu_one); 

  if (two) {
    TH1D* cu_two = (TH1D*)two->GetCumulative(); 
    cu_two->SetName(TString::Format("%s_cumulative", two->GetName()).Data()); 
    cu_two->SetLineColor(kAzure-4);
    cu_two->SetMarkerColor(kAzure-4);   
    cu_two->SetLineStyle(2);
    cu_two->SetTitle(two->GetTitle()); 
    cu_two->Draw("same"); 
    HarryPlotter::CheckAndStore(out, cu_two); 
  }
  if (three) { 
    TH1D* cu_three = (TH1D*)three->GetCumulative(); 
    cu_three->SetLineColor(kOrange+7); 
    cu_three->SetMarkerColor(kOrange+7); 
    cu_three->SetTitle(three->GetTitle()); 
    cu_three->Draw("same"); 
    HarryPlotter::CheckAndStore(out, cu_three); 
  }

  if (four) { 
    TH1D* cu_four = (TH1D*)four->GetCumulative(); 
    cu_four->SetLineColor(kBlue); 
    cu_four->SetMarkerColor(kBlue); 
    cu_four->SetLineStyle(2);
    cu_four->SetTitle(four->GetTitle());
    cu_four->Draw("same"); 
    HarryPlotter::CheckAndStore(out, cu_four); 
  }
  c->BuildLegend();
  out->cd(); 
  c->Write(); 
}

template <typename T> void HarryPlotter::AverageBackground(TFile* out, T signal, T background) { 
  TString signalName = signal->GetName(); 
  TString backgroundName = background->GetName(); 
    
  TH1D* avg_background = (TH1D*)background->Clone(backgroundName+"avg"); 
  avg_background->Reset(); 
  double d_avg_background = background->Integral()/background->GetNbinsX(); 
  for (int ix = 1; ix < avg_background->GetNbinsX()+1; ++ix) { 
    avg_background->SetBinContent(ix, d_avg_background); 
    //avg_background->SetBinError(ix,1e-8);
  }
  TH1D* signal_with_bkg = (TH1D*)signal->Clone(signalName+"_bkg"); 
  signal_with_bkg->Add(avg_background); 
  TH1D* signal_with_bkg_ratio = (TH1D*)signal_with_bkg->Clone(signalName+"_bkg_ratio"); 
  signal_with_bkg_ratio->Divide(avg_background); 
  
  out->cd(); 
  avg_background->Write(); 
  signal_with_bkg->Write(); 
  signal_with_bkg_ratio->Write(); 
  
}
