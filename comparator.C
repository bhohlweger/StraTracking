void comparator(TString addon) { 
  double xiccMass = 3.596; //3.621; 
  double xiccWindow = 0.08;
  
  TFile* xi = TFile::Open(    TString::Format("outxiccSelector_xi%s.root"  , addon.Data()), "read"); 
  TFile* xic = TFile::Open(   TString::Format("outxiccSelector_xic%s.root" , addon.Data()), "read"); 
  TFile* xicc = TFile::Open(  TString::Format("outxiccSelector_xicc%s.root", addon.Data()), "read"); 
  TFile* output = TFile::Open(TString::Format("outComp_cut%s.root"         , addon.Data()), "recreate"); 
  TList* inList = xi->GetListOfKeys(); 
  TIter next(inList); 
  TObject* obj = nullptr; 
  
  while ((obj = next())) {
    TString objName = TString::Format("%s",obj->GetName()); 
    if (objName.Contains("Counter")  || objName.Contains("_vs_") ) { 
      continue; 
    }

    TH1D* xiHist = (TH1D*)xi->Get(obj->GetName()); 
    TH1D* xicHist = (TH1D*)xic->Get(obj->GetName()); 
    TH1D* xiccHist = (TH1D*)xicc->Get(obj->GetName()); 
    
    xiHist->SetLineColor(kPink+7); 
    xiHist->SetMarkerColor(kPink+7); 
    xicHist->SetLineColor(38); 
    xicHist->SetMarkerColor(38); 
    xiccHist->SetLineColor(kGreen+3); 
    xiccHist->SetMarkerColor(kGreen+3); 
    
    if (!(objName.Contains("mass") && objName.Contains("xi_cc"))) { 
      xiHist->GetYaxis()->SetTitleOffset(1.2); 
      xiHist->GetXaxis()->SetNdivisions(506); 
      xiHist->GetYaxis()->SetNdivisions(504); 
      /*
      xiHist->Rebin(4); 
      xicHist->Rebin(4); 
      xiccHist->Rebin(4); 
      */
      xiHist->SetTitle("#Xi^{#minus} + 3#times#pi_{Pythia}"); 
      xicHist->SetTitle("#Xi^{+}_{c} + #pi_{Pythia}"); 
      xiccHist->SetTitle("#Xi_{cc}^{++}"); 

      xiHist->Scale(1./xiHist->GetMaximum()); 
      xicHist->Scale(1./xicHist->GetMaximum()); 
      xiccHist->Scale(1./xiccHist->GetMaximum()); 
    
      xiHist->GetXaxis()->SetTitle(obj->GetName()); 
      xiHist->GetYaxis()->SetRangeUser(0.,1.7); 
      xiHist->GetYaxis()->SetTitle("Count Normalized to Maximum"); 
    } else { 
      continue; 
    }
    auto c = c11(obj->GetName()); 
    auto p = (TPad*)gROOT->FindObject(TString::Format("p%s",obj->GetName()).Data()); 
    p->cd(); 
    
    auto leg = new TLegend(0.16, 0.45, 0.56, 0.8, "#splitline{ALICE 3 (Layout v1) Full Simluation}{Pythia pp #sqrt{s} = 13 TeV + GEANT3}");
    leg->SetFillStyle(0); 

    p->cd(); 
    leg->AddEntry(xiHist, xiHist->GetTitle(), "l");
    leg->AddEntry(xicHist, xicHist->GetTitle(), "l");
    leg->AddEntry(xiccHist, xiccHist->GetTitle(), "l");
    xiHist->Draw("hist"); 
    xicHist->Draw("histsame"); 
    xiccHist->Draw("histsame"); 
    leg->Draw("same"); 
    output->cd();     
    c->Write();
    c->Close(); 
  }
  output->Close(); 
  xi->Close(); 
  xic->Close(); 
  xicc->Close(); 
}
