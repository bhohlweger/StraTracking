void comparator(TString addon) { 
  double xiccMass = 3.596; //3.621; 
  double xiccWindow = 0.08;
  
  TFile* xi = TFile::Open(    TString::Format("outxiccSelector_xi%s.root"  , addon.Data()), "read"); 
  TFile* xic = TFile::Open(   TString::Format("outxiccSelector_xic%s.root" , addon.Data()), "read"); 
  TFile* xicc = TFile::Open(  TString::Format("outxiccSelector_xicc%s.root", addon.Data()), "read"); 
  TFile* output = TFile::Open(TString::Format("outComp_cut%s.root"         , addon.Data()), "recreate"); 
  
  TH1D* h_xiCounter = (TH1D*)xi->Get("df_xi_c_candCounter");
  double normXi = 1/(9.25e-7*h_xiCounter->GetBinContent(1)); //4 for the production cross section, 20 for the BR
  
  TH1D* h_xicCounter = (TH1D*)xic->Get("df_xi_c_candCounter");
  double normXic = 1/(1.24e-4*h_xicCounter->GetBinContent(1)); //30 for production cross section, 20 for the BR
  
  TH1D* h_xiccCounter = (TH1D*)xicc->Get("df_xi_c_candCounter");
  double normXicc = (1.)/(h_xiccCounter->GetBinContent(1)); 
  
  TList* inList = xi->GetListOfKeys(); 
  TIter next(inList); 
  TObject* obj = nullptr; 
  TH1D* sumHist = nullptr;

  while ((obj = next())) {
    TString objName = TString::Format("%s",obj->GetName()); 
    if (objName.Contains("Counter")) { 
      continue; 
    }
    if (sumHist) { 
      sumHist = nullptr;
    }
    TH1D* xiHist = (TH1D*)xi->Get(obj->GetName()); 
    xiHist->GetXaxis()->SetTitle(obj->GetName()); 
    xiHist->SetLineColor(kPink+7); 
    xiHist->SetMarkerColor(kPink+7); 
    TH1D* xicHist = (TH1D*)xic->Get(obj->GetName()); 
    xicHist->SetLineColor(38); 
    xicHist->SetMarkerColor(38); 
    TH1D* xiccHist = (TH1D*)xicc->Get(obj->GetName()); 
    xiccHist->SetLineColor(kGreen+3); 
    xiccHist->SetMarkerColor(kGreen+3); 
    
    auto c = c11(obj->GetName()); 
    auto p = (TPad*)gROOT->FindObject(TString::Format("p%s",obj->GetName()).Data()); 
    p->cd(); 
    
    if (!(objName.Contains("mass") && objName.Contains("xi_cc"))) { 
      xiHist->SetTitle("#Xi + Pythia #pi"); 
      xicHist->SetTitle("#Xi_{c} + Pythia #pi"); 
      xiccHist->SetTitle("#Xi_{cc}"); 

      xiHist->Scale(1./xiHist->GetMaximum()); 
      xicHist->Scale(1./xicHist->GetMaximum()); 
      xiccHist->Scale(1./xiccHist->GetMaximum()); 
      
      xiHist->GetYaxis()->SetTitle("Count Normalized to Maximum"); 
    } else { 
      xiHist->SetTitle("#Xi + Pythia #pi"); 
      xiHist->Rebin(4); 
      xicHist->SetTitle("#Xi_{c} + Pythia #pi"); 
      xicHist->Rebin(4); 
      xiccHist->SetTitle("#Xi_{cc}"); 
      xiccHist->Rebin(4); 
      
      xiHist->Scale(normXi); 
      xicHist->Scale(normXic); 
      xiccHist->Scale(normXicc); 
      
      sumHist = (TH1D*)xiHist->Clone(TString::Format("sumHist_%s", xiHist->GetName())); 
      sumHist->Add(xicHist); 
      sumHist->Add(xiccHist); 
      sumHist->SetTitle("Sum"); 
      sumHist->GetYaxis()->SetTitle("Counts relative to #Xi_{cc}");       
      //get some number bby 
      double nxi = xiHist->Integral(xiHist->FindBin(xiccMass-xiccWindow),xiHist->FindBin(xiccMass+xiccWindow)); 
      double nxic = xicHist->Integral(xicHist->FindBin(xiccMass-xiccWindow),xicHist->FindBin(xiccMass+xiccWindow)); 
      double nxicc = xiccHist->Integral(xiccHist->FindBin(xiccMass-xiccWindow),xiccHist->FindBin(xiccMass+xiccWindow)); 
      std::cout <<"For " << obj->GetName() << ": \n nxi=" << nxi << " nxic="<< nxic << " nxicc=" << nxicc << " S/B=" << nxicc/(nxi+nxic) << std::endl;
    }

    if (sumHist) {
      sumHist->GetYaxis()->SetTitleOffset(1.2); 
      sumHist->GetXaxis()->SetNdivisions(506); 
      sumHist->GetYaxis()->SetNdivisions(504); 
      
      sumHist->SetLineColor(kOrange+7); 
      
      sumHist->Draw("hist"); 
      xiHist->GetYaxis()->SetRangeUser(0, sumHist->GetMaximum()*1.2); 
      xiHist->Draw("histsame");
      xicHist->Draw("histsame"); 
      xiccHist->Draw("histsame"); 

      TLine one; 
      one.SetLineColor(kBlack); 
      one.SetLineStyle(4);
      one.DrawLine(xiccMass-xiccWindow, 0, xiccMass-xiccWindow, sumHist->GetMaximum()); 
      one.DrawLine(xiccMass+xiccWindow, 0, xiccMass+xiccWindow, sumHist->GetMaximum()); 
      
    } else { 
      xiHist->GetYaxis()->SetTitleOffset(1.2); 
      xiHist->GetXaxis()->SetNdivisions(506); 
      xiHist->GetYaxis()->SetNdivisions(504); 
      xiHist->Draw("hist"); 
      xicHist->Draw("histsame"); 
      xiccHist->Draw("histsame"); 
    }
    auto leg = p->BuildLegend(0.59, 0.4, 0.9, 0.82, "Injected:", "l");
    leg->SetFillStyle(0); 
    leg->Draw("same"); 
    output->cd();     
    c->Write();
  }
  output->Close(); 
  xi->Close(); 
  xic->Close(); 
  xicc->Close(); 
}
