void xiperator(TString addon) { 
  double xiMass = 1.322; //3.621; 
  double xiWindow = 0.012;

  TFile* mb = TFile::Open(    TString::Format("outxiccSelector_mb%s.root"  , addon.Data()), "read"); 
  TFile* xi = TFile::Open(    TString::Format("outxiccSelector_xi%s.root"  , addon.Data()), "read"); 
  TFile* output = TFile::Open(TString::Format("outComp_xiperator%s.root"   , addon.Data()), "recreate"); 
  
  TH1D* h_mbCounter = (TH1D*)mb->Get("df_xi_c_candCounter");
  double normMb = 1/(9.91e-3*h_mbCounter->GetBinContent(1)); //4 for the production cross section, 20 for the BR
    
  TH1D* h_xiCounter = (TH1D*)xi->Get("df_xi_c_candCounter");
  double normXi = 1/(h_xiCounter->GetBinContent(1)); //4 for the production cross section, 20 for the BR
  
  TList* inList = xi->GetListOfKeys(); 
  TIter next(inList); 
  TObject* obj = nullptr; 
  TH1D* sumHist = nullptr;

  double nxi; 
  double nmb; 
  
  while ((obj = next())) {
    TString objName = TString::Format("%s",obj->GetName()); 
    if (objName.Contains("Counter") || objName.Contains("_vs_")) { 
      continue; 
    }
    if (sumHist) { 
      sumHist = nullptr;
    }

    TH1D* xiHist = (TH1D*)xi->Get(obj->GetName()); 
    xiHist->GetXaxis()->SetTitle(obj->GetName()); 
    xiHist->SetLineColor(kPink+7); 
    xiHist->SetMarkerColor(kPink+7); 
    TH1D* mbHist = (TH1D*)mb->Get(obj->GetName()); 
    mbHist->SetLineColor(42); 
    mbHist->SetMarkerColor(42); 
    
    xiHist->SetTitle("Pure #Xi^{#minus}"); 
    mbHist->SetTitle("Pythia MB w/o #Xi^{#minus}"); 

    if (!(objName.Contains("mass"))) { 
      xiHist->GetYaxis()->SetTitleOffset(1.2); 
      xiHist->GetXaxis()->SetNdivisions(506); 
      xiHist->GetYaxis()->SetNdivisions(504); 

      xiHist->Scale(1./xiHist->GetMaximum()); 
      mbHist->Scale(1./mbHist->GetMaximum()); 
      
      xiHist->GetYaxis()->SetTitle("Count Normalized to Maximum"); 
    } else { 
      if (objName.Contains("xi_cc") || objName.Contains("xi_c")) { 
	continue;
      }

      xiHist->Scale(normXi); 
      mbHist->Scale(normMb); 
      
      sumHist = (TH1D*)xiHist->Clone(TString::Format("sumHist_%s", xiHist->GetName())); 
      sumHist->Add(mbHist); 
      sumHist->SetTitle("Sum scaled #Xi^{#minus} and Pythia MB"); 
      sumHist->GetYaxis()->SetTitle("Counts relative to #Xi");       
      sumHist->GetYaxis()->SetRangeUser(0, sumHist->GetMaximum()*3); 
      sumHist->GetXaxis()->SetRangeUser(1.15, 1.5); 
      sumHist->GetXaxis()->SetTitle("IM(#Lambda,#pi^{-}) (GeV/#it{c}^{2})");
            
      sumHist->GetYaxis()->SetTitleOffset(1.4); 
      sumHist->GetXaxis()->SetNdivisions(506); 
      sumHist->GetYaxis()->SetNdivisions(504); 
      
      xiHist->SetFillColor(kPink+7); 
      xiHist->SetFillStyle(1001);
      
      sumHist->SetLineColor(kOrange+7); 
      sumHist->SetLineStyle(2);

      //get some number bby 
      nxi = xiHist->Integral(xiHist->FindBin(xiMass-xiWindow),xiHist->FindBin(xiMass+xiWindow)); 
      nmb = mbHist->Integral(mbHist->FindBin(xiMass-xiWindow),mbHist->FindBin(xiMass+xiWindow)); 
      std::cout <<"For " << obj->GetName() << ": \n nxi=" << nxi << " nmb="<< nmb << " S/B=" << nxi/(nmb) << std::endl;
    }
    auto c = c11(obj->GetName()); 
    auto p = (TPad*)gROOT->FindObject(TString::Format("p%s",obj->GetName()).Data()); 
    p->cd(); 
    
    auto leg = new TLegend(0.15, 0.51, 0.56, 0.8, "#splitline{ALICE 3 Full Simluation}{Pythia pp #sqrt{s} = 13 TeV + GEANT3}");
    leg->SetFillStyle(0); 
        
    if (sumHist) {
      sumHist->Draw("hist"); 
      xiHist->Draw("histsame");
      // mbHist->Draw("histsame"); 
      leg->AddEntry(xiHist, xiHist->GetTitle(), "l"); 
      leg->AddEntry(sumHist, sumHist->GetTitle(), "l"); 
      TLatex* myTex = GenTex(); 
      myTex->DrawLatex(0.17, 0.45, TString::Format("#Xi^{#minus} S/B = %.1f",nxi/nmb).Data());
    } else { 
      xiHist->Draw("hist"); 
      mbHist->Draw("histsame"); 
      leg->AddEntry(xiHist, xiHist->GetTitle(), "l"); 
      leg->AddEntry(mbHist, mbHist->GetTitle(), "l"); 
    }
    leg->Draw("same"); 
    output->cd();     
    c->Write();
    c->Close();
  }
  output->Close(); 
  xi->Close(); 
  mb->Close(); 
}
