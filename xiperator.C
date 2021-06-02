void xiperator(TString addon) { 
  double xiMass = 1.322; //3.621; 
  double xiWindow = 0.008;

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
    TH1D* mbHist = (TH1D*)mb->Get(obj->GetName()); 
    mbHist->SetLineColor(42); 
    mbHist->SetMarkerColor(42); 
    
    auto c = c11(obj->GetName()); 
    auto p = (TPad*)gROOT->FindObject(TString::Format("p%s",obj->GetName()).Data()); 
    p->cd(); 
    
    xiHist->SetTitle("#Xi injected in Pythia"); 
    mbHist->SetTitle("Pythia MB"); 

    if (!(objName.Contains("mass"))) { 
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
      sumHist->SetTitle("Sum"); 
      sumHist->GetYaxis()->SetTitle("Counts relative to #Xi");       
      //get some number bby 
      double nxi = xiHist->Integral(xiHist->FindBin(xiMass-xiWindow),xiHist->FindBin(xiMass+xiWindow)); 
      double nmb = mbHist->Integral(mbHist->FindBin(xiMass-xiWindow),mbHist->FindBin(xiMass+xiWindow)); 
      std::cout <<"For " << obj->GetName() << ": \n nxi=" << nxi << " nmb="<< nmb << " S/B=" << nxi/(nmb) << std::endl;
    }

    if (sumHist) {
      sumHist->GetYaxis()->SetTitleOffset(1.2); 
      sumHist->GetXaxis()->SetNdivisions(506); 
      sumHist->GetYaxis()->SetNdivisions(504); 
      
      sumHist->SetLineColor(kOrange+7); 
      
      sumHist->Draw("hist"); 
      xiHist->GetYaxis()->SetRangeUser(0, sumHist->GetMaximum()*1.2); 
      xiHist->Draw("histsame");
      mbHist->Draw("histsame"); 

      TLine one; 
      one.SetLineColor(kBlack); 
      one.SetLineStyle(4);
      one.DrawLine(xiMass-xiWindow, 0, xiMass-xiWindow, sumHist->GetMaximum()); 
      one.DrawLine(xiMass+xiWindow, 0, xiMass+xiWindow, sumHist->GetMaximum()); 
      
    } else { 
      xiHist->GetYaxis()->SetTitleOffset(1.2); 
      xiHist->GetXaxis()->SetNdivisions(506); 
      xiHist->GetYaxis()->SetNdivisions(504); 
      xiHist->Draw("hist"); 
      mbHist->Draw("histsame"); 
    }
    auto leg = p->BuildLegend(0.59, 0.4, 0.9, 0.82, "Injected:", "l");
    leg->SetFillStyle(0); 
    leg->Draw("same"); 
    output->cd();     
    c->Write();
  }
  output->Close(); 
  xi->Close(); 
  mb->Close(); 
}
