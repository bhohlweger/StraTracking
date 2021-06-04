void comparator(TString addon) { 
  double xiccMass = 3.596; //3.621; 
  double xiccWindow = 0.08;
  
  TFile* xi = TFile::Open(    TString::Format("outxiccSelector_xi%s.root"  , addon.Data()), "read"); 
  TFile* xic = TFile::Open(   TString::Format("outxiccSelector_xic%s.root" , addon.Data()), "read"); 
  TFile* xicc = TFile::Open(  TString::Format("outxiccSelector_xicc%s.root", addon.Data()), "read"); 
  TFile* output = TFile::Open(TString::Format("outComp_cut%s.root"         , addon.Data()), "recreate"); 
  
  TH1D* h_xiRedCounter = (TH1D*)xi->Get("cutCounter");
  TH1D* h_xiCounter = (TH1D*)xi->Get("df_xi_c_candCounter");
  double normXi = 1/(9.25e-7*h_xiCounter->GetBinContent(1)); //4 for the production cross section, 20 for the BR

  TH1D* h_xicRedCounter = (TH1D*)xic->Get("cutCounter");  
  TH1D* h_xicCounter = (TH1D*)xic->Get("df_xi_c_candCounter");
  double normXic = 1/(1.24e-4*h_xicCounter->GetBinContent(1)); //30 for production cross section, 20 for the BR
  
  TH1D* h_xiccRedCounter = (TH1D*)xicc->Get("cutCounter");
  TH1D* h_xiccCounter = (TH1D*)xicc->Get("df_xi_c_candCounter");
  double normXicc = (1.)/(h_xiccCounter->GetBinContent(1)); 
  
  TList* inList = xi->GetListOfKeys(); 
  TIter next(inList); 
  TObject* obj = nullptr; 
  TH1D* sumHist = nullptr;
  TH1D* sumHistBkg = nullptr;
  int counter = 0; 
  double redXi; 
  double redXic;
  double redXicc; 

  while ((obj = next())) {
    TString objName = TString::Format("%s",obj->GetName()); 
    if (objName.Contains("Counter")  || objName.Contains("_vs_") ) { 
      continue; 
    }
    
    if (sumHist || sumHistBkg) { 
      sumHist = nullptr;
      sumHistBkg = nullptr;
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
      continue; 
      
      xiHist->GetYaxis()->SetTitleOffset(1.2); 
      xiHist->GetXaxis()->SetNdivisions(506); 
      xiHist->GetYaxis()->SetNdivisions(504); 
      
      xiHist->SetTitle("#Xi^{#minus} + #pi_{Pythia} + #pi_{Pythia} + #pi_{Pythia}"); 
      xicHist->SetTitle("#Xi^{+}_{c} + #pi_{Pythia}"); 
      xiccHist->SetTitle("#Xi_{cc}^{++}"); 

      xiHist->Scale(1./xiHist->GetMaximum()); 
      xicHist->Scale(1./xicHist->GetMaximum()); 
      xiccHist->Scale(1./xiccHist->GetMaximum()); 
    
      xiHist->GetXaxis()->SetTitle(obj->GetName()); 
      xiHist->GetYaxis()->SetRangeUser(0.,1.4); 
      xiHist->GetYaxis()->SetTitle("Count Normalized to Maximum"); 
    } else { 
      
      xiHist->Rebin(4);       
      xicHist->Rebin(4); 
      xiccHist->Rebin(4); 
      
      xiHist->Scale(normXi); 
      xicHist->Scale(normXic); 
      xiccHist->Scale(normXicc); 
      
      redXi = h_xiRedCounter->GetBinContent(counter); 
      redXi = TMath::Abs(redXi) > 1e-30?1./redXi:-99; 
      redXic = h_xicRedCounter->GetBinContent(counter);
      redXic = TMath::Abs(redXic) > 1e-30?1./redXic:-99; 
      redXicc = h_xiccRedCounter->GetBinContent(counter); 
      redXicc = TMath::Abs(redXicc) > 1e-30?1./redXicc:-99; 

      sumHist = (TH1D*)xiHist->Clone(TString::Format("sumHist_%s", xiHist->GetName())); 
      sumHist->Add(xicHist); 
      sumHist->Add(xiccHist); 
      sumHist->SetTitle("Sum"); 
      sumHist->GetXaxis()->SetTitle("IM(#Xi_{c}^{+},#pi^{+}) (GeV/#it{c}^{2})");
      sumHist->GetYaxis()->SetTitle("Counts relative to #Xi_{cc}^{++}");       
      
      sumHist->GetYaxis()->SetRangeUser(0, sumHist->GetMaximum()*6); 
      sumHist->GetYaxis()->SetTitleOffset(1.4); 
      sumHist->GetXaxis()->SetNdivisions(506); 
      sumHist->GetYaxis()->SetNdivisions(504); 
      
      sumHist->SetLineColor(kOrange+7); 
      sumHist->SetLineStyle(2);

      sumHistBkg = (TH1D*)xiHist->Clone(TString::Format("sumHist_%s", xiHist->GetName())); 
      sumHistBkg->Add(xicHist); 
      sumHistBkg->SetTitle("Sum Background"); 
      sumHistBkg->GetYaxis()->SetTitle("Counts relative to #Xi_{cc}^{++}");       
      //TString::Format("#splitline{Sum Background}{(Back. Red: #Xi^{#minus} = %f, #Xi_{c}^{+} = %f)}", redXi, redXic)
      sumHistBkg->SetLineColor(kAzure-3); 
      
      
      xiHist->SetTitle("#Xi^{#minus} + #pi_{Pythia} + #pi_{Pythia} + #pi_{Pythia}"); 
      xicHist->SetTitle("#Xi^{+}_{c} + #pi_{Pythia}"); 
      //TString::Format("#Xi_{cc}^{++} Signal Red. by %f",redXicc)); 
      xiccHist->SetTitle("#Xi_{cc}^{++}");
      xiccHist->SetFillColor(kGreen+3);
      xiccHist->SetFillStyle(1001); 

      //get some number bby 
      double nxi = xiHist->Integral(xiHist->FindBin(xiccMass-xiccWindow),xiHist->FindBin(xiccMass+xiccWindow)); 
      double nxic = xicHist->Integral(xicHist->FindBin(xiccMass-xiccWindow),xicHist->FindBin(xiccMass+xiccWindow)); 
      double nxicc = xiccHist->Integral(xiccHist->FindBin(xiccMass-xiccWindow),xiccHist->FindBin(xiccMass+xiccWindow)); 
      std::cout <<"For " << obj->GetName() << ": \n nxi=" << nxi << " nxic="<< nxic << " nxicc=" << nxicc << " S/B=" << nxicc/(nxi+nxic) << std::endl;
      counter++; 
    }
    auto c = c11(obj->GetName()); 
    auto p = (TPad*)gROOT->FindObject(TString::Format("p%s",obj->GetName()).Data()); 
    p->cd(); 
    
    auto leg = new TLegend(0.16, 0.45, 0.56, 0.8, "#splitline{ALICE 3 Full Simluation}{Pythia pp #sqrt{s} = 13 TeV + GEANT3}");
    leg->SetFillStyle(0); 

    if (sumHist) {
      sumHist->Draw("hist"); 
      sumHistBkg->Draw("histsame"); 
      xiccHist->Draw("histsame"); 
      sumHist->Draw("samehist"); 
      leg->AddEntry(xiccHist, xiccHist->GetTitle(), "l");
      leg->AddEntry(sumHistBkg, sumHistBkg->GetTitle(), "l");
      leg->AddEntry(sumHist, sumHist->GetTitle(), "l");      
      
      TLatex* myTex = GenTex(); 
      myTex->DrawLatex(0.63, 0.55, TString::Format("Cut Variation %d", counter).Data());
      myTex->DrawLatex(0.63, 0.48, TString::Format("Reduction Factors:").Data());
      if (redXi > 0) {
	myTex->DrawLatex(0.63, 0.41, TString::Format("#Xi^{-}:%1.2e", redXi).Data());
      } else { 
	myTex->DrawLatex(0.63, 0.41, TString::Format("#Xi^{-}: Complete Red.").Data()); 
      }
      if (redXic > 0) {
	myTex->DrawLatex(0.63, 0.34, TString::Format("#Xi_{c}^{+}:%1.2e", redXic).Data());
      } else { 
	myTex->DrawLatex(0.63, 0.34, TString::Format("#Xi_{c}^{+}: Complete Red.").Data()); 
      }
      if (redXicc > 0) {
	myTex->DrawLatex(0.63, 0.26, TString::Format("#Xi_{cc}^{++}:%.1f", redXicc).Data());
      } else { 
	myTex->DrawLatex(0.63, 0.26, TString::Format("#Xi_{cc}^{++}: Complete Red.").Data()); 
      }
    } else { 
      leg->AddEntry(xiHist, xiHist->GetTitle(), "l");
      leg->AddEntry(xicHist, xicHist->GetTitle(), "l");
      leg->AddEntry(xiccHist, xiccHist->GetTitle(), "l");
      xiHist->Draw("hist"); 
      xicHist->Draw("histsame"); 
      xiccHist->Draw("histsame"); 
    }
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
