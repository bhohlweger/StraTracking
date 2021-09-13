void invariantMassXicc(TString addon) { 
  double xiccMass = 3.621; 
  double xiccWindow = 0.030;
  
  TFile* xi = TFile::Open(    TString::Format("outxiccSelector_xi%s.root"  , addon.Data()), "read"); 
  TFile* xic = TFile::Open(   TString::Format("outxiccSelector_xic%s.root" , addon.Data()), "read"); 
  TFile* xicc = TFile::Open(  TString::Format("outxiccSelector_xicc%s.root", addon.Data()), "read"); 
  TFile* output = TFile::Open(TString::Format("outIMXicc%s.root"         , addon.Data()), "recreate"); 
  
  TH1D* h_xiRedCounter = (TH1D*)xi->Get("cutCounter");

  TH1D* h_xiCounter = (TH1D*)xi->Get("df_xi_c_candCounter");
  std::cout << h_xiCounter->GetBinContent(1) << std::endl;
  TH1D* h_xicCounter = (TH1D*)xic->Get("df_xi_c_candCounter");
  std::cout << h_xicCounter->GetBinContent(1) << std::endl;
  TH1D* h_xiccCounter = (TH1D*)xicc->Get("df_xi_c_candCounter");
  std::cout << h_xiccCounter->GetBinContent(1) << std::endl;

  TH2D* h_xipteta = (TH2D*)xi->Get("ptetaXiGen"); 
  double etaCut = 1.5; 
  int xi_eta_min = h_xipteta->GetYaxis()->FindBin(-etaCut); 
  int xi_eta_max = h_xipteta->GetYaxis()->FindBin(etaCut); 
  int nEvtsXi = h_xipteta->ProjectionY()->Integral(xi_eta_min, xi_eta_max); 
  std::cout <<" nEvts xi: " << nEvtsXi << std::endl;
  double normXi = 1/(1.63e-6*nEvtsXi); 

  TH1D* h_xicRedCounter = (TH1D*)xic->Get("cutCounter");  
  TH2D* h_xicpteta = (TH2D*)xic->Get("ptetaXicGen"); 

  int xic_eta_min = h_xicpteta->GetYaxis()->FindBin(-etaCut); 
  int xic_eta_max = h_xicpteta->GetYaxis()->FindBin(etaCut); 
  int nEvtsXic = h_xicpteta->ProjectionY()->Integral(xic_eta_min, xic_eta_max); 
  std::cout <<" nEvts xic: " << nEvtsXic << std::endl;
  double normXic = 1/(2.1e-4*nEvtsXic); 
  
  TH1D* h_xiccRedCounter = (TH1D*)xicc->Get("cutCounter");
  TH2D* h_xiccpteta = (TH2D*)xicc->Get("ptetaXiccGen"); 
    
  int xicc_eta_min = h_xiccpteta->GetYaxis()->FindBin(-etaCut); 
  int xicc_eta_max = h_xiccpteta->GetYaxis()->FindBin(etaCut); 
  int nEvtsXicc = h_xiccpteta->ProjectionY()->Integral(xicc_eta_min, xicc_eta_max); 
  std::cout <<" nEvts xicc: " << nEvtsXicc << std::endl;
    
  TH1D* xiccGenPt = (TH1D*)xicc->Get("pTXiccGenerados"); 
  xiccGenPt->Sumw2(); 
  //xiccGenPt->Rebin(4);
  output->cd(); 
  xiccGenPt->Write(); 

  double normXicc = (1.)/(nEvtsXicc); 
  
  TList* inList = xi->GetListOfKeys(); 
  TIter next(inList); 
  TObject* obj = nullptr; 
  TH1D* sumHist = nullptr;
  TH1D* sumHistBkg = nullptr;
  TH1D* avgBkg; 
  int counter = -4;  // four non tracked histos
  //double redXi; 
  //double redXic;
  //double redXicc; 

  while ((obj = next())) {
    TString objName = TString::Format("%s",obj->GetName()); 
    if (objName.Contains("Counter")  || objName.Contains("_vs_")|| !(objName.Contains("mass") && objName.Contains("xi_cc"))) { 
      continue; 
    }
    
    if (sumHist || sumHistBkg) { 
      sumHist = nullptr;
      sumHistBkg = nullptr;
      avgBkg = nullptr;
    }
    std::cout <<"=========================================\n" << obj->GetName() << "\n=========================================\n"; 
    TH1D* xiHist = (TH1D*)xi->Get(obj->GetName()); 
    TH1D* xicHist = (TH1D*)xic->Get(obj->GetName()); 
    TH1D* xiccHist = (TH1D*)xicc->Get(obj->GetName()); 

    if (!xiHist) { 
      std::cout << "xi hist missing \n"; 
    }
    if (!xicHist) { 
      std::cout << "xic hist missing \n"; 
    }
    if (!xiccHist) { 
      std::cout << "xicc hist missing \n"; 
    }
    
    TString toReplace = "mass"; 
    TString replacement = "pt_vs_eta"; 
    
    TString pTIdent = objName.ReplaceAll(toReplace, toReplace.Length(), replacement, replacement.Length()); 
    std::cout << "pTIdent: " << pTIdent.Data() << std::endl; 
    TH2D* pT_vs_eta = (TH2D*)xicc->Get(pTIdent.Data()); 
    TH1D* pT_Identified = nullptr; 
    if (pT_vs_eta && pTIdent.Contains("df_xi_cc_im_xi_cc_pt_vs_eta_")) { 
      pT_Identified = pT_vs_eta->ProjectionX(TString::Format("Efficiency_%s", objName.Data())); 
      pT_Identified->Sumw2();
      pT_Identified->Rebin(4);
      pT_Identified->Divide(xiccGenPt); 
    }
    xiHist->Rebin(4);       
    xicHist->Rebin(4); 
    xiccHist->Rebin(4); 

    double signalCounts = xiccHist->Integral(xiccHist->FindBin(xiccMass-xiccWindow), xiccHist->FindBin(xiccMass+xiccWindow)); 
    double bkgCounts = xicHist->Integral(xicHist->FindBin(xiccMass-xiccWindow), xicHist->FindBin(xiccMass+xiccWindow)) + xiHist->Integral(xiHist->FindBin(xiccMass-xiccWindow), xiHist->FindBin(xiccMass+xiccWindow)); 
    
    double relErrSignal = TMath::Sqrt(signalCounts)/signalCounts; 
    double relErrBkg = TMath::Sqrt(bkgCounts)/bkgCounts; 
    
    xiHist->Scale(normXi); 
    xicHist->Scale(normXic); 
    xiccHist->Scale(normXicc); 

    double signalCountsScaled = xiccHist->Integral(xiccHist->FindBin(xiccMass-xiccWindow), xiccHist->FindBin(xiccMass+xiccWindow)); 
    double bkgCountsScaled = xicHist->Integral(xicHist->FindBin(xiccMass-xiccWindow), xicHist->FindBin(xiccMass+xiccWindow)) + xiHist->Integral(xiHist->FindBin(xiccMass-xiccWindow), xiHist->FindBin(xiccMass+xiccWindow)); 
    std::cout << "Signal Efficiency pT integrated = " << TString::Format("%.3e",signalCounts*normXicc).Data() << std::endl; 
    std::cout << "Signal counts = " << signalCounts << " rel. Err. = " << relErrSignal << "\nSignal counts after scaling = " << signalCountsScaled << " Abs. Err. = " << relErrSignal*signalCountsScaled << std::endl; 
    std::cout << "Bkg counts = " << bkgCounts << " rel. Err. = " << relErrBkg << "\nBkg counts after scaling = " << bkgCountsScaled << " Abs Err. = " << relErrBkg*bkgCountsScaled << std::endl; 
    
    /*
    redXi = h_xiRedCounter->GetBinContent(counter); 
    redXi = TMath::Abs(redXi) > 1e-30?1./redXi:-99; 
    redXic = h_xicRedCounter->GetBinContent(counter);
    redXic = TMath::Abs(redXic) > 1e-30?1./redXic:-99; 
    redXicc = h_xiccRedCounter->GetBinContent(counter); 
    redXicc = TMath::Abs(redXicc) > 1e-30?1./redXicc:-99; 
    */
    counter++; 
    
    sumHist = (TH1D*)xiHist->Clone(TString::Format("sumHist_%s", xiHist->GetName())); 
    sumHist->Add(xicHist); 
    sumHist->Add(xiccHist); 
    
    sumHistBkg = (TH1D*)xiHist->Clone(TString::Format("sumHist_%s", xiHist->GetName())); 
    sumHistBkg->Add(xicHist); 
    
    avgBkg = (TH1D*)sumHistBkg->Clone(TString::Format("avg%s", sumHistBkg->GetName()).Data()); 
    
    int whiteBins = 0; 
    TH1D* myavgBkg = (TH1D*)sumHistBkg->Clone(TString::Format("myavg%s", sumHistBkg->GetName()).Data()); 
    for (int iBin = 1; iBin < myavgBkg->GetNbinsX(); ++iBin) { 
      double cent = myavgBkg->GetBinCenter(iBin); 
      if (cent < 3 || cent > 4) { 
	myavgBkg->SetBinContent(iBin,0); 
      } else { 
	whiteBins++; 
      }
    }
    myavgBkg->Rebin(myavgBkg->GetNbinsX()); 
    myavgBkg->Scale(1./whiteBins); 
      
    for (int iBin = 1; iBin < avgBkg->GetNbinsX(); ++iBin) { 
      avgBkg->SetBinContent(iBin,myavgBkg->GetBinContent(1)); 
    }

    auto leg = new TLegend(0.16, 0.38, 0.56, 0.65);//, "#splitline{ALICE 3 Full Simulation}{Pythia pp #sqrt{s} = 13 TeV + GEANT3}");
    leg->SetFillStyle(0); 
    
    
    auto leg2 = new TLegend(0.16, 0.38, 0.56, 0.65);//, "#splitline{ALICE 3 Full Simulation}{Pythia pp #sqrt{s} = 13 TeV + GEANT3}");

    auto c1 = c11(TString::Format("%s1",obj->GetName()).Data()); 
    auto p1 = (TPad*)gROOT->FindObject(TString::Format("p%s1",obj->GetName()).Data()); 
    
    auto c2 = c11(TString::Format("%s2",obj->GetName()).Data()); 
    auto p2 = (TPad*)gROOT->FindObject(TString::Format("p%s2",obj->GetName()).Data()); 
    
    TLine *HoldTheLine = new TLine();
    HoldTheLine->SetLineColor(kBlack); 
    HoldTheLine->SetLineWidth(4); 
    HoldTheLine->SetLineStyle(2); 
        
    xiHist->SetLineColor(kPink+7); 
    xiHist->SetMarkerColor(kPink+7); 
    xicHist->SetLineColor(38); 
    xicHist->SetMarkerColor(38); 
    xiccHist->SetLineColor(kGreen+3); 
    xiccHist->SetMarkerColor(kGreen+3); 
    
    sumHist->SetTitle("Sum");     
    sumHist->GetXaxis()->SetTitle("IM(#Xi_{c}^{+},#pi^{+}) (GeV/#it{c}^{2})");
    sumHist->GetYaxis()->SetTitle("Counts relative to #Xi_{cc}^{++}");       
    sumHist->SetLineColor(kOrange+7); 
    sumHist->SetLineStyle(2);
    
    sumHist->GetXaxis()->SetRangeUser(3., 4.); 
    sumHist->GetYaxis()->SetRangeUser(0, sumHist->GetMaximum()*6); 
    sumHist->GetYaxis()->SetTitleOffset(1.4); 
    sumHist->GetXaxis()->SetNdivisions(506); 
    sumHist->GetYaxis()->SetMaxDigits(3); 
    sumHist->GetYaxis()->SetNdivisions(504); 
      
    sumHistBkg->SetTitle("Sum Background"); 
    sumHistBkg->GetYaxis()->SetTitle("Counts relative to #Xi_{cc}^{++}");       
    sumHistBkg->SetLineColor(kAzure-3); 
          
    avgBkg->SetTitle("Average Background"); 
    avgBkg->GetXaxis()->SetTitle("IM(#Xi_{c}^{+},#pi^{+}) (GeV/#it{c}^{2})");
    avgBkg->GetYaxis()->SetTitle("Counts relative to #Xi_{cc}^{++}");       
      
    avgBkg->GetXaxis()->SetRangeUser(3., 4.); 
    avgBkg->GetYaxis()->SetRangeUser(0, avgBkg->GetMaximum()*6); 
    avgBkg->GetYaxis()->SetTitleOffset(1.4); 
    avgBkg->GetXaxis()->SetNdivisions(506); 
    avgBkg->GetYaxis()->SetNdivisions(504); 
      
    avgBkg->SetLineColor(kOrange+7); 
    avgBkg->SetLineStyle(2);

    xiHist->SetTitle("#Xi^{#minus} + 3#times#pi_{Pythia}"); 
    xicHist->SetTitle("#Xi^{+}_{c} + #pi_{Pythia}"); 
    xiccHist->SetTitle("#Xi_{cc}^{++}");
    xiccHist->SetFillColor(kGreen+3);
    xiccHist->SetFillStyle(1001); 

    leg->AddEntry(xiccHist, xiccHist->GetTitle(), "l");
    leg->AddEntry(xicHist, xicHist->GetTitle(), "l");
    leg->AddEntry(xiHist, xiHist->GetTitle(), "l");
    leg->AddEntry(sumHist, sumHist->GetTitle(), "l");      

    leg2->SetFillStyle(0); 
    leg2->AddEntry(xiccHist, xiccHist->GetTitle(), "l");
    leg2->AddEntry(avgBkg, avgBkg->GetTitle(), "l");

    //get some number bby 
    double nxi = xiHist->Integral(xiHist->FindBin(xiccMass-xiccWindow),xiHist->FindBin(xiccMass+xiccWindow)); 
    double nxic = xicHist->Integral(xicHist->FindBin(xiccMass-xiccWindow),xicHist->FindBin(xiccMass+xiccWindow)); 
    double nxicc = xiccHist->Integral(xiccHist->FindBin(xiccMass-xiccWindow),xiccHist->FindBin(xiccMass+xiccWindow)); 
    double nbkgavg = avgBkg->Integral(avgBkg->FindBin(xiccMass-xiccWindow),avgBkg->FindBin(xiccMass+xiccWindow)); 
    std::cout <<"nxi = " << nxi << " nxic = "<< nxic << " nxicc = " << nxicc << " S/B = " << nxicc/(nxi+nxic) << std::endl;
    std::cout <<"Average Bakcground values: \n" << "nbkgavg = " << nbkgavg << "S/B = " << nxicc/nbkgavg << std::endl; 

    p1->cd(); 
    
    sumHist->Draw("hist"); 
    //sumHistBkg->Draw("histsame"); 
    xiccHist->Draw("same"); 
    xiccHist->Draw("histsame"); 
    xicHist->Draw("same"); 
    xiHist->Draw("same"); 
    sumHist->Draw("samehist"); 
    leg->Draw("same"); 
    HoldTheLine->DrawLine(xiccMass-xiccWindow, 0, xiccMass-xiccWindow, xiccHist->GetMaximum()*1.5); 
    HoldTheLine->DrawLine(xiccMass+xiccWindow, 0, xiccMass+xiccWindow, xiccHist->GetMaximum()*1.5); 
    
    TLatex* myTex = GenTex(); 
    //myTex->DrawLatex(0.61, 0.55, TString::Format("Cut Variation %d", counter).Data());
    /*
    myTex->DrawLatex(0.63, 0.55, TString::Format("Reduction Factors:").Data());
    if (redXi > 0) {
      myTex->DrawLatex(0.63, 0.48, TString::Format("#Xi^{-}:  %1.2e", redXi).Data());
    } else { 
      myTex->DrawLatex(0.63, 0.48, TString::Format("#Xi^{-}:  Complete Red.").Data()); 
    }
    if (redXic > 0) {
      myTex->DrawLatex(0.63, 0.41, TString::Format("#Xi_{c}^{+}: %1.2e", redXic).Data());
    } else { 
      myTex->DrawLatex(0.63, 0.41, TString::Format("#Xi_{c}^{+}: Complete Red.").Data()); 
    }
    if (redXicc > 0) {
      myTex->DrawLatex(0.63, 0.34, TString::Format("#Xi_{cc}^{++}: %.1f", redXicc).Data());
    } else { 
      myTex->DrawLatex(0.63, 0.34, TString::Format("#Xi_{cc}^{++}: Complete Red.").Data()); 
    }
    */
    myTex->DrawLatex(0.18,0.75,"#splitline{ALICE 3 Study (Layout v1)}{#splitline{Full Simulation}{Pythia pp #sqrt{s} = 13 TeV + GEANT3}}");
  
    p2->cd(); 
    avgBkg->Draw("hist"); 
    xiccHist->Draw("sameHist"); 
    leg2->Draw("same"); 
    myTex->DrawLatex(0.18,0.75,"#splitline{ALICE 3 Study (Layout v1)}{#splitline{Full Simulation}{Pythia pp #sqrt{s} = 13 TeV + GEANT3}}");
    

    c1->Write();
    c1->Close(); 
    c2->Write();      
    c2->Close();
    output->cd();     
    myavgBkg->Write();
    avgBkg->Write(); 
    if (pT_vs_eta) pT_vs_eta->Write(); 
    if (pT_Identified) pT_Identified->Write();
    xiccHist->Write(TString::Format("xicc_%s",xiccHist->GetName())); 
    xicHist->Write(TString::Format("xic_%s",xiccHist->GetName())); 
    xiHist->Write(TString::Format("xi_%s",xiccHist->GetName())); 

  }
  output->Close(); 
  xi->Close(); 
  xic->Close(); 
  xicc->Close(); 
}
