void invariantMassXiccPbPb(TString addon) { 
  double xiccMass = 3.621; 
  double xiccWindow = 0.12;

  TFile* mb = TFile::Open(    TString::Format("outxiccSelector_mb%s.root"  , addon.Data()), "read"); 
  TFile* xicc = TFile::Open(  TString::Format("outxiccSelector_xicc%s.root", addon.Data()), "read"); 
  TFile* output = TFile::Open(TString::Format("outIMXiccPBPB%s.root"         , addon.Data()), "recreate"); 

  double etaCut = 1.5; 


  TH1D* h_mbRedCounter = (TH1D*)mb->Get("cutCounter");
  TH1D* h_mbCounter = (TH1D*)mb->Get("df_xi_c_candCounter");
  int nEvtsMb = h_mbCounter->GetBinContent(1); 
  std::cout <<" nEvts mb: " << nEvtsMb << std::endl;
  double normMb = 137/(0.05*nEvtsMb); //missing branching ratio ...
  
  TH1D* h_xiccCounter = (TH1D*)xicc->Get("df_xi_c_candCounter");
  std::cout << h_xiccCounter->GetBinContent(1) << std::endl;
  
  TH1D* h_xiccRedCounter = (TH1D*)xicc->Get("cutCounter");
  TH2D* h_xiccpteta = (TH2D*)xicc->Get("ptetaXiccGen"); 
    
  int xicc_eta_min = h_xiccpteta->GetYaxis()->FindBin(-etaCut); 
  int xicc_eta_max = h_xiccpteta->GetYaxis()->FindBin(etaCut); 
  int nEvtsXicc = h_xiccpteta->ProjectionY()->Integral(xicc_eta_min, xicc_eta_max); 
  std::cout <<" nEvts xicc: " << nEvtsXicc << std::endl;
    
  TH1D* xiccGenPt = (TH1D*)h_xiccpteta->ProjectionX("pTXiccGenerados",xicc_eta_min, xicc_eta_max);
  xiccGenPt->Sumw2(); 
  output->cd(); 
  xiccGenPt->Rebin(4); 
  xiccGenPt->Write(); 
  
  TH1D* xiccidentbefcuts = (TH1D*)xicc->Get("df_xi_c_qa_xi_cc_pt"); 
  xiccidentbefcuts->Write("identficados");
  xiccidentbefcuts->Divide(xiccGenPt);
  xiccidentbefcuts->Write("effciendados");

  

  double normXicc = (1.)/(nEvtsXicc); 
  
  TList* inList = mb->GetListOfKeys(); 
  TIter next(inList); 
  TObject* obj = nullptr; 
  TH1D* sumHist = nullptr;
  TH1D* mbHist = nullptr;
  TH1D* avgBkg; 
  int counter = -4;  // four non tracked histos
  /*
  double redMB;
  double redXicc; 
  */
  while ((obj = next())) {
    TString objName = TString::Format("%s",obj->GetName()); 
    if (objName.Contains("Counter")  || objName.Contains("_vs_")|| !(objName.Contains("mass") && objName.Contains("xi_cc"))) { 
      continue; 
    }
    
    if (sumHist || mbHist) { 
      sumHist = nullptr;
      mbHist = nullptr;
      avgBkg = nullptr;
    }
    std::cout <<"=========================================\n" << obj->GetName() << "\n=========================================\n"; 
    TH1D* mbHist = (TH1D*)mb->Get(obj->GetName()); 
    TH1D* xiccHist = (TH1D*)xicc->Get(obj->GetName()); 

    TString toReplace = "mass_stra"; 
    TString replacement = "pt_vs_y"; 
    
    TString pTIdent = objName.ReplaceAll(toReplace, toReplace.Length(), replacement, replacement.Length()); 
    
    TH2D* pT_vs_eta = (TH2D*)xicc->Get(pTIdent.Data()); 
    TH1D* pT_Identified = nullptr; 
    if (pT_vs_eta) { 
      pT_Identified = pT_vs_eta->ProjectionX(TString::Format("Efficiency_%s", objName.Data())); 
      pT_Identified->Sumw2();
      output->cd();
      pT_Identified->Rebin(4); 
      pT_Identified->Write(TString::Format("pTSpectrum_%s", objName.Data())); 
      pT_Identified->Divide(xiccGenPt); 
    }
    mbHist->Rebin(4); 
    xiccHist->Rebin(4); 

    double signalCounts = xiccHist->Integral(xiccHist->FindBin(xiccMass-xiccWindow), xiccHist->FindBin(xiccMass+xiccWindow)); 
    double bkgCounts = mbHist->Integral(mbHist->FindBin(xiccMass-xiccWindow), mbHist->FindBin(xiccMass+xiccWindow)); 
    
    double relErrSignal = TMath::Sqrt(signalCounts)/signalCounts; 
    double relErrBkg = TMath::Sqrt(bkgCounts)/bkgCounts; 
    
    mbHist->Scale(normMb); 
    xiccHist->Scale(normXicc); 

    double signalCountsScaled = xiccHist->Integral(xiccHist->FindBin(xiccMass-xiccWindow), xiccHist->FindBin(xiccMass+xiccWindow)); 
    double bkgCountsScaled = mbHist->Integral(mbHist->FindBin(xiccMass-xiccWindow), mbHist->FindBin(xiccMass+xiccWindow));

    std::cout << "Signal Efficiency pT integrated = " << TString::Format("%.3e",signalCounts*normXicc).Data() << std::endl; 
    std::cout << "Signal counts = " << signalCounts << " rel. Err. = " << relErrSignal << "\nSignal counts after scaling = " << signalCountsScaled << " Abs. Err. = " << relErrSignal*signalCountsScaled << std::endl; 
    std::cout << "Bkg counts = " << bkgCounts << " rel. Err. = " << relErrBkg << "\nBkg counts after scaling = " << bkgCountsScaled << " Abs Err. = " << relErrBkg*bkgCountsScaled << std::endl; 
    /*
    redMB = h_mbRedCounter->GetBinContent(counter);
    redMB = TMath::Abs(redMB) > 1e-30?1./redMB:-99; 
    redXicc = h_xiccRedCounter->GetBinContent(counter); 
    redXicc = TMath::Abs(redXicc) > 1e-30?1./redXicc:-99; 
    counter++; 
    */
    sumHist = (TH1D*)mbHist->Clone(TString::Format("sumHist_%s", mbHist->GetName())); 
    sumHist->Add(xiccHist); 
    
    avgBkg = (TH1D*)mbHist->Clone(TString::Format("avg%s", mbHist->GetName()).Data()); 
    
    int whiteBins = 0; 
    TH1D* myavgBkg = (TH1D*)mbHist->Clone(TString::Format("myavg%s", mbHist->GetName()).Data()); 
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
        
    mbHist->SetLineColor(38); 
    mbHist->SetMarkerColor(38); 
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
      
    mbHist->SetTitle("MB"); 
    mbHist->GetYaxis()->SetTitle("Counts relative to #Xi_{cc}^{++}");       
    mbHist->SetLineColor(kAzure-3); 
          
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

    xiccHist->SetTitle("#Xi_{cc}^{++}");
    xiccHist->SetFillColor(kGreen+3);
    xiccHist->SetFillStyle(1001); 

    leg->AddEntry(xiccHist, xiccHist->GetTitle(), "l");
    leg->AddEntry(mbHist, mbHist->GetTitle(), "l");
    //leg->AddEntry(sumHist, sumHist->GetTitle(), "l");      

    leg2->SetFillStyle(0); 
    leg2->AddEntry(xiccHist, xiccHist->GetTitle(), "l");
    leg2->AddEntry(avgBkg, avgBkg->GetTitle(), "l");

    //get some number bby 
    double nmb = mbHist->Integral(mbHist->FindBin(xiccMass-xiccWindow),mbHist->FindBin(xiccMass+xiccWindow)); 
    double nxicc = xiccHist->Integral(xiccHist->FindBin(xiccMass-xiccWindow),xiccHist->FindBin(xiccMass+xiccWindow)); 
    double nbkgavg = avgBkg->Integral(avgBkg->FindBin(xiccMass-xiccWindow),avgBkg->FindBin(xiccMass+xiccWindow)); 
    std::cout <<" nmb = "<< nmb << " nxicc = " << nxicc << " S/B = " << nxicc/(double)nmb << std::endl;
    std::cout <<"Average Background values: \n" << "nbkgavg = " << nbkgavg << "S/B = " << nxicc/nbkgavg << std::endl; 

    p1->cd(); 
    sumHist->Draw("hist");      
    xiccHist->Draw("histsame"); 
    mbHist->Draw("same"); 
    leg->Draw("same"); 
    HoldTheLine->DrawLine(xiccMass-xiccWindow, 0, xiccMass-xiccWindow, xiccHist->GetMaximum()*1.5); 
    HoldTheLine->DrawLine(xiccMass+xiccWindow, 0, xiccMass+xiccWindow, xiccHist->GetMaximum()*1.5); 
    
    TLatex* myTex = GenTex(); 
    /*
    if (counter > 0) { 
      //myTex->DrawLatex(0.61, 0.55, TString::Format("Cut Variation %d", counter).Data());
      myTex->DrawLatex(0.63, 0.55, TString::Format("Reduction Factors:").Data());
      if (redMB > 0) {
	myTex->DrawLatex(0.63, 0.48, TString::Format("MB: %1.2e", redMB).Data());
      } else { 
	myTex->DrawLatex(0.63, 0.48, TString::Format("MB: Complete Red.").Data()); 
      }
      if (redXicc > 0) {
	myTex->DrawLatex(0.63, 0.41, TString::Format("#Xi_{cc}^{++}: %.1f", redXicc).Data());
      } else { 
	myTex->DrawLatex(0.63, 0.41, TString::Format("#Xi_{cc}^{++}: Complete Red.").Data()); 
      }
    } 
    */
    myTex->DrawLatex(0.18,0.75,"#splitline{ALICE 3 Study (Layout v1)}{#splitline{Full Simulation}{Pythia pp #sqrt{s} = 13 TeV + GEANT3}}");
  
    p2->cd(); 
    avgBkg->Draw("hist"); 
    xiccHist->Draw("sameHist"); 
    leg2->Draw("same"); 
    myTex->DrawLatex(0.18,0.75,"#splitline{ALICE 3 Study (Layout v1)}{#splitline{Full Simulation}{Pythia pp #sqrt{s} = 13 TeV + GEANT3}}");
    
    output->cd();
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
    mbHist->Write(TString::Format("xic_%s",xiccHist->GetName())); 

  }
  output->Close(); 
  mb->Close(); 
  xicc->Close(); 
}
