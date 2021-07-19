void invariantMassXiccBefCuts(TString addon) { 
  double xiccMass = 3.596; //3.621; 
  double xiccWindow = 0.08;

  TFile* mb = TFile::Open(    TString::Format("outxiccSelector_mb%s.root"  , addon.Data()), "read"); 
  TFile* xi = TFile::Open(    TString::Format("outxiccSelector_xi%s.root"  , addon.Data()), "read"); 
  TFile* xic = TFile::Open(   TString::Format("outxiccSelector_xic%s.root" , addon.Data()), "read"); 
  TFile* xicc = TFile::Open(  TString::Format("outxiccSelector_xicc%s.root", addon.Data()), "read"); 
  TFile* output = TFile::Open(TString::Format("outIMXiccBefCuts%s.root"         , addon.Data()), "recreate"); 
  
  TH1D* h_mbRedCounter = (TH1D*)mb->Get("cutCounter");

  TH1D* h_mbCounter = (TH1D*)mb->Get("df_xi_c_candCounter");
  int nEvtsMb = h_mbCounter->GetBinContent(1); 
  std::cout <<" nEvts mb: " << nEvtsMb << std::endl;
  double normMb = 1/(1.614e-8*nEvtsMb); 

  TH2D* h_xipteta = (TH2D*)xi->Get("ptetaXiGen"); 
  TH1D* h_xiRedCounter = (TH1D*)xi->Get("cutCounter");
  
  int xi_eta_min = h_xipteta->GetYaxis()->FindBin(-0.5); 
  int xi_eta_max = h_xipteta->GetYaxis()->FindBin(0.5); 
  int nEvtsXi = h_xipteta->ProjectionY()->Integral(xi_eta_min, xi_eta_max); 
  std::cout <<" nEvts xi: " << nEvtsXi << std::endl;
  double normXi = 1/(1.63e-6*nEvtsXi); 

  TH1D* h_xicRedCounter = (TH1D*)xic->Get("cutCounter");  
  TH2D* h_xicpteta = (TH2D*)xic->Get("ptetaXicGen"); 

  int xic_eta_min = h_xicpteta->GetYaxis()->FindBin(-0.5); 
  int xic_eta_max = h_xicpteta->GetYaxis()->FindBin(0.5); 
  int nEvtsXic = h_xicpteta->ProjectionY()->Integral(xic_eta_min, xic_eta_max); 
  std::cout <<" nEvts xic: " << nEvtsXic << std::endl;
  double normXic = 1/(2.1e-4*nEvtsXic); 
  
  TH1D* h_xiccRedCounter = (TH1D*)xicc->Get("cutCounter");
  TH2D* h_xiccpteta = (TH2D*)xicc->Get("ptetaXiccGen"); 
    
  int xicc_eta_min = h_xiccpteta->GetYaxis()->FindBin(-0.5); 
  int xicc_eta_max = h_xiccpteta->GetYaxis()->FindBin(0.5); 
  int nEvtsXicc = h_xiccpteta->ProjectionY()->Integral(xicc_eta_min, xicc_eta_max); 
  std::cout <<" nEvts xicc: " << nEvtsXicc << std::endl;
    
  TH1D* xiccGenPt = (TH1D*)h_xiccpteta->ProjectionX("pTXiccGenerados",xicc_eta_min, xicc_eta_max);
  xiccGenPt->Sumw2(); 
  output->cd(); 
  xiccGenPt->Write(); 
  
  double normXicc = (1.)/(nEvtsXicc); 
 
  TH1D* sumHist = nullptr;
  TH1D* sumHistBkg = nullptr;
  TH1D* avgBkg; 
  int counter = 0; 
  std::vector<TString> histNames = {"h_df_in_im_xi_cc_mass_stra", "h_df_lmb_xi_cc_mass_stra", "h_df_xi_xi_cc_mass_stra", "df_xi_c_xi_cc_mass_stra", "df_xi_cc_im_xi_cc_mass_stra_c5"}; 
  //TString histName = "h_df_in_im_xi_cc_mass_stra";
  for ( auto histName : histNames) { 
    //TString histName = "df_xi_c_xi_cc_mass_stra"; 
    TH1D* mbHist = (TH1D*)mb->Get(histName.Data()); 
    TH1D* xiHist = (TH1D*)xi->Get(histName.Data()); 
    TH1D* xicHist = (TH1D*)xic->Get(histName.Data()); 
    TH1D* xiccHist = (TH1D*)xicc->Get(histName.Data()); 
  
    mbHist->Rebin(4);
    xiHist->Rebin(4);       
    xicHist->Rebin(4); 
    xiccHist->Rebin(4); 
  
    mbHist->Scale(normMb); 
    xiHist->Scale(normXi); 
    xicHist->Scale(normXic); 
    xiccHist->Scale(normXicc); 
  
    sumHist = (TH1D*)xiHist->Clone(TString::Format("sumHist_%s", xiHist->GetName())); 
    sumHist->Add(xicHist); 
    sumHist->Add(xiccHist); 
    sumHist->Add(mbHist); 

    sumHistBkg = (TH1D*)xiHist->Clone(TString::Format("sumHist_%s", xiHist->GetName())); 
    sumHistBkg->Add(xicHist); 
  
  
    auto leg = new TLegend(0., 0.29, 0.56, 0.68);
    leg->SetFillStyle(0); 
   
    auto c1 = c11Leg(histName.Data()); 
    auto p1 = (TPad*)gROOT->FindObject(TString::Format("p1%s", histName.Data())); 
    auto p2 = (TPad*)gROOT->FindObject(TString::Format("p2%s", histName.Data())); 
    p2->SetLeftMargin(0); 
    p2->SetRightMargin(p2->GetRightMargin()*2); 
    c1->cd(); 
    p1->Draw(); 
    p2->Draw(); 
  
    mbHist->SetTitle("p + 5#times#pi_{Pythia}");   
    mbHist->SetLineColor(kViolet-3); 
    mbHist->SetMarkerColor(kViolet-3); 

    xiHist->SetTitle("#Xi^{#minus} + 3#times#pi_{Pythia}"); 
    xiHist->SetLineColor(kPink+7); 
    xiHist->SetMarkerColor(kPink+7); 
  
    xicHist->SetTitle("#Xi^{+}_{c} + #pi_{Pythia}"); 
    xicHist->SetLineColor(38); 
    xicHist->SetMarkerColor(38); 

    xiccHist->SetTitle("#Xi_{cc}^{++}");
    xiccHist->SetLineColor(kGreen+3); 
    xiccHist->SetMarkerColor(kGreen+3); 
    xiccHist->SetFillColor(kGreen+3);
    xiccHist->SetFillStyle(1001); 
  
    sumHist->SetTitle("Sum");     
    sumHist->GetXaxis()->SetTitle("IM(#Xi_{c}^{+},#pi^{+}) (GeV/#it{c}^{2})");
    sumHist->GetYaxis()->SetTitle("Counts relative to #Xi_{cc}^{++}");       
    sumHist->SetLineColor(kOrange+7); 
    sumHist->SetLineStyle(2);
    
    sumHist->GetXaxis()->SetRangeUser(3., 4.); 
    sumHist->GetYaxis()->SetRangeUser(1e-5, 1e6); 
    sumHist->GetYaxis()->SetTitleOffset(1.4); 
    sumHist->GetXaxis()->SetNdivisions(506); 
    sumHist->GetXaxis()->SetMaxDigits(3); 
    //  sumHist->GetYaxis()->SetNdivisions(504); 
      
    sumHistBkg->SetTitle("Sum Background"); 
    sumHistBkg->GetYaxis()->SetTitle("Counts relative to #Xi_{cc}^{++}");       
    sumHistBkg->SetLineColor(kAzure-3); 
 
    leg->AddEntry(xiccHist, xiccHist->GetTitle(), "l");
    leg->AddEntry(xicHist, xicHist->GetTitle(), "l");
    leg->AddEntry(xiHist, xiHist->GetTitle(), "l");
    leg->AddEntry(mbHist, mbHist->GetTitle(), "l");
    leg->AddEntry(sumHist, sumHist->GetTitle(), "l");      
   
    p1->cd(); 
    p1->SetLogy();
    sumHist->Draw("hist"); 
    xiccHist->Draw("same"); 
    xiccHist->Draw("histsame"); 
    xicHist->Draw("same"); 
    xiHist->Draw("same"); 
    mbHist->Draw("same"); 
    sumHist->Draw("samehist"); 
  
    p2->cd();
    leg->Draw("same"); 
    auto myTex = GenTex(); 
    myTex->DrawLatex(0.02,0.78,"#splitline{ALICE 3 Study (Layout v1)}{#splitline{Full Simulation}{Pythia pp #sqrt{s} = 13 TeV + GEANT3}}");
    
    c1->SaveAs(TString::Format("%s.pdf", c1->GetName()).Data());
    c1->Write();
    c1->Close(); 
  
    output->cd();     
    xiccHist->Write(TString::Format("xicc_%s",xiccHist->GetName())); 
    xicHist->Write(TString::Format("xic_%s",xicHist->GetName())); 
    xiHist->Write(TString::Format("xi_%s",xiHist->GetName())); 
    mbHist->Write(TString::Format("mb_%s",mbHist->GetName())); 
  }
  output->Close(); 
  mb->Close(); 
  xi->Close(); 
  xic->Close(); 
  xicc->Close(); 
}
