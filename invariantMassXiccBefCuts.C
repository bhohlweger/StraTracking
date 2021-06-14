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
  double normMb = 1/(1.614e-8*h_mbCounter->GetBinContent(1)); 
  
  TH1D* h_xiRedCounter = (TH1D*)xi->Get("cutCounter");
  TH1D* h_xiCounter = (TH1D*)xi->Get("df_xi_c_candCounter");
  double normXi = 1/(1.63e-6*h_xiCounter->GetBinContent(1)); 

  TH1D* h_xicRedCounter = (TH1D*)xic->Get("cutCounter");  
  TH1D* h_xicCounter = (TH1D*)xic->Get("df_xi_c_candCounter");
  double normXic = 1/(2.1e-4*h_xicCounter->GetBinContent(1)); 
  
  TH1D* h_xiccRedCounter = (TH1D*)xicc->Get("cutCounter");
  TH1D* h_xiccCounter = (TH1D*)xicc->Get("df_xi_c_candCounter");
  double normXicc = (1.)/(h_xiccCounter->GetBinContent(1)); 
  
  TH1D* sumHist = nullptr;
  TH1D* sumHistBkg = nullptr;
  TH1D* avgBkg; 
  int counter = 0; 
  
  TString histName = "h_df_in_im_xi_cc_mass_stra";

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
   
  auto c1 = c11Leg("1"); 
  auto p1 = (TPad*)gROOT->FindObject("p11"); 
  auto p2 = (TPad*)gROOT->FindObject("p21"); 
  p2->SetLeftMargin(0); 
  p2->SetRightMargin(p2->GetRightMargin()*2); 
  c1->cd(); 
  p1->Draw(); 
  p2->Draw(); 
  
  mbHist->SetTitle("Pythia MB");   
  mbHist->SetLineColor(kAzure-3); 
  mbHist->SetMarkerColor(kAzure-3); 

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

  c1->Write();
  c1->Close(); 
  
  output->cd();     
  xiccHist->Write(TString::Format("xicc_%s",xiccHist->GetName())); 
  xicHist->Write(TString::Format("xic_%s",xicHist->GetName())); 
  xiHist->Write(TString::Format("xi_%s",xiHist->GetName())); 
  mbHist->Write(TString::Format("mb_%s",mbHist->GetName())); 
  
  output->Close(); 
  mb->Close(); 
  xi->Close(); 
  xic->Close(); 
  xicc->Close(); 
}
