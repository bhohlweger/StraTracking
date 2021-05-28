double treshold_bkg = 0.01;  
double treshold_sgn = 0.99;  


void StupidPlot() {
  
  TFile* corre = TFile::Open("outDCAProductCuts_ccc.root","read");
  //Cut Plots
  TH1D* om_c_topo  = (TH1D*)corre->Get("om_c_dca_om_dca_pi_topo_cumulative");  
  TH1D* om_cc_topo  = (TH1D*)corre->Get("om_cc_dca_om_dca_pi_topo_cumulative");  
  TH1D* om_ccc_topo  = (TH1D*)corre->Get("om_ccc_dca_om_dca_pi_topo_cumulative");  
  
  TH1D* om_c_stra  = (TH1D*)corre->Get("om_c_dca_om_dca_pi_stra_cumulative");  
  TH1D* om_cc_stra  = (TH1D*)corre->Get("om_cc_dca_om_dca_pi_stra_cumulative");  
  TH1D* om_ccc_stra  = (TH1D*)corre->Get("om_ccc_dca_om_dca_pi_stra_cumulative");  
  
  TH1D* ca_c_topo  = (TH1D*)corre->Get("ca_c_dca_om_dca_pi_topo_cumulative");  
  TH1D* ca_cc_topo  = (TH1D*)corre->Get("ca_cc_dca_om_dca_pi_topo_cumulative");  
  TH1D* ca_ccc_topo  = (TH1D*)corre->Get("ca_ccc_dca_om_dca_pi_topo_cumulative");  
  
  TH1D* ca_c_stra  = (TH1D*)corre->Get("ca_c_dca_om_dca_pi_stra_cumulative");  
  TH1D* ca_cc_stra  = (TH1D*)corre->Get("ca_cc_dca_om_dca_pi_stra_cumulative");  
  TH1D* ca_ccc_stra  = (TH1D*)corre->Get("ca_ccc_dca_om_dca_pi_stra_cumulative");  
  
    
  auto c4 = c11("4"); 
  auto p4 = (TPad*)gROOT->FindObject("p4"); 
  p4->cd(); 
  om_c_topo->SetTitle("#Omega_{c}^{0}#rightarrow#Omega^{-}+#pi^{+};DCA_{#Omega}#timesDCA_{#pi} (#mum^{2});Cumulative Counts"); 
  om_c_topo->GetYaxis()->SetRangeUser(0,1.4); 
  om_c_topo->GetXaxis()->SetNdivisions(506);
  om_c_topo->GetYaxis()->SetNdivisions(504);
  om_c_topo->GetYaxis()->SetTitleOffset(1.2); 
  om_c_topo->SetLineColor(kPink+7); 
  om_c_topo->SetMarkerColor(kPink+7); 
  om_c_topo->Draw("hist"); 
  
  ca_c_topo->SetTitle("Background #Omega^{-}+#pi^{+}"); 
  ca_c_topo->SetLineColor(38); 
  ca_c_topo->SetMarkerColor(38); 
  ca_c_topo->Draw("samehist"); 

  p4->BuildLegend(0.173, 0.522, 0.472, 0.793, "#splitline{1.0 < p_{T} (#Omega_{ccc}^{++}) (GeV/#it{c}) < 3.0}{(Topo. Track.)}"); 
  std::cout << "Chomega \n"; 
  std::cout << "To suppress the background down to " << treshold_bkg*100 << "% Topological Reconstruction requires a cut below " << ca_c_topo->GetBinCenter(ca_c_topo->FindFirstBinAbove(treshold_bkg)) << " and retains " << om_c_topo->GetBinContent(ca_c_topo->FindFirstBinAbove(treshold_bkg))*100 << "% of the signal \n"; 
  
  std::cout << "To keep " << treshold_sgn*100 << "% of the signal with topological Reconstruction, select up to " << om_c_topo->GetBinCenter(om_c_topo->FindFirstBinAbove(treshold_sgn)) << " where you keep  " << ca_c_topo->GetBinContent(om_c_topo->FindFirstBinAbove(treshold_sgn))*100 << "% of the background \n"; 
  
  
  auto c8 = c11("8"); 
  auto p8 = (TPad*)gROOT->FindObject("p8"); 
  p8->cd(); 
  om_cc_topo->SetTitle("#Omega_{cc}^{+}#rightarrow#Omega^{0}_{c}+#pi^{+};DCA_{#Omega}#timesDCA_{#pi} (#mum^{2});Cumulative Counts"); 
  om_cc_topo->GetYaxis()->SetRangeUser(0,1.4); 
  om_cc_topo->GetXaxis()->SetNdivisions(506);
  om_cc_topo->GetYaxis()->SetNdivisions(504);
  om_cc_topo->GetYaxis()->SetTitleOffset(1.2); 
  om_cc_topo->SetLineColor(kPink+7); 
  om_cc_topo->SetMarkerColor(kPink+7); 
  om_cc_topo->Draw("hist"); 
  
  ca_cc_topo->SetTitle("Background #Omega^{0}_{c}+#pi^{+}"); 
  ca_cc_topo->SetLineColor(38); 
  ca_cc_topo->SetMarkerColor(38); 
  ca_cc_topo->Draw("samehist"); 

  p8->BuildLegend(0.173, 0.522, 0.472, 0.793, "#splitline{1.0 < p_{T} (#Omega_{ccc}^{++}) (GeV/#it{c}) < 3.0}{(Topo. Track.)}"); 
  std::cout << "CChomega\n"; 
  std::cout << "To suppress the background down to " << treshold_bkg*100 << "% Topological Reconstruction requires a cut below " << ca_cc_topo->GetBinCenter(ca_cc_topo->FindFirstBinAbove(treshold_bkg)) << " and retains " << om_cc_topo->GetBinContent(ca_cc_topo->FindFirstBinAbove(treshold_bkg))*100 << "% of the signal \n"; 

  std::cout << "To keep " << treshold_sgn*100 << "% of the signal with topological Reconstruction, select up to " << om_cc_topo->GetBinCenter(om_cc_topo->FindFirstBinAbove(treshold_sgn)) << " where you keep  " << ca_cc_topo->GetBinContent(om_cc_topo->FindFirstBinAbove(treshold_sgn))*100 << "% of the background \n"; 
  
  auto c9 = c11("9"); 
  auto p9 = (TPad*)gROOT->FindObject("p9"); 
  p9->cd(); 
  om_ccc_topo->SetTitle("#Omega_{ccc}^{++}#rightarrow#Omega^{+}_{cc}+#pi^{+};DCA_{#Omega}#timesDCA_{#pi} (#mum^{2});Cumulative Counts"); 
  om_ccc_topo->GetYaxis()->SetRangeUser(0,1.4); 
  om_ccc_topo->GetXaxis()->SetNdivisions(506);
  om_ccc_topo->GetYaxis()->SetNdivisions(504);
  om_ccc_topo->GetYaxis()->SetTitleOffset(1.2); 
  om_ccc_topo->SetLineColor(kPink+7); 
  om_ccc_topo->SetMarkerColor(kPink+7); 
  om_ccc_topo->Draw("hist"); 
  
  ca_ccc_topo->SetTitle("Background #Omega^{++}_{c}+#pi^{+}"); 
  ca_ccc_topo->SetLineColor(38); 
  ca_ccc_topo->SetMarkerColor(38); 
  ca_ccc_topo->Draw("samehist"); 

  p9->BuildLegend(0.173, 0.522, 0.472, 0.793, "#splitline{1.0 < p_{T} (#Omega_{ccc}^{++}) (GeV/#it{c}) < 3.0}{(Topo. Track.)}"); 
  std::cout << "CCChomega\n"; 
  std::cout << "To suppress the background down to " << treshold_bkg*100 << "% Topological Reconstruction requires a cut below " << ca_ccc_topo->GetBinCenter(ca_ccc_topo->FindFirstBinAbove(treshold_bkg)) << " and retains " << om_ccc_topo->GetBinContent(ca_ccc_topo->FindFirstBinAbove(treshold_bkg))*100 << "% of the signal \n"; 
  
  std::cout << "To keep " << treshold_sgn*100 << "% of the signal with topological Reconstruction, select up to " << om_ccc_topo->GetBinCenter(om_ccc_topo->FindFirstBinAbove(treshold_sgn)) << " where you keep  " << ca_ccc_topo->GetBinContent(om_ccc_topo->FindFirstBinAbove(treshold_sgn))*100 << "% of the background \n"; 


  auto c5 = c11("5"); 
  auto p5 = (TPad*)gROOT->FindObject("p5"); 
  p5->cd(); 
  om_c_stra->SetTitle("#Omega_{c}^{0}#rightarrow#Omega^{-}+#pi^{+};DCA_{#Omega}#timesDCA_{#pi} (#mum^{2});Cumulative Counts"); 
  om_c_stra->GetYaxis()->SetRangeUser(0,1.4); 
  om_c_stra->GetXaxis()->SetNdivisions(506);
  om_c_stra->GetYaxis()->SetNdivisions(504);
  om_c_stra->GetYaxis()->SetTitleOffset(1.2); 
  om_c_stra->SetLineColor(kPink+7); 
  om_c_stra->SetMarkerColor(kPink+7); 
  om_c_stra->Draw("hist"); 
  
  ca_c_stra->SetTitle("Background #Omega^{-}+#pi^{+}"); 
  ca_c_stra->SetLineColor(38); 
  ca_c_stra->SetMarkerColor(38); 
  ca_c_stra->Draw("samehist"); 

  p5->BuildLegend(0.173, 0.522, 0.472, 0.793, "#splitline{1.0 < p_{T} (#Omega_{ccc}^{++}) (GeV/#it{c}) < 3.0}{(Stra. Track. Min. 3 Hits Added)}"); 
  std::cout << "Chomega\n"; 
  std::cout << "To suppress the background down to " << treshold_bkg*100 << "% Strangeness Tracking requires a cut below " << ca_c_stra->GetBinCenter(ca_c_stra->FindFirstBinAbove(treshold_bkg)) << " and retains " << om_c_stra->GetBinContent(ca_c_stra->FindFirstBinAbove(treshold_bkg))*100 << "% of the signal \n"; 
  
  std::cout << "To keep " << treshold_sgn*100 << "% of the signal with strageness tracking, select up to " << om_c_stra->GetBinCenter(om_c_stra->FindFirstBinAbove(treshold_sgn)) << " where you keep  " << ca_c_stra->GetBinContent(om_c_stra->FindFirstBinAbove(treshold_sgn))*100 << "% of the background \n"; 
  

  auto c6 = c11("6"); 
  auto p6 = (TPad*)gROOT->FindObject("p6"); 
  p6->cd(); 
  om_cc_stra->SetTitle("#Omega_{cc}^{+}#rightarrow#Omega^{0}_{c}+#pi^{+};DCA_{#Omega}#timesDCA_{#pi} (#mum^{2});Cumulative Counts"); 
  om_cc_stra->GetYaxis()->SetRangeUser(0,1.4); 
  om_cc_stra->GetXaxis()->SetNdivisions(506);
  om_cc_stra->GetYaxis()->SetNdivisions(504);
  om_cc_stra->GetYaxis()->SetTitleOffset(1.2); 
  om_cc_stra->SetLineColor(kPink+7); 
  om_cc_stra->SetMarkerColor(kPink+7); 
  om_cc_stra->Draw("hist"); 
  
  ca_cc_stra->SetTitle("Background #Omega^{0}_{c}+#pi^{+}"); 
  ca_cc_stra->SetLineColor(38); 
  ca_cc_stra->SetMarkerColor(38); 
  ca_cc_stra->Draw("samehist"); 

  p6->BuildLegend(0.173, 0.522, 0.472, 0.793, "#splitline{1.0 < p_{T} (#Omega_{ccc}^{++}) (GeV/#it{c}) < 3.0}{(Stra. Track. Min. 3 Hits Added)}"); 
  std::cout << "CChomega\n"; 
  std::cout << "To suppress the background down to " << treshold_bkg*100 << "% Strangeness Tracking requires a cut below " << ca_cc_stra->GetBinCenter(ca_cc_stra->FindFirstBinAbove(treshold_bkg)) << " and retains " << om_cc_stra->GetBinContent(ca_cc_stra->FindFirstBinAbove(treshold_bkg))*100 << "% of the signal \n"; 

  std::cout << "To keep " << treshold_sgn*100 << "% of the signal with strageness tracking, select up to " << om_cc_stra->GetBinCenter(om_cc_stra->FindFirstBinAbove(treshold_sgn)) << " where you keep  " << ca_cc_stra->GetBinContent(om_cc_stra->FindFirstBinAbove(treshold_sgn))*100 << "% of the background \n";  

  auto c7 = c11("7"); 
  auto p7 = (TPad*)gROOT->FindObject("p7"); 
  p7->cd(); 
  om_ccc_stra->SetTitle("#Omega_{ccc}^{++}#rightarrow#Omega^{+}_{cc}+#pi^{+};DCA_{#Omega}#timesDCA_{#pi} (#mum^{2});Cumulative Counts"); 
  om_ccc_stra->GetYaxis()->SetRangeUser(0,1.4); 
  om_ccc_stra->GetXaxis()->SetNdivisions(506);
  om_ccc_stra->GetYaxis()->SetNdivisions(504);
  om_ccc_stra->GetYaxis()->SetTitleOffset(1.2); 
  om_ccc_stra->SetLineColor(kPink+7); 
  om_ccc_stra->SetMarkerColor(kPink+7); 
  om_ccc_stra->Draw("hist"); 
  
  ca_ccc_stra->SetTitle("Background #Omega^{++}_{c}+#pi^{+}"); 
  ca_ccc_stra->SetLineColor(38); 
  ca_ccc_stra->SetMarkerColor(38); 
  ca_ccc_stra->Draw("samehist"); 

  p7->BuildLegend(0.173, 0.522, 0.472, 0.793, "#splitline{1.0 < p_{T} (#Omega_{ccc}^{++}) (GeV/#it{c}) < 3.0}{(Stra. Track. Min. 3 Hits Added)}"); 
  std::cout << "CCChomega\n"; 
  std::cout << "To suppress the background down to " << treshold_bkg*100 << "% Strangeness Tracking requires a cut below " << ca_ccc_stra->GetBinCenter(ca_ccc_stra->FindFirstBinAbove(treshold_bkg)) << " and retains " << om_ccc_stra->GetBinContent(ca_ccc_stra->FindFirstBinAbove(treshold_bkg))*100 << "% of the signal \n"; 
 
  std::cout << "To keep " << treshold_sgn*100 << "% of the signal with strageness tracking, select up to " << om_ccc_stra->GetBinCenter(om_ccc_stra->FindFirstBinAbove(treshold_sgn)) << " where you keep  " << ca_ccc_stra->GetBinContent(om_ccc_stra->FindFirstBinAbove(treshold_sgn))*100 << "% of the background \n";  

  c4->SaveAs("DCAProd_Topo_C.pdf"); 
  c8->SaveAs("DCAProd_Topo_CC.pdf"); 
  c9->SaveAs("DCAProd_Topo_CCC.pdf"); 
  c5->SaveAs("DCAProd_StraTrack_C.pdf"); 
  c6->SaveAs("DCAProd_StraTrack_CC.pdf"); 
  c7->SaveAs("DCAProd_StraTrack_CCC.pdf"); 
  
  //Mass Plots
  TH1D* signal = (TH1D*)corre->Get("om_c_mass_topo");
  TH1D* bkg_corre = (TH1D*)corre->Get("ca_mass_topo");

  TH1D* signal_topo = (TH1D*)corre->Get("om_c_mass_topo_cut");
  TH1D* bkg_corre_topo = (TH1D*)corre->Get("ca_mass_topo_cut");
  
  TH1D* signal_stra = (TH1D*)corre->Get("om_c_mass_stra_cut");
  TH1D* bkg_corre_stra = (TH1D*)corre->Get("ca_mass_stra_cut");
  
  TH1D* InvSum = (TH1D*)signal->Clone("InvMasSum");
  InvSum->Add(bkg_corre);
  
  TH1D* InvSum_topo = (TH1D*)signal_topo->Clone("InvMasSum_topo");
  InvSum_topo->Add(bkg_corre_topo);
  
  TH1D* InvSum_stra = (TH1D*)signal_stra->Clone("InvMasSum_stra");
  InvSum_stra->Add(bkg_corre_stra);
    
  InvSum->GetXaxis()->SetTitle("IM(#Omega_{cc}, #pi) (GeV/#it{c}^{ 2})");
  InvSum->GetYaxis()->SetRangeUser(0,30000.);
  InvSum->GetYaxis()->SetTitle("Count Density"); 
  InvSum->GetYaxis()->SetNdivisions(504);
  InvSum->GetYaxis()->SetTitleOffset(1.2);
  InvSum->GetXaxis()->SetNdivisions(506); 
  InvSum->SetTitle("All Candidates");
  
  auto c1 = c11("1");
  auto p1 = (TPad*)gROOT->FindObject("p1"); 
  p1->cd();
  InvSum->DrawCopy("Hist");
  
  signal->SetTitle("#Omega_{ccc} signal"); 
  signal->SetFillColorAlpha(kWhite, 1.0);
  signal->SetLineColor(kWhite); 
  signal->DrawCopy("sameHist");
  
  signal->SetFillColorAlpha(kOrange+7, .2);
  signal->SetLineColor(kOrange+7); 
  signal->DrawCopy("sameHist");

  auto leg1 = p1->BuildLegend(0.18, 0.5, 0.5, 0.8, "1.0 < p_{T} (#Omega_{ccc}) (GeV/#it{c}) < 3.0", "l"); 
  
  InvSum->SetLineColorAlpha(kBlack, 0.2); 
  auto c2 = c11("2");
  auto p2 = (TPad*)gROOT->FindObject("p2"); 
  p2->cd();
  InvSum_topo->SetTitle("All Candidates"); 
  signal_topo->SetTitle("#Omega_{ccc} signal"); 
  
  InvSum->DrawCopy("Hist");
  InvSum_topo->DrawCopy("sameHist");
  
  signal_topo->SetFillColorAlpha(kWhite, 1.0);
  signal_topo->SetLineColor(kWhite); 
  signal_topo->DrawCopy("sameHist");
   
  signal_topo->SetFillColorAlpha(kOrange+7, .2);
  signal_topo->SetLineColor(kOrange+7); 
  signal_topo->DrawCopy("sameHist");

  auto leg2 = p2->BuildLegend(0.18, 0.5, 0.5, 0.8, "1.0 < p_{T} (#Omega_{ccc}) (GeV/#it{c}) < 3.0 Topological Reco", "l");
  
  auto c3 = c11("3");
  auto p3 = (TPad*)gROOT->FindObject("p3"); 
  p3->cd();
  InvSum_stra->SetTitle("All Candidates"); 
  signal_stra->SetTitle("#Omega_{ccc} signal"); 
  
  InvSum->DrawCopy("Hist");
  InvSum_stra->DrawCopy("sameHist");

  signal_stra->SetFillColorAlpha(kWhite, 1.0);
  signal_stra->SetLineColor(kWhite); 
  signal_stra->DrawCopy("sameHist");
  
  signal_stra->SetFillColorAlpha(kOrange+7, .2);
  signal_stra->SetLineColor(kOrange+7); 
  signal_stra->DrawCopy("sameHist");

  auto leg3 = p3->BuildLegend(0.18, 0.5, 0.5, 0.8, "1.0 < p_{T} (#Omega_{ccc}) (GeV/#it{c}) < 3.0 Strangeness Tracking", "l");
  
}
