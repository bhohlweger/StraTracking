std::vector<float> ptbins = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};

void DCAOmega() { 
    
  TFile* input = TFile::Open("outResolutionsDCA_xi.root", "read"); 
  TH1D* dca_xy_topo = (TH1D*)input->Get("dca_xy_topo"); 
  TH1D* dca_z_topo = (TH1D*)input->Get("dca_z_topo"); 
  
  TH2F* dca_xy_stra_vs_ah = (TH2F*)input->Get("dca_xy_vs_AddedHits_stra"); 
  TH2F* dca_z_stra_vs_ah = (TH2F*)input->Get("dca_z_vs_AddedHits_stra"); 
  
  TH1D* dca_xy_stra = (TH1D*)dca_xy_stra_vs_ah->ProjectionY("dca_xy_stra", 2, 6);
  TH1D* dca_z_stra = (TH1D*)dca_z_stra_vs_ah->ProjectionY("dca_z_stra", 2, 6); 
  
  TH1D* width_dca_xy = (TH1D*)input->Get("width_dca_xy_integrated"); 
  TH1D* width_dca_z = (TH1D*)input->Get("width_dca_z_integrated"); 

  TH1D* width_dca_xy_1_2 = (TH1D*)input->Get("width_dca_xy_pT_1.0_2.0"); 
  TH1D* width_dca_xy_2_3 = (TH1D*)input->Get("width_dca_xy_pT_2.0_3.0"); 
  TH1D* width_dca_xy_3_4 = (TH1D*)input->Get("width_dca_xy_pT_3.0_4.0"); 
  TH1D* width_dca_xy_4_5 = (TH1D*)input->Get("width_dca_xy_pT_4.0_5.0"); 
  TH1D* width_dca_xy_5_6 = (TH1D*)input->Get("width_dca_xy_pT_5.0_6.0"); 
  
  
  
  dca_xy_topo->Scale(1./dca_xy_topo->GetMaximum()); 
  dca_z_topo->Scale(1./dca_z_topo->GetMaximum()); 
  
  dca_xy_stra->Scale(1./dca_xy_stra->GetMaximum()); 
  dca_z_stra->Scale(1./dca_z_stra->GetMaximum()); 
  
  auto c1 = c21_SharedAxis("1"); 
  auto p1_1 = (TPad*)gROOT->FindObject("p11"); 
  p1_1->cd(); 
  ScaleToPad(dca_xy_topo, 1.0, .5); 
  dca_xy_topo->SetTitle("Topological Tracking; DCA_{xy} (#mum); Density w.r.t. Maximum"); 
  dca_xy_topo->GetXaxis()->SetRangeUser(-199,199); 
  dca_xy_topo->GetXaxis()->SetNdivisions(506);

  dca_xy_topo->GetYaxis()->SetRangeUser(0,1.8); 
  dca_xy_topo->GetYaxis()->SetNdivisions(504);
  dca_xy_topo->SetLineColor(kPink+7); 
  dca_xy_topo->SetMarkerColor(kPink+7); 
  dca_xy_topo->Draw("hist"); 
  
  dca_xy_stra->SetTitle("Strangeness Tracking"); 
  dca_xy_stra->SetLineColor(38); 
  dca_xy_stra->SetMarkerColor(38); 
  dca_xy_stra->Draw("samehist"); 

  p1_1->BuildLegend(0.365, 0.541, 0.826, 0.812, "1.0 < p_{T} (#Xi^{-}) (GeV/#it{c}) < 6.0"); 
  
  
  auto p1_2 = (TPad*)gROOT->FindObject("p21"); 
  p1_2->cd(); 
  ScaleToPad(dca_z_topo, 1.0, .5); 
  dca_z_topo->SetTitle("Topological Tracking; DCA_{z} (#mum); Density w.r.t. Maximum");

  dca_z_topo->GetXaxis()->SetRangeUser(-199,199); 
  dca_z_topo->GetXaxis()->SetNdivisions(506);

  dca_z_topo->GetYaxis()->SetRangeUser(0,1.8); 
  dca_z_topo->GetYaxis()->SetNdivisions(504);
  //dca_z_topo->GetYaxis()->SetTitleOffset(1.2); 
  dca_z_topo->SetLineColor(kPink+7); 
  dca_z_topo->SetMarkerColor(kPink+7); 
  dca_z_topo->Draw("hist"); 
  
  dca_z_stra->SetTitle("Strangeness Tracking"); 
  dca_z_stra->SetLineColor(38); 
  dca_z_stra->SetMarkerColor(38); 
  dca_z_stra->Draw("samehist"); 

  //p1_2->BuildLegend(0.173, 0.572, 0.472, 0.843, "1.0 < p_{T} (#Xi^{-}) (GeV/#it{c}) < 6.0"); 


  auto c2 = c11("2"); 
  auto p2_1 = (TPad*)gROOT->FindObject("p2"); 
  p2_1->cd(); 
  p2_1->SetRightMargin(0.15);
  p2_1->SetBottomMargin(0.18);
  
  //ScaleToPad(width_dca_xy, 1.0, .5); 
  width_dca_xy->SetTitle("DCA_{xy}");
  width_dca_xy->GetYaxis()->SetTitle("#sigma_{DCA} (#mum)"); 
  width_dca_xy->GetYaxis()->SetRangeUser(0,width_dca_xy->GetMaximum()*1.2); 
  width_dca_xy->GetXaxis()->SetNdivisions(506);
  width_dca_xy->GetXaxis()->SetRange(1, 7); 
  width_dca_xy->GetYaxis()->SetNdivisions(504);
  width_dca_xy->GetYaxis()->SetTitleOffset(1.6); 
  width_dca_xy->SetLineColor(kOrange+7); 
  width_dca_xy->SetMarkerColor(kOrange+7); 
  for (int ib = 2; ib <=width_dca_xy->GetNbinsX(); ++ib) { 
    width_dca_xy->GetXaxis()->SetBinLabel(ib, TString::Format("%d add. Hits", ib-1)); 
  }

  width_dca_xy->Draw("hist"); 
  
  width_dca_z->SetTitle("DCA_{z}"); 
  width_dca_z->SetLineColor(kGreen+3); 
  width_dca_z->SetMarkerColor(kGreen+3); 
  width_dca_z->SetLineStyle(9); 
  width_dca_z->Draw("samehist"); 

  p2_1->BuildLegend(0.45, 0.550, 0.75, 0.82, "1.0 < p_{T} (#Xi^{-}) (GeV/#it{c}) < 6.0"); 
  

  auto c3 = c11("3"); 
  auto p3_1 = (TPad*)gROOT->FindObject("p3"); 
  p3_1->cd(); 
  p3_1->SetRightMargin(0.15);
  p3_1->SetBottomMargin(0.18);
  
  //ScaleToPad(width_dca_xy, 1.0, .5); 
  width_dca_xy_1_2->SetTitle("1.0 < p_{T} < 2.0");
  width_dca_xy_1_2->GetYaxis()->SetTitle("#sigma_{DCA} (#mum)"); 
  width_dca_xy_1_2->GetYaxis()->SetRangeUser(0,width_dca_xy_1_2->GetMaximum()*1.2); 
  width_dca_xy_1_2->GetXaxis()->SetNdivisions(506);
  width_dca_xy_1_2->GetXaxis()->SetRange(1, 7); 
  width_dca_xy_1_2->GetYaxis()->SetNdivisions(504);
  width_dca_xy_1_2->GetYaxis()->SetTitleOffset(1.6); 
  width_dca_xy_1_2->SetLineColor(kOrange+7); 
  width_dca_xy_1_2->SetMarkerColor(kOrange+7); 
  for (int ib = 2; ib <=width_dca_xy->GetNbinsX(); ++ib) { 
    width_dca_xy_1_2->GetXaxis()->SetBinLabel(ib, TString::Format("%d add. Hits", ib-1)); 
  }

  width_dca_xy_1_2->Draw("hist"); 
  

  width_dca_xy_2_3->SetTitle("2.0 < p_{T} < 3.0");
  width_dca_xy_2_3->SetLineColor(38); 
  width_dca_xy_2_3->SetMarkerColor(38); 
  width_dca_xy_2_3->Draw("samehist"); 
  
  width_dca_xy_3_4->SetTitle("3.0 < p_{T} < 4.0");
  width_dca_xy_3_4->SetLineColor(kPink+7); 
  width_dca_xy_3_4->SetMarkerColor(kPink+7);   
  width_dca_xy_3_4->Draw("samehist"); 
  
  width_dca_xy_4_5->SetTitle("4.0 < p_{T} < 5.0");
  width_dca_xy_4_5->SetLineColor(kGreen+3); 
  width_dca_xy_4_5->SetMarkerColor(kGreen+3); 
  width_dca_xy_4_5->Draw("samehist"); 
  
  width_dca_xy_5_6->SetTitle("5.0 < p_{T} < 6.0");
  width_dca_xy_5_6->SetLineColor(kBlack); 
  width_dca_xy_5_6->SetMarkerColor(kBlack); 
  width_dca_xy_5_6->Draw("samehist"); 

  p3_1->BuildLegend(0.35, 0.550, 0.65, 0.82, "#Xi^{-} p_{T} (GeV/#it{c})"); 
}
