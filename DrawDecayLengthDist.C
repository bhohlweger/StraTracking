/* 

The math ... 

dN/dt = - lmb*N = - N / tau -> Decay Law ... (Eq. 1)

Now do the variable substitution 
-> v = beta * c = L / t 
-> t = L / beta * c 
-> dt / dL = 1 / beta * c 
-> dt = dL / beta * c 

Insert into (Eq. 1)

dN/N = - dL/(beta*c*tau)

Integration on each side gives 

ln(N) = - L/(beta*c*tau) + C ( < -- C related to the starting population)

-> N(L) = exp(-L/(beta*c*tau))*exp(C) ( exp(C) is the Starting Population N_0)

-> N(L)/N_0 = exp(-L/(beta*c*tau)) 

beta comes from momentum: 

p = m * v = m_0 * gamma * beta * c 

in the following: m = m_0, gamma = 1/sqrt(1-beta*beta), solve for beta 

-> beta = 1/sqrt(1+(c*m/p)*(c*m/p))

-> 1/beta = sqrt(1+(c*m/p)*(c*m/p))

*/  


double decayLengthDist(double* x, double* par) { 
  double L = x[0]; // in m 
  double tau = par[0]; // in s 
  double mass = par[1]; // in GeV/cc (ignore one c and tf don't multiply by c)
  double p = par[2]; // in GeV/c
  double c = 299792458; // m/s 
  return TMath::Exp(-TMath::Sqrt(1+(mass/p)*(mass/p))/(c*tau)*L); 
}

/////////////////////////
//helper functions ... //
/////////////////////////
 
TCanvas* cc21(const char* c) {
  //assumes canvas is used by 100% scaling                                                                                                                                                                                                    
  TString cName = TString::Format("c%s", c);
  TString p1Name = TString::Format("p1%s", c);
  TString p2Name = TString::Format("p2%s", c);
  auto c1 = new TCanvas(cName,cName,0,0,1143,600);
  TPad* p1 = new TPad(p1Name,p1Name, 0, 0, 0.5, 1.0);
  p1->SetLeftMargin(c1->GetLeftMargin()/0.5);
  p1->SetRightMargin(c1->GetRightMargin()/0.5);
  TPad* p2 = new TPad(p2Name,p2Name, 0.5, 0, 1.0, 1.0);
  p2->SetLeftMargin(c1->GetLeftMargin()/0.5);
  p2->SetRightMargin(c1->GetRightMargin()/0.5);
  c1->cd();
  p1->Draw();
  p2->Draw();
  return c1;
}

void StyleBox() {
  const int NCont = 255;
  //gROOT->ForceStyle();                                                                                                                                                                                                                      
  gStyle->Reset("Plain");
  gStyle->SetNumberContours(NCont);

  gStyle->SetTitleBorderSize(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPalette(8, 0);
  gStyle->SetCanvasBorderMode(0);

  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.1);
  gStyle->SetPadTopMargin(0.15);
  gStyle->SetPadRightMargin(0.05);

  gStyle->SetFrameLineWidth(1);
  gStyle->SetHistLineWidth(1);
  gStyle->SetFuncWidth(2);
  gStyle->SetLineWidth(2);
  gStyle->SetHatchesLineWidth(1);
  gStyle->SetHatchesSpacing(0.5);

  //gStyle->SetErrorX(0.005);                                                                                                                                                                                                                 

  gStyle->SetCanvasColor(0);
  gStyle->SetPadColor(0);
  gStyle->SetHistLineColor(kBlack);
  gStyle->SetFuncColor(kGreen);
  gStyle->SetLabelColor(kBlack, "xyz");
  gStyle->SetPalette(kCividis);

  gStyle->SetTextFont(43);
  gStyle->SetLabelFont(43, "xyz");
  gStyle->SetTitleFont(43, "xyz");
  gStyle->SetLegendFont(43);

  gStyle->SetTextSizePixels(28);
  gStyle->SetLabelSize(28, "xyz");
  gStyle->SetTitleSize(28, "xyz");
  gStyle->SetLegendTextSize(28);

  gStyle->SetLabelOffset(0.01, "xy");
  gStyle->SetTitleOffset(1.2, "y");
  gStyle->SetTitleOffset(1.25, "x");

  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  gStyle->SetLegendBorderSize(0);


}


/////////////////////////
//end helper functions //
/////////////////////////


void DrawDecayLengthDist() { 
  StyleBox(); 
  double tau_xi = 1.6695e-10;//s -> decay length of ~ 4.9 cm 
  double mass_xi = 1.322; //GeV/cc
  double tau_omega = 8.21e-11;//s -> decay length of ~ 2.7 cm 
  double mass_omega = 1.672; //GeV/cc
  double momentum = 1.0; //GeV/c
  
  std::vector<float> layerPos = {0, 0.005, 0.012, 0.025, 0.0375, 0.070, 0.12, 0.20};
  TF1* f_dcLngth = new TF1("DecLngth", decayLengthDist, 0, 0.5, 3); 
  //draw the xi 

  
  auto c1 = cc21("1"); 
  TPad* p1 = (TPad*)gROOT->FindObject("p11") ; 
  TPad* p2 = (TPad*)gROOT->FindObject("p21") ; 
  p1->cd(); 
  auto d1 = p1->DrawFrame(0, 0, 0.1, 1.2); 
  d1->SetTitle(";L (m);N(L)/N_{0} (Surv. Prob)");
  d1->GetXaxis()->SetNdivisions(504); 
  d1->GetYaxis()->SetNdivisions(506); 
  d1->GetYaxis()->SetTitleOffset(1.4); 
  f_dcLngth->SetParameters(tau_xi, mass_xi, momentum); 
  f_dcLngth->SetLineColor(kPink+7); 
  f_dcLngth->SetLineWidth(3);
  f_dcLngth->DrawCopy("same");
  
  TH1D* ProbHits_Xi = new TH1D("ProbHitsXi", "#Xi (#it{c}#tau = 4.9 cm)", layerPos.size()-1, 0, layerPos.size()); 
  //Evaluate this at Layerpos to get the probability of a candidates surviving ... 
  for (int iLy = 1; iLy < layerPos.size(); ++iLy) { 
    ProbHits_Xi->SetBinContent(iLy, f_dcLngth->Eval(layerPos[iLy])); 
    ProbHits_Xi->GetXaxis()->SetBinLabel(iLy, TString::Format("Layer %d", iLy).Data()); 
  }
  p2->cd(); 
  ProbHits_Xi->SetLineColor(kPink+7); 
  ProbHits_Xi->SetMarkerColor(kPink+7); 
  ProbHits_Xi->SetLineWidth(3);
  ProbHits_Xi->GetXaxis()->SetNdivisions(504); 
  ProbHits_Xi->GetXaxis()->SetTitle("");   
  ProbHits_Xi->GetYaxis()->SetNdivisions(506); 
  ProbHits_Xi->GetYaxis()->SetTitle("Surv. Prob at each Layer");   
  ProbHits_Xi->GetYaxis()->SetTitleOffset(1.4); 
  ProbHits_Xi->DrawCopy("");
  
  p1->cd(); 
  f_dcLngth->SetParameters(tau_omega, mass_omega, momentum); 
  f_dcLngth->SetLineColor(38); 
  f_dcLngth->DrawCopy("same");
  
  TH1D* ProbHits_Omega = new TH1D("ProbHitsOmega", "#Omega (#it{c}#tau = 2.5 cm)", layerPos.size()-1, 0, layerPos.size()); 
  //Evaluate this at Layerpos to get the probability of a candidates surviving ... 
  for (int iLy = 1; iLy < layerPos.size(); ++iLy) { 
    ProbHits_Omega->SetBinContent(iLy, f_dcLngth->Eval(layerPos[iLy])); 
    ProbHits_Omega->GetXaxis()->SetBinLabel(iLy, TString::Format("Layer %d", iLy).Data()); 
  }
  ProbHits_Omega->SetLineColor(38); 
  ProbHits_Omega->SetMarkerColor(38); 
  ProbHits_Omega->SetLineWidth(3);
  p2->cd(); 
  ProbHits_Omega->DrawCopy("same");
  
  p2->BuildLegend( 0.453,0.636,0.753,0.846, TString::Format("p = %.1f GeV/#it{c}", momentum).Data()); 
  
  c1->SaveAs("DecayLengthDist.pdf");

}
