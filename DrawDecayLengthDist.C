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

void DrawDecayLengthDist() {  
  double tau_xi = 1.6695e-10;//s -> decay length of ~ 4.9 cm 
  double mass_xi = 1.322; //GeV/cc
  double tau_omega = 8.21e-11;//s -> decay length of ~ 2.7 cm 
  double mass_omega = 1.672; //GeV/cc
  double momentum = 1.0; //GeV/c
  
  std::vector<float> layerPos = {0, 0.005, 0.012, 0.025, 0.0375, 0.070, 0.12};
  TF1* f_dcLngth_Xi = new TF1("DecLngthXi", decayLengthDist, 0, 0.5, 3); 
  TF1* f_dcLngth_Om = new TF1("DecLngthOmega", decayLengthDist, 0, 0.5, 3); 

  
  auto c1 = c11("1"); 
  TPad* p = (TPad*)gROOT->FindObject("p1") ; 
  p->cd(); 
  p->SetLogx(); 
  auto d1 = p->DrawFrame(2e-3, 0, 0.2, 1.7); 
  d1->SetTitle(";L (m);N(L)/N_{0} (Surv. Prob.)");
  d1->GetXaxis()->SetNdivisions(504); 
  d1->GetYaxis()->SetNdivisions(506); 
  d1->GetYaxis()->SetTitleOffset(1.4); 
  f_dcLngth_Xi->SetTitle("#Xi^{#minus} (#it{c}#tau = 4.9cm)");
  f_dcLngth_Xi->SetParameters(tau_xi, mass_xi, momentum); 
  f_dcLngth_Xi->SetLineColor(kPink+7); 
  f_dcLngth_Xi->SetLineWidth(2);
  f_dcLngth_Xi->DrawCopy("same");
  
  p->cd(); 
  f_dcLngth_Om->SetTitle("#Omega^{#minus} (#it{c}#tau = 2.5cm)");
  f_dcLngth_Om->SetParameters(tau_omega, mass_omega, momentum); 
  f_dcLngth_Om->SetLineColor(38); 
  f_dcLngth_Om->SetLineWidth(2);
  f_dcLngth_Om->DrawCopy("same");
  
  auto leg = p->BuildLegend( 0.57, 0.59, 0.87, 0.82, TString::Format("p = %.1f GeV/#it{c}", momentum).Data(), "l"); 
  leg->SetFillStyle(0); 

  TArrow layers; 
  layers.SetLineColorAlpha(kBlack, 0.7); 
  //layers.SetLineStyle(2);
  layers.SetLineWidth(2);
  
  TLatex text; 
  text.SetTextFont(43);
  text.SetTextSizePixels(18); 
  text.SetTextColor(1);

  TArrow xiArr; 
  xiArr.SetLineColorAlpha(kPink+7, 0.7);
  xiArr.SetLineWidth(2); 
  xiArr.SetLineStyle(2); 

  TArrow omArr; 
  omArr.SetLineColorAlpha(38, 0.7);
  omArr.SetLineWidth(2); 
  omArr.SetLineStyle(2); 

  for (int iLy = 1; iLy < layerPos.size(); ++iLy) { 
    if (iLy < 4) { 
      layers.DrawArrow(layerPos[iLy], f_dcLngth_Xi->Eval(layerPos[iLy]), layerPos[iLy], f_dcLngth_Xi->Eval(layerPos[iLy])+0.12, 0.01, "<|"); 
      layers.DrawArrow(layerPos[iLy], f_dcLngth_Om->Eval(layerPos[iLy]) -0.12, layerPos[iLy], f_dcLngth_Om->Eval(layerPos[iLy]), 0.01, "|>"); 
      
      text.DrawLatex(layerPos[iLy]*0.84, 
		     (0.93*f_dcLngth_Xi->Eval(layerPos[iLy])+f_dcLngth_Om->Eval(layerPos[iLy]))/2., 
		     TString::Format("Layer %d", iLy-1));
      
      xiArr.DrawArrow(layerPos[iLy], f_dcLngth_Xi->Eval(layerPos[iLy]), 0.2, f_dcLngth_Xi->Eval(layerPos[iLy]), 0.01, "|>"); 
      omArr.DrawArrow(0, f_dcLngth_Om->Eval(layerPos[iLy]), layerPos[iLy], f_dcLngth_Om->Eval(layerPos[iLy]), 0.01, "<|"); 
    }
    
  }

  c1->SaveAs("DecayLengthDist.pdf");

}
