double mydGauss (double *x, double *par) { 
  return par[0]*TMath::Gaus(x[0], par[1], par[2], true) + par[3]*TMath::Gaus(x[0], par[4], par[5], true);
}

double myxdGauss (double *x, double *par) { 
  return x[0]*(par[0]*TMath::Gaus(x[0], par[1], par[2], true) + par[3]*TMath::Gaus(x[0], par[4], par[5], true));
}

double myxsqdGauss (double *x, double *par) { 
  return (x[0]-par[6])*(x[0]-par[6])*(par[0]*TMath::Gaus(x[0], par[1], par[2], true) + par[3]*TMath::Gaus(x[0], par[4], par[5], true));
}

void Fitme(TH1* hist, double mass) { 
  
  TF1* doubleGauss = new TF1("dGauss", mydGauss, hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax(), 6); 
  doubleGauss->SetNpx(10000); 
  
  TF1* doubleGauss_x = new TF1("xdGauss", myxdGauss, hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax(), 6); 
  doubleGauss->SetNpx(10000); 
  
  TF1* doubleGauss_xsq = new TF1("xsqdGauss", myxsqdGauss, hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax(), 7); 
  doubleGauss->SetNpx(10000); 
  

  std::cout << "For " << hist->GetName() << std::endl; 
  
  hist->Scale(1./hist->Integral()); 
  
  doubleGauss->SetParameter(0, 0.9*hist->GetMaximum()); //norm 1
  doubleGauss->SetParLimits(0, 1e-6, hist->GetMaximum()); //norm 1

  doubleGauss->SetParameter(1, mass); //mu 1 
  doubleGauss->SetParLimits(1, 0.98*mass, 1.02*mass); //mu 1 
  
  doubleGauss->SetParameter(2, 0.5*hist->GetStdDev(1)); //sig 1 
  doubleGauss->SetParLimits(2, hist->GetStdDev(1)*1e-1, hist->GetStdDev(1)*10); //sig 1 

  doubleGauss->SetParameter(3, 0.1*hist->GetMaximum()); //norm 2
  doubleGauss->SetParLimits(3, 1e-6, hist->GetMaximum()); //norm 2

  doubleGauss->SetParameter(4, mass); //mu 2 
  doubleGauss->SetParLimits(4, 0.98*mass, 1.02*mass); //mu 2 

  doubleGauss->SetParameter(5,1.5*hist->GetStdDev(1)); //sig 2  
  doubleGauss->SetParLimits(5,hist->GetStdDev(1)*1e-1, hist->GetStdDev(1)*10); //sig 2 
  
  hist->Fit(doubleGauss, "MBR+");
  
  std::cout << "doubleGauss->GetParameter(2): " << doubleGauss->GetParameter(2) << std::endl; 
  std::cout << "doubleGauss->GetParameter(5): " << doubleGauss->GetParameter(5) << std::endl; 
  
  doubleGauss_x->SetParameter(0, doubleGauss->GetParameter(0)); 
  doubleGauss_x->SetParameter(1, doubleGauss->GetParameter(1)); 
  doubleGauss_x->SetParameter(2, doubleGauss->GetParameter(2)); 
  doubleGauss_x->SetParameter(3, doubleGauss->GetParameter(3)); 
  doubleGauss_x->SetParameter(4, doubleGauss->GetParameter(4)); 
  doubleGauss_x->SetParameter(5, doubleGauss->GetParameter(5)); 
  
  double norm = doubleGauss->Integral(hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax()); 
  double mean = doubleGauss_x->Integral(hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax())/norm; 
  
  std::cout << "mean: " << mean << std::endl; 

  doubleGauss_xsq->SetParameter(0, doubleGauss->GetParameter(0)); 
  doubleGauss_xsq->SetParameter(1, doubleGauss->GetParameter(1)); 
  doubleGauss_xsq->SetParameter(2, doubleGauss->GetParameter(2)); 
  doubleGauss_xsq->SetParameter(3, doubleGauss->GetParameter(3)); 
  doubleGauss_xsq->SetParameter(4, doubleGauss->GetParameter(4)); 
  doubleGauss_xsq->SetParameter(5, doubleGauss->GetParameter(5)); 
  doubleGauss_xsq->SetParameter(6, mean); 
  
  double sigma = doubleGauss_xsq->Integral(hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax())/norm; 
  std::cout << "sigma: " << sigma << " rms: " << TMath::Sqrt(sigma) << std::endl; 
  
  delete doubleGauss; 
  return; 

}



void massFitter() { 
  double xiccMass = 3.596; //3.621; 
  double xicMass = 2.468;

  double xiccWindow = 0.08;


  TFile* xicc = TFile::Open("outxiccSelector_xicc-fulltracker.root", "read"); 
  
  TH1D* lambda_mass = (TH1D*)xicc->Get("df_xi_im_lmb_mass"); 
  TH1D* xi_mass = (TH1D*)xicc->Get("df_xi_im_xi_mass"); 
  TH1D* xic_mass = (TH1D*)xicc->Get("df_xi_c_im_xi_c_mass_stra"); 
  TH1D* xicc_mass = (TH1D*)xicc->Get("df_xi_cc_im_xi_cc_mass_stra_c1"); 
  
  TFile* output = TFile::Open("outMassFit.root", "recreate"); 
  output->cd(); 
  
  Fitme(lambda_mass, 1.116); 
  Fitme(xi_mass, 1.322); 
  Fitme(xic_mass, xicMass); 
  Fitme(xicc_mass, xiccMass); 
  
  lambda_mass->Write();
  xi_mass->Write();
  xic_mass->Write(); 
  xicc_mass->Write(); 
  
  output->Close();
  xicc->Close(); 
}
