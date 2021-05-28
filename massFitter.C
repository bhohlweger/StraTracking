void Fitme(TH1* hist, double mass) { 
  
  TF1* doubleGauss = new TF1("dGauss", "gaus(0)+gaus(3)", hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax()); 
  
  doubleGauss->SetNpx(10000); 
  std::cout << "For " << hist->GetName() << std::endl; 
  hist->Scale(1./hist->Integral()); 
  
  doubleGauss->SetParameter(0, 0.9*hist->GetMaximum()); //norm 1
  doubleGauss->SetParLimits(0, 1e-6, hist->GetMaximum()); //norm 1

  doubleGauss->SetParameter(1, mass); //mu 1 
  doubleGauss->SetParLimits(1, 0.99*mass, 1.01*mass); //mu 1 
  
  doubleGauss->SetParameter(2, 0.5*hist->GetStdDev(1)); //sig 1 
  doubleGauss->SetParLimits(2, hist->GetStdDev(1)*5e-2, hist->GetStdDev(1)*5e2); //sig 1 

  doubleGauss->SetParameter(3, 0.1*hist->GetMaximum()); //norm 2
  doubleGauss->SetParLimits(3, 1e-6, hist->GetMaximum()); //norm 2

  doubleGauss->SetParameter(4, mass); //mu 2 
  doubleGauss->SetParLimits(4, 0.99*mass, 1.01*mass); //mu 2 

  doubleGauss->SetParameter(5,5*hist->GetStdDev(1)); //sig 2  
  doubleGauss->SetParLimits(5,hist->GetStdDev(1)*5e-2, hist->GetStdDev(1)*5e2); //sig 2 
  hist->Fit(doubleGauss, "MBR+");

  
  TF1* singleGauss = new TF1("sGauss", "gaus(0)", hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax()); 
  singleGauss->SetParameters(doubleGauss->GetParameter(0),doubleGauss->GetParameter(1),doubleGauss->GetParameter(2));
  double norm1 = singleGauss->Integral(hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax()); 
  double sig1 = doubleGauss->GetParameter(2); 
  
  singleGauss->SetParameters(doubleGauss->GetParameter(3),doubleGauss->GetParameter(4),doubleGauss->GetParameter(5));
  double norm2 = singleGauss->Integral(hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax()); 
  double sig2 = doubleGauss->GetParameter(5); 
  double norm = norm1+norm2; 
  
  if (norm < 1e-6) { 
    std::cout << "For " << hist->GetName() << " the fit failed ... Norm1: " << norm1 << " Norm2: " << norm2 << std::endl; 
    norm = 1e266;
  }
  std::cout << "norm1: " << norm1 << std::endl; 
  std::cout << "doubleGauss->GetParameter(2): " << doubleGauss->GetParameter(2) << std::endl; 
  std::cout << "norm2: " << norm2 << std::endl; 
  std::cout << "doubleGauss->GetParameter(5): " << doubleGauss->GetParameter(5) << std::endl; 
  std::cout << "Efffective width: " << (sig1*norm1+sig2*norm2)/norm << std::endl;
    
  delete doubleGauss; 
  return; 

}



void massFitter() { 
  double xiccMass = 3.596; //3.621; 
  double xicMass = 2.468;

  double xiccWindow = 0.08;


  TFile* xicc = TFile::Open("outxiccSelector_xicc.root", "read"); 
  
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
