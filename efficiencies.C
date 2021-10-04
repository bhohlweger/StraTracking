Double_t myLevyPt(const Double_t *pt, const Double_t *par)
{
  //Levy Fit Function
  Double_t lMass  = par[3]; //pion Mass
  Double_t ldNdy  = par[0];
  Double_t lTemp = par[1];
  Double_t lPower = par[2];

  Double_t lBigCoef = ((lPower-1)*(lPower-2)) / (lPower*lTemp*(lPower*lTemp+lMass*(lPower-2)));
  Double_t lInPower = 1 + (TMath::Sqrt(pt[0]*pt[0]+lMass*lMass)-lMass) / (lPower*lTemp);

  return ldNdy * pt[0] * lBigCoef * TMath::Power(lInPower,(-1)*lPower);
}

void efficiencies(TString infile) { 
  TFile* xicc = TFile::Open(  TString::Format("%s", infile.Data()), "read"); 
  TF1* fSpectra = new TF1("fSpectra",myLevyPt, 0.0,20., 4);
  TFile* out = TFile::Open("outEff.root", "recreate"); 
  fSpectra ->SetNpx( 1000 );   
  float m = 3.6212;
  fSpectra->SetParameters(0.689877, 0.981593, 8.71805,m);
  
  TH1D* xiccGenPt = (TH1D*)xicc->Get("pTXiccGenerados"); 
  xiccGenPt->Sumw2(); 
  
  TH1D* efficiencyReco = (TH1D*)xicc->Get("EfficiencyNoCutting"); 
  
  TList* inList = xicc->GetListOfKeys(); 
  TIter next(inList); 
  TObject* obj = nullptr; 
  
  while ((obj = next())) {
    TString objName = TString::Format("%s",obj->GetName()); 
    if(!(objName.Contains("pt_vs_eta"))) { 
      continue; 
    }
    TH2D* pT_vs_eta = (TH2D*)xicc->Get(obj->GetName()); 
    TH1D* pT_Identified = (TH1D*)pT_vs_eta->ProjectionX(TString::Format("Eff%s", objName.Data())); 
    pT_Identified->Sumw2();
    pT_Identified->Rebin(4);
    pT_Identified->Divide(xiccGenPt); 
    
    double SumWeight = 0; 
    double WeightedAverage = 0; 
    for (int iBin = 1 ; iBin <= pT_Identified->GetNbinsX(); ++iBin )  {
      double binWeight = fSpectra->Integral(pT_Identified->GetXaxis()->GetBinLowEdge(iBin), pT_Identified->GetXaxis()->GetBinUpEdge(iBin)); 
      WeightedAverage += binWeight*pT_Identified->GetBinContent(iBin); 
      SumWeight += binWeight; 
    }
    std::cout << "For " << objName.Data() << " the average efficiency is " << WeightedAverage/(double)SumWeight << " including xicc and xic BR: " << WeightedAverage/(double)SumWeight*0.05*0.0286  << std::endl;
    out->cd(); 
    pT_Identified->Write(); 
  }
  out->cd();
  fSpectra->Write("Spectrum"); 
}
