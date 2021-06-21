#include <iostream> 
#include "HarryPlotter.hh"
#include "TF1.h"

void HarryPlotter::StyleBox() { 

  const int NCont = 255;
  //gROOT->ForceStyle();                                                                                                                                                                                                   
  gStyle->Reset("Plain");
  gStyle->SetNumberContours(NCont);

  gStyle->SetTitleBorderSize(0);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(8, 0);
  gStyle->SetCanvasBorderMode(0);

  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.1);
  gStyle->SetPadTopMargin(0.15);
  gStyle->SetPadRightMargin(0.05);

  gStyle->SetFrameLineWidth(1);
  gStyle->SetHistLineWidth(2);
  gStyle->SetFuncWidth(2);
  gStyle->SetLineWidth(2);
  gStyle->SetHatchesLineWidth(1);
  gStyle->SetHatchesSpacing(0.5);
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
  gStyle->SetTitleOffset(0.8, "y");
  gStyle->SetTitleOffset(1.25, "x");

  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  gStyle->SetLegendBorderSize(0);

}

void HarryPlotter::Normalize2DBinByBin(TH2D* hist) {
  for (int iBinX = 1; iBinX <= hist->GetNbinsX(); ++iBinX) { 
    auto projY = hist->ProjectionY("_py", iBinX, iBinX);
    projY->Scale(1./projY->Integral()); 
    for (int iBinY = 1; iBinY <= projY->GetNbinsX(); ++iBinY) { 
      hist->SetBinContent(iBinX, iBinY, projY->GetBinContent(iBinY)); 
    }
    delete projY; 
  }      
};



double HarryPlotter::FitDCA(TH1* hist) { 
  
  auto mydGauss = [] (double *x, double *par) { 
    return par[0]*TMath::Gaus(x[0], par[1], par[2], false) + par[3]*TMath::Gaus(x[0], par[4], par[5], false);
  };
  
  auto myxdGauss = [] (double *x, double *par) { 
    return x[0]*(par[0]*TMath::Gaus(x[0], par[1], par[2], false) + par[3]*TMath::Gaus(x[0], par[4], par[5], false));
  };
  
  auto myxsqdGauss = [] (double *x, double *par) { 
    return (x[0]-par[6])*(x[0]-par[6])*(par[0]*TMath::Gaus(x[0], par[1], par[2], false) + par[3]*TMath::Gaus(x[0], par[4], par[5], false));
  };
  
  
  TF1* doubleGauss = new TF1("dGauss", mydGauss, hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax(), 6); 
  doubleGauss->SetNpx(10000); 
  
  TF1* doubleGauss_x = new TF1("xdGauss", myxdGauss, hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax(), 6); 
  doubleGauss->SetNpx(10000); 
  
  TF1* doubleGauss_xsq = new TF1("xsqdGauss", myxsqdGauss, hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax(), 7); 
  doubleGauss->SetNpx(10000); 
  

  std::cout << "For " << hist->GetName() << std::endl; 
  std::cout << " maximum: " << hist->GetBinContent(hist->FindBin(0)) << std::endl;

  hist->Scale(1./hist->Integral()); 
  
  doubleGauss->SetParameter(0, 0.9*hist->GetBinContent(hist->FindBin(0))); //norm 1
  doubleGauss->SetParLimits(0, 1e-6, 100*hist->GetBinContent(hist->FindBin(0))); //norm 1
			    
  doubleGauss->SetParameter(1, hist->GetMean()); //mu 1 
  doubleGauss->SetParLimits(1, -10, 10); //mu 1 
  
  doubleGauss->SetParameter(2, 0.5*hist->GetStdDev(1)); //sig 1 
  doubleGauss->SetParLimits(2, hist->GetStdDev(1)*1e-2, hist->GetStdDev(1)*1e2); //sig 1 
  
  doubleGauss->SetParameter(3, 0.9*hist->GetBinContent(hist->FindBin(0))); //norm 2
  doubleGauss->SetParLimits(3, 1e-6, 100*hist->GetBinContent(hist->FindBin(0))); //norm 2
			    
  doubleGauss->SetParameter(4, hist->GetMean()); //mu 2 
  doubleGauss->SetParLimits(4, -10, 10); //mu 2 

  doubleGauss->SetParameter(5,1.5*hist->GetStdDev(1)); //sig 2  
  doubleGauss->SetParLimits(5,hist->GetStdDev(1)*1e-2, hist->GetStdDev(1)*1e2); //sig 2 
  
  hist->Fit(doubleGauss, "MBR+");
  
  //std::cout << "doubleGauss->GetParameter(2): " << doubleGauss->GetParameter(2) << std::endl; 
  //std::cout << "doubleGauss->GetParameter(5): " << doubleGauss->GetParameter(5) << std::endl; 
  
  doubleGauss_x->SetParameter(0, doubleGauss->GetParameter(0)); 
  doubleGauss_x->SetParameter(1, doubleGauss->GetParameter(1)); 
  doubleGauss_x->SetParameter(2, doubleGauss->GetParameter(2)); 
  doubleGauss_x->SetParameter(3, doubleGauss->GetParameter(3)); 
  doubleGauss_x->SetParameter(4, doubleGauss->GetParameter(4)); 
  doubleGauss_x->SetParameter(5, doubleGauss->GetParameter(5)); 
  
  double norm = doubleGauss->Integral(hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax()); 
  double mean = doubleGauss_x->Integral(hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax())/norm; 
  
  //std::cout << "mean: " << mean << std::endl; 

  doubleGauss_xsq->SetParameter(0, doubleGauss->GetParameter(0)); 
  doubleGauss_xsq->SetParameter(1, doubleGauss->GetParameter(1)); 
  doubleGauss_xsq->SetParameter(2, doubleGauss->GetParameter(2)); 
  doubleGauss_xsq->SetParameter(3, doubleGauss->GetParameter(3)); 
  doubleGauss_xsq->SetParameter(4, doubleGauss->GetParameter(4)); 
  doubleGauss_xsq->SetParameter(5, doubleGauss->GetParameter(5)); 
  doubleGauss_xsq->SetParameter(6, mean); 
  
  double sigma = doubleGauss_xsq->Integral(hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax())/norm; 
  //std::cout << "sigma: " << sigma << " rms: " << TMath::Sqrt(sigma) << std::endl; 
  
  delete doubleGauss; 
  return TMath::Sqrt(sigma); 
}

float HarryPlotter::YFromMomentum(float ptot, float pT, float mass) { 
  //ptot and pT in GeV/c, mass in GeV/c2
  float pL = TMath::Sqrt(ptot*ptot-pT*pT); 
  float e = TMath::Sqrt(ptot*ptot + mass*mass); 
  return 0.5*TMath::Log((e+pL)/(e-pL));
}; 
