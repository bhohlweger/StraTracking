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
  TF1* doubleGauss = new TF1("dGauss", "gaus(0)+gaus(3)", hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax()); 
  doubleGauss->SetNpx(10000); 
  std::cout << "For " << hist->GetName() << std::endl; 
  hist->Scale(1./hist->Integral()); 
  
  doubleGauss->SetParameter(0, 0.9*hist->GetMaximum()); //norm 1
  doubleGauss->SetParLimits(0, 1e-6, hist->GetMaximum()); //norm 1

  doubleGauss->SetParameter(1, hist->GetMean()); //mu 1 
  doubleGauss->SetParLimits(1, -1, 1); //mu 1 
  
  doubleGauss->SetParameter(2, 0.5*hist->GetStdDev(1)); //sig 1 
  doubleGauss->SetParLimits(2, hist->GetStdDev(1)*5e-1, hist->GetStdDev(1)*5e1); //sig 1 

  doubleGauss->SetParameter(3, 0.1*hist->GetMaximum()); //norm 2
  doubleGauss->SetParLimits(3, 1e-6, hist->GetMaximum()); //norm 2

  doubleGauss->SetParameter(4, hist->GetMean()); //mu 2 
  doubleGauss->SetParLimits(4, -1, 1); //mu 2 

  doubleGauss->SetParameter(5,5*hist->GetStdDev(1)); //sig 2  
  doubleGauss->SetParLimits(5,hist->GetStdDev(1)*5e-1, hist->GetStdDev(1)*5e1); //sig 2 
  
  hist->Fit(doubleGauss, "MBR+");

  TF1* singleGauss = new TF1("sGauss", "gaus(0)", hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax()); 
  singleGauss->SetParameters(doubleGauss->GetParameter(0),doubleGauss->GetParameter(1),doubleGauss->GetParameter(2));
  double norm1 = singleGauss->Integral(hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax()); 
  double sig1 = doubleGauss->GetParameter(2); 
  singleGauss->SetParameters(doubleGauss->GetParameter(3),doubleGauss->GetParameter(4),doubleGauss->GetParameter(5));
  
  double norm2 = singleGauss->Integral(hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax()); 
  double sig2 = doubleGauss->GetParameter(5); 
  double norm = norm1+norm2; 
    
  /* 
  if (norm < 1e-6) { 
    std::cout << "For " << hist->GetName() << " the fit failed ... Norm1: " << norm1 << " Norm2: " << norm2 << std::endl; 
    norm = 1e266;
  }
  */
  std::cout << "norm1: " << norm1 << std::endl; 
  std::cout << "doubleGauss->GetParameter(2): " << doubleGauss->GetParameter(2) << std::endl; 
  std::cout << "norm2: " << norm2 << std::endl; 
  std::cout << "doubleGauss->GetParameter(5): " << doubleGauss->GetParameter(5) << std::endl; 
  delete doubleGauss; 
  return (sig1*norm1+sig2*norm2)/norm; 
}
