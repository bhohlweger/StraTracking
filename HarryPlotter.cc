#include "HarryPlotter.hh"

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



