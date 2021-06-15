{
  printf("Loading GentleFemtoStyle Box \n");
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
  gStyle->SetHistFillStyle(0); 
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


TCanvas* c11(const char* c) {
  //assumes canvas is used by 70% scaling 
  TString cName = TString::Format("c%s", c);
  TString pName = TString::Format("p%s", c); 
  auto c1 = new TCanvas(cName,cName,0,0,1143,600);
  c1->SetFillStyle(0);
  TPad* p1 = new TPad(pName,pName, 0.15, 0, .85, 1.0);
  p1->SetFillStyle(0);
  p1->SetLeftMargin(c1->GetLeftMargin()/0.7);
  p1->SetRightMargin(c1->GetRightMargin()/0.7); 
  c1->cd();
  p1->Draw();
  p1->cd();
  return c1; 
}

TCanvas* c11Leg(const char* c) {
  //assumes canvas is used by 70% scaling 
  TString cName = TString::Format("c%s", c);
  TString p1Name = TString::Format("p1%s", c);
  TString p2Name = TString::Format("p2%s", c);
  
  auto c1 = new TCanvas(cName,cName,0,0,1143,600);
  c1->SetFillStyle(0);
  TPad* p1 = new TPad(p1Name,p1Name, 0., 0, .6, 1.0);
  p1->SetFillStyle(0);
  p1->SetLeftMargin(c1->GetLeftMargin()/0.6);
  p1->SetRightMargin(c1->GetRightMargin()/0.6);
  TPad* p2 = new TPad(p2Name,p2Name, 0.6, 0, 1., 1.0);
  p2->SetFillStyle(0);
  p2->SetLeftMargin(c1->GetLeftMargin()/0.4);
  p2->SetRightMargin(c1->GetRightMargin()/0.4);
  
  c1->cd();
  p1->Draw();
  p2->Draw();
  return c1; 
}

TCanvas* c12(const char* c) {
  //assumes canvas is used by 100% scaling 
  TString cName = TString::Format("c%s", c);
  TString pName = TString::Format("p%s", c);
  auto c1 = new TCanvas(cName,cName,0,0,1143,600);
  TPad* p1 = new TPad(pName,pName, 0.25, 0, 0.85, 1.0);
  p1->SetLeftMargin(c1->GetLeftMargin()/0.6);
  p1->SetRightMargin(c1->GetRightMargin()/0.6); 
  c1->cd();
  p1->Draw();
  return c1; 
}

TCanvas* c21(const char* c) {
  //assumes canvas is used by 100% scaling 
  TString cName = TString::Format("c%s", c);
  TString p1Name = TString::Format("p1%s", c);
  TString p2Name = TString::Format("p2%s", c); 
  auto c1 = new TCanvas(cName,cName,0,0,1143,600);
  c1->SetFillStyle(0);
  TPad* p1 = new TPad(p1Name,p1Name, 0, 0, 0.5, 1.0);
  p1->SetFillStyle(0);
  p1->SetLeftMargin(c1->GetLeftMargin()/0.5);
  p1->SetRightMargin(c1->GetRightMargin()/0.5); 
  TPad* p2 = new TPad(p2Name,p2Name, 0.5, 0, 1.0, 1.0);
  p2->SetFillStyle(0);
  p2->SetLeftMargin(c1->GetLeftMargin()/0.5);
  p2->SetRightMargin(c1->GetRightMargin()/0.5); 
  c1->cd();
  p1->Draw();
  p2->Draw(); 
  return c1; 
}

TCanvas* c31(const char* c) {
  //assumes canvas is used by 100% scaling 
  TString cName = TString::Format("c%s", c);
  TString p1Name = TString::Format("p1%s", c);
  TString p2Name = TString::Format("p2%s", c);
  TString p3Name = TString::Format("p3%s", c); 
  auto c1 = new TCanvas(cName,cName,0,0,1143,600);
  c1->SetFillStyle(0);
  TPad* p1 = new TPad(p1Name,p1Name, 0., 0, .6, 1.0);
  p1->SetFillStyle(0);
  p1->SetLeftMargin(c1->GetLeftMargin()/0.6);
  p1->SetRightMargin(c1->GetRightMargin()/0.6); 

  TPad* p2 = new TPad(p2Name,p2Name, 0.6, 0.5, 1.0, 1.0);
  p2->SetFillStyle(0);
  p2->SetLeftMargin(c1->GetLeftMargin()/0.4);
  p2->SetRightMargin(c1->GetRightMargin()/0.4);

  p2->SetTopMargin(c1->GetTopMargin()/0.5);
  p2->SetBottomMargin(c1->GetBottomMargin()/0.5);
  

  TPad* p3 = new TPad(p3Name,p3Name, 0.6, 0., 1., 0.5);
  p3->SetFillStyle(0);
  p3->SetLeftMargin(c1->GetLeftMargin()/0.4);
  p3->SetRightMargin(c1->GetRightMargin()/0.4);
  
  p3->SetTopMargin(c1->GetTopMargin()/0.5);
  p3->SetBottomMargin(c1->GetBottomMargin()/0.6);
  
  c1->cd();
  p1->Draw();
  p2->Draw();
  p3->Draw(); 
  return c1; 
}


TCanvas* c21_SharedAxis(const char* c) {
  //assumes canvas is used by 100% scaling 
  TString cName = TString::Format("c%s", c);
  TString p1Name = TString::Format("p1%s", c);
  TString p2Name = TString::Format("p2%s", c); 
  auto c1 = new TCanvas(cName,cName,0,0,1143,600);
  c1->SetFillStyle(0);
  TPad* p1 = new TPad(p1Name,p1Name, 0, 0, 0.5, 1.0);
  p1->SetFillStyle(0);
  p1->SetLeftMargin(c1->GetLeftMargin()/0.5+c1->GetRightMargin()/0.5);
  p1->SetRightMargin(0); 
  TPad* p2 = new TPad(p2Name,p2Name, 0.5, 0, 1.0, 1.0);
  p2->SetFillStyle(0);
  p2->SetLeftMargin(0);
  p2->SetRightMargin(c1->GetLeftMargin()/0.5+c1->GetRightMargin()/0.5); 
  c1->cd();
  p1->Draw();
  p2->Draw(); 
  return c1; 
}


void ScaleToPad(TH1* h, float scaleX, float scaleY) {
  h->GetXaxis()->SetTitleOffset(gStyle->GetTitleOffset("x")/scaleX);
  h->GetXaxis()->SetLabelOffset(gStyle->GetLabelOffset("x")/scaleX);

  h->GetYaxis()->SetTitleOffset(gStyle->GetTitleOffset("y")/scaleY);
  h->GetYaxis()->SetLabelOffset(gStyle->GetLabelOffset("y")/scaleY);
  
  return; 
} 

void ScaleToPad(TH2* h, float scaleX, float scaleY) {
  h->GetXaxis()->SetTitleOffset(gStyle->GetTitleOffset("x")/scaleX);
  h->GetXaxis()->SetLabelOffset(gStyle->GetLabelOffset("x")/scaleX);

  h->GetYaxis()->SetTitleOffset(gStyle->GetTitleOffset("y")/scaleY);
  h->GetYaxis()->SetLabelOffset(gStyle->GetLabelOffset("y")/scaleY);
  
  return; 
} 


TLatex* GenTex() {
  TLatex *text = new TLatex();
  text->SetTextFont(43);
  text->SetTextSizePixels(28); 
  text->SetNDC();
  text->SetTextColor(1);
  return text; 
}


TBrowser* tb() {
  return new TBrowser();
} 
