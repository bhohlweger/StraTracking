#ifndef HARRYPLOTTER_HH
#define HARRYPLOTTER_HH
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TString.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TFile.h" 
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RResultPtr.hxx"

class HarryPlotter { 
public: 
  HarryPlotter() {};
  static void StyleBox(); 
  template <typename T> static void CheckAndStore(TFile *out, T key); 
  template <typename T> static void NormToMaxCheckAndStore(TFile *out, T key); 
  template <typename T> static void CumulateAndStore(TFile *out, T one); 	   
  template <typename T, typename F> static void PlotAndStore(TString name, TFile* out, T one, T two, T three, T four);
  template <typename T,typename F> static void CumulatePlotAndStore(TString name, TFile* out, T one, T two, T three, T four); 
  template <typename T> static void AverageBackground(TFile* out, T signal, T background); 
  
  static void Normalize2DBinByBin(TH2D* hist); 
  
  static bool TopoCuts_om(float pT) {return ((1. < pT)&&(pT < 6.0));}; 
  static bool StraCuts_om(float pT, int iAdded) {return ((1. < pT)&&(pT < 6.0)&&(iAdded>=2));}; 
  static bool TopoCuts_om_c(float pT, float decayRad) {return ((1.0 < pT)&&(pT < 3.0)&&(decayRad>1e-20));}; 
  static bool StraCuts_om_c(float pT, int iAdded, float decayRad) {return ((1.0 < pT)&&(pT < 3.0)&&(iAdded>=2)&&(decayRad>1e-20));}; 
  static bool TopoCuts_ca_c(float pT, float decayRad, bool isOmegac) {return ((1.0 < pT)&&(pT < 3.0)&&(decayRad>1e-20)&&(!isOmegac));}; 
  static bool StraCuts_ca_c(float pT, int iAdded, float decayRad, bool isOmegac) {return ((1.0 < pT)&&(pT < 3.0)&&(iAdded>=2)&&(decayRad>1e-20)&&(!isOmegac));}; 
  static bool AddedHitsCut(int iAdded) { return (iAdded>=2);}; 
  
  static double FitDCA(TH1* hist); 
  
  static bool StraCut_dca(float dcaprod){return (dcaprod <  -41);}; 
  static bool TopoCut_dca(float dcaprod){return (dcaprod <  -509);}; 
  
  static bool StraCut_c_dca(float dcaprod){return (dcaprod <  353);}; 
  static bool TopoCut_c_dca(float dcaprod){return (dcaprod <  3526);}; 

  static bool StraCut_cc_dca(float dcaprod){return (dcaprod <  126);}; 
  static bool TopoCut_cc_dca(float dcaprod){return (dcaprod <  2873);}; 

  static bool StraCut_ccc_dca(float dcaprod){return (dcaprod <  100);}; 
  static bool TopoCut_ccc_dca(float dcaprod){return (dcaprod <  1940);}; 

  static float Distance(float x1, float y1, float z1) {return (float)TMath::Sqrt(x1*x1+y1*y1+z1*z1);};  
  static float Distance_xy(float x1, float y1, float x2, float y2) {
    double dx = x1-x2; 
    double dy = y1-y2; 
    return (float)(dx>0?1:-1)*TMath::Sqrt(dx*dx+dy*dy);
  };  
  static float Distance_z(float z1, float z2) {return (float)(z1-z2);};   
  static float YFromMomentum(float ptot, float pT, float mass, float eta); 
  static float RelErr(float generated, float measured) {return (float)(generated-measured)/generated;}; 
  static const char* FilePath() { return "/localstore/alice/hohlweger/analysis/StrangnessTracking/210319";}; 
  static const char* FilePathPythia() { return "/localstore/alice/hohlweger/analysis/StrangnessTracking/210406/omegacccpythia_new.root";}; 
  
  static std::vector<float> Getptbins() { std::vector<float> ptbins = {0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 20.0}; return ptbins;}; 
  static std::vector<float> GetposLayers() { std::vector<float> layerPos = {0, 0.5, 1.2, 2.5, 3.75, 7.0, 12, 20, 30, 45, 60, 80, 100}; ; return layerPos;}; 
  
}; 
#include "HarryPlotter.tpp"
#endif
