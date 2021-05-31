#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RResultPtr.hxx"

#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TString.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TFile.h" 
#include "HarryPlotter.hh"

std::vector<float> ptbins = {1.,6.};//HarryPlotter::Getptbins(); 
std::vector<float> layerPos = HarryPlotter::GetposLayers();

/* 
   workflow: 
   1) take the xis for kind of granted -> QA cuts
   -> additional hits + pT cut 
   -> transverse Radius (difference between truth and reco, causality) + Daughter DCA (?)
   -> invariant Mass cut (?)
	  
   2) Ultimately: Select the Xi_c from Xi_cc 
   -> QA Cuts: transverse radius, dca among daughters, dca to PV ... invariant mass
   -> Pion DCA distribution... 
	  
   3) Select the Xi_cc 
   -> Xi_c from Xi_cc Cuts: DCA to PV 
   What this needs to be able: 
   1) process signal + background files 

   Some random documentation: 

   Things to think about: 
   1) Do we compare Topo & Strangeness Tracking? 
   2) Switch between signal and background
*/ 


int main(int argc, char **argv) {
  
  const char* fileName = argv[1]; 
  const char* outAddon = (argv[2])?argv[2]:""; 
  HarryPlotter::StyleBox(); 

  ROOT::EnableImplicitMT(); // Tell ROOT you want to go parallel          
  TString filePath = TString::Format("%s", fileName); 

  bool forceXic = false; 
  if ( filePath.Contains("xic.root")) {
    forceXic = true;
  }

  TFile *file = new TFile(filePath, "READ");
  
  ROOT::RDataFrame df("fTreeCandidates", file);
  
  auto df_in = forceXic?df.Filter("fTrueXic"):df.Filter("fTrueXic||!fTrueXic"); 
  
  //Towards the Xi for Xic->Xi+2pi 
  //Fill some Histograms
  auto h_df_in_lmb_ddca = df_in.Histo1D({"df_in_lmb_ddca", "lmb prong dca", 500, 0, 500}, "fXiV0DauDCA"); 
  auto h_df_in_lmb_trad = df_in.Histo1D({"df_in_lmb_trad", "lmb trad", 500, 0, 10}, "fV0DecayRadius"); 
  auto h_df_in_lmb_trad_mc = df_in.Histo1D({"df_in_lmb_trad_mc", "lmb trad", 500, 0, 10}, "fV0DecayRadiusMC"); 
  auto h_df_in_lmb_trad_diff = df_in.Define("XiV0DecayRadDiff", "TMath::Abs(fV0DecayRadiusMC-fV0DecayRadius)").Histo1D({"df_in_lmb_trad_diff", "lmb trad", 500, 0, 2}, "XiV0DecayRadDiff"); 
  
  auto h_df_in_xi_ddca = df_in.Histo1D({"df_in_xi_ddca", "xi prong dca", 500, 0, 2000}, "fXiCascDauDCA"); 
  auto h_df_in_xi_trad = df_in.Histo1D({"df_in_xi_trad", "xi trad", 250, 0, 10}, "fXiDecayRadius"); 
  auto h_df_in_xi_trad_mc = df_in.Histo1D({"df_in_xi_trad_mc", "xi trad", 250, 0, 10}, "fXiDecayRadiusMC"); 
  auto h_df_in_xi_trad_diff = df_in.Define("XiDecayRadDiff", "TMath::Abs(fXiDecayRadiusMC-fXiDecayRadius)").Histo1D({"df_in_xi_trad_diff", "xi trad", 250, 0, 2}, "XiDecayRadDiff"); 
    
  auto h_df_in_trad_diff_lmb_xi = df_in.Define("XiLmbDecayRadDiff", "fV0DecayRadius-fXiDecayRadius").Histo1D({"df_in_trad_diff_lmb_xi", "lmb-xi trad", 500, -100, 150}, "XiLmbDecayRadDiff"); 

  auto h_df_in_lmb_mass = df_in.Histo1D({"df_in_lmb_mass", "lmb inv mass", 500, 1., 2.3}, "fLambdaMass"); 
  auto h_df_in_xi_mass = df_in.Histo1D({"df_in_xi_mass", "xi inv mass", 500, 1.2, 2.5}, "fXiMass"); 
  
  //Select Xi
  float radDiffMax = 1.; //cm 
  float invMassDiff = 0.008; //8 MeV/c2 mass window 
  float lmbMass = 1.116; 
  float xiMass = 1.322; 
  float xipTmin = 1.0; 
  float xipTmax = 6.0; 
  int addedHitsMin = 1; 
  
  float xicMass = 2.468; 
  float xiccMass = 3.621; 
  
  auto radCut = [&radDiffMax](float radDiff) { return (radDiff < radDiffMax);};
  auto invMassLmbCut = [&invMassDiff, &lmbMass](float invMass) { return (TMath::Abs(invMass-lmbMass) < invMassDiff); }; 
  auto invMassXiCut = [&invMassDiff, &xiMass](float invMass) { return (TMath::Abs(invMass-xiMass) < invMassDiff); }; 
  auto pTCut = [&xipTmin, &xipTmax](float pT) { return (xipTmin < pT)&&(pT < xipTmax); }; 
  auto hitsCut = [&addedHitsMin](int AddedHits) { return (AddedHits >= addedHitsMin); }; 
  
  auto decLengthXic = [&xicMass](float len, float mom){ return TMath::Abs(mom)>1e-4?len*xicMass/mom:-999; }; 
  auto decLengthXicc = [&xiccMass](float len, float mom){ return TMath::Abs(mom)>1e-4?len*xiccMass/mom:-999; }; 
  
  auto df_xi_im = df_in
    .Define("XiV0DecayRadDiff", "TMath::Abs(fV0DecayRadiusMC-fV0DecayRadius)")
    .Define("XiDecayRadDiff", "TMath::Abs(fXiDecayRadiusMC-fXiDecayRadius)")
    .Define("XiLmbDecayRadDiff", "fV0DecayRadius-fXiDecayRadius")
    .Define("XiXicDecayRadDiffTopo", "fXiDecayRadius-fXicDecayRadiusTopo")
    .Define("XiXicDecayRadDiffStra", "fXiDecayRadius-fXicDecayRadiusStraTrack")
    .Define("DCAxyProd"  , "fXiDCAxyToPVStraTrack*fXicPionDCAxyToPV1*fXicPionDCAxyToPV2*fXicDCAxyToPVStraTrack*fPicDCAxyToPVStraTrack") //O(1e7)
    .Define("DCAxySqSum"  , "fXiDCAxyToPVStraTrack*fXiDCAxyToPVStraTrack+fXicPionDCAxyToPV1*fXicPionDCAxyToPV1+fXicPionDCAxyToPV2*fXicPionDCAxyToPV2+fXicDCAxyToPVStraTrack*fXicDCAxyToPVStraTrack+fPicDCAxyToPVStraTrack*fPicDCAxyToPVStraTrack") //O(1e7)
    .Define("DCAzProd"   , "fXiDCAzToPVStraTrack*fXicPionDCAzToPV1*fXicPionDCAzToPV2*fXicDCAzToPVStraTrack*fPicDCAzToPVStraTrack")
    .Define("DCAzSqSum"  , "fXiDCAzToPVStraTrack*fXiDCAzToPVStraTrack+fXicPionDCAzToPV1*fXicPionDCAzToPV1+fXicPionDCAzToPV2*fXicPionDCAzToPV2+fXicDCAzToPVStraTrack*fXicDCAzToPVStraTrack+fPicDCAzToPVStraTrack*fPicDCAzToPVStraTrack") //O(1e7)
    .Define("DCAProdAdd" , "DCAxyProd*DCAxyProd+DCAzProd*DCAzProd")
    .Define("DCAProdMult", "DCAxyProd*DCAzProd")
    .Define("DCAxyProdXicPi", "fPicDCAxyToPVStraTrack*fXicDCAxyToPVStraTrack")
    .Define("DCAzProdXicPi", "fPicDCAzToPVStraTrack*fXicDCAzToPVStraTrack")
    .Define("fXicInvDecayLengthToPVTopo", decLengthXic, {"fXicDecayDistanceFromPVTopo", "lPXiCTopo"})
    .Define("fXiccInvDecayLengthToPVTopo", decLengthXicc, {"fXiccDecayDistanceFromPVTopo", "lPXiCCTopo"})
    .Define("fXicInvDecayLengthToPVStra", decLengthXic, {"fXicDecayDistanceFromPVStraTrack", "lPXiCStraTrack"})
    .Define("fXiccInvDecayLengthToPVStra", decLengthXicc, {"fXiccDecayDistanceFromPVStraTrack", "lPXiCCStraTrack"})
    .Define("fXicInvDecayLengthToDVTopo", decLengthXic, {"fXiCtoXiLengthTopo", "lPXiCTopo"})
    .Define("fXiccInvDecayLengthToDVTopo", decLengthXicc, {"fXiCCtoXiCLengthTopo", "lPXiCCTopo"})
    .Define("fXicInvDecayLengthToDVStra", decLengthXic, {"fXiCtoXiLengthStraTrack", "lPXiCStraTrack"})
    .Define("fXiccInvDecayLengthToDVStra", decLengthXicc, {"fXiCCtoXiCLengthStraTrack", "lPXiCCStraTrack"})    
    .Filter(hitsCut, {"fXiHitsAdded"})
    .Filter(pTCut, {"fXiPtMC"})
    .Filter(radCut, {"XiV0DecayRadDiff"})
    .Filter(radCut, {"XiDecayRadDiff"})
    .Filter("XiLmbDecayRadDiff > 0")
    ; 
  
  auto h_df_xi_im_lmb_mass = df_xi_im.Histo1D({"df_xi_im_lmb_mass", "lmb inv mass", 500, 1., 2.3}, "fLambdaMass"); 
  auto h_df_xi_im_xi_mass = df_xi_im.Histo1D({"df_xi_im_xi_mass", "xi inv mass", 500, 1.2, 2.5}, "fXiMass"); 
  
  auto df_xi = df_xi_im
    .Filter(invMassLmbCut, {"fLambdaMass"})
    .Filter(invMassXiCut, {"fXiMass"})
    ;
  
  //Towards the Xi_c for Xi_cc->Xi_c+pi: 
  //Fill some histograms 
  //Topoor 
  
  auto h_df_xi_trad_diff_xi_xi_c_topo = df_xi.Histo1D({"df_xi_trad_diff_xi_xi_c_topo", "xi-xi_c trad", 500, -50, 100}, "XiXicDecayRadDiffTopo") ;

  auto h_df_xi_xi_c_mass_topo = df_xi.Histo1D({"df_xi_xi_c_mass_topo", "xi_c inv mass", 700, 1.6, 3.2}, "fXicMassTopo"); 
  
  auto h_df_xi_xi_c_ddca_topo = df_xi.Histo1D({"df_xi_xi_c_ddca_topo", "xi_c prong dca", 500, 0, 700}, "fXicDaughterDCATopo"); 
  auto h_df_xi_xi_c_ddist_pv_topo = df_xi.Histo1D({"df_xi_xi_c_dist_pv_topo", "xi_c decay dist", 1500, 0, 0.60}, "fXicInvDecayLengthToPVTopo"); 
  auto h_df_xi_xi_c_ddist_dv_topo = df_xi.Histo1D({"df_xi_xi_c_dist_dv_topo", "xi_c decay dist", 1500, 0, 0.60}, "fXicInvDecayLengthToDVTopo"); 
  auto h_df_xi_xi_c_trad_topo = df_xi.Histo1D({"df_xi_xi_c_trad_topo", "xi_c trad", 2000, 0, 0.4}, "fXicDecayRadiusTopo"); 

  auto h_df_xi_xi_dca_xy_topo = df_xi.Histo1D({"df_xi_xi_dca_xy_topo", "xi dca xy topo", 1000, -500, 500}, "fXiDCAxyToPVTopo");  
  auto h_df_xi_xi_dca_z_topo = df_xi.Histo1D({"df_xi_xi_dca_z_topo", "xi dca z topo", 1000, -500, 500}, "fXiDCAzToPVTopo");  

  //Stra

  auto h_df_xi_trad_diff_xi_xi_c_stra = df_xi.Histo1D({"df_xi_trad_diff_xi_xi_c_stra", "xi-xi_c trad", 500, -50, 200}, "XiXicDecayRadDiffStra") ;

  auto h_df_xi_xi_c_mass_stra = df_xi.Histo1D({"df_xi_xi_c_mass_stra", "xi_c inv mass", 700, 1.6, 3.2}, "fXicMassStraTrack"); 
   
  auto h_df_xi_xi_c_ddca_stra = df_xi.Histo1D({"df_xi_xi_c_ddca_stra", "xi_c prong dca", 500, 0, 700}, "fXicDaughterDCAStraTrack"); 
  auto h_df_xi_xi_c_ddist_pv_stra = df_xi.Histo1D({"df_xi_xi_c_dist_pv_stra", "xi_c decay dist", 1500, 0, 0.60}, "fXicInvDecayLengthToPVStra"); //redefine as mL/p 
  auto h_df_xi_xi_c_ddist_dv_stra = df_xi.Histo1D({"df_xi_xi_c_dist_dv_stra", "xi_c decay dist", 1500, 0, 0.60}, "fXicInvDecayLengthToDVStra"); //redefine as mL/p 
  auto h_df_xi_xi_c_trad_stra = df_xi.Histo1D({"df_xi_xi_c_trad_stra", "xi_c trad", 2000, 0, 0.4}, "fXicDecayRadiusStraTrack"); 
  
  auto h_df_xi_xi_dca_xy_stra = df_xi.Histo1D({"df_xi_xi_dca_xy_stra", "xi dca xy stra", 1000, -500, 500}, "fXiDCAxyToPVStraTrack");  
  auto h_df_xi_xi_dca_z_stra = df_xi.Histo1D({"df_xi_xi_dca_z_stra", "xi dca z stra", 1000, -500, 500}, "fXiDCAzToPVStraTrack");  
  
  //Pions (from the xic) 
  
  auto h_df_xi_pi_one_dca_xy = df_xi.Histo1D({"df_xi_pi_one_dca_xy", "pi1 dca xy stra", 1000, -500, 500}, "fXicPionDCAxyToPV1");  
  auto h_df_xi_pi_one_dca_z = df_xi.Histo1D({"df_xi_pi_one_dca_z", "pi1 dca z stra", 1000, -500, 500}, "fXicPionDCAzToPV1");  
  
  auto h_df_xi_pi_two_dca_xy = df_xi.Histo1D({"df_xi_pi_two_dca_xy", "pi2 dca xy stra", 1000, -500, 500}, "fXicPionDCAxyToPV2");  
  auto h_df_xi_pi_two_dca_z = df_xi.Histo1D({"df_xi_pi_two_dca_z", "pi2 dca z stra", 1000, -500, 500}, "fXicPionDCAzToPV2");  
  
  //Select the Xi_c 
  float invMassDiffXic = 0.08; //8 MeV/c2 mass window 
   
  auto invMassXicCut = [&invMassDiffXic, &xicMass](float invMass) { return (TMath::Abs(invMass-xicMass) < invMassDiffXic); }; 

  auto df_xi_c_im = df_xi
    .Define("XicXiccDecayRadDiffTopo", "fXicDecayRadiusTopo-fXiccDecayRadiusTopo")
    .Define("XicXiccDecayRadDiffStra", "fXicDecayRadiusStraTrack-fXiccDecayRadiusStraTrack")
    .Filter("XiXicDecayRadDiffStra > 0")
    ;
  
  auto h_df_xi_c_im_xi_c_mass_stra = df_xi_c_im.Histo1D({"df_xi_c_im_xi_c_mass_stra", "xi_c inv mass", 700, 1.6, 3.2}, "fXicMassStraTrack"); 

  auto df_xi_c = df_xi_c_im.Filter(invMassXicCut, {"fXicMassStraTrack"});

  //Towards the actual Xi_cc
  //Fill some histograms 
  
  //Topo
  
  auto h_df_xi_c_trad_diff_xi_xi_c_topo = df_xi_c.Histo1D({"df_xi_c_trad_diff_xi_xi_c_topo", "xi_c-xi_cc trad", 500, -100, 150}, "XicXiccDecayRadDiffTopo") ;

  auto h_df_xi_c_xi_cc_ddca_topo = df_xi_c.Histo1D({"df_xi_c_xi_cc_ddca_topo", "xi_cc prong dca", 500, 0, 500}, "fXiccDaughterDCATopo"); 
  auto h_df_xi_c_xi_cc_ddist_pv_topo = df_xi_c.Histo1D({"df_xi_c_xi_cc_dist_pv_topo", "xi_cc decay dist", 1500, 0, 0.5}, "fXiccInvDecayLengthToPVTopo"); 
  auto h_df_xi_c_xi_cc_ddist_dv_topo = df_xi_c.Histo1D({"df_xi_c_xi_cc_dist_dv_topo", "xi_cc decay dist", 1500, 0, 0.5}, "fXiccInvDecayLengthToDVTopo"); 
  auto h_df_xi_c_xi_cc_trad_topo = df_xi_c.Histo1D({"df_xi_c_xi_cc_trad_topo", "xi_cc trad", 2000, 0, 0.5}, "fXiccDecayRadiusTopo"); 

  auto h_df_xi_c_xi_c_dca_xy_topo = df_xi_c.Histo1D({"df_xi_c_xi_c_dca_xy_topo", "xi_c dca xy topo", 1000, -500, 500}, "fXicDCAxyToPVTopo");  
  auto h_df_xi_c_xi_c_dca_z_topo  = df_xi_c.Histo1D({"df_xi_c_xi_c_dca_z_topo", "xi_c dca z topo", 1000, -500, 500}, "fXicDCAzToPVTopo");  

  auto h_df_xi_c_xi_cc_dca_xy_topo = df_xi_c.Histo1D({"df_xi_c_xi_cc_dca_xy_topo", "xi_cc dca xy topo", 1000, -500, 500}, "fXiccDCAxyToPVTopo");  
  auto h_df_xi_c_xi_cc_dca_z_topo  = df_xi_c.Histo1D({"df_xi_c_xi_cc_dca_z_topo", "xi_cc dca z topo", 1000, -500, 500}, "fXiccDCAzToPVTopo");  
  
  //Stra

  auto h_df_xi_c_trad_diff_xi_xi_c_stra = df_xi_c.Histo1D({"df_xi_c_trad_diff_xi_xi_c_stra", "xi_c-xi_cc trad", 500, -100, 150}, "XicXiccDecayRadDiffStra") ;

  auto h_df_xi_c_xi_cc_ddca_stra = df_xi_c.Histo1D({"df_xi_c_xi_cc_ddca_stra", "xi_cc prong dca", 500, 0, 500}, "fXiccDaughterDCAStraTrack"); 
  auto h_df_xi_c_xi_cc_ddist_pv_stra = df_xi_c.Histo1D({"df_xi_c_xi_cc_dist_pv_stra", "xi_cc decay dist", 1500, 0, 0.50}, "fXiccInvDecayLengthToPVStra"); 
  auto h_df_xi_c_xi_cc_ddist_dv_stra = df_xi_c.Histo1D({"df_xi_c_xi_cc_dist_dv_stra", "xi_cc decay dist", 1500, 0, 0.50}, "fXiccInvDecayLengthToDVStra"); 
  auto h_df_xi_c_xi_cc_trad_stra = df_xi_c.Histo1D({"df_xi_c_xi_cc_trad_stra", "xi_cc trad", 2000, 0, 0.5}, "fXiccDecayRadiusStraTrack"); 
  
  auto h_df_xi_c_xi_c_dca_xy_stra = df_xi_c.Histo1D({"df_xi_c_xi_c_dca_xy_stra", "xi_c dca xy stra", 1000, -500, 500}, "fXicDCAxyToPVStraTrack");  
  auto h_df_xi_c_xi_c_dca_z_stra  = df_xi_c.Histo1D({"df_xi_c_xi_c_dca_z_stra", "xi_c dca z stra", 1000, -500, 500}, "fXicDCAzToPVStraTrack");  
  
  auto h_df_xi_c_xi_cc_dca_xy_stra = df_xi_c.Histo1D({"df_xi_c_xi_cc_dca_xy_stra", "xi_cc dca xy stra", 1000, -500, 500}, "fXiccDCAxyToPVStraTrack");  
  auto h_df_xi_c_xi_cc_dca_z_stra  = df_xi_c.Histo1D({"df_xi_c_xi_cc_dca_z_stra", "xi_cc dca z stra", 1000, -500, 500}, "fXiccDCAzToPVStraTrack");  

  //Pions (from the xicc)
 
  auto h_df_xi_c_pi_dca_xy = df_xi_c.Histo1D({"df_xi_c_pi_dca_xy", "xi_c dca xy stra", 1000, -500, 500}, "fPicDCAxyToPVStraTrack");  
  auto h_df_xi_c_pi_dca_z  = df_xi_c.Histo1D({"df_xi_c_pi_dca_z", "xi_c dca z stra", 1000, -500, 500}, "fPicDCAzToPVStraTrack");  
  

  //some product magic .. 

  auto h_df_xi_c_dcaxyProd   = df_xi_c.Histo1D({"df_xi_c_dca_xy_prod", "dca xy prod", 1000, -1e8, 1e8},"DCAxyProd"  ); 
  auto h_df_xi_c_dcazProd    = df_xi_c.Histo1D({"df_xi_c_dca_z_prod", "dca z prod", 1000, -1e8, 1e8},"DCAzProd"   ); 
  auto h_df_xi_c_dcaProdAdd  = df_xi_c.Histo1D({"df_xi_c_dca_prod_add", "dca prod add", 1000, 0, 1e15},"DCAProdAdd" ); 
  auto h_df_xi_c_dcaProdMult = df_xi_c.Histo1D({"df_xi_c_dca_prod_mult", "dca prod mult", 1000, -1e9, 1e9},"DCAProdMult"); 
  auto h_df_xi_c_dcaxySqSum   = df_xi_c.Histo1D({"df_xi_c_dca_xy_sqsum", "dca xy sqsum", 1000, 0, 1000000},"DCAxySqSum"  ); 
  auto h_df_xi_c_dcazSqSum   = df_xi_c.Histo1D({"df_xi_c_dca_z_sqsum", "dca z sqsum", 1000, 0, 1000000},"DCAzSqSum"  ); 
  
  auto h_df_xi_c_dcaxy_prod_xic_pi   = df_xi_c.Histo1D({"df_xi_c_dca_xy_prod_xic_pi", "dca xy prod", 1000, -1000, 1000},"DCAxyProdXicPi"  ); 
  auto h_df_xi_c_dcaz_prod_xic_pi   = df_xi_c.Histo1D({"df_xi_c_dca_z_prod_xic_pi", "dca z prod", 1000, -1000, 1000},"DCAzProdXicPi"  ); 
  
  auto h_df_xi_c_xi_cc_mass_stra = df_xi_c.Filter("XicXiccDecayRadDiffStra > 0").Histo1D({"df_xi_c_xi_cc_mass_stra", "xi_cc inv mass", 700, 2.6, 4.6}, "fXiccMassStraTrack"); 
  
  //Select the Xi_cc
   
  auto invMassXiccCut = [&invMassDiffXic, &xiccMass](float invMass) { return (TMath::Abs(invMass-xiccMass) < invMassDiffXic); }; 
  
  
  auto df_xi_cc_im_c1 = df_xi_c
    .Filter("XicXiccDecayRadDiffStra > 0")
    .Filter("TMath::Abs(fXicPionDCAxyToPV1) > 30")
    .Filter("TMath::Abs(fXicPionDCAzToPV1) > 35")
    .Filter("TMath::Abs(fXicPionDCAxyToPV2) > 30")
    .Filter("TMath::Abs(fXicPionDCAzToPV2) > 35")
    .Filter("TMath::Abs(fPicDCAxyToPVTopo) > 30")
    .Filter("TMath::Abs(fPicDCAzToPVTopo) > 35")
    ;
  
  //Fill some final histograms  
  auto h_df_xi_cc_im_xi_cc_mass_stra_c1 = df_xi_cc_im_c1.Histo1D({"df_xi_cc_im_xi_cc_mass_stra_c1", "xi_cc inv mass", 700, 2.6, 4.6}, "fXiccMassStraTrack"); 

  auto df_xi_cc_im_c2 = df_xi_c
    .Filter("XicXiccDecayRadDiffStra > 0")
    .Filter("fXicDecayRadiusStraTrack>0.005")
    .Filter("fXiccDecayRadiusStraTrack>0.005")
    .Filter("TMath::Abs(fXicPionDCAxyToPV1) > 30")
    .Filter("TMath::Abs(fXicPionDCAzToPV1) > 35")
    .Filter("TMath::Abs(fXicPionDCAxyToPV2) > 30")
    .Filter("TMath::Abs(fXicPionDCAzToPV2) > 35")
    .Filter("TMath::Abs(fPicDCAxyToPVTopo) > 30")
    .Filter("TMath::Abs(fPicDCAzToPVTopo) > 35")
    ;
  
  //Fill some final histograms  
  auto h_df_xi_cc_im_xi_cc_mass_stra_c2 = df_xi_cc_im_c2.Histo1D({"df_xi_cc_im_xi_cc_mass_stra_c2", "xi_cc inv mass", 700, 2.6, 4.6}, "fXiccMassStraTrack"); 
  
  auto df_xi_cc_im_c3 = df_xi_c
    .Filter("XicXiccDecayRadDiffStra > 0")
    .Filter("fXicDecayRadiusStraTrack>0.005")
    .Filter("fXiccDecayRadiusStraTrack>0.005")
    .Filter("TMath::Abs(fXicDCAxyToPVStraTrack)>10")
    .Filter("TMath::Abs(fXicPionDCAxyToPV1) > 30")
    .Filter("TMath::Abs(fXicPionDCAzToPV1) > 35")
    .Filter("TMath::Abs(fXicPionDCAxyToPV2) > 30")
    .Filter("TMath::Abs(fXicPionDCAzToPV2) > 35")
    .Filter("TMath::Abs(fPicDCAxyToPVTopo) > 30")
    .Filter("TMath::Abs(fPicDCAzToPVTopo) > 35")
    ;
  
  //Fill some final histograms  
  auto h_df_xi_cc_im_xi_cc_mass_stra_c3 = df_xi_cc_im_c3.Histo1D({"df_xi_cc_im_xi_cc_mass_stra_c3", "xi_cc inv mass", 700, 2.6, 4.6}, "fXiccMassStraTrack"); 
  
  auto df_xi_cc_im_c4 = df_xi_c
    .Filter("XicXiccDecayRadDiffStra > 0")
    .Filter("fXicDecayRadiusStraTrack>0.005")
    .Filter("fXiccDecayRadiusStraTrack>0.005")
    .Filter("TMath::Abs(fXiDCAxyToPVStraTrack)>10")
    .Filter("TMath::Abs(fXiDCAzToPVStraTrack)>10")
    .Filter("TMath::Abs(fXicDCAxyToPVStraTrack)>10")
    .Filter("TMath::Abs(fXicDCAzToPVStraTrack)>10")
    .Filter("TMath::Abs(fXicPionDCAxyToPV1) > 20")
    .Filter("TMath::Abs(fXicPionDCAzToPV1) > 20")
    .Filter("TMath::Abs(fXicPionDCAxyToPV2) > 20")
    .Filter("TMath::Abs(fXicPionDCAzToPV2) > 20")
    .Filter("TMath::Abs(fPicDCAxyToPVTopo) > 20")
    .Filter("TMath::Abs(fPicDCAzToPVTopo) > 20")
    ;
  
  //Fill some final histograms  
  auto h_df_xi_cc_im_xi_cc_mass_stra_c4 = df_xi_cc_im_c4.Histo1D({"df_xi_cc_im_xi_cc_mass_stra_c4", "xi_cc inv mass", 700, 2.6, 4.6}, "fXiccMassStraTrack"); 

  auto counter = df_xi_c.Filter("XicXiccDecayRadDiffStra > 0").Count(); 
  
  auto h_cand_counter = new TH1D("df_xi_c_candCounter", "candCounter", 1, 0, 1); 
  h_cand_counter->SetBinContent(1, *counter); 

  TString outName = TString::Format("outxiccSelector_%s.root",outAddon )  ; 
  TFile* out = TFile::Open(outName.Data(), "recreate");
  
  out->cd(); 
  //to xi 
  HarryPlotter::CheckAndStore(out, h_df_in_lmb_ddca);
  HarryPlotter::CheckAndStore(out, h_df_in_lmb_trad); 
  HarryPlotter::CheckAndStore(out, h_df_in_lmb_trad_mc);
  HarryPlotter::CheckAndStore(out, h_df_in_lmb_trad_diff);
  HarryPlotter::CheckAndStore(out, h_df_in_lmb_mass);
  
  HarryPlotter::CheckAndStore(out, h_df_in_xi_ddca);
  HarryPlotter::CheckAndStore(out, h_df_in_xi_trad); 
  HarryPlotter::CheckAndStore(out, h_df_in_xi_trad_mc);
  HarryPlotter::CheckAndStore(out, h_df_in_xi_trad_diff);
  HarryPlotter::CheckAndStore(out, h_df_in_xi_mass);

  HarryPlotter::CheckAndStore(out, h_df_in_trad_diff_lmb_xi); 
  

  //to xi_c
  HarryPlotter::CheckAndStore(out, h_df_xi_im_lmb_mass);
  HarryPlotter::CheckAndStore(out, h_df_xi_im_xi_mass);

  HarryPlotter::CheckAndStore(out, h_df_xi_xi_c_mass_topo); 
  HarryPlotter::CheckAndStore(out, h_df_xi_trad_diff_xi_xi_c_topo);
 
  HarryPlotter::CheckAndStore(out, h_df_xi_xi_c_ddca_topo); 
  HarryPlotter::CheckAndStore(out, h_df_xi_xi_c_ddist_pv_topo); 
  HarryPlotter::CheckAndStore(out, h_df_xi_xi_c_ddist_dv_topo); 
  HarryPlotter::CheckAndStore(out, h_df_xi_xi_c_trad_topo); 
    
  HarryPlotter::CheckAndStore(out, h_df_xi_xi_dca_xy_topo);
  HarryPlotter::CheckAndStore(out, h_df_xi_xi_dca_z_topo);
   
  HarryPlotter::CheckAndStore(out, h_df_xi_xi_c_mass_stra); 
  
  HarryPlotter::CheckAndStore(out, h_df_xi_trad_diff_xi_xi_c_stra);
  HarryPlotter::CheckAndStore(out, h_df_xi_xi_c_ddca_stra); 
  HarryPlotter::CheckAndStore(out, h_df_xi_xi_c_ddist_pv_stra); 
  HarryPlotter::CheckAndStore(out, h_df_xi_xi_c_ddist_dv_stra); 
  HarryPlotter::CheckAndStore(out, h_df_xi_xi_c_trad_stra); 

  HarryPlotter::CheckAndStore(out, h_df_xi_xi_dca_xy_stra);
  HarryPlotter::CheckAndStore(out, h_df_xi_xi_dca_z_stra);
    
  HarryPlotter::CheckAndStore(out, h_df_xi_pi_one_dca_xy);
  HarryPlotter::CheckAndStore(out, h_df_xi_pi_one_dca_z);
  HarryPlotter::CheckAndStore(out, h_df_xi_pi_two_dca_xy);
  HarryPlotter::CheckAndStore(out, h_df_xi_pi_two_dca_z);
  
  HarryPlotter::CheckAndStore(out, h_df_xi_c_im_xi_c_mass_stra); 
  
  //to xi_cc
  HarryPlotter::CheckAndStore(out,h_df_xi_c_trad_diff_xi_xi_c_topo);
  HarryPlotter::CheckAndStore(out,h_df_xi_c_xi_cc_ddca_topo);
  HarryPlotter::CheckAndStore(out,h_df_xi_c_xi_cc_ddist_pv_topo);
  HarryPlotter::CheckAndStore(out,h_df_xi_c_xi_cc_ddist_dv_topo);
  HarryPlotter::CheckAndStore(out,h_df_xi_c_xi_cc_trad_topo);
  HarryPlotter::CheckAndStore(out,h_df_xi_c_xi_c_dca_xy_topo);
  HarryPlotter::CheckAndStore(out,h_df_xi_c_xi_c_dca_z_topo);
  HarryPlotter::CheckAndStore(out, h_df_xi_c_xi_cc_dca_xy_topo);
  HarryPlotter::CheckAndStore(out, h_df_xi_c_xi_cc_dca_z_topo);


  HarryPlotter::CheckAndStore(out,h_df_xi_c_trad_diff_xi_xi_c_stra);
  HarryPlotter::CheckAndStore(out,h_df_xi_c_xi_cc_ddca_stra);
  HarryPlotter::CheckAndStore(out,h_df_xi_c_xi_cc_ddist_pv_stra);
  HarryPlotter::CheckAndStore(out,h_df_xi_c_xi_cc_ddist_dv_stra);
  HarryPlotter::CheckAndStore(out,h_df_xi_c_xi_cc_trad_stra);
  HarryPlotter::CheckAndStore(out,h_df_xi_c_xi_c_dca_xy_stra);
  HarryPlotter::CheckAndStore(out,h_df_xi_c_xi_c_dca_z_stra);
  HarryPlotter::CheckAndStore(out, h_df_xi_c_xi_cc_dca_xy_stra);
  HarryPlotter::CheckAndStore(out, h_df_xi_c_xi_cc_dca_z_stra);
    
  HarryPlotter::CheckAndStore(out,h_df_xi_c_pi_dca_xy);
  HarryPlotter::CheckAndStore(out,h_df_xi_c_pi_dca_z);
  
  HarryPlotter::CheckAndStore(out,h_df_xi_c_dcaxyProd  );
  HarryPlotter::CheckAndStore(out,h_df_xi_c_dcazProd   );
  HarryPlotter::CheckAndStore(out,h_df_xi_c_dcaProdAdd );
  HarryPlotter::CheckAndStore(out,h_df_xi_c_dcaProdMult);
  
  HarryPlotter::CheckAndStore(out,h_df_xi_c_dcaxySqSum);
  HarryPlotter::CheckAndStore(out,h_df_xi_c_dcazSqSum);

  HarryPlotter::CheckAndStore(out,h_df_xi_c_dcaxy_prod_xic_pi);
  HarryPlotter::CheckAndStore(out,h_df_xi_c_dcaz_prod_xic_pi);
  
  HarryPlotter::CheckAndStore(out,h_df_xi_c_xi_cc_mass_stra); 
  
  //xi_cc selected
  HarryPlotter::CheckAndStore(out, h_df_xi_cc_im_xi_cc_mass_stra_c1); 
  HarryPlotter::CheckAndStore(out, h_df_xi_cc_im_xi_cc_mass_stra_c2); 
  HarryPlotter::CheckAndStore(out, h_df_xi_cc_im_xi_cc_mass_stra_c3); 
  HarryPlotter::CheckAndStore(out, h_df_xi_cc_im_xi_cc_mass_stra_c4); 

  HarryPlotter::CheckAndStore(out, h_cand_counter); 

  out->Close(); 
  return 0; 
}

