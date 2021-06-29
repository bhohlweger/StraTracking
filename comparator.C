void comparator(TString addon, TString mode) { 
  double xiccMass = 3.596; //3.621; 
  double xiccWindow = 0.12;
  
  std::vector<TFile*> inFiles; 

  TSystemDirectory dir("ZeDirectory", gSystem->pwd());
  auto files = dir.GetListOfFiles();
  for (auto fileObj : *files)  {
    auto file = (TSystemFile*) fileObj;
    TString fileName = TString::Format("%s", file->GetName()); 
    if (fileName.Contains(addon) && !fileName.Contains("outComp_cut")) { 
      std::cout << "Reading file " << file->GetName() << std::endl;
      inFiles.emplace_back(TFile::Open(file->GetName(), "read"));
    }
  }
  std::vector<int> colors = {kPink+7, 38, kGreen+3, kOrange+7, kBlack, kViolet+1, kCyan+1}; 
  
  std::vector<TString> histnames;
  if (mode.Contains("ptcomp")) { 
    std::cout <<"You chose mode pT comparison!\n"; 
    histnames =  {"0-2 GeV/#it{c}", "2-4 GeV/#it{c}", "4-7 GeV/#it{c}", "7-10 GeV/#it{c}"};
  } else if (mode.Contains("sgnbkgcomp")) {
    std::cout <<"You chose mode Signal/Background comparison!\n"; 
    histnames = {"#Xi^{#minus} + 3#times#pi_{Pythia}", "#Xi^{+}_{c} + #pi_{Pythia}", "#Xi_{cc}^{++}"}; 
  } else {
    std::cout << "Mode " << mode << " not supported. Choose from (1) ptcomp or (2) sgnbkgcomp \n"; 
  }
  
  TFile* output = TFile::Open(TString::Format("outComp_cut%s.root", addon.Data()), "recreate"); 
  
  TList* inList = inFiles[0]->GetListOfKeys(); 
  TIter next(inList); 
  TObject* obj = nullptr; 
  
  while ((obj = next())) {
    TString objName = TString::Format("%s",obj->GetName()); 
    std::cout << objName.Data() << std::endl; 
    if (objName.Contains("Counter")  || objName.Contains("_vs_") ) { 
      continue; 
    }
    if ((objName.Contains("mass") && objName.Contains("df_xi_cc"))) { 
      continue; 
    }

    auto c = c11(obj->GetName()); 
    auto p = (TPad*)gROOT->FindObject(TString::Format("p%s",obj->GetName()).Data()); 
    p->cd(); 
    
    auto leg = new TLegend(0.16, 0.35, 0.56, 0.8, "#splitline{ALICE 3 (Layout v1) Full Simulation}{Pythia pp #sqrt{s} = 13 TeV + GEANT3}");
    leg->SetFillStyle(0); 

    p->cd(); 
    
    int fileCounter = 0; 
    
    for (auto it : inFiles) { 
      TH1D* Hist = (TH1D*)it->Get(obj->GetName()); 
      if (!Hist) { 
	continue;
      }
      Hist->SetLineColor(colors[fileCounter]); 
      Hist->SetMarkerColor(colors[fileCounter]); 
 
      Hist->GetYaxis()->SetTitleOffset(1.2); 
      Hist->GetXaxis()->SetNdivisions(506); 
      Hist->GetYaxis()->SetNdivisions(504); 
	
      Hist->SetTitle(histnames[fileCounter].Data()); 

      Hist->Scale(1./Hist->GetMaximum()); 
    
      Hist->GetXaxis()->SetTitle(obj->GetName()); 
      Hist->GetYaxis()->SetRangeUser(0.,2.5); 
      Hist->GetYaxis()->SetTitle("Count Normalized to Maximum"); 
      fileCounter==0?Hist->Draw("hist"):Hist->Draw("histsame"); 
      leg->AddEntry(Hist, Hist->GetTitle(), "l");
      fileCounter++;
    }
    leg->Draw("same"); 
    output->cd();     
    c->Write();
    c->Close(); 
  }
  output->Close(); 
  for (auto it : inFiles) { 
    it->Close();
  }
}
