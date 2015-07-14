#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TRandom.h"



void CanvasPartition(TCanvas *C, const Int_t Nx = 3, const Int_t Ny = 3,
		     Float_t lMargin = 0.15, Float_t rMargin = 0.05,
		     Float_t bMargin = 0.15, Float_t tMargin = 0.15);

void Canvas3()
{

  gStyle->SetOptStat(0);

  TCanvas *C = (TCanvas*) gROOT->FinObject("C");
  if (C) delete C;
  C = new TCanvas("C","canvas",1024,640);
  C->SetFillStyle(4000);

  //Number of Pads

  const Int_t Nx = 3;
  const Int_t Ny = 3;


  //Margins
   Float_t lMargin = 0.12;
  Float_t rMargin = 0.05;
  Float_t bMargin = 0.15;
  Float_t tMargin = 0.05;
  
  //Canvas Setup
  CanvasPartition(C,Nx,Ny,lMargin,rMargin,bMargin,tMargin);

  TFile *f = new TFile("DYToMuMu_M%d_13TeV_pythia8_GEN.root");
  TTree *Zmass = (TTree*)f->get("h_Zmass");
  Float_t m;
  Zmass->SetBranchAddress("mass",&m);

  //Create Histogram

  TH1F *hZm = new TH1F("hZm","Invariant mass", 1000, 0, 10000);
  Long64_t nentries = Zmass->GetEntries();
  for(Long64_t i = 0; i < nentries; i++){
    Zmass->GetEntry(i);
    hZm->Fill(m);
    hZm->GetXaxis()->SetTitle("number of bosons");
    hZm->GetYaxis()->SetTitle("mass (GeV)");
  }

 TFile *f = new TFile("DYToMuMu_M%d_13TeV_pythia8_GEN.root");
  TTree *Zpt = (TTree*)f->get("h_Zpt");
  Float_t pt;
  Zmass->SetBranchAddress("pt",&pt);

  //Create Histogram

  TH1F *hZpt = new TH1F("hZpt","Transverse momentum", 500, 0, 2500);
  Long64_t nentries = Zpt->GetEntries();
  for(Long64_t i = 0; i < nentries; i++){
    Zpt->GetEntry(i);
    hZpt->Fill(pt);
    hZpt->GetXaxis()->SetTitle("number of bosons");
    hZpt->GetYaxis()->SetTitle("p_T");
  }


  //Create Histogram

  TH2F *hZm_vs_Zpt = new TH2F("hZm_vs_Zpt","Invariant mass vs. pt", 1000, 0, 10000);
  Long64_t nentries = hZm_vs_Zpt->GetEntries();
  for(Long64_t i = 0; i < nentries; i++){
    hZm_vs_Zpt->GetEntry(i);
    hZm_vs_Zpt->Fill(m,pt);
    hZm_vs_Zpt->GetXaxis()->SetTitle("boson mass");
    hZm_vs_Zpt->GetYaxis()->SetTitle("boson p_T");
  }

  TPad *pad[Nx][Ny];
