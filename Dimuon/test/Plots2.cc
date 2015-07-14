#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TRandom.h"
#include "TTreeReader.h"
#include "TCanvas.h"
#include "TTreeReaderValue.h"

void Plots2()
{
  TCanvas *c1 =new TCanvas("c1","Boson",1000,1000);
 TCanvas *c2 = new TCanvas("c2","Muon 1");
 TCanvas *c3 = new TCanvas("c3","Muon 2");


 c1->Divide(1,3);
 c1->cd(1);

  TH1F *h_zMass = new TH1F("h_Zmass", "h_Zmass", 1000, 0., 3500.);
  TH1F *h_zPt = new TH1F("h_Zpt", "h_Zpt", 1000, 0., 2500.);
  TH1F *h_mu1Pt = new TH1F("h_Mu1Pt","h_Mu1Pt", 1000, 0., 2500.);
  TH1F *h_mu1Eta = new TH1F("h_Mu1Eta", "h_Mu1Eta", 100, -10., 10.);
  TH1F *h_mu2Pt = new TH1F("h_Mu2Pt", "h_Mu2Pt", 1000, 0., 2500.);
  TH1F *h_mu2Eta = new TH1F("h_Mu2Eta","h_Mu2Eta", 100, -10., 10.); 
  TH2F *h_zPt_vs_zMass = new TH2F("h_zPt_vs_zMass", "h_zPt_vs_zMass", 1000, 0., 3500., 1000, 0., 2500.);
  TH2F *h_mu1Pt_vs_mu1Eta = new TH2F("h_mu1Pt_vs_mu1Eta", "h_mu1Pt_vs_mu1Eta", 100, -10., 10., 1000, 0., 2500.);
  TH2F *h_mu2Pt_vs_mu2Eta = new TH2F("h_mu2Pt_vs_mu2Eta", "h_mu2Pt_vs_mu2Eta", 100, -10., 10., 1000, 0., 2500.);


  TFile *f = TFile::Open("/afs/cern.ch/user/s/szaleski/CMSSW_744_MCGen/src/GenStudy/Dimuon/DYoutput.root");
  if (f == 0){
    std::cout<<"Error: cannot open DYoutput.root \n";
 return;
  }

  TTreeReader myReader("Dimuon/pdfTree", f);

  std::cout << "\n\nmade it this far\nTTreeReader="<<myReader.GetEntries(true)<<"\n";

TTreeReaderValue<Float_t> fZ(myReader, "bosonP4.mass");
TTreeReaderValue<Float_t> Zpt(myReader, "bosonP4.pt");
TTreeReaderValue<Float_t> mu1pt(myReader, "decay1P4.pt");
TTreeReaderValue<Float_t> mu1eta(myReader, "decay1P4.eta");
TTreeReaderValue<Float_t> mu2pt(myReader, "decay2P4.pt");
TTreeReaderValue<Float_t> mu2eta(myReader, "decay2P4.eta");

 std::cout << "\n\nmade it this far\n\n";
 int counter = 0;


 while (myReader.Next()) {
  std::cout<<"\n\nMMade it into the loop\n\n";
// Do the analysis...
      
//for (int iParticle = 0; iParticle < fZ.GetSize(); ++iParticle) {

            h_zMass->Fill(*fZ);
	    h_zPt->Fill(*Zpt);
	    h_zPt_vs_zMass->Fill(*Zpt, *fZ);	
	    h_mu1Pt->Fill(*mu1pt);
	    h_mu1Eta->Fill(*mu1eta);
	    //h_mu1Pt_vs_mu1Eta->Fill(*mu1eta, *mu1pt);
	    h_mu2Pt->Fill(*mu2pt);
	    h_mu2Eta->Fill(*mu2eta);
	    // h_mu2Pt_vs_mu2Eta->Fill(*mu2eta, *mu2pt);

     std::cout << counter << " Zmass "<< *fZ<<"\n";
	    counter++;
	    //    	        }
 }
 h_zMass->Draw();

 c1->cd(2);
 h_zPt->Draw();

 c1->cd(3);
 h_zPt_vs_zMass->Draw();
 
  c1->Print("Boson.pdf", "pdf");

  c2->Divide(1,3);
  c2->cd(1);

  h_mu1Pt->Draw();

  c2->cd(2);
  h_mu1Eta->Draw();

  c2->Print("Muon1.pdf", "pdf");
 
  c3->Divide(1,3);
  c3->cd(1);

    h_mu2Pt->Draw();

  c3->cd(2);

  h_mu2Eta->Draw();

  c3->Print("Muon2.pdf", "pdf");

  
  //TTree *Zmass = (TTree*)f->Get("h_Zmass");
  /*TTree *Zpt = (TTree*)f->Get("h_Zpt");
  TTree *Mu1pt = (TTree*)f->Get("h_mu1pt");
  TTree *Mu1eta = (TTree*)f->Get("h_mu1eta");
  TTree *Mu2pt = (TTree*)f->Get("h_mu2pt");
  TTree *Mu2eta = (TTree*)f->Get("h_mu2eta");*/
  //Float_t m;
  /*  Float_t zpt;
  Float_t mu1pt;
  Float_t mu1eta;
  Float_t mu2pt;
  Float_t mu2eta;*/
  // Zmass->SetBranchAddress("mass",&m);
  /*Zpt->SetBranchAddress("Zpt",&zpt);
 Mu1pt->SetBranchAddress("Mu1pt",&mu1pt);
 Mu1eta->SetBranchAddress("Mu1eta",&mu1eta);
 Mu2pt->SetBranchAddress("Mu2pt",&mu2pt);
 Mu2eta->SetBranchAddress("Mu2eta",&mu2eta);

 // std::cout << endl << endl << "made it this far" << endl << endl;
 */
  //Create Histogram
  // c1->Divide(1,3);
  // c1->cd(1);
 /* TH1F *hZm = new TH1F("hZm","Invariant mass", 1000, 0., 10000.);
 Long64_t nentries = Zmass->GetEntries();
  for(Long64_t i = 1; i <= nentries; i++){
    Zmass->GetEntry(i);
    hZm->Fill(m);
  }
    hZm->GetXaxis()->SetTitle("number of bosons");
    hZm->GetYaxis()->SetTitle("mass (GeV)");
  hZm->Draw();
 */

  /*  c1->cd(2);
  TH1F *hZpt = new TH1F("hZpt","Transverse momentum", 500, 0., 2500.);
  //  Long64_t nentries = Zpt->GetEntries();
  for(Long64_t i = 1; i <= nentries; i++){
    Zpt->GetEntry(i);
    hZpt->Fill(zpt);
  }

    hZpt->GetXaxis()->SetTitle("number of bosons");
    hZpt->GetYaxis()->SetTitle("p_T");


  hZpt->Draw();

  c1->Print("Boson.pdf", "pdf");

  c2->Divide(1,3);

  c2->cd(1);
  //Create Histogram

  TH1F *hMu1pt = new TH1F("hMu1pt","P_t", 1000, 0., 2500.);
  // Long64_t nentries = Mu1pt->GetEntries();
  for(Long64_t i = 1; i <= nentries; i++){
    Mu1pt->GetEntry(i);
    hMu1pt->Fill(mu1pt);
  }

 hMu1pt->GetXaxis()->SetTitle("number of muons");
    hMu1pt->GetYaxis()->SetTitle("p_T");

  hMu1pt->Draw();

  c2->cd(2);

  //Create Histogram

  TH1F *hMu1eta = new TH1F("hMu1eta","Pseudorapidity", 500, -10., 10.);
  // Long64_t nentries = Mu1eta->GetEntries();
  for(Long64_t i = 1; i <= nentries; i++){
    Mu1eta->GetEntry(i);
    hMu1eta->Fill(mu1eta);
  }

    hMu1eta->GetXaxis()->SetTitle("number of muons");
    hMu1eta->GetYaxis()->SetTitle("Pseudorapidity");


  hMu1eta->Draw();

  c2->Print("muon1.pdf", "pdf");


  c3->Divide(1,3);
 c3->cd(1);
  //Create Histogram

   TH1F *hMu2pt = new TH1F("hMu2pt","Invariant mass", 1000, 0., 2500.);
   // Long64_t nentries = Mu2pt->GetEntries();
  for(Long64_t i = 1; i <= nentries; i++){
//    Mu2pt->GetEntry(i);
    hMu2pt->Fill(mu2pt);
  }

    hMu2pt->GetXaxis()->SetTitle("number of muons");
    hMu2pt->GetYaxis()->SetTitle("p_T");

  hMu2pt->Draw();

  c3->cd(2);

  //Create Histogram

  TH1F *hMu2eta = new TH1F("hMu2eta","Pseudorapidity", 500,-10., 10.);
  // Long64_t nentries = Mu2eta->GetEntries();
  for(Long64_t i = 1; i <= nentries; i++){
    Mu2eta->GetEntry(i);
    hMu2eta->Fill(mu2eta);
  }

    hMu2eta->GetXaxis()->SetTitle("number of muons");
    hMu2eta->GetYaxis()->SetTitle("Pseudorapidity");


  hMu2eta->Draw();

c3->Print("muon2.pdf", "pdf");

  */
  return;
}
