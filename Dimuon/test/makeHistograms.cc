#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TRandom.h"
#include "TH1.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TCanvas.h"
#include "TVector2.h"
#include "TPad.h"
#include "TPaveStats.h"
#include "TString.h"


// function prototype     
// the name of formal parameter can ba anything in a function protoype
void Plot( TCanvas *c, TH1F * hist, TFile *f, std::string const & quantity, std::string const & identifier, std::string const & units, int isLogy);


// genStudyHistogram is equivalent to fucntin main in c++
void makeHistograms(std::string const & inputfile, std::string const & particleType, std::string const& outName) {
  
  TCanvas *c = new TCanvas("c","Invariant Mass",1000,1000);
  TFile *outfile = new TFile(TString(outName+".root"), "RECREATE");
  // Create a histogram for the values we read.
  //TH1F *h_InvariantMass = new TH1F("hInvariantMass", "demoInvarinatMass", 100, 0, 1000);
  
  // Hisgrogram  for Invariant mass of Dimuon pair
  TH1F *h_InvariantMass = new TH1F("hInvariantMass", TString("di"+particleType+"_InvarinatMass"), 500, 0, 5000);
  
  TH1F *h_cosTheta = new TH1F("h_cosTheta", TString("di"+particleType+"_cosTheta"), 100, -1.1, 1.1);
  TH1F *h_tanPhi = new TH1F("h_tanPhi", TString("di"+particleType+"_tanPhi"), 100, -5.1, 5.1);
  TH1F *h_thetaCS = new TH1F("h_thetaCS", TString("di"+particleType+"_thetaCS"), 100, -0.2, 3.2);
  TH1F *h_phiCS = new TH1F("h_phiCS", TString("di"+particleType+"_phiCS"), 100,-3.15, 3.15);

  // Histogram for kinematic variables of Muon1
  TH1F *h_Pt1 = new TH1F("hPt1", TString(particleType+"Pt1"), 100, 0, 1000);
  TH1F *h_Eta1 = new TH1F("hEta1", TString(particleType+"Eta1"), 100, -5.0, 5.0);
  TH1F *h_Phi1 = new TH1F("hPhi1", TString(particleType+"Phi1"), 100, -4.0, 4.0);
  TH1F *h_Energy1 = new TH1F("hEnergy1", TString(particleType+"Energy1"), 100, 0, 1000);
  TH1F *h_Et1 = new TH1F("hEt1", TString(particleType+"Et1"), 100, 0, 1000);
  TH1F *h_Mass1 = new TH1F("hMass1", TString(particleType+"Mass1"), 100, 0, 1000);
  TH1F *h_Theta1 = new TH1F("hTheta1", TString(particleType+"Theta1"), 100, -5.0, 5.0);
  
  // Histograms for kinematic variable of particleType2
  TH1F *h_Pt2 = new TH1F("hPt2", TString(particleType+"Pt2"), 100, 0, 1000);
  TH1F *h_Eta2 = new TH1F("hEta2", TString(particleType+"Eta2"), 100, -5.0, 5.0);
  TH1F *h_Phi2 = new TH1F("hPhi2", TString(particleType+"Phi2"), 100, -4.0, 4.0);
  TH1F *h_Energy2 = new TH1F("hEnergy2", TString(particleType+"Energy2"), 100, 0, 1000);
  TH1F *h_Et2 = new TH1F("hEt2", TString(particleType+"Et2"), 100, 0, 1000);
  TH1F *h_Mass2 = new TH1F("hMass2", TString(particleType+"Mass2"), 100, 0, 1000);
  TH1F *h_Theta2 = new TH1F("hTheta2", TString(particleType+"Theta2"), 100, -4.0, 4.0);
  
  
  // Histograms for kinematic variable of boson                                                                                                                                     
  TH1F *h_boson_Pt = new TH1F("hBosonPt", "BosonPt", 100, 0, 1000);
  TH1F *h_boson_Eta = new TH1F("hBosonEta", "BosonEta", 100, -5.0, 5.0);
  TH1F *h_boson_Phi = new TH1F("hBosonPhi", "BosonPhi", 100, -4.0, 4.0);
  TH1F *h_boson_Energy = new TH1F("hBosonEnergy", "BosonEnergy", 100, 0, 1000);
  TH1F *h_boson_Et = new TH1F("hBosonEt", "BosonEt", 100, 0, 1000);
  TH1F *h_boson_Mass = new TH1F("hBosonMass", "BosonMass", 100, 0, 1000);
  TH1F *h_boson_Theta = new TH1F("hBosonTheta", "BosonTheta", 100, -4.0, 4.0);
  
  
  // Histograms for PID
  /* TH1F *h_Decay1_PID = new TH1F("hDecay1PID", "Decay1PID", 100, 0, 1000);
     TH1F *h_Decay2_PID = new TH1F("hDecay2PID", "Decay1PID", 100, 0, 1000);
     TH1F *h_Boson_PID = new TH1F("hBosonPID", "BosonPID", 100, 0, 1000);
     
     // Histogram for Cross section
     TH1F *h_cross_section = new TH1F("hcrossSection", "CrossSection", 100, 0, 1000);
  */
  // Open the file containing the tree. inputfile is given from command  line.
  // Syntax for compiling and passing the file name inputfile from command line:
  // root -x -b -q genStudyHistogram.cc++(\"Pythia8_Sep21_CI_L13000_13TeV_Sep21.root\")
  
  //   TFile *myFile = TFile::Open(TString("/afs/cern.ch/work/s/szaleski/private/CMSSW_7_6_3_patch2/src/WSUDiLeptons/CosmicEndpoint/test/" + inputfile) );
  TFile *myFile = TFile::Open(TString("$CMSSW_BASE/src/GenStudy/Dimuon/test/" + inputfile) );
  
  // Suceeds if a pointer myFile is a null
  // Null pointer is assumed to point nothing
  
  if (!myFile)
    {

      std::cout << "Failed to open the file" << std::endl;
      return;
    }
  else 
    {
      std::cout << "File opened sucessfully" << std::endl;
      
    }// end of checkig input file 
  
  TTreeReader massReader ("Dimuon/pdfTree", myFile);

  // kinematic varibales for decay1 
  TTreeReaderValue<Float_t> decay1Pt(massReader, "decay1P4.pt");
  TTreeReaderValue<Float_t> decay1Eta(massReader, "decay1P4.eta");
  TTreeReaderValue<Float_t> decay1Phi(massReader, "decay1P4.phi");
  // Added items
  TTreeReaderValue<Float_t> decay1Energy(massReader, "decay1P4.energy");
  TTreeReaderValue<Float_t> decay1Et(massReader, "decay1P4.et");
  TTreeReaderValue<Float_t> decay1Mass(massReader, "decay1P4.mass");
  TTreeReaderValue<Float_t> decay1Theta(massReader, "decay1P4.theta");
  
  
  // kinematic variables for Decay2
  TTreeReaderValue<Float_t> decay2Pt(massReader, "decay2P4.pt");
  TTreeReaderValue<Float_t> decay2Eta(massReader, "decay2P4.eta");
  TTreeReaderValue<Float_t> decay2Phi(massReader, "decay2P4.phi");
  //Added items
  TTreeReaderValue<Float_t> decay2Energy(massReader, "decay2P4.energy");
  TTreeReaderValue<Float_t> decay2Et(massReader, "decay2P4.et");
  TTreeReaderValue<Float_t> decay2Mass(massReader, "decay2P4.mass");
  TTreeReaderValue<Float_t> decay2Theta(massReader, "decay2P4.theta");
  
  // kinematic variables for boson
  TTreeReaderValue<Float_t> bosonPt(massReader, "bosonP4.pt");
  TTreeReaderValue<Float_t> bosonEta(massReader, "bosonP4.eta");
  TTreeReaderValue<Float_t> bosonPhi(massReader, "bosonP4.phi");                                                                                                              
  TTreeReaderValue<Float_t> bosonEnergy(massReader, "bosonP4.energy");
  TTreeReaderValue<Float_t> bosonEt(massReader, "bosonP4.et");
  TTreeReaderValue<Float_t> bosonMass(massReader, "bosonP4.mass");
  TTreeReaderValue<Float_t> bosonTheta(massReader, "bosonP4.theta");
  
  // TTreeReaderValues for PID
  /*TTreeReaderValue<Int_t> decay1_PID(massReader, "decay1PID");
    TTreeReaderValue<Int_t> decay2_PID(massReader, "decay2PID");
    TTreeReaderValue<Int_t> boson_PID(massReader, "bosonPID");
  */
  
  TTreeReaderValue<Double_t> cosTheta(massReader, "cosTheta");
  TTreeReaderValue<Double_t> tanPhi(massReader, "tanPhi");
  TTreeReaderValue<Double_t> thetaCS(massReader, "csTheta");
  TTreeReaderValue<Double_t> phiCS(massReader, "csPhi");
  
  // TTreeReaerValye for cross section
  // TTreeReaderValue<Double_t> cross_Section(massReader, "crossSec");
  
  while (massReader.Next())
    { 
      // if((*decay1Eta < 2.4 && *decay2Eta < 2.1) || (*decay1Eta < 2.1 && *decay2Eta < 2.4)){
      
      // Fills histogram for Invariant Mass of Dimuon Pair
      h_InvariantMass->Fill(sqrt( 2 * (*decay1Pt) * (*decay2Pt)*(cosh ((* decay1Eta) - (*decay2Eta)  ) - cos(TVector2::Phi_mpi_pi((* decay1Phi)- (*decay2Phi)) ))));
      
      h_cosTheta->Fill(*cosTheta);
      h_tanPhi->Fill(*tanPhi);
      h_thetaCS->Fill(*thetaCS);
      h_phiCS->Fill(*phiCS);
      
      
      // Fills histogram for kinematic variables of Decay1
      h_Pt1->Fill(*decay1Pt);
      h_Eta1->Fill(*decay1Eta);
      h_Phi1->Fill(*decay1Phi);
      h_Energy1->Fill(*decay1Energy);
      h_Et1->Fill(*decay1Et);
      h_Mass1->Fill(*decay1Mass);
      h_Theta1->Fill(*decay1Theta);
      
      // Fills histogram for Kinematic variales of Decay2
      h_Pt2->Fill(*decay2Pt);
      h_Eta2->Fill(*decay2Eta);
      h_Phi2->Fill(*decay2Phi);
      h_Energy2->Fill(*decay2Energy);
      h_Et2->Fill(*decay2Et);
      h_Mass2->Fill(*decay2Mass);
      h_Theta2->Fill(*decay2Theta);
      
      // Fills histogram for Kinematic variales of boson                                                                                                                        
      h_boson_Pt->Fill(*bosonPt);
      h_boson_Eta->Fill(*bosonEta);
      h_boson_Phi->Fill(*bosonPhi);
      h_boson_Energy->Fill(*bosonEnergy);
      h_boson_Et->Fill(*bosonEt);
      h_boson_Mass->Fill(*bosonMass);
      h_boson_Theta->Fill(*bosonTheta);
      
      
      //Fills histogram for PIDs
      /*h_Decay1_PID->Fill(*decay1_PID);
	h_Decay2_PID->Fill(*decay2_PID);
	h_Boson_PID->Fill(*boson_PID);
      */
      // Fills histogram for cross section
      //h_cross_section->Fill(*cross_Section);
      //}
    }
  
  // function call for Invariant Masss of Dileptn pair  
  
  std::string label;
  
  if (particleType == "electron" || particleType == "Electron") label = "ee";
  else label = "#mu#mu";

  int isLogy = 1;
  Plot(c, h_InvariantMass, outfile, "m", label, "[GeV]",isLogy);
  //function call for kinematic variables of Muon1
  /*Plot(c, h_Pt1, outfile, "p_T1", label, "[GeV]", isLogy);
  Plot(c, h_Mass1, outfile, "m_01", label, "[GeV]", isLogy);
  Plot(c, h_Pt2, outfile, "p_T2", label, "[GeV]",isLogy);
  Plot(c, h_Mass2, outfile, "m_02", label,"[GeV]", isLogy);
  Plot(c, h_boson_Pt, outfile, "Boson_pT", "Z","[GeV]", isLogy);
  Plot(c, h_boson_Mass, outfile, "Boson_m", "Z" ,"[GeV]", isLogy);
   
  isLogy =0;
 //Plot call for muon 1
  Plot(c, h_Eta1, outfile, "#eta_1", label, "", isLogy);
  Plot(c, h_Phi1, outfile, "#phi_1", label, "[rad]", isLogy);
  Plot(c, h_Energy1, outfile, "E_1", label, "[GeV]", isLogy );
  Plot(c, h_Et1, outfile, "E_T1", label, "[GeV]", isLogy);
  Plot(c, h_Theta1, outfile, "#theta_1", label, "[rad]", isLogy);
  
  //Plot call for muon 2
  Plot(c, h_Eta2, outfile, "#eta_2", label, "", isLogy);
  Plot(c, h_Phi2, outfile, "#phi_2", label, "[rad]", isLogy);
  Plot(c, h_Energy2, outfile, "E_2", label, "[GeV]", isLogy);
  Plot(c, h_Et2, outfile, "E_T2", label, "[GeV]", isLogy);
  Plot(c, h_Theta2, outfile, "#theta_2", label, "[rad]", isLogy);
  
  //Plot call for boson
  Plot(c, h_boson_Eta, outfile, "Boson_eta", label, "", isLogy);
  Plot(c, h_boson_Phi, outfile, "Boson_phi", label, "[rad]", isLogy);
  Plot(c, h_boson_Energy, outfile, "Boson_E", label, "[GeV]", isLogy);
  Plot(c, h_boson_Et, outfile, "Boson_ET", label, "[GeV]", isLogy);
  Plot(c, h_boson_Theta, outfile,"Boson_theta", label, "[rad]", isLogy);
  
  // function call for PIDs
  Plot(c, h_Decay1_PID, outfile, "Decay1 PID in numbers","Decay1_PID", isLogy);
  Plot(c, h_Decay2_PID, outfile, "Decay2 PID in numbers","Decay2_PID", isLogy);
  Plot(c, h_Boson_PID, outfile, "Boson  PID in numbers","Boson_PID", isLogy);
  */
  // function call  for cross section
  //   Plot(c, h_cross_section, outfile, "Cross section in pico barn,","cross_Section", isLogy);
  
  
  outfile->Write();
  outfile->Close();
  return;
  
}

// function definition Plot
void Plot( TCanvas *c, TH1F * hist, TFile *f, std::string const & quantity, std::string const & identifier, std::string const & units, int isLogy)
{
  

  f->cd();
  c->cd();
  // h_InvariantMass->SetMarkerColor(kRed + 2);                                                                                                                                   
  hist->SetLineColor(kAzure + 4);
  hist->SetLineWidth(2);
  hist->SetLineStyle(1);
  //hist->SetStats(0);                       // Remove Stats box
  //     h_InvariantMass->SetMarkerStyle(33);                                                                                                                                     
  hist->Draw();
  
  // sets Y axis into LogY if we pass the argument 1 in a function call
  // Y axis stays same if we pass the argument 0 in a functin call
  if(isLogy == 1)      c->SetLogy(1);
  else c->SetLogy(0);
  // end of if statement to set Y into LogY

  
  hist->GetXaxis()->SetTitle(TString(quantity + "_{" + identifier + "} " + units));
  hist->GetYaxis()->SetTitle(TString("N_{"+identifier+"}"));
  c->Update();
  c->SaveAs(TString(identifier+"_"+quantity + ".png"));
  //  f->cd();
  c->Write();
  return;

}//end of function definition Plot
