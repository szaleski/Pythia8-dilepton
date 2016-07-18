// -*- C++ -*-
//
// Package:    GenStudy/Dimuon
// Class:      Dimuon
// 
/**\class Dimuon Dimuon.cc GenStudy/Dimuon/plugins/Dimuon.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Shawn Gregory Zaleski
//         Created:  Tue, 30 Jun 2015 14:01:36 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include <FWCore/ServiceRegistry/interface/Service.h>
#include <CommonTools/UtilAlgos/interface/TFileService.h>
#include <TTree.h>
#include <TVector2.h>
#include <TH1F.h>
#include <TH2F.h>
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include <vector>
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
//
// class declaration
//

class Dimuon : public edm::EDAnalyzer {
public:
  explicit Dimuon(const edm::ParameterSet&);
  ~Dimuon();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  const reco::Candidate* getDaughter(const reco::Candidate* part,int pid);
  const reco::Candidate* getLastDaughter(const reco::Candidate* part,int pid);
  const reco::Candidate* getBoson( const reco::GenParticleCollection& genParts);
  const reco::Candidate* getMother(const reco::Candidate* part, int pid);
  //const reco::Candidate* getDYBoson(const reco::Candidate* part int pid)
  bool isBoson(int pid);
  bool isMuon(int pid);
  bool checkBosonStatus(const reco::GenParticleCollection& genParts);


  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;



  struct P4Struct {
    float energy,et,eta,phi,pt,mass,theta;
    void fill(const math::XYZTLorentzVector& p4){
      if(p4.Pt()!=0 && p4.Et()!=0){
	energy = p4.E();
	et = p4.Et();
	eta = p4.Eta();
	phi = p4.Phi();
	pt =p4.Pt();
	mass = p4.mag();
	theta = p4.Theta();
      }else clear();
    }
    void clear(){energy=-999;et=-999;eta=-999;phi=-0;pt=-999;mass=-999;theta=-999;}
    static std::string contents(){return "energy/F:et:eta:phi:pt:mass:theta";}
  }; 

  TTree* tree_;



  TH1F * h_Zmass, *h_Zpt,*h_Zeta,*h_Zphi,*h_Zcharge;
  TH1F *h_muMinusmass,*h_muMinuspt,*h_muMinuseta,*h_muMinusphi,*h_muMinuscharge;
  TH1F *h_muPlusmass,*h_muPluspt,*h_muPluseta,*h_muPlusphi,*h_muPluscharge;
  TH1F *h_dphi,*h_dtheta, *h_dr, *h_thetaMuMinus,*h_thetaMuPlus;
  TH1F *h_massInvar, *h_dimuonPt, *h_dimuonEta, *h_dimuonPhi;

  TH2F *h2_pt1_vs_pt2,*h2_eta1_vs_eta2,*h2_phi1_vs_phi2;

  P4Struct bosonP4_; // as a sanity check we have the right event...
  P4Struct muMinusP4_;
  P4Struct muPlusP4_;
  int muMinusPID_;
  int muPlusPID_;
  int bosonId_;
  double crossSec;

  bool isCI_, debug_;
  edm::InputTag genPartsTag_;
  int decayParticlePID_, status_, leptonId_;


  // ----------member data ---------------------------
  
};

void Dimuon::beginJob()
{
  edm::Service<TFileService> fs;
  h_Zmass = fs->make<TH1F>("Zmass" , "m", 1000, 0., 600);
  h_Zpt  = fs->make<TH1F>( "Zpt"  , "p_{t}", 500,  0., 2500. );
  h_Zeta = fs->make<TH1F>( "Zeta" , "#eta" , 100, -10., 10.    );
  h_Zphi = fs->make<TH1F>( "Zphi" , "#phi" , 100,  -3.20, 3.20   );
  h_Zcharge = fs->make<TH1F>( "Zcharge" , "Q" ,3,  -1.5, 1.5    );
  h_muMinusmass = fs->make<TH1F>("muMinusmass" , "m", 1000, 0., 500);
  h_muMinuspt  = fs->make<TH1F>( "muMinuspt"  , "p_{t}", 500,  0., 2500. );
  h_muMinuseta = fs->make<TH1F>( "muMinuseta" , "#eta" , 100, -5., 5.    );
  h_muMinusphi = fs->make<TH1F>( "muMinusphi" , "#phi" , 100,  -3.15, 3.15   );
  h_muMinuscharge = fs->make<TH1F>( "muMinuscharge" , "Q" ,3,  -1.5, 1.5    );

  h_muPlusmass = fs->make<TH1F>("muPlusmass" , "m", 1000, 0., 500);
  h_muPluspt  = fs->make<TH1F>( "muPluspt"  , "p_{t}", 500,  0., 2500. );
  h_muPluseta = fs->make<TH1F>( "muPluseta" , "#eta" , 100, -5., 5.    );
  h_muPlusphi = fs->make<TH1F>( "muPlusphi" , "#phi" , 100,  -3.15, 3.15   );
  h_muPluscharge = fs->make<TH1F>( "muPluscharge" , "Q" ,3,  -1.5, 1.5    );

  h_dphi = fs->make<TH1F>("delta phi", "#delta #phi", 100, -3.15, 3.15 );       
  h_dtheta = fs->make<TH1F>("delta theta", "#delta #theta", 100, -3.15, 3.15); 
  h_dr = fs->make<TH1F>("delta r", "#delta r", 100, 0, 10);
  h_thetaMuMinus = fs->make<TH1F>("theta muMinus", "#theta", 100, -3.15, 3.15);      
  h_thetaMuPlus = fs->make<TH1F>("theta muPlus", "#theta", 100, -3.15, 3.15); 
  h_massInvar = fs->make<TH1F>("Invariant mass", "Invariant mass", 1000, 0., 600.);
  h_dimuonPt = fs->make<TH1F>("Dimuon Pt", "Dimuon Pt", 500, 0, 2500);
  h_dimuonEta = fs->make<TH1F>("Dimuon eta", "Dimuon #eta", 100, -5, 5);
  h_dimuonPhi = fs->make<TH1F>("Dimuon Phi", "Dimuon #phi", 100, -3.15, 3.15);

  h2_pt1_vs_pt2   = fs->make<TH2F>( "pt1_vs_pt2"   , "p_{t,1} vs. p_{t,2}"   , 500,  0., 2500., 500,  0., 2500.);
  h2_eta1_vs_eta2 = fs->make<TH2F>( "eta1_vs_eta2" , "#eta_{1} vs. #eta_{2}" , 100, -5., 5.   , 100, -5., 5.   );
  h2_phi1_vs_phi2 = fs->make<TH2F>( "phi1_vs_phi2" , "#phi_{1} vs. #phi_{2}" , 100,  -3.15, 3.15  , 100,  -3.15, 3.15  );

  tree_= fs->make<TTree>("pdfTree","PDF Tree");
  // tree_->Branch("evtId",&evtId_,EventId::contents().c_str());
  tree_->Branch("bosonP4",&bosonP4_,P4Struct::contents().c_str());
  tree_->Branch("decay1P4",&muMinusP4_,P4Struct::contents().c_str());
  tree_->Branch("decay2P4",&muPlusP4_,P4Struct::contents().c_str());
  tree_->Branch("decay1PID",&muMinusPID_,"decay1PID/I");
  tree_->Branch("decay2PID",&muPlusPID_,"decay2PID/I");
  tree_->Branch("bosonPID",&bosonId_,"bosonPID/I");
  tree_->Branch("crossSec", &crossSec, "crossSec/D");
  // tree_->Branch("pdfInfo",&pdfInfo_,PDFInfo::contents().c_str());
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
Dimuon::Dimuon(const edm::ParameterSet& iConfig)

{
  debug_=iConfig.getParameter<bool>("debug");
  isCI_=iConfig.getParameter<bool>("isCI");
  status_=iConfig.getParameter<int>("status");
  genPartsTag_=iConfig.getParameter<edm::InputTag>("genPartsTag");
  decayParticlePID_ = iConfig.getParameter<int>("decayParticlePID");
  leptonId_ = iConfig.getParameter<int>("leptonId");
  //now do what ever initialization is needed

}


Dimuon::~Dimuon()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
Dimuon::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //  isCI_ = true;
  using namespace edm;
  edm::Handle<reco::GenParticleCollection> genPartsHandle;
  iEvent.getByLabel(genPartsTag_,genPartsHandle);
  const reco::GenParticleCollection& genParts = *genPartsHandle;


  bosonId_=0;
  bosonP4_.clear();
  muMinusP4_.clear();
  muPlusP4_.clear();
  muMinusPID_=0;
  muPlusPID_=0;

  const reco::Candidate* boson;
  const reco::Candidate* daughter1;
  const reco::Candidate* daughter2;
  const reco::Candidate* mother1;
  //  const reco::Candidate* mother2;
  const reco::Candidate* muMinus;
  const reco::Candidate* muPlus;
  math::XYZTLorentzVectorD dimuon;
  double dimuonPx, dimuonPz, dimuonPt, pseudorapidity, Phi;
  
      for(auto &part : genParts){
	if((part.pdgId() == 1 || part.pdgId() == 2 || part.pdgId() == 3 || part.pdgId() == 4 || part.pdgId() == 5 || part.pdgId() == 6) && 
	   (abs(part.daughter(0)->pdgId()) == 11 || abs(part.daughter(0)->pdgId()) == 13)){
	  std::cout << "\nFound the quark! " << "\nQuark is: " << part.pdgId() << "\tStatus is: " << part.status() << "\tNumber of daughters are: " <<
	    part.numberOfDaughters() << "\tFirst daughter is:"  << part.daughter(0)->pdgId() << "\tSecond daughter is: " << part.daughter(1)->pdgId() << std::endl;
	  //      if(part.status() < -20 && part.status() > -30){ std::cout << "\nFound the Z boson!";
	  std::cout << "\nkinematic properties of the particles are: " << std::endl;
	  std::cout << "\npT1: " << part.daughter(0)->pt() << "\tpT2: " << part.daughter(1)->pt() << std::endl;
	  std::cout << "\neta1: " << part.daughter(0)->eta() << "\teta2: " << part.daughter(1)->eta() << std::endl;
	  std::cout << "\nphi1: " << part.daughter(0)->phi() << "\tphi2: " << part.daughter(1)->phi() << std::endl;
	  
	  daughter1 = getLastDaughter(part.daughter(0), part.daughter(0)->pdgId());
	  daughter2 = getLastDaughter(part.daughter(1), part.daughter(1)->pdgId());
	  std::cout << "\nDaughter particle is: " << daughter1->pdgId() << "tStatus is: " << daughter1->status()
		    << "\tDaughter2 is: " << daughter2->pdgId() << "\tStatus is: " << daughter2->status() << std::endl;
	  boson = nullptr;
	  if(!daughter1 || !daughter2){
	    std::cout<<"daughter1::0x"<<std::hex<<daughter1<<std::dec<<std::endl;
	    std::cout<<"daughter2::0x"<<std::hex<<daughter2<<std::dec<<std::endl;
	  }
	}

	if(part.pdgId() == 23 && (abs(part.daughter(0)->pdgId()) == 11 || abs(part.daughter(0)->pdgId()) == 13)){
	  std::cout << "\nFound the Z boson! " << "\tStatus is: " << part.status() << "\tNumber of daughters are: " <<
	    part.numberOfDaughters() << "\tFirst daughter is:"  << part.daughter(0)->pdgId() << "\tSecond daughter is: " << part.daughter(1)->pdgId() << std::endl;
	  //      if(part.status() < -20 && part.status() > -30){ std::cout << "\nFound the Z boson!";
	  std::cout << "\nkinematic properties of the particles are: " << std::endl;
	  std::cout << "\npT1: " << part.daughter(0)->pt() << "\tpT2: " << part.daughter(1)->pt() << std::endl;
	  std::cout << "\neta1: " << part.daughter(0)->eta() << "\teta2: " << part.daughter(1)->eta() << std::endl;
	  std::cout << "\nphi1: " << part.daughter(0)->phi() << "\tphi2: " << part.daughter(1)->phi() << std::endl;
	  
	  daughter1 = getLastDaughter(part.daughter(0), part.daughter(0)->pdgId());
	  daughter2 = getLastDaughter(part.daughter(1), part.daughter(1)->pdgId());
	  std::cout << "\nDaughter particle is: " << daughter1->pdgId() << "tStatus is: " << daughter1->status()
		    << "\tDaughter2 is: " << daughter2->pdgId() << "\tStatus is: " << daughter2->status() << std::endl;
	  mother1 = &part;
	  boson = mother1;
	  if(!boson || !daughter1 || !daughter2){
	    std::cout<<"boson::0x"<<std::hex<<boson<<std::dec<<std::endl;
	    std::cout<<"daughter1::0x"<<std::hex<<daughter1<<std::dec<<std::endl;
	    std::cout<<"daughter2::0x"<<std::hex<<daughter2<<std::dec<<std::endl;
	  }

	}
      }


    if(debug_){
      std::cout << "Eta of daughter1 is: " << daughter1->eta() << "\n";
      std::cout << "Eta of daughter2 is: " << daughter2->eta() << "\n";
    }

    //    if(!isCI_){
      if(boson){
	bosonId_=boson->pdgId();
	bosonP4_.fill(boson->p4());
	
	h_Zmass->Fill(boson->mass());
	h_Zpt->Fill(boson->pt());
	h_Zeta ->Fill(boson->eta());
	h_Zphi ->Fill(boson->phi());
	h_Zcharge->Fill(boson->charge());
      }
      //    }

    if(daughter1->charge() > 0 && daughter2->charge() < 0){
      muMinus = daughter2;
      muPlus = daughter1;
    }
    else if(daughter1->charge() < 0 && daughter2->charge() > 0){
      muMinus = daughter1;
      muPlus = daughter2;
    }

    else return;

    if(debug_){  
      std::cout<< "\n\nDaughter1: pId = " << muMinus->pdgId() << "\tpT = " << muMinus->pt() << "\teta = " 
	       << muMinus->eta() << "\tphi = " << muMinus->phi() << "\tq = " << muMinus->charge();
      std::cout<< "\nDaughter2: pId = " << muPlus->pdgId() << "\tpT = " << muPlus->pt() << "\teta = " << muPlus->eta() << "\tphi = " << 
	muPlus->phi() << "\tq = " << muPlus->charge();
    }
    
  muMinusP4_.fill(muMinus->p4());
    muMinusPID_=muMinus->pdgId();
    if(debug_){  
      std::cout<< "\n\nDaughter1: pId = " << muMinus->pdgId() << "\tpT = " << muMinus->pt() << "\teta = " 
	       << muMinus->eta() << "\tphi = " << muMinus->phi() << "\tq = " << muMinus->charge();
    std::cout<< "\nDaughter2: pId = " << muPlus->pdgId() << "\tpT = " << muPlus->pt() << "\teta = " << muPlus->eta() << "\tphi = " << 
      muPlus->phi() << "\tq = " << muPlus->charge();
    }


    h_muMinusmass->Fill(muMinus->mass());
    h_muMinuspt->Fill(muMinus->pt());
    h_muMinuseta->Fill(muMinus->eta());
    h_muMinusphi->Fill(muMinus->phi());
    h_muMinuscharge->Fill(muMinus->charge());
    h_thetaMuMinus->Fill(muMinus->theta());  
  
    muPlusP4_.fill(muPlus->p4());
    muPlusPID_=muPlus->pdgId();

    h_muPlusmass->Fill(muPlus->mass());
    h_muPluspt->Fill(muPlus->pt());
    h_muPluseta->Fill(muPlus->eta());
    h_muPlusphi->Fill(muPlus->phi());
    h_muPluscharge->Fill(muPlus->charge());
    h_thetaMuPlus->Fill(muPlus->theta());    

    dimuon = muMinus->p4() + muPlus->p4();

    dimuonPt =dimuon.pt();
    dimuonPz = dimuon.pz();
    pseudorapidity = asinh(dimuonPz/dimuonPt);
    dimuonPx = dimuon.px();
    Phi = acos(dimuonPx/dimuonPt);

    h_dphi->Fill(TVector2::Phi_mpi_pi(muMinus->phi()- muPlus->phi()));
    h_dtheta->Fill(TVector2::Phi_mpi_pi(muMinus->theta()- muPlus->theta()));
    h_dr->Fill(reco::deltaR(muMinus->p4(),muPlus->p4()));
    h_massInvar->Fill(sqrt(2 * daughter1->pt() * daughter2->pt() *( cosh(daughter1->eta() - daughter2->eta()) - cos(TVector2::Phi_mpi_pi(daughter1->phi() - daughter2->phi())))));
    h_dimuonPt->Fill(dimuonPt);
    h_dimuonEta->Fill(pseudorapidity);
    h_dimuonPhi->Fill(Phi);
    h2_phi1_vs_phi2->Fill(muMinus->phi(),muPlus->phi());  
    h2_eta1_vs_eta2->Fill(muMinus->eta(),muPlus->eta());
    h2_pt1_vs_pt2->Fill(muMinus->pt(),muPlus->pt());


    //  }

  std::cout << "\n\n===========================================================================================================" << std::endl;
  tree_->Fill();  
}

bool Dimuon::isBoson(int pid)
{
  if(pid==23 || abs(pid)==22 || pid==32){
    if(debug_) std::cout << "\n\nFound Boson\n";
    return true;
  }
  else return false;
}

bool Dimuon::isMuon(int pid){
  if(abs(pid)==11 || abs(pid) ==13){
    if(debug_) std::cout << "\n\nFound A Muon!\n";
    return true;
  }
  else return false;
}

bool Dimuon::checkBosonStatus( const reco::GenParticleCollection& genParts){
  const reco::Candidate* boson = getBoson(genParts);
  if(boson == nullptr){
    if(debug_) std::cout << "\nBoson is: "  << boson;
    return false;
  }

  else if( boson->status() != 22){
    if(debug_)  std::cout <<"\nBoson Status is: "<< boson->status();
    return false;
  }
 
    return true;
 }

const reco::Candidate* Dimuon::getBoson( const reco::GenParticleCollection& genParts)
{
  for(auto &part : genParts){
    if(isBoson(part.pdgId())){
      if(debug_){
      std::cout << "\npId is: " << part.pdgId();
      std::cout << "\nStatus is: " << part.status();
      }
      return getLastDaughter(&part,part.pdgId());
    }
  }
  return nullptr;
}


const reco::Candidate* Dimuon::getLastDaughter(const reco::Candidate* part,int pid)
{
  for(size_t partNr =0; part && partNr<part->numberOfDaughters();partNr++){
    if(part->daughter(partNr)->pdgId()==pid) return getLastDaughter(part->daughter(partNr),pid);
  }
  return part;
}
       
const reco::Candidate* Dimuon::getDaughter(const reco::Candidate* part,int pid)
{  
  for(size_t partNr =0; part && partNr<part->numberOfDaughters();partNr++){
    if(part->daughter(partNr)->pdgId()==pid) return part->daughter(partNr);
  }
  return nullptr;
}

 const reco::Candidate* Dimuon::getMother(const reco::Candidate* part, int pid)
{
  for(size_t partNr = 0; part && partNr < part->numberOfMothers(); partNr++){
    if(part->mother(partNr)->pdgId() == pid) return getMother(part->mother(partNr),pid);
  
    else if(abs(part->mother(partNr)->pdgId()) == 1 || abs(part->mother(partNr)->pdgId()) == 2 ||
	    abs(part->mother(partNr)->pdgId()) == 3 || abs(part->mother(partNr)->pdgId()) == 4 ||
		 abs(part->mother(partNr)->pdgId()) == 5 || abs(part->mother(partNr)->pdgId()) == 6 ||
	    abs(part->mother(partNr)->pdgId()) == 7 || abs(part->mother(partNr)->pdgId()) == 8 ||
	    abs(part->mother(partNr)->pdgId()) == 23 || abs(part->mother(partNr)->pdgId()) == 32  || 
	    abs(part->mother(partNr)->pdgId()) == 22) return part->mother(partNr);
  }  
  return nullptr;
  
}


// ------------ method called once each job just before starting event loop  ------------

// ------------ method called once each job just after ending the event loop  ------------
void 
Dimuon::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
  
void 
Dimuon::beginRun(edm::Run const& iRun, edm::EventSetup const& iEventSetup)
{
  edm::Handle< GenRunInfoProduct > genInfoProduct;
  iRun.getByLabel("generator", genInfoProduct );
  crossSec = genInfoProduct->internalXSec().value();

  std::cout<< "Cross Section is: "  << crossSec << std::endl;  

}
  

// ------------ method called when ending the processing of a run  ------------
/*
  void 
  Dimuon::endRun(edm::Run const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
  void 
  Dimuon::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
  void 
  Dimuon::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Dimuon::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Dimuon);
