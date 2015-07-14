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


#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

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

 static const reco::Candidate* getDaughter(const reco::Candidate* part,int pid);
  static const reco::Candidate* getLastDaughter(const reco::Candidate* part,int pid);
  static const reco::Candidate* getBoson( const reco::GenParticleCollection& genParts);
  static bool isBoson(int pid);

  //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  edm::InputTag genPartsTag;

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
  TH1F *h_mu1mass,*h_mu1pt,*h_mu1eta,*h_mu1phi,*h_mu1charge;
  TH1F *h_mu2mass,*h_mu2pt,*h_mu2eta,*h_mu2phi,*h_mu2charge;
  TH1F *h_dphi,*h_dtheta,*h_thetaMu1,*h_thetaMu2;

  TH2F *h2_pt1_vs_pt2,*h2_eta1_vs_eta2,*h2_phi1_vs_phi2;

  P4Struct bosonP4_; // as a sanity check we have the right event...
  P4Struct mu1P4_;
  P4Struct mu2P4_;
  int mu1PID_;
  int mu2PID_;
  int bosonId_;


  edm::InputTag genPartsTag_;
  int decayParticlePID_;


  // ----------member data ---------------------------
  
};

  void Dimuon::beginJob()
  {
    edm::Service<TFileService> fs;
    h_Zmass = fs->make<TH1F>("Zmass" , "m", 1000, 0., 500);
    h_Zpt  = fs->make<TH1F>( "Zpt"  , "p_{t}", 500,  0., 1000. );
    h_Zeta = fs->make<TH1F>( "Zeta" , "#eta" , 100, -10., 10.    );
    h_Zphi = fs->make<TH1F>( "Zphi" , "#phi" , 100,  -3.20, 3.20   );
    h_Zcharge = fs->make<TH1F>( "Zcharge" , "Q" ,3,  -1.5, 1.5    );
    h_mu1mass = fs->make<TH1F>("mu1mass" , "m", 1000, 0., 500);
    h_mu1pt  = fs->make<TH1F>( "mu1pt"  , "p_{t}", 500,  0., 1000. );
    h_mu1eta = fs->make<TH1F>( "mu1eta" , "#eta" , 100, -5., 5.    );
    h_mu1phi = fs->make<TH1F>( "mu1phi" , "#phi" , 100,  -3.15, 3.15   );
    h_mu1charge = fs->make<TH1F>( "mu1charge" , "Q" ,3,  -1.5, 1.5    );

    h_mu2mass = fs->make<TH1F>("mu2mass" , "m", 1000, 0., 500);
    h_mu2pt  = fs->make<TH1F>( "mu2pt"  , "p_{t}", 500,  0., 1000. );
    h_mu2eta = fs->make<TH1F>( "mu2eta" , "#eta" , 100, -5., 5.    );
    h_mu2phi = fs->make<TH1F>( "mu2phi" , "#phi" , 100,  -3.15, 3.15   );
    h_mu2charge = fs->make<TH1F>( "mu2charge" , "Q" ,3,  -1.5, 1.5    );

    h_dphi = fs->make<TH1F>("delta phi", "#delta #phi", 100, -3.15, 3.15 );       
    h_dtheta = fs->make<TH1F>("delta theta", "#delta #theta", 100, -3.15, 3.15);  
    h_thetaMu1 = fs->make<TH1F>("theta muon_1", "#theta", 100, -3.15, 3.15);      
    h_thetaMu2 = fs->make<TH1F>("theta muon_2", "#theta", 100, -3.15, 3.15);      

    h2_pt1_vs_pt2   = fs->make<TH2F>( "pt1_vs_pt2"   , "p_{t,1} vs. p_{t,2}"   , 500,  0., 2500., 500,  0., 2500.);
    h2_eta1_vs_eta2 = fs->make<TH2F>( "eta1_vs_eta2" , "#eta_{1} vs. #eta_{2}" , 100, -5., 5.   , 100, -5., 5.   );
    h2_phi1_vs_phi2 = fs->make<TH2F>( "phi1_vs_phi2" , "#phi_{1} vs. #phi_{2}" , 100,  -3.15, 3.15  , 100,  -3.15, 3.15  );

    tree_= fs->make<TTree>("pdfTree","PDF Tree");
    // tree_->Branch("evtId",&evtId_,EventId::contents().c_str());
    tree_->Branch("bosonP4",&bosonP4_,P4Struct::contents().c_str());
    tree_->Branch("decay1P4",&mu1P4_,P4Struct::contents().c_str());
    tree_->Branch("decay2P4",&mu2P4_,P4Struct::contents().c_str());
    tree_->Branch("decay1PID",&mu1PID_,"decay1PID/I");
    tree_->Branch("decay2PID",&mu2PID_,"decay2PID/I");
    tree_->Branch("bosonPID",&bosonId_,"bosonPID/I");
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
    genPartsTag_=iConfig.getParameter<edm::InputTag>("genPartsTag");
    decayParticlePID_ = iConfig.getParameter<int>("decayParticlePID");
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
    using namespace edm;
 edm::Handle<reco::GenParticleCollection> genPartsHandle;
  iEvent.getByLabel(genPartsTag_,genPartsHandle);
  const reco::GenParticleCollection& genParts = *genPartsHandle;


 bosonId_=0;
  bosonP4_.clear();
  mu1P4_.clear();
  mu2P4_.clear();
  mu1PID_=0;
  mu2PID_=0;



  const reco::Candidate* boson = getBoson(genParts);
  const reco::Candidate* daughter1 = getDaughter(boson,decayParticlePID_);
  const reco::Candidate* daughter2 = getDaughter(boson,decayParticlePID_*-1); 

  if(boson){
    bosonId_=boson->pdgId();
    bosonP4_.fill(boson->p4());

    h_Zmass->Fill(boson->mass());
    h_Zpt->Fill(boson->pt());
    h_Zeta ->Fill(boson->eta());
    h_Zphi ->Fill(boson->phi());
    h_Zcharge->Fill(boson->charge());
  }

  if(daughter1){
    mu1P4_.fill(daughter1->p4());
    mu1PID_=daughter1->pdgId();

    h_mu1mass->Fill(daughter1->mass());
    h_mu1pt->Fill(daughter1->pt());
    h_mu1eta->Fill(daughter1->eta());
    h_mu1phi->Fill(daughter1->phi());
    h_mu1charge->Fill(daughter1->charge());
    h_thetaMu1->Fill(daughter1->theta());  
  }

  if(daughter2){
    mu2P4_.fill(daughter2->p4());
    mu2PID_=daughter2->pdgId();

    h_mu2mass->Fill(daughter2->mass());
    h_mu2pt->Fill(daughter2->pt());
    h_mu2eta->Fill(daughter2->eta());
    h_mu2phi->Fill(daughter2->phi());
    h_mu2charge->Fill(daughter2->charge());
    h_thetaMu2->Fill(daughter2->theta());    
}

  h_dphi->Fill(TVector2::Phi_mpi_pi(daughter1->phi()- daughter2->phi()));
  h_dtheta->Fill(TVector2::Phi_mpi_pi(daughter1->theta()- daughter2->theta()));

  h2_phi1_vs_phi2->Fill(daughter1->phi(),daughter2->phi());  
  h2_eta1_vs_eta2->Fill(daughter1->eta(),daughter2->eta());
  h2_pt1_vs_pt2->Fill(daughter1->pt(),daughter2->pt());

  tree_->Fill();  





#ifdef THIS_IS_AN_EVENT_EXAMPLE
    Handle<ExampleData> pIn;
    iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
    ESHandle<SetupData> pSetup;
    iSetup.get<SetupRecord>().get(pSetup);
#endif
  }


bool Dimuon::isBoson(int pid)
{
  if(pid==23 || abs(pid)==22 || pid==32) return true;
  else return false;
}

const reco::Candidate* Dimuon::getBoson( const reco::GenParticleCollection& genParts)
{
  for(auto &part : genParts){
    if(isBoson(part.pdgId())){
      return getLastDaughter(&part,part.pdgId());
    }
  }
  return nullptr;
}
 
const reco::Candidate* Dimuon::getLastDaughter(const reco::Candidate* part,int pid)
{
 std::cout << "daughters" << std::endl;
  for(size_t partNr =0; part && partNr<part->numberOfDaughters();partNr++){
    std::cout << pid<< "  "   << partNr << "\t" << part->daughter(partNr)->pdgId() << std::endl;
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

  // ------------ method called once each job just before starting event loop  ------------

  // ------------ method called once each job just after ending the event loop  ------------
  void 
  Dimuon::endJob() 
  {
  }

  // ------------ method called when starting to processes a run  ------------
  /*
    void 
    Dimuon::beginRun(edm::Run const&, edm::EventSetup const&)
    {
    }
  */

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
