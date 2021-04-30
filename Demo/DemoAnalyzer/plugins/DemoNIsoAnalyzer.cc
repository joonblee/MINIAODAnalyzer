// -*- C++ -*-
//
// Package:    Demo/DemoNIsoAnalyzer
// Class:      DemoNIsoAnalyzer
// 
/**\class DemoNIsoAnalyzer DemoNIsoAnalyzer.cc Demo/DemoNIsoAnalyzer/plugins/DemoNIsoAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Jhovanny Andres Mejia Guisao
//         Created:  Tue, 05 Sep 2017 22:37:14 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"


#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"


#include <TH2D.h>
#include <TMath.h>
#include <TLorentzVector.h>


/////////////////////
// -- For Muons -- //
/////////////////////
#include "DataFormats/MuonReco/interface/MuonFwd.h" // -- Forward declarations of muon variables
#include "DataFormats/MuonReco/interface/Muon.h" // -- A reconstructed Muon. (tracker alone, muon detector alone, combined muon plus tracker)
#include "DataFormats/PatCandidates/interface/Muon.h" // -- Analysis-level muon class, pat::Muon implements the analysis-level muon class within the 'pat' namespace. PAT(Physics Analysis Toolkit)
#include "DataFormats/MuonReco/interface/MuonTimeExtra.h"
#include "DataFormats/MuonReco/interface/MuonTimeExtraMap.h"
#include "DataFormats/MuonReco/interface/MuonCocktails.h"
#include "DataFormats/MuonReco/interface/MuonCosmicCompatibility.h" // -- Variables shows how muon is cosmic-like

///////////////////////////////////
// -- For Electrons & Photons -- //
///////////////////////////////////
#include "DataFormats/PatCandidates/interface/Photon.h" // -- Analysis-level Photon class, pat::Photon implements the analysis-level photon class within the 'pat' namespace
#include "DataFormats/PatCandidates/interface/Electron.h" // -- Analysis-level electron class,  pat::Electron implements the analysis-level electron class within the 'pat' namespace
#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h" // -- It seems as a function refer effective areas of ECAL
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"


///////////////////
// -- For MET -- //
///////////////////
#include "DataFormats/PatCandidates/interface/MET.h" // -- Analysis-level MET class, pat::MET implements an analysis-level missing energy class as a 4-vector within the 'pat' namespace
#include "DataFormats/METReco/interface/PFMETFwd.h"
#include "DataFormats/METReco/interface/PFMET.h" // --  MET made from Particle Flow candidates, with various fraction functions ex)photonEtFraction()
#include "DataFormats/METReco/interface/PFMETCollection.h"

//////////////////////////
// -- Track & Vertex -- //
//////////////////////////
#include "DataFormats/TrackReco/interface/TrackFwd.h" // -- Forward definitions for tracker variables
#include "DataFormats/TrackReco/interface/Track.h" // -- reconstructed tracks that are stored in the AOD and RECO. also contains a reference to more detailed information(RECO)
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h" // -- A composite Candidate  with error matrix and other vertex fix information.
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h" // -- Foward declaration of VertexCompositeCandidate class's variables
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h" // -- Least-squares vertex fitter implemented in the Kalman Filter formalism
#include "PhysicsTools/RecoUtils/interface/CandCommonVertexFitter.h" // --
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h" // -- Helper class to build TransientTrack from the persistent Track.
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "DataFormats/PatCandidates/interface/Vertexing.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexTools/interface/InvariantMassFromVertex.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

////////////////////
// -- For Jets -- //
////////////////////
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/PatCandidates/interface/Jet.h" // -- Analysis-level calorimeter jet class, Jet implements the analysis-level calorimeter jet class within the 'pat' namespace.
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"

////////////////////////////////
// -- For PackedCandidates -- //
////////////////////////////////
#include "DataFormats/PatCandidates/interface/PackedCandidate.h" // muy importante para MiniAOD

//////////////////////////
// -- Track & Vertex -- //
//////////////////////////
#include "DataFormats/TrackReco/interface/TrackFwd.h" // -- Forward definitions for tracker variables
#include "DataFormats/TrackReco/interface/Track.h" // -- reconstructed tracks that are stored in the AOD and RECO. also contains a reference to more detailed information(RECO)
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h" // -- A composite Candidate  with error matrix and other vertex fix information.
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h" // -- Foward declaration of VertexCompositeCandidate class's variables
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h" // -- Least-squares vertex fitter implemented in the Kalman Filter formalism
#include "PhysicsTools/RecoUtils/interface/CandCommonVertexFitter.h" // --
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h" // -- Helper class to build TransientTrack from the persistent Track.
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "DataFormats/PatCandidates/interface/Vertexing.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexTools/interface/InvariantMassFromVertex.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"


// -- calculations -- //
#include "DataFormats/Math/interface/deltaR.h"

//
// class declaration
//
class DemoNIsoAnalyzer : public edm::EDAnalyzer {
   public:
      explicit DemoNIsoAnalyzer(const edm::ParameterSet&);
      ~DemoNIsoAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;


      // ----------member data ---------------------------
      //edm::EDGetTokenT<edm::View<pat::PackedCandidate>> trackTags_;
      edm::EDGetTokenT<pat::PackedCandidateCollection> pcToken;
      edm::EDGetTokenT<std::vector<pat::Jet>> JetToken;
      edm::EDGetTokenT<std::vector<pat::Muon>> MuonToken;
      edm::EDGetTokenT<edm::View<pat::Electron>> ElectronToken;
      edm::EDGetTokenT<edm::View<pat::Photon>> PhotonToken;
      edm::EDGetTokenT<std::vector<pat::MET>> MetToken;
      edm::EDGetTokenT<reco::VertexCollection> pvToken;
      //edm::EDGetTokenT< edm::View<reco::Track> > TrackToken;

      TH1D *jetPt, *jetEta, *jetPhi;
      TH1D *muonPt, *muonEta, *muonPhi;
      TH1D *metPt, *metPhi, *metUncorPt, *metUncorPhi;
      TH1D *nJets;
      TH1D *MASS;
      TH1D *dR_JetMu, *dR_MuMu, *dR_JetZ, *ZPt;

      // ------ functions ----- //
      std::vector<pat::Muon> W_selection(edm::Handle<std::vector<pat::Muon>>, edm::Handle<std::vector<pat::MET>>, edm::Handle<reco::VertexCollection>);
      std::vector<pat::Muon> Z_selection(edm::Handle<std::vector<pat::Muon>>, edm::Handle<reco::VertexCollection>);
      std::vector<pat::Jet> jet_selection(edm::Handle<std::vector<pat::Jet>>, double);
      std::vector<pat::Muon> muon_selection(edm::Handle<std::vector<pat::Muon>>, edm::Handle<reco::VertexCollection>, double);
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
DemoNIsoAnalyzer::DemoNIsoAnalyzer(const edm::ParameterSet& iConfig) :  

//trackTags_(consumes<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("tracks")))
pcToken       ( consumes< pat::PackedCandidateCollection > (iConfig.getUntrackedParameter<edm::InputTag>("pfCands"))   ),
JetToken      ( consumes< std::vector<pat::Jet> >          (iConfig.getUntrackedParameter<edm::InputTag>("Jet")) ),
MuonToken     ( consumes< std::vector<pat::Muon> >         (iConfig.getUntrackedParameter<edm::InputTag>("Muon")) ),
ElectronToken ( consumes< edm::View<pat::Electron> >       (iConfig.getUntrackedParameter<edm::InputTag>("Electron")) ),
PhotonToken   ( consumes< edm::View<pat::Photon> >         (iConfig.getUntrackedParameter<edm::InputTag>("Photon")) ),
MetToken      ( consumes< std::vector<pat::MET> >          (iConfig.getUntrackedParameter<edm::InputTag>("MET")) ),
pvToken       ( consumes< reco::VertexCollection >         (iConfig.getUntrackedParameter<edm::InputTag>("PrimaryVertex")) )
{
  //now do what ever initialization is needed
  edm::Service<TFileService> fs;
  jetPt               = fs->make<TH1D>("jet pt" , "jet pt" , 200 , 0 , 200 );
  jetEta              = fs->make<TH1D>("jet eta" , "jet eta" , 100 , -5. , 5. );
  jetPhi              = fs->make<TH1D>("jet phi" , "jet phi" , 64 , -3.2 , 3.2 );
  muonPt              = fs->make<TH1D>("muon pt" , "muon pt" , 200 , 0 , 200 );
  muonEta             = fs->make<TH1D>("muon eta" , "muon eta" , 100 , -5. , 5. );
  muonPhi             = fs->make<TH1D>("muon phi" , "muon phi" , 64 , -3.2 , 3.2 );
  metUncorPt          = fs->make<TH1D>("met uncor pt" , "met uncor pt" , 200 , 0 , 200 );
  metUncorPhi         = fs->make<TH1D>("met uncor phi" , "met uncor phi" , 64 , -3.2 , 3.2 );
  metPt               = fs->make<TH1D>("met pt" , "met pt" , 200 , 0 , 200 );
  metPhi              = fs->make<TH1D>("met phi" , "met phi" , 64 , -3.2 , 3.2 );

  nJets               = fs->make<TH1D>("nJets", "nJets", 10,0,10);

  MASS                = fs->make<TH1D>("MASS","MASS", 2000, 0, 200);
  dR_JetMu            = fs->make<TH1D>("dR_JetMu","dR_JetMu", 1000, 0., 1.);
  dR_JetZ             = fs->make<TH1D>("dR_JetZ","dR_JetZ",1000,0.,1.);
  dR_MuMu             = fs->make<TH1D>("dR_MuMu","dR_MuMu",1000,0.,1.);
  ZPt                 = fs->make<TH1D>("ZPt","ZPt",1000,0,1000);
}


DemoNIsoAnalyzer::~DemoNIsoAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void DemoNIsoAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  using namespace edm;
  using namespace std;
    
  edm::Handle< pat::PackedCandidateCollection > pcHandle;
  iEvent.getByToken(pcToken, pcHandle);
  edm::Handle< std::vector<pat::Jet> > JetHandle;
  iEvent.getByToken(JetToken, JetHandle);
  edm::Handle<std::vector<pat::Muon>> MuonHandle;
  iEvent.getByToken(MuonToken, MuonHandle);
  edm::Handle<edm::View<pat::Electron>> ElectronHandle;
  iEvent.getByToken(ElectronToken, ElectronHandle);
  edm::Handle<edm::View<pat::Photon>> PhotonHandle;
  iEvent.getByToken(PhotonToken, PhotonHandle);
  edm::Handle<std::vector<pat::MET>> MetHandle;
  iEvent.getByToken(MetToken, MetHandle);
  edm::Handle<reco::VertexCollection> pvHandle;
  iEvent.getByToken(pvToken, pvHandle);



  std::cout<<std::endl;
  std::cout<<"Event "<<std::endl;
  cout<<"[Event id] run # = "<<iEvent.id().run()<<", evt # = "<<iEvent.id().event()<<endl;
  std::cout<<"nMuon = "<<MuonHandle->size()<<"  |  nJet = "<<JetHandle->size()<<std::endl;


  if( MuonHandle->size() < 2 ) return;
  if( JetHandle->size() == 0 ) return;


  //vector<pat::Muon> muons = Z_selection(MuonHandle, pvHandle);
  double SubMuonPtCut = 13.; double LeadMuonPtCut = 32.;
  vector<pat::Muon> muons = muon_selection(MuonHandle, pvHandle, SubMuonPtCut);
  if( muons.size() == 0 ) return;

  double jetPtCut = 70.;
  vector<pat::Jet> jets = jet_selection(JetHandle, jetPtCut); 
  if( jets.size() == 0 ) return;
 
  vector<pat::Muon> niso;
  for(unsigned iMuon = 0; iMuon < muons.size(); iMuon++) {
    const pat::Muon muon = muons.at(iMuon);
    for( auto jet : jets ) {
      if( deltaR( jet, muon ) < 0.3 ) { niso.push_back(muon); break; }
    }
  }
  if( niso.size() < 2 ) return;


  vector<pat::Muon> dimuon;
  pat::Jet this_jet;
  for( auto jet : jets ) {
    for(unsigned iMuon = 0; iMuon < niso.size(); iMuon++) {
      const pat::Muon muon_1 = niso.at(iMuon);
      if( !(muon_1.pt() > LeadMuonPtCut) ) break;
      if( !( deltaR( jet, muon_1 ) < 0.3 ) ) continue;

      for(unsigned jMuon = iMuon+1; jMuon < niso.size(); jMuon++) {
        const pat::Muon muon_2 = niso.at(jMuon);
        if( !( deltaR( jet, muon_2 ) < 0.3 ) ) continue;

        dimuon.push_back( muon_1 );
        dimuon.push_back( muon_2 );
        this_jet = jet;
        break;
      }
      if( dimuon.size() == 2 ) break;
    }
    if( dimuon.size() == 2 ) break;
  }
  if( dimuon.size() != 2 ) return;


  for( auto muon : muons ) { muonPt->Fill(muon.pt()); muonEta->Fill(muon.eta()); muonPhi->Fill(muon.phi()); }
  for( auto jet : jets ) { jetPt->Fill(jet.pt()); jetEta->Fill(jet.eta()); jetPhi->Fill(jet.phi()); }
  nJets->Fill(jets.size());

  const pat::MET  met  = MetHandle->front();
  metUncorPt->Fill(met.pt()); metUncorPhi->Fill(met.phi());
  metPt->Fill(met.corPt(pat::MET::Type1XY)); metPhi->Fill(met.corPhi(pat::MET::Type1XY));


  TLorentzVector muon_1_4vec, muon_2_4vec, jet_4vec;
  muon_1_4vec.SetPtEtaPhiM(dimuon[0].pt(), dimuon[0].eta(), dimuon[0].phi(), dimuon[0].mass());
  muon_2_4vec.SetPtEtaPhiM(dimuon[1].pt(), dimuon[1].eta(), dimuon[1].phi(), dimuon[1].mass());
  jet_4vec.SetPtEtaPhiE(this_jet.pt(), this_jet.eta(), this_jet.phi(), this_jet.energy());
  double mass = (muon_1_4vec + muon_2_4vec).M();

  MASS->Fill(mass);

  dR_JetMu->Fill( jet_4vec.DeltaR( muon_1_4vec ) );
  dR_JetMu->Fill( jet_4vec.DeltaR( muon_2_4vec ) );
  dR_JetZ->Fill( jet_4vec.DeltaR( muon_1_4vec + muon_2_4vec ) );
  dR_MuMu->Fill( muon_1_4vec.DeltaR( muon_2_4vec ) );
  ZPt->Fill( (muon_1_4vec+muon_2_4vec).Pt() );


  /*
  cout<<endl;
  cout<<"============================================================="<<endl;
  cout<<"[Event id] run # = "<<iEvent.id().run()<<", evt # = "<<iEvent.id().event()<<endl;
  cout<<"[Jet Info] #(AllJets) = "<<JetHandle->size()<<", #(SelectedJets) = "<<jets.size()<<endl;
  */


#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
   
}

std::vector<pat::Muon> DemoNIsoAnalyzer::Z_selection(edm::Handle<std::vector<pat::Muon>> MuonHandle, edm::Handle<reco::VertexCollection> pvHandle)
{
  std::vector<pat::Muon> muons;

  const reco::Vertex &vtx = pvHandle->front();
  for(unsigned iMuon = 0; iMuon < MuonHandle->size(); iMuon++) {
    const pat::Muon muon_1 = MuonHandle->at(iMuon);

    if( !(muon_1.pt() > 32) ) return muons;
    if( !(fabs(muon_1.eta()) < 2.4) ) continue;
    //if( !(muon_1.isolationR03().sumPt/muon_1.pt() < 0.1) ) continue;
    if( !(muon_1.isTightMuon(vtx)) ) continue;

    for(unsigned jMuon = iMuon+1; jMuon < MuonHandle->size(); jMuon++) {
      const pat::Muon muon_2 = MuonHandle->at(jMuon);

      if( !(muon_2.pt() > 13) ) break;
      if( !(fabs(muon_2.eta()) < 2.4) ) continue;
      //if( !(muon_2.isolationR03().sumPt/muon_2.pt() < 0.1) ) continue;
      if( !(muon_2.isTightMuon(vtx)) ) continue;

      TLorentzVector muon_1_4vec, muon_2_4vec;
      muon_1_4vec.SetPtEtaPhiM(muon_1.pt(), muon_1.eta(), muon_1.phi(), muon_1.mass());
      muon_2_4vec.SetPtEtaPhiM(muon_2.pt(), muon_2.eta(), muon_2.phi(), muon_2.mass());
      double mass = (muon_1_4vec + muon_2_4vec).M();
      //if( 15 < mass && mass < 25 ) { 
        MASS->Fill(mass);
        muons.push_back(muon_1); muons.push_back(muon_2); 
        return muons;
      //}
    }
  }

  return muons;
}

std::vector<pat::Muon> DemoNIsoAnalyzer::muon_selection(edm::Handle<std::vector<pat::Muon>> MuonHandle, edm::Handle<reco::VertexCollection> pvHandle, double muonPt)
{
  std::vector<pat::Muon> muons;

  const reco::Vertex &vtx = pvHandle->front();
  for(unsigned iMuon = 0; iMuon < MuonHandle->size(); iMuon++) {
    const pat::Muon muon = MuonHandle->at(iMuon);

    if( !(muon.pt() > muonPt) ) return muons;
    if( !(fabs(muon.eta()) < 2.4) ) continue;
    //if( !(muon.isolationR03().sumPt/muon.pt() < 0.1) ) continue;
    if( !(muon.isTightMuon(vtx)) ) continue;

    TLorentzVector muon_4vec;
    muon_4vec.SetPtEtaPhiM(muon.pt(), muon.eta(), muon.phi(), muon.mass());
    muons.push_back(muon); 
  }

  return muons;
}


std::vector<pat::Jet> DemoNIsoAnalyzer::jet_selection(edm::Handle<std::vector<pat::Jet>> JetHandle, double jetPt)
{
  std::vector<pat::Jet> jets;
  for(unsigned iJet = 0; iJet < JetHandle->size(); iJet++) {
    const pat::Jet jet = JetHandle->at(iJet);
    if( !(jetPt < jet.pt()) ) continue;
    if( !(fabs(jet.eta()) < 2.4) ) continue;

    double NHF      = jet.neutralHadronEnergyFraction();
    double NEMF     = jet.neutralEmEnergyFraction();
    double NumConst = jet.chargedMultiplicity()+jet.neutralMultiplicity();
    double MUF      = jet.muonEnergyFraction();
    double CHF      = jet.chargedHadronEnergyFraction();
    double CHM      = jet.chargedMultiplicity();
    double CEMF     = jet.chargedEmEnergyFraction();
    double NumNeutralParticles = jet.neutralMultiplicity();
    double eta      = jet.eta();

    bool tightJetID = (NHF<0.90 && NEMF<0.90 && NumConst>1) && ((fabs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF < 0.99) || fabs(eta)>2.4) && fabs(eta)<=2.7;

    if( !tightJetID ) continue;

    jets.push_back( jet );

  }
  return jets;
}


// ------------ method called once each job just before starting event loop  ------------
void DemoNIsoAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void DemoNIsoAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
DemoNIsoAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
DemoNIsoAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
DemoNIsoAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
DemoNIsoAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DemoNIsoAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DemoNIsoAnalyzer);
