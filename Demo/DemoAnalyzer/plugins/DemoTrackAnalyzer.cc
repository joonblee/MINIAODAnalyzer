// -*- C++ -*-
//
// Package:    Demo/DemoTrackAnalyzer
// Class:      DemoTrackAnalyzer
// 
/**\class DemoTrackAnalyzer DemoTrackAnalyzer.cc Demo/DemoTrackAnalyzer/plugins/DemoTrackAnalyzer.cc

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
class DemoTrackAnalyzer : public edm::EDAnalyzer {
   public:
      explicit DemoTrackAnalyzer(const edm::ParameterSet&);
      ~DemoTrackAnalyzer();

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
      TH1D *nTrks, *nGoodTrks, *nGoodPVTrks, *trkPt, *goodTrkPt, *goodPVTrkPt, *trkPtSumRatio, *goodTrkPtSumRatio, *goodPVTrkPtSumRatio, *pcPtRatio;
      TH1D *nGoodPVTrks_pt30, *nGoodPVTrks_pt100, *goodPVTrkPt_pt30, *goodPVTrkPt_pt100;
      TH1D *nJets;
      TH1D *MASS, *M_T;
      TH2D *nTrksPerJetPt, *nGoodTrksPerJetPt, *nGoodPVTrksPerJetPt;


      std::vector< std::vector<int> > vec_Ntrks, vec_Ngoodtrks, vec_NgoodPVtrks;
      int MaxJetPt;


      // ------ functions ----- //
      std::vector<pat::Muon> W_selection(edm::Handle<std::vector<pat::Muon>>, edm::Handle<std::vector<pat::MET>>, edm::Handle<reco::VertexCollection>);
      std::vector<pat::Muon> Z_selection(edm::Handle<std::vector<pat::Muon>>, edm::Handle<reco::VertexCollection>);
      std::vector<pat::Jet> jet_selection(edm::Handle<std::vector<pat::Jet>>, std::vector<pat::Muon>);
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
DemoTrackAnalyzer::DemoTrackAnalyzer(const edm::ParameterSet& iConfig) :  

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

  nTrks               = fs->make<TH1D>("ntracks" , "ntracks" , 30 , 0 , 30 );
  nGoodTrks           = fs->make<TH1D>("ntracks_GoodTrk" , "ntracks_GoodTrk" , 30 , 0 , 30 );
  nGoodPVTrks         = fs->make<TH1D>("ntracks_GoodPVTrk" , "ntracks_GoodPVTrk" , 30 , 0 , 30 );

  nGoodPVTrks_pt30    = fs->make<TH1D>("ntracks_GoodPVTrk_pt30" , "ntracks_GoodPVTrk_pt30" , 30 , 0 , 30 );
  nGoodPVTrks_pt100   = fs->make<TH1D>("ntracks_GoodPVTrk_pt100" , "ntracks_GoodPVTrk_pt100" , 30 , 0 , 30 );


  trkPt               = fs->make<TH1D>("track_pt" , "track_pt" , 8000 , 0 , 80 );
  goodTrkPt           = fs->make<TH1D>("track_pt_GoodTrk" , "track_pt_GoodTrk" , 8000 , 0 , 80 );
  goodPVTrkPt         = fs->make<TH1D>("track_pt_GoodPVTrk" , "track_pt_GoodPVTrk" , 80000 , 0 , 80 );

  goodPVTrkPt_pt30    = fs->make<TH1D>("track_pt_GoodPVTrk_pt30" , "track_pt_GoodPVTrk_pt30" , 8000 , 0 , 80 );
  goodPVTrkPt_pt100   = fs->make<TH1D>("track_pt_GoodPVTrk_pt100" , "track_pt_GoodPVTrk_pt100" , 8000 , 0 , 80 );


  trkPtSumRatio       = fs->make<TH1D>("trkPtSumRatio" , "trkPtSumRatio" , 100 , 0 , 1 );
  goodTrkPtSumRatio   = fs->make<TH1D>("trkPtSumRatio_GoodTrk" , "trkPtSumRatio_GoodTrk" , 100 , 0 , 1 );
  goodPVTrkPtSumRatio = fs->make<TH1D>("trkPtSumRatio_GoodPVTrk" , "trkPtSumRatio_GoodPVTrk" , 100 , 0 , 1 );
  pcPtRatio           = fs->make<TH1D>("pcPtRatio", "pcPtRatio", 100, 0.5, 1.5);
  nTrksPerJetPt       = fs->make<TH2D>("nTrksPerJetPt", "ntracks/JetPt", 200, 0, 200, 30, 0, 30);
  nGoodTrksPerJetPt   = fs->make<TH2D>("nGoodTrksPerJetPt", "ntracks/JetPt_GoodTrk", 200, 0, 200, 30, 0, 30);
  nGoodPVTrksPerJetPt = fs->make<TH2D>("nGoodPVTrksPerJetPt", "ntracks/JetPt_GoodPVTrk", 200, 0, 200, 30, 0, 30);

  nJets               = fs->make<TH1D>("nJets", "nJets", 10,0,10);

  MASS               = fs->make<TH1D>("MASS","MASS", 2000, 0, 200);
  M_T                 = fs->make<TH1D>("M_T","M_T", 2000, 0, 200);
}


DemoTrackAnalyzer::~DemoTrackAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void DemoTrackAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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


  if( MuonHandle->size() == 0 ) return;
  if( JetHandle->size() == 0 ) return;


  vector<pat::Muon> muons = W_selection(MuonHandle, MetHandle, pvHandle);
  //vector<pat::Muon> muons = Z_selection(MuonHandle, pvHandle);
  if( muons.size() == 0 ) return;


  for( auto muon : muons ) { muonPt->Fill(muon.pt()); muonEta->Fill(muon.eta()); muonPhi->Fill(muon.phi()); }
  const pat::MET  met  = MetHandle->front();
  metUncorPt->Fill(met.pt()); metUncorPhi->Fill(met.phi());
  metPt->Fill(met.corPt(pat::MET::Type1XY)); metPhi->Fill(met.corPhi(pat::MET::Type1XY));


  vector<pat::Jet> jets = jet_selection(JetHandle, muons);
  
  for( auto jet : jets ) { jetPt->Fill(jet.pt()); jetEta->Fill(jet.eta()); jetPhi->Fill(jet.phi()); }
  nJets->Fill(jets.size());
  if( jets.size() == 0 ) return;
  
  /*
  cout<<endl;
  cout<<"============================================================="<<endl;
  cout<<"[Event id] run # = "<<iEvent.id().run()<<", evt # = "<<iEvent.id().event()<<endl;
  cout<<"[Jet Info] #(AllJets) = "<<JetHandle->size()<<", #(SelectedJets) = "<<jets.size()<<endl;
  */


  MaxJetPt = 500;
  for(int i=0; i<MaxJetPt; i++) { 
    vec_Ntrks.push_back({0,0}); vec_Ngoodtrks.push_back({0,0}); vec_NgoodPVtrks.push_back({0,0});
  }


  for(unsigned iJet=0; iJet < jets.size(); iJet++) {
    const pat::Jet jet = jets[iJet];

    vector<pat::PackedCandidate> jetContents;
    bool isMuon=false;
    for(unsigned iPC=0; iPC < pcHandle->size(); iPC++) {
      const pat::PackedCandidate pc = pcHandle->at(iPC);
      if( deltaR( jet , pc ) < 0.4 ) {
        //if( (abs(pc.pdgId())==11 || abs(pc.pdgId())==13 || abs(pc.pdgId())==15) && pc.pt() > 0.5*jet.pt() ) { isMuon=true; break; }
        jetContents.push_back(pc);
      }
    }

    if( jetContents.size() == 0 ) continue;
    if( isMuon ) continue;

    /*
    cout<<" - "<<endl;
    cout<<"[Jet Info] jet"<<iJet<<": ("<<jet.pt()<<", "<<jet.eta()<<", "<<jet.phi()<<")"<<endl;
    */

    double pt=0, trkpt=0, goodtrkpt=0, PVtrkpt=0;
    int nTrks_=0, nGoodTrks_=0, nGoodPVTrks_=0;
    for(unsigned iCont=0; iCont<jetContents.size(); iCont++) {
      pat::PackedCandidate jc = jetContents[iCont];

      pt += jc.pt();

      //bool fromPV = (jc.fromPV()>1 || fabs(jc.dz()) < 0.1);
      //cout<<"[Packed Candidate Info] jc"<<iCont<<": ("<<jc.pt()<<", "<<jc.eta()<<", "<<jc.phi()<<") id: "<<jc.pdgId()<<", track? "<<jc.hasTrackDetails()<<", PV"<<fromPV<<", purity"<<jc.trackHighPurity()<<", nPixHit "<<jc.numberOfPixelHits()<<", nHits "<<jc.numberOfHits()<<", jcLayer "<<jc.trackerLayersWithMeasurement()<<endl;
      //cout<<"[Packed Candidate Info] jc"<<iCont<<": ("<<jc.pt()<<", "<<jc.eta()<<", "<<jc.phi()<<") id: "<<jc.pdgId()<<", track? "<<jc.hasTrackDetails()<<", PV"<<jc.fromPV()<<", PV(dz)"<<fromPV<<", purity"<<jc.trackHighPurity()<<endl;

      if( jc.hasTrackDetails() ) {
        trkpt += jc.pt(); nTrks_++; trkPt->Fill(jc.pt());
        if( jc.trackHighPurity() ) {
          goodtrkpt += jc.pt(); nGoodTrks_++; goodTrkPt->Fill(jc.pt());
          if( jc.fromPV() > 0 || jc.fromPV() == -1 ) {
            PVtrkpt += jc.pt(); nGoodPVTrks_++; goodPVTrkPt->Fill(jc.pt());
            if( 30 < jet.pt() && jet.pt() < 40 ) goodPVTrkPt_pt30->Fill(jc.pt());
            else if( 100 < jet.pt() && jet.pt() < 120 ) goodPVTrkPt_pt100->Fill(jc.pt());
          }
        }
      }

    }

    nTrks->Fill(nTrks_);
    nGoodTrks->Fill(nGoodTrks_);
    nGoodPVTrks->Fill(nGoodPVTrks_);
    if( 30 < jet.pt() && jet.pt() < 40 ) nGoodPVTrks_pt30->Fill(nGoodPVTrks_);
    else if( 100 < jet.pt() && jet.pt() < 120 ) nGoodPVTrks_pt100->Fill(nGoodPVTrks_);

    nTrksPerJetPt->Fill(jet.pt(), nTrks_);
    nGoodTrksPerJetPt->Fill(jet.pt(), nGoodTrks_);
    nGoodPVTrksPerJetPt->Fill(jet.pt(), nGoodPVTrks_);

    trkPtSumRatio->Fill(trkpt/jet.pt());
    goodTrkPtSumRatio->Fill(goodtrkpt/jet.pt());
    goodPVTrkPtSumRatio->Fill(PVtrkpt/jet.pt());

    pcPtRatio->Fill(pt/jet.pt());

    if( jet.pt() < MaxJetPt ) {
      vec_Ntrks[(int)jet.pt()][0]++; vec_Ntrks[(int)jet.pt()][1]+=nTrks_;
      vec_Ngoodtrks[(int)jet.pt()][0]++; vec_Ngoodtrks[(int)jet.pt()][1]+=nGoodTrks_;
      vec_NgoodPVtrks[(int)jet.pt()][0]++; vec_NgoodPVtrks[(int)jet.pt()][1]+=nGoodPVTrks_;
    }

//    cout<<"[Packed Candidate Info] total pt       = "<<pt       <<endl;
//    cout<<"[Packed Candidate Info] track pt       = "<<trkpt    <<", ratio = "<<trkpt/pt<<endl;
//    cout<<"[Packed Candidate Info] good trk pt    = "<<goodtrkpt<<", ratio = "<<goodtrkpt/pt<<endl;
//    cout<<"[Packed Candidate Info] good PV trk pt = "<<PVtrkpt  <<", ratio = "<<PVtrkpt/pt<<endl;
  }


#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
   
}

std::vector<pat::Muon> DemoTrackAnalyzer::W_selection(edm::Handle<std::vector<pat::Muon>> MuonHandle, edm::Handle<std::vector<pat::MET>> MetHandle, edm::Handle<reco::VertexCollection> pvHandle)
{
  std::vector<pat::Muon> muons;

  const pat::MET  met  = MetHandle->front();
  if( !(met.pt() > 20) ) return muons;
  TLorentzVector met_4vec;
  //met_4vec.SetPtEtaPhiE(met.pt(), 0., met.phi(), met.pt());
  met_4vec.SetPtEtaPhiE(met.corPt(pat::MET::Type1XY), 0., met.corPhi(pat::MET::Type1XY), met.corPt(pat::MET::Type1XY));

  const reco::Vertex &vtx = pvHandle->front();
  for(unsigned iMuon = 0; iMuon < MuonHandle->size(); iMuon++) {
    const pat::Muon muon = MuonHandle->at(iMuon);

    if( !(muon.pt() > 10 && fabs(muon.eta()) < 2.4) ) continue;
    if( !(muon.isGlobalMuon()) ) continue;

    muons.push_back( muon );
  }

  if( muons.size() != 1 ) { muons.clear(); return muons; }

  pat::Muon muon = muons[0];
  if( muon.pt() > 26 && muon.isolationR03().sumPt/muon.pt() < 0.05 && muon.isTightMuon(vtx) ) {
    TLorentzVector muon_4vec;
    muon_4vec.SetPtEtaPhiM(muon.pt(), muon.eta(), muon.phi(), muon.mass());
    TLorentzVector W_4vec = muon_4vec + met_4vec;
    double mt = sqrt(2*muon.pt()*met.pt()*(1-cos(muon.phi()-met.phi()))); //W_4vec.Mt();
    MASS->Fill(W_4vec.M()); M_T->Fill(mt);
    if( 60 < mt && mt < 100 ) return muons;
  }

  muons.clear();
  return muons;
}


std::vector<pat::Muon> DemoTrackAnalyzer::Z_selection(edm::Handle<std::vector<pat::Muon>> MuonHandle, edm::Handle<reco::VertexCollection> pvHandle)
{
  std::vector<pat::Muon> muons;

  const reco::Vertex &vtx = pvHandle->front();
  for(unsigned iMuon = 0; iMuon < MuonHandle->size(); iMuon++) {
    const pat::Muon muon_1 = MuonHandle->at(iMuon);

    if( !(muon_1.pt() > 26) ) return muons;
    if( !(fabs(muon_1.eta()) < 2.4) ) continue;
    if( !(muon_1.isolationR03().sumPt/muon_1.pt() < 0.1) ) continue;
    if( !(muon_1.isTightMuon(vtx)) ) continue;

    for(unsigned jMuon = iMuon+1; jMuon < MuonHandle->size(); jMuon++) {
      const pat::Muon muon_2 = MuonHandle->at(jMuon);

      if( !(muon_2.pt() > 10) ) break;
      if( !(fabs(muon_2.eta()) < 2.4) ) continue;
      if( !(muon_2.isolationR03().sumPt/muon_2.pt() < 0.1) ) continue;
      if( !(muon_2.isTightMuon(vtx)) ) continue;

      TLorentzVector muon_1_4vec, muon_2_4vec;
      muon_1_4vec.SetPtEtaPhiM(muon_1.pt(), muon_1.eta(), muon_1.phi(), muon_1.mass());
      muon_2_4vec.SetPtEtaPhiM(muon_2.pt(), muon_2.eta(), muon_2.phi(), muon_2.mass());
      double mass = (muon_1_4vec + muon_2_4vec).M();
      MASS->Fill(mass);
      if( 70 < mass && mass < 110 ) { muons.push_back(muon_1); muons.push_back(muon_2); return muons; }
    }
  }

  return muons;
}


std::vector<pat::Jet> DemoTrackAnalyzer::jet_selection(edm::Handle<std::vector<pat::Jet>> JetHandle, std::vector<pat::Muon> muons)
{
  std::vector<pat::Jet> jets;
  for(unsigned iJet = 0; iJet < JetHandle->size(); iJet++) {
    const pat::Jet jet = JetHandle->at(iJet);
    bool muonVeto = false;
    for( auto muon : muons ) {
      if( deltaR( jet , muon ) < 0.5 ) { 
        muonVeto = true; break;
      }
    }
    if( muonVeto ) continue;
    if( !(30 < jet.pt()) ) continue;
    if( !(fabs(jet.eta()) < 2.4) ) continue;
    jets.push_back( jet );
  }
  return jets;
}


// ------------ method called once each job just before starting event loop  ------------
void DemoTrackAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void DemoTrackAnalyzer::endJob() 
{

}

// ------------ method called when starting to processes a run  ------------
/*
void 
DemoTrackAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
DemoTrackAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
DemoTrackAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
DemoTrackAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DemoTrackAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DemoTrackAnalyzer);
