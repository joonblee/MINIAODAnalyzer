import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                            # replace 'myfile.root' with the source file you want to use
                            # fileNames = cms.untracked.vstring('file:/afs/cern.ch/cms/Tutorials/TWIKI_DATA/TTJets_8TeV_53X.root')
                            fileNames = cms.untracked.vstring('file:/data6/Users/joonblee/MINIAODAnalyzer/CMSSW_10_4_0_patch1/src/Demo/DemoAnalyzer/Inputs/___INPUTFILES___')
                            )

#process.demo = cms.EDAnalyzer('DemoAnalyzer'
#                              )

process.TFileService = cms.Service("TFileService", 
                                   fileName = cms.string("/data6/Users/joonblee/MINIAODAnalyzer/CMSSW_10_4_0_patch1/src/Demo/DemoAnalyzer/___OUTPUTFILES___"),
                                   closeFileFast = cms.untracked.bool(True)
)

process.demo = cms.EDAnalyzer('DemoTrackAnalyzer',
                              pfCands = cms.untracked.InputTag("packedPFCandidates"),
                              Jet = cms.untracked.InputTag("slimmedJets"),
                              Muon = cms.untracked.InputTag("slimmedMuons"),
                              Electron = cms.untracked.InputTag("slimmedElectrons"),
                              Photon = cms.untracked.InputTag("slimmedPhotons"),
                              MET = cms.untracked.InputTag("slimmedMETs"),
                              PrimaryVertex = cms.untracked.InputTag("offlineSlimmedPrimaryVertices")
                              )

process.p = cms.Path(process.demo)
