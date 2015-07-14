import FWCore.ParameterSet.Config as cms

dimuon = cms.EDAnalyzer('Dimuon',
                      genPartsTag=cms.InputTag("genParticles"),
                      decayParticlePID=cms.int32(13)                        
)
