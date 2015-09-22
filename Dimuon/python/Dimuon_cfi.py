import FWCore.ParameterSet.Config as cms

dimuon = cms.EDAnalyzer('Dimuon',
                      genPartsTag=cms.InputTag("genParticles"),
                      decayParticlePID=cms.int32(13),
                      isCI=cms.bool(False),
                      debug=cms.bool(False),
                      status=cms.int32(23),
)
