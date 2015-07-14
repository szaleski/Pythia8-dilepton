import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:$CMSSW_BASE/src/GenStudy/Dimuon/test/DYToMuMu_M-1_13TeV_pythia8_GEN.root'
    )
)
from GenStudy.Dimuon.Dimuon_cfi import *
process.Dimuon=dimuon.clone()

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("DYoutput_cut.root")
)



process.p = cms.Path(process.Dimuon)
