import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing ('analysis')
from GenStudy.Dimuon.dileptonMcCmndLineOptions import registerDefaultMCOptions
registerDefaultMCOptions(options)
options.parseArguments()

import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.MessageLogger.cerr.FwkReport = cms.untracked.PSet(
    reportEvery = cms.untracked.int32(1),
    limit = cms.untracked.int32(10000000)
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(200000) )
#filedir = '/DYToMuMu/szaleski-M200_13TeV_pythia8_GEN-c932f08e13e94abe45c501c7b6a032e1/USER'
process.source = cms.Source("PoolSource",
            # replace 'myfile.root' with the source file you want to use
                        fileNames = cms.untracked.vstring(
                                  'file:$CMSSW_BASE/src/GenStudy/Dimuon/test/%s'%(options.inFile)
                                )
                            )

from GenStudy.Dimuon.Dimuon_cfi import *
process.Dimuon=dimuon.clone()
process.Dimuon.debug=3

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("file:%s.root"%(options.filename)),  
#                                   fileName = cms.string("Pythia8_Jul22_CIEE_13TeV_DYM1300.root")  
#                                   fileName = cms.string("Pythia8_Jul6_CIEE_13TeV_CIM1500LLCon.root")
#                                   fileName = cms.string("Pythia8_Jul6_CIEE_13TeV_CIM1800LLConTest.root")
)



process.p = cms.Path(process.Dimuon)
