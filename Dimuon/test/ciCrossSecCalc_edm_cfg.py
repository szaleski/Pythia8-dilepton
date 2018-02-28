
import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing ('analysis')
from GenStudy.Dimuon.mcCmdLineOptions_cfi import registerDefaultMCOptions
registerDefaultMCOptions(options)
options.register ('zPrimeModel',
                  "zPrimeSSM",
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.string,
                  "which Z' model to use")
options.register ('interferenceMode',
                  3,
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.int,          
                  "Z/gamma/Z' interference setting")
options.parseArguments()


import FWCore.ParameterSet.Config as cms


process = cms.Process('GEN')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedNominalCollision2015_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.MessageLogger.cerr.FwkReport = cms.untracked.PSet(
    reportEvery = cms.untracked.int32(500),
    limit = cms.untracked.int32(10000000)
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

# Input source
process.source = cms.Source("EmptySource")

process.options = cms.untracked.PSet(

)
process.source = cms.Source("EmptySource",
			    firstLuminosityBlock = cms.untracked.uint32(options.seed),
			    numberEventsInLuminosityBlock = cms.untracked.uint32(100)
)

process.RandomNumberGeneratorService.generator.initialSeed = options.seed*10000

# Other statements
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '93X_upgrade2023_realistic_v3')
isDY = "on"
isCImumu = "off"
isCIee = "off"

if options.ciGen:
	isDY = "off"
	if options.pdgId == 11:
		isCIee = "on"
		isCImumu = "off"
	else:
		isCIee = "off"
		isCImumu = "on"


process.generator = cms.EDFilter("Pythia8GeneratorFilter",
        comEnergy =  cms.double(options.comEnergy*1000),
        crossSection = cms.untracked.double(1e10),
        filterEfficiency = cms.untracked.double(1),
        maxEventsToPrint = cms.untracked.int32(0),
        pythiaHepMCVerbosity = cms.untracked.bool(False),
        pythiaPylistVerbosity = cms.untracked.int32(1),
        PythiaParameters = cms.PSet(
                processParameters = cms.vstring(
                        'Main:timesAllowErrors    = 10000',
                        'ParticleDecays:limitTau0 = on',
                        'ParticleDecays:tauMax = 10',
                        'Tune:pp 5',
                        'Tune:ee 3',
 			'WeakSingleBoson:ffbar2ffbar(s:gmZ)= %s'%(isDY),
                        'ContactInteractions:QCffbar2eebar = %s'%(isCIee),
                        'ContactInteractions:QCffbar2mumubar = %s'%(isCImumu),
                        '23:onMode = off',
                        '23:onIfAny = '+str(options.pdgId),
			'PartonLevel:MPI = '+str(options.ULE),
			'PartonLevel:ISR = '+str(options.ISR),
			'PartonLevel:FSR = '+str(options.FSR),
			'PhaseSpace:pTHatMin = '+str(options.pTMin),
                        'PhaseSpace:mHatMin = '+str(options.minMass),
                        'PhaseSpace:mHatMax = '+str(options.maxMass),
			'PhaseSpace:pTHatMax = '+str(options.pTMax),
			'ContactInteractions:Lambda = '+str(options.Lambda),
			'ContactInteractions:etaLL = '+str(options.helicityLL),
			'ContactInteractions:etaLR = '+str(options.helicityLR),
			'ContactInteractions:etaRR = '+str(options.helicityRR),

                ),
                parameterSets = cms.vstring('processParameters')
        )
)


process.ProductionFilterSequence = cms.Sequence(process.generator)


process.AODSIMoutput = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    ),
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(4),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM'),
        filterName = cms.untracked.string('')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(15728640),
    fileName =
cms.untracked.string('file:CITo_PID%d_M%d_D%d_L%d_LL%d_LR_%d_RR_%d_13TeV_pythia8_GEN_Sep9.root'%(options.pdgId,options.minMass,options.pTMin,options.Lambda,options.helicityLL,options.helicityLR,options.helicityRR)),
    outputCommands = process.AODSIMEventContent.outputCommands
)

process.AODSIMoutput_step = cms.EndPath(process.AODSIMoutput)










# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)#*process.crossSecTreeMaker*process.pdfTreeMaker)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.endjob_step,process.AODSIMoutput_step)

# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.ProductionFilterSequence * getattr(process,path)._seq


###-- Dump config -----------------------------------------------
file = open('GenTest.txt','w')
file.write(str(process.dumpPython()))
file.close()

