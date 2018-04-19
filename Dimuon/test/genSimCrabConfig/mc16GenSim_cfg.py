# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: Configuration/GenProduction/python/ThirteenTeV/MC16_CITo2E_M300_CUETP8M1_Lam10TeVConLL_13TeV_Pythia8.py --fileout file:EXO-RunIISummer15GS-09252.root --mc --eventcontent RAWSIM --customise SLHCUpgradeSimulations/Configuration/postLS1Customs.customisePostLS1,Configuration/DataProcessing/Utils.addMonitoring --datatier GEN-SIM --conditions MCRUN2_71_V1::All --beamspot Realistic50ns13TeVCollision --step GEN,SIM --magField 38T_PostLS1 --python_filename /afs/cern.ch/cms/PPD/PdmV/work/McM/submit/EXO-RunIISummer15GS-09252/EXO-RunIISummer15GS-09252_1_cfg.py --no_exec -n 329
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing ('analysis')
from GenStudy.Dimuon.genSimCmdLineOptions_cfi import registerDefaultMCOptions
registerDefaultMCOptions(options)
options.parseArguments()

process = cms.Process('SIM')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.Geometry.GeometrySimDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic50ns13TeVCollision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

# Input source
process.source = cms.Source("EmptySource")


process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.19 $'),
    annotation = cms.untracked.string('Configuration/GenProduction/python/ThirteenTeV/MC16_CITo2E_M300_CUETP8M1_Lam10TeVConLL_13TeV_Pythia8.py nevts:329'),
    name = cms.untracked.string('Applications')
)

# Output definition

process.RAWSIMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    fileName = cms.untracked.string('file:MC16output.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    )
)

# Additional output definition

# Other statements
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_mcRun2_asymptotic_2016_TrancheIV_v6', '')

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
    pythiaPylistVerbosity = cms.untracked.int32(1),
    filterEfficiency = cms.untracked.double(1.0),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    comEnergy = cms.double(13000.0),
    maxEventsToPrint = cms.untracked.int32(5),
    PythiaParameters = cms.PSet(
        pythia8CommonSettings = cms.vstring('Tune:preferLHAPDF = 2', 
            'Main:timesAllowErrors = 10000', 
            'Check:epTolErr = 0.01', 
            'Beams:setProductionScalesFromLHEF = off', 
            'SLHA:keepSM = on', 
            'SLHA:minMassSM = 1000.', 
            'ParticleDecays:limitTau0 = on', 
            'ParticleDecays:tau0Max = 10', 
            'ParticleDecays:allowPhotonRadiation = on'),
        pythia8CUEP8M1Settings = cms.vstring('Tune:pp 14', 
            'Tune:ee 7', 
            'MultipartonInteractions:pT0Ref=2.4024', 
            'MultipartonInteractions:ecmPow=0.25208', 
            'MultipartonInteractions:expPow=1.6'),
        processParameters = cms.vstring(
#            'ParticleDecays:limitTau0 = on',
#            'ParticleDecays:tauMax = 10',
#            'Tune:pp 5',
#            'Tune:ee 3',
#            'WeakSingleBoson:ffbar2ffbar(s:gmZ)= on',
#            '23:onMode = off',
#            '23:onIfAny = 13',
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
            'ContactInteractions:etaLL = '+str(options.helicityLL),
            'ContactInteractions:etaLR = '+str(options.helicityLR),
            'ContactInteractions:etaRR = '+str(options.helicityRR),
            'ContactInteractions:Lambda = '+str(options.Lambda),

#            'ContactInteractions:QCffbar2mumubar = on',
#            'PhaseSpace:mHatMin = 300',

#            'ContactInteractions:Lambda = 10000',
#            'ContactInteractions:etaLL = -1',
#            'ContactInteractions:etaRR = 0',
#            'ContactInteractions:etaLR = 0',


       ),
        parameterSets = cms.vstring('pythia8CommonSettings', 
            'pythia8CUEP8M1Settings', 
            'processParameters')
    )
)


# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RAWSIMoutput_step = cms.EndPath(process.RAWSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.endjob_step,process.RAWSIMoutput_step)
# filter all path with the production filter sequence
for path in process.paths:
    getattr(process,path)._seq = process.generator * getattr(process,path)._seq 

# customisation of the process.

# Automatic addition of the customisation function from Configuration.DataProcessing.Utils
from Configuration.DataProcessing.Utils import addMonitoring 

#call to customisation function addMonitoring imported from Configuration.DataProcessing.Utils
process = addMonitoring(process)

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.postLS1Customs
from SLHCUpgradeSimulations.Configuration.postLS1Customs import customisePostLS1 

#call to customisation function customisePostLS1 imported from SLHCUpgradeSimulations.Configuration.postLS1Customs
process = customisePostLS1(process)

# End of customisation functions
