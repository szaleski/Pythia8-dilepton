from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'Dimuon_GENSIM17_M300to800_CI_L16000_LL_Con_13TeV_Feb5'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'PrivateMC'
config.JobType.psetName = 'mc17genSimcfg.py'
config.JobType.numCores = 8
config.JobType.pyCfgParams = ['minMass=300','maxMass=800','Lambda=16000','helicityLL=-1','ciGen=1','pdgId=13']

#config.Data.inputDataset = '/GenericTTbar/HC-CMSSW_5_3_1_START53_V5-v1/GEN-SIM-RECO'
config.Data.outputPrimaryDataset = 'CITo2Mu_L16TeV_GENSIM17_Test'
#config.Data.inputDBS = 'global'
config.Data.splitting = 'EventBased'
config.Data.unitsPerJob = 500
NJOBS = 100
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.outLFNDirBase = '/store/user/szaleski/'
config.Data.publication = True
config.Data.outputDatasetTag = 'MuMu_16TeV_GENSIM17_LLConM300'

config.Site.storageSite = 'T3_US_FNALLPC' 
