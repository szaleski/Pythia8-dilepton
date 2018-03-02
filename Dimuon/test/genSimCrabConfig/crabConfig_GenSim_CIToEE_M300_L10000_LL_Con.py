from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'Dielectron_GENSIM_M300to800_CI_L10000_LL_Con_13TeV_Aug14'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'PrivateMC'
config.JobType.psetName = 'MC16GenSim_cfg.py'
config.JobType.pyCfgParams = ['minMass=300','maxMass=800','ciGen=1','helicityLL=-1','Lambda=10000','pdgId=11']

#config.Data.inputDataset = '/GenericTTbar/HC-CMSSW_5_3_1_START53_V5-v1/GEN-SIM-RECO'
config.Data.outputPrimaryDataset = 'CITo2E_GENSIM_Lam10'
#config.Data.inputDBS = 'global'
config.Data.splitting = 'EventBased'
config.Data.unitsPerJob = 500
NJOBS = 100
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.outLFNDirBase = '/store/user/szaleski/' 
config.Data.publication = True
config.Data.outputDatasetTag = 'EE_GENSIM_LLConM300'

config.Site.storageSite = 'T3_US_FNALLPC' 
