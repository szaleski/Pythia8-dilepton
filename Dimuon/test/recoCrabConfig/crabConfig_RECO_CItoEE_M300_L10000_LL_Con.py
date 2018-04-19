from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'Dielectron_DIGIRECO_M300to800_CI_L10000_LL_Con_13TeV_Aug22'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'mc16RECO_cfg.py'
#config.JobType.numCores = 4

config.Data.inputDataset = '/CITo2E_GENSIM_Lam10/szaleski-EE_GENSIM_LLConM300-917d3a542e5db16211424ace778223b3/USER'
#config.Data.outputPrimaryDataset = 'CIToEEDigiRaw'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
#NJOBS = 100
#config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = True
config.Data.outputDatasetTag = 'EE_DIGIRAW_Lam10_LLConM300'

config.Site.storageSite = 'T3_US_FNALLPC' 
