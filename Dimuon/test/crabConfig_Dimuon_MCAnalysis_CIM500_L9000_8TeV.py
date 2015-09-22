from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'Dimuon_MCAnalysis_M500_CI_L9000_Con_8TeV_Sep1'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'PrivateMC'
config.JobType.psetName = 'ciCrossSecCalc_edm_cfg.py'
config.JobType.pyCfgParams = ['pTMin=0','minMass=500','Lambda=9000','helicityLL=-1','ciGen=1','comEnergy=8']

#config.Data.inputDataset = '/GenericTTbar/HC-CMSSW_5_3_1_START53_V5-v1/GEN-SIM-RECO'
config.Data.primaryDataset = 'DYToMuMu'
#config.Data.inputDBS = 'global'
config.Data.splitting = 'EventBased'
config.Data.unitsPerJob = 50000
NJOBS = 100
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
#config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = True
config.Data.publishDataName = 'CI_M500_L9000_8TeV_pythia8_GEN'

config.Site.storageSite = 'T3_US_FNALLPC' 
