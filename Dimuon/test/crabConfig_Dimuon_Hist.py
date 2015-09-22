from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'Dimuon_Histogram_M500_CI_L13000_Eta_13TeV_Sep7'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'Dimuon_cfg.py'
#config.JobType.pyCfgParams = ['pTMin=0','minMass=0','maxMass=100','ciGen=0','comEnergy=13']

config.Data.inputDataset = 	'/DYToMuMu/szaleski-M200_13TeV_pythia8_GEN-8adfb1edc07237acc3774d76f94543e7/USER'

#config.Data.primaryDataset = 'DYToMuMu'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 20
#NJOBS = 100
#config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.outLFNDirBase = '/store/user/szaleski/DimuonHistograms/CI/M_500/L13000/Sep7/'
#config.Data.publication = True
#config.Data.publishDataName = 'M200_13TeV_pythia8_GEN'

config.Site.storageSite = 'T3_US_FNALLPC' 
