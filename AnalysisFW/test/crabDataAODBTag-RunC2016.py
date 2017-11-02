from WMCore.Configuration import Configuration
config = Configuration()

config.section_("User")
config.User.voGroup = 'dcms'

config.section_("General")
config.General.requestName = 'JetData-June2016-MagnField-Run2016C_25ns-BTag-80X-PromptReco-v2'
config.General.workArea = 'JetDataWithJsonFile-MagnField-Run2016C_25ns-BTag-80X-PromptReco-v2'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'ProcessedTreeProducer_dataPAT25ns_BTag_cfg.py'
config.JobType.allowUndistributedCMSSW = True


config.section_("Data")
config.Data.inputDataset = '/JetHT/Run2016C-PromptReco-v2/AOD'
#config.Data.inputDataset = '/JetHT/Run2015D-16Dec2015-v2/AOD'
#config.Data.inputDataset = '/JetHT/Run2015C_25ns-05Oct2015-v2/AOD '
config.Data.splitting = 'LumiBased' #LumiBased'
config.Data.unitsPerJob = 10
config.Data.lumiMask = 'Cert_271036-276097_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt'
config.Data.runRange = '271036-276097' # '193093-194075'
config.Data.outputDatasetTag = 'CRAB3_JetData-25ns-Run2016C_v2-BTag-80X-PromptReco-v2'

config.section_("Site")
config.Site.storageSite = "T2_DE_DESY"
#config.Site.whitelist = ['T2_DE_DESY']
