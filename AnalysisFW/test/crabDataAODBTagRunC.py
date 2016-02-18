from WMCore.Configuration import Configuration
config = Configuration()

config.section_("User")
config.User.voGroup = 'dcms'

config.section_("General")
config.General.requestName = 'JetData-July2015-MagnField-Run2015C_25ns-BTag-76X-16Dec-v3'
config.General.workArea = 'JetDataWithJsonFile-MagnField-Run2015C_25ns-BTag-76X-16Dec-v3'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'ProcessedTreeProducer_dataPAT25ns_BTag_RunC_cfg.py'

config.section_("Data")
#config.Data.inputDataset = '/JetHT/Run2015D-PromptReco-v3/AOD'
#config.Data.inputDataset = '/JetHT/Run2015D-PromptReco-v4/AOD'
config.Data.inputDataset = '/JetHT/Run2015C_25ns-16Dec2015-v1/AOD'
config.Data.splitting = 'LumiBased' #LumiBased'
config.Data.unitsPerJob = 10
#config.Data.lumiMask = 'Cert_246908-251883_13TeV_PromptReco_Collisions15_JSON_v3.txt'
#config.Data.lumiMask = 'Cert_254833_13TeV_PromptReco_Collisions15_JSON.txt'
config.Data.lumiMask = 'Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_v2.txt'
config.Data.runRange = '246908-260627' # '193093-194075'
config.Data.outputDatasetTag = 'CRAB3_JetData-December2015-25ns-Run2015C-BTag-76X-16Dec-v3'

config.section_("Site")
config.Site.storageSite = "T2_DE_DESY"
