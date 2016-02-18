from WMCore.Configuration import Configuration
config = Configuration()

config.section_("User")
config.User.voGroup = 'dcms'

config.section_("General")
config.General.requestName = 'JetData-July2015-MagnField-Run2015D_25ns-BTag-76X-ReReco-v3'
config.General.workArea = 'JetDataWithJsonFile-MagnField-Run2015D_25ns-BTag-76X-ReReco-v3'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'ProcessedTreeProducer_dataPAT25ns_BTag_cfg.py'

config.section_("Data")
config.Data.inputDataset = '/JetHT/Run2015D-16Dec2015-v1/AOD'
#config.Data.inputDataset = '/JetHT/Run2015D-16Dec2015-v1/AOD'
#config.Data.inputDataset = '/JetHT/Run2015C_25ns-05Oct2015-v1/AOD '
config.Data.splitting = 'LumiBased' #LumiBased'
config.Data.unitsPerJob = 10
config.Data.lumiMask = 'Cert_13TeV_16Dec2015ReReco_Collisions15_25ns_JSON.txt'
config.Data.runRange = '246908-260627' # '193093-194075'
config.Data.outputDatasetTag = 'CRAB3_JetData-25ns-Run2015D_v3-BTag-76X-ReReco-v3'

config.section_("Site")
config.Site.storageSite = "T2_DE_DESY"
#config.Site.whitelist = ['T2_DE_DESY']
