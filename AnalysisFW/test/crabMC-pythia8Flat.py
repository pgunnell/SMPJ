from WMCore.Configuration import Configuration
config = Configuration()

config.section_("User")
config.User.voGroup = 'dcms'

config.section_("General")
config.General.requestName = 'Herwig-25ns-Ntuples-July2016-Flat-v2'
config.General.workArea = 'Herwig-25ns-Ntuples-July2016-TuneCUETM1-Flat-v2'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'ProcessedTreeProducer_MCPAT25ns_BTag_cfg.py'

config.section_("Data")
#config.Data.inputDataset = '/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/RunIISpring16DR80-PUSpring16_magnetOn_80X_mcRun2_asymptotic_2016_v3-v1/AODSIM'
config.Data.inputDataset = '/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/AODSIM'
config.Data.splitting = 'EventAwareLumiBased' #LumiBased'
config.Data.unitsPerJob = 12000
config.Data.outputDatasetTag = 'CRAB3_Herwig-25ns-NtuplesTuneCUETS1-July2016-Flat-v2'

config.section_("Site")
config.Site.storageSite = "T2_DE_DESY"

#/QCD_Pt_15to30_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2/AODSIM
#/QCD_Pt_30to50_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2/AODSIM
#/QCD_Pt_50to80_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2/AODSIM

#/QCD_Pt_80to120_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM
#/QCD_Pt_120to170_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM
#/QCD_Pt_170to300_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2/AODSIM
#/QCD_Pt_300to470_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM
#/QCD_Pt_470to600_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2/AODSIM
#/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v3/AODSIM
#/QCD_Pt_800to1000_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2/AODSIM
#/QCD_Pt_1000to1400_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM
#/QCD_Pt_1400to1800_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM
#/QCD_Pt_1800to2400_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM
#/QCD_Pt_2400to3200_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM
#/QCD_Pt_3200toInf_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM
