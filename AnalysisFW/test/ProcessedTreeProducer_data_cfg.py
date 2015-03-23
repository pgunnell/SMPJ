import FWCore.ParameterSet.Config as cms

process = cms.Process("Ana")
process.load('FWCore.MessageService.MessageLogger_cfi')
##-------------------- Communicate with the DB -----------------------
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'FT_53_V21_AN6::All'
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.Geometry.GeometryIdeal_cff')
process.load('RecoJets.Configuration.RecoPFJets_cff')
process.load('RecoJets.Configuration.RecoJets_cff')
process.load('CommonTools/RecoAlgos/HBHENoiseFilterResultProducer_cfi')
##-------------------- Import the JEC services -----------------------
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
#-------- HCAL Laser Filter ------#
process.load("RecoMET.METFilters.hcalLaserEventFilter_cfi")
process.load('EventFilter.HcalRawToDigi.hcallasereventfilter2012_cfi')
#process.load('EventFilter.HcalRawToDigi.hcallaserhbhehffilter2012_cfi')
## The ECAL dead cell trigger primitive filter _______________________________||
process.load('RecoMET.METFilters.EcalDeadCellTriggerPrimitiveFilter_cfi')
## The EE bad SuperCrystal filter ____________________________________________||
process.load('RecoMET.METFilters.eeBadScFilter_cfi')
## TRACKING Filter 
process.load('RecoMET.METFilters.trackingFailureFilter_cfi')
## The ECAL laser correction filter
process.load('RecoMET.METFilters.ecalLaserCorrFilter_cfi')

#############   Set the number of events #############
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
#############   Format MessageLogger #################
process.MessageLogger.cerr.FwkReport.reportEvery = 10
#############   Define the source file ###############
process.source = cms.Source("PoolSource",
#    fileNames = cms.untracked.vstring('/store/data/Run2012C/JetHT/AOD/PromptReco-v2/000/198/941/221BEAC8-7ACF-E111-8A3D-E0CB4E4408E7.root')
   fileNames = cms.untracked.vstring('root://eoscms//eos/cms//store/data/Run2012D/SingleMu/AOD/22Jan2013-v1/20000/F893D87B-1585-E211-B143-001EC9D8D973.root')
)
############# processed tree producer ##################
process.TFileService = cms.Service("TFileService",fileName = cms.string('ProcessedTree_data.root'))

process.ak7 = cms.EDAnalyzer('ProcessedTreeProducer',
    ## jet collections ###########################
    pfjets          = cms.InputTag('ak7PFJets'),
    calojets        = cms.InputTag('ak7CaloJets'),
    ## database entry for the uncertainties ######
    PFPayloadName   = cms.string('AK7PF'),
    CaloPayloadName = cms.string('AK7Calo'),
    jecUncSrc       = cms.string('Summer13_V4_DATA_UncertaintySources_AK7PF.txt'),
    jecUncSrcNames  = cms.vstring('Absolute','HighPtExtra','SinglePionECAL','SinglePionHCAL','FlavorQCD',
                                  'Time','RelativeJEREC1','RelativeJEREC2','RelativeJERHF','RelativePtBB',
                                  'RelativePtEC1','RelativePtEC2','RelativePtHF','RelativeFSR','RelativeStatEC2',
                                  'RelativeStatHF','PileUpDataMC','PileUpPtBB','PileUpPtEC','PileUpPtHF','PileUpBias',
                                  'SubTotalPileUp','SubTotalRelative','SubTotalPt','SubTotalMC','Total','TotalNoFlavor',
                                  'FlavorZJet','FlavorPhotonJet','FlavorPureGluon','FlavorPureQuark','FlavorPureCharm',
                                  'FlavorPureBottom'),
    
    ## calojet ID and extender for the JTA #######
    calojetID       = cms.InputTag('ak7JetID'),
    calojetExtender = cms.InputTag('ak7JetExtender'),
    ## set the conditions for good Vtx counting ##
    offlineVertices = cms.InputTag('offlinePrimaryVertices'),
    goodVtxNdof     = cms.double(4), 
    goodVtxZ        = cms.double(24),
    ## rho #######################################
    srcCaloRho      = cms.InputTag('kt6CaloJets','rho'),
    srcPFRho        = cms.InputTag('kt6PFJets','rho'),
    ## preselection cuts #########################
    maxY            = cms.double(5.0), 
    minPFPt         = cms.double(20),
    minPFFatPt      = cms.double(10),
    maxPFFatEta     = cms.double(2.5),
    minCaloPt       = cms.double(20),
    minNPFJets      = cms.int32(1),
    minNCaloJets    = cms.int32(1), 
    minJJMass       = cms.double(-1),
    ## trigger ###################################
    printTriggerMenu = cms.untracked.bool(True),
    processName     = cms.string('HLT'),
    triggerName     = cms.vstring('HLT_PFJet40_v3','HLT_PFJet40_v4','HLT_PFJet40_v5','HLT_PFJet40_v6','HLT_PFJet40_v7','HLT_PFJet40_v8','HLT_PFJet40_v9',
                                  'HLT_PFJet80_v3','HLT_PFJet80_v4','HLT_PFJet80_v5','HLT_PFJet80_v6','HLT_PFJet80_v7','HLT_PFJet80_v8','HLT_PFJet80_v9',
                                  'HLT_PFJet140_v3','HLT_PFJet140_v4','HLT_PFJet140_v5','HLT_PFJet140_v6','HLT_PFJet140_v7','HLT_PFJet140_v8','HLT_PFJet140_v9',
                                  'HLT_PFJet200_v3','HLT_PFJet200_v4','HLT_PFJet200_v5','HLT_PFJet200_v6','HLT_PFJet200_v7','HLT_PFJet200_v8','HLT_PFJet200_v9',
                                  'HLT_PFJet260_v3','HLT_PFJet260_v4','HLT_PFJet260_v5','HLT_PFJet260_v6','HLT_PFJet260_v7','HLT_PFJet260_v8','HLT_PFJet260_v9',
                                  'HLT_PFJet320_v3','HLT_PFJet320_v4','HLT_PFJet320_v5','HLT_PFJet320_v6','HLT_PFJet320_v7','HLT_PFJet320_v8','HLT_PFJet320_v9',
                                  'HLT_PFJet400_v3','HLT_PFJet400_v4','HLT_PFJet400_v5','HLT_PFJet400_v6','HLT_PFJet400_v7','HLT_PFJet400_v8','HLT_PFJet400_v9',
                                  'HLT_IsoMu24_eta2p1_v11', 'HLT_IsoMu24_eta2p1_v12', 'HLT_IsoMu24_eta2p1_v13', 'HLT_IsoMu24_eta2p1_v14', 'HLT_IsoMu24_eta2p1_v15'
    ),
    triggerResults  = cms.InputTag("TriggerResults","","HLT"),
    triggerEvent    = cms.InputTag("hltTriggerSummaryAOD","","HLT"),
    ## jec services ##############################
    pfjecService    = cms.string('ak7PFL1FastL2L3Residual'),
    calojecService  = cms.string('ak7CaloL1FastL2L3Residual')
)

process.ak5 = process.ak7.clone(
    pfjets           = 'ak5PFJets',
    calojets         = 'ak5CaloJets',
    PFPayloadName    = 'AK5PF',
    CaloPayloadName  = 'AK5Calo',
    jecUncSrc        = 'Summer13_V4_DATA_UncertaintySources_AK5PF.txt',
    calojetID        = 'ak5JetID',
    calojetExtender  = 'ak5JetExtender',
    pfjecService     = 'ak5PFL1FastL2L3Residual',
    calojecService   = 'ak5CaloL1FastL2L3Residual',
    printTriggerMenu = False 
)

############# hlt filter #########################
process.hltFilter = cms.EDFilter('HLTHighLevel',
    TriggerResultsTag  = cms.InputTag('TriggerResults','','HLT'),
    HLTPaths           = cms.vstring(
    'HLT_PFJet40_v*','HLT_PFJet80_v*','HLT_PFJet140_v*','HLT_PFJet200_v*','HLT_PFJet260_v*','HLT_PFJet320_v*','HLT_PFJet400_v*', 'HLT_IsoMu24_eta2p1_v*'
    ),
    eventSetupPathsKey = cms.string(''),
    andOr              = cms.bool(True), #----- True = OR, False = AND between the HLTPaths
    throw              = cms.bool(False)
)

process.hltMuFilter = cms.EDFilter('HLTHighLevel',
    TriggerResultsTag  = cms.InputTag('TriggerResults','','HLT'),
    HLTPaths           = cms.vstring(
    'HLT_IsoMu24_eta2p1_v*'
    ),
    eventSetupPathsKey = cms.string(''),
    andOr              = cms.bool(True), #----- True = OR, False = AND between the HLTPaths
    throw              = cms.bool(False)
)



process.path = cms.Path(process.hltFilter * process.HBHENoiseFilterResultProducer * process.hcalLaserEventFilter *
process.hcallasereventfilter2012 * #process.hcallaserhbhehffilter2012 *
process.EcalDeadCellTriggerPrimitiveFilter * process.eeBadScFilter * process.ak5 * process.ak7)


