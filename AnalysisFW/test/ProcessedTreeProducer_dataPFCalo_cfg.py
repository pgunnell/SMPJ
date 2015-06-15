# -*- coding: utf-8 -*-
import FWCore.ParameterSet.Config as cms

process = cms.Process("Ntuplizer")
process.load('FWCore.MessageService.MessageLogger_cfi')
##-------------------- Communicate with the DB -----------------------
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'MCRUN2_74_V9A::All'
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
#rocess.load('EventFilter.HcalRawToDigi.hcallaserhbhehffilter2012_cfi')
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
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
#############   Define the source file ###############
process.source = cms.Source("PoolSource",
#    fileNames = cms.untracked.vstring('/store/data/Run2012C/JetHT/AOD/PromptReco-v2/000/198/941/221BEAC8-7ACF-E111-8A3D-E0CB4E4408E7.root')
   #fileNames = cms.untracked.vstring('root://eoscms//eos/cms//store/data/Run2012D/SingleMu/AOD/22Jan2013-v1/20000/F893D87B-1585-E211-B143-001EC9D8D973.root')
    fileNames = cms.untracked.vstring('file://./CE5DBB3C-170B-E511-A908-02163E01473C.root')


)
############# processed tree producer ##################
process.TFileService = cms.Service("TFileService",fileName = cms.string('DATA_ProcessedTreePFCalo.root'))

process.ak4 = cms.EDAnalyzer('ProcessedTreeProducerPFCalo',
    ## jet collections ###########################
    pfjets          = cms.InputTag('ak4PFJets'),
    calojets        = cms.InputTag('ak4CaloJets'),
    ## database entry for the uncertainties ######
    #PFPayloadName   = cms.string('AK4PF'),
    #CaloPayloadName = cms.string('AK4Calo'),
    PFPayloadName   = cms.string(''),
    CaloPayloadName = cms.string(''),
    jecUncSrc       = cms.string(''),
    jecUncSrcNames  = cms.vstring(''),
    
    ## calojet ID and extender for the JTA #######
    calojetID       = cms.InputTag('ak4JetID'),
    calojetExtender = cms.InputTag('ak4JetExtender'),
    ## set the conditions for good Vtx counting ##
    offlineVertices = cms.InputTag('offlinePrimaryVertices'),
    goodVtxNdof     = cms.double(4), 
    goodVtxZ        = cms.double(24),
    ## rho #######################################
   srcCaloRho      = cms.InputTag('fixedGridRhoFastjetAllCalo'),
   srcPFRho        = cms.InputTag('fixedGridRhoFastjetAll'),
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
    triggerName     = cms.vstring('HLT_L1SingleJet36_v1','HLT_L1SingleJet68_v1','HLT_ZeroBias_part0_v1','HLT_ZeroBias_part1_v1','HLT_ZeroBias_part3_v1','HLT_ZeroBias_part4_v1','HLT_ZeroBias_part5_v1'),
    triggerResults  = cms.InputTag("TriggerResults","","HLT"),
    triggerEvent    = cms.InputTag("hltTriggerSummaryAOD","","HLT"),
    ## jec services ##############################
    pfjecService    = cms.string('ak4PFL1FastL2L3Residual'),
    calojecService  = cms.string('ak4CaloL1FastL2L3Residual')
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



#process.path = cms.Path(process.hltFilter * process.HBHENoiseFilterResultProducer * process.hcalLaserEventFilter *
#process.hcallasereventfilter2012 * #process.hcallaserhbhehffilter2012 *
#process.EcalDeadCellTriggerPrimitiveFilter * process.eeBadScFilter  process.ak4)
process.path = cms.Path(process.HBHENoiseFilterResultProducer*process.EcalDeadCellTriggerPrimitiveFilter * process.eeBadScFilter *process.ak4)

