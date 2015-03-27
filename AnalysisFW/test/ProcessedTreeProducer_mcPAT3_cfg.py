# -*- coding: utf-8 -*-
import FWCore.ParameterSet.Config as cms

process = cms.Process("Ntuplizer")
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#! Conditions
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.load('Configuration.StandardSequences.Geometry_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('RecoJets.Configuration.GenJetParticles_cff')
process.load('RecoJets.Configuration.RecoGenJets_cff')

process.GlobalTag.globaltag = "PHYS14_25_V2::All"


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#! Input
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
                                    
inFiles = cms.untracked.vstring(
'file:///mnt/storage/gflouris/08C07BB6-376F-E411-BE9F-C4346BC7EE18.root' 

   )
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(20))
process.source = cms.Source("PoolSource", fileNames = inFiles )

from jetToolbox_cff import *
jetToolbox( process, 'ak4', 'ak4JetSubs', 'out', JETCorrLevels = ['None']) 
jetToolbox( process, 'ak5', 'ak5JetSubs', 'out', JETCorrLevels = ['None']) 
jetToolbox( process, 'ak7', 'ak7JetSubs', 'out', JETCorrLevels = ['None']) 


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#! Services
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.load('CommonTools.UtilAlgos.TFileService_cfi')
process.TFileService.fileName=cms.string('MC_ProcessedTreeProducer.root')


#validator_sequence = cms.Sequence()
#setattr(process,"validator_sequence",validator_sequence)

process.ak4 =  cms.EDAnalyzer('ProcessedTreeProducer',
	## jet collections ###########################
	pfjets          = cms.InputTag('selectedPatJetsAK4PFCHS'),
	pfjetschs       = cms.InputTag('selectedPatJetsAK4PFCHS'),
	## MET collection ####
	pfmet           = cms.InputTag('pfMet'),
	genjets         = cms.untracked.InputTag('ak4GenJets'),
	## database entry for the uncertainties ######
	PFPayloadName   = cms.string(''),
	PFPayloadNameCHS= cms.string(''),
	CaloPayloadName = cms.string(''),
	jecUncSrc       = cms.string(''),
	jecUncSrcCHS    = cms.string(''), 
	jecUncSrcNames  = cms.vstring(''),
	## set the conditions for good Vtx counting ##
	offlineVertices = cms.InputTag('goodOfflinePrimaryVertices'),
	goodVtxNdof     = cms.double(4), 
	goodVtxZ        = cms.double(24),
	## rho #######################################
	srcCaloRho      = cms.InputTag('ak4CaloJets','rho'),
	srcPFRho        = cms.InputTag('ak4CaloJets','rho'),
	srcPU           = cms.untracked.InputTag('addPileupInfo'),
	## preselection cuts #########################
	maxY            = cms.double(5.0), 
	minPFPt         = cms.double(20),
	minPFFatPt      = cms.double(10),
	maxPFFatEta     = cms.double(2.5),
	minNPFJets      = cms.int32(1),
	minGenPt        = cms.untracked.double(20),
	minJJMass       = cms.double(-1),
	isMCarlo        = cms.untracked.bool(True),
	useGenInfo      = cms.untracked.bool(True),
	## trigger ###################################
	printTriggerMenu = cms.untracked.bool(True),
	processName     = cms.string('HLT'),
	triggerName     = cms.vstring('HLT_PFJet260_v1'),
	triggerResults  = cms.InputTag("TriggerResults","","HLT"),
	triggerEvent    = cms.InputTag("hltTriggerSummaryAOD","","HLT"),
	## jec services ##############################
	#pfjecService    = cms.string('ak7PFL1FastL2L3Residual'),
	#calojecService  = cms.string('ak7CaloL1FastL2L3Residual')
)

process.ak7GenJets = process.ak5GenJets.clone()
process.ak7GenJets.rParam = cms.double(0.7)

process.ak7 = process.ak4.clone(
	pfjets          = cms.InputTag('selectedPatJetsAK7PFCHS'),
	pfjetschs       = cms.InputTag('selectedPatJetsAK7PFCHS'),
	## MET collection ####
	pfmet           = cms.InputTag('pfMet'),
	genjets         = cms.untracked.InputTag('ak7GenJets'),
)

process.ak5GenJets = process.ak5GenJets.clone()
process.ak5GenJets.rParam = cms.double(0.5)

process.ak5 = process.ak4.clone(
	pfjets          = cms.InputTag('selectedPatJetsAK5PFCHS'),
	pfjetschs       = cms.InputTag('selectedPatJetsAK5PFCHS'),
	## MET collection ####
	pfmet           = cms.InputTag('pfMet'),
	genjets         = cms.untracked.InputTag('ak5GenJets'),
)


process.p = cms.Path( process.ak4*process.ak5*process.ak7 )


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#! Output and Log
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )
process.options.allowUnscheduled = cms.untracked.bool(True)

process.output = cms.OutputModule("PoolOutputModule",                                                                                                                                                     
                                  #outputCommands = cms.untracked.vstring('drop *','keep *_puppi_*_*'),
                                  outputCommands = cms.untracked.vstring('keep *'),
                                  fileName       = cms.untracked.string ("Output.root")                                                                                                                   
)
# schedule definition                                                                                                       
process.outpath  = cms.EndPath(process.out) 

#!
#! THAT'S ALL! CAN YOU BELIEVE IT? :-D
#!
