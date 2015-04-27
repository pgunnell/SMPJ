# -*- coding: utf-8 -*-

import FWCore.ParameterSet.Config as cms

import FWCore.ParameterSet.Config as cms

from RecoJets.Configuration.RecoPFJets_cff import ak4PFJets, ak4PFJetsCHS
from RecoJets.Configuration.RecoGenJets_cff import ak4GenJets
from RecoJets.JetProducers.SubJetParameters_cfi import SubJetParameters
from RecoJets.JetProducers.PFJetParameters_cfi import *
from RecoJets.JetProducers.GenJetParameters_cfi import *
from RecoJets.JetProducers.AnomalousCellParameters_cfi import *
from RecoJets.JetProducers.CATopJetParameters_cfi import *
from PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff import *
from PhysicsTools.PatAlgos.selectionLayer1.jetSelector_cfi import selectedPatJets
from PhysicsTools.PatAlgos.tools.jetTools import addJetCollection


## Modified version of jetToolBox from https://github.com/cms-jet/jetToolbox
## Options for PUMethod: Puppi, CS, SK, CHS
def jetToolbox( proc, jetType, jetSequence,PUMethod=''):
	JETCorrPayload='None'
	JETCorrLevels = [ 'None' ]
	#JECLevels = [ 'L1Offset', 'L1FastJet', 'L1JPTOffset', 'L2Relative', 'L3Absolute', 'L5Falvour', 'L7Parton' ]

	algorithm='AntiKt' # CambridgeAachen' , 'Kt'
	size = jetType[-1:] #[-1:] takes the last char from string 'akX'
	jetSize = float('0.'+jetType[-1:])
	jetALGO = jetType.upper()
	jetalgo = jetType.lower()

	print 'Running processes with: '+str(jetALGO)+' PF '+PUMethod+' jet algorithm with radius parameter '+str(jetSize)

	JETCorrPayload = 'AK'+size+'PFchs'
	#JETCorrLevels = ['L1FastJet', 'L2Relative', 'L3Absolute']
	JETCorrLevels = [] #No JEC corrections
	JEC = ( JETCorrPayload, JETCorrLevels , 'None')

	#################################################################################
	####### Toolbox start 
	#################################################################################

	elemToKeep = []
	jetSeq = cms.Sequence()
	genParticlesLabel = ''
	pvLabel = ''
	tvLabel = ''
	toolsUsed = []

	#### For MiniAOD
	genParticlesLabel = 'prunedGenParticles'
	pvLabel = 'offlineSlimmedPrimaryVertices'
	svLabel = 'slimmedSecondaryVertices'
	tvLabel = 'unpackedTracksAndVertices'
	pfCand = 'packedPFCandidates'

	setattr( proc, 'chs', cms.EDFilter('CandPtrSelector', src = cms.InputTag('packedPFCandidates'), cut = cms.string('fromPV')) )
	jetSeq += getattr(proc, 'chs')


	### Filter out neutrinos from packed GenParticles
	#setattr( proc, 'packedGenParticlesForJets', 
			#cms.EDFilter("CandPtrSelector", 
				#src = cms.InputTag("packedGenParticles"), 
				#cut = cms.string("abs(pdgId) != 12 && abs(pdgId) != 14 && abs(pdgId) != 16")
				#))
	#jetSeq += getattr(proc, 'packedGenParticlesForJets' )
	    
	setattr( proc, jetalgo+'GenJets', 
			ak4GenJets.clone( src = 'packedGenParticles', #'packedGenParticlesForJets', 
				rParam = jetSize, 
				jetAlgorithm = algorithm ) ) 
	jetSeq += getattr(proc, jetalgo+'GenJets' )
	#fixedGridRhoFastjetAll.pfCandidatesTag = 'packedPFCandidates'

	#for Inclusive Vertex Finder
	proc.load('PhysicsTools.PatAlgos.slimming.unpackedTracksAndVertices_cfi')


	####  Creating PATjets
	if 'CHS' in PUMethod: 
		setattr( proc, jetalgo+'PFJetsCHS', 
				ak4PFJetsCHS.clone( 
					doAreaFastjet = True, 
					rParam = jetSize, 
					jetAlgorithm = algorithm ) ) 
		getattr( proc, jetalgo+'PFJetsCHS').src = 'chs'
		jetSeq += getattr(proc, jetalgo+'PFJetsCHS' )
	else: 
		PUMethod = ''
		setattr( proc, jetalgo+'PFJets', 
				ak4PFJets.clone( 
					doAreaFastjet = True, 
					rParam = jetSize, 
					jetAlgorithm = algorithm ) ) 
		getattr( proc, jetalgo+'PFJets').src = 'packedPFCandidates'
		jetSeq += getattr(proc, jetalgo+'PFJets' )

	addJetCollection(
			proc,
			labelName = jetALGO+'PF'+PUMethod,
			jetSource = cms.InputTag( jetalgo+'PFJets'+PUMethod),
			algo = jetalgo,
			rParam = jetSize,
			jetCorrections = JEC, #( 'AK'+size+'PFchs', cms.vstring( ['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
			pfCandidates = cms.InputTag( pfCand ),  #'packedPFCandidates'),
			svSource = cms.InputTag( svLabel ),   #'slimmedSecondaryVertices'),
			genJetCollection = cms.InputTag( jetalgo+'GenJets'),
			pvSource = cms.InputTag( pvLabel )#'offlineSlimmedPrimaryVertices'),
			) 

	getattr( proc, 'patJetCorrFactors'+jetALGO+'PF'+PUMethod ).primaryVertices = pvLabel  #'offlineSlimmedPrimaryVertices' 
	#getattr( proc, 'jetTracksAssociatorAtVertex'+jetALGO+'PF'+PUMethod ).tracks = tvLabel  # 'unpackedTracksAndVertices'

	getattr(proc,'patJetPartons').particles = cms.InputTag( genParticlesLabel ) #'prunedGenParticles')
	getattr(proc,'patJetPartonMatch'+jetALGO+'PF'+PUMethod).matched = cms.InputTag( genParticlesLabel ) #'prunedGenParticles')
	if hasattr(proc,'pfInclusiveSecondaryVertexFinderTagInfos'+jetALGO+'PF'+PUMethod):
		    getattr(proc,'pfInclusiveSecondaryVertexFinderTagInfos'+jetALGO+'PF'+PUMethod).extSVCollection = cms.InputTag( svLabel ) #'slimmedSecondaryVertices')
	getattr(proc,'patJets'+jetALGO+'PF'+PUMethod).addAssociatedTracks = cms.bool(False) # needs to be disabled since there is no track collection present in MiniAOD
	getattr(proc,'patJets'+jetALGO+'PF'+PUMethod).addJetCharge = cms.bool(False)        # needs to be disabled since there is no track collection present in MiniAOD
		
	elemToKeep += [ 'keep *_selectedPatJets'+jetALGO+'PF'+PUMethod+'_*_*' ]
	#elemToKeep += [ 'drop *_selectedPatJets'+jetALGO+'PF'+PUMethod+'_calo*_*' ]
	#elemToKeep += [ 'drop *_selectedPatJets'+jetALGO+'PF'+PUMethod+'_tagInfos_*' ]

	### "return"
	setattr(proc, jetSequence, jetSeq)
	#getattr(proc, outputFile).outputCommands += elemToKeep

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
'file:///mnt/storage/gflouris/12FCC45C-5D78-E411-BEB2-00261894388B.root'
#'file:///afs/cern.ch/work/g/gflouris/public/SMPJ_AnalysisFW/08C07BB6-376F-E411-BE9F-C4346BC7EE18.root'
   )
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))
process.source = cms.Source("PoolSource", fileNames = inFiles )

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#! Services
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.load('CommonTools.UtilAlgos.TFileService_cfi')
process.TFileService.fileName=cms.string('MC_ProcessedTreeProducer_miniAOD_2.root')

jetToolbox( process, 'ak4', 'ak4JetSubs', PUMethod='') 
jetToolbox( process, 'ak4', 'ak4JetSubs', PUMethod='CHS') 

process.ak4 =  cms.EDAnalyzer('ProcessedTreeProducer_miniAOD',
	## jet collections ###########################
	pfjets          = cms.InputTag('selectedPatJetsAK7PF'),
	pfjetschs       = cms.InputTag('selectedPatJetsAK7PFCHS'),
	## MET collection ####
	pfmet           = cms.InputTag('slimmedMETs'),
	genjets         = cms.untracked.InputTag('ak4GenJets'),
	## database entry for the uncertainties ######
	PFPayloadName   = cms.string(''),
	PFPayloadNameCHS= cms.string(''),
	CaloPayloadName = cms.string(''),
	jecUncSrc       = cms.string(''),
	jecUncSrcCHS    = cms.string(''),
	jecUncSrcNames  = cms.vstring(''),
	## set the conditions for good Vtx counting ##
	offlineVertices = cms.InputTag('offlineSlimmedPrimaryVertices'),
	goodVtxNdof     = cms.double(4),
	goodVtxZ        = cms.double(24),
	## rho #######################################
	srcCaloRho      = cms.InputTag('fixedGridRhoFastjetAll'),
	srcPFRho        = cms.InputTag('fixedGridRhoFastjetAllCalo'),
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


jetToolbox( process, 'ak8', 'ak8JetSubs','CHS')
jetToolbox( process, 'ak8', 'ak8JetSubs')

process.ak8GenJets = process.ak4GenJets.clone()
process.ak8GenJets.rParam = cms.double(0.8)

process.ak8 = process.ak4.clone(
	pfjets          = cms.InputTag('selectedPatJetsAK8PF'),
	pfjetschs       = cms.InputTag('selectedPatJetsAK8PFCHS'),
	genjets         = cms.untracked.InputTag('ak8GenJets'),
)
jetToolbox( process, 'ak7', 'ak7JetSubs', PUMethod='CHS') 
jetToolbox( process, 'ak7', 'ak7JetSubs', PUMethod='') 

process.ak7GenJets = process.ak4GenJets.clone()
process.ak7GenJets.rParam = cms.double(0.7)

process.ak7 = process.ak4.clone(
	pfjets          = cms.InputTag('selectedPatJetsAK7PF'),
	pfjetschs       = cms.InputTag('selectedPatJetsAK7PFCHS'),
	genjets         = cms.untracked.InputTag('ak7GenJets'),
)

jetToolbox( process, 'ak5', 'ak5JetSubs','CHS')
jetToolbox( process, 'ak5', 'ak5JetSubs')

process.ak5GenJets = process.ak5GenJets.clone()
process.ak5GenJets.rParam = cms.double(0.5)

process.ak5 = process.ak4.clone(
	pfjets          = cms.InputTag('selectedPatJetsAK5PFCHS'),
	pfjetschs       = cms.InputTag('selectedPatJetsAK5PFCHS'),
	genjets         = cms.untracked.InputTag('ak5GenJets'),
)


process.p = cms.Path( process.ak4*process.ak7*process.ak5*process.ak7*process.ak8 )


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#! Output and Log
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )
process.options.allowUnscheduled = cms.untracked.bool(True)
