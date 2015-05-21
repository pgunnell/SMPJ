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
from PhysicsTools.PatAlgos.tools.jetTools import *
from PhysicsTools.PatAlgos.producersLayer1.metProducer_cfi import patMETs
from PhysicsTools.PatAlgos.patSequences_cff import *
from PhysicsTools.PatAlgos.tools.metTools import *
from RecoJets.JetProducers.PileupJetID_cfi import *




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
	elemToKeep = []

	print 'Running processes with: '+str(jetALGO)+' PF '+PUMethod+' jet algorithm with radius parameter '+str(jetSize)

	JETCorrPayload = 'AK'+size+'PFchs'
	JETCorrLevels = ['L1FastJet', 'L2Relative', 'L3Absolute']
	#JETCorrLevels = [] #No JEC corrections
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

	proc.load('RecoJets.Configuration.GenJetParticles_cff')
	setattr( proc, jetalgo+'GenJetsNoNu', ak4GenJets.clone( src = 'genParticlesForJetsNoNu', rParam = jetSize, jetAlgorithm = algorithm ) ) 
	jetSeq += getattr(proc, jetalgo+'GenJetsNoNu' )



	proc.load('CommonTools.ParticleFlow.pfNoPileUpJME_cff')
	#setattr( proc, jetalgo+'GenJetsNoNu', ak4GenJets.clone( src = 'genParticlesForJetsNoNu', rParam = jetSize, jetAlgorithm = algorithm ) )
	#jetSeq += getattr(proc, jetalgo+'GenJetsNoNu' )
	####  Creating PATjets
	if( PUMethod=='CHS') :
	  setattr( proc, jetalgo+'PFJetsCHS', ak4PFJets.clone( rParam = jetSize, jetAlgorithm = algorithm ) )
	  jetSeq += getattr(proc, jetalgo+'PFJetsCHS' )

	  setattr( proc, jetalgo+'PFJetsCHS',
			  ak4PFJetsCHS.clone(
				  doAreaFastjet = True,
				  rParam = jetSize,
				  jetAlgorithm = algorithm ) )
	  jetSeq += getattr(proc, jetalgo+'PFJetsCHS' )

	else :
	  setattr( proc, jetalgo+'PFJets', ak4PFJets.clone( rParam = jetSize, jetAlgorithm = algorithm ) )
	  jetSeq += getattr(proc, jetalgo+'PFJets' )

	  setattr( proc, jetalgo+'PFJets',
			  ak4PFJets.clone(
				  doAreaFastjet = True,
				  rParam = jetSize,
				  jetAlgorithm = algorithm ) )
	  jetSeq += getattr(proc, jetalgo+'PFJets' )
	  PUMethod=''


	addJetCollection(
			proc,
			labelName = jetALGO+'PF'+PUMethod,
			jetSource = cms.InputTag( jetalgo+'PFJets'+PUMethod),
			algo = jetalgo,
			rParam = jetSize,
			jetCorrections =  JEC, #( 'AK'+size+'PFchs', cms.vstring( ['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
			pfCandidates = cms.InputTag( 'particleFlow' ),  #'packedPFCandidates'),
			svSource = cms.InputTag('slimmedSecondaryVertices'),
			genJetCollection = cms.InputTag( jetalgo+'GenJetsNoNu'),
			pvSource = cms.InputTag( 'offlinePrimaryVertices' ), #'offlineSlimmedPrimaryVertices'),
			jetTrackAssociation = True,

			#userFloats = cms.PSet( src = cms.VInputTag(cms.InputTag(jetALGO+"pileupJetId","fullDiscriminant"))),
			)


	process.load('RecoJets.JetProducers.pileupjetidproducer_cfi')
	proc.pileupJetIdCalculator.jets = cms.InputTag(jetalgo+"PFJets"+PUMethod)
	proc.pileupJetIdEvaluator.jets  = cms.InputTag(jetalgo+"PFJets"+PUMethod)
	#proc.pileupJetIdCalculator.jec  = cms.string(jetalgo+"PFJets"+PUMethod.lower())
	#proc.pileupJetIdEvaluator.jec   = cms.string(jetalgo+"PFJets"+PUMethod.lower())
	
	proc.patJetsAK4PFCHS.userData.userFloats.src += ['pileupJetId:fullDiscriminant']
	#getattr(proc, 'patJets'+jetALGO+'PF'+PUMethod).userData.userFloats.src += ['pileupJetId:fullDiscriminant']  

	elemToKeep += [ 'keep *_selectedPatJets'+jetALGO+'PF'+PUMethod+'_*_*' ]
	elemToKeep += [ 'keep *_patMetsPF*_*_*' ]
	elemToKeep +=  ['keep *_pileupJetId_*_*']

	getattr(proc,'patJetPartons').particles = cms.InputTag( 'genParticles' ) #'prunedGenParticles')
	setattr(proc, 'selectedPatJets'+jetALGO+'PF'+PUMethod, selectedPatJets.clone( src = 'patJets'+jetALGO+'PF'+PUMethod ) )
	setattr(proc, jetSequence, jetSeq)




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
process.load('RecoJets.JetProducers.TrackJetParameters_cfi')
process.load('RecoJets.JetProducers.PileupJetIDParams_cfi')

process.load("PhysicsTools.PatAlgos.slimming.pileupJetId_cfi")

#process.GlobalTag.globaltag = "PHYS14_25_V2::All" #Mad

process.GlobalTag.globaltag = "MCRUN2_74_V9A::All" #Pythia



#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#! Input
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

inFiles = cms.untracked.vstring(
#'file:///mnt/storage/gflouris/08C07BB6-376F-E411-BE9F-C4346BC7EE18.root' #Madgraph PHYS14
'file:///mnt/storage/gflouris/2430A1EC-00FA-E411-8641-0025905A7786.root' #Pythia
#'file:///afs/cern.ch/work/g/gflouris/public/SMPJ_AnalysisFW/08C07BB6-376F-E411-BE9F-C4346BC7EE18.root'
   )
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(200))
process.source = cms.Source("PoolSource", fileNames = inFiles )

jetToolbox( process, 'ak4', 'ak4JetSubs','CHS')
jetToolbox( process, 'ak4', 'ak4JetSubs')


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#! Services
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.load('CommonTools.UtilAlgos.TFileService_cfi')
process.TFileService.fileName=cms.string('MC_ProcessedTreeProducer_2.root')


process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
process.load("PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff")

addMETCollection(process,'patMETPF','pfMetT1')

process.ak4 =  cms.EDAnalyzer('ProcessedTreeProducer',
	## jet collections ###########################
	pfjets          = cms.InputTag('selectedPatJetsAK4PF'),
	pfjetschs       = cms.InputTag('selectedPatJetsAK4PFCHS'),
	## MET collection ####
	pfmet           = cms.InputTag('patMETs'),
	genjets         = cms.untracked.InputTag('ak4GenJetsNoNu'),
	## database entry for the uncertainties ######
	PFPayloadName   = cms.string(''),
	PFPayloadNameCHS= cms.string(''),
	jecUncSrc       = cms.string(''),
	jecUncSrcCHS    = cms.string(''),
	jecUncSrcNames  = cms.vstring(''),
	## set the conditions for good Vtx counting ##
	offlineVertices = cms.InputTag('offlinePrimaryVertices'),
	goodVtxNdof     = cms.double(4),
	goodVtxZ        = cms.double(24),
	## rho #######################################
	srcCaloRho      = cms.InputTag('fixedGridRhoFastjetAllCalo'),
	srcPFRho        = cms.InputTag('fixedGridRhoFastjetAll'),
	srcPU           = cms.untracked.InputTag('addPileupInfo'),
	## preselection cuts #########################
	maxY            = cms.double(5.0),
	minPFPt         = cms.double(20),
	minNPFJets      = cms.int32(1),
	minGenPt        = cms.untracked.double(20),
	minJJMass       = cms.double(-1),
	isMCarlo        = cms.untracked.bool(True),
	useGenInfo      = cms.untracked.bool(True),
	## trigger ###################################
	printTriggerMenu = cms.untracked.bool(True),
	processName     = cms.string('HLT'),
	triggerName     = cms.vstring('HLT_PFJet40_v1','HLT_PFJet60_v1', 'HLT_PFJet80_v1', 'HLT_PFJet140_v1', 'HLT_PFJet200_v1', 'HLT_PFJet260_v1', 
				      'HLT_PFJet320_v1', 'HLT_PFJet400_v1', 'HLT_PFJet450_v1', 'HLT_PFJet500_v1'
				      ),
	triggerResults  = cms.InputTag("TriggerResults","","HLT"),
	triggerEvent    = cms.InputTag("hltTriggerSummaryAOD","","HLT"),
	## jec services ##############################
	#pfjecService    = cms.string('ak7PFL1FastL2L3Residual'),
)


jetToolbox( process, 'ak8', 'ak8JetSubs','CHS')
jetToolbox( process, 'ak8', 'ak8JetSubs')

process.ak8 = process.ak4.clone(
	pfjets          = cms.InputTag('selectedPatJetsAK8PF'),
	pfjetschs       = cms.InputTag('selectedPatJetsAK8PFCHS'),
	## MET collection ####
	pfmet           = cms.InputTag('patMETs'),
	genjets         = cms.untracked.InputTag('ak8GenJets'),
)

jetToolbox( process, 'ak7', 'ak5JetSubs','CHS')
jetToolbox( process, 'ak7', 'ak7JetSubs')

process.ak7GenJets = process.ak5GenJets.clone()
process.ak7GenJets.rParam = cms.double(0.7)

process.ak7 = process.ak4.clone(
	pfjets          = cms.InputTag('selectedPatJetsAK7PF'),
	pfjetschs       = cms.InputTag('selectedPatJetsAK7PFCHS'),
	## MET collection ####
	pfmet           = cms.InputTag('patMETs'),
	genjets         = cms.untracked.InputTag('ak7GenJetsNoNu'),
)

jetToolbox( process, 'ak5', 'ak5JetSubs','CHS')
jetToolbox( process, 'ak5', 'ak5JetSubs')

process.ak5GenJets = process.ak5GenJets.clone()
process.ak5GenJets.rParam = cms.double(0.5)

process.ak5 = process.ak4.clone(
	pfjets          = cms.InputTag('selectedPatJetsAK5PF'),
	pfjetschs       = cms.InputTag('selectedPatJetsAK5PFCHS'),
	## MET collection ####
	pfmet           = cms.InputTag('patMETs'),
	genjets         = cms.untracked.InputTag('ak5GenJetsNoNu'),
)


process.p = cms.Path( process.ak4*process.ak5*process.ak7*process.ak8 )


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#! Output and Log
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )
process.options.allowUnscheduled = cms.untracked.bool(True)
