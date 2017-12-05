# -*- coding: utf-8 -*-

import FWCore.ParameterSet.Config as cms

from RecoJets.Configuration.RecoPFJets_cff import ak4PFJets, ak4PFJetsCHS
from RecoJets.Configuration.RecoGenJets_cff import ak4GenJets
from RecoJets.JetProducers.SubJetParameters_cfi import SubJetParameters
from RecoJets.JetProducers.PFJetParameters_cfi import *
from RecoJets.JetProducers.GenJetParameters_cfi import *
from RecoJets.JetProducers.AnomalousCellParameters_cfi import *
from PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff import *
from PhysicsTools.PatAlgos.selectionLayer1.jetSelector_cfi import selectedPatJets
from PhysicsTools.PatAlgos.tools.jetTools import *
#from PhysicsTools.PatAlgos.producersLayer1.metProducer_cfi import patMETs
from PhysicsTools.PatAlgos.patSequences_cff import *
#nnfrom PhysicsTools.PatAlgos.tools.metTools import addMETCollection
#from RecoJets.JetProducers.pileupjetidproducer_cfi import *
from PhysicsTools.PatAlgos.patTemplate_cfg import *
from PhysicsTools.PatAlgos.tools.jetTools import *
from RecoJets.JetProducers.QGTagger_cfi import QGTagger

## Modified version of jetToolBox from https://github.com/cms-jet/jetToolbox
## Options for PUMethod: Puppi, CS, SK, CHS

# -*- coding: utf-8 -*-
import FWCore.ParameterSet.Config as cms

process = cms.Process("Ntuplizer")
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#! Conditions
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('RecoJets.Configuration.GenJetParticles_cff')
process.load('RecoJets.Configuration.RecoGenJets_cff')
process.load('RecoJets.JetProducers.TrackJetParameters_cfi')
process.load('RecoJets.JetProducers.PileupJetIDParams_cfi')
process.load("PhysicsTools.PatAlgos.patSequences_cff")

#process.load("PhysicsTools.PatAlgos.slimming.pileupJetId_cfi")
process.load('RecoJets.JetProducers.QGTagger_cfi')
process.QGTagger.srcJets = cms.InputTag('ak4PFJetsCHS')    # Could be reco::PFJetCollection or pat::JetCollection (both AOD and miniAOD)
process.QGTagger.jetsLabel  = cms.string('QGL_AK4PFchs')        # Other options: see https://twiki.cern.ch/twiki/bin/viewauth/CMS/QGDataBaseVersion

#process.GlobalTag.globaltag = "80X_dataRun2_ICHEP16_repro_v0"
process.GlobalTag.globaltag = "80X_dataRun2_2016SeptRepro_v7" #ReReco RunB-G
#process.GlobalTag.globaltag = "80X_dataRun2_Prompt_v16"

process.GlobalTag.toGet = cms.VPSet(
        cms.PSet(record = cms.string("BTagTrackProbability3DRcd"),
            tag = cms.string("JPcalib_Data80X_2016_v3"),
            connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS")
            )
        )



##-------------------- Import the JEC services -----------------------
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#! Input
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

inFiles = cms.untracked.vstring(
#'root://cms-xrd-global.cern.ch//store/data/Run2016G/JetHT/AOD/23Sep2016-v1/100000/0645BD20-F486-E611-A724-002590D0B054.root'
'root://cms-xrd-global.cern.ch//store/data/Run2016G/JetHT/AOD/23Sep2016-v1/100000/0006CE1E-9986-E611-8DFB-6C3BE5B5C0B0.root'
#'root://cms-xrd-global.cern.ch//store/data/Run2016G/JetHT/AOD/18Apr2017-v1/100000/001D282A-9134-E711-9003-0090FAA57780.root'
#'0002DAEF-97EA-E611-8DEC-001E67E69DEC.root'
)

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))
process.source = cms.Source("PoolSource", fileNames = inFiles )

def jetToolbox( proc, jetType, jetSequence,PUMethod='', bTagDiscriminators = None):

	JETCorrPayload='None'
	JETCorrLevels = [ 'None' ]
	bTagDiscriminators = [#'trackCountingHighEffBJetTags',
						  #'trackCountingHighPurBJetTags',
						  #'pfTrackCountingHighEffBJetTags',
						  #'pfTrackCountingHighPurBJetTags',
						  #'softPFMuonByIP3dBJetTags',
						  #'softPFElectronByIP3dBJetTags',
						  #'softPFMuonBJetTags',
						  #'softPFElectronBJetTags',
						  #'simpleSecondaryVertexHighEffBJetTags',
						  #'simpleSecondaryVertexHighPurBJetTags',
						  #'pfSimpleSecondaryVertexHighEffBJetTags',
						  #'pfSimpleSecondaryVertexHighPurBJetTags',
						  #'combinedSecondaryVertexV2BJetTags',
						  'deepFlavourJetTags:probb',
						  'deepFlavourJetTags:probc',
						  'deepFlavourJetTags:probudsg',
						  'deepFlavourJetTags:probbb',
						  'deepFlavourJetTags:probcc',
						  'negativeDeepFlavourJetTags:probb',
						  'negativeDeepFlavourJetTags:probc',
						  'negativeDeepFlavourJetTags:probudsg',
						  'negativeDeepFlavourJetTags:probbb',
						  'negativeDeepFlavourJetTags:probcc',
						  'positiveDeepFlavourJetTags:probb',
						  'positiveDeepFlavourJetTags:probc',
						  'positiveDeepFlavourJetTags:probudsg',
						  'positiveDeepFlavourJetTags:probbb',
						  'positiveDeepFlavourJetTags:probcc',
						  'pfCombinedCvsLJetTags',
						  'pfCombinedCvsBJetTags',
						  'pfBoostedDoubleSecondaryVertexAK8BJetTags',
						  'pfCombinedSecondaryVertexV2BJetTags',
						  'pfPositiveCombinedSecondaryVertexV2BJetTags',  #implemented
						  'pfNegativeCombinedSecondaryVertexV2BJetTags',  #implemented
						  'pfCombinedInclusiveSecondaryVertexV2BJetTags', #implemented
						  'pfCombinedMVAV2BJetTags',                      #implemented
						  'pfJetProbabilityBJetTags']                     #implemented

	algorithm='AntiKt' # CambridgeAachen' , 'Kt'
	size = jetType[-1:] #[-1:] takes the last char from string 'akX'
	jetSize = float('0.'+jetType[-1:])
	jetALGO = jetType.upper()
	jetalgo = jetType.lower()
	elemToKeep = []

	print 'Running processes with: '+str(jetALGO)+' PF '+PUMethod+' jet algorithm with radius parameter '+str(jetSize)

	JETCorrPayload = 'AK'+size+'PF'+PUMethod.lower()
	JETCorrLevels = ['L1FastJet', 'L2Relative', 'L3Absolute','L2L3Residual']
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
			btagDiscriminators = bTagDiscriminators,
			rParam = jetSize,
			jetCorrections =  JEC, #( 'AK'+size+'PFchs', cms.vstring( ['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
			pfCandidates = cms.InputTag( 'particleFlow' ),  #'packedPFCandidates'),
			svSource = cms.InputTag('inclusiveCandidateSecondaryVertices'),
			genJetCollection = cms.InputTag( jetalgo+'GenJetsNoNu'),
			pvSource = cms.InputTag( 'offlinePrimaryVertices' ), #'offlineSlimmedPrimaryVertices'),
			jetTrackAssociation = True,
			)

	getattr(proc,'patJetPartons').particles = cms.InputTag( 'genParticles' ) #'prunedGenParticles')

	QGjetsLabel='chs'

	setattr( proc, 'QGTagger'+jetALGO+'PF'+PUMethod,
		 QGTagger.clone(
			srcJets = cms.InputTag(jetalgo+'PFJets'+PUMethod),    # Could be reco::PFJetCollection or pat::JetCollection (both AOD and miniAOD)
			jetsLabel = cms.string('QGL_AK4PF'+QGjetsLabel)        # Other options (might need to add an ESSource for it): see https://twiki.cern.ch/twiki/bin/viewauth/CMS/QGDataBaseVersion
			)
		 )
	elemToKeep += [ 'keep *_QGTagger'+jetALGO+'PF'+PUMethod+'_*_*' ]
	getattr( proc, 'patJets'+jetALGO+'PF'+PUMethod).userData.userFloats.src += ['QGTagger'+jetALGO+'PF'+PUMethod+':qgLikelihood']
	jetSeq += getattr(proc, 'QGTagger'+jetALGO+'PF'+PUMethod )

	toolsUsed.append( 'QGTagger'+jetALGO+'PF'+PUMethod )

	#getattr( proc, 'selectedPatJets').userData.userFloats.src += ['QGTagger:qgLikelihood']
	#setattr(proc, 'selectedPatJets'+jetALGO+'PF'+PUMethod, selectedPatJets.clone( src = 'patJets'+jetALGO+'PF'+PUMethod ) )
	#jetSeq += getattr(proc, 'QGTagger')

	setattr(proc, jetSequence, jetSeq)




#jetToolbox( process, 'ak4', 'ak4JetSubs')
jetToolbox( process, 'ak4', 'ak4JetSubs','CHS')

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#! Services
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.load('CommonTools.UtilAlgos.TFileService_cfi')
process.TFileService.fileName=cms.string('DATA_ProcessedTreeProducer_2.root')

# PAT Layer 1
#process.load("PhysicsTools.PatAlgos.patLayer0_cff") # need to load this
#process.load("PhysicsTools.PatAlgos.patLayer1_cff") # even if we run only layer 1

process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
process.load("PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff")
#from PhysicsTools.PatAlgos.tools.metTools import addMETCollection
#addMETCollection(process,'patMETPF','pfMetT1')
#patMETs.addGenMet = cms.bool(False)
#patMETs.genMETSource = cms.InputTag("")
process.patJets.addTagInfos = True
process.patJets.addAssociatedTracks = True

#process.out.outputCommands += ['keep *_QGTagger_*_*']
process.load("PhysicsTools.PatAlgos.patSequences_cff")

process.ak4 =  cms.EDAnalyzer('ProcessedTreeProducerBTag',
	## jet collections ###########################
	pfjetschs       = cms.InputTag('selectedPatJetsAK4PFCHS'),
	pfpujetid       = cms.string('AK4PFpileupJetIdEvaluator:fullDiscriminant'),
	pfchsjetpuid    = cms.string('AK4PFCHSpileupJetIdEvaluator:fullDiscriminant'),
	## MET collection ####
	pfmet           = cms.InputTag('patMETs'),
	genjets         = cms.untracked.InputTag('ak4GenJetsNoNu'),
	## database entry for the uncertainties ######
	PFPayloadNameCHS= cms.string('AK4PFchs'),
	jecUncSrcCHS    = cms.string(''),
	jecUncSrcNames  = cms.vstring(''),
	## set the conditions for good Vtx counting ##
	offlineVertices = cms.InputTag('offlinePrimaryVertices'),
	beamSpot        = cms.InputTag('offlineBeamSpot'),
	goodVtxNdof     = cms.double(4),
	goodVtxZ        = cms.double(24),
	## rho #######################################
	srcCaloRho      = cms.InputTag('fixedGridRhoFastjetAllCalo'),
	srcPFRho        = cms.InputTag('fixedGridRhoFastjetAll'),
	srcPULabel      = cms.untracked.InputTag('addPileupInfo','addPileupInfo'),
	## preselection cuts #########################
	maxY            = cms.double(5.0),
	minPFPt         = cms.double(64.0),
	minNPFJets      = cms.int32(1),
	minGenPt        = cms.untracked.double(20),
	minJJMass       = cms.double(-1),
	isMCarlo        = cms.untracked.bool(False),
	useGenInfo      = cms.untracked.bool(False),
	AK4             = cms.untracked.bool(True),
	## trigger ###################################
	printTriggerMenu = cms.untracked.bool(True),
	processName     = cms.string('HLT'),
	triggerName     =  cms.vstring(
		'HLT_PFJet40_v9','HLT_PFJet60_v9','HLT_PFJet80_v9','HLT_PFJet140_v9','HLT_PFJet200_v9','HLT_PFJet260_v9','HLT_PFJet320_v9','HLT_PFJet400_v9','HLT_PFJet450_v9','HLT_PFJet500_v9',
		'HLT_PFJet40_v8','HLT_PFJet60_v8','HLT_PFJet80_v8','HLT_PFJet140_v8','HLT_PFJet200_v8','HLT_PFJet260_v8','HLT_PFJet320_v8','HLT_PFJet400_v8','HLT_PFJet450_v8','HLT_PFJet500_v8',
		'HLT_PFJet40_v7','HLT_PFJet60_v7','HLT_PFJet80_v7','HLT_PFJet140_v7','HLT_PFJet200_v7','HLT_PFJet260_v7','HLT_PFJet320_v7','HLT_PFJet400_v7','HLT_PFJet450_v7','HLT_PFJet500_v7',
		'HLT_PFJet40_v6','HLT_PFJet60_v6','HLT_PFJet80_v6','HLT_PFJet140_v6','HLT_PFJet200_v6','HLT_PFJet260_v6','HLT_PFJet320_v6','HLT_PFJet400_v6','HLT_PFJet450_v6','HLT_PFJet500_v6',
		'HLT_PFJet40_v5','HLT_PFJet60_v5','HLT_PFJet80_v5','HLT_PFJet140_v5','HLT_PFJet200_v5','HLT_PFJet260_v5','HLT_PFJet320_v5','HLT_PFJet400_v5','HLT_PFJet450_v5','HLT_PFJet500_v5',
		'HLT_PFJet40_v4','HLT_PFJet60_v4', 'HLT_PFJet80_v4', 'HLT_PFJet140_v4','HLT_PFJet200_v4','HLT_PFJet260_v4','HLT_PFJet320_v4','HLT_PFJet400_v4','HLT_PFJet450_v4','HLT_PFJet500_v4',
		'HLT_PFHT125_v2','HLT_PFHT200_v2','HLT_PFHT250_v2','HLT_PFHT300_v2','HLT_PFHT350_v2','HLT_PFHT400_v2','HLT_PFHT475_v2','HLT_PFHT600_v2','HLT_PFHT650_v2','HLT_PFHT800_v2','HLT_PFHT900_v2',
		'HLT_PFHT125_v3','HLT_PFHT200_v3','HLT_PFHT250_v3','HLT_PFHT300_v3','HLT_PFHT350_v3','HLT_PFHT400_v3','HLT_PFHT475_v3','HLT_PFHT600_v3','HLT_PFHT650_v3','HLT_PFHT800_v3','HLT_PFHT900_v3',
		'HLT_PFHT125_v4','HLT_PFHT200_v4','HLT_PFHT250_v4','HLT_PFHT300_v4','HLT_PFHT350_v4','HLT_PFHT400_v4','HLT_PFHT475_v4','HLT_PFHT600_v4','HLT_PFHT650_v4','HLT_PFHT800_v4','HLT_PFHT900_v4',
		'HLT_PFHT125_v5','HLT_PFHT200_v5','HLT_PFHT250_v5','HLT_PFHT300_v5','HLT_PFHT350_v5','HLT_PFHT400_v5','HLT_PFHT475_v5','HLT_PFHT600_v5','HLT_PFHT650_v5','HLT_PFHT800_v5','HLT_PFHT900_v5',
		'HLT_PFHT125_v1','HLT_PFHT200_v1','HLT_PFHT250_v1','HLT_PFHT300_v1','HLT_PFHT350_v1','HLT_PFHT400_v1','HLT_PFHT475_v1','HLT_PFHT600_v1','HLT_PFHT650_v1','HLT_PFHT800_v1','HLT_PFHT900_v1',
		'HLT_AK8PFJet40_v1','HLT_AK8PFJet60_v1','HLT_AK8PFJet80_v1','HLT_AK8PFJet140_v1','HLT_AK8PFJet200_v1','HLT_AK8PFJet260_v1','HLT_AK8PFJet320_v1','HLT_AK8PFJet400_v1','HLT_AK8PFJet450_v1','HLT_AK8PFJet500_v1',
		'HLT_AK8PFJet40_v2','HLT_AK8PFJet60_v2','HLT_AK8PFJet80_v2','HLT_AK8PFJet140_v2','HLT_AK8PFJet200_v2','HLT_AK8PFJet260_v2','HLT_AK8PFJet320_v2','HLT_AK8PFJet400_v2','HLT_AK8PFJet450_v2','HLT_AK8PFJet500_v2',
		'HLT_AK8PFJet40_v3','HLT_AK8PFJet60_v3','HLT_AK8PFJet80_v3','HLT_AK8PFJet140_v3','HLT_AK8PFJet200_v3','HLT_AK8PFJet260_v3','HLT_AK8PFJet320_v3','HLT_AK8PFJet400_v3','HLT_AK8PFJet450_v3','HLT_AK8PFJet500_v3',
		'HLT_AK8PFJet40_v4','HLT_AK8PFJet60_v4','HLT_AK8PFJet80_v4','HLT_AK8PFJet140_v4','HLT_AK8PFJet200_v4','HLT_AK8PFJet260_v4','HLT_AK8PFJet320_v4','HLT_AK8PFJet400_v4','HLT_AK8PFJet450_v4','HLT_AK8PFJet500_v4',
		'HLT_AK8PFJet40_v5','HLT_AK8PFJet60_v5','HLT_AK8PFJet80_v5','HLT_AK8PFJet140_v5','HLT_AK8PFJet200_v5','HLT_AK8PFJet260_v5','HLT_AK8PFJet320_v5','HLT_AK8PFJet400_v5','HLT_AK8PFJet450_v5','HLT_AK8PFJet500_v5'),
	triggerResults  = cms.InputTag("TriggerResults","","HLT"),
	triggerEvent    = cms.InputTag("hltTriggerSummaryAOD","","HLT"),
	## jec services ##############################
	#new tokens
		EventInfo       = cms.InputTag("generator"),
		GenParticles    = cms.InputTag("genparticles"),
		HBHENoiseFilterResultLabel = cms.InputTag("HBHENoiseFilterResultProducer", "HBHENoiseFilterResult"),
		HBHENoiseFilterResultNoMinZLabel = cms.InputTag("HBHENoiseFilterResultProducerNoMinZ", "HBHENoiseFilterResult"),
				  jetFlavourInfos = cms.InputTag("genJetFlavourInfos"),
	#pfjecService    = cms.string('ak7PFL1FastL2L3Residual'),
)

#jetToolbox( process, 'ak8', 'ak8JetSubs','CHS')

#process.ak8 = process.ak4.clone(
#	pfjetschs       = cms.InputTag('selectedPatJetsAK8PFCHS'),
#	pfpujetid       = cms.string('AK8PFpileupJetIdEvaluator:fullDiscriminant'),
#	pfchsjetpuid    = cms.string('AK8PFCHSpileupJetIdEvaluator:fullDiscriminant'),
#    PFPayloadNameCHS= cms.string('AK8PFchs'),
#	AK4             = cms.untracked.bool(False),
#)

jetToolbox( process, 'ak7', 'ak7JetSubs','CHS')

process.ak7 = process.ak4.clone(
	pfjetschs       = cms.InputTag('selectedPatJetsAK7PFCHS'),
	pfpujetid       = cms.string('AK7PFpileupJetIdEvaluator:fullDiscriminant'),
	pfchsjetpuid    = cms.string('AK7PFCHSpileupJetIdEvaluator:fullDiscriminant'),
	PFPayloadNameCHS= cms.string('AK7PFchs'),
	AK4             = cms.untracked.bool(False),
)

#jetToolbox( process, 'ak2', 'ak2JetSubs','CHS')

#process.ak2 = process.ak4.clone(
#	pfjetschs       = cms.InputTag('selectedPatJetsAK4PFCHS'),
#	pfpujetid       = cms.string('AK4PFpileupJetIdEvaluator:fullDiscriminant'),
#	pfchsjetpuid    = cms.string('AK4PFCHSpileupJetIdEvaluator:fullDiscriminant'),
#    PFPayloadNameCHS= cms.string('AK4PFchs'),
#	AK4             = cms.untracked.bool(False),
#)

############# hlt filter #########################
process.hltFilter = cms.EDFilter('HLTHighLevel',
	TriggerResultsTag  = cms.InputTag('TriggerResults','','HLT'),
	HLTPaths  =  cms.vstring(
		'HLT_PFJet40_v9','HLT_PFJet60_v9','HLT_PFJet80_v9','HLT_PFJet140_v9','HLT_PFJet200_v9','HLT_PFJet260_v9','HLT_PFJet320_v9','HLT_PFJet400_v9','HLT_PFJet450_v9','HLT_PFJet500_v9',
		'HLT_PFJet40_v8','HLT_PFJet60_v8','HLT_PFJet80_v8','HLT_PFJet140_v8','HLT_PFJet200_v8','HLT_PFJet260_v8','HLT_PFJet320_v8','HLT_PFJet400_v8','HLT_PFJet450_v8','HLT_PFJet500_v8',
		'HLT_PFJet40_v7','HLT_PFJet60_v7','HLT_PFJet80_v7','HLT_PFJet140_v7','HLT_PFJet200_v7','HLT_PFJet260_v7','HLT_PFJet320_v7','HLT_PFJet400_v7','HLT_PFJet450_v7','HLT_PFJet500_v7',
		'HLT_PFJet40_v6','HLT_PFJet60_v6','HLT_PFJet80_v6','HLT_PFJet140_v6','HLT_PFJet200_v6','HLT_PFJet260_v6','HLT_PFJet320_v6','HLT_PFJet400_v6','HLT_PFJet450_v6','HLT_PFJet500_v6',
		'HLT_PFJet40_v5','HLT_PFJet60_v5','HLT_PFJet80_v5','HLT_PFJet140_v5','HLT_PFJet200_v5','HLT_PFJet260_v5','HLT_PFJet320_v5','HLT_PFJet400_v5','HLT_PFJet450_v5','HLT_PFJet500_v5',
		'HLT_PFJet40_v4','HLT_PFJet60_v4', 'HLT_PFJet80_v4', 'HLT_PFJet140_v4','HLT_PFJet200_v4','HLT_PFJet260_v4','HLT_PFJet320_v4','HLT_PFJet400_v4','HLT_PFJet450_v4','HLT_PFJet500_v4',
		'HLT_PFHT125_v2','HLT_PFHT200_v2','HLT_PFHT250_v2','HLT_PFHT300_v2','HLT_PFHT350_v2','HLT_PFHT400_v2','HLT_PFHT475_v2','HLT_PFHT600_v2','HLT_PFHT650_v2','HLT_PFHT800_v2','HLT_PFHT900_v2',
		'HLT_PFHT125_v3','HLT_PFHT200_v3','HLT_PFHT250_v3','HLT_PFHT300_v3','HLT_PFHT350_v3','HLT_PFHT400_v3','HLT_PFHT475_v3','HLT_PFHT600_v3','HLT_PFHT650_v3','HLT_PFHT800_v3','HLT_PFHT900_v3',
		'HLT_PFHT125_v4','HLT_PFHT200_v4','HLT_PFHT250_v4','HLT_PFHT300_v4','HLT_PFHT350_v4','HLT_PFHT400_v4','HLT_PFHT475_v4','HLT_PFHT600_v4','HLT_PFHT650_v4','HLT_PFHT800_v4','HLT_PFHT900_v4',
		'HLT_PFHT125_v5','HLT_PFHT200_v5','HLT_PFHT250_v5','HLT_PFHT300_v5','HLT_PFHT350_v5','HLT_PFHT400_v5','HLT_PFHT475_v5','HLT_PFHT600_v5','HLT_PFHT650_v5','HLT_PFHT800_v5','HLT_PFHT900_v5',
		'HLT_PFHT125_v1','HLT_PFHT200_v1','HLT_PFHT250_v1','HLT_PFHT300_v1','HLT_PFHT350_v1','HLT_PFHT400_v1','HLT_PFHT475_v1','HLT_PFHT600_v1','HLT_PFHT650_v1','HLT_PFHT800_v1','HLT_PFHT900_v1',
		'HLT_AK8PFJet40_v1','HLT_AK8PFJet60_v1','HLT_AK8PFJet80_v1','HLT_AK8PFJet140_v1','HLT_AK8PFJet200_v1','HLT_AK8PFJet260_v1','HLT_AK8PFJet320_v1','HLT_AK8PFJet400_v1','HLT_AK8PFJet450_v1','HLT_AK8PFJet500_v1',
		'HLT_AK8PFJet40_v2','HLT_AK8PFJet60_v2','HLT_AK8PFJet80_v2','HLT_AK8PFJet140_v2','HLT_AK8PFJet200_v2','HLT_AK8PFJet260_v2','HLT_AK8PFJet320_v2','HLT_AK8PFJet400_v2','HLT_AK8PFJet450_v2','HLT_AK8PFJet500_v2',
		'HLT_AK8PFJet40_v3','HLT_AK8PFJet60_v3','HLT_AK8PFJet80_v3','HLT_AK8PFJet140_v3','HLT_AK8PFJet200_v3','HLT_AK8PFJet260_v3','HLT_AK8PFJet320_v3','HLT_AK8PFJet400_v3','HLT_AK8PFJet450_v3','HLT_AK8PFJet500_v3',
		'HLT_AK8PFJet40_v4','HLT_AK8PFJet60_v4','HLT_AK8PFJet80_v4','HLT_AK8PFJet140_v4','HLT_AK8PFJet200_v4','HLT_AK8PFJet260_v4','HLT_AK8PFJet320_v4','HLT_AK8PFJet400_v4','HLT_AK8PFJet450_v4','HLT_AK8PFJet500_v4',
		'HLT_AK8PFJet40_v5','HLT_AK8PFJet60_v5','HLT_AK8PFJet80_v5','HLT_AK8PFJet140_v5','HLT_AK8PFJet200_v5','HLT_AK8PFJet260_v5','HLT_AK8PFJet320_v5','HLT_AK8PFJet400_v5','HLT_AK8PFJet450_v5','HLT_AK8PFJet500_v5'),
	eventSetupPathsKey = cms.string(''),
	andOr              = cms.bool(True), #----- True = OR, False = AND between the HLTPaths
	throw              = cms.bool(False)
)

##MET Filters
process.load('RecoMET.METFilters.CSCTightHaloFilter_cfi')
process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')
process.load('CommonTools.RecoAlgos.HBHENoiseFilter_cfi')

process.HBHENoiseFilterResultProducerNoMinZ = process.HBHENoiseFilterResultProducer.clone(minZeros = cms.int32(99999))

process.goodVertices = cms.EDFilter("VertexSelector",
   filter = cms.bool(False),
   src = cms.InputTag("offlinePrimaryVertices"),
   cut = cms.string("!isFake && ndof >= 4 && abs(z) <= 24 && position.rho <= 2"),
)

process.load('RecoMET.METFilters.trackingFailureFilter_cfi')
process.load('RecoMET.METFilters.EcalDeadCellTriggerPrimitiveFilter_cfi')
process.load('RecoMET.METFilters.eeBadScFilter_cfi')
process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')
process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
process.load('RecoMET.METFilters.globalTightHalo2016Filter_cfi')

process.allMetFilterPaths=cms.Sequence(process.HBHENoiseFilter*process.HBHENoiseIsoFilter*process.EcalDeadCellTriggerPrimitiveFilter*process.eeBadScFilter*process.BadPFMuonFilter*process.BadChargedCandidateFilter*process.globalTightHalo2016Filter)

##Type1 patMET Producer
process.load('PhysicsTools.PatAlgos.recoLayer0.metCorrections_cff')
process.load('PhysicsTools.PatAlgos.producersLayer1.metProducer_cfi')
process.patMETs.addGenMET = cms.bool(False)
#process.patMETs.metSource = cms.InputTag("pfChMet")

#Try scheduled processs
process.path = cms.Path(process.goodVertices*process.trackingFailureFilter*
			process.hltFilter*
			process.HBHENoiseFilterResultProducer*
				process.HBHENoiseFilterResultProducerNoMinZ*process.allMetFilterPaths*
			process.patMETCorrections*process.patMETs
			*process.QGTagger*process.ak4*process.ak7)


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#! Output and Log
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )
process.options.allowUnscheduled = cms.untracked.bool(True)
