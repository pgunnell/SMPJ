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
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:mc', '')
process.GlobalTag.globaltag = "80X_mcRun2_asymptotic_2016_TrancheIV_v6"
#process.GlobalTag.globaltag = "80X_mcRun2_asymptotic_2016_miniAODv2_v1"

process.load('RecoJets.JetProducers.QGTagger_cfi')
process.QGTagger.srcJets = cms.InputTag('ak4PFJetsCHS')    # Could be reco::PFJetCollection or pat::JetCollection (both AOD and miniAOD)
process.QGTagger.jetsLabel  = cms.string('QGL_AK4PFchs')        # Other options: see https://twiki.cern.ch/twiki/bin/viewauth/CMS/QGDataBaseVersion

##-------------------- Import the JEC services -----------------------
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#! Input
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

inFiles = cms.untracked.vstring(
#'file:00430D8B-320F-E611-BC64-02163E00EA8C.root'
#/afs/cern.ch/work/e/eeren/public/04E758F7-34F9-E411-AFF0-002618FDA259.root'
#'root://eoscms.cern.ch//eos/cms/store/mc/RunIISpring16DR80/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/AODSIM/PUSpring16_magnetOn_80X_mcRun2_asymptotic_2016_v3-v1/00000/1276C35E-EA0F-E611-8336-02163E00F718.root'
#'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/QCD_Pt_1400to1800_TuneCUETP8M1_13TeV_pythia8/AODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/60000/381409B5-9B08-E611-98F3-0025905B85A0.root'
'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16DR80Premix/QCD_Pt_300to470_TuneCUETP8M1_13TeV_pythia8/AODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/110000/02CD261D-10A9-E611-8521-FA163E5CAE75.root'
   )

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(700))
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

	GetJetMCFlavour = ['True']
        #JECLevels = [ 'L1Offset', 'L1FastJet', 'L1JPTOffset', 'L2Relative', 'L3Absolute', 'L5Falvour', 'L7Parton' ]

	algorithm='AntiKt' # CambridgeAachen' , 'Kt'
	size = jetType[-1:] #[-1:] takes the last char from string 'akX'
	jetSize = float('0.'+jetType[-1:])
	jetALGO = jetType.upper()
	jetalgo = jetType.lower()
	elemToKeep = []

	print 'Running processes with: '+str(jetALGO)+' PF '+PUMethod+' jet algorithm with radius parameter '+str(jetSize)

	JETCorrPayload = 'AK'+size+'PF'+PUMethod.lower()
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

	#setattr( proc, jetALGO+'PF'+PUMethod+'pileupJetIdCalculator',
#			pileupJetIdCalculator.clone(
#				jets = cms.InputTag(jetalgo+'PFJets'+PUMethod),
#				rho = cms.InputTag("fixedGridRhoFastjetAll"),
#				vertexes = cms.InputTag('offlinePrimaryVertices'),
#				applyJec = cms.bool(True),
#				inputIsCorrected = cms.bool(False)
#				))

#	setattr( proc, jetALGO+'PF'+PUMethod+'pileupJetIdEvaluator',
#			pileupJetIdEvaluator.clone(
#				jetids = cms.InputTag(jetALGO+'PF'+PUMethod+'pileupJetIdCalculator'),
#				jets = cms.InputTag(jetalgo+'PFJets'+PUMethod),
#				rho = cms.InputTag("fixedGridRhoFastjetAll"),
#				vertexes = cms.InputTag('offlinePrimaryVertices'),
#                               applyJec = cms.bool(True),
#                                inputIsCorrected = cms.bool(False)
#
#				)
#			)

	#getattr( proc, 'patJets'+jetALGO+'PF'+PUMethod).userData.userFloats.src += [jetALGO+'PF'+PUMethod+'pileupJetIdEvaluator:fullDiscriminant']
	#getattr( proc, 'patJets'+jetALGO+'PF'+PUMethod).userData.userInts.src += [jetALGO+'PF'+PUMethod+'pileupJetIdEvaluator:cutbasedId',jetALGO+'PF'+PUMethod+'pileupJetIdEvaluator:fullId']

	QGjetsLabel='chs'

	setattr( proc, 'QGTagger'+jetALGO+'PF'+PUMethod,
		 QGTagger.clone(
			srcJets = cms.InputTag(jetalgo+'PFJets'+PUMethod),    # Could be reco::PFJetCollection or pat::JetCollection (both AOD and miniAOD)
			jetsLabel = cms.string('QGL_AK4PF'+QGjetsLabel)        # Other options (might need to add an ESSource for it): see https://twiki.cern.ch/twiki/bi
			)
		 )
	elemToKeep += [ 'keep *_QGTagger'+jetALGO+'PF'+PUMethod+'_*_*' ]
	getattr( proc, 'patJets'+jetALGO+'PF'+PUMethod).userData.userFloats.src += ['QGTagger'+jetALGO+'PF'+PUMethod+':qgLikelihood']
	jetSeq += getattr(proc, 'QGTagger'+jetALGO+'PF'+PUMethod )

	toolsUsed.append( 'QGTagger'+jetALGO+'PF'+PUMethod )

	getattr(proc,'patJetPartons').particles = cms.InputTag( 'genParticles' ) #'prunedGenParticles')
	setattr(proc, 'selectedPatJets'+jetALGO+'PF'+PUMethod, selectedPatJets.clone( src = 'patJets'+jetALGO+'PF'+PUMethod ) )
	setattr(proc, jetSequence, jetSeq)


#jetToolbox( process, 'ak4', 'ak4JetSubs')
jetToolbox( process, 'ak4', 'ak4JetSubs','CHS')

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#! Services
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.load('CommonTools.UtilAlgos.TFileService_cfi')
process.TFileService.fileName=cms.string('MC_ProcessedTreeProducer_without_L2L3.root')

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
process.patJets.getJetMCFlavour = True

genJetCollection = 'ak4GenJetsNoNu'
genParticleCollection = 'genParticles'

from PhysicsTools.JetMCAlgos.HadronAndPartonSelector_cfi import selectedHadronsAndPartons
process.selectedHadronsAndPartons = selectedHadronsAndPartons.clone(
    particles = genParticleCollection
)

from PhysicsTools.JetMCAlgos.AK4PFJetsMCFlavourInfos_cfi import ak4JetFlavourInfos

process.genJetFlavourInfos = ak4JetFlavourInfos.clone(
    jets = genJetCollection,
)

process.ak4 =  cms.EDAnalyzer('ProcessedTreeProducerBTag',
	## jet collections ###########################
	#pfjets          = cms.InputTag('selectedPatJetsAK4PF'),
	pfjetschs       = cms.InputTag('selectedPatJetsAK4PFCHS'),
	#pfpujetid       = cms.string('AK4PFpileupJetIdEvaluator:fullDiscriminant'),
	pfchsjetpuid    = cms.string('AK4PFCHSpileupJetIdEvaluator:fullDiscriminant'),
	## MET collection ####
	pfmet           = cms.InputTag('patMETs'),
	genjets         = cms.untracked.InputTag('ak4GenJetsNoNu'),
	## database entry for the uncertainties ######
	#PFPayloadName   = cms.string('AK4PF'),
	PFPayloadNameCHS= cms.string('AK4PFchs'),
	#jecUncSrc       = cms.string(''),
	jecUncSrcCHS    = cms.string(''),
        jecUncSrcNames  = cms.vstring(''),
        AK4             = cms.untracked.bool(True),
	## set the conditions for good Vtx counting ##
	offlineVertices = cms.InputTag('offlinePrimaryVertices'),
	goodVtxNdof     = cms.double(4),
	goodVtxZ        = cms.double(24),
	## rho #######################################
	srcCaloRho      = cms.InputTag('fixedGridRhoFastjetAllCalo'),
	srcPFRho        = cms.InputTag('fixedGridRhoFastjetAll'),
        srcPULabel      = cms.untracked.InputTag('addPileupInfo'),
	## preselection cuts #########################
	maxY            = cms.double(5.0),
	minPFPt         = cms.double(20),
	minNPFJets      = cms.int32(1),
	minGenPt        = cms.untracked.double(20),
	minJJMass       = cms.double(-1),
        isMCarlo        = cms.untracked.bool(True),
        useGenInfo      = cms.untracked.bool(True),
	## trigger ###################################
	printTriggerMenu = cms.untracked.bool(False),
	processName     = cms.string('HLT'),
	triggerName     = cms.vstring(''),
	triggerResults  = cms.InputTag("TriggerResults","","HLT"),
	triggerEvent    = cms.InputTag("hltTriggerSummaryAOD","","HLT"),
	## jec services ##############################
	HBHENoiseFilterResultLabel = cms.InputTag("HBHENoiseFilterResultProducer", "HBHENoiseFilterResult"),
	HBHENoiseFilterResultNoMinZLabel = cms.InputTag("HBHENoiseFilterResultProducerNoMinZ", "HBHENoiseFilterResult"),
        EventInfo       = cms.InputTag("generator"),
	GenParticles    = cms.InputTag("genParticles"),
	beamSpot        = cms.InputTag('offlineBeamSpot'),
	jetFlavourInfos = cms.InputTag("genJetFlavourInfos")

	#pfjecService    = cms.string('ak7PFL1FastL2L3Residual'),
)

#jetToolbox( process, 'ak8', 'ak8JetSubs','CHS')
#jetToolbox( process, 'ak8', 'ak8JetSubs')

#process.ak8 = process.ak4.clone(
	#pfjets          = cms.InputTag('selectedPatJetsAK8PF'),
#	pfjetschs       = cms.InputTag('selectedPatJetsAK8PFCHS'),
#	pfpujetid       = cms.string('AK8PFpileupJetIdEvaluator:fullDiscriminant'),
#	pfchsjetpuid    = cms.string('AK8PFCHSpileupJetIdEvaluator:fullDiscriminant'),
#        #PFPayloadName   = cms.string('AK8PF'),
#        PFPayloadNameCHS= cms.string('AK8PFchs'),
#        AK4             = cms.untracked.bool(False),
#)

jetToolbox( process, 'ak7', 'ak7JetSubs','CHS')
#jetToolbox( process, 'ak7', 'ak7JetSubs')

process.ak7 = process.ak4.clone(
	pfjets          = cms.InputTag('selectedPatJetsAK7PF'),
    pfjetschs       = cms.InputTag('selectedPatJetsAK7PFCHS'),
    pfpujetid       = cms.string('AK7PFpileupJetIdEvaluator:fullDiscriminant'),
    pfchsjetpuid    = cms.string('AK7PFCHSpileupJetIdEvaluator:fullDiscriminant'),
    PFPayloadNameCHS= cms.string('AK7PFchs'),
	AK4             = cms.untracked.bool(False),
    genjets         = cms.untracked.InputTag('ak7GenJetsNoNu')
)

#jetToolbox( process, 'ak5', 'ak5JetSubs','CHS')
#jetToolbox( process, 'ak5', 'ak5JetSubs')

#process.ak5 = process.ak4.clone(
#	pfjets          = cms.InputTag('selectedPatJetsAK5PF'),
#	pfjetschs       = cms.InputTag('selectedPatJetsAK5PFCHS'),
#	pfpujetid       = cms.string('AK5PFpileupJetIdEvaluator:fullDiscriminant'),
#	pfchsjetpuid    = cms.string('AK5PFCHSpileupJetIdEvaluator:fullDiscriminant'),
        #PFPayloadName   = cms.string('AK5PF'),
#        PFPayloadNameCHS= cms.string('AK5PFchs'),
#	AK4             = cms.untracked.bool(False),
#)

process.goodVertices = cms.EDFilter("VertexSelector",
   filter = cms.bool(False),
   src = cms.InputTag("offlinePrimaryVertices"),
   cut = cms.string("!isFake && ndof >= 4 && abs(z) <= 24 && position.rho <= 2"),
)
process.load('RecoMET.METFilters.trackingFailureFilter_cfi')

##Type1 patMET Producer
process.load('PhysicsTools.PatAlgos.recoLayer0.metCorrections_cff')
process.load('PhysicsTools.PatAlgos.producersLayer1.metProducer_cfi')
process.patMETs.addGenMET = cms.bool(False)
#Try scheduled processs
process.path = cms.Path(process.ak4*process.ak7)

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#! Output and Log
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )
process.options.allowUnscheduled = cms.untracked.bool(True)
