# -*- coding: utf-8 -*-
###############################################
####
####   Jet Substructure Toolbox (jetToolBox)
####   Python function for easy access of jet substructure tools implemented in CMS
####
####   Alejandro Gomez Espinosa (gomez@physics.rutgers.edu)
####   First version: 2015 - 01 - 28
####
###############################################
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


def jetToolbox( proc, jetType, jetSequence, outputFile, 
		PUMethod='CHS',                    #### Options: Puppi, CS, SK
		JETCorrPayload='None', JETCorrLevels = [ 'None' ],
		subJETCorrPayload='None', subJETCorrLevels = [ 'None' ],
		miniAOD=False,
		Cut = '', 
		addPruning=False, zCut=0.1, rCut=0.5, addPrunedSubjets=False,
		addSoftDrop=False, betaCut=0.0,  zCutSD=0.1, addSoftDropSubjets=False,
		addTrimming=False, rFiltTrim=0.2, ptFrac=0.03,
		addFiltering=False, rfilt=0.3, nfilt=3,
		addCMSTopTagger=False,
		addMassDrop=False,
		addHEPTopTagger=False,
		addNsub=False, maxTau=4, 
		addQJets=False 
		):
	
	###############################################################################
	#######  Verifying some inputs and defining variables
	###############################################################################
	print '|---- jetToolbox: Initialyzing...'
	supportedJetAlgos = { 'ak': 'AntiKt', 'ca' : 'CambridgeAachen', 'kt' : 'Kt' }
	recommendedJetAlgos = [ 'ak4', 'ak8', 'ca4', 'ca8', 'ca10' ]
	payloadList = [ 'AK1PFchs', 'AK2PFchs', 'AK3PFchs', 'AK4PFchs', 'AK5PFchs', 'AK6PFchs', 'AK7PFchs', 'AK8PFchs', 'AK9PFchs', 'AK10PFchs' ]
	JECLevels = [ 'L1Offset', 'L1FastJet', 'L1JPTOffset', 'L2Relative', 'L3Absolute', 'L5Falvour', 'L7Parton' ]
	jetAlgo = ''
	algorithm = ''
	size = ''
	for type, tmpAlgo in supportedJetAlgos.iteritems(): 
		if type in jetType.lower():
			jetAlgo = type
			algorithm = tmpAlgo
			size = jetType.replace( type, '' )

	jetSize = 0.
	if int(size) in range(0, 20): jetSize = int(size)/10.
	else: print '|---- jetToolBox: jetSize has not a valid value. Insert a number between 1 and 20 after algorithm, like: AK8'
	### Trick for uppercase/lowercase algo name
	jetALGO = jetAlgo.upper()+size
	jetalgo = jetAlgo.lower()+size
	if jetalgo not in recommendedJetAlgos: print '|---- jetToolBox: CMS recommends the following jet algoritms: '+' '.join(recommendedJetAlgos)+'. You are using', jetalgo,'.'

	if JETCorrPayload not in payloadList:
		if( int(size) > 10 ): 
			size = '10' 
			print '|---- jetToolbox: For jets bigger than 1.0, the jet corrections are AK10PFchs.'
		if not 'None' in JETCorrPayload: print '|---- jetToolBox: Payload given for Jet corrections ('+JETCorrPayload+') is not correct. Using a default AK'+size+'PFchs instead.'
		JETCorrPayload = 'AK'+size+'PFchs'
	else: print '|---- jetToolBox: Using '+JETCorrPayload+' payload for jet corrections.'

	if not set(JETCorrLevels).issubset(set(JECLevels)) :
		if not 'None' in JETCorrLevels: print '|---- jetToolbox: JEC levels given ( '+' '.join(JETCorrLevels)+' ) are incorrect. Using the default levels: L1FastJet, L2Relative, L3Absolute.'
		JETCorrLevels = ['L1FastJet', 'L2Relative', 'L3Absolute']
	
	print '|---- jetToolBox: Applying these jet corrections: ( '+JETCorrPayload+', '+' '.join(JETCorrLevels)+' )'
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


	print '|---- jetToolBox: JETTOOLBOX RUNNING ON AOD FOR '+jetALGO+' JETS USING '+PUMethod

	genParticlesLabel = 'genParticles'
	pvLabel = 'offlinePrimaryVertices'
	tvLabel = 'generalTracks'
	pfCand = 'particleFlow'
	svLabel = 'inclusiveSecondaryVertices'

	proc.load('RecoJets.Configuration.GenJetParticles_cff')
	proc.load('CommonTools.ParticleFlow.pfNoPileUpJME_cff')
	setattr( proc, jetalgo+'GenJetsNoNu', ak4GenJets.clone( src = 'genParticlesForJetsNoNu', rParam = jetSize, jetAlgorithm = algorithm ) ) 
	jetSeq += getattr(proc, jetalgo+'GenJetsNoNu' )
	setattr( proc, jetalgo+'PFJetsCHS', ak4PFJets.clone( rParam = jetSize, jetAlgorithm = algorithm ) ) 
	jetSeq += getattr(proc, jetalgo+'PFJetsCHS' )

	

	## b-tag discriminators
	#bTagDiscriminators = [
			#'pfTrackCountingHighEffBJetTags',
			#'pfTrackCountingHighPurBJetTags',
			#'pfJetProbabilityBJetTags',
			#'pfJetBProbabilityBJetTags',
			#'pfSimpleSecondaryVertexHighEffBJetTags',
			#'pfSimpleSecondaryVertexHighPurBJetTags',
			#'pfCombinedSecondaryVertexBJetTags',
			#'pfCombinedInclusiveSecondaryVertexV2BJetTags'
	    #]

	bTagDiscriminators = [
			'positiveOnlyJetProbabilityJetTags',
	    ]

	####  Creating PATjets
	setattr( proc, jetalgo+'PFJetsCHS', 
			ak4PFJetsCHS.clone( 
				doAreaFastjet = True, 
				rParam = jetSize, 
				jetAlgorithm = algorithm ) ) 
	if miniAOD: getattr( proc, jetalgo+'PFJetsCHS').src = 'chs'
	jetSeq += getattr(proc, jetalgo+'PFJetsCHS' )
	#elemToKeep += [ 'keep *_'+jetalgo+'PFJetsCHS_*_*' ]

	#if miniAOD: setattr( proc, jetalgo+'PFJets'+PUMethod+'Constituents', cms.EDFilter("MiniAODJetConstituentSelector", src = cms.InputTag( jetalgo+'PFJets'+PUMethod ), cut = cms.string( 'pt > 100.0 ' ) ))
	#else: setattr( proc, jetalgo+'PFJets'+PUMethod+'Constituents', cms.EDFilter("PFJetConstituentSelector", src = cms.InputTag( jetalgo+'PFJets'+PUMethod ), cut = cms.string( Cut ) ))
	#jetSeq += getattr(proc, jetalgo+'PFJets'+PUMethod+'Constituents' )

	addJetCollection(
			proc,
			labelName = jetALGO+'PF'+PUMethod,
			jetSource = cms.InputTag( jetalgo+'PFJets'+PUMethod),
			algo = jetalgo,
			rParam = jetSize,
			jetCorrections = JEC, #( 'AK'+size+'PFchs', cms.vstring( ['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
			pfCandidates = cms.InputTag( pfCand ),  #'packedPFCandidates'),
			svSource = cms.InputTag( svLabel ),   #'slimmedSecondaryVertices'),
			genJetCollection = cms.InputTag( jetalgo+'GenJetsNoNu'),
			pvSource = cms.InputTag( pvLabel ), #'offlineSlimmedPrimaryVertices'),
			btagDiscriminators = bTagDiscriminators,
			outputModules = ['outputFile']
			) 

	getattr( proc, 'patJetCorrFactors'+jetALGO+'PF'+PUMethod ).primaryVertices = pvLabel  #'offlineSlimmedPrimaryVertices' 
	getattr( proc, 'jetTracksAssociatorAtVertex'+jetALGO+'PF'+PUMethod ).tracks = tvLabel  # 'unpackedTracksAndVertices'

	getattr(proc,'patJetPartons').particles = cms.InputTag( genParticlesLabel ) #'prunedGenParticles')
	getattr(proc,'patJetPartonMatch'+jetALGO+'PF'+PUMethod).matched = cms.InputTag( genParticlesLabel ) #'prunedGenParticles')
	if hasattr(proc,'pfInclusiveSecondaryVertexFinderTagInfos'+jetALGO+'PF'+PUMethod):
		    getattr(proc,'pfInclusiveSecondaryVertexFinderTagInfos'+jetALGO+'PF'+PUMethod).extSVCollection = cms.InputTag( svLabel ) #'slimmedSecondaryVertices')

	setattr( proc, 'selectedPatJets'+jetALGO+'PF'+PUMethod, selectedPatJets.clone( src = 'patJets'+jetALGO+'PF'+PUMethod, cut = Cut ) )
	elemToKeep += [ 'keep *_selectedPatJets'+jetALGO+'PF'+PUMethod+'_*_*' ]
	elemToKeep += [ 'drop *_selectedPatJets'+jetALGO+'PF'+PUMethod+'_calo*_*' ]
	elemToKeep += [ 'drop *_selectedPatJets'+jetALGO+'PF'+PUMethod+'_tagInfos_*' ]

	print '|---- jetToolBox: Running '+', '.join(toolsUsed)+'.'

	### "return"
	setattr(proc, jetSequence, jetSeq)
	if	hasattr(proc, outputFile): getattr(proc, outputFile).outputCommands += elemToKeep
	else: setattr( proc, outputFile, 
			cms.OutputModule('PoolOutputModule', 
				fileName = cms.untracked.string('jettoolbox.root'), 
				outputCommands = cms.untracked.vstring( elemToKeep ) ) )


