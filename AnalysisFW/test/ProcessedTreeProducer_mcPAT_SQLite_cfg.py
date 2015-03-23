## import skeleton process
from PhysicsTools.PatAlgos.patTemplate_cfg import *

runOnMC = True

process.source.fileNames = [
 "root://eoscms//eos/cms/store/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/TestGED_START53_V6-v1/0000/CC17C457-89AC-E111-A4EC-002481E15070.root"
 ]


if runOnMC:
        process.GlobalTag.globaltag = cms.string('START53_V27::All')
else:
        process.GlobalTag.globaltag = cms.string('FT53_V21A_AN6::All')


### ----------------------------- loading for the sqlite file ------------------------- ####
process.load("CondCore.DBCommon.CondDBCommon_cfi")
from CondCore.DBCommon.CondDBSetup_cfi import *
process.jec = cms.ESSource("PoolDBESSource",
      DBParameters = cms.PSet(
        messageLevel = cms.untracked.int32(0)
        ),
      timetype = cms.string('runnumber'),
      toGet = cms.VPSet(
      cms.PSet(
            record = cms.string('JetCorrectionsRecord'),
            tag    = cms.string('JetCorrectorParametersCollection_Winter14_V5_MC_AK5PFchs'),
            # tag    = cms.string('JetCorrectorParametersCollection_Summer12_V7_MC_AK5PF'),
            label  = cms.untracked.string('AK5PFchs')
            ),
      cms.PSet(
            record = cms.string('JetCorrectionsRecord'),
            tag    = cms.string('JetCorrectorParametersCollection_Winter14_V5_MC_AK7PFchs'),
            # tag    = cms.string('JetCorrectorParametersCollection_Summer12_V7_MC_AK7PF'),
            label  = cms.untracked.string('AK7PFchs')
            ),
      cms.PSet(
            record = cms.string('JetCorrectionsRecord'),
            tag    = cms.string('JetCorrectorParametersCollection_Winter14_V5_MC_AK5PF'),
            # tag    = cms.string('JetCorrectorParametersCollection_Summer12_V7_MC_AK5PF'),
            label  = cms.untracked.string('AK5PF')
            ),
      cms.PSet(
            record = cms.string('JetCorrectionsRecord'),
            tag    = cms.string('JetCorrectorParametersCollection_Winter14_V5_MC_AK7PF'),
            # tag    = cms.string('JetCorrectorParametersCollection_Summer12_V7_MC_AK7PF'),
            label  = cms.untracked.string('AK7PF')
            ),
      cms.PSet(
            record = cms.string('JetCorrectionsRecord'),
            tag    = cms.string('JetCorrectorParametersCollection_Winter14_V5_MC_AK5Calo'),
            # tag    = cms.string('JetCorrectorParametersCollection_Summer12_V7_MC_AK5Calo'),
            label  = cms.untracked.string('AK5Calo')
            ),
      cms.PSet(
            record = cms.string('JetCorrectionsRecord'),
            tag    = cms.string('JetCorrectorParametersCollection_Winter14_V5_MC_AK7Calo'),
            # tag    = cms.string('JetCorrectorParametersCollection_Summer12_V7_MC_AK7Calo'),
            label  = cms.untracked.string('AK7Calo')
            ),
      #..................................................
      ## here you add as many jet types as you need
      ## note that the tag name is specific for the particular sqlite file 
      ), 
      connect = cms.string('sqlite:Winter14_V5_MC.db')
     # uncomment above tag lines and this comment to use MC JEC
     # connect = cms.string('sqlite:Summer12_V7_MC.db')
)
## add an es_prefer statement to resolve a possible conflict from simultaneous connection to a global tag
process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')


##--------------------- ******************************************************************** --------------------------- ##


# load the PAT config
process.load("RecoTauTag.Configuration.RecoPFTauTag_cff") # re-run tau discriminators (new version)
process.load("PhysicsTools.PatAlgos.patSequences_cff")
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('RecoJets.Configuration.RecoPFJets_cff')
process.load('RecoJets.Configuration.RecoJets_cff')
process.load('RecoJets.JetProducers.TrackJetParameters_cfi')
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')


# Configure PAT to use PF2PAT instead of AOD sources
# this function will modify the PAT sequences.
from PhysicsTools.PatAlgos.tools.pfTools import *
from PhysicsTools.PatAlgos.tools.coreTools import *
from PhysicsTools.PatAlgos.tools.metTools import *
from PhysicsTools.PatAlgos.tools.jetTools import *
from PhysicsTools.PatAlgos.tools.coreTools import *
from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector


##--------- good primary vertices ---------------
#process.goodOfflinePrimaryVertices = cms.EDFilter(
#    "PrimaryVertexObjectFilter",
#    filterParams = pvSelector.clone( minNdof = cms.double(4.0), maxZ = cms.double(24.0) ),
#    src=cms.InputTag('offlinePrimaryVertices')
#    )

process.goodOfflinePrimaryVertices = cms.EDFilter(
    "VertexSelector",
    filter = cms.bool(False),
    src = cms.InputTag("offlinePrimaryVertices"),
    cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.rho < 2")
    )


process.ak5PFJets.doAreaFastjet = True
process.ak7PFJets.doAreaFastjet = True
process.kt6PFJets.doRhoFastjet = True



# Configure PAT to use PF2PAT instead of AOD sources
# this function will modify the PAT sequences.
from PhysicsTools.PatAlgos.tools.pfTools import *

# An empty postfix means that only PF2PAT is run,
# otherwise both standard PAT and PF2PAT are run. In the latter case PF2PAT
# collections have standard names + postfix (e.g. patElectronPFlow)
postfix = "CHS"
jetAlgo = "AK5"
#usePF2PAT(process,runPF2PAT=True, jetAlgo=jetAlgo, runOnMC=runOnMC, postfix=postfix)
usePF2PAT(process,runPF2PAT=True, jetAlgo=jetAlgo, runOnMC=runOnMC, postfix=postfix,
         jetCorrections=('AK5PFchs', ['L1FastJet','L2Relative','L3Absolute']),
         pvCollection=cms.InputTag('goodOfflinePrimaryVertices'),
         typeIMetCorrections=True
         )
process.pfPileUpCHS.checkClosestZVertex = False
process.pfPileUpCHS.Enable = True
process.pfPileUpCHS.Vertices = cms.InputTag('goodOfflinePrimaryVertices')
process.pfJetsCHS.doAreaFastjet = True
process.pfJetsCHS.doRhoFastjet = False



process.kt6PFJetsCHS = process.kt6PFJets.clone(
    src = cms.InputTag('pfNoElectron'+postfix),
    doAreaFastjet = cms.bool(True),
    doRhoFastjet = cms.bool(True)
    )

process.kt6PFJetsISO = process.kt6PFJets.clone(
    Rho_EtaMax = cms.double(2.4)
    )

getattr(process,"patPF2PATSequence"+postfix).replace(
    getattr(process,"pfNoElectron"+postfix),
    getattr(process,"pfNoElectron"+postfix)*process.kt6PFJetsCHS)


#----- recommendation from JES: use the standard rho for CHS ------
#process.patJetCorrFactorsCHS.rho = cms.InputTag("kt6PFJetsCHS", "rho")


#------------- Running a second PF2PAT for Ak7chs jets ----------------#
postfix2 = "CHS7"
jetAlgo2 = "AK7"
usePF2PAT(process,runPF2PAT=True, jetAlgo=jetAlgo2, runOnMC=runOnMC, postfix=postfix2,
         jetCorrections=('AK7PFchs', ['L1FastJet','L2Relative','L3Absolute']),
         pvCollection=cms.InputTag('goodOfflinePrimaryVertices'),
         typeIMetCorrections=True
         )

# --- modifying the pfpileup+postfix2 -------------------- #
process.pfPileUpCHS7.checkClosestZVertex = False
process.pfPileUpCHS7.Enable = True
process.pfPileUpCHS7.Vertices = cms.InputTag('goodOfflinePrimaryVertices')
process.pfJetsCHS7.doAreaFastjet = True
process.pfJetsCHS7.doRhoFastjet = False

process.ak5PFJets.doAreaFastjet = True
process.ak7PFJets.doAreaFastjet = True
process.kt6PFJets.doRhoFastjet = True

#----- recommendation from JES: use the standard rho for CHS ------
#process.patJetCorrFactorsCHS7.rho = cms.InputTag("kt6PFJetsCHS", "rho")

# to use tau-cleaned jet collection uncomment the following:
getattr(process,"pfNoTau"+postfix).enable = True
getattr(process,"pfNoTau"+postfix2).enable = True

#------ removing the MC matching and let it run -------#
if not runOnMC:
        removeMCMatchingPF2PAT( process, '' )
        runOnData(process)
  

# ------------- Adding Ak5 and Ak7 jet collection to process ------------- #
addPfMET(process, 'PF')
# ------- Adding non CHS jets to process ------//
addJetCollection(process,cms.InputTag('ak5PFJets'),
                 'AK5', 'PFCorr',
                 doJTA        = True,
                 doBTagging   = False,
                 jetCorrLabel = ('AK5PF', cms.vstring(['L1FastJet','L2Relative','L3Absolute'])),
                 doType1MET   = True,
                 doL1Cleaning = True,
                 doL1Counters = False,
                 genJetCollection=cms.InputTag("ak5GenJets"),
                 doJetID      = True,
                 jetIdLabel   = "ak5"
                 )

addJetCollection(process,cms.InputTag('ak7PFJets'),
                 'AK7', 'PFCorr',
                 doJTA        = True,
                 doBTagging   = False,
                 jetCorrLabel = ('AK7PF', cms.vstring(['L1FastJet','L2Relative','L3Absolute'])),
                 doType1MET   = True,
                 doL1Cleaning = True,
                 doL1Counters = False,
                 genJetCollection=cms.InputTag("ak7GenJets"),
                 doJetID      = True,
                 jetIdLabel   = "ak7"
                 )

#process.patJetCorrFactorsAK5PFcorr.primaryVertices = cms.InputTag('goodOfflinePrimaryVertices')
#process.patJetCorrFactorsAK5PFcorr.rho = cms.InputTag("kt6PFJetsCHS", "rho")
#for icorr in [process.patJetCorrFactorsAK5PFcorr,
#              process.patJetCorrFactorsAK7PFcorr
#              ] :
#    icorr.rho = cms.InputTag("kt6PFJets", "rho")
#    icorr.primaryVertices = cms.InputTag('goodOfflinePrimaryVertices')





################### declaring the EDAnalyzer ##############################3

process.ak7 = cms.EDAnalyzer('ProcessedTreeProducer',
    ## jet collections ###########################
    pfjets          = cms.InputTag('selectedPatJetsAK7PFCorr'),
    pfjetschs       = cms.InputTag('selectedPatJetsCHS7'),
    genjets         = cms.untracked.InputTag('ak7GenJets'),
    ## MET collection ####
    pfmet           = cms.InputTag('pfMETCHS7'),
    ## database entry for the uncertainties ######
    PFPayloadName   = cms.string(''),
    PFPayloadNameCHS= cms.string(''),
    CaloPayloadName = cms.string(''),
    jecUncSrc       = cms.string(''),
    jecUncSrcCHS    = cms.string(''), 
    jecUncSrcNames  = cms.vstring(''),
    srcPU           = cms.untracked.InputTag('addPileupInfo'),
    ## set the conditions for good Vtx counting ##
    offlineVertices = cms.InputTag('goodOfflinePrimaryVertices'),
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
    minNPFJets      = cms.int32(1),
    minJJMass       = cms.double(-1),
    isMCarlo        = cms.untracked.bool(runOnMC),
    useGenInfo      = cms.untracked.bool(True),
    ## trigger ###################################
    printTriggerMenu = cms.untracked.bool(True),
    processName     = cms.string('HLT'),
    triggerName     = cms.vstring('HLT_PFJet40_v3','HLT_PFJet40_v4','HLT_PFJet40_v5','HLT_PFJet40_v6','HLT_PFJet40_v7','HLT_PFJet40_v8','HLT_PFJet40_v9'
                            #      'HLT_IsoMu24_eta2p1_v11', 'HLT_IsoMu24_eta2p1_v12', 'HLT_IsoMu24_eta2p1_v13', 'HLT_IsoMu24_eta2p1_v14', 'HLT_IsoMu24_eta2p1_v15'
    ),
    triggerResults  = cms.InputTag("TriggerResults","","HLT"),
    triggerEvent    = cms.InputTag("hltTriggerSummaryAOD","","HLT"),
    ## jec services ##############################
    #pfjecService    = cms.string('ak7PFL1FastL2L3Residual'),
    #calojecService  = cms.string('ak7CaloL1FastL2L3Residual')
)

process.ak5 = process.ak7.clone(

    ## jet collections ###########################
    pfjets          = cms.InputTag('selectedPatJetsAK5PFCorr'),
    pfjetschs       = cms.InputTag('selectedPatJetsCHS'),
    genjets         = cms.untracked.InputTag('ak5GenJets'),
    ## MET collection ####
    pfmet           = cms.InputTag('pfMETCHS'),
    ## database entry for the uncertainties ######
    printTriggerMenu = False
)




# Let it run
process.p = cms.Path(
    process.goodOfflinePrimaryVertices*
    process.patDefaultSequence *
    getattr(process,"patPF2PATSequence"+postfix) *
    getattr(process,"patPF2PATSequence"+postfix2) * 
    process.ak5 *
    process.ak7
#    second PF2PAT
#    + getattr(process,"patPF2PATSequence"+postfix2)
)
#if not postfix == "":
#    process.p += process.recoTauClassicHPSSequence # re-run tau discriminators (new version)
#    process.p += process.patDefaultSequence

# Add PF2PAT output to the created file
#from PhysicsTools.PatAlgos.patEventContent_cff import patEventContentNoCleaning
#process.out.outputCommands = cms.untracked.vstring('drop *',
#                                                   'keep recoPFCandidates_particleFlow_*_*',
#                                                   *patEventContentNoCleaning )


# top projections in PF2PAT:
getattr(process,"pfNoPileUp"+postfix).enable = True
getattr(process,"pfNoMuon"+postfix).enable = True
getattr(process,"pfNoElectron"+postfix).enable = True
getattr(process,"pfNoTau"+postfix).enable = False
getattr(process,"pfNoJet"+postfix).enable = True

# verbose flags for the PF2PAT modules
getattr(process,"pfNoMuon"+postfix).verbose = False

## ------------------------------------------------------
#  In addition you usually want to change the following
#  parameters:
## ------------------------------------------------------
#
#   process.GlobalTag.globaltag =  ...    ##  (according to https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideFrontierConditions)
#                                         ##
#   process.source.fileNames =  ...       ##  (e.g. 'file:AOD.root')
#                                         ##
process.maxEvents.input = 1000
#
process.MessageLogger.cerr.FwkReport.reportEvery = 100                                         ##
#   process.out.outputCommands = [ ... ]  ##  (e.g. taken from PhysicsTools/PatAlgos/python/patEventContent_cff.py)
#
############# processed tree producer ##################
process.TFileService = cms.Service("TFileService",fileName = cms.string('ProcessedTree_mc.root'))
           
#process.out.fileName = 'patTuple_PATandPF2PAT.root'
#                                         ##
process.options.wantSummary = False   ##  (to suppress the long output at the end of the job)

del process.outpath
