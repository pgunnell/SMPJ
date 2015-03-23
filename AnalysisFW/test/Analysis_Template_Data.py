# -*- coding: utf-8 -*-
import FWCore.ParameterSet.Config as cms

process = cms.Process("myprocess")
process.TFileService=cms.Service("TFileService",fileName=cms.string('My_azi_data_Run2012A_test.root'))

##-------------------- Define the source  ----------------------------
process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(1)
        )
process.source = cms.Source("EmptySource")

##-------------------- User analyzer  --------------------------------
process.efficiency  = cms.EDAnalyzer('Analysis_Template_Data',

     filename        = cms.string('/mnt/storage/Azimuthal/Data_19thNovember2014/19thNovember2014_Run2012A.root'),

     treename        = cms.string('ProcessedTree'),
     dirname         = cms.string('ak7'),
     
     minPt           = cms.double(100.0),
     ymax            = cms.double(2.5),
     JetID           = cms.int32(2),

     printOk         = cms.int32(1),

          
     hcalNoiseFilter = cms.int32(0),
     isMCarlo        = cms.untracked.bool(False),

     jecUncSrcNames  = cms.vstring('AbsoluteStat','AbsoluteScale','AbsoluteFlavMap', 'AbsoluteMPFBias', 'Fragmentation', 'SinglePionECAL',
				    'SinglePionHCAL', 'FlavorQCD', 'TimeEta', 'TimePt', 'RelativeJEREC1', 'RelativeJEREC2', 'RelativeJERHF',
				    'RelativePtBB', 'RelativePtEC1', 'RelativePtEC2', 'RelativePtHF', 'RelativeFSR', 'RelativeStatFSR',
				    'RelativeStatEC2', 'RelativeStatHF', 'PileUpDataMC', 'PileUpPtRef', 'PileUpPtBB', 'PileUpPtEC1', 
				    'PileUpPtEC2', 'PileUpPtHF', 'PileUpMuZero', 'PileUpEnvelope', 
				    'SubTotalPileUp', 'SubTotalRelative', 'SubTotalPt', 'SubTotalScale', 'SubTotalAbsolute', 'SubTotalMC', 'Total','TotalNoFlavor',
				    'TotalNoTime','TotalNoFlavorNoTime',
				    'FlavorZJet','FlavorPhotonJet','FlavorPureGluon','FlavorPureQuark','FlavorPureCharm',
				    'FlavorPureBottom','TimeRunA','TimeRunB','TimeRunC','TimeRunD','CorrelationGroupMPFInSitu',
				    'CorrelationGroupIntercalibration','CorrelationGroupbJES','CorrelationGroupFlavor',
				    'CorrelationGroupUncorrelated'),
		  )

process.p = cms.Path(process.efficiency)

