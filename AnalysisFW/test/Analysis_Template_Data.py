# -*- coding: utf-8 -*-
import FWCore.ParameterSet.Config as cms

process = cms.Process("myprocess")
process.TFileService=cms.Service("TFileService",fileName=cms.string('TriggerTagAndProbe-ak4.root'))

##-------------------- Define the source  ----------------------------
process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(1)
        )
process.source = cms.Source("EmptySource")

##-------------------- User analyzer  --------------------------------
process.efficiency  = cms.EDAnalyzer('Analysis_Template_MC',
                                     filename        = cms.vstring('','','','','root://eoscms//eos/cms/store/group/phys_smp/Multijet/13TeV/Data/Ntuples-Data-MagneticField-JetHt-JsonFile-25-Run2015C_v2.root','root://eoscms//eos/cms/store/group/phys_smp/Multijet/13TeV/Data/Ntuples-Data-MagneticField-JetHt-JsonFile-25-Run2015B_v2.root'),
                                     #filename        = cms.vstring('','','','','root://eoscms//eos/cms/store/group/phys_smp/Multijet/13TeV/Data/Ntuples-Data-MagneticField-JetHt-JsonFile-25-Run2015B_v2.root'),

     treename        = cms.string('ProcessedTree'),
     dirname         = cms.string('ak4'),
     
     minPt           = cms.double(20.0),
     ymax            = cms.double(4.7),
     JetID           = cms.int32(2),

     printOk         = cms.int32(1),
     pseudoglobaltag = cms.string('Summer15_50nsV5'), #Make sure that the name is the same as in data direcotry                                              
     jettype         = cms.string('AK4PFchs'),     #AKXPF or AKXPFchs  -  this is only for JECs  
                                     
          
     hcalNoiseFilter = cms.int32(0),
     isMCarlo        = cms.untracked.bool(False),
     PUReweighting        = cms.untracked.bool(False),
     LowPileUp  = cms.untracked.bool(False),                                   
     MCSlice = cms.int32(0),                                     
     jecUncSrc  = cms.string('../data/Summer15_50nsV5/Summer15_50nsV5_DATA_Uncertainty_AK4PFchs.txt'),#test.txt'),#Summer15_50nsV5_DATA_Uncertainty_AK8PFchs.txt'),
                                     #jecUncSrc  = cms.string('../data/Summer15_50nsV5/test.txt'),#Summer15_50nsV5_DATA_UncertaintySources_AK8PFchs.txt'),
     jecUncSrcNames  = cms.vstring(''),#'AbsoluteStat,AbsoluteFlavMap', 'AbsoluteMPFBias', 'Fragmentation', 'SinglePionECAL',
                                   #'SinglePionHCAL', 'FlavorQCD', 'TimeEta', 'TimePt', 'RelativeJEREC1', 'RelativeJEREC2', 'RelativeJERHF',
				   #'RelativePtBB', 'RelativePtEC1', 'RelativePtEC2', 'RelativePtHF', 'RelativeFSR', 'RelativeStatFSR',
				   #'RelativeStatEC2', 'RelativeStatHF', 'PileUpDataMC', 'PileUpPtRef', 'PileUpPtBB', 'PileUpPtEC1', 
				   #'PileUpPtEC2', 'PileUpPtHF', 'PileUpMuZero', 'PileUpEnvelope', 
				   #'SubTotalPileUp', 'SubTotalRelative', 'SubTotalPt', 'SubTotalScale', 'SubTotalAbsolute', 'SubTotalMC', 'Total','TotalNoFlavor',
				   #'TotalNoTime','TotalNoFlavorNoTime',
				   #'FlavorZJet','FlavorPhotonJet','FlavorPureGluon','FlavorPureQuark','FlavorPureCharm',
				   #'FlavorPureBottom','TimeRunA','TimeRunB','TimeRunC','TimeRunD','CorrelationGroupMPFInSitu',
				   #'CorrelationGroupIntercalibration','CorrelationGroupbJES','CorrelationGroupFlavor',
				   #'CorrelationGroupUncorrelated'),
		  )

process.p = cms.Path(process.efficiency)

