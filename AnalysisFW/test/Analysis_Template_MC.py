# -*- coding: utf-8 -*-
import FWCore.ParameterSet.Config as cms

process = cms.Process("myprocess")
process.TFileService=cms.Service("TFileService",fileName=cms.string("out_Analysis_Template_MC.root"))

##-------------------- Define the source  ----------------------------
process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(1)
        )
process.source = cms.Source("EmptySource")


##-------------------- User analyzer  --------------------------------
process.efficiency  = cms.EDAnalyzer('Analysis_Template_MC',
    filename        = cms.string("file://./MC_ProcessedTreeProducer_2.root"),

    treename        = cms.string('ProcessedTree'),
    dirname         = cms.string('ak7'),

    #Variables for JECs on the fly
    pseudoglobaltag = cms.string('PHYS14_V4'), #Make sure that the name is the same as in data direcotry
    jettype         = cms.string('AK7PF'),     #AKXPF or AKXPFchs  -  this is only for JECs
    #

    minPt           = cms.double(100.0),
    ymax            = cms.double(2.5),
    JetID           = cms.int32(2),

    printOk         = cms.int32(0),

    #xsection        = cms.double(10360000.0),
    #mcevents        = cms.double(52000000),

    isMCarlo        = cms.untracked.bool(True),
    jecUncSrcNames  = cms.vstring('Absolute','HighPtExtra','SinglePionECAL','SinglePionHCAL','FlavorQCD',
                               'Time','RelativeJEREC1','RelativeJEREC2','RelativeJERHF','RelativePtBB',
                               'RelativePtEC1','RelativePtEC2','RelativePtHF','RelativeFSR','RelativeStatEC2',
                               'RelativeStatHF','PileUpDataMC','PileUpPtBB','PileUpPtEC','PileUpPtHF','PileUpBias',
                               'TotalNoFlavor', 'FlavorZJet','FlavorPhotonJet','FlavorPureGluon','FlavorPureQuark',
                               'FlavorPureCharm','FlavorPureBottom',
                               'SubTotalPileUp','SubTotalRelative','SubTotalPt','SubTotalMC','Total'),

)

process.p = cms.Path(process.efficiency)
