# -*- coding: utf-8 -*-
import FWCore.ParameterSet.Config as cms

process = cms.Process("myprocess")
process.TFileService=cms.Service("TFileService",fileName=cms.string("Pythia8-MagneticField-Sliced-AK4-Summer15_50nsV5-HighStat-ReweightRunC.root"))

##-------------------- Define the source  ----------------------------
process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(1)
        )
process.source = cms.Source("EmptySource")


##-------------------- User analyzer  --------------------------------
process.efficiency  = cms.EDAnalyzer('Analysis_Template_MC',
                                                     filename        = cms.vstring(
#filename        = cms.vstring("root://eoscms//eos/cms/store/group/phys_smp/Multijet/13TeV/MC/Pythia8-CUETP8M1-Ntuples-PFJets.root"),
"root://eoscms//eos/cms/store/group/phys_smp/Multijet/13TeV/MC/Pythia8_CUETP8M1_PtHatSliced/Pythia8-CUETP8M1-Ntuples-5-10GeVPtHat.root","root://eoscms//eos/cms/store/group/phys_smp/Multijet/13TeV/MC/Pythia8_CUETP8M1_PtHatSliced/Pythia8-CUETP8M1-Ntuples-10-15GeVPtHat.root","root://eoscms//eos/cms/store/group/phys_smp/Multijet/13TeV/MC/Pythia8_CUETP8M1_PtHatSliced/Pythia8-CUETP8M1-Ntuples-15-30GeVPtHat.root","root://eoscms//eos/cms/store/group/phys_smp/Multijet/13TeV/MC/Pythia8_CUETP8M1_PtHatSliced/Pythia8-CUETP8M1-Ntuples-30-50GeVPtHat.root","root://eoscms//eos/cms/store/group/phys_smp/Multijet/13TeV/MC/Pythia8_CUETP8M1_PtHatSliced/Pythia8-CUETP8M1-Ntuples-50-80GeVPtHat.root","root://eoscms//eos/cms/store/group/phys_smp/Multijet/13TeV/MC/Pythia8_CUETP8M1_PtHatSliced/Pythia8-CUETP8M1-Ntuples-80-120GeVPtHat.root","root://eoscms//eos/cms/store/group/phys_smp/Multijet/13TeV/MC/Pythia8_CUETP8M1_PtHatSliced/Pythia8-CUETP8M1-Ntuples-120-170GeVPtHat.root","root://eoscms//eos/cms/store/group/phys_smp/Multijet/13TeV/MC/Pythia8_CUETP8M1_PtHatSliced/Pythia8-CUETP8M1-Ntuples-170-300GeVPtHat.root","root://eoscms//eos/cms/store/group/phys_smp/Multijet/13TeV/MC/Pythia8_CUETP8M1_PtHatSliced/Pythia8-CUETP8M1-Ntuples-300-470GeVPtHat.root","root://eoscms//eos/cms/store/group/phys_smp/Multijet/13TeV/MC/Pythia8_CUETP8M1_PtHatSliced/Pythia8-CUETP8M1-Ntuples-470-600GeVPtHat.root","root://eoscms//eos/cms/store/group/phys_smp/Multijet/13TeV/MC/Pythia8_CUETP8M1_PtHatSliced/Pythia8-CUETP8M1-Ntuples-600-800GeVPtHat.root","root://eoscms//eos/cms/store/group/phys_smp/Multijet/13TeV/MC/Pythia8_CUETP8M1_PtHatSliced/Pythia8-CUETP8M1-Ntuples-800-1000GeVPtHat.root","root://eoscms//eos/cms/store/group/phys_smp/Multijet/13TeV/MC/Pythia8_CUETP8M1_PtHatSliced/Pythia8-CUETP8M1-Ntuples-1000-1400GeVPtHat.root","root://eoscms//eos/cms/store/group/phys_smp/Multijet/13TeV/MC/Pythia8_CUETP8M1_PtHatSliced/Pythia8-CUETP8M1-Ntuples-1400-1800GeVPtHat.root","root://eoscms//eos/cms/store/group/phys_smp/Multijet/13TeV/MC/Pythia8_CUETP8M1_PtHatSliced/Pythia8-CUETP8M1-Ntuples-1800-2400GeVPtHat.root","root://eoscms//eos/cms/store/group/phys_smp/Multijet/13TeV/MC/Pythia8_CUETP8M1_PtHatSliced/Pythia8-CUETP8M1-Ntuples-2400-3200GeVPtHat.root","root://eoscms//eos/cms/store/group/phys_smp/Multijet/13TeV/MC/Pythia8_CUETP8M1_PtHatSliced/Pythia8-CUETP8M1-Ntuples-3200-InfGeVPtHat.root"),
                                     treename        = cms.string('ProcessedTree'),
    dirname         = cms.string('ak4'),

    #Variables for JECs on the fly
    pseudoglobaltag = cms.string('Summer15_50nsV5'), #Make sure that the name is the same as in data direcotry
    jettype         = cms.string('AK4PFchs'),     #AKXPF or AKXPFchs  -  this is only for JECs
    #
    minPt           = cms.double(50.0),
    ymax            = cms.double(4.7),
    JetID           = cms.int32(2),

    printOk         = cms.int32(18),

    #xsection        = cms.double(10360000.0),
    #mcevents        = cms.double(52000000),
    MCSlice         = cms.int32(1),                                 
    PUReweighting   = cms.untracked.bool(True),                                 
    LowPileUp  = cms.untracked.bool(False),                                 
    isMCarlo        = cms.untracked.bool(True),
                                     jecUncSrc  = cms.string('../data/Summer15_50nsV5/Summer15_50nsV5_MC_Uncertainty_AK4PFchs.txt'),    
    #jecUncSrc  = cms.string('../data/Summer15_50nsV5/test.txt'),    
    jecUncSrcNames  = cms.vstring(''),
     #jecUncSrcNames  = cms.vstring('Absolute','HighPtExtra','SinglePionECAL','SinglePionHCAL','FlavorQCD',
                               #'Time','RelativeJEREC1','RelativeJEREC2','RelativeJERHF','RelativePtBB',
                               #'RelativePtEC1','RelativePtEC2','RelativePtHF','RelativeFSR','RelativeStatEC2',
                               #'RelativeStatHF','PileUpDataMC','PileUpPtBB','PileUpPtEC','PileUpPtHF','PileUpBias',
                               #'TotalNoFlavor', 'FlavorZJet','FlavorPhotonJet','FlavorPureGluon','FlavorPureQuark',
                               #'FlavorPureCharm','FlavorPureBottom',
                               #'SubTotalPileUp','SubTotalRelative','SubTotalPt','SubTotalMC','Total'),

)

process.p = cms.Path(process.efficiency)
