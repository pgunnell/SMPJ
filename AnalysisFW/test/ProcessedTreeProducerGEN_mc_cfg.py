# -*- coding: utf-8 -*-

import FWCore.ParameterSet.Config as cms

process = cms.Process("Ntuplizer")
process.load("Configuration.EventContent.EventContent_cff")

inFiles = cms.untracked.vstring(
'file:///mnt/storage/gflouris/MC_13TeV/Pythia8_CUETM1_Flat/QCD_Pt_15to5000_CUETM1_Flat_13TeV_pythia8_cffGEN_1.root' #Pythia
)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.source = cms.Source("PoolSource", fileNames = inFiles )

process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 10000
process.load('CommonTools.UtilAlgos.TFileService_cfi')
process.TFileService.fileName=cms.string('MC_ProcessedTreeProducerGEN.root')

process.ak5 =  cms.EDAnalyzer('ProcessedTreeProducerGEN',
	genjets         = cms.untracked.InputTag('ak5GenJets'),
	## preselection cuts #########################
	maxY            = cms.double(5.0),
	minGenPt        = cms.untracked.double(20),
)


process.ak7 = process.ak5.clone(
	genjets         = cms.untracked.InputTag('ak7GenJets'),
)

process.p = cms.Path( process.ak5*process.ak7)

