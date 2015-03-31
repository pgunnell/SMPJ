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
filename        = cms.string("file://./MC_ProcessedTreeProducer.root"),

treename        = cms.string('ProcessedTree'),
dirname         = cms.string('ak7'),
   
minPt           = cms.double(100.0),
ymax            = cms.double(2.5),
JetID           = cms.int32(2),

printOk         = cms.int32(1),
	
#xsection        = cms.double(10360000.0),
#mcevents        = cms.double(52000000),
	          
)	
process.p = cms.Path(process.efficiency)
