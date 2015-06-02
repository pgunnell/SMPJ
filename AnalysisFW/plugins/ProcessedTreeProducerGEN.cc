#include <iostream>
#include <sstream>
#include <istream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <functional>
#include "TTree.h"
#include <vector>
#include <cassert>
#include <TLorentzVector.h>

#include "SMPJ/AnalysisFW/plugins/ProcessedTreeProducerGEN.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/JetExtendedAssociation.h"
#include "DataFormats/JetReco/interface/JetID.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

ProcessedTreeProducerGEN::ProcessedTreeProducerGEN(edm::ParameterSet const& cfg)
{
  mMaxY              = cfg.getParameter<double>                    ("maxY");
  mGenJetsName       = cfg.getUntrackedParameter<edm::InputTag>    ("genjets",edm::InputTag(""));
  mMinGenPt          = cfg.getUntrackedParameter<double>           ("minGenPt",30);
}
//////////////////////////////////////////////////////////////////////////////////////////
void ProcessedTreeProducerGEN::beginJob()
{
  mTree = fs->make<TTree>("ProcessedTree","ProcessedTree");
  mEvent = new QCDEvent();
  mTree->Branch("events","QCDEvent",&mEvent);
}
//////////////////////////////////////////////////////////////////////////////////////////
void ProcessedTreeProducerGEN::endJob()
{
}
//////////////////////////////////////////////////////////////////////////////////////////
void ProcessedTreeProducerGEN::beginRun(edm::Run const & iRun, edm::EventSetup const& iSetup)
{

}
//////////////////////////////////////////////////////////////////////////////////////////
void ProcessedTreeProducerGEN::analyze(edm::Event const& event, edm::EventSetup const& iSetup)
{
  vector<LorentzVector> mGenJets;
  QCDEventHdr mEvtHdr;
  //-------------- Basic Event Info ------------------------------
  mEvtHdr.setRun(event.id().run());
  mEvtHdr.setEvt(event.id().event());
  mEvtHdr.setLumi(event.luminosityBlock());
  mEvtHdr.setBunch(event.bunchCrossing());
  //-------------- Generator Info -------------------------------------
  Handle<GenEventInfoProduct> hEventInfo;
  event.getByLabel("generator", hEventInfo);
  mEvtHdr.setPthat(hEventInfo->binningValues()[0]);
  mEvtHdr.setWeight(hEventInfo->weight());

  Handle<GenJetCollection>  genjets;
  event.getByLabel(mGenJetsName,genjets);
    for(GenJetCollection::const_iterator i_gen = genjets->begin(); i_gen != genjets->end(); i_gen++) {
      if (i_gen->pt() > mMinGenPt && fabs(i_gen->y()) < mMaxY) {
        mGenJets.push_back(i_gen->p4());
      }
    }

  //-------------- fill the tree -------------------------------------
  mEvent->setEvtHdr(mEvtHdr);
  mEvent->setGenJets(mGenJets);
  if ( genjets->size() > 0 )
    mTree->Fill();
}
//////////////////////////////////////////////////////////////////////////////////////////
ProcessedTreeProducerGEN::~ProcessedTreeProducerGEN()
{
}

DEFINE_FWK_MODULE(ProcessedTreeProducerGEN);
