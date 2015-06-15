#ifndef ProcessedTreeProducerPFCalo_h
#define ProcessedTreeProducerPFCalo_h

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "SMPJ/AnalysisFW/interface/QCDJet.h"
#include "SMPJ/AnalysisFW/interface/QCDEvent.h"
#include "SMPJ/AnalysisFW/interface/QCDEventHdr.h"
#include "SMPJ/AnalysisFW/interface/QCDPFJet.h"
#include "SMPJ/AnalysisFW/interface/QCDCaloJet.h"
#include "SMPJ/AnalysisFW/interface/QCDMET.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

using namespace edm;
using namespace reco;
using namespace std;
using namespace trigger;

class ProcessedTreeProducerPFCalo : public edm::EDAnalyzer 
{
  public:
    typedef reco::Particle::LorentzVector LorentzVector;
    explicit ProcessedTreeProducerPFCalo(edm::ParameterSet const& cfg);
    virtual void beginJob();
    virtual void beginRun(edm::Run const &, edm::EventSetup const& iSetup);
    virtual void analyze(edm::Event const& evt, edm::EventSetup const& iSetup);
    virtual void endJob();
    virtual ~ProcessedTreeProducerPFCalo();
  private:  
    void buildTree();
    static bool sort_calojets(QCDCaloJet j1, QCDCaloJet j2) {
      return j1.ptCor() > j2.ptCor();
    }
    static bool sort_pfjets(QCDPFJet j1, QCDPFJet j2) {
      return j1.ptCor() > j2.ptCor();
    }
    //---- configurable parameters --------  
    bool   mIsMCarlo;
    bool   mUseGenInfo;
    bool   mPrintTriggerMenu;
    bool   isPFJecUncSet_,isCaloJecUncSet_;
    int    mGoodVtxNdof,mMinNCaloJets,mMinNPFJets;
    double mGoodVtxZ; 
    double mMinCaloPt,mMinPFPt,mMinGenPt,mMaxY,mMinJJMass;
    std::string mCaloJECservice;
    std::string mPFJECservice;
    std::string mPFPayloadName;
    std::string mCaloPayloadName;
    std::string mPFJECUncSrc;
    std::vector<std::string> mPFJECUncSrcNames;
    edm::InputTag mCaloJetsName;
    edm::InputTag mPFJetsName;
    edm::InputTag mGenJetsName;
    edm::InputTag mCaloJetID;
    edm::InputTag mCaloJetExtender;
    edm::InputTag mOfflineVertices;
    edm::InputTag mSrcCaloRho;
    edm::InputTag mSrcPFRho;
    edm::InputTag mSrcPU;
    //---- TRIGGER -------------------------
    std::string   processName_;
    std::vector<std::string> triggerNames_;
    std::vector<unsigned int> triggerIndex_;
    edm::InputTag triggerResultsTag_;
    edm::InputTag triggerEventTag_;
    edm::Handle<edm::TriggerResults>   triggerResultsHandle_;
    edm::Handle<trigger::TriggerEvent> triggerEventHandle_;
    HLTConfigProvider hltConfig_;
    //---- CORRECTORS ----------------------
    const JetCorrector *mPFJEC;
    const JetCorrector *mCALOJEC;
    JetCorrectionUncertainty *mCALOUnc;
    JetCorrectionUncertainty *mPFUnc;
    std::vector<JetCorrectionUncertainty*> mPFUncSrc;
    
    edm::Service<TFileService> fs;
    TTree *mTree;
    TH1F *mTriggerPassHisto,*mTriggerNamesHisto; 
    //---- TREE variables --------
    QCDEvent *mEvent;
};

#endif
