#ifndef ProcessedTreeProducerBTag_h
#define ProcessedTreeProducerBTag_h

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/JetCorrFactors.h"
#include "SMPJ/AnalysisFW/interface/QCDJet.h"
#include "SMPJ/AnalysisFW/interface/QCDEvent.h"
#include "SMPJ/AnalysisFW/interface/QCDEventHdr.h"
#include "SMPJ/AnalysisFW/interface/QCDPFJet.h"
#include "SMPJ/AnalysisFW/interface/QCDMET.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

using namespace edm;
using namespace reco;
using namespace std;
using namespace trigger;

class ProcessedTreeProducerBTag : public edm::EDAnalyzer
{
  public:
    typedef reco::Particle::LorentzVector LorentzVector;
    explicit ProcessedTreeProducerBTag(edm::ParameterSet const& cfg);
    virtual void beginJob();
    virtual void beginRun(edm::Run const &, edm::EventSetup const& iSetup);
    virtual void analyze(edm::Event const& evt, edm::EventSetup const& iSetup);
    virtual void endJob();
    virtual ~ProcessedTreeProducerBTag();
  private:
    void buildTree();
    static bool sort_pfjets(QCDPFJet j1, QCDPFJet j2) {
      return j1.ptCor() > j2.ptCor();
    }
    //---- configurable parameters --------
    bool   mIsMCarlo;
    bool   mUseGenInfo;
    bool   mPrintTriggerMenu;
    bool   isPFJecUncSet_,isPFJecUncSetCHS_;
    int    mGoodVtxNdof,mMinNPFJets;
    double mGoodVtxZ;
    double mMinPFPt,mMinPFFatPt,mMaxPFFatEta,mMinGenPt,mMaxY,mMinJJMass;
    std::string mPFJECservice;
    std::string mPFPayloadName;
    std::string mPFPayloadNameCHS;
    std::string pfpujetid;
	std::string pfchsjetpuid;

    // unc file for non CHS jet ---- //
    std::string mPFJECUncSrc;
    // unc file for CHS jet ---- //
    std::string mPFJECUncSrcCHS;
    std::vector<std::string> mPFJECUncSrcNames;
    // ---- non CHS jet input tag ----- //
    edm::InputTag mPFJetsName;
    // ----CHS jet input tag ----- //
    edm::InputTag mPFJetsNameCHS;
    edm::InputTag mGenJetsName;
    edm::InputTag mOfflineVertices;
    edm::InputTag mSrcCaloRho;
    edm::InputTag mSrcPFRho;
    edm::InputTag mSrcPU;
    //edm::InputTag mPFMET;
    edm::EDGetTokenT<pat::METCollection> mPFMET;
    //edm::InputTag mHBHENoiseFilter;
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
    // ---- non CHS jet uncertainty ------ //
    JetCorrectionUncertainty *mPFUnc;
    // ---- non CHS jet uncertainty ------ //
    JetCorrectionUncertainty *mPFUncCHS;
    //------- non CHS jet uncertainty sources -------- //
    std::vector<JetCorrectionUncertainty*> mPFUncSrc;
    // -------- CHS jet uncertainty sources -------- //
    std::vector<JetCorrectionUncertainty*> mPFUncSrcCHS;

    edm::Service<TFileService> fs;
    TTree *mTree;
    TH1F *mTriggerPassHisto,*mTriggerNamesHisto;
    //---- TREE variables --------
    QCDEvent *mEvent;
};

#endif
