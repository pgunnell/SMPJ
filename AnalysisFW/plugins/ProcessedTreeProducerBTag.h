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
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"

#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

//Hadron level definition
#include "SimDataFormats/JetMatching/interface/JetFlavourInfo.h"              
#include "SimDataFormats/JetMatching/interface/JetFlavourInfoMatching.h"
#include "DataFormats/Math/interface/deltaR.h"

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
    bool   mAK4;
    bool   mUseGenInfo;
    bool   mPrintTriggerMenu;
    int    mGoodVtxNdof,mMinNPFJets;
    double mGoodVtxZ;
    double mMinPFPt,mMinGenPt,mMaxY;
    std::string pfchsjetpuid;

    // ---- non CHS jet input tag ----- //
    edm::EDGetTokenT<reco::VertexCollection> mOfflineVertices;
    edm::EDGetTokenT<reco::BeamSpot> mBeamSpot;
    edm::EDGetTokenT<edm::View<pat::Jet> >mPFJetsNameCHS;
    edm::EDGetTokenT<GenEventInfoProduct> mhEventInfo;
    edm::EDGetTokenT<edm::ValueMap<float>> qgToken;
    // ----CHS jet input tag ----- //
    edm::EDGetTokenT<double> mSrcCaloRho;
    edm::EDGetTokenT<double> mSrcPFRho;
    edm::EDGetTokenT<pat::METCollection> mPFMET;
    edm::EDGetTokenT<GenJetCollection> mGenJetsName;
    edm::EDGetTokenT<reco::GenParticleCollection> mgenParticles;
    //edm::InputTag mHBHENoiseFilter;
    //---- TRIGGER -------------------------
    std::string   processName_;
    std::vector<std::string> triggerNames_;
    std::vector<unsigned int> triggerIndex_;
    //edm::InputTag mSrcPU;
    edm::EDGetTokenT<edm::TriggerResults> triggerResultsTag_;
    edm::EDGetTokenT<trigger::TriggerEvent> triggerEventTag_;
    edm::EDGetTokenT<std::vector<PileupSummaryInfo> > mSrcPU;
    edm::Handle<edm::TriggerResults>   triggerResultsHandle_;
    edm::Handle<trigger::TriggerEvent> triggerEventHandle_;

    edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
    edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
    edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;
    edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescalesL1Min_;
    edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescalesL1Max_;

    edm::EDGetTokenT<GenEventInfoProduct> genEvtInfoToken;
    edm::EDGetTokenT<LHEEventProduct> lheEvtInfoToken;

    edm::Handle<GenEventInfoProduct> genEvtInfo;
    edm::Handle<LHEEventProduct> lheEvtInfo;

    //hadron jet definition
    edm::EDGetTokenT<reco::JetFlavourInfoMatchingCollection> jetFlavourInfosToken_;

    std::string mPFJECUncSrcCHS;
    std::vector<std::string> mPFJECUncSrcNames;
    std::string mPFPayloadNameCHS;

    bool saveWeights_;

    HLTConfigProvider hltConfig_;
    //---- CORRECTORS ----------------------
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

    int getMatchedPartonGen(edm::Event const& event, GenJetCollection::const_iterator i_gen);
    int getMatchedHadronGen(edm::Event const& event, GenJetCollection::const_iterator i_gen);

    HLTPrescaleProvider hltPrescale_;

};

#endif
