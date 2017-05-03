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

#include "SMPJ/AnalysisFW/plugins/ProcessedTreeProducerBTag.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/JetExtendedAssociation.h"
#include "DataFormats/JetReco/interface/JetID.h"
#include "DataFormats/METReco/interface/HcalNoiseSummary.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"

#include "PhysicsTools/PatUtils/interface/bJetSelector.h"
#include "PhysicsTools/PatExamples/interface/BTagPerformance.h"
#include "PhysicsTools/PatExamples/interface/PatBTagCommonHistos.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Common/interface/ValueMap.h"

//hadron-level definition
#include "SimDataFormats/JetMatching/interface/JetFlavourInfo.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourInfoMatching.h"
#include "DataFormats/Math/interface/deltaR.h"

ProcessedTreeProducerBTag::ProcessedTreeProducerBTag(edm::ParameterSet const& cfg):
  mOfflineVertices(mayConsume<reco::VertexCollection>(cfg.getParameter<edm::InputTag>("offlineVertices"))),
  mBeamSpot(mayConsume<reco::BeamSpot>(cfg.getParameter<edm::InputTag>("beamSpot"))),
  //mPFJetsName(mayConsume<pat::Jet>(cfg.getParameter<edm::InputTag>("pfjets"))),
  //mPFJetsNameCHS(mayConsume<pat::Jet>(cfg.getParameter<edm::InputTag>("pfjetschs"))),
  mSrcCaloRho(mayConsume<double>(cfg.getParameter<edm::InputTag>("srcCaloRho"))),
  mSrcPFRho(mayConsume<double>(cfg.getParameter<edm::InputTag>("srcPFRho"))),
  mPFMET(mayConsume<pat::METCollection>(cfg.getParameter<edm::InputTag>("pfmet"))),
  mGenJetsName(mayConsume<GenJetCollection>(cfg.getUntrackedParameter<edm::InputTag>("genjets",edm::InputTag("")))),
  triggerResultsTag_(mayConsume<edm::TriggerResults>(cfg.getParameter<edm::InputTag>("triggerResults"))),
  triggerEventTag_(mayConsume<trigger::TriggerEvent>(cfg.getParameter<edm::InputTag>("triggerEvent"))),
  mHBHENoiseFilterResultLabel(mayConsume<bool>(cfg.getParameter<edm::InputTag>("HBHENoiseFilterResultLabel"))),
  mHBHENoiseFilterResultNoMinZLabel(mayConsume<bool>(cfg.getParameter<edm::InputTag>("HBHENoiseFilterResultNoMinZLabel"))),
  mSrcPU(mayConsume<std::vector<PileupSummaryInfo> >(cfg.getUntrackedParameter<edm::InputTag>("srcPULabel"))),
  hltPrescale_(cfg, consumesCollector(), *this)// ",edm::InputTag("addPileupInfo"))))
{
//  mPFJECservice      = cfg.getParameter<std::string>               ("pfjecService");
  //mPFPayloadName = cfg.getParameter<std::string>("PFPayloadName");
  mPFPayloadNameCHS  = cfg.getParameter<std::string>               ("PFPayloadNameCHS");
  //pfpujetid       = cfg.getParameter<std::string>               ("pfpujetid");
  pfchsjetpuid    = cfg.getParameter<std::string>               ("pfchsjetpuid");
  mGoodVtxNdof       = cfg.getParameter<double>                    ("goodVtxNdof");
  mGoodVtxZ          = cfg.getParameter<double>                    ("goodVtxZ");
  mMinPFPt           = cfg.getParameter<double>                    ("minPFPt");
  mMinJJMass         = cfg.getParameter<double>                    ("minJJMass");
  mMaxY              = cfg.getParameter<double>                    ("maxY");
  mMinNPFJets        = cfg.getParameter<int>                       ("minNPFJets");
  //mOfflineVertices   = cfg.getParameter<edm::InputTag>             ("offlineVertices");
  mPrintTriggerMenu  = cfg.getUntrackedParameter<bool>             ("printTriggerMenu",false);
  mIsMCarlo          = cfg.getUntrackedParameter<bool>             ("isMCarlo",false);
  mAK4          = cfg.getUntrackedParameter<bool>             ("AK4",false);
  mUseGenInfo        = cfg.getUntrackedParameter<bool>             ("useGenInfo",false);
  mMinGenPt          = cfg.getUntrackedParameter<double>           ("minGenPt",30);
  processName_       = cfg.getParameter<std::string>               ("processName");
  triggerNames_      = cfg.getParameter<std::vector<std::string> > ("triggerName");
  //mPFJECUncSrc       = cfg.getParameter<std::string>               ("jecUncSrc");
  mPFJECUncSrcCHS    = cfg.getParameter<std::string>               ("jecUncSrcCHS");
  mPFJECUncSrcNames  = cfg.getParameter<std::vector<std::string> > ("jecUncSrcNames");
  //mPFJetsName = consumes<edm::View<pat::Jet> >(cfg.getParameter<edm::InputTag>("pfjets"));
  mPFJetsNameCHS = consumes<edm::View<pat::Jet> >(cfg.getParameter<edm::InputTag>("pfjetschs"));
  mhEventInfo = consumes<GenEventInfoProduct>(cfg.getParameter<edm::InputTag>("EventInfo"));
  mgenParticles = consumes<reco::GenParticleCollection>(cfg.getParameter<edm::InputTag>("GenParticles"));
  qgToken = consumes<edm::ValueMap<float>>(edm::InputTag("QGTagger", "qgLikelihood"));
  jetFlavourInfosToken_ = consumes<reco::JetFlavourInfoMatchingCollection>( cfg.getParameter<edm::InputTag>("jetFlavourInfos"));
  //mPFJetsName        = cfg.getParameter<edm::InputTag>             ("pfjets");
  //mPFJetsNameCHS     = cfg.getParameter<edm::InputTag>             ("pfjetschs");
  //mSrcPU = cfg.getUntrackedParameter<edm::InputTag> ("srcPU",edm::InputTag("addPileupInfo"));
  //triggerResultsTag_ = cfg.getParameter<edm::InputTag>             ("triggerResults");
  //triggerEventTag_   = cfg.getParameter<edm::InputTag>             ("triggerEvent");
  //New additions
  //beamSpot_(consumes<reco::BeamSpot>(cfg.getParameter<std::string>("offlineBeamSpot")));
}
//////////////////////////////////////////////////////////////////////////////////////////
void ProcessedTreeProducerBTag::beginJob()
{
  mTree = fs->make<TTree>("ProcessedTree","ProcessedTree");
  mEvent = new QCDEvent();
  mTree->Branch("events","QCDEvent",&mEvent);
  mTriggerNamesHisto = fs->make<TH1F>("TriggerNames","TriggerNames",1,0,1);
  mTriggerNamesHisto->SetBit(TH1::kUserContour);
  for(unsigned i=0;i<triggerNames_.size();i++)
    mTriggerNamesHisto->Fill(triggerNames_[i].c_str(),1);
  mTriggerPassHisto = fs->make<TH1F>("TriggerPass","TriggerPass",1,0,1);
  mTriggerPassHisto->SetBit(TH1::kUserContour);
  //isPFJecUncSet_ = false;
  isPFJecUncSetCHS_ = false;
}
//////////////////////////////////////////////////////////////////////////////////////////
void ProcessedTreeProducerBTag::endJob()
{
}
//////////////////////////////////////////////////////////////////////////////////////////
void ProcessedTreeProducerBTag::beginRun(edm::Run const & iRun, edm::EventSetup const& iSetup)
{
  bool changed(true);
  if (hltConfig_.init(iRun,iSetup,processName_,changed) && hltPrescale_.init(iRun, iSetup, processName_, changed) ) {
    if (changed) {
      // check if trigger names in (new) config
      cout<<"New trigger menu found !!!"<<endl;
      triggerIndex_.clear();
      const unsigned int n(hltConfig_.size());
      for(unsigned itrig=0;itrig<triggerNames_.size();itrig++) {
        triggerIndex_.push_back(hltConfig_.triggerIndex(triggerNames_[itrig]));
        cout<<triggerNames_[itrig]<<" "<<triggerIndex_[itrig]<<" ";
        if (triggerIndex_[itrig] >= n)
          cout<<"does not exist in the current menu"<<endl;
        else
          cout<<"exists"<<endl;
      }
      cout << "Available TriggerNames are: " << endl;
      if (mPrintTriggerMenu)
        hltConfig_.dump("Triggers");
    }
  }
  else {
    cout << "ProcessedTreeProducerBTag::analyze:"
         << " config extraction failure with process name "
         << processName_ << endl;
  }
}
//////////////////////////////////////////////////////////////////////////////////////////
void ProcessedTreeProducerBTag::analyze(edm::Event const& event, edm::EventSetup const& iSetup)
{
  vector<QCDPFJet>      mPFJets;
  vector<QCDPFJet>      mPFJetsCHS;
  vector<LorentzVector> mGenJets;
  vector<float> GenFlavour;
  vector<float> GenHadronFlavour;
  QCDEventHdr mEvtHdr;
  QCDMET mPFMet;

  //-------------- Basic Event Info ------------------------------
  mEvtHdr.setRun(event.id().run());
  mEvtHdr.setEvt(event.id().event());
  mEvtHdr.setLumi(event.luminosityBlock());
  mEvtHdr.setBunch(event.bunchCrossing());
  //-------------- Beam Spot --------------------------------------
  Handle<reco::BeamSpot> beamSpot;
  event.getByToken(mBeamSpot,beamSpot);
  if (beamSpot.isValid())
    mEvtHdr.setBS(beamSpot->x0(),beamSpot->y0(),beamSpot->z0());
  else
    mEvtHdr.setBS(-999,-999,-999);


  //-------------- HCAL Noise Summary -----------------------------
  Handle<bool> noiseSummary;
  Handle<bool> noiseSummary_NoMinZ;

  if (!mIsMCarlo) {
    // event.getByLabel(mHBHENoiseFilter,noiseSummary);
    event.getByToken(mHBHENoiseFilterResultLabel, noiseSummary);
    //event.getByToken(mHBHENoiseFilterResultProducer, noiseSummary);
    mEvtHdr.setHCALNoise(*noiseSummary);

    event.getByToken(mHBHENoiseFilterResultNoMinZLabel, noiseSummary_NoMinZ);
    mEvtHdr.setHCALNoiseNoMinZ(*noiseSummary_NoMinZ);
    
  }
  else{ 
   mEvtHdr.setHCALNoise(true);
   mEvtHdr.setHCALNoiseNoMinZ(true);
  }
  //-------------- Trigger Info -----------------------------------
  event.getByToken(triggerResultsTag_,triggerResultsHandle_);
  if (!triggerResultsHandle_.isValid()) {
    cout << "ProcessedTreeProducerBTag::analyze: Error in getting TriggerResults product from Event!" << endl;
    return;
  }
  event.getByToken(triggerEventTag_,triggerEventHandle_);
  if (!triggerEventHandle_.isValid()) {
    cout << "ProcessedTreeProducerBTag::analyze: Error in getting TriggerEvent product from Event!" << endl;
    return;
  }
  vector<int> L1Prescales,HLTPrescales,Fired;
  vector<vector<LorentzVector> > mL1Objects,mHLTObjects;
  // sanity check
  assert(triggerResultsHandle_->size() == hltConfig_.size());
  //------ loop over all trigger names ---------
  for(unsigned itrig=0;itrig<triggerNames_.size() && !mIsMCarlo;itrig++) {
    bool accept(false);
    int preL1(-1);
    int preHLT(-1);
    int tmpFired(-1);
    vector<LorentzVector> vvL1,vvHLT;
    if (triggerIndex_[itrig] < hltConfig_.size()) {
      accept = triggerResultsHandle_->accept(triggerIndex_[itrig]);
      //      const std::pair<int,int> prescales(hltConfig_.prescaleValues(event,iSetup,triggerNames_[itrig]));

      ///In detail
      //get prescale info from hltConfig_
      
      std::pair<std::vector<std::pair<std::string,int> >,int> detailedPrescaleInfo = hltPrescale_.prescaleValuesInDetail(event, iSetup, triggerNames_[itrig]);	 
      preHLT = detailedPrescaleInfo.second ;

      // save l1 prescale values in standalone vector
      std::vector <int> l1prescalevals;
      for( size_t varind = 0; varind < detailedPrescaleInfo.first.size(); varind++ ){
	l1prescalevals.push_back(detailedPrescaleInfo.first.at(varind).second);
      }
   
      //find and save minimum l1 prescale of any ORed L1 that seeds the HLT
      std::vector<int>::iterator result = std::min_element(std::begin(l1prescalevals), std::end(l1prescalevals));
      size_t minind = std::distance(std::begin(l1prescalevals), result);
      // sometimes there are no L1s associated with a HLT. In that case, this branch stores -1 for the l1prescale
      preL1 = minind < l1prescalevals.size() ? l1prescalevals.at(minind) : -1 ;//commented for 76X
      
      ///end in detail
      if (!accept)
        tmpFired = 0;
      else {
        mTriggerPassHisto->Fill(triggerNames_[itrig].c_str(),1);
        tmpFired = 1;
      }

      //--------- modules on this trigger path--------------
      const vector<string>& moduleLabels(hltConfig_.moduleLabels(triggerIndex_[itrig]));
      const unsigned int moduleIndex(triggerResultsHandle_->index(triggerIndex_[itrig]));
      bool foundL1(false);
      for(unsigned int j=0; j<=moduleIndex; ++j) {
        const string& moduleLabel(moduleLabels[j]);
        const string  moduleType(hltConfig_.moduleType(moduleLabel));
        //--------check whether the module is packed up in TriggerEvent product
        const unsigned int filterIndex(triggerEventHandle_->filterIndex(InputTag(moduleLabel,"",processName_)));
        if (filterIndex<triggerEventHandle_->sizeFilters()) {
          const Vids& VIDS (triggerEventHandle_->filterIds(filterIndex));
          const Keys& KEYS(triggerEventHandle_->filterKeys(filterIndex));
          const size_type nI(VIDS.size());
          const size_type nK(KEYS.size());
          assert(nI==nK);
          const size_type n(max(nI,nK));
          const TriggerObjectCollection& TOC(triggerEventHandle_->getObjects());
          if (foundL1) {
            for(size_type i=0; i!=n; ++i) {
              const TriggerObject& TO(TOC[KEYS[i]]);
              TLorentzVector P4;
              P4.SetPtEtaPhiM(TO.pt(),TO.eta(),TO.phi(),TO.mass());
              LorentzVector qcdhltobj(P4.Px(),P4.Py(),P4.Pz(),P4.E());
              vvHLT.push_back(qcdhltobj);
              //cout<<TO.pt()<<endl;
            }
          }
          else {
            for(size_type i=0; i!=n; ++i) {
              const TriggerObject& TO(TOC[KEYS[i]]);
              TLorentzVector P4;
              P4.SetPtEtaPhiM(TO.pt(),TO.eta(),TO.phi(),TO.mass());
              LorentzVector qcdl1obj(P4.Px(),P4.Py(),P4.Pz(),P4.E());
              vvL1.push_back(qcdl1obj);
              //cout<<TO.pt()<<endl;
            }
            foundL1 = true;
          }
        }
      }// loop over modules
    }// if the trigger exists in the menu
    //cout<<triggerNames_[itrig]<<" "<<triggerIndex_[itrig]<<" "<<accept<<" "<<tmpFired<<endl;
    Fired.push_back(tmpFired);
    L1Prescales.push_back(preL1);
    HLTPrescales.push_back(preHLT);
    mL1Objects.push_back(vvL1);
    mHLTObjects.push_back(vvHLT);
  }// loop over trigger names
  mEvent->setTrigDecision(Fired);
  mEvent->setPrescales(L1Prescales,HLTPrescales);
  mEvent->setL1Obj(mL1Objects);
  mEvent->setHLTObj(mHLTObjects);

  //-------------- Vertex Info -----------------------------------
  Handle<reco::VertexCollection> recVtxs;
  event.getByToken(mOfflineVertices,recVtxs);
  //------------- reject events without reco vertices ------------
  int VtxGood(0);
  bool isPVgood(false);
  float PVx(0),PVy(0),PVz(0),PVndof(0);
  for(VertexCollection::const_iterator i_vtx = recVtxs->begin(); i_vtx != recVtxs->end(); i_vtx++) {
    int index = i_vtx-recVtxs->begin();
    if (index == 0) {
      PVx    = i_vtx->x();
      PVy    = i_vtx->y();
      PVz    = i_vtx->z();
      PVndof = i_vtx->ndof();
    }
    if (!(i_vtx->isFake()) && i_vtx->ndof() >= mGoodVtxNdof && fabs(i_vtx->z()) <= mGoodVtxZ) {
      if (index == 0) {
        isPVgood = true;
      }
      VtxGood++;
    }
  }
  mEvtHdr.setVertices(recVtxs->size(),VtxGood);
  mEvtHdr.setPV(isPVgood,PVndof,PVx,PVy,PVz);
  //-------------- Rho ------------------------------------------------
  Handle<double> rhoCalo;
  event.getByToken(mSrcCaloRho,rhoCalo);
  Handle<double> rhoPF;
  event.getByToken(mSrcPFRho,rhoPF);
  mEvtHdr.setRho(*rhoCalo,*rhoPF);
  //-------------- Generator Info -------------------------------------
  Handle<GenEventInfoProduct> hEventInfo;
  //-------------- Simulated PU Info ----------------------------------
  Handle<std::vector<PileupSummaryInfo> > PupInfo;
  if (mIsMCarlo && mUseGenInfo) {
    event.getByToken(mhEventInfo, hEventInfo);
    if(hEventInfo->hasBinningValues())
    mEvtHdr.setPthat(hEventInfo->binningValues()[0]);
    else
    mEvtHdr.setPthat(0);

    mEvtHdr.setWeight(hEventInfo->weight());
    event.getByToken(mSrcPU, PupInfo);
    std::vector<PileupSummaryInfo>::const_iterator PUI;
    int nbx = PupInfo->size();
    int ootpuEarly(0),ootpuLate(0),intpu(0);
    float Tnpv = -1.; // new variable for computing pileup weight factor for the event

    for(PUI = PupInfo->begin(); PUI != PupInfo->end(); ++PUI) {
      if (PUI->getBunchCrossing() < 0)
        ootpuEarly += PUI->getPU_NumInteractions();      
      else if (PUI->getBunchCrossing() > 0)
        ootpuLate += PUI->getPU_NumInteractions();
      else {
        intpu += PUI->getPU_NumInteractions();
        Tnpv = PUI->getTrueNumInteractions();
       }

    }

    mEvtHdr.setPU(nbx,ootpuEarly,ootpuLate,intpu);
    mEvtHdr.setTrPu(Tnpv);
  }
  else {
    mEvtHdr.setPthat(0);
    mEvtHdr.setWeight(0);
    mEvtHdr.setPU(0,0,0,0);
    mEvtHdr.setTrPu(0);
  }

  //---------------- Jets ---------------------------------------------
  //mPFJEC   = JetCorrector::getJetCorrector(mPFJECservice,iSetup);
 //event.getByToken(mvaFullPUDiscriminantToken_ ,puJetIdMva);
  
  /*   LET's remove PFjets // Engin 
  edm::ESHandle<JetCorrectorParametersCollection> PFJetCorParColl;
  if (mPFPayloadName != "" && !isPFJecUncSet_){
    iSetup.get<JetCorrectionsRecord>().get(mPFPayloadName,PFJetCorParColl);
    JetCorrectorParameters const& PFJetCorPar = (*PFJetCorParColl)["Uncertainty"];
    mPFUnc = new JetCorrectionUncertainty(PFJetCorPar);
    if (mPFJECUncSrc != "") {
      for(unsigned isrc=0;isrc<mPFJECUncSrcNames.size();isrc++) {
        JetCorrectorParameters *par = new JetCorrectorParameters(mPFJECUncSrc,mPFJECUncSrcNames[isrc]);
        JetCorrectionUncertainty *tmpUnc = new JetCorrectionUncertainty(*par);
        mPFUncSrc.push_back(tmpUnc);
      } // for(unsigned isrc=0;isrc<mPFJECUncSrcNames.size();isrc++)
    } // if (mPFJECUncSrc != "")
    isPFJecUncSet_ = true;
  } // if (mPFPayloadName != "" && !isPFJecUncSet_)
  */

  Handle<GenJetCollection>  genjets;
  if (mIsMCarlo) {
    event.getByToken(mGenJetsName,genjets);
    for(GenJetCollection::const_iterator i_gen = genjets->begin(); i_gen != genjets->end(); i_gen++) {
      if (i_gen->pt() > mMinGenPt && fabs(i_gen->y()) < mMaxY) {
        mGenJets.push_back(i_gen->p4());

	//ADD FLAVOUR AT GEN LEVEL
	int FlavourGen = getMatchedPartonGen(event,i_gen);
	//if(FlavourGen<-100) cout<<FlavourGen<<" "<<i_gen->pt()<<" "<<i_gen->eta()<<" "<<i_gen->phi()<<endl;
	GenFlavour.push_back(FlavourGen);

	//OLD DEFINITION STRONG MODEL DEPENDENCE
	//int FlavourGenHadron = getMatchedHadronGen(event,i_gen);
	//if(FlavourGenHadron==5) cout<<FlavourGen<<" "<<i_gen->pt()<<" "<<i_gen->eta()<<" "<<i_gen->phi()<<endl;
	//GenHadronFlavour.push_back(FlavourGenHadron);
	

      }
    }

    edm::Handle<reco::JetFlavourInfoMatchingCollection> theJetFlavourInfos;
    event.getByToken(jetFlavourInfosToken_, theJetFlavourInfos );
    
    for ( reco::JetFlavourInfoMatchingCollection::const_iterator j  = theJetFlavourInfos->begin();j != theJetFlavourInfos->end();++j ) {
      //std::cout << "-------------------- Jet Flavour Info --------------------" << std::endl;
      
      //const reco::Jet *aJet = (*j).first.get();
      reco::JetFlavourInfo aInfo = (*j).second;
      //std::cout << std::setprecision(2) << std::setw(6) << std::fixed
      //<< "[printJetFlavourInfo] Jet " << (j - theJetFlavourInfos->begin()) << " pt, eta, rapidity, phi = " << aJet->pt() << ", "
      //<< aJet->eta() << ", "
      //<< aJet->rapidity() << ", "
      //<< aJet->phi()
      //<< std::endl;
      // ----------------------- Hadrons -------------------------------
      //std::cout << " Hadron-based flavour: " << aInfo.getHadronFlavour() << std::endl;
      
      int FlavourGenHadron = aInfo.getHadronFlavour();
      //if(FlavourGenHadron==5) cout<<FlavourGenHadron<<" "<<aJet->pt()<<" "<<aJet->eta()<<" "<<aJet->phi()<<" HADRONFLAV"<<endl;
      GenHadronFlavour.push_back(FlavourGenHadron);
    }
  }
  
  //----------- PFJets non CHS part -------------------------
  //Let's remove NON-CHS part  Engin 
  
  //edm::Handle<edm::View<pat::Jet> > patjets;
  //event.getByToken(mPFJetsName,patjets);

  /*edm::Handle<reco::JetTagCollection> btagDiscriminators;
    event.getByLabel("pfCombinedInclusiveSecondaryVertexV2BJetTags", btagDiscriminators);  */
  /*
  for(edm::View<pat::Jet>::const_iterator i_pfjet=patjets->begin(); i_pfjet!=patjets->end(); ++i_pfjet)
    {
      QCDPFJet qcdpfjet;
      
      if(i_pfjet->isPFJet() ){
	
	double scale = 1./i_pfjet->jecFactor(0); // --- the value of the JEC factor

	//---- preselection -----------------
	if (fabs(i_pfjet->y()) > mMaxY) continue;
	
	//---- vertex association -----------
	//---- get the vector of tracks -----
	reco::TrackRefVector vTrks(i_pfjet->associatedTracks());
	float sumTrkPt(0.0),sumTrkPtBeta(0.0),sumTrkPtBetaStar(0.0),beta(0.0),betaStar(0.0);
	int mpuTrk(0), mlvTrk(0); // # of pile-up tracks & lead-vertex tracks ## Juska
	int mjtTrk(0); // multiplicity of _all_ tracks in jet (also vtx-unassociated!) ## Juska
	
       //---- loop over the tracks of the jet ----
       //std::cout << "starting the loop yo!" << std::endl; // debug
       //std::cout << "vTrks.size()" << vTrks.size() << std::endl;
	for(reco::TrackRefVector::const_iterator i_trk = vTrks.begin(); i_trk != vTrks.end(); i_trk++) {
	  if (recVtxs->size() == 0) break;
	  sumTrkPt += (*i_trk)->pt();
	  mjtTrk++; //Juska
	  //---- loop over all vertices ----------------------------
	  for(unsigned ivtx = 0;ivtx < recVtxs->size();ivtx++) {
	    //---- loop over the tracks associated with the vertex ---
	    if (!((*recVtxs)[ivtx].isFake()) && (*recVtxs)[ivtx].ndof() >= mGoodVtxNdof && fabs((*recVtxs)[ivtx].z()) <= mGoodVtxZ) {
	      for(reco::Vertex::trackRef_iterator i_vtxTrk = (*recVtxs)[ivtx].tracks_begin(); i_vtxTrk != (*recVtxs)[ivtx].tracks_end(); ++i_vtxTrk) {
		//---- match the jet track to the track from the vertex ----
		reco::TrackRef trkRef(i_vtxTrk->castTo<reco::TrackRef>());
		//---- check if the tracks match -------------------------
		if (trkRef == (*i_trk)) {
		  if (ivtx == 0) {
		    sumTrkPtBeta += (*i_trk)->pt();
		    mlvTrk++; //Juska
		  }
		  else {
		    sumTrkPtBetaStar += (*i_trk)->pt();
		    mpuTrk++; //Juska
		  }
		  break;
		} // if (trkRef == (*i_trk))
	      } // for(reco::Vertex::trackRef_iterator i_vtxTrk = (*recVtxs)[ivtx].tracks_begin(); i_vtxTrk != (*recVtxs)[ivtx].tracks_end(); ++i_vtxTrk)
	    } // if (!((*recVtxs)[ivtx].isFake()) && (*recVtxs)[ivtx].ndof() >= mGoodVtxNdof && fabs((*recVtxs)[ivtx].z()) <= mGoodVtxZ)
          } // for(unsigned ivtx = 0;ivtx < recVtxs->size();ivtx++)
	} // for(reco::TrackRefVector::const_iterator i_trk = vTrks.begin(); i_trk != vTrks.end(); i_trk++)
	if (sumTrkPt > 0) {
	  beta     = sumTrkPtBeta/sumTrkPt;
	  betaStar = sumTrkPtBetaStar/sumTrkPt;
	} //if (sumTrkPt > 0)
        qcdpfjet.setBeta(beta);
        qcdpfjet.setBetaStar(betaStar);
	
	//---- jec uncertainty --------------
	double unc(0.0);
	vector<float> uncSrc(0);
	if (mPFPayloadName != "") {
	  mPFUnc->setJetEta(i_pfjet->eta());
	  mPFUnc->setJetPt(i_pfjet->pt());
	  unc = mPFUnc->getUncertainty(true);
	} // if (mPFPayloadName != "")
	if (mPFJECUncSrc != "") {
	  for(unsigned isrc=0;isrc<mPFJECUncSrcNames.size();isrc++) {
	    mPFUncSrc[isrc]->setJetEta(i_pfjet->eta());
	    mPFUncSrc[isrc]->setJetPt(i_pfjet->pt());
	    float unc1 = mPFUncSrc[isrc]->getUncertainty(true);
	    uncSrc.push_back(unc1);
	  } // for(unsigned isrc=0;isrc<mPFJECUncSrcNames.size();isrc++)
	} // if (mPFJECUncSrc != "")
	
	qcdpfjet.setP4(i_pfjet->p4());
	qcdpfjet.setCor(scale);
	qcdpfjet.setUnc(unc);
	qcdpfjet.setUncSrc(uncSrc);
	qcdpfjet.setArea(i_pfjet->jetArea());
      
	double chf   = i_pfjet->chargedHadronEnergyFraction();
	double nhf   = i_pfjet->neutralHadronEnergyFraction(); //+ i_pfjet->HFHadronEnergyFraction();
	double nemf   = i_pfjet->neutralEmEnergyFraction(); // equals to old phf with HF info included
	//double elf   = i_pfjet->electronEnergyFraction(); equals to cemf
	double cemf  = i_pfjet->chargedEmEnergyFraction();
	double muf   = i_pfjet->muonEnergyFraction();
	double hf_hf = i_pfjet->HFHadronEnergyFraction();
	double hf_phf= i_pfjet->HFEMEnergyFraction();
      
	int hf_hm    = i_pfjet->HFHadronMultiplicity();
	int hf_phm   = i_pfjet->HFEMMultiplicity();
	int chm      = i_pfjet->chargedHadronMultiplicity();
	int nhm      = i_pfjet->neutralHadronMultiplicity();
	int phm      = i_pfjet->photonMultiplicity();
	int elm      = i_pfjet->electronMultiplicity();
	int mum      = i_pfjet->muonMultiplicity();
	int npr      = i_pfjet->chargedMultiplicity() + i_pfjet->neutralMultiplicity();
	//bool looseID  = (npr>1 && phf<0.99 && nhf<0.99 && ((fabs(i_pfjet->eta())<=2.4 && elf<0.99 && chf>0 && chm>0) || fabs(i_pfjet->eta())>2.4));
	//bool tightID  = (npr>1 && phf<0.99 && nhf<0.99 && ((fabs(i_pfjet->eta())<=2.4 && nhf<0.9 && phf<0.9 && elf<0.99 && chf>0 && chm>0) || fabs(i_pfjet->eta())>2.4));
        // https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID
	float eta    = i_pfjet->eta();
	int cm       = i_pfjet->chargedMultiplicity();
	bool looseID = (nhf<0.99 && nemf<0.99 && npr>1 && muf<0.8) && ((fabs(eta) <= 2.4 && chf>0 && cm>0 && cemf<0.99) || fabs(eta)>2.4);
	bool tightID = (nhf<0.90 && nemf<0.90 && npr>1 && muf<0.8) && ((fabs(eta)<=2.4 && chf>0 && cm>0 && cemf<0.90) || fabs(eta)>2.4);
	
	double TCHE = i_pfjet->bDiscriminator("trackCountingHighEffBJetTags");
	double TCHP = i_pfjet->bDiscriminator("trackCountingHighPurBJetTags");
	double TCHEpf = i_pfjet->bDiscriminator("pfTrackCountingHighEffBJetTags");
	double TCHPpf = i_pfjet->bDiscriminator("pfTrackCountingHighPurBJetTags");
	
	double SoftMuonTagByIP = i_pfjet->bDiscriminator("softPFMuonByIP3dBJetTags");
	double SoftElectronTagByIP = i_pfjet->bDiscriminator("softPFElectronByIP3dBJetTags");
	double SoftMuonTag = i_pfjet->bDiscriminator("softPFMuonBJetTags");
	double SoftElectronTag = i_pfjet->bDiscriminator("softPFElectronBJetTags");
	
	double SimpleSecVertexHE = i_pfjet->bDiscriminator("simpleSecondaryVertexHighEffBJetTags");
	double SimpleSecVertexHP = i_pfjet->bDiscriminator("simpleSecondaryVertexHighPurBJetTags");
	double SimpleSecVertexHEpf = i_pfjet->bDiscriminator("pfSimpleSecondaryVertexHighEffBJetTags");
	double SimpleSecVertexHPpf = i_pfjet->bDiscriminator("pfSimpleSecondaryVertexHighPurBJetTags");
	
	double CSV = i_pfjet->bDiscriminator("combinedSecondaryVertexBJetTags");
	double CSVpf = i_pfjet->bDiscriminator("pfCombinedSecondaryVertexBJetTags");
	double CinclSVpf = i_pfjet->bDiscriminator("pfCombinedInclusiveSecondaryVertexBJetTags");
	double CMVApf = i_pfjet->bDiscriminator("pfCombinedMVABJetTags");
	double CSVSoftLeptonpf = i_pfjet->bDiscriminator("pfCombinedSecondaryVertexSoftLeptonBJetTags");
	
	double CSVpfPositive = i_pfjet->bDiscriminator("pfPositiveCombinedSecondaryVertexBJetTags");
	double CSVpfNegative = i_pfjet->bDiscriminator("pfNegativeCombinedSecondaryVertexBJetTags");

	//the three recommended
	double pfJetProbabilityBJetTags=i_pfjet->bDiscriminator("pfJetProbabilityBJetTags");
	double pfCombinedInclusiveSecondaryVertexV2BJetTags= i_pfjet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
	double pfCombinedMVAV2BJetTags=i_pfjet->bDiscriminator("pfCombinedMVAV2BJetTags");

	qcdpfjet.setLooseID(looseID);
	qcdpfjet.setTightID(tightID);
	qcdpfjet.setFrac(chf,nhf,nemf,cemf,muf);
	qcdpfjet.setMulti(npr,chm,nhm,phm,elm,mum);
	qcdpfjet.setHFFrac(hf_hf,hf_phf);
	qcdpfjet.setHFMulti(hf_hm,hf_phm);

	float partonFlavour=-100;
	float hadronFlavour=-100;
	
	if (mIsMCarlo && mUseGenInfo) {
	  partonFlavour = i_pfjet->partonFlavour();
	  hadronFlavour = i_pfjet->hadronFlavour();
	}

	//if(i_pfjet->pt()>114) cout<<partonFlavour<<" "<<hadronFlavour<<" Reco "<<i_pfjet->pt()<<" "<<i_pfjet->eta()<<" "<<i_pfjet->phi()<<" "<<CSV<<endl;

	qcdpfjet.setFlavour(partonFlavour,hadronFlavour);
	
	double hof   = i_pfjet->hoEnergyFraction(); // Juska
	qcdpfjet.setVtxInfo(mpuTrk,mlvTrk,mjtTrk);
	qcdpfjet.setHO(hof);
	
	//Filling B-tag infos
	qcdpfjet.setTCHETag(TCHE,TCHP,TCHEpf,TCHPpf);
	qcdpfjet.setSoftLeptonTag(SoftMuonTagByIP,SoftElectronTagByIP,SoftMuonTag,SoftElectronTag);
	qcdpfjet.setSimpleSecondaryVertexTag(SimpleSecVertexHE,SimpleSecVertexHP,SimpleSecVertexHEpf,SimpleSecVertexHPpf);
	qcdpfjet.setCombinedSecondaryVertexTag(CSV,CSVpf,CinclSVpf,CSVSoftLeptonpf,CMVApf);
	qcdpfjet.setPositiveNegativeCSV(CSVpfPositive,CSVpfNegative);
	qcdpfjet.setTagRecommended(pfJetProbabilityBJetTags,pfCombinedInclusiveSecondaryVertexV2BJetTags,pfCombinedMVAV2BJetTags);
	
	float pileupJetId = -999;
	if ( i_pfjet->hasUserFloat(pfpujetid) )    pileupJetId = i_pfjet->userFloat(pfpujetid);
	qcdpfjet.SetPUJetId(pileupJetId);
	
	if (mIsMCarlo) {
	  GenJetCollection::const_iterator i_matched;
	  float rmin(999);
	  for(GenJetCollection::const_iterator i_gen = genjets->begin(); i_gen != genjets->end(); i_gen++) {
	    double deltaR = reco::deltaR(*i_pfjet,*i_gen);
	    if (deltaR < rmin) {
	      rmin = deltaR;
	      i_matched = i_gen;
	    }
	  }
	  if (genjets->size() == 0) {
	    LorentzVector tmpP4(0.0,0.0,0.0,0.0);
	    qcdpfjet.setGen(tmpP4,0);
	  }
	  else
	    qcdpfjet.setGen(i_matched->p4(),rmin);
	} // if (mIsMCarlo)
	else {
	  LorentzVector tmpP4(0.0,0.0,0.0,0.0);
	  qcdpfjet.setGen(tmpP4,0);
	}
	if (qcdpfjet.pt() >= mMinPFPt && qcdpfjet.ptCor() >= mMinPFPt/2.)
	  mPFJets.push_back(qcdpfjet);
	

      } // if(iJet->isPFJet() )
    } // --- end of non chs patjet iterator loop -------------------- //
  */

  // ========================******************************************===================== //
  
  // -------- CHS Uncertainty part ----------------//
  edm::ESHandle<JetCorrectorParametersCollection> PFJetCorParCollCHS;
  if (mPFPayloadNameCHS != "" && !isPFJecUncSetCHS_){
    iSetup.get<JetCorrectionsRecord>().get(mPFPayloadNameCHS,PFJetCorParCollCHS);
    JetCorrectorParameters const& PFJetCorParCHS = (*PFJetCorParCollCHS)["Uncertainty"];
    mPFUncCHS = new JetCorrectionUncertainty(PFJetCorParCHS);
    if (mPFJECUncSrcCHS != "") {
      for(unsigned isrc=0;isrc<mPFJECUncSrcNames.size();isrc++) {
        JetCorrectorParameters *parchs = new JetCorrectorParameters(mPFJECUncSrcCHS,mPFJECUncSrcNames[isrc]);
        JetCorrectionUncertainty *tmpUncCHS = new JetCorrectionUncertainty(*parchs);
        mPFUncSrcCHS.push_back(tmpUncCHS);
      } // for(unsigned isrc=0;isrc<mPFJECUncSrcNames.size();isrc++)
    } // if (mPFJECUncSrcCHS != "")
    isPFJecUncSetCHS_ = true;
  } // if (mPFPayloadNameCHS != "" && !isPFJecUncSetCHS_)


  //----------- PFJets  CHS part -------------------------

  edm::Handle<edm::ValueMap<float>> qgHandle; 
  event.getByToken(qgToken, qgHandle);
 
  edm::Handle<edm::View<pat::Jet> > patjetschs;
  event.getByToken(mPFJetsNameCHS,patjetschs);

  /*for(edm::View<pat::Jet>::const_iterator i_pfjet=patjetschs->begin(); i_pfjet!=patjetschs->end(); ++i_pfjet){
    edm::RefToBase<pat::Jet> jetRef(edm::Ref<edm::View<pat::Jet> >(patjetschs, i_pfjet - patjetschs->begin()));
    float qgLikelihood = (*qgHandle)[jetRef];
    
    cout<<qgLikelihood<<endl;
    }*/
    
  for(edm::View<pat::Jet>::const_iterator i_pfjetchs=patjetschs->begin(); i_pfjetchs!=patjetschs->end(); ++i_pfjetchs)
     {
       QCDPFJet qcdpfjetchs;
       
       if(i_pfjetchs->isPFJet() ){
	 
	 double scaleCHS = 1./i_pfjetchs->jecFactor(0); // --- the value of the JEC factor
	 
	 //---- preselection -----------------
	 if (fabs(i_pfjetchs->y()) > mMaxY) continue;
	 
	 //---- vertex association -----------
	 //---- get the vector of tracks -----
	 reco::TrackRefVector vTrksCHS(i_pfjetchs->associatedTracks());
	 float sumTrkPtCHS(0.0),sumTrkPtBetaCHS(0.0),sumTrkPtBetaStarCHS(0.0),betaCHS(0.0),betaStarCHS(0.0);
	 
	 // Dunno how useful these are in chs jets...
	 int mpuTrk(0), mlvTrk(0); // # of pile-up tracks & lead-vertex tracks ## Juska
	 int mjtTrk(0); // multiplicity of _all_ tracks in jet (also vtx-unassociated!) ## Juska
	 
	 //---- loop over the tracks of the jet ----
	 
	 for(reco::TrackRefVector::const_iterator i_trkchs = vTrksCHS.begin(); i_trkchs != vTrksCHS.end(); i_trkchs++) {
	   if (recVtxs->size() == 0) break;
	   sumTrkPtCHS += (*i_trkchs)->pt();
	   mjtTrk++; //Juska
	   //---- loop over all vertices ----------------------------
	   for(unsigned ivtx = 0;ivtx < recVtxs->size();ivtx++) {
	     //---- loop over the tracks associated with the vertex ---
	     if (!((*recVtxs)[ivtx].isFake()) && (*recVtxs)[ivtx].ndof() >= mGoodVtxNdof && fabs((*recVtxs)[ivtx].z()) <= mGoodVtxZ) {
	       for(reco::Vertex::trackRef_iterator i_vtxTrk = (*recVtxs)[ivtx].tracks_begin(); i_vtxTrk != (*recVtxs)[ivtx].tracks_end(); ++i_vtxTrk) {
		 //---- match the chsjet track to the track from the vertex ----
		 reco::TrackRef trkRef(i_vtxTrk->castTo<reco::TrackRef>());
		 //---- check if the tracks match -------------------------
		 if (trkRef == (*i_trkchs)) {
		   if (ivtx == 0) {
		     sumTrkPtBetaCHS += (*i_trkchs)->pt();
		     mlvTrk++; //Juska
		   }
		   else {
		     sumTrkPtBetaStarCHS += (*i_trkchs)->pt();
		     mpuTrk++; //Juska
		   }
		   break;
		 } // if (trkRef == (*i_trk))
	       } // for(reco::Vertex::trackRef_iterator i_vtxTrk = (*recVtxs)[ivtx].tracks_begin(); i_vtxTrk != (*recVtxs)[ivtx].tracks_end(); ++i_vtxTrk)
	     } // if (!((*recVtxs)[ivtx].isFake()) && (*recVtxs)[ivtx].ndof() >= mGoodVtxNdof && fabs((*recVtxs)[ivtx].z()) <= mGoodVtxZ)
	   } // for(unsigned ivtx = 0;ivtx < recVtxs->size();ivtx++)
         } // for(reco::TrackRefVector::const_iterator i_trk = vTrks.begin(); i_trk != vTrks.end(); i_trk++)
         if (sumTrkPtCHS > 0) {
	   betaCHS     = sumTrkPtBetaCHS/sumTrkPtCHS;
	   betaStarCHS = sumTrkPtBetaStarCHS/sumTrkPtCHS;
         } //if (sumTrkPt > 0)
	 qcdpfjetchs.setBeta(betaCHS);
	 qcdpfjetchs.setBetaStar(betaStarCHS);
	 
	 //---- jec uncertainty --------------
	 double uncCHS(0.0);
	 vector<float> uncSrcCHS(0);
	 if (mPFPayloadNameCHS != "") {
	   mPFUncCHS->setJetEta(i_pfjetchs->eta());
	   mPFUncCHS->setJetPt(i_pfjetchs->pt());
	   uncCHS = mPFUncCHS->getUncertainty(true);
	 } // if (mPFPayloadName != "")
	 if (mPFJECUncSrcCHS != "") {
	   for(unsigned isrc=0;isrc<mPFJECUncSrcNames.size();isrc++) {
	     mPFUncSrcCHS[isrc]->setJetEta(i_pfjetchs->eta());
	     mPFUncSrcCHS[isrc]->setJetPt(i_pfjetchs->pt());
	     float unc1 = mPFUncSrcCHS[isrc]->getUncertainty(true);
	     uncSrcCHS.push_back(unc1);
	   } // for(unsigned isrc=0;isrc<mPFJECUncSrcNames.size();isrc++)
	 } // if (mPFJECUncSrc != "")
	 
	 qcdpfjetchs.setP4(i_pfjetchs->p4());
	 qcdpfjetchs.setCor(scaleCHS);
	 qcdpfjetchs.setUnc(uncCHS);
	 qcdpfjetchs.setUncSrc(uncSrcCHS);
	 qcdpfjetchs.setArea(i_pfjetchs->jetArea());
	 
	 double chf   = i_pfjetchs->chargedHadronEnergyFraction();
	 double nhf   = i_pfjetchs->neutralHadronEnergyFraction();// + i_pfjetchs->HFHadronEnergyFraction();
	 double nemf  = i_pfjetchs->neutralEmEnergyFraction(); // equals to deprecated phf but has HF info too
	 double cemf  = i_pfjetchs->chargedEmEnergyFraction(); // equals to deprecated elf
	 double muf   = i_pfjetchs->muonEnergyFraction();
	 double hf_hf = i_pfjetchs->HFHadronEnergyFraction();
	 double hf_phf= i_pfjetchs->HFEMEnergyFraction();
	 int hf_hm    = i_pfjetchs->HFHadronMultiplicity();
	 int hf_phm   = i_pfjetchs->HFEMMultiplicity();
	 int chm      = i_pfjetchs->chargedHadronMultiplicity();
	 int nhm      = i_pfjetchs->neutralHadronMultiplicity();
	 int phm      = i_pfjetchs->photonMultiplicity();
	 int elm      = i_pfjetchs->electronMultiplicity();
	 int mum      = i_pfjetchs->muonMultiplicity();
	 int npr      = i_pfjetchs->chargedMultiplicity() + i_pfjetchs->neutralMultiplicity();
	 
	 
	 float eta    = i_pfjetchs->eta();
	 int cm       = i_pfjetchs->chargedMultiplicity();
     
     bool looseID, tightID;
     if ( fabs(eta) <= 2.7){
        looseID = (nhf<0.99 && nemf<0.99 && npr > 1) && ((fabs(eta)<=2.4 && chf>0 && chm>0 && cemf <0.99) || fabs(eta)>2.4);
        tightID = (nhf<0.90 && nemf<0.90 && npr > 1) && ((fabs(eta)<=2.4 && chf>0 && chm>0 && cemf <0.99) || fabs(eta)>2.4);
     }
     else if ( fabs(eta) <= 3.0){
        looseID = (nemf<0.90 && (npr - cm)>2 ); 
        tightID = (nemf<0.90 && (npr - cm)>2 ); 
     }
     else {
        looseID = (nemf<0.90 && (npr - cm)>10 ); 
        tightID = (nemf<0.90 && (npr - cm)>10 ); 
     }

        // https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID
   
   //OLD JET IDs, 
  //	 bool looseID = (nhf<0.99 && nemf<0.99 && npr>1 && muf<0.8) && ((fabs(eta) <= 2.4 && chf>0 && cm>0 && cemf<0.99) || fabs(eta)>2.4);
  //	 bool tightID = (nhf<0.90 && nemf<0.90 && npr>1 && muf<0.8) && ((fabs(eta)<=2.4 && chf>0 && cm>0 && cemf<0.90) || fabs(eta)>2.4);
	 
	 qcdpfjetchs.setLooseID(looseID);
	 qcdpfjetchs.setTightID(tightID);
	 qcdpfjetchs.setFrac(chf,nhf,nemf,cemf,muf);
	 qcdpfjetchs.setMulti(npr,chm,nhm,phm,elm,mum,cm);
	 qcdpfjetchs.setHFFrac(hf_hf,hf_phf);
	 qcdpfjetchs.setHFMulti(hf_hm,hf_phm);
	 
	 
	 double hof   = i_pfjetchs->hoEnergyFraction(); // Juska
	 qcdpfjetchs.setVtxInfo(mpuTrk,mlvTrk,mjtTrk);
	 qcdpfjetchs.setHO(hof);

	 /*
     double TCHE = i_pfjetchs->bDiscriminator("trackCountingHighEffBJetTags");
	 double TCHP = i_pfjetchs->bDiscriminator("trackCountingHighPurBJetTags");
	 double TCHEpf = i_pfjetchs->bDiscriminator("pfTrackCountingHighEffBJetTags");
	 double TCHPpf = i_pfjetchs->bDiscriminator("pfTrackCountingHighPurBJetTags");
	 
	 double SoftMuonTagByIP = i_pfjetchs->bDiscriminator("softPFMuonByIP3dBJetTags");
	 double SoftElectronTagByIP = i_pfjetchs->bDiscriminator("softPFElectronByIP3dBJetTags");
	 double SoftMuonTag = i_pfjetchs->bDiscriminator("softPFMuonBJetTags");
	 double SoftElectronTag = i_pfjetchs->bDiscriminator("softPFElectronBJetTags");
	 
	 double SimpleSecVertexHE = i_pfjetchs->bDiscriminator("simpleSecondaryVertexHighEffBJetTags");
	 double SimpleSecVertexHP = i_pfjetchs->bDiscriminator("simpleSecondaryVertexHighPurBJetTags");
	 double SimpleSecVertexHEpf = i_pfjetchs->bDiscriminator("pfSimpleSecondaryVertexHighEffBJetTags");
	 double SimpleSecVertexHPpf = i_pfjetchs->bDiscriminator("pfSimpleSecondaryVertexHighPurBJetTags");
	 
	 double CSV = i_pfjetchs->bDiscriminator("combinedSecondaryVertexBJetTags");
	 double CSVpf = i_pfjetchs->bDiscriminator("pfCombinedSecondaryVertexBJetTags");
	 double CinclSVpf = i_pfjetchs->bDiscriminator("pfCombinedInclusiveSecondaryVertexBJetTags");
	 double CMVApf = i_pfjetchs->bDiscriminator("pfCombinedMVABJetTags");
	 double CSVSoftLeptonpf = i_pfjetchs->bDiscriminator("pfCombinedSecondaryVertexSoftLeptonBJetTags");
	 */
	 double CSVpfPositive = i_pfjetchs->bDiscriminator("pfPositiveCombinedSecondaryVertexV2BJetTags");
	 double CSVpfNegative = i_pfjetchs->bDiscriminator("pfNegativeCombinedSecondaryVertexV2BJetTags");

     double pfBoostedDoubleSecondaryVertex = i_pfjetchs->bDiscriminator("pfBoostedDoubleSecondaryVertexAK8BJetTags");
     //C taggers
     //
     double pfCombinedCvsL = i_pfjetchs->bDiscriminator("pfCombinedCvsLJetTags");
     double pfCombinedCvsB = i_pfjetchs->bDiscriminator("pfCombinedCvsBJetTags");


	 //the three recommended                                                                                                                                        
	 
	 double pfJetProbabilityBJetTags=i_pfjetchs->bDiscriminator("pfJetProbabilityBJetTags");
	 double pfCombinedInclusiveSecondaryVertexV2BJetTags= i_pfjetchs->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
	 double pfCombinedMVAV2BJetTags=i_pfjetchs->bDiscriminator("pfCombinedMVAV2BJetTags");

	 float partonFlavour=-100;
	 float hadronFlavour=-100;
	 
	 if (mIsMCarlo && mUseGenInfo) {
	   partonFlavour = i_pfjetchs->partonFlavour();
	   hadronFlavour = i_pfjetchs->hadronFlavour();
	 }
	 
	 qcdpfjetchs.setFlavour(partonFlavour,hadronFlavour);

	 float QGTagger=-100;

	 if(mAK4){
	   QGTagger = i_pfjetchs->userFloat("QGTaggerAK4PFCHS:qgLikelihood");
	 }
     


	 //Filling B-tag infos
	 //qcdpfjetchs.setTCHETag(TCHE,TCHP,TCHEpf,TCHPpf);
	 //qcdpfjetchs.setSoftLeptonTag(SoftMuonTagByIP,SoftElectronTagByIP,SoftMuonTag,SoftElectronTag);
	 //qcdpfjetchs.setSimpleSecondaryVertexTag(SimpleSecVertexHE,SimpleSecVertexHP,SimpleSecVertexHEpf,SimpleSecVertexHPpf);
	 //qcdpfjetchs.setCombinedSecondaryVertexTag(CSV,CSVpf,CinclSVpf,CSVSoftLeptonpf,CMVApf);
	 qcdpfjetchs.setPositiveNegativeCSV(CSVpfPositive,CSVpfNegative);
	 qcdpfjetchs.setTagRecommended(pfJetProbabilityBJetTags,pfCombinedInclusiveSecondaryVertexV2BJetTags,pfCombinedMVAV2BJetTags);
	 qcdpfjetchs.setQGTagger(QGTagger);	 
     qcdpfjetchs.setBoosted(pfBoostedDoubleSecondaryVertex);
     qcdpfjetchs.setCTagger(pfCombinedCvsL,pfCombinedCvsB);

	 float pileupJetId = -999;
	 if ( i_pfjetchs->hasUserFloat(pfchsjetpuid) )   {  pileupJetId = i_pfjetchs->userFloat(pfchsjetpuid);}
	 qcdpfjetchs.SetPUJetId(pileupJetId);
	 
	 if (mIsMCarlo) {
	   GenJetCollection::const_iterator i_matchedchs;
	   float rmin(999);
	   for(GenJetCollection::const_iterator i_gen = genjets->begin(); i_gen != genjets->end(); i_gen++) {
	     double deltaR = reco::deltaR(*i_pfjetchs,*i_gen);
	     if (deltaR < rmin) {
	       rmin = deltaR;
	       i_matchedchs = i_gen;
	     }
	   }
	   if (genjets->size() == 0) {
	     LorentzVector tmpP4(0.0,0.0,0.0,0.0);
	     qcdpfjetchs.setGen(tmpP4,0);
	   }
	   else
	     qcdpfjetchs.setGen(i_matchedchs->p4(),rmin);
	 } // if (mIsMCarlo)
	 else {
	   LorentzVector tmpP4(0.0,0.0,0.0,0.0);
	   qcdpfjetchs.setGen(tmpP4,0);
	 }
	 if (qcdpfjetchs.pt() >= mMinPFPt)
	   mPFJetsCHS.push_back(qcdpfjetchs);
	 

       } // if(i_pfjetchs->isPFJet() )
     } // for(edm::View<pat::Jet>::const_iterator i_pfjetchs=patjetschs->begin(); i_pfjetchs!=patjetschs->end(); ++i_pfjetchs)
  
  
  //---------------- met ---------------------------------------------
  Handle<pat::METCollection> pfmet;
  event.getByToken(mPFMET, pfmet);
  const pat::MET &met = pfmet->front();
  mPFMet.setVar(met.et(),met.sumEt(),met.phi());
  

  //-------------- fill the tree -------------------------------------
  sort(mPFJets.begin(),mPFJets.end(),sort_pfjets);
  mEvent->setEvtHdr(mEvtHdr);
  //mEvent->setPFJets(mPFJets);
  mEvent->setPFJetsCHS(mPFJetsCHS); // -- later substitute chs jets
  mEvent->setGenJets(mGenJets);
  mEvent->setGenFlavour(GenFlavour);
  mEvent->setGenHadronFlavour(GenHadronFlavour);
  mEvent->setPFMET(mPFMet);
  mEvent->setL1Obj(mL1Objects);
  mEvent->setHLTObj(mHLTObjects);
  if ((mEvent->nPFJetsCHS() >= (unsigned)mMinNPFJets) ) {
    if ((mEvent->pfchsmjjcor(0) >= mMinJJMass) ) {
      //    cout<<"Feeling tree ----"<<endl;
      mTree->Fill();
    }
  }
  //if (mPFPayloadName != "") {
  //delete mPFUnc;
  //delete mPFUncSrc;
  //}

}


/////////////// Matching Flavour ///////////////////////////////

int ProcessedTreeProducerBTag::getMatchedPartonGen(edm::Event const& event,GenJetCollection::const_iterator i_gen)
{

  int jetFlavour=-100;
  bool switchB=0;
  bool switchC=0;

  edm::Handle<reco::GenParticleCollection> genParticles;
  event.getByToken(mgenParticles, genParticles);

  for (size_t i = 0; i < genParticles->size (); ++i) {
      const GenParticle & genIt = (*genParticles)[i];
      int pdgId = genIt.pdgId();
      double DeltaR=deltaR(genIt.p4().eta(),genIt.p4().phi(),i_gen->eta(),i_gen->phi());
      double DeltaRmin=0.3;
      if (DeltaR < DeltaRmin ){

	DeltaRmin=DeltaR;
	if(abs(pdgId)==5){ jetFlavour=5; switchB=true;}
	if(abs(pdgId)==4){ jetFlavour=4; switchC=true;}
	if(abs(pdgId)<=3 && abs(pdgId)>=1){ jetFlavour=1; }
	if(abs(pdgId)==21){ jetFlavour=21; }
      }
      
      if (switchB) {jetFlavour=5;}
      if (switchC && !switchB) {jetFlavour=4;}

  }

  return jetFlavour;

}

int ProcessedTreeProducerBTag::getMatchedHadronGen(edm::Event const& event,GenJetCollection::const_iterator i_gen)
{

  int jetFlavour=-100;

  edm::Handle<reco::GenParticleCollection> genParticles;
  event.getByToken(mgenParticles, genParticles);

  for (size_t i = 0; i < genParticles->size (); ++i) {
    const GenParticle & genIt = (*genParticles)[i];
    
    int aid = abs(genIt.pdgId());
    if (aid/100 == 5 || aid/1000==5) {
      // 2J+1 == 1 (mesons) or 2 (baryons)
      if (aid%10 == 1 || aid%10 == 2) {
	// No B decaying to B
	if (aid != 5222 && aid != 5112 && aid != 5212 && aid != 5322) {
	  double DeltaR=deltaR(genIt.p4().eta(),genIt.p4().phi(),i_gen->eta(),i_gen->phi());
	  if(sqrt(DeltaR)<0.5){
	    jetFlavour=5;
	  }
	  else jetFlavour=21;
	}
      }
    }
  }

  return jetFlavour;
}


//////////////////////////////////////////////////////////////////////////////////////////
ProcessedTreeProducerBTag::~ProcessedTreeProducerBTag()
{

}

DEFINE_FWK_MODULE(ProcessedTreeProducerBTag);
