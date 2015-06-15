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

#include "SMPJ/AnalysisFW/plugins/ProcessedTreeProducerPFCalo.h"
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
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/JetExtendedAssociation.h"
#include "DataFormats/JetReco/interface/JetID.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/HcalNoiseSummary.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"

ProcessedTreeProducerPFCalo::ProcessedTreeProducerPFCalo(edm::ParameterSet const& cfg)
{
  mPFJECservice      = cfg.getParameter<std::string>               ("pfjecService");
  mCaloJECservice    = cfg.getParameter<std::string>               ("calojecService");
  mPFPayloadName     = cfg.getParameter<std::string>               ("PFPayloadName");
  mCaloPayloadName   = cfg.getParameter<std::string>               ("CaloPayloadName");
  mGoodVtxNdof       = cfg.getParameter<double>                    ("goodVtxNdof");
  mGoodVtxZ          = cfg.getParameter<double>                    ("goodVtxZ");
  mMinCaloPt         = cfg.getParameter<double>                    ("minCaloPt");
  mMinPFPt           = cfg.getParameter<double>                    ("minPFPt");
  mMinJJMass         = cfg.getParameter<double>                    ("minJJMass");
  mMaxY              = cfg.getParameter<double>                    ("maxY");
  mMinNCaloJets      = cfg.getParameter<int>                       ("minNCaloJets");
  mMinNPFJets        = cfg.getParameter<int>                       ("minNPFJets");
  mCaloJetID         = cfg.getParameter<edm::InputTag>             ("calojetID");
  mCaloJetExtender   = cfg.getParameter<edm::InputTag>             ("calojetExtender");
  mOfflineVertices   = cfg.getParameter<edm::InputTag>             ("offlineVertices");
  mPFJetsName        = cfg.getParameter<edm::InputTag>             ("pfjets");
  mCaloJetsName      = cfg.getParameter<edm::InputTag>             ("calojets");
  mSrcCaloRho        = cfg.getParameter<edm::InputTag>             ("srcCaloRho");
  mSrcPFRho          = cfg.getParameter<edm::InputTag>             ("srcPFRho");
  mSrcPU             = cfg.getUntrackedParameter<edm::InputTag>    ("srcPU",edm::InputTag("addPileupInfo"));
  mGenJetsName       = cfg.getUntrackedParameter<edm::InputTag>    ("genjets",edm::InputTag(""));
  mPrintTriggerMenu  = cfg.getUntrackedParameter<bool>             ("printTriggerMenu",false);
  mIsMCarlo          = cfg.getUntrackedParameter<bool>             ("isMCarlo",false);
  mUseGenInfo        = cfg.getUntrackedParameter<bool>             ("useGenInfo",false);
  mMinGenPt          = cfg.getUntrackedParameter<double>           ("minGenPt",30);
  processName_       = cfg.getParameter<std::string>               ("processName");
  triggerNames_      = cfg.getParameter<std::vector<std::string> > ("triggerName");
  triggerResultsTag_ = cfg.getParameter<edm::InputTag>             ("triggerResults");
  triggerEventTag_   = cfg.getParameter<edm::InputTag>             ("triggerEvent");
  mPFJECUncSrc       = cfg.getParameter<std::string>               ("jecUncSrc");
  mPFJECUncSrcNames  = cfg.getParameter<std::vector<std::string> > ("jecUncSrcNames");
}
//////////////////////////////////////////////////////////////////////////////////////////
void ProcessedTreeProducerPFCalo::beginJob()
{
  mTree = fs->make<TTree>("ProcessedTree","ProcessedTree");
  mEvent = new QCDEvent();
  mTree->Branch("events","QCDEvent",&mEvent);
  mTriggerNamesHisto = fs->make<TH1F>("TriggerNames","TriggerNames",1,0,1);
  mTriggerNamesHisto->SetBit(TH1::kCanRebin);
  for(unsigned i=0;i<triggerNames_.size();i++)
    mTriggerNamesHisto->Fill(triggerNames_[i].c_str(),1);
  mTriggerPassHisto = fs->make<TH1F>("TriggerPass","TriggerPass",1,0,1);
  mTriggerPassHisto->SetBit(TH1::kCanRebin);
  isPFJecUncSet_ = false;
  isCaloJecUncSet_ = false;
}
//////////////////////////////////////////////////////////////////////////////////////////
void ProcessedTreeProducerPFCalo::endJob()
{
}
//////////////////////////////////////////////////////////////////////////////////////////
void ProcessedTreeProducerPFCalo::beginRun(edm::Run const & iRun, edm::EventSetup const& iSetup)
{
  bool changed(true);
  if (hltConfig_.init(iRun,iSetup,processName_,changed)) {
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
    cout << "ProcessedTreeProducerPFCalo::analyze:"
         << " config extraction failure with process name "
         << processName_ << endl;
  }
}
//////////////////////////////////////////////////////////////////////////////////////////
void ProcessedTreeProducerPFCalo::analyze(edm::Event const& event, edm::EventSetup const& iSetup)
{
  vector<QCDCaloJet>    mCaloJets;
  vector<QCDPFJet>      mPFJets;
  vector<QCDPFJet>      tmpPFJets;
  vector<LorentzVector> mGenJets;
  QCDEventHdr mEvtHdr;
  QCDMET mCaloMet,mPFMet;
  //-------------- Basic Event Info ------------------------------
  mEvtHdr.setRun(event.id().run());
  mEvtHdr.setEvt(event.id().event());
  mEvtHdr.setLumi(event.luminosityBlock());
  mEvtHdr.setBunch(event.bunchCrossing());
  //-------------- Beam Spot --------------------------------------
  Handle<reco::BeamSpot> beamSpot;
  event.getByLabel("offlineBeamSpot", beamSpot);
  if (beamSpot.isValid())
    mEvtHdr.setBS(beamSpot->x0(),beamSpot->y0(),beamSpot->z0());
  else
    mEvtHdr.setBS(-999,-999,-999);

  //-------------- HCAL Noise Summary -----------------------------
  Handle<bool> noiseSummary;
  if (!mIsMCarlo) {
    event.getByLabel(edm::InputTag("HBHENoiseFilterResultProducer","HBHENoiseFilterResult"), noiseSummary);
    mEvtHdr.setHCALNoise(*noiseSummary);
  }
  else
   mEvtHdr.setHCALNoise(true);
  //-------------- Trigger Info -----------------------------------
  event.getByLabel(triggerResultsTag_,triggerResultsHandle_);
  if (!triggerResultsHandle_.isValid()) {
    cout << "ProcessedTreeProducerPFCalo::analyze: Error in getting TriggerResults product from Event!" << endl;
    return;
  }
  event.getByLabel(triggerEventTag_,triggerEventHandle_);
  if (!triggerEventHandle_.isValid()) {
    cout << "ProcessedTreeProducerPFCalo::analyze: Error in getting TriggerEvent product from Event!" << endl;
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
      const std::pair<int,int> prescales(hltConfig_.prescaleValues(event,iSetup,triggerNames_[itrig]));
      preL1    = prescales.first;
      preHLT   = prescales.second;
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
  event.getByLabel(mOfflineVertices,recVtxs);
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
  event.getByLabel(mSrcCaloRho,rhoCalo);
  Handle<double> rhoPF;
  event.getByLabel(mSrcPFRho,rhoPF);
  mEvtHdr.setRho(*rhoCalo,*rhoPF);
  //-------------- Generator Info -------------------------------------
  Handle<GenEventInfoProduct> hEventInfo;
  //-------------- Simulated PU Info ----------------------------------
  Handle<std::vector<PileupSummaryInfo> > PupInfo;
  if (mIsMCarlo && mUseGenInfo) {
    event.getByLabel("generator", hEventInfo);
    mEvtHdr.setPthat(hEventInfo->binningValues()[0]);
    mEvtHdr.setWeight(hEventInfo->weight());
//    event.getByLabel(mSrcPU, PupInfo);
//    std::vector<PileupSummaryInfo>::const_iterator PUI;
//    int nbx = PupInfo->size();
//    int ootpuEarly(0),ootpuLate(0),intpu(0);
//    float Tnpv = -1.; // new variable for computing pileup weight factor for the event
//    for(PUI = PupInfo->begin(); PUI != PupInfo->end(); ++PUI) {
//      if (PUI->getBunchCrossing() < 0)
//        ootpuEarly += PUI->getPU_NumInteractions();
//      else if (PUI->getBunchCrossing() > 0)
//        ootpuLate += PUI->getPU_NumInteractions();
//      else {
//        intpu += PUI->getPU_NumInteractions();
//        Tnpv = PUI->getTrueNumInteractions();
//       }
//    }
//
//    mEvtHdr.setPU(nbx,ootpuEarly,ootpuLate,intpu);
//    mEvtHdr.setTrPu(Tnpv);
  }
  else {
    mEvtHdr.setPthat(0);
    mEvtHdr.setWeight(0);
//    mEvtHdr.setPU(0,0,0,0);
//    mEvtHdr.setTrPu(0);
  }
  //---------------- Jets ---------------------------------------------
  mPFJEC   = JetCorrector::getJetCorrector(mPFJECservice,iSetup);
  mCALOJEC = JetCorrector::getJetCorrector(mCaloJECservice,iSetup);
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
      }
    }
    isPFJecUncSet_ = true;
  }
  edm::ESHandle<JetCorrectorParametersCollection> CaloJetCorParColl;
  if (mCaloPayloadName != "" && !isCaloJecUncSet_){
    iSetup.get<JetCorrectionsRecord>().get(mCaloPayloadName,CaloJetCorParColl);
    JetCorrectorParameters const& CaloJetCorPar = (*CaloJetCorParColl)["Uncertainty"];
    mCALOUnc = new JetCorrectionUncertainty(CaloJetCorPar);
    isCaloJecUncSet_ = true;
  }
  Handle<GenJetCollection>  genjets;
  Handle<PFJetCollection>   pfjets;
  Handle<CaloJetCollection> calojets;
  Handle<JetExtendedAssociation::Container> calojetExtender;
  Handle<ValueMap<reco::JetID> > calojetID;
  event.getByLabel(mPFJetsName,pfjets);
  event.getByLabel(mCaloJetsName,calojets);
  event.getByLabel(mCaloJetExtender,calojetExtender);
  event.getByLabel(mCaloJetID,calojetID);
  if (mIsMCarlo) {
    event.getByLabel(mGenJetsName,genjets);
    for(GenJetCollection::const_iterator i_gen = genjets->begin(); i_gen != genjets->end(); i_gen++) {
      if (i_gen->pt() > mMinGenPt && fabs(i_gen->y()) < mMaxY) {
        mGenJets.push_back(i_gen->p4());
      }
    }
  }
  //----------- PFJets -------------------------
  for(PFJetCollection::const_iterator i_pfjet = pfjets->begin(); i_pfjet != pfjets->end(); i_pfjet++) {
    QCDPFJet qcdpfjet;
    //int index = i_pfjet-pfjets->begin();
    //edm::RefToBase<reco::Jet> pfjetRef(edm::Ref<PFJetCollection>(pfjets,index));
    double scale = mPFJEC->correction(*i_pfjet,event,iSetup);
    //---- preselection -----------------
    if (fabs(i_pfjet->y()) > mMaxY) continue;
    //---- vertex association -----------
    //---- get the vector of tracks -----
    reco::TrackRefVector vTrks(i_pfjet->getTrackRefs());
    float sumTrkPt(0.0),sumTrkPtBeta(0.0),sumTrkPtBetaStar(0.0),beta(0.0),betaStar(0.0);
    //---- loop over the tracks of the jet ----
    for(reco::TrackRefVector::const_iterator i_trk = vTrks.begin(); i_trk != vTrks.end(); i_trk++) {
      if (recVtxs->size() == 0) break;
      sumTrkPt += (*i_trk)->pt();
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
              }
              else {
                sumTrkPtBetaStar += (*i_trk)->pt();
              }
              break;
            }
          }
        }
      }
    }
    if (sumTrkPt > 0) {
      beta     = sumTrkPtBeta/sumTrkPt;
      betaStar = sumTrkPtBetaStar/sumTrkPt;
    }
    qcdpfjet.setBeta(beta);
    qcdpfjet.setBetaStar(betaStar);
    //---- jec uncertainty --------------
    double unc(0.0);
    vector<float> uncSrc(0);
    if (mPFPayloadName != "") {
      mPFUnc->setJetEta(i_pfjet->eta());
      mPFUnc->setJetPt(scale * i_pfjet->pt());
      unc = mPFUnc->getUncertainty(true);
    }
    if (mPFJECUncSrc != "") {
      for(unsigned isrc=0;isrc<mPFJECUncSrcNames.size();isrc++) {
        mPFUncSrc[isrc]->setJetEta(i_pfjet->eta());
        mPFUncSrc[isrc]->setJetPt(scale * i_pfjet->pt());
        float unc1 = mPFUncSrc[isrc]->getUncertainty(true);
        uncSrc.push_back(unc1);
      }
    }
    qcdpfjet.setP4(i_pfjet->p4());
    qcdpfjet.setCor(scale);
    qcdpfjet.setUnc(unc);
    qcdpfjet.setUncSrc(uncSrc);
    qcdpfjet.setArea(i_pfjet->jetArea());
    double chf   = i_pfjet->chargedHadronEnergyFraction();
    double nhf   = (i_pfjet->neutralHadronEnergy() + i_pfjet->HFHadronEnergy())/i_pfjet->energy();
    double phf   = i_pfjet->photonEnergyFraction();
    double elf   = i_pfjet->electronEnergyFraction();
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
    bool looseID  = (npr>1 && phf<0.99 && nhf<0.99 && ((fabs(i_pfjet->eta())<=2.4 && elf<0.99 && chf>0 && chm>0) || fabs(i_pfjet->eta())>2.4));
    bool tightID  = (npr>1 && phf<0.99 && nhf<0.99 && ((fabs(i_pfjet->eta())<=2.4 && nhf<0.9 && phf<0.9 && elf<0.99 && chf>0 && chm>0) || fabs(i_pfjet->eta())>2.4));
    qcdpfjet.setLooseID(looseID);
    qcdpfjet.setTightID(tightID);
    qcdpfjet.setFrac(chf,nhf,phf,elf,muf);
    qcdpfjet.setMulti(npr,chm,nhm,phm,elm,mum);
    qcdpfjet.setHFFrac(hf_hf,hf_phf);
    qcdpfjet.setHFMulti(hf_hm,hf_phm);
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
    }
    else {
      LorentzVector tmpP4(0.0,0.0,0.0,0.0);
      qcdpfjet.setGen(tmpP4,0);
    }
    if (qcdpfjet.ptCor() >= mMinPFPt)
      mPFJets.push_back(qcdpfjet);
    if (  qcdpfjet.looseID())
      tmpPFJets.push_back(qcdpfjet);
  }
  //----------- CaloJets -----------------------
  for(CaloJetCollection::const_iterator i_calojet = calojets->begin(); i_calojet != calojets->end(); i_calojet++) {
    int index = i_calojet-calojets->begin();
    edm::RefToBase<reco::Jet> calojetRef(edm::Ref<CaloJetCollection>(calojets,index));
    double scale = mCALOJEC->correction(*i_calojet,event,iSetup);
    //---- preselection -----------------
    if (fabs(i_calojet->y()) > mMaxY) continue;
    double unc(0.0);
    vector<float> uncSrc(0);
    if (mCaloPayloadName != "") {
      mCALOUnc->setJetEta(i_calojet->eta());
      mCALOUnc->setJetPt(scale * i_calojet->pt());
      unc = mCALOUnc->getUncertainty(true);
    }
    QCDCaloJet qcdcalojet;
    qcdcalojet.setP4(i_calojet->p4());
    qcdcalojet.setCor(scale);
    qcdcalojet.setUnc(unc);
    qcdcalojet.setUncSrc(uncSrc);
    qcdcalojet.setArea(i_calojet->jetArea());
    double emf    = i_calojet->emEnergyFraction();
    int n90hits   = int((*calojetID)[calojetRef].n90Hits);
    double fHPD   = (*calojetID)[calojetRef].fHPD;
    double fRBX   = (*calojetID)[calojetRef].fRBX;
    int nTrkVtx   = JetExtendedAssociation::tracksAtVertexNumber(*calojetExtender,*i_calojet);
    int nTrkCalo  = JetExtendedAssociation::tracksAtCaloNumber(*calojetExtender,*i_calojet);
    bool looseID  = ((emf>0.01 || fabs(i_calojet->eta())>2.6) && (n90hits>1) && (fHPD<0.98));
    bool tightID  = ((emf>0.01 || fabs(i_calojet->eta())>2.6) && (n90hits>1) && ((fHPD<0.98 && i_calojet->pt()<=25) || (fHPD<0.95 && i_calojet->pt()>25)));
    qcdcalojet.setVar(emf,fHPD,fRBX,n90hits,nTrkCalo,nTrkVtx);
    qcdcalojet.setLooseID(looseID);
    qcdcalojet.setTightID(tightID);
    if (mIsMCarlo) {
      GenJetCollection::const_iterator i_matched;
      float rmin(999);
      for(GenJetCollection::const_iterator i_gen = genjets->begin(); i_gen != genjets->end(); i_gen++) {
        double deltaR = reco::deltaR(*i_calojet,*i_gen);
        if (deltaR < rmin) {
          rmin = deltaR;
          i_matched = i_gen;
        }
      }
      if (genjets->size() == 0) {
        LorentzVector tmpP4(0.0,0.0,0.0,0.0);
        qcdcalojet.setGen(tmpP4,0);
      }
      else
        qcdcalojet.setGen(i_matched->p4(),rmin);
    }
    else {
      LorentzVector tmpP4(0.0,0.0,0.0,0.0);
      qcdcalojet.setGen(tmpP4,0);
    }
    if (qcdcalojet.ptCor() >= mMinCaloPt)
      mCaloJets.push_back(qcdcalojet);
  }

  //---------------- met ---------------------------------------------
  Handle<PFMETCollection> pfmet;
  Handle<CaloMETCollection> calomet;
  event.getByLabel("pfMet",pfmet);
  event.getByLabel("caloMet",calomet);
  mPFMet.setVar((*pfmet)[0].et(),(*pfmet)[0].sumEt(),(*pfmet)[0].phi());
  mCaloMet.setVar((*calomet)[0].et(),(*calomet)[0].sumEt(),(*calomet)[0].phi());
  //-------------- fill the tree -------------------------------------
  sort(mCaloJets.begin(),mCaloJets.end(),sort_calojets);
  sort(mPFJets.begin(),mPFJets.end(),sort_pfjets);
  mEvent->setEvtHdr(mEvtHdr);
  mEvent->setCaloJets(mCaloJets);
  mEvent->setPFJets(mPFJets);
  mEvent->setGenJets(mGenJets);
  mEvent->setCaloMET(mCaloMet);
  mEvent->setPFMET(mPFMet);
  mEvent->setL1Obj(mL1Objects);
  mEvent->setHLTObj(mHLTObjects);
  if ((mEvent->nPFJets() >= (unsigned)mMinNPFJets) && (mEvent->nCaloJets() >= (unsigned)mMinNCaloJets)) {
    if ((mEvent->pfmjjcor(0) >= mMinJJMass) || (mEvent->calomjjcor(0) >= mMinJJMass) ) {
      mTree->Fill();
    }
  }
  //if (mPFPayloadName != "") {
    //delete mPFUnc;
    //delete mPFUncSrc;
 //}
  //if (mCaloPayloadName != "")
    //delete mCALOUnc;
}
//////////////////////////////////////////////////////////////////////////////////////////
ProcessedTreeProducerPFCalo::~ProcessedTreeProducerPFCalo()
{
}

DEFINE_FWK_MODULE(ProcessedTreeProducerPFCalo);
