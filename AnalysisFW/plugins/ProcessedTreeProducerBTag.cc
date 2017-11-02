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
  mSrcCaloRho(mayConsume<double>(cfg.getParameter<edm::InputTag>("srcCaloRho"))),
  mSrcPFRho(mayConsume<double>(cfg.getParameter<edm::InputTag>("srcPFRho"))),
  mPFMET(mayConsume<pat::METCollection>(cfg.getParameter<edm::InputTag>("pfmet"))),
  mGenJetsName(mayConsume<GenJetCollection>(cfg.getUntrackedParameter<edm::InputTag>("genjets",edm::InputTag("")))),
  triggerResultsTag_(mayConsume<edm::TriggerResults>(cfg.getParameter<edm::InputTag>("triggerResults"))),
  triggerEventTag_(mayConsume<trigger::TriggerEvent>(cfg.getParameter<edm::InputTag>("triggerEvent"))),
  mSrcPU(mayConsume<std::vector<PileupSummaryInfo> >(cfg.getUntrackedParameter<edm::InputTag>("srcPULabel"))),
  hltPrescale_(cfg, consumesCollector(), *this)// ",edm::InputTag("addPileupInfo"))))
{
  mPFPayloadNameCHS = cfg.getParameter<std::string> ("PFPayloadNameCHS");
  mPFJECUncSrcCHS    = cfg.getParameter<std::string>               ("jecUncSrcCHS");
  mPFJECUncSrcNames = cfg.getParameter<std::vector<std::string> > ("jecUncSrcNames");
  pfchsjetpuid    = cfg.getParameter<std::string>               ("pfchsjetpuid");
  mGoodVtxNdof       = cfg.getParameter<double>                    ("goodVtxNdof");
  mGoodVtxZ          = cfg.getParameter<double>                    ("goodVtxZ");
  mMinPFPt           = cfg.getParameter<double>                    ("minPFPt");
  mMaxY              = cfg.getParameter<double>                    ("maxY");
  mMinNPFJets        = cfg.getParameter<int>                       ("minNPFJets");
  mPrintTriggerMenu  = cfg.getUntrackedParameter<bool>             ("printTriggerMenu",false);
  mIsMCarlo          = cfg.getUntrackedParameter<bool>             ("isMCarlo",false);
  mAK4               = cfg.getUntrackedParameter<bool>             ("AK4",false);
  mUseGenInfo        = cfg.getUntrackedParameter<bool>             ("useGenInfo",false);
  mMinGenPt          = cfg.getUntrackedParameter<double>           ("minGenPt",30);
  processName_       = cfg.getParameter<std::string>               ("processName");
  triggerNames_      = cfg.getParameter<std::vector<std::string> > ("triggerName");
  mPFJetsNameCHS = consumes<edm::View<pat::Jet> >(cfg.getParameter<edm::InputTag>("pfjetschs"));
  mhEventInfo = consumes<GenEventInfoProduct>(cfg.getParameter<edm::InputTag>("EventInfo"));
  mgenParticles = consumes<reco::GenParticleCollection>(cfg.getParameter<edm::InputTag>("GenParticles"));
  qgToken = consumes<edm::ValueMap<float>>(edm::InputTag("QGTagger", "qgLikelihood"));
  jetFlavourInfosToken_ = consumes<reco::JetFlavourInfoMatchingCollection>( cfg.getParameter<edm::InputTag>("jetFlavourInfos"));
  triggerBits_ = consumes<edm::TriggerResults>(cfg.getParameter<edm::InputTag>("triggerResults"));
  triggerObjects_  = consumes<pat::TriggerObjectStandAloneCollection>(cfg.getParameter<edm::InputTag> ("triggerObjects"));
  triggerPrescales_ = consumes<pat::PackedTriggerPrescales>(cfg.getParameter<edm::InputTag>("prescales"));
  triggerPrescalesL1Max_ = consumes<pat::PackedTriggerPrescales>(cfg.getParameter<edm::InputTag>("prescalesL1Max"));
  triggerPrescalesL1Min_ = consumes<pat::PackedTriggerPrescales>(cfg.getParameter<edm::InputTag>("prescalesL1Min"));
  genEvtInfoToken = consumes<GenEventInfoProduct>(edm::InputTag("generator"));
  lheEvtInfoToken = consumes<LHEEventProduct>(edm::InputTag("externalLHEProducer"));  
  saveWeights_ = cfg.getParameter<bool>("saveWeights");
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
  vector<QCDPFJet>      mPFJetsCHS;
  vector<LorentzVector> mGenJets;
  vector<float> GenFlavour;
  vector<float> GenHadronFlavour;
  QCDEventHdr mEvtHdr;
  QCDMET mPFMet;
  vector<double> mPDFWeights;
  vector<double> mScaleWeights;

  bool save_event=false;

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

  //-------------- Trigger Info -----------------------------------
  
  edm::Handle<edm::TriggerResults> triggerBits;
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
  edm::Handle<pat::PackedTriggerPrescales> triggerPrescalesL1Max;
  edm::Handle<pat::PackedTriggerPrescales> triggerPrescalesL1Min;

  vector<int> L1Prescales,HLTPrescales,Fired;                                                                                                                          
  vector<vector<LorentzVector> > mHLTObjects;  
  vector<LorentzVector> vvL1,vvHLT;
    
  if(!mIsMCarlo){
    event.getByToken(triggerBits_, triggerBits);
    event.getByToken(triggerObjects_, triggerObjects);
    event.getByToken(triggerPrescales_, triggerPrescales);
    event.getByToken(triggerPrescalesL1Min_, triggerPrescalesL1Min);
    event.getByToken(triggerPrescalesL1Max_, triggerPrescalesL1Max);
    
    //Variables
    
    const edm::TriggerNames &names = event.triggerNames(*triggerBits);
    
    for(unsigned int k=0;k<triggerNames_.size();k++) {
      
      for(unsigned int itrig=0;itrig<triggerBits->size();itrig++) {
	//if(triggerBits->accept(itrig)) std::cout << "Trigger " << names.triggerName(itrig) <<                                                                                                         
	//  ", prescale " << triggerPrescales->getPrescaleForIndex(itrig) <<                                                                                            
	//  ": " << (triggerBits->accept(itrig) ? "PASS" : "fail (or not run)")                                                                                          
	//	  << std::endl;       
	
	string trigger_name = string(names.triggerName(itrig));
	//--- erase the last character, i.e. the version number----                                                                                                                              
	trigger_name.pop_back();
	if (trigger_name == triggerNames_[k]) {

	  HLTPrescales.push_back(triggerPrescales->getPrescaleForIndex(itrig));
	  if(triggerPrescalesL1Max->getPrescaleForIndex(itrig)>triggerPrescalesL1Min->getPrescaleForIndex(itrig))
	    L1Prescales.push_back(triggerPrescalesL1Max->getPrescaleForIndex(itrig));
	  else L1Prescales.push_back(triggerPrescalesL1Min->getPrescaleForIndex(itrig));

	  if(triggerBits->accept(itrig)) Fired.push_back(1);
	  if(!triggerBits->accept(itrig)) Fired.push_back(0);
	  
	  if(triggerBits->accept(itrig)) mTriggerPassHisto->Fill(1,1);
	}
      }
    }
    
    //std::cout << "\n === TRIGGER OBJECTS === " << std::endl;
    for (pat::TriggerObjectStandAlone obj : *triggerObjects) { // note: not "const &" since we want to call unpackPathNames
      
      obj.unpackPathNames(names);
      
      TLorentzVector P4;                                                                                                                                         
      P4.SetPtEtaPhiM(obj.pt(),obj.eta(),obj.phi(),obj.mass());                                                                                                      
      LorentzVector qcdhltobj(P4.Px(),P4.Py(),P4.Pz(),P4.E());                                                                                                   
      vvHLT.push_back(qcdhltobj);
      
      std::vector<std::string> pathNamesAll  = obj.pathNames(false);
      std::vector<std::string> pathNamesLast = obj.pathNames(true);
      
      mHLTObjects.push_back(vvHLT); 
    }
  
    mEvent->setTrigDecision(Fired);
    mEvent->setPrescales(L1Prescales,HLTPrescales);
    mEvent->setHLTObj(mHLTObjects);
  }

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
  
    event.getByToken(genEvtInfoToken,genEvtInfo);
    event.getByToken(lheEvtInfoToken,lheEvtInfo);

    float lheOriginalXWGTUP_ = lheEvtInfo->originalXWGTUP();
    
    if (saveWeights_) {
      for(unsigned i=0;i<lheEvtInfo->weights().size();i++) {
	string wtid(lheEvtInfo->weights()[i].id);
	float wgt(lheEvtInfo->weights()[i].wgt);
	if (wtid == "1002" || wtid == "2") mScaleWeights.push_back(wgt/lheOriginalXWGTUP_);
	if (wtid == "1003" || wtid == "3") mScaleWeights.push_back(wgt/lheOriginalXWGTUP_);
	if (wtid == "1004" || wtid == "4") mScaleWeights.push_back(wgt/lheOriginalXWGTUP_);
	if (wtid == "1005" || wtid == "5") mScaleWeights.push_back(wgt/lheOriginalXWGTUP_);
	if (wtid == "1007" || wtid == "7") mScaleWeights.push_back(wgt/lheOriginalXWGTUP_);
	if (wtid == "1009" || wtid == "9") mScaleWeights.push_back(wgt/lheOriginalXWGTUP_); 
	
	if ((stoi(wtid) > 2000 && stoi(wtid) <= 2102) || (stoi(wtid) > 10 && stoi(wtid) <= 110)) {
	  mPDFWeights.push_back(wgt/lheOriginalXWGTUP_);
	}
      }
    }

    mEvtHdr.setScaleWeight(mScaleWeights);
    mEvtHdr.setPDFWeight(mPDFWeights);
  }
  
  else {
    mEvtHdr.setPthat(0);
    mEvtHdr.setWeight(0);
    mEvtHdr.setPU(0,0,0,0);
    mEvtHdr.setTrPu(0);
  }

  //---------------- Jets ---------------------------------------------

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

      }
    }

    edm::Handle<reco::JetFlavourInfoMatchingCollection> theJetFlavourInfos;
    event.getByToken(jetFlavourInfosToken_, theJetFlavourInfos );
    
    for ( reco::JetFlavourInfoMatchingCollection::const_iterator j  = theJetFlavourInfos->begin();j != theJetFlavourInfos->end();++j ) {
      reco::JetFlavourInfo aInfo = (*j).second;
      int FlavourGenHadron = aInfo.getHadronFlavour();
      //if(FlavourGenHadron==5) cout<<FlavourGenHadron<<" "<<aJet->pt()<<" "<<aJet->eta()<<" "<<aJet->phi()<<" HADRONFLAV"<<endl;
      GenHadronFlavour.push_back(FlavourGenHadron);
    }
  }
  
  //----------- PFJets  CHS part -------------------------

  //uncertainties + jet pT corrected or not (otherwise it is done)

  edm::Handle<edm::ValueMap<float>> qgHandle; 
  event.getByToken(qgToken, qgHandle);
 
  edm::Handle<edm::View<pat::Jet> > patjetschs;
  event.getByToken(mPFJetsNameCHS,patjetschs);

  edm::ESHandle<JetCorrectorParametersCollection> PFJetCorParCollCHS;

  //mPFUncCHS = new JetCorrectionUncertainty(PFJetCorParCHS);//"Summer16_23Sep2016V4_MC_Uncertainty_AK8PFchs.txt");     
  
  for(edm::View<pat::Jet>::const_iterator i_pfjetchs=patjetschs->begin(); i_pfjetchs!=patjetschs->end(); ++i_pfjetchs)
    {
      QCDPFJet qcdpfjetchs;
      
      if(i_pfjetchs->isPFJet() ){
	
	double scaleCHS = 1./i_pfjetchs->jecFactor(0); // --- the value of the JEC factor
	
	//---- preselection -----------------
	if (fabs(i_pfjetchs->y()) > mMaxY) continue;
	if (fabs(i_pfjetchs->pt()) < mMinPFPt) continue;
	
	save_event=true;

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
	  iSetup.get<JetCorrectionsRecord>().get(mPFPayloadNameCHS,PFJetCorParCollCHS);
	  JetCorrectorParameters const& PFJetCorParCHS = (*PFJetCorParCollCHS)["Uncertainty"];
	  
	  mPFUncCHS = new JetCorrectionUncertainty(PFJetCorParCHS);//"Summer16_23Sep2016V4_MC_Uncertainty_AK8PFchs.txt");     
	  
	  mPFUncCHS->setJetEta(i_pfjetchs->eta());
	  mPFUncCHS->setJetPt(i_pfjetchs->pt());
	  uncCHS = mPFUncCHS->getUncertainty(true);
	  cout<<uncCHS<<endl;
	} // if (mPFPayloadName != "")
	/*if (mPFJECUncSrcCHS != "") {
	  for(unsigned isrc=0;isrc<mPFJECUncSrcNames.size();isrc++) {
	    
	    JetCorrectorParameters const& PFJetCorParCHS = (*PFJetCorParCollCHS)["Uncertainty"];
	    mPFUncSrcCHS[isrc] = new JetCorrectionUncertainty(PFJetCorParCHS);//"Summer16_23Sep2016V4_MC_Uncertainty_AK8PFchs.txt");     
	    
	    mPFUncSrcCHS[isrc]->setJetEta(i_pfjetchs->eta());
	    mPFUncSrcCHS[isrc]->setJetPt(i_pfjetchs->pt());
	    float unc1 = mPFUncSrcCHS[isrc]->getUncertainty(true);
	    uncSrcCHS.push_back(unc1);
	  } // for(unsigned isrc=0;isrc<mPFJECUncSrcNames.size();isrc++)
	} // if (mPFJECUncSrc != "")
	*/
	
	qcdpfjetchs.setP4(i_pfjetchs->p4());
	qcdpfjetchs.setCor(scaleCHS);
	qcdpfjetchs.setUnc(uncCHS);
	qcdpfjetchs.setUncSrc(uncSrcCHS);
	qcdpfjetchs.setArea(i_pfjetchs->jetArea());

	double chf   = i_pfjetchs->chargedHadronEnergyFraction();
	double nhf   = i_pfjetchs->neutralHadronEnergyFraction();// + i_pfjetchs->HFHadronEnergyFraction();
	double nemf   = i_pfjetchs->neutralEmEnergyFraction(); // equals to deprecated phf but has HF info too
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
	// https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID
	float eta    = i_pfjetchs->eta();
	int cm       = i_pfjetchs->chargedMultiplicity();
	bool looseID = (nhf<0.99 && nemf<0.99 && npr>1 && muf<0.8) && ((fabs(eta) <= 2.4 && chf>0 && cm>0 && cemf<0.99) || fabs(eta)>2.4);
	bool tightID = (nhf<0.90 && nemf<0.90 && npr>1 && muf<0.8) && ((fabs(eta)<=2.4 && chf>0 && cm>0 && cemf<0.90) || fabs(eta)>2.4);
	 
	qcdpfjetchs.setLooseID(looseID);
	qcdpfjetchs.setTightID(tightID);
	qcdpfjetchs.setFrac(chf,nhf,nemf,cemf,muf);
	qcdpfjetchs.setMulti(npr,chm,nhm,phm,elm,mum);
	qcdpfjetchs.setHFFrac(hf_hf,hf_phf);
	qcdpfjetchs.setHFMulti(hf_hm,hf_phm);
	 
	double hof   = i_pfjetchs->hoEnergyFraction(); // Juska
	qcdpfjetchs.setVtxInfo(mpuTrk,mlvTrk,mjtTrk);
	qcdpfjetchs.setHO(hof);
	
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
	  QGTagger = i_pfjetchs->userFloat("QGTagger:qgLikelihood");
	}
	
	qcdpfjetchs.setQGTagger(QGTagger);
	
	//Filling B-tag infos
	qcdpfjetchs.setTagRecommended(pfJetProbabilityBJetTags,pfCombinedInclusiveSecondaryVertexV2BJetTags,pfCombinedMVAV2BJetTags);
	
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
  sort(mPFJetsCHS.begin(),mPFJetsCHS.end(),sort_pfjets);
  mEvent->setEvtHdr(mEvtHdr);
  mEvent->setPFJetsCHS(mPFJetsCHS);
  if (mIsMCarlo) {
    mEvent->setGenJets(mGenJets);
    //mEvent->setGenFlavour(GenFlavour);
    //mEvent->setGenHadronFlavour(GenHadronFlavour);
  }

  mEvent->setPFMET(mPFMet);
  mEvent->setHLTObj(mHLTObjects);

  if (save_event) {
    mTree->Fill();
  }
  
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
