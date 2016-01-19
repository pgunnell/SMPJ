//Author K. Kousouris
//Modified by: S. Ganguly

#ifndef QCDEvent_h
#define QCDEvent_h
#include "SMPJ/AnalysisFW/interface/QCDJet.h"
#include "SMPJ/AnalysisFW/interface/QCDMET.h"
#include "SMPJ/AnalysisFW/interface/QCDCaloJet.h"
#include "SMPJ/AnalysisFW/interface/QCDPFJet.h"
#include "SMPJ/AnalysisFW/interface/QCDEventHdr.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include <vector>

class QCDEvent
{
    public:
      typedef reco::Particle::LorentzVector LorentzVector;
      //------------ Constructor ------------------------------
      QCDEvent();
      //------------ Destructor -------------------------------
      ~QCDEvent();
      //------------ Set methods ------------------------------
      void setCaloMET(const QCDMET& fCaloMET)                     {CaloMet_ = fCaloMET;}
      void setPFMET(const QCDMET& fPFMET)                         {PFMet_ = fPFMET;}
      void setEvtHdr(const QCDEventHdr& fEvtHdr)                  {EvtHdr_ = fEvtHdr;}
      void setCaloJets(const std::vector<QCDCaloJet>& fCaloJets);
      void setPFJets(const std::vector<QCDPFJet>& fPFJets);
      void setPFJetsCHS(const std::vector<QCDPFJet>& fPFJetsCHS);
      //void setFatJets(const std::vector<QCDJet>& fFatJets);
      void setGenJets(const std::vector<LorentzVector>& fGenJets);
      void setL1Obj(const std::vector<std::vector<LorentzVector> >& fL1Obj);
      void setHLTObj(const std::vector<std::vector<LorentzVector> >& fHLTObj);
      void setFilterId(const std::vector<std::vector<int> >& filterIdList) {filterIdList_ = filterIdList;}
      void setPrescales(const std::vector<int>& fPreL1, const std::vector<int>& fPreHLT) {L1Prescale_ = fPreL1; HLTPrescale_ = fPreHLT;}
      void setTrigDecision(const std::vector<int>& fTrigDecision) {TriggerDecision_ = fTrigDecision;}
      void setTrigPathList(const std::vector<std::string>& trigPathList) {triggerList_ = trigPathList;}

      //------------ Get methods -------------------------------
      unsigned int nTriggers()                         const {return TriggerDecision_.size();}
      unsigned int nL1Obj(int i)                       const {return L1Obj_[i].size();}
      unsigned int nHLTObj(int i)                      const {return HLTObj_[i].size();}
      unsigned int nPFJets()                           const {return PFJets_.size();}
      unsigned int nPFJetsCHS()                        const {return PFJetsCHS_.size();}
      //unsigned int nFatJets()                          const {return FatJets_.size();}
      unsigned int nCaloJets()                         const {return CaloJets_.size();}
      unsigned int nGenJets()                          const {return GenJets_.size();}
      int nGoodJets(int unc, int id, float ymax, float ptmin, std::vector<QCDJet> jets);
      int fired(int i)                                 const {return TriggerDecision_[i];}
      int preL1(int i)                                 const {return L1Prescale_[i];}
      int preHLT(int i)                                const {return HLTPrescale_[i];}
      float pfmjj();

      float calomjj();
      float genmjj();
      float pfchsmjjcor(int unc);
      float pfchsmjjcor(int unc,int src);
      float pfmjjcor(int unc);
      float pfmjjcor(int unc,int src);
      //float fatmjjcor(int unc);
      float calomjjcor(int unc);
      float pfmjjgen();
      float calomjjgen();
      const QCDMET&        calomet()                   const {return CaloMet_;}
      const QCDMET&        pfmet()                     const {return PFMet_;}
      const LorentzVector& hltobj(int itrig, int iobj) const {return (HLTObj_[itrig])[iobj];}
      const LorentzVector& l1obj(int itrig, int iobj)  const {return (L1Obj_[itrig])[iobj];}
      const LorentzVector& genjet(int i)               const {return GenJets_[i];}
      const QCDPFJet&      pfjet(int i)                const {return PFJets_[i];}
      const QCDPFJet&      pfjetchs(int i)             const {return PFJetsCHS_[i];}
      //const QCDJet&        fatjet(int i)               const {return FatJets_[i];}
      const QCDCaloJet&    calojet(int i)              const {return CaloJets_[i];}
      const QCDEventHdr&   evtHdr()                    const {return EvtHdr_;}

      const std::vector<std::vector<LorentzVector>>& HLTObj() const {return HLTObj_;}
      const std::vector<std::string>& trigPathList() const {return triggerList_;}
      const std::vector<std::vector<int>>& filterIdList() const {return filterIdList_;}
      const std::vector<LorentzVector>& hltObjsForPath(int i) const {return HLTObj_[i];}
      const std::vector<int>& filterIdsForPath(int i) const {return filterIdList_[i];}

    private:
      std::vector<std::vector<int> > filterIdList_;

      //---- event header (contains all the event info) --------------
      QCDEventHdr                              EvtHdr_;
      //---- CALO met object -----------------------------------------
      QCDMET                                   CaloMet_;
      //---- PF met object -------------------------------------------
      QCDMET                                   PFMet_;
      //---- trigger decision vector ---------------------------------
      std::vector<int>                         TriggerDecision_;
      std::vector<std::string>                 triggerList_;
      //---- L1 prescale vector --------------------------------------
      std::vector<int>                         L1Prescale_;
      //---- HLT prescale vector -------------------------------------
      std::vector<int>                         HLTPrescale_;
      //---- HLT objects ---------------------------------------------
      std::vector<std::vector<LorentzVector> > HLTObj_;
      //---- L1 objects ----------------------------------------------
      std::vector<std::vector<LorentzVector> > L1Obj_;
      //---- Genjets -------------------------------------------------
      std::vector<LorentzVector>               GenJets_;
      //---- CaloJets ------------------------------------------------
      std::vector<QCDCaloJet>                  CaloJets_;
      //---- PFJets --------------------------------------------------
      std::vector<QCDPFJet>                    PFJets_;
      //---- PFJetsCHS -----------------------------------------------
      std::vector<QCDPFJet>                    PFJetsCHS_;
      //---- FatJets -------------------------------------------------
      //std::vector<QCDJet>                      FatJets_;
};
#endif
