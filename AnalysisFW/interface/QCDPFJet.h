//Author K. Kousouris
//Modified by: S. Ganguly

#ifndef QCDPFJet_h
#define QCDPFJet_h
#include "SMPJ/AnalysisFW/interface/QCDJet.h"
#include "SMPJ/AnalysisFW/interface/QCDPFJetBTag.h"
#include "TLorentzVector.h"
class QCDPFJet : public QCDPFJetBTag {
   public:
     //------------ Constructor ------------------------------
     QCDPFJet() {chf_=0;nhf_=0;nemf_=0;cemf_=0;muf_=0;chm_=0;nhm_=0;phm_=0;elm_=0;mum_=0,cm_=0;
//     pfParticles_.clear();
     }
     //------------ Destructor -------------------------------
     ~QCDPFJet() {}
     //------------ Set methods ------------------------------
     void setFrac(float fchf, float fnhf, float fnemf, float fcemf, float fmuf)  {chf_ = fchf; nhf_ = fnhf; nemf_ = fnemf; cemf_ = fcemf; muf_ = fmuf;}
     void setMulti(int fncand, int fchm, int fnhm, int fphm, int felm, int fmum, int fcm) {ncand_ = fncand; chm_ = fchm; nhm_ = fnhm; phm_ = fphm; elm_ = felm; mum_ = fmum; cm_ = fcm; }
     void setBeta(float fbeta) {beta_ = fbeta;}
     void setBetaStar(float fbetaStar) {betaStar_ = fbetaStar;}
     void setHFFrac(float fhf_hf, float fhf_phf) {hf_hf_ = fhf_hf; hf_phf_ = fhf_phf;}
     void setHFMulti(int fhf_hm, int fhf_phm) {hf_hm_ = fhf_hm; hf_phm_ = fhf_phm;}
     void setVtxInfo(int mpuTrk, int mlvTrk, int mjtTrk) { mpuTrk_ = mpuTrk; mlvTrk_ = mlvTrk; mjtTrk_ = mjtTrk;} // Juska
     void setHO(float hof) {hof_ = hof;} // Juska
     void SetPUJetId(float pujid) { pujid_ = pujid; }
     void SetCaloJetPt(float calojetpt) { calojetpt_ = calojetpt; }
     void SetCaloJetEf(float calojetef) { calojetef_ = calojetef; }   

    /*
     void setPFParticles(std::vector<LorentzVector>& fpfFParticles) {
      for(unsigned i=0; i<fpfFParticles.size(); i++)
       {
        pfParticles_.push_back(fpfFParticles[i]);
      }
     } // setPFParticles
     */

     //------------ Get methods ------------------------------
     float beta()     const {return beta_;}
     float betaStar() const {return betaStar_;}
     float chf()      const {return chf_;}
     float nhf()      const {return nhf_;}
     float nemf()      const {return nemf_;}
     float cemf()      const {return cemf_;}
     float muf()      const {return muf_;}
     float hf_hf()    const {return hf_hf_;}
     float hf_phf()   const {return hf_phf_;}
     int chm()        const {return chm_;}
     int nhm()        const {return nhm_;}
     int phm()        const {return phm_;}
     int elm()        const {return elm_;}
     int mum()        const {return mum_;}
     int hf_hm()      const {return hf_hm_;}
     int hf_phm()     const {return hf_phm_;}
     int ncand()      const {return ncand_;}
     float CaloJetPt() const {return calojetpt_;}
     float CaloJetEf() const {return calojetef_;}
 //    int nParticles() const {return pfParticles_.size();}
 //    const LorentzVector& getPFParticles(int i) const {return pfParticles_[i];}

     int mpuTrk()     const {return mpuTrk_;} // Juska
     int mlvTrk()     const {return mlvTrk_;} //
     int mjtTrk()     const {return mjtTrk_;} //
     float hof()      const {return hof_;}    //

     float PUJetId()  const { return pujid_ ;}


   private:
     //---- charged hadron energy fraction ----
     float chf_;
     //---- neutral hadron energy fraction ----
     float nhf_;
     //---- photon energy fraction ------------
     float nemf_;
     //---- electron energy fraction ----------
     float cemf_;
     //---- muon energy fraction --------------
     float muf_;
     //-----HF Hadron Fraction ---------------
     float hf_hf_;
     //-----HF Photon Fraction ------------
     float hf_phf_;
     //-----HF Hadron Multiplicity---------
     int hf_hm_;
     //-----HF Photon Multiplicity --------
     int hf_phm_;
     //---- charged hadron multiplicity -------
     int chm_;
     //---- neutral hadron multiplicity -------
     int nhm_;
     //---- photon multiplicity ---------------
     int phm_;
     //---- electron multiplicity -------------
     int elm_;
     //---- muon multiplicity -----------------
     int mum_;
     //---- number of PF candidates -----------
     int ncand_;
     // --- charged multiplicity ------
     int cm_;
     //---- fraction of track pt coming from the signal vertex ---
     float beta_;
     //---- fraction of track pt NOT coming from the signal vertex ---
     float betaStar_;
     //----- PF Particles ----//
//     std::vector<LorentzVector> pfParticles_;

     // Juska:
     int mpuTrk_; // PU-tracks in jet
     int mlvTrk_; // lead vertex tracks in jet
     int mjtTrk_; // all tracks in jet (also unassociated)
     float hof_; // Hadronic Outer energy fraction

     float pujid_;

     float calojetpt_;
     float calojetef_;

    };
#endif
