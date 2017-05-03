 //Author P. Gunnellini

#ifndef MyMuon_h
#define MyMuon_h
#include "DataFormats/PatCandidates/interface/Muon.h"
//-------- Generic Jet class for QCD analyses ---------------
class MyMuon
{
   public:
     typedef reco::Particle::LorentzVector LorentzVector;
     //------------ Constructor ------------------------------
     MyMuon() {}
     //------------ Destructor -------------------------------
     ~MyMuon() {}
     //------------ Sett methods -----------------------------
     void setP4(LorentzVector fP4) {P4_ = fP4;}
     void setGen(LorentzVector fP4, float fgenR) {genP4_ = fP4;genR_ = fgenR;}
     void setCor(float fCor)                     {cor_  = fCor;} 
     void setUnc(float fUnc)                     {unc_  = fUnc;} 
     void setUncSrc(std::vector<float> fUncSrc)  {uncSrc_ = fUncSrc;}
     //New Instructions
     void setMuonDxyVertex(double fDxyVertex)    {DxyVertex_ = fDxyVertex;}
     void setMuonDzVertex(double fDzVertex)    {DzVertex_ = fDzVertex;}
     void setPDGId(double fPDGID) {PDGID_ = fPDGID;}
     void setPfIso(double fPfIso) {PfIso_ = fPfIso;}

     void setChargedHadronIso(double fChargedHadronIso) {fChargedHadronIso = ChargedHadronIso_;}
     void setNeutralHadronIso(double fNeutralHadronIso) {fNeutralHadronIso = NeutralHadronIso_;}
     void setPhotonIso(double fPhotonIso) {fPhotonIso = PhotonIso_;} 
     void setPuChargedHadronIso(double fPuChargedHadronIso) {fPuChargedHadronIso = PuChargedHadronIso_;}

     //------------ Get methods ------------------------------
     const LorentzVector& p4()    const {return P4_;}
     const LorentzVector& genp4() const {return genP4_;}
     float pt()                   const {return P4_.pt()/cor_;}
     float genpt()                const {return genP4_.pt();}
     float geneta()               const {return genP4_.eta();} 
     float genR()                 const {return genR_;} 
     float ptCor()                const {return P4_.pt();}
     float e()                    const {return P4_.energy()/cor_;}
     float eCor()                 const {return P4_.energy();}
     float eta()                  const {return P4_.eta();}
     float y()                    const {return P4_.Rapidity();}
     float phi()                  const {return P4_.phi();}
     float mass()                 const {return P4_.mass();}
     float cor()                  const {return cor_;}
     float unc()                  const {return unc_;} 
     float uncSrc(int i)          const {return uncSrc_[i];}

     float Iso()                  const {return PfIso_;}
     float pdgId()                const {return PDGID_;}
     float MuonDZVertex()          const {return DzVertex_;}
     float MuonDxyVertex()          const {return DxyVertex_;}

     //   int  nParticles()            const {return pfParticles_.size();}
     //   const LorentzVector& getPFParticles(int i) const {return pfParticles_[i];}
 private:
     //------ jet 4-momentum vector------------------
     LorentzVector P4_;
     //------ matched genjet 4-momentum vector-------
     LorentzVector genP4_;
     //------ matching radius -----------------------
     float genR_;
     //------ jec factor ----------------------------
     float cor_;
     // ----- All components of JEC Factor ----------
     std::vector<double> jecLabels_;
     //------ jec uncertainty -----------------------
     float unc_;
     //------ jec uncertainty sources ---------------
     std::vector<float> uncSrc_;
     //------ jet area ------------------------------
     float area_;

     //Muon Vertices
     double DxyVertex_;
     double DzVertex_;
     double PDGID_;
     double PfIso_;

     double PuChargedHadronIso_;
     double ChargedHadronIso_;
     double NeutralHadronIso_;
     double PhotonIso_;

 };
#endif
