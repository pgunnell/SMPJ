//Author K. Kousouris
//Modified by: S. Ganguly

#ifndef QCDPFJetBTag_h
#define QCDPFJetBTag_h
#include "SMPJ/AnalysisFW/interface/QCDJet.h"
#include "TLorentzVector.h"
class QCDPFJetBTag : public QCDJet {
   public:
     //------------ Constructor ------------------------------
     QCDPFJetBTag() {}
     //------------ Destructor -------------------------------
     ~QCDPFJetBTag() {}
     //------------ Set methods ------------------------------
     /*
     void setTCHETag(float ftche, float ftchp, float ftchepf, float ftchppf)  {TCHE_ = ftche; TCHP_ = ftchp; TCHEpf_ = ftchepf; TCHPpf_ = ftchppf; }
     void setSoftLeptonTag(float fsoftmuonbyip, float fsoftelectronbyip, float fsoftmuon, float fsoftelectron)  {SoftMuonTagByIP_ = fsoftmuonbyip; SoftElectronTagByIP_ = fsoftelectronbyip; SoftMuonTag_ = fsoftmuon; SoftElectronTag_ = fsoftelectron; }
     void setSimpleSecondaryVertexTag(float fsimplesecvertexhe, float fsimplesecvertexhp, float fsimplesecvertexhepf, float fsimplesecvertexhppf)  {SimpleSecVertexHE_ = fsimplesecvertexhe; SimpleSecVertexHP_ = fsimplesecvertexhp; SimpleSecVertexHEpf_ = fsimplesecvertexhepf; SimpleSecVertexHPpf_ = fsimplesecvertexhppf; }
     void setCombinedSecondaryVertexTag(float fcsv, float fcsvpf, float fcinclsvpf, float fcsvsoftleptonpf, float fcmvapf)  {CSV_ = fcsv; CSVpf_ = fcsvpf; CinclSVpf_ = fcinclsvpf; CSVSoftLeptonpf_ = fcsvsoftleptonpf; CMVApf_= fcmvapf;}
     */
     void setPositiveNegativeCSV(float fcsvpfpositive, float fcsvpfnegative) { CSVpfPositive_ = fcsvpfpositive; CSVpfNegative_ = fcsvpfnegative;}

     void setTagRecommended(float recommend1, float recommend2, float recommend3) { recommend1_ = recommend1; recommend2_ = recommend2; recommend3_ = recommend3; } 

     void setDeepCSV (float DeepCSVb, float
             DeepCSVc, float DeepCSVl, float DeepCSVbb , float DeepCSVcc , float DeepCSVbN ,
             float DeepCSVcN , float DeepCSVlN , float DeepCSVbbN, float DeepCSVccN, float
             DeepCSVbP , float DeepCSVcP , float DeepCSVlP , float DeepCSVbbP, float
             DeepCSVccP) {
         DeepCSVb_ = DeepCSVb;
         DeepCSVc_ = DeepCSVc;
         DeepCSVl_ = DeepCSVl;
         DeepCSVbb_ = DeepCSVbb;
         DeepCSVcc_ = DeepCSVcc;
         DeepCSVbN_ = DeepCSVbN;
         DeepCSVcN_ = DeepCSVcN;
         DeepCSVlN_ = DeepCSVlN;
         DeepCSVbbN_ = DeepCSVbbN;
         DeepCSVccN_ = DeepCSVccN;
         DeepCSVbP_ = DeepCSVbP;
         DeepCSVcP_ = DeepCSVcP;
         DeepCSVlP_ = DeepCSVlP;
         DeepCSVbbP_ = DeepCSVbbP;
         DeepCSVccP_ = DeepCSVccP;}


     void setFlavour(float fpartonflavour, float fhadronflavour) {partonFlavour_ = fpartonflavour; hadronFlavour_ = fhadronflavour;}

     void setQGTagger(float fQGTagger) {QGtagger_ = fQGTagger;}

     void setBoosted(float fboosted) {boosted_ = fboosted;}
     void setCTagger(float fpfCombinedCvsL, float fpfCombinedCvsB) {pfCombinedCvsL_ = fpfCombinedCvsL; pfCombinedCvsB_ = fpfCombinedCvsB;}

     //------------ Get methods ------------------------------
     /*
     float tche()     const {return TCHE_;}
     float tchp()     const {return TCHP_;}
     float tchepf()      const {return TCHEpf_;}
     float tchppf()      const {return TCHPpf_;}
     float softmuonbyip()      const {return SoftMuonTagByIP_;}
     float softelectronbyip()      const {return SoftElectronTagByIP_;}
     float softmuon()      const {return SoftMuonTag_;}
     float softelectron()      const {return SoftElectronTag_;}
     float simplesecvertexhe()      const {return SimpleSecVertexHE_;}
     float simplesecvertexhp()      const {return SimpleSecVertexHP_;}
     float simplesecvertexhepf()      const {return SimpleSecVertexHEpf_;}
     float simplesecvertexhppf()      const {return SimpleSecVertexHPpf_;}
     float csv()      const {return CSV_;}
     float csvpf()      const {return CSVpf_;}
     float cinclsvpf()      const {return CinclSVpf_;}
     float cmvapf()      const {return CMVApf_;}
     float csvsoftleptonpf()      const {return CSVSoftLeptonpf_;}
     */
     float csvpfpositive()      const {return CSVpfPositive_;}
     float csvpfnegative()      const {return CSVpfNegative_;}
    
     float pfBoostedDouble()  const {return boosted_;} 

     float partonflavour()      const {return partonFlavour_;}
     float hadronflavour()      const {return hadronFlavour_;}

     float qgtagger()      const {return QGtagger_;}

     float pfJetProbabilityBJetTags() const {return recommend1_;}
     float pfCombinedInclusiveSecondaryVertexV2BJetTags() const {return recommend2_;}
     float pfCombinedMVAV2BJetTags() const {return recommend2_;}
     float pfCombinedCvsL() const {return pfCombinedCvsL_;}
     float pfCombinedCvsB() const {return pfCombinedCvsB_;}

   private:
     /*
     float TCHE_;
     float TCHP_;
     float TCHEpf_;
     float TCHPpf_;
     float SoftMuonTagByIP_;
     float SoftElectronTagByIP_;
     float SoftMuonTag_;
     float SoftElectronTag_;
     float SimpleSecVertexHE_;
     float SimpleSecVertexHP_;
     float SimpleSecVertexHEpf_;
     float SimpleSecVertexHPpf_;
     float CSV_;
     float CSVpf_;
     float CinclSVpf_;
     float CMVApf_;
     float CSVSoftLeptonpf_;
     */
     float CSVpfPositive_;
     float CSVpfNegative_;
     
     float boosted_;

     float QGtagger_;

     float partonFlavour_;
     float hadronFlavour_;
     float recommend1_;
     float recommend2_;
     float recommend3_;
     // ctaggers
     float pfCombinedCvsL_;
     float pfCombinedCvsB_;
     float DeepCSVb_, DeepCSVc_, DeepCSVl_, DeepCSVbb_, DeepCSVcc_, DeepCSVbN_,
           DeepCSVcN_, DeepCSVlN_, DeepCSVbbN_, DeepCSVccN_,  DeepCSVbP_, DeepCSVcP_,
           DeepCSVlP_, DeepCSVbbP_, DeepCSVccP_;

     

    };
#endif
