#ifndef NewTreeAnalyzer_h
#define NewTreeAnalyzer_h

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SMPJ/AnalysisFW/interface/QCDJet.h"
#include "SMPJ/AnalysisFW/interface/QCDEvent.h"
#include "SMPJ/AnalysisFW/interface/QCDEventHdr.h"
#include "SMPJ/AnalysisFW/interface/QCDCaloJet.h"
#include "SMPJ/AnalysisFW/interface/QCDPFJet.h"
#include "SMPJ/AnalysisFW/interface/QCDMET.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "DataFormats/JetReco/interface/JetExtendedAssociation.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"

//#include "PhysicsTools/Utilities/interface/LumiReweighting.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TProfile.h"
#include "TLorentzVector.h"
#include "TRandom.h"
#include <TMath.h>
#include <boost/shared_ptr.hpp>
using namespace edm;
using namespace std;

 int  HLTJetPtN[7] = {40,80,140,200,260,320,400};
 double HLTJetPtS[7] = {56,133,220,300,395,507,5000};
 double HLTJetPtT[7] = {80,130,200,260,330,400,520};
  int  ATrig[7] = {16,36,68,92,128,128,128} ;
  double EtabinN[7] = {0.0,0.5,1.0,1.5,2.0,2.5,3.0};
  float metHT[4] = {0.0, 0.3, 0.9, 1.0};
  int PUBin[8] = {0,5,10,15,20,25,30,40};

  unsigned  NEVENTSN = 100000;
  double dR = 0.5;
  const int nbin1 = 50, nbin2 = 20, nbin3 = 39;
  double ptMin = 20;
 const double Ptbin1[nbin3+1] = {30,45,60,75,90,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507,548,592,638,686,737,790,846,905,967,1032,1101,1172,1248,1327,1410,1497,1588,1784,2000,2116};

// const int nx[6] = {39,36,33,30,25,19};
 const int nx[6] = {39,39,38,34,27,23};
 const int nGenx[6]  = {45,45,44,40,33,26};
 const int n1x[6] = {47,47,46,42,35,28};

const double x[6][50]= {
 {56,64,74,84,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507,548,592,638,686,737,790,846,905,967,1032,1101,1172,1248,1327,1410,1497,1588,1784,2116,2500},

 {56,64,74,84,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507,548,592,638,686,737,790,846,905,967,1032,1101,1172,1248,1327,1410,1497,1588,1784,2116,2500},

 {56,64,74,84,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507,548,592,638,686,737,790,846,905,967,1032,1101,1172,1248,1327,1410,1497,1588,1784,2116},

 {56,64,74,84,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507,548,592,638,686,737,790,846,905,967,1032,1101,1172,1248,1327,1410},

 {56,64,74,84,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507,548,592,638,686,737,790,846,905},

 {56,64,74,84,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507,548,592,638,686}
 };

const double x1[6][50]= {
 {18,21,24,28,32,37,43,49,56,64,74,84,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507,548,592,638,686,737,790,846,905,967,1032,1101,1172,1248,1327,1410,1497,1588,1784,2116,2500},

 {18,21,24,28,32,37,43,49,56,64,74,84,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507,548,592,638,686,737,790,846,905,967,1032,1101,1172,1248,1327,1410,1497,1588,1784,2116,2500},

 {18,21,24,28,32,37,43,49,56,64,74,84,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507,548,592,638,686,737,790,846,905,967,1032,1101,1172,1248,1327,1410,1497,1588,1784,2116},

 {18,21,24,28,32,37,43,49,56,64,74,84,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507,548,592,638,686,737,790,846,905,967,1032,1101,1172,1248,1327,1410},

 {18,21,24,28,32,37,43,49,56,64,74,84,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507,548,592,638,686,737,790,846,905},

 {18,21,24,28,32,37,43,49,56,64,74,84,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507,548}
 };

const double xGen[6][50]= {
 {28,32,37,43,49,56,64,74,84,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507,548,592,638,686,737,790,846,
 905,967,1032,1101,1172,1248,1327,1410,1497,1588,1784,2116,2500, 3000},

 {28,32,37,43,49,56,64,74,84,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507,548,592,638,686,737,790,846,
 905,967,1032,1101,1172,1248,1327,1410,1497,1588,1784,2116,2500, 3000},

 {28,32,37,43,49,56,64,74,84,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507,548,592,638,686,737,790,846,
 905,967,1032,1101,1172,1248,1327,1410,1497,1588,1784,2116, 2500},

 {28,32,37,43,49,56,64,74,84,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507,548,592,638,686,737,790,846,
 905,967,1032,1101,1172,1248,1327,1410, 1497},

 {28,32,37,43,49,56,64,74,84,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507,548,592,638,686,737,790,846,
 905, 967},

 {28,32,37,43,49,56,64,74,84,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507,548, 592}
 };

int  VertexN[4] = {0,10,13,10000};

  double DeltaPhi(double phi1, double phi2)
     {
          if (fabs(phi1 - phi2) <= M_PI)
           {
             return ( fabs(phi1 - phi2) ) ;
        }

       else
        {
         return (2*M_PI - fabs(phi1 - phi2));
       }

   }



class NewTreeAnalyzer : public edm::EDAnalyzer
 {
  public:
    explicit NewTreeAnalyzer(edm::ParameterSet const& cfg);
    virtual void beginJob();
    virtual void analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup);
    virtual void endJob();
    virtual ~NewTreeAnalyzer();
    typedef reco::Particle::LorentzVector LorentzVector;
  private:
     int mJetID,mHCALNoise,mNEvents;
     static const int nJetTrig = 7, nJetVersn = 7, nMuonTrig = 10 , nMuonVersn = 10, nEtabin = 6, nvtx = 3,
     nMetBin = 6, nUnc = 52, nPU = 7;
     double mMinPt;
     std::vector<int> mRefTrigIndex;
     std::vector<double> mEtabin,mL1Pt,mHLTPt,mVertex;
     std::vector<std::string> mRefTriggerJet,mRefTriggerMuon;
     std::string mFileName,mTreeName,mDirName,mPuFileName,mPuTrigName;
     std::string mPFPayloadNameA, mPFPayloadNameCHSA, mPFJECUncSrcA, mPFJECUncSrcCHSA;
     std::string mJECL1FastFile, mJECL1FastFileCHS, mJECL2RelativeFile, mJECL2RelativeFileCHS,
                 mJECL3AbsoluteFile, mJECL3AbsoluteFileCHS, mJECL2L3ResidualFile, mJECL2L3ResidualFileCHS;
     static bool sort_pfjets(QCDPFJet j1, QCDPFJet j2) {
      return j1.ptCor() > j2.ptCor();
    }

 //    char mPuTrigName;


    // ---- Jet Corrector Parameter ---- //
    JetCorrectorParameters *L1Fast, *L2Relative, *L3Absolute, *L2L3Residual;
    vector<JetCorrectorParameters> vecL1Fast, vecL2Relative, vecL3Absolute, vecL2L3Residual;
    FactorizedJetCorrector *jecL1Fast, *jecL2Relative, *jecL3Absolute, *jecL2L3Residual;

    // ======== ***************** ======= //
     edm::Service<TFileService> fs;
     TTree *mTree;
     TFile *mInf, *mPuf;
     TDirectoryFile *mDir;
     bool mIsMCarlo, mIsCHS, isPFJecUncSet_;
     std::vector<std::string> mJECUncSrcNames;
     // - jec unc sources ----//
     JetCorrectionUncertainty *mPFUnc;
     std::vector<JetCorrectionUncertainty*> mPFUncSrc;

   //---- TREE variable --------
     QCDEvent *Event;

    TString HLTJet[nJetTrig][nJetVersn],  HLTMuon[nMuonTrig][nMuonVersn];

     char name[200],title[200],trigtitle[200];

     int  ihltj[nJetTrig][nJetVersn] ,prescalej[nJetTrig][nJetVersn];

     bool hltPassj[nJetTrig][nJetVersn],  hltPassm[nMuonTrig][nMuonVersn];

  //------- PU Reweighting Function ------------//
   void LumiReWeighting( std::vector< float > MC_distr, std::vector< float > Lumi_distr);
   TH1F *weights_, *MC_distr_, *Data_distr_, *den;
   double puweight(float npv);
   std::vector<float> WSummer2012;
   std::vector<float> WData2012;


   //--------- Histogram Declaration --------------------//

 TH1F *PFPt[nEtabin][nJetTrig][nJetVersn], *PFPtUp[nEtabin][nUnc][nJetTrig][nJetVersn], *PFPtDown[nEtabin][nUnc][nJetTrig][nJetVersn];

 TH1F *PFPt_Tot[nEtabin], *PFPtUp_Tot[nEtabin][nUnc], *PFPtDown_Tot[nEtabin][nUnc], *PFPtUp_Ratio_Tot[nEtabin][nUnc], *PFPtDown_Ratio_Tot[nEtabin][nUnc];


    TH1F *PFJet0[nJetTrig][nJetVersn][nEtabin][nvtx] , *PFJet1[nJetTrig][nJetVersn][nEtabin][nvtx], *PFJet2[nJetTrig][nJetVersn][nEtabin], *PFJet2Eta[nJetTrig][nJetVersn][nEtabin], *PFJet2Phi[nJetTrig][nJetVersn][nEtabin];
    TH1F *PFJetMJJ0[nJetTrig][nJetVersn][nEtabin][nvtx] , *PFJetMJJ1[nJetTrig][nJetVersn][nEtabin][nvtx];

    TH1F *PFJet[nJetTrig][nJetVersn][nEtabin], *PFJetYield[nJetTrig][nJetVersn][nEtabin], *PFJetL[nJetTrig][nJetVersn], *PFJetGen[nJetTrig][nJetVersn][nEtabin];

    TH1F *PFJetPU[nJetTrig][nJetVersn][nEtabin][nPU];

    TH1F *ChHadFr[nJetTrig][nJetVersn][nEtabin], *NuHadFr[nJetTrig][nJetVersn][nEtabin], *PhFr[nJetTrig][nJetVersn][nEtabin], *ElFr[nJetTrig][nJetVersn][nEtabin], *MuFr[nJetTrig][nJetVersn][nEtabin] ;
   TH1F *ChHadFrPT[nJetTrig][nJetVersn][nEtabin], *NuHadFrPT[nJetTrig][nJetVersn][nEtabin], *PhFrPT[nJetTrig][nJetVersn][nEtabin], *ElFrPT[nJetTrig][nJetVersn][nEtabin], *MuFrPT[nJetTrig][nJetVersn][nEtabin] ;


  TH2F* PrescalevsRun[nJetTrig][nJetVersn], *BetavsPT[nJetTrig][nJetVersn][nEtabin]; TH1F *MetHT[nJetTrig][nJetVersn], *NPV[nJetTrig][nJetVersn];

  TH1F *TagJet[nJetTrig][nJetVersn][nEtabin], *ProbeJet[nJetTrig][nJetVersn][nEtabin], *PV, *GenJet[nEtabin];
  TH1F *Unfolding_GenPt[nEtabin], *Unfolding_RecoPt[nEtabin];
  TH2F *Unfolding_RecoGenPt[nEtabin];
  TH1F *PtRes[nEtabin][nbin1], *hAux[nJetTrig][nJetVersn][nEtabin], *hAux1[nJetTrig][nJetVersn][nEtabin]; // Auxiliary histos fo jet fraction
  TH1F *MPFEventPT[nJetTrig][nJetVersn][nEtabin], *MPFEventCountPT[nJetTrig][nJetVersn][nEtabin][100];
  TProfile *AvgX_Pt[nEtabin][nbin1], *AvgCHF_NPV[nJetTrig][nJetVersn][nEtabin], *AvgNHF_NPV[nJetTrig][nJetVersn][nEtabin], *AvgPHF_NPV[nJetTrig][nJetVersn][nEtabin];
  TProfile *PFJetPt_L1Fast[nJetTrig][nJetVersn][nEtabin], *PFJetPt_L2Relative[nJetTrig][nJetVersn][nEtabin], *PFJetPt_L3Absolute[nJetTrig][nJetVersn][nEtabin],
           *PFJetPt_L2L3Residual[nJetTrig][nJetVersn][nEtabin], *PFJetPt_TotalJEC[nJetTrig][nJetVersn][nEtabin];

 };

#endif
