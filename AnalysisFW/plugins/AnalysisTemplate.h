#ifndef My_azimuthal_MC_h
#define My_azimuthal_MC_h
 
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
//#include "PhysicsTools/Utilities/interface/LumiReweighting.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TProfile.h"
#include <TMath.h>
#include <boost/shared_ptr.hpp>
using namespace edm;
using namespace std;

#include "SMPJ/AnalysisFW/plugins/JECs.h"
//#include "SMPJ/AnalysisFW/plugins/JECs_S15.h"

//#include "../Unfolding/RooUnfold-1.1.1/src/RooUnfold.h"
//#include "../Unfolding/RooUnfold-1.1.1/src/RooUnfoldResponse.h"

  
class Analysis_Template_MC : public edm::EDAnalyzer
 {
  public:
    explicit Analysis_Template_MC(edm::ParameterSet const& cfg);
    virtual void beginJob();
    virtual void analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup);
    virtual void endJob();
    virtual ~Analysis_Template_MC();

  private:
    //---- configurable parameters --------    
    std::string mjettype,mGlobalTag,mTreeName,mDirName;
    std::vector<std::string> mFileName;
    double mMinPt, mYMax;
    int    mJetID;        // looseID==1 tightID==2
    int    mprintOk;       // noPrint=0  Print=1
    int    mMCSlice;       // noPrint=0  Print=1
    bool mIsMCarlo;
    bool mPUReweighting;
    bool mLowPileUp;
    std::vector<std::string> mJECUncSrcNames;
    string mJECUncSrc;    

    edm::Service<TFileService> fs;
    //std::vector<TTree> *mTree;
    //std::vector<TFile> *mInf;
    TFile *mPuf,*mInf;
    TTree *mTree;
    //std::vector<TDirectoryFile> *mDir;
    TDirectoryFile *mDir;

    //---- TREE variable --------
    QCDEvent *Event;

    JECs *jecs;

    //--------- Histogram Declaration --------------------//
    // Vertices
    TH1F *num_of_Vtx;
    TH1F *num_of_VtxGood;
    
    TH1F *mc_pthat;
    TH1F *mc_pthat_weighted;   

    ///Measurement Gen jets
    TH1F *pt0_GENJet;
    TH1F *pt1_GENJet;
    TH1F *y0_GENJet;
    TH1F *y1_GENJet;
    TH1F *phi0_GENJet;
    TH1F *phi1_GENJet;

    TH1F *Multiplicity_GENJet;

    ///Measurement Det jets                                                                                                                                              
    TH1F *pt0_DETJet;
    TH1F *pt1_DETJet;
    TH1F *pt0_DETJetUncor;
    TH1F *pt1_DETJetUncor;
    TH1F *y0_DETJet;
    TH1F *y1_DETJet;
    TH1F *phi0_DETJet;
    TH1F *phi1_DETJet; 

    TH1F *Multiplicity_DETJet;

    TH1F *pt0_DETInclJet;
    TH1F *pt0_DETInclJetUncor;
    TH1F *y0_DETInclJet;
    TH1F *phi0_DETInclJet;
  
    // Purity and stability
    TH1F *PF_MatchedInclusiveJets;
    TH1F *Gen_MatchedInclusiveJets;
    TH2F *TwoD_MatchedInclusiveJets;
    TH1F *PF_FakeInclusiveJets;

    TH1F *Gen_MissInclusiveJets;
    TH1F *DeltaR_Jets;

    TH2F *PileUpVSVertex;

    TH1F *AcceptancePtJets;
    TH1F *PurityPtJets;
    TH1F *BackgroundPtJets;
    TH1F *StabilityPtJets;

    TH1F *AcceptancePtJets_1bin;
    TH1F *PurityPtJets_1bin;
    TH1F *AcceptancePtJets_2bin;
    TH1F *PurityPtJets_2bin;
    TH1F *AcceptancePtJets_3bin;
    TH1F *PurityPtJets_3bin;
    TH1F *AcceptancePtJets_4bin;
    TH1F *PurityPtJets_4bin;
    TH1F *AcceptancePtJets_5bin;
    TH1F *PurityPtJets_5bin;
    TH1F *AcceptancePtJets_6bin;
    TH1F *PurityPtJets_6bin;
    TH1F *AcceptancePtJets_7bin;
    TH1F *PurityPtJets_7bin;

    TH1F *BackgroundPtJets_1bin;
    TH1F *StabilityPtJets_1bin;
    TH1F *BackgroundPtJets_2bin;
    TH1F *StabilityPtJets_2bin;
    TH1F *BackgroundPtJets_3bin;
    TH1F *StabilityPtJets_3bin;
    TH1F *BackgroundPtJets_4bin;
    TH1F *StabilityPtJets_4bin;
    TH1F *BackgroundPtJets_5bin;
    TH1F *StabilityPtJets_5bin;
    TH1F *BackgroundPtJets_6bin;
    TH1F *StabilityPtJets_6bin;
    TH1F *BackgroundPtJets_7bin;
    TH1F *StabilityPtJets_7bin;
 
    TH1F* pt_DETInclJet_1bin;
    TH1F* pt_DETInclJet_2bin;
    TH1F* pt_DETInclJet_3bin;
    TH1F* pt_DETInclJet_4bin;
    TH1F* pt_DETInclJet_5bin;
    TH1F* pt_DETInclJet_6bin;
    TH1F* pt_DETInclJet_7bin;

    TH1F* pt_DETInclJetCrossSectNorm_1bin;
    TH1F* pt_DETInclJetCrossSectNorm_2bin;
    TH1F* pt_DETInclJetCrossSectNorm_3bin;
    TH1F* pt_DETInclJetCrossSectNorm_4bin;
    TH1F* pt_DETInclJetCrossSectNorm_5bin;
    TH1F* pt_DETInclJetCrossSectNorm_6bin;
    TH1F* pt_DETInclJetCrossSectNorm_7bin;

    TH1F* pt_GENInclJetCrossSectNorm_1bin;
    TH1F* pt_GENInclJetCrossSectNorm_2bin;
    TH1F* pt_GENInclJetCrossSectNorm_3bin;
    TH1F* pt_GENInclJetCrossSectNorm_4bin;
    TH1F* pt_GENInclJetCrossSectNorm_5bin;
    TH1F* pt_GENInclJetCrossSectNorm_6bin;
    TH1F* pt_GENInclJetCrossSectNorm_7bin;

    TH1F* pt_DETInclJet60_1bin;
    TH1F* pt_DETInclJet60_2bin;
    TH1F* pt_DETInclJet60_3bin;
    TH1F* pt_DETInclJet60_4bin;
    TH1F* pt_DETInclJet60_5bin;
    TH1F* pt_DETInclJet60_6bin;
    TH1F* pt_DETInclJet60_7bin;

    TH1F* pt_DETInclJet80_1bin;
    TH1F* pt_DETInclJet80_2bin;
    TH1F* pt_DETInclJet80_3bin;
    TH1F* pt_DETInclJet80_4bin;
    TH1F* pt_DETInclJet80_5bin;
    TH1F* pt_DETInclJet80_6bin;
    TH1F* pt_DETInclJet80_7bin;

    TH1F* pt_DETInclJet140_1bin;
    TH1F* pt_DETInclJet140_2bin;
    TH1F* pt_DETInclJet140_3bin;
    TH1F* pt_DETInclJet140_4bin;
    TH1F* pt_DETInclJet140_5bin;
    TH1F* pt_DETInclJet140_6bin;
    TH1F* pt_DETInclJet140_7bin;

    TH1F* pt_DETInclJet200_1bin;
    TH1F* pt_DETInclJet200_2bin;
    TH1F* pt_DETInclJet200_3bin;
    TH1F* pt_DETInclJet200_4bin;
    TH1F* pt_DETInclJet200_5bin;
    TH1F* pt_DETInclJet200_6bin;
    TH1F* pt_DETInclJet200_7bin;
    
    TH1F* pt_DETInclJet260_1bin;
    TH1F* pt_DETInclJet260_2bin;
    TH1F* pt_DETInclJet260_3bin;
    TH1F* pt_DETInclJet260_4bin;
    TH1F* pt_DETInclJet260_5bin;
    TH1F* pt_DETInclJet260_6bin;
    TH1F* pt_DETInclJet260_7bin;

    TH1F* pt_DETInclJet320_1bin;
    TH1F* pt_DETInclJet320_2bin;
    TH1F* pt_DETInclJet320_3bin;
    TH1F* pt_DETInclJet320_4bin;
    TH1F* pt_DETInclJet320_5bin;
    TH1F* pt_DETInclJet320_6bin;
    TH1F* pt_DETInclJet320_7bin;

    TH1F* pt_DETInclJet400_1bin;
    TH1F* pt_DETInclJet400_2bin;
    TH1F* pt_DETInclJet400_3bin;
    TH1F* pt_DETInclJet400_4bin;
    TH1F* pt_DETInclJet400_5bin;
    TH1F* pt_DETInclJet400_6bin;
    TH1F* pt_DETInclJet400_7bin;

    TH1F* pt_DETInclJet450_1bin;
    TH1F* pt_DETInclJet450_2bin;
    TH1F* pt_DETInclJet450_3bin;
    TH1F* pt_DETInclJet450_4bin;
    TH1F* pt_DETInclJet450_5bin;
    TH1F* pt_DETInclJet450_6bin;
    TH1F* pt_DETInclJet450_7bin;

    TH1F* pt_GENInclJet_1bin;
    TH1F* pt_GENInclJet_2bin;
    TH1F* pt_GENInclJet_3bin;
    TH1F* pt_GENInclJet_4bin;
    TH1F* pt_GENInclJet_5bin;
    TH1F* pt_GENInclJet_6bin;
    TH1F* pt_GENInclJet_7bin;

    TH1F* PF_MatchedInclusiveJets_1bin;
    TH1F* PF_MatchedInclusiveJets_2bin;
    TH1F* PF_MatchedInclusiveJets_3bin;
    TH1F* PF_MatchedInclusiveJets_4bin;
    TH1F* PF_MatchedInclusiveJets_5bin;
    TH1F* PF_MatchedInclusiveJets_6bin;
    TH1F* PF_MatchedInclusiveJets_7bin;

    TH1F* Gen_MatchedInclusiveJets_1bin;
    TH1F* Gen_MatchedInclusiveJets_2bin;
    TH1F* Gen_MatchedInclusiveJets_3bin;
    TH1F* Gen_MatchedInclusiveJets_4bin;
    TH1F* Gen_MatchedInclusiveJets_5bin;
    TH1F* Gen_MatchedInclusiveJets_6bin;
    TH1F* Gen_MatchedInclusiveJets_7bin;

    TH1F* PF_FakeInclusiveJets_1bin;
    TH1F* PF_FakeInclusiveJets_2bin;
    TH1F* PF_FakeInclusiveJets_3bin;
    TH1F* PF_FakeInclusiveJets_4bin;
    TH1F* PF_FakeInclusiveJets_5bin;
    TH1F* PF_FakeInclusiveJets_6bin;
    TH1F* PF_FakeInclusiveJets_7bin;

    TH1F* Gen_MissInclusiveJets_1bin;
    TH1F* Gen_MissInclusiveJets_2bin;
    TH1F* Gen_MissInclusiveJets_3bin;
    TH1F* Gen_MissInclusiveJets_4bin;
    TH1F* Gen_MissInclusiveJets_5bin;
    TH1F* Gen_MissInclusiveJets_6bin;
    TH1F* Gen_MissInclusiveJets_7bin;

    TH2F *TwoD_MatchedInclusiveJets_1bin;
    TH2F *TwoD_MatchedInclusiveJets_2bin;
    TH2F *TwoD_MatchedInclusiveJets_3bin;
    TH2F *TwoD_MatchedInclusiveJets_4bin;
    TH2F *TwoD_MatchedInclusiveJets_5bin;
    TH2F *TwoD_MatchedInclusiveJets_6bin;
    TH2F *TwoD_MatchedInclusiveJets_7bin;

    TH1F* Resolution1D;
    TH1F* ResolutionForward1D;

    TProfile *ResolutionTagAndProbe;

    TProfile *ResolutionInclusiveJets;
    TProfile *ResolutionInclusiveJets_1bin;
    TProfile *ResolutionInclusiveJets_2bin;
    TProfile *ResolutionInclusiveJets_3bin;
    TProfile *ResolutionInclusiveJets_4bin;
    TProfile *ResolutionInclusiveJets_5bin;
    TProfile *ResolutionInclusiveJets_6bin;
    TProfile *ResolutionInclusiveJets_7bin;

    TH1F *hist_leading_pt_emulated_Jet60;
    TH1F *hist_leading_pt_all_Jet60;
    TH1F *hist_leading_eta_emulated_Jet60;
    TH1F *hist_leading_eta_all_Jet60;    
    TH1F *hist_leading_pt_HLT_Jet60U_eff;
    TH1F *hist_leading_eta_HLT_Jet60U_eff;

    TH1F *hist_leading_pt_emulated_Jet80;
    TH1F *hist_leading_pt_all_Jet80;
    TH1F *hist_leading_eta_emulated_Jet80;
    TH1F *hist_leading_eta_all_Jet80;
    TH1F *hist_leading_pt_HLT_Jet80U_eff;
    TH1F *hist_leading_eta_HLT_Jet80U_eff;

    TH1F *hist_leading_pt_emulated_Jet140;
    TH1F *hist_leading_pt_all_Jet140;
    TH1F *hist_leading_eta_emulated_Jet140;
    TH1F *hist_leading_eta_all_Jet140;
    TH1F *hist_leading_pt_HLT_Jet140U_eff;
    TH1F *hist_leading_eta_HLT_Jet140U_eff;

    TH1F *hist_leading_pt_emulated_Jet200;
    TH1F *hist_leading_pt_all_Jet200;
    TH1F *hist_leading_eta_emulated_Jet200;
    TH1F *hist_leading_eta_all_Jet200;
    TH1F *hist_leading_pt_HLT_Jet200U_eff;
    TH1F *hist_leading_eta_HLT_Jet200U_eff;

    TH1F *hist_leading_pt_emulated_Jet260;
    TH1F *hist_leading_pt_all_Jet260;
    TH1F *hist_leading_eta_emulated_Jet260;
    TH1F *hist_leading_eta_all_Jet260;
    TH1F *hist_leading_pt_HLT_Jet260U_eff;
    TH1F *hist_leading_eta_HLT_Jet260U_eff;

    TH1F *hist_leading_pt_emulated_Jet320;
    TH1F *hist_leading_pt_all_Jet320;
    TH1F *hist_leading_eta_emulated_Jet320;
    TH1F *hist_leading_eta_all_Jet320;
    TH1F *hist_leading_pt_HLT_Jet320U_eff;
    TH1F *hist_leading_eta_HLT_Jet320U_eff;

    TH1F *hist_leading_pt_emulated_Jet400;
    TH1F *hist_leading_pt_all_Jet400;
    TH1F *hist_leading_eta_emulated_Jet400;
    TH1F *hist_leading_eta_all_Jet400;
    TH1F *hist_leading_pt_HLT_Jet400U_eff;
    TH1F *hist_leading_eta_HLT_Jet400U_eff;

    TH1F *hist_leading_pt_emulated_Jet450;
    TH1F *hist_leading_pt_all_Jet450;
    TH1F *hist_leading_eta_emulated_Jet450;
    TH1F *hist_leading_eta_all_Jet450;
    TH1F *hist_leading_pt_HLT_Jet450U_eff;
    TH1F *hist_leading_eta_HLT_Jet450U_eff;

    TH1F *hist_leading_pt_emulated_Jet500;
    TH1F *hist_leading_pt_all_Jet500;
    TH1F *hist_leading_eta_emulated_Jet500;
    TH1F *hist_leading_eta_all_Jet500;
    TH1F *hist_leading_pt_HLT_Jet500U_eff;
    TH1F *hist_leading_eta_HLT_Jet500U_eff;

    TH1F* TagAndProbeEff;
    TH1F* TagAndProbeNum;
    TH1F* TagAndProbeDen; 

    TH1F* pt_DETInclJetUP_1bin;
    TH1F* pt_DETInclJetUP_2bin;
    TH1F* pt_DETInclJetUP_3bin;
    TH1F* pt_DETInclJetUP_4bin;
    TH1F* pt_DETInclJetUP_5bin;
    TH1F* pt_DETInclJetUP_6bin;
    TH1F* pt_DETInclJetUP_7bin;

    TH1F* pt_DETInclJetDOWN_1bin;
    TH1F* pt_DETInclJetDOWN_2bin;
    TH1F* pt_DETInclJetDOWN_3bin;
    TH1F* pt_DETInclJetDOWN_4bin;
    TH1F* pt_DETInclJetDOWN_5bin;
    TH1F* pt_DETInclJetDOWN_6bin;
    TH1F* pt_DETInclJetDOWN_7bin;

    TH1F* Chargedhf0_DETJet;
    TH1F* Chargedef0_DETJet;
    TH1F* Neutralhf0_DETJet;
    TH1F* Photonef0_DETJet;
    TH1F* Hadronef0_DETJet;
    TH1F* Muonef0_DETJet;
    TH1F* Electromagneticef0_DETJet;
   
    TH1F* ChargedhMultiplicity0_DETJet;
    TH1F* ChargedeMultiplicity0_DETJet;
    TH1F* NeutralhMultiplicity0_DETJet;
    TH1F* PhotoneMultiplicity0_DETJet;
    TH1F* HadroneMultiplicity0_DETJet;
    TH1F* MuoneMultiplicity0_DETJet;
    TH1F* ElectromagneticeMultiplicity0_DETJet;

    TH1F* Emf0_DETJet;
    TH1F* Hpd0_DETJet;
    TH1F* Chf0_DETJet;
    TH1F* Nhf0_DETJet;
    TH1F* Pef0_DETJet;
    TH1F* Eef0_DETJet;
    TH1F* Mef0_DETJet;

    TH1F* TruePileUpMC;
    TH1F* TruePileUpMCInteger;
    TH1F* TruePileUpDataInteger;

    TH1F* MET_DET;
    TH1F* METPhi_DET;
    TH1F* FractionMET_DET;

    /*RooUnfoldResponse resp_jetpt1etabin;
    RooUnfoldResponse resp_jetpt2etabin;
    RooUnfoldResponse resp_jetpt3etabin;
    RooUnfoldResponse resp_jetpt4etabin;
    RooUnfoldResponse resp_jetpt5etabin;
    RooUnfoldResponse resp_jetpt6etabin;
    RooUnfoldResponse resp_jetpt7etabin;*/

 };

#endif

