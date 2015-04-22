#include <iostream>
#include <sstream>
#include <istream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <functional>
#include <vector>
#include <cassert>
#include "TMath.h"
#include "TRandom.h"
#include "TFile.h"

#include "SMPJ/AnalysisFW/plugins/Analysis_Template_MC.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"

using namespace std;

//---------------------------- Constructor Of The Class TriggerTurnOn -------------------------- //
Analysis_Template_MC::Analysis_Template_MC(edm::ParameterSet const& cfg)
{
     mFileName       = cfg.getParameter<std::string>               ("filename");
     mTreeName       = cfg.getParameter<std::string>               ("treename");
     mDirName        = cfg.getParameter<std::string>               ("dirname");

     mGlobalTag        = cfg.getParameter<std::string>               ("pseudoglobaltag");
     mjettype        = cfg.getParameter<std::string>               ("jettype");

     mMinPt     = cfg.getParameter<double> ("minPt");
     mYMax      = cfg.getParameter<double> ("ymax");
     mJetID     = cfg.getParameter<int>  ("JetID");

     mprintOk   = cfg.getParameter<int>  ("printOk");

     mIsMCarlo       = cfg.getUntrackedParameter<bool>             ("isMCarlo");
     mJECUncSrcNames = cfg.getParameter<std::vector<std::string> > ("jecUncSrcNames");
}

//------------------------------ Declaration Of The Function beginjob() ------------------------//
void Analysis_Template_MC::beginJob()
 {

     mInf = TFile::Open(mFileName.c_str());
     mDir = (TDirectoryFile*)mInf->Get(mDirName.c_str());
     mTree = (TTree*)mDir->Get(mTreeName.c_str());
     Event = new QCDEvent();
     TBranch *branch = mTree->GetBranch("events");
     branch->SetAddress(&Event);

     //------------------ Init jec constructor --------------------------- //
     jecs = new JECs(mIsMCarlo, mGlobalTag, mjettype);

     //------------------ Histogram Booking --------------------------- //
     num_of_Vtx     = fs->make<TH1F>("num_of_Vtx","num_of_Vtx",100,0.,100.);
     num_of_VtxGood = fs->make<TH1F>("num_of_VtxGood","num_of_VtxGood",100,0.,100.);

     mc_pthat = fs->make<TH1F>("mc_pthat","mc_pthat",200,0.,2000.);
     mc_pthat_weighted = fs->make<TH1F>("mc_pthat_weighted","mc_pthat_weighted",200,0.,2000.);
     mc_pthat_weighted->Sumw2();


     pt0_GENJet  = fs->make<TH1F>("pt0_GENJet","pt0_GENJet",200,0.,2000.); pt0_GENJet->Sumw2();
     pt1_GENJet  = fs->make<TH1F>("pt1_GENJet","pt1_GENJet",200,0.,2000.); pt1_GENJet->Sumw2();
     y0_GENJet = fs->make<TH1F>("y0_GENJet","y0_GENJet",60,-3.,3.); y0_GENJet->Sumw2();
     y1_GENJet = fs->make<TH1F>("y1_GENJet","y1_GENJet",60,-3.,3.); y1_GENJet->Sumw2();
     phi0_GENJet = fs->make<TH1F>("phi0_GENJet","phi0_GENJet",60, -TMath::Pi(),TMath::Pi()); phi0_GENJet->Sumw2();
     phi1_GENJet = fs->make<TH1F>("phi1_GENJet","phi1_GENJet",60, -TMath::Pi(),TMath::Pi()); phi1_GENJet->Sumw2();

 } // end of function beginJob()





 //------------------------ endjob() function declaration ---------------------- //
 void Analysis_Template_MC::endJob()
 {
   mInf->Close();




 } // closing endJob()





 //--------------------------- analyze() fuction declaration ------------------ //
void Analysis_Template_MC::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup)
 {

  unsigned NEntries = mTree->GetEntries();
  cout<<"Reading TREE: "<<NEntries<<" events"<<endl;

   int decade = 0 ;

   float hweight=1.;  ///Initial value to one

   for(unsigned  l=0; l<NEntries; l++) {
   //for(unsigned  l=0; l<5; l++) {


    //----------- progress report -------------
    double progress = 10.0*l/(1.0*NEntries);
    int k = TMath::FloorNint(progress);
    if (k > decade)
      cout<<10*k<<" %"<<endl;
    decade = k;
    //----------- read the event --------------
    mTree->GetEntry(l);

    ///pthat & mc_weight
    float pthat = Event->evtHdr().pthat();
    float mc_weight = Event->evtHdr().weight();
    if(mprintOk==1) printf("\npthat=%f  mc_weight=%e\n",pthat,mc_weight);
    hweight=mc_weight;
    mc_pthat->Fill(pthat);
    mc_pthat_weighted->Fill(pthat,hweight);

    /////////////////////////////////////////////////////////////////////////////////////////////
    ///Examine GenJets
    unsigned n_genJets = Event->nGenJets();


    /// Dump all GEN jets
    if(mprintOk==1){
       printf("Number of GENJets=%d\n",n_genJets);
       for(unsigned j=0; j<n_genJets; ++j){
          printf("j=%2d  pt=%8.3f  y=%6.3f  phi=%6.3f\n",j,Event->genjet(j).pt(),Event->genjet(j).Rapidity(),Event->genjet(j).phi());
       }
    }

    ///Apply Jet cuts.  Very General to all existing GEN Jets
    int GENjet_ok[100]; for(int ii=0;ii<100;++ii){GENjet_ok[ii]=0;}

    for(unsigned j=0; j< n_genJets; ++j){
	if(Event->genjet(j).pt()<mMinPt) continue;
	if(fabs(Event->genjet(j).Rapidity())>mYMax) continue;
	GENjet_ok[j]=1;
    }



    /// Keep events where leading Jet[0] and Jet[1] survived cuts
    if((GENjet_ok[0]==1)&&(GENjet_ok[1]==1)) {

       ///////////////////////////////////// Measurement with Gen Jets ///////////////////////////////////////////////////

            float ptmax_gen=Event->genjet(0).pt();

            if((ptmax_gen>=200.)){
	      pt0_GENJet->Fill(Event->genjet(0).pt(),hweight);
	      pt1_GENJet->Fill(Event->genjet(1).pt(),hweight);
	      y0_GENJet->Fill(Event->genjet(0).Rapidity(),hweight);
	      y1_GENJet->Fill(Event->genjet(1).Rapidity(),hweight);
	      phi0_GENJet->Fill(Event->genjet(0).phi(),hweight);
	      phi1_GENJet->Fill(Event->genjet(1).phi(),hweight);

	    }


    } //end of GEN Jets


    /////////////////////////////////////////////////////////////////////////////////////////////
    /// PFJets
    /////////////////////////////////////////////Vertex!!!!/////////////////////////////////////
    unsigned n_PFJets = Event->nPFJets();

    jecs->JEC_corrections(Event, n_PFJets, mIsMCarlo);


    /// Vertex selection
    if(mprintOk==1) cout<<"Vertex info: numVtx="<<Event->evtHdr().nVtx()<<"  numVtxGood="<<Event->evtHdr().nVtxGood()<<"  isPVgood()="<<Event->evtHdr().isPVgood()<<"   pfRho="<<Event->evtHdr().pfRho()<<endl;

    num_of_Vtx->Fill(Event->evtHdr().nVtx());
    /// Keep events with PVgood
    if (Event->evtHdr().isPVgood() != 1) continue;

    num_of_VtxGood->Fill(Event->evtHdr().nVtxGood());

    /// Dump all jets No cuts only isPVgood()
    if(mprintOk==1){
       printf("Number of PFJets=%d\n",n_PFJets);
       for(unsigned j=0; j<n_PFJets; ++j){
          printf("j=%2d  pt=%8.3f  y=%6.3f  phi=%6.3f   cor=%6.3f   tightID=%d\n",j,Event->pfjet(j).ptCor(),Event->pfjet(j).y(),Event->pfjet(j).phi(),Event->pfjet(j).cor(),Event->pfjet(j).tightID());
       }
    }

} // end of event loop


} // closing analyze() function



Analysis_Template_MC::~Analysis_Template_MC()
{
}


DEFINE_FWK_MODULE(Analysis_Template_MC);

