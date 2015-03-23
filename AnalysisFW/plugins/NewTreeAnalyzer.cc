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

#include "SMPJ/AnalysisFW/plugins/NewTreeAnalyzer.h"




//---------------------------- Constructor Of The Class NewTreeAnalyzer -------------------------- //
NewTreeAnalyzer::NewTreeAnalyzer(edm::ParameterSet const& cfg)
{
     mFileName               = cfg.getParameter<std::string>               ("filename");
     mPuFileName             = cfg.getParameter<std::string>               ("pileupfile");
     mPuTrigName             = cfg.getParameter<std::string>               ("pileuptrig");
     mTreeName               = cfg.getParameter<std::string>               ("treename");
     mDirName                = cfg.getParameter<std::string>               ("dirname");
     mIsMCarlo               = cfg.getUntrackedParameter<bool>             ("isMCarlo");
     mIsCHS                  = cfg.getUntrackedParameter<bool>             ("isCHS");
     mPFPayloadNameA         = cfg.getParameter<std::string>               ("PFPayloadName");
     mPFPayloadNameCHSA      = cfg.getParameter<std::string>               ("PFPayloadNameCHS");
     mPFJECUncSrcA           = cfg.getParameter<std::string>               ("jecUncSrc");
     mPFJECUncSrcCHSA        = cfg.getParameter<std::string>               ("jecUncSrcCHS");
     mJECL1FastFile          = cfg.getParameter<std::string>               ("jecL1FastFile");
     mJECL1FastFileCHS       = cfg.getParameter<std::string>               ("jecL1FastFileCHS");
     mJECL2RelativeFile      = cfg.getParameter<std::string>               ("jecL2RelativeFile");
     mJECL2RelativeFileCHS   = cfg.getParameter<std::string>               ("jecL2RelativeFileCHS");
     mJECL3AbsoluteFile      = cfg.getParameter<std::string>               ("jecL3AbsoluteFile");
     mJECL3AbsoluteFileCHS   = cfg.getParameter<std::string>               ("jecL3AbsoluteFileCHS");
     mJECL2L3ResidualFile    = cfg.getParameter<std::string>               ("jecL2L3ResidualFile");
     mJECL2L3ResidualFileCHS = cfg.getParameter<std::string>               ("jecL2L3ResidualFileCHS");
     mJECUncSrcNames         = cfg.getParameter<std::vector<std::string> > ("jecUncSrcNames");
}

//------------------------------ Declaration Of The Function beginjob() ------------------------//
void NewTreeAnalyzer::beginJob()
 {

     mInf = TFile::Open(mFileName.c_str());
     mDir = (TDirectoryFile*)mInf->Get(mDirName.c_str());
     mTree = (TTree*)mDir->Get(mTreeName.c_str());
     Event = new QCDEvent();
     TBranch *branch = mTree->GetBranch("events");
     branch->SetAddress(&Event);
     isPFJecUncSet_ = false;

     // --- Initializing the jet correctors ---- //
     if(mIsCHS) {
     L1Fast       = new JetCorrectorParameters(mJECL1FastFileCHS.c_str());
     L2Relative   = new JetCorrectorParameters(mJECL2RelativeFileCHS.c_str());
     L3Absolute   = new JetCorrectorParameters(mJECL3AbsoluteFileCHS.c_str());
     if(!mIsMCarlo)
     L2L3Residual = new JetCorrectorParameters(mJECL2L3ResidualFileCHS.c_str());
     } // if(mIsCHS)
     else {
     L1Fast       = new JetCorrectorParameters(mJECL1FastFile.c_str());
     L2Relative   = new JetCorrectorParameters(mJECL2RelativeFile.c_str());
     L3Absolute   = new JetCorrectorParameters(mJECL3AbsoluteFile.c_str());
     if(!mIsMCarlo)
     L2L3Residual = new JetCorrectorParameters(mJECL2L3ResidualFile.c_str());
     } // else -- non-chs

     vecL1Fast.push_back(*L1Fast);
     vecL2Relative.push_back(*L2Relative);
     vecL3Absolute.push_back(*L3Absolute);
     if(!mIsMCarlo)
     vecL2L3Residual.push_back(*L2L3Residual);
    
     jecL1Fast       = new FactorizedJetCorrector(vecL1Fast);
     jecL2Relative   = new FactorizedJetCorrector(vecL2Relative);
     jecL3Absolute   = new FactorizedJetCorrector(vecL3Absolute);
     if(!mIsMCarlo)
     jecL2L3Residual = new FactorizedJetCorrector(vecL2L3Residual);
     

    if(!mIsMCarlo) {
   //------------------------- Initializing The Trigger Variables -------------------------- //
      for (int i=0; i<nJetTrig; i++)
         {
           for (int j=0; j<nJetVersn; j++)
              {
                sprintf(trigtitle, "HLT_PFJet%i_v%i",HLTJetPtN[i], j+3);
                HLTJet[i][j] = trigtitle;
             }
         }

   for (int i=0; i<nJetTrig; i++)
         {
           for (int j=0; j<nJetVersn; j++)
              {
                  ihltj[i][j] = -1 ;
            }
        }


   TH1F *hTrigNames = (TH1F*)mDir->Get("TriggerNames");
   // ------------ Assigning An Integer To Each Trigger Path---------------//
    for(int ibin=0;ibin<hTrigNames->GetNbinsX();ibin++) {
       TString ss(hTrigNames->GetXaxis()->GetBinLabel(ibin+1));

       for (int ii=0; ii<nJetTrig; ii++)
          {
            for (int ij=0; ij<nJetVersn; ij++)
               {
                 if (ss == HLTJet[ii][ij]) {
                 ihltj[ii][ij] = ibin;
                 continue;
               }
            } // for (int ij=0; ij<nJetVersn; ij++)
         } // for (int ii=0; ii<nJetTrig; ii++)
      } // for(int ibin=0;ibin<55;ibin++)
  // ----------------------- Checking For The Trigger Assignment ------------------- // 
    for (int ij=0; ij<nJetTrig; ij++)
        {
         for (int ik=0; ik<nJetVersn; ik++)
            {
              if (ihltj[ij][ik] == -1) {
              cout<<"The requested trigger ("<<HLTJet[ij][ik]<<") is not found "<<endl;
   //         break;
               }
              else {
              cout<<HLTJet[ij][ik]<<" --> "<<ihltj[ij][ik]<<endl;
             }

          }
      }
  }//if(!mIsMCarlo)

     //------------------ Histogram Booking --------------------------- //

       if(!mIsMCarlo) {
      //---------- For JES UNCERTAINTY -----------------//
      for (int k=0; k<nEtabin; k++)
        {
         for(int i=0; i<nUnc; i++)
           {
            for (int m=0; m<nJetTrig; m++)
               {
                for (int n=0; n<nJetVersn; n++)
                  {
                   sprintf(name, "PtUp__eta%i_unc%s_%i_v%i",k+1,mJECUncSrcNames[i].c_str(),HLTJetPtN[m],n+3);
                   PFPtUp[k][i][m][n] = fs->make<TH1F>(name,name,nx[k],&x[k][0]);

                   sprintf(name, "PtDown__eta%i_unc%s_%i_v%i",k+1,mJECUncSrcNames[i].c_str(),HLTJetPtN[m],n+3);
                   PFPtDown[k][i][m][n] = fs->make<TH1F>(name,name,nx[k],&x[k][0]);
               } // for (int n=0; n<nJetVersn; n++)              

            } // for (int m=0; m<nJetTrig; m++)

         } // for(int i=0; i<mJECUncSrcNames.size(); i++)

             for (int m=0; m<nJetTrig; m++)
               {
                for (int n=0; n<nJetVersn; n++)
                  {
                    sprintf(name, "PtUp__eta%i_%i_v%i",k+1,HLTJetPtN[m],n+3);
                    PFPt[k][m][n] = fs->make<TH1F>(name,name,nx[k],&x[k][0]);

                    sprintf(name, "BetavsPT_Trig%i_v%i_eta%i", HLTJetPtN[m],n+3, k+1);
                    BetavsPT[m][n][k] = fs->make<TH2F>(name, name, nx[k],&x[k][0], 100,0,1);
                 } // for (int n=0; n<nJetVersn; n++)
                } // for (int m=0; m<nJetTrig; m++)


      } // for (int k=0; k<nEtabin; k++) 
    }//if(!mIsMCarlo) 

     //----------- For X-Section and Triggers .. etc....  -------//
      for (int m=0; m<nJetTrig; m++)
        {
         for (int n=0; n<nJetVersn; n++)
           {
            for (int k=0; k<nEtabin; k++)
              {
               for (int r=0; r<nvtx; r++)
                 {
                   if(!mIsMCarlo) {
                   sprintf(name, "PFJet0_%i_v%i_eta%i_PU%i", HLTJetPtN[m],n+3,k+1,r);
                   PFJet0[m][n][k][r] = fs->make<TH1F>(name,name,nx[k],&x[k][0]);
                   PFJet0[m][n][k][r]->Sumw2();

                   sprintf(name, "PFJet1_%i_v%i_eta%i_PU%i", HLTJetPtN[m],n+3,k+1,r);
                   PFJet1[m][n][k][r] = fs->make<TH1F>(name,name,nx[k],&x[k][0]);
                   PFJet1[m][n][k][r]->Sumw2();
                   }
                  } // for (int r=0; r<nvtx; r++)
                   //----------------- Histogram For JetId Variables ---------------- //
                   sprintf(name,"ChHadFrac_Jet%i_v%i_eta%i",HLTJetPtN[m],n+3,k+1);
                   ChHadFr[m][n][k] = fs->make<TH1F>(name,name,30,0,1);
                   ChHadFr[m][n][k]->Sumw2();

                   sprintf(name,"NuHadFrac_Jet%i_v%i_eta%i",HLTJetPtN[m],n+3,k+1);
                   NuHadFr[m][n][k] = fs->make<TH1F>(name,name,30,0,1);
                   NuHadFr[m][n][k]->Sumw2();

                   sprintf(name,"PhFrac_Jet%i_v%i_eta%i",HLTJetPtN[m],n+3,k+1);
                   PhFr[m][n][k] = fs->make<TH1F>(name,name,30,0,1);
                   PhFr[m][n][k]->Sumw2();

                   sprintf(name,"ElFrac_Jet%i_v%i_eta%i",HLTJetPtN[m],n+3,k+1);
                   ElFr[m][n][k] = fs->make<TH1F>(name,name,30,0,1);
                   ElFr[m][n][k]->Sumw2();

                   sprintf(name,"MuFrac_Jet%i_v%i_eta%i",HLTJetPtN[m],n+3,k+1);
                   MuFr[m][n][k] = fs->make<TH1F>(name,name,30,0,1);
                   MuFr[m][n][k]->Sumw2();

                 //-------------------- Histogram For Spectrum Construction ---------------//
                   sprintf(name, "PFJet_%i_v%i_eta%i", HLTJetPtN[m],n+3,k+1);
                   PFJet[m][n][k] = fs->make<TH1F>(name,name,nx[k],&x[k][0]);
                   PFJet[m][n][k]->Sumw2();

                   sprintf(name, "PFJetYield_%i_v%i_eta%i", HLTJetPtN[m],n+3,k+1);
                   PFJetYield[m][n][k] = fs->make<TH1F>(name,name,nx[k],&x[k][0]);
                   PFJetYield[m][n][k]->Sumw2();

                   sprintf(name, "PFJet_L1Fast_%i_v%i_eta%i", HLTJetPtN[m],n+3,k+1);
                   PFJetPt_L1Fast[m][n][k] = fs->make<TProfile>(name,name,nx[k],&x[k][0],0,5);
                   PFJetPt_L1Fast[m][n][k]->Sumw2();

                   sprintf(name, "PFJet_L2Relative_%i_v%i_eta%i", HLTJetPtN[m],n+3,k+1);
                   PFJetPt_L2Relative[m][n][k] = fs->make<TProfile>(name,name,nx[k],&x[k][0],0,5);
                   PFJetPt_L2Relative[m][n][k]->Sumw2();

                   sprintf(name, "PFJet_L3Absolute_%i_v%i_eta%i", HLTJetPtN[m],n+3,k+1);
                   PFJetPt_L3Absolute[m][n][k] = fs->make<TProfile>(name,name,nx[k],&x[k][0],0,5);
                   PFJetPt_L3Absolute[m][n][k]->Sumw2();

                   sprintf(name, "PFJet_L2L3Residual_%i_v%i_eta%i", HLTJetPtN[m],n+3,k+1);
                   PFJetPt_L2L3Residual[m][n][k] = fs->make<TProfile>(name,name,nx[k],&x[k][0],0,5);
                   PFJetPt_L2L3Residual[m][n][k]->Sumw2();

                   sprintf(name, "PFJet_TotalJEC_%i_v%i_eta%i", HLTJetPtN[m],n+3,k+1);
                   PFJetPt_TotalJEC[m][n][k] = fs->make<TProfile>(name,name,nx[k],&x[k][0],0,5);
                   PFJetPt_TotalJEC[m][n][k]->Sumw2();


                   if(mIsMCarlo) {
                   sprintf(name, "PFJetGen_%i_v%i_eta%i", HLTJetPtN[m],n+3,k+1);
                   PFJetGen[m][n][k] = fs->make<TH1F>(name,name,nx[k],&x[k][0]);
                   PFJetGen[m][n][k]->Sumw2();
                   }

                    if(!mIsMCarlo) {
                    for(int r=0; r<nPU; r++)
                      {
                       sprintf(name, "PFJet_%i_v%i_eta%i_NPU%i", HLTJetPtN[m],n+3,k+1,r);
                       PFJetPU[m][n][k][r] = fs->make<TH1F>(name,name,nx[k],&x[k][0]);
                       PFJetPU[m][n][k][r]->Sumw2();
                     } // for(int r=0; r<nPU; r++)
                    } // if(!mIsMCarlo)

                   sprintf(name, "TagJet_%i_v%i_eta%i", HLTJetPtN[m],n+3,k+1);
                   TagJet[m][n][k] = fs->make<TH1F>(name,name,nx[k],&x[k][0]);
                   TagJet[m][n][k]->Sumw2();

                   sprintf(name, "ProbeJet_%i_v%i_eta%i", HLTJetPtN[m],n+3,k+1);
                   ProbeJet[m][n][k] = fs->make<TH1F>(name,name,nx[k],&x[k][0]);
                   ProbeJet[m][n][k]->Sumw2();
              } // for (int k=0; k<nEtabin; k++)
             } // for (int n=0; n<nJetVersn; n++)
            } // for (int m=0; m<nJetTrig; m++) 


             for (int k=0; k<nEtabin; k++)
              {
              if(mIsMCarlo) {
                for(int i=0; i<n1x[k]; i++)
                  {
                    sprintf(name, "PtRes_eta%i_bin%i", k+1, i+1);
                    PtRes[k][i] = fs->make<TH1F>(name,name,400,0.2,1.8);
                    PtRes[k][i]->Sumw2();

                    sprintf(name, "AvgX_eta%i_bin%i", k+1, i+1);
                    AvgX_Pt[k][i] = fs->make<TProfile>(name,name,n1x[k],0,n1x[k],0,3000);
                    AvgX_Pt[k][i]->Sumw2();
                 } // for(int i=0; i<nbin1; i++)

                   sprintf(name, "GenJet_eta%i", k+1);
                   GenJet[k] = fs->make<TH1F>(name,name,nx[k],&x[k][0]);
                   GenJet[k]->Sumw2();

                   sprintf(name, "Unfolding_GenPt_YBin%i", k+1);
                   Unfolding_GenPt[k] = fs->make<TH1F>(name,name,nGenx[k],&xGen[k][0]);
                   Unfolding_GenPt[k]->Sumw2();

                   sprintf(name, "Unfolding_RecoPt_YBin%i", k+1);
                   Unfolding_RecoPt[k] = fs->make<TH1F>(name,name,nGenx[k],&xGen[k][0]);
                   Unfolding_RecoPt[k]->Sumw2();

                   sprintf(name, "Unfolding_RecoGenPt_YBin%i", k+1);
                   Unfolding_RecoGenPt[k] = fs->make<TH2F>(name,name,nGenx[k],&xGen[k][0], nGenx[k],&xGen[k][0]);
                   Unfolding_RecoGenPt[k]->Sumw2();
               } // if(mIsMCarlo)
            } //for (int k=0; k<nEtabin; k++)

    //---------- Array For MC Weight PU -----------------//

 float Summer2012[600] = {0,0,0,0,0,3,3,1,2,3,2,4,2,1,1,4,6,6,7,8,6,7,6,8,8,10,11,13,4,12,9,13,13,6,9,37,58,50,48,51,45,47,48,52,56,102,82,101,118,95,108,97,112,98,83,273,275,252,272,266,269,284,267,241,258,1937,1907,1903,1989,2026,1973,1986,1966,1939,1969,5891,5986,6154,5988,6169,5949,5957,6007,6015,6004,10448,10329,10307,10547,10268,10307,10179,10278,10400,10311,13674,13331,13719,13645,13704,13547,13696,13862,13617,13741,16621,16528,16680,16772,16778,16638,16867,16723,16530,16507,20353,20254,20218,20348,20622,20212,20519,20530,20290,20410,25453,25546,25278,25440,25274,25667,25453,25311,25624,25673,32341,32562,32175,32250,32379,32673,32189,32152,32166,32138,40510,40639,40894,40730,40480,40630,40911,40621,40894,40461,49231,49589,49296,49156,49323,49440,49325,49313,49176,49067,54863,54835,54741,54882,54851,55019,54384,54602,54635,54751,56616,56096,56520,56937,56651,56939,56599,56675,56363,57077,55646,55470,55209,55770,55702,55862,55232,55965,55378,55334,52529,52526,52843,52714,52798,52667,52241,52206,52524,52798,49657,49510,49605,49411,49353,49385,49648,49660,49726,49542,47139,46760,47337,46690,47492,46819,47277,47333,47184,46952,44951,45061,45338,44962,45262,44916,45287,45253,44877,44921,43538,43336,43107,43052,43122,42700,43469,42982,43629,43142,41220,41052,41393,41199,41199,41193,41156,41046,41001,41367,38811,38923,38783,38851,38488,38808,38878,38841,38764,38559,36691,36497,36389,36056,36297,36221,36306,36285,36277,36201,33871,33580,33690,33749,33420,33562,33947,33496,33680,33736,31110,30668,30614,30695,30818,30720,30835,30631,30682,30722,27544,27587,27849,27735,27783,27593,27954,27770,28280,27590,24877,24821,24943,24879,25242,24766,24866,24763,25031,24598,22253,21803,21954,22038,21962,22010,21875,22147,21909,21974,19226,19147,19227,18986,19306,19210,19105,19349,19516,19221,16868,16603,16941,16628,16677,16725,16568,16784,16439,16760,14165,14308,14441,14170,14275,14332,14175,14431,14283,14311,12065,12035,12180,11760,12184,12003,11753,12075,12045,11893,10069,9947,10237,9997,10079,9977,10055,9883,10165,10078,8318,8392,8370,8250,8257,8375,8455,8419,8184,8396,6771,6784,6797,6698,6827,6651,6698,6871,6910,6763,5450,5489,5475,5456,5509,5523,5548,5512,5507,5410,4432,4348,4436,4590,4424,4366,4432,4331,4346,4361,3444,3555,3609,3509,3518,3453,3460,3502,3414,3549,2776,2693,2733,2675,2799,2800,2732,2794,2765,2748,2103,2139,2090,2184,2112,2069,2119,2093,2213,2160,1647,1716,1677,1640,1699,1590,1636,1661,1615,1655,1306,1182,1287,1209,1231,1243,1209,1261,1236,1198,954,947,944,953,961,969,976,952,914,912,693,739,741,672,691,710,680,763,698,692,505,548,473,538,519,509,495,577,542,527,371,404,401,393,399,335,387,403,402,411,280,294,275,305,283,275,286,297,291,266,191,198,203,232,191,209,199,198,186,222,148,142,136,131,142,133,148,159,153,113,116,104,114,102,100,96,112,89,90,111,65,58,58,67,60,72,53,60,62,65,57,46,41,47,70,39,51,50,61,45,27,34,31,22,37,47,41,43,28,34,16,24,17,26,24,15,27,18,24,19,21,14,13,8,15,8,13,18,10,9,7,4,6,3,7};


   mPuf= TFile::Open(mPuFileName.c_str());
   TH1F *dataPU;
   dataPU = (TH1F*)gDirectory-> Get(mPuTrigName.c_str());

/*   for(int i=20; i<65; i++)
    {
      dataPU->SetBinContent(i,0); 
   } */

   for(int i=0; i<600; i++)
     {
      WSummer2012.push_back(Summer2012[i]);
      WData2012.push_back(dataPU->GetBinContent(i+1));
   } // for(int i=0; i<60; i++)

 LumiReWeighting(WSummer2012,WData2012);


} // void NewTreeAnalyzer::beginJob()


void NewTreeAnalyzer::endJob()
 {
   mInf->Close();

} // void NewTreeAnalyzer::endJob()

//------------------------ PU Weight Computation ------------------//
double NewTreeAnalyzer::puweight(float npv)
 {
  int bin = weights_->GetXaxis()->FindBin( npv );
  return weights_->GetBinContent( bin );
 }

// --------------------- PU Reweighting Function -------------------- //
void NewTreeAnalyzer::LumiReWeighting( std::vector< float > MC_distr, std::vector< float > Lumi_distr)
  {
    // no histograms for input: use vectors
        // now, make histograms out of them:

        // first, check they are the same size...

        if( MC_distr.size() != Lumi_distr.size() ){

          std::cerr <<"ERROR: LumiReWeighting: input vectors have different sizes. Quitting... \n";
          return;

        }

        Int_t NBins = MC_distr.size();

 /*         MC_distr_ = fs->make<TH1F>("MC_distr","MC dist",NBins,-0.5, float(NBins)-0.5);
          Data_distr_ = fs->make<TH1F>("Data_distr","Data dist",NBins,-0.5, float(NBins)-0.5);

         weights_ = fs->make<TH1F>("luminumer","luminumer",NBins,-0.5, float(NBins)-0.5);
          den = fs->make<TH1F>("lumidenom","lumidenom",NBins,-0.5, float(NBins)-0.5); */

        MC_distr_ = fs->make<TH1F>("MC_distr","MC dist",NBins,-0.5, 59.5);
          Data_distr_ = fs->make<TH1F>("Data_distr","Data dist",NBins,-0.5, 59.5);

         weights_ = fs->make<TH1F>("luminumer","luminumer",NBins,-0.5, 59.5);
          den = fs->make<TH1F>("lumidenom","lumidenom",NBins,-0.5, 59.5);

        for(int ibin = 1; ibin<NBins+1; ++ibin ) {
          weights_->SetBinContent(ibin, Lumi_distr[ibin-1]);
          Data_distr_->SetBinContent(ibin, Lumi_distr[ibin-1]);
          den->SetBinContent(ibin,MC_distr[ibin-1]);
          MC_distr_->SetBinContent(ibin,MC_distr[ibin-1]);
        }

        // check integrals, make sure things are normalized
        float deltaH = weights_->Integral();
        if(fabs(1.0 - deltaH) > 0.02 ) { //*OOPS*...
          weights_->Scale( 1.0/ deltaH );
          Data_distr_->Scale( 1.0/ deltaH );
        }
        float deltaMC = den->Integral();
        if(fabs(1.0 - deltaMC) > 0.02 ) {
          den->Scale(1.0/ deltaMC );
          MC_distr_->Scale(1.0/ deltaMC );
        }

        weights_->Divide( den );  // so now the average weight should be 1.0    

        if(mIsMCarlo) {
        std::cout << " Lumi/Pileup Reweighting: Computed Weights per In-Time Nint " << std::endl;

        for(int ibin = 1; ibin<NBins+1; ++ibin){
          std::cout << "   " << ibin-1 << " " << weights_->GetBinContent(ibin) << std::endl;
        }
       } // if(mIsMCarlo)  

 } //void NewTreeAnalyzer::LumiReWeighting

 //--------------------------- analyze() fuction declaration ------------------ //
void NewTreeAnalyzer::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup)
 {
  unsigned NEntries = mTree->GetEntries();
  cout<<"Reading TREE: "<<NEntries<<" events"<<endl;

   int decade = 0 ;

  //    for(unsigned  l=0; l<NEVENTSN; l++) {
    for(unsigned  l=0; l<NEntries; l++) {

    //----------- progress report -------------
//    double progress = 10.0*l/(1.0*NEVENTSN);
    double progress = 10.0*l/(1.0*NEntries);
    int k = TMath::FloorNint(progress);
    if (k > decade)
      cout<<10*k<<" %"<<endl;
    decade = k;
    //----------- read the event --------------
    mTree->GetEntry(l);
    double wt1 = 1.0, wtpu = 1.0 ;
    if(mIsMCarlo) {
    wt1 = Event->evtHdr().weight();
    wtpu = puweight(Event->evtHdr().trpu());
    } // if(mIsMCarlo)
   double wt = wt1*wtpu;
//   cout<<"wt = "<<wt<<endl;
//   cout<<"PU  = "<<Event->evtHdr().trpu()<<endl;

// ------- for MC the trigger decision is always true ------ //
 if(mIsMCarlo) {
      for (int i=0; i<nJetTrig; i++)
         {
           for (int j=0; j<nJetVersn; j++)
              {
                   hltPassj[i][j] = true ;
                   prescalej[i][j] = 1;
            } // for (int j=0; j<nJetVersn; j++)
        } // for (int i=0; i<nJetTrig; i++)


 } // if(mIsMCarlo)

 if(!mIsMCarlo){
  //-------------------------- Initializing the Boolians and the Prescale Values----------------------//
        for (int i=0; i<nJetTrig; i++)
         {
           for (int j=0; j<nJetVersn; j++)
              {
                   hltPassj[i][j] = false ;
                   prescalej[i][j] = 1;
            } // for (int j=0; j<nJetVersn; j++)
        } // for (int i=0; i<nJetTrig; i++)


  //----------------------- Computng The Prescale Values For a given event ----------------------- //
    for (int i=0; i<nJetTrig; i++)
         {
           for (int j=0; j<nJetVersn; j++)
              {
                 if (ihltj[i][j] == -1)
                 hltPassj[i][j] = false; // no trigger set
                 else {
                 if (Event->fired(ihltj[i][j]) > 0) {
                 hltPassj[i][j] = true;
                 prescalej[i][j] = Event->preL1(ihltj[i][j]) * Event->preHLT(ihltj[i][j]);
        //         if(i==4)
        //         cout<<HLTJet[i][j]<<" has prescale = "<<prescalej[i][j]<<endl;
                 } // if (Event->fired(ihltj[i][j]) > 0)
               } // else
           } // for (int j=0; j<nJetVersn; j++)
      } // for (int i=0; i<nJetTrig; i++)

 } // if(!mIsMCarlo)


    //--------------------------- Filling Up The Histograms ------------------------------ //
      for (int i=0; i<nJetTrig-1; i++)
       {
        for (int j=0; j<nJetVersn; j++)
         {
          if (hltPassj[i][j])
           { 
             int NJets;
             if(mIsCHS) NJets = (int)Event->nPFJetsCHS();
             else       NJets = (int)Event->nPFJets();

             // ---- declaration of the vector of jets ---- //
             vector<QCDPFJet> mPFJets; mPFJets.clear();
    
            // ----- looping over the number of jets in the event ---- //        
             for(int iJet = 0; iJet < NJets; iJet++)
              {
               vector<double> JecFactors; JecFactors.clear(); 
               QCDPFJet pfjet;
               if(mIsCHS)
               pfjet = Event->pfjetchs(iJet); // ---- accessing the already corrected chs jet in the event -- //
               else
               pfjet = Event->pfjet(iJet); // ----- accessing the uncorrected jet in the event for non-chs jet----- //

               // ---- Old JEC factor ----- //
               double oldJecFactor;
               if(mIsCHS)          
               oldJecFactor = pfjet.cor(); // --- the chs jet 4-vector are corrected by old JEC factor ---- //
               else
               oldJecFactor = 1; // ---- the non-chs jet 4-vector are not corrected by old JEC factor ---- //

               LorentzVector oldJetP4 = pfjet.p4(); // ---- accessing the 4-vector of the jet ---- //
               TLorentzVector tmpJet;
               tmpJet.SetPxPyPzE(oldJetP4.px(), oldJetP4.py(), oldJetP4.pz(), oldJetP4.energy());

               TLorentzVector UnCorrectedJet = tmpJet * (1./oldJecFactor); // --- obtaining the uncorrected jet 4-vector by scaling down by the old jec factor --- //
                     
               // ---- Evaluating the L1Fast correction factor ---- //
               jecL1Fast->setJetPt(UnCorrectedJet.Pt()); 
               jecL1Fast->setJetA(pfjet.area());
               jecL1Fast->setRho(Event->evtHdr().pfRho());
               jecL1Fast->setJetEta(UnCorrectedJet.Eta());

               double corFactorL1Fast = jecL1Fast->getCorrection();
               //cout<<"L1Fast Cor Factor = "<<corFactorL1Fast<<endl;
               TLorentzVector tmpJetL1FastCorrected =
                              UnCorrectedJet * corFactorL1Fast; // ---- getting the jet corrected at L1Fast level ----- //

               // ---- Evaluating the L2Relative correction factor ---- //
               jecL2Relative->setJetPt(tmpJetL1FastCorrected.Pt());
               jecL2Relative->setJetEta(tmpJetL1FastCorrected.Eta());
               
               double corFactorL2Relative = jecL2Relative->getCorrection();
               //cout<<"L2Relative Cor Factor"<<corFactorL2Relative<<endl;
               TLorentzVector tmpJetL1FastL2RelativeCorrected =
                              tmpJetL1FastCorrected * corFactorL2Relative; //  ---- getting the jet corrected at L1Fast*L2Relative level ----- //

               // ---- Evaluating the L3Absolute correction factor ---- // 
               jecL3Absolute->setJetPt(tmpJetL1FastL2RelativeCorrected.Pt()); 
               jecL3Absolute->setJetEta(tmpJetL1FastL2RelativeCorrected.Eta());

               double corFactorL3Absolute = jecL3Absolute->getCorrection();
               //cout<<"L3Absolute Cor Factor"<<corFactorL3Absolute<<endl;
               TLorentzVector tmpJetL1FastL2RelativeL3AbsoluteCorrected =
                              tmpJetL1FastL2RelativeCorrected * corFactorL3Absolute; // -- getting the jet corrected at L1Fast*L2Relative*L3Absolute level -- //          

               // ---- Evaluating the L2L3Rsidual correction factor ---- //
               TLorentzVector tmpJetL1FastL2RelativeL3AbsoluteL2L3ResidualCorrected;
               double corFactorL2L3Residual(-1.);
               if(!mIsMCarlo) {
               jecL2L3Residual->setJetPt(tmpJetL1FastL2RelativeL3AbsoluteCorrected.Pt());
               jecL2L3Residual->setJetEta(tmpJetL1FastL2RelativeL3AbsoluteCorrected.Eta());

               corFactorL2L3Residual = jecL2L3Residual->getCorrection();
               //cout<<"L2L3Rsidual Cor Factor"<<corFactorL2L3Residual<<endl;
                              tmpJetL1FastL2RelativeL3AbsoluteL2L3ResidualCorrected = 
                              tmpJetL1FastL2RelativeL3AbsoluteCorrected * corFactorL2L3Residual; //  -- getting the jet corrected at L1Fast*L2Relative*L3Absolute*L2L3Residual level -- /
               } // if(!mIsMCarlo)  


               LorentzVector correctedJetP4; double CorFactor;
               if(!mIsMCarlo) { // -- if data, then take the full L1FastL2RelativeL3AbsoluteL2L3Residual corrected jet - //
               correctedJetP4 = LorentzVector(tmpJetL1FastL2RelativeL3AbsoluteL2L3ResidualCorrected.Px(),
                                              tmpJetL1FastL2RelativeL3AbsoluteL2L3ResidualCorrected.Py(),
                                              tmpJetL1FastL2RelativeL3AbsoluteL2L3ResidualCorrected.Pz(),
                                              tmpJetL1FastL2RelativeL3AbsoluteL2L3ResidualCorrected.E());

               CorFactor = tmpJetL1FastL2RelativeL3AbsoluteL2L3ResidualCorrected.E()/UnCorrectedJet.E(); 
               } 
               else { // -- if mc, the take the L1FastL2RelativeL3Absolute corrected jet --//
               correctedJetP4 = LorentzVector(tmpJetL1FastL2RelativeL3AbsoluteCorrected.Px(),
                                              tmpJetL1FastL2RelativeL3AbsoluteCorrected.Py(),
                                              tmpJetL1FastL2RelativeL3AbsoluteCorrected.Pz(),
                                              tmpJetL1FastL2RelativeL3AbsoluteCorrected.E());


               CorFactor = tmpJetL1FastL2RelativeL3AbsoluteCorrected.E()/UnCorrectedJet.E();
               }

               // ----- Storing the JEC correction factor in each label at a vector --- //
               JecFactors.push_back(corFactorL1Fast); JecFactors.push_back(corFactorL2Relative);
               JecFactors.push_back(corFactorL3Absolute); JecFactors.push_back(corFactorL2L3Residual);
               JecFactors.push_back(CorFactor);
               // --- declaring the new jet with corrected P4 and unaltered energy fractions ---//
               QCDPFJet newpfjet = pfjet; // -- initializing as the old corrected pfjet ----//
               newpfjet.setP4(correctedJetP4); // -- replacing the old P4 by newly corrected P4 -- //
               newpfjet.setCor(CorFactor);  // -- replacing the old corfactor by the new corfactor -- //
               newpfjet.setJecLabels(JecFactors); // ---- setting each label of JEC  ------ //
  
               if(newpfjet.ptCor() >= 20) 
               mPFJets.push_back(newpfjet);
             } // for(int iJet = 0; iJet < NJets; iJet++)  
     
             // ----- pT sorting the PF jets ------- //
             sort(mPFJets.begin(),mPFJets.end(),sort_pfjets);
       
 // ----- Initializing for JEC uncertainty sources -------------- //
    
    std::string mPFPayloadName, mPFJECUncSrc ;
    if(mIsCHS) {
    mPFPayloadName = mPFPayloadNameCHSA;
    mPFJECUncSrc   = mPFJECUncSrcCHSA; 
    }
    else {
    mPFPayloadName = mPFPayloadNameA;
    mPFJECUncSrc   = mPFJECUncSrcA;
    } // else
  
  //  cout<<"pfpayload name = "<<mPFPayloadName.c_str()<<", jec unc size = "<<mJECUncSrcNames.size()<<endl;
    
    edm::ESHandle<JetCorrectorParametersCollection> PFJetCorParColl; 
    if (mPFPayloadName != "" && !isPFJecUncSet_){
    iSetup.get<JetCorrectionsRecord>().get(mPFPayloadName,PFJetCorParColl);
    JetCorrectorParameters const& PFJetCorPar = (*PFJetCorParColl)["Uncertainty"];
     mPFUnc = new JetCorrectionUncertainty(PFJetCorPar);
     if (mPFJECUncSrc != "") {
       for(unsigned isrc=0;isrc<mJECUncSrcNames.size();isrc++) {
         JetCorrectorParameters *par = new JetCorrectorParameters(mPFJECUncSrc,mJECUncSrcNames[isrc]);
         JetCorrectionUncertainty *tmpUnc = new JetCorrectionUncertainty(*par);
         mPFUncSrc.push_back(tmpUnc);
       } // for(unsigned isrc=0;isrc<mJECUncSrcNames.size();isrc++)
     } // if (mPFJECUncSrc != "")
     isPFJecUncSet_ = true;
   } // if (mPFPayloadName != "" && !isPFJecUncSet_) 
    

          // ------- Emulating the higher trigger path from a lower one ------- //
            bool hltcut = false, l1cut = false;
            if(!mIsMCarlo) {
               //----------------- L1 Theshold Checking --------------------------//
               for ( unsigned l1iobj=0; l1iobj<Event->nL1Obj(ihltj[i][j]); l1iobj++ )
                 {
                   if (Event->l1obj(ihltj[i][j],l1iobj).pt() > ATrig[i+1])
                    {
                      l1cut = true ;
                    } // if (Event->l1obj(ihltj[i][j],l1iobj).pt() > ATrig[i+1]) 
                  } // for ( unsigned l1iobj=0; l1iobj<Event->nL1Obj(ihltj[i][j]); l1iobj++ )


               //------------ HLT Threshold Checking For HLTObj ------------------ //
               for ( unsigned hltiobj=0; hltiobj<Event->nHLTObj(ihltj[i][j]); hltiobj++ )
                 {
                  if (Event->hltobj(ihltj[i][j],hltiobj).pt() >  HLTJetPtN[i+1])
                   {
                     hltcut = true ;
                  } // if (Event->hltobj(ihltj[i][j],hltiobj).pt() >  HLTJetPtN[i+1])

                } // for ( unsigned hltiobj=0; hltiobj<Event->nHLTObj(ihltj[i][j]); hltiobj++ ) 

               } // if(!mIsMCarlo)
  

     // ---- starting the loop over newly corrected jets ----- //
     for(int iJet = 0; iJet < (int)mPFJets.size(); iJet++)
      {
       QCDPFJet thisJet = mPFJets[iJet];

        bool cutID =(thisJet.tightID() && (thisJet.elf() < 0.9) && (thisJet.muf() < 0.9) && (thisJet.nhf() < 0.9)  && (thisJet.phf() < 0.9) &&
        (Event->pfmet().met_o_sumet()<0.3) && (Event->evtHdr().pfRho()<100) && (Event->evtHdr().isPVgood()));

        if(cutID) {

       for (int k=0; k<nEtabin; k++)
        {
         if ((fabs(thisJet.y())>= EtabinN[k]) && (fabs(thisJet.y())<EtabinN[k+1])) {

               if(!mIsMCarlo) {
               for (int r=0; r<nvtx; r++)
                    {
                    if((Event->evtHdr().nVtxGood()>VertexN[r])&&(Event->evtHdr().nVtxGood()<=VertexN[r+1])) {
                       PFJet0[i][j][k][r] -> Fill(thisJet.ptCor(),wt);
                   //    cout<<"weight = "<<wt<<endl;
                       if(l1cut && hltcut) {
                       PFJet1[i][j][k][r] -> Fill(thisJet.ptCor(),wt);
                       } // if(l1cut && hltcut)

                  //    cout<<"Uncertainty factor = "<<Event->pfjet(p).uncSrc(16)<<endl;
                    } // if((Event->evtHdr().nVtxGood()>VertexN[r])&&(Event->evtHdr().nVtxGood()<=VertexN[r+1]))
                   } // for (int r=0; r<nvtx; r++)
                  } // if(!mIsMCarlo)

                   if(thisJet.ptCor()>HLTJetPtT[i]) {
                   //------------------------- Fill Of JetID variable -------------------//
                    ChHadFr[i][j][k] ->Fill(thisJet.chf(),wt);
                    NuHadFr[i][j][k] ->Fill(thisJet.nhf(),wt);
                    PhFr[i][j][k] ->Fill(thisJet.phf(),wt);
                    ElFr[i][j][k] ->Fill(thisJet.elf(),wt);
                    MuFr[i][j][k] ->Fill(thisJet.muf(),wt);
                   //----------------------- Fill For Data-MC --------------------------//
                   //PFJet2[i][j][k] -> Fill(Event->pfjet(p).ptCor(),wt);
                   //PFJet2Eta[i][j][k] -> Fill(Event->pfjet(p).eta(),wt);
                   //PFJet2Phi[i][j][k] -> Fill(Event->pfjet(p).phi(),wt);
                   } // if(Event->pfjet(p).ptCor()>HLTJetPtT[i])

                  //-------------- Filling Histogram For Spectrum --------------------//
                  if(thisJet.ptCor()>=HLTJetPtS[i] && thisJet.ptCor()<HLTJetPtS[i+1]) {
                   PFJet[i][j][k] -> Fill(thisJet.ptCor(),prescalej[i][j]);
                   PFJetYield[i][j][k] -> Fill(thisJet.ptCor(),wt1);
                   // --- Filling the JEC factors ----- //
                   PFJetPt_L1Fast[i][j][k] -> Fill(thisJet.ptCor(),thisJet.jecLabels(0));
                   PFJetPt_L2Relative[i][j][k] -> Fill(thisJet.ptCor(),thisJet.jecLabels(1));
                   PFJetPt_L3Absolute[i][j][k] -> Fill(thisJet.ptCor(),thisJet.jecLabels(2));
                   PFJetPt_L2L3Residual[i][j][k] -> Fill(thisJet.ptCor(),thisJet.jecLabels(3));
                   PFJetPt_TotalJEC[i][j][k] -> Fill(thisJet.ptCor(),thisJet.jecLabels(4)); 
 
                   if(!mIsMCarlo)
                   BetavsPT[i][j][k]->Fill(thisJet.ptCor(),thisJet.beta(),prescalej[i][j]);
                   } // if(Event->pfjet(p).ptCor()>=HLTJetPtS[i] && Event->pfjet(p).ptCor()<HLTJetPtS[i+1])    

                  if(mIsMCarlo){
                   if(thisJet.genpt()>=HLTJetPtS[i] && thisJet.genpt()<HLTJetPtS[i+1]) {
                   PFJetGen[i][j][k] -> Fill(thisJet.genpt(),wt1);
                   }
                  }

                  if(!mIsMCarlo){
                  for(int r=0; r<nPU; r++)
                   {
                    if((Event->evtHdr().nVtxGood()>=PUBin[r]) && (Event->evtHdr().nVtxGood()<PUBin[r+1])) {
                    if(thisJet.ptCor()>=HLTJetPtS[i] && thisJet.ptCor()<HLTJetPtS[i+1]) {
                     PFJetPU[i][j][k][r] -> Fill(thisJet.ptCor(),prescalej[i][j]);
                    } // if(Event->pfjet(p).ptCor()>=HLTJetPtS[i] && Event->pfjet(p).ptCor()<HLTJetPtS[i+1])  
                   } // if((Event->evtHdr().nVtxGood()>=PUBin[r]) && (Event->evtHdr().nVtxGood()<PUBin[r+1]))
                  } // for(int r=0; r<nPU; r++)  
                  } // if(!mIsMCarlo)


                  //--------------Filling The Profile Plots For Stability Check ------//
                  if(!mIsMCarlo) {

                  if(thisJet.ptCor()>=HLTJetPtS[i] && thisJet.ptCor()<HLTJetPtS[i+1]) {
                  PFPt[k][i][j] -> Fill(thisJet.ptCor(),wt);
                  }


                 //---- jec uncertainty --------------
                 double unc1(0.0);
                 vector<float> uncSrc(0);
                 if (mPFPayloadName != "") {
                 mPFUnc->setJetEta(thisJet.eta());
                 mPFUnc->setJetPt(thisJet.ptCor());
                 unc1 = mPFUnc->getUncertainty(true);
                 }

                  if (mPFJECUncSrc != "") {
                  for(unsigned isrc=0;isrc<mJECUncSrcNames.size();isrc++)
                    {  
                     mPFUncSrc[isrc]->setJetEta(thisJet.eta());
                     mPFUncSrc[isrc]->setJetPt(thisJet.ptCor());
                     float unc = mPFUncSrc[isrc]->getUncertainty(true);
                     uncSrc.push_back(unc);
                     double ptUp = (1+unc)*thisJet.ptCor();
                     double ptDown = (1-unc)*thisJet.ptCor();
           
                     if(ptUp>=HLTJetPtS[i] && ptUp<HLTJetPtS[i+1]) {
                     PFPtUp[k][isrc][i][j] ->Fill(ptUp,wt);
                     } // if(ptUp>=HLTJetPtS[i] && ptUp<HLTJetPtS[i+1])

                     if(ptDown>=HLTJetPtS[i] && ptDown<HLTJetPtS[i+1]) {
                     PFPtDown[k][isrc][i][j] ->Fill(ptDown,wt);
                      } 
                   } // for(unsigned isrc=0;isrc<mJECUncSrcNames.size();isrc++) 
                  } // if (mPFJECUncSrc != "") 
              
                  thisJet.setUnc(unc1);
                  thisJet.setUncSrc(uncSrc);

                  } // if(!mIsMCarlo)


        } // if ((fabs(thisJet.y())>= EtabinN[k]) && (fabs(thisJet.y())<EtabinN[k+1]))  
       } // for (int k=0; k<nEtabin; k++) 
      } // if(cutID)

     } // for(int iJet = 0; iJet < (int)mPFJets.size(); iJet++)

       
          if(!mIsMCarlo && mPFJets.size()>=2){
        //----------------------- JetID Efficiency Tag and Probe ---------------------//
        int irand = (gRandom->Uniform()>0.5) ? 0 : 1;
        int irand2 = (irand+1)%2;

        QCDPFJet jet1 = mPFJets[irand] ;
        QCDPFJet jet2 = mPFJets[irand2];

       bool cutID2=((jet1.elf() < 0.9) && (jet1.muf() < 0.9) && (jet1.nhf() < 0.9) && (jet1.phf() < 0.9) &&  (Event->pfmet().met_o_sumet()<0.3) );

      double dphi = DeltaPhi(jet1.phi(),jet2.phi());

        if(cutID2 && (dphi>2.7)) {
            
        for (int k=0; k<nEtabin; k++)
           {
            if ((fabs(jet2.y())>= EtabinN[k]) && (fabs(jet2.y())<EtabinN[k+1]))  {
                if(jet1.tightID() && (fabs(jet1.y())<1.3)) {

             if(jet2.ptCor()>=HLTJetPtS[i] && jet2.ptCor()<HLTJetPtS[i+1]) {
                 TagJet[i][j][k]->Fill(jet2.ptCor(),prescalej[i][j]);

                  if(jet2.tightID()  &&  (jet2.elf() < 0.9) && (jet2.muf() < 0.9) && (jet2.nhf() < 0.9) && (jet2.phf() < 0.9)) {
                  ProbeJet[i][j][k]->Fill(jet2.ptCor(),prescalej[i][j]);
                   } // if(Event->pfjet(irand2).tightID())

                 } // if(Event->pfjet(irand2).ptCor()>=HLTJetPtS[i] && Event->pfjet(irand2).ptCor()<HLTJetPtS[i+1])
                } // if(Event->pfjet(irand).tightID() && (fabs(Event->pfjet(irand).y())<1.3)) 

              } // if ((fabs(Event->pfjet(irand2).y())>= EtabinN[k]) && (fabs(Event->pfjet(irand2).y())<EtabinN[k+1]))
           } // for (int k=0; k<nEtabin; k++) 
        } // if(cutID2 && (dphi>2.7)) 
      } // if(!mIsMCarlo)



      //---------------- Jet Resolution for MC ---------------//
     if(i==0 && j==0) { // --- run the tag and probe part only once in MC, i.e. 1-trigger-1-version --- //
     if(mIsMCarlo) {  

     for (int iJet = 0; iJet < (int)mPFJets.size(); iJet++)
       {

        QCDPFJet thisJet = mPFJets[iJet];

  bool cutID=(thisJet.tightID() && (thisJet.elf() < 0.9) && (thisJet.muf() < 0.9) && (thisJet.nhf() < 0.9) && (thisJet.phf() < 0.9) && (Event->pfmet().met_o_sumet()<0.3));

      if(cutID) {
      for (int k=0; k<nEtabin; k++)
         {
          if ((fabs(thisJet.y())>= EtabinN[k]) && (fabs(thisJet.y())<EtabinN[k+1])) {
          for(int m=0; m<n1x[k]; m++)
            {
             if((thisJet.genpt()>x1[k][m]) && (thisJet.genpt()<x1[k][m+1]) && (thisJet.genR() < 0.35)) {
              PtRes[k][m] ->Fill(thisJet.ptCor()/thisJet.genpt());
              AvgX_Pt[k][m] ->Fill(m,thisJet.genpt());
           } // if((Event->pfjet(p).genpt()>Ptbin1[m]) && (Event->pfjet(p).genpt()<Ptbin1[m+1]))
         } // for(int m=0; m<nbin1; m++)

          Unfolding_GenPt[k]->Fill(thisJet.genpt());
          Unfolding_RecoPt[k]->Fill(thisJet.ptCor());
          Unfolding_RecoGenPt[k]->Fill(thisJet.ptCor(), thisJet.genpt());

        } // if ((fabs(Event->pfjet(p).y())>= EtabinN[k]) && (fabs(Event->pfjet(p).y())<EtabinN[k+1]))
       } // for (int k=0; k<nEtabin; k++)
      } // if(cutID)
     } // for (unsigned p=0;p<Event->nPFJets();p++)    

 // ------------------------------- Unbiased GenJet Pt Distribution ------------------// 
     for (unsigned p=0;p<Event->nGenJets();p++)
        {
         for (int k=0; k<nEtabin; k++)
          {
            if ((fabs(Event->genjet(p).y())>= EtabinN[k]) && (fabs(Event->genjet(p).y())<EtabinN[k+1])) {
              GenJet[k]->Fill(Event->genjet(p).pt());
           }  // if ((fabs(Event->genjet(p).y())>= EtabinN[k]) && (fabs(Event->genjet(p).y())<EtabinN[k+1]))
        } // for (int k=0; k<nEtabin; k++)
     } // for (unsigned p=0;p<Event->nGenJets();p++)

  } // if(mIsMCarlo)
   
 } // if(i==0 && j==0)
 
 
    } //if (hltPassj[i][j])
   } // for (int j=0; j<nJetVersn; j++)
  } // for (int i=0; i<nJetTrig-1; i++) 



  

 } // for(unsigned  l=0; l<NEntries; l++)

} // void NewTreeAnalyzer::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup)

// ----- Declaring the destructor of the class ---- //
NewTreeAnalyzer::~NewTreeAnalyzer()
{
}

DEFINE_FWK_MODULE(NewTreeAnalyzer);

