#ifndef JECs_h
#define JECs_h
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SMPJ/AnalysisFW/interface/QCDJet.h"

#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

using namespace edm;
using namespace std;



class JECs  {

  typedef reco::Particle::LorentzVector LorentzVector;

  public:
  JECs(bool IsMCarlo, string GlobalTag, string JETTYPE, string jecUncSrc, vector<string> jecUncSrcNames);

    static bool sort_pfjets(QCDPFJet j1, QCDPFJet j2) {
        return j1.ptCor() > j2.ptCor();
    }

    static bool sort_calojets(QCDCaloJet j1, QCDCaloJet j2) {
        return j1.ptCor() > j2.ptCor();
    }

    JetCorrectorParameters *L1Fast, *L2Relative, *L3Absolute, *L2L3Residual;
    vector<JetCorrectorParameters> vecL1Fast, vecL2Relative, vecL3Absolute, vecL2L3Residual;
    FactorizedJetCorrector *jecL1Fast, *jecL2Relative, *jecL3Absolute, *jecL2L3Residual;
    JetCorrectionUncertainty *mPFUnc;
    JetCorrectionUncertainty *mCaloUnc;

    JetCorrectorParameters *par;
    JetCorrectionUncertainty *tmpUnc;
    std::vector<JetCorrectionUncertainty*> mPFUncSrc;
    std::vector<JetCorrectionUncertainty*> mCaloUncSrc;
    std::vector<string> mJECUncSrcNames;
    
    virtual void JEC_corrections(QCDEvent *Event, unsigned n_PFJets, bool IsMCarlo, vector<string> jecUncSrcNames){


             // ---- declaration of the vector of jets ---- //
             vector<QCDPFJet> mPFJets; mPFJets.clear();

            // ----- looping over the number of jets in the event ---- //
             for(unsigned iJet = 0; iJet < n_PFJets; iJet++)
              {
               vector<double> JecFactors; JecFactors.clear();
               QCDPFJet pfjet;
               pfjet = Event->pfjet(iJet); // ----- accessing the uncorrected jet in the event for non-chs jet----- //

               // ---- Old JEC factor ----- //
               double oldJecFactor = pfjet.cor();

               LorentzVector oldJetP4= pfjet.p4(); // ---- accessing the 4-vector of the jet ---- //
               TLorentzVector tmpJet;
               tmpJet.SetPxPyPzE(oldJetP4.px(), oldJetP4.py(), oldJetP4.pz(), oldJetP4.energy());

               TLorentzVector UnCorrectedJet = tmpJet * (1./oldJecFactor); // --- obtaining the uncorrected jet 4-vector by scaling down by the old jec factor --- //
	       //TLorentzVector UnCorrectedJet = tmpJet; ///Jets are already uncorrected
               // ---- Evaluating the L1Fast correction factor ---- //
               jecL1Fast->setJetPt(UnCorrectedJet.Pt());
               jecL1Fast->setJetA(pfjet.area());
               jecL1Fast->setRho(Event->evtHdr().pfRho());
               jecL1Fast->setJetEta(UnCorrectedJet.Eta());
	       //std::cout<<std::endl<<pfjet.area()<<std::endl;

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
               if(!IsMCarlo) {
               jecL2L3Residual->setJetPt(tmpJetL1FastL2RelativeL3AbsoluteCorrected.Pt());
               jecL2L3Residual->setJetEta(tmpJetL1FastL2RelativeL3AbsoluteCorrected.Eta());

               corFactorL2L3Residual = jecL2L3Residual->getCorrection();
               //cout<<"L2L3Rsidual Cor Factor"<<corFactorL2L3Residual<<endl;
                              tmpJetL1FastL2RelativeL3AbsoluteL2L3ResidualCorrected =
                              tmpJetL1FastL2RelativeL3AbsoluteCorrected * corFactorL2L3Residual; //  -- getting the jet corrected at L1Fast*L2Relative*L3Absolute*L2L3Residual level -- /
               } // if(!IsMCarlo)


               LorentzVector correctedJetP4; double CorFactor;
               if(!IsMCarlo) { // -- if data, then take the full L1FastL2RelativeL3AbsoluteL2L3Residual corrected jet - //
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
               pfjet.setP4(correctedJetP4); // -- replacing the old P4 by newly corrected P4 -- //

	       //JEC Uncertainties
	       double unc(0.0);
	       vector<float> uncSrc(0);
	       mPFUnc->setJetEta(pfjet.eta());
	       mPFUnc->setJetPt(pfjet.pt());
	       unc = mPFUnc->getUncertainty(true);
	       
	       for(unsigned isrc=0;isrc<jecUncSrcNames.size();isrc++) {
		 mPFUncSrc[isrc]->setJetEta(pfjet.eta());
	       	 mPFUncSrc[isrc]->setJetPt(pfjet.pt());
	       	 float unc1 = mPFUncSrc[isrc]->getUncertainty(true);
	       	 uncSrc.push_back(unc1);
	       } // for(unsigned isrc=0;isrc<mPFJECUncSrcNames.size();isrc++)
	       ///end of JEC Uncertainties
	       
	       pfjet.setCor(CorFactor);  // -- replacing the old corfactor by the new corfactor -- //
	       pfjet.setJecLabels(JecFactors); // ---- setting each label of JEC  ------ //
	       pfjet.setUnc(unc);
	       pfjet.setUncSrc(uncSrc);
	       mPFJets.push_back(pfjet);
	       
        } // for(int iJet = 0; iJet < NJets; iJet++)

    // ----- pT sorting the PF jets ------- //
    sort(mPFJets.begin(),mPFJets.end(),sort_pfjets);
    Event->setPFJets(mPFJets);
    }

    virtual void JEC_CHScorrections(QCDEvent *Event, unsigned n_PFJetsCHS, bool IsMCarlo, vector<string> jecUncSrcNames){

             // ---- declaration of the vector of jets ---- //
             vector<QCDPFJet> mPFJetsCHS; mPFJetsCHS.clear();

            // ----- looping over the number of jets in the event ---- //
             for(unsigned iJet = 0; iJet < n_PFJetsCHS; iJet++)
              {
               vector<double> JecFactors; JecFactors.clear();
               QCDPFJet pfjetchs;
               pfjetchs = Event->pfjetchs(iJet); // ----- accessing the uncorrected jet in the event for non-chs jet----- //

               // ---- Old JEC factor ----- //
               double oldJecFactor = pfjetchs.cor();

               LorentzVector oldJetP4 = pfjetchs.p4(); // ---- accessing the 4-vector of the jet ---- //
               TLorentzVector tmpJet;
               tmpJet.SetPxPyPzE(oldJetP4.px(), oldJetP4.py(), oldJetP4.pz(), oldJetP4.energy());

               TLorentzVector UnCorrectedJet = tmpJet * (1./oldJecFactor); // --- obtaining the uncorrected jet 4-vector by scaling down by the old jec factor --- //
	       //TLorentzVector UnCorrectedJet = tmpJet; ///Jets are already uncorrected
               // ---- Evaluating the L1Fast correction factor ---- //
               jecL1Fast->setJetPt(UnCorrectedJet.Pt());
               jecL1Fast->setJetA(pfjetchs.area());
               jecL1Fast->setRho(Event->evtHdr().pfRho());
               jecL1Fast->setJetEta(UnCorrectedJet.Eta());
	       //std::cout<<std::endl<<pfjetchs.area()<<std::endl;

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
	       //TLorentzVector tmpJetL1FastL2RelativeL3AbsoluteCorrected = UnCorrectedJet;

               TLorentzVector tmpJetL1FastL2RelativeL3AbsoluteL2L3ResidualCorrected;
               double corFactorL2L3Residual(-1.);
               if(!IsMCarlo) {
		 jecL2L3Residual->setJetPt(tmpJetL1FastL2RelativeL3AbsoluteCorrected.Pt());
		 jecL2L3Residual->setJetEta(tmpJetL1FastL2RelativeL3AbsoluteCorrected.Eta());

		 corFactorL2L3Residual = jecL2L3Residual->getCorrection();
               //cout<<"L2L3Rsidual Cor Factor"<<corFactorL2L3Residual<<endl;
		 tmpJetL1FastL2RelativeL3AbsoluteL2L3ResidualCorrected =
		   tmpJetL1FastL2RelativeL3AbsoluteCorrected * corFactorL2L3Residual; //  -- getting the jet corrected at L1Fast*L2Relative*L3Absolute*L2L3Residual level -- /
               } // if(!IsMCarlo)


               LorentzVector correctedJetP4; double CorFactor;
               if(!IsMCarlo) { // -- if data, then take the full L1FastL2RelativeL3AbsoluteL2L3Residual corrected jet - //
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
               pfjetchs.setP4(correctedJetP4); // -- replacing the old P4 by newly corrected P4 -- //

                ///JEC Uncertainties
	       double unc(0.0);
	       vector<float> uncSrc(0);
	       mPFUnc->setJetEta(pfjetchs.eta());
	       mPFUnc->setJetPt(pfjetchs.pt());
	       unc = mPFUnc->getUncertainty(true);

//
	       for(unsigned isrc=0;isrc<jecUncSrcNames.size();isrc++) {
		 //for(unsigned isrc=0;isrc<2;isrc++) {
		 mPFUncSrc[isrc]->setJetEta(pfjetchs.eta());
		 mPFUncSrc[isrc]->setJetPt(pfjetchs.pt());
		 float unc1 = mPFUncSrc[isrc]->getUncertainty(true);
		 uncSrc.push_back(unc1);
	       } // for(unsigned isrc=0;isrc<mPFJECUncSrcNames.size();isrc++)
	       ///end of JEC Uncertainties

	       pfjetchs.setCor(CorFactor);  // -- replacing the old corfactor by the new corfactor -- //
	       pfjetchs.setJecLabels(JecFactors); // ---- setting each label of JEC  ------ //
	       pfjetchs.setUnc(unc);
               pfjetchs.setUncSrc(uncSrc);
               mPFJetsCHS.push_back(pfjetchs);

        } // for(int iJet = 0; iJet < NJets; iJet++)

    // ----- pT sorting the PF jets ------- //
    sort(mPFJetsCHS.begin(),mPFJetsCHS.end(),sort_pfjets);
    Event->setPFJetsCHS(mPFJetsCHS);
    }

    virtual void JEC_Calocorrections(QCDEvent *Event, unsigned n_CaloJets, bool IsMCarlo){

             // ---- declaration of the vector of jets ---- //
             vector<QCDCaloJet> mCaloJets; mCaloJets.clear();

            // ----- looping over the number of jets in the event ---- //
             for(unsigned iJet = 0; iJet < n_CaloJets; iJet++)
              {
               vector<double> JecFactors; JecFactors.clear();
               QCDCaloJet calojet;
               calojet = Event->calojet(iJet); // ----- accessing the uncorrected jet in the event for non-chs jet----- //

               // ---- Old JEC factor ----- //
               //double oldJecFactor = calojet.cor();

               LorentzVector oldJetP4 = calojet.p4(); // ---- accessing the 4-vector of the jet ---- //
               TLorentzVector tmpJet;
               tmpJet.SetPxPyPzE(oldJetP4.px(), oldJetP4.py(), oldJetP4.pz(), oldJetP4.energy());

               //TLorentzVector UnCorrectedJet = tmpJet * (1./oldJecFactor); // --- obtaining the uncorrected jet 4-vector by scaling down by the old jec factor --- //
                TLorentzVector UnCorrectedJet = tmpJet; ///Jets are already uncorrected
               // ---- Evaluating the L1Fast correction factor ---- //
               jecL1Fast->setJetPt(UnCorrectedJet.Pt());
               jecL1Fast->setJetA(calojet.area());
               jecL1Fast->setRho(Event->evtHdr().caloRho());
               jecL1Fast->setJetEta(UnCorrectedJet.Eta());
	       //std::cout<<std::endl<<calojet.area()<<std::endl;

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
               if(!IsMCarlo) {
               jecL2L3Residual->setJetPt(tmpJetL1FastL2RelativeL3AbsoluteCorrected.Pt());
               jecL2L3Residual->setJetEta(tmpJetL1FastL2RelativeL3AbsoluteCorrected.Eta());

               corFactorL2L3Residual = jecL2L3Residual->getCorrection();
               //cout<<"L2L3Rsidual Cor Factor"<<corFactorL2L3Residual<<endl;
                              tmpJetL1FastL2RelativeL3AbsoluteL2L3ResidualCorrected =
                              tmpJetL1FastL2RelativeL3AbsoluteCorrected * corFactorL2L3Residual; //  -- getting the jet corrected at L1Fast*L2Relative*L3Absolute*L2L3Residual level -- /
               } // if(!IsMCarlo)


               LorentzVector correctedJetP4; double CorFactor;
               if(!IsMCarlo) { // -- if data, then take the full L1FastL2RelativeL3AbsoluteL2L3Residual corrected jet - //
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
               calojet.setP4(correctedJetP4); // -- replacing the old P4 by newly corrected P4 -- //

                ///JEC Uncertainties
//                double unc(0.0);
//                vector<float> uncSrc(0);
//                mCaloUnc->setJetEta(calojet.eta());
//                mCaloUnc->setJetPt(calojet.pt());
//                unc = mCaloUnc->getUncertainty(true);
//
//                for(unsigned isrc=0;isrc<mJECUncSrcNames.size();isrc++) {
//                    mCaloUncSrc[isrc]->setJetEta(calojet.eta());
//                    mCaloUncSrc[isrc]->setJetPt(calojet.pt());
//                    float unc1 = mCaloUncSrc[isrc]->getUncertainty(true);
//                    uncSrc.push_back(unc1);
//                } // for(unsigned isrc=0;isrc<mCaloJECUncSrcNames.size();isrc++)
                ///end of JEC Uncertainties

                calojet.setCor(CorFactor);  // -- replacing the old corfactor by the new corfactor -- //
                calojet.setJecLabels(JecFactors); // ---- setting each label of JEC  ------ //
                //calojet.setUnc(unc);
                //calojet.setUncSrc(uncSrc);
                mCaloJets.push_back(calojet);

        } // for(int iJet = 0; iJet < NJets; iJet++)

    // ----- pT sorting the CALO jets ------- //
    sort(mCaloJets.begin(),mCaloJets.end(),sort_calojets);
    Event->setCaloJets(mCaloJets);
    }
};

JECs::JECs(bool IsMCarlo, string GlobalTag, string JETTYPE, string jecUncSrc, vector<string> mJECUncSrcNames){

  string file_data_mc = "DATA";
  if(IsMCarlo)    file_data_mc = "MC";

     L1Fast       = new JetCorrectorParameters("../data/"+GlobalTag+"/"+GlobalTag+"_"+file_data_mc+"_L1FastJet_"+JETTYPE+".txt");
     L2Relative   = new JetCorrectorParameters("../data/"+GlobalTag+"/"+GlobalTag+"_"+file_data_mc+"_L2Relative_"+JETTYPE+".txt");
     L3Absolute   = new JetCorrectorParameters("../data/"+GlobalTag+"/"+GlobalTag+"_"+file_data_mc+"_L3Absolute_"+JETTYPE+".txt");
     if(!IsMCarlo)
        L2L3Residual = new JetCorrectorParameters("../data/"+GlobalTag+"/"+GlobalTag+"_DATA_L2L3Residual_"+JETTYPE+".txt");

     vecL1Fast.push_back(*L1Fast);
     vecL2Relative.push_back(*L2Relative);
     vecL3Absolute.push_back(*L3Absolute);
     if(!IsMCarlo)
        vecL2L3Residual.push_back(*L2L3Residual);


     jecL1Fast       = new FactorizedJetCorrector(vecL1Fast);
     jecL2Relative   = new FactorizedJetCorrector(vecL2Relative);
     jecL3Absolute   = new FactorizedJetCorrector(vecL3Absolute);
     if(!IsMCarlo)
       jecL2L3Residual = new FactorizedJetCorrector(vecL2L3Residual);

     ///Read Uncertainty txt files
     mPFUnc = new JetCorrectionUncertainty("../data/"+GlobalTag+"/"+GlobalTag+"_"+file_data_mc+"_Uncertainty_"+JETTYPE+".txt");    

     for(unsigned isrc=0;isrc<mJECUncSrcNames.size();isrc++) {
       cout<<mJECUncSrcNames[isrc]<<endl;

       par = new JetCorrectorParameters(jecUncSrc,mJECUncSrcNames[isrc]);
       tmpUnc = new JetCorrectionUncertainty(*par);

       mPFUncSrc.push_back(tmpUnc);
     }
}



#endif

