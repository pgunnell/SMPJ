#include "SMPJ/AnalysisFW/interface/QCDEvent.h"
//---------------------------------------------------
QCDEvent::QCDEvent()
{
  PFJetsCHS_.clear();
  GenJets_.clear();
}
//---------------------------------------------------
QCDEvent::~QCDEvent()
{
}
//---------------------------------------------------
void QCDEvent::setPFJetsCHS(const std::vector<QCDPFJet>& fPFJetsCHS)
{
  PFJetsCHS_.clear();
  for(unsigned i=0;i<fPFJetsCHS.size();i++) {
    PFJetsCHS_.push_back(fPFJetsCHS[i]);
  }
}
//---------------------------------------------------
void QCDEvent::setGenJets(const std::vector<LorentzVector>& fGenJets)
{
  GenJets_.clear();
  for(unsigned i=0;i<fGenJets.size();i++) {
    GenJets_.push_back(fGenJets[i]);
  }
}
//---------------------------------------------------
void QCDEvent::setHLTObj(const std::vector<std::vector<LorentzVector> >& fHLTObj)
{
  HLTObj_.clear();
  for(unsigned i=0;i<fHLTObj.size();i++) {
    std::vector<LorentzVector> vv;
    for(unsigned j=0;j<fHLTObj[i].size();j++) {
      vv.push_back(fHLTObj[i][j]);
    }
    HLTObj_.push_back(vv);
  }
}
//---------------------------------------------------
