#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/TauReco/interface/Tau.h"
#include "DataFormats/TauReco/interface/TauDiscriminatorByIsolation.h"

#include "Math/GenVector/VectorUtil.h"
#include "Math/GenVector/PxPyPzE4D.h"

#include <memory>
#include <string>
#include <iostream>

#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1.h>

using namespace edm;
using namespace reco; 
using namespace std;

class PFTauTest : public EDAnalyzer {
public:
  explicit PFTauTest(const ParameterSet&);
  ~PFTauTest() {}
  virtual void analyze(const Event& iEvent,const EventSetup& iSetup);
  virtual void beginJob();
  virtual void endJob();
private:
  string PFTaus_;
  string TauDiscriminatorByIsolationProducer_;
  int nEvent;
  vector<float> nEventsUsed;
  vector<float> nEventsRiso;
  int nEventTaggedJets;
};

PFTauTest::PFTauTest(const ParameterSet& iConfig){
  PFTaus_                              = iConfig.getParameter<string>("PFTaus");
  TauDiscriminatorByIsolationProducer_ = iConfig.getParameter<string>("TauDiscriminatorByIsolationProducer");
  nEvent=0;
  nEventTaggedJets=0;
  nEventsRiso.reserve(6);
  nEventsUsed.reserve(6);
  for(int i=0;i<6;i++){
    nEventsRiso[i]=0.;
    nEventsUsed[i]=0.;
  }
}

void PFTauTest::beginJob(){}

void PFTauTest::analyze(const Event& iEvent, const EventSetup& iSetup){
  cout<<"********"<<endl;
  cout<<"Event number "<<nEvent++<<endl;
  
  Handle<TauCollection> tauHandle;
  iEvent.getByLabel(PFTaus_,tauHandle);
  const TauCollection& myTauCollection=*(tauHandle.product()); 

  Handle<TauDiscriminatorByIsolation> theTauDiscriminatorByIsolation;
  iEvent.getByLabel(TauDiscriminatorByIsolationProducer_,theTauDiscriminatorByIsolation);

  cout<<"***"<<endl;
  cout<<"Found "<<myTauCollection.size()<<" had. tau-jet candidates"<<endl;
  int i_Tau=0;
  for (TauCollection::size_type iTau=0;iTau<tauHandle->size();iTau++) {
    TauRef theTau(tauHandle,iTau);
    //Prints out some quantities
    cout<<"Jet Number "<<i_Tau<<endl;
    cout<<"DiscriminatorByIsolation value "<<(*theTauDiscriminatorByIsolation)[theTau]<<endl;
    cout<<"Pt of the Tau "<<(*theTau).pt()<<endl;
    PFCandidateRef theLeadPFCand = (*theTau).getleadPFChargedHadrCand();
    if(!theLeadPFCand){
      cout<<"No Lead PFCand "<<endl;
    }else{
      cout<<"Lead PFCand pt "<<(*theLeadPFCand).pt()<<endl;
      cout<<"InvariantMass of the Tau "<<(*theTau).getInvariantMass()<<endl;
      cout<<"Vertex of the Tau "<<(*theTau).vz()<<endl;
      cout<<"Charge of the Tau "<<(*theTau).charge()<<endl;
      cout<<"Em energy fraction "<<(*theTau).getEmEnergyFraction()<<endl;
      cout<<"Max Hadron energy "<<(*theTau).getMaximumHcalEnergy()<<endl;
      cout<<"# PF charged hadr. cand's "<<(*theTau).getSelectedPFChargedHadrCands().size()<<endl;
      cout<<"# PF neutral hadr. cand's "<<(*theTau).getSelectedPFNeutrHadrCands().size()<<endl;
      cout<<"# PF gamma cand's "<<(*theTau).getSelectedPFGammaCands().size()<<endl;
      cout<<"Number of SignalPFChargedHadrCands = "<<(*theTau).getSignalPFChargedHadrCands().size()<<endl;
      cout<<"Number of IsolationPFChargedHadrCands = "<<(*theTau).getIsolationPFChargedHadrCands().size()<<endl;
      cout<<"Number of SignalPFGammaCands = "<<(*theTau).getSignalPFGammaCands().size()<<endl;
      cout<<"Number of IsolationPFGammaCands = "<<(*theTau).getIsolationPFGammaCands().size()<<endl;
      cout<<"Sum pT of IsolationPFChargedHadrCands  = "<<(*theTau).getSumPtIsolation()<<endl;
      cout<<"Sum ET of IsolationPFGammaCands = "<<(*theTau).getEMIsolation()<<endl;	
    }
    i_Tau++;    
  }    
}
void PFTauTest::endJob(){}

DEFINE_SEAL_MODULE();
DEFINE_ANOTHER_FWK_MODULE(PFTauTest);
