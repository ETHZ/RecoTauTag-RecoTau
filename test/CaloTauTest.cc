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

class CaloTauTest : public EDAnalyzer {
public:
  explicit CaloTauTest(const ParameterSet&);
  ~CaloTauTest() {}
  virtual void analyze(const Event& iEvent,const EventSetup& iSetup);
  virtual void beginJob();
  virtual void endJob();
private:
  string CaloTaus_;
  string TauDiscriminatorByIsolationProducer_;
  int nEvent;
  vector<float> nEventsUsed;
  vector<float> nEventsRiso;
  int nEventTaggedJets;
};

CaloTauTest::CaloTauTest(const ParameterSet& iConfig){
  CaloTaus_                            = iConfig.getParameter<string>("CaloTaus");
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

void CaloTauTest::beginJob(){}

void CaloTauTest::analyze(const Event& iEvent, const EventSetup& iSetup){
  cout<<endl;
  cout<<"********"<<endl;
  cout<<"Event number "<<nEvent++<<endl;
  
  Handle<TauCollection> tauHandle;
  iEvent.getByLabel(CaloTaus_,tauHandle);
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
    TrackRef theLeadTk = (*theTau).getleadTrack();
    if(!theLeadTk) {
      cout<<"No Lead Tk "<<endl;
    }else{
      cout<<"Lead Tk pt "<<(*theLeadTk).pt()<<endl;
      cout<<"InvariantMass of the Tau "<<(*theTau).getInvariantMass()<<endl;
      cout<<"Vertex of the Tau "<<(*theTau).vz()<<endl;
      cout<<"Charge of the Tau "<<(*theTau).charge()<<endl;
      cout<<"Em energy fraction "<<(*theTau).getEmEnergyFraction()<<endl;
      cout<<"Max Hadron energy "<<(*theTau).getMaximumHcalEnergy()<<endl;
      cout<<"# Tracks "<<(*theTau).getSelectedTracks().size()<<endl;
      cout<<"Number of Signal Tracks = "<<(*theTau).getSignalTracks().size()<<endl;
      cout<<"Number of Isolation Tracks = "<<(*theTau).getIsolationTracks().size()<<endl;
      cout<<"Sum pT of Isolation Tracks = "<<(*theTau).getSumPtIsolation()<<endl;
    }
    i_Tau++;    
  }    
}
void CaloTauTest::endJob() { }

DEFINE_SEAL_MODULE();
DEFINE_ANOTHER_FWK_MODULE(CaloTauTest);
