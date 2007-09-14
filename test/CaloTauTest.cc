#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/TauReco/interface/CaloTau.h"
#include "DataFormats/TauReco/interface/CaloTauDiscriminatorByIsolation.h"

#include "Math/GenVector/VectorUtil.h"
#include "Math/GenVector/PxPyPzE4D.h"
#include "DataFormats/TrackReco/interface/Track.h"

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
  string CaloTauProducer_;
  string CaloTauDiscriminatorByIsolationProducer_;
  int nEvent;
  vector<float> nEventsUsed;
  vector<float> nEventsRiso;
  int nEventTaggedJets;
};

CaloTauTest::CaloTauTest(const ParameterSet& iConfig){
  CaloTauProducer_                            = iConfig.getParameter<string>("CaloTauProducer");
  CaloTauDiscriminatorByIsolationProducer_    = iConfig.getParameter<string>("CaloTauDiscriminatorByIsolationProducer");
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
  
  Handle<CaloTauCollection> theCaloTauHandle;
  iEvent.getByLabel(CaloTauProducer_,theCaloTauHandle);
  const CaloTauCollection& theCaloTauCollection=*(theCaloTauHandle.product()); 

  Handle<CaloTauDiscriminatorByIsolation> theCaloTauDiscriminatorByIsolation;
  iEvent.getByLabel(CaloTauDiscriminatorByIsolationProducer_,theCaloTauDiscriminatorByIsolation);

  cout<<"***"<<endl;
  cout<<"Found "<<theCaloTauCollection.size()<<" had. tau-jet candidates"<<endl;
  int i_CaloTau=0;
  for (CaloTauCollection::size_type iCaloTau=0;iCaloTau<theCaloTauHandle->size();iCaloTau++) {
    CaloTauRef theCaloTau(theCaloTauHandle,iCaloTau);
    //Prints out some quantities
    cout<<"Jet Number "<<i_CaloTau<<endl;
    cout<<"CaloDiscriminatorByIsolation value "<<(*theCaloTauDiscriminatorByIsolation)[theCaloTau]<<endl;
    cout<<"Pt of the CaloTau "<<(*theCaloTau).pt()<<endl;
    TrackRef theLeadTk=(*theCaloTau).leadTrack();
    if(!theLeadTk){
      cout<<"No Lead Tk "<<endl;
    }else{
      cout<<"Lead Tk pt "<<(*theLeadTk).pt()<<endl;
      cout<<"InvariantMass of the Tracks system "<<(*theCaloTau).TracksInvariantMass()<<endl;
      cout<<"InvariantMass of the signal Tracks system "<<(*theCaloTau).signalTracksInvariantMass()<<endl;
      cout<<"Inner point position (x,y,z) of the CaloTau ("<<(*theCaloTau).vx()<<","<<(*theCaloTau).vy()<<","<<(*theCaloTau).vz()<<")"<<endl;
      cout<<"Charge of the CaloTau "<<(*theCaloTau).charge()<<endl;
      cout<<"Et of the highest Et HCAL hit "<<(*theCaloTau).maximumHCALhitEt()<<endl;
      cout<<"# Tracks "<<(*theCaloTau).caloTauTagInfoRef()->Tracks().size()<<endl;
      cout<<"# Signal Tracks = "<<(*theCaloTau).signalTracks().size()<<endl;
      cout<<"# Isolation Tracks = "<<(*theCaloTau).isolationTracks().size()<<endl;
      cout<<"Sum of Pt of the Tracks in isolation annulus = "<<(*theCaloTau).isolationTracksPtSum()<<endl;
      cout<<"Sum of Et of the ECAL RecHits in other isolation annulus = "<<(*theCaloTau).isolationECALhitsEtSum()<<endl;
    }
    i_CaloTau++;    
  }    
}
void CaloTauTest::endJob() { }

DEFINE_SEAL_MODULE();
DEFINE_ANOTHER_FWK_MODULE(CaloTauTest);
