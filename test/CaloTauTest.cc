#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/TauReco/interface/Tau.h"

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
  int nEvent;
  vector<float> nEventsUsed;
  vector<float> nEventsRiso;
  int nEventTaggedJets;
};

CaloTauTest::CaloTauTest(const ParameterSet& iConfig){
  CaloTaus_        = iConfig.getParameter<string>("Taus");
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

  cout<<"***"<<endl;
  cout<<"Found "<<myTauCollection.size()<<" had. tau-jet candidates"<<endl;
  int it=0;
  for(TauCollection::const_iterator iT =myTauCollection.begin();iT !=myTauCollection.end();iT++)
    {  
      //Prints out some quantities
      cout<<"Jet Number "<<it<<endl;
      it++;
      cout<<"Pt of the Tau "<<iT->pt()<<endl;
      TrackRef theLeadTk = iT->getLeadingTrack();
      if(!theLeadTk) {
	cout<<"No Lead Tk "<<endl;
      }else{
	cout<<"Lead Tk pt "<<(*theLeadTk).pt()<<endl;
	cout<<"InvariantMass of the Tau "<<iT->getInvariantMass()<<endl;
	cout<<"Vertex of the Tau "<<iT->vz()<<endl;
	cout<<"Charge of the Tau "<<iT->charge()<<endl;
	cout<<"Em Over Hadron energy "<<iT->getEmOverChargedEnergy()<<endl;
	cout<<"Max Hadron energy "<<iT->getMaximumHcalTowerEnergy()<<endl;
	cout<<"# Tracks "<<iT->getSelectedTracks().size()<<endl;
	cout<<"Number of Signal Tracks = "<<iT->getSignalTracks().size()<<endl;
	cout<<"Number of Isolation Tracks = "<<iT->getIsolationTracks().size()<<endl;
	cout<<"Sum pT of Isolation Tracks = "<<iT->getSumPtIsolation()<<endl;
	//	cout<<"Sum E_T of Isolation Gamma Candidates = "<<iT->getEMIsolation()<<endl;
	
      }
    }    
}
void CaloTauTest::endJob() { }

DEFINE_SEAL_MODULE();
DEFINE_ANOTHER_FWK_MODULE(CaloTauTest);
