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

class PFTauTest : public EDAnalyzer {
public:
  explicit PFTauTest(const ParameterSet&);
  ~PFTauTest() {}
  virtual void analyze(const Event& iEvent,const EventSetup& iSetup);
  virtual void beginJob();
  virtual void endJob();
private:
  string PFTaus_;
  int nEvent;
  vector<float> nEventsUsed;
  vector<float> nEventsRiso;
  int nEventTaggedJets;
};

PFTauTest::PFTauTest(const ParameterSet& iConfig){
  PFTaus_        = iConfig.getParameter<string>("PFTaus");
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
  cout<<endl;
  cout<<"********"<<endl;
  cout<<"Event number "<<nEvent++<<endl;
  
  Handle<TauCollection> tauHandle;
  iEvent.getByLabel(PFTaus_,tauHandle);
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
      cout<<"InvariantMass of the Tau "<<iT->getInvariantMass()<<endl;
      cout<<"Vertex of the Tau "<<iT->vz()<<endl;
      cout<<"Charge of the Tau "<<iT->charge()<<endl;
      cout<<"Em Over Hadron energy "<<iT->getEmOverHadronEnergy()<<endl;
      cout<<"Max Hadron energy "<<iT->getMaximumHcalTowerEnergy()<<endl;
      cout<<"# PF charged hadr. cand's "<<iT->getSelectedChargedHadrons().size()<<endl;
      cout<<"# PF neutral hadr. cand's "<<iT->getSelectedNeutralHadrons().size()<<endl;
      cout<<"# PF gamma cand's "<<iT->getSelectedGammaCandidates().size()<<endl;
      PFCandidateRef theLeadPFCand = iT->getLeadingChargedHadron();
      if(!theLeadPFCand) cout<<"No Lead PFCand "<<endl;
      else{
	cout<<"Lead PFCand pt "<<(*theLeadPFCand).pt()<<endl;
	cout<<"Number of SignalPFChargedHadrons = "<<iT->getSignalChargedHadrons().size()<<endl;
	cout<<"Number of IsolationPFChargedHadrons = "<<iT->getIsolationChargedHadrons().size()<<endl;
	cout<<"Number of SignalPFGammaCandidate = "<<iT->getSignalGammaCandidates().size()<<endl;
	cout<<"Number of IsolationPFGammaCandidate = "<<iT->getIsolationGammaCandidates().size()<<endl;
	cout<<"Sum pT of Isolation Charged Hadrons = "<<iT->getSumPtIsolation()<<endl;
	cout<<"Sum E_T of Isolation Gamma Candidates = "<<iT->getEMIsolation()<<endl;
	
      }
    }    
}
void PFTauTest::endJob() { }

DEFINE_SEAL_MODULE();
DEFINE_ANOTHER_FWK_MODULE(PFTauTest);
