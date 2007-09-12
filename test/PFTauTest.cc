#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminatorByIsolation.h"

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
  string PFTauProducer_;
  string PFTauDiscriminatorByIsolationProducer_;
  int nEvent;
  vector<float> nEventsUsed;
  vector<float> nEventsRiso;
  int nEventTaggedJets;
};

PFTauTest::PFTauTest(const ParameterSet& iConfig){
  PFTauProducer_                         = iConfig.getParameter<string>("PFTauProducer");
  PFTauDiscriminatorByIsolationProducer_ = iConfig.getParameter<string>("PFTauDiscriminatorByIsolationProducer");
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
  
  Handle<PFTauCollection> thePFTauHandle;
  iEvent.getByLabel(PFTauProducer_,thePFTauHandle);
  const PFTauCollection& thePFTauCollection=*(thePFTauHandle.product()); 

  Handle<PFTauDiscriminatorByIsolation> thePFTauDiscriminatorByIsolation;
  iEvent.getByLabel(PFTauDiscriminatorByIsolationProducer_,thePFTauDiscriminatorByIsolation);

  cout<<"***"<<endl;
  cout<<"Found "<<thePFTauCollection.size()<<" had. tau-jet candidates"<<endl;
  int i_PFTau=0;
  for (PFTauCollection::size_type iPFTau=0;iPFTau<thePFTauHandle->size();iPFTau++) {
    PFTauRef thePFTau(thePFTauHandle,iPFTau);
    //Prints out some quantities
    cout<<"Jet Number "<<i_PFTau<<endl;
    cout<<"PFDiscriminatorByIsolation value "<<(*thePFTauDiscriminatorByIsolation)[thePFTau]<<endl;
    cout<<"Pt of the PFTau "<<(*thePFTau).pt()<<endl;
    PFCandidateRef theLeadPFCand = (*thePFTau).leadPFChargedHadrCand();
    if(!theLeadPFCand){
      cout<<"No Lead PFCand "<<endl;
    }else{
      cout<<"Lead PFCand Pt "<<(*theLeadPFCand).pt()<<endl;
      cout<<"Vertex of the PFTau "<<(*thePFTau).vz()<<endl;
      cout<<"Charge of the PFTau "<<(*thePFTau).charge()<<endl;
      cout<<"Et of the highest Et HCAL PFCluster "<<(*thePFTau).maximumHCALPFClusterEt()<<endl;
      cout<<"# PF charged hadr. cand's "<<(*thePFTau).pfTauTagInfoRef()->PFChargedHadrCands().size()<<endl;
      cout<<"# PF neutral hadr. cand's "<<(*thePFTau).pfTauTagInfoRef()->PFNeutrHadrCands().size()<<endl;
      cout<<"# PF gamma cand's "<<(*thePFTau).pfTauTagInfoRef()->PFGammaCands().size()<<endl;
      cout<<"Number of SignalPFChargedHadrCands = "<<(*thePFTau).signalPFChargedHadrCands().size()<<endl;
      cout<<"Number of IsolationPFChargedHadrCands = "<<(*thePFTau).isolationPFChargedHadrCands().size()<<endl;
      cout<<"Number of SignalPFGammaCands = "<<(*thePFTau).signalPFGammaCands().size()<<endl;
      cout<<"Number of IsolationPFGammaCands = "<<(*thePFTau).isolationPFGammaCands().size()<<endl;
      cout<<"Sum of Pt of charged hadr. PFCandidates in isolation annulus  = "<<(*thePFTau).isolationPFChargedHadrCandsPtSum()<<endl;
      cout<<"Sum of Et of gamma PFCandidates in other isolation annulus = "<<(*thePFTau).isolationPFGammaCandsEtSum()<<endl;	
    }
    i_PFTau++;    
  }    
}
void PFTauTest::endJob(){}

DEFINE_SEAL_MODULE();
DEFINE_ANOTHER_FWK_MODULE(PFTauTest);
