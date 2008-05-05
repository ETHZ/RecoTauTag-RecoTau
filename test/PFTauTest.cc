#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminatorByIsolation.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"

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
  string PFTauDiscriminatorAgainstElectronProducer_;
  int nEvent;
};

PFTauTest::PFTauTest(const ParameterSet& iConfig){
  PFTauProducer_                         = iConfig.getParameter<string>("PFTauProducer");
  PFTauDiscriminatorByIsolationProducer_ = iConfig.getParameter<string>("PFTauDiscriminatorByIsolationProducer");
  PFTauDiscriminatorAgainstElectronProducer_    = iConfig.getParameter<string>("PFTauDiscriminatorAgainstElectronProducer");
  nEvent=0;
}

void PFTauTest::beginJob(){}

void PFTauTest::analyze(const Event& iEvent, const EventSetup& iSetup){
  cout<<"********"<<endl;
  cout<<"Event number "<<nEvent++<<endl;
  
  Handle<PFTauCollection> thePFTauHandle;
  iEvent.getByLabel(PFTauProducer_,thePFTauHandle);
  
  Handle<PFTauDiscriminatorByIsolation> thePFTauDiscriminatorByIsolation;
  iEvent.getByLabel(PFTauDiscriminatorByIsolationProducer_,thePFTauDiscriminatorByIsolation);

  Handle<PFTauDiscriminator> thePFTauDiscriminatorAgainstElectron;
  iEvent.getByLabel(PFTauDiscriminatorAgainstElectronProducer_,thePFTauDiscriminatorAgainstElectron);

  cout<<"***"<<endl;
  cout<<"Found "<<thePFTauHandle->size()<<" hadr. tau-jet candidates ->"<<endl;
  cout<<endl;
  int i_PFTau=0;
  //  for (PFTauCollection::size_type iPFTau=0;iPFTau<thePFTauHandle->size();iPFTau++) {
  const PFTauCollection & myTau  = *(thePFTauHandle.product()); 
  for(PFTauCollection::const_iterator pippo = myTau.begin();pippo!=myTau.end();pippo++){
    //PFTauRef thePFTau(thePFTauHandle,iPFTau);
    //Prints out some quantities
    cout<<"PFTau object number "<<i_PFTau<<endl;
    cout<<"*** check initial PFJet object ***"<<endl;
    cout<<"Its constituents :"<<endl;
    CandidateBaseRefVector theCandidateBaseRefVector=(*pippo).pfTauTagInfoRef()->pfjetRef()->getJetConstituents();
    for(unsigned int i_Constit=0;i_Constit!=theCandidateBaseRefVector.size();i_Constit++) { 
      const PFCandidate* thePFCand=dynamic_cast<const PFCandidate*>(&*(theCandidateBaseRefVector[i_Constit]));
      cout<<*(thePFCand)<<endl;
      
    }
    cout<<"*** check intermediate PFTauTagInfo object ***"<<endl;
    cout <<*(pippo);
    cout<<"***"<<endl;
    i_PFTau++;    
  }    
}
void PFTauTest::endJob(){}

DEFINE_SEAL_MODULE();
DEFINE_ANOTHER_FWK_MODULE(PFTauTest);
