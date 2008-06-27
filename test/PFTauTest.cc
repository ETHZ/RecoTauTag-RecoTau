#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminatorByIsolation.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
#include "DataFormats/TauReco/interface/PFTauTagInfo.h"

#include "RecoTauTag/TauTagTools/interface/PFTauElementsOperators.h"
#include "RecoTauTag/TauTagTools/interface/TauTagTools.h"

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
#include <TH2.h>
#include <TDirectory.h>

#include "FWCore/ServiceRegistry/interface/Service.h" // Framework services
#include "PhysicsTools/UtilAlgos/interface/TFileService.h" // Framework service for histograms

#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include <Math/GenVector/VectorUtil.h>
#include "TLorentzVector.h"

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
  string PFTauDiscriminatorAgainstMuonProducer_;
  string GenParticles_;

  
  int nEvent;
  int nTauMatchPFTau;
  int nElecMatchPFTau;
  int nTauElecPreID;
  int nElecElecPreID;
  int nTauNonElecPreID;
  int nElecNonElecPreID;

  // files
  std::string _fileName;
  TFile *_file;
  TDirectory *_dir;

  // Gen Objects
  std::vector<TLorentzVector> _GenElecs;
  std::vector<TLorentzVector> _GenMuons;
  std::vector<TLorentzVector> _GenTauElecs;
  std::vector<TLorentzVector> _GenTauMuons;
  std::vector<TLorentzVector> _GenTauHads;

  // Efficiency plots
  TH1F* h_TauEff_ElecPreID;
  TH1F* h_ElecEff_ElecPreID;
  TH1F* h_TauEff_PFTauD;
  TH1F* h_ElecEff_PFTauD;

  TH1F* h_TauEff_EoP;
  TH1F* h_TauEff_HoP;
  TH2F* h_TauEff_HvsEoP;

  TH1F* h_TauEff_eid0_EoP;
  TH1F* h_TauEff_eid0_HoP;
  TH2F* h_TauEff_eid0_HvsEoP;

  TH1F* h_TauEff_eid1_EoP;
  TH1F* h_TauEff_eid1_HoP;
  TH2F* h_TauEff_eid1_HvsEoP;

  TH1F* h_ElecEff_EoP;
  TH1F* h_ElecEff_HoP;
  TH2F* h_ElecEff_HvsEoP;

  TH1F* h_ElecEff_eid0_EoP;
  TH1F* h_ElecEff_eid0_HoP;
  TH2F* h_ElecEff_eid0_HvsEoP;

  TH1F* h_ElecEff_eid1_EoP;
  TH1F* h_ElecEff_eid1_HoP;
  TH2F* h_ElecEff_eid1_HvsEoP;
};

PFTauTest::PFTauTest(const ParameterSet& iConfig){
  PFTauProducer_                         = iConfig.getParameter<string>("PFTauProducer");
  PFTauDiscriminatorByIsolationProducer_ = iConfig.getParameter<string>("PFTauDiscriminatorByIsolationProducer");
  PFTauDiscriminatorAgainstElectronProducer_    = iConfig.getParameter<string>("PFTauDiscriminatorAgainstElectronProducer");
  PFTauDiscriminatorAgainstMuonProducer_    = iConfig.getParameter<string>("PFTauDiscriminatorAgainstMuonProducer");
  GenParticles_    = iConfig.getParameter<string>("GenParticles");

  nEvent=0;

  nTauMatchPFTau=0;
  nElecMatchPFTau=0;
  nTauElecPreID=0;
  nElecElecPreID=0;
  nTauNonElecPreID=0;
  nElecNonElecPreID=0;

  edm::Service<TFileService> fs;
  TFileDirectory dir = fs->mkdir("histos");

  // Book histograms
  h_TauEff_ElecPreID = dir.make<TH1F>("TauEff_ElecPreID","TauEff_ElecPreID",3,0.,1.01);
  h_TauEff_ElecPreID->Sumw2();
  h_ElecEff_ElecPreID = dir.make<TH1F>("ElecEff_ElecPreID","ElecEff_ElecPreID",3,0.,1.01);
  h_ElecEff_ElecPreID->Sumw2();
  h_TauEff_PFTauD = dir.make<TH1F>("TauEff_PFTauD","TauEff_PFTauD",3,0.,1.01);
  h_TauEff_PFTauD->Sumw2();
  h_ElecEff_PFTauD = dir.make<TH1F>("ElecEff_PFTauD","ElecEff_PFTauD",3,0.,1.01);
  h_ElecEff_PFTauD->Sumw2();

  h_TauEff_EoP = dir.make<TH1F>("TauEff_EoP","TauEff_EoP",100,0.,5.);
  h_TauEff_EoP->Sumw2();
  h_TauEff_HoP = dir.make<TH1F>("TauEff_HoP","TauEff_HoP",100,0.,5.);
  h_TauEff_HoP->Sumw2();
  h_TauEff_HvsEoP = dir.make<TH2F>("TauEff_HvsEoP","TauEff_HvsEoP",60,0.,4.,60,0.,4.);
  h_TauEff_HvsEoP->Sumw2();

  h_TauEff_eid0_EoP = dir.make<TH1F>("TauEff_eid0_EoP","TauEff_eid0_EoP",100,0.,5.);
  h_TauEff_eid0_EoP->Sumw2();
  h_TauEff_eid0_HoP = dir.make<TH1F>("TauEff_eid0_HoP","TauEff_eid0_HoP",100,0.,5.);
  h_TauEff_eid0_HoP->Sumw2();
  h_TauEff_eid0_HvsEoP = dir.make<TH2F>("TauEff_eid0_HvsEoP","TauEff_eid0_HvsEoP",60,0.,4.,60,0.,4.);
  h_TauEff_eid0_HvsEoP->Sumw2();

  h_TauEff_eid1_EoP = dir.make<TH1F>("TauEff_eid1_EoP","TauEff_eid1_EoP",100,0.,5.);
  h_TauEff_eid1_EoP->Sumw2();
  h_TauEff_eid1_HoP = dir.make<TH1F>("TauEff_eid1_HoP","TauEff_eid1_HoP",100,0.,5.);
  h_TauEff_eid1_HoP->Sumw2();
  h_TauEff_eid1_HvsEoP = dir.make<TH2F>("TauEff_eid1_HvsEoP","TauEff_eid1_HvsEoP",60,0.,4.,60,0.,4.);
  h_TauEff_eid1_HvsEoP->Sumw2();

  
  h_ElecEff_EoP = dir.make<TH1F>("ElecEff_EoP","ElecEff_EoP",100,0.,5.);
  h_ElecEff_EoP->Sumw2();
  h_ElecEff_HoP = dir.make<TH1F>("ElecEff_HoP","ElecEff_HoP",100,0.,5.);
  h_ElecEff_HoP->Sumw2();
  h_ElecEff_HvsEoP = dir.make<TH2F>("ElecEff_HvsEoP","ElecEff_HvsEoP",60,0.,4.,60,0.,4.);
  h_ElecEff_HvsEoP->Sumw2();

  h_ElecEff_eid0_EoP = dir.make<TH1F>("ElecEff_eid0_EoP","ElecEff_eid0_EoP",100,0.,5.);
  h_ElecEff_eid0_EoP->Sumw2();
  h_ElecEff_eid0_HoP = dir.make<TH1F>("ElecEff_eid0_HoP","ElecEff_eid0_HoP",100,0.,5.);
  h_ElecEff_eid0_HoP->Sumw2();
  h_ElecEff_eid0_HvsEoP = dir.make<TH2F>("ElecEff_eid0_HvsEoP","ElecEff_eid0_HvsEoP",60,0.,4.,60,0.,4.);
  h_ElecEff_eid0_HvsEoP->Sumw2();

  h_ElecEff_eid1_EoP = dir.make<TH1F>("ElecEff_eid1_EoP","ElecEff_eid1_EoP",100,0.,5.);
  h_ElecEff_eid1_EoP->Sumw2();
  h_ElecEff_eid1_HoP = dir.make<TH1F>("ElecEff_eid1_HoP","ElecEff_eid1_HoP",100,0.,5.);
  h_ElecEff_eid1_HoP->Sumw2();
  h_ElecEff_eid1_HvsEoP = dir.make<TH2F>("ElecEff_eid1_HvsEoP","ElecEff_eid1_HvsEoP",60,0.,4.,60,0.,4.);
  h_ElecEff_eid1_HvsEoP->Sumw2();

}

void PFTauTest::beginJob(){
}

void PFTauTest::analyze(const Event& iEvent, const EventSetup& iSetup){
  //cout<<"********"<<endl;
  //cout<<"Event number "<<nEvent++<<endl;



  ////////////////////////////////////////////////////////  

  edm::Handle<HepMCProduct> evt;
  iEvent.getByLabel(GenParticles_, evt);
  HepMC::GenEvent * generated_event = new HepMC::GenEvent(*(evt->GetEvent()));

  //init
  _GenTauHads.clear();;
  _GenTauMuons.clear();;
  _GenTauElecs.clear();;
  _GenMuons.clear();;
  _GenElecs.clear();;
  
  int ntaujet=0,nelec=0,nmuon=0;  
  TLorentzVector taunet,tauelec,taumuon;
  HepMC::GenEvent::particle_iterator p;
  for (p = generated_event->particles_begin(); p != generated_event->particles_end(); p++) {
    if(abs((*p)->pdg_id()) == 15&&(*p)->status()==2) { 
      bool lept_decay = false;
      
      TLorentzVector tau((*p)->momentum().px(),(*p)->momentum().py(),(*p)->momentum().pz(),(*p)->momentum().e());
      HepMC::GenVertex::particle_iterator z = (*p)->end_vertex()->particles_begin(HepMC::descendants);
      for(; z != (*p)->end_vertex()->particles_end(HepMC::descendants); z++) {
	if(abs((*z)->pdg_id()) == 11 || abs((*z)->pdg_id()) == 13) lept_decay=true;
	if(abs((*z)->pdg_id()) == 11) {
	  nelec++;
	  tauelec.SetPxPyPzE((*z)->momentum().px(),(*z)->momentum().py(),(*z)->momentum().pz(),(*z)->momentum().e());
	}
	if(abs((*z)->pdg_id()) == 13) {
	  nmuon++;
	  taumuon.SetPxPyPzE((*z)->momentum().px(),(*z)->momentum().py(),(*z)->momentum().pz(),(*z)->momentum().e());
	  _GenTauMuons.push_back(taumuon);
	}
	if(abs((*z)->pdg_id()) == 16)
	  taunet.SetPxPyPzE((*z)->momentum().px(),(*z)->momentum().py(),(*z)->momentum().pz(),(*z)->momentum().e());
	
      }
      if(lept_decay==false) {
	ntaujet++;
	TLorentzVector jetMom=tau-taunet;
	_GenTauHads.push_back(jetMom);
      }
    } else if(abs((*p)->pdg_id()) == 11&&(*p)->status()==1) { 
      TLorentzVector elec((*p)->momentum().px(),(*p)->momentum().py(),(*p)->momentum().pz(),(*p)->momentum().e());
      _GenElecs.push_back(elec);
    } else if(abs((*p)->pdg_id()) == 13&&(*p)->status()==1) { 
      TLorentzVector mu((*p)->momentum().px(),(*p)->momentum().py(),(*p)->momentum().pz(),(*p)->momentum().e());
      _GenMuons.push_back(mu);
    }

  }

  ////////////////////////////////////////////////////////  

  
  
  Handle<PFTauCollection> thePFTauHandle;
  iEvent.getByLabel(PFTauProducer_,thePFTauHandle);
  
  Handle<PFTauDiscriminatorByIsolation> thePFTauDiscriminatorByIsolation;
  iEvent.getByLabel(PFTauDiscriminatorByIsolationProducer_,thePFTauDiscriminatorByIsolation);

  Handle<PFTauDiscriminator> thePFTauDiscriminatorAgainstElectron;
  iEvent.getByLabel(PFTauDiscriminatorAgainstElectronProducer_,thePFTauDiscriminatorAgainstElectron);
  
  //Handle<PFTauDiscriminator> thePFTauDiscriminatorAgainstMuon;
  //iEvent.getByLabel(PFTauDiscriminatorAgainstMuonProducer_,thePFTauDiscriminatorAgainstMuon);


  Handle<PFJetCollection> thePFJetHandle;
  iEvent.getByLabel("iterativeCone5PFJets",thePFJetHandle);

  /*
  int n=0;
  for (PFJetCollection::size_type iPFJet=0;iPFJet<thePFJetHandle->size();iPFJet++) { 
    PFJetRef thePFJet(thePFJetHandle,iPFJet); 

    std::vector <const reco::PFCandidate*> PFJetCands = (*thePFJet).getPFConstituents();

    for (unsigned iPFJetCand=0;iPFJetCand<PFJetCands.size();iPFJetCand++) { 
      if(PFJetCands[iPFJetCand]!=NULL){
	cout<<"****************************"<<endl;
	cout<<"mva: "<<PFJetCands[iPFJetCand]->mva_e_pi()<<endl;
	cout<<PFJetCands[iPFJetCand]->et()<<" "<<PFJetCands[iPFJetCand]->eta()<<" "<<PFJetCands[iPFJetCand]->phi()<<endl;
      } else {
	cout<<++n<<" ============================="<<endl;
      }
    }
  }
  */
  
  //int n = 0;
  // Tau Loop
  for (PFTauCollection::size_type iPFTau=0;iPFTau<thePFTauHandle->size();iPFTau++) { 
    PFTauRef thePFTau(thePFTauHandle,iPFTau); 

    /*
    PFTau myPFTau(*thePFTau); 
    PFTauElementsOperators myPFTauElementsOperators(myPFTau);
    TFormula myMatchingConeSizeTFormula=myPFTauElementsOperators.computeConeSizeTFormula("0.1","Matching cone size");
    double myMatchingConeSize=myPFTauElementsOperators.computeConeSize(myMatchingConeSizeTFormula,0.,0.6);
    PFCandidateRef myleadPFCand=myPFTauElementsOperators.leadPFChargedHadrCand("DR",myMatchingConeSize,5.);

    if(myleadPFCand.isNonnull()){
      cout<<"****************************"<<endl;
      cout<<myleadPFCand->mva_e_pi()<<endl;
      cout<<myleadPFCand->et()<<" "<<myleadPFCand->eta()<<" "<<myleadPFCand->phi()<<endl;
      //}
      cout<<(*thePFTau).electronPreIDDecision()<<endl;
      cout<<(*thePFTau).emFraction()<<" "<<(*thePFTau).hcal3x3OverPLead()<<endl;
      cout<<(*thePFTau).et()<<" "<<(*thePFTau).eta()<<" "<<(*thePFTau).phi()<<endl;
      cout<<"DiscriminatonByIsolation: "<<(*thePFTauDiscriminatorByIsolation)[thePFTau]<<endl;
    } else {
      cout<<++n<<" ============================="<<endl;
    }
    */
    
    /*
    if ((*thePFTau).et() > 15. && abs((*thePFTau).eta()) < 2.5) {
      if ((*thePFTau).electronPreIDDecision()) {
	cout<<"****************************"<<endl;
	cout<<(*thePFTau).et()<<" "<<(*thePFTau).eta()<<" "<<(*thePFTau).phi()<<endl;
	cout<<(*thePFTau).electronPreIDDecision()<<" "<<(*thePFTau).emFraction()<<" "<<(*thePFTau).hcal3x3OverPLead()<<endl;
	cout<<"DiscriminatonByIsolation: "<<(*thePFTauDiscriminatorByIsolation)[thePFTau]<<endl;
      }
      if ((*thePFTau).electronPreIDDecision()) {
	nTauNonElecPreID++;		
      } else {
	nTauElecPreID++;		
      }      
    }
    */
    if ((*thePFTauDiscriminatorByIsolation)[thePFTau] == 1) {
      if ((*thePFTau).et() > 15. && abs((*thePFTau).eta()) < 2.5) {
	for (unsigned int i = 0; i<_GenTauHads.size();i++) {
	  if (_GenTauHads[i].Et() >= 10. && abs(_GenTauHads[i].Eta()) <= 2.5 ) {
	    TLorentzVector pftau((*thePFTau).px(),(*thePFTau).py(),(*thePFTau).pz(),(*thePFTau).energy());
	    double deltaR = ROOT::Math::VectorUtil::DeltaR(_GenTauHads[i],pftau);
	    if (deltaR<0.25) {
	      nTauMatchPFTau++;
	      
	      /*
		if ((*thePFTau).electronPreIDDecision()) {
		cout<<"****************************"<<endl;
		cout<<(*thePFTau).et()<<" "<<(*thePFTau).eta()<<" "<<(*thePFTau).phi()<<endl;
		cout<<(*thePFTau).electronPreIDDecision()<<" "<<(*thePFTau).emFraction()<<" "<<(*thePFTau).hcal3x3OverPLead()<<endl;
		}
	      */
	      
	      h_TauEff_EoP->Fill((*thePFTau).ecalStripSumEOverPLead());
	      h_TauEff_HoP->Fill((*thePFTau).hcalTotOverPLead());
	      h_TauEff_HvsEoP->Fill((*thePFTau).ecalStripSumEOverPLead(),(*thePFTau).hcalTotOverPLead());
	      
	      if ((*thePFTau).electronPreIDDecision()) {
		h_TauEff_ElecPreID->Fill(1.);
		h_TauEff_eid1_EoP->Fill((*thePFTau).ecalStripSumEOverPLead());
		h_TauEff_eid1_HoP->Fill((*thePFTau).hcalTotOverPLead());
		h_TauEff_eid1_HvsEoP->Fill((*thePFTau).ecalStripSumEOverPLead(),(*thePFTau).hcalTotOverPLead());
		////nTauElecPreID++;		
	      } else {
		h_TauEff_ElecPreID->Fill(0.);
		h_TauEff_ElecPreID->Fill(1.);
		h_TauEff_eid0_EoP->Fill((*thePFTau).ecalStripSumEOverPLead());
		h_TauEff_eid0_HoP->Fill((*thePFTau).hcalTotOverPLead());
		h_TauEff_eid0_HvsEoP->Fill((*thePFTau).ecalStripSumEOverPLead(),(*thePFTau).hcalTotOverPLead());
		////nTauNonElecPreID++;		
	      }
	      //cout<<"DiscriminatorAgainstElectron: "<<(*thePFTauDiscriminatorAgainstElectron)[thePFTau]<<endl;
	      if ((*thePFTauDiscriminatorAgainstElectron)[thePFTau] == 1) {
		h_TauEff_PFTauD->Fill(1.);
		nTauNonElecPreID++;
	      } else {
		h_TauEff_PFTauD->Fill(0.);
		nTauElecPreID++;
	      }
	    }
	  }
	}
      }
    }
  
  }
  
  // Electron Loop
  for (PFTauCollection::size_type iPFTau=0;iPFTau<thePFTauHandle->size();iPFTau++) { 
    PFTauRef thePFTau(thePFTauHandle,iPFTau); 
    /*
    if ((*thePFTau).et() > 15. && abs((*thePFTau).eta()) < 2.5) {
      if ((*thePFTau).electronPreIDDecision()) {
	cout<<"****************************"<<endl;
	cout<<(*thePFTau).et()<<" "<<(*thePFTau).eta()<<" "<<(*thePFTau).phi()<<endl;
	cout<<(*thePFTau).electronPreIDDecision()<<" "<<(*thePFTau).emFraction()<<" "<<(*thePFTau).hcal3x3OverPLead()<<endl;
	cout<<"DiscriminatonByIsolation: "<<(*thePFTauDiscriminatorByIsolation)[thePFTau]<<endl;
      }
      if ((*thePFTau).electronPreIDDecision()) {
	nElecNonElecPreID++;		
      } else {
	nElecElecPreID++;		
      }      
    }
    */
    if ((*thePFTauDiscriminatorByIsolation)[thePFTau] == 1) {
      if ((*thePFTau).et() > 15. && abs((*thePFTau).eta()) < 2.5) {
	for (unsigned int i = 0; i<_GenElecs.size();i++) {
	  if (_GenElecs[i].Et() >= 10. && abs(_GenElecs[i].Eta()) <= 2.5 ) {
	    TLorentzVector pftau((*thePFTau).px(),(*thePFTau).py(),(*thePFTau).pz(),(*thePFTau).energy());
	    double deltaR = ROOT::Math::VectorUtil::DeltaR(_GenElecs[i],pftau);
	    if (deltaR<0.25) {
	      nElecMatchPFTau++;
	      
	      /*
		if ((*thePFTau).electronPreIDDecision()) {
		cout<<"****************************"<<endl;
		cout<<(*thePFTau).et()<<" "<<(*thePFTau).eta()<<" "<<(*thePFTau).phi()<<endl;
		cout<<(*thePFTau).electronPreIDDecision()<<" "<<(*thePFTau).emFraction()<<" "<<(*thePFTau).hcal3x3OverPLead()<<endl;
		}
	      */

	      h_ElecEff_EoP->Fill((*thePFTau).ecalStripSumEOverPLead());
	      h_ElecEff_HoP->Fill((*thePFTau).hcalTotOverPLead());
	      h_ElecEff_HvsEoP->Fill((*thePFTau).ecalStripSumEOverPLead(),(*thePFTau).hcalTotOverPLead());
	      
	      if ((*thePFTau).electronPreIDDecision()) {
		h_ElecEff_ElecPreID->Fill(1.);
		h_ElecEff_eid1_EoP->Fill((*thePFTau).ecalStripSumEOverPLead());
		h_ElecEff_eid1_HoP->Fill((*thePFTau).hcalTotOverPLead());
		h_ElecEff_eid1_HvsEoP->Fill((*thePFTau).ecalStripSumEOverPLead(),(*thePFTau).hcalTotOverPLead());
		////nElecElecPreID++;		
	      } else {
		h_ElecEff_ElecPreID->Fill(0.);
		h_ElecEff_eid0_EoP->Fill((*thePFTau).ecalStripSumEOverPLead());
		h_ElecEff_eid0_HoP->Fill((*thePFTau).hcalTotOverPLead());
		h_ElecEff_eid0_HvsEoP->Fill((*thePFTau).ecalStripSumEOverPLead(),(*thePFTau).hcalTotOverPLead());
		////nElecNonElecPreID++;		
	      }
	      //cout<<"DiscriminatorAgainstElectron: "<<(*thePFTauDiscriminatorAgainstElectron)[thePFTau]<<endl;
	      if ((*thePFTauDiscriminatorAgainstElectron)[thePFTau] == 1) {
		h_ElecEff_PFTauD->Fill(1.);
		nElecNonElecPreID++;
		//cout<<"****************************"<<endl;
		//cout<<(*thePFTau).et()<<" "<<(*thePFTau).eta()<<" "<<(*thePFTau).phi()<<endl;
		//cout<<(*thePFTau).electronPreIDDecision()<<" "<<(*thePFTau).emFraction()<<" "<<(*thePFTau).hcal3x3OverPLead()<<endl;
		
	      } else {
		h_ElecEff_PFTauD->Fill(0.);
		nElecElecPreID++;		
	      }
	    }
	  }
	}
      }
    }
      
  }
  

  /*
  cout<<"***"<<endl;
  cout<<"Found "<<thePFTauHandle->size()<<" hadr. tau-jet candidates ->"<<endl;
  cout<<endl;
  int i_PFTau=0;
  //  for (PFTauCollection::size_type iPFTau=0;iPFTau<thePFTauHandle->size();iPFTau++) {
  const PFTauCollection & myTau  = *(thePFTauHandle.product()); 
  for(PFTauCollection::const_iterator pippo = myTau.begin();pippo!=myTau.end();pippo++){
    PFTauRef thePFTau(thePFTauHandle,i_PFTau);
    //Prints out some quantities
    cout<<"PFTau object number "<<i_PFTau<<endl;
    cout<<"*** check initial PFJet object ***"<<endl;
    cout<<"Its constituents :"<<endl;
    vector<CandidatePtr> theCandidateBasePtrVector=(*pippo).pfTauTagInfoRef()->pfjetRef()->getJetConstituents();
    for(unsigned int i_Constit=0;i_Constit!=theCandidateBasePtrVector.size();i_Constit++) { 
      const PFCandidate* thePFCand=dynamic_cast<const PFCandidate*>(&*(theCandidateBasePtrVector[i_Constit]));
      cout<<*(thePFCand)<<endl;
      
    }
    cout<<"*** check intermediate PFTauTagInfo object ***"<<endl;
    cout <<*(pippo);
    cout<<"***"<<endl;
    cout<<"**** Discriminators "<<endl;
    float discriminator = (*thePFTauDiscriminatorByIsolation)[thePFTau];
    cout<<"     Isolation "<< discriminator <<endl;
    discriminator = (*thePFTauDiscriminatorAgainstElectron)[thePFTau];
      cout<<"     Against Electrons "<< discriminator<<endl;
    discriminator = (*thePFTauDiscriminatorAgainstMuon)[thePFTau];
    cout<<"     Against Muons "<<discriminator <<endl;

    i_PFTau++;  

  }    
  */
  
}

void PFTauTest::endJob(){
  cout<<"**********************************"<<endl;
  cout<<"Electron rejection efficiencies"<<endl;
  cout<<"**********************************"<<endl;
  if (nTauMatchPFTau>0) {
    cout<<"Taus ("<<nTauMatchPFTau<<"):"<<endl;
    cout<<"---------------------"<<endl;
    cout<<"  Passed: "<< (float)nTauNonElecPreID/(float)nTauMatchPFTau <<endl;
    cout<<"Rejected: "<< (float)nTauElecPreID/(float)nTauMatchPFTau <<endl;
  }
  if (nElecMatchPFTau>0) {
    cout<<"\nElectrons ("<<nElecMatchPFTau<<"):"<<endl;
    cout<<"---------------------"<<endl;
    cout<<"  Passed: "<< (float)nElecNonElecPreID/(float)nElecMatchPFTau <<endl;
    cout<<"Rejected: "<< (float)nElecElecPreID/(float)nElecMatchPFTau <<endl;
  }
}

DEFINE_SEAL_MODULE();
DEFINE_ANOTHER_FWK_MODULE(PFTauTest);
