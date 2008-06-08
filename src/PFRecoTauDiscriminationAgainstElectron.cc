#include "RecoTauTag/RecoTau/interface/PFRecoTauDiscriminationAgainstElectron.h"

void PFRecoTauDiscriminationAgainstElectron::produce(Event& iEvent,const EventSetup& iEventSetup){
  Handle<PFTauCollection> thePFTauCollection;
  iEvent.getByLabel(PFTauProducer_,thePFTauCollection);

  ESHandle<MagneticField> myMF;
  iEventSetup.get<IdealMagneticFieldRecord>().get(myMF);
  const MagneticField* MagneticField = myMF.product(); 

  // fill the AssociationVector object
  auto_ptr<PFTauDiscriminator> thePFTauDiscriminatorAgainstElectron(new PFTauDiscriminator(PFTauRefProd(thePFTauCollection)));

  for(size_t iPFTau=0;iPFTau<thePFTauCollection->size();++iPFTau) {
    PFTauRef thePFTauRef(thePFTauCollection,iPFTau);

    // Check if track goes to Ecal crack
    TrackRef myleadTk=(*thePFTauRef).leadPFChargedHadrCand()->trackRef();
    math::XYZPoint myleadTkEcalPos;

    if(myleadTk.isNonnull()){ 
      if(MagneticField!=0){ 
	myleadTkEcalPos = TauTagTools::propagTrackECALSurfContactPoint(MagneticField,myleadTk);
      } else {
	// temporary: outer position is not correct!
	myleadTkEcalPos = myleadTk->outerPosition();
      }
      if (applyCut_ecalCrack_ && isInEcalCrack(abs((double)myleadTkEcalPos.eta()))) {
	thePFTauDiscriminatorAgainstElectron->setValue(iPFTau,0);
	continue;
      }
    }

    bool decision = false;
    bool emfPass = true, htotPass = true, hmaxPass = true; 
    bool h3x3Pass = true, estripPass = true, erecovPass = true, epreidPass = true;

    if (applyCut_emFraction_) {
      if ((*thePFTauRef).emFraction() > emFraction_maxValue_) {
	emfPass = false;
      }
    }
    if (applyCut_hcalTotOverPLead_) {
      if ((*thePFTauRef).hcalTotOverPLead() < hcalTotOverPLead_minValue_) {
	htotPass = false;
      }
    }
    if (applyCut_hcalMaxOverPLead_) {
      if ((*thePFTauRef).hcalMaxOverPLead() < hcalMaxOverPLead_minValue_) {
	hmaxPass = false;
      }
    }
    if (applyCut_hcal3x3OverPLead_) {
      if ((*thePFTauRef).hcal3x3OverPLead() < hcal3x3OverPLead_minValue_) {
	h3x3Pass = false;
      }
    }
    if (applyCut_ecalStripSumEOverPLead_) {
      if ((*thePFTauRef).ecalStripSumEOverPLead() > ecalStripSumEOverPLead_minValue_ &&
	  (*thePFTauRef).ecalStripSumEOverPLead() < ecalStripSumEOverPLead_maxValue_) {
	estripPass = false;
      } else {
	estripPass = true;
      }
    }
    if (applyCut_bremsRecoveryEOverPLead_) {
      if ((*thePFTauRef).bremsRecoveryEOverPLead() > bremsRecoveryEOverPLead_minValue_ &&
	  (*thePFTauRef).bremsRecoveryEOverPLead() < bremsRecoveryEOverPLead_maxValue_) {
	erecovPass = false;
      } else {
	erecovPass = true;
      } 
    }
    if (applyCut_electronPreID_) {
      if (
	  ((*thePFTauRef).electronPreIDDecision() &&
	   ((*thePFTauRef).ecalStripSumEOverPLead() < elecPreID1_SumEOverPLead_maxValue ||
	    (*thePFTauRef).hcal3x3OverPLead() > elecPreID1_Hcal3x3_minValue))
	  ||
	  (!(*thePFTauRef).electronPreIDDecision() &&
	   ((*thePFTauRef).ecalStripSumEOverPLead() < elecPreID0_SumEOverPLead_maxValue ||
	    (*thePFTauRef).hcal3x3OverPLead() > elecPreID0_Hcal3x3_minValue))
	  ){
	epreidPass = true;
      }  else {
	epreidPass = false;
      }
    }

    decision = emfPass && htotPass && hmaxPass && 
      h3x3Pass && estripPass && erecovPass && epreidPass;
    if (decision) {
      thePFTauDiscriminatorAgainstElectron->setValue(iPFTau,1);
    } else {
      thePFTauDiscriminatorAgainstElectron->setValue(iPFTau,0);
    }
  }


  iEvent.put(thePFTauDiscriminatorAgainstElectron);
}

bool
PFRecoTauDiscriminationAgainstElectron::isInEcalCrack(double eta) const{  
  return (eta < 0.018 || 
	  (eta>0.423 && eta<0.461) ||
	  (eta>0.770 && eta<0.806) ||
	  (eta>1.127 && eta<1.163) ||
	  (eta>1.460 && eta<1.558));
}
