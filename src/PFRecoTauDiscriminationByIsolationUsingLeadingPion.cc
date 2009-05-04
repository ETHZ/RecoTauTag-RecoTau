#include "RecoTauTag/RecoTau/interface/PFRecoTauDiscriminationByIsolationUsingLeadingPion.h"

void PFRecoTauDiscriminationByIsolationUsingLeadingPion::produce(Event& iEvent,const EventSetup& iEventSetup){
  Handle<PFTauCollection> thePFTauCollection;
  iEvent.getByLabel(PFTauProducer_,thePFTauCollection);
  
  // fill the AssociationVector object
  auto_ptr<PFTauDiscriminator> thePFTauDiscriminatorByIsolationUsingLeadingPion(new PFTauDiscriminator(PFTauRefProd(thePFTauCollection)));
  
  for(size_t iPFTau=0;iPFTau<thePFTauCollection->size();++iPFTau) {
    PFTauRef thePFTauRef(thePFTauCollection,iPFTau);
    PFTau thePFTau=*thePFTauRef;
    math::XYZVector thePFTau_XYZVector=thePFTau.momentum();   
    PFTauElementsOperators thePFTauElementsOperators(thePFTau);
    
    if (ApplyDiscriminationByTrackerIsolation_){  
      // optional selection by a tracker isolation : ask for 0 charged hadron PFCand / reco::Track in an isolation annulus around a leading PFCand / reco::Track axis
      double theTrackerIsolationDiscriminator = 1.;
      if (ManipulateTracks_insteadofChargedHadrCands_){
         const TrackRefVector& isolationTracks = thePFTau.isolationTracks();
         unsigned int tracksAboveThreshold = 0;
         for(size_t iTrack = 0; iTrack < isolationTracks.size(); ++iTrack)
         {
            if(isolationTracks[iTrack]->pt() > maxChargedPt_) {
               if(++tracksAboveThreshold > TrackerIsolAnnulus_Tracksmaxn_)
               {
                  theTrackerIsolationDiscriminator = 0.;
                  break;
               }
            }
         }
      } else { //use pf candidates instead
         const PFCandidateRefVector& pfIsoChargedCands = thePFTau.isolationPFChargedHadrCands();
         unsigned int tracksAboveThreshold = 0;
         for(size_t iIsoCand = 0; iIsoCand < pfIsoChargedCands.size(); ++iIsoCand)
         {
            if(pfIsoChargedCands[iIsoCand]->pt() > maxChargedPt_) {
               if(++tracksAboveThreshold > TrackerIsolAnnulus_Candsmaxn_) {
                  theTrackerIsolationDiscriminator = 0.;
                  break;
               }
            }
         }
      }

      if (theTrackerIsolationDiscriminator == 0.){
	thePFTauDiscriminatorByIsolationUsingLeadingPion->setValue(iPFTau,0.);
        continue;
      }
    }    
    
    if (ApplyDiscriminationByECALIsolation_){
      // optional selection by an ECAL isolation : ask for 0 gamma PFCand in an isolation annulus around a leading PFCand
      double theECALIsolationDiscriminator =1.;
      const PFCandidateRefVector& pfIsoGammaCands = thePFTau.isolationPFGammaCands();
      unsigned int gammasAboveThreshold = 0;
      for(size_t iIsoGamma = 0; iIsoGamma < pfIsoGammaCands.size(); ++iIsoGamma)
      {
         if(pfIsoGammaCands[iIsoGamma]->pt() > maxGammaPt_) {
            if(++gammasAboveThreshold > ECALIsolAnnulus_Candsmaxn_) {
               theECALIsolationDiscriminator = 0;
               break;
            }
         }
      }
      if (theECALIsolationDiscriminator==0.){
	thePFTauDiscriminatorByIsolationUsingLeadingPion->setValue(iPFTau,0);
	continue;
      }
    }
    
    // not optional selection : ask for a leading (Pt>minPt) PFCand / reco::Track in a matching cone around the PFJet axis
    double theleadElementDiscriminator=1.;
    if (ManipulateTracks_insteadofChargedHadrCands_){
       if (!thePFTau.leadTrack()) theleadElementDiscriminator=0.;
    }else{
       if (!thePFTau.leadPFCand() || (!thePFTau.leadPFChargedHadrCand())) {
          theleadElementDiscriminator=0.;
       }
    }
    if (theleadElementDiscriminator < 0.5) thePFTauDiscriminatorByIsolationUsingLeadingPion->setValue(iPFTau,0);
    else thePFTauDiscriminatorByIsolationUsingLeadingPion->setValue(iPFTau,1); //passes everything
  }    
  
  iEvent.put(thePFTauDiscriminatorByIsolationUsingLeadingPion);
  
}


