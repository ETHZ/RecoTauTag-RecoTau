import FWCore.ParameterSet.Config as cms

pfRecoTauDiscriminationByLeadingTrackFinding = cms.EDFilter("PFRecoTauDiscriminationByLeadingTrackFinding",
    PFTauProducer = cms.InputTag('pfRecoTauProducer'),
)
