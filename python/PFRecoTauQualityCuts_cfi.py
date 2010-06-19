import FWCore.ParameterSet.Config as cms

# A set of quality cuts used for the PFTaus.  Note that the quality cuts are
# different for the signal and isolation regions.  (Currently, only in Nhits)

PFTauQualityCuts = cms.PSet(
    signalQualityCuts = cms.PSet(
        minTrackPt                   = cms.double(0.5),  # filter PFChargedHadrons below given pt
        maxTrackChi2                 = cms.double(100.), # require track Chi2
        maxTransverseImpactParameter = cms.double(0.1),  # w.r.t. PV
        maxDeltaZ                    = cms.double(1.0),  # w.r.t. PV
        minTrackPixelHits            = cms.uint32(0),    # pixel-only hits (note that these cuts are turned off, 
                                                          # the tracking cuts might be higher)
        minTrackHits                 = cms.uint32(3),    # total track hits
        minGammaEt                   = cms.double(0.5),  # filter PFgammas below given Pt
        useTracksInsteadOfPFHadrons  = cms.bool(False),  # if true, use generalTracks, instead of PFChargedHadrons
    ),
    isolationQualityCuts = cms.PSet(
        minTrackPt                   = cms.double(1.0), 
        maxTrackChi2                 = cms.double(100.),
        maxTransverseImpactParameter = cms.double(0.1),
        maxDeltaZ                    = cms.double(1.0),
        minTrackPixelHits            = cms.uint32(0),
        minTrackHits                 = cms.uint32(5),
        minGammaEt                   = cms.double(1.5),
        useTracksInsteadOfPFHadrons  = cms.bool(False),
    )
)


