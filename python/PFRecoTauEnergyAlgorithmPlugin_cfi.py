import FWCore.ParameterSet.Config as cms

pfTauEnergyAlgorithmPlugin = cms.PSet(
    dRaddNeutralHadron = cms.double(-1.), # CV: disabled adding PFNeutralHadrons
        minNeutralHadronEt = cms.double(20.),
    dRaddPhoton = cms.double(-1.), # CV: disabled adding PFGammas
    minGammaEt = cms.double(5.)
)
