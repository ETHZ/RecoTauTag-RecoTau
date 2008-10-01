import FWCore.ParameterSet.Config as cms


process = cms.Process("TAU")


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(500)
)
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
'rfio:/castor/cern.ch/user/p/pjanot/CMSSW219/reco_QCDpt30_50_Full.root'

    )
)

process.load("Configuration.StandardSequences.GeometryPilot2_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")


# Conditions: fake or frontier
# process.load("Configuration.StandardSequences.FakeConditions_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'IDEAL_V9::All'

process.DQMStore = cms.Service("DQMStore")

process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")

process.load("RecoTauTag.RecoTau.PFRecoTauDiscriminationAgainstElectron_cfi")

process.load("RecoTauTag.RecoTau.PFRecoTauDiscriminationAgainstMuon_cfi")

process.pfTaus = cms.EDFilter("PFTauSelector",
    src = cms.InputTag("pfRecoTauProducer"),
    discriminator = cms.InputTag("pfRecoTauDiscriminationByIsolation")
)
process.pfTauTag = cms.EDAnalyzer("PFTauTagVal",
    OutPutFile = cms.string('pftautag.root'), ## This name is modified to reflect releaseversion and histograms stored
    PFTauProducer = cms.string('pfRecoTauProducer'),
    DataType = cms.string('QCD'),
    OutPutHistograms = cms.string('OneProngAndThreeProng'),
    ExtensionName = cms.InputTag("PFTauIsolationValidation"),
    PFTauDiscriminatorAgainstElectronProducer = cms.string('pfRecoTauDiscriminationAgainstElectron'),
    PFTauDiscriminatorAgainstMuonProducer = cms.string('pfRecoTauDiscriminationAgainstMuon'),
    GenJetProd = cms.InputTag("iterativeCone5GenJets")
)
process.p1 = cms.Path(
    process.PFTau +
    process.pfRecoTauDiscriminationAgainstElectron +
    process.pfRecoTauDiscriminationAgainstMuon +
    process.pfTauTag+
    process.pfTaus
    )


process.load("Configuration.EventContent.EventContent_cff")
process.aod = cms.OutputModule("PoolOutputModule",
    process.AODSIMEventContent,
    fileName = cms.untracked.string('/tmp/gennai/aodFullSimBJets.root')
)
process.aod.outputCommands.append('keep edmHepMCProduct_*_*_*')
process.aod.outputCommands.append('keep recoPFTaus_*_*_*')

#process.outpath = cms.EndPath(process.aod)

#
