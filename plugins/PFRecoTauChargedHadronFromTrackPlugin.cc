/*
 * PFRecoTauChargedHadronFromTrackPlugin
 *
 * Build PFRecoTauChargedHadron objects
 * using charged PFCandidates as input
 *
 * Author: Christian Veelken, LLR
 *
 * $Id $
 */

#include "RecoTauTag/RecoTau/interface/PFRecoTauChargedHadronPlugins.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TauReco/interface/PFRecoTauChargedHadron.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "FastSimulation/BaseParticlePropagator/interface/BaseParticlePropagator.h"
#include "FastSimulation/Particle/interface/RawParticle.h"

#include "RecoTauTag/RecoTau/interface/RecoTauCommonUtilities.h"
#include "RecoTauTag/RecoTau/interface/RecoTauQualityCuts.h"
#include "RecoTauTag/RecoTau/interface/RecoTauVertexAssociator.h"

#include <memory>

namespace reco { 

namespace tau {

class PFRecoTauChargedHadronFromTrackPlugin : public PFRecoTauChargedHadronBuilderPlugin 
{
 public:
  explicit PFRecoTauChargedHadronFromTrackPlugin(const edm::ParameterSet&);
  virtual ~PFRecoTauChargedHadronFromTrackPlugin();
  // Return type is auto_ptr<ChargedHadronVector>
  return_type operator()(const reco::PFJet&) const;
  // Hook to update PV information
  virtual void beginEvent();
  
 private:
  typedef std::vector<reco::PFCandidatePtr> PFCandPtrs;

  RecoTauVertexAssociator vertexAssociator_;

  RecoTauQualityCuts* qcuts_;

  edm::InputTag srcTracks_;
  double dRcone_;

  math::XYZVector magneticFieldStrength_;
};

PFRecoTauChargedHadronFromTrackPlugin::PFRecoTauChargedHadronFromTrackPlugin(const edm::ParameterSet& pset)
  : PFRecoTauChargedHadronBuilderPlugin(pset),
    vertexAssociator_(pset.getParameter<edm::ParameterSet>("qualityCuts")),
    qcuts_(0)
{
  edm::ParameterSet qcuts_pset = pset.getParameterSet("qualityCuts").getParameterSet("signalQualityCuts");
  qcuts_ = new RecoTauQualityCuts(qcuts_pset);

  srcTracks_ = pset.getParameter<edm::InputTag>("srcTracks");
  dRcone_ = pset.getParameter<double>("dRcone");
}
  
PFRecoTauChargedHadronFromTrackPlugin::~PFRecoTauChargedHadronFromTrackPlugin()
{
  delete qcuts_;
}

// Update the primary vertex
void PFRecoTauChargedHadronFromTrackPlugin::beginEvent() 
{
  vertexAssociator_.setEvent(*this->evt());

  edm::ESHandle<MagneticField> magneticField;
  evtSetup()->get<IdealMagneticFieldRecord>().get(magneticField);
  magneticFieldStrength_ = magneticField->inTesla(GlobalPoint(0.,0.,0.));
}

PFRecoTauChargedHadronFromTrackPlugin::return_type PFRecoTauChargedHadronFromTrackPlugin::operator()(const reco::PFJet& jet) const 
{
  ChargedHadronVector output;

  const edm::Event& evt = (*this->evt());

  edm::Handle<reco::TrackCollection> tracks;
  evt.getByLabel(srcTracks_, tracks);

  qcuts_->setPV(vertexAssociator_.associatedVertex(jet));

  size_t numTracks = tracks->size();
  for ( size_t iTrack = 0; iTrack < numTracks; ++iTrack ) {
    reco::TrackRef track(tracks, iTrack);
    
    // consider tracks in vicinity of tau-jet candidate only
    double dR = deltaR(track->eta(), track->phi(), jet.eta(), jet.phi());
    if ( dR > dRcone_ ) continue;
    
    // ignore tracks which fail quality cuts
    if ( !qcuts_->filterTrack(track) ) continue;

    reco::Candidate::Charge trackCharge_int = 0;
    if ( track->charge() > 0. ) trackCharge_int = +1;
    else if ( track->charge() < 0. ) trackCharge_int = -1;

    const double chargedPionMass = 0.13957; // GeV
    reco::Candidate::PolarLorentzVector p4_polar(track->pt(), track->eta(), track->phi(), chargedPionMass);
    reco::Candidate::LorentzVector p4(p4_polar.px(), p4_polar.py(), p4_polar.pz(), p4_polar.E());

    reco::Vertex::Point vtx(0.,0.,0.);
    if ( vertexAssociator_.associatedVertex(jet).isNonnull() ) vtx = vertexAssociator_.associatedVertex(jet)->position();

    std::auto_ptr<PFRecoTauChargedHadron> chargedHadron(new PFRecoTauChargedHadron(trackCharge_int, p4, vtx, 0, true, PFRecoTauChargedHadron::kTrack));
    chargedHadron->track_ = edm::Ptr<reco::Track>(tracks, iTrack);

    // CV: take code for propagating track to ECAL entrance 
    //     from RecoParticleFlow/PFTracking/src/PFTrackTransformer.cc
    //     to make sure propagation is done in the same way as for charged PFCandidates
    double chargedPionEn_out = sqrt(chargedPionMass*chargedPionMass + track->outerMomentum().Mag2());
    reco::Candidate::LorentzVector chargedPionP4_out(track->outerMomentum().x(), track->outerMomentum().y(), track->outerMomentum().z(), chargedPionEn_out);
    XYZTLorentzVector chargedPionPos_out(track->outerPosition().x(), track->outerPosition().y(), track->outerPosition().z(), 0.);
    BaseParticlePropagator trackPropagator(RawParticle(chargedPionP4_out, chargedPionPos_out), 0., 0., magneticFieldStrength_.z());
    trackPropagator.setCharge(track->charge());
    trackPropagator.propagateToEcalEntrance(false);
    if ( trackPropagator.getSuccess() != 0 ) { 
      chargedHadron->positionAtECALEntrance_ = trackPropagator.vertex();
    } else {
      edm::LogWarning("PFRecoTauChargedHadronFromTrackPlugin::operator()") 
	<< "Failed to propagate track: Pt = " << track->pt() << ", eta = " << track->eta() << ", phi = " << track->phi() << " to ECAL entrance !!" << std::endl;
      chargedHadron->positionAtECALEntrance_ = math::XYZPointF(0.,0.,0.);
    }

    output.push_back(chargedHadron);
  }

  return output.release();
}

}} // end namespace reco::tau

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(PFRecoTauChargedHadronBuilderPluginFactory, reco::tau::PFRecoTauChargedHadronFromTrackPlugin, "PFRecoTauChargedHadronFromTrackPlugin");
