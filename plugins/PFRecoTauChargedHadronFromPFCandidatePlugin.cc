/*
 * PFRecoTauChargedHadronFromPFCandidatePlugin
 *
 * Build PFRecoTauChargedHadron objects
 * using charged PFCandidates and/or PFNeutralHadrons as input
 *
 * Author: Christian Veelken, LLR
 *
 * $Id $
 */

#include "RecoTauTag/RecoTau/interface/PFRecoTauChargedHadronPlugins.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TauReco/interface/PFRecoTauChargedHadron.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "RecoTauTag/RecoTau/interface/RecoTauCommonUtilities.h"
#include "RecoTauTag/RecoTau/interface/RecoTauQualityCuts.h"
#include "RecoTauTag/RecoTau/interface/RecoTauVertexAssociator.h"

#include <memory>

namespace reco 
{ 

namespace tau
{ 

class PFRecoTauChargedHadronFromPFCandidatePlugin : public PFRecoTauChargedHadronBuilderPlugin 
{
 public:
  explicit PFRecoTauChargedHadronFromPFCandidatePlugin(const edm::ParameterSet&);
  virtual ~PFRecoTauChargedHadronFromPFCandidatePlugin();
  // Return type is auto_ptr<ChargedHadronVector>
  return_type operator()(const reco::PFJet&) const;
  // Hook to update PV information
  virtual void beginEvent();
  
 private:
  typedef std::vector<reco::PFCandidatePtr> PFCandPtrs;

  RecoTauVertexAssociator vertexAssociator_;

  RecoTauQualityCuts* qcuts_;

  std::vector<int> inputPdgIds_;  // type of candidates to clusterize

  double dRmergeNeutralHadron_;
  double dRmergePhoton_;
};

PFRecoTauChargedHadronFromPFCandidatePlugin::PFRecoTauChargedHadronFromPFCandidatePlugin(const edm::ParameterSet& pset)
  : PFRecoTauChargedHadronBuilderPlugin(pset),
    vertexAssociator_(pset.getParameter<edm::ParameterSet>("qualityCuts")),
    qcuts_(0)
{
  edm::ParameterSet qcuts_pset = pset.getParameterSet("qualityCuts").getParameterSet("signalQualityCuts");
  qcuts_ = new RecoTauQualityCuts(qcuts_pset);

  inputPdgIds_ = pset.getParameter<std::vector<int> >("chargedHadronCandidatesParticleIds");

  dRmergeNeutralHadron_ = pset.getParameter<double>("dRmergeNeutralHadron");
  dRmergePhoton_ = pset.getParameter<double>("dRmergePhoton");
}
  
PFRecoTauChargedHadronFromPFCandidatePlugin::~PFRecoTauChargedHadronFromPFCandidatePlugin()
{
  delete qcuts_;
}

// Update the primary vertex
void PFRecoTauChargedHadronFromPFCandidatePlugin::beginEvent() 
{
  vertexAssociator_.setEvent(*evt());
}

PFRecoTauChargedHadronFromPFCandidatePlugin::return_type PFRecoTauChargedHadronFromPFCandidatePlugin::operator()(const reco::PFJet& jet) const 
{
  ChargedHadronVector output;

  // Get the candidates passing our quality cuts
  qcuts_->setPV(vertexAssociator_.associatedVertex(jet));
  PFCandPtrs candsVector = qcuts_->filterCandRefs(pfCandidates(jet, inputPdgIds_));

  for ( PFCandPtrs::iterator cand = candsVector.begin();
	cand != candsVector.end(); ++cand ) {
    PFRecoTauChargedHadron::PFRecoTauChargedHadronAlgorithm algo = PFRecoTauChargedHadron::kUndefined;
    if ( fabs((*cand)->charge()) > 0.5 ) algo = PFRecoTauChargedHadron::kChargedPFCandidate;
    else algo = PFRecoTauChargedHadron::kPFNeutralHadron;

    std::auto_ptr<PFRecoTauChargedHadron> chargedHadron(new PFRecoTauChargedHadron(**cand, algo));

    chargedHadron->chargedPFCandidate_ = (*cand);

    chargedHadron->positionAtECALEntrance_ = (*cand)->positionAtECALEntrance();

    std::vector<reco::PFCandidatePtr> jetConstituents = jet.getPFConstituents();
    for ( std::vector<reco::PFCandidatePtr>::const_iterator jetConstituent = jetConstituents.begin();
	  jetConstituent != jetConstituents.end(); ++jetConstituent ) {
      // CV: take care of not double-counting energy in case "charged" PFCandidate is in fact a PFNeutralHadron
      if ( (*jetConstituent) == chargedHadron->chargedPFCandidate_ ) continue;

      reco::PFCandidate::ParticleType jetConstituentType = (*jetConstituent)->particleId();
      if ( !(jetConstituentType == reco::PFCandidate::h0 || jetConstituentType == reco::PFCandidate::gamma) ) continue;

      double dR = deltaR((*jetConstituent)->positionAtECALEntrance(), chargedHadron->positionAtECALEntrance_);
      double dRmerge = -1.;      
      if      ( jetConstituentType == reco::PFCandidate::h0    ) dRmerge = dRmergeNeutralHadron_;
      else if ( jetConstituentType == reco::PFCandidate::gamma ) dRmerge = dRmergePhoton_;
      if ( dR < dRmerge ) chargedHadron->neutralPFCandidates_.push_back(*jetConstituent);
    }

    output.push_back(chargedHadron);
  }

  return output.release();
}

}} // end namespace reco::tau

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(PFRecoTauChargedHadronBuilderPluginFactory, reco::tau::PFRecoTauChargedHadronFromPFCandidatePlugin, "PFRecoTauChargedHadronFromPFCandidatePlugin");
