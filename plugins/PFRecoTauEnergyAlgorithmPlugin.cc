/*
 * =============================================================================
 *       Filename:  RecoTauEnergyAlgorithmPlugin.cc
 *
 *    Description:  Determine best estimate for tau energy
 *                  for tau candidates reconstructed in different decay modes
 *
 *        Created:  04/09/2013 11:40:00
 *
 *         Authors:  Christian Veelken (LLR)
 *
 * =============================================================================
 */

#include "RecoTauTag/RecoTau/interface/RecoTauBuilderPlugins.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/TauReco/interface/PFRecoTauChargedHadron.h"
#include "DataFormats/TauReco/interface/PFRecoTauChargedHadronFwd.h"
#include "RecoTauTag/RecoTau/interface/pfRecoTauChargedHadronAuxFunctions.h"

#include "DataFormats/Math/interface/deltaR.h"

#include <vector>

namespace reco { namespace tau {

class PFRecoTauEnergyAlgorithmPlugin : public RecoTauModifierPlugin
{
 public:

  explicit PFRecoTauEnergyAlgorithmPlugin(const edm::ParameterSet&);
  virtual ~PFRecoTauEnergyAlgorithmPlugin();
  void operator()(PFTau&) const;
  virtual void beginJob(edm::EDProducer*);
  virtual void beginEvent();
  virtual void endEvent();

 private:
  
  double dRaddNeutralHadron_;
  double minNeutralHadronEt_;
  double dRaddPhoton_;
  double minGammaEt_;

  int verbosity_;
};

PFRecoTauEnergyAlgorithmPlugin::PFRecoTauEnergyAlgorithmPlugin(const edm::ParameterSet& cfg)
  : RecoTauModifierPlugin(cfg),
    dRaddNeutralHadron_(cfg.getParameter<double>("dRaddNeutralHadron")),
    minNeutralHadronEt_(cfg.getParameter<double>("minNeutralHadronEt")),
    dRaddPhoton_(cfg.getParameter<double>("dRaddPhoton")),
    minGammaEt_(cfg.getParameter<double>("minGammaEt"))
{
  verbosity_ = ( cfg.exists("verbosity") ) ?
    cfg.getParameter<int>("verbosity") : 0;
}

PFRecoTauEnergyAlgorithmPlugin::~PFRecoTauEnergyAlgorithmPlugin()
{}

void PFRecoTauEnergyAlgorithmPlugin::beginJob(edm::EDProducer* producer)
{}

void PFRecoTauEnergyAlgorithmPlugin::beginEvent()
{}

namespace
{
  double getTrackPerr2(const reco::Track& track)
  {
    double trackPerr = track.p()*(track.ptError()/track.pt());
    return trackPerr*trackPerr;
  }

  void updateTauP4(PFTau& tau, double sf, const reco::Candidate::LorentzVector& addP4)
  {
    // preserve tau candidate mass when adding extra neutral energy
    double tauPx_modified = tau.px() + sf*addP4.px();
    double tauPy_modified = tau.py() + sf*addP4.py();
    double tauPz_modified = tau.pz() + sf*addP4.pz();
    double tauMass = tau.mass();
    double tauEn_modified = sqrt(tauPx_modified*tauPx_modified + tauPy_modified*tauPy_modified + tauPz_modified*tauPz_modified + tauMass*tauMass);
    reco::Candidate::LorentzVector tauP4_modified(tauPx_modified, tauPy_modified, tauPz_modified, tauEn_modified);
    tau.setP4(tauP4_modified);
  }
}

void PFRecoTauEnergyAlgorithmPlugin::operator()(PFTau& tau) const
{
  if ( verbosity_ ) {
    std::cout << "<PFRecoTauEnergyAlgorithmPlugin::operator()>:" << std::endl;
  }

  // Add high Pt PFNeutralHadrons and PFGammas that are not "used" by tau decay mode object
  std::vector<reco::PFCandidatePtr> addNeutrals;
  reco::Candidate::LorentzVector addNeutralsSumP4;
  std::vector<reco::PFCandidatePtr> jetConstituents = tau.jetRef()->getPFConstituents();
  for ( std::vector<reco::PFCandidatePtr>::const_iterator jetConstituent = jetConstituents.begin();
	jetConstituent != jetConstituents.end(); ++jetConstituent ) {
    reco::PFCandidate::ParticleType jetConstituentType = (*jetConstituent)->particleId();
    if ( !((jetConstituentType == reco::PFCandidate::h0    && (*jetConstituent)->et() > minNeutralHadronEt_) ||
	   (jetConstituentType == reco::PFCandidate::gamma && (*jetConstituent)->et() > minGammaEt_        )) ) continue;

    bool isSignalPFCand = false;
    const std::vector<reco::PFCandidatePtr>& signalPFCands = tau.signalPFCands();
    for ( std::vector<reco::PFCandidatePtr>::const_iterator signalPFCand = signalPFCands.begin();
	  signalPFCand != signalPFCands.end(); ++signalPFCand ) {
      if ( (*jetConstituent) == (*signalPFCand) ) isSignalPFCand = true;
    }
    if ( isSignalPFCand ) continue;
    
    double dR = deltaR((*jetConstituent)->p4(), tau.p4());
    double dRadd = -1.;      
    if      ( jetConstituentType == reco::PFCandidate::h0    ) dRadd = dRaddNeutralHadron_;
    else if ( jetConstituentType == reco::PFCandidate::gamma ) dRadd = dRaddPhoton_;
    if ( dR < dRadd ) {
      addNeutrals.push_back(*jetConstituent);
      addNeutralsSumP4 += (*jetConstituent)->p4();
    }
  }
  
  unsigned numChargedHadronTracks = 0;
  double chargedHadronTracksSumP = 0.;
  double chargedHadronTracksSumPerr2 = 0.;
  const std::vector<PFRecoTauChargedHadron> chargedHadrons = tau.signalTauChargedHadronCandidates();
  for ( std::vector<PFRecoTauChargedHadron>::const_iterator chargedHadron = chargedHadrons.begin();
	chargedHadron != chargedHadrons.end(); ++chargedHadron ) {
    if ( chargedHadron->algoIs(PFRecoTauChargedHadron::kTrack) ) {
      ++numChargedHadronTracks;
      const edm::Ptr<Track>& chargedHadronTrack = chargedHadron->getTrack();
      chargedHadronTracksSumP += chargedHadronTrack->p();
      chargedHadronTracksSumPerr2 += getTrackPerr2(*chargedHadronTrack);
    }
  }
  if ( verbosity_ ) {
    std::cout << "numChargedHadronTracks = " << numChargedHadronTracks 
	      << " (chargedHadronTracksSumP = " << chargedHadronTracksSumP << " +/- " << sqrt(chargedHadronTracksSumPerr2) << ")" << std::endl;
    std::cout << "addNeutralsSumP4: En = " << addNeutralsSumP4.energy() << std::endl;
  }
 
  if ( numChargedHadronTracks == 0 ) {
    // This is the easy case: 
    // All tau energy is taken from PFCandidates reconstructed by PFlow algorithm
    // and there is no issue with double-counting of energy.
    updateTauP4(tau, 1., addNeutralsSumP4);
    return;
  } else {
    // This is the difficult case: 
    // The tau energy needs to be computed for an arbitrary mix of charged and neutral PFCandidates plus reco::Tracks.
    // We need to make sure not to double-count energy deposited by reco::Track in ECAL and/or HCAL as neutral PFCandidates.
    
    // Check if we have enough energy in collection of PFNeutralHadrons and PFGammas that are not "used" by tau decay mode object
    // to balance track momenta:
    if ( chargedHadronTracksSumP < addNeutralsSumP4.energy() ) {
      double scaleFactor = 1. - chargedHadronTracksSumP/addNeutralsSumP4.energy();
      assert(scaleFactor >= 0. && scaleFactor <= 1.);
      updateTauP4(tau, scaleFactor, addNeutralsSumP4);
      return;
    }

    // Determine which neutral PFCandidates are close to PFChargedHadrons
    // and have been merged into ChargedHadrons
    std::vector<reco::PFCandidatePtr> mergedNeutrals;
    reco::Candidate::LorentzVector mergedNeutralsSumP4;
    for ( std::vector<PFRecoTauChargedHadron>::const_iterator chargedHadron = chargedHadrons.begin();
	  chargedHadron != chargedHadrons.end(); ++chargedHadron ) {
      const std::vector<reco::PFCandidatePtr>& neutralPFCands = chargedHadron->getNeutralPFCandidates();
      for ( std::vector<reco::PFCandidatePtr>::const_iterator neutralPFCand = neutralPFCands.begin();
	    neutralPFCand != neutralPFCands.end(); ++neutralPFCand ) {
	mergedNeutrals.push_back(*neutralPFCand);
	mergedNeutralsSumP4 += (*neutralPFCand)->p4();
      }
    }
    if ( verbosity_ ) {
      std::cout << "mergedNeutralsSumP4: En = " << mergedNeutralsSumP4.energy() << std::endl;
    }

    // Check if track momenta are balanced by sum of PFNeutralHadrons and PFGammas that are not "used" by tau decay mode object
    // plus neutral PFCandidates close to PFChargedHadrons:
    if ( chargedHadronTracksSumP < (addNeutralsSumP4.energy() + mergedNeutralsSumP4.energy()) ) {
      double scaleFactor = ((addNeutralsSumP4.energy() + mergedNeutralsSumP4.energy()) - chargedHadronTracksSumP)/mergedNeutralsSumP4.energy();
      assert(scaleFactor >= 0. && scaleFactor <= 1.);
      reco::Candidate::LorentzVector diffP4;
      for ( std::vector<PFRecoTauChargedHadron>::const_iterator chargedHadron = chargedHadrons.begin();
	    chargedHadron != chargedHadrons.end(); ++chargedHadron ) {
	if ( chargedHadron->getNeutralPFCandidates().size() >= 1 ) {
	  const PFRecoTauChargedHadron* chargedHadron_modified = &(*chargedHadron);
	  setChargedHadronP4(*(const_cast<PFRecoTauChargedHadron*>(chargedHadron_modified)), scaleFactor);
	  diffP4 += (chargedHadron->p4() - chargedHadron_modified->p4());
	}
      }
      updateTauP4(tau, -1., diffP4);
      return;
    }

    // Determine energy sum of all PFNeutralHadrons interpreted as ChargedHadrons with missing track
    unsigned numChargedHadronNeutrals = 0;
    std::vector<reco::PFCandidatePtr> chargedHadronNeutrals;
    reco::Candidate::LorentzVector chargedHadronNeutralsSumP4;
    for ( std::vector<PFRecoTauChargedHadron>::const_iterator chargedHadron = chargedHadrons.begin();
	  chargedHadron != chargedHadrons.end(); ++chargedHadron ) {
      if ( chargedHadron->algoIs(PFRecoTauChargedHadron::kPFNeutralHadron) ) {
	++numChargedHadronNeutrals;
	chargedHadronNeutrals.push_back(chargedHadron->getChargedPFCandidate());
	chargedHadronNeutralsSumP4 += chargedHadron->getChargedPFCandidate()->p4();
      }
    }
    if ( verbosity_ ) {
      std::cout << "numChargedHadronNeutrals = " << numChargedHadronNeutrals << std::endl;
      std::cout << "chargedHadronNeutralsSumP4: En = " << chargedHadronNeutralsSumP4.energy() << std::endl;
    }
    
    // Check if sum of PFNeutralHadrons and PFGammas that are not "used" by tau decay mode object
    // plus neutral PFCandidates close to PFChargedHadrons plus PFNeutralHadrons interpreted as ChargedHadrons with missing track balances track momenta
    if ( chargedHadronTracksSumP < (addNeutralsSumP4.energy() + mergedNeutralsSumP4.energy() + chargedHadronNeutralsSumP4.energy()) ) {
      double scaleFactor = ((addNeutralsSumP4.energy() + mergedNeutralsSumP4.energy() + chargedHadronNeutralsSumP4.energy()) - chargedHadronTracksSumP)/chargedHadronNeutralsSumP4.energy();
      assert(scaleFactor >= 0. && scaleFactor <= 1.);
      reco::Candidate::LorentzVector diffP4;
      for ( std::vector<PFRecoTauChargedHadron>::const_iterator chargedHadron = chargedHadrons.begin();
	    chargedHadron != chargedHadrons.end(); ++chargedHadron ) {
	 const PFRecoTauChargedHadron* chargedHadron_modified = &(*chargedHadron);
	(const_cast<PFRecoTauChargedHadron*>(chargedHadron_modified))->neutralPFCandidates_.clear();
	if ( chargedHadron->algoIs(PFRecoTauChargedHadron::kPFNeutralHadron) ) {
	  const PFCandidatePtr& chargedPFCand = chargedHadron->getChargedPFCandidate();
	  double chargedHadronPx_modified = scaleFactor*chargedPFCand->px();
	  double chargedHadronPy_modified = scaleFactor*chargedPFCand->py();
	  double chargedHadronPz_modified = scaleFactor*chargedPFCand->pz();	
	  reco::Candidate::LorentzVector chargedHadronP4_modified = compChargedHadronP4(chargedHadronPx_modified, chargedHadronPy_modified, chargedHadronPz_modified);
	  (const_cast<PFRecoTauChargedHadron*>(chargedHadron_modified))->setP4(chargedHadronP4_modified);
	  diffP4 += (chargedHadron->p4() - chargedHadron_modified->p4());
	}
      }
      updateTauP4(tau, -1., diffP4);
      return;
    } else {
      if ( numChargedHadronNeutrals == 0 ) {
	// Adjust momenta of ChargedHadrons build from reco::Tracks to match sum of energy deposits in ECAL + HCAL
	reco::Candidate::LorentzVector diffP4;
	for ( std::vector<PFRecoTauChargedHadron>::const_iterator chargedHadron = chargedHadrons.begin();
	      chargedHadron != chargedHadrons.end(); ++chargedHadron ) {
	  if ( chargedHadron->algoIs(PFRecoTauChargedHadron::kTrack) ) {
	    const PFRecoTauChargedHadron* chargedHadron_modified = &(*chargedHadron);
	    (const_cast<PFRecoTauChargedHadron*>(chargedHadron_modified))->neutralPFCandidates_.clear();
	    const edm::Ptr<Track>& track = chargedHadron->getTrack();
	    double trackP = track->p();
	    double trackPerr2 = getTrackPerr2(*track);	  
	    //std::cout << "trackP = " << trackP << " +/- " << sqrt(trackPerr2) << std::endl;
	    //std::cout << "chargedHadronTracksSumPerr2 = " << chargedHadronTracksSumPerr2 << std::endl;
	    //std::cout << "trackPerr2 = " << trackPerr2 << std::endl;
	    //std::cout << "addNeutralsSumP4.energy = " << addNeutralsSumP4.energy() << std::endl;
	    //std::cout << "mergedNeutralsSumP4.energy = " << mergedNeutralsSumP4.energy() << std::endl;
	    //std::cout << "chargedHadronNeutralsSumP4.energy = " << chargedHadronNeutralsSumP4.energy() << std::endl;
	    //std::cout << "chargedHadronTracksSumP = " << chargedHadronTracksSumP << std::endl;	    
	    // CV: adjust track momenta such that difference beeen (measuredTrackP - adjustedTrackP)/sigmaMeasuredTrackP is minimal
	    //    (expression derived using Mathematica)
	    double trackP_modified = 
              (trackP*(chargedHadronTracksSumPerr2 - trackPerr2) 
	     + trackPerr2*(addNeutralsSumP4.energy() + mergedNeutralsSumP4.energy() + chargedHadronNeutralsSumP4.energy() - (chargedHadronTracksSumP - trackP)))/
              chargedHadronTracksSumPerr2;
	    // CV: trackP_modified may actually become negative in case sum of energy deposits in ECAL + HCAL is small
	    //     and one of the tracks has a significantly larger momentum uncertainty than the other tracks.
	    //     In this case set track momentum to small positive value.
	    if ( trackP_modified < 1.e-1 ) trackP_modified = 1.e-1;
	    //std::cout << "trackP (modified) = " << trackP_modified << std::endl;
	    double scaleFactor = trackP_modified/trackP;
	    assert(scaleFactor >= 0. && scaleFactor <= 1.);
	    double chargedHadronPx_modified = scaleFactor*track->px();
	    double chargedHadronPy_modified = scaleFactor*track->py();
	    double chargedHadronPz_modified = scaleFactor*track->pz();
	    reco::Candidate::LorentzVector chargedHadronP4_modified = compChargedHadronP4(chargedHadronPx_modified, chargedHadronPy_modified, chargedHadronPz_modified);
	    (const_cast<PFRecoTauChargedHadron*>(chargedHadron_modified))->setP4(chargedHadronP4_modified);
	    diffP4 += (chargedHadron->p4() - chargedHadron_modified->p4());
	  }
	}
	updateTauP4(tau, -1., diffP4);
	return;
      } else {
	// Interpretation of PFNeutralHadrons as ChargedHadrons with missing track is not compatible 
	// with the fact that sum of reco::Track momenta exceeds sum of energy deposits in ECAL + HCAL:
	// kill tau candidate (by setting its four-vector to zero)
	reco::Candidate::LorentzVector tauP4_modified(0.,0.,0.,0.);
	tau.setP4(tauP4_modified);
	tau.setStatus(-1);
	return;
      }
    }
  }

  // CV: You should never come here.
  assert(0);
}

void PFRecoTauEnergyAlgorithmPlugin::endEvent()
{}

}} // end namespace reco::tau

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(RecoTauModifierPluginFactory, reco::tau::PFRecoTauEnergyAlgorithmPlugin, "PFRecoTauEnergyAlgorithmPlugin");
