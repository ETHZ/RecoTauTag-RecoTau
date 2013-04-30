#include <vector>

#include "RecoTauTag/RecoTau/interface/RecoTauBuilderPlugins.h"
#include "RecoTauTag/RecoTau/interface/RecoTauCommonUtilities.h"

#include "RecoTauTag/RecoTau/interface/CombinatoricGenerator.h"
#include "RecoTauTag/RecoTau/interface/RecoTauCrossCleaning.h"
#include "RecoTauTag/RecoTau/interface/ConeTools.h"

#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFRecoTauChargedHadron.h"
#include "DataFormats/TauReco/interface/RecoTauPiZero.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "RecoTauTag/RecoTau/interface/RecoTauConstructor.h"
#include "RecoTauTag/RecoTau/interface/RecoTauQualityCuts.h"

namespace reco { namespace tau {

typedef std::vector<reco::PFRecoTauChargedHadron> ChargedHadronList;
typedef std::vector<RecoTauPiZero> PiZeroList;
typedef tau::CombinatoricGenerator<ChargedHadronList> PFCombo;

class RecoTauBuilderCombinatoricPlugin : public RecoTauBuilderPlugin 
{
 public:
  explicit RecoTauBuilderCombinatoricPlugin(const edm::ParameterSet&);
  virtual ~RecoTauBuilderCombinatoricPlugin() {}

  return_type operator()(
      const reco::PFJetRef&, 
      const std::vector<reco::PFRecoTauChargedHadron>&, 
      const std::vector<RecoTauPiZero>&, 
      const std::vector<PFCandidatePtr>&) const;

 private:
  RecoTauQualityCuts qcuts_;

  double isolationConeSize_;

  struct decayModeInfo 
  {
    uint32_t maxPiZeros_;
    uint32_t maxPFCHs_;
    uint32_t nCharged_;
    uint32_t nPiZeros_;
  };
  std::vector<decayModeInfo> decayModesToBuild_;

  int verbosity_;
};

RecoTauBuilderCombinatoricPlugin::RecoTauBuilderCombinatoricPlugin(const edm::ParameterSet& pset)
  : RecoTauBuilderPlugin(pset),
    qcuts_(pset.getParameterSet("qualityCuts").getParameterSet("signalQualityCuts")),
    isolationConeSize_(pset.getParameter<double>("isolationConeSize")) 
{
  typedef std::vector<edm::ParameterSet> VPSet;
  const VPSet& decayModes = pset.getParameter<VPSet>("decayModes");
  for ( VPSet::const_iterator decayMode = decayModes.begin();
	decayMode != decayModes.end(); ++decayMode ) {
    decayModeInfo info;
    info.nCharged_ = decayMode->getParameter<uint32_t>("nCharged");
    info.nPiZeros_ = decayMode->getParameter<uint32_t>("nPiZeros");
    info.maxPFCHs_ = decayMode->getParameter<uint32_t>("maxTracks");
    info.maxPiZeros_ = decayMode->getParameter<uint32_t>("maxPiZeros");
    decayModesToBuild_.push_back(info);
    
  }
  
  verbosity_ = ( pset.exists("verbosity") ) ?
    pset.getParameter<int>("verbosity") : 0;
}

// define template specialization for cross-cleaning 
namespace xclean
{
  template<>
  inline void CrossCleanPiZeros<PFCombo::combo_iterator>::initialize(PFCombo::combo_iterator signalTracksBegin, PFCombo::combo_iterator signalTracksEnd) 
  {
    // Get the list of objects we need to clean
    for ( PFCombo::combo_iterator i = signalTracksBegin; i != signalTracksEnd; ++i ) {
      const reco::CompositePtrCandidate::daughters& daughters = i->daughterPtrVector();
      for ( reco::CompositePtrCandidate::daughters::const_iterator daughter = daughters.begin();
	    daughter != daughters.end(); ++daughter ) {
	toRemove_.insert(reco::CandidatePtr(*daughter));
      }
    }
  }

  template<>
  inline void CrossCleanPtrs<PiZeroList>::initialize(const PiZeroList& piZeros) 
  {
    BOOST_FOREACH( const PFCandidatePtr &ptr, flattenPiZeros(piZeros) ) {
      toRemove_.insert(CandidatePtr(ptr));
    }
  }

  template<>
  inline void CrossCleanPtrs<ChargedHadronList>::initialize(const ChargedHadronList& chargedHadrons) 
  {
    for ( ChargedHadronList::const_iterator i = chargedHadrons.begin(); i != chargedHadrons.end(); ++i ) {
      const reco::CompositePtrCandidate::daughters& daughters = i->daughterPtrVector();
      for ( reco::CompositePtrCandidate::daughters::const_iterator daughter = daughters.begin();
	    daughter != daughters.end(); ++daughter ) {
	toRemove_.insert(reco::CandidatePtr(*daughter));
      }
    }
  }
}

namespace
{
  // auxiliary class for sorting pizeros by descending transverse momentum
  class SortPi0sDescendingPt 
  {
   public:
    bool operator()(const RecoTauPiZero& a, const RecoTauPiZero& b) const 
    {
      return a.pt() > b.pt();
    }  
  };
}

RecoTauBuilderCombinatoricPlugin::return_type
RecoTauBuilderCombinatoricPlugin::operator()(
    const reco::PFJetRef& jet, 
    const std::vector<reco::PFRecoTauChargedHadron>& chargedHadrons, 
    const std::vector<RecoTauPiZero>& piZeros, 
    const std::vector<PFCandidatePtr>& regionalExtras) const 
{
  if ( verbosity_ ) {
    std::cout << "<RecoTauBuilderCombinatoricPlugin::operator()>:" << std::endl;
  }
  
  // Define output.  
  output_type output;
  
  reco::VertexRef primaryVertexRef = primaryVertex(jet);
  
  // Update the primary vertex used by the quality cuts.  The PV is supplied by
  // the base class.
  qcuts_.setPV(primaryVertexRef);
  
  typedef std::vector<PFCandidatePtr> PFCandPtrs;
  
  if ( verbosity_ ) {
    std::cout << "#chargedHadrons = " << chargedHadrons.size() << std::endl;
    int idx = 0;
    for ( ChargedHadronList::const_iterator chargedHadron = chargedHadrons.begin();
	  chargedHadron != chargedHadrons.end(); ++chargedHadron ) {
      std::cout << "chargedHadron #" << idx << ": Pt = " << chargedHadron->pt() << ", eta = " << chargedHadron->eta() << ", phi = " << chargedHadron->phi() << std::endl;
      ++idx;
    }
    std::cout << "#piZeros = " << piZeros.size() << std::endl;
    idx = 0;
    for ( PiZeroList::const_iterator piZero = piZeros.begin();
   	  piZero != piZeros.end(); ++piZero ) {
      std::cout << "piZero #" << idx << ": Pt = " << piZero->pt() << ", eta = " << piZero->eta() << ", phi = " << piZero->phi() << std::endl;
      ++idx;
    }
  }

  PFCandPtrs pfchs = qcuts_.filterCandRefs(pfChargedCands(*jet));
  PFCandPtrs pfnhs = qcuts_.filterCandRefs(pfCandidates(*jet, reco::PFCandidate::h0));
  
  /// Apply quality cuts to the regional junk around the jet.  Note that the
  /// particle contents of the junk is exclusive to the jet content.
  PFCandPtrs regionalJunk = qcuts_.filterCandRefs(regionalExtras);
    
  // Loop over the decay modes we want to build
  for ( std::vector<decayModeInfo>::const_iterator decayMode = decayModesToBuild_.begin();
	decayMode != decayModesToBuild_.end(); ++decayMode ) {
    // Find how many piZeros are in this decay mode
    size_t piZerosToBuild = decayMode->nPiZeros_;
    // Find how many tracks are in this decay mode
    size_t tracksToBuild = decayMode->nCharged_;
    if ( verbosity_ ) {
      std::cout << "piZerosToBuild = " << piZerosToBuild << std::endl;
      std::cout << "tracksToBuild = " << tracksToBuild << std::endl;
    }
    
    // Skip decay mode if jet doesn't have the multiplicity to support it
    if ( chargedHadrons.size() < tracksToBuild ) continue;

    // Find the start and end of potential signal tracks
    ChargedHadronList::const_iterator chargedHadron_begin = chargedHadrons.begin();
    ChargedHadronList::const_iterator chargedHadron_end = chargedHadrons.end();
    chargedHadron_end = takeNElements(chargedHadron_begin, chargedHadron_end, decayMode->maxPFCHs_);

    // Build our track combo generator
    PFCombo trackCombos(chargedHadron_begin, chargedHadron_end, tracksToBuild);

    PFCandPtrs::iterator pfch_end = pfchs.end();
    pfch_end = takeNElements(pfchs.begin(), pfch_end, decayMode->maxPFCHs_);

    //-------------------------------------------------------
    // Begin combinatoric loop for this decay mode
    //-------------------------------------------------------
    
    // Loop over the different combinations of tracks
    for ( PFCombo::iterator trackCombo = trackCombos.begin();
	  trackCombo != trackCombos.end(); ++trackCombo ) {
      xclean::CrossCleanPiZeros<PFCombo::combo_iterator> xCleaner(trackCombo->combo_begin(), trackCombo->combo_end());
      
      PiZeroList cleanPiZeros = xCleaner(piZeros);
      
      // CV: sort collection of cross-cleaned pi0s by descending Pt
      std::sort(cleanPiZeros.begin(), cleanPiZeros.end(), SortPi0sDescendingPt());
      
      // Skip decay mode if we don't have enough remaining clean pizeros to
      // build it.
      if ( cleanPiZeros.size() < piZerosToBuild ) continue;
      
      // Find the start and end of potential signal tracks
      PiZeroList::iterator piZero_begin = cleanPiZeros.begin();
      PiZeroList::iterator piZero_end = cleanPiZeros.end();
      piZero_end = takeNElements(piZero_begin, piZero_end, decayMode->maxPiZeros_);
      
      // Build our piZero combo generator
      typedef tau::CombinatoricGenerator<PiZeroList> PiZeroCombo;
      PiZeroCombo piZeroCombos(piZero_begin, piZero_end, piZerosToBuild);
      // Loop over the different combinations of PiZeros
      for ( PiZeroCombo::iterator piZeroCombo = piZeroCombos.begin();
            piZeroCombo != piZeroCombos.end(); ++piZeroCombo ) {
        // Output tau
        RecoTauConstructor tau(jet, getPFCands(), true);
        // Reserve space in our collections
        tau.reserve(
	    RecoTauConstructor::kSignal,
	    RecoTauConstructor::kChargedHadron, tracksToBuild);
        tau.reserve(
            RecoTauConstructor::kSignal,
            RecoTauConstructor::kGamma, 2*piZerosToBuild); // k-factor = 2
        tau.reservePiZero(RecoTauConstructor::kSignal, piZerosToBuild);

        // FIXME - are all these reserves okay?  will they get propagated to the
        // dataformat size if they are wrong?
        tau.reserve(
            RecoTauConstructor::kIsolation,
            RecoTauConstructor::kChargedHadron, chargedHadrons.size() - tracksToBuild);
        tau.reserve(
            RecoTauConstructor::kIsolation,
	    RecoTauConstructor::kGamma,
	    (cleanPiZeros.size() - piZerosToBuild)*2);
        tau.reservePiZero(
	    RecoTauConstructor::kIsolation,
	    (cleanPiZeros.size() - piZerosToBuild));

        // Get signal PiZero constituents and add them to the tau.
        // The sub-gammas are automatically added.
        tau.addPiZeros(
            RecoTauConstructor::kSignal,
            piZeroCombo->combo_begin(), piZeroCombo->combo_end());

	// Set signal and isolation components for charged hadrons, after
        // converting them to a PFCandidateRefVector
	//
	// NOTE: signal ChargedHadrons need to be added **after** signal PiZeros
	//       to avoid double-counting PFGammas as part of PiZero and merged with ChargedHadron
	//
        tau.addTauChargedHadrons(
            RecoTauConstructor::kSignal, 
            trackCombo->combo_begin(), trackCombo->combo_end());

        // Now build isolation collections
        // Load our isolation tools
        using namespace reco::tau::cone;
        PFCandPtrDRFilter isolationConeFilter(tau.p4(), 0, isolationConeSize_);

        // Cross cleaning predicate.  Remove any PFCandidatePtrs that are
        // contained within existing ChargedHadrons or PiZeros.  This predicate will return false
        // for any object that overlaps with chargedHadrons or cleanPiZeros.
	xclean::CrossCleanPtrs<PiZeroList> pfCandXCleaner_pizeros(cleanPiZeros);
	xclean::CrossCleanPtrs<ChargedHadronList> pfCandXCleaner_chargedHadrons(chargedHadrons);
	typedef xclean::PredicateAND<xclean::CrossCleanPtrs<PiZeroList>, xclean::CrossCleanPtrs<ChargedHadronList> > pfCandXCleanerType;
        pfCandXCleanerType pfCandXCleaner(pfCandXCleaner_pizeros, pfCandXCleaner_chargedHadrons);
        // And this cleaning filter predicate with our Iso cone filter
        xclean::PredicateAND<PFCandPtrDRFilter, pfCandXCleanerType> pfCandFilter(isolationConeFilter, pfCandXCleaner);

	ChargedHadronDRFilter isolationConeFilterChargedHadron(tau.p4(), 0, isolationConeSize_);
        PiZeroDRFilter isolationConeFilterPiZero(tau.p4(), 0, isolationConeSize_);

        // Additionally make predicates to select the different PF object types
        // of the regional junk objects to add
        typedef xclean::PredicateAND<xclean::FilterPFCandByParticleId,
	    PFCandPtrDRFilter> RegionalJunkConeAndIdFilter;

        xclean::FilterPFCandByParticleId
          pfchCandSelector(reco::PFCandidate::h);
        xclean::FilterPFCandByParticleId
          pfgammaCandSelector(reco::PFCandidate::gamma);
        xclean::FilterPFCandByParticleId
          pfnhCandSelector(reco::PFCandidate::h0);

        RegionalJunkConeAndIdFilter pfChargedJunk(
            pfchCandSelector, // select charged stuff from junk
            isolationConeFilter); // only take those in iso cone

        RegionalJunkConeAndIdFilter pfGammaJunk(
            pfgammaCandSelector, // select gammas from junk
            isolationConeFilter); // only take those in iso cone

        RegionalJunkConeAndIdFilter pfNeutralJunk(
            pfnhCandSelector, // select neutral stuff from junk
            isolationConeFilter); // select stuff in iso cone
 
	tau.addPiZeros(
            RecoTauConstructor::kIsolation,
            boost::make_filter_iterator(
                isolationConeFilterPiZero,
                piZeroCombo->remainder_begin(), piZeroCombo->remainder_end()),
            boost::make_filter_iterator(
                isolationConeFilterPiZero,
                piZeroCombo->remainder_end(), piZeroCombo->remainder_end()));
	
        tau.addPiZeros(
            RecoTauConstructor::kIsolation,
            boost::make_filter_iterator(
                isolationConeFilterPiZero,
                piZero_end, cleanPiZeros.end()),
            boost::make_filter_iterator(
                isolationConeFilterPiZero,
                cleanPiZeros.end(), cleanPiZeros.end()));

        // Filter the isolation candidates in a DR cone
	//
	// NOTE: isolation ChargedHadrons need to be added **after** signal and isolation PiZeros
	//       to avoid double-counting PFGammas as part of PiZero and merged with ChargedHadron
	//
        tau.addTauChargedHadrons(
            RecoTauConstructor::kIsolation,
            boost::make_filter_iterator(
                isolationConeFilterChargedHadron,
                trackCombo->remainder_begin(), trackCombo->remainder_end()),
            boost::make_filter_iterator(
                isolationConeFilterChargedHadron,
                trackCombo->remainder_end(), trackCombo->remainder_end()));

        // Add all the candidates that weren't included in the combinatoric
        // generation
        tau.addPFCands(
            RecoTauConstructor::kIsolation, RecoTauConstructor::kChargedHadron,
            boost::make_filter_iterator(
                pfCandFilter,
                pfch_end, pfchs.end()),
            boost::make_filter_iterator(
                pfCandFilter,
                pfchs.end(), pfchs.end()));
        // Add all charged candidates that are in the iso cone but weren't in the
        // original PFJet
        tau.addPFCands(
            RecoTauConstructor::kIsolation, RecoTauConstructor::kChargedHadron,
            boost::make_filter_iterator(
                pfChargedJunk, regionalJunk.begin(), regionalJunk.end()),
            boost::make_filter_iterator(
                pfChargedJunk, regionalJunk.end(), regionalJunk.end()));
	
        // Add all gammas that are in the iso cone but weren't in the
        // orginal PFJet
        tau.addPFCands(
            RecoTauConstructor::kIsolation, RecoTauConstructor::kGamma,
            boost::make_filter_iterator(
                pfGammaJunk, regionalJunk.begin(), regionalJunk.end()),
            boost::make_filter_iterator(
                pfGammaJunk, regionalJunk.end(), regionalJunk.end()));

        // Add all the neutral hadron candidates to the isolation collection
        tau.addPFCands(
            RecoTauConstructor::kIsolation, RecoTauConstructor::kNeutralHadron,
            boost::make_filter_iterator(
                pfCandFilter,
                pfnhs.begin(), pfnhs.end()),
            boost::make_filter_iterator(
                pfCandFilter,
                pfnhs.end(), pfnhs.end()));
        // Add all the neutral hadrons from the region collection that are in
        // the iso cone to the tau
        tau.addPFCands(
            RecoTauConstructor::kIsolation,  RecoTauConstructor::kNeutralHadron,
            boost::make_filter_iterator(
              pfNeutralJunk, regionalJunk.begin(), regionalJunk.end()),
            boost::make_filter_iterator(
              pfNeutralJunk, regionalJunk.end(), regionalJunk.end()));

        std::auto_ptr<reco::PFTau> tauPtr = tau.get(true);
	
	if ( primaryVertexRef.isNonnull() ) tauPtr->setVertex(primaryVertexRef->position());

        output.push_back(tauPtr);
      }
    }
  }

  return output.release();
}

}}  // end namespace reco::tau

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_EDM_PLUGIN(RecoTauBuilderPluginFactory,
                  reco::tau::RecoTauBuilderCombinatoricPlugin,
                  "RecoTauBuilderCombinatoricPlugin");
