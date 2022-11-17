// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
///
/// \Marta said: use the V0s as starting point for jet clustering, not their decay products (pions) and tro to get angular distance from leading jet momentum to V0 candidate
/// \author
//	Johanna Lömker
//	I have to embed this either in analysis-je-corelationV0Jet or in analysistutorial-... !!
//	
//	to run (step vzerotemplateexample):  o2-analysis-pid-tpc --configuration json://config.json | o2-analysis-multiplicity-table --configuration json://config.json | o2-analysistutorial-correlationV0Jet --configuration json://config.json | o2-analysis-track-propagation --configuration json://config.json | o2-analysis-timestamp --configuration json://config.json | o2-analysis-lf-lambdakzerobuilder --configuration json://config.json | o2-analysis-event-selection --configuration json://config.json -b
//
//	\since 2022
///
/// I assume that ZDCVZeroIteration.cxx is not the V0 particles, but rather detector part !
///
/// \brief jet analysis tasks (subscribing to jet finder task).
///        o2-analysis-jetfinder --aod-file AO2D.root | o2-analysistutorial-jet-analysis
/// \author Jochen Klein
/// \since
/// \brief in this tutorial the information contained in table TrackSelection is
///        used to retrieve tracks of a given type.
///        This requires tables TracksDCA and TrackSelection. Therefor use
///        o2-analysis-trackextension --aod-file AO2D.root | o2-analysis-trackselection | o2-analysistutorial-histogram-track-selection --select n
///        with n: 0 = no selection, 1 = globalTracks, 2 = globalTracksSDD
/// \author
/// \since
#include "Framework/runDataProcessing.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"

#include <cmath>
#include <string>
#include <TMath.h>
#include <TVector2.h>

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "EventFiltering/filterTables.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

//#include "PWGJE/Core/JetFinder.h"// this gives me troubles bc the fastjet/*.hh are not existing ?!

//#include "fastjet/PseudoJet.h"
//#include "fastjet/ClusterSequenceArea.h"

#include "PWGJE/DataModel/EMCALClusters.h"
#include "PWGJE/DataModel/Jet.h"

#include "Framework/HistogramRegistry.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct correlationvzerojets{
  //configurables for collision and track filter ! Crosscheck with configurables applied when using the jet finder
  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};
  Configurable<float> trackPtCut{"trackPtCut", 0.1, "minimum constituent pT"};
  Configurable<float> trackEtaCut{"trackEtaCut", 0.9, "constituent eta cut"};
  Configurable<int> nBins{"nBins", 100, "N bins in all histos"};

  Filter collisionFilter = nabs(aod::collision::posZ) < vertexZCut;
  Filter trackFilter = (nabs(aod::track::eta) < trackEtaCut) && (requireGlobalTrackInFilter()) && (aod::track::pt > trackPtCut);
  // using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra,  aod::TracksDCA, aod::TrackSelection>; //->change in process function track part to: soa::Filtered<TrackCandidates> const& tracks
  HistogramRegistry registry{
    "registry",
    {
      {"hVtxZ", "hVtxZ", {HistType::kTH1F, {{nBins, -15., 15.}}}},
      {"hMK0Short", "hMK0Short", {HistType::kTH1F, {{200,0.450f,0.550f}}}},//why these exact numbers - around mass peak i assume :P
      {"hMLambda", "hMLambda", {HistType::kTH1F, {{200, 0, 10}}}},
      {"hMAntiLambda", "hMAntiLambda", {HistType::kTH1F, {{200, 0, 10}}}},
      {"hMassTrueK0Short", "hMassTrueK0Short", {HistType::kTH1F, {{200,0.450f,0.550f}}}},
      {"hMassTrueLambda", "hMassTrueLambda", {HistType::kTH1F, {{200, 0, 10}}}},
      {"hMassTrueAntiLambda", "hMassTrueAntiLambda", {HistType::kTH1F, {{200, 0, 10}}}},

      {"jetPt", "jet p_{T};p_{T} (GeV/#it{c})", {HistType::kTH1F, {{100, 0, 100}}}},
      {"constPt", "constituent p_{T};p_{T} (GeV/#it{c})", {HistType::kTH1F, {{100, 0, 100}}}}

    }
  };
  
  void Jet(aod::Jet const& jet, aod::JetTrackConstituents const& constituents, aod::Tracks const& tracks)
  {
    registry.fill(HIST("jetPt"),jet.pt());
    for (const auto& c : constituents) {
      LOGF(debug, "jet %d: track id %d, track pt %g", jet.index(), c.trackId(), c.track().pt());
      registry.fill(HIST("constPt"), c.track().pt());
    }
  }
  PROCESS_SWITCH(correlationvzerojets, Jet, "process jets", true);
//now write the jet processing nd then add the V0 tables and blah on same iteration over collision
  void V0(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& collision, soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection>> const& tracks, aod::V0Datas const& V0s){
    if(!collision.sel7()){// sel8 is event selection based on TVX mainly for run3 data, but sel7=bool -> Event selection decision based on V0A & V0C -> https://aliceo2group.github.io/analysis-framework/docs/datamodel/helperTaskTables.html
      return;
     }
    // registry.fill(HIST("hVtxZ"),collision.posZ());

     // if (collision.hasJetChHighPt() >= bTriggerDecision) {
      registry.fill(HIST("hVtxZ"),collision.posZ()); // Inclusive Track Cross TPC Rows
      //}
      for(auto& v0 : V0s){
        registry.fill(HIST("hMK0Short"),v0.mK0Short());
      }
  }
  PROCESS_SWITCH(correlationvzerojets, V0, "process v0", true);
};
  
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{//run all via , separation
  // adaptAnalysisTask<V0>(cfgc),
  // adaptAnalysisTask<Jets>(cfgc),
   adaptAnalysisTask<correlationvzerojets>(cfgc)
  };
}
