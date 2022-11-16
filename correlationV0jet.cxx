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

#include "Common/Core/TrackSelection.h"
#include "PWGJE/DataModel/Jet.h"

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct Jets{
  Configurable<int> nBins{"nBins", 100, "N bins in all histos"};
  HistogramRegistry registry{
    "registry",
    {
      {"jetPt", "jet p_{T};p_{T} (GeV/#it{c})", {HistType::kTH1F, {{100, 0, 100}}}},
      {"constPt", "constituent p_{T};p_{T} (GeV/#it{c})", {HistType::kTH1F, {{100, 0, 100}}}}
    }
  };

  // Jet and JetConstituents are tables which are filled by executable o2-analysis-jetfinder
  void process(aod::Jet const& jet, aod::JetTrackConstituents const& constituents, aod::Tracks const& tracks)
  {
    registry.fill(HIST("jetPt"), jet.pt());
    for (const auto& c : constituents) {
      LOGF(debug, "jet %d: track id %d, track pt %g", jet.index(), c.trackId(), c.track().pt());
      registry.fill(HIST("constPt"), c.track().pt());
    }
    for(auto& track : tracks ){
      if(track.has_collision()){
       // V0(track.collision(), V0s);
        

      }
    }
  }
};
//or with process switch ?
struct V0{
  Configurable<int> nBins{"nBins", 100, "N bins in all histos"};
  HistogramRegistry registry{
    "registry",
    {
      {"hVtxZ", "hVtxZ", {HistType::kTH1F, {{nBins, -15., 15.}}}},
      {"hMK0Short", "hMK0Short", {HistType::kTH1F, {{200,0.450f,0.550f}}}},//why these exact numbers - around mass peak i assume :P
      {"hMLambda", "hMLambda", {HistType::kTH1F, {{200, 0, 10}}}},
      {"hMAntiLambda", "hMAntiLambda", {HistType::kTH1F, {{200, 0, 10}}}},
      {"hMassTrueK0Short", "hMassTrueK0Short", {HistType::kTH1F, {{200,0.450f,0.550f}}}},
      {"hMassTrueLambda", "hMassTrueLambda", {HistType::kTH1F, {{200, 0, 10}}}},
      {"hMassTrueAntiLambda", "hMassTrueAntiLambda", {HistType::kTH1F, {{200, 0, 10}}}}
    }
  };
      
  void process(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, aod::V0Datas const& V0s){
     //now getting V0
    if(!collision.sel8()){//event selection based on TVX -> https://aliceo2group.github.io/analysis-framework/docs/datamodel/helperTaskTables.html
      return;
     }
     registry.fill(HIST("hVtxZ"),collision.posZ());
     for(auto& v0 : V0s){
       registry.fill(HIST("hMassK0Short"), v0.mK0Short());
     }
   }
};



WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{//run all via , separation
   adaptAnalysisTask<V0>(cfgc),
   adaptAnalysisTask<Jets>(cfgc)
  };
}
