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
//		this means I would need to produce my own jet table in the main process - buuuut : I cannot include the jet.h bc I get some fastjet/blah.hh does not exist problems
//		-> should ask the experts, but first I pepare some more kinematics on the V0 and then on the decay products - pion and proton - tutorial for it
//		-> should ask how the track tables from the jet process and v0 process are now related bc I assume the jet finder might mere things ?
//		otoh: the track pt sprectrum from all tracks is enourmous and the constituent track and full jet pt are small
/// \author
//	Johanna Lömker
//	to run (step vzerotemplateexample):  o2-analysis-pid-tpc --configuration json://config.json | o2-analysis-multiplicity-table --configuration json://config.json | o2-analysistutorial-correlationV0Jet --configuration json://config.json | o2-analysis-track-propagation --configuration json://config.json | o2-analysis-timestamp --configuration json://config.json | o2-analysis-lf-lambdakzerobuilder --configuration json://config.json | o2-analysis-event-selection --configuration json://config.json -b
//	\since 2022
///        with n: 0 = no selection, 1 = globalTracks, 2 = globalTracksSDD
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

using MyTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTPCPi, aod::pidTPCPr>;
//maybe something like MyJets for the different set of joined tables ?
//- a global process function over all traks in this thing and then we give a set of selected tracks into the jet and V0 iterations ?
//
//
struct correlationvzerojets{
  //configurables for collision and track filter ! Crosscheck with configurables applied when using the jet finder
  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};
  Configurable<float> trackPtCut{"trackPtCut", 0.1, "minimum constituent pT"};
  Configurable<float> trackEtaCut{"trackEtaCut", 0.9, "constituent eta cut"};

  //configurables for the V0 identification of lambda and K0short through  selection criteria on daughters
  Configurable<double> v0cospa{"v0cospa", 0.97, "V0 CosPA" };// double -> N.B. dcos(x)/dx = 0 at x=0
  Configurable<float> dcav0dau{"dcav0dau", 1.0, "DCA V0 Daughters"};
  Configurable<float> dcanegtopv{"dcanegtopv", -1, "DCA Neg To PV"};
  Configurable<float> dcapostopv{"dcapostopv", -1, "DCA Pos To PV"};//why would they also put -1: bc this is the initiaization, updated in process du dummy 
  Configurable<float> v0radius{"v0radius", 0.5, "v0radius"};
  Configurable<int> nBins{"nBins", 100, "N bins in all pT histos"};

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
      {"hTrackPt", "track p_{T}; p_{T} (GeV/#it{c})", {HistType::kTH1F, {{nBins, 0, 100}}}},
      {"hPtLambda", "v0.ptlambda p_{T}; p_{T} (GeV/#it{c})", {HistType::kTH1F, {{nBins, 0, 100}}}},
      {"hPtTrackV0", "v0.pt p_{T}; p_{T} (GeV/#it{c})", {HistType::kTH1F, {{nBins, 0, 100}}}},

      {"jetPt", "jet p_{T};p_{T} (GeV/#it{c})", {HistType::kTH1F, {{100, 0, 100}}}},
      {"constTrackPt", "constituent track p_{T};p_{T} (GeV/#it{c})", {HistType::kTH1F, {{100, 0, 100}}}}

    }
  };
  
  void Jet(aod::Jet const& jet, aod::JetTrackConstituents const& constituents, aod::Tracks const& tracks)
  {//after the fastjet mystery is solved I can produce my own jet table as in the ChJetTriggerQA.cxx and then the process switch should be used fro MC/Run2/Run3 data 
    registry.fill(HIST("jetPt"),jet.pt());
    for (const auto& c : constituents) {
      LOGF(debug, "jet %d: track id %d, track pt %g", jet.index(), c.trackId(), c.track().pt());
      registry.fill(HIST("constTrackPt"), c.track().pt());
    }
  }
  PROCESS_SWITCH(correlationvzerojets, Jet, "process jets", true);
//now write the jet processing nd then add the V0 tables and blah on same iteration over collision
// void V0(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& collision, soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection>> const& tracks, aod::V0Datas const& V0s){ // for process funtion over all tracks within event selection - so this should be used and then subproccess over MyJetTracks and MyV0Tracks 
  void V0(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& collision, MyTracks const& tracks, aod::V0Datas const& V0s){
    if(!collision.sel7()){// sel8 is event selection based on TVX mainly for run3 data, but sel7=bool -> Event selection decision based on V0A & V0C -> https://aliceo2group.github.io/analysis-framework/docs/datamodel/helperTaskTables.html
      return;
     }

     registry.fill(HIST("hVtxZ"),collision.posZ()); // Inclusive Tracks from sel7 selections and the aod::pidTPCPi, aod::pidTPCPr
	
      for(auto& v0 : V0s){
	// particle identification by using the tpcNSigma track variable
        float nsigma_pos_proton = TMath::Abs(v0.posTrack_as<MyTracks>().tpcNSigmaPr());
        float nsigma_neg_proton = TMath::Abs(v0.negTrack_as<MyTracks>().tpcNSigmaPr());
        float nsigma_pos_pion = TMath::Abs(v0.posTrack_as<MyTracks>().tpcNSigmaPi());
        float nsigma_neg_pion = TMath::Abs(v0.negTrack_as<MyTracks>().tpcNSigmaPi());

        if(v0.v0radius() > v0radius && v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) > v0cospa ){
          if( nsigma_pos_pion < 4 && nsigma_neg_pion < 4 ){
	    registry.fill(HIST("hMK0Short"), v0.mK0Short());
	  }
	  if( nsigma_pos_proton < 4 && nsigma_neg_proton < 4){
            registry.fill(HIST("hMLambda"), v0.mLambda());
	  }
	  if( nsigma_pos_pion < 4 && nsigma_neg_proton < 4){
            registry.fill(HIST("hMAntiLambda"), v0.mAntiLambda());	
	  }
        // registry.fill(HIST("hPtLambda"), v0.Lambda());
        registry.fill(HIST("hPtTrackV0"), v0.pt());
	
      }

      }	
      for(auto& track : tracks){
        registry.fill(HIST("hTrackPt"), track.pt());
      }
  }
  PROCESS_SWITCH(correlationvzerojets, V0, "process v0", true);
};
  
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{//run all via , separation
   adaptAnalysisTask<correlationvzerojets>(cfgc)
  };
}
