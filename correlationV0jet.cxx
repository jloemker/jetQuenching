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

using MyTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTPCPi, aod::pidTPCPr>;
using MyTracksRun2 = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov, aod::TracksDCA, aod::pidTPCPi, aod::pidTPCPr>;//run 2 data stores Tracks in AO2D
using MyTracksRun3 = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA, aod::pidTPCPi, aod::pidTPCPr>;//run 3 data stores TracksIU (inner most update) in AO2D
//using LabeledV0s = soa::Join<aod::V0Datas, aod::McV0Labels>;

struct vzerotemplateexample {
  //Configurable for number of bins
  Configurable<int> nBins{"nBins", 100, "N bins in all histos"};

  //Selection criteria: 5 basic V0 selection criteria!
  Configurable<double> v0cospa{"v0cospa", 0.97, "V0 CosPA" };// double -> N.B. dcos(x)/dx = 0 at x=0
  Configurable<float> dcav0dau{"dcav0dau", 1.0, "DCA V0 Daughters"};
  Configurable<float> dcanegtopv{"dcanegtopv", -1, "DCA Neg To PV"};
  Configurable<float> dcapostopv{"dcapostopv", -1, "DCA Pos To PV"};//why would they also put -1 here ? 
  Configurable<float> v0radius{"v0radius", 0.5, "v0radius"};

  //Cannot filter on dynamic columns, so we cut on DCA to PV and DCA between daus only
  Filter preFilterV0 = nabs(aod::v0data::dcapostopv) > dcapostopv&&
                       nabs(aod::v0data::dcanegtopv) > dcanegtopv&&
                       aod::v0data::dcaV0daughters < dcav0dau;
  // histogram defined with HistogramRegistry
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
      {"constituentPt", "constituent p_{T};p_{T} (GeV/#it{c})", {HistType::kTH1F, {{100, 0, 100}}}}
    }
  };



  /* Jet and JetConstituents are tables which are filled by executable o2-analysis-jetfinder
  void getJets(aod::Jet const& jet, aod::JetTrackConstituents const& constituents, aod::Tracks const& tracks)
  {
    registry.fill(HIST("jetPt"),jet.pt());
    for (const auto& c : constituents) {
      LOGF(debug, "jet %d: track id %d, track pt %g", jet.index(), c.trackId(), c.track().pt());
      registry.fill(HIST("constituentPt"), c.track().pt());
    }
  }*/


  template <class TMyTracks, typename TV0>
  void processV0Candidate(TV0 const& v0, float const& pvx, float const& pvy, float const& pvz)
  {//function to process a vzero candidate freely, actually with right track type !
   auto posTrackCast = v0.template posTrack_as<TMyTracks>();
   auto negTrackCast = v0.template negTrack_as<TMyTracks>();

   float nsigma_pos_proton = TMath::Abs(posTrackCast.tpcNSigmaPr());
   float nsigma_neg_proton = TMath::Abs(posTrackCast.tpcNSigmaPr());//shouldn't it be negTrackCast ??
   float nsigma_pos_pion = TMath::Abs(negTrackCast.tpcNSigmaPi());//shouldn't it be posTrackCast ??
   float nsigma_neg_pion = TMath::Abs(negTrackCast.tpcNSigmaPi());
   
   if( v0.v0radius() > v0radius && v0.v0cosPA(pvx, pvy, pvz) > v0cospa ){
     if( nsigma_pos_pion < 4 && nsigma_neg_pion < 4){// we want to use V0 for jet clustering, but we need decay product for V0 reconsturction, right ?
       registry.fill(HIST("hMK0Short"), v0.mK0Short());
     }
     if( nsigma_pos_proton < 4 && nsigma_neg_proton < 4){
       registry.fill(HIST("hMLambda"), v0.mLambda());
     }
     if( nsigma_pos_pion < 4 && nsigma_neg_proton < 4){
       registry.fill(HIST("hMAntiLambda"), v0.mAntiLambda());
     }   
     // for this i need an additional .h file and another subscription in my pipeline !
     /*An adding the final STEP *5* - lambdakzerobuilder is able to do MC association for me #wow !
     if( v0.has_mcParticle()){//soem association was made !
       auto v0mcparticle = v0.mcParticle();
       //Check particle PDG code to see if this is the one you want
       if( v0mcparticle.pdgCode() == 310 ) registry.fill(HIST("hMassTrueK0Short"), v0.mK0Short());
       if( v0mcparticle.pdgCode() == 3122 ) registry.fill(HIST("hMassTrueLambda"), v0.mLambda());// wtf why the same odg code for lamda and anti lambda ??
       if( v0mcparticle.pdgCode() == 3122 ) registry.fill(HIST("hMassTrueAntiLambda"), v0.mAntiLambda());
     }*/
   }
  }
  //*
    //Basic event selection (all helper tasks are now included!)
    //https://aliceo2group.github.io/analysis-framework/docs/datamodel/helperTaskTables.html
    //o2::aod::evsel::Sel7              sel7    bool    Event selection decision based on V0A & V0C
    //o2::aod::evsel::Sel8              sel8    bool    Event selection decision based on TVX)
    // (...)
    //check getter here: https://aliceo2group.github.io/analysis-framework/docs/datamodel/ao2dTables.html
  //*
  //Process function for old data (run2)//change from LabeldV0s to Vo0Datas and removed Mcparticle
  void processRun2(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, soa::Filtered<aod::V0Datas> const& V0s, MyTracksRun2 const& tracks)
  {
    if( !collision.sel7()){//sel7=bool -> Event selection decision based on V0A & V0C + all helper tasks required for this !
      return;
    }
    registry.get<TH1>(HIST("hVtxZ"))->Fill(collision.posZ());
    for( auto& v0 : V0s){
      processV0Candidate<MyTracksRun2>(v0, collision.posX(), collision.posY(), collision.posZ() );
    }
  }
  PROCESS_SWITCH(vzerotemplateexample, processRun2, "Process Run 2 data", true);
  //Define process function for run3 data
  void processRun3(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, soa::Filtered<aod::V0Datas> const& V0s, MyTracksRun3 const& tracks)
  {
    if( !collision.sel8()){//sel8=bool -> Event selection decision based on TVX) ! Event selection is based on a different detector part for run2 / run3 !
      return;
    }
    registry.get<TH1>(HIST("hVtxZ"))->Fill(collision.posZ());
    for( auto& v0 : V0s){
      processV0Candidate<MyTracksRun3>(v0, collision.posX(), collision.posY(), collision.posZ());
    }
  }
  PROCESS_SWITCH(vzerotemplateexample, processRun3, "Process Run 3 data", false);

};

struct Jets{
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
      {"hMassTrueAntiLambda", "hMassTrueAntiLambda", {HistType::kTH1F, {{200, 0, 10}}}},

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

  void V0(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, aod::V0Datas const& V0s){
     //now getting V0
    if(!collision.sel8()){//event selection based on TVX -> https://aliceo2group.github.io/analysis-framework/docs/datamodel/helperTaskTables.html
      return;
     }
     registry.fill(HIST("hVertexZ"),collision.posZ());
     for(auto& v0 : V0s){
       registry.fill(HIST("hMassK0Short"), v0.mK0Short());
     }
   }
};



WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{//run all via , separation
   // adaptAnalysisTask<vzerotemplateexample>(cfgc)
   adaptAnalysisTask<Jets>(cfgc)
  };
}
