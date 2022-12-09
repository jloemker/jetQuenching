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
/// \Marta said: use the V0s as starting point for jet clustering, not their decay products (pions) and try to get angular distance from leading jet momentum to V0 candidate
//		-> should ask how the track tables from the jet process and v0 process are now related bc I assume the jet finder might mere things ?
//		! the track pt sprectrum from all tracks is enourmous and the constituent track and full jet pt are small
/// \author
//	Johanna Lömker
//	to run (step vzerotemplateexample): check the run.sh and config.json
//
//	o2-analysis-pid-tpc --configuration json://config.json | o2-analysis-multiplicity-table --configuration json://config.json | o2-analysistutorial-correlationV0Jet --configuration json://config.json | o2-analysis-track-propagation --configuration json://config.json | o2-analysis-timestamp --configuration json://config.json | o2-analysis-lf-lambdakzerobuilder --configuration json://config.json | o2-analysis-event-selection --configuration json://config.json -b
//	\since 2022
///        with n: 0 = no selection, 1 = globalTracks, 2 = globalTracksSDD

#include <cmath>
#include <string>
#include <TMath.h>
#include <TVector2.h>

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"

#include "EventFiltering/filterTables.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"

#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGJE/Core/JetFinder.h"// this gives me troubles bc the fastjet/*.hh are not existing ?!
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/EMCALClusters.h"
//#include "fastjet/PseudoJet.h"
//#include "fastjet/ClusterSequenceArea.h"



using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using MyTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTPCPi, aod::pidTPCPr>;//This was used for the V0 invariant mass tutorial
using JetTrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>;//This was used in the ChJetTriggerQA
using CombinedTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::pidTPCPi, aod::pidTPCPr, aod::TrackSelection>;//This is what I use for my analysis (for now)

struct correlationvzerojets{
  //configurables for collision and track filter ! Crosscheck with configurables applied when using the jet finder
  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};
  Configurable<float> trackPtCut{"trackPtCut", 0.1, "minimum constituent pT"};
  Configurable<float> trackEtaCut{"trackEtaCut", 0.9, "constituent eta cut"};//corresponds to TPCVolume in triggerQA
 
  //configurables for jets
  Configurable<float> JetR{"JetR", 0.4, "jet resolution parameter"}; // jet cone radius
  Configurable<float> JetPtMin{ "JetPtMin", 0.1,"minimum jet pT constituent cut"}; // minimum jet constituent pT
  Configurable<int> bTriggerDecision{"bTriggerDecision", 0,"Charged Jet Trigger Decision Selection"}; // 0=MB Event, 1=Event selected by EPN

  //configurables for the V0 identification of lambda and K0short through  selection criteria on daughters
  Configurable<double> v0cospa{"v0cospa", 0.97, "V0 CosPA" };// double -> N.B. dcos(x)/dx = 0 at x=0  what is this ??
  Configurable<float> dcav0dau{"dcav0dau", 1.0, "DCA V0 Daughters"};
  Configurable<float> dcanegtopv{"dcanegtopv", -1, "DCA Neg To PV"};
  Configurable<float> dcapostopv{"dcapostopv", -1, "DCA Pos To PV"};//why would they also put -1: bc this is the initiaization, updated in process du dummy 
  Configurable<float> v0radius{"v0radius", 0.5, "v0radius"};//what is this ? decay cone ? no idea

  Configurable<int> nBins{"nBins", 200, "N bins in histos"};
  Configurable<int> nBinsPt{"nBinsPt", 200, "N bins in pT histos"};
  Configurable<int> nBinsEta{"nBinsEta", 200, "N bins in Eta histos"};
  Configurable<int> nBinsPhi{"nBinsPhi", 200, "N bins in Phi histos"};

  float fiducialVolume;                        // 0.9 - jetR
  std::vector<fastjet::PseudoJet> jetConstituents;
  std::vector<fastjet::PseudoJet> jetReclustered;
  JetFinder jetReclusterer;

  Filter collisionFilter = nabs(aod::collision::posZ) < vertexZCut;
  Filter trackFilter = (nabs(aod::track::eta) < trackEtaCut) && (requireGlobalTrackInFilter()) && (aod::track::pt > trackPtCut);
  // TrackSelection globalTracks;
  void init(o2::framework::InitContext&)
  {
    jetReclusterer.isReclustering = true;
    jetReclusterer.algorithm = fastjet::JetAlgorithm::antikt_algorithm;
    jetReclusterer.jetR = JetR;
    jetReclusterer.jetEtaMin = -trackEtaCut;
    jetReclusterer.jetEtaMax = trackEtaCut;
    jetReclusterer.jetPtMax = 1e9;

    fiducialVolume = trackEtaCut - JetR;
  }


  // using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra,  aod::TracksDCA, aod::TrackSelection>; //->change in process function track part to: soa::Filtered<TrackCandidates> const& tracks
  HistogramRegistry registry{
    "registry",
    {
      {"hCollVtxZ", "hCollisionVtxZ", {HistType::kTH1F, {{nBins, -15., 15.}}}},

      {"hMK0Short", "hMK0Short; M (GeV/#it{c})", {HistType::kTH1F, {{200,0.450f,0.550f}}}},//why these exact numbers - around mass peak i assume :P
      {"hPtK0Short", "if K0Short : v0.pt()  p_{T}; p_{T} (GeV/#it{c})", {HistType::kTH1F, {{nBinsPt, 0, 100}}}},
      {"hEtaK0Short", "if K0Short : v0.eta() #eta; #eta", {HistType::kTH1F, {{nBinsEta, -0.9, 0.9}}}},
      {"hPhiK0Short", "if K0Short : v0.phi() #phi; #phi", {HistType::kTH1F, {{nBinsPhi, 0, 6.3}}}},

      {"hMLambda", "hMLambda; M (GeV/#it{c})", {HistType::kTH1F, {{200, 0, 10}}}},
      {"hPtLambda", "if Lambda : v0.pt()  p_{T}; p_{T} (GeV/#it{c})", {HistType::kTH1F, {{nBinsPt, 0, 100}}}},
      {"hEtaLambda", "if Lambda : v0.eta() #eta; #eta", {HistType::kTH1F, {{nBinsEta, -0.9, 0.9}}}},
      {"hPhiLambda", "if Lambda : v0.phi() #phi; #phi", {HistType::kTH1F, {{nBinsPhi, 0, 6.3}}}},

      {"hMAntiLambda", "hMAntiLambda; M (GeV/#it{c})", {HistType::kTH1F, {{200, 0, 10}}}},
      {"hPtAntiLambda", "if AntiLambda : v0.pt()  p_{T}; p_{T} (GeV/#it{c})", {HistType::kTH1F, {{nBinsPt, 0, 100}}}},
      {"hEtaAntiLambda", "if AntiLambda : v0.eta() #eta; #eta", {HistType::kTH1F, {{nBinsEta, -0.9, 0.9}}}},
      {"hPhiAntiLambda", "if AntiLambda : v0.phi() #phi; #phi", {HistType::kTH1F, {{nBinsPhi, 0, 6.3}}}},

      //Analysis plots - adjust binning and normalize that stuff ... something with integral ...
      {"LambdaOverKaonPt", "(Lambda + AntiLambda)/2K p_{T}; p_{T} (GeV/#it{c})", {HistType::kTH1F, {{nBinsPt, 0, 60}}}},

      //control plots	
      {"hPtTrackV0inRadius", "V0 from V0Datas in V0 radius  p_{T}; p_{T} (GeV/#it{c})", {HistType::kTH1F, {{nBins, 0, 100}}}},
      {"hEtaTrackV0inRadius", "V0 from V0Datas in V0 radius #eta; #eta", {HistType::kTH1F, {{nBinsEta, -0.9, 0.9}}}},
      {"hPhiTrackV0inRadius", "V0 from V0Datas in V0 radius  #phi; #phi", {HistType::kTH1F, {{nBinsPhi, 0, 6.3}}}},

      {"hPtV0", "V0 in V0Datas p_{T}; p_{T} (GeV/#it{c})", {HistType::kTH1F, {{nBinsPt, 0, 100}}}},
      {"hEtaV0", "V0 in V0Datas #eta; #eta", {HistType::kTH1F, {{nBinsEta, -0.9, 0.9}}}},
      {"hPhiV0", "V0 in V0Datas  #phi; #phi", {HistType::kTH1F, {{nBinsPhi, 0, 6.3}}}},

      {"hTrackPt", "MyTracks in collision track p_{T}; p_{T} (GeV/#it{c})", {HistType::kTH1F, {{nBinsPt, 0, 100}}}},
      {"hTrackEta", "MyTracks in collision: track #eta; #eta", {HistType::kTH1F, {{nBinsEta, -0.9, 0.9}}}},
      {"hTrackPhi", "MyTracks in collision: track #phi; #phi", {HistType::kTH1F, {{nBinsPhi, 0, 6.3}}}},

      //jets as next 
      {"jetVtx", "jet vtxZ; vtxZ [cm] ", {HistType::kTH1F, {{nBins, -15, 15}}}},
      {"JetTrackPt", "inclusive track p_{T}; p_{T} (GeV/#it{c})", {HistType::kTH1F, {{nBinsPt, 0, 100}}}},
      {"JetTrackEta", "inclusive track #eta; #eta", {HistType::kTH1F, {{nBinsEta, -0.9, 0.9}}}},
      {"JetTrackPhi", "inclusive track #phi; #phi", {HistType::kTH1F, {{nBinsPhi, 0, 6.3}}}},

      {"JetLeadTrackPt", "leading track p_{T}; p_{T} (GeV/#it{c})", {HistType::kTH1F, {{nBinsPt, 0, 100}}}},
      {"JetLeadTrackEta", "leading track #eta; #eta", {HistType::kTH1F, {{nBinsEta, -0.9, 0.9}}}},
      {"JetLeadTrackPhi", "leading track #phi; #phi", {HistType::kTH1F, {{nBinsPhi, 0, 6.3}}}},

      {"JetLeadJetPt", "leading jet p_{T}; p_{T} (GeV/#it{c})", {HistType::kTH1F, {{nBinsPt, 0, 100}}}},
      {"JetLeadJetEta", "leading jet #eta; #eta", {HistType::kTH1F, {{nBinsEta, -0.9, 0.9}}}},
      {"JetLeadJetPhi", "leading jet #phi; #phi", {HistType::kTH1F, {{nBinsPhi, 0, 6.3}}}},
      
      {"jetV0Pt", "V0 in jet code p_{T};p_{T} (GeV/#it{c})", {HistType::kTH1F, {{nBinsPt, 0, 100}}}},
      {"jetV0Eta", "V0 in jet code #eta; #eta", {HistType::kTH1F, {{nBinsEta,-0.9, 0.9}}}},
      {"jetV0Phi", "V0 in jet code #phi; #phi ", {HistType::kTH1F, {{nBinsPhi, 0, 6.3}}}},
      {"jetWithV0Pt", "V0 in jet with collId requ. p_{T};p_{T} (GeV/#it{c})", {HistType::kTH1F, {{nBinsPt, 0, 100}}}},
      {"jetWithV0Eta", "V0 in jet with collId requ. #eta; #eta", {HistType::kTH1F, {{nBinsEta,-0.9, 0.9}}}},
      {"jetWithV0Phi", "V0 in jet with collId requ.  #phi; #phi ", {HistType::kTH1F, {{nBinsPhi, 0, 6.3}}}}, 

      {"AngularDistance", "Angular distance(leading J - V0); #Delta R", {HistType::kTH1F, {{nBins, 0, 10}}}}
     // {"jetPt", "jet p_{T};p_{T} (GeV/#it{c})", {HistType::kTH1F, {{nBinsPt, 0, 100}}}},
     // {"jetconstTrackPt", "constituent track p_{T};p_{T} (GeV/#it{c})", {HistType::kTH1F, {{nBinsPt, 0, 100}}}},
     // {"jetTrackPt",  "jetTrack in tracks p_{T};p_{T} (GeV/#it{c})", {HistType::kTH1F, {{nBinsPt, 0, 100}}}},
     // {"jetMyTrackPt",  "jetMyTrack in mytracks p_{T};p_{T} (GeV/#it{c})", {HistType::kTH1F, {{nBinsPt, 0, 100}}}}
    }
  };

  void Jet(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::JetFilters>>::iterator const& collision, soa::Filtered<CombinedTracks> const& tracks, aod::V0Datas const& V0s)//figure out why only aod::Tracks or MyTracks and/or if they are th same
  { 
   if (collision.hasJetChHighPt() >= bTriggerDecision) {
      jetConstituents.clear();
      jetReclustered.clear();

      float leadingJetPt = -1.0;
      float leadingJetEta = -2.0;
      float leadingJetPhi = -1.0;
      float leadingTrackPt = -1.0;
      float leadingTrackEta = -2.0;
      float leadingTrackPhi = -1.0;
      
      float angularDistance = -99;

      registry.fill(HIST("jetVtx"),collision.posZ());
      for (auto& trk : tracks) { //loop over filtered tracks in full TPC volume having pT > 100 MeV
        if (!trk.isQualityTrack()) {
          continue; // skip bad quality tracks
        }

       registry.fill(HIST("JetTrackPt"), trk.pt());
       registry.fill(HIST("JetTrackEta"),trk.eta());
       registry.fill(HIST("JetTrackPhi"),trk.phi());
       
       fillConstituents(trk,jetConstituents); // ./PWGJE/Core/JetFinder.h
                            // recombination scheme is assumed
                            // to be Escheme with pion mass
        if (trk.pt() > leadingTrackPt) { // Find leading track pT in full TPC volume
          leadingTrackPt = trk.pt();
          leadingTrackEta = trk.eta();
          leadingTrackPhi = trk.phi();
        }
      }//end tracks

      for(auto& v0 : V0s){
        registry.fill(HIST("jetV0Pt"), v0.pt());
	registry.fill(HIST("jetV0Eta"), v0.eta());
	registry.fill(HIST("jetV0Phi"), v0.phi());
        
       // fillConstituents(v0,jetConstituents); // if this works, then we need two - 1 from tracks one from v0 or even from v0 daughters ?
      }

      if (leadingTrackPt > -1.) {
        registry.fill(HIST("JetLeadTrackPt"), leadingTrackPt);
	registry.fill(HIST("JetLeadTrackPhi"), leadingTrackPhi);
        registry.fill(HIST("JetLeadTrackEta"), leadingTrackEta);
      }

      // Reconstruct jet from tracks
      fastjet::ClusterSequenceArea clusterSeq(jetReclusterer.findJets(jetConstituents, jetReclustered));

      // Find leading jet pT in full TPC volume
      for (auto& jet : jetReclustered) {
        if (fabs(jet.eta()) < 2 * trackEtaCut){ // AK:  cfgTPCVolume) {

          if (jet.perp() > leadingJetPt) {
            leadingJetPt = jet.perp();
            leadingJetEta = jet.eta();
            leadingJetPhi = jet.phi();
          }
        }
      }

      if (leadingJetPt > -1.) {
        registry.fill(HIST("JetLeadJetPt"), leadingJetPt);
        registry.fill(HIST("JetLeadJetPhi"), leadingJetPhi);
        registry.fill(HIST("JetLeadJetEta"), leadingJetEta);
        for(auto& v0 : V0s){
	 if(v0.collisionId() == collision.globalIndex()){
           registry.fill(HIST("jetWithV0Pt"), v0.pt());
           registry.fill(HIST("jetWithV0Eta"), v0.eta());
           registry.fill(HIST("jetWithV0Phi"), v0.phi());

	   angularDistance = sqrt(pow(leadingJetEta-v0.eta(),2) + pow(leadingJetPhi-v0.phi(),2) );
	   registry.fill(HIST("AngularDistance"), angularDistance);
  	  }
        }//end V0s    
       }//end if leading jet
      
	
    }//end collision

  }
  PROCESS_SWITCH(correlationvzerojets, Jet, "process jets", true);

// void V0(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& collision, soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection>> const& tracks, aod::V0Datas const& V0s){ // for process funtion over all tracks within event selection - so this should be used and then subproccess over MyJetTracks and MyV0Tracks 
  void V0(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& collision,  soa::Filtered<CombinedTracks> const& tracks, aod::V0Datas const& V0s){//no more my tracks
    if(!collision.sel8()){// sel8 is event selection based on TVX mainly for run3 data, but sel7=bool -> Event selection decision based on V0A & V0C -> https://aliceo2group.github.io/analysis-framework/docs/datamodel/helperTaskTables.html
      return;
    }

    registry.fill(HIST("hCollVtxZ"),collision.posZ()); // Inclusive Tracks from sel7 selections and the aod::pidTPCPi, aod::pidTPCPr
	
    for(auto& v0 : V0s){
      // particle identification by using the tpcNSigma track variable
      float nsigma_pos_proton = TMath::Abs(v0.posTrack_as<CombinedTracks>().tpcNSigmaPr());//what does tpcNSigmaPr / Pi exactly ??
      float nsigma_neg_proton = TMath::Abs(v0.negTrack_as<CombinedTracks>().tpcNSigmaPr());
      float nsigma_pos_pion = TMath::Abs(v0.posTrack_as<CombinedTracks>().tpcNSigmaPi());
      float nsigma_neg_pion = TMath::Abs(v0.negTrack_as<CombinedTracks>().tpcNSigmaPi());
      
      if(v0.v0radius() > v0radius && v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) > v0cospa ){//i dont understand this criteria 
        if( nsigma_pos_pion < 4 && nsigma_neg_pion < 4 ){
          registry.fill(HIST("hMK0Short"), v0.mK0Short());
	  registry.fill(HIST("hPtK0Short"), v0.pt());
	  registry.fill(HIST("hEtaK0Short"), v0.eta());
          registry.fill(HIST("hPhiK0Short"), v0.phi());
       	}
	if( nsigma_pos_proton < 4 && nsigma_neg_proton < 4){
          registry.fill(HIST("hMLambda"), v0.mLambda());
	  registry.fill(HIST("hPtLambda"), v0.pt());
          registry.fill(HIST("hEtaLambda"), v0.eta());
          registry.fill(HIST("hPhiLambda"), v0.phi());
	}
	if( nsigma_pos_pion < 4 && nsigma_neg_proton < 4){
          registry.fill(HIST("hMAntiLambda"), v0.mAntiLambda());	
	  registry.fill(HIST("hPtAntiLambda"), v0.pt());
          registry.fill(HIST("hEtaAntiLambda"), v0.eta());
          registry.fill(HIST("hPhiAntiLambda"), v0.phi());
         }
        // Get V0 within radius and v0cosPA
        registry.fill(HIST("hPtTrackV0inRadius"), v0.pt());
        registry.fill(HIST("hEtaTrackV0inRadius"), v0.eta());
        registry.fill(HIST("hPhiTrackV0inRadius"), v0.phi());
	}//to have control plots for the V0 in the jet iteration above
	registry.fill(HIST("hPtV0"), v0.pt());
	registry.fill(HIST("hEtaV0"), v0.eta());
        registry.fill(HIST("hPhiV0"), v0.phi());
      }	
 
      for(auto& track : tracks){//Check kinematics of MyTracks within the sel7 process selection !
        registry.fill(HIST("hTrackPt"), track.pt());
        registry.fill(HIST("hTrackEta"), track.eta());
        registry.fill(HIST("hTrackPhi"), track.phi());
      }

      // Baryon over Meson ratio
      for(int i = 0; i<nBinsPt; i++){
        double lamb = registry.get<TH1>(HIST("hPtLambda"))->GetBinContent(i);
        double antl = registry.get<TH1>(HIST("hPtAntiLambda"))->GetBinContent(i);
        double kaon = registry.get<TH1>(HIST("hPtK0Short"))->GetBinContent(i);

        if(kaon != 0){
          double ratio = (lamb+antl)/(2*kaon);
          registry.fill(HIST("LambdaOverKaonPt"), registry.get<TH1>(HIST("hPtLambda"))->GetBinCenter(i), ratio);
        }
      }

   }
  PROCESS_SWITCH(correlationvzerojets, V0, "process v0 and their track QA", true);

};
  
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{//run all via , separation
  // adaptAnalysisTask<JetQATask>(cfgc)
     adaptAnalysisTask<correlationvzerojets>(cfgc)
  //  adaptAnalysisTask<jets>
  };
}
