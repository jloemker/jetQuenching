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
//  \author
//	Johanna Lömker
//	to run (step vzerotemplateexample): check the run.sh and config.json
//
//	\since 2022
/// with n: 0 = no selection, 1 = globalTracks, 2 = globalTracksSDD

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
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/EMCALClusters.h"



using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using MyTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTPCPi, aod::pidTPCPr>;//This was used for the V0 invariant mass tutorial
using JetTrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>;//This was used in the ChJetTriggerQA
using CombinedTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::pidTPCPi, aod::pidTPCPr, aod::TrackSelection>;//This is what I use for my analysis (for now)
//there are also pidTOFFullPr ... but this did not work, bc of the lambda k zero builder
struct correlationvzerojets{
  //configurables for collision and track filter ! Crosscheck with configurables applied when using the jet finder
  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};
  Configurable<float> trackPtCutMin{"trackPtCut", 0.1, "minimum constituent pT"};
  Configurable<float> trackEtaCut{"trackEtaCut", 0.8, "constituent eta cut"};//corresponds to TPCVolume in triggerQA
 
  //configurables for jets
  Configurable<float> JetR{"JetR", 0.4, "jet resolution parameter"}; // jet cone radius
  Configurable<float> JetPtMin{ "JetPtMin", 0.1,"minimum jet pT constituent cut"}; // minimum jet constituent pT
  Configurable<int> bTriggerDecision{"bTriggerDecision", 0,"Charged Jet Trigger Decision Selection"}; // 0=MB Event, 1=Event selected by EPN

  //configurables for the V0 identification of lambda and K0short through  selection criteria on daughters
  //Configurable<double> v0cospa{"v0cospa", 0.9, "V0 CosPA in prefilter" }; 
  Configurable<float> V0daugPtMin{"V0daugPtMin", 1, "Minimum pt of daughter tracks"};
  Configurable<float> V0ptMin{"V0ptMin", 0.6, "MinimumV0pt"};
  Configurable<float> V0ptMax{"V0pTMax", 20, "MaximumV0pt"}; 
  Configurable<double> v0cospaLamb{"v0cospaLamb", 0.995, "V0 CosPA Lambda" };// double -> N.B. dcos(x)/dx = 0 at x=0 
  Configurable<double> v0cospaK0s{"v0cospaK0s", 0.97, "V0 CosPA Kaons" };
  Configurable<float> dcav0dau{"dcav0dau", 1.0, "DCA V0 Daughters"};//dca between daughters in sigma units
  Configurable<float> dcanegtopv{"dcanegtopv", 0.06, "DCA Neg To PV"};
  Configurable<float> dcapostopv{"dcapostopv", 0.06, "DCA Pos To PV"};//why would they also put -1: bc this is the initiaization, updated in process du dummy 
  Configurable<float> MinV0radius{"MinV0radius", 5, "minimum cut on v0radius"};//minium decay length - irrelevant bc i set 4 in pre-selection
  Configurable<float> MaxV0radius{"MaxV0radius", 200, "maximum cut on v0radius"};

  Configurable<int> nBins{"nBins", 200, "N bins in histos"};
  Configurable<int> nBinsPt{"nBinsPt", 200, "N bins in pT histos"};
  Configurable<int> nBinsEta{"nBinsEta", 200, "N bins in Eta histos"};
  Configurable<int> nBinsPhi{"nBinsPhi", 200, "N bins in Phi histos"};

  float fiducialVolume;                        // 0.9 - jetR
  std::vector<fastjet::PseudoJet> jetConstituents;
  std::vector<fastjet::PseudoJet> jetReclustered;
  JetFinder jetReclusterer;

  Filter collisionFilter = nabs(aod::collision::posZ) < vertexZCut;
  Filter trackFilter = (nabs(aod::track::eta) < trackEtaCut) && (requireGlobalTrackInFilter()) && (aod::track::pt > trackPtCutMin);//here also a max and then i can evaluate in different pt bins - see literature for proper ranges
  Filter preFilterV0 = nabs(aod::v0data::dcapostopv) > dcapostopv&& nabs(aod::v0data::dcanegtopv) > dcanegtopv&& aod::v0data::dcaV0daughters < dcav0dau;
  
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

  HistogramRegistry registry{
    "registry",
    {
      {"hCollVtxZ", "hCollisionVtxZ", {HistType::kTH1F, {{nBins, -15., 15.}}}},
      {"hV0radius", "hV0radius", {HistType::kTH1F, {{nBins, 0., 200.}}}},
      {"hV0cospa", "hV0cospa", {HistType::kTH1F, {{nBins, 0.8, 1.1}}}},

      {"hMK0Short", "hMK0Short; M (GeV/#it{c})", {HistType::kTH1F, {{nBins, 0, 4}}}},
      {"hPtK0Short", "if K0Short : v0.pt()  p_{T}; p_{T} (GeV/#it{c})", {HistType::kTH1F, {{nBinsPt, 0, 20}}}},
      {"hEtaK0Short", "if K0Short : v0.eta() #eta; #eta", {HistType::kTH1F, {{nBinsEta, -0.9, 0.9}}}},
      {"hPhiK0Short", "if K0Short : v0.phi() #phi; #phi", {HistType::kTH1F, {{nBinsPhi, 0, 6.3}}}},

      {"hMLambda", "hMLambda; M (GeV/#it{c})", {HistType::kTH1F, {{nBins, 0, 4}}}},
      {"hPtLambda", "if Lambda : v0.pt()  p_{T}; p_{T} (GeV/#it{c})", {HistType::kTH1F, {{nBinsPt, 0, 20}}}},
      {"hEtaLambda", "if Lambda : v0.eta() #eta; #eta", {HistType::kTH1F, {{nBinsEta, -0.9, 0.9}}}},
      {"hPhiLambda", "if Lambda : v0.phi() #phi; #phi", {HistType::kTH1F, {{nBinsPhi, 0, 6.3}}}},

      {"hMAntiLambda", "hMAntiLambda; M (GeV/#it{c})", {HistType::kTH1F, {{nBins, 0, 4}}}},
      {"hPtAntiLambda", "if AntiLambda : v0.pt()  p_{T}; p_{T} (GeV/#it{c})", {HistType::kTH1F, {{nBinsPt, 0, 20}}}},
      {"hEtaAntiLambda", "if AntiLambda : v0.eta() #eta; #eta", {HistType::kTH1F, {{nBinsEta, -0.9, 0.9}}}},
      {"hPhiAntiLambda", "if AntiLambda : v0.phi() #phi; #phi", {HistType::kTH1F, {{nBinsPhi, 0, 6.3}}}},

      //Daughter QA
      {"hKPtPosPion", "Koan PosPion ; p_{T} (GeV/#it{c})", {HistType::kTH1F, {{nBinsPt, 0, 100}}}},
      {"hKPtNegPion", "Kaon NegPion ; p_{T} (GeV/#it{c})", {HistType::kTH1F, {{nBinsPt, 0, 100}}}},
      {"hLPtPosPr", "Lambda PosPr ; p_{T} (GeV/#it{c})", {HistType::kTH1F, {{nBinsPt, 0, 100}}}},
      {"hLPtNegPi", "Lambda NegPr ; p_{T} (GeV/#it{c})", {HistType::kTH1F, {{nBinsPt, 0, 100}}}},//ofc. there is no neg proton, pos/neg refers to the tracks
      {"hALPtPosPion", "AntiLambda PosPion ; p_{T} (GeV/#it{c})", {HistType::kTH1F, {{nBinsPt, 0, 100}}}},
      {"hALPtNegPr", "AntiLambda NegPr ; p_{T} (GeV/#it{c})", {HistType::kTH1F, {{nBinsPt, 0, 100}}}},
      {"hKEtaPosPion", "Kaon PosPion; #eta", {HistType::kTH1F, {{nBinsEta, -0.9, 0.9}}}},
      {"hKEtaNegPion", "Kaon NegPion; #eta", {HistType::kTH1F, {{nBinsEta, -0.9, 0.9}}}},
      {"hLEtaPosPr", "Lambda PosPr; #eta", {HistType::kTH1F, {{nBinsEta, -0.9, 0.9}}}},
      {"hLEtaNegPi", "Lambda NegPr; #eta", {HistType::kTH1F, {{nBinsEta, -0.9, 0.9}}}},
      {"hALEtaPosPion", "AntiLambda PosPion; #eta", {HistType::kTH1F, {{nBinsEta, -0.9, 0.9}}}},
      {"hALEtaNegPr", "AntiLambda NegPr; #eta", {HistType::kTH1F, {{nBinsEta, -0.9, 0.9}}}},
      {"hKPhiPosPion", "Kaon PosPion; #phi", {HistType::kTH1F, {{nBinsPhi, 0, 6.3}}}},
      {"hKPhiNegPion", "Kaon NegPion; #phi", {HistType::kTH1F, {{nBinsPhi, 0, 6.3}}}},
      {"hLPhiPosPr", "Lambda PosPr; #phi", {HistType::kTH1F, {{nBinsPhi, 0, 6.3}}}},
      {"hLPhiNegPi", "Lambda NegPr; #phi", {HistType::kTH1F, {{nBinsPhi, 0, 6.3}}}},
      {"hALPhiPosPion", "AntiLambda PosPion; #phi", {HistType::kTH1F, {{nBinsPhi, 0, 6.3}}}},
      {"hALPhiNegPr", "AntiLambda NegPr; #phi", {HistType::kTH1F, {{nBinsPhi, 0, 6.3}}}},
      //Analysis plots - adjust binning and normalize that stuff ... something with integral ... move this whole shabang to the plotting script
      {"LambdaOverKaonPt", "(Lambda + AntiLambda)/2K p_{T}; p_{T} (GeV/#it{c})", {HistType::kTH1F, {{nBinsPt, 0, 20}}}},

      //control plots	
      {"MdeltaPhi", "Martas #Delta #phi; #phi", {HistType::kTH1F, {{nBinsPhi, -3.2, 6.3}}}},
      {"JdeltaPhi", "Johannas #Delta #phi; #phi", {HistType::kTH1F, {{nBinsPhi, -3.2, 6.3}}}},
      {"hdeltaEta", "#Delta #eta; #eta", {HistType::kTH1F, {{nBinsEta, -2, 2}}}},
      {"V0CollIdVsGlobIndex", " V0.CollId vs. coll.GlobalIndex; V0.CollId; coll.GlobalIndex", {HistType::kTH2F, {{100, 0, 100}, {100, 0, 100}}}},
      {"V0CollIdVsV0GlobIndex", " V0.CollId vs. V0.GlobalIndex; V0.CollId; V0.GlobalIndex", {HistType::kTH2F, {{100, 0, 100}, {100, 0, 100}}}},
      {"V0GlobVsCollGlobIndex", " V0.GlobalIndex vs. coll.GlobalIndex; V0.GlobalIndex; coll.GlobalIndex", {HistType::kTH2F, {{100, 0, 100}, {100, 0, 100}}}},

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
    }
  };

  //I think a second struct where i process jets withit the dca filter is required.. then I could compare jets with and without ?..hmm..
  void Jet(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::JetFilters>>::iterator const& collision, soa::Filtered<CombinedTracks> const& tracks, aod::V0Datas const& V0s)//figure out why only aod::Tracks or MyTracks and/or if they are th same
  { 
  //same initial collision requirement as in V0 process -> Event selection decision based on V0A & V0C
   if(!collision.sel7()){return;}
   if (collision.hasJetChHighPt() >= bTriggerDecision) {//0=MB Event, 1=Event selected by EPN
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
          continue; // skip bad quality tracks in even selection helper task -> o2::aod::track::IsQualityTrack |	D |	isQualityTrack |	bool
        }

       registry.fill(HIST("JetTrackPt"), trk.pt());
       registry.fill(HIST("JetTrackEta"),trk.eta());
       registry.fill(HIST("JetTrackPhi"),trk.phi());
       
       fillConstituents(trk,jetConstituents); // ./PWGJE/Core/JetFinder.h - maybe, if I add the V0 shabang, I could fill from tracks of V0 within v0 radius ?
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
      if (leadingJetPt > -1.) {// at some point i should add the pid too !
        registry.fill(HIST("JetLeadJetPt"), leadingJetPt);
        registry.fill(HIST("JetLeadJetPhi"), leadingJetPhi);
        registry.fill(HIST("JetLeadJetEta"), leadingJetEta);
        for(auto& v0 : V0s){
	      //if(v0.collisionId() == collision.globalIndex()){
          registry.fill(HIST("V0CollIdVsGlobIndex"), v0.collisionId(), collision.globalIndex());
          registry.fill(HIST("V0CollIdVsV0GlobIndex"), v0.globalIndex(), collision.globalIndex());
          registry.fill(HIST("V0GlobVsCollGlobIndex"), v0.collisionId(), v0.globalIndex());
          //if(v0.globalIndex() == collision.globalIndex()){
          // from V0 table o2::aod::v0::CollisionId |	I |	collisionId |	int32 |	Collision index oooor: o2::soa::Index |	GI |	globalIndex |	int64_t
          // from collision table o2::soa::Index |	GI |	globalIndex |	int64_t
          if(v0.collisionId() == collision.globalIndex()){
            registry.fill(HIST("jetWithV0Pt"), v0.pt());
            registry.fill(HIST("jetWithV0Eta"), v0.eta());
            registry.fill(HIST("jetWithV0Phi"), v0.phi());
            float dPhi = leadingJetPhi - v0.phi();//fast jet uses -pi to pi azimuthal range while Alice uses 0 to 2pi 
            if(dPhi > fastjet::pi ){ dPhi -= 2*fastjet::pi;}
            if(dPhi < -fastjet::pi ){ dPhi += 2*fastjet::pi;}
            registry.fill(HIST("MdeltaPhi"), dPhi);
            registry.fill(HIST("hdeltaEta"), leadingJetEta-v0.eta() );
	          angularDistance = sqrt(pow(leadingJetEta-v0.eta(),2) + pow(dPhi,2) );//this is the correct dPhi that I use from now on
	          registry.fill(HIST("AngularDistance"), angularDistance);
  	     }
        }//end V0s    
       }//end if leading jet
      
    }//end collision

  }
  PROCESS_SWITCH(correlationvzerojets, Jet, "process jets", true);

  void V0(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& collision,  soa::Filtered<CombinedTracks> const& tracks, aod::V0Datas const& V0s){//no more my tracks
    if(!collision.sel7()){// sel8 is event selection based on TVX mainly for run3 data, but sel7=bool -> Event selection decision based on V0A & V0C -> https://aliceo2group.github.io/analysis-framework/docs/datamodel/helperTaskTables.html
      return;
    }

    registry.fill(HIST("hCollVtxZ"),collision.posZ()); // Inclusive Tracks from sel7 selections and the aod::pidTPCPi, aod::pidTPCPr
	
    for(auto& v0 : V0s){
      //particle identification by using the tpcNSigma track variable
      //in principle i could try the same with tofNSigma - if I include tof in my tracks ... good or bad -> gonna ask ...
      float nsigma_pos_proton = TMath::Abs(v0.posTrack_as<CombinedTracks>().tpcNSigmaPr());// o2::aod::pidtpc::TPCNSigmaPr 	|	tpcNSigmaPr |	float |	Nsigma separation with the TPC detector for proton
      float nsigma_neg_proton = TMath::Abs(v0.negTrack_as<CombinedTracks>().tpcNSigmaPr());// this is the TPC dE/dx, for 22 < 4; for 2017/2018 < 5 !
      float nsigma_pos_pion = TMath::Abs(v0.posTrack_as<CombinedTracks>().tpcNSigmaPi());// o2::aod::pidtpc::TPCNSigmaPi 	|	tpcNSigmaPi |	float |	Nsigma separation with the TPC detector for pion
      float nsigma_neg_pion = TMath::Abs(v0.negTrack_as<CombinedTracks>().tpcNSigmaPi());

      if( v0.posTrack_as<CombinedTracks>().pt() > V0daugPtMin && v0.negTrack_as<CombinedTracks>().pt() > V0daugPtMin){
        if(v0.v0radius() > MinV0radius && v0.v0radius() < MaxV0radius && v0.pt() < V0ptMax && v0.pt() > V0ptMin) {//v0 radius and pT bin cuts
        //get citted v0.radius distribution
        registry.fill(HIST("hV0radius"), v0.v0radius());
        registry.fill(HIST("hV0cospa"), v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()));

        //for invariant mass rejection as function of pT -- there are sources that select a fixed mass rejection but pt dependent cospa
        float upperLambda = 1.13688 + 0.00527838*v0.pt() + 0.084222*exp(-3.80595*v0.pt());
        float lowerLambda = 1.09501 - 0.00523272*v0.pt() - 0.075269*exp(-3.46339*v0.pt()); 
        float upperKaon = 0.563707 + 0.0114979*v0.pt();
        float lowerKaon = 0.43006 - 0.0110029*v0.pt();

        if( nsigma_pos_pion < 5 && nsigma_neg_pion < 5 && v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) > v0cospaK0s ){//topological daughter cut
          if(v0.mK0Short() > upperKaon || v0.mK0Short() < lowerKaon ){continue;}
          registry.fill(HIST("hMK0Short"), v0.mK0Short());
	        registry.fill(HIST("hPtK0Short"), v0.pt());
	        registry.fill(HIST("hEtaK0Short"), v0.eta());
          registry.fill(HIST("hPhiK0Short"), v0.phi());
          //for QA of daughters
          registry.fill(HIST("hKPtPosPion"), v0.posTrack_as<CombinedTracks>().pt()); 
          registry.fill(HIST("hKPtNegPion"), v0.negTrack_as<CombinedTracks>().pt()); 
          registry.fill(HIST("hKEtaPosPion"), v0.posTrack_as<CombinedTracks>().eta()); 
          registry.fill(HIST("hKEtaNegPion"), v0.negTrack_as<CombinedTracks>().eta()); 
          registry.fill(HIST("hKPhiPosPion"), v0.posTrack_as<CombinedTracks>().phi()); 
          registry.fill(HIST("hKPhiNegPion"), v0.negTrack_as<CombinedTracks>().phi()); 
       	}
	      if( nsigma_pos_proton < 5 && nsigma_neg_pion < 5 && v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) > v0cospaLamb ){
          if(v0.mLambda() > upperLambda || v0.mLambda() < lowerLambda ){continue;}
          registry.fill(HIST("hMLambda"), v0.mLambda());
	        registry.fill(HIST("hPtLambda"), v0.pt());
          registry.fill(HIST("hEtaLambda"), v0.eta());
          registry.fill(HIST("hPhiLambda"), v0.phi());
          //for QA of daughters -- this has to be corrected !
          registry.fill(HIST("hLPtPosPr"), v0.posTrack_as<CombinedTracks>().pt()); 
          registry.fill(HIST("hLPtNegPi"), v0.negTrack_as<CombinedTracks>().pt()); 
          registry.fill(HIST("hLEtaPosPr"), v0.posTrack_as<CombinedTracks>().eta()); 
          registry.fill(HIST("hLEtaNegPi"), v0.negTrack_as<CombinedTracks>().eta()); 
          registry.fill(HIST("hLPhiPosPr"), v0.posTrack_as<CombinedTracks>().phi()); 
          registry.fill(HIST("hLPhiNegPi"), v0.negTrack_as<CombinedTracks>().phi()); 
	      }
	      if( nsigma_pos_pion < 5 && nsigma_neg_proton < 5 && v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) > v0cospaLamb ){
          if(v0.mAntiLambda() > upperLambda || v0.mAntiLambda() < lowerLambda ){continue;}
          registry.fill(HIST("hMAntiLambda"), v0.mAntiLambda());	
	        registry.fill(HIST("hPtAntiLambda"), v0.pt());
          registry.fill(HIST("hEtaAntiLambda"), v0.eta());
          registry.fill(HIST("hPhiAntiLambda"), v0.phi());
          //for QA of daughters - maybe V0 specific plots
          registry.fill(HIST("hALPtPosPion"), v0.posTrack_as<CombinedTracks>().pt()); 
          registry.fill(HIST("hALPtNegPr"), v0.negTrack_as<CombinedTracks>().pt()); 
          registry.fill(HIST("hALEtaPosPion"), v0.posTrack_as<CombinedTracks>().eta()); 
          registry.fill(HIST("hALEtaNegPr"), v0.negTrack_as<CombinedTracks>().eta()); 
          registry.fill(HIST("hALPhiPosPion"), v0.posTrack_as<CombinedTracks>().phi()); 
          registry.fill(HIST("hALPhiNegPr"), v0.negTrack_as<CombinedTracks>().phi()); 
        }
        // Get V0 within radius and v0cosPA
        registry.fill(HIST("hPtTrackV0inRadius"), v0.pt());
        registry.fill(HIST("hEtaTrackV0inRadius"), v0.eta());
        registry.fill(HIST("hPhiTrackV0inRadius"), v0.phi());
	    }//if in v0 radius
	    registry.fill(HIST("hPtV0"), v0.pt());
	    registry.fill(HIST("hEtaV0"), v0.eta());
      registry.fill(HIST("hPhiV0"), v0.phi());
    }//if in pT daughter range
    }//end of V0's	
 
    for(auto& track : tracks){//Check kinematics of tracks
      if (!track.isQualityTrack()) {
        continue; // skip bad quality tracks just like in the beginning of the jet process
      }
      registry.fill(HIST("hTrackPt"), track.pt());
      registry.fill(HIST("hTrackEta"), track.eta());
      registry.fill(HIST("hTrackPhi"), track.phi());
    }
   }
  PROCESS_SWITCH(correlationvzerojets, V0, "process v0 and their track QA", true);
};
  
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{//run all via , separation
  // adaptAnalysisTask<JetQATask>(cfgc)
     adaptAnalysisTask<correlationvzerojets>(cfgc)
  };
}
