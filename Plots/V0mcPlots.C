void EffAndRes2D(TH2F *pt, TH2F *eta, TH2F *phi, TString varname, TString objname){
    TCanvas *can = new TCanvas("ResEff", "ResEff", 900, 600);
    TPad *pad = new TPad("pad", "pad", 0., 0., 0.3, 1);
    pad->SetBottomMargin(0.1); 
    pad->SetLeftMargin(0.15);
	pad->SetRightMargin(0.15);
    pad->Draw();            
    pad->cd(); 
    pt->Draw("COLZ");
    can->cd();
    TPad *pad2 = new TPad("pad", "pad", 0.3, 0., 0.6, 1);
    pad2->SetBottomMargin(0.1); 
    pad2->SetLeftMargin(0.15);
	pad2->SetRightMargin(0.15);
    pad2->Draw();            
    pad2->cd(); 
    eta->Draw("COLZ");
    can->cd();
    TPad *pad3 = new TPad("pad", "pad", 0.6, 0., 0.9, 1);
    pad3->SetBottomMargin(0.1); 
    pad3->SetLeftMargin(0.15);
	pad3->SetRightMargin(0.15);
    pad3->Draw();            
    pad3->cd(); 
    phi->Draw("COLZ");
    can->SaveAs("results/"+varname+"/"+objname+"_2D.pdf");
}

void plotPerKinBin(TH2F *VarVsPt, int nBins, TString kin, TString varname){
    for(int i = 0; i < nBins; i++){
        TH1D* hx = VarVsPt->ProjectionX("hx",i);
        TH1D* hy = VarVsPt->ProjectionY("hy",i);
        TCanvas *can = new TCanvas("can", "projection", 800, 400);
        can->Divide(2,1);    
        can->cd(1);
        hx->SetTitle(Form("X projection of "+varname+" in "+kin+" bin No.%d",i));
        hx->GetYaxis()->SetTitle("Number of entries");
        hx->Draw();
        can->cd(2);
        hy->SetTitle(Form("Bin No.%d Y projection of "+varname+" in "+kin,i));
        hy->GetYaxis()->SetTitle("Number of entries");
        hy->Draw();
        can->SaveAs(Form("results/"+varname+"/"+kin+"/Bin_%d.pdf",i));
    }
}


void plotProjection(TH2F *VarVsKin, TString kin, TString varname){
    TH1D* hx = VarVsKin->ProjectionX("hx");
    TH1D* hy = VarVsKin->ProjectionY("hy");
    TCanvas *can = new TCanvas("can", "projection", 800, 400);
    can->Divide(2,1);    
    can->cd(1);
    hx->SetTitle("X projection of "+varname);
    hx->GetYaxis()->SetTitle("Number of entries");
    hx->Draw();
    can->cd(2);
    hy->SetTitle("Y projection of "+varname);
    hy->GetYaxis()->SetTitle("Number of entries");
    hy->Draw();
    can->SaveAs("results/"+varname+"/"+kin+"_projection.pdf");
        
}

//Oh well..nach muede kommt doof - moin johanna 
//trotzdem sollte hier noch rungekabeled werden und dann noch baryon over meson fuer den vergleich mit pp
void plotVariablePerObservable(TH2F *VarVsKin, TString kin, TString varname){
    TH1D* hx = VarVsKin->ProjectionX("hx");
    TH1D* hy = VarVsKin->ProjectionY("hy");
    TH1D* h = VarVsKin->ProjectionX("h");

    for(int i = 0; i < hx->GetNbinsX(); i++){
        float val = 0;
        for(int j = 0; j < hy->GetNbinsX(); j++){
            float cont = VarVsKin->GetBinContent(i,j);
            if(hy->GetBinContent(j) == 0){continue;}
            val = val+cont/hy->Integral();
            if(val>0){cout<<val<<endl;}   
        }
        h->SetBinContent(i ,val);
    }
    TCanvas *can = new TCanvas("can", "projection", 800, 400);
    //h->SetTitle(Form("X projection of "+varname+" in "+kin+" bin No.%d",i));
    h->GetYaxis()->SetTitle("varname");
    h->Draw();
    can->SaveAs("results/"+varname+"/"+kin+".pdf");
}


//and we can also look at the baryon over meson ratio for the kinematic selection and after adding pT dependent true information also for this
//-> these results should be consitent with 1/2 bc this simulated dataset was pp run 3.
void V0mcPlots(){
    TFile *Result = new TFile("/home/johannalomker/alice/analysis/jetQuenching/McRun3setup/AnalysisResults.root");
    bool perBin = false;
    /*      
      {"","#Delta p_{T} = MC truth V0 - V0; MC truth V0 p_{T} (GeV/#it{c}); #Delta p_{T}", {HistType::kTH2F, {{nBinsPt, 0, 20}, {nBins, -3, 3}}}},
      {"","#Delta #eta = MC truth V0 - V0; MC truth V0 #eta; #Delta #eta", {HistType::kTH2F, {{nBinsPt, -0.9, 0.9}, {nBins, -1, 1}}}},
      {"","#Delta #phi = MC truth V0 - V0; MC truth V0 #phi (GeV/#it{c}); #Delta #phi", {HistType::kTH2F, {{nBinsPt, 0, 6.32}, {nBins, -5, 5}}}},
      
      {"","Eff.(p_{T}) = V0 / MC truth V0; MC truth V0 p_{T} (GeV/#it{c}); Eff.(p_{T})", {HistType::kTH2F, {{nBinsPt, 0, 20}, {nBins, 0, 2}}}},
      {"","Eff.(#eta) = V0 / MC truth V0; MC truth V0 #eta; Eff.(#eta)", {HistType::kTH2F, {{nBinsEta, -0.9, 0.9}, {nBins, 0, 2}}}},
      {"","Eff.(#phi) = V0 / MC truth V0; MC truth V0 #phi; Eff.(#phi)", {HistType::kTH2F, {{nBinsPhi, 0, 6.32}, {nBins, 0, 2}}}},

      
      {"","#Delta p_{T} = MC Track - Track; MC truth track p_{T} (GeV/#it{c}); #Delta p_{T}", {HistType::kTH2F, {{nBinsPt, 0, 20}, {nBins, -10, 10}}}},
      {"","#Delta #eta = MC Track - Track; MC truth track #eta; #Delta #eta", {HistType::kTH2F, {{nBinsEta, -0.9, 0.9}, {nBins, -1, 1}}}},
      {"","#Delta #phi = MC Track - Track; MC truth track #phi; #Delta #phi", {HistType::kTH2F, {{nBinsPhi, 0, 6.32}, {nBins, -5, 5}}}},

      {"effTpt"," Eff.(p_{T}) = Track / MC truth track; MC truth track p_{T} (GeV/#it{c}); Eff.(p_{T})", {HistType::kTH2F, {{nBinsPt, 0, 20}, {nBins, 0, 2}}}},
      {"effTeta"," Eff.(#eta) = Track / MC truth track; MC truth track #eta; Eff.(#eta)", {HistType::kTH2F, {{nBinsEta, -0.9, 0.9}, {nBins, 0, 2}}}},
      {"effTphi"," Eff.(#phi) = Track / MC truth track; MC truth track #phi; Eff.(#phi)", {HistType::kTH2F, {{nBinsPhi, 0, 6.32}, {nBins, 0, 2}}}},
    */
    //Resolution plots V0
    TH2F *resoV0pt = (TH2F*) Result->Get("correlationvzerojets/resoV0pt");
    TH2F *resoV0eta = (TH2F*) Result->Get("correlationvzerojets/resoV0eta");
    TH2F *resoV0phi = (TH2F*) Result->Get("correlationvzerojets/resoV0phi");
    //Resolution plots Tracks
    TH2F *resoTpt = (TH2F*) Result->Get("correlationvzerojets/resoTpt");
    TH2F *resoTeta = (TH2F*) Result->Get("correlationvzerojets/resoTeta");
    TH2F *resoTphi = (TH2F*) Result->Get("correlationvzerojets/resoTphi");
    //Efficiency plots V0
    TH2F *effV0pt = (TH2F*) Result->Get("correlationvzerojets/effV0pt");
    TH2F *effV0eta = (TH2F*) Result->Get("correlationvzerojets/effV0eta");
    TH2F *effV0phi = (TH2F*) Result->Get("correlationvzerojets/effV0phi");
    //Efficiency plots Tracks
    TH2F *effTpt = (TH2F*) Result->Get("correlationvzerojets/effTpt");
    TH2F *effTeta = (TH2F*) Result->Get("correlationvzerojets/effTeta");
    TH2F *effTphi = (TH2F*) Result->Get("correlationvzerojets/effTphi");
    //1) 2D histos as a function 
    //V0
    EffAndRes2D(resoV0pt, resoV0eta , resoV0phi, "resolution", "V0");
    EffAndRes2D(effV0pt, effV0eta , effV0phi, "efficiency", "V0");
    //Track
    EffAndRes2D(resoTpt, resoTeta , resoTphi, "resolution", "track");
    EffAndRes2D(effTpt, effTeta , effTphi, "efficiency", "track");

    int Nbins = 50;
    //2) projections per kinematic observable bin
    //V0
    if(perBin == true){
    plotPerKinBin(resoV0pt, Nbins, "V0pT", "resolution");
    plotPerKinBin(resoV0eta, Nbins, "V0eta", "resolution");
    plotPerKinBin(resoV0phi, Nbins, "V0phi",  "resolution");

    plotPerKinBin(effV0pt, Nbins, "V0pT", "efficiency");
    plotPerKinBin(effV0eta, Nbins, "V0eta", "efficiency");
    plotPerKinBin(effV0phi, Nbins, "V0phi",  "efficiency");
    //Tracks
    plotPerKinBin(resoTpt, Nbins, "TpT", "resolution");
    plotPerKinBin(resoTeta, Nbins, "Teta", "resolution");
    plotPerKinBin(resoTphi, Nbins, "Tphi",  "resolution");

    plotPerKinBin(effTpt, Nbins, "TpT", "efficiency");
    plotPerKinBin(effTeta, Nbins, "Teta", "efficiency");
    plotPerKinBin(effTphi, Nbins, "Tphi",  "efficiency");
    }
    //3) pure projections 
    //V0
    plotProjection(resoV0pt, "V0pT", "resolution");
    plotProjection(resoV0eta, "V0eta", "resolution");
    plotProjection(resoV0phi, "V0phi", "resolution");

    plotProjection(effV0pt, "V0pT", "efficiency");
    plotProjection(effV0eta, "V0eta", "efficiency");
    plotProjection(effV0phi, "V0phi",  "efficiency"); 
    //Tracks
    plotProjection(resoTpt, "TpT", "resolution");
    plotProjection(resoTeta, "Teta", "resolution");
    plotProjection(resoTphi, "Tphi", "resolution");

    plotProjection(effTpt, "TpT", "efficiency");
    plotProjection(effTeta, "Teta", "efficiency");
    plotProjection(effTphi, "Tphi",  "efficiency");  

    //4) projection as function of pT (integrated -- fill histos in loop over pt)
    plotVariablePerObservable(resoV0pt, "V0pT", "resolution");
    plotVariablePerObservable(effV0pt, "V0pT", "efficiency");

}
