//Add the Lamber over kaon - normalize the ratio ! newPtb->Scale(1/Nb->Integral()) -> ok, not required anymore .. but why is the error so large ?;
void fullQaPlots(){
        TString input;
        //Output the Histos from correlationV0jet - the lambda over kaon, angular distance are not in here, but all qa's 
        TFile *AResult = new TFile("/home/johannalomker/alice/analysis/jetQuenching/AnalysisResultsWithJets.root");
        //TFile *AResult = new TFile("/home/johannalomker/alice/analysis/jetQuenching/AnalysisV0CollisionId.root");
        //TFile *AResult = new TFile("/home/johannalomker/alice/analysis/jetQuenching/AnalysisV0GlobalIndex.root");
        
        //Output the Histos from jet-filter
        TH2F *PtEtaGood2T = (TH2F*) AResult->Get("jet-filter/spectra/ptetaGoodTracks");
        TH2F *PtEtaRej2T = (TH2F*) AResult->Get("jet-filter/spectra/ptetaRejectedTracks");
        TH2F *PtEtaSel2CJ = (TH2F*) AResult->Get("jet-filter/spectra/ptetaJetChSelected");

        TH2F *PtPhiGood2T = (TH2F*) AResult->Get("jet-filter/spectra/ptphiGoodTracks");
        TH2F *PtPhiRej2T = (TH2F*) AResult->Get("jet-filter/spectra/ptphiRejectedTracks");
        TH2F *PtPhiSel2CJ = (TH2F*) AResult->Get("jet-filter/spectra/ptphiJetChSelected");

        TH1F *CollZpos = (TH1F*) AResult->Get("jet-filter/spectra/fCollZpos");
        TH1F *NProcessed = (TH1F*) AResult->Get("jet-filter/spectra/fProcessedEvents");

        TCanvas *c = new TCanvas("c0", "canvas", 800, 400);
        c->Divide(2);
        c->cd(1);
        CollZpos->Draw("E");
        c->cd(2);
        NProcessed->Draw("E");
        c->SaveAs("plots/NPosZ.pdf");

        //correlationvzerojets (no subdir)
        TH1F *VtxZ = (TH1F*) AResult->Get("correlationvzerojets/hCollVtxZ");
        TH1F *VtxZjets = (TH1F*) AResult->Get("correlationvzerojets/jetVtx");

        TCanvas *can = new TCanvas("can", "vtx", 800, 400);
        can->Divide(2,1);
        can->cd(1);
        VtxZjets->SetTitle("Vtx Z Jets process");
        VtxZjets->Draw("E");
        can->cd(2);
        VtxZ->SetTitle("Vtx Z V0 process");
        VtxZ->Draw("E");
        can->SaveAs("plots/VtxJetsV0.pdf");

        TH1F *V0radius = (TH1F*) AResult->Get("correlationvzerojets/hV0radius");
        TH1F *cosPA = (TH1F*) AResult->Get("correlationvzerojets/hV0cospa");
        TCanvas *can3 = new TCanvas("can3", "V0 cut variables", 800, 400);
        can3->Divide(2,1);
        can3->cd(1);
        V0radius->Draw("E");
        can3->cd(2);
        cosPA->Draw("E");
        can3->SaveAs("plots/V0radiusAndCospa.pdf");

        //to check the dPhi and dEta calculation in jets!
        TH1F *J = (TH1F*) AResult->Get("correlationvzerojets/hdeltaEta");
        TH1F *M = (TH1F*) AResult->Get("correlationvzerojets/MdeltaPhi");
        TCanvas *can4 = new TCanvas("can4", "dPhi", 800, 400);
        can4->Divide(2,1);
        can4->cd(1);
        J->SetTitle("#Delta #eta");
        J->Draw("E");
        can4->cd(2);
        M->SetTitle("#Delta #Phi");
        M->Draw("E");
        can4->SaveAs("plots/DeltaPhiEta.pdf");

        //to understand what these indeces in jets are
        TH2F *VColGlob = (TH2F*) AResult->Get("correlationvzerojets/V0CollIdVsGlobIndex");
        TH2F *VColVGlob = (TH2F*) AResult->Get("correlationvzerojets/V0CollIdVsV0GlobIndex");
        TH2F *VGlobGlob = (TH2F*) AResult->Get("correlationvzerojets/V0GlobVsCollGlobIndex");
        TCanvas *V = new TCanvas("V", "Index", 900, 400);
        V->Divide(3,1);
        V->cd(1);
        VColGlob->Draw("COLZ");
        V->cd(2);
        VColVGlob->Draw("COLZ");
        V->cd(3);
        VGlobGlob->Draw("COLZ");
        V->SaveAs("plots/IndexMystery.pdf");        

        gStyle -> SetOptStat(0);
        //V0 observables (Mass, pT, phi, eta)
        TH1F *MK = (TH1F*) AResult->Get("correlationvzerojets/hMK0Short");
        TH1F *ML = (TH1F*) AResult->Get("correlationvzerojets/hMLambda");
        TH1F *MAL = (TH1F*) AResult->Get("correlationvzerojets/hMAntiLambda");

        TH1F *KPt = (TH1F*) AResult->Get("correlationvzerojets/hPtK0Short");
        TH1F *LPt = (TH1F*) AResult->Get("correlationvzerojets/hPtLambda");
        TH1F *ALPt = (TH1F*) AResult->Get("correlationvzerojets/hPtAntiLambda");

        TH1F *KPhi = (TH1F*) AResult->Get("correlationvzerojets/hPhiK0Short");
        TH1F *LPhi = (TH1F*) AResult->Get("correlationvzerojets/hPhiLambda");
        TH1F *ALPhi = (TH1F*) AResult->Get("correlationvzerojets/hPhiAntiLambda");

        TH1F *KEta = (TH1F*) AResult->Get("correlationvzerojets/hEtaK0Short");
        TH1F *LEta = (TH1F*) AResult->Get("correlationvzerojets/hEtaLambda");
        TH1F *ALEta = (TH1F*) AResult->Get("correlationvzerojets/hEtaAntiLambda");

        auto L1 = new TLegend(0.2,0.49, 0.8,0.51);
        L1->SetHeader("","C");
        L1->SetNColumns(3);
        L1->SetTextSize(0.026);
        L1->AddEntry(ML, "Lambda", "lep");
        L1->AddEntry(MAL, "AntiLambda", "lep");
        L1->AddEntry(MK, "K0Short", "lep");
        L1->SetBorderSize(0);
        L1->SetFillStyle(0);

        TCanvas *c1 = new TCanvas("c1", "canvas", 800, 800);
        c1->Divide(2,2);
        c1->cd(1);
        MK->SetTitle("");
        MK->SetLineColor(3);
        MK->GetXaxis()->SetRangeUser(0,3);
        MK->Draw();
        MAL->SetLineColor(2);
        MAL->Draw("same");
        ML->SetLineColor(1);
        ML->Draw("same");
        c1->cd(2);
        KPt->SetTitle("");
        KPt->GetXaxis()->SetRangeUser(0,10);
        KPt->SetLineColor(3);
        KPt->Draw();
        LPt->SetLineColor(1);
        LPt->Draw("same");
        ALPt->SetLineColor(2);
        ALPt->Draw("same");
        c1->cd(3); 
        c1->SetLogy();
        KPhi->GetYaxis()->SetRangeUser(0, 2500);
        KPhi->SetTitle("");
        KPhi->SetLineColor(3);
        KPhi->Draw();
        c1->SetLogy();
        LPhi->SetLineColor(1);
        LPhi->Draw("same");
        ALPhi->SetLineColor(2);
        ALPhi->Draw("same");
        c1->cd(4);
        KEta->SetTitle("");
        KEta->SetLineColor(3);
        KEta->Draw();
        LEta->SetLineColor(1);
        LEta->Draw("same");
        ALEta->SetLineColor(2);
        ALEta->Draw("same");
        c1->cd();
        L1->Draw(" ");
        c1->SaveAs("plots/V0Candidates.pdf");

        //here I want to add my daughter track QA
        TH1F *KPiPospt = (TH1F*) AResult->Get("correlationvzerojets/hKPtPosPion");
        TH1F *KPiNegpt = (TH1F*) AResult->Get("correlationvzerojets/hKPtNegPion");
        TH1F *KPiPoseta = (TH1F*) AResult->Get("correlationvzerojets/hKEtaPosPion");
        TH1F *KPiNegeta = (TH1F*) AResult->Get("correlationvzerojets/hKEtaNegPion");
        TH1F *KPiPosphi = (TH1F*) AResult->Get("correlationvzerojets/hKPhiPosPion");
        TH1F *KPiNegphi = (TH1F*) AResult->Get("correlationvzerojets/hKPhiNegPion");

        auto L4 = new TLegend(0.2,0.9, 0.8,1);//for each mother one legend
        L4->SetHeader("Kaon Daughters","C");
        L4->SetNColumns(4);
        L4->SetTextSize(0.036);
        L4->AddEntry(KPiPospt, "Positive #pi tracks", "lep");
        L4->AddEntry(KPiNegpt, "Negative #pi tracks", "lep");
        L4->SetBorderSize(0);
        L4->SetFillStyle(0);

        TH1F *LPrPospt = (TH1F*) AResult->Get("correlationvzerojets/hLPtPosPr");
        TH1F *LPiNegpt = (TH1F*) AResult->Get("correlationvzerojets/hLPtNegPi");
        TH1F *LPrPoseta = (TH1F*) AResult->Get("correlationvzerojets/hLEtaPosPr");
        TH1F *LPiNegeta = (TH1F*) AResult->Get("correlationvzerojets/hLEtaNegPi");
        TH1F *LPrPosphi = (TH1F*) AResult->Get("correlationvzerojets/hLPhiPosPr");
        TH1F *LPiNegphi = (TH1F*) AResult->Get("correlationvzerojets/hLPhiNegPi");

        auto L41 = new TLegend(0.2,0.9, 0.8,1);
        L41->SetHeader("Lambda Daughters","C");
        L41->SetNColumns(4);
        L41->SetTextSize(0.036);
        L41->AddEntry(LPrPospt, "Positive p tracks", "lep");
        L41->AddEntry(LPiNegpt, "Negative #pi tracks", "lep");
        L41->SetBorderSize(0);
        L41->SetFillStyle(0);

        TH1F *ALPiPospt = (TH1F*) AResult->Get("correlationvzerojets/hALPtPosPion");
        TH1F *ALPrNegpt = (TH1F*) AResult->Get("correlationvzerojets/hALPtNegPr");
        TH1F *ALPiPoseta = (TH1F*) AResult->Get("correlationvzerojets/hALEtaPosPion");
        TH1F *ALPrNegeta = (TH1F*) AResult->Get("correlationvzerojets/hALEtaNegPr");
        TH1F *ALPiPosphi = (TH1F*) AResult->Get("correlationvzerojets/hALPhiPosPion");
        TH1F *ALPrNegphi = (TH1F*) AResult->Get("correlationvzerojets/hALPhiNegPr");

        auto L42 = new TLegend(0.2,0.9, 0.8,1);
        L42->SetHeader("Anti Lambda Daughters","C");
        L42->SetNColumns(4);
        L42->SetTextSize(0.036);
        L42->AddEntry(ALPiPospt, "Positive #pi tracks", "lep");
        L42->AddEntry(ALPrNegpt, "Negative p tracks", "lep");
        L42->SetBorderSize(0);
        L42->SetFillStyle(0);

        TCanvas *c4 = new TCanvas("c4", "V0 daughter tracks", 1200, 1200);
        c4->Divide(1,3);
        c4->cd(1);
        TPad *pad4 = new TPad("pad4", "pad4", 0., 0., 0.3, 0.9);
        pad4->SetBottomMargin(0.15); 
        pad4->SetLeftMargin(0.15);
	pad4->SetRightMargin(0.1);
        pad4->Draw();             // Draw the upper pad: pad1
        pad4->cd();               // pad1 becomes the current pad
	KPiPospt->GetXaxis()->SetRangeUser(0,15);
        //PiPospt->GetYaxis()->SetRangeUser(0,120000);
	KPiPospt->SetLineColor(3);
        KPiPospt->SetMarkerColor(3);
        KPiPospt->SetMarkerStyle(23);
        KPiPospt->SetTitle(" ");
        KPiPospt->DrawCopy();
        KPiNegpt->SetLineColor(8);
        KPiNegpt->SetMarkerStyle(26);
        KPiNegpt->SetMarkerColor(8);
        KPiNegpt->DrawCopy("Esame");
        c4->cd(1);
        TPad *p4 = new TPad("pa4", "pa4", 0.3, 0., 0.6, 0.9);
        p4->SetBottomMargin(0.15); 
        p4->SetLeftMargin(0.15);
	p4->SetRightMargin(0.1);
        p4->Draw();             
        p4->cd();               
	KPiPoseta->SetLineColor(3);
        KPiPoseta->SetMarkerStyle(23);
        KPiPoseta->SetMarkerColor(3);
        KPiPoseta->SetTitle(" ");
        KPiPoseta->DrawCopy();
        KPiNegeta->SetLineColor(8);
        KPiNegeta->SetMarkerStyle(26);
        KPiNegeta->SetMarkerColor(8);
        KPiNegeta->DrawCopy("Esame");
        c4->cd(1);
        TPad *pa4 = new TPad("pa4", "pa4", 0.6, 0., 0.9, 0.9);
        pa4->SetBottomMargin(0.15); 
        pa4->SetLeftMargin(0.15);
	pa4->SetRightMargin(0.1);
        pa4->Draw();             
        pa4->cd();               
        KPiPosphi->SetLineColor(3);
        KPiPosphi->SetMarkerStyle(23);
        KPiPosphi->SetMarkerColor(3);
        KPiPosphi->SetTitle(" ");
        KPiPosphi->DrawCopy();
        KPiNegphi->SetLineColor(8);
        KPiNegphi->SetMarkerStyle(26);
        KPiNegphi->SetMarkerColor(8);
        KPiNegphi->DrawCopy("Esame");
        c4->cd(1);
        L4->Draw();
        //                           Lambdas
        c4->cd(2);
        TPad *pad41 = new TPad("pad41", "pad41", 0., 0., 0.3, 0.9);
        pad41->SetBottomMargin(0.15); 
        pad41->SetLeftMargin(0.15);
	pad41->SetRightMargin(0.1);
        pad41->Draw();             
        pad41->cd();               
        LPrPospt->GetXaxis()->SetRangeUser(0,15);
	LPrPospt->SetLineColor(1);
        LPrPospt->SetMarkerColor(1);
        LPrPospt->SetMarkerStyle(23);
        LPrPospt->SetTitle(" ");
        LPrPospt->DrawCopy();
        LPiNegpt->SetLineColor(13);
        LPiNegpt->SetMarkerStyle(26);
        LPiNegpt->SetMarkerColor(13);
        LPiNegpt->DrawCopy("Esame");
        c4->cd(2);
        TPad *p41 = new TPad("pa41", "pa41", 0.3, 0., 0.6, 0.9);
        p41->SetBottomMargin(0.15); 
        p41->SetLeftMargin(0.15);
	p41->SetRightMargin(0.1);
        p41->Draw();            
        p41->cd();              
	LPrPoseta->SetLineColor(1);
        LPrPoseta->SetMarkerColor(1);
        LPrPoseta->SetMarkerStyle(23);
        LPrPoseta->SetTitle(" ");
        LPrPoseta->DrawCopy();
        LPiNegeta->SetLineColor(13);
        LPiNegeta->SetMarkerStyle(26);
        LPiNegeta->SetMarkerColor(13);
        LPiNegeta->DrawCopy("Esame");
        c4->cd(2);
        TPad *pa41 = new TPad("pa41", "pa41", 0.6, 0., 0.9, 0.9);
        pa41->SetBottomMargin(0.15); 
        pa41->SetLeftMargin(0.15);
	pa41->SetRightMargin(0.1);
        pa41->Draw();             
        pa41->cd();                
	LPrPosphi->SetLineColor(1);
        LPrPosphi->SetMarkerColor(1);
        LPrPosphi->SetMarkerStyle(23);
        LPrPosphi->SetTitle(" ");
        LPrPosphi->DrawCopy();
        LPiNegphi->SetLineColor(13);
        LPiNegphi->SetMarkerStyle(26);
        LPiNegphi->SetMarkerColor(13);
        LPiNegphi->DrawCopy("Esame");
        c4->cd(2);
        L41->Draw();
        c4->cd(3);
        //                           Anti Lambdas
        TPad *pad42 = new TPad("pad42", "pad42", 0., 0., 0.3, 0.9);
        pad42->SetBottomMargin(0.15); 
        pad42->SetLeftMargin(0.15);
	pad42->SetRightMargin(0.1);
        pad42->Draw();             
        pad42->cd();               
        ALPiPospt->GetXaxis()->SetRangeUser(0,15);
	ALPiPospt->SetLineColor(2);
        ALPiPospt->SetMarkerColor(2);
        ALPiPospt->SetMarkerStyle(23);
        ALPiPospt->SetTitle(" ");
        ALPiPospt->DrawCopy();
        ALPrNegpt->SetLineColor(46);
        ALPrNegpt->SetMarkerStyle(26);
        ALPrNegpt->SetMarkerColor(46);
        ALPrNegpt->DrawCopy("Esame");
        c4->cd(3);
        TPad *p42 = new TPad("pa42", "pa42", 0.3, 0., 0.6, 0.9);
        p42->SetBottomMargin(0.15); 
        p42->SetLeftMargin(0.15);
	p42->SetRightMargin(0.1);
        p42->Draw();            
        p42->cd();              
	ALPiPoseta->SetLineColor(2);
        ALPiPoseta->SetMarkerColor(2);
        ALPiPoseta->SetMarkerStyle(23);
        ALPiPoseta->SetTitle(" ");
        ALPiPoseta->DrawCopy();
        ALPrNegeta->SetLineColor(46);
        ALPrNegeta->SetMarkerStyle(26);
        ALPrNegeta->SetMarkerColor(46);
        ALPrNegeta->DrawCopy("Esame");
        c4->cd(3);
        TPad *pa42 = new TPad("pa42", "pa42", 0.6, 0., 0.9, 0.9);
        pa42->SetBottomMargin(0.15); 
        pa42->SetLeftMargin(0.15);
	pa42->SetRightMargin(0.1);
        pa42->Draw();            
        pa42->cd();               
	ALPiPosphi->SetLineColor(2);
        ALPiPosphi->SetMarkerColor(2);
        ALPiPosphi->SetMarkerStyle(23);
        ALPiPosphi->SetTitle(" ");
        ALPiPosphi->DrawCopy();
        ALPrNegphi->SetLineColor(46);
        ALPrNegphi->SetMarkerStyle(26);
        ALPrNegphi->SetMarkerColor(46);
        ALPrNegphi->DrawCopy("Esame");
        c4->cd(3);
        L42->Draw();
        c4->SaveAs("plots/qaV0Daughters.pdf");

        TH1F *VRpt = (TH1F*) AResult->Get("correlationvzerojets/hPtTrackV0inRadius");
        TH1F *Vpt = (TH1F*) AResult->Get("correlationvzerojets/hPtV0");
        TH1F *Tpt = (TH1F*) AResult->Get("correlationvzerojets/hTrackPt");

        TH1F *VRe = (TH1F*) AResult->Get("correlationvzerojets/hEtaTrackV0inRadius");
        TH1F *Ve = (TH1F*) AResult->Get("correlationvzerojets/hEtaV0");
        TH1F *Te = (TH1F*) AResult->Get("correlationvzerojets/hTrackEta");

        TH1F *VRphi = (TH1F*) AResult->Get("correlationvzerojets/hPhiTrackV0inRadius");
        TH1F *Vphi = (TH1F*) AResult->Get("correlationvzerojets/hPhiV0");
        TH1F *Tphi = (TH1F*) AResult->Get("correlationvzerojets/hTrackPhi");

        //jets as next 
        TH1F *Jpt = (TH1F*) AResult->Get("correlationvzerojets/JetTrackPt");
        TH1F *LTpt = (TH1F*) AResult->Get("correlationvzerojets/JetLeadTrackPt");
        TH1F *LJpt = (TH1F*) AResult->Get("correlationvzerojets/JetLeadJetPt");
        TH1F *JVpt = (TH1F*) AResult->Get("correlationvzerojets/jetV0Pt");
        TH1F *JVCpt = (TH1F*) AResult->Get("correlationvzerojets/jetWithV0Pt"); 

        TH1F *Je = (TH1F*) AResult->Get("correlationvzerojets/JetTrackEta");
        TH1F *LTe = (TH1F*) AResult->Get("correlationvzerojets/JetLeadTrackEta");
        TH1F *LJe = (TH1F*) AResult->Get("correlationvzerojets/JetLeadJetEta");
        TH1F *JVe = (TH1F*) AResult->Get("correlationvzerojets/jetV0Eta");
        TH1F *JVCe = (TH1F*) AResult->Get("correlationvzerojets/jetWithV0Eta"); 

        TH1F *Jphi = (TH1F*) AResult->Get("correlationvzerojets/JetTrackPhi");
        TH1F *LTphi = (TH1F*) AResult->Get("correlationvzerojets/JetLeadTrackPhi");
        TH1F *LJphi = (TH1F*) AResult->Get("correlationvzerojets/JetLeadJetPhi");
        TH1F *JVphi = (TH1F*) AResult->Get("correlationvzerojets/jetV0Phi");
        TH1F *JVCphi = (TH1F*) AResult->Get("correlationvzerojets/jetWithV0Phi"); 
        
        auto L2 = new TLegend(0.2,0.49, 0.8,0.51);
        L2->SetHeader("","C");
        L2->SetNColumns(4);
        L2->SetTextSize(0.026);
        L2->AddEntry(VRpt, "hPtTrackV0inRadius", "lep");//from V0 process
        L2->AddEntry(Vpt, "hPtV0", "lep");//from V0 process
        L2->AddEntry(Tpt, "hTrackPt", "lep");//from V0 process
        L2->AddEntry(Jpt, "JetTrackPt", "lep");//from jet process
        L2->SetBorderSize(0);
        L2->SetFillStyle(0);
        
        TCanvas *c2 = new TCanvas("c2", "V0 process", 1200, 800); 
        c2->Divide(3,2);
        c2->cd(1);
        c2->SetTopMargin(0.1);
        Vpt->SetTitle("");
        Vpt->GetXaxis()->SetRangeUser(0,10);
        Vpt->SetLineColor(2);
        Vpt->SetMarkerStyle(26);
        Vpt->SetMarkerColor(2);
        Vpt->Draw("E");
        VRpt->SetLineColor(1);
        VRpt->SetMarkerStyle(32);
        VRpt->SetMarkerColor(1);
        VRpt->Draw("same");
        c2->cd(2);
        c2->SetTopMargin(0.1);
        Ve->SetTitle("");
        Ve->SetLineColor(1);
        Ve->SetMarkerStyle(26);
        Ve->SetMarkerColor(1);
        Ve->Draw("E");
        VRe->SetLineColor(2);
        VRe->SetMarkerStyle(32);
        VRe->SetMarkerColor(2);
        VRe->Draw("same");
        c2->cd(3);
        c2->SetTopMargin(0.1);
        Vphi->SetTitle("");
        Vphi->SetLineColor(1);
        Vphi->SetMarkerStyle(26);
        Vphi->SetMarkerColor(1);
        Vphi->Draw("E");
        VRphi->SetLineColor(2);
        VRphi->SetMarkerStyle(32);
        VRphi->SetMarkerColor(2);
        VRphi->Draw("same");
        c2->cd(4);
        Jpt->SetTitle("");
        Jpt->GetXaxis()->SetRangeUser(0,10);
        Jpt->SetLineColor(3);
        Jpt->SetMarkerStyle(24);
        Jpt->SetMarkerColor(3);
        Jpt->Draw();
        Tpt->SetLineColor(2);
        Tpt->SetMarkerStyle(24);
        Tpt->SetMarkerColor(2);
        Tpt->Draw("Esame");
        c2->cd(5);
        Je->SetTitle("");
        Je->SetLineColor(3);
        Je->SetMarkerStyle(24);
        Je->SetMarkerColor(3);
        Je->Draw();
        Te->SetLineColor(2);
        Te->SetMarkerStyle(24);
        Te->SetMarkerColor(2);
        Te->Draw("Esame");
        c2->cd(6); 
        Jphi->SetTitle("");       
        Jphi->SetLineColor(3);
        Jphi->SetMarkerStyle(24);
        Jphi->SetMarkerColor(3);
        Jphi->Draw();
        Tphi->SetTitle("");
        Tphi->SetLineColor(2);
        Tphi->SetMarkerStyle(24);
        Tphi->SetMarkerColor(2);
        Tphi->Draw("Esame");
        c2->cd();
        L2->Draw();
        c2->SaveAs("plots/QA.pdf");

        auto L3 = new TLegend(0.2,0.9, 0.8,1);
        L3->SetHeader("","C");
        L3->SetNColumns(4);
        L3->SetTextSize(0.036);
        L3->AddEntry(LTpt, "JetLeadTrack", "lep");//from jet process
        L3->AddEntry(LJpt, "JetLeadJet", "lep");//from jet process
        L3->AddEntry(JVpt, "jetV0", "lep");//from jet process
        L3->AddEntry(JVCpt, "jetWithV0", "lep");//from jet process
        L3->SetBorderSize(0);
        L3->SetFillStyle(0);

        TCanvas *c3 = new TCanvas("c3", "Jets and V0 in jetprocess", 1200, 400);
        c3->cd();
        TPad *pa = new TPad("pa", "pa", 0., 0., 0.3, 0.9);
        pa->SetBottomMargin(0.15); 
        pa->SetLeftMargin(0.15);
	pa->SetRightMargin(0.1);
        pa->Draw();             // Draw the upper pad: pad1
        pa->cd();               // pad1 becomes the current pad
	JVpt->GetXaxis()->SetRangeUser(0,10);
        //JVpt->GetYaxis()->SetRangeUser(0,120000);
	JVpt->SetLineColor(2);
        JVpt->SetMarkerStyle(26);
        JVpt->SetMarkerColor(2);
        JVpt->SetTitle(" ");
        JVpt->DrawCopy();
        LJpt->SetLineColor(3);
        LJpt->SetMarkerStyle(26);
        LJpt->SetMarkerColor(3);
        LJpt->DrawCopy("Esame");
        LTpt->SetLineColor(1);
        LTpt->SetMarkerStyle(23);
        LTpt->SetMarkerColor(1);
        LTpt->DrawCopy("Esame");
        JVCpt->SetLineColor(2);
        JVCpt->SetMarkerStyle(32);
        JVCpt->SetMarkerColor(2);
        JVCpt->DrawCopy("Esame");
        c3->cd();
        TPad *pa2 = new TPad("pa2", "pa2", 0.3, 0., 0.6, 0.9);
        pa2->SetBottomMargin(0.15); 
        pa2->SetLeftMargin(0.15);
	pa2->SetRightMargin(0.1);
        pa2->Draw();             
        pa2->cd();               
	JVe->SetLineColor(2);
        JVe->SetMarkerStyle(26);
        JVe->SetMarkerColor(2);
        JVe->SetTitle(" ");
        JVe->DrawCopy();
        LJe->SetLineColor(3);
        LJe->SetMarkerStyle(26);
        LJe->SetMarkerColor(3);
        LJe->DrawCopy("Esame");
        LTe->SetLineColor(1);
        LTe->SetMarkerStyle(23);
        LTe->SetMarkerColor(1);
        LTe->DrawCopy("Esame");
        JVCe->SetLineColor(2);
        JVCe->SetMarkerStyle(32);
        JVCe->SetMarkerColor(2);
        JVCe->DrawCopy("Esame");
        c3->cd();

        TPad *pa3 = new TPad("pa3", "pa3", 0.6, 0., 0.9, 0.9);//finish this plot in the train
        pa3->SetBottomMargin(0.15); 
        pa3->SetLeftMargin(0.15);
	pa3->SetRightMargin(0.1);
        pa3->Draw();             
        pa3->cd();               
        JVphi->GetYaxis()->SetRangeUser(0,5500);
        JVphi->SetLineColor(2);
        JVphi->SetMarkerStyle(26);
        JVphi->SetMarkerColor(2);
        JVphi->SetTitle(" ");
        JVphi->DrawCopy();
        LJphi->SetLineColor(3);
        LJphi->SetMarkerStyle(26);
        LJphi->SetMarkerColor(3);
        LJphi->DrawCopy("Esame");
        LTphi->SetLineColor(1);
        LTphi->SetMarkerStyle(23);
        LTphi->SetMarkerColor(1);
        LTphi->DrawCopy("Esame");
        JVCphi->SetLineColor(2);
        JVCphi->SetMarkerStyle(32);
        JVCphi->SetMarkerColor(2);
        JVCphi->DrawCopy("Esame");
        c3->cd();
        L3->Draw();

        c3->SaveAs("plots/qaJets.pdf");
}
