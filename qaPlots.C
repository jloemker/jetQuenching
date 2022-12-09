//Add the vertex/counting plots
//Add the angular distance
//Add the Lamber over kaon - normalize the ratio ! newPtb->Scale(1/Nb->Integral());
void qaPlots(){
        gStyle -> SetOptStat(0);
        TString input;
        //Output the Histos from correlationV0jet - the lambda over kaon, angular distance and vtx collision are not in here, but all qa's 
        TFile *AResult = new TFile("/home/johannalomker/alice/analysis/jetQuenching/ResultsSel8.root");
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
        MAL->SetTitle("");
        MAL->SetLineColor(2);
        MAL->GetXaxis()->SetRangeUser(0,3);
        MAL->Draw();
        MK->SetLineColor(3);
        MK->Draw("same");
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
        c1->cd(3);                              //here we need to modify the y scale bc, the other lines are there but rather low...
        c1->SetLogy();
        KPhi->GetYaxis()->SetRangeUser(0,1000);
        KPhi->SetTitle("");
        KPhi->SetLineColor(3);
        KPhi->GetYaxis()->SetMoreLogLabels(5);
        KPhi->Draw();
        c1->SetLogy();
        LPhi->GetYaxis()->SetMoreLogLabels();
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
        c1->SaveAs("Sel8_V0Candidates.pdf");

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
        
        TCanvas *c2 = new TCanvas("c2", "V0 process", 1200, 800); //they are on top of each other
        c2->Divide(3,2);
        c2->cd(1);
        c2->SetTopMargin(0.1);
        VRpt->SetTitle("");
        VRpt->GetXaxis()->SetRangeUser(0,10);
        VRpt->SetLineColor(1);
        VRpt->SetMarkerStyle(32);
        VRpt->SetMarkerColor(1);
        VRpt->Draw("E");
        Vpt->SetLineColor(2);
        Vpt->SetMarkerStyle(26);
        Vpt->SetMarkerColor(2);
        Vpt->Draw("same");
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
        c2->SaveAs("Sel8_QA.pdf");

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

        //here i need pads otherwise the legend will not work - once this is done: rerun the whole shabang with sel8 !
        TCanvas *c3 = new TCanvas("c3", "Jets and V0 in jetprocess", 1200, 400);//fix colors markers and add legends
        c3->cd();
        TPad *pa = new TPad("pa", "pa", 0., 0., 0.3, 0.9);
        pa->SetBottomMargin(0.15); 
        pa->SetLeftMargin(0.15);
	pa->SetRightMargin(0.1);
        pa->Draw();             // Draw the upper pad: pad1
        pa->cd();               // pad1 becomes the current pad
	LTpt->GetXaxis()->SetRangeUser(0,10);
        LTpt->GetYaxis()->SetRangeUser(0,120000);
	LTpt->SetLineColor(1);
        LTpt->SetMarkerStyle(23);
        LTpt->SetTitle(" ");
        LTpt->DrawCopy();
        LJpt->SetLineColor(3);
        LJpt->SetMarkerStyle(26);
        LJpt->SetMarkerColor(3);
        LJpt->DrawCopy("Esame");
        JVpt->SetLineColor(2);
        JVpt->SetMarkerStyle(26);
        JVpt->SetMarkerColor(2);
        JVpt->DrawCopy("Esame");
        JVCpt->SetLineColor(2);
        JVCpt->SetMarkerStyle(32);
        JVCpt->SetMarkerColor(2);
        JVCpt->DrawCopy("Esame");
 
        c3->cd();
        TPad *pa2 = new TPad("pa2", "pa2", 0.3, 0., 0.6, 0.9);
        pa2->SetBottomMargin(0.15); 
        pa2->SetLeftMargin(0.15);
	pa2->SetRightMargin(0.1);
        pa2->Draw();             // Draw the upper pad: pad1
        pa2->cd();               // pad1 becomes the current pad
        LTe->GetYaxis()->SetRangeUser(0,2200);
	LTe->SetLineColor(1);
        LTe->SetMarkerStyle(23);
        LTe->SetTitle(" ");
        LTe->DrawCopy();
        LJe->SetLineColor(3);
        LJe->SetMarkerStyle(26);
        LJe->SetMarkerColor(3);
        LJe->DrawCopy("Esame");
        JVe->SetLineColor(2);
        JVe->SetMarkerStyle(26);
        JVe->SetMarkerColor(2);
        JVe->DrawCopy("Esame");
        JVCe->SetLineColor(2);
        JVCe->SetMarkerStyle(32);
        JVCe->SetMarkerColor(2);
        JVCe->DrawCopy("Esame");
        c3->cd();

        TPad *pa3 = new TPad("pa3", "pa3", 0.6, 0., 0.9, 0.9);//finish this plot in the train
        pa3->SetBottomMargin(0.15); 
        pa3->SetLeftMargin(0.15);
	pa3->SetRightMargin(0.1);
        pa3->Draw();             // Draw the upper pad: pad1
        pa3->cd();               // pad1 becomes the current pad
        LTphi->GetYaxis()->SetRangeUser(0,1800);
        LTphi->SetLineColor(1);
        LTphi->SetMarkerStyle(23);
        LTphi->SetTitle(" ");
        LTphi->DrawCopy();
        LJphi->SetLineColor(3);
        LJphi->SetMarkerStyle(26);
        LJphi->SetMarkerColor(3);
        LJphi->DrawCopy("Esame");
        JVphi->SetLineColor(2);
        JVphi->SetMarkerStyle(26);
        JVphi->SetMarkerColor(2);
        JVphi->DrawCopy("Esame");
        JVCphi->SetLineColor(2);
        JVCphi->SetMarkerStyle(32);
        JVCphi->SetMarkerColor(2);
        JVCphi->DrawCopy("Esame");
        c3->cd();
        L3->Draw();

        c3->SaveAs("Sel8_qaJets.pdf");

}
