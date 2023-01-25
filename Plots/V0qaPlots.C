//Add the Lamber over kaon - normalize the ratio ! newPtb->Scale(1/Nb->Integral()) -> ok, not required anymore .. but why is the error so large ?;
void plotPerPt(TH2F *VarVsPt, int nBins, TString particle, TString varname){
        for(int i = 0; i < nBins; i++){
                //TH1F h = CosPaVsMAntiLambda->SetShowProjectionX(i);
                TCanvas *can = new TCanvas("can", "projection", 800, 400);
                TH1D* h = VarVsPt->ProjectionX("h",i);
                h->SetTitle(Form(particle+" "+varname+" in pT bin %d",i));
                h->Draw();
                can->SaveAs(Form("plots/perPtBin/"+particle+"/"+varname+"Pt_%d.pdf",i));
        }
}

//this is for now the main function
void V0qaPlots(){
        TFile *Result = new TFile("/home/johannalomker/alice/analysis/jetQuenching/AnalysisResults.root");

        TH2F *InvMvsPtK0Short = (TH2F*) Result->Get("correlationvzerojets/InvMvsPtK0Short");
        TH2F *InvMvsPtLambda = (TH2F*) Result->Get("correlationvzerojets/InvMvsPtLambda");
        TH2F *InvMvsPtAntiLambda = (TH2F*) Result->Get("correlationvzerojets/InvMvsPtAntiLambda");

        TH2F *CosPaVspTK0Short = (TH2F*) Result->Get("correlationvzerojets/CosPaVspTK0Short");
        TH2F *CosPaVspTLambda = (TH2F*) Result->Get("correlationvzerojets/CosPaVspTLambda");
        TH2F *CosPaVspTAntiLambda = (TH2F*) Result->Get("correlationvzerojets/CosPaVspTAntiLambda");

        TH2F *CosPaVsMK0Short = (TH2F*) Result->Get("correlationvzerojets/CosPaVsMK0Short");
        TH2F *CosPaVsMLambda = (TH2F*) Result->Get("correlationvzerojets/CosPaVsMLambda");
        TH2F *CosPaVsMAntiLambda = (TH2F*) Result->Get("correlationvzerojets/CosPaVsMAntiLambda");
        //one 2D disrivution for each plot
        TCanvas *V = new TCanvas("var vs pT or m", "var vs pT or m", 900, 1000);
        V->Divide(3,3);
        V->cd(1);
        InvMvsPtK0Short->GetXaxis()->SetRangeUser(0.2,0.8);
        InvMvsPtK0Short->Draw("COLZ");
        V->cd(2);
        InvMvsPtLambda->GetXaxis()->SetRangeUser(1,1.3);
        InvMvsPtLambda->Draw("COLZ");
        V->cd(3);
        InvMvsPtAntiLambda->GetXaxis()->SetRangeUser(1,1.3);
        InvMvsPtAntiLambda->Draw("COLZ");
        V->cd(4);
        CosPaVspTK0Short->Draw("COLZ");
        V->cd(5);
        CosPaVspTLambda->Draw("COLZ");
        V->cd(6);
        CosPaVspTAntiLambda->Draw("COLZ");
        V->cd(7);
        CosPaVsMK0Short->GetYaxis()->SetRangeUser(0.2,0.8);
        CosPaVsMK0Short->Draw("COLZ");
        V->cd(8);
        CosPaVsMLambda->GetYaxis()->SetRangeUser(1,1.3);
        CosPaVsMLambda->Draw("COLZ");
        V->cd(9);
        CosPaVsMAntiLambda->GetYaxis()->SetRangeUser(1,1.3);
        CosPaVsMAntiLambda->Draw("COLZ");
        V->SaveAs("plots/VarVsPt.pdf");  

        //loop over pT bins for distinct invariant mass and cospa distribution  
        plotPerPt(InvMvsPtK0Short, 100, "Kaon", "InvMass");
        plotPerPt(InvMvsPtLambda, 100, "Lambda", "InvMass");
        plotPerPt(InvMvsPtAntiLambda, 100, "AntiLambda", "InvMass");

        plotPerPt(CosPaVspTK0Short, 100, "Kaon", "CosPa");
        plotPerPt(CosPaVspTLambda, 100, "Lambda", "CosPa");
        plotPerPt(CosPaVspTAntiLambda, 100, "AntiLambda", "CosPa");

        //create differences and ratios for a selection criteria on cospa
}


void V0qaPlotsPtCuts(){
        TString input;
        //Output the Histos from correlationV0jet - the lambda over kaon, angular distance are not in here, but all qa's 
        //TFile *AResult = new TFile("/home/johannalomker/alice/analysis/jetQuenching/AnalysisResultsNoPtcut.root");
        TFile *Result1 = new TFile("/home/johannalomker/alice/analysis/jetQuenching/AnalysisResultsNoPtcut.root");
        TString input1 = "/home/johannalomker/alice/analysis/jetQuenching/AnalysisResultsNoPtcut.root";
        TFile *Result2 = new TFile("/home/johannalomker/alice/analysis/jetQuenching/AnalysisResultsNoDaughPtMin02.root");
        TString input2 = "/home/johannalomker/alice/analysis/jetQuenching/AnalysisResultsNoDaughPtMin02.root";
        TFile *Result3 = new TFile("/home/johannalomker/alice/analysis/jetQuenching/AnalysisResultsNoDaughPtMin04.root");
        TString input3 = "/home/johannalomker/alice/analysis/jetQuenching/AnalysisResultsNoDaughPtMin04.root";
        TFile *Result4 = new TFile("/home/johannalomker/alice/analysis/jetQuenching/AnalysisResultsNoDaughPtMin06.root");
        TString input4 = "/home/johannalomker/alice/analysis/jetQuenching/AnalysisResultsNoDaughPtMin06.root";
        //TFile *AResult = new TFile("/home/johannalomker/alice/analysis/jetQuenching/AnalysisV0CollisionId.root");
        //TFile *AResult = new TFile("/home/johannalomker/alice/analysis/jetQuenching/AnalysisV0GlobalIndex.root");

        TFile *AResult = new TFile(input1);; 

        //correlationvzerojets (no subdir)
        TH1F *VtxZ = (TH1F*) AResult->Get("correlationvzerojets/hCollVtxZ");

        TCanvas *can = new TCanvas("can", "vtx", 800, 400);
        can->Divide(2,1);
        can->cd(2);
        VtxZ->SetTitle("Vtx Z V0 process");
        VtxZ->Draw("E");
        can->SaveAs("plots/VtxJetsV0.pdf");

        TH1F *V0radius = (TH1F*) AResult->Get("correlationvzerojets/hV0radius");
        TH1F *cosPA = (TH1F*) AResult->Get("correlationvzerojets/hV0cospa");
        TCanvas *can3 = new TCanvas("can3", "V0 cut variables", 800, 400);
        can3->Divide(2,1);
        can3->cd(1);
        V0radius->GetYaxis()->SetTitle("Number of entries");
        V0radius->Draw("E");
        can3->cd(2);
        cosPA->GetYaxis()->SetTitle("Number of entries");
        cosPA->Draw("E");
        can3->SaveAs("plots/V0radiusAndCospa.pdf");

        //to check the dPhi and dEta calculation in jets!
        TH1F *J = (TH1F*) AResult->Get("correlationvzerojets/hdeltaEta");
        TH1F *M = (TH1F*) AResult->Get("correlationvzerojets/MdeltaPhi");
        TCanvas *can4 = new TCanvas("can4", "dPhi", 800, 400);
        can4->Divide(2,1);
        can4->cd(1);
        J->SetTitle("#Delta #eta");
        J->GetYaxis()->SetTitle("Number of entries");
        J->Draw("E");
        can4->cd(2);
        M->GetYaxis()->SetTitle("Number of entries");
        M->SetTitle("#Delta #Phi");
        M->Draw("E");
        can4->SaveAs("plots/DeltaPhiEta.pdf");

        gStyle -> SetOptStat(0);
        for(int i = 0; i<4; i++){
                TString result;
                if(i == 0){AResult = new TFile(input1);}
                if(i == 1){AResult = new TFile(input2);}
                if(i == 2){AResult = new TFile(input3);}
                if(i == 3){AResult = new TFile(input4);}
                
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
                TPad *t = new TPad("t", "t", 0., 0., 1, 1);
                t->SetBottomMargin(0.15); 
                t->SetLeftMargin(0.15);
	        t->SetRightMargin(0.1);
                t->Draw();          
                t->cd();               
                MK->SetTitle("");
                MK->GetYaxis()->SetTitleOffset(2);
                MK->GetYaxis()->SetTitle("Number of entries");
                MK->SetLineColor(3);
                MK->GetXaxis()->SetRangeUser(0,3);
                MK->Draw();
                MAL->SetLineColor(2);
                MAL->Draw("same");
                ML->SetLineColor(1);
                ML->Draw("same");
                c1->cd(2);
                TPad *t2 = new TPad("t2", "t2", 0., 0., 1, 1);
                t2->SetBottomMargin(0.15); 
                t2->SetLeftMargin(0.15);
        	t2->SetRightMargin(0.1);
                t2->SetLogy();
                t2->Draw();          
                t2->cd();  
                KPt->GetYaxis()->SetTitleOffset(2);
                KPt->GetYaxis()->SetTitle("Number of entries");
                KPt->SetTitle("");
                KPt->GetXaxis()->SetRangeUser(0,10);
                KPt->SetLineColor(3);
                KPt->Draw();
                LPt->SetLineColor(1);
                LPt->Draw("same");
                ALPt->SetLineColor(2);
                ALPt->Draw("same");
                c1->cd(3); 
                TPad *t3 = new TPad("t3", "t3", 0., 0., 1, 1);
                t3->SetBottomMargin(0.15); 
                t3->SetLeftMargin(0.15);
        	t3->SetRightMargin(0.1);
                t3->Draw();          
                t3->cd(); 
                KPhi->GetYaxis()->SetTitleOffset(2);
                KPhi->GetYaxis()->SetTitle("Number of entries");
                //KPhi->GetYaxis()->SetRangeUser(0, 2500);
                KPhi->SetTitle("");
                KPhi->SetLineColor(3);
                KPhi->Draw();
                c1->SetLogy();
                LPhi->SetLineColor(1);
                LPhi->Draw("same");
                ALPhi->SetLineColor(2);
                ALPhi->Draw("same");
                c1->cd(4);
                TPad *t4 = new TPad("t4", "t4", 0., 0., 1, 1);
                t4->SetBottomMargin(0.15); 
                t4->SetLeftMargin(0.15);
        	t4->SetRightMargin(0.1);
                t4->Draw();          
                t4->cd();  
                KEta->GetYaxis()->SetTitleOffset(2);
                KEta->GetYaxis()->SetTitle("Number of entries");
                KEta->SetTitle("");
                KEta->SetLineColor(3);
                KEta->Draw();
                LEta->SetLineColor(1);
                LEta->Draw("same");
                ALEta->SetLineColor(2);
                ALEta->Draw("same");
                c1->cd();
                L1->Draw(" ");
                c1->SaveAs(Form("plots/V0Candidates_input%i.pdf",i));

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

                int N = 2;//rebin variable
                int n = 2;//rebin phi
                TCanvas *c4 = new TCanvas("c4", "V0 daughter tracks", 1200, 1200);//phi and eta need all rebin by N, with N as free choice here
                c4->Divide(1,3);
                c4->cd(1);
                TPad *pad4 = new TPad("pad4", "pad4", 0., 0., 0.3, 0.9);
                pad4->SetBottomMargin(0.15); 
                pad4->SetLeftMargin(0.15);
	        pad4->SetRightMargin(0.1);
                pad4->SetLogy();
                pad4->Draw();             // Draw the upper pad: pad1
                pad4->cd();               // pad1 becomes the current pad
	        KPiPospt->GetXaxis()->SetRangeUser(0,15);
                KPiPospt->Sumw2();
                //KPiPospt->Rebin(n);
                //PiPospt->GetYaxis()->SetRangeUser(0,120000);
                KPiPospt->GetYaxis()->SetTitle("Number of entries");
	        KPiPospt->SetLineColor(3);
                KPiPospt->SetMarkerColor(3);
                KPiPospt->SetMarkerStyle(23);
                KPiPospt->SetTitle(" ");
                KPiPospt->DrawCopy();
                KPiNegpt->Sumw2();
                //KPiNegpt->Rebin(n);
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
                KPiPoseta->Sumw2();
                KPiPoseta->Rebin(N);
                KPiPoseta->GetYaxis()->SetTitle("Number of entries");            
	        KPiPoseta->SetLineColor(3);
                KPiPoseta->SetMarkerStyle(23);
                KPiPoseta->SetMarkerColor(3);
                KPiPoseta->SetTitle(" ");
                KPiPoseta->DrawCopy();
                KPiNegeta->Sumw2();
                KPiNegeta->Rebin(N);
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
                KPiPosphi->Sumw2();
                KPiPosphi->Rebin(n);
                KPiPosphi->GetYaxis()->SetTitle("Number of entries");            
                KPiPosphi->SetLineColor(3);
                KPiPosphi->SetMarkerStyle(23);
                KPiPosphi->SetMarkerColor(3);
                KPiPosphi->SetTitle(" ");
                KPiPosphi->DrawCopy();
                KPiNegphi->Sumw2();
                KPiNegphi->Rebin(n);
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
                pad41->SetLogy();
                pad41->Draw();             
                pad41->cd();   
                LPrPospt->Sumw2();
                //LPrPospt->Rebin(n);
                LPrPospt->GetYaxis()->SetTitle("Number of entries");            
                LPrPospt->GetXaxis()->SetRangeUser(0,15);
	        LPrPospt->SetLineColor(1);
                LPrPospt->SetMarkerColor(1);
                LPrPospt->SetMarkerStyle(23);
                LPrPospt->SetTitle(" ");
                LPrPospt->DrawCopy();
                LPiNegpt->Sumw2();
                //LPiNegpt->Rebin(n);
                LPiNegpt->SetLineColor(13);
                LPiNegpt->SetMarkerStyle(26);
                LPiNegpt->SetMarkerColor(13);
                LPiNegpt->DrawCopy("Esame");
                c4->cd(2);
                TPad *p41 = new TPad("pa41", "pa41", 0.3, 0., 0.6, 0.9);
                p41->SetBottomMargin(0.15); 
                p41->SetLeftMargin(0.15);
        	p41->SetRightMargin(0.1);
                p41->SetLogy();
                p41->Draw();            
                p41->cd();     
                LPrPoseta->Sumw2();
                LPrPoseta->Rebin(N);
                LPrPoseta->GetYaxis()->SetTitle("Number of entries");         
	        LPrPoseta->SetLineColor(1);
                LPrPoseta->SetMarkerColor(1);
                LPrPoseta->SetMarkerStyle(23);
                LPrPoseta->SetTitle(" ");
                LPrPoseta->DrawCopy();
                LPiNegeta->Sumw2();
                LPiNegeta->Rebin(N);
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
                LPrPosphi->Sumw2();
                LPrPosphi->Rebin(n);
                LPrPosphi->GetYaxis()->SetTitle("Number of entries");             
	        LPrPosphi->SetLineColor(1);
                LPrPosphi->SetMarkerColor(1);
                LPrPosphi->SetMarkerStyle(23);
                LPrPosphi->SetTitle(" ");
                LPrPosphi->DrawCopy();
                LPiNegphi->Sumw2();
                LPiNegphi->Rebin(n);
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
                pad42->SetLogy();
                pad42->Draw();             
                pad42->cd();   
                ALPiPospt->Sumw2();
                //ALPiPospt->Rebin(n);
                ALPiPospt->GetYaxis()->SetTitle("Number of entries");            
                ALPiPospt->GetXaxis()->SetRangeUser(0,15);
	        ALPiPospt->SetLineColor(2);
                ALPiPospt->SetMarkerColor(2);
                ALPiPospt->SetMarkerStyle(23);
                ALPiPospt->SetTitle(" ");
                ALPiPospt->DrawCopy();
                ALPrNegpt->Sumw2();
                //ALPrNegpt->Rebin(n);
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
                ALPiPoseta->Sumw2();
                ALPiPoseta->Rebin(N);
                ALPiPoseta->GetYaxis()->SetTitle("Number of entries");            
	        ALPiPoseta->SetLineColor(2);
                ALPiPoseta->SetMarkerColor(2);
                ALPiPoseta->SetMarkerStyle(23);
                ALPiPoseta->SetTitle(" ");
                ALPiPoseta->DrawCopy();
                ALPrNegeta->Sumw2();
                ALPrNegeta->Rebin(N);
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
                ALPiPosphi->Sumw2();
                ALPiPosphi->Rebin(n);
                ALPiPosphi->GetYaxis()->SetTitle("Number of entries");          
	        ALPiPosphi->SetLineColor(2);
                ALPiPosphi->SetMarkerColor(2);
                ALPiPosphi->SetMarkerStyle(23);
                ALPiPosphi->SetTitle(" ");
                ALPiPosphi->DrawCopy();
                ALPrNegphi->Sumw2();
                ALPrNegphi->Rebin(n);
                ALPrNegphi->SetLineColor(46);
                ALPrNegphi->SetMarkerStyle(26);
                ALPrNegphi->SetMarkerColor(46);
                ALPrNegphi->DrawCopy("Esame");
                c4->cd(3);
                L42->Draw();
                //c4->SaveAs("plots/qaV0Daughters.pdf");
                c4->SaveAs(Form("plots/qaV0Daughters_input%d.pdf",i));
        }//For loop over test files

        //Not relevant for group update

        TH1F *VRpt = (TH1F*) AResult->Get("correlationvzerojets/hPtTrackV0inRadius");
        TH1F *Vpt = (TH1F*) AResult->Get("correlationvzerojets/hPtV0");
        TH1F *Tpt = (TH1F*) AResult->Get("correlationvzerojets/hTrackPt");

        TH1F *VRe = (TH1F*) AResult->Get("correlationvzerojets/hEtaTrackV0inRadius");
        TH1F *Ve = (TH1F*) AResult->Get("correlationvzerojets/hEtaV0");
        TH1F *Te = (TH1F*) AResult->Get("correlationvzerojets/hTrackEta");

        TH1F *VRphi = (TH1F*) AResult->Get("correlationvzerojets/hPhiTrackV0inRadius");
        TH1F *Vphi = (TH1F*) AResult->Get("correlationvzerojets/hPhiV0");
        TH1F *Tphi = (TH1F*) AResult->Get("correlationvzerojets/hTrackPhi");

        auto L2 = new TLegend(0.2,0.49, 0.8,0.51);
        L2->SetHeader("","C");
        L2->SetNColumns(4);
        L2->SetTextSize(0.026);
        L2->AddEntry(VRpt, "hPtTrackV0inRadius", "lep");//from V0 process
        L2->AddEntry(Vpt, "hPtV0", "lep");//from V0 process
        L2->AddEntry(Tpt, "hTrackPt", "lep");//from V0 process
        L2->SetBorderSize(0);
        L2->SetFillStyle(0);
        
        TCanvas *c2 = new TCanvas("c2", "V0 process", 1200, 800); 
        c2->Divide(3,2);
        c2->cd(1);
        c2->SetTopMargin(0.1);
        Vpt->GetYaxis()->SetTitle("Number of entries");
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
        Ve->GetYaxis()->SetTitle("Number of entries");
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
        Vphi->GetYaxis()->SetTitle("Number of entries");
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
        Tpt->GetYaxis()->SetTitle("Number of entries");
        Tpt->SetLineColor(2);
        Tpt->SetMarkerStyle(24);
        Tpt->SetMarkerColor(2);
        Tpt->Draw("Esame");
        c2->cd(5);
        Te->GetYaxis()->SetTitle("Number of entries");
        Te->SetLineColor(2);
        Te->SetMarkerStyle(24);
        Te->SetMarkerColor(2);
        Te->Draw("Esame");
        c2->cd(6); 
        Tphi->GetYaxis()->SetTitle("Number of entries");
        Tphi->SetTitle("");
        Tphi->SetLineColor(2);
        Tphi->SetMarkerStyle(24);
        Tphi->SetMarkerColor(2);
        Tphi->Draw("Esame");
        c2->cd();
        L2->Draw();
        c2->SaveAs("plots/V0QA.pdf");

}
