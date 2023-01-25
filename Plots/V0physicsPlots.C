void V0physicsPlots(){
        TString input;
        //Output the Histos from correlationV0jet - the lambda over kaon, angular distance and vtx collision are not in here, but all qa's 
        TFile *Result1 = new TFile("/home/johannalomker/alice/analysis/jetQuenching/AnalysisResultsNoPtcut.root");
        TFile *Result2 = new TFile("/home/johannalomker/alice/analysis/jetQuenching/AnalysisResultsNoDaughPtMin02.root");
        TFile *Result3 = new TFile("/home/johannalomker/alice/analysis/jetQuenching/AnalysisResultsNoDaughPtMin04.root");
        TFile *Result4 = new TFile("/home/johannalomker/alice/analysis/jetQuenching/AnalysisResultsNoDaughPtMin06.root");

        gStyle -> SetOptStat(0);
        //correlationvzerojets (no subdir)
        TH1F *Lamb1 = (TH1F*) Result1->Get("correlationvzerojets/hPtLambda");
        TH1F *ALamb1 = (TH1F*) Result1->Get("correlationvzerojets/hPtAntiLambda");
        TH1F *Kaon1 = (TH1F*) Result1->Get("correlationvzerojets/hPtK0Short");
        TH1F *Lamb2 = (TH1F*) Result2->Get("correlationvzerojets/hPtLambda");
        TH1F *ALamb2 = (TH1F*) Result2->Get("correlationvzerojets/hPtAntiLambda");
        TH1F *Kaon2 = (TH1F*) Result2->Get("correlationvzerojets/hPtK0Short");
        TH1F *Lamb3 = (TH1F*) Result3->Get("correlationvzerojets/hPtLambda");
        TH1F *ALamb3 = (TH1F*) Result3->Get("correlationvzerojets/hPtAntiLambda");
        TH1F *Kaon3 = (TH1F*) Result3->Get("correlationvzerojets/hPtK0Short");
        TH1F *Lamb4 = (TH1F*) Result4->Get("correlationvzerojets/hPtLambda");
        TH1F *ALamb4 = (TH1F*) Result4->Get("correlationvzerojets/hPtAntiLambda");
        TH1F *Kaon4 = (TH1F*) Result4->Get("correlationvzerojets/hPtK0Short");
        
	Lamb1->Sumw2();
	Lamb2->Sumw2();
        Lamb3->Sumw2();
        Lamb4->Sumw2();

        TH1F *Sum1 = (TH1F*) Lamb1->Clone();
        TH1F *Sum2 = (TH1F*) Lamb2->Clone();
	TH1F *Sum3 = (TH1F*) Lamb3->Clone();
        TH1F *Sum4 = (TH1F*) Lamb4->Clone();

        Sum1->Add(ALamb1);
        Sum2->Add(ALamb2);
        Sum3->Add(ALamb3);
        Sum4->Add(ALamb4);
	
        Sum1->Divide(Sum1, Kaon1, 1.0, 2.0);
        Sum2->Divide(Sum2, Kaon2, 1.0, 2.0);
        Sum3->Divide(Sum3, Kaon3, 1.0, 2.0);
        Sum4->Divide(Sum4, Kaon4, 1.0, 2.0);

        auto L = new TLegend(0.2,0.65, 0.6,0.9);
        L->SetHeader("PbPb with minimum daughter p_{T}","C");
        //L->SetNColumns(4);
        L->SetTextSize(0.036);
        L->AddEntry(Sum1, "No cut ", "lep");
        L->AddEntry(Sum2 , " > 0.2 GeV/c ", "lep");
        L->AddEntry(Sum3 , " > 0.4 GeV/c ", "lep");
        L->AddEntry(Sum4 , " > 0.6 GeV/c ", "lep");
        L->SetBorderSize(0);
        L->SetFillStyle(0);

        Sum1->Rebin(2);
        Sum2->Rebin(2);
        Sum3->Rebin(2);
        Sum4->Rebin(2);

        TCanvas *can = new TCanvas("can", "ratio", 800, 500);
        TPad *pad = new TPad("pad", "pad", 0., 0.3, 1, 1);
        pad->SetBottomMargin(0.001); 
        pad->SetLeftMargin(0.15);
	pad->SetRightMargin(0.1);
        pad->Draw();             // Draw the upper pad: pad1
        pad->cd();               // pad1 becomes the current pad
        Sum1->GetXaxis()->SetRangeUser(0,12);//make more bins in low pT and then rebin properly !
        Sum1->GetYaxis()->SetRangeUser(0,1.1);
        Sum1->GetYaxis()->SetTitle("#frac{(#Lambda + #hat{#Lambda})}{2 K_{S}^{0}}");
        Sum1->SetTitle(" ");
        Sum1->SetLineColor(kRed+1);
        Sum1->DrawCopy("E");
        Sum2->SetLineColor(kBlue+1);
        //Sum2->Draw("Esame");
        Sum3->SetLineColor(kGreen+1);
        //Sum3->Draw("Esame");
        Sum4->SetLineColor(kBlack);
        //Sum4->Draw("Esame");
        L->Draw();
        can->cd();

        TH1F *div2 = (TH1F*) Sum2->Clone();
        TH1F *div3 = (TH1F*) Sum3->Clone();
        TH1F *div4 = (TH1F*) Sum4->Clone();

        div2->Sumw2();
        div3->Sumw2();
        div4->Sumw2();
        
        div2->Add(Sum1,-1);
        div3->Add(Sum1,-1);
        div4->Add(Sum1,-1);
        
        TPad *pad2 = new TPad("pad", "pad", 0., 0., 1, 0.3);
        pad2->SetBottomMargin(0.25); 
        pad2->SetLeftMargin(0.15);
	pad2->SetRightMargin(0.1);
        pad2->Draw();             // Draw the upper pad: pad1
        pad2->cd();               // pad1 becomes the current pad
        div2->SetTitle("");
        div2->GetYaxis()->SetTitle("Cut - NoCut");
        div2->GetYaxis()->SetTitleSize(0.099);
        div2->GetYaxis()->SetTitleOffset(0.4);
        div2->GetYaxis()->SetLabelSize(0.099);
        div2->GetXaxis()->SetTitleSize(0.099);
        div2->GetXaxis()->SetLabelSize(0.099);
        div2->GetXaxis()->SetRangeUser(0,12);//make more bins in low pT and then rebin properly !
        div2->GetYaxis()->SetRangeUser(-0.8,0.8);
        div2->GetYaxis()->SetNdivisions(4);
        div2->Draw("E");
        //div3->Draw("Esame");
        //div4->Draw("Esame");
        can->SaveAs("plots/Group_BaryonOverMeson.pdf");

}
