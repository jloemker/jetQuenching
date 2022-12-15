void physicsPlots(){
        TString input;
        //Output the Histos from correlationV0jet - the lambda over kaon, angular distance and vtx collision are not in here, but all qa's 
        TFile *AResult = new TFile("/home/johannalomker/alice/analysis/jetQuenching/AnalysisV0CollisionId.root");
        //TFile *AResult = new TFile("/home/johannalomker/alice/analysis/jetQuenching/AnalysisV0GlobalIndex.root");

        //correlationvzerojets (no subdir)
        TH1F *Lamb = (TH1F*) AResult->Get("correlationvzerojets/hPtLambda");
        TH1F *ALamb = (TH1F*) AResult->Get("correlationvzerojets/hPtAntiLambda");
        TH1F *Kaon = (TH1F*) AResult->Get("correlationvzerojets/hPtK0Short");
        TH1F *BaryonMeson = (TH1F*) AResult->Get("correlationvzerojets/LambdaOverKaonPt");

        TH1F *Sum = (TH1F*) Lamb->Clone();
        Sum->Add(ALamb);
        Sum->Divide(Sum, Kaon, 1.0, 2.0);

        for(int i = 0; i<Lamb->GetNbinsX(); i++){
          double lamb = Lamb->GetBinContent(i);
          double errLamb = Lamb->GetBinError(i);
          double alamb = ALamb->GetBinContent(i);
          double errALamb = ALamb->GetBinError(i);
          double kaon = Kaon->GetBinContent(i);
          double errKaon = Kaon->GetBinError(i);
          if(kaon != 0 && errKaon != 0){
            double ratio = (lamb+alamb)/(2*kaon);
            double err =  pow( pow( (errLamb/(2*kaon)) ,2) + pow( (errALamb/(2*kaon)) ,2) + pow( (-(lamb+alamb)*errKaon/ (2*pow(kaon,2))) , 2 ), 1/2 );
            BaryonMeson->SetBinContent(i, ratio);
            BaryonMeson->SetBinError(i, err);
          }
        }

        TCanvas *can = new TCanvas("can", "ratio", 800, 400);
        can->Divide(2,1);
        can->cd(1);
        BaryonMeson->SetTitle("Error Prop.");
        BaryonMeson->GetXaxis()->SetRangeUser(0,10);
        BaryonMeson->Draw("E");
        can->cd(2);
        Sum->GetXaxis()->SetRangeUser(0,10);
        Sum->SetTitle("ROOT add/divide");
        Sum->Draw("E");
        can->SaveAs("plots/BaryonOverMeson.pdf");

        TH1F *AngularDistance = (TH1F*) AResult->Get("correlationvzerojets/AngularDistance");
        TCanvas *can2 = new TCanvas("can2", "AngularDistance", 800, 400);  
        can2->cd(1);
        AngularDistance->Draw("E");
        can2->SaveAs("plots/AngularDistance.pdf");
}
