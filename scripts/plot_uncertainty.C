//#include "plot.h"

void plot_uncertainty(){
  //loadStyle();
    gROOT->SetBatch(0);
    gStyle->SetOptStat(0000);
    gStyle->SetOptFit(1111);
    gStyle->SetOptTitle(0);
    //gStyle->SetPadTickX(1);
    //gStyle->SetPadTickY(1);
    gStyle->SetPadColor(kWhite);
    gStyle->SetStatY(0.90);
    gStyle->SetStatX(0.90);
    gStyle->SetStatW(0.15);
    gStyle->SetStatH(0.2);
    gStyle->SetLabelSize(0.04,"X");
    gStyle->SetLabelFont(62,"X");
    gStyle->SetTitleSize(0.05,"X");
    gStyle->SetTitleFont(62,"X");
    gStyle->SetTitleOffset(0.8,"X");

    gStyle->SetLabelSize(0.04,"Y");
    gStyle->SetLabelFont(62,"Y");
    gStyle->SetTitleSize(0.06,"Y");
    gStyle->SetTitleFont(62,"Y");
    gStyle->SetTitleOffset(0.7,"Y");
    gStyle->SetTitleX(0.22);
    gStyle->SetTitleY(0.98);
    gStyle->SetTitleSize(0.04,"t");
    //gStyle->SetTitleTextColor(kRed);
    gStyle->SetTitleBorderSize(0);
    gStyle->SetCanvasBorderSize(0);
    //gStyle->SetTitleFontSize(0);
    //gStyle->SetCanvasColor(kWhite);
    //gStyle->SetPadGridX(kTRUE);
    //gStyle->SetPadGridY(kTRUE);
    //gStyle->SetGridStyle();
    //gStyle->SetErrorX(0.0001);

    std::vector<std::string> syst_unc_type;
    syst_unc_type.push_back("GENIE");
    syst_unc_type.push_back("Flux");
    syst_unc_type.push_back("Reinteraction");
    syst_unc_type.push_back("CCMEC+CCQE");
    syst_unc_type.push_back("Detector_Syst");
    
    TFile * inputfile[5];
    TH1D* Total_Muon_Momentum_9999[5];
    TH1D* Total_Proton_Momentum_9999[5];
    TH1D* Total_Muon_CosTheta_9999[5];
    TH1D* Total_Proton_CosTheta_9999[5];
    TH1D* Total_Thetamup_9999[5];

    TH1D* total_mumom_1;
    TH1D* total_pmom_1;
    TH1D* total_muangle_1;
    TH1D* total_pangle_1;
    TH1D* total_thetamup_1;

    double bins_mumom[7] = {0.00, 0.18, 0.30, 0.48, 0.75, 1.14, 2.50};
    double bins_mucostheta[13] = {-1.00, -0.82, -0.66, -0.39, -0.16, 0.05, 0.25, 0.43, 0.59, 0.73, 0.83, 0.91, 1.00};
    int n_bins_mumom = 6;
    int n_bins_mucostheta = 12;


    double bins_pmom[11] = {0.30, 0.41, 0.495, 0.56, 0.62, 0.68, 0.74, 0.80, 0.87, 0.93, 1.50};
    double bins_pcostheta[10] = {-1.00, -0.50, 0.00, 0.27, 0.45, 0.62, 0.76, 0.86, 0.94, 1.00};

    int n_bins_pmom = 10; 
    int n_bins_pcostheta = 9;

    double bins_muptheta[7] = {0.00, 1.00, 1.31, 1.5, 1.7, 2.04, 3.14};
    int n_bins_muptheta = 6;

    int dd=8888; 
    total_mumom_1=new TH1D(Form("Total_Muon_Momentum_%d",dd), Form("Total_Muon_Momentum_%d",dd), n_bins_mumom, bins_mumom);
    total_mumom_1->GetXaxis()->SetTitle("P_{#mu}[GeV]");
    total_mumom_1->GetYaxis()->SetTitle("Relative Uncertainty");
    total_mumom_1->SetLineColor(7);
    total_mumom_1->SetLineWidth(4);
    total_mumom_1->SetFillStyle(0);
    total_mumom_1->SetMinimum(0); 
    total_mumom_1->SetMaximum(1.0);
 
    total_pmom_1=new TH1D(Form("Total_Proton_Momentum_%d",dd), Form("Total_Proon_Momentum_%d",dd), n_bins_pmom, bins_pmom);
    total_pmom_1->GetXaxis()->SetTitle("P_{proton}[GeV]");
    total_pmom_1->GetYaxis()->SetTitle("Relative Uncertainty");
    total_pmom_1->SetLineColor(7);
    total_pmom_1->SetLineWidth(4);
    total_pmom_1->SetFillStyle(0);
    total_pmom_1->SetMinimum(0); 
    total_pmom_1->SetMaximum(1.0);
 
    total_muangle_1=new TH1D(Form("Total_Muon_CosTheta_%d",dd), Form("Total_Muon_CosTheta_%d",dd), n_bins_mucostheta, bins_mucostheta);
    total_muangle_1->GetXaxis()->SetTitle("Cos#theta_{#mu}");
    total_muangle_1->GetYaxis()->SetTitle("Relative Uncertainty");
    total_muangle_1->SetLineColor(7);
    total_muangle_1->SetLineWidth(4);
    total_muangle_1->SetFillStyle(0);
    total_muangle_1->SetMinimum(0); 
    total_muangle_1->SetMaximum(1.0);
 
    total_pangle_1=new TH1D(Form("Total_Proton_CosTheta_%d",dd), Form("Total_Proton_CosTheta_%d",dd), n_bins_pcostheta, bins_pcostheta);
    total_pangle_1->GetXaxis()->SetTitle("Cos#theta_{proton}");
    total_pangle_1->GetYaxis()->SetTitle("Relative Uncertainty");
    total_pangle_1->SetLineColor(7);
    total_pangle_1->SetLineWidth(4);
    total_pangle_1->SetFillStyle(0);
    total_pangle_1->SetMinimum(0); 
    total_pangle_1->SetMaximum(1.0);
 
    total_thetamup_1=new TH1D(Form("Total_Thetamup_%d",dd), Form("Total_Thetamup_%d",dd), n_bins_muptheta, bins_muptheta); 
    total_thetamup_1->GetXaxis()->SetTitle("#theta_{#mu, proton}");
    total_thetamup_1->GetYaxis()->SetTitle("Relative Uncertainty");
    total_thetamup_1->SetLineColor(7);
    total_thetamup_1->SetLineWidth(4);
    total_thetamup_1->SetFillStyle(0);
    total_thetamup_1->SetMinimum(0); 
    total_thetamup_1->SetMaximum(1.0);

    TLegend * leg=new TLegend(0.6, 0.7, 0.9, 0.9); 
    for(int i=0; i<syst_unc_type.size(); i++){
       string inputfilename=("output_"+syst_unc_type[i]+".root");
       inputfile[i]=new TFile(inputfilename.c_str());
       Total_Muon_Momentum_9999[i]=(TH1D*)inputfile[i]->Get("Total_Muon_Momentum_9999");
       Total_Proton_Momentum_9999[i]=(TH1D*)inputfile[i]->Get("Total_Proton_Momentum_9999");
       Total_Muon_CosTheta_9999[i]=(TH1D*)inputfile[i]->Get("Total_Muon_CosTheta_9999");
       Total_Proton_CosTheta_9999[i]=(TH1D*)inputfile[i]->Get("Total_Proton_CosTheta_9999");
       Total_Thetamup_9999[i]=(TH1D*)inputfile[i]->Get("Total_Thetamup_9999");

       int nbins;
       nbins=Total_Muon_Momentum_9999[i]->GetNbinsX();
       for(int j=0; j<nbins; j++){
         total_mumom_1->SetBinContent(j+1, total_mumom_1->GetBinContent(j+1) + Total_Muon_Momentum_9999[i]->GetBinContent(j+1)*Total_Muon_Momentum_9999[i]->GetBinContent(j+1)); 
       }
       nbins=Total_Proton_Momentum_9999[i]->GetNbinsX();
       for(int j=0; j<nbins; j++){
         total_pmom_1->SetBinContent(j+1, total_pmom_1->GetBinContent(j+1) + Total_Proton_Momentum_9999[i]->GetBinContent(j+1)*Total_Proton_Momentum_9999[i]->GetBinContent(j+1)); 
       }
       nbins=Total_Muon_CosTheta_9999[i]->GetNbinsX();
       for(int j=0; j<nbins; j++){
         total_muangle_1->SetBinContent(j+1, total_muangle_1->GetBinContent(j+1) + Total_Muon_CosTheta_9999[i]->GetBinContent(j+1)*Total_Muon_CosTheta_9999[i]->GetBinContent(j+1)); 
       }
       nbins=Total_Proton_CosTheta_9999[i]->GetNbinsX();
       for(int j=0; j<nbins; j++){
         total_pangle_1->SetBinContent(j+1, total_pangle_1->GetBinContent(j+1) + Total_Proton_CosTheta_9999[i]->GetBinContent(j+1)*Total_Proton_CosTheta_9999[i]->GetBinContent(j+1)); 
       }
       nbins=Total_Thetamup_9999[i]->GetNbinsX();
       for(int j=0; j<nbins; j++){
         total_thetamup_1->SetBinContent(j+1, total_thetamup_1->GetBinContent(j+1) + Total_Proton_CosTheta_9999[i]->GetBinContent(j+1)*Total_Proton_CosTheta_9999[i]->GetBinContent(j+1)); 
       }




       Total_Muon_Momentum_9999[i]->SetMaximum(1.0*Total_Muon_Momentum_9999[i]->GetMaximum());      
       Total_Muon_Momentum_9999[i]->SetLineColor(i+2);
       if(i+2==5) Total_Muon_Momentum_9999[i]->SetLineColor(28); 
 
       Total_Proton_Momentum_9999[i]->SetMaximum(1.0*Total_Proton_Momentum_9999[i]->GetMaximum());      
       Total_Proton_Momentum_9999[i]->SetLineColor(i+2);
       if(i+2==5) Total_Proton_Momentum_9999[i]->SetLineColor(28); 
 
       Total_Muon_CosTheta_9999[i]->SetMaximum(1.0*Total_Muon_CosTheta_9999[i]->GetMaximum());      
       Total_Muon_CosTheta_9999[i]->SetLineColor(i+2);
       if(i+2==5) Total_Muon_CosTheta_9999[i]->SetLineColor(28); 
 
 
       Total_Proton_CosTheta_9999[i]->SetMaximum(1.0*Total_Proton_CosTheta_9999[i]->GetMaximum());      
       Total_Proton_CosTheta_9999[i]->SetLineColor(i+2);
       if(i+2==5) Total_Proton_CosTheta_9999[i]->SetLineColor(28); 
 
       Total_Thetamup_9999[i]->SetMaximum(1.0*Total_Thetamup_9999[i]->GetMaximum());      
       Total_Thetamup_9999[i]->SetLineColor(i+2);
       if(i+2==5) Total_Thetamup_9999[i]->SetLineColor(28); 

       leg->AddEntry(Total_Muon_Momentum_9999[i], syst_unc_type[i].c_str(), "l");
    }

    int k=0;

    TCanvas * c_mumom=new TCanvas("c_mumom", "c_mumom");

    for(k=0; k<total_mumom_1->GetNbinsX(); k++){
        total_mumom_1->SetBinContent(k+1, sqrt(total_mumom_1->GetBinContent(k+1)));
    }
    int ind=0; 
    for(ind=0; ind<syst_unc_type.size(); ind++){
        if(ind==0){
           Total_Muon_Momentum_9999[ind]->Draw();
        }else{
           Total_Muon_Momentum_9999[ind]->Draw("same");
        }
          
    }
    leg->AddEntry(total_mumom_1, "Total Unc.");
    total_mumom_1->Draw("same");
    leg->Draw("same");
    c_mumom->SaveAs("Relative_unc_mumom.png");
    //=====================================================================================
    TCanvas * c_pmom=new TCanvas("c_pmom", "c_pmom");

    for(k=0; k<total_pmom_1->GetNbinsX(); k++){
        total_pmom_1->SetBinContent(k+1, sqrt(total_pmom_1->GetBinContent(k+1)));
    }
    for(ind=0; ind<syst_unc_type.size(); ind++){
        if(ind==0){
           Total_Proton_Momentum_9999[ind]->Draw();
        }else{
           Total_Proton_Momentum_9999[ind]->Draw("same");
        }
          
    }
    total_pmom_1->Draw("same");
    leg->Draw("same");
    c_pmom->SaveAs("Relative_unc_pmom.png");
    //=================================================================================
    TCanvas * c_muangle=new TCanvas("c_muangle", "c_muangle");

    for(k=0; k<total_muangle_1->GetNbinsX(); k++){
        total_muangle_1->SetBinContent(k+1, sqrt(total_muangle_1->GetBinContent(k+1)));
    }
    for(ind=0; ind<syst_unc_type.size(); ind++){
        if(ind==0){
           Total_Muon_CosTheta_9999[ind]->Draw();
        }else{
           Total_Muon_CosTheta_9999[ind]->Draw("same");
        }
          
    }
    total_muangle_1->Draw("same");
    leg->Draw("same");
    c_muangle->SaveAs("Relative_unc_muangle.png");
    //===================================================================================== 
    TCanvas * c_pangle=new TCanvas("c_pangle", "c_pangle");

    for(k=0; k<total_pangle_1->GetNbinsX(); k++){
        total_pangle_1->SetBinContent(k+1, sqrt(total_pangle_1->GetBinContent(k+1)));
    }
    for(ind=0; ind<syst_unc_type.size(); ind++){
        if(ind==0){
           Total_Proton_CosTheta_9999[ind]->Draw();
        }else{
           Total_Proton_CosTheta_9999[ind]->Draw("same");
        }
          
    }
    total_pangle_1->Draw("same");
    leg->Draw("same");
    c_pangle->SaveAs("Relative_unc_pangle.png");
    //=============================================================================================
    TCanvas * c_thetamup=new TCanvas("c_thetamup", "c_thetamup");

    for(k=0; k<total_thetamup_1->GetNbinsX(); k++){
        total_thetamup_1->SetBinContent(k+1, sqrt(total_thetamup_1->GetBinContent(k+1)));
    }
    for(ind=0; ind<syst_unc_type.size(); ind++){
        if(ind==0){
           Total_Thetamup_9999[ind]->Draw();
        }else{
           Total_Thetamup_9999[ind]->Draw("same");
        }
          
    }
    total_thetamup_1->Draw("same");
    leg->Draw("same");
    c_thetamup->SaveAs("Relative_unc_thetamup.png");
    //================================================================================================== 
  //open the xsec root files and get the central value of the cross sections
  TFile *xsec_inputfile= new TFile("xsec_file_CV.root");
  TH1D * xsec_mumom_mc;
  TH1D * xsec_mumom_CV;
  TH1D * xsec_pmom_mc;
  TH1D * xsec_pmom_CV;
  TH1D * xsec_muangle_mc;
  TH1D * xsec_muangle_CV;
  TH1D * xsec_pangle_mc;
  TH1D * xsec_pangle_CV;
  TH1D * xsec_thetamup_mc;
  TH1D * xsec_thetamup_CV;
  
  std::cout<<"libo test 0"<<std::endl; 
  xsec_mumom_mc=(TH1D*)xsec_inputfile->Get("xsec_mumom_mc_CV");
  xsec_mumom_CV=(TH1D*)xsec_inputfile->Get("xsec_mumom_CV");

  xsec_pmom_mc=(TH1D*)xsec_inputfile->Get("xsec_pmom_mc_CV");
  xsec_pmom_CV=(TH1D*)xsec_inputfile->Get("xsec_pmom_CV");

  xsec_muangle_mc=(TH1D*)xsec_inputfile->Get("xsec_muangle_mc_CV");
  xsec_muangle_CV=(TH1D*)xsec_inputfile->Get("xsec_muangle_CV");

  xsec_pangle_mc=(TH1D*)xsec_inputfile->Get("xsec_pangle_mc_CV");
  xsec_pangle_CV=(TH1D*)xsec_inputfile->Get("xsec_pangle_CV");

  xsec_thetamup_mc=(TH1D*)xsec_inputfile->Get("xsec_thetamup_mc_CV");
  xsec_thetamup_CV=(TH1D*)xsec_inputfile->Get("xsec_thetamup_CV");

  TLegend *legnew=new TLegend(0.65, 0.7, 0.9, 0.9);
  legnew->AddEntry(xsec_mumom_mc, "GENIE-default");
  legnew->AddEntry(xsec_mumom_CV, "Measured(Stat+Syst Unc.");
 
  //==============================================================================   
   TCanvas *c1=new TCanvas("c1", "c1");
   xsec_mumom_mc->SetMaximum(0.5*xsec_mumom_mc->GetMaximum());
   xsec_mumom_CV->SetMaximum(0.5*xsec_mumom_mc->GetMaximum()); 
   xsec_mumom_mc->SetTitleOffset(0.8, "X");
   xsec_mumom_CV->SetTitleOffset(0.8, "X");
   for(int jj=0; jj<xsec_mumom_CV->GetNbinsX(); jj++){
         
         xsec_mumom_CV->SetBinError(jj+1, xsec_mumom_CV->GetBinContent(jj+1)*TMath::Sqrt(TMath::Power((xsec_mumom_CV->GetBinError(jj+1)/xsec_mumom_CV->GetBinContent(jj+1)), 2.0)+TMath::Power(total_mumom_1->GetBinContent(jj+1),2.0)));
   }


   xsec_mumom_mc->Draw();
   xsec_mumom_CV->Draw("E1, same");
   legnew->Draw("same");

   c1->SaveAs("xsec_mumom.png");


   TCanvas *c2=new TCanvas("c2", "c2");
   xsec_pmom_mc->SetMaximum(0.8*xsec_pmom_mc->GetMaximum());
   xsec_pmom_CV->SetMaximum(0.8*xsec_pmom_mc->GetMaximum()); 
   xsec_pmom_mc->SetTitleOffset(0.8, "X");
   xsec_pmom_CV->SetTitleOffset(0.8, "X");
   for(int jj=0; jj<xsec_pmom_CV->GetNbinsX(); jj++){
         
         xsec_pmom_CV->SetBinError(jj+1, xsec_pmom_CV->GetBinContent(jj+1)*TMath::Sqrt(TMath::Power((xsec_pmom_CV->GetBinError(jj+1)/xsec_pmom_CV->GetBinContent(jj+1)), 2.0)+TMath::Power(total_pmom_1->GetBinContent(jj+1),2.0)));
   }


   xsec_pmom_mc->Draw();
   xsec_pmom_CV->Draw("E1, same");
   legnew->Draw("same");

   c2->SaveAs("xsec_pmom.png");

   TCanvas *c3=new TCanvas("c3", "c3");
   xsec_muangle_mc->SetMaximum(0.5*xsec_muangle_mc->GetMaximum());
   xsec_muangle_CV->SetMaximum(0.5*xsec_muangle_mc->GetMaximum()); 
   xsec_muangle_mc->SetTitleOffset(0.8, "X");
   xsec_muangle_CV->SetTitleOffset(0.8, "X");
   for(int jj=0; jj<xsec_muangle_CV->GetNbinsX(); jj++){
         
         xsec_muangle_CV->SetBinError(jj+1, xsec_muangle_CV->GetBinContent(jj+1)*TMath::Sqrt(TMath::Power((xsec_muangle_CV->GetBinError(jj+1)/xsec_muangle_CV->GetBinContent(jj+1)), 2.0)+TMath::Power(total_muangle_1->GetBinContent(jj+1),2.0)));
   }


   xsec_muangle_mc->Draw();
   xsec_muangle_CV->Draw("E1, same");
   legnew->Draw("same");
   
   c3->SaveAs("xsec_muangle.png");

   TCanvas *c4=new TCanvas("c4", "c4");
   xsec_pangle_mc->SetMaximum(0.5*xsec_pangle_mc->GetMaximum());
   xsec_pangle_CV->SetMaximum(0.5*xsec_pangle_mc->GetMaximum()); 
   xsec_pangle_mc->SetTitleOffset(0.8, "X");
   xsec_pangle_CV->SetTitleOffset(0.8, "X");
   for(int jj=0; jj<xsec_pangle_CV->GetNbinsX(); jj++){
         
         xsec_pangle_CV->SetBinError(jj+1, xsec_pangle_CV->GetBinContent(jj+1)*TMath::Sqrt(TMath::Power((xsec_pangle_CV->GetBinError(jj+1)/xsec_pangle_CV->GetBinContent(jj+1)), 2.0)+TMath::Power(total_pangle_1->GetBinContent(jj+1),2.0)));
   }


   xsec_pangle_mc->Draw();
   xsec_pangle_CV->Draw("E1, same");
   legnew->Draw("same");

   c4->SaveAs("xsec_pangle.png");

   TCanvas *c5=new TCanvas("c5", "c5");
   xsec_thetamup_mc->SetMaximum(0.2*xsec_thetamup_mc->GetMaximum());
   xsec_thetamup_CV->SetMaximum(0.2*xsec_thetamup_mc->GetMaximum()); 
   xsec_thetamup_mc->SetTitleOffset(0.8, "X");
   xsec_thetamup_CV->SetTitleOffset(0.8, "X");
   for(int jj=0; jj<xsec_thetamup_CV->GetNbinsX(); jj++){
         
         xsec_thetamup_CV->SetBinError(jj+1, xsec_thetamup_CV->GetBinContent(jj+1)*TMath::Sqrt(TMath::Power((xsec_thetamup_CV->GetBinError(jj+1)/xsec_thetamup_CV->GetBinContent(jj+1)), 2.0)+TMath::Power(total_thetamup_1->GetBinContent(jj+1),2.0)));
   }


   xsec_thetamup_mc->Draw();
   xsec_thetamup_CV->Draw("E1, same");
   legnew->Draw("same");
   c5->SaveAs("xsec_thetamup.png");















   //========================================================================================================



}
