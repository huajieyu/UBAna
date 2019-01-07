void stackHists(THStack *stack, TH1D *histarray_sig[], TH1D *histarray_bac[], TH1D *histarray_data[], const double &normfac, const double &scale_onoffbeam, const double &scale_dirt_MC){
  histarray_data[0]->SetLineColor(kBlack);
  histarray_data[0]->SetLineWidth(2);
  histarray_data[0]->SetLineStyle(1);

  histarray_sig[0]-> SetFillColor(2); 
  histarray_sig[0]-> SetLineWidth(1); 
  histarray_sig[0]->Scale(normfac);
  histarray_bac[0]-> SetFillColor(kBlue+1);
  histarray_bac[0]-> SetLineWidth(1);
  histarray_bac[0]->Scale(normfac);

  histarray_bac[1]-> SetFillColor(kGreen+2);
  histarray_bac[1]-> SetLineWidth(1);
  histarray_bac[1]->Scale(normfac);
  histarray_bac[2]-> SetFillColor(kGray+1);
  histarray_bac[2]-> SetLineWidth(1);
  histarray_bac[2]->Scale(normfac);
  histarray_bac[3]-> SetFillColor(kOrange+1);
  histarray_bac[3]-> SetLineWidth(1);
  histarray_bac[3]->Scale(normfac);

  histarray_bac[4]-> SetFillColor(kMagenta);
  histarray_bac[4]-> SetLineWidth(1);
  histarray_bac[4]->Scale(normfac);
  histarray_bac[5]->SetFillColor(46);
  histarray_bac[5]->SetLineWidth(1);
  histarray_bac[5]->Scale(normfac);
 
  histarray_data[1]->SetFillStyle(3005);
  histarray_data[1]->SetFillColor(28);
  histarray_data[1]->SetLineWidth(1);
  histarray_data[1]->Scale(scale_onoffbeam);

  histarray_data[3]->SetFillStyle(3004);
  histarray_data[3]->SetFillColor(11);
  histarray_data[3]->SetLineWidth(1);
  histarray_data[3]->Scale(scale_dirt_MC*normfac);

  // Merge background histograms as needed
  //histarray_bac[1]->Add(histarray_bac[2]); // CC0p0pi add CC0pNpi
  //histarray_bac[1]->Add(histarray_bac[3]); // CC0p0pi add CCNpNpi
  // Stack histograms into stack
  stack->Add(histarray_sig[0]); // signal
  stack->Add(histarray_bac[0]); // 
  stack->Add(histarray_bac[1]); // 
  stack->Add(histarray_bac[2]); // 
  stack->Add(histarray_bac[3]); // 
  stack->Add(histarray_bac[4]); // 
  stack->Add(histarray_bac[5]);
  stack->Add(histarray_data[3]); // Dirt MC
  stack->Add(histarray_data[1]); // EXT data
}
float Chi2Calc(TH1D *histo_MC, TH1D *histo_bnb, TH1D *histo_extbnb, float scale_offbeam, float norm_MC){
   int nbins=histo_MC->GetNbinsX();
   float N_MC[nbins], N_BNB[nbins], N_EXTBNB[nbins];

   float fac1=histo_MC->Integral()*norm_MC+histo_extbnb->Integral()*scale_offbeam;
   float fac2=histo_bnb->Integral();
   

   float Chi2=0.0;
   for(int ii=0; ii<nbins-1; ii++){
      N_MC[ii]=histo_MC->GetBinContent(ii+1);
      N_BNB[ii]=histo_bnb->GetBinContent(ii+1);
      N_EXTBNB[ii]=histo_extbnb->GetBinContent(ii+1);
      N_BNB[ii]*=fac1/fac2;

      if(N_MC[ii]!=0){
      Chi2 +=pow((N_MC[ii]*norm_MC+N_EXTBNB[ii]*scale_offbeam-N_BNB[ii]),2)/N_MC[ii];
      }
   }
   //Chi2=Chi2/nbins;

   return Chi2;
}

void plot_com_bacsep(){
  //loadStyle();
  int tune=1;
  int cosmicCut=1;
  
  TFile *input0; // on-beam
  TFile *input1; // off-beam
  TFile *input2; // MC
  TFile *input3; // dirt
  std::cout<<"Setup input root files "<<std::endl;
  if (cosmicCut){
    input0 = new TFile("/uboone/data/users/jiangl/ubxsec_static/v06_26_01_22_Dec20/ubxsecana_output_data_onbeam_ubcodev06_26_01_22.root");
    input1 = new TFile("/uboone/data/users/jiangl/ubxsec_static/v06_26_01_22_Dec20/ubxsecana_output_data_offbeam_ubcodev06_26_01_22.root");
  }
  else{
    input0 = new TFile("/uboone/data/users/jiangl/ubxsec_static/v06_26_01_22_Dec20/ubxsecana_output_data_onbeam_ubcodev06_26_01_22.root");
    input1 = new TFile("/uboone/data/users/jiangl/ubxsec_static/v06_26_01_22_Dec20/ubxsecana_output_data_offbeam_ubcodev06_26_01_22.root");
  }
  
  if (tune==3){
    if (cosmicCut){input2= new TFile("/uboone/data/users/jiangl/ubxsec_static/v06_26_01_22_Dec20/ubxsecana_output_mc_bnbcosmic_ubcodev06_26_01_22.root");}
    else{input2= new TFile("/uboone/data/users/jiangl/ubxsec_static/v06_26_01_22_Dec20/ubxsecana_output_mc_bnbcosmic_ubcodev06_26_01_22.root");}
  }
  else{
    if (cosmicCut){input2= new TFile("/uboone/data/users/jiangl/ubxsec_static/v06_26_01_22_Dec20/ubxsecana_output_mc_bnbcosmic_ubcodev06_26_01_22.root");}
    else{input2= new TFile("/uboone/data/users/jiangl/ubxsec_static/v06_26_01_22_Dec20/ubxsecana_output_mc_bnbcosmic_ubcodev06_26_01_22.root");}
  }
  input3 = new TFile("/uboone/data/users/jiangl/ubxsec_static/v06_26_01_22_Dec20/ubxsecana_output_mc_bnbdirt_ubcodev06_26_01_22.root");

  gROOT->SetBatch();
/*  
*/
/*
*/

  //----------------------------------------------------------------------
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
    gStyle->SetTitleOffset(0.85,"X");

    gStyle->SetLabelSize(0.04,"Y");
    gStyle->SetLabelFont(62,"Y");
    gStyle->SetTitleSize(0.06,"Y");
    gStyle->SetTitleFont(62,"Y");
    gStyle->SetTitleOffset(1.0,"Y");
    gStyle->SetTitleX(0.22);
    gStyle->SetTitleY(0.98);
    gStyle->SetTitleSize(0.04,"t");
    //gStyle->SetTitleTextColor(kRed);
    gStyle->SetTitleBorderSize(0);
    gStyle->SetCanvasBorderSize(0);
    //gStyle->SetTitleFontSize(0);
    //gStyle->SetCanvasColor(kWhite);
    //gStyle->SetPadGridX(kTRUE);
    gStyle->SetPadGridY(kTRUE);
    //gStyle->SetGridStyle();
  //----------------------------------------------------------------------
/* 
*/   //==========================================================
  std::cout<<"Setup the global enviromment "<<std::endl;

  TH1D                  *h_range_allsel[4];

  h_range_allsel[0]=(TH1D*)input0->Get("h_trklen_total");
  //h_range_allsel[0]->Rebin(4); 
  h_range_allsel[0]->Sumw2();

  h_range_allsel[1]=(TH1D*)input1->Get("h_trklen_total");
  //h_range_allsel[1]->Rebin(4);
  h_range_allsel[1]->Sumw2();

  h_range_allsel[2]=(TH1D*)input2->Get("h_trklen_total");
  //h_range_allsel[2]->Rebin(4);
  h_range_allsel[2]->Sumw2();

  h_range_allsel[3]=(TH1D*)input3->Get("h_trklen_total");
  //h_range_allsel[2]->Rebin(4);
  h_range_allsel[3]->Sumw2();


  TH1D               *h_range_sig[2];
  TH1D               *h_range_bac[8];
  h_range_sig[0]=(TH1D*)input2->Get("h_trklen_signal");
  //h_range_sig[0]->Rebin(4);
  h_range_sig[0]->Sumw2();

  h_range_bac[0]=(TH1D*)input2->Get("h_trklen_cosmic");
  //h_range_bac[0]->Rebin(4);
  h_range_bac[0]->Sumw2();
 
  h_range_bac[1]=(TH1D*)input2->Get("h_trklen_outfv");
  //h_range_bac[1]->rebin(4);
  h_range_bac[1]->Sumw2();

  h_range_bac[2]=(TH1D*)input2->Get("h_trklen_nc");
  //h_range_bac[2]->rebin(4);
  h_range_bac[2]->Sumw2();

  h_range_bac[3]=(TH1D*)input2->Get("h_trklen_anumu");
  //h_range_bac[3]->Rebin(4);
  h_range_bac[3]->Sumw2();
 
  h_range_bac[4]=(TH1D*)input2->Get("h_trklen_nue");
  //h_range_bac[4]->Rebin(4);
  h_range_bac[4]->Sumw2();

  h_range_bac[5]=(TH1D*)input2->Get("h_trklen_ccother");
  //h_range_bac[5]->Rebin(4);
  h_range_bac[5]->Sumw2();
 
  //===============================================================
  cout<<"get the histogram of track mucand range<<<<<<<<<<<"<<endl; 
  //===================================================================
  TH1D                  *h_prange_allsel[4];

  h_prange_allsel[0]=(TH1D*)input0->Get("h_trkplen_total");
  //h_prange_allsel[0]->Rebin(4); 
  h_prange_allsel[0]->Sumw2();

  h_prange_allsel[1]=(TH1D*)input1->Get("h_trkplen_total");
  //h_prange_allsel[1]->Rebin(4);
  h_prange_allsel[1]->Sumw2();

  h_prange_allsel[2]=(TH1D*)input2->Get("h_trkplen_total");
  //h_prange_allsel[2]->Rebin(4);
  h_prange_allsel[2]->Sumw2();
  
  h_prange_allsel[3]=(TH1D*)input3->Get("h_trkplen_total");
  //h_prange_allsel[2]->Rebin(4);
  h_prange_allsel[3]->Sumw2();

  TH1D               *h_prange_sig[2];
  TH1D               *h_prange_bac[8];
  h_prange_sig[0]=(TH1D*)input2->Get("h_trkplen_signal");
  //h_prange_sig[0]->Rebin(4);
  h_prange_sig[0]->Sumw2();

  h_prange_bac[0]=(TH1D*)input2->Get("h_trkplen_cosmic");
  //h_prange_bac[0]->Rebin(4);
  h_prange_bac[0]->Sumw2();
 
  h_prange_bac[1]=(TH1D*)input2->Get("h_trkplen_outfv");
  //h_prange_bac[1]->rebin(4);
  h_prange_bac[1]->Sumw2();

  h_prange_bac[2]=(TH1D*)input2->Get("h_trkplen_nc");
  //h_prange_bac[2]->rebin(4);
  h_prange_bac[2]->Sumw2();

  h_prange_bac[3]=(TH1D*)input2->Get("h_trkplen_anumu");
  //h_prange_bac[3]->Rebin(4);
  h_prange_bac[3]->Sumw2();
 
  h_prange_bac[4]=(TH1D*)input2->Get("h_trkplen_nue");
  //h_prange_bac[4]->Rebin(4);
  h_prange_bac[4]->Sumw2();

  h_prange_bac[5]=(TH1D*)input2->Get("h_trkplen_ccother");
  //h_prange_bac[5]->Rebin(4);
  h_prange_bac[5]->Sumw2();
   //===============================================================
  cout<<"get the histogram of track pcand prange<<<<<<<<<<<"<<endl; 




  //=============================================================
  TH1D                  *h_phi_allsel[4];

  h_phi_allsel[0]=(TH1D*)input0->Get("h_trkphi_total");
  //h_phi_allsel[0]->Rebin(4); 
  h_phi_allsel[0]->Sumw2();
 

  h_phi_allsel[1]=(TH1D*)input1->Get("h_trkphi_total");
  //h_phi_allsel[1]->Rebin(4);
  h_phi_allsel[1]->Sumw2();

  h_phi_allsel[2]=(TH1D*)input2->Get("h_trkphi_total");
  //h_phi_allsel[2]->Rebin(4);
  h_phi_allsel[2]->Sumw2();

  h_phi_allsel[3]=(TH1D*)input3->Get("h_trkphi_total");
  //h_phi_allsel[2]->Rebin(4);
  h_phi_allsel[3]->Sumw2();
  //
  TH1D               *h_phi_sig[2];
  TH1D               *h_phi_bac[8];
  h_phi_sig[0]=(TH1D*)input2->Get("h_trkphi_signal");
  //h_phi_sig[0]->Rebin(4);
  //h_phi_sig[0]->Sumw2();

  h_phi_bac[0]=(TH1D*)input2->Get("h_trkphi_cosmic");
  //h_phi_bac[0]->Rebin(4);
  //h_phi_bac[0]->Sumw2();
 
  h_phi_bac[1]=(TH1D*)input2->Get("h_trkphi_outfv");
  //h_phi_bac[1]->Rebin(4);
  //h_phi_bac[1]->Sumw2();

  h_phi_bac[2]=(TH1D*)input2->Get("h_trkphi_nc");
  //h_phi_bac[2]->Rebin(4);
  //h_phi_bac[2]->Sumw2();


  h_phi_bac[3]=(TH1D*)input2->Get("h_trkphi_anumu");
  //h_phi_bac[3]->Rebin(4);
  //h_phi_bac[3]->Sumw2();
 
  h_phi_bac[4]=(TH1D*)input2->Get("h_trkphi_nue");
  //h_phi_bac[4]->Rebin(4);
  //h_phi_bac[4]->Sumw2();

  h_phi_bac[5]=(TH1D*)input2->Get("h_trkphi_ccother");
  //h_phi_bac[5]->Rebin(4);
  //h_phi_bac[5]->Sumw2();
 
   //==================================================================== 
  cout<<"get the histogram of track phi <<<<<<<<<<<"<<endl; 
  //========================================================================
  TH1D                  *h_pphi_allsel[4];

  h_pphi_allsel[0]=(TH1D*)input0->Get("h_trkpphi_total");
  //h_pphi_allsel[0]->Rebin(4); 
  h_pphi_allsel[0]->Sumw2();
 

  h_pphi_allsel[1]=(TH1D*)input1->Get("h_trkpphi_total");
  //h_pphi_allsel[1]->Rebin(4);
  h_pphi_allsel[1]->Sumw2();

  h_pphi_allsel[2]=(TH1D*)input2->Get("h_trkpphi_total");
  //h_pphi_allsel[2]->Rebin(4);
  h_pphi_allsel[2]->Sumw2();

  h_pphi_allsel[3]=(TH1D*)input3->Get("h_trkpphi_total");
  //h_pphi_allsel[2]->Rebin(4);
  h_pphi_allsel[3]->Sumw2();

  TH1D               *h_pphi_sig[2];
  TH1D               *h_pphi_bac[8];
  h_pphi_sig[0]=(TH1D*)input2->Get("h_trkpphi_signal");
  //h_pphi_sig[0]->Rebin(4);
  h_pphi_sig[0]->Sumw2();

  h_pphi_bac[0]=(TH1D*)input2->Get("h_trkpphi_cosmic");
  //h_pphi_bac[0]->Rebin(4);
  h_pphi_bac[0]->Sumw2();
 
  h_pphi_bac[1]=(TH1D*)input2->Get("h_trkpphi_outfv");
  //h_pphi_bac[1]->Rebin(4);
  h_pphi_bac[1]->Sumw2();

  h_pphi_bac[2]=(TH1D*)input2->Get("h_trkpphi_nc");
  //h_pphi_bac[2]->Rebin(4);
  h_pphi_bac[2]->Sumw2();


  h_pphi_bac[3]=(TH1D*)input2->Get("h_trkpphi_anumu");
  //h_pphi_bac[3]->Rebin(4);
  h_pphi_bac[3]->Sumw2();
 
  h_pphi_bac[4]=(TH1D*)input2->Get("h_trkpphi_nue");
  //h_pphi_bac[4]->Rebin(4);
  h_pphi_bac[4]->Sumw2();

  h_pphi_bac[5]=(TH1D*)input2->Get("h_trkpphi_ccother");
  //h_pphi_bac[5]->Rebin(4);
  h_pphi_bac[5]->Sumw2();
 
   
  //========================================================================
  TH1D                  *h_costheta_allsel[4];

  h_costheta_allsel[0]=(TH1D*)input0->Get("h_trktheta_classic_total");
  //h_costheta_allsel[0]->Rebin(5); 
  h_costheta_allsel[0]->Sumw2();
 

  h_costheta_allsel[1]=(TH1D*)input1->Get("h_trktheta_classic_total");
  //h_costheta_allsel[1]->Rebin(5);
  h_costheta_allsel[1]->Sumw2();

  h_costheta_allsel[2]=(TH1D*)input2->Get("h_trktheta_classic_total");
  //h_costheta_allsel[2]->Rebin(5);
  h_costheta_allsel[2]->Sumw2();

  h_costheta_allsel[3]=(TH1D*)input3->Get("h_trktheta_classic_total");
  //h_costheta_allsel[3]->Rebin(5);
  h_costheta_allsel[3]->Sumw2();

  TH1D               *h_costheta_sig[2];
  TH1D               *h_costheta_bac[8];
  h_costheta_sig[0]=(TH1D*)input2->Get("h_trktheta_classic_signal");
  //h_costheta_sig[0]->Rebin(5);
  //h_costheta_sig[0]->Sumw2();

  h_costheta_bac[0]=(TH1D*)input2->Get("h_trktheta_classic_cosmic");
  //h_costheta_bac[0]->Rebin(5);
  //h_costheta_bac[0]->Sumw2();
 
  h_costheta_bac[1]=(TH1D*)input2->Get("h_trktheta_classic_outfv");
  //h_costheta_bac[1]->Rebin(5);
  //h_costheta_bac[1]->Sumw2();

  h_costheta_bac[2]=(TH1D*)input2->Get("h_trktheta_classic_nc");
  //h_costheta_bac[2]->Rebin(5);
  //h_costheta_bac[2]->Sumw2();

  h_costheta_bac[3]=(TH1D*)input2->Get("h_trktheta_classic_anumu");
  //h_costheta_bac[3]->Rebin(5);
  //h_costheta_bac[3]->Sumw2();
 
  h_costheta_bac[4]=(TH1D*)input2->Get("h_trktheta_classic_nue");
  //h_costheta_bac[4]->Rebin(5);
  //h_costheta_bac[4]->Sumw2();

  h_costheta_bac[5]=(TH1D*)input2->Get("h_trktheta_classic_ccother");
  //h_costheta_bac[5]->Rebin(5);
  //h_costheta_bac[5]->Sumw2();

   //========================================================
  cout<<"get histograms of muon candidate angular distirbution<<<<<<"<<endl;
   //===============================================================
  TH1D                  *h_pcostheta_allsel[4];

  h_pcostheta_allsel[0]=(TH1D*)input0->Get("h_trkptheta_classic_total");
  //h_pcostheta_allsel[0]->Rebin(5); 
  h_pcostheta_allsel[0]->Sumw2();
 

  h_pcostheta_allsel[1]=(TH1D*)input1->Get("h_trkptheta_classic_total");
  //h_pcostheta_allsel[1]->Rebin(5);
  h_pcostheta_allsel[1]->Sumw2();

  h_pcostheta_allsel[2]=(TH1D*)input2->Get("h_trkptheta_classic_total");
  //h_pcostheta_allsel[2]->Rebin(5);
  h_pcostheta_allsel[2]->Sumw2();

  h_pcostheta_allsel[3]=(TH1D*)input3->Get("h_trkptheta_classic_total");
  //h_pcostheta_allsel[3]->Rebin(5);
  h_pcostheta_allsel[3]->Sumw2();

  TH1D               *h_pcostheta_sig[2];
  TH1D               *h_pcostheta_bac[8];
  h_pcostheta_sig[0]=(TH1D*)input2->Get("h_trkptheta_classic_signal");
  //h_pcostheta_sig[0]->Rebin(5);
  h_pcostheta_sig[0]->Sumw2();

  h_pcostheta_bac[0]=(TH1D*)input2->Get("h_trkptheta_classic_cosmic");
  //h_pcostheta_bac[0]->Rebin(5);
  h_pcostheta_bac[0]->Sumw2();
 
  h_pcostheta_bac[1]=(TH1D*)input2->Get("h_trkptheta_classic_outfv");
  //h_pcostheta_bac[1]->Rebin(5);
  h_pcostheta_bac[1]->Sumw2();

  h_pcostheta_bac[2]=(TH1D*)input2->Get("h_trkptheta_classic_nc");
  //h_pcostheta_bac[2]->Rebin(5);
  h_pcostheta_bac[2]->Sumw2();

  h_pcostheta_bac[3]=(TH1D*)input2->Get("h_trkptheta_classic_anumu");
  //h_pcostheta_bac[3]->Rebin(5);
  h_pcostheta_bac[3]->Sumw2();
 
  h_pcostheta_bac[4]=(TH1D*)input2->Get("h_trkptheta_classic_nue");
  //h_pcostheta_bac[4]->Rebin(5);
  h_pcostheta_bac[4]->Sumw2();

  h_pcostheta_bac[5]=(TH1D*)input2->Get("h_trkptheta_classic_ccother");
  //h_pcostheta_bac[5]->Rebin(5);
  h_pcostheta_bac[5]->Sumw2();

   //========================================================
 
  //========================================================================
  TH1D                  *h_plep_allsel[4];

  h_plep_allsel[0]=(TH1D*)input0->Get("h_trkmom_classic_total");
  //h_plep_allsel[0]->Rebin(5); 
  h_plep_allsel[0]->Sumw2();
 

  h_plep_allsel[1]=(TH1D*)input1->Get("h_trkmom_classic_total");
  //h_plep_allsel[1]->Rebin(5);
  h_plep_allsel[1]->Sumw2();

  h_plep_allsel[2]=(TH1D*)input2->Get("h_trkmom_classic_total");
  //h_plep_allsel[2]->Rebin(5);
  h_plep_allsel[2]->Sumw2();

  h_plep_allsel[3]=(TH1D*)input3->Get("h_trkmom_classic_total");
  //h_plep_allsel[3]->Rebin(5);
  h_plep_allsel[3]->Sumw2();

  TH1D               *h_plep_sig[2];
  TH1D               *h_plep_bac[8];
  h_plep_sig[0]=(TH1D*)input2->Get("h_trkmom_classic_signal");
  //h_plep_sig[0]->Rebin(5);
  h_plep_sig[0]->Sumw2();

  h_plep_bac[0]=(TH1D*)input2->Get("h_trkmom_classic_cosmic");
  //h_plep_bac[0]->Rebin(5);
  h_plep_bac[0]->Sumw2();
 
  h_plep_bac[1]=(TH1D*)input2->Get("h_trkmom_classic_outfv");
  //h_plep_bac[1]->Rebin(5);
  h_plep_bac[1]->Sumw2();

  h_plep_bac[2]=(TH1D*)input2->Get("h_trkmom_classic_nc");
  //h_plep_bac[2]->Rebin(5);
  h_plep_bac[2]->Sumw2();

  h_plep_bac[3]=(TH1D*)input2->Get("h_trkmom_classic_anumu");
  //h_plep_bac[3]->Rebin(5);
  h_plep_bac[3]->Sumw2();
 
  h_plep_bac[4]=(TH1D*)input2->Get("h_trkmom_classic_nue");
  //h_plep_bac[4]->Rebin(5);
  h_plep_bac[4]->Sumw2();

  h_plep_bac[5]=(TH1D*)input2->Get("h_trkmom_classic_ccother");
  //h_plep_bac[5]->Rebin(5);
  h_plep_bac[5]->Sumw2();
 
  std::cout<<"get histograms of the muon momentum "<<std::endl;
   //=============================================================== 
  TH1D                  *h_phad_allsel[4];

  h_phad_allsel[0]=(TH1D*)input0->Get("h_trkpmom_classic_total");
  h_phad_allsel[0]->Rebin(5); 
  //h_phad_allsel[0]->GetXaxis()->SetRange(1,15);
  //h_phad_allsel[0]->Sumw2();
 

  h_phad_allsel[1]=(TH1D*)input1->Get("h_trkpmom_classic_total");
  h_phad_allsel[1]->Rebin(5);
  //h_phad_allsel[1]->GetXaxis()->SetRange(1,15);
  //h_phad_allsel[1]->Sumw2();

  h_phad_allsel[2]=(TH1D*)input2->Get("h_trkpmom_classic_total");
  h_phad_allsel[2]->Rebin(5);
  //h_phad_allsel[2]->GetXaxis()->SetRange(1,15);
  //h_phad_allsel[2]->Sumw2();
  
  h_phad_allsel[3]=(TH1D*)input3->Get("h_trkpmom_classic_total");
  h_phad_allsel[3]->Rebin(5);
  //h_phad_allsel[3]->GetXaxis()->SetRange(1,15);
  //h_phad_allsel[3]->Sumw2();
  
  
  TH1D               *h_phad_sig[2];
  TH1D               *h_phad_bac[8];
  h_phad_sig[0]=(TH1D*)input2->Get("h_trkpmom_classic_signal");
  h_phad_sig[0]->Rebin(5);
  //h_phad_sig[0]->GetXaxis()->SetRange(1,15);
  //h_phad_sig[0]->Sumw2();

  h_phad_bac[0]=(TH1D*)input2->Get("h_trkpmom_classic_cosmic");
  h_phad_bac[0]->Rebin(5);
  //h_phad_bac[0]->GetXaxis()->SetRange(1,15);
  //h_phad_bac[0]->Sumw2();

  h_phad_bac[1]=(TH1D*)input2->Get("h_trkpmom_classic_outfv");
  h_phad_bac[1]->Rebin(5);
  //h_phad_bac[1]->GetXaxis()->SetRange(1,15);
  //h_phad_bac[1]->Sumw2();

  h_phad_bac[2]=(TH1D*)input2->Get("h_trkpmom_classic_nc");
  h_phad_bac[2]->Rebin(5);
  //h_phad_bac[2]->GetXaxis()->SetRange(1,15);
  //h_phad_bac[2]->Sumw2();

  h_phad_bac[3]=(TH1D*)input2->Get("h_trkpmom_classic_anumu");
  h_phad_bac[3]->Rebin(5);
  //h_phad_bac[3]->GetXaxis()->SetRange(1,15);
  //h_phad_bac[3]->Sumw2();
 
  h_phad_bac[4]=(TH1D*)input2->Get("h_trkpmom_classic_nue");
  h_phad_bac[4]->Rebin(5);
  //h_phad_bac[4]->GetXaxis()->SetRange(1,15);
  //h_phad_bac[4]->Sumw2();

  h_phad_bac[5]=(TH1D*)input2->Get("h_trkpmom_classic_ccother");
  h_phad_bac[5]->Rebin(5);
  //h_phad_bac[5]->GetXaxis()->SetRange(1,15);
  //h_phad_bac[5]->Sumw2();
 
   //=============================================================== 
  cout<<"get the histogram of proton momentum phad<<<<<<<<<<<<"<<endl;
  //==============================================================
 /* 
  */

  //======================================================================
  TH1D                  *h_thetamup_allsel[4];

  h_thetamup_allsel[0]=(TH1D*)input0->Get("h_thetamup_total");
  //h_thetamup_allsel[0]->Rebin(4); 
  h_thetamup_allsel[0]->Sumw2();
 

  h_thetamup_allsel[1]=(TH1D*)input1->Get("h_thetamup_total");
  //h_thetamup_allsel[1]->Rebin(4);
  h_thetamup_allsel[1]->Sumw2();

  h_thetamup_allsel[2]=(TH1D*)input2->Get("h_thetamup_total");
  //h_thetamup_allsel[2]->Rebin(4);
  h_thetamup_allsel[2]->Sumw2();

  h_thetamup_allsel[3]=(TH1D*)input3->Get("h_thetamup_total");
  //h_thetamup_allsel[3]->Rebin(4);
  h_thetamup_allsel[3]->Sumw2();

  TH1D               *h_thetamup_sig[2];
  TH1D               *h_thetamup_bac[8];
  h_thetamup_sig[0]=(TH1D*)input2->Get("h_thetamup_signal");
  //h_thetamup_sig[0]->Rebin(4);
  h_thetamup_sig[0]->Sumw2();

  h_thetamup_bac[0]=(TH1D*)input2->Get("h_thetamup_cosmic");
  //h_thetamup_bac[0]->Rebin(4);
  h_thetamup_bac[0]->Sumw2();
 
  h_thetamup_bac[1]=(TH1D*)input2->Get("h_thetamup_outfv");
  //h_thetamup_bac[1]->Rebin(4);
  h_thetamup_bac[1]->Sumw2();

  h_thetamup_bac[2]=(TH1D*)input2->Get("h_thetamup_nc");
  //h_thetamup_bac[2]->Rebin(4);
  h_thetamup_bac[2]->Sumw2();
 
  h_thetamup_bac[3]=(TH1D*)input2->Get("h_thetamup_anumu");
  //h_thetamup_bac[3]->Rebin(4);
  h_thetamup_bac[3]->Sumw2();

  h_thetamup_bac[4]=(TH1D*)input2->Get("h_thetamup_nue");
  //h_thetamup_bac[4]->Rebin(4);
  h_thetamup_bac[4]->Sumw2();
 
  h_thetamup_bac[5]=(TH1D*)input2->Get("h_thetamup_ccother");
  //h_thetamup_bac[5]->Rebin(4);
  h_thetamup_bac[5]->Sumw2();

  //====================================================

  cout<<"get the histogram of thetamup<<<<<<<<<<<<<<<"<<endl;
  //===================================================
  TH1D                  *h_ptmis_allsel[4];

  h_ptmis_allsel[0]=(TH1D*)input0->Get("h_ptmis_total");
  //h_ptmis_allsel[0]->Rebin(4); 
  h_ptmis_allsel[0]->Sumw2();
 

  h_ptmis_allsel[1]=(TH1D*)input1->Get("h_ptmis_total");
  //h_ptmis_allsel[1]->Rebin(4);
  h_ptmis_allsel[1]->Sumw2();

  h_ptmis_allsel[2]=(TH1D*)input2->Get("h_ptmis_total");
  //h_ptmis_allsel[2]->Rebin(4);
  h_ptmis_allsel[2]->Sumw2();
  
  h_ptmis_allsel[3]=(TH1D*)input3->Get("h_ptmis_total");
  //h_ptmis_allsel[3]->Rebin(4);
  h_ptmis_allsel[3]->Sumw2();

  TH1D               *h_ptmis_sig[2];
  TH1D               *h_ptmis_bac[8];
  h_ptmis_sig[0]=(TH1D*)input2->Get("h_ptmis_signal");
  //h_ptmis_sig[0]->Rebin(4);
  h_ptmis_sig[0]->Sumw2();

  h_ptmis_bac[0]=(TH1D*)input2->Get("h_ptmis_cosmic");
  //h_ptmis_bac[0]->Rebin(4);
  h_ptmis_bac[0]->Sumw2();
 
  h_ptmis_bac[1]=(TH1D*)input2->Get("h_ptmis_outfv");
  //h_ptmis_bac[1]->Rebin(4);
  h_ptmis_bac[1]->Sumw2();

  h_ptmis_bac[2]=(TH1D*)input2->Get("h_ptmis_nc");
  //h_ptmis_bac[2]->Rebin(4);
  h_ptmis_bac[2]->Sumw2();
 
  h_ptmis_bac[3]=(TH1D*)input2->Get("h_ptmis_anumu");
  //h_ptmis_bac[3]->Rebin(4);
  h_ptmis_bac[3]->Sumw2();

  h_ptmis_bac[4]=(TH1D*)input2->Get("h_ptmis_nue");
  //h_ptmis_bac[4]->Rebin(4);
  h_ptmis_bac[4]->Sumw2();
 
  h_ptmis_bac[5]=(TH1D*)input2->Get("h_ptmis_ccother");
  //h_ptmis_bac[5]->Rebin(4);
  h_ptmis_bac[5]->Sumw2();
  //===================================================
  cout<<"get the histogram of ptmis<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"<<endl;
  //====================================================
  TH1D                  *h_etatest_allsel[4];

  h_etatest_allsel[0]=(TH1D*)input0->Get("h_etatest_total");
  //h_etatest_allsel[0]->Rebin(4); 
  h_etatest_allsel[0]->Sumw2();
 

  h_etatest_allsel[1]=(TH1D*)input1->Get("h_etatest_total");
  //h_etatest_allsel[1]->Rebin(4);
  h_etatest_allsel[1]->Sumw2();

  h_etatest_allsel[2]=(TH1D*)input2->Get("h_etatest_total");
  //h_etatest_allsel[2]->Rebin(4);
  h_etatest_allsel[2]->Sumw2();

  h_etatest_allsel[3]=(TH1D*)input3->Get("h_etatest_total");
  //h_etatest_allsel[3]->Rebin(4);
  h_etatest_allsel[3]->Sumw2();

  TH1D               *h_etatest_sig[2];
  TH1D               *h_etatest_bac[8];
  h_etatest_sig[0]=(TH1D*)input2->Get("h_etatest_signal");
  //h_etatest_sig[0]->Rebin(4);
  h_etatest_sig[0]->Sumw2();

  h_etatest_bac[0]=(TH1D*)input2->Get("h_etatest_cosmic");
  //h_etatest_bac[0]->Rebin(4);
  h_etatest_bac[0]->Sumw2();
 
  h_etatest_bac[1]=(TH1D*)input2->Get("h_etatest_outfv");
  //h_etatest_bac[1]->Rebin(4);
  h_etatest_bac[1]->Sumw2();

  h_etatest_bac[2]=(TH1D*)input2->Get("h_etatest_nc");
  //h_etatest_bac[2]->Rebin(4);
  h_etatest_bac[2]->Sumw2();
 
  h_etatest_bac[3]=(TH1D*)input2->Get("h_etatest_anumu");
  //h_etatest_bac[3]->Rebin(4);
  h_etatest_bac[3]->Sumw2();

  h_etatest_bac[4]=(TH1D*)input2->Get("h_etatest_nue");
  //h_etatest_bac[4]->Rebin(4);
  h_etatest_bac[4]->Sumw2();
 
  h_etatest_bac[5]=(TH1D*)input2->Get("h_etatest_ccother");
  //h_etatest_bac[5]->Rebin(4);
  h_etatest_bac[5]->Sumw2();
  //===================================================
  cout<<"get the histogram of etatest<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"<<endl;
  //====================================================
 



  //======================================================================
  TH1D                  *h_enucal_allsel[4];

  h_enucal_allsel[0]=(TH1D*)input0->Get("h_enucal_total");
  //h_enucal_allsel[0]->Rebin(4); 
  h_enucal_allsel[0]->Sumw2();
 

  h_enucal_allsel[1]=(TH1D*)input1->Get("h_enucal_total");
  //h_enucal_allsel[1]->Rebin(4);
  h_enucal_allsel[1]->Sumw2();

  h_enucal_allsel[2]=(TH1D*)input2->Get("h_enucal_total");
  //h_enucal_allsel[2]->Rebin(4);
  h_enucal_allsel[2]->Sumw2();

  h_enucal_allsel[3]=(TH1D*)input3->Get("h_enucal_total");
  //h_enucal_allsel[3]->Rebin(4);
  h_enucal_allsel[3]->Sumw2();

  TH1D               *h_enucal_sig[2];
  TH1D               *h_enucal_bac[8];
  h_enucal_sig[0]=(TH1D*)input2->Get("h_enucal_signal");
  //h_enucal_sig[0]->Rebin(4);
  h_enucal_sig[0]->Sumw2();

  h_enucal_bac[0]=(TH1D*)input2->Get("h_enucal_cosmic");
  //h_enucal_bac[0]->Rebin(4);
  h_enucal_bac[0]->Sumw2();
 
  h_enucal_bac[1]=(TH1D*)input2->Get("h_enucal_outfv");
  //h_enucal_bac[1]->Rebin(4);
  h_enucal_bac[1]->Sumw2();

  h_enucal_bac[2]=(TH1D*)input2->Get("h_enucal_nc");
  //h_enucal_bac[2]->Rebin(4);
  h_enucal_bac[2]->Sumw2();
 
  h_enucal_bac[3]=(TH1D*)input2->Get("h_enucal_anumu");
  //h_enucal_bac[3]->Rebin(4);
  h_enucal_bac[3]->Sumw2();

  h_enucal_bac[4]=(TH1D*)input2->Get("h_enucal_nue");
  //h_enucal_bac[4]->Rebin(4);
  h_enucal_bac[4]->Sumw2();
 
  h_enucal_bac[5]=(TH1D*)input2->Get("h_enucal_ccother");
  //h_enucal_bac[5]->Rebin(4);
  h_enucal_bac[5]->Sumw2();
  //===================================================
  cout<<"get the histogram of enucal<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"<<endl;
  //==========================================================
  TH1D                  *h_alphat_allsel[4];

  h_alphat_allsel[0]=(TH1D*)input0->Get("h_alphat_total");
  //h_alphat_allsel[0]->Rebin(4); 
  h_alphat_allsel[0]->Sumw2();
 

  h_alphat_allsel[1]=(TH1D*)input1->Get("h_alphat_total");
  //h_alphat_allsel[1]->Rebin(4);
  h_alphat_allsel[1]->Sumw2();

  h_alphat_allsel[2]=(TH1D*)input2->Get("h_alphat_total");
  //h_alphat_allsel[2]->Rebin(4);
  h_alphat_allsel[2]->Sumw2();

  h_alphat_allsel[3]=(TH1D*)input3->Get("h_alphat_total");
  //h_alphat_allsel[3]->Rebin(4);
  h_alphat_allsel[3]->Sumw2();

  TH1D               *h_alphat_sig[2];
  TH1D               *h_alphat_bac[8];
  h_alphat_sig[0]=(TH1D*)input2->Get("h_alphat_signal");
  //h_alphat_sig[0]->Rebin(4);
  h_alphat_sig[0]->Sumw2();

  h_alphat_bac[0]=(TH1D*)input2->Get("h_alphat_cosmic");
  //h_alphat_bac[0]->Rebin(4);
  h_alphat_bac[0]->Sumw2();
 
  h_alphat_bac[1]=(TH1D*)input2->Get("h_alphat_outfv");
  //h_alphat_bac[1]->Rebin(4);
  h_alphat_bac[1]->Sumw2();

  h_alphat_bac[2]=(TH1D*)input2->Get("h_alphat_nc");
  //h_alphat_bac[2]->Rebin(4);
  h_alphat_bac[2]->Sumw2();
 
  h_alphat_bac[3]=(TH1D*)input2->Get("h_alphat_anumu");
  //h_alphat_bac[3]->Rebin(4);
  h_alphat_bac[3]->Sumw2();

  h_alphat_bac[4]=(TH1D*)input2->Get("h_alphat_nue");
  //h_alphat_bac[4]->Rebin(4);
  h_alphat_bac[4]->Sumw2();
 
  h_alphat_bac[5]=(TH1D*)input2->Get("h_alphat_ccother");
  //h_alphat_bac[5]->Rebin(4);
  h_alphat_bac[5]->Sumw2();
  //===================================================
  cout<<"get the histogram of alphat<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"<<endl;
  //======================================================
  TH1D                  *h_pmult_allsel[4];

  h_pmult_allsel[0]=(TH1D*)input0->Get("h_pmult_total");
  //h_pmult_allsel[0]->Rebin(4); 
  h_pmult_allsel[0]->Sumw2();
 

  h_pmult_allsel[1]=(TH1D*)input1->Get("h_pmult_total");
  //h_pmult_allsel[1]->Rebin(4);
  h_pmult_allsel[1]->Sumw2();

  h_pmult_allsel[2]=(TH1D*)input2->Get("h_pmult_total");
  //h_pmult_allsel[2]->Rebin(4);
  h_pmult_allsel[2]->Sumw2();
  
  h_pmult_allsel[3]=(TH1D*)input3->Get("h_pmult_total");
  //h_pmult_allsel[3]->Rebin(4);
  h_pmult_allsel[3]->Sumw2();

  TH1D               *h_pmult_sig[2];
  TH1D               *h_pmult_bac[8];
  h_pmult_sig[0]=(TH1D*)input2->Get("h_pmult_signal");
  //h_pmult_sig[0]->Rebin(4);
  h_pmult_sig[0]->Sumw2();

  h_pmult_bac[0]=(TH1D*)input2->Get("h_pmult_cosmic");
  //h_pmult_bac[0]->Rebin(4);
  h_pmult_bac[0]->Sumw2();
 
  h_pmult_bac[1]=(TH1D*)input2->Get("h_pmult_outfv");
  //h_pmult_bac[1]->Rebin(4);
  h_pmult_bac[1]->Sumw2();

  h_pmult_bac[2]=(TH1D*)input2->Get("h_pmult_nc");
  //h_pmult_bac[2]->Rebin(4);
  h_pmult_bac[2]->Sumw2();
 
  h_pmult_bac[3]=(TH1D*)input2->Get("h_pmult_anumu");
  //h_pmult_bac[3]->Rebin(4);
  h_pmult_bac[3]->Sumw2();

  h_pmult_bac[4]=(TH1D*)input2->Get("h_pmult_nue");
  //h_pmult_bac[4]->Rebin(4);
  h_pmult_bac[4]->Sumw2();
 
  h_pmult_bac[5]=(TH1D*)input2->Get("h_pmult_ccother");
  //h_pmult_bac[5]->Rebin(4);
  h_pmult_bac[5]->Sumw2();
  //===================================================
  cout<<"get the histogram of pmult<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"<<endl;
  
  /* 
 * TH1D                  *h_vtxx_allsel[3];

  h_vtxx_allsel[0]=(TH1D*)input0->Get("fvertex_x");
  h_vtxx_allsel[0]->Rebin(4); 
  h_vtxx_allsel[0]->Sumw2();
 

  h_vtxx_allsel[1]=(TH1D*)input1->Get("fvertex_x");
  h_vtxx_allsel[1]->Rebin(4);
  h_vtxx_allsel[1]->Sumw2();

  h_vtxx_allsel[2]=(TH1D*)input2->Get("fvertex_x");
  h_vtxx_allsel[2]->Rebin(4);
  h_vtxx_allsel[2]->Sumw2();

  TH1D               *h_vtxx_sig[2];
  TH1D               *h_vtxx_bac[8];
  h_vtxx_sig[0]=(TH1D*)input2->Get("fsigvertex_x_0");
  h_vtxx_sig[0]->Rebin(4);
  h_vtxx_sig[0]->Sumw2();

  h_vtxx_bac[0]=(TH1D*)input2->Get("fbacvertex_x_0");
  h_vtxx_bac[0]->Rebin(4);
  h_vtxx_bac[0]->Sumw2();
 
  h_vtxx_bac[1]=(TH1D*)input2->Get("fbacvertex_x_1");
  h_vtxx_bac[1]->Rebin(4);
  h_vtxx_bac[1]->Sumw2();

  h_vtxx_bac[2]=(TH1D*)input2->Get("fbacvertex_x_2");
  h_vtxx_bac[2]->Rebin(4);
  h_vtxx_bac[2]->Sumw2();
 
  h_vtxx_bac[3]=(TH1D*)input2->Get("fbacvertex_x_3");
  h_vtxx_bac[3]->Rebin(4);
  h_vtxx_bac[3]->Sumw2();

  h_vtxx_bac[4]=(TH1D*)input2->Get("fbacvertex_x_4");
  h_vtxx_bac[4]->Rebin(4);
  h_vtxx_bac[4]->Sumw2();
 
  h_vtxx_bac[5]=(TH1D*)input2->Get("fbacvertex_x_5");
  h_vtxx_bac[5]->Rebin(4);
  h_vtxx_bac[5]->Sumw2();

  h_vtxx_bac[6]=(TH1D*)input2->Get("fbacvertex_x_6");
  h_vtxx_bac[6]->Rebin(4);
  h_vtxx_bac[6]->Sumw2();

  h_vtxx_bac[7]=(TH1D*)input2->Get("fbacvertex_x_7");
  h_vtxx_bac[7]->Rebin(4);
  h_vtxx_bac[7]->Sumw2();
   //===============================================================
  TH1D                  *h_vtxy_allsel[3];

  h_vtxy_allsel[0]=(TH1D*)input0->Get("fvertex_y");
  h_vtxy_allsel[0]->Rebin(4); 
  h_vtxy_allsel[0]->Sumw2();
 

  h_vtxy_allsel[1]=(TH1D*)input1->Get("fvertex_y");
  h_vtxy_allsel[1]->Rebin(4);
  h_vtxy_allsel[1]->Sumw2();

  h_vtxy_allsel[2]=(TH1D*)input2->Get("fvertex_y");
  h_vtxy_allsel[2]->Rebin(4);
  h_vtxy_allsel[2]->Sumw2();

  TH1D               *h_vtxy_sig[2];
  TH1D               *h_vtxy_bac[8];
  h_vtxy_sig[0]=(TH1D*)input2->Get("fsigvertex_y_0");
  h_vtxy_sig[0]->Rebin(4);
  h_vtxy_sig[0]->Sumw2();

  h_vtxy_bac[0]=(TH1D*)input2->Get("fbacvertex_y_0");
  h_vtxy_bac[0]->Rebin(4);
  h_vtxy_bac[0]->Sumw2();
 
  h_vtxy_bac[1]=(TH1D*)input2->Get("fbacvertex_y_1");
  h_vtxy_bac[1]->Rebin(4);
  h_vtxy_bac[1]->Sumw2();

  h_vtxy_bac[2]=(TH1D*)input2->Get("fbacvertex_y_2");
  h_vtxy_bac[2]->Rebin(4);
  h_vtxy_bac[2]->Sumw2();
 
  h_vtxy_bac[3]=(TH1D*)input2->Get("fbacvertex_y_3");
  h_vtxy_bac[3]->Rebin(4);
  h_vtxy_bac[3]->Sumw2();

  h_vtxy_bac[4]=(TH1D*)input2->Get("fbacvertex_y_4");
  h_vtxy_bac[4]->Rebin(4);
  h_vtxy_bac[4]->Sumw2();
 
  h_vtxy_bac[5]=(TH1D*)input2->Get("fbacvertex_y_5");
  h_vtxy_bac[5]->Rebin(4);
  h_vtxy_bac[5]->Sumw2();

  h_vtxy_bac[6]=(TH1D*)input2->Get("fbacvertex_y_6");
  h_vtxy_bac[6]->Rebin(4);
  h_vtxy_bac[6]->Sumw2();

  h_vtxy_bac[7]=(TH1D*)input2->Get("fbacvertex_y_7");
  h_vtxy_bac[7]->Rebin(4);
  h_vtxy_bac[7]->Sumw2();
   //===============================================================
  TH1D                  *h_vtxz_allsel[3];

  h_vtxz_allsel[0]=(TH1D*)input0->Get("fvertex_z");
  h_vtxz_allsel[0]->Rebin(4); 
  h_vtxz_allsel[0]->Sumw2();
 

  h_vtxz_allsel[1]=(TH1D*)input1->Get("fvertex_z");
  h_vtxz_allsel[1]->Rebin(4);
  h_vtxz_allsel[1]->Sumw2();

  h_vtxz_allsel[2]=(TH1D*)input2->Get("fvertex_z");
  h_vtxz_allsel[2]->Rebin(4);
  h_vtxz_allsel[2]->Sumw2();

  TH1D               *h_vtxz_sig[2];
  TH1D               *h_vtxz_bac[8];
  h_vtxz_sig[0]=(TH1D*)input2->Get("fsigvertex_z_0");
  h_vtxz_sig[0]->Rebin(4);
  h_vtxz_sig[0]->Sumw2();

  h_vtxz_bac[0]=(TH1D*)input2->Get("fbacvertex_z_0");
  h_vtxz_bac[0]->Rebin(4);
  h_vtxz_bac[0]->Sumw2();
 
  h_vtxz_bac[1]=(TH1D*)input2->Get("fbacvertex_z_1");
  h_vtxz_bac[1]->Rebin(4);
  h_vtxz_bac[1]->Sumw2();

  h_vtxz_bac[2]=(TH1D*)input2->Get("fbacvertex_z_2");
  h_vtxz_bac[2]->Rebin(4);
  h_vtxz_bac[2]->Sumw2();
 
  h_vtxz_bac[3]=(TH1D*)input2->Get("fbacvertex_z_3");
  h_vtxz_bac[3]->Rebin(4);
  h_vtxz_bac[3]->Sumw2();

  h_vtxz_bac[4]=(TH1D*)input2->Get("fbacvertex_z_4");
  h_vtxz_bac[4]->Rebin(4);
  h_vtxz_bac[4]->Sumw2();
 
  h_vtxz_bac[5]=(TH1D*)input2->Get("fbacvertex_z_5");
  h_vtxz_bac[5]->Rebin(4);
  h_vtxz_bac[5]->Sumw2();

  h_vtxz_bac[6]=(TH1D*)input2->Get("fbacvertex_z_6");
  h_vtxz_bac[6]->Rebin(4);
  h_vtxz_bac[6]->Sumw2();
 
  h_vtxz_bac[7]=(TH1D*)input2->Get("fbacvertex_z_7");
  h_vtxz_bac[7]->Rebin(4);
  h_vtxz_bac[7]->Sumw2();
 
  //====================================================================== 
  TH1D                  *h_mustartx_allsel[3];

  h_mustartx_allsel[0]=(TH1D*)input0->Get("ftrkstartx_mucand");
  h_mustartx_allsel[0]->Rebin(4); 
  h_mustartx_allsel[0]->Sumw2();
 

  h_mustartx_allsel[1]=(TH1D*)input1->Get("ftrkstartx_mucand");
  h_mustartx_allsel[1]->Rebin(4);
  h_mustartx_allsel[1]->Sumw2();

  h_mustartx_allsel[2]=(TH1D*)input2->Get("ftrkstartx_mucand");
  h_mustartx_allsel[2]->Rebin(4);
  h_mustartx_allsel[2]->Sumw2();

  TH1D               *h_mustartx_sig[2];
  TH1D               *h_mustartx_bac[8];
  h_mustartx_sig[0]=(TH1D*)input2->Get("fsig_trkstart_mucand_x_0");
  h_mustartx_sig[0]->Rebin(4);
  h_mustartx_sig[0]->Sumw2();

  h_mustartx_bac[0]=(TH1D*)input2->Get("fbac_trkstart_mucand_x_0");
  h_mustartx_bac[0]->Rebin(4);
  h_mustartx_bac[0]->Sumw2();
 
  h_mustartx_bac[1]=(TH1D*)input2->Get("fbac_trkstart_mucand_x_1");
  h_mustartx_bac[1]->Rebin(4);
  h_mustartx_bac[1]->Sumw2();

  h_mustartx_bac[2]=(TH1D*)input2->Get("fbac_trkstart_mucand_x_2");
  h_mustartx_bac[2]->Rebin(4);
  h_mustartx_bac[2]->Sumw2();
 
  h_mustartx_bac[3]=(TH1D*)input2->Get("fbac_trkstart_mucand_x_3");
  h_mustartx_bac[3]->Rebin(4);
  h_mustartx_bac[3]->Sumw2();

  h_mustartx_bac[4]=(TH1D*)input2->Get("fbac_trkstart_mucand_x_4");
  h_mustartx_bac[4]->Rebin(4);
  h_mustartx_bac[4]->Sumw2();
 
  h_mustartx_bac[5]=(TH1D*)input2->Get("fbac_trkstart_mucand_x_5");
  h_mustartx_bac[5]->Rebin(4);
  h_mustartx_bac[5]->Sumw2();

  h_mustartx_bac[6]=(TH1D*)input2->Get("fbac_trkstart_mucand_x_6");
  h_mustartx_bac[6]->Rebin(4);
  h_mustartx_bac[6]->Sumw2();

  h_mustartx_bac[7]=(TH1D*)input2->Get("fbac_trkstart_mucand_x_7");
  h_mustartx_bac[7]->Rebin(4);
  h_mustartx_bac[7]->Sumw2();
   //=============================================================================== 
  TH1D                  *h_mustarty_allsel[3];

  h_mustarty_allsel[0]=(TH1D*)input0->Get("ftrkstarty_mucand");
  h_mustarty_allsel[0]->Rebin(4); 
  h_mustarty_allsel[0]->Sumw2();
 

  h_mustarty_allsel[1]=(TH1D*)input1->Get("ftrkstarty_mucand");
  h_mustarty_allsel[1]->Rebin(4);
  h_mustarty_allsel[1]->Sumw2();

  h_mustarty_allsel[2]=(TH1D*)input2->Get("ftrkstarty_mucand");
  h_mustarty_allsel[2]->Rebin(4);
  h_mustarty_allsel[2]->Sumw2();

  TH1D               *h_mustarty_sig[2];
  TH1D               *h_mustarty_bac[8];
  h_mustarty_sig[0]=(TH1D*)input2->Get("fsig_trkstart_mucand_y_0");
  h_mustarty_sig[0]->Rebin(4);
  h_mustarty_sig[0]->Sumw2();

  h_mustarty_bac[0]=(TH1D*)input2->Get("fbac_trkstart_mucand_y_0");
  h_mustarty_bac[0]->Rebin(4);
  h_mustarty_bac[0]->Sumw2();
 
  h_mustarty_bac[1]=(TH1D*)input2->Get("fbac_trkstart_mucand_y_1");
  h_mustarty_bac[1]->Rebin(4);
  h_mustarty_bac[1]->Sumw2();

  h_mustarty_bac[2]=(TH1D*)input2->Get("fbac_trkstart_mucand_y_2");
  h_mustarty_bac[2]->Rebin(4);
  h_mustarty_bac[2]->Sumw2();
 
  h_mustarty_bac[3]=(TH1D*)input2->Get("fbac_trkstart_mucand_y_3");
  h_mustarty_bac[3]->Rebin(4);
  h_mustarty_bac[3]->Sumw2();

  h_mustarty_bac[4]=(TH1D*)input2->Get("fbac_trkstart_mucand_y_4");
  h_mustarty_bac[4]->Rebin(4);
  h_mustarty_bac[4]->Sumw2();
 
  h_mustarty_bac[5]=(TH1D*)input2->Get("fbac_trkstart_mucand_y_5");
  h_mustarty_bac[5]->Rebin(4);
  h_mustarty_bac[5]->Sumw2();

  h_mustarty_bac[6]=(TH1D*)input2->Get("fbac_trkstart_mucand_y_6");
  h_mustarty_bac[6]->Rebin(4);
  h_mustarty_bac[6]->Sumw2();

  h_mustarty_bac[7]=(TH1D*)input2->Get("fbac_trkstart_mucand_y_7");
  h_mustarty_bac[7]->Rebin(4);
  h_mustarty_bac[7]->Sumw2();
   //===================================================================
  TH1D                  *h_mustartz_allsel[3];

  h_mustartz_allsel[0]=(TH1D*)input0->Get("ftrkstartz_mucand");
  h_mustartz_allsel[0]->Rebin(4); 
  h_mustartz_allsel[0]->Sumw2();
 

  h_mustartz_allsel[1]=(TH1D*)input1->Get("ftrkstartz_mucand");
  h_mustartz_allsel[1]->Rebin(4);
  h_mustartz_allsel[1]->Sumw2();

  h_mustartz_allsel[2]=(TH1D*)input2->Get("ftrkstartz_mucand");
  h_mustartz_allsel[2]->Rebin(4);
  h_mustartz_allsel[2]->Sumw2();

  TH1D               *h_mustartz_sig[2];
  TH1D               *h_mustartz_bac[8];
  h_mustartz_sig[0]=(TH1D*)input2->Get("fsig_trkstart_mucand_z_0");
  h_mustartz_sig[0]->Rebin(4);
  h_mustartz_sig[0]->Sumw2();

  h_mustartz_bac[0]=(TH1D*)input2->Get("fbac_trkstart_mucand_z_0");
  h_mustartz_bac[0]->Rebin(4);
  h_mustartz_bac[0]->Sumw2();
 
  h_mustartz_bac[1]=(TH1D*)input2->Get("fbac_trkstart_mucand_z_1");
  h_mustartz_bac[1]->Rebin(4);
  h_mustartz_bac[1]->Sumw2();

  h_mustartz_bac[2]=(TH1D*)input2->Get("fbac_trkstart_mucand_z_2");
  h_mustartz_bac[2]->Rebin(4);
  h_mustartz_bac[2]->Sumw2();
 
  h_mustartz_bac[3]=(TH1D*)input2->Get("fbac_trkstart_mucand_z_3");
  h_mustartz_bac[3]->Rebin(4);
  h_mustartz_bac[3]->Sumw2();

  h_mustartz_bac[4]=(TH1D*)input2->Get("fbac_trkstart_mucand_z_4");
  h_mustartz_bac[4]->Rebin(4);
  h_mustartz_bac[4]->Sumw2();
 
  h_mustartz_bac[5]=(TH1D*)input2->Get("fbac_trkstart_mucand_z_5");
  h_mustartz_bac[5]->Rebin(4);
  h_mustartz_bac[5]->Sumw2();

  h_mustartz_bac[6]=(TH1D*)input2->Get("fbac_trkstart_mucand_z_6");
  h_mustartz_bac[6]->Rebin(4);
  h_mustartz_bac[6]->Sumw2();

  h_mustartz_bac[7]=(TH1D*)input2->Get("fbac_trkstart_mucand_z_7");
  h_mustartz_bac[7]->Rebin(4);
  h_mustartz_bac[7]->Sumw2();


  //==========================================================================
  TH1D                  *h_muendx_allsel[3];

  h_muendx_allsel[0]=(TH1D*)input0->Get("ftrkendx_mucand");
  h_muendx_allsel[0]->Rebin(4); 
  h_muendx_allsel[0]->Sumw2();
 

  h_muendx_allsel[1]=(TH1D*)input1->Get("ftrkendx_mucand");
  h_muendx_allsel[1]->Rebin(4);
  h_muendx_allsel[1]->Sumw2();

  h_muendx_allsel[2]=(TH1D*)input2->Get("ftrkendx_mucand");
  h_muendx_allsel[2]->Rebin(4);
  h_muendx_allsel[2]->Sumw2();

  TH1D               *h_muendx_sig[2];
  TH1D               *h_muendx_bac[8];
  h_muendx_sig[0]=(TH1D*)input2->Get("fsig_trkend_mucand_x_0");
  h_muendx_sig[0]->Rebin(4);
  h_muendx_sig[0]->Sumw2();

  h_muendx_bac[0]=(TH1D*)input2->Get("fbac_trkend_mucand_x_0");
  h_muendx_bac[0]->Rebin(4);
  h_muendx_bac[0]->Sumw2();
 
  h_muendx_bac[1]=(TH1D*)input2->Get("fbac_trkend_mucand_x_1");
  h_muendx_bac[1]->Rebin(4);
  h_muendx_bac[1]->Sumw2();

  h_muendx_bac[2]=(TH1D*)input2->Get("fbac_trkend_mucand_x_2");
  h_muendx_bac[2]->Rebin(4);
  h_muendx_bac[2]->Sumw2();
 
  h_muendx_bac[3]=(TH1D*)input2->Get("fbac_trkend_mucand_x_3");
  h_muendx_bac[3]->Rebin(4);
  h_muendx_bac[3]->Sumw2();

  h_muendx_bac[4]=(TH1D*)input2->Get("fbac_trkend_mucand_x_4");
  h_muendx_bac[4]->Rebin(4);
  h_muendx_bac[4]->Sumw2();
 
  h_muendx_bac[5]=(TH1D*)input2->Get("fbac_trkend_mucand_x_5");
  h_muendx_bac[5]->Rebin(4);
  h_muendx_bac[5]->Sumw2();

  h_muendx_bac[6]=(TH1D*)input2->Get("fbac_trkend_mucand_x_6");
  h_muendx_bac[6]->Rebin(4);
  h_muendx_bac[6]->Sumw2();

  h_muendx_bac[7]=(TH1D*)input2->Get("fbac_trkend_mucand_x_7");
  h_muendx_bac[7]->Rebin(4);
  h_muendx_bac[7]->Sumw2();
   //=============================================================================== 
  TH1D                  *h_muendy_allsel[3];

  h_muendy_allsel[0]=(TH1D*)input0->Get("ftrkendy_mucand");
  h_muendy_allsel[0]->Rebin(4); 
  h_muendy_allsel[0]->Sumw2();
 

  h_muendy_allsel[1]=(TH1D*)input1->Get("ftrkendy_mucand");
  h_muendy_allsel[1]->Rebin(4);
  h_muendy_allsel[1]->Sumw2();

  h_muendy_allsel[2]=(TH1D*)input2->Get("ftrkendy_mucand");
  h_muendy_allsel[2]->Rebin(4);
  h_muendy_allsel[2]->Sumw2();

  TH1D               *h_muendy_sig[2];
  TH1D               *h_muendy_bac[8];
  h_muendy_sig[0]=(TH1D*)input2->Get("fsig_trkend_mucand_y_0");
  h_muendy_sig[0]->Rebin(4);
  h_muendy_sig[0]->Sumw2();

  h_muendy_bac[0]=(TH1D*)input2->Get("fbac_trkend_mucand_y_0");
  h_muendy_bac[0]->Rebin(4);
  h_muendy_bac[0]->Sumw2();
 
  h_muendy_bac[1]=(TH1D*)input2->Get("fbac_trkend_mucand_y_1");
  h_muendy_bac[1]->Rebin(4);
  h_muendy_bac[1]->Sumw2();

  h_muendy_bac[2]=(TH1D*)input2->Get("fbac_trkend_mucand_y_2");
  h_muendy_bac[2]->Rebin(4);
  h_muendy_bac[2]->Sumw2();
 
  h_muendy_bac[3]=(TH1D*)input2->Get("fbac_trkend_mucand_y_3");
  h_muendy_bac[3]->Rebin(4);
  h_muendy_bac[3]->Sumw2();

  h_muendy_bac[4]=(TH1D*)input2->Get("fbac_trkend_mucand_y_4");
  h_muendy_bac[4]->Rebin(4);
  h_muendy_bac[4]->Sumw2();
 
  h_muendy_bac[5]=(TH1D*)input2->Get("fbac_trkend_mucand_y_5");
  h_muendy_bac[5]->Rebin(4);
  h_muendy_bac[5]->Sumw2();

  h_muendy_bac[6]=(TH1D*)input2->Get("fbac_trkend_mucand_y_6");
  h_muendy_bac[6]->Rebin(4);
  h_muendy_bac[6]->Sumw2();

  h_muendy_bac[7]=(TH1D*)input2->Get("fbac_trkend_mucand_y_7");
  h_muendy_bac[7]->Rebin(4);
  h_muendy_bac[7]->Sumw2();
   //===================================================================
  TH1D                  *h_muendz_allsel[3];

  h_muendz_allsel[0]=(TH1D*)input0->Get("ftrkendz_mucand");
  h_muendz_allsel[0]->Rebin(4); 
  h_muendz_allsel[0]->Sumw2();
 

  h_muendz_allsel[1]=(TH1D*)input1->Get("ftrkendz_mucand");
  h_muendz_allsel[1]->Rebin(4);
  h_muendz_allsel[1]->Sumw2();

  h_muendz_allsel[2]=(TH1D*)input2->Get("ftrkendz_mucand");
  h_muendz_allsel[2]->Rebin(4);
  h_muendz_allsel[2]->Sumw2();

  TH1D               *h_muendz_sig[2];
  TH1D               *h_muendz_bac[8];
  h_muendz_sig[0]=(TH1D*)input2->Get("fsig_trkend_mucand_z_0");
  h_muendz_sig[0]->Rebin(4);
  h_muendz_sig[0]->Sumw2();

  h_muendz_bac[0]=(TH1D*)input2->Get("fbac_trkend_mucand_z_0");
  h_muendz_bac[0]->Rebin(4);
  h_muendz_bac[0]->Sumw2();
 
  h_muendz_bac[1]=(TH1D*)input2->Get("fbac_trkend_mucand_z_1");
  h_muendz_bac[1]->Rebin(4);
  h_muendz_bac[1]->Sumw2();

  h_muendz_bac[2]=(TH1D*)input2->Get("fbac_trkend_mucand_z_2");
  h_muendz_bac[2]->Rebin(4);
  h_muendz_bac[2]->Sumw2();
 
  h_muendz_bac[3]=(TH1D*)input2->Get("fbac_trkend_mucand_z_3");
  h_muendz_bac[3]->Rebin(4);
  h_muendz_bac[3]->Sumw2();

  h_muendz_bac[4]=(TH1D*)input2->Get("fbac_trkend_mucand_z_4");
  h_muendz_bac[4]->Rebin(4);
  h_muendz_bac[4]->Sumw2();
 
  h_muendz_bac[5]=(TH1D*)input2->Get("fbac_trkend_mucand_z_5");
  h_muendz_bac[5]->Rebin(4);
  h_muendz_bac[5]->Sumw2();

  h_muendz_bac[6]=(TH1D*)input2->Get("fbac_trkend_mucand_z_6");
  h_muendz_bac[6]->Rebin(4);
  h_muendz_bac[6]->Sumw2();

  h_muendz_bac[7]=(TH1D*)input2->Get("fbac_trkend_mucand_z_7");
  h_muendz_bac[7]->Rebin(4);
  h_muendz_bac[7]->Sumw2();
  //=========================================================================
  TH1D                  *h_trunmean_muon_allsel[3];

  h_trunmean_muon_allsel[0]=(TH1D*)input0->Get("trunmeandqdx_mucand");
  h_trunmean_muon_allsel[0]->Rebin(4); 
  h_trunmean_muon_allsel[0]->Sumw2();
 

  h_trunmean_muon_allsel[1]=(TH1D*)input1->Get("trunmeandqdx_mucand");
  h_trunmean_muon_allsel[1]->Rebin(4);
  h_trunmean_muon_allsel[1]->Sumw2();

  h_trunmean_muon_allsel[2]=(TH1D*)input2->Get("trunmeandqdx_mucand");
  h_trunmean_muon_allsel[2]->Rebin(4);
  h_trunmean_muon_allsel[2]->Sumw2();

  TH1D               *h_trunmean_muon_sig[2];
  TH1D               *h_trunmean_muon_bac[8];
  h_trunmean_muon_sig[0]=(TH1D*)input2->Get("trunmeandqdx_mucand_sig_0");
  h_trunmean_muon_sig[0]->Rebin(4);
  h_trunmean_muon_sig[0]->Sumw2();

  h_trunmean_muon_bac[0]=(TH1D*)input2->Get("trunmeandqdx_mucand_bac_0");
  h_trunmean_muon_bac[0]->Rebin(4);
  h_trunmean_muon_bac[0]->Sumw2();
 
  h_trunmean_muon_bac[1]=(TH1D*)input2->Get("trunmeandqdx_mucand_bac_1");
  h_trunmean_muon_bac[1]->Rebin(4);
  h_trunmean_muon_bac[1]->Sumw2();

  h_trunmean_muon_bac[2]=(TH1D*)input2->Get("trunmeandqdx_mucand_bac_2");
  h_trunmean_muon_bac[2]->Rebin(4);
  h_trunmean_muon_bac[2]->Sumw2();
 
  h_trunmean_muon_bac[3]=(TH1D*)input2->Get("trunmeandqdx_mucand_bac_3");
  h_trunmean_muon_bac[3]->Rebin(4);
  h_trunmean_muon_bac[3]->Sumw2();

  h_trunmean_muon_bac[4]=(TH1D*)input2->Get("trunmeandqdx_mucand_bac_4");
  h_trunmean_muon_bac[4]->Rebin(4);
  h_trunmean_muon_bac[4]->Sumw2();
 
  h_trunmean_muon_bac[5]=(TH1D*)input2->Get("trunmeandqdx_mucand_bac_5");
  h_trunmean_muon_bac[5]->Rebin(4);
  h_trunmean_muon_bac[5]->Sumw2();

  h_trunmean_muon_bac[6]=(TH1D*)input2->Get("trunmeandqdx_mucand_bac_6");
  h_trunmean_muon_bac[6]->Rebin(4);
  h_trunmean_muon_bac[6]->Sumw2();

  h_trunmean_muon_bac[7]=(TH1D*)input2->Get("trunmeandqdx_mucand_bac_7");
  h_trunmean_muon_bac[7]->Rebin(4);
  h_trunmean_muon_bac[7]->Sumw2();
  //===================================================================
  TH1D                  *h_trunmean_proton_allsel[3];

  h_trunmean_proton_allsel[0]=(TH1D*)input0->Get("trunmeandqdx_pcand");
  h_trunmean_proton_allsel[0]->Rebin(5); 
  h_trunmean_proton_allsel[0]->Sumw2();
 

  h_trunmean_proton_allsel[1]=(TH1D*)input1->Get("trunmeandqdx_pcand");
  h_trunmean_proton_allsel[1]->Rebin(5);
  h_trunmean_proton_allsel[1]->Sumw2();

  h_trunmean_proton_allsel[2]=(TH1D*)input2->Get("trunmeandqdx_pcand");
  h_trunmean_proton_allsel[2]->Rebin(5);
  h_trunmean_proton_allsel[2]->Sumw2();

  TH1D               *h_trunmean_proton_sig[2];
  TH1D               *h_trunmean_proton_bac[8];
  h_trunmean_proton_sig[0]=(TH1D*)input2->Get("trunmeandqdx_pcand_sig_0");
  h_trunmean_proton_sig[0]->Rebin(5);
  h_trunmean_proton_sig[0]->Sumw2();

  h_trunmean_proton_bac[0]=(TH1D*)input2->Get("trunmeandqdx_pcand_bac_0");
  h_trunmean_proton_bac[0]->Rebin(5);
  h_trunmean_proton_bac[0]->Sumw2();
 
  h_trunmean_proton_bac[1]=(TH1D*)input2->Get("trunmeandqdx_pcand_bac_1");
  h_trunmean_proton_bac[1]->Rebin(5);
  h_trunmean_proton_bac[1]->Sumw2();

  h_trunmean_proton_bac[2]=(TH1D*)input2->Get("trunmeandqdx_pcand_bac_2");
  h_trunmean_proton_bac[2]->Rebin(5);
  h_trunmean_proton_bac[2]->Sumw2();
 
  h_trunmean_proton_bac[3]=(TH1D*)input2->Get("trunmeandqdx_pcand_bac_3");
  h_trunmean_proton_bac[3]->Rebin(5);
  h_trunmean_proton_bac[3]->Sumw2();

  h_trunmean_proton_bac[4]=(TH1D*)input2->Get("trunmeandqdx_pcand_bac_4");
  h_trunmean_proton_bac[4]->Rebin(5);
  h_trunmean_proton_bac[4]->Sumw2();
 
  h_trunmean_proton_bac[5]=(TH1D*)input2->Get("trunmeandqdx_pcand_bac_5");
  h_trunmean_proton_bac[5]->Rebin(5);
  h_trunmean_proton_bac[5]->Sumw2();

  h_trunmean_proton_bac[6]=(TH1D*)input2->Get("trunmeandqdx_pcand_bac_6");
  h_trunmean_proton_bac[6]->Rebin(5);
  h_trunmean_proton_bac[6]->Sumw2();

  h_trunmean_proton_bac[7]=(TH1D*)input2->Get("trunmeandqdx_pcand_bac_7");
  h_trunmean_proton_bac[7]->Rebin(5);
  h_trunmean_proton_bac[7]->Sumw2();
 
  //=============================================================================
  TH1D      *h_neutron_mom_allsel[3];
  h_neutron_mom_allsel[0]=(TH1D*)input0->Get("neutron_mom");
  h_neutron_mom_allsel[0]->Rebin(4);
  h_neutron_mom_allsel[0]->Sumw2();

  h_neutron_mom_allsel[1]=(TH1D*)input1->Get("neutron_mom");
  h_neutron_mom_allsel[1]->Rebin(4);
  h_neutron_mom_allsel[1]->Sumw2();

  h_neutron_mom_allsel[2]=(TH1D*)input2->Get("neutron_mom");
  h_neutron_mom_allsel[2]->Rebin(4);
  h_neutron_mom_allsel[2]->Sumw2();
 
  TH1D      *h_thetax_mu_allsel[3];
  h_thetax_mu_allsel[0]=(TH1D*)input0->Get("thetax_mu");
  h_thetax_mu_allsel[0]->Rebin(4);
  h_thetax_mu_allsel[0]->Sumw2();

  h_thetax_mu_allsel[1]=(TH1D*)input1->Get("thetax_mu");
  h_thetax_mu_allsel[1]->Rebin(4);
  h_thetax_mu_allsel[1]->Sumw2();

  h_thetax_mu_allsel[2]=(TH1D*)input2->Get("thetax_mu");
  h_thetax_mu_allsel[2]->Rebin(4);
  h_thetax_mu_allsel[2]->Sumw2();
 
  TH1D      *h_angletest_0505_mucand_allsel[3];
  h_angletest_0505_mucand_allsel[0]=(TH1D*)input0->Get("h_angletest_0505_mucand");
  h_angletest_0505_mucand_allsel[0]->Rebin(4);
  h_angletest_0505_mucand_allsel[0]->Sumw2();

  h_angletest_0505_mucand_allsel[1]=(TH1D*)input1->Get("h_angletest_0505_mucand");
  h_angletest_0505_mucand_allsel[1]->Rebin(4);
  h_angletest_0505_mucand_allsel[1]->Sumw2();

  h_angletest_0505_mucand_allsel[2]=(TH1D*)input2->Get("h_angletest_0505_mucand");
  h_angletest_0505_mucand_allsel[2]->Rebin(4);
  h_angletest_0505_mucand_allsel[2]->Sumw2();
 
  TH1D      *h_angletest_pi05_mucand_allsel[3];
  h_angletest_pi05_mucand_allsel[0]=(TH1D*)input0->Get("h_angletest_pi05_mucand");
  h_angletest_pi05_mucand_allsel[0]->Rebin(4);
  h_angletest_pi05_mucand_allsel[0]->Sumw2();

  h_angletest_pi05_mucand_allsel[1]=(TH1D*)input1->Get("h_angletest_pi05_mucand");
  h_angletest_pi05_mucand_allsel[1]->Rebin(4);
  h_angletest_pi05_mucand_allsel[1]->Sumw2();

  h_angletest_pi05_mucand_allsel[2]=(TH1D*)input2->Get("h_angletest_pi05_mucand");
  h_angletest_pi05_mucand_allsel[2]->Rebin(4);
  h_angletest_pi05_mucand_allsel[2]->Sumw2();
 
  TH1D      *h_angletest_0505_pcand_allsel[3];
  h_angletest_0505_pcand_allsel[0]=(TH1D*)input0->Get("h_angletest_0505_pcand");
  h_angletest_0505_pcand_allsel[0]->Rebin(4);
  h_angletest_0505_pcand_allsel[0]->Sumw2();

  h_angletest_0505_pcand_allsel[1]=(TH1D*)input1->Get("h_angletest_0505_pcand");
  h_angletest_0505_pcand_allsel[1]->Rebin(4);
  h_angletest_0505_pcand_allsel[1]->Sumw2();

  h_angletest_0505_pcand_allsel[2]=(TH1D*)input2->Get("h_angletest_0505_pcand");
  h_angletest_0505_pcand_allsel[2]->Rebin(4);
  h_angletest_0505_pcand_allsel[2]->Sumw2();
 
  TH1D      *h_angletest_pi05_pcand_allsel[3];
  h_angletest_pi05_pcand_allsel[0]=(TH1D*)input0->Get("h_angletest_pi05_pcand");
  h_angletest_pi05_pcand_allsel[0]->Rebin(4);
  h_angletest_pi05_pcand_allsel[0]->Sumw2();

  h_angletest_pi05_pcand_allsel[1]=(TH1D*)input1->Get("h_angletest_pi05_pcand");
  h_angletest_pi05_pcand_allsel[1]->Rebin(4);
  h_angletest_pi05_pcand_allsel[1]->Sumw2();

  h_angletest_pi05_pcand_allsel[2]=(TH1D*)input2->Get("h_angletest_pi05_pcand");
  h_angletest_pi05_pcand_allsel[2]->Rebin(4);
  h_angletest_pi05_pcand_allsel[2]->Sumw2();
  
  //=============================================================================  
*/
  std::cout<<"Start making plots<<<<<<<<<<<<<<"<<std::endl; 
  ////~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~

  std::cout<<"Calculating POT normalization factors "<<h_range_allsel[0]->GetEntries()<<std::endl;

 
  Float_t E1DCNT_wcut_bnb=35388924.0;
  Float_t E1DCNT_wcut_extbnb=72299264.0;


 

  Float_t mcbnbcos_POT;
  Float_t dirt_POT;

  if (tune==3){mcbnbcos_POT=1.8e20;} // Tune1
  else{mcbnbcos_POT=1.8e20;} // Tune1

  dirt_POT=3.70426e+20;//Dirt


  Float_t dataPOT=1.592e20;// ??????????????/
 
  
  Double_t scalefac=E1DCNT_wcut_bnb/E1DCNT_wcut_extbnb;

  Double_t normfac=dataPOT/mcbnbcos_POT;
  cout<<"normalization factor for monte carlo sample is: "<<normfac<<endl;

  Double_t scale_onoffbeam=0.0;
  scale_onoffbeam=scalefac;
  cout<<"scale factor for extbnb is : "<<scale_onoffbeam<<endl;

  Double_t scale_dirt_MC = mcbnbcos_POT/dirt_POT;
  cout<<"scale factor for dirt to MC is : " << scale_dirt_MC << endl;
  cout<<"scale factor for dirt to data is : " << scale_dirt_MC*normfac << endl;

  bool areanorm=false;

  double areanorm_fac=1.0;
  //if(areanorm==true){
   //  if(tune==1){
     //areanorm_fac=(67.+1413.)/1143.;}
     //if(tune==3){
     //areanorm_fac=(67.+1265.)/1143.;}
  //} else if(areanorm==false){areanorm_fac=1.0;}
  cout<<"areanorm = "<<areanorm_fac<<endl;




//--------------------------------------------------------------------------
TLine *line= new TLine(0, 1, 700, 1);
line->SetLineColor(kRed);
line->SetLineStyle(9);


TLatex* prelim = new TLatex(0.94,0.93, "MicroBooNE Preliminary");
prelim->SetTextFont(62);
prelim->SetTextColor(kGray+2); 
prelim->SetNDC(); 
prelim->SetTextSize(1/30.); 
prelim->SetTextAlign(32); 
//prelim->Draw();


 
//------------------------------------------------------------------
  TH1D *h_range_MCandOff=(TH1D*)h_range_allsel[2]->Clone(Form("%s_MC+Off", h_range_allsel[2]->GetName()));
  cout<<"get all the histograms!"<<endl;

  TCanvas *c1 = new TCanvas("c1"); 
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  pad1->SetBottomMargin(0); // Upper and lower plot are joined
  pad1->SetGridx();         // Vertical grid
  pad1->Draw();             // Draw the upper pad: pad1
  pad1->cd();               // pad1 becomes the current pad
  
  std::cout<<"Draw the histogram of on beam data "<<std::endl;
  h_range_allsel[0]->GetXaxis()->SetTitle("Track Length of Muon Candidate[cm]");
  h_range_allsel[0]->GetYaxis()->SetTitle("No. of Tracks");
  h_range_allsel[0]->GetYaxis()->SetTitleSize(0.06);
  h_range_allsel[0]->GetYaxis()->SetTitleOffset(0.6);
  h_range_allsel[0]->SetMaximum(1000);
  h_range_allsel[0]->Draw();
  std::cout<<"start stacking the mc and off beam data "<<std::endl;  
  //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
  THStack *hs_range = new THStack("hs_range","");
  stackHists(hs_range, h_range_sig, h_range_bac, h_range_allsel, normfac, scale_onoffbeam, scale_dirt_MC);
  hs_range -> Draw("HIST,SAME");
  std::cout<<"stacked the backgrounds to signal"<<std::endl;
  
  //h_range_allsel[0]->Scale(areanorm_fac);
  //h_range_allsel[0]->SetMaximum(400);
  //h_range_allsel[0]->Draw("same");
  
  //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
  TLegend *legendR = new TLegend(.60, .55, .90, .90); // right-aligned
  TLegend *legendL = new TLegend(.10, .55, .35, .90); // left-aligned

  legendR->AddEntry(h_range_allsel[0], "on_beam");
  legendR->AddEntry(h_range_allsel[1], "off_beam");
  legendR -> AddEntry(h_range_sig[0], "Signal", "f");
  legendR -> AddEntry(h_range_bac[0], "Cosmic", "f");
  legendR -> AddEntry(h_range_bac[1], "OUTFV", "f");
  legendR -> AddEntry(h_range_bac[2], "NC", "f");
  legendR -> AddEntry(h_range_bac[3], "#bar{#nu}_{#mu} CC", "f");
  legendR -> AddEntry(h_range_bac[4], "#nu_{e}, #bar{#nu}_{e} CC", "f");
  legendR -> AddEntry(h_range_bac[5], "CC0P, CCpion", "f");
  legendR -> AddEntry(h_range_allsel[3], "dirt");

  legendL->AddEntry(h_range_allsel[0], "on_beam");
  legendL->AddEntry(h_range_allsel[1], "off_beam");
  legendL -> AddEntry(h_range_sig[0], "Signal", "f");
  legendL -> AddEntry(h_range_bac[0], "Cosmic", "f");
  legendL -> AddEntry(h_range_bac[1], "OUTFV", "f");
  legendL -> AddEntry(h_range_bac[2], "NC", "f");
  legendL -> AddEntry(h_range_bac[3], "#bar{#nu}_{#mu} CC", "f");
  legendL -> AddEntry(h_range_bac[4], "#nu_{e}, #bar{#nu}_{e} CC", "f");
  legendL -> AddEntry(h_range_bac[5], "CC0P, CCpion", "f");
  legendL -> AddEntry(h_range_allsel[3], "dirt");


  
  legendR->Draw("same"); prelim->Draw("same");
  
  //Draw ratio between two histograms and save the plot
  // lower plot will be in pad2
  c1->cd();          // Go back to the main canvas before defining pad2
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.2);
  pad2->SetGridx(); // vertical grid
  pad2->Draw();
  pad2->cd();       // pad2 becomes the current pad
   
  TH1D *h_range_temp_0=(TH1D*)h_range_allsel[0]->Clone("h_range_temp_0"); 
  TH1D *h_range_temp_1=(TH1D*)h_range_allsel[1]->Clone("h_range_temp_1"); 
  TH1D *h_range_temp_2=(TH1D*)h_range_allsel[2]->Clone("h_range_temp_2"); 
  h_range_temp_2->Scale(normfac);
  h_range_MCandOff->Add(h_range_temp_2,h_range_temp_1,1,1); 

  h_range_temp_0->Divide(h_range_MCandOff);
  h_range_temp_0->SetMinimum(0);
  h_range_temp_0->SetMaximum(2);
  h_range_temp_0->GetXaxis()->SetTitleSize(0.15);
  h_range_temp_0->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  h_range_temp_0->GetXaxis()->SetLabelSize(11);
 
  
  h_range_temp_0->GetYaxis()->SetTitle("Data/MC");
  h_range_temp_0->GetYaxis()->SetTitleSize(18);
  h_range_temp_0->GetYaxis()->SetTitleFont(43);
  h_range_temp_0->GetYaxis()->SetTitleOffset(0.9);
  h_range_temp_0->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  h_range_temp_0->GetYaxis()->SetLabelSize(11);
  h_range_temp_0->Draw();
  line->Draw("same");
  c1->cd();
  
  std::cout<<"libo test 0"<<std::endl;
  if (cosmicCut){
    if(tune==3){c1->Print("figures/Tune3/BackSep/h_range_allsel.png");}
    else{c1->Print("figures/Tune1/BackSep/h_range_allsel.png");}
  }
  else{
    if(tune==3){c1->Print("figures_noCosmicCut/Tune3/BackSep/h_range_allsel.png");}
    else{c1->Print("figures_noCosmicCut/Tune1/BackSep/h_range_allsel.png");}
  }
  std::cout<<"libo test 1"<<std::endl; 
  int nbins_range=h_range_allsel[2]->GetNbinsX();
  float chi2_range=Chi2Calc(h_range_allsel[2], h_range_allsel[0], h_range_allsel[1],scalefac, normfac);
  cout<<"chi2_range = "<<chi2_range<<"  number of bins: "<<nbins_range<<"  chi2/dof is : "<<chi2_range/nbins_range<<endl;

 //==============================================================================
 line= new TLine(0, 1, 150, 1);
 line->SetLineColor(kRed);
 line->SetLineStyle(9);

  TH1D *h_prange_MCandOff=(TH1D*)h_prange_allsel[2]->Clone(Form("%s_MC+Off", h_prange_allsel[2]->GetName()));
  cout<<"get all the histograms!"<<endl;
  pad1->cd();               // pad1 becomes the current pad
  
  std::cout<<"Draw the histogram of on beam data "<<std::endl;
  h_prange_allsel[0]->GetXaxis()->SetTitle("Track Length of Proton Candidate[cm]");
  h_prange_allsel[0]->GetYaxis()->SetTitle("No. of Tracks");
  h_prange_allsel[0]->GetYaxis()->SetTitleSize(0.06);
  h_prange_allsel[0]->GetYaxis()->SetTitleOffset(0.6);
  h_prange_allsel[0]->Draw();
  std::cout<<"start stacking the mc and off beam data "<<std::endl;  
  //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
  THStack *hs_prange = new THStack("hs_prange","");
  stackHists(hs_prange, h_prange_sig, h_prange_bac, h_prange_allsel, normfac, scale_onoffbeam, scale_dirt_MC);
  hs_prange -> Draw("HIST,SAME");
  std::cout<<"stacked the backgrounds to signal"<<std::endl;
  
  h_prange_allsel[0]->Scale(areanorm_fac);
  h_prange_allsel[0]->SetMaximum(1400);
  h_prange_allsel[0]->Draw("same");
  
  //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
  legendR->Draw("same"); prelim->Draw("same");
  
  //Draw ratio between two histograms and save the plot
  // lower plot will be in pad2
  c1->cd();          // Go back to the main canvas before defining pad2
  pad2->cd();       // pad2 becomes the current pad
   
  TH1D *h_prange_temp_0=(TH1D*)h_prange_allsel[0]->Clone("h_prange_temp_0"); 
  TH1D *h_prange_temp_1=(TH1D*)h_prange_allsel[1]->Clone("h_prange_temp_1"); 
  TH1D *h_prange_temp_2=(TH1D*)h_prange_allsel[2]->Clone("h_prange_temp_2"); 
  h_prange_temp_2->Scale(normfac);
  h_prange_MCandOff->Add(h_prange_temp_2,h_prange_temp_1,1,1); 

  h_prange_temp_0->Divide(h_prange_MCandOff);
  h_prange_temp_0->SetMinimum(0);
  h_prange_temp_0->SetMaximum(2);
  h_prange_temp_0->GetXaxis()->SetTitleSize(0.15);
  h_prange_temp_0->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  h_prange_temp_0->GetXaxis()->SetLabelSize(11);
 
  
  h_prange_temp_0->GetYaxis()->SetTitle("Data/MC");
  h_prange_temp_0->GetYaxis()->SetTitleSize(18);
  h_prange_temp_0->GetYaxis()->SetTitleFont(43);
  h_prange_temp_0->GetYaxis()->SetTitleOffset(0.9);
  h_prange_temp_0->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  h_prange_temp_0->GetYaxis()->SetLabelSize(11);
  h_prange_temp_0->Draw();
  line->Draw("same");
  c1->cd();
  
  std::cout<<"libo test 0"<<std::endl;
  if (cosmicCut){
    if(tune==3){c1->Print("figures/Tune3/BackSep/h_prange_allsel.png");}
    else{c1->Print("figures/Tune1/BackSep/h_prange_allsel.png");}
  }
  else{
    if(tune==3){c1->Print("figures_noCosmicCut/Tune3/BackSep/h_prange_allsel.png");}
    else{c1->Print("figures_noCosmicCut/Tune1/BackSep/h_prange_allsel.png");}
  }
  std::cout<<"libo test 1"<<std::endl; 
  int nbins_prange=h_prange_allsel[2]->GetNbinsX();
  float chi2_prange=Chi2Calc(h_prange_allsel[2], h_prange_allsel[0], h_prange_allsel[1],scalefac, normfac);
  cout<<"chi2_prange = "<<chi2_prange<<"  number of bins: "<<nbins_prange<<"  chi2/dof is : "<<chi2_prange/nbins_prange<<endl;




   
  //=============================================================================================
  line= new TLine(-3.14, 1, 3.14, 1);
  line->SetLineColor(kRed);
  line->SetLineStyle(9);
  TCanvas *c2 = new TCanvas("c2"); 

  TH1D *h_onoff_phi=(TH1D*)h_phi_allsel[1]->Clone(Form("%s_on-off", h_phi_allsel[1]->GetName()));
  TH1D *h_phi_MCandOff=(TH1D*)h_phi_allsel[2]->Clone(Form("%s_MC+Off", h_phi_allsel[2]->GetName()));
  cout<<"get all the histograms!"<<endl;

  pad1->SetBottomMargin(0); // Upper and lower plot are joined
  pad1->SetGridx();         // Vertical grid
  pad1->Draw();             // Draw the upper pad: pad1
  pad1->cd();               // pad1 becomes the current pad

  h_phi_allsel[0]->GetXaxis()->SetTitle("#phi of Muon Candidate [rad]");
  h_phi_allsel[0]->GetYaxis()->SetTitle("No. of Tracks");
  h_phi_allsel[0]->GetYaxis()->SetTitleSize(0.06);
  h_phi_allsel[0]->GetYaxis()->SetTitleOffset(0.6);
  //h_phi_allsel[0]->SetMaximum(500);
  h_phi_allsel[0]->SetMinimum(0.0);
  h_phi_allsel[0]->Draw();
  
  THStack *hs_phi = new THStack("hs_phi","");
  stackHists(hs_phi, h_phi_sig, h_phi_bac, h_phi_allsel, normfac, scale_onoffbeam, scale_dirt_MC);
  hs_phi -> Draw("HIST,SAME");

  h_phi_allsel[0]->Scale(areanorm_fac);
  h_phi_allsel[0]->SetMaximum(800);
  h_phi_allsel[0]->SetMinimum(0);
  h_phi_allsel[0]->Draw("same");
  //h_onoff_phi->Draw("E1CSAME");
  //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
  
  legendL->Draw("same"); prelim->Draw("same");
  
  c2->cd();          // Go back to the main canvas before defining pad2
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.2);
  pad2->SetGridx(); // vertical grid
  pad2->Draw();
  pad2->cd();       // pad2 becomes the current pad

  TH1D *h_phi_temp_0=(TH1D*)h_phi_allsel[0]->Clone("h_phi_temp_0"); 
  TH1D *h_phi_temp_1=(TH1D*)h_phi_allsel[1]->Clone("h_phi_temp_1"); 
  TH1D *h_phi_temp_2=(TH1D*)h_phi_allsel[2]->Clone("h_phi_temp_2"); 
  h_phi_temp_2->Scale(normfac);
  h_phi_MCandOff->Add(h_phi_temp_2,h_phi_temp_1,1,1); 

  h_phi_temp_0->Divide(h_phi_MCandOff);
  h_phi_temp_0->SetMinimum(0);
  h_phi_temp_0->SetMaximum(2);
  h_phi_temp_0->GetXaxis()->SetTitleSize(0.15);
  h_phi_temp_0->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  h_phi_temp_0->GetXaxis()->SetLabelSize(11);
 

  h_phi_temp_0->GetYaxis()->SetTitle("Data/MC");
  h_phi_temp_0->GetYaxis()->SetTitleSize(18);
  h_phi_temp_0->GetYaxis()->SetTitleFont(43);
  h_phi_temp_0->GetYaxis()->SetTitleOffset(0.9);
  h_phi_temp_0->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  h_phi_temp_0->GetYaxis()->SetLabelSize(11);
  h_phi_temp_0->Draw();
  line->Draw("same");
  c2->cd();
  
  if (cosmicCut){
    if(tune==3){c2->Print("figures/Tune3/BackSep/h_phi_allsel.png");}
    else{c2->Print("figures/Tune1/BackSep/h_phi_allsel.png");}
  }
  else{
    if(tune==3){c2->Print("figures_noCosmicCut/Tune3/BackSep/h_phi_allsel.png");}
    else{c2->Print("figures_noCosmicCut/Tune1/BackSep/h_phi_allsel.png");}
  }
  int nbins_phi=h_phi_allsel[2]->GetNbinsX();
  float chi2_phi=Chi2Calc(h_phi_allsel[2], h_phi_allsel[0], h_phi_allsel[1],scalefac, normfac);
  cout<<"chi2_phi = "<<chi2_phi<<"  number of bins: "<<nbins_phi<<"  chi2/dof is : "<<chi2_phi/nbins_phi<<endl;
 
  c2->Update();   
  TH1D *h_onoff_pphi=(TH1D*)h_pphi_allsel[1]->Clone(Form("%s_on-off", h_pphi_allsel[1]->GetName()));
  TH1D *h_pphi_MCandOff=(TH1D*)h_pphi_allsel[2]->Clone(Form("%s_MC+Off", h_pphi_allsel[2]->GetName()));
  cout<<"get all the histograms!"<<endl;

  pad1->SetBottomMargin(0); // Upper and lower plot are joined
  pad1->SetGridx();         // Vertical grid
  pad1->Draw();             // Draw the upper pad: pad1
  pad1->cd();               // pad1 becomes the current pad

  h_pphi_allsel[0]->GetXaxis()->SetTitle("#phi of Proton Candidate [rad]");
  h_pphi_allsel[0]->GetYaxis()->SetTitle("No. of Tracks");
  h_pphi_allsel[0]->GetYaxis()->SetTitleSize(0.06);
  h_pphi_allsel[0]->GetYaxis()->SetTitleOffset(0.6);
  //h_pphi_allsel[0]->SetMaximum(500);
  h_pphi_allsel[0]->SetMinimum(0.0);
  h_pphi_allsel[0]->Draw();
  
  THStack *hs_pphi = new THStack("hs_pphi","");
  stackHists(hs_pphi, h_pphi_sig, h_pphi_bac, h_pphi_allsel, normfac, scale_onoffbeam, scale_dirt_MC);
  hs_pphi -> Draw("HIST,SAME");

  h_pphi_allsel[0]->Scale(areanorm_fac);
  h_pphi_allsel[0]->SetMaximum(800);
  h_pphi_allsel[0]->SetMinimum(0);
  h_pphi_allsel[0]->Draw("same");
  //h_onoff_pphi->Draw("E1CSAME");
  //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
  
  legendL->Draw("same"); prelim->Draw("same");
  
  c2->cd();          // Go back to the main canvas before defining pad2
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.2);
  pad2->SetGridx(); // vertical grid
  pad2->Draw();
  pad2->cd();       // pad2 becomes the current pad

  TH1D *h_pphi_temp_0=(TH1D*)h_pphi_allsel[0]->Clone("h_pphi_temp_0"); 
  TH1D *h_pphi_temp_1=(TH1D*)h_pphi_allsel[1]->Clone("h_pphi_temp_1"); 
  TH1D *h_pphi_temp_2=(TH1D*)h_pphi_allsel[2]->Clone("h_pphi_temp_2"); 
  h_pphi_temp_2->Scale(normfac);
  h_pphi_MCandOff->Add(h_pphi_temp_2,h_pphi_temp_1,1,1); 

  h_pphi_temp_0->Divide(h_pphi_MCandOff);
  h_pphi_temp_0->SetMinimum(0);
  h_pphi_temp_0->SetMaximum(2);
  h_pphi_temp_0->GetXaxis()->SetTitleSize(0.15);
  h_pphi_temp_0->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  h_pphi_temp_0->GetXaxis()->SetLabelSize(11);
 

  h_pphi_temp_0->GetYaxis()->SetTitle("Data/MC");
  h_pphi_temp_0->GetYaxis()->SetTitleSize(18);
  h_pphi_temp_0->GetYaxis()->SetTitleFont(43);
  h_pphi_temp_0->GetYaxis()->SetTitleOffset(0.9);
  h_pphi_temp_0->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  h_pphi_temp_0->GetYaxis()->SetLabelSize(11);
  h_pphi_temp_0->Draw();
  line->Draw("same");
  c2->cd();
  
  if (cosmicCut){
    if(tune==3){c2->Print("figures/Tune3/BackSep/h_pphi_allsel.png");}
    else{c2->Print("figures/Tune1/BackSep/h_pphi_allsel.png");}
  }
  else{
    if(tune==3){c2->Print("figures_noCosmicCut/Tune3/BackSep/h_pphi_allsel.png");}
    else{c2->Print("figures_noCosmicCut/Tune1/BackSep/h_pphi_allsel.png");}
  }
  int nbins_pphi=h_pphi_allsel[2]->GetNbinsX();
  float chi2_pphi=Chi2Calc(h_pphi_allsel[2], h_pphi_allsel[0], h_pphi_allsel[1],scalefac, normfac);
  cout<<"chi2_pphi = "<<chi2_pphi<<"  number of bins: "<<nbins_pphi<<"  chi2/dof is : "<<chi2_pphi/nbins_pphi<<endl;
 



  //==========================================================================================
  line= new TLine(-1., 1, 1., 1);
  line->SetLineColor(kRed);
  line->SetLineStyle(9);
  TCanvas *c3 = new TCanvas("c3"); 

  TH1D *h_onoff_costheta=(TH1D*)h_costheta_allsel[1]->Clone(Form("%s_on-off", h_costheta_allsel[1]->GetName()));
  TH1D *h_costheta_MCandOff=(TH1D*)h_costheta_allsel[2]->Clone(Form("%s_MC+Off", h_costheta_allsel[2]->GetName()));
  cout<<"get all the histograms!"<<endl;

  pad1->SetBottomMargin(0); // Upper and lower plot are joined
  pad1->SetGridx();         // Vertical grid
  pad1->Draw();             // Draw the upper pad: pad1
  pad1->cd();               // pad1 becomes the current pad

  h_costheta_allsel[0]->GetXaxis()->SetTitle("Track CosTheta of Muon Candidate");
  h_costheta_allsel[0]->GetYaxis()->SetTitle("No. of Tracks");
  h_costheta_allsel[0]->GetYaxis()->SetTitleSize(0.06);
  h_costheta_allsel[0]->GetYaxis()->SetTitleOffset(0.6);
  //h_costheta_allsel[0]->SetMaximum(500);
  h_costheta_allsel[0]->SetMinimum(0.0);
  h_costheta_allsel[0]->Draw();
  
  THStack *hs_costheta = new THStack("hs_costheta","");
  stackHists(hs_costheta, h_costheta_sig, h_costheta_bac, h_costheta_allsel, normfac, scale_onoffbeam, scale_dirt_MC);
  hs_costheta -> Draw("HIST,SAME");

  h_costheta_allsel[0]->Scale(areanorm_fac);
  h_costheta_allsel[0]->SetMaximum(1200);
  h_costheta_allsel[0]->Draw("same");
  //h_onoff_costheta->Draw("E1CSAME");
  //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
  
  legendL->Draw("same"); prelim->Draw("same");
  
  c3->cd();          // Go back to the main canvas before defining pad2
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.2);
  pad2->SetGridx(); // vertical grid
  pad2->Draw();
  pad2->cd();       // pad2 becomes the current pad

  TH1D *h_costheta_temp_0=(TH1D*)h_costheta_allsel[0]->Clone("h_costheta_temp_0"); 
  TH1D *h_costheta_temp_1=(TH1D*)h_costheta_allsel[1]->Clone("h_costheta_temp_1"); 
  TH1D *h_costheta_temp_2=(TH1D*)h_costheta_allsel[2]->Clone("h_costheta_temp_2"); 
  h_costheta_temp_2->Scale(normfac);
  h_costheta_MCandOff->Add(h_costheta_temp_2,h_costheta_temp_1,1,1); 

  h_costheta_temp_0->Divide(h_costheta_MCandOff);
  h_costheta_temp_0->SetMinimum(0);
  h_costheta_temp_0->SetMaximum(2);
  h_costheta_temp_0->GetXaxis()->SetTitleSize(0.15);
  h_costheta_temp_0->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  h_costheta_temp_0->GetXaxis()->SetLabelSize(11);
 

  h_costheta_temp_0->GetYaxis()->SetTitle("Data/MC");
  h_costheta_temp_0->GetYaxis()->SetTitleSize(18);
  h_costheta_temp_0->GetYaxis()->SetTitleFont(43);
  h_costheta_temp_0->GetYaxis()->SetTitleOffset(0.9);
  h_costheta_temp_0->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  h_costheta_temp_0->GetYaxis()->SetLabelSize(11);
  h_costheta_temp_0->Draw();
  line->Draw("same");
  c3->cd();
  
  if (cosmicCut){
    if(tune==3){c3->Print("figures/Tune3/BackSep/h_costheta_allsel.png");}
    else{c3->Print("figures/Tune1/BackSep/h_costheta_allsel.png");}
  }
  else{
    if(tune==3){c3->Print("figures_noCosmicCut/Tune3/BackSep/h_costheta_allsel.png");}
    else{c3->Print("figures_noCosmicCut/Tune1/BackSep/h_costheta_allsel.png");}
  }
  int nbins_costheta=h_costheta_allsel[2]->GetNbinsX();
  float chi2_costheta=Chi2Calc(h_costheta_allsel[2], h_costheta_allsel[0], h_costheta_allsel[1],scalefac, normfac);
  cout<<"chi2_costheta = "<<chi2_costheta<<"  number of bins: "<<nbins_costheta<<"  chi2/dof is : "<<chi2_costheta/nbins_costheta<<endl;
 
 //=========================================================================================  

  TH1D *h_onoff_pcostheta=(TH1D*)h_pcostheta_allsel[1]->Clone(Form("%s_on-off", h_pcostheta_allsel[1]->GetName()));
  TH1D *h_pcostheta_MCandOff=(TH1D*)h_pcostheta_allsel[2]->Clone(Form("%s_MC+Off", h_pcostheta_allsel[2]->GetName()));
  cout<<"get all the histograms!"<<endl;

  pad1->SetBottomMargin(0); // Upper and lower plot are joined
  pad1->SetGridx();         // Vertical grid
  pad1->Draw();             // Draw the upper pad: pad1
  pad1->cd();               // pad1 becomes the current pad

  h_pcostheta_allsel[0]->GetXaxis()->SetTitle("Track CosTheta of Proton Candidate");
  h_pcostheta_allsel[0]->GetYaxis()->SetTitle("No. of Tracks");
  h_pcostheta_allsel[0]->GetYaxis()->SetTitleSize(0.06);
  h_pcostheta_allsel[0]->GetYaxis()->SetTitleOffset(0.6);
  //h_pcostheta_allsel[0]->SetMaximum(500);
  h_pcostheta_allsel[0]->SetMinimum(0.0);
  h_pcostheta_allsel[0]->Draw();
  
  THStack *hs_pcostheta = new THStack("hs_pcostheta","");
  stackHists(hs_pcostheta, h_pcostheta_sig, h_pcostheta_bac, h_pcostheta_allsel, normfac, scale_onoffbeam, scale_dirt_MC);
  hs_pcostheta -> Draw("HIST,SAME");

  h_pcostheta_allsel[0]->Scale(areanorm_fac);
  h_pcostheta_allsel[0]->SetMaximum(800);
  h_pcostheta_allsel[0]->Draw("same");
  //h_onoff_pcostheta->Draw("E1CSAME");
  //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
  
  legendL->Draw("same"); prelim->Draw("same");
  
  c3->cd();          // Go back to the main canvas before defining pad2
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.2);
  pad2->SetGridx(); // vertical grid
  pad2->Draw();
  pad2->cd();       // pad2 becomes the current pad

  TH1D *h_pcostheta_temp_0=(TH1D*)h_pcostheta_allsel[0]->Clone("h_pcostheta_temp_0"); 
  TH1D *h_pcostheta_temp_1=(TH1D*)h_pcostheta_allsel[1]->Clone("h_pcostheta_temp_1"); 
  TH1D *h_pcostheta_temp_2=(TH1D*)h_pcostheta_allsel[2]->Clone("h_pcostheta_temp_2"); 
  h_pcostheta_temp_2->Scale(normfac);
  h_pcostheta_MCandOff->Add(h_pcostheta_temp_2,h_pcostheta_temp_1,1,1); 

  h_pcostheta_temp_0->Divide(h_pcostheta_MCandOff);
  h_pcostheta_temp_0->SetMinimum(0);
  h_pcostheta_temp_0->SetMaximum(2);
  h_pcostheta_temp_0->GetXaxis()->SetTitleSize(0.15);
  h_pcostheta_temp_0->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  h_pcostheta_temp_0->GetXaxis()->SetLabelSize(11);
 

  h_pcostheta_temp_0->GetYaxis()->SetTitle("Data/MC");
  h_pcostheta_temp_0->GetYaxis()->SetTitleSize(18);
  h_pcostheta_temp_0->GetYaxis()->SetTitleFont(43);
  h_pcostheta_temp_0->GetYaxis()->SetTitleOffset(0.9);
  h_pcostheta_temp_0->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  h_pcostheta_temp_0->GetYaxis()->SetLabelSize(11);
  h_pcostheta_temp_0->Draw();
  line->Draw("same");
  c3->cd();
  
  if (cosmicCut){
    if(tune==3){c3->Print("figures/Tune3/BackSep/h_pcostheta_allsel.png");}
    else{c3->Print("figures/Tune1/BackSep/h_pcostheta_allsel.png");}
  }
  else{
    if(tune==3){c3->Print("figures_noCosmicCut/Tune3/BackSep/h_pcostheta_allsel.png");}
    else{c3->Print("figures_noCosmicCut/Tune1/BackSep/h_pcostheta_allsel.png");}
  }
  int nbins_pcostheta=h_pcostheta_allsel[2]->GetNbinsX();
  float chi2_pcostheta=Chi2Calc(h_pcostheta_allsel[2], h_pcostheta_allsel[0], h_pcostheta_allsel[1],scalefac, normfac);
  cout<<"chi2_pcostheta = "<<chi2_pcostheta<<"  number of bins: "<<nbins_pcostheta<<"  chi2/dof is : "<<chi2_pcostheta/nbins_pcostheta<<endl;
   
 //================================================================
  cout<<"start using canvas c5<<<<<<<<<<<<<<<<<<<<<<<<"<<endl;
  //===============================================================
  TH1D *h_onoff_plep=(TH1D*)h_plep_allsel[1]->Clone(Form("%s_on-off", h_plep_allsel[1]->GetName()));
  TH1D *h_plep_MCandOff=(TH1D*)h_plep_allsel[2]->Clone(Form("%s_MC+Off", h_plep_allsel[2]->GetName()));
  cout<<"get all the histograms!"<<endl;
  TCanvas *c5 = new TCanvas("c5"); 
  pad1->SetBottomMargin(0); // Upper and lower plot are joined
  pad1->SetGridx();         // Vertical grid
  pad1->Draw();             // Draw the upper pad: pad1
  pad1->cd();               // pad1 becomes the current pad


  h_plep_allsel[0]->SetLineColor(kBlack);
  h_plep_allsel[0]->SetLineWidth(2);
  h_plep_allsel[0]->SetLineStyle(1);
  h_plep_allsel[0]->GetXaxis()->SetTitle("Momentum of the muon candidate[GeV/c]");
  h_plep_allsel[0]->GetYaxis()->SetTitle("No. of Tracks");
  h_plep_allsel[0]->GetYaxis()->SetTitleSize(0.06);
  h_plep_allsel[0]->GetYaxis()->SetTitleOffset(0.6);
  //h_plep_allsel[0]->SetMaximum(1.5*h_plep_allsel[0]->GetMaximum());
  //h_plep_allsel[0]->SetMaximum(1200);
  h_plep_allsel[0]->Draw();  

  THStack *hs_plep = new THStack("hs_plep","");
  stackHists(hs_plep, h_plep_sig, h_plep_bac, h_plep_allsel, normfac, scale_onoffbeam, scale_dirt_MC);
  hs_plep -> Draw("HIST,SAME");

  h_plep_allsel[0]->Scale(areanorm_fac);
  h_plep_allsel[0]->SetMaximum(1400);
  //h_plep_allsel[0]->Draw("same");  
  //h_onoff_plep->Draw("E1CSAME");
  //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
 
  legendR->Draw("same"); prelim->Draw("same");

  //Draw ratio between two histograms and save the plot
  // lower plot will be in pad2
  c5->cd();          // Go back to the main canvas before defining pad2
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.2);
  pad2->SetGridx(); // vertical grid
  pad2->Draw();
  pad2->cd();       // pad2 becomes the current pad

  TH1D *h_plep_temp_0=(TH1D*)h_plep_allsel[0]->Clone("h_plep_temp_0"); 
  TH1D *h_plep_temp_1=(TH1D*)h_plep_allsel[1]->Clone("h_plep_temp_1"); 
  TH1D *h_plep_temp_2=(TH1D*)h_plep_allsel[2]->Clone("h_plep_temp_2"); 
  h_plep_temp_2->Scale(normfac);
  h_plep_MCandOff->Add(h_plep_temp_2,h_plep_temp_1,1,1); 

  h_plep_temp_0->Divide(h_plep_MCandOff);
  h_plep_temp_0->SetMinimum(0);
  h_plep_temp_0->SetMaximum(2);
  h_plep_temp_0->GetXaxis()->SetTitleSize(0.15);
  h_plep_temp_0->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  h_plep_temp_0->GetXaxis()->SetLabelSize(11);
 

  h_plep_temp_0->GetYaxis()->SetTitle("Data/MC");
  h_plep_temp_0->GetYaxis()->SetTitleSize(18);
  h_plep_temp_0->GetYaxis()->SetTitleFont(43);
  h_plep_temp_0->GetYaxis()->SetTitleOffset(0.9);
  h_plep_temp_0->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  h_plep_temp_0->GetYaxis()->SetLabelSize(11);
  h_plep_temp_0->Draw();
  line= new TLine(0, 1, 2.5, 1);
  line->SetLineColor(kRed);
  line->SetLineStyle(9);
   line->Draw("same");
 
  c5->cd();


   if (cosmicCut){
     if(tune==3){c5->Print("figures/Tune3/BackSep/h_plep_allsel.png");}
     else{c5->Print("figures/Tune1/BackSep/h_plep_allsel.png");}
   }
   else{
     if(tune==3){c5->Print("figures_noCosmicCut/Tune3/BackSep/h_plep_allsel.png");}
     else{c5->Print("figures_noCosmicCut/Tune1/BackSep/h_plep_allsel.png");}
   }
  
  int nbins_plep=h_plep_allsel[2]->GetNbinsX();
  float chi2_plep=Chi2Calc(h_plep_allsel[2], h_plep_allsel[0], h_plep_allsel[1],scalefac, normfac);
  cout<<"chi2_plep = "<<chi2_plep<<"  number of bins: "<<nbins_plep<<"  chi2/dof is : "<<chi2_plep/nbins_plep<<endl;

 //==============================================================

  c5->Update();
  TH1D *h_onoff_phad=(TH1D*)h_phad_allsel[1]->Clone(Form("%s_on-off", h_phad_allsel[1]->GetName()));
  TH1D *h_phad_MCandOff=(TH1D*)h_phad_allsel[2]->Clone(Form("%s_MC+Off", h_phad_allsel[2]->GetName()));
  cout<<"get all the histograms!"<<endl;

  

  pad1->SetBottomMargin(0); // Upper and lower plot are joined
  pad1->SetGridx();         // Vertical grid
  pad1->Draw();             // Draw the upper pad: pad1
  pad1->cd();               // pad1 becomes the current pad


  h_phad_allsel[0]->SetLineColor(kBlack);
  h_phad_allsel[0]->SetLineWidth(2);
  h_phad_allsel[0]->SetLineStyle(1);
  h_phad_allsel[0]->GetXaxis()->SetTitle("Momentum of the proton candidate[GeV/c]");
  h_phad_allsel[0]->GetYaxis()->SetTitle("No. of Tracks");
  h_phad_allsel[0]->GetYaxis()->SetTitleSize(0.06);
  h_phad_allsel[0]->GetYaxis()->SetTitleOffset(0.6);
  //h_phad_allsel[0]->SetMaximum(1.4*h_phad_allsel[0]->GetMaximum());
  h_phad_allsel[0]->SetMaximum(800);
  h_phad_allsel[0]->Draw();  

  THStack *hs_phad = new THStack("hs_phad","");
  stackHists(hs_phad, h_phad_sig, h_phad_bac, h_phad_allsel, normfac, scale_onoffbeam, scale_dirt_MC);
  hs_phad -> Draw("HIST,SAME");

  h_phad_allsel[0]->Scale(areanorm_fac);
  h_phad_allsel[0]->SetMaximum(1000);
  h_phad_allsel[0]->Draw("same");  
   //h_onoff_phad->Draw("E1CSAME");
  //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
 
  legendR->Draw("same"); prelim->Draw("same");

  c5->cd();          // Go back to the main canvas before defining pad2
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.2);
  pad2->SetGridx(); // vertical grid
  pad2->Draw();
  pad2->cd();       // pad2 becomes the current pad

  TH1D *h_phad_temp_0=(TH1D*)h_phad_allsel[0]->Clone("h_phad_temp_0"); 
  TH1D *h_phad_temp_1=(TH1D*)h_phad_allsel[1]->Clone("h_phad_temp_1"); 
  TH1D *h_phad_temp_2=(TH1D*)h_phad_allsel[2]->Clone("h_phad_temp_2"); 
  h_phad_temp_2->Scale(normfac);
  h_phad_MCandOff->Add(h_phad_temp_2,h_phad_temp_1,1,1); 

  h_phad_temp_0->Divide(h_phad_MCandOff);
  h_phad_temp_0->SetMinimum(0);
  h_phad_temp_0->SetMaximum(2);
  h_phad_temp_0->GetXaxis()->SetTitleSize(0.15);
  h_phad_temp_0->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  h_phad_temp_0->GetXaxis()->SetLabelSize(11);
 

  h_phad_temp_0->GetYaxis()->SetTitle("Data/MC");
  h_phad_temp_0->GetYaxis()->SetTitleSize(18);
  h_phad_temp_0->GetYaxis()->SetTitleFont(43);
  h_phad_temp_0->GetYaxis()->SetTitleOffset(0.9);
  h_phad_temp_0->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  h_phad_temp_0->GetYaxis()->SetLabelSize(11);
  h_phad_temp_0->Draw();
  line= new TLine(0, 1, 1.5, 1);
  line->SetLineColor(kRed);
  line->SetLineStyle(9);
  line->Draw("same");
 
  c5->cd();


   if (cosmicCut){
     if(tune==3){c5->Print("figures/Tune3/BackSep/h_phad_allsel.png");}
     else{c5->Print("figures/Tune1/BackSep/h_phad_allsel.png");}
   }
   else{
     if(tune==3){c5->Print("figures_noCosmicCut/Tune3/BackSep/h_phad_allsel.png");}
     else{c5->Print("figures_noCosmicCut/Tune1/BackSep/h_phad_allsel.png");}
   }

  int nbins_phad=h_phad_allsel[2]->GetNbinsX();
  float chi2_phad=Chi2Calc(h_phad_allsel[2], h_phad_allsel[0], h_phad_allsel[1],scalefac, normfac);
  cout<<"chi2_phad = "<<chi2_phad<<"  number of bins: "<<nbins_phad<<"  chi2/dof is : "<<chi2_phad/nbins_phad<<endl;

  //=============================================================================================================
  //================================================================================
  TH1D *h_onoff_thetamup=(TH1D*)h_thetamup_allsel[1]->Clone(Form("%s_on-off", h_thetamup_allsel[1]->GetName()));
  TH1D *h_thetamup_MCandOff=(TH1D*)h_thetamup_allsel[2]->Clone(Form("%s_MC+Off", h_thetamup_allsel[2]->GetName()));
  cout<<"get all the histograms!"<<endl;

  
  TCanvas *c6 = new TCanvas("c6"); 

  pad1->SetBottomMargin(0); // Upper and lower plot are joined
  pad1->SetGridx();         // Vertical grid
  pad1->Draw();             // Draw the upper pad: pad1
  pad1->cd();               // pad1 becomes the current pad


  h_thetamup_allsel[0]->SetLineColor(kBlack);
  h_thetamup_allsel[0]->SetLineWidth(2);
  h_thetamup_allsel[0]->SetLineStyle(1);
  h_thetamup_allsel[0]->GetXaxis()->SetTitle("#theta_{#mu p}[rad]");
  h_thetamup_allsel[0]->GetYaxis()->SetTitle("Selected Events");
  h_thetamup_allsel[0]->GetYaxis()->SetTitleSize(0.06);
  h_thetamup_allsel[0]->GetYaxis()->SetTitleOffset(0.6);
  //h_thetamup_allsel[0]->SetMaximum(1.5*h_thetamup_allsel[0]->GetMaximum());
  //h_thetamup_allsel[0]->SetMaximum(1200);
  h_thetamup_allsel[0]->Draw();  

  THStack *hs_thetamup = new THStack("hs_thetamup","");
  stackHists(hs_thetamup, h_thetamup_sig, h_thetamup_bac, h_thetamup_allsel, normfac, scale_onoffbeam, scale_dirt_MC);
  hs_thetamup -> Draw("HIST,SAME");

  h_thetamup_allsel[0]->Scale(areanorm_fac);
  h_thetamup_allsel[0]->SetMaximum(600);
  //h_thetamup_allsel[0]->Draw("same");  
  //h_onoff_thetamup->Draw("E1CSAME");
  //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
 
  legendR->Draw("same"); prelim->Draw("same");

  //Draw ratio between two histograms and save the plot
  // lower plot will be in pad2
  c6->cd();          // Go back to the main canvas before defining pad2
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.2);
  pad2->SetGridx(); // vertical grid
  pad2->Draw();
  pad2->cd();       // pad2 becomes the current pad

  TH1D *h_thetamup_temp_0=(TH1D*)h_thetamup_allsel[0]->Clone("h_thetamup_temp_0"); 
  TH1D *h_thetamup_temp_1=(TH1D*)h_thetamup_allsel[1]->Clone("h_thetamup_temp_1"); 
  TH1D *h_thetamup_temp_2=(TH1D*)h_thetamup_allsel[2]->Clone("h_thetamup_temp_2"); 
  h_thetamup_temp_2->Scale(normfac);
  h_thetamup_MCandOff->Add(h_thetamup_temp_2,h_thetamup_temp_1,1,1); 

  h_thetamup_temp_0->Divide(h_thetamup_MCandOff);
  h_thetamup_temp_0->SetMinimum(0);
  h_thetamup_temp_0->SetMaximum(2);
  h_thetamup_temp_0->GetXaxis()->SetTitleSize(0.15);
  h_thetamup_temp_0->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  h_thetamup_temp_0->GetXaxis()->SetLabelSize(11);
 

  h_thetamup_temp_0->GetYaxis()->SetTitle("Data/MC");
  h_thetamup_temp_0->GetYaxis()->SetTitleSize(18);
  h_thetamup_temp_0->GetYaxis()->SetTitleFont(43);
  h_thetamup_temp_0->GetYaxis()->SetTitleOffset(0.9);
  h_thetamup_temp_0->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  h_thetamup_temp_0->GetYaxis()->SetLabelSize(11);
  h_thetamup_temp_0->Draw();
  line= new TLine(0, 1, 3.14, 1);
  line->SetLineColor(kRed);
  line->SetLineStyle(9);
   line->Draw("same");
 
  c6->cd();


   if (cosmicCut){
     if(tune==3){c6->Print("figures/Tune3/BackSep/h_thetamup_allsel.png");}
     else{c6->Print("figures/Tune1/BackSep/h_thetamup_allsel.png");}
   }
   else{
     if(tune==3){c6->Print("figures_noCosmicCut/Tune3/BackSep/h_thetamup_allsel.png");}
     else{c6->Print("figures_noCosmicCut/Tune1/BackSep/h_thetamup_allsel.png");}
   }
  
  int nbins_thetamup=h_thetamup_allsel[2]->GetNbinsX();
  float chi2_thetamup=Chi2Calc(h_thetamup_allsel[2], h_thetamup_allsel[0], h_thetamup_allsel[1],scalefac, normfac);
  cout<<"chi2_thetamup = "<<chi2_thetamup<<"  number of bins: "<<nbins_thetamup<<"  chi2/dof is : "<<chi2_thetamup/nbins_thetamup<<endl;
 //======================================================================================
  //================================================================================
  c6->Update();

  TH1D *h_onoff_ptmis=(TH1D*)h_ptmis_allsel[1]->Clone(Form("%s_on-off", h_ptmis_allsel[1]->GetName()));
  TH1D *h_ptmis_MCandOff=(TH1D*)h_ptmis_allsel[2]->Clone(Form("%s_MC+Off", h_ptmis_allsel[2]->GetName()));
  cout<<"get all the histograms!"<<endl;
  
  pad1->SetBottomMargin(0); // Upper and lower plot are joined
  pad1->SetGridx();         // Vertical grid
  pad1->Draw();             // Draw the upper pad: pad1
  pad1->cd();               // pad1 becomes the current pad


  h_ptmis_allsel[0]->SetLineColor(kBlack);
  h_ptmis_allsel[0]->SetLineWidth(2);
  h_ptmis_allsel[0]->SetLineStyle(1);
  h_ptmis_allsel[0]->GetXaxis()->SetTitle("Missing-Pt[GeV/c]");
  h_ptmis_allsel[0]->GetYaxis()->SetTitle("Selected Events");
  h_ptmis_allsel[0]->GetYaxis()->SetTitleSize(0.06);
  h_ptmis_allsel[0]->GetYaxis()->SetTitleOffset(0.6);
  //h_ptmis_allsel[0]->SetMaximum(1.5*h_ptmis_allsel[0]->GetMaximum());
  //h_ptmis_allsel[0]->SetMaximum(1200);
  h_ptmis_allsel[0]->Draw();  

  THStack *hs_ptmis = new THStack("hs_ptmis","");
  stackHists(hs_ptmis, h_ptmis_sig, h_ptmis_bac, h_ptmis_allsel, normfac, scale_onoffbeam, scale_dirt_MC);
  hs_ptmis -> Draw("HIST,SAME");

  h_ptmis_allsel[0]->Scale(areanorm_fac);
  h_ptmis_allsel[0]->SetMaximum(600);
  //h_ptmis_allsel[0]->Draw("same");  
  //h_onoff_ptmis->Draw("E1CSAME");
  //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
 
  legendR->Draw("same"); prelim->Draw("same");

  //Draw ratio between two histograms and save the plot
  // lower plot will be in pad2
  c6->cd();          // Go back to the main canvas before defining pad2
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.2);
  pad2->SetGridx(); // vertical grid
  pad2->Draw();
  pad2->cd();       // pad2 becomes the current pad

  TH1D *h_ptmis_temp_0=(TH1D*)h_ptmis_allsel[0]->Clone("h_ptmis_temp_0"); 
  TH1D *h_ptmis_temp_1=(TH1D*)h_ptmis_allsel[1]->Clone("h_ptmis_temp_1"); 
  TH1D *h_ptmis_temp_2=(TH1D*)h_ptmis_allsel[2]->Clone("h_ptmis_temp_2"); 
  h_ptmis_temp_2->Scale(normfac);
  h_ptmis_MCandOff->Add(h_ptmis_temp_2,h_ptmis_temp_1,1,1); 

  h_ptmis_temp_0->Divide(h_ptmis_MCandOff);
  h_ptmis_temp_0->SetMinimum(0);
  h_ptmis_temp_0->SetMaximum(2);
  h_ptmis_temp_0->GetXaxis()->SetTitleSize(0.15);
  h_ptmis_temp_0->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  h_ptmis_temp_0->GetXaxis()->SetLabelSize(11);
 

  h_ptmis_temp_0->GetYaxis()->SetTitle("Data/MC");
  h_ptmis_temp_0->GetYaxis()->SetTitleSize(18);
  h_ptmis_temp_0->GetYaxis()->SetTitleFont(43);
  h_ptmis_temp_0->GetYaxis()->SetTitleOffset(0.9);
  h_ptmis_temp_0->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  h_ptmis_temp_0->GetYaxis()->SetLabelSize(11);
  h_ptmis_temp_0->Draw();
  line= new TLine(0, 1, 1.0, 1);
  line->SetLineColor(kRed);
  line->SetLineStyle(9);
   line->Draw("same");
 
  c6->cd();


   if (cosmicCut){
     if(tune==3){c6->Print("figures/Tune3/BackSep/h_ptmis_allsel.png");}
     else{c6->Print("figures/Tune1/BackSep/h_ptmis_allsel.png");}
   }
   else{
     if(tune==3){c6->Print("figures_noCosmicCut/Tune3/BackSep/h_ptmis_allsel.png");}
     else{c6->Print("figures_noCosmicCut/Tune1/BackSep/h_ptmis_allsel.png");}
   }
  
  int nbins_ptmis=h_ptmis_allsel[2]->GetNbinsX();
  float chi2_ptmis=Chi2Calc(h_ptmis_allsel[2], h_ptmis_allsel[0], h_ptmis_allsel[1],scalefac, normfac);
  cout<<"chi2_ptmis = "<<chi2_ptmis<<"  number of bins: "<<nbins_ptmis<<"  chi2/dof is : "<<chi2_ptmis/nbins_ptmis<<endl;

   //=============================================================================================
  c6->Update();

  TH1D *h_onoff_etatest=(TH1D*)h_etatest_allsel[1]->Clone(Form("%s_on-off", h_etatest_allsel[1]->GetName()));
  TH1D *h_etatest_MCandOff=(TH1D*)h_etatest_allsel[2]->Clone(Form("%s_MC+Off", h_etatest_allsel[2]->GetName()));
  cout<<"get all the histograms!"<<endl;
  
  pad1->SetBottomMargin(0); // Upper and lower plot are joined
  pad1->SetGridx();         // Vertical grid
  pad1->Draw();             // Draw the upper pad: pad1
  pad1->cd();               // pad1 becomes the current pad


  h_etatest_allsel[0]->SetLineColor(kBlack);
  h_etatest_allsel[0]->SetLineWidth(2);
  h_etatest_allsel[0]->SetLineStyle(1);
  h_etatest_allsel[0]->GetXaxis()->SetTitle("#eta");
  h_etatest_allsel[0]->GetYaxis()->SetTitle("Selected Events");
  h_etatest_allsel[0]->GetYaxis()->SetTitleSize(0.06);
  h_etatest_allsel[0]->GetYaxis()->SetTitleOffset(0.6);
  //h_etatest_allsel[0]->SetMaximum(1.5*h_etatest_allsel[0]->GetMaximum());
  //h_etatest_allsel[0]->SetMaximum(1200);
  h_etatest_allsel[0]->Draw();  

  THStack *hs_etatest = new THStack("hs_etatest","");
  stackHists(hs_etatest, h_etatest_sig, h_etatest_bac, h_etatest_allsel, normfac, scale_onoffbeam, scale_dirt_MC);
  hs_etatest -> Draw("HIST,SAME");

  h_etatest_allsel[0]->Scale(areanorm_fac);
  h_etatest_allsel[0]->SetMaximum(1400);
  //h_etatest_allsel[0]->Draw("same");  
  //h_onoff_etatest->Draw("E1CSAME");
  //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
 
  legendL->Draw("same"); prelim->Draw("same");

  //Draw ratio between two histograms and save the plot
  // lower plot will be in pad2
  c6->cd();          // Go back to the main canvas before defining pad2
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.2);
  pad2->SetGridx(); // vertical grid
  pad2->Draw();
  pad2->cd();       // pad2 becomes the current pad

  TH1D *h_etatest_temp_0=(TH1D*)h_etatest_allsel[0]->Clone("h_etatest_temp_0"); 
  TH1D *h_etatest_temp_1=(TH1D*)h_etatest_allsel[1]->Clone("h_etatest_temp_1"); 
  TH1D *h_etatest_temp_2=(TH1D*)h_etatest_allsel[2]->Clone("h_etatest_temp_2"); 
  h_etatest_temp_2->Scale(normfac);
  h_etatest_MCandOff->Add(h_etatest_temp_2,h_etatest_temp_1,1,1); 

  h_etatest_temp_0->Divide(h_etatest_MCandOff);
  h_etatest_temp_0->SetMinimum(0);
  h_etatest_temp_0->SetMaximum(2);
  h_etatest_temp_0->GetXaxis()->SetTitleSize(0.15);
  h_etatest_temp_0->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  h_etatest_temp_0->GetXaxis()->SetLabelSize(11);
 

  h_etatest_temp_0->GetYaxis()->SetTitle("Data/MC");
  h_etatest_temp_0->GetYaxis()->SetTitleSize(18);
  h_etatest_temp_0->GetYaxis()->SetTitleFont(43);
  h_etatest_temp_0->GetYaxis()->SetTitleOffset(0.9);
  h_etatest_temp_0->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  h_etatest_temp_0->GetYaxis()->SetLabelSize(11);
  h_etatest_temp_0->Draw();
  line= new TLine(0, 1, 2000, 1);
  line->SetLineColor(kRed);
  line->SetLineStyle(9);
  line->Draw("same");
 
  c6->cd();


   if (cosmicCut){
     if(tune==3){c6->Print("figures/Tune3/BackSep/h_etatest_allsel.png");}
     else{c6->Print("figures/Tune1/BackSep/h_etatest_allsel.png");}
   }
   else{
     if(tune==3){c6->Print("figures_noCosmicCut/Tune3/BackSep/h_etatest_allsel.png");}
     else{c6->Print("figures_noCosmicCut/Tune1/BackSep/h_etatest_allsel.png");}
   }
  
  int nbins_etatest=h_etatest_allsel[2]->GetNbinsX();
  float chi2_etatest=Chi2Calc(h_etatest_allsel[2], h_etatest_allsel[0], h_etatest_allsel[1],scalefac, normfac);
  cout<<"chi2_etatest = "<<chi2_etatest<<"  number of bins: "<<nbins_etatest<<"  chi2/dof is : "<<chi2_etatest/nbins_etatest<<endl;
  //=============================================================================================== 
  c6->Update();

  TH1D *h_onoff_enucal=(TH1D*)h_enucal_allsel[1]->Clone(Form("%s_on-off", h_enucal_allsel[1]->GetName()));
  TH1D *h_enucal_MCandOff=(TH1D*)h_enucal_allsel[2]->Clone(Form("%s_MC+Off", h_enucal_allsel[2]->GetName()));
  cout<<"get all the histograms!"<<endl;
  
  pad1->SetBottomMargin(0); // Upper and lower plot are joined
  pad1->SetGridx();         // Vertical grid
  pad1->Draw();             // Draw the upper pad: pad1
  pad1->cd();               // pad1 becomes the current pad


  h_enucal_allsel[0]->SetLineColor(kBlack);
  h_enucal_allsel[0]->SetLineWidth(2);
  h_enucal_allsel[0]->SetLineStyle(1);
  h_enucal_allsel[0]->GetXaxis()->SetTitle("Calculated Neutrino Energy[MeV]");
  h_enucal_allsel[0]->GetYaxis()->SetTitle("Selected Events");
  h_enucal_allsel[0]->GetYaxis()->SetTitleSize(0.06);
  h_enucal_allsel[0]->GetYaxis()->SetTitleOffset(0.6);
  //h_enucal_allsel[0]->SetMaximum(1.5*h_enucal_allsel[0]->GetMaximum());
  //h_enucal_allsel[0]->SetMaximum(1200);
  h_enucal_allsel[0]->Draw();  

  THStack *hs_enucal = new THStack("hs_enucal","");
  stackHists(hs_enucal, h_enucal_sig, h_enucal_bac, h_enucal_allsel, normfac, scale_onoffbeam, scale_dirt_MC);
  hs_enucal -> Draw("HIST,SAME");

  h_enucal_allsel[0]->Scale(areanorm_fac);
  h_enucal_allsel[0]->SetMaximum(600);
  //h_enucal_allsel[0]->Draw("same");  
  //h_onoff_enucal->Draw("E1CSAME");
  //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
 
  legendR->Draw("same"); prelim->Draw("same");

  //Draw ratio between two histograms and save the plot
  // lower plot will be in pad2
  c6->cd();          // Go back to the main canvas before defining pad2
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.2);
  pad2->SetGridx(); // vertical grid
  pad2->Draw();
  pad2->cd();       // pad2 becomes the current pad

  TH1D *h_enucal_temp_0=(TH1D*)h_enucal_allsel[0]->Clone("h_enucal_temp_0"); 
  TH1D *h_enucal_temp_1=(TH1D*)h_enucal_allsel[1]->Clone("h_enucal_temp_1"); 
  TH1D *h_enucal_temp_2=(TH1D*)h_enucal_allsel[2]->Clone("h_enucal_temp_2"); 
  h_enucal_temp_2->Scale(normfac);
  h_enucal_MCandOff->Add(h_enucal_temp_2,h_enucal_temp_1,1,1); 

  h_enucal_temp_0->Divide(h_enucal_MCandOff);
  h_enucal_temp_0->SetMinimum(0);
  h_enucal_temp_0->SetMaximum(2);
  h_enucal_temp_0->GetXaxis()->SetTitleSize(0.15);
  h_enucal_temp_0->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  h_enucal_temp_0->GetXaxis()->SetLabelSize(11);
 

  h_enucal_temp_0->GetYaxis()->SetTitle("Data/MC");
  h_enucal_temp_0->GetYaxis()->SetTitleSize(18);
  h_enucal_temp_0->GetYaxis()->SetTitleFont(43);
  h_enucal_temp_0->GetYaxis()->SetTitleOffset(0.9);
  h_enucal_temp_0->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  h_enucal_temp_0->GetYaxis()->SetLabelSize(11);
  h_enucal_temp_0->Draw();
  line= new TLine(0, 1, 2000, 1);
  line->SetLineColor(kRed);
  line->SetLineStyle(9);
  line->Draw("same");
 
  c6->cd();


   if (cosmicCut){
     if(tune==3){c6->Print("figures/Tune3/BackSep/h_enucal_allsel.png");}
     else{c6->Print("figures/Tune1/BackSep/h_enucal_allsel.png");}
   }
   else{
     if(tune==3){c6->Print("figures_noCosmicCut/Tune3/BackSep/h_enucal_allsel.png");}
     else{c6->Print("figures_noCosmicCut/Tune1/BackSep/h_enucal_allsel.png");}
   }
  
  int nbins_enucal=h_enucal_allsel[2]->GetNbinsX();
  float chi2_enucal=Chi2Calc(h_enucal_allsel[2], h_enucal_allsel[0], h_enucal_allsel[1],scalefac, normfac);
  cout<<"chi2_enucal = "<<chi2_enucal<<"  number of bins: "<<nbins_enucal<<"  chi2/dof is : "<<chi2_enucal/nbins_enucal<<endl;
 //=======================================================================
  TH1D *h_onoff_pmult=(TH1D*)h_pmult_allsel[1]->Clone(Form("%s_on-off", h_pmult_allsel[1]->GetName()));
  TH1D *h_pmult_MCandOff=(TH1D*)h_pmult_allsel[2]->Clone(Form("%s_MC+Off", h_pmult_allsel[2]->GetName()));
  cout<<"get all the histograms!"<<endl;

  
  TCanvas *c7 = new TCanvas("c7"); 

  pad1->SetBottomMargin(0); // Upper and lower plot are joined
  pad1->SetGridx();         // Vertical grid
  pad1->Draw();             // Draw the upper pad: pad1
  pad1->cd();               // pad1 becomes the current pad


  h_pmult_allsel[0]->SetLineColor(kBlack);
  h_pmult_allsel[0]->SetLineWidth(2);
  h_pmult_allsel[0]->SetLineStyle(1);
  h_pmult_allsel[0]->GetXaxis()->SetTitle("Proton Multiplicity");
  h_pmult_allsel[0]->GetYaxis()->SetTitle("Selected Events");
  h_pmult_allsel[0]->GetYaxis()->SetTitleSize(0.06);
  h_pmult_allsel[0]->GetYaxis()->SetTitleOffset(0.6);
  //h_pmult_allsel[0]->SetMaximum(1.5*h_pmult_allsel[0]->GetMaximum());
  //h_pmult_allsel[0]->SetMaximum(1200);
  h_pmult_allsel[0]->Draw();  

  THStack *hs_pmult = new THStack("hs_pmult","");
  stackHists(hs_pmult, h_pmult_sig, h_pmult_bac, h_pmult_allsel, normfac, scale_onoffbeam, scale_dirt_MC);
  hs_pmult -> Draw("HIST,SAME");

  h_pmult_allsel[0]->Scale(areanorm_fac);
  h_pmult_allsel[0]->SetMaximum(6000);
  h_pmult_allsel[0]->Draw("same");  
  //h_onoff_pmult->Draw("E1CSAME");
  //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
 
  legendR->Draw("same"); prelim->Draw("same");

  //Draw ratio between two histograms and save the plot
  // lower plot will be in pad2
  c7->cd();          // Go back to the main canvas before defining pad2
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.2);
  pad2->SetGridx(); // vertical grid
  pad2->Draw();
  pad2->cd();       // pad2 becomes the current pad

  TH1D *h_pmult_temp_0=(TH1D*)h_pmult_allsel[0]->Clone("h_pmult_temp_0"); 
  TH1D *h_pmult_temp_1=(TH1D*)h_pmult_allsel[1]->Clone("h_pmult_temp_1"); 
  TH1D *h_pmult_temp_2=(TH1D*)h_pmult_allsel[2]->Clone("h_pmult_temp_2"); 
  h_pmult_temp_2->Scale(normfac);
  h_pmult_MCandOff->Add(h_pmult_temp_2,h_pmult_temp_1,1,1); 

  h_pmult_temp_0->Divide(h_pmult_MCandOff);
  h_pmult_temp_0->SetMinimum(0);
  h_pmult_temp_0->SetMaximum(2);
  h_pmult_temp_0->GetXaxis()->SetTitleSize(0.15);
  h_pmult_temp_0->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  h_pmult_temp_0->GetXaxis()->SetLabelSize(11);
 

  h_pmult_temp_0->GetYaxis()->SetTitle("Data/MC");
  h_pmult_temp_0->GetYaxis()->SetTitleSize(18);
  h_pmult_temp_0->GetYaxis()->SetTitleFont(43);
  h_pmult_temp_0->GetYaxis()->SetTitleOffset(0.9);
  h_pmult_temp_0->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  h_pmult_temp_0->GetYaxis()->SetLabelSize(11);
  h_pmult_temp_0->Draw();
  line= new TLine(0, 1, 9.5, 1);
  line->SetLineColor(kRed);
  line->SetLineStyle(9);
   line->Draw("same");
 
  c7->cd();


   if (cosmicCut){
     if(tune==3){c7->Print("figures/Tune3/BackSep/h_pmult_allsel.png");}
     else{c7->Print("figures/Tune1/BackSep/h_pmult_allsel.png");}
   }
   else{
     if(tune==3){c7->Print("figures_noCosmicCut/Tune3/BackSep/h_pmult_allsel.png");}
     else{c7->Print("figures_noCosmicCut/Tune1/BackSep/h_pmult_allsel.png");}
   }
  
  int nbins_pmult=h_pmult_allsel[2]->GetNbinsX();
  float chi2_pmult=Chi2Calc(h_pmult_allsel[2], h_pmult_allsel[0], h_pmult_allsel[1],scalefac, normfac);
  cout<<"chi2_pmult = "<<chi2_pmult<<"  number of bins: "<<nbins_pmult<<"  chi2/dof is : "<<chi2_pmult/nbins_pmult<<endl;
  //========================================================================
  c7->Update();

  TH1D *h_onoff_alphat=(TH1D*)h_alphat_allsel[1]->Clone(Form("%s_on-off", h_alphat_allsel[1]->GetName()));
  TH1D *h_alphat_MCandOff=(TH1D*)h_alphat_allsel[2]->Clone(Form("%s_MC+Off", h_alphat_allsel[2]->GetName()));
  cout<<"get all the histograms!"<<endl;
  
  pad1->SetBottomMargin(0); // Upper and lower plot are joined
  pad1->SetGridx();         // Vertical grid
  pad1->Draw();             // Draw the upper pad: pad1
  pad1->cd();               // pad1 becomes the current pad


  h_alphat_allsel[0]->SetLineColor(kBlack);
  h_alphat_allsel[0]->SetLineWidth(2);
  h_alphat_allsel[0]->SetLineStyle(1);
  h_alphat_allsel[0]->GetXaxis()->SetTitle("#delta#alpha_{T}");
  h_alphat_allsel[0]->GetYaxis()->SetTitle("Selected Events");
  h_alphat_allsel[0]->GetYaxis()->SetTitleSize(0.06);
  h_alphat_allsel[0]->GetYaxis()->SetTitleOffset(0.6);
  //h_alphat_allsel[0]->SetMaximum(1.5*h_alphat_allsel[0]->GetMaximum());
  //h_alphat_allsel[0]->SetMaximum(1200);
  h_alphat_allsel[0]->Draw();  

  THStack *hs_alphat = new THStack("hs_alphat","");
  stackHists(hs_alphat, h_alphat_sig, h_alphat_bac, h_alphat_allsel, normfac, scale_onoffbeam, scale_dirt_MC);
  hs_alphat -> Draw("HIST,SAME");

  h_alphat_allsel[0]->Scale(areanorm_fac);
  h_alphat_allsel[0]->SetMaximum(1800);
  //h_alphat_allsel[0]->Draw("same");  
  //h_onoff_alphat->Draw("E1CSAME");
  //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
 
  legendR->Draw("same"); prelim->Draw("same");

  //Draw ratio between two histograms and save the plot
  // lower plot will be in pad2
  c7->cd();          // Go back to the main canvas before defining pad2
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.2);
  pad2->SetGridx(); // vertical grid
  pad2->Draw();
  pad2->cd();       // pad2 becomes the current pad

  TH1D *h_alphat_temp_0=(TH1D*)h_alphat_allsel[0]->Clone("h_alphat_temp_0"); 
  TH1D *h_alphat_temp_1=(TH1D*)h_alphat_allsel[1]->Clone("h_alphat_temp_1"); 
  TH1D *h_alphat_temp_2=(TH1D*)h_alphat_allsel[2]->Clone("h_alphat_temp_2"); 
  h_alphat_temp_2->Scale(normfac);
  h_alphat_MCandOff->Add(h_alphat_temp_2,h_alphat_temp_1,1,1); 

  h_alphat_temp_0->Divide(h_alphat_MCandOff);
  h_alphat_temp_0->SetMinimum(0);
  h_alphat_temp_0->SetMaximum(2);
  h_alphat_temp_0->GetXaxis()->SetTitleSize(0.15);
  h_alphat_temp_0->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  h_alphat_temp_0->GetXaxis()->SetLabelSize(11);
 

  h_alphat_temp_0->GetYaxis()->SetTitle("Data/MC");
  h_alphat_temp_0->GetYaxis()->SetTitleSize(18);
  h_alphat_temp_0->GetYaxis()->SetTitleFont(43);
  h_alphat_temp_0->GetYaxis()->SetTitleOffset(0.9);
  h_alphat_temp_0->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  h_alphat_temp_0->GetYaxis()->SetLabelSize(11);
  h_alphat_temp_0->Draw();
  line= new TLine(0, 1, 3.14, 1);
  line->SetLineColor(kRed);
  line->SetLineStyle(9);
   line->Draw("same");
 
  c7->cd();


   if (cosmicCut){
     if(tune==3){c7->Print("figures/Tune3/BackSep/h_alphat_allsel.png");}
     else{c7->Print("figures/Tune1/BackSep/h_alphat_allsel.png");}
   }
   else{
     if(tune==3){c7->Print("figures_noCosmicCut/Tune3/BackSep/h_alphat_allsel.png");}
     else{c7->Print("figures_noCosmicCut/Tune1/BackSep/h_alphat_allsel.png");}
   }
  
  int nbins_alphat=h_alphat_allsel[2]->GetNbinsX();
  float chi2_alphat=Chi2Calc(h_alphat_allsel[2], h_alphat_allsel[0], h_alphat_allsel[1],scalefac, normfac);
  cout<<"chi2_alphat = "<<chi2_alphat<<"  number of bins: "<<nbins_alphat<<"  chi2/dof is : "<<chi2_alphat/nbins_alphat<<endl;
  
  
  //==========================================================================
  // Dirt plots
  //==========================================================================
  
  // Dirt-only distributions
  c1->cd();
  c1->Clear();
  h_range_allsel[3]->Draw();
  c1->Print("figures/Dirt/h_range_dirt.png");
  
  c1->cd();
  c1->Clear();
  h_prange_allsel[3]->Draw();
  c1->Print("figures/Dirt/h_prange_dirt.png");
  
  c1->cd();
  c1->Clear();
  h_phi_allsel[3]->Draw();
  c1->Print("figures/Dirt/h_phi_dirt.png");
  
  c1->cd();
  c1->Clear();
  h_pphi_allsel[3]->Draw();
  c1->Print("figures/Dirt/h_pphi_dirt.png");
  
  c1->cd();
  c1->Clear();
  h_costheta_allsel[3]->Draw();
  c1->Print("figures/Dirt/h_costheta_dirt.png");
  
  c1->cd();
  c1->Clear();
  h_pcostheta_allsel[3]->Draw();
  c1->Print("figures/Dirt/h_pcostheta_dirt.png");
  
  c1->cd();
  c1->Clear();
  h_plep_allsel[3]->Draw();
  c1->Print("figures/Dirt/h_plep_dirt.png");
  
  c1->cd();
  c1->Clear();
  h_phad_allsel[3]->Draw();
  c1->Print("figures/Dirt/h_phad_dirt.png");
  
  c1->cd();
  c1->Clear();
  h_thetamup_allsel[3]->Draw();
  c1->Print("figures/Dirt/h_thetamup_dirt.png");
  
  c1->cd();
  c1->Clear();
  h_ptmis_allsel[3]->Draw();
  c1->Print("figures/Dirt/h_ptmis_dirt.png");
  
  c1->cd();
  c1->Clear();
  h_etatest_allsel[3]->Draw();
  c1->Print("figures/Dirt/h_etatest_dirt.png");
  
  c1->cd();
  c1->Clear();
  h_enucal_allsel[3]->Draw();
  c1->Print("figures/Dirt/h_enucal_dirt.png");
  
  c1->cd();
  c1->Clear();
  h_alphat_allsel[3]->Draw();
  c1->Print("figures/Dirt/h_alphat_dirt.png");
  
  c1->cd();
  c1->Clear();
  h_pmult_allsel[3]->Draw();
  c1->Print("figures/Dirt/h_pmult_dirt.png");
  
  // Dirt as a fraction of MC stack (with off-beam)
  TH1 *data_cl;
  TH1 *dirt_cl;

  // range
  c1->cd();
  c1->Clear();
  data_cl = (TH1*)h_range_allsel[0]->Clone();
  dirt_cl = (TH1*)h_range_allsel[3]->Clone();
  dirt_cl->Divide(data_cl);
  c1->Clear();
  dirt_cl->Draw();
  c1->Print("figures/Dirt/h_range_dirtfrac.png");

  // prange
  c1->cd();
  c1->Clear();
  data_cl = (TH1*)h_prange_allsel[0]->Clone();
  dirt_cl = (TH1*)h_prange_allsel[3]->Clone();
  dirt_cl->Divide(data_cl);
  c1->Clear();
  dirt_cl->Draw();
  c1->Print("figures/Dirt/h_prange_dirtfrac.png");

  // phi
  c1->cd();
  c1->Clear();
  data_cl = (TH1*)h_phi_allsel[0]->Clone();
  dirt_cl = (TH1*)h_phi_allsel[3]->Clone();
  dirt_cl->Divide(data_cl);
  c1->Clear();
  dirt_cl->Draw();
  c1->Print("figures/Dirt/h_phi_dirtfrac.png");

  // pphi
  c1->cd();
  c1->Clear();
  data_cl = (TH1*)h_pphi_allsel[0]->Clone();
  dirt_cl = (TH1*)h_pphi_allsel[3]->Clone();
  dirt_cl->Divide(data_cl);
  c1->Clear();
  dirt_cl->Draw();
  c1->Print("figures/Dirt/h_pphi_dirtfrac.png");

  // costheta
  c1->cd();
  c1->Clear();
  data_cl = (TH1*)h_costheta_allsel[0]->Clone();
  dirt_cl = (TH1*)h_costheta_allsel[3]->Clone();
  dirt_cl->Divide(data_cl);
  c1->Clear();
  dirt_cl->Draw();
  c1->Print("figures/Dirt/h_costheta_dirtfrac.png");

  // pcostheta
  c1->cd();
  c1->Clear();
  data_cl = (TH1*)h_pcostheta_allsel[0]->Clone();
  dirt_cl = (TH1*)h_pcostheta_allsel[3]->Clone();
  dirt_cl->Divide(data_cl);
  c1->Clear();
  dirt_cl->Draw();
  c1->Print("figures/Dirt/h_pcostheta_dirtfrac.png");

  // plep
  c1->cd();
  c1->Clear();
  data_cl = (TH1*)h_plep_allsel[0]->Clone();
  dirt_cl = (TH1*)h_plep_allsel[3]->Clone();
  dirt_cl->Divide(data_cl);
  c1->Clear();
  dirt_cl->Draw();
  c1->Print("figures/Dirt/h_plep_dirtfrac.png");

  // phad
  c1->cd();
  c1->Clear();
  data_cl = (TH1*)h_phad_allsel[0]->Clone();
  dirt_cl = (TH1*)h_phad_allsel[3]->Clone();
  dirt_cl->Divide(data_cl);
  c1->Clear();
  dirt_cl->Draw();
  c1->Print("figures/Dirt/h_phad_dirtfrac.png");

  // ptmis
  c1->cd();
  c1->Clear();
  data_cl = (TH1*)h_ptmis_allsel[0]->Clone();
  dirt_cl = (TH1*)h_ptmis_allsel[3]->Clone();
  dirt_cl->Divide(data_cl);
  c1->Clear();
  dirt_cl->Draw();
  c1->Print("figures/Dirt/h_ptmis_dirtfrac.png");

  // etatest
  c1->cd();
  c1->Clear();
  data_cl = (TH1*)h_etatest_allsel[0]->Clone();
  dirt_cl = (TH1*)h_etatest_allsel[3]->Clone();
  dirt_cl->Divide(data_cl);
  c1->Clear();
  dirt_cl->Draw();
  c1->Print("figures/Dirt/h_etatest_dirtfrac.png");

  // enucal
  c1->cd();
  c1->Clear();
  data_cl = (TH1*)h_enucal_allsel[0]->Clone();
  dirt_cl = (TH1*)h_enucal_allsel[3]->Clone();
  dirt_cl->Divide(data_cl);
  c1->Clear();
  dirt_cl->Draw();
  c1->Print("figures/Dirt/h_enucal_dirtfrac.png");

  // alphat
  c1->cd();
  c1->Clear();
  data_cl = (TH1*)h_alphat_allsel[0]->Clone();
  dirt_cl = (TH1*)h_alphat_allsel[3]->Clone();
  dirt_cl->Divide(data_cl);
  c1->Clear();
  dirt_cl->Draw();
  c1->Print("figures/Dirt/h_alphat_dirtfrac.png");

  // pmult
  c1->cd();
  c1->Clear();
  data_cl = (TH1*)h_pmult_allsel[0]->Clone();
  dirt_cl = (TH1*)h_pmult_allsel[3]->Clone();
  dirt_cl->Divide(data_cl);
  c1->Clear();
  dirt_cl->Draw();
  c1->Print("figures/Dirt/h_pmult_dirtfrac.png");

 }
