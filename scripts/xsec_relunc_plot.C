void xsec_relunc_plot(){

    std::cout<<"beginning of the analysis"<<std::endl;
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

  
  TFile *output_file=new TFile("output.root", "recreate");
  std::cout<<"created output file "<<std::endl;
  //open the xsec root files and get the central value of the cross sections
  TFile *xsec_inputfile= new TFile("/pnfs/uboone/persistent/users/jiangl/CC1uNPSelection/xsec_covmatrix_withadditionalcuts/xsec_file_CV_Tune1.root");
  //TFile *xsec_inputfile= new TFile("xsec_file_CV.root");
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



  TLegend *leg=new TLegend(0.65, 0.7, 0.9, 0.9);
  leg->AddEntry(xsec_mumom_mc, "GENIE-default");
  leg->AddEntry(xsec_mumom_CV, "Measured");
  


  std::cout<<"libo test 1"<<std::endl;

  std::cout<<"start set the detector systematic list "<<std::endl;
  TFile* inputfile1;
  TFile* inputfile2;

  double bins_mumom[7] = {0.1, 0.18, 0.30, 0.48, 0.75, 1.14, 2.50};
  double bins_mucostheta[13] = {-1.00, -0.82, -0.66, -0.39, -0.16, 0.05, 0.25, 0.43, 0.59, 0.73, 0.83, 0.91, 1.00};
  int n_bins_mumom = 6;
  int n_bins_mucostheta = 12;

  double bins_testmu[8]={0.0, 0.1, 0.18, 0.30, 0.48, 0.75, 1.14, 2.50};
  double bins_testp[12]={0.0, 0.30, 0.41, 0.495, 0.56, 0.62, 0.68, 0.74, 0.80, 0.87, 0.93, 1.20};
  int n_bins_testmu = 7;
  int n_bins_testp = 11;


  double bins_pmom[11] = {0.30, 0.41, 0.495, 0.56, 0.62, 0.68, 0.74, 0.80, 0.87, 0.93, 1.20};
  double bins_pcostheta[10] = {-1.00, -0.50, 0.00, 0.27, 0.45, 0.62, 0.76, 0.86, 0.94, 1.00};
  int n_bins_pmom = 10; 
  int n_bins_pcostheta = 9;

  double bins_muptheta[7] = {0.00, 0.8, 1.2, 1.57, 1.94, 2.34, 3.14};
  int n_bins_muptheta = 6;
 
int sys_sel=0;
string syst_unc_name="";
std::vector<std::string> syst_unc_list;

bool flux_syst = false;
bool genie_syst = false;
bool extra_syst = false;
bool total_flux_syst = false;
bool total_extra_syst = false;
bool det_syst = false;
bool total_det_syst = false;
bool tune3_syst = false;

if(sys_sel==0){
genie_syst=true;
syst_unc_name="genie";
syst_unc_list.push_back("genie");
}
if(sys_sel==1){
total_flux_syst=true;
syst_unc_name="total_flux";
syst_unc_list.push_back("flux");
}
if(sys_sel==2){
total_extra_syst=true;
syst_unc_name="total_extra";
syst_unc_list.push_back("extra_syst");
}

if(sys_sel==3){
flux_syst=true;
syst_unc_name="flux";
syst_unc_list.push_back("flux_FluxUnisim");
syst_unc_list.push_back("flux_kminus");
syst_unc_list.push_back("flux_kplus");
syst_unc_list.push_back("flux_kzero");
syst_unc_list.push_back("flux_piminus");
syst_unc_list.push_back("flux_piplus");
}
if(sys_sel==4){
extra_syst=true;
syst_unc_name="extra_syst:CCMEC, QE";
syst_unc_list.push_back("extra_syst_ccmec");
syst_unc_list.push_back("extra_syst_ccqe");
//syst_unc_list.push_back("extra_syst_reint_proton");
//syst_unc_list.push_back("extra_syst_reint_piplus");
//syst_unc_list.push_back("extra_syst_reint_piminus");
}

if(sys_sel==5){
extra_syst=true;
syst_unc_name="extra_syst:Reinteraction";
//syst_unc_list.push_back("extra_syst_ccmec");
//syst_unc_list.push_back("extra_syst_ccqe");
syst_unc_list.push_back("extra_syst_reint_proton");
syst_unc_list.push_back("extra_syst_reint_piplus");
syst_unc_list.push_back("extra_syst_reint_piminus");
}



if(sys_sel==6){
  det_syst=true;
  syst_unc_name="det_syst";
  //syst_unc_list.push_back("DLdown"); //smaller than DLup so turning it off
  syst_unc_list.push_back("DLup");
  //syst_unc_list.push_back("DTup"); //smaller than DTdown so turning it off
  syst_unc_list.push_back("DTdown");
  syst_unc_list.push_back("noiseAmpUp");
  //syst_unc_list.push_back("noiseAmpDown"); //smaller than noiseAmpUp so turning it off
  //syst_unc_list.push_back("squeezeResp"); //smaller than stretchResp so turning it off
  syst_unc_list.push_back("stretchResp");
  //syst_unc_list.push_back("downPEnoise"); //smaller than upPEnoise so turning it off
  syst_unc_list.push_back("upPEnoise");
  syst_unc_list.push_back("LArG4BugFix");
  syst_unc_list.push_back("altDeadChannels");
  syst_unc_list.push_back("dataSCE");
  syst_unc_list.push_back("enhancedexttpcvis");
  syst_unc_list.push_back("lifetime10ms");
  syst_unc_list.push_back("withDIC");
  syst_unc_list.push_back("recombination");
 }

 if (sys_sel==7) {
   tune3_syst=true;
   syst_unc_name="Tune3";
   syst_unc_list.push_back("Tune3");
 }
 if (sys_sel==8) {
   syst_unc_name="GENIE,Flux,Reinteraction";
   syst_unc_list.push_back("genie");
   syst_unc_list.push_back("flux");
   syst_unc_list.push_back("extra_syst");
 }

auto mg=new TMultiGraph;

TH1D *gr_testmu;
TH1D *gr_testp;

TH1D *total_mumom_1;
TH1D *gr_mumom_1[20];
//TH1D *gr_mumom_2[20];
TH1D *total_pmom_1;
TH1D *gr_pmom_1[20];
//TH1D *gr_pmom_2[20];
TH1D *total_muangle_1;
TH1D *gr_muangle_1[20];
//TH1D *gr_muangle_2[20];
TH1D *total_pangle_1;
TH1D *gr_pangle_1[20];
//TH1D *gr_pangle_2[20];
TH1D *total_thetamup_1;
TH1D *gr_thetamup_1[20];
//TH1D *gr_thetamup_2[20];
 int dd = 9999;

 total_mumom_1=new TH1D(Form("Total_Muon_Momentum_%d",dd), Form("Total_Muon_Momentum_%d",dd), n_bins_mumom, bins_mumom);
 total_mumom_1->GetXaxis()->SetTitle("P_{#mu}[GeV]");
 total_mumom_1->GetYaxis()->SetTitle("Relative Uncertainty");
 total_mumom_1->SetLineColor(2);
 total_mumom_1->SetLineWidth(4);
 total_mumom_1->SetFillStyle(0);
 total_mumom_1->SetMinimum(0); 
 total_mumom_1->SetMaximum(1.0);
 
 total_pmom_1=new TH1D(Form("Total_Proton_Momentum_%d",dd), Form("Total_Proton_Momentum_%d",dd), n_bins_pmom, bins_pmom);
 total_pmom_1->GetXaxis()->SetTitle("P_{proton}[GeV]");
 total_pmom_1->GetYaxis()->SetTitle("Relative Uncertainty");
 total_pmom_1->SetLineColor(2);
 total_pmom_1->SetLineWidth(4);
 total_pmom_1->SetFillStyle(0);
 total_pmom_1->SetMinimum(0); 
 total_pmom_1->SetMaximum(1.0);
 
 total_muangle_1=new TH1D(Form("Total_Muon_CosTheta_%d",dd), Form("Total_Muon_CosTheta_%d",dd), n_bins_mucostheta, bins_mucostheta);
 total_muangle_1->GetXaxis()->SetTitle("Cos#theta_{#mu}");
 total_muangle_1->GetYaxis()->SetTitle("Relative Uncertainty");
 total_muangle_1->SetLineColor(2);
 total_muangle_1->SetLineWidth(4);
 total_muangle_1->SetFillStyle(0);
 total_muangle_1->SetMinimum(0); 
 total_muangle_1->SetMaximum(1.0);
 
 total_pangle_1=new TH1D(Form("Total_Proton_CosTheta_%d",dd), Form("Total_Proton_CosTheta_%d",dd), n_bins_pcostheta, bins_pcostheta);
 total_pangle_1->GetXaxis()->SetTitle("Cos#theta_{proton}");
 total_pangle_1->GetYaxis()->SetTitle("Relative Uncertainty");
 total_pangle_1->SetLineColor(2);
 total_pangle_1->SetLineWidth(4);
 total_pangle_1->SetFillStyle(0);
 total_pangle_1->SetMinimum(0); 
 total_pangle_1->SetMaximum(1.0);
 
 total_thetamup_1=new TH1D(Form("Total_Thetamup_%d",dd), Form("Total_Thetamup_%d",dd), n_bins_muptheta, bins_muptheta); 
 total_thetamup_1->GetXaxis()->SetTitle("#theta_{#mu, proton}");
 total_thetamup_1->GetYaxis()->SetTitle("Relative Uncertainty");
 total_thetamup_1->SetLineColor(2);
 total_thetamup_1->SetLineWidth(4);
 total_thetamup_1->SetFillStyle(0);
 total_thetamup_1->SetMinimum(0); 
 total_thetamup_1->SetMaximum(1.0);

gr_testmu= new TH1D("gr_testmu", "gr_testmu", n_bins_testmu, bins_testmu);
gr_testmu->GetXaxis()->SetTitle("P_{#mu}[GeV]");
gr_testmu->GetYaxis()->SetTitle("Relative Uncertainty");
gr_testp= new TH1D("gr_testp", "gr_testp", n_bins_testp, bins_testp);
gr_testp->GetXaxis()->SetTitle("P_{proton}[GeV]");
gr_testp->GetYaxis()->SetTitle("Relative Uncertainty");

for(int tt=0; tt<n_bins_testmu; tt++){
  gr_testmu->SetBinContent(0, tt);
}

for(int tt=0; tt<n_bins_testp; tt++){
  gr_testp->SetBinContent(0, tt);
}


for(int k=0; k<syst_unc_list.size(); k++){
   gr_mumom_1[k]=new TH1D(Form("Muon_Momentum_%d",k), Form("Muon_Momentum_%d",k), n_bins_mumom, bins_mumom); 
   gr_mumom_1[k]->GetXaxis()->SetTitle("P_{#mu}[GeV]");
   gr_mumom_1[k]->GetYaxis()->SetTitle("Relative Uncertainty");

   //gr_mumom_2[k]=new TH1D(Form("Muon_Momentum_showerastrack_%d",k), Form("Muon_Momentum_showerastrack_%d",k), n_bins_mumom, bins_mumom);
   //gr_mumom_2[k]->GetXaxis()->SetTitle("P_{#mu}[GeV]");
   //gr_mumom_2[k]->GetYaxis()->SetTitle("Relative Uncertainty");

   gr_pmom_1[k]=new TH1D(Form("Proton_Momentum_%d",k), Form("Proton_Momentum_%d",k), n_bins_pmom, bins_pmom); 
   gr_pmom_1[k]->GetXaxis()->SetTitle("P_{proton}[GeV]");
   gr_pmom_1[k]->GetYaxis()->SetTitle("Relative Uncertainty");

   //gr_pmom_2[k]=new TH1D(Form("Proton_Momentum_showerastrack_%d",k), Form("Proton_Momentum_showerastrack_%d",k), n_bins_pmom, bins_pmom);
   //gr_pmom_2[k]->GetXaxis()->SetTitle("P_{proton}[GeV]");
   //gr_pmom_2[k]->GetYaxis()->SetTitle("Relative Uncertainty");

   gr_muangle_1[k]=new TH1D(Form("Muon_CosTheta_%d",k), Form("Muon_CosTheta_%d",k), n_bins_mucostheta, bins_mucostheta); 
   gr_muangle_1[k]->GetXaxis()->SetTitle("Cos#theta_{#mu}");
   gr_muangle_1[k]->GetYaxis()->SetTitle("Relative Uncertainty");

   //gr_muangle_2[k]=new TH1D(Form("Muon_CosTheta_showerastrack_%d",k), Form("Muon_CosTheta_showerastrack_%d",k), n_bins_mucostheta, bins_mucostheta);
   //gr_muangle_2[k]->GetXaxis()->SetTitle("Cos#theta_{#mu}");
   //gr_muangle_2[k]->GetYaxis()->SetTitle("Relative Uncertainty");

   gr_pangle_1[k]=new TH1D(Form("Proton_CosTheta_%d",k), Form("Proton_CosTheta_%d",k), n_bins_pcostheta, bins_pcostheta); 
   gr_pangle_1[k]->GetXaxis()->SetTitle("Cos#theta_{proton}");
   gr_pangle_1[k]->GetYaxis()->SetTitle("Relative Uncertainty");

   //gr_pangle_2[k]=new TH1D(Form("Proton_CosTheta_showerastrack_%d",k), Form("Proton_CosTheta_showerastrack_%d",k), n_bins_pcostheta, bins_pcostheta);
   //gr_pangle_2[k]->GetXaxis()->SetTitle("Cos#theta_{proton}");
   //gr_pangle_2[k]->GetYaxis()->SetTitle("Relative Uncertainty");

   gr_thetamup_1[k]=new TH1D(Form("Thetamup_%d",k), Form("Thetamup_%d",k), n_bins_muptheta, bins_muptheta); 
   gr_thetamup_1[k]->GetXaxis()->SetTitle("#theta_{#mu, proton}");
   gr_thetamup_1[k]->GetYaxis()->SetTitle("Relative Uncertainty");

   //gr_thetamup_2[k]=new TH1D(Form("Thetamup_showerastrack_%d",k), Form("Thetamup_showerastrack_%d",k), n_bins_muptheta, bins_muptheta);
   //gr_thetamup_2[k]->GetXaxis()->SetTitle("#theta_{#mu, proton}");
   //gr_thetamup_2[k]->GetYaxis()->SetTitle("Relative Uncertainty");
}
for(int j=0; j<syst_unc_list.size(); j++){

   string inputfilename1;
   //string inputfilename2;
   if(sys_sel==0 || sys_sel ==1 || sys_sel==2) {
   inputfilename1="/pnfs/uboone/persistent/users/jiangl/CC1uNPSelection/xsec_covmatrix_withadditionalcuts/covariance_"+syst_unc_list[j]+".root";
   //inputfilename2="/pnfs/uboone/persistent/users/jiangl/CC1uNPSelection/xsec_covmatrix_withadditionalcuts/covariance_"+syst_unc_list[j]+".root";
   }

   if(sys_sel==3){
   inputfilename1="/pnfs/uboone/persistent/users/jiangl/CC1uNPSelection/xsec_covmatrix_withadditionalcuts/covariance_"+syst_unc_list[j]+".root";
   //inputfilename2="/pnfs/uboone/persistent/users/jiangl/CC1uNPSelection/showerastrack_folder/covariance_"+syst_unc_list[j]+".root";
   }
   else if(sys_sel==4){
   inputfilename1="/pnfs/uboone/persistent/users/jiangl/CC1uNPSelection/xsec_covmatrix_withadditionalcuts/covariance_"+syst_unc_list[j]+".root";
   //inputfilename2="/pnfs/uboone/persistent/users/jiangl/CC1uNPSelection/extra_showerbuiltastrack/covariance_"+syst_unc_list[j]+".root";
   }
   else if(sys_sel==5){
   inputfilename1="/pnfs/uboone/persistent/users/jiangl/CC1uNPSelection/xsec_covmatrix_withadditionalcuts/covariance_"+syst_unc_list[j]+".root";
   //inputfilename2="/pnfs/uboone/persistent/users/jiangl/CC1uNPSelection/extra_showerbuiltastrack/covariance_"+syst_unc_list[j]+".root";
   }
   else if(sys_sel==6){
   inputfilename1="/pnfs/uboone/persistent/users/jiangl/CC1uNPSelection/xsec_covmatrix_withadditionalcuts/covariance_detector_"+syst_unc_list[j]+".root";
   //inputfilename2="/pnfs/uboone/persistent/users/jiangl/CC1uNPSelection/extra_showerbuiltastrack/covariance_"+syst_unc_list[j]+".root";
   }
   std::cout<<inputfilename1<<std::endl;
   //std::cout<<inputfilename2<<std::endl;
   inputfile1=new TFile(inputfilename1.c_str());
   //inputfile2=new TFile(inputfilename2.c_str());


   TH2D *covariance_matrix_mumom_1;
   //TH2D *covariance_matrix_mumom_2;
   if(genie_syst==true){
   covariance_matrix_mumom_1=(TH2D*)inputfile1->Get("frac_covariance_matrix_genie_mumom");
   //covariance_matrix_mumom_2=(TH2D*)inputfile2->Get("frac_covariance_matrix_genie_mumom"); 
   }else if(flux_syst==true || total_flux_syst == true){
   covariance_matrix_mumom_1=(TH2D*)inputfile1->Get("frac_covariance_matrix_flux_mumom");
   //covariance_matrix_mumom_2=(TH2D*)inputfile2->Get("frac_covariance_matrix_flux_mumom");     
   }else if(extra_syst==true || total_extra_syst == true){
   covariance_matrix_mumom_1=(TH2D*)inputfile1->Get("frac_covariance_matrix_extra_syst_mumom");
   //covariance_matrix_mumom_2=(TH2D*)inputfile2->Get("frac_covariance_matrix_extra_syst_mumom");     
   }else if(det_syst==true || total_det_syst == true){
   covariance_matrix_mumom_1=(TH2D*)inputfile1->Get("frac_covariance_matrix_detector_mumom");
   //covariance_matrix_mumom_2=(TH2D*)inputfile2->Get("frac_covariance_matrix_extra_syst_mumom");     
   }else if(tune3_syst==true){
     covariance_matrix_mumom_1=(TH2D*)inputfile1->Get("frac_covariance_matrix_detector_mumom");
     //covariance_matrix_mumom_2=(TH2D*)inputfile2->Get("frac_covariance_matrix_extra_syst_mumom");     
   }


   
   TH2D *covariance_matrix_pmom_1;
   //TH2D *covariance_matrix_pmom_2;
   if(genie_syst ==true){
   covariance_matrix_pmom_1=(TH2D*)inputfile1->Get("frac_covariance_matrix_genie_pmom");
   //covariance_matrix_pmom_2=(TH2D*)inputfile2->Get("frac_covariance_matrix_genie_pmom"); 
   }else if(flux_syst==true || total_flux_syst == true){
   covariance_matrix_pmom_1=(TH2D*)inputfile1->Get("frac_covariance_matrix_flux_pmom");
   //covariance_matrix_pmom_2=(TH2D*)inputfile2->Get("frac_covariance_matrix_flux_pmom");     
   }else if(extra_syst==true || total_extra_syst == true){
   covariance_matrix_pmom_1=(TH2D*)inputfile1->Get("frac_covariance_matrix_extra_syst_pmom");
   //covariance_matrix_pmom_2=(TH2D*)inputfile2->Get("frac_covariance_matrix_extra_syst_pmom");     
   }else if(det_syst==true || total_det_syst == true){
   covariance_matrix_pmom_1=(TH2D*)inputfile1->Get("frac_covariance_matrix_detector_pmom");
   //covariance_matrix_pmom_2=(TH2D*)inputfile2->Get("frac_covariance_matrix_extra_syst_pmom");
   }else if(tune3_syst==true){
     covariance_matrix_pmom_1=(TH2D*)inputfile1->Get("frac_covariance_matrix_detector_pmom");
     //covariance_matrix_mumom_2=(TH2D*)inputfile2->Get("frac_covariance_matrix_extra_syst_mumom");     
   }
 
 
   TH2D *covariance_matrix_muangle_1;
   //TH2D *covariance_matrix_muangle_2;
   if(genie_syst==true){
   covariance_matrix_muangle_1=(TH2D*)inputfile1->Get("frac_covariance_matrix_genie_muangle");
   //covariance_matrix_muangle_2=(TH2D*)inputfile2->Get("frac_covariance_matrix_genie_muangle"); 
   }else if(flux_syst==true || total_flux_syst == true){
   covariance_matrix_muangle_1=(TH2D*)inputfile1->Get("frac_covariance_matrix_flux_muangle");
   //covariance_matrix_muangle_2=(TH2D*)inputfile2->Get("frac_covariance_matrix_flux_muangle");     
   }else if(extra_syst==true || total_extra_syst == true){
   covariance_matrix_muangle_1=(TH2D*)inputfile1->Get("frac_covariance_matrix_extra_syst_muangle");
   //covariance_matrix_muangle_2=(TH2D*)inputfile2->Get("frac_covariance_matrix_extra_syst_muangle");     
   }else if(det_syst==true || total_det_syst == true){
   covariance_matrix_muangle_1=(TH2D*)inputfile1->Get("frac_covariance_matrix_detector_muangle");
   //covariance_matrix_muangle_2=(TH2D*)inputfile2->Get("frac_covariance_matrix_extra_syst_muangle");     
   }else if(tune3_syst==true){
     covariance_matrix_muangle_1=(TH2D*)inputfile1->Get("frac_covariance_matrix_detector_muangle");
     //covariance_matrix_mumom_2=(TH2D*)inputfile2->Get("frac_covariance_matrix_extra_syst_mumom");     
   }

 
   TH2D *covariance_matrix_pangle_1;
   //TH2D *covariance_matrix_pangle_2;
   if(genie_syst==true){
   covariance_matrix_pangle_1=(TH2D*)inputfile1->Get("frac_covariance_matrix_genie_pangle");
   //covariance_matrix_pangle_2=(TH2D*)inputfile2->Get("frac_covariance_matrix_genie_pangle"); 
   }else if(flux_syst==true || total_flux_syst == true){
   covariance_matrix_pangle_1=(TH2D*)inputfile1->Get("frac_covariance_matrix_flux_pangle");
   //covariance_matrix_pangle_2=(TH2D*)inputfile2->Get("frac_covariance_matrix_flux_pangle");     
   }else if(extra_syst==true || total_extra_syst == true){
   covariance_matrix_pangle_1=(TH2D*)inputfile1->Get("frac_covariance_matrix_extra_syst_pangle");
   //covariance_matrix_pangle_2=(TH2D*)inputfile2->Get("frac_covariance_matrix_extra_syst_pangle");     
   }else if(det_syst==true || total_det_syst == true){
   covariance_matrix_pangle_1=(TH2D*)inputfile1->Get("frac_covariance_matrix_detector_pangle");
   //covariance_matrix_pangle_2=(TH2D*)inputfile2->Get("frac_covariance_matrix_extra_syst_pangle");     
   }else if(tune3_syst==true){
     covariance_matrix_pangle_1=(TH2D*)inputfile1->Get("frac_covariance_matrix_detector_pangle");
     //covariance_matrix_mumom_2=(TH2D*)inputfile2->Get("frac_covariance_matrix_extra_syst_mumom");     
   }


   TH2D *covariance_matrix_thetamup_1;
   //TH2D *covariance_matrix_thetamup_2;
   if(genie_syst==true){
   covariance_matrix_thetamup_1=(TH2D*)inputfile1->Get("frac_covariance_matrix_genie_thetamup");
   //covariance_matrix_thetamup_2=(TH2D*)inputfile2->Get("frac_covariance_matrix_genie_thetamup"); 
   }else if(flux_syst==true || total_flux_syst == true){
   covariance_matrix_thetamup_1=(TH2D*)inputfile1->Get("frac_covariance_matrix_flux_thetamup");
   //covariance_matrix_thetamup_2=(TH2D*)inputfile2->Get("frac_covariance_matrix_flux_thetamup");     
   }else if(extra_syst==true || total_extra_syst == true){
   covariance_matrix_thetamup_1=(TH2D*)inputfile1->Get("frac_covariance_matrix_extra_syst_thetamup");
   //covariance_matrix_thetamup_2=(TH2D*)inputfile2->Get("frac_covariance_matrix_extra_syst_thetamup");     
   }else if(det_syst==true || total_det_syst == true){
   covariance_matrix_thetamup_1=(TH2D*)inputfile1->Get("frac_covariance_matrix_detector_thetamup");
   //covariance_matrix_thetamup_2=(TH2D*)inputfile2->Get("frac_covariance_matrix_extra_syst_thetamup");     
   }else if(tune3_syst==true){
     covariance_matrix_thetamup_1=(TH2D*)inputfile1->Get("frac_covariance_matrix_detector_thetamup");
     //covariance_matrix_mumom_2=(TH2D*)inputfile2->Get("frac_covariance_matrix_extra_syst_mumom");     
   }

   std::cout<<"start creating a canvas to make plots"<<std::endl;

    
   std::cout<<"total number of bins of mumom is "<<covariance_matrix_mumom_1->GetNbinsX()<<std::endl;

   Int_t nbins=covariance_matrix_mumom_1->GetNbinsX();
   Double_t xsec_unc_mumom_1[nbins];
   //Double_t xsec_unc_mumom_2[nbins];
   Double_t xbinind_mumom[nbins];
   
   std::cout<<"libo test 2"<<std::endl;
  
   for(int i=0; i<nbins; i++){
     xsec_unc_mumom_1[i]=covariance_matrix_mumom_1->GetBinContent(i+1,i+1);
     //***********************************************************
     if (syst_unc_list[j]=="lifetime10ms") { xsec_unc_mumom_1[i]=sqrt(xsec_unc_mumom_1[i])*0.1; }
     //THIS IS HERE BECAUSE WE NEED TO SCALE THE LIFETIME UNCERTAINTY BY 0.10 SINCE ONLY 10% OF THE DATA MIGHT BE AFFECTED
     else { xsec_unc_mumom_1[i]=sqrt(xsec_unc_mumom_1[i]); }

     //xsec_unc_mumom_2[i]=covariance_matrix_mumom_2->GetBinContent(i+1,i+1);
     //xsec_unc_mumom_2[i]=sqrt(xsec_unc_mumom_2[i]);

     xbinind_mumom[i]=bins_mumom[i]+(bins_mumom[i+1]-bins_mumom[i])/2.0;
     gr_mumom_1[j]->SetBinContent(i+1, xsec_unc_mumom_1[i]);
     total_mumom_1->SetBinContent(i+1, total_mumom_1->GetBinContent(i+1) + xsec_unc_mumom_1[i]*xsec_unc_mumom_1[i]);
     //gr_mumom_2[j]->SetBinContent(i+1, xsec_unc_mumom_2[i]);

     //std::cout<<"j=  "<<j<<" "<<i<<" th"<<xsec_unc_mumom_1[i]<<std::endl; 
   }
   //gr_mumom_1[j]=new TGraph(nbins, xbinind_mumom, xsec_unc_mumom_1);
   //gr_mumom_1[j]->SetMarkerStyle(21);
   //gr_mumom_1[j]->SetDrawOption("AP");
   if (j<4) gr_mumom_1[j]->SetLineColor(kCyan+j);
   else if (j<8) gr_mumom_1[j]->SetLineColor(kGreen+j-4);
   else if (j<12) gr_mumom_1[j]->SetLineColor(kMagenta+j-8);
   else if (j<16) gr_mumom_1[j]->SetLineColor(kGray+j-12);
   else gr_mumom_1[j]->SetLineColor(kYellow+j-16);
   //   if(j+2==2) {gr_mumom_1[j]->SetLineColor(28);}
   gr_mumom_1[j]->SetLineWidth(4);
   gr_mumom_1[j]->SetFillStyle(0);
   gr_mumom_1[j]->SetMinimum(0); 
   gr_mumom_1[j]->SetMaximum(1.0);

   if(sys_sel==4 || sys_sel ==5) {gr_mumom_1[j]->SetLineColor(6+j);}
   //gr_mumom_2[j]->SetMarkerStyle(22);
   //gr_mumom_2[j]->SetDrawOption("P");
   //gr_mumom_2[j]->SetLineColor(2+j);
   //if(j+2==5) {gr_mumom_2[j]->SetLineColor(28);}
   //gr_mumom_2[j]->SetLineStyle(2);
   //gr_mumom_2[j]->SetLineWidth(4);
   //gr_mumom_2[j]->SetFillStyle(0);
   //gr_mumom_2[j]->SetMinimum(0);
   //gr_mumom_2[j]->SetMaximum(0.3);   

   
   //================================================================= 
   std::cout<<"total number of bins of pmom is "<<covariance_matrix_pmom_1->GetNbinsX()<<std::endl;

   nbins=covariance_matrix_pmom_1->GetNbinsX();
   Double_t xsec_unc_pmom_1[nbins];
   //Double_t xsec_unc_pmom_2[nbins];
   Double_t xbinind_pmom[nbins];
   //mg->SetTitle("pmom; pmom(nbins); Relative Uncertainty");
   //std::cout<<"libo test 0"<<std::endl;
   for(int i=0; i<nbins; i++){
     std::cout<<"The value of the frac cov matrix is: " << covariance_matrix_pmom_1->GetBinContent(i+1,i+1) << std::endl;
     std::cout<<"The sqrt of the frac cov matrix is: " << sqrt(covariance_matrix_pmom_1->GetBinContent(i+1,i+1)) << std::endl;
     xsec_unc_pmom_1[i]=covariance_matrix_pmom_1->GetBinContent(i+1,i+1);
     //***********************************************************
     if (syst_unc_list[j]=="lifetime10ms") { xsec_unc_pmom_1[i]=sqrt(xsec_unc_pmom_1[i])*0.1; }
     //THIS IS HERE BECAUSE WE NEED TO SCALE THE LIFETIME UNCERTAINTY BY 0.10 SINCE ONLY 10% OF THE DATA MIGHT BE AFFECTED
     else { xsec_unc_pmom_1[i]=sqrt(xsec_unc_pmom_1[i]); }

     //xsec_unc_pmom_2[i]=covariance_matrix_pmom_2->GetBinContent(i+1,i+1);
     //xsec_unc_pmom_2[i]=sqrt(xsec_unc_pmom_2[i]);

     xbinind_pmom[i]=bins_pmom[i]+(bins_pmom[i+1]-bins_pmom[i])/2.0;
     std::cout<<"The value filled into the plot is: " << xsec_unc_pmom_1[i] << std::endl; 
     gr_pmom_1[j]->SetBinContent(i+1, xsec_unc_pmom_1[i]);
     total_pmom_1->SetBinContent(i+1, total_pmom_1->GetBinContent(i+1) + xsec_unc_pmom_1[i]*xsec_unc_pmom_1[i]);
     //gr_pmom_2[j]->SetBinContent(i+1, xsec_unc_pmom_2[i]);

     //std::cout<<"j=  "<<j<<" "<<i<<" th"<<xsec_unc_pmom_1[i]<<std::endl; 
   }
   //gr_pmom_1[j]=new TGraph(nbins, xbinind_pmom, xsec_unc_pmom_1);
   //gr_pmom_1[j]->SetMarkerStyle(21);
   //gr_pmom_1[j]->SetDrawOption("AP");
   if (j<4) gr_pmom_1[j]->SetLineColor(kCyan+j);
   else if (j<8) gr_pmom_1[j]->SetLineColor(kGreen+j-4);
   else if (j<12) gr_pmom_1[j]->SetLineColor(kMagenta+j-8);
   else if (j<16) gr_pmom_1[j]->SetLineColor(kGray+j-12);
   else gr_pmom_1[j]->SetLineColor(kYellow+j-16);
   //gr_pmom_1[j]->SetLineColor(2+j);
   //if(j+2==5) {gr_pmom_1[j]->SetLineColor(28);}
   gr_pmom_1[j]->SetLineWidth(4);
   gr_pmom_1[j]->SetFillStyle(0);
   gr_pmom_1[j]->SetMinimum(0); 
   gr_pmom_1[j]->SetMaximum(1.0);

   if(sys_sel==4 || sys_sel ==5) {gr_pmom_1[j]->SetLineColor(6+j);}
   //gr_pmom_2[j]->SetMarkerStyle(22);
   //gr_pmom_2[j]->SetDrawOption("P");
   //gr_pmom_2[j]->SetLineColor(2+j);
   //if(j+2==5) {gr_pmom_2[j]->SetLineColor(28);}
   //gr_pmom_2[j]->SetLineStyle(2);
   //gr_pmom_2[j]->SetLineWidth(4);
   //gr_pmom_2[j]->SetFillStyle(0);
   //gr_pmom_2[j]->SetMinimum(0);
   //gr_pmom_2[j]->SetMaximum(0.3);   
   //=======================================================================
   std::cout<<"total number of bins of muangle is "<<covariance_matrix_muangle_1->GetNbinsX()<<std::endl;

   nbins=covariance_matrix_muangle_1->GetNbinsX();
   Double_t xsec_unc_muangle_1[nbins];
   //Double_t xsec_unc_muangle_2[nbins];
   Double_t xbinind_muangle[nbins];
   //mg->SetTitle("muangle; muangle(nbins); Relative Uncertainty");
   //std::cout<<"libo test 0"<<std::endl;
   for(int i=0; i<nbins; i++){
     xsec_unc_muangle_1[i]=covariance_matrix_muangle_1->GetBinContent(i+1,i+1);
     //***********************************************************
     if (syst_unc_list[j]=="lifetime10ms") { xsec_unc_muangle_1[i]=sqrt(xsec_unc_muangle_1[i])*0.1; }
     //THIS IS HERE BECAUSE WE NEED TO SCALE THE LIFETIME UNCERTAINTY BY 0.10 SINCE ONLY 10% OF THE DATA MIGHT BE AFFECTED
     else { xsec_unc_muangle_1[i]=sqrt(xsec_unc_muangle_1[i]); }

     //xsec_unc_muangle_2[i]=covariance_matrix_muangle_2->GetBinContent(i+1,i+1);
     //xsec_unc_muangle_2[i]=sqrt(xsec_unc_muangle_2[i]);

     //xbinind_muangle[i]=bins_muangle[i]+(bins_muangle[i+1]-bins_muangle[i])/2.0;
     gr_muangle_1[j]->SetBinContent(i+1, xsec_unc_muangle_1[i]);
     total_muangle_1->SetBinContent(i+1, total_muangle_1->GetBinContent(i+1) + xsec_unc_muangle_1[i]*xsec_unc_muangle_1[i]);
     //gr_muangle_2[j]->SetBinContent(i+1, xsec_unc_muangle_2[i]);

     //std::cout<<"j=  "<<j<<" "<<i<<" th"<<xsec_unc_muangle_1[i]<<std::endl; 
   }
   //gr_muangle_1[j]=new TGraph(nbins, xbinind_muangle, xsec_unc_muangle_1);
   //gr_muangle_1[j]->SetMarkerStyle(21);
   //gr_muangle_1[j]->SetDrawOption("AP");
   if (j<4) gr_muangle_1[j]->SetLineColor(kCyan+j);
   else if (j<8) gr_muangle_1[j]->SetLineColor(kGreen+j-4);
   else if (j<12) gr_muangle_1[j]->SetLineColor(kMagenta+j-8);
   else if (j<16) gr_muangle_1[j]->SetLineColor(kGray+j-12);
   else gr_muangle_1[j]->SetLineColor(kYellow+j-16);
   //gr_muangle_1[j]->SetLineColor(2+j);
   //if(j+2==5) {gr_muangle_1[j]->SetLineColor(28);}
   gr_muangle_1[j]->SetLineWidth(4);
   gr_muangle_1[j]->SetFillStyle(0);
   gr_muangle_1[j]->SetMinimum(0); 
   gr_muangle_1[j]->SetMaximum(1.0);

   if(sys_sel==4 || sys_sel ==5) {gr_muangle_1[j]->SetLineColor(6+j);}
   //gr_muangle_2[j]->SetMarkerStyle(22);
   //gr_muangle_2[j]->SetDrawOption("P");
   //gr_muangle_2[j]->SetLineColor(2+j);
   //if(j+2==5) {gr_muangle_2[j]->SetLineColor(28);}
   //gr_muangle_2[j]->SetLineStyle(2);
   //gr_muangle_2[j]->SetLineWidth(4);
   //gr_muangle_2[j]->SetFillStyle(0);
   //gr_muangle_2[j]->SetMinimum(0);
   //gr_muangle_2[j]->SetMaximum(0.3);   
   //==============================================================
   std::cout<<"total number of bins of pangle is "<<covariance_matrix_pangle_1->GetNbinsX()<<std::endl;

   nbins=covariance_matrix_pangle_1->GetNbinsX();
   Double_t xsec_unc_pangle_1[nbins];
   //Double_t xsec_unc_pangle_2[nbins];
   Double_t xbinind_pangle[nbins];
   //mg->SetTitle("pangle; pangle(nbins); Relative Uncertainty");
   //std::cout<<"libo test 0"<<std::endl;
   for(int i=0; i<nbins; i++){
     xsec_unc_pangle_1[i]=covariance_matrix_pangle_1->GetBinContent(i+1,i+1);
     if (syst_unc_list[j]=="lifetime10ms") { xsec_unc_pangle_1[i]=sqrt(xsec_unc_pangle_1[i])*0.1; }
     //THIS IS HERE BECAUSE WE NEED TO SCALE THE LIFETIME UNCERTAINTY BY 0.10 SINCE ONLY 10% OF THE DATA MIGHT BE AFFECTED
     else { xsec_unc_pangle_1[i]=sqrt(xsec_unc_pangle_1[i]); }

     //xsec_unc_pangle_2[i]=covariance_matrix_pangle_2->GetBinContent(i+1,i+1);
     //xsec_unc_pangle_2[i]=sqrt(xsec_unc_pangle_2[i]);

     //xbinind_pangle[i]=bins_pangle[i]+(bins_pangle[i+1]-bins_pangle[i])/2.0;
     gr_pangle_1[j]->SetBinContent(i+1, xsec_unc_pangle_1[i]);
     total_pangle_1->SetBinContent(i+1, total_pangle_1->GetBinContent(i+1) + xsec_unc_pangle_1[i]*xsec_unc_pangle_1[i]);
     //gr_pangle_2[j]->SetBinContent(i+1, xsec_unc_pangle_2[i]);

     //std::cout<<"j=  "<<j<<" "<<i<<" th"<<xsec_unc_pangle_1[i]<<std::endl; 
   }
   //gr_pangle_1[j]=new TGraph(nbins, xbinind_pangle, xsec_unc_pangle_1);
   //gr_pangle_1[j]->SetMarkerStyle(21);
   //gr_pangle_1[j]->SetDrawOption("AP");
   if (j<4) gr_pangle_1[j]->SetLineColor(kCyan+j);
   else if (j<8) gr_pangle_1[j]->SetLineColor(kGreen+j-4);
   else if (j<12) gr_pangle_1[j]->SetLineColor(kMagenta+j-8);
   else if (j<16) gr_pangle_1[j]->SetLineColor(kGray+j-12);
   else gr_pangle_1[j]->SetLineColor(kYellow+j-16);
   //gr_pangle_1[j]->SetLineColor(2+j);
   //if(j+2==5) {gr_pangle_1[j]->SetLineColor(28);}
   gr_pangle_1[j]->SetLineWidth(4);
   gr_pangle_1[j]->SetFillStyle(0);
   gr_pangle_1[j]->SetMinimum(0); 
   gr_pangle_1[j]->SetMaximum(1.0);

   if(sys_sel==4 || sys_sel ==5) {gr_pangle_1[j]->SetLineColor(6+j);}
   //gr_pangle_2[j]->SetMarkerStyle(22);
   //gr_pangle_2[j]->SetDrawOption("P");
   //gr_pangle_2[j]->SetLineColor(2+j);
   //if(j+2==5) {gr_pangle_2[j]->SetLineColor(28);}
   //gr_pangle_2[j]->SetLineStyle(2);
   //gr_pangle_2[j]->SetLineWidth(4);
   //gr_pangle_2[j]->SetFillStyle(0);
   //gr_pangle_2[j]->SetMinimum(0);
   //gr_pangle_2[j]->SetMaximum(0.3);   
   //====================================================================
   std::cout<<"total number of bins of thetamup is "<<covariance_matrix_thetamup_1->GetNbinsX()<<std::endl;

   nbins=covariance_matrix_thetamup_1->GetNbinsX();
   Double_t xsec_unc_thetamup_1[nbins];
   //Double_t xsec_unc_thetamup_2[nbins];
   Double_t xbinind_thetamup[nbins];
   //mg->SetTitle("thetamup; thetamup(nbins); Relative Uncertainty");
   //std::cout<<"libo test 0"<<std::endl;
   for(int i=0; i<nbins; i++){
     xsec_unc_thetamup_1[i]=covariance_matrix_thetamup_1->GetBinContent(i+1,i+1);
     if (syst_unc_list[j]=="lifetime10ms") { xsec_unc_thetamup_1[i]=sqrt(xsec_unc_thetamup_1[i])*0.1; }
     //THIS IS HERE BECAUSE WE NEED TO SCALE THE LIFETIME UNCERTAINTY BY 0.10 SINCE ONLY 10% OF THE DATA MIGHT BE AFFECTED
     else { xsec_unc_thetamup_1[i]=sqrt(xsec_unc_thetamup_1[i]); }

     //xsec_unc_thetamup_2[i]=covariance_matrix_thetamup_2->GetBinContent(i+1,i+1);
     //xsec_unc_thetamup_2[i]=sqrt(xsec_unc_thetamup_2[i]);

     //xbinind_thetamup[i]=bins_thetamup[i]+(bins_thetamup[i+1]-bins_thetamup[i])/2.0;
     gr_thetamup_1[j]->SetBinContent(i+1, xsec_unc_thetamup_1[i]);
     total_thetamup_1->SetBinContent(i+1, total_thetamup_1->GetBinContent(i+1) + xsec_unc_thetamup_1[i]*xsec_unc_thetamup_1[i]);
     //gr_thetamup_2[j]->SetBinContent(i+1, xsec_unc_thetamup_2[i]);

     //std::cout<<"j=  "<<j<<" "<<i<<" th"<<xsec_unc_thetamup_1[i]<<std::endl; 
   }
   //gr_thetamup_1[j]=new TGraph(nbins, xbinind_thetamup, xsec_unc_thetamup_1);
   //gr_thetamup_1[j]->SetMarkerStyle(21);
   //gr_thetamup_1[j]->SetDrawOption("AP");
   if (j<4) gr_thetamup_1[j]->SetLineColor(kCyan+j);
   else if (j<8) gr_thetamup_1[j]->SetLineColor(kGreen+j-4);
   else if (j<12) gr_thetamup_1[j]->SetLineColor(kMagenta+j-8);
   else if (j<16) gr_thetamup_1[j]->SetLineColor(kGray+j-12);
   else gr_thetamup_1[j]->SetLineColor(kYellow+j-16);
   //gr_thetamup_1[j]->SetLineColor(2+j);
   //if(j+2==5) {gr_thetamup_1[j]->SetLineColor(28);}
   gr_thetamup_1[j]->SetLineWidth(4);
   gr_thetamup_1[j]->SetFillStyle(0);
   gr_thetamup_1[j]->SetMinimum(0); 
   gr_thetamup_1[j]->SetMaximum(1.0);

   if(sys_sel==4 || sys_sel ==5) {gr_thetamup_1[j]->SetLineColor(6+j);}
   //gr_thetamup_2[j]->SetMarkerStyle(22);
   //gr_thetamup_2[j]->SetDrawOption("P");
   //gr_thetamup_2[j]->SetLineColor(2+j);
   //if(j+2==5) {gr_thetamup_2[j]->SetLineColor(28);}
   //gr_thetamup_2[j]->SetLineStyle(2);
   //gr_thetamup_2[j]->SetLineWidth(4);
   //gr_thetamup_2[j]->SetFillStyle(0);
   //gr_thetamup_2[j]->SetMinimum(0);
   //gr_thetamup_2[j]->SetMaximum(0.3);   

   /*
   */
 }
 
for(int i=0; i<total_mumom_1->GetNbinsX(); i++){  total_mumom_1->SetBinContent(i+1, sqrt(total_mumom_1->GetBinContent(i+1))); }
for(int i=0; i<total_muangle_1->GetNbinsX(); i++){  total_muangle_1->SetBinContent(i+1, sqrt(total_muangle_1->GetBinContent(i+1))); }
for(int i=0; i<total_pmom_1->GetNbinsX(); i++){  total_pmom_1->SetBinContent(i+1, sqrt(total_pmom_1->GetBinContent(i+1))); }
for(int i=0; i<total_pangle_1->GetNbinsX(); i++){  total_pangle_1->SetBinContent(i+1, sqrt(total_pangle_1->GetBinContent(i+1))); }
for(int i=0; i<total_thetamup_1->GetNbinsX(); i++){  total_thetamup_1->SetBinContent(i+1, sqrt(total_thetamup_1->GetBinContent(i+1))); }

 output_file->cd();
 total_mumom_1->Write();
 total_pmom_1->Write();
 total_muangle_1->Write();
 total_pangle_1->Write();
 total_thetamup_1->Write();

  //=============================================================================================================

  TCanvas *c_mumom = new TCanvas("c_mumom"); 
  TLegend* legend = new TLegend(0.6, 0.6, 0.9, 0.9);


  std::cout<<"size of the syst unc vector is "<<syst_unc_list.size()<<std::endl;

  if(sys_sel==0 || sys_sel==4 || sys_sel ==3){gr_testmu->SetMaximum(0.3);}
  if(sys_sel==5){gr_testmu->SetMaximum(0.1);}
  gr_testmu->Draw();
  for(int kk=0; kk<syst_unc_list.size(); kk++){
    if(sys_sel==0 ||sys_sel==4 || sys_sel==3) gr_mumom_1[kk]->SetMaximum(0.5);
    if(sys_sel==5) gr_mumom_1[kk]->SetMaximum(0.1);
    if(kk==0){
         gr_mumom_1[kk]->Draw("same"); }
    else{
         gr_mumom_1[kk]->Draw("same"); }
   legend->AddEntry(gr_mumom_1[kk], syst_unc_list[kk].c_str());
  }
  total_mumom_1->SetLineStyle(10);  
  total_mumom_1->Draw("same");
  legend->AddEntry(total_mumom_1,"Total");
  legend->Draw("same");
  
  string outputfilename="relative_unc_mumom_"+syst_unc_name +".png"; 
  c_mumom->SaveAs(outputfilename.c_str());
  //c_mumom->Delete();

  
  TCanvas *c_pmom = new TCanvas("c_pmom"); 
  if(sys_sel==0 || sys_sel==4 || sys_sel==3) {gr_testp->SetMaximum(0.3);}
  if(sys_sel==5) {gr_testp->SetMaximum(0.2);}
  gr_testp->Draw();
  for(int kk=0; kk<syst_unc_list.size(); kk++){
    if(sys_sel==0) gr_pmom_1[kk]->SetMaximum(0.15);
    if(sys_sel==4) gr_pmom_1[kk]->SetMaximum(0.5);
    if(sys_sel==3) gr_pmom_1[kk]->SetMaximum(0.5);
    if(sys_sel==5) gr_pmom_1[kk]->SetMaximum(0.2);
    if(kk==0){
         gr_pmom_1[kk]->Draw("same"); }
    else{
         gr_pmom_1[kk]->Draw("same"); }
  }  
  total_pmom_1->SetLineStyle(10);  
  total_pmom_1->Draw("same");
  //legend->AddEntry(total_pmom_1,"Total");
  legend->Draw("same");
  
  outputfilename="relative_unc_pmom_"+syst_unc_name +".png"; 
  c_pmom->SaveAs(outputfilename.c_str());
  //c_pmom->Delete();

  TCanvas *c_muangle = new TCanvas("c_muangle"); 
  for(int kk=0; kk<syst_unc_list.size(); kk++){
    if(sys_sel==0 ||sys_sel==4 || sys_sel==3) gr_muangle_1[kk]->SetMaximum(0.5);
    if(sys_sel==5) gr_muangle_1[kk]->SetMaximum(0.1);
    if(kk==0){
         gr_muangle_1[kk]->Draw(); }
    else{
         gr_muangle_1[kk]->Draw("same"); }
  }  
  total_muangle_1->SetLineStyle(10);  
  total_muangle_1->Draw("same");
  //legend->AddEntry(total_muangle_1,"Total");
  legend->Draw("same");
  
  outputfilename="relative_unc_muangle_"+syst_unc_name +".png"; 
  c_muangle->SaveAs(outputfilename.c_str());
  //c_muangle->Delete();

 
  TCanvas *c_pangle = new TCanvas("c_pangle"); 
  for(int kk=0; kk<syst_unc_list.size(); kk++){
    if(sys_sel==0 ||sys_sel==4 || sys_sel==3) gr_pangle_1[kk]->SetMaximum(0.5);
    if(sys_sel==5) gr_pangle_1[kk]->SetMaximum(0.1);
    if(kk==0){
         gr_pangle_1[kk]->Draw(); }
    else{
         gr_pangle_1[kk]->Draw("same"); }
  }  
  total_pangle_1->SetLineStyle(10);  
  total_pangle_1->Draw("same");
  //legend->AddEntry(total_pangle_1,"Total");
  legend->Draw("same");
  
  outputfilename="relative_unc_pangle_"+syst_unc_name +".png"; 
  c_pangle->SaveAs(outputfilename.c_str());
  //c_pangle->Delete();

  TCanvas *c_thetamup = new TCanvas("c_thetamup"); 
  for(int kk=0; kk<syst_unc_list.size(); kk++){
    if(sys_sel==0 ||sys_sel==4 || sys_sel==3) gr_thetamup_1[kk]->SetMaximum(0.5);
    if(sys_sel==5) gr_thetamup_1[kk]->SetMaximum(0.1);
    if(kk==0){
         gr_thetamup_1[kk]->Draw(); }
    else{
         gr_thetamup_1[kk]->Draw("same"); }
  }  
  total_thetamup_1->SetLineStyle(10);  
  total_thetamup_1->Draw("same");
  //legend->AddEntry(total_thetamup_1,"Total");
  legend->Draw("same");
  
  outputfilename="relative_unc_thetamup_"+syst_unc_name +".png"; 
  c_thetamup->SaveAs(outputfilename.c_str());
  //c_thetamup->Delete();

//========================================================================
 //  TCanvas *cc_mumom = new TCanvas("cc_mumom"); 
//   TLegend* legend2 = new TLegend(0.6, 0.7, 0.9, 0.9);


//   std::cout<<"size of the syst unc vector is "<<syst_unc_list.size()<<std::endl;

//   for(int kk=0; kk<syst_unc_list.size(); kk++){
//     if(kk==0){
//          //gr_mumom_2[kk]->Draw(); }
//     }else{
//          //gr_mumom_2[kk]->Draw("same"); 
//     }
//     //legend2->AddEntry(gr_mumom_2[kk], syst_unc_list[kk].c_str());
//   }  
//   legend2->Draw("same");
  
//   outputfilename="relative_unc_mumom_showerastrack_"+syst_unc_name +".png"; 
//   cc_mumom->SaveAs(outputfilename.c_str());
//   //cc_mumom->Delete();

  
//   TCanvas *cc_pmom = new TCanvas("cc_pmom"); 
//   for(int kk=0; kk<syst_unc_list.size(); kk++){
//     if(kk==0){
//          //gr_pmom_2[kk]->Draw(); }
//     } else{
//          //gr_pmom_2[kk]->Draw("same"); 
//     }
//   }  
//   legend2->Draw("same");
  
//   outputfilename="relative_unc_pmom_showerastrack_"+syst_unc_name +".png"; 
//   cc_pmom->SaveAs(outputfilename.c_str());
//   //c_pmom->Delete();

//   TCanvas *cc_muangle = new TCanvas("cc_muangle"); 
//   for(int kk=0; kk<syst_unc_list.size(); kk++){
//     if(kk==0){
//          //gr_muangle_2[kk]->Draw(); }
//     } else{
//          //gr_muangle_2[kk]->Draw("same"); 
//     }
//   }  
//   legend2->Draw("same");
  
//   outputfilename="relative_unc_muangle_showerastrack_"+syst_unc_name +".png"; 
//   cc_muangle->SaveAs(outputfilename.c_str());
//   //cc_muangle->Delete();

 
//   TCanvas *cc_pangle = new TCanvas("cc_pangle"); 
//   for(int kk=0; kk<syst_unc_list.size(); kk++){
//     if(kk==0){
//          //gr_pangle_2[kk]->Draw(); }
//     } else{
//          //gr_pangle_2[kk]->Draw("same"); 
//     }
//   }  
//   legend2->Draw("same");
  
//   outputfilename="relative_unc_pangle_showerastrack_"+syst_unc_name +".png"; 
//   cc_pangle->SaveAs(outputfilename.c_str());
//   //c_pangle->Delete();

//   TCanvas *cc_thetamup = new TCanvas("cc_thetamup"); 
//   for(int kk=0; kk<syst_unc_list.size(); kk++){
//     if(kk==0){
//          //gr_thetamup_2[kk]->Draw(); }
//     }    else{
//          //gr_thetamup_2[kk]->Draw("same"); }
//     }
//   }  
//   legend2->Draw("same");
  
 //  outputfilename="relative_unc_thetamup_showerastrack"+syst_unc_name +".png"; 
//   cc_thetamup->SaveAs(outputfilename.c_str());
//   //cc_thetamup->Delete();
//   TCanvas *c1=new TCanvas("c1", "c1");
//   xsec_mumom_mc->SetMaximum(0.5*xsec_mumom_mc->GetMaximum());
//   xsec_mumom_CV->SetMaximum(0.5*xsec_mumom_mc->GetMaximum()); 
//   xsec_mumom_mc->SetTitleOffset(0.8, "X");
//   xsec_mumom_CV->SetTitleOffset(0.8, "X");

//   xsec_mumom_mc->Draw();
//   xsec_mumom_CV->Draw("same");
//   leg->Draw("same");

//   TCanvas *c2=new TCanvas("c2", "c2");
//   xsec_pmom_mc->SetMaximum(0.5*xsec_pmom_mc->GetMaximum());
//   xsec_pmom_CV->SetMaximum(0.5*xsec_pmom_mc->GetMaximum()); 
//   xsec_pmom_mc->SetTitleOffset(0.8, "X");
//   xsec_pmom_CV->SetTitleOffset(0.8, "X");

//   xsec_pmom_mc->Draw();
//   xsec_pmom_CV->Draw("same");
//   leg->Draw("same");

//   TCanvas *c3=new TCanvas("c3", "c3");
//   xsec_muangle_mc->SetMaximum(0.5*xsec_muangle_mc->GetMaximum());
//   xsec_muangle_CV->SetMaximum(0.5*xsec_muangle_mc->GetMaximum()); 
//   xsec_muangle_mc->SetTitleOffset(0.8, "X");
//   xsec_muangle_CV->SetTitleOffset(0.8, "X");

//   xsec_muangle_mc->Draw();
//   xsec_muangle_CV->Draw("same");
//   leg->Draw("same");

//   TCanvas *c4=new TCanvas("c4", "c4");
//   xsec_pangle_mc->SetMaximum(0.5*xsec_pangle_mc->GetMaximum());
//   xsec_pangle_CV->SetMaximum(0.5*xsec_pangle_mc->GetMaximum()); 
//   xsec_pangle_mc->SetTitleOffset(0.8, "X");
//   xsec_pangle_CV->SetTitleOffset(0.8, "X");

//   xsec_pangle_mc->Draw();
//   xsec_pangle_CV->Draw("same");
//   leg->Draw("same");

//   TCanvas *c5=new TCanvas("c5", "c5");
//   xsec_thetamup_mc->SetMaximum(0.2*xsec_thetamup_mc->GetMaximum());
//   xsec_thetamup_CV->SetMaximum(0.2*xsec_thetamup_mc->GetMaximum()); 
//   xsec_thetamup_mc->SetTitleOffset(0.8, "X");
//   xsec_thetamup_CV->SetTitleOffset(0.8, "X");

//   xsec_thetamup_mc->Draw();
//   xsec_thetamup_CV->Draw("same");
//   leg->Draw("same");


  output_file->Write();
  
//========================================================================
}
