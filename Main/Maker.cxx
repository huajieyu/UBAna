#ifndef __MAIN_MAKER_CXX__
#define __MAIN_MAKER_CXX__

#include "Maker.h"

using namespace Base;

void Main::Maker::SetInputFile(std::string in)
{
  filen = in;
}

void Main::Maker::SetOutputFile(std::string in)
{
  fileoutn = in;
}

void Main::Maker::SetEntries(int e)
{
  maxEntries = e;
}

void Main::Maker::SetInitialEntry(int e)
{
  _initial_entry = e;
}

void Main::Maker::SetBeamSpillStart(double v)
{
  _beamSpillStarts = v;
}

void Main::Maker::SetBeamSpillEnd(double v)
{
  _beamSpillEnds = v;
}

void Main::Maker::SetFlashShift(double v)
{
  _flashShift = v;
}

void Main::Maker::SetGainCalibration(double v)
{
  _gainCalib = v;
}

void Main::Maker::SetCalculatePOT(bool v)
{
  evalPOT = v;
}


void Main::Maker::SetIsData(bool v)
{
  isdata = v;
}

float temp_pmom;
float temp_pangle;
float temp_pphi;
float temp_thetamup;

float temp_geantpmom;
float temp_geantpenergy;


int muind=-999;
int pind=-999;

bool pion_reco=false;
bool pion_reint=false;
float thetamup = -999.0;
float ptmis = -999.0;
float alphat= -999.0;
float phit = -999.0;
float enucal=-999.0;
float etatest=-999.0;
float FVx = 256.35;
float FVy = 233;
float FVz = 1036.8;
float borderx = 10.;
float bordery = 20.;
float borderz = 10.;


//This function returns if a 3D point is within the fiducial volume
bool Main::Maker::inCV(float x, float y, float z) {
    if(x < (FVx - borderx) && (x > borderx) && (y < (FVy/2. - bordery)) && (y > (-FVy/2. + bordery)) && (z < (FVz - borderz)) && (z > borderz)) return true;
    //if(x < (FVx - borderx) && (x > borderx) && (y < (FVy/2. - bordery)) && (y > (-FVy/2. + bordery)) && (z < (FVz - 85)) && (z > borderz)) return true;
    else return false;
}
float Main::Maker::getEta(vector<vector<double>> canddQdx, vector<vector<double>> trkRR, vector<float> trklen,  int muind, int pind){
      int nhitmu=0; int nhitp=0; double deltaEmu=0; double deltaEp=0;
      //loop over all the dQdx of the muon candidate
      for(unsigned int kk=0; kk<canddQdx[muind].size(); kk++){
	if(trkRR[muind][kk] > (trklen[muind]-5) ) continue;
	deltaEmu = deltaEmu+ canddQdx[muind][kk];     
	nhitmu=nhitmu+1;	
      } 
      //loop over all the dQdx of the proton candidate
      for(unsigned int kk=0; kk<canddQdx[pind].size(); kk++){
	if(trkRR[pind][kk] > (trklen[pind]-5) ) continue;
	deltaEp = deltaEp+canddQdx[pind][kk];      
	nhitp=nhitp+1;    
      }
      //double nuttest=(deltaEp-deltaEmu)/(deltaEmu+deltaEp);
      double nuttest2=(deltaEp/nhitp-deltaEmu/nhitmu)/(deltaEmu/nhitmu+deltaEp/nhitp);
       
  
    return nuttest2;
}
float Main::Maker::getAngle(float mom1, float theta1, float phi1, float mom2, float theta2, float phi2) {
      float px1=mom1*TMath::Sin(theta1)*TMath::Cos(phi1);
      float py1=mom1*TMath::Sin(theta1)*TMath::Sin(phi1);
      float pz1=mom1*TMath::Cos(theta1);


      float px2=mom2*TMath::Sin(theta2)*TMath::Cos(phi2);
      float py2=mom2*TMath::Sin(theta2)*TMath::Sin(phi2);
      float pz2=mom2*TMath::Cos(theta2);    

      TVector3 p4_mom1(px1, py1, pz1);
      TVector3 p4_mom2(px2, py2, pz2);
      return p4_mom1.Angle(p4_mom2);

}
float Main::Maker::Ecalomiss(float Esum, float PTmiss, int np) {
   Esum *= 1000; //convert to MeV
   PTmiss *= 1000; //convert to MeV
   float Eexcit = 30.4; //in MeV
   float Mass = 0; // in MeV
   if(np == 0) Mass = 37.2050e3; //Ar40
   else if(np == 1) Mass = 36.2758e3; //Ar39
   else if(np == 2) Mass = 35.3669e3; //Cl38
   else if(np == 3) Mass = 34.4201e3; //S37
   else if(np == 4) Mass = 33.4957e3; //P36
   else if(np == 5) Mass = 32.5706e3; //Si35
   else if(np == 6) Mass = 31.6539e3; //Al34
   else if(np == 7) Mass = 30.7279e3; //Mg33
   else if(np == 8) Mass = 29.8111e3; //Na32
   else if(np == 9) Mass = 28.8918e3; //Ne31
   else if(np >= 10) Mass = 27.9789e3; //F30

   float Ekinrecoil = sqrt(PTmiss*PTmiss + Mass*Mass) - Mass;
   return Esum + Eexcit + Ekinrecoil; // return result in MeV
}




void Main::Maker::PrintConfig()
{
  LOG_INFO() << "--- Main::Maker::PrintConfig" << std::endl;

  LOG_INFO() << "--- _breakdownPlots: " << _breakdownPlots << std::endl;
  LOG_INFO() << "--- _makePlots " << _makePlots << std::endl;
  LOG_INFO() << "--- _fill_bootstrap_flux " << _fill_bootstrap_flux << std::endl;
  LOG_INFO() << "--- _fill_bootstrap_genie " << _fill_bootstrap_genie << std::endl;
  LOG_INFO() << "--- _target_flux_syst " << _target_flux_syst << std::endl;
  LOG_INFO() << "--- _check_duplicate_events " << _check_duplicate_events << std::endl;

  LOG_INFO() << "--- _beamSpillStarts " << _beamSpillStarts << std::endl;
  LOG_INFO() << "--- _beamSpillEnds " << _beamSpillEnds << std::endl;
  LOG_INFO() << "--- _flashShift " << _flashShift << std::endl;
  LOG_INFO() << "--- _gainCalib " << _gainCalib << std::endl;

  LOG_INFO() << "--- filen " << filen << std::endl;
  LOG_INFO() << "--- evalPOT " << evalPOT << std::endl;
  LOG_INFO() << "--- maxEntries " << maxEntries << std::endl;
  LOG_INFO() << "--- isdata " << isdata << std::endl;

  LOG_INFO() << "--- _pe_cut " << _pe_cut << std::endl;

  LOG_INFO() << "--- targetPOT " << targetPOT << std::endl;
}

void Main::Maker::PrintMaUpMECOff()
{
  for (int i = 0; i < 10; i++) {
    std::cout << "**************************** RUNNING WITH MA+1SIGMA AND MEC OFF ****************************" << std::endl;
  }
}

void Main::Maker::PrintReweighKaons()
{
  for (int i = 0; i < 10; i++) {
    std::cout << "**************************** RUNNING WITH KAON FLUX SCALED BY " << _kaon_reweigh_factor << " ****************************" << std::endl;
  }
}








//____________________________________________________________________________________________________
void Main::Maker::DrawProgressBar(double progress, double barWidth) {
  
  std::cout << "[";
  int pos = barWidth * progress;
  for (int i = 0; i < barWidth; ++i) {
    if (i < pos) std::cout << "=";
    else if (i == pos) std::cout << ">";
    else std::cout << " ";
  }
  std::cout << "] " << int(progress * 100.0) << " %\r";
  std::cout.flush();
}

//____________________________________________________________________________________________________
void Main::Maker::DrawPOT2(double pot, double target)
{
  //std::string str = "Simulated POT:" + std::to_string(pot);
  
  std::stringstream sstm;
  sstm << "Simulated POT: " << pot;
  std::string str = sstm.str();
  
  TLatex* pot_latex = new TLatex(.10, .96, str.c_str());
  pot_latex->SetTextColor(kGray+2);
  pot_latex->SetNDC();
  pot_latex->SetTextSize(1/30.);
  pot_latex->SetTextAlign(10); //left adjusted
  pot_latex->Draw();
  
  
  std::stringstream sstm2;
  sstm2 << "Scaled to POT: " << target;
  str = sstm2.str();
  
  TLatex* pot_latex_2 = new TLatex(.10, .92, str.c_str());
  pot_latex_2->SetTextFont(62);
  pot_latex_2->SetTextColor(kGray+2);
  pot_latex_2->SetNDC();
  pot_latex_2->SetTextSize(1/30.);
  pot_latex_2->SetTextAlign(10);//left adjusted
  pot_latex_2->Draw();
  
}



//____________________________________________________________________________________________________
double Main::Maker::eff_uncertainty(int _n, int _N) {

  double n = (double) _n;
  double N = (double) _N;

  double unc = 1/std::sqrt(N) * std::sqrt((n/N)*(1-n/N));

  return unc;

}


//___________________________________________________________________________________________________
void Main::Maker::FillBootstrap(double fill_value,
                                double evt_wgt,
                                std::map<std::string,std::map<std::string,TH1D*>> hmap_trkmom_genie_pm1_bs, 
                                std::string channel_namel, 
                                std::vector<std::string> fname, 
                                std::vector<double> wgts_genie) {

  auto iter = hmap_trkmom_genie_pm1_bs.find(channel_namel);
  if (iter == hmap_trkmom_genie_pm1_bs.end()) {
    LOG_CRITICAL() << "Can't find " << channel_namel << std::endl;
    throw std::exception();
  }
  std::map<std::string,TH1D*> this_map = iter->second;

  this_map["nominal"]->Fill(fill_value, evt_wgt);

  for (size_t i = 0; i < fname.size(); i++) {

    //std::cout << i << ": Fill value: " << fill_value << ", weight: " << wgts_genie.at(i) << ", name is " << fname.at(i) << std::endl;

    this_map[fname.at(i)]->Fill(fill_value, wgts_genie.at(i) * evt_wgt);

  }


} 

//___________________________________________________________________________________________________
void Main::Maker::FillBootstrap(double fill_value1,
                                double fill_value2,
                                double evt_wgt,
                                std::map<std::string,std::map<std::string,TH2D*>> hmap, 
                                std::string channel_namel, 
                                std::vector<std::string> fname, 
                                std::vector<double> wgts_genie) {

  auto iter = hmap.find(channel_namel);
  if (iter == hmap.end()) {
    LOG_CRITICAL() << "Can't find " << channel_namel << std::endl;
    throw std::exception();
  }
  std::map<std::string,TH2D*> this_map = iter->second;

  this_map["nominal"]->Fill(fill_value1, fill_value2, evt_wgt);

  for (size_t i = 0; i < fname.size(); i++) {

    //std::cout << i << ": Fill value: " << fill_value << ", weight: " << wgts_genie.at(i) << ", name is " << fname.at(i) << std::endl;

    this_map[fname.at(i)]->Fill(fill_value1, fill_value2, wgts_genie.at(i) * evt_wgt);

  }


} 


//___________________________________________________________________________________________________
void Main::Maker::FillBootstrap(double fill_value1,
                                double fill_value2,
                                double evt_wgt,
                                std::map<std::string,std::map<std::string,UBTH2Poly*>> hmap, 
                                std::string channel_namel, 
                                std::vector<std::string> fname, 
                                std::vector<double> wgts_genie) {

  auto iter = hmap.find(channel_namel);
  if (iter == hmap.end()) {
    LOG_CRITICAL() << "Can't find " << channel_namel << std::endl;
    throw std::exception();
  }
  std::map<std::string,UBTH2Poly*> this_map = iter->second;

  this_map["nominal"]->Fill(fill_value1, fill_value2, evt_wgt);

  for (size_t i = 0; i < fname.size(); i++) {

    //std::cout << i << ": Fill value: " << fill_value << ", weight: " << wgts_genie.at(i) << ", name is " << fname.at(i) << std::endl;

    this_map[fname.at(i)]->Fill(fill_value1, fill_value2, wgts_genie.at(i) * evt_wgt);

  }


} 

//___________________________________________________________________________________________________
void Main::Maker::FillBootstrap(double fill_value1,
                                double fill_value2,
                                double evt_wgt,
                                std::map<std::string,TH2D*> hmap_trkmom_genie_pm1_bs, 
                                std::vector<std::string> fname, 
                                std::vector<double> wgts_genie) {


  hmap_trkmom_genie_pm1_bs["nominal"]->Fill(fill_value1, fill_value2, evt_wgt);

  for (size_t i = 0; i < fname.size(); i++) {

    hmap_trkmom_genie_pm1_bs[fname.at(i)]->Fill(fill_value1, fill_value2, wgts_genie.at(i) * evt_wgt);

    //std::cout << "Fill value: " << fill_value << ", weight: " << wgts_genie_pm1.at(i) << std::endl;

  }

}




//______________________________________________________________________________________________________
void Main::Maker::FillBootstrap_test(double fill_value1,
                               double fill_value2,
                               double evt_wgt, 
                               std::map<std::string,TH2D*> bs_genie_pm1_reco_true_mom,
                               std::vector<std::string> fname,
                               std::vector<double> wgts_genie){
    
    bs_genie_pm1_reco_true_mom["nominal"]->Fill(fill_value1, fill_value2, evt_wgt);
    //std::cout<<"start printing out the genie parameters"<<std::endl;
    for(size_t i= 0; i< fname.size(); i++){
        //std::cout<<"GENIE_par.push_back(\""<<fname.at(i)<<"\");"<<std::endl;
        //std::cout<<"i= "<<i<<" fname "<<fname.at(i)<<" weight fac "<<wgts_genie.at(i)<<std::endl;
        bs_genie_pm1_reco_true_mom[fname.at(i)]->Fill(fill_value1, fill_value2, wgts_genie.at(i) * evt_wgt);
        //bs_genie_pm1_reco_true_mom[fname.at(i)]->Fill(fill_value1, fill_value2, wgts_genie.at(i) * evt_wgt);
    }
    //std::cout<<"end of print out genie parameters"<<std::endl;
}









//___________________________________________________________________________________________________
void Main::Maker::FillBootstrap(double fill_value1, // reco value x (costheta)
                                double fill_value2, // reco value y (momentum)
                                int m, // true bin m (costheta)
                                int n, // true bin n (momentum)
                                double evt_wgt,
                                std::map<std::string,std::vector<std::vector<TH2D*>>> bs_reco_per_true, 
                                std::vector<std::string> fname, 
                                std::vector<double> wgts) {


  bs_reco_per_true["nominal"][m][n]->Fill(fill_value1, fill_value2, evt_wgt);

  for (size_t i = 0; i < fname.size(); i++) {

    bs_reco_per_true[fname.at(i)][m][n]->Fill(fill_value1, fill_value2, wgts.at(i) * evt_wgt);

    //std::cout << "Fill value: " << fill_value << ", weight: " << wgts_genie_pm1.at(i) << std::endl;

  }

}


//___________________________________________________________________________________________________
void Main::Maker::FillBootstrap(double fill_value1, // reco value x (costheta)
                                double fill_value2, // reco value y (momentum)
                                int m, // true bin m (1 number, unrolled)
                                double evt_wgt,
                                std::map<std::string,std::vector<UBTH2Poly*>> bs_poly_reco_per_true, 
                                std::vector<std::string> fname, 
                                std::vector<double> wgts) {


  bs_poly_reco_per_true["nominal"][m]->Fill(fill_value1, fill_value2, evt_wgt);

  for (size_t i = 0; i < fname.size(); i++) {

    bs_poly_reco_per_true[fname.at(i)][m]->Fill(fill_value1, fill_value2, wgts.at(i) * evt_wgt);

    //std::cout << "Fill value: " << fill_value << ", weight: " << wgts_genie_pm1.at(i) << std::endl;

  }

}

//___________________________________________________________________________________________________
void Main::Maker::FillBootstrap(int m, // true bin m (1 number, unrolled)
                                int j, // reco bin i (1 number, unrolled)
                                double evt_wgt,
                                std::map<std::string,std::vector<std::vector<double>>> & bs_poly_reco_per_true, 
                                std::vector<std::string> fname, 
                                std::vector<double> wgts) {


  if (j < 0) j = 0; // Negative bins are overflows, and are all added to entry 0 of the vector

  bs_poly_reco_per_true["nominal"][m][j] += evt_wgt;

  for (size_t i = 0; i < fname.size(); i++) {

    bs_poly_reco_per_true[fname.at(i)][m][j] += wgts.at(i) * evt_wgt;

    //std::cout << "Fill value: " << fill_value << ", weight: " << wgts_genie_pm1.at(i) << std::endl;

  }

}




//___________________________________________________________________________________________________
void Main::Maker::AddPolyBins(UBTH2Poly * h) {

  // std::map<int, std::pair<int, int>> _exclusion_map;
  // _exclusion_map[0] = std::make_pair(2, 3);

  h->SetNBinsX(n_bins_double_mucostheta);

  for (int y = 0; y < n_bins_double_mumom; y++) {
    for (int x = 0; x < n_bins_double_mucostheta; x++) {

      auto it = _exclusion_map.find(x);
      if (it != _exclusion_map.end()) {
        if (y == it->second.first) {
          h->AddBin(bins_double_mucostheta[it->first], bins_double_mumom[it->second.first], bins_double_mucostheta[it->first+1], bins_double_mumom[it->second.second+1]);
          continue;
        } else if (y == it->second.second) {
          continue;
        }
      }
      h->AddBin(bins_double_mucostheta[x], bins_double_mumom[y], bins_double_mucostheta[x+1], bins_double_mumom[y+1]);
    }
  }
}



//___________________________________________________________________________________________________
void Main::Maker::AddPolyBins(BootstrapTH2DPoly h) {

  // std::map<int, std::pair<int, int>> _exclusion_map;
  // _exclusion_map[0] = std::make_pair(2, 3);

  // h->SetNBinsX(n_bins_double_mucostheta);

  for (int y = 0; y < n_bins_double_mumom; y++) {
    for (int x = 0; x < n_bins_double_mucostheta; x++) {

      auto it = _exclusion_map.find(x);
      if (it != _exclusion_map.end()) {
        if (y == it->second.first) {
          h.AddBin(bins_double_mucostheta[it->first], bins_double_mumom[it->second.first], bins_double_mucostheta[it->first+1], bins_double_mumom[it->second.second+1]);
          continue;
        } else if (y == it->second.second) {
          continue;
        }
      }
      h.AddBin(bins_double_mucostheta[x], bins_double_mumom[y], bins_double_mucostheta[x+1], bins_double_mumom[y+1]);
    }
  }
}
















void Main::Maker::MakeFile() 
{

	clock_t begin = clock();

  double n_signal = 0;


  system("mkdir -p output/");
  
  // CSV file for dqdx and track lenght values
  std::ofstream _csvfile;
  _csvfile.open ("./dqdx_trklen.csv", std::ofstream::out | std::ofstream::trunc);
  _csvfile << "dqdx,trklen,y" << std::endl;

   
  if (isdata) {
    LOG_NORMAL() << "Running on a data file." << std::endl;
  }
  else {
    LOG_NORMAL() << "Running on a MC file." << std::endl;
  }

  

  if (_scale_cosmics) {
    for (int i = 0; i < 10; i++) {
      std::cout << "****** Scaling cosmic background by " << _scale_factor_cosmic << "******" << std::endl;
    }
  }
  

  //*************************
  //* Starting ROOT application
  //*************************

  
  //TApplication* rootapp = new TApplication("ROOT Application",&argc, argv);
  gROOT->SetBatch(kTRUE);
  gROOT->ProcessLine("gErrorIgnoreLevel = 2001;"); // 1001: INFO, 2001: WARNINGS, 3001: ERRORS
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();


  //gROOT->ProcessLine(".x rootlogon.C");



  LOG_NORMAL() << "Opening output file with name " << fileoutn << std::endl;
  TFile *file_out = new TFile(fileoutn.c_str(),"RECREATE");
  if ( file_out->IsOpen() ) {
    LOG_NORMAL() << "File opened successfully." << std::endl;
  } else {
    LOG_CRITICAL() << "File not opened (maybe not found?). File: " << fileoutn << std::endl;
    exit(0);
  }
  
  string pattern = filen;
  
  
  
  //*************************
  //* Getting POTs
  //*************************
  
  double totalPOT = 0.;
  
  if (maxEntries > 0) evalPOT = false;
  
  if (evalPOT) {
    
    LOG_NORMAL() << " ~~~~~~~~~~~~~~ " << endl;
    LOG_NORMAL() << " |   Calculating POT" << endl;
    LOG_NORMAL() << " |" << endl;
    TChain *cpot;
    cpot = new TChain("UBXSec/pottree");
    cpot->Add(pattern.c_str());
    LOG_NORMAL() << " | Number of entries in the pot tree: " << cpot->GetEntries() << endl;
    Double_t pot;
    cpot->SetBranchAddress("pot", &pot);
    for (int potEntry = 0; potEntry < cpot->GetEntries(); potEntry++) {
      cpot->GetEntry(potEntry);
      totalPOT += pot;
    } // end loop entries
    LOG_NORMAL() << " | Total POT: " << totalPOT << endl;
    LOG_NORMAL() << " ~~~~~~~~~~~~~~ " << endl << endl;
  } // end if evalPOT
  else
    totalPOT = -1.;
  
  double pot_scaling = 1.;
  if (evalPOT) pot_scaling = targetPOT/totalPOT;
  
  
  TChain *chain_ubxsec;
  chain_ubxsec = new TChain("UBXSec/tree");
  chain_ubxsec->Add(pattern.c_str());
  
  LOG_NORMAL() << "Using file: " << pattern << endl;
  
  int Nfiles = chain_ubxsec->GetNtrees();
  LOG_NORMAL() << "Number of files: " << Nfiles << endl;
  
  int evts = chain_ubxsec -> GetEntries();
  LOG_NORMAL() << "Number of events used is: " << evts << endl;
  
  UBXSecEvent * t = new UBXSecEvent(chain_ubxsec);
  //ActivateBranches(t);

  _event_histo_1d = new UBXSecEventHisto1D();
  _event_histo_1d->InitializeBootstraps();

  _event_histo = new UBXSecEventHisto();
  _event_histo->InitializeBootstraps();
  // _event_histo->OpenFile(fileoutn + "mytest");

  _n_poly_bins = _event_histo->_n_poly_bins;
  LOG_NORMAL() << "Number of polybins: " << _n_poly_bins << std::endl;

  double nsignal = 0;


  double nsignal_qe = 0;
  double nsignal_res = 0;
  double nsignal_dis = 0;
  double nsignal_coh = 0;
  double nsignal_mec = 0;
  
  double signal_sel = 0;
  double bkg_anumu_sel = 0;
  double bkg_nue_sel = 0;
  double bkg_nc_sel = 0;
  double bkg_outfv_sel = 0;
  double bkg_cosmic_sel = 0;
  double bkg_cosmic_top_sel = 0;
  double bkg_ccother_sel = 0;
  double signal_sel_qe = 0;
  double signal_sel_res = 0;
  double signal_sel_dis = 0;
  double signal_sel_coh = 0;
  double signal_sel_mec = 0;
  
  int nEvtsWFlashInBeamSpill = 0;
  int nNumuCC = 0;
  
  double nue_cc_fv = 0;
  double nue_cc_selected = 0;
  double nue_cc_selected_total = 0;
  double nue_cc_selected_total_energy_range = 0;
  double nue_selected_total_energy_range = 0;
  double n_nue_electron = 0;
  double n_nue_proton = 0;
  double n_nue_pion = 0;

  int nSignalWMuonReco = 0;
  int nSignalMuonRecoVtxOk = 0;
  
  int nSignalFlashMatched = 0;
  
  int n_slc_nu_origin = 0;

  std::map<std::string, double> selected_events_percut;
  selected_events_percut["initial"] = 0.;
  selected_events_percut["beamflash"] = 0.;
  selected_events_percut["flash_match"] = 0.;
  selected_events_percut["flash_match_deltax"] = 0.;
  selected_events_percut["flash_match_deltaz"] = 0.;
  selected_events_percut["fiducial_volume"] = 0.;
  selected_events_percut["quality"] = 0.;
  selected_events_percut["mcs_length_quality"] = 0.;
  selected_events_percut["mip_consistency"] = 0.;
  selected_events_percut["ntrk2"] = 0.;
  selected_events_percut["pinCV"] = 0. ;
  selected_events_percut["minCol"] = 0.;
  selected_events_percut["chi2"] = 0. ;


  std::map<std::string, double> selected_signal_events_percut;
  selected_signal_events_percut["initial"] = 0.;
  selected_signal_events_percut["beamflash"] = 0.;
  selected_signal_events_percut["flash_match"] = 0.;
  selected_signal_events_percut["flash_match_deltax"] = 0.;
  selected_signal_events_percut["flash_match_deltaz"] = 0.;
  selected_signal_events_percut["fiducial_volume"] = 0.;
  selected_signal_events_percut["quality"] = 0.;
  selected_signal_events_percut["mcs_length_quality"] = 0.;
  selected_signal_events_percut["mip_consistency"] = 0.;
  selected_signal_events_percut["ntrk2"] = 0.;
  selected_signal_events_percut["pinCV"] = 0. ;
  selected_signal_events_percut["minCol"] = 0.;
  selected_signal_events_percut["chi2"] = 0. ;


  
  TTree* _shower_tree = new TTree("shower_tree", "shower_tree");
  double _s_nupdg, _s_track_pdg, _s_tpcobj_origin, _s_shower_length, _s_shower_phi, _s_shower_theta, _s_shower_openangle, _s_shower_startx, _s_shower_starty, _s_shower_startz, _s_flash_z;
  _shower_tree->Branch("s_nupdg", &_s_nupdg, "s_nupdg/D");
  _shower_tree->Branch("s_track_pdg", &_s_track_pdg, "s_track_pdg/D");
  _shower_tree->Branch("s_tpcobj_origin", &_s_tpcobj_origin, "s_tpcobj_origin/D");
  _shower_tree->Branch("s_shower_length", &_s_shower_length, "s_shower_length/D");
  _shower_tree->Branch("s_shower_phi", &_s_shower_phi, "s_shower_phi/D");
  _shower_tree->Branch("s_shower_theta", &_s_shower_theta, "s_shower_theta/D");
  _shower_tree->Branch("s_shower_openangle", &_s_shower_openangle, "s_shower_openangle/D");
  _shower_tree->Branch("s_shower_startx", &_s_shower_startx, "s_shower_startx/D");
  _shower_tree->Branch("s_shower_starty", &_s_shower_starty, "s_shower_starty/D");
  _shower_tree->Branch("s_shower_startz", &_s_shower_startz, "s_shower_startz/D");
  _shower_tree->Branch("s_flash_z", &_s_flash_z, "s_flash_z/D");


  TTree* _true_reco_tree = new TTree("true_reco_tree", "true_reco_tree");
  _true_reco_tree->Branch("mom_true", &_mom_true, "mom_true/D");
  _true_reco_tree->Branch("mom_mcs", &_mom_mcs, "mom_mcs/D");
  _true_reco_tree->Branch("pmom_true",&_pmom_true, "pmom_true/D");
  _true_reco_tree->Branch("pmom_reco",&_pmom_reco, "pmom_reco/D");
  _true_reco_tree->Branch("contained", &_contained, "contained/O");
  _true_reco_tree->Branch("selected", &_selected, "selected/O");
  _true_reco_tree->Branch("angle_true", &_angle_true, "angle_true/D");
  _true_reco_tree->Branch("angle_reco", &_angle_reco, "angle_reco/D");
  _true_reco_tree->Branch("event_weight", &_event_weight_fortree, "event_weight/D");
  _true_reco_tree->Branch("wgtsnames_genie_multisim", "std::vector<std::string>", &_wgtsnames_genie_multisim);
  _true_reco_tree->Branch("wgts_genie_multisim", "std::vector<double>", &_wgts_genie_multisim);
  _true_reco_tree->Branch("wgtsnames_extra_syst", "std::vector<std::string>", &_wgtsnames_extra_syst);
  _true_reco_tree->Branch("wgts_extra_syst", "std::vector<double>", &_wgts_extra_syst);
  _true_reco_tree->Branch("wgtsnames_flux_multisim", "std::vector<std::string>", &_wgtsnames_flux_multisim);
  _true_reco_tree->Branch("wgts_flux_multisim", "std::vector<double>", &_wgts_flux_multisim);


  

  //
  // Truth histograms stacked in interaction type - Selected
  //
  std::map<std::string,TH1D*> hmap_mctruth_nuenergy;
  hmap_mctruth_nuenergy["total"] = new TH1D("h_mctruth_nuenergy_total", ";True Neutrino Energy [GeV];Selected Events", 20, 0, 3);
  hmap_mctruth_nuenergy["qe"] = new TH1D("h_mctruth_nuenergy_qe", ";True Neutrino Energy [GeV];Selected Events", 20, 0, 3);
  hmap_mctruth_nuenergy["res"] = new TH1D("h_mctruth_nuenergy_res", ";True Neutrino Energy [GeV];Selected Events", 20, 0, 3);
  hmap_mctruth_nuenergy["dis"] = new TH1D("h_mctruth_nuenergy_dis", ";True Neutrino Energy [GeV];Selected Events", 20, 0, 3);
  hmap_mctruth_nuenergy["coh"] = new TH1D("h_mctruth_nuenergy_coh", ";True Neutrino Energy [GeV];Selected Events", 20, 0, 3);
  hmap_mctruth_nuenergy["mec"] = new TH1D("h_mctruth_nuenergy_mec", ";True Neutrino Energy [GeV];Selected Events", 20, 0, 3);
  hmap_mctruth_nuenergy["other"] = new TH1D("h_mctruth_nuenergy_other", ";True Neutrino Energy [GeV];Selected Events", 20, 0, 3);

  std::map<std::string,TH1D*> hmap_mctruth_mumom;
  hmap_mctruth_mumom["total"] = new TH1D("h_mctruth_mumom_total", ";True Muon Momentum [GeV];Selected Events", 20, 0.1, 2.5);
  hmap_mctruth_mumom["qe"] = new TH1D("h_mctruth_mumom_qe", ";True Muon Momentum [GeV];Selected Events", 20, 0.1, 2.5);
  hmap_mctruth_mumom["res"] = new TH1D("h_mctruth_mumom_res", ";True Muon Momentum [GeV];Selected Events", 20, 0.1, 2.5);
  hmap_mctruth_mumom["dis"] = new TH1D("h_mctruth_mumom_dis", ";True Muon Momentum [GeV];Selected Events", 20, 0.1, 2.5);
  hmap_mctruth_mumom["coh"] = new TH1D("h_mctruth_mumom_coh", ";True Muon Momentum [GeV];Selected Events", 20, 0.1, 2.5);
  hmap_mctruth_mumom["mec"] = new TH1D("h_mctruth_mumom_mec", ";True Muon Momentum [GeV];Selected Events", 20, 0.1, 2.5);
  hmap_mctruth_mumom["other"] = new TH1D("h_mctruth_mumom_other", ";True Muon Momentum [GeV];Selected Events", 20, 0.1, 2.5);

  std::map<std::string,TH1D*> hmap_mctruth_mucostheta;
  hmap_mctruth_mucostheta["total"] = new TH1D("h_mctruth_mucostheta_total", ";True Muon cos(#theta);Selected Events", 25, -1, 1);
  hmap_mctruth_mucostheta["qe"] = new TH1D("h_mctruth_mucostheta_qe", ";True Muon cos(#theta);Selected Events", 25, -1, 1);
  hmap_mctruth_mucostheta["res"] = new TH1D("h_mctruth_mucostheta_res", ";True Muon cos(#theta);Selected Events", 25, -1, 1);
  hmap_mctruth_mucostheta["dis"] = new TH1D("h_mctruth_mucostheta_dis", ";True Muon cos(#theta);Selected Events", 25, -1, 1);
  hmap_mctruth_mucostheta["coh"] = new TH1D("h_mctruth_mucostheta_coh", ";True Muon cos(#theta);Selected Events", 25, -1, 1);
  hmap_mctruth_mucostheta["mec"] = new TH1D("h_mctruth_mucostheta_mec", ";True Muon cos(#theta);Selected Events", 25, -1, 1);
  hmap_mctruth_mucostheta["other"] = new TH1D("h_mctruth_mucostheta_other", ";True Muon cos(#theta);Selected Events", 25, -1, 1);

  std::map<std::string,TH1D*> hmap_mctruth_muphi;
  hmap_mctruth_muphi["total"] = new TH1D("h_mctruth_muphi_total", ";True Muon #phi;Selected Events", 20, -3.15, 3.15);
  hmap_mctruth_muphi["qe"] = new TH1D("h_mctruth_muphi_qe", ";True Muon #phi;Selected Events", 20, -3.15, 3.15);
  hmap_mctruth_muphi["res"] = new TH1D("h_mctruth_muphi_res", ";True Muon #phi;Selected Events", 20, -3.15, 3.15);
  hmap_mctruth_muphi["dis"] = new TH1D("h_mctruth_muphi_dis", ";True Muon #phi;Selected Events", 20, -3.15, 3.15);
  hmap_mctruth_muphi["coh"] = new TH1D("h_mctruth_muphi_coh", ";True Muon #phi;Selected Events", 20, -3.15, 3.15);
  hmap_mctruth_muphi["mec"] = new TH1D("h_mctruth_muphi_mec", ";True Muon #phi;Selected Events", 20, -3.15, 3.15);
  hmap_mctruth_muphi["other"] = new TH1D("h_mctruth_muphi_other", ";True Muon #phi;Selected Events", 20, -3.15, 3.15);

  std::map<std::string,TH1D*> hmap_mctruth_chargedmult;
  hmap_mctruth_chargedmult["total"] = new TH1D("h_mctruth_chargedmult_total", ";True Charged Particle Multiplicity;Selected Events", 10, 0, 10);
  hmap_mctruth_chargedmult["qe"] = new TH1D("h_mctruth_chargedmult_qe", ";True Charged Particle Multiplicity;Selected Events", 10, 0, 10);
  hmap_mctruth_chargedmult["res"] = new TH1D("h_mctruth_chargedmult_res", ";True Charged Particle Multiplicity;Selected Events", 10, 0, 10);
  hmap_mctruth_chargedmult["dis"] = new TH1D("h_mctruth_chargedmult_dis", ";True Charged Particle Multiplicity;Selected Events", 10, 0, 10);
  hmap_mctruth_chargedmult["coh"] = new TH1D("h_mctruth_chargedmult_coh", ";True Charged Particle Multiplicity;Selected Events", 10, 0, 10);
  hmap_mctruth_chargedmult["mec"] = new TH1D("h_mctruth_chargedmult_mec", ";True Charged Particle Multiplicity;Selected Events", 10, 0, 10);
  hmap_mctruth_chargedmult["other"] = new TH1D("h_mctruth_chargedmult_other", ";True Charged Particle Multiplicity;Selected Events", 10, 0, 10);

  std::map<std::string,TH2D*> hmap_mctruth_mucostheta_mumom;
  hmap_mctruth_mucostheta_mumom["total"] = new TH2D("hmap_mctruth_mucostheta_mumom_total", ";True Muon cos(#theta);True Muon Momentum", 25, -1, 1, 20, 0, 2.5);
  hmap_mctruth_mucostheta_mumom["qe"] = new TH2D("hmap_mctruth_mucostheta_mumom_qe", ";True Muon cos(#theta);True Muon Momentum", 25, -1, 1, 20, 0, 2.5);
  hmap_mctruth_mucostheta_mumom["res"] = new TH2D("hmap_mctruth_mucostheta_mumom_res", ";True Muon cos(#theta);True Muon Momentum", 25, -1, 1, 20, 0, 2.5);
  hmap_mctruth_mucostheta_mumom["dis"] = new TH2D("hmap_mctruth_mucostheta_mumom_dis", ";True Muon cos(#theta);True Muon Momentum", 25, -1, 1, 20, 0, 2.5);
  hmap_mctruth_mucostheta_mumom["coh"] = new TH2D("hmap_mctruth_mucostheta_mumom_coh", ";True Muon cos(#theta);True Muon Momentum", 25, -1, 1, 20, 0, 2.5);
  hmap_mctruth_mucostheta_mumom["mec"] = new TH2D("hmap_mctruth_mucostheta_mumom_mec", ";True Muon cos(#theta);True Muon Momentum", 25, -1, 1, 20, 0, 2.5);
  hmap_mctruth_mucostheta_mumom["other"] = new TH2D("hmap_mctruth_mucostheta_mumom_other", ";True Muon cos(#theta);True Muon Momentum", 25, -1, 1, 20, 0, 2.5);

  std::map<std::string,TH1D*> hmap_mctruth_pmom;
  hmap_mctruth_pmom["total"] = new TH1D("h_mctruth_pmom_total", ";True Proton Momentum [GeV];Selected Events", 24, 0, 1.2);
  hmap_mctruth_pmom["qe"] = new TH1D("h_mctruth_pmom_qe", ";True Proton Momentum [GeV];Selected Events", 24, 0, 1.2);
  hmap_mctruth_pmom["res"] = new TH1D("h_mctruth_pmom_res", ";True Proton Momentum [GeV];Selected Events", 24, 0, 1.2);
  hmap_mctruth_pmom["dis"] = new TH1D("h_mctruth_pmom_dis", ";True Proton Momentum [GeV];Selected Events", 24, 0, 1.2);
  hmap_mctruth_pmom["coh"] = new TH1D("h_mctruth_pmom_coh", ";True Proton Momentum [GeV];Selected Events", 24, 0, 1.2);
  hmap_mctruth_pmom["mec"] = new TH1D("h_mctruth_pmom_mec", ";True Proton Momentum [GeV];Selected Events", 24, 0, 1.2);
  hmap_mctruth_pmom["other"] = new TH1D("h_mctruth_pmom_other", ";True Proton Momentum [GeV];Selected Events", 24, 0, 1.2);

  std::map<std::string,TH1D*> hmap_mctruth_pcostheta;
  hmap_mctruth_pcostheta["total"] = new TH1D("h_mctruth_pcostheta_total", ";True Proton cos(#theta);Selected Events", 25, -1, 1);
  hmap_mctruth_pcostheta["qe"] = new TH1D("h_mctruth_pcostheta_qe", ";True Proton cos(#theta);Selected Events", 25, -1, 1);
  hmap_mctruth_pcostheta["res"] = new TH1D("h_mctruth_pcostheta_res", ";True Proton cos(#theta);Selected Events", 25, -1, 1);
  hmap_mctruth_pcostheta["dis"] = new TH1D("h_mctruth_pcostheta_dis", ";True Proton cos(#theta);Selected Events", 25, -1, 1);
  hmap_mctruth_pcostheta["coh"] = new TH1D("h_mctruth_pcostheta_coh", ";True Proton cos(#theta);Selected Events", 25, -1, 1);
  hmap_mctruth_pcostheta["mec"] = new TH1D("h_mctruth_pcostheta_mec", ";True Proton cos(#theta);Selected Events", 25, -1, 1);
  hmap_mctruth_pcostheta["other"] = new TH1D("h_mctruth_pcostheta_other", ";True Proton cos(#theta);Selected Events", 25, -1, 1);

  std::map<std::string,TH1D*> hmap_mctruth_pphi;
  hmap_mctruth_pphi["total"] = new TH1D("h_mctruth_pphi_total", ";True Proton #phi [Rad];Selected Events", 20, -3.15, 3.15);
  hmap_mctruth_pphi["qe"] = new TH1D("h_mctruth_pphi_qe", ";True Proton #phi [Rad;Selected Events", 20, -3.15, 3.15);
  hmap_mctruth_pphi["res"] = new TH1D("h_mctruth_pphi_res", ";True Proton #phi [Rad];Selected Events", 20, -3.15, 3.15);
  hmap_mctruth_pphi["dis"] = new TH1D("h_mctruth_pphi_dis", ";True Proton #phi [Rad];Selected Events", 20, -3.15, 3.15);
  hmap_mctruth_pphi["coh"] = new TH1D("h_mctruth_pphi_coh", ";True Proton #phi [Rad];Selected Events", 20, -3.15, 3.15);
  hmap_mctruth_pphi["mec"] = new TH1D("h_mctruth_pphi_mec", ";True Proton #phi [Rad];Selected Events", 20, -3.15, 3.15);
  hmap_mctruth_pphi["other"] = new TH1D("h_mctruth_pphi_other", ";True Proton #phi [Rad];Selected Events", 20, -3.15, 3.15);

  std::map<std::string,TH1D*> hmap_mctruth_thetamup;
  hmap_mctruth_thetamup["total"] = new TH1D("h_mctruth_thetamup_total", ";#theta_{#mu p};Selected Events", 20, 0, 3.15);
  hmap_mctruth_thetamup["qe"] = new TH1D("h_mctruth_thetamup_qe", ";#theta_{#mu p};Selected Events", 20, 0, 3.15);
  hmap_mctruth_thetamup["res"] = new TH1D("h_mctruth_thetamup_res", ";#theta_{#mu p};Selected Events", 20, 0, 3.15);
  hmap_mctruth_thetamup["dis"] = new TH1D("h_mctruth_thetamup_dis", ";#theta_{#mu p};Selected Events", 20, 0, 3.15);
  hmap_mctruth_thetamup["coh"] = new TH1D("h_mctruth_thetamup_coh", ";#theta_{#mu p};Selected Events", 20, 0, 3.15);
  hmap_mctruth_thetamup["mec"] = new TH1D("h_mctruth_thetamup_mec", ";#theta_{#mu p};Selected Events", 20, 0, 3.15);
  hmap_mctruth_thetamup["other"] = new TH1D("h_mctruth_thetamup_other", ";#theta_{#mu p};Selected Events", 20, 0, 3.15);


  //
  // Truth histograms stacked in interaction type - Generated
  //
  std::map<std::string,TH1D*> hmap_mctruth_nuenergy_gen;
  hmap_mctruth_nuenergy_gen["total"] = new TH1D("h_mctruth_nuenergy_gen_total", ";True Neutrino Energy [GeV];Selected Events", 20, 0, 3);
  hmap_mctruth_nuenergy_gen["qe"] = new TH1D("h_mctruth_nuenergy_gen_qe", ";True Neutrino Energy [GeV];Selected Events", 20, 0, 3);
  hmap_mctruth_nuenergy_gen["res"] = new TH1D("h_mctruth_nuenergy_gen_res", ";True Neutrino Energy [GeV];Selected Events", 20, 0, 3);
  hmap_mctruth_nuenergy_gen["dis"] = new TH1D("h_mctruth_nuenergy_gen_dis", ";True Neutrino Energy [GeV];Selected Events", 20, 0, 3);
  hmap_mctruth_nuenergy_gen["coh"] = new TH1D("h_mctruth_nuenergy_gen_coh", ";True Neutrino Energy [GeV];Selected Events", 20, 0, 3);
  hmap_mctruth_nuenergy_gen["mec"] = new TH1D("h_mctruth_nuenergy_gen_mec", ";True Neutrino Energy [GeV];Selected Events", 20, 0, 3);
  hmap_mctruth_nuenergy_gen["other"] = new TH1D("h_mctruth_nuenergy_gen_other", ";True Neutrino Energy [GeV];Selected Events", 20, 0, 3);

  std::map<std::string,TH1D*> hmap_mctruth_mumom_gen;
  hmap_mctruth_mumom_gen["total"] = new TH1D("h_mctruth_mumom_gen_total", ";True Muon Momentum [GeV];Selected Events", 20, 0.1, 2.5);
  hmap_mctruth_mumom_gen["qe"] = new TH1D("h_mctruth_mumom_gen_qe", ";True Muon Momentum [GeV];Selected Events", 20, 0.1, 2.5);
  hmap_mctruth_mumom_gen["res"] = new TH1D("h_mctruth_mumom_gen_res", ";True Muon Momentum [GeV];Selected Events", 20, 0.1, 2.5);
  hmap_mctruth_mumom_gen["dis"] = new TH1D("h_mctruth_mumom_gen_dis", ";True Muon Momentum [GeV];Selected Events", 20, 0.1, 2.5);
  hmap_mctruth_mumom_gen["coh"] = new TH1D("h_mctruth_mumom_gen_coh", ";True Muon Momentum [GeV];Selected Events", 20, 0.1, 2.5);
  hmap_mctruth_mumom_gen["mec"] = new TH1D("h_mctruth_mumom_gen_mec", ";True Muon Momentum [GeV];Selected Events", 20, 0.1, 2.5);
  hmap_mctruth_mumom_gen["other"] = new TH1D("h_mctruth_mumom_gen_other", ";True Muon Momentum [GeV];Selected Events", 20, 0.1, 2.5);

  std::map<std::string,TH1D*> hmap_mctruth_mucostheta_gen;
  hmap_mctruth_mucostheta_gen["total"] = new TH1D("h_mctruth_mucostheta_gen_total", ";True Muon cos(#theta);Selected Events", 25, -1, 1);
  hmap_mctruth_mucostheta_gen["qe"] = new TH1D("h_mctruth_mucostheta_gen_qe", ";True Muon cos(#theta);Selected Events", 25, -1, 1);
  hmap_mctruth_mucostheta_gen["res"] = new TH1D("h_mctruth_mucostheta_gen_res", ";True Muon cos(#theta);Selected Events", 25, -1, 1);
  hmap_mctruth_mucostheta_gen["dis"] = new TH1D("h_mctruth_mucostheta_gen_dis", ";True Muon cos(#theta);Selected Events", 25, -1, 1);
  hmap_mctruth_mucostheta_gen["coh"] = new TH1D("h_mctruth_mucostheta_gen_coh", ";True Muon cos(#theta);Selected Events", 25, -1, 1);
  hmap_mctruth_mucostheta_gen["mec"] = new TH1D("h_mctruth_mucostheta_gen_mec", ";True Muon cos(#theta);Selected Events", 25, -1, 1);
  hmap_mctruth_mucostheta_gen["other"] = new TH1D("h_mctruth_mucostheta_gen_other", ";True Muon cos(#theta);Selected Events", 25, -1, 1);

  std::map<std::string,TH1D*> hmap_mctruth_muphi_gen;
  hmap_mctruth_muphi_gen["total"] = new TH1D("h_mctruth_muphi_gen_total", ";True Muon #phi;Selected Events", 20, -3.15, 3.15);
  hmap_mctruth_muphi_gen["qe"] = new TH1D("h_mctruth_muphi_gen_qe", ";True Muon #phi;Selected Events", 20, -3.15, 3.15);
  hmap_mctruth_muphi_gen["res"] = new TH1D("h_mctruth_muphi_gen_res", ";True Muon #phi;Selected Events", 20, -3.15, 3.15);
  hmap_mctruth_muphi_gen["dis"] = new TH1D("h_mctruth_muphi_gen_dis", ";True Muon #phi;Selected Events", 20, -3.15, 3.15);
  hmap_mctruth_muphi_gen["coh"] = new TH1D("h_mctruth_muphi_gen_coh", ";True Muon #phi;Selected Events", 20, -3.15, 3.15);
  hmap_mctruth_muphi_gen["mec"] = new TH1D("h_mctruth_muphi_gen_mec", ";True Muon #phi;Selected Events", 20, -3.15, 3.15);
  hmap_mctruth_muphi_gen["other"] = new TH1D("h_mctruth_muphi_gen_other", ";True Muon #phi;Selected Events", 20, -3.15, 3.15);

  std::map<std::string,TH1D*> hmap_mctruth_chargedmult_gen;
  hmap_mctruth_chargedmult_gen["total"] = new TH1D("h_mctruth_chargedmult_gen_total", ";True Charged Particle Multiplicity;Selected Events", 10, 0, 10);
  hmap_mctruth_chargedmult_gen["qe"] = new TH1D("h_mctruth_chargedmult_gen_qe", ";True Charged Particle Multiplicity;Selected Events", 10, 0, 10);
  hmap_mctruth_chargedmult_gen["res"] = new TH1D("h_mctruth_chargedmult_gen_res", ";True Charged Particle Multiplicity;Selected Events", 10, 0, 10);
  hmap_mctruth_chargedmult_gen["dis"] = new TH1D("h_mctruth_chargedmult_gen_dis", ";True Charged Particle Multiplicity;Selected Events", 10, 0, 10);
  hmap_mctruth_chargedmult_gen["coh"] = new TH1D("h_mctruth_chargedmult_gen_coh", ";True Charged Particle Multiplicity;Selected Events", 10, 0, 10);
  hmap_mctruth_chargedmult_gen["mec"] = new TH1D("h_mctruth_chargedmult_gen_mec", ";True Charged Particle Multiplicity;Selected Events", 10, 0, 10);
  hmap_mctruth_chargedmult_gen["other"] = new TH1D("h_mctruth_chargedmult_gen_other", ";True Charged Particle Multiplicity;Selected Events", 10, 0, 10);

  std::map<std::string,TH2D*> hmap_mctruth_mucostheta_mumom_gen;
  hmap_mctruth_mucostheta_mumom_gen["total"] = new TH2D("hmap_mctruth_mucostheta_mumom_gen_total", ";True Muon cos(#theta);True Muon Momentum [GeV]", 25, -1, 1, 20, 0, 2.5);
  hmap_mctruth_mucostheta_mumom_gen["qe"] = new TH2D("hmap_mctruth_mucostheta_mumom_gen_qe", ";True Muon cos(#theta);True Muon Momentum [GeV]", 25, -1, 1, 20, 0, 2.5);
  hmap_mctruth_mucostheta_mumom_gen["res"] = new TH2D("hmap_mctruth_mucostheta_mumom_gen_res", ";True Muon cos(#theta);True Muon Momentum [GeV]", 25, -1, 1, 20, 0, 2.5);
  hmap_mctruth_mucostheta_mumom_gen["dis"] = new TH2D("hmap_mctruth_mucostheta_mumom_gen_dis", ";True Muon cos(#theta);True Muon Momentum [GeV]", 25, -1, 1, 20, 0, 2.5);
  hmap_mctruth_mucostheta_mumom_gen["coh"] = new TH2D("hmap_mctruth_mucostheta_mumom_gen_coh", ";True Muon cos(#theta);True Muon Momentum [GeV]", 25, -1, 1, 20, 0, 2.5);
  hmap_mctruth_mucostheta_mumom_gen["mec"] = new TH2D("hmap_mctruth_mucostheta_mumom_gen_mec", ";True Muon cos(#theta);True Muon Momentum [GeV]", 25, -1, 1, 20, 0, 2.5);
  hmap_mctruth_mucostheta_mumom_gen["other"] = new TH2D("hmap_mctruth_mucostheta_mumom_gen_other", ";True Muon cos(#theta);True Muon Momentum [GeV]", 25, -1, 1, 20, 0, 2.5);

  //-----------------------------------------------------------------------------------------------------------------------------
  //
  // True v.s. reco histograms for constructing smearing matrices->These histograms have been redefined as Bootstrap
  //
  /*std::map<std::string,TH2D*> bs_genie_pm1_true_reco_mom;
  bs_genie_pm1_true_reco_mom["nominal"] = new TH2D("bs_genie_pm1_true_reco_mom_nominal", ";Muon Momentum (Truth) [GeV]; Muon Momentum [GeV]", n_bins_mumom, bins_mumom, n_bins_mumom, bins_mumom);
  
  std::map<std::string,TH2D*> bs_genie_pm1_true_reco_pmom;
  bs_genie_pm1_true_reco_pmom["nominal"] = new TH2D("bs_genie_pm1_true_reco_pmom_nominal", ";Proton Momentum (Truth) [GeV]; Proton Momentum [GeV]", n_bins_pmom, bins_pmom, n_bins_pmom, bins_pmom);
  
  std::map<std::string,TH2D*> bs_genie_pm1_true_reco_muangle;
  bs_genie_pm1_true_reco_muangle["nominal"] = new TH2D("bs_genie_pm1_true_reco_muangle_nominal", ";Muon Angle (Truth); Muon Angle", n_bins_mucostheta, bins_mucostheta, n_bins_mucostheta, bins_mucostheta);
  
  std::map<std::string,TH2D*> bs_genie_pm1_true_reco_pangle;
  bs_genie_pm1_true_reco_pangle["nominal"] = new TH2D("bs_genie_pm1_true_reco_pangle_nominal", ";Proton Angle (Truth); Proton Angle", n_bins_pcostheta, bins_pcostheta, n_bins_pcostheta, bins_pcostheta);
  

  std::map<std::string,TH2D*> bs_genie_pm1_true_reco_thetamup;
  bs_genie_pm1_true_reco_thetamup["nominal"] = new TH2D("bs_genie_pm1_true_reco_thetamup_nominal", ";#theta_{#mu p} (Truth); #theta_{#mu p}", n_bins_muptheta, bins_muptheta, n_bins_muptheta, bins_muptheta);
  */




  
  TH1D* h_eff_num = new TH1D("h_eff_num", "h_eff_num", 15, 0, 3);
  TH1D* h_eff_den = new TH1D("h_eff_den", "h_eff_den", 15, 0, 3);
  TEfficiency* pEff = new TEfficiency("eff",";Neutrino Energy (truth) [GeV];Efficiency",6, 0, 4);


  // Efficiency - GENIE pm1sigma
  BootstrapTH1D bs_genie_pm1_eff_mumom_num("bs_genie_pm1_eff_mumom_num", "bs_genie_pm1_eff_mumom_num_title", n_bins_mumom, bins_mumom);
  BootstrapTH1D bs_genie_pm1_eff_mumom_den("bs_genie_pm1_eff_mumom_den", "bs_genie_pm1_eff_mumom_den_title", n_bins_mumom, bins_mumom);

  BootstrapTH1D bs_genie_pm1_eff_pmom_num("bs_genie_pm1_eff_pmom_num", "bs_genie_pm1_eff_pmom_num_title", n_bins_pmom, bins_pmom);
  BootstrapTH1D bs_genie_pm1_eff_pmom_den("bs_genie_pm1_eff_pmom_den", "bs_genie_pm1_eff_pmom_den_title", n_bins_pmom, bins_pmom);

  BootstrapTH1D bs_genie_pm1_eff_muangle_num("bs_genie_pm1_eff_muangle_num", "bs_genie_pm1_eff_muangle_num_title", n_bins_mucostheta, bins_mucostheta);
  BootstrapTH1D bs_genie_pm1_eff_muangle_den("bs_genie_pm1_eff_muangle_den", "bs_genie_pm1_eff_muangle_den_title", n_bins_mucostheta, bins_mucostheta);

  BootstrapTH1D bs_genie_pm1_eff_pangle_num("bs_genie_pm1_eff_pangle_num", "bs_genie_pm1_eff_pangle_num_title", n_bins_pcostheta, bins_pcostheta);
  BootstrapTH1D bs_genie_pm1_eff_pangle_den("bs_genie_pm1_eff_pangle_den", "bs_genie_pm1_eff_pangle_den_title", n_bins_pcostheta, bins_pcostheta);

  BootstrapTH1D bs_genie_pm1_eff_thetamup_num("bs_genie_pm1_eff_thetamup_num", "bs_genie_pm1_eff_thetamup_num_title", n_bins_muptheta, bins_muptheta);
  BootstrapTH1D bs_genie_pm1_eff_thetamup_den("bs_genie_pm1_eff_thetamup_den", "bs_genie_pm1_eff_thetamup_den_title", n_bins_muptheta, bins_muptheta);

  //-------------------------------------------------------------------------------------------------------------------------------------
  

  TH1D* h_eff_mult_num = new TH1D("h_eff_mult_num", "h_eff_mult_num", 20, 0, 20);
  TH1D* h_eff_mult_den = new TH1D("h_eff_mult_den", "h_eff_mult_den", 20, 0, 20);
  TH1D* h_eff_mult_ch_num = new TH1D("h_eff_mult_ch_num", "h_eff_mult_ch_num", 10, 0, 15);
  TH1D* h_eff_mult_ch_den = new TH1D("h_eff_mult_ch_den", "h_eff_mult_ch_den", 10, 0, 15);

  TH1D* h_eff_allprotons_den=new TH1D("h_eff_allprotons_den", "h_eff_allprotons_den", 30, 0.3, 1.5); 
  TH1D* h_eff_reintprotons_den=new TH1D("h_eff_reintprotons_den", "h_eff_reintprotons_den", 30, 0.3, 1.5); 
  TH1D* h_eff_containedprotons_den=new TH1D("h_eff_containedprotons_den", "h_eff_containedprotons_den", 30, 0.3, 1.5); 
  
  TH1D* h_eff_allprotons_num=new TH1D("h_eff_allprotons_num", "h_eff_allprotons_num", 30, 0.3, 1.5); 
  TH1D* h_eff_reintprotons_num=new TH1D("h_eff_reintprotons_num", "h_eff_reintprotons_num", 30, 0.3, 1.5); 
  TH1D* h_eff_containedprotons_num=new TH1D("h_eff_containedprotons_num", "h_eff_containedprotons_num", 30, 0.3, 1.5); 

  TH1D* h_nhits_lowmu=new TH1D("h_nhit_lowmu", "h_nhit_lowmu", 200, -0.5, 199.5);
  TH1D* h_nhits_lowmu_uplane=new TH1D("h_nhit_lowmu_uplane", "h_nhit_lowmu_uplane", 200, -0.5, 199.5);
  TH1D* h_nhits_lowmu_vplane=new TH1D("h_nhit_lowmu_vplane", "h_nhit_lowmu_vplane", 200, -0.5, 199.5);
  TH2D* h_trklen_truevsreco_lowmu=new TH2D("h_trklen_truevsreco_lowmu", "h_trklen_truevsreco_lowmu", 100, 0.0, 100.0, 100, 0.0, 100.0);
  h_trklen_truevsreco_lowmu->GetXaxis()->SetTitle("True length of Muon candidate [cm]");
  h_trklen_truevsreco_lowmu->GetYaxis()->SetTitle("Reco length of Muon candidate [cm]");
  TH1D* h_vtx_resolution_lowmu= new TH1D("h_vtx_resolution_lowmu", "h_vtx_resolution_lowmu", 100, 0.0, 100.0);
  TH1D* h_costhetay_lowmu=new TH1D("h_costhetay_lowmu", "h_costhetay_lowmu", 100, -1.0, 1.0); 
  
  TH1D* h_nhits_firstbin=new TH1D("h_nhit_firstbin", "h_nhit_firstbin", 200, -0.5, 199.5);
  TH1D* h_nhits_firstbin_uplane=new TH1D("h_nhit_firstbin_uplane", "h_nhit_firstbin_uplane", 100, -0.5, 199.5);
  TH1D* h_nhits_firstbin_vplane=new TH1D("h_nhit_firstbin_vplane", "h_nhit_firstbin_vplane", 100, -0.5, 199.5);
  TH2D* h_trklen_truevsreco_firstbin=new TH2D("h_trklen_truevsreco_firstbin", "h_trklen_truevsreco_firstbin", 100, 0.0, 100.0, 100, 0.0, 100.0);
  h_trklen_truevsreco_firstbin->GetXaxis()->SetTitle("True length of Muon candidate [cm]");
  h_trklen_truevsreco_firstbin->GetYaxis()->SetTitle("Reco length of Muon candidate [cm]");
  TH1D* h_vtx_resolution_firstbin= new TH1D("h_vtx_resolution_firstbin", "h_vtx_resolution_firstbin", 100, 0.0, 100.0);
  TH1D* h_costhetay_firstbin=new TH1D("h_costhetay_firstbin", "h_costhetay_firstbin", 100, -1.0, 1.0); 
 
  TH1D* h_eff_qe_num = new TH1D("h_eff_qe_num", "h_eff_qe_num", 15, 0, 3);
  TH1D* h_eff_qe_den = new TH1D("h_eff_qe_den", "h_eff_qe_den", 15, 0, 3);
  TH1D* h_eff_res_num = new TH1D("h_eff_res_num", "h_eff_res_num", 15, 0, 3);
  TH1D* h_eff_res_den = new TH1D("h_eff_res_den", "h_eff_res_den", 15, 0, 3);
  TH1D* h_eff_dis_num = new TH1D("h_eff_dis_num", "h_eff_dis_num", 15, 0, 3);
  TH1D* h_eff_dis_den = new TH1D("h_eff_dis_den", "h_eff_dis_den", 15, 0, 3);
  TH1D* h_eff_coh_num = new TH1D("h_eff_coh_num", "h_eff_coh_num", 15, 0, 3);
  TH1D* h_eff_coh_den = new TH1D("h_eff_coh_den", "h_eff_coh_den", 15, 0, 3);
  TH1D* h_eff_mec_num = new TH1D("h_eff_mec_num", "h_eff_mec_num", 15, 0, 3);
  TH1D* h_eff_mec_den = new TH1D("h_eff_mec_den", "h_eff_mec_den", 15, 0, 3);

  TH1D* h_truth_xsec_mumom = new TH1D("h_truth_xsec_mumom", "h_truth_xsec_mumom", n_bins_mumom, bins_mumom);
  TH1D* h_truth_xsec_muangle = new TH1D("h_truth_xsec_muangle", "h_truth_xsec_muangle", n_bins_mucostheta, bins_mucostheta);

  TH1D* h_truth_xsec_pmom  = new TH1D("h_truth_xsec_pmom ", "h_truth_xsec_pmom ", n_bins_pmom, bins_pmom);
  TH1D* h_truth_xsec_pangle = new TH1D("h_truth_xsec_pangle ", "h_truth_xsec_pangle ", n_bins_pcostheta, bins_pcostheta);

  TH1D* h_truth_xsec_thetamup = new TH1D("h_truth_xsec_thetamup", "h_truth_xsec_thetamup", n_bins_muptheta, bins_muptheta);

  TH1D* h_nue_selected_energy = new TH1D("h_nue_selected_energy", ";True Neutrino Energy [GeV];#nu_{e} Selected Events", 100, 0, 1.5);


  TH1D* h_true_nu_eng_beforesel = new TH1D("h_true_nu_eng_beforesel", ";True Neutrino Energy [GeV];Events", 200, 0, 3);
  TH1D* h_true_nu_eng_afterflash = new TH1D("h_true_nu_eng_afterflash", ";True Neutrino Energy [GeV];Events", 200, 0, 3);
  TH1D* h_true_nu_eng_aftersel = new TH1D("h_true_nu_eng_aftersel", ";True Neutrino Energy [GeV];Events", 200, 0, 3);

  
  
  TH1D* h_chi2 = new TH1D("h_chi2", "h_chi2", 50, 0, 50);
  TH1D* h_flsTime = new TH1D("h_flsTime", ";Flash time w.r.t. trigger [#mus];Flashes", 125, 0, 25);
  TH1D* h_flsTime_wcut = new TH1D("h_flsTime_wcut", ";Flash time w.r.t. trigger [#mus];Flashes (> 50PE)", 500, 0, 25);
  TH1D* h_flsTime_wcut_2 = new TH1D("h_flsTime_wcut_2", "(2);Flash time w.r.t. trigger [#mus];Flashes (> 50PE)", 125, 0, 25);
  TH1D* h_flsTime_wcut_3 = new TH1D("h_flsTime_wcut_3", "(3);Flash time w.r.t. trigger [#mus];Flashes (> 50PE)", 125, 0, 25);
  TH1D* h_flsTime_wcut_4 = new TH1D("h_flsTime_wcut_4", "(4);Flash time w.r.t. trigger [#mus];Flashes (> 50PE)", 125, 0, 25);
  TH1D* h_flsTime_wcut_5 = new TH1D("h_flsTime_wcut_5", "(5);Flash time w.r.t. trigger [#mus];Flashes (> 50PE)", 125, 0, 25);
  TH1D* h_flsTime_wcut_6 = new TH1D("h_flsTime_wcut_6", "(6);Flash time w.r.t. trigger [#mus];Flashes (> 50PE)", 125, 0, 25);
  TH1D* h_flsTime_wcut_7 = new TH1D("h_flsTime_wcut_7", "(7);Flash time w.r.t. trigger [#mus];Flashes (> 50PE)", 125, 0, 25);
  TH1D* h_flsTime_wcut_8 = new TH1D("h_flsTime_wcut_8", "(8);Flash time w.r.t. trigger [#mus];Flashes (> 50PE)", 125, 0, 25);
  h_flsTime->Sumw2(); h_flsTime_wcut->Sumw2(); h_flsTime_wcut_2->Sumw2(); h_flsTime_wcut_3->Sumw2(); h_flsTime_wcut_4->Sumw2(); h_flsTime_wcut_5->Sumw2(); h_flsTime_wcut_6->Sumw2(); h_flsTime_wcut_7->Sumw2(); h_flsTime_wcut_8->Sumw2();

  TH1D* h_flsPe_wcut = new TH1D("h_flsPe_wcut", ";Flash PE;Flashes (> 50PE)", 700, 0, 50000);
  TH2D* h_flsTime_flsPe_wcut = new TH2D("h_flsTime_flsPe_wcut", "Flashes (> 50PE);Flash time w.r.t. trigger [#mus];Flash PE", 125, 0, 25, 700, 0, 50000);
  
  TH1D* h_deltax = new TH1D("h_deltax", "(4);QLL X - TPC X [cm];;", 500, -200,200);
  TH2D* h_deltax_2d = new TH2D("h_deltax_2d", "(4);QLL X [cm];TPC X [cm]", 70, -100,350, 70, -100,350);
  TH1D* h_deltaz_4 = new TH1D("h_deltaz_4", "(4);Delta z [cm];", 100, -200,200);
  TH1D* h_deltaz_6 = new TH1D("h_deltaz_6", "(6);Delta z [cm];", 100, -200,200);

  TH1D* h_nslices = new TH1D("h_nslices", ";Number of slices per event;Entries per bin", 15, 0, 15);
  TH1D* h_vtx_resolution = new TH1D("h_nslh_vtx_resolutionices", ";Vertex resolution (2D) [cm];Entries per bin", 300, 0, 500);
  
  TH2D* h_frac_diff = new TH2D("h_frac_diff", ";PMT ID; Fractional difference", 32, 0, 32, 80, -2, 2);
  TH2D* h_frac_diff_others = new TH2D("h_frac_diff", ";PMT ID; Fractional difference", 32, 0, 32, 80, -2, 2);
  double hypo_spec_x[32], hypo_spec_y[32];
  double meas_spec_x[32], meas_spec_y[32];
  double numc_spec_x[32], numc_spec_y[32];
  
  TH1D* h_xdiff = new TH1D("h_xdiff", "h_xdiff", 1000, -100,100);
  TH1D* h_xdiff_others = new TH1D("h_xdiff_others", "h_xdiff_others", 1000, -100,100);
  TH1D* h_zdiff = new TH1D("h_zdiff", "h_zdiff", 1000, 0,1000);
  TH1D* h_zdiff_others = new TH1D("h_zdiff_others", "h_zdiff_others", 1000, 0,1000);
  // Before selection
  std::map<std::string,TH1D*> hmap_xdiff_b;
  hmap_xdiff_b["total"] = new TH1D("h_xdiff_total_b", ";QLL x - TPC x [cm];", 80, -200,200);//  100, 0,22
  hmap_xdiff_b["signal"] = new TH1D("h_xdiff_signal_b", ";QLL x - TPC x [cm];", 80, -200,200);
  hmap_xdiff_b["background"] = new TH1D("h_xdiff_background_b", ";QLL x - TPC x [cm];", 80, -200,200);
  std::map<std::string,TH1D*> hmap_zdiff_b;
  hmap_zdiff_b["total"] = new TH1D("h_zdiff_total_b", ";Hypo z - Flash z [cm];", 160, -400,400);
  hmap_zdiff_b["signal"] = new TH1D("h_zdiff_signal_b", ";Hypo z - Flash z [cm];", 160, -400,400);
  hmap_zdiff_b["background"] = new TH1D("h_zdiff_background_b", ";Hypo z - Flash z [cm];", 160, -400,400);
  // After selection
  std::map<std::string,TH1D*> hmap_xdiff;
  hmap_xdiff["total"] = new TH1D("h_xdiff_total", ";QLL x - TPC x [cm];", 80, -200,200);
  hmap_xdiff["signal"] = new TH1D("h_xdiff_signal", ";QLL x - TPC x [cm];", 80, -200,200);
  hmap_xdiff["background"] = new TH1D("h_xdiff_background", ";QLL x - TPC x [cm];", 80, -200,200);
  std::map<std::string,TH1D*> hmap_zdiff;
  hmap_zdiff["total"] = new TH1D("h_zdiff_total", ";Hypo z - Flash z [cm];", 160, -400,400);
  hmap_zdiff["signal"] = new TH1D("h_zdiff_signal", ";Hypo z - Flash z [cm];", 160, -400,400);
  hmap_zdiff["background"] = new TH1D("h_zdiff_background", ";Hypo z - Flash z [cm];", 160, -400,400);
  std::map<std::string,TH1D*> hmap_pediff;
  hmap_pediff["total"] = new TH1D("h_pediff_total", ";Hypo PE - Flash PE [p.e.];", 160, -400,400);
  hmap_pediff["signal"] = new TH1D("h_pediff_signal", ";Hypo PE - Flash PE [p.e.];", 160, -400,400);
  hmap_pediff["background"] = new TH1D("h_pediff_background", ";Hypo PE - Flash PE [p.e.];", 160, -600,600);

  std::map<std::string,TH1D*> hmap_vtxcheck_angle;
  hmap_vtxcheck_angle["total"] = new TH1D("h_vtxcheck_angle_total", ";Angle [rad];Entries per bin", 80, 0, 4);
  hmap_vtxcheck_angle["signal"] = new TH1D("h_vtxcheck_angle_signal", ";Angle [rad];Entries per bin", 80, 0, 4);
  hmap_vtxcheck_angle["background"] = new TH1D("h_vtxcheck_angle_background", ";Angle [rad];Entries per bin", 80, 0, 4);
  
  TH1D* h_vtxcheck_angle_good = new TH1D("h_vtxcheck_angle_good", ";Angle [rad];Entries per bin", 100, 0, 4);
  TH1D* h_vtxcheck_angle_bad  = new TH1D("h_vtxcheck_angle_bad",  ";Angle [rad];Entries per bin",  100, 0, 4);
  
  std::map<std::string,TH1D*> hmap_residuals_std;
  hmap_residuals_std["total"] = new TH1D("h_residuals_std_total", ";#sigma_{r_{i}};Entries per bin", 40, 0, 10);
  hmap_residuals_std["signal"] = new TH1D("h_residuals_std_signal", ";Angle [rad];Entries per bin", 40, 0, 10);
  hmap_residuals_std["background"] = new TH1D("h_residuals_std_background", ";Angle [rad];Entries per bin", 40, 0, 10);
  std::map<std::string,TH1D*> hmap_residuals_mean;
  hmap_residuals_mean["total"] = new TH1D("h_residuals_mean_total", ";<r_{i}>;Entries per bin", 40, -5, 5);
  hmap_residuals_mean["signal"] = new TH1D("h_residuals_mean_signal", ";<r_{i}>;Entries per bin", 40, -5, 5);
  hmap_residuals_mean["background"] = new TH1D("h_residuals_mean_background", ";<r_{i}>;Entries per bin", 40, -5, 5);
  std::map<std::string,TH1D*> hmap_perc_used_hits;
  hmap_perc_used_hits["total"] = new TH1D("h_perc_used_hits_total", ";Fraction of used hits in cluster;Entries per bin", 30, 0, 1);
  hmap_perc_used_hits["signal"] = new TH1D("h_perc_used_hits_signal", ";Fraction of used hits in cluster;Entries per bin", 30, 0, 1);
  hmap_perc_used_hits["background"] = new TH1D("h_perc_used_hits_background", ";Fraction of used hits in cluster;Entries per bin", 30, 0, 1);

  std::map<std::string,TH1D*> hmap_mom_mcs_length;
  hmap_mom_mcs_length["total"] = new TH1D("h_mom_mcs_length_total", ";(MCS - Length) Reconstructed Momentum [GeV];Entries per bin", 20, -0.5, 1.5);
  hmap_mom_mcs_length["signal"] = new TH1D("h_mom_mcs_length_signal", ";(MCS - Length) Reconstructed Momentum [GeV];Entries per bin", 20, -0.5, 1.5);
  hmap_mom_mcs_length["background"] = new TH1D("h_mom_mcs_length_background", ";(MCS - Length) Reconstructed Momentum [GeV];Entries per bin", 20, -0.5, 1.5);

  TH1D* h_muon_track_eff  = new TH1D("h_muon_track_eff",  ";Muon track efficiency;Entries per bin",  100, 0, 1);
  TH1D* h_muon_track_pur  = new TH1D("h_muon_track_pur",  ";Muon track purity;Entries per bin",  100, 0, 1);
  
  TH1D* h_mueff_num = new TH1D("h_mueff_num", "h_mueff_num", 30, 0, 2);
  TH1D* h_mueff_2_num = new TH1D("h_mueff_2_num", "h_mueff_2_num", 30, 0, 2);
  TH1D* h_mueff_den = new TH1D("h_mueff_den", "h_mueff_den", 30, 0, 2);
  TH1D* h_mueff_angle_num = new TH1D("h_mueff_angle_num", "h_mueff_num", 15, -1, 1);
  TH1D* h_mueff_angle_den = new TH1D("h_mueff_angle_den", "h_mueff_den", 15, -1, 1);

  TH2D* h_mu_eff_mom = new TH2D("h_mu_eff_mom", ";True Muon Momentum [GeV]; Efficiency", 50, 0, 2, 20, 0, 1);
  TH2D* h_mu_pur_mom = new TH2D("h_mu_pur_mom", ";True Muon Momentum [GeV]; Purity", 50, 0, 2, 20, 0, 1);
  TH2D* h_mu_eff_mom_sel = new TH2D("h_mu_eff_mom_sel", "After Selection;True Muon Momentum [GeV]; Efficiency", 50, 0, 2, 20, 0, 1);
  
  TH2D* h_mumom_nue = new TH2D("h_mumom_nue", ";True Neutrino Energy [GeV]; True Muon Momentum [GeV]", 50, 0, 2, 50, 0, 4);
  
  TH1D* h_acpt_tagged  = new TH1D("h_acpt_tagged",  ";Tagged TPC Objects;Entries per bin",  10, 0, 10);
  
  TH1D* h_slice_origin = new TH1D("h_slice_origin",  ";;",  3, -0.5, 2.5);
  
  TH1D* h_slice_npfp = new TH1D("h_slice_npfp",  ";npfp;",  10, 0, 10);
  TH1D* h_slice_npfp_others = new TH1D("h_slice_npfp_others",  ";npfp;",  10, 0, 10);
  
  TH1D* h_slice_ntrack = new TH1D("h_slice_ntrack",  ";npfp;",  10, 0, 10);
  TH1D* h_slice_ntrack_others = new TH1D("h_slice_ntrack_others",  ";npfp;",  10, 0, 10);
  
  TH1D* h_fm_score = new TH1D("h_fm_score",  ";fm score;",  500, 0, 10);
  TH1D* h_fm_score_others = new TH1D("h_fm_score_other",  ";fm score;",  500, 0, 10);
  TH2D* h_fm_score_pe = new TH2D("h_fm_score_pe",  ";fm score;Reco PE",  500, 0, 10, 500, 0, 2000);
  
  TH1D* h_n_slc_flsmatch = new TH1D("h_n_slc_flsmatch",  ";n slices flash matched per event;",  10, 0, 10);



  std::map<std::string,TH1D*> hmap_trklen;
  hmap_trklen["total"] = new TH1D("h_trklen_total", "; Track length;", 30, 0, 700);
  hmap_trklen["signal"] = new TH1D("h_trklen_signal", "; Track length;", 30, 0, 700);
  hmap_trklen["cosmic"] = new TH1D("h_trklen_cosmic", "; Track length;", 30, 0, 700);
  hmap_trklen["cosmic_stopmu"] = new TH1D("h_trklen_cosmic_stopmu", "; Track length;", 30, 0, 700);
  hmap_trklen["cosmic_nostopmu"] = new TH1D("h_trklen_cosmic_nostopmu", "; Track length;", 30, 0, 700);
  hmap_trklen["outfv"] = new TH1D("h_trklen_outfv", "; Track length;", 30, 0, 700);
  hmap_trklen["outfv_stopmu"] = new TH1D("h_trklen_outfv_stopmu", "; Track length;", 30, 0, 700);
  hmap_trklen["outfv_nostopmu"] = new TH1D("h_trklen_outfv_nostopmu", "; Track length;", 30, 0, 700);
  hmap_trklen["nc"] = new TH1D("h_trklen_nc", "; Track length;", 30, 0, 700);
  hmap_trklen["nc_proton"] = new TH1D("h_trklen_nc_proton", "; Track length;", 30, 0, 700);
  hmap_trklen["nc_pion"] = new TH1D("h_trklen_nc_pion", "; Track length;", 30, 0, 700);
  hmap_trklen["nc_other"] = new TH1D("h_trklen_nc_other", "; Track length;", 30, 0, 700);
  hmap_trklen["anumu"] = new TH1D("h_trklen_anumu", "; Track length;", 30, 0, 700);
  hmap_trklen["nue"] = new TH1D("h_trklen_nue", "; Track length;", 30, 0, 700);
  hmap_trklen["signal_stopmu"] = new TH1D("h_trklen_signal_stopmu", "; Track length;", 30, 0, 700);
  hmap_trklen["signal_nostopmu"] = new TH1D("h_trklen_signal_nostopmu", "; Track length;", 30, 0, 700);
  hmap_trklen["cc_other"] = new TH1D("h_trklen_ccother", "; Track length;", 30, 0, 700);  
  hmap_trklen["cc_pion"] = new TH1D("h_trklen_ccpion", "; Track length;", 30, 0, 700);
  hmap_trklen["cc_0proton"] = new TH1D("h_trklen_cc0proton", "; Track length;", 30, 0, 700);

 
  std::map<std::string,TH1D*> hmap_trkmom_classic;
  hmap_trkmom_classic["total"] = new TH1D("h_trkmom_classic_total", "; Track momentum;", 25, 0, 2.5); // 20, 0, 2.5
  hmap_trkmom_classic["signal"] = new TH1D("h_trkmom_classic_signal", "; Track momentum;", 25, 0, 2.5);
  hmap_trkmom_classic["cosmic"] = new TH1D("h_trkmom_classic_cosmic", "; Track momentum;", 25, 0, 2.5);
  hmap_trkmom_classic["cosmic_stopmu"] = new TH1D("h_trkmom_classic_cosmic_stopmu", "; Track momentum;", 25, 0, 2.5);
  hmap_trkmom_classic["cosmic_nostopmu"] = new TH1D("h_trkmom_classic_cosmic_nostopmu", "; Track momentum;", 25, 0, 2.5);
  hmap_trkmom_classic["outfv"] = new TH1D("h_trkmom_classic_outfv", "; Track momentum;", 25, 0, 2.5);
  hmap_trkmom_classic["outfv_stopmu"] = new TH1D("h_trkmom_classic_outfv_stopmu", "; Track momentum;", 25, 0, 2.5);
  hmap_trkmom_classic["outfv_nostopmu"] = new TH1D("h_trkmom_classic_outfv_nostopmu", "; Track momentum;", 25, 0, 2.5);
  hmap_trkmom_classic["nc"] = new TH1D("h_trkmom_classic_nc", "; Track momentum;", 25, 0, 2.5);
  hmap_trkmom_classic["nc_proton"] = new TH1D("h_trkmom_classic_nc_proton", "; Track momentum;", 25, 0, 2.5);
  hmap_trkmom_classic["nc_pion"] = new TH1D("h_trkmom_classic_nc_pion", "; Track momentum;", 25, 0, 2.5);
  hmap_trkmom_classic["nc_other"] = new TH1D("h_trkmom_classic_nc_other", "; Track momentum;", 25, 0, 2.5);
  hmap_trkmom_classic["anumu"] = new TH1D("h_trkmom_classic_anumu", "; Track momentum;", 25, 0, 2.5);
  hmap_trkmom_classic["nue"] = new TH1D("h_trkmom_classic_nue", "; Track momentum;", 25, 0, 2.5);
  hmap_trkmom_classic["signal_stopmu"] = new TH1D("h_trkmom_classic_signal_stopmu", "; Track momentum;", 25, 0, 2.5);
  hmap_trkmom_classic["signal_nostopmu"] = new TH1D("h_trkmom_classic_signal_nostopmu", "; Track momentum;", 25, 0, 2.5);
  hmap_trkmom_classic["cc_other"] = new TH1D("h_trkmom_classic_ccother", "; Track momentum;", 25, 0, 2.5);
  hmap_trkmom_classic["cc_pion"] = new TH1D("h_trkmom_classic_ccpion", "; Track momentum;", 25, 0, 2.5);
  hmap_trkmom_classic["cc_0proton"] = new TH1D("h_trkmom_classic_cc0proton", "; Track momentum;", 25, 0, 2.5);
  // Number of events histograms - Cross Section Muon Momentum - GENIE pm1sigma

 

  std::map<std::string,TH1D*> hmap_trkphi;
  hmap_trkphi["total"] = new TH1D("h_trkphi_total", "; Track #phi;", 20, -3.15, 3.15);
  hmap_trkphi["signal"] = new TH1D("h_trkphi_signal", "; Track #phi;", 20, -3.15, 3.15);
  hmap_trkphi["cosmic"] = new TH1D("h_trkphi_cosmic", "; Track #phi;", 20, -3.15, 3.15);
  hmap_trkphi["outfv"] = new TH1D("h_trkphi_outfv", "; Track #phi;", 20, -3.15, 3.15);
  hmap_trkphi["nc"] = new TH1D("h_trkphi_nc", "; Track #phi;", 20, -3.15, 3.15);
  hmap_trkphi["anumu"] = new TH1D("h_trkphi_anumu", "; Track #phi;", 20, -3.15, 3.15);
  hmap_trkphi["nue"] = new TH1D("h_trkphi_nue", "; Track #phi;", 20, -3.15, 3.15);
  hmap_trkphi["cosmic_stopmu"] = new TH1D("h_trkphi_cosmic_stopmu", "; Track #phi;", 20, -3.15, 3.15);
  hmap_trkphi["cosmic_nostopmu"] = new TH1D("h_trkphi_cosmic_nostopmu", "; Track #phi;", 20, -3.15, 3.15);
  hmap_trkphi["outfv_stopmu"] = new TH1D("h_trkphi_outfv_stopmu", "; Track #phi;", 20, -3.15, 3.15);
  hmap_trkphi["outfv_nostopmu"] = new TH1D("h_trkphi_outfv_nostopmu", "; Track #phi;", 20, -3.15, 3.15);
  hmap_trkphi["nc_proton"] = new TH1D("h_trkphi_nc_proton", "; Track #phi;", 20, -3.15, 3.15);
  hmap_trkphi["nc_pion"] = new TH1D("h_trkphi_nc_pion", "; Track #phi;", 20, -3.15, 3.15);
  hmap_trkphi["nc_other"] = new TH1D("h_trkphi_nc_other", "; Track #phi;", 20, -3.15, 3.15);
  hmap_trkphi["signal_stopmu"] = new TH1D("h_trkphi_signal_stopmu", "; Track #phi;", 20, -3.15, 3.15);
  hmap_trkphi["signal_nostopmu"] = new TH1D("h_trkphi_signal_nostopmu", "; Track #phi;", 20, -3.15, 3.15);
  hmap_trkphi["cc_other"] = new TH1D("h_trkphi_ccother", "; Track #phi;", 20, -3.15,3.15);
  hmap_trkphi["cc_pion"] = new TH1D("h_trkphi_ccpion", "; Track #phi;", 20, -3.15, 3.15);
  hmap_trkphi["cc_0proton"] = new TH1D("h_trkphi_cc0proton", "; Track #phi;", 20, -3.15, 3.15);
  std::map<std::string,TH1D*> hmap_trktheta_classic;
  hmap_trktheta_classic["total"] = new TH1D("h_trktheta_classic_total", "; Track cos(#theta);", 30, -1, 1); // 30, -1, 1
  hmap_trktheta_classic["signal"] = new TH1D("h_trktheta_classic_signal", "; Track cos(#theta);", 30, -1, 1);
  hmap_trktheta_classic["cosmic"] = new TH1D("h_trktheta_classic_cosmic", "; Track cos(#theta);", 30, -1, 1);
  hmap_trktheta_classic["outfv"] = new TH1D("h_trktheta_classic_outfv", "; Track cos(#theta);", 30, -1, 1);
  hmap_trktheta_classic["nc"] = new TH1D("h_trktheta_classic_nc", "; Track cos(#theta);", 30, -1, 1);
  hmap_trktheta_classic["anumu"] = new TH1D("h_trktheta_classic_anumu", "; Track cos(#theta);", 30, -1, 1);
  hmap_trktheta_classic["nue"] = new TH1D("h_trktheta_classic_nue", "; Track cos(#theta);", 30, -1, 1);
  hmap_trktheta_classic["cosmic_stopmu"] = new TH1D("h_trktheta_classic_cosmic_stopmu", "; Track cos(#theta);", 30, -1, 1);
  hmap_trktheta_classic["cosmic_nostopmu"] = new TH1D("h_trktheta_classic_cosmic_nostopmu", "; Track cos(#theta);", 30, -1, 1);
  hmap_trktheta_classic["outfv_stopmu"] = new TH1D("h_trktheta_classic_outfv_stopmu", "; Track cos(#theta);", 30, -1, 1);
  hmap_trktheta_classic["outfv_nostopmu"] = new TH1D("h_trktheta_classic_outfv_nostopmu", "; Track cos(#theta);", 30, -1, 1);
  hmap_trktheta_classic["nc_proton"] = new TH1D("h_trktheta_classic_nc_proton", "; Track cos(#theta);", 30, -1, 1);
  hmap_trktheta_classic["nc_pion"] = new TH1D("h_trktheta_classic_nc_pion", "; Track cos(#theta);", 30, -1, 1);
  hmap_trktheta_classic["nc_other"] = new TH1D("h_trktheta_classic_nc_other", "; Track cos(#theta);", 30, -1, 1);
  hmap_trktheta_classic["signal_stopmu"] = new TH1D("h_trktheta_classic_signal_stopmu", "; Track cos(#theta);", 30, -1, 1);
  hmap_trktheta_classic["signal_nostopmu"] = new TH1D("h_trktheta_classic_signal_nostopmu", "; Track cos(#theta);", 30, -1, 1);
  hmap_trktheta_classic["cc_other"] = new TH1D("h_trktheta_classic_ccother", "; Track cos(#theta);", 30, -1, 1);
  hmap_trktheta_classic["cc_pion"] = new TH1D("h_trktheta_classic_ccpion", "; Track cos(#theta);", 30, -1, 1);
  hmap_trktheta_classic["cc_0proton"] = new TH1D("h_trktheta_classic_cc0proton", "; Track cos(#theta);", 30, -1, 1);
  std::map<std::string,TH1D*> hmap_multpfp;
  hmap_multpfp["total"] = new TH1D("h_multpfp_total", "; PFP Multiplicity", 10, 0, 10);
  hmap_multpfp["signal"] = new TH1D("h_multpfp_signal", "; PFP Multiplicity;", 10, 0, 10);
  hmap_multpfp["cosmic"] = new TH1D("h_multpfp_cosmic", "; PFP Multiplicity;", 10, 0, 10);
  hmap_multpfp["outfv"] = new TH1D("h_multpfp_outfv", "; PFP Multiplicity;", 10, 0, 10);
  hmap_multpfp["nc"] = new TH1D("h_multpfp_nc", "; PFP Multiplicity;", 10, 0, 10);
  hmap_multpfp["anumu"] = new TH1D("h_multpfp_anumu", "; PFP Multiplicity;", 10, 0, 10);
  hmap_multpfp["nue"] = new TH1D("h_multpfp_nue", "; PFP Multiplicity;", 10, 0, 10);
  hmap_multpfp["cosmic_stopmu"] = new TH1D("h_multpfp_cosmic_stopmu", "; PFP Multiplicity;", 10, 0, 10);
  hmap_multpfp["cosmic_nostopmu"] = new TH1D("h_multpfp_cosmic_nostopmu", "; PFP Multiplicity;", 10, 0, 10);
  hmap_multpfp["outfv_stopmu"] = new TH1D("h_multpfp_outfv_stopmu", "; PFP Multiplicity;", 10, 0, 10);
  hmap_multpfp["outfv_nostopmu"] = new TH1D("h_multpfp_outfv_nostopmu", "; PFP Multiplicity;", 10, 0, 10);
  hmap_multpfp["nc_proton"] = new TH1D("h_multpfp_nc_proton", "; PFP Multiplicity;", 10, 0, 10);
  hmap_multpfp["nc_pion"] = new TH1D("h_multpfp_nc_pion", "; PFP Multiplicity;", 10, 0, 10);
  hmap_multpfp["nc_other"] = new TH1D("h_multpfp_nc_other", "; PFP Multiplicity;", 10, 0, 10);
  hmap_multpfp["signal_stopmu"] = new TH1D("h_multpfp_signal_stopmu", "; PFP Multiplicity;", 10, 0, 10);
  hmap_multpfp["signal_nostopmu"] = new TH1D("h_multpfp_signal_nostopmu", "; PFP Multiplicity;", 10, 0, 10);
  hmap_multpfp["cc_other"] = new TH1D("h_multpfp_ccother", "; PFP Multiplicity;", 10, 0, 10);
  
  std::map<std::string,TH1D*> hmap_multtracktol;
  hmap_multtracktol["total"] = new TH1D("h_multtracktol_total", "; Track Multiplicity (5 cm)", 10, 0, 10);
  hmap_multtracktol["signal"] = new TH1D("h_multtracktol_signal", "; Track Multiplicity (5 cm);", 10, 0, 10);
  hmap_multtracktol["cosmic"] = new TH1D("h_multtracktol_cosmic", "; Track Multiplicity (5 cm);", 10, 0, 10);
  hmap_multtracktol["outfv"] = new TH1D("h_multtracktol_outfv", "; Track Multiplicity (5 cm);", 10, 0, 10);
  hmap_multtracktol["nc"] = new TH1D("h_multtracktol_nc", "; Track Multiplicity (5 cm);", 10, 0, 10);
  hmap_multtracktol["anumu"] = new TH1D("h_multtracktol_anumu", "; Track Multiplicity (5 cm);", 10, 0, 10);
  hmap_multtracktol["nue"] = new TH1D("h_multtracktol_nue", "; Track Multiplicity (5 cm);", 10, 0, 10);
  hmap_multtracktol["cosmic_stopmu"] = new TH1D("h_multtracktol_cosmic_stopmu", "; Track Multiplicity (5 cm);", 10, 0, 10);
  hmap_multtracktol["cosmic_nostopmu"] = new TH1D("h_multtracktol_cosmic_nostopmu", "; Track Multiplicity (5 cm);", 10, 0, 10);
  hmap_multtracktol["outfv_stopmu"] = new TH1D("h_multtracktol_outfv_stopmu", "; Track Multiplicity (5 cm);", 10, 0, 10);
  hmap_multtracktol["outfv_nostopmu"] = new TH1D("h_multtracktol_outfv_nostopmu", "; Track Multiplicity (5 cm);", 10, 0, 10);
  hmap_multtracktol["nc_proton"] = new TH1D("h_multtracktol_nc_proton", "; Track Multiplicity (5 cm);", 10, 0, 10);
  hmap_multtracktol["nc_pion"] = new TH1D("h_multtracktol_nc_pion", "; Track Multiplicity (5 cm);", 10, 0, 10);
  hmap_multtracktol["nc_other"] = new TH1D("h_multtracktol_nc_other", "; Track Multiplicity (5 cm);", 10, 0, 10);
  hmap_multtracktol["signal_stopmu"] = new TH1D("h_multtracktol_signal_stopmu", "; Track Multiplicity (5 cm);", 10, 0, 10);
  hmap_multtracktol["signal_nostopmu"] = new TH1D("h_multtracktol_signal_nostopmu", "; Track Multiplicity (5 cm);", 10, 0, 10);
  hmap_multtracktol["cc_other"] = new TH1D("h_multtracktol_ccother", "; Track Multiplicity (5 cm);", 10, 0, 10);
  std::map<std::string,TH1D*> hmap_dqdx_trunc;
  hmap_dqdx_trunc["total"] = new TH1D("h_dqdx_trunc_total", ";<dQ/dx>_{trunc};", 40, 0, 200000);//40, 0, 800);
  hmap_dqdx_trunc["muon"] = new TH1D("h_dqdx_trunc_muon", ";<dQ/dx>_{trunc};", 40, 0, 200000);
  hmap_dqdx_trunc["proton"] = new TH1D("h_dqdx_trunc_proton", ";<dQ/dx>_{trunc};", 40, 0, 200000);
  hmap_dqdx_trunc["pion"] = new TH1D("h_dqdx_trunc_pions", ";<dQ/dx>_{trunc};", 40, 0, 200000);
  hmap_dqdx_trunc["photon"] = new TH1D("h_dqdx_trunc_photon", ";<dQ/dx>_{trunc};", 40, 0, 200000);
  hmap_dqdx_trunc["electron"] = new TH1D("h_dqdx_trunc_electron", ";<dQ/dx>_{trunc};", 40, 0, 200000);
  hmap_dqdx_trunc["else"] = new TH1D("h_dqdx_trunc_else", ";<dQ/dx>_{trunc};", 40, 0, 200000);
  
  TH2D *h_dqdx_trunc_length = new TH2D("h_dqdx_trunc_length", ";Candidate Track <dQ/dx>_{trunc};Track Length [cm]", 40, 0, 200000,40, 0, 700);
  TH2D *h_dqdx_trunc_length_muon = new TH2D("h_dqdx_trunc_length_muon", ";Candidate Track <dQ/dx>_{trunc};Track Length [cm]", 40, 0, 200000,40, 0, 700);
  TH2D *h_dqdx_trunc_length_proton = new TH2D("h_dqdx_trunc_length_proton", ";Candidate Track <dQ/dx>_{trunc};Track Length [cm]", 40, 0, 200000,40, 0, 700);

  std::map<std::string,TH1D*> hmap_vtxx;
  hmap_vtxx["total"] = new TH1D("h_vtxx_total", ";Candidate Neutrino Vertex X [cm];", 40, 0, 275);
  hmap_vtxx["signal"] = new TH1D("h_vtxx_signal", ";Candidate Neutrino Vertex X [cm];", 40, 0,275);
  hmap_vtxx["background"] = new TH1D("h_vtxx_background", ";Candidate Neutrino Vertex X [cm];", 40, 0,275);
  
  std::map<std::string,TH1D*> hmap_vtxx_b;
  hmap_vtxx_b["total"] = new TH1D("h_vtxx_total_b", ";Candidate Neutrino Vertex X [cm];", 40, 0, 275);
  hmap_vtxx_b["signal"] = new TH1D("h_vtxx_signal_b", ";Candidate Neutrino Vertex X [cm];", 40, 0,275);
  hmap_vtxx_b["background"] = new TH1D("h_vtxx_background_b", ";Candidate Neutrino Vertex X [cm];", 40, 0,275);
  
  std::map<std::string,TH1D*> hmap_vtxy;
  hmap_vtxy["total"] = new TH1D("h_vtxy_total", ";Candidate Neutrino Vertex Y [cm];", 40, -125,125);
  hmap_vtxy["signal"] = new TH1D("h_vtxy_signal", ";Candidate Neutrino Vertex Y [cm];", 40, -125,125);
  hmap_vtxy["background"] = new TH1D("h_vtxy_background", ";Candidate Neutrino Vertex Y [cm];", 40, -125,125);
  
  std::map<std::string,TH1D*> hmap_vtxy_b;
  hmap_vtxy_b["total"] = new TH1D("h_vtxy_total_b", ";Candidate Neutrino Vertex Y [cm];", 40, -125,125);
  hmap_vtxy_b["signal"] = new TH1D("h_vtxy_signal_b", ";Candidate Neutrino Vertex Y [cm];", 40, -125,125);
  hmap_vtxy_b["background"] = new TH1D("h_vtxy_background_b", ";Candidate Neutrino Vertex Y [cm];", 40, -125,125);
  
  std::map<std::string,TH1D*> hmap_vtxz;
  hmap_vtxz["total"] = new TH1D("h_vtxz_total", ";Candidate Neutrino Vertex Z [cm];", 50, 0,1050);
  hmap_vtxz["signal"] = new TH1D("h_vtxz_signal", ";Candidate Neutrino Vertex Z [cm];", 50, 0,1050);
  hmap_vtxz["background"] = new TH1D("h_vtxz_background", ";Candidate Neutrino Vertex Z [cm];", 50, 0,1050);
  
  std::map<std::string,TH1D*> hmap_vtxz_b;
  hmap_vtxz_b["total"] = new TH1D("h_vtxz_total_b", ";Candidate Neutrino Vertex Z [cm];", 50, 0,1050);
  hmap_vtxz_b["signal"] = new TH1D("h_vtxz_signal_b", ";Candidate Neutrino Vertex Z [cm];", 50, 0,1050);
  hmap_vtxz_b["background"] = new TH1D("h_vtxz_background_b", ";Candidate Neutrino Vertex Z [cm];", 50, 0,1050);
  
  TH2D * h_vtx_xz = new TH2D("h_vtx_xz", ";X;Z", 40, 0, 275, 50, 0,1050);
  TH2D * h_vtx_xy = new TH2D("h_vtx_xy", ";X;Y", 40, 0, 275, 40, -125,125);

  std::map<std::string,TH1D*> hmap_vtxz_upborder;
  hmap_vtxz_upborder["total"] = new TH1D("h_vtxz_upborder_total", ";Candidate Neutrino Vertex Z [cm];", 3, 0,1050);
  hmap_vtxz_upborder["signal"] = new TH1D("h_vtxz_upborder_signal", ";Candidate Neutrino Vertex Z [cm];", 3, 0,1050);
  hmap_vtxz_upborder["background"] = new TH1D("h_vtxz_upborder_background", ";Candidate Neutrino Vertex Z [cm];", 3, 0,1050);

  std::map<std::string,TH1D*> hmap_vtxx_upborder;
  hmap_vtxx_upborder["total"] = new TH1D("h_vtxx_upborder_total", ";Candidate Neutrino Vertex X [cm];", 5, 0, 275);
  hmap_vtxx_upborder["signal"] = new TH1D("h_vtxx_upborder_signal", ";Candidate Neutrino Vertex X [cm];", 5, 0, 275);
  hmap_vtxx_upborder["background"] = new TH1D("h_vtxx_upborder_background", ";Candidate Neutrino Vertex X [cm];", 5, 0, 275);

  std::map<std::string,TH1D*> hmap_flsmatch_score;
  hmap_flsmatch_score["total"] = new TH1D("h_flsmatch_score_total", ";1/(-log(L));", 80, 0, 1.5);
  hmap_flsmatch_score["signal"] = new TH1D("h_flsmatch_score_signal", ";1/(-log(L));", 80, 0, 1.5);
  hmap_flsmatch_score["background"] = new TH1D("h_flsmatch_score_background", ";1/(-log(L));", 80, 0, 1.5);

  std::map<std::string,TH1D*> hmap_flsmatch_score_second;
  hmap_flsmatch_score_second["total"] = new TH1D("h_flsmatch_score_second_total", ";1/(-log(L));", 80, 0, 1.5);
  hmap_flsmatch_score_second["signal"] = new TH1D("h_flsmatch_score_second_signal", ";1/(-log(L));", 80, 0, 1.5);
  hmap_flsmatch_score_second["background"] = new TH1D("h_flsmatch_score_second_background", ";1/(-log(L));", 80, 0, 1.5);

  std::map<std::string,TH1D*> hmap_flsmatch_score_difference;
  hmap_flsmatch_score_difference["total"] = new TH1D("h_flsmatch_score_difference_total", ";1/(-log(L));", 80, 0, 0.2);
  hmap_flsmatch_score_difference["signal"] = new TH1D("h_flsmatch_score_difference_signal", ";1/(-log(L));", 80, 0, 0.2);
  hmap_flsmatch_score_difference["background"] = new TH1D("h_flsmatch_score_difference_background", ";1/(-log(L));", 80, 0, 0.2);

  std::map<std::string,TH1D*> hmap_ntpcobj;
  hmap_ntpcobj["total"] = new TH1D("h_ntpcobj_total", ";1/(-log(L));", 10, 0, 10);
  hmap_ntpcobj["signal"] = new TH1D("h_ntpcobj_signal", ";1/(-log(L));", 10, 0, 10);
  hmap_ntpcobj["background"] = new TH1D("h_ntpcobj_background", ";1/(-log(L));", 10, 0, 10);
  
  TH1D* h_pot = new TH1D("h_pot", "First bin contains number of POT (not valid on data)", 1, 0, 1);
  TH1D* h_nevts = new TH1D("h_nevts", "First bin contains number of events", 1, 0, 1);

  
  TH1D * h_deltall_cosmic_stop = new TH1D("h_deltall_cosmic_stop", "Cosmic stopping Muons;MCS Delta LL;", 400, -30, 30);
  TH2D * h_deltall_length_cosmic_stop = new TH2D("h_deltall_length_cosmic_stop", "Cosmic stopping Muons;MCS Delta LL;Track Length [cm]", 70, -30, 30, 70, 0, 700);
  TH1D * h_deltall_cosmic_nostop = new TH1D("h_deltall_cosmic_nostop", "Cosmic non-stopping Muons;MCS Delta LL;", 400, -30, 30);
  TH2D * h_deltall_length_cosmic_nostop = new TH2D("h_deltall_length_cosmic_nostop", "Cosmic non-stopping Muons;MCS Delta LL;Track Length [cm]", 70, -30, 30, 70, 0, 700);
  TH1D * h_deltall_nu = new TH1D("h_deltall_nu", "Neutrino origin;MCS Delta LL;", 400, -30, 30);
  TH2D * h_deltall_length_nu = new TH2D("h_deltall_length_nu", "Neutrino origin;MCS Delta LL;Track Length [cm]", 70, -30, 30, 70, 0, 700);

  TH1D* h_trklen_first = new TH1D("h_trklen_first", "h_trklen_first", 60, 0, 700);
  TH1D* h_trklen_second = new TH1D("h_trklen_second", "h_trklen_second", 60, 0, 700);
  //=========================================================================================================
  TH1D* chi2_proton_hypothesis_pcand=new TH1D("proton_hypothesis_pcand", "proton_hypothesis_pcand", 100, 0, 400);

  std::map<std::string, TH1D*> chi2_proton_hypothesis;
  chi2_proton_hypothesis["proton"] = new TH1D("proton_proton_chi2", "proton_proton_chi2", 100, 0, 400);
  chi2_proton_hypothesis["muon"] = new TH1D("proton_muon_chi2", "proton_muon_chi2", 100, 0, 400);
  chi2_proton_hypothesis["pion"] = new TH1D("proton_pion_chi2", "proton_pion_chi2", 100, 0, 400);
  chi2_proton_hypothesis["kaon"] = new TH1D("proton_kaon_chi2", "proton_kaon_chi2", 100, 0, 400);
  chi2_proton_hypothesis["other"] = new TH1D("proton_other_chi2", "proton_other_chi2", 100, 0, 400);
  std::map<std::string, TH1D*> chi2_muon_hypothesis;
  chi2_muon_hypothesis["proton"] = new TH1D("muon_proton_chi2", "muon_proton_chi2", 100, 0, 400);
  chi2_muon_hypothesis["muon"] = new TH1D("muon_muon_chi2", "muon_muon_chi2", 100, 0, 400);
  chi2_muon_hypothesis["pion"] = new TH1D("muon_pion_chi2", "muon_pion_chi2", 100, 0, 400);
  chi2_muon_hypothesis["kaon"] = new TH1D("muon_kaon_chi2", "muon_kaon_chi2", 100, 0, 400);
  chi2_muon_hypothesis["other"] = new TH1D("muon_other_chi2", "muon_other_chi2", 100, 0, 400);
  std::map<std::string, TH1D*> chi2_pion_hypothesis;
  chi2_pion_hypothesis["proton"] = new TH1D("pion_proton_chi2", "pion_proton_chi2", 100, 0, 400);
  chi2_pion_hypothesis["muon"] = new TH1D("pion_muon_chi2", "pion_muon_chi2", 100, 0, 400);
  chi2_pion_hypothesis["pion"] = new TH1D("pion_pion_chi2", "pion_pion_chi2", 100, 0, 400);
  chi2_pion_hypothesis["kaon"] = new TH1D("pion_kaon_chi2", "pion_kaon_chi2", 100, 0, 400);
  chi2_pion_hypothesis["other"] = new TH1D("pion_other_chi2", "pion_other_chi2", 100, 0, 400);
  std::map<std::string, TH1D*> chi2_kaon_hypothesis;
  chi2_kaon_hypothesis["proton"] = new TH1D("kaon_proton_chi2", "kaon_proton_chi2", 100, 0, 400);
  chi2_kaon_hypothesis["muon"] = new TH1D("kaon_muon_chi2", "kaon_muon_chi2", 100, 0, 400);
  chi2_kaon_hypothesis["pion"] = new TH1D("kaon_pion_chi2", "kaon_pion_chi2", 100, 0, 400);
  chi2_kaon_hypothesis["kaon"] = new TH1D("kaon_kaon_chi2", "kaon_kaon_chi2", 100, 0, 400);
  chi2_kaon_hypothesis["other"] = new TH1D("kaon_other_chi2", "kaon_other_chi2", 100, 0, 400);

  //proton resolution
  //======================================================================================
  //check how many percent of protons of all the tracks with number of hits less than 5
  TH1D* h_nhits_lt5=new TH1D("h_nhits_lt5", "h_nhits_lt5", 5, -0.5, 4.5);
  TH1D* h_nhits_lt5_proton=new TH1D("h_nhits_lt5_proton", "h_nhits_lt5_proton", 5, -0.5, 4.5);



  TH1D* h_pmom_reint=new TH1D("h_pmom_reint",";Proton Momentum [GeV/c]; Nevts", 30, 0.3, 1.5);
  TH1D* h_pmom_total=new TH1D("h_pmom_total",";Proton MOmentum [GeV/c]; Nevts", 30, 0.3, 1.5);

  TH1D* h_ngenie_proton=new TH1D("h_ngenie_proton", ";Number of Protons < 300 MeV/c;", 20, -0.5, 19.5);
  // histograms for broken tracks study  
  TH1D* h_true_photonmom_cc0p= new TH1D("h_true_photonmom_cc0p", ";Photon's Momentum(True)[GeV/c];Ntracks", 30, 0.0, 1.0);
  TH1D* h_reco_photonmom_cc0p= new TH1D("h_reco_photonmom_cc0p", ";Photon's Momentum(True)[GeV/c];Ntracks", 30, 0.0, 1.0);
  TH1D* h_true_thetamup_cc0p=new TH1D("h_true_thetamup_cc0p", "h_true_thetamup_cc0p", 60, 0.0, 3.14);
  TH1D* h_reco_thetamup_cc0p=new TH1D("h_reco_thetamup_cc0p", "h_reco_thetamup_cc0p", 60, 0.0, 3.14);

  TH2D* h_true_reco_protonmom_cc0p= new TH2D("h_true_reco_protonmom_cc0p", ";Proton's Momentum(True)[GeV/c];Proton's Momentum(Reco) [GeV]", 30, 0.0, 1.0, 30, 0.0, 1.0);
  TH2D* h_true_reco_protonlen_cc0p= new TH2D("h_true_reco_protonlen_cc0p", ";Proton's Track Length(True)[cm];Proton's Track Length(Reco) [cm]", 40, 0.0, 10.0, 40, 0.0, 10.0);
 
  //***************************************************************************************
  TH1D* h_true_pmom_cc0p=new TH1D("h_true_pmom_cc0p", "h_true_pmom_cc0p", 80, 0, 0.8);
  //check if the pions reinteracted or reconstructed
  TH1D* h_pion_reco= new TH1D("h_pion_reco", ";pion_reco ; Nevts", 2, -0.5, 1.5);
  TH1D* h_pion_reint=new TH1D("h_pion_reint",";pion_reint; Nevts", 2, -0.5, 1.5);

 
  TH1D* h_pimom_reint=new TH1D("h_pimom_reint", ";Pion Momentum [GeV/c]; Ntracks", 30, 0.0,1.2);
  TH1D* h_pimom_total=new TH1D("h_pimom_total", ";Pion Momentum [GeV/c]; Ntracks", 30, 0.0,1.2);

  TH1D* h_true_pimom_noreco=new TH1D("h_true_pimom_noreco", ";Pion Momentum [GeV/c]; Ntracks", 30, 0.0, 1.2);
  TH1D* h_true_pilen_noreco=new TH1D("h_true_pilen_noreco", ";Pion Length [cm]; Ntracks", 40, 0.0, 20.0);
  TH1D* h_true_pizlen_noreco=new TH1D("h_true_pizlen_noreco", ";Pion Length [cm]; Ntracks", 40, 0.0, 20.0);


  TH2D* h_dEdx_vs_rr_pcand=new TH2D("h_dEdx_vs_rr_pcand", "h_dEdx_vs_rr_pcand",100,0,100,100, 0, 20);
  h_dEdx_vs_rr_pcand->GetXaxis()->SetTitle("Residual Range [cm]");
  h_dEdx_vs_rr_pcand->GetYaxis()->SetTitle("dEdx [MeV/cm]"); 
  TH2D* h_dEdx_vs_rr_pcand_2=new TH2D("h_dEdx_vs_rr_pcand_2", "h_dEdx_vs_rr_pcand",100,0,100,100, 0, 20);
  h_dEdx_vs_rr_pcand_2->GetXaxis()->SetTitle("Residual Range [cm]");
  h_dEdx_vs_rr_pcand_2->GetYaxis()->SetTitle("dEdx [MeV/cm]"); 
 
  TH2D* h_dEdx_vs_rr_pshower=new TH2D("h_dEdx_vs_rr_pshower", "h_dEdx_vs_rr_pshower", 100, 0, 100, 100, 0, 20);
  h_dEdx_vs_rr_pshower->GetXaxis()->SetTitle("Residual Range [cm]");
  h_dEdx_vs_rr_pshower->GetYaxis()->SetTitle("dEdx [MeV/cm]"); 
 
   

 
  TH2D* h_tmdqdx_vs_rr_pcand=new TH2D("h_tmdqdx_vs_rr_pcand", "h_tmdqdx_vs_rr_pcand", 100, 0, 100, 100, 0, 2000);
  
  //==========================================================================================================
  std::map<std::string,TH1D*> hmap_trkplen;
  hmap_trkplen["total"] = new TH1D("h_trkplen_total", "; Track length;", 30, 0, 150); // 20, 0, 2.5
  hmap_trkplen["signal"] = new TH1D("h_trkplen_signal", "; Track length;", 30, 0, 150);
  hmap_trkplen["cosmic"] = new TH1D("h_trkplen_cosmic", "; Track length;", 30, 0, 150);
  hmap_trkplen["outfv"] = new TH1D("h_trkplen_outfv", "; Track length;", 30, 0, 150);
  hmap_trkplen["nc"] = new TH1D("h_trkplen_nc", "; Track length;", 30, 0, 150);
  hmap_trkplen["anumu"] = new TH1D("h_trkplen_anumu", "; Track length;", 30, 0, 150);
  hmap_trkplen["nue"] = new TH1D("h_trkplen_nue", "; Track length;", 30, 0, 150);
  hmap_trkplen["cc_other"] = new TH1D("h_trkplen_ccother", "; Track length;", 30, 0, 150);
  hmap_trkplen["cc_pion"] = new TH1D("h_trkplen_ccpion", "; Track length;", 30, 0, 150);
  hmap_trkplen["cc_0proton"] = new TH1D("h_trkplen_cc0proton", "; Track length;", 30, 0, 150);
  std::map<std::string,TH1D*> hmap_trkpmom_classic;
  hmap_trkpmom_classic["total"] = new TH1D("h_trkpmom_classic_total", "; Track momentum;", 150, 0, 1.5); // 20, 0, 2.5
  hmap_trkpmom_classic["signal"] = new TH1D("h_trkpmom_classic_signal", "; Track momentum;", 150, 0, 1.5);
  hmap_trkpmom_classic["cosmic"] = new TH1D("h_trkpmom_classic_cosmic", "; Track momentum;", 150, 0, 1.5);
  hmap_trkpmom_classic["outfv"] = new TH1D("h_trkpmom_classic_outfv", "; Track momentum;", 150, 0, 1.5);
  hmap_trkpmom_classic["nc"] = new TH1D("h_trkpmom_classic_nc", "; Track momentum;", 150, 0, 1.5);
  hmap_trkpmom_classic["anumu"] = new TH1D("h_trkpmom_classic_anumu", "; Track momentum;", 150, 0, 1.5);
  hmap_trkpmom_classic["nue"] = new TH1D("h_trkpmom_classic_nue", "; Track momentum;", 150, 0, 1.5);
  hmap_trkpmom_classic["cc_other"] = new TH1D("h_trkpmom_classic_ccother", "; Track momentum;", 150, 0, 1.5);
  hmap_trkpmom_classic["cc_pion"] = new TH1D("h_trkpmom_classic_ccpion", "; Track momentum;", 150, 0, 1.5);
  hmap_trkpmom_classic["cc_0proton"]= new TH1D("h_trkpmom_classic_cc0proton", "; Track momentum;", 150, 0, 1.5);
  std::map<std::string,TH1D*> hmap_trkptheta_classic;
  hmap_trkptheta_classic["total"] = new TH1D("h_trkptheta_classic_total", "; Track cos(#theta);", 30, -1, 1); // 30, -1, 1
  hmap_trkptheta_classic["signal"] = new TH1D("h_trkptheta_classic_signal", "; Track cos(#theta);", 30, -1, 1);
  hmap_trkptheta_classic["cosmic"] = new TH1D("h_trkptheta_classic_cosmic", "; Track cos(#theta);", 30, -1, 1);
  hmap_trkptheta_classic["outfv"] = new TH1D("h_trkptheta_classic_outfv", "; Track cos(#theta);", 30, -1, 1);
  hmap_trkptheta_classic["nc"] = new TH1D("h_trkptheta_classic_nc", "; Track cos(#theta);", 30, -1, 1);
  hmap_trkptheta_classic["anumu"] = new TH1D("h_trkptheta_classic_anumu", "; Track cos(#theta);", 30, -1, 1);
  hmap_trkptheta_classic["nue"] = new TH1D("h_trkptheta_classic_nue", "; Track cos(#theta);", 30, -1, 1);
  hmap_trkptheta_classic["cc_other"] = new TH1D("h_trkptheta_classic_ccother", "; Track cos(#theta);", 30, -1, 1);
  hmap_trkptheta_classic["cc_pion"] = new TH1D("h_trkptheta_classic_ccpion", "; Track cos(#theta);", 30, -1, 1);
  hmap_trkptheta_classic["cc_0proton"] = new TH1D("h_trkptheta_classic_cc0proton", "; Track cos(#theta);", 30, -1, 1);
  std::map<std::string,TH1D*> hmap_trkpphi;
  hmap_trkpphi["total"] = new TH1D("h_trkpphi_total", "; Track #phi;", 20, -3.15, 3.15);
  hmap_trkpphi["signal"] = new TH1D("h_trkpphi_signal", "; Track #phi;", 20, -3.15, 3.15);
  hmap_trkpphi["cosmic"] = new TH1D("h_trkpphi_cosmic", "; Track #phi;", 20, -3.15, 3.15);
  hmap_trkpphi["outfv"] = new TH1D("h_trkpphi_outfv", "; Track #phi;", 20, -3.15, 3.15);
  hmap_trkpphi["nc"] = new TH1D("h_trkpphi_nc", "; Track #phi;", 20, -3.15, 3.15);
  hmap_trkpphi["anumu"] = new TH1D("h_trkpphi_anumu", "; Track #phi;", 20, -3.15, 3.15);
  hmap_trkpphi["nue"] = new TH1D("h_trkpphi_nue", "; Track #phi;", 20, -3.15, 3.15);
  hmap_trkpphi["cc_other"] = new TH1D("h_trkpphi_ccother", "; Track #phi;", 20, -3.15,3.15);
  hmap_trkpphi["cc_pion"] = new TH1D("h_trkpphi_ccpion", "; Track #phi;", 20,  -3.15, 3.15);
  hmap_trkpphi["cc_0proton"] = new TH1D("h_trkpphi_cc0proton", "; Track #phi;", 20, -3.15, 3.15);
  std::map<std::string,TH1D*> hmap_thetamup;
  hmap_thetamup["total"] = new TH1D("h_thetamup_total", "; #theta_{#mu P};", 30, 0, 3.15);
  hmap_thetamup["signal"] = new TH1D("h_thetamup_signal", "; #theta_{#mu P};", 30, 0, 3.15);
  hmap_thetamup["cosmic"] = new TH1D("h_thetamup_cosmic", "; #theta_{#mu P};", 30, 0, 3.15);
  hmap_thetamup["outfv"] = new TH1D("h_thetamup_outfv", "; #theta_{#mu P};", 30, 0, 3.15);
  hmap_thetamup["nc"] = new TH1D("h_thetamup_nc", "; #theta_{#mu P};", 30, 0, 3.15);
  hmap_thetamup["anumu"] = new TH1D("h_thetamup_anumu", "; #theta_{#mu P};", 30, 0, 3.15);
  hmap_thetamup["nue"] = new TH1D("h_thetamup_nue", "; #theta_{#mu P};", 30, 0, 3.15);
  hmap_thetamup["cc_other"] = new TH1D("h_thetamup_ccother", "; #theta_{#mu P};", 30, 0,3.15);
  std::map<std::string,TH1D*> hmap_ptmis;
  hmap_ptmis["total"] = new TH1D("h_ptmis_total", "; Missing-Pt;", 30, 0, 1.0);
  hmap_ptmis["signal"] = new TH1D("h_ptmis_signal", "; Missing-Pt;", 30, 0, 1.0);
  hmap_ptmis["cosmic"] = new TH1D("h_ptmis_cosmic", "; Missing-Pt;", 30, 0, 1.0);
  hmap_ptmis["outfv"] = new TH1D("h_ptmis_outfv", "; Missing-Pt;", 30, 0, 1.0);
  hmap_ptmis["nc"] = new TH1D("h_ptmis_nc", "; Missing-Pt;", 30, 0, 1.0);
  hmap_ptmis["anumu"] = new TH1D("h_ptmis_anumu", "; Missing-Pt;", 30, 0, 1.0);
  hmap_ptmis["nue"] = new TH1D("h_ptmis_nue", "; Missing-Pt;", 30, 0, 1.0);
  hmap_ptmis["cc_other"] = new TH1D("h_ptmis_ccother", "; Missing-Pt;", 30, 0,1.0);
  std::map<std::string,TH1D*> hmap_etatest;
  hmap_etatest["total"] = new TH1D("h_etatest_total", "; #eta;", 30, -1.0, 1.0);
  hmap_etatest["signal"] = new TH1D("h_etatest_signal", "; #eta;", 30, -1.0, 1.0);
  hmap_etatest["cosmic"] = new TH1D("h_etatest_cosmic", "; #eta;", 30, -1.0, 1.0);
  hmap_etatest["outfv"] = new TH1D("h_etatest_outfv", "; #eta;", 30, -1.0, 1.0);
  hmap_etatest["nc"] = new TH1D("h_etatest_nc", "; #eta;", 30, -1.0, 1.0);
  hmap_etatest["anumu"] = new TH1D("h_etatest_anumu", "; #eta;", 30, -1.0, 1.0);
  hmap_etatest["nue"] = new TH1D("h_etatest_nue", "; #eta;", 30, -1.0, 1.0);
  hmap_etatest["cc_other"] = new TH1D("h_etatest_ccother", "; #eta;", 30, -1.0,1.0);
  
  std::map<std::string,TH1D*> hmap_alphat;
  hmap_alphat["total"] = new TH1D("h_alphat_total", "; #alpha_{T};", 30, 0, 3.14);
  hmap_alphat["signal"] = new TH1D("h_alphat_signal", "; #alpha_{T};", 30, 0, 3.14);
  hmap_alphat["cosmic"] = new TH1D("h_alphat_cosmic", "; #alpha_{T};", 30, 0, 3.14);
  hmap_alphat["outfv"] = new TH1D("h_alphat_outfv", "; #alpha_{T};", 30, 0, 3.14);
  hmap_alphat["nc"] = new TH1D("h_alphat_nc", "; #alpha_{T};", 30, 0, 3.14);
  hmap_alphat["anumu"] = new TH1D("h_alphat_anumu", "; #alpha_{T};", 30, 0, 3.14);
  hmap_alphat["nue"] = new TH1D("h_alphat_nue", "; #alpha_{T};", 30, 0, 3.14);
  hmap_alphat["cc_other"] = new TH1D("h_alphat_ccother", "; #alpha_{T};", 30, 0,3.14);
  std::map<std::string,TH1D*> hmap_phit;
  hmap_phit["total"] = new TH1D("h_phit_total", "; #alpha_{T};", 30, 0, 3.14);
  hmap_phit["signal"] = new TH1D("h_phit_signal", "; #alpha_{T};", 30, 0, 3.14);
  hmap_phit["cosmic"] = new TH1D("h_phit_cosmic", "; #alpha_{T};", 30, 0, 3.14);
  hmap_phit["outfv"] = new TH1D("h_phit_outfv", "; #alpha_{T};", 30, 0, 3.14);
  hmap_phit["nc"] = new TH1D("h_phit_nc", "; #alpha_{T};", 30, 0, 3.14);
  hmap_phit["anumu"] = new TH1D("h_phit_anumu", "; #alpha_{T};", 30, 0, 3.14);
  hmap_phit["nue"] = new TH1D("h_phit_nue", "; #alpha_{T};", 30, 0, 3.14);
  hmap_phit["cc_other"] = new TH1D("h_phit_ccother", "; #alpha_{T};", 30, 0,3.14);
  std::map<std::string,TH1D*> hmap_enucal;
  hmap_enucal["total"] = new TH1D("h_enucal_total", "; Calculated Enu;", 30, 0, 2000.0);
  hmap_enucal["signal"] = new TH1D("h_enucal_signal", "; Calculated Enu;", 30, 0, 2000.0);
  hmap_enucal["cosmic"] = new TH1D("h_enucal_cosmic", "; Calculated Enu;", 30, 0, 2000.0);
  hmap_enucal["outfv"] = new TH1D("h_enucal_outfv", "; Calculated Enu;", 30, 0, 2000.0);
  hmap_enucal["nc"] = new TH1D("h_enucal_nc", "; Calculated Enu;", 30, 0, 2000.0);
  hmap_enucal["anumu"] = new TH1D("h_enucal_anumu", "; Calculated Enu;", 30, 0, 2000.0);
  hmap_enucal["nue"] = new TH1D("h_enucal_nue", "; Calculated Enu;", 30, 0, 2000.0);
  hmap_enucal["cc_other"] = new TH1D("h_enucal_ccother", "; Calculated Enu;", 30, 0,2000.0);
  std::map<std::string,TH1D*> hmap_pmult;
  hmap_pmult["total"] = new TH1D("h_pmult_total", "; Proton Multiplicity;", 10, -0.5, 9.5);
  hmap_pmult["signal"] = new TH1D("h_pmult_signal", "; Proton Multiplicity;", 10, -0.5, 9.5);
  hmap_pmult["cosmic"] = new TH1D("h_pmult_cosmic", "; Proton Multiplicity;", 10, -0.5, 9.5);
  hmap_pmult["outfv"] = new TH1D("h_pmult_outfv", "; Proton Multiplicity;", 10, -0.5, 9.5);
  hmap_pmult["nc"] = new TH1D("h_pmult_nc", "; Proton Multiplicity;", 10, -0.5, 9.5);
  hmap_pmult["anumu"] = new TH1D("h_pmult_anumu", "; Proton Multiplicity;", 10, -0.5, 9.5);
  hmap_pmult["nue"] = new TH1D("h_pmult_nue", "; Proton Multiplicity;", 10, -0.5, 9.5);
  hmap_pmult["cc_other"] = new TH1D("h_pmult_ccother", "; Proton Multiplicity;", 10, -0.5,9.5);

  std::map<std::string,TH1D*> hmap_nhits_leadingp;
  hmap_nhits_leadingp["total"] = new TH1D("h_nhits_leadingp_total", "; nhits_leadingp;", 30, 0, 300.0);
  hmap_nhits_leadingp["signal"] = new TH1D("h_nhits_leadingp_signal", "; nhits_leadingp;", 30, 0, 300.0);
  hmap_nhits_leadingp["cosmic"] = new TH1D("h_nhits_leadingp_cosmic", "; nhits_leadingp;", 30, 0, 300.0);
  hmap_nhits_leadingp["outfv"] = new TH1D("h_nhits_leadingp_outfv", "; nhits_leadingp;", 30, 0, 300.0);
  hmap_nhits_leadingp["nc"] = new TH1D("h_nhits_leadingp_nc", "; nhits_leadingp;", 30, 0, 300.0);
  hmap_nhits_leadingp["anumu"] = new TH1D("h_nhits_leadingp_anumu", "; nhits_leadingp;", 30, 0, 300.0);
  hmap_nhits_leadingp["nue"] = new TH1D("h_nhits_leadingp_nue", "; nhits_leadingp;", 30, 0, 300.0);
  hmap_nhits_leadingp["cc_other"] = new TH1D("h_nhits_leadingp_ccother", "; nhits_leadingp;", 30, 0,300.0);
 


   //===========================================================================================================


  std::vector<std::string> fname_genie_pm1;
  std::vector<std::string> fname_genie_multisim;
  std::vector<std::string> fname_mc_stat_multisim;
  std::vector<std::string> fname_extra_syst;
  std::vector<std::string> fname_flux_multisim;


  if (_maup_mecoff && !isdata) {
    PrintMaUpMECOff();
  }

  if (_reweigh_kaons) {
    PrintReweighKaons();
  }
    
  int barWidth = 70;
  
  if(maxEntries > 0.) evts = maxEntries;

  evts += _initial_entry;

  LOG_NORMAL() << "Looping over " << evts - _initial_entry << " events starting with entry " << _initial_entry << std::endl;
  
  int total_events = 0;

  std::vector<int> run_numbers, subrun_numbers, event_numbers;
  run_numbers.resize(evts); subrun_numbers.resize(evts); event_numbers.resize(evts);
  
  //ofstream outfile;
  //outfile.open("rewght_genie_rwfac.txt");
   
  for(int i = _initial_entry; i < evts; i++) {
    
    if (i != 0) DrawProgressBar((double)i/(double)evts, barWidth);
    
    chain_ubxsec->GetEntry(i);
    
    total_events ++;
    
    //cout << "***** Event " << i << endl;
    //cout << "***** Event Number " << t->event << endl;
    run_numbers.at(i) = t->run;
    subrun_numbers.at(i) = t->subrun;
    event_numbers.at(i) = t->event;

    // Check for duplicate MC events
    if (_check_duplicate_events){
      
      if (std::count (event_numbers.begin(), event_numbers.end(), t->event) > 1) {

        // Now check the subrun
        for (size_t i_ev = 0; i_ev < event_numbers.size(); i_ev++) {
          if (event_numbers.at(i_ev) == t->event) {

            if (run_numbers.at(i_ev) == t->run && subrun_numbers.at(i_ev) == t->subrun) {
              std::cout << "Found duplicate event: " << t->event << std::endl;
            }
            break;
          }
        }
      }
    }


    // ************************
    //
    // Total event weight (BNB Correction)
    //
    // ************************

    double event_weight = t->bnb_weight;
    event_weight *= _extra_weight;
    if (isdata) event_weight = 1.;

    if (t->file_type == "dirt") event_weight /= _extra_weight;


    bool is_from_kaon = false;

    // ************************
    //
    // Check if running with Ma+1sigma and MEC off
    //
    // ************************

    if(_maup_mecoff && !isdata) {

      // Remove MEC events
      if (t->mode == 10) {
        continue;
      }

      // Scale up Ma CCQE
      for (size_t i = 0; i < t->evtwgt_genie_pm1_weight.size(); i++) {
        if (t->evtwgt_genie_pm1_funcname.at(i) == "genie_qema_Genie") {
          event_weight *= t->evtwgt_genie_pm1_weight.at(i).at(0);
        }
      }
    }

    if (!isdata && false) {
      LOG_CRITICAL() << "SPECIAL WEIGHTS APPLIED!!! MODEL 0" << std::endl;
      if (t->mode == 0) { // QE
        event_weight *= 0.95; 
      }
      if (t->mode == 1) { // RES
        event_weight *= 0.75; 
      }
      if (t->mode == 2) { // DIS
        event_weight *= 0.85; 
      }
      if (t->mode == 3) { // COH
        event_weight *= 1.00; 
      }
      if (t->mode == 10) { // MEC
        event_weight *= 0.85; 
      }
    }

    if (!isdata && false) {
      LOG_CRITICAL() << "SPECIAL WEIGHTS APPLIED!!! MODEL 1" << std::endl;
      if (t->mode == 0) { // QE
        event_weight *= 0.90; 
      }
      if (t->mode == 1) { // RES
        event_weight *= 0.00; 
      }
      if (t->mode == 2) { // DIS
        event_weight *= 3.00; 
      }
      if (t->mode == 3) { // COH
        event_weight *= 1.00; 
      }
      if (t->mode == 10) { // MEC
        event_weight *= 1.10; 
      }
    }

    if (!isdata && false) {
      LOG_CRITICAL() << "SPECIAL WEIGHTS APPLIED!!! MODEL 2" << std::endl;
      if (t->mode == 0) { // QE
        event_weight *= 1.00; 
      }
      if (t->mode == 1) { // RES
        event_weight *= 1.50; 
      }
      if (t->mode == 2) { // DIS
        event_weight *= 1.00; 
      }
      if (t->mode == 3) { // COH
        event_weight *= 1.00; 
      }
      if (t->mode == 10) { // MEC
        event_weight *= 0.00; 
      }
    }




    // ************************
    //
    // Set weight names, prepare bootstraps -- PM1SIGMA
    //
    // ************************

    // Set the weight names, just do it once (first event only)
    //
    if (i == _initial_entry && !isdata && _fill_bootstrap_genie) {
      ofstream myfile;
      myfile.open("rewght_genie_pm1.txt");
      for (auto name : t->evtwgt_genie_pm1_funcname) {
        fname_genie_pm1.push_back(name + "_p1");
        fname_genie_pm1.push_back(name + "_m1");

        //open a txt file and save the reweight factors
        myfile<<std::setw(30)<<std::setprecision(20)<<name + "_p1 "<<std::setw(30)<<std::setprecision(20)<<name + "_m1"<<std::endl;  
        
        
        //==============================================================
      }
      myfile.close();
      //std::cout<<"libo test 1"<<std::endl;
      // Number of events
      for (auto iter : _event_histo_1d->hmap_trkmom_genie_pm1_bs) {

        std::string this_name = iter.first;
        std::map<std::string, TH1D*> bs_map = iter.second;

        // Now emplace the histograms for the variations
        for (size_t i = 0; i < fname_genie_pm1.size(); i++) {

          std::string histo_name = "h_trkmom_" + this_name + "_" + fname_genie_pm1.at(i);
          double this_bins_mumom[7] = {0.1, 0.18, 0.30, 0.48, 0.75, 1.14, 2.50};
          _event_histo_1d->hmap_trkmom_genie_pm1_bs[this_name][fname_genie_pm1.at(i)] = new TH1D(histo_name.c_str(), "; Track momentum;", 6, this_bins_mumom); 

        }

      }
      //std::cout<<"libo test 2"<<std::endl;
      for (auto iter : _event_histo_1d->hmap_trkpmom_genie_pm1_bs) {

        std::string this_name = iter.first;
        std::map<std::string, TH1D*> bs_map = iter.second;

        // Now emplace the histograms for the variations
        for (size_t i = 0; i < fname_genie_pm1.size(); i++) {

          std::string histo_name = "h_trkpmom_" + this_name + "_" + fname_genie_pm1.at(i);
          double this_bins_pmom[11] = {0.30, 0.41, 0.495, 0.56, 0.62, 0.68, 0.74, 0.80, 0.87, 0.93, 1.20};
          _event_histo_1d->hmap_trkpmom_genie_pm1_bs[this_name][fname_genie_pm1.at(i)] = new TH1D(histo_name.c_str(), "; Track momentum;", 10, this_bins_pmom); 

        }

      }
      //std::cout<<"libo test 3"<<std::endl;
      for (auto iter : _event_histo_1d->hmap_trktheta_genie_pm1_bs) {

        std::string this_name = iter.first;
        std::map<std::string, TH1D*> bs_map = iter.second;

        // Now emplace the histograms for the variations
        for (size_t i = 0; i < fname_genie_pm1.size(); i++) {

          std::string histo_name = "h_trktheta_" + this_name + "_" + fname_genie_pm1.at(i);
          double this_bins_mucostheta[13] = {-1.00, -0.82, -0.66, -0.39, -0.16, 0.05, 0.25, 0.43, 0.59, 0.73, 0.83, 0.91, 1.00};
          _event_histo_1d->hmap_trktheta_genie_pm1_bs[this_name][fname_genie_pm1.at(i)] = new TH1D(histo_name.c_str(), "; Track Angle;", 12, this_bins_mucostheta); 

        }

      }
      //std::cout<<"libo test 4"<<std::endl;
      for (auto iter : _event_histo_1d->hmap_trkptheta_genie_pm1_bs) {

        std::string this_name = iter.first;
        std::map<std::string, TH1D*> bs_map = iter.second;

        // Now emplace the histograms for the variations
        for (size_t i = 0; i < fname_genie_pm1.size(); i++) {

          std::string histo_name = "h_trkptheta_" + this_name + "_" + fname_genie_pm1.at(i);
          double this_bins_pcostheta[10] = {-1.00, -0.50, 0.00, 0.27, 0.45, 0.62, 0.76, 0.86, 0.94, 1.00};
          _event_histo_1d->hmap_trkptheta_genie_pm1_bs[this_name][fname_genie_pm1.at(i)] = new TH1D(histo_name.c_str(), "; Track Angle;", 9, this_bins_pcostheta); 
      
        }

      }
      //std::cout<<"libo test 5"<<std::endl;
      for (auto iter : _event_histo_1d->hmap_thetamup_genie_pm1_bs) {

        std::string this_name = iter.first;
        std::map<std::string, TH1D*> bs_map = iter.second;

        // Now emplace the histograms for the variations
        for (size_t i = 0; i < fname_genie_pm1.size(); i++) {

          std::string histo_name = "h_thetamup_" + this_name + "_" + fname_genie_pm1.at(i);
          double this_bins_muptheta[7] = {0.00, 0.8, 1.2, 1.57, 1.94, 2.34, 3.14};
          _event_histo_1d->hmap_thetamup_genie_pm1_bs[this_name][fname_genie_pm1.at(i)] = new TH1D(histo_name.c_str(), "; Track Angle (MuP);", 6, this_bins_muptheta); 

        }

      }
       // Efficiency
      /*for (size_t i = 0; i < fname_genie_pm1.size(); i++) {
        double this_bins_mumom[7] = {0.1, 0.18, 0.30, 0.48, 0.75, 1.14, 2.50};
        double this_bins_mucostheta[13] = {-1.00, -0.82, -0.66, -0.39, -0.16, 0.05, 0.25, 0.43, 0.59, 0.73, 0.83, 0.91, 1.00};
        double this_bins_pmom[11] = {0.30, 0.41, 0.495, 0.56, 0.62, 0.68, 0.74, 0.80, 0.87, 0.93, 1.20};
        double this_bins_pcostheta[10] = {-1.00, -0.50, 0.00, 0.27, 0.45, 0.62, 0.76, 0.86, 0.94, 1.00};
        double this_bins_muptheta[7] = {0.00, 0.8, 1.2, 1.57, 1.94, 2.34, 3.14};
 
  
        std::string histo_name;// = "bs_genie_pm1_eff_mumom_num_" + fname_genie_pm1.at(i);

        histo_name = "bs_genie_pm1_true_reco_mom_" + fname_genie_pm1.at(i);
        _event_histo_1d->bs_genie_pm1_true_reco_mom[fname_genie_pm1.at(i)] = new TH2D(histo_name.c_str(), ";Muon Momentum (Truth) [GeV]; Muon Momentum (MCS) [GeV]", 6, this_bins_mumom, 6, this_bins_mumom);

        histo_name = "bs_genie_pm1_true_reco_pmom_" + fname_genie_pm1.at(i);
        _event_histo_1d->bs_genie_pm1_true_reco_pmom[fname_genie_pm1.at(i)] = new TH2D(histo_name.c_str(), ";Proton Momentum (Truth) [GeV]; Proton Momentum (MCS) [GeV]", 10, this_bins_pmom, 10, this_bins_pmom);

        histo_name = "bs_genie_pm1_true_reco_muangle_" + fname_genie_pm1.at(i);
        _event_histo_1d->bs_genie_pm1_true_reco_muangle[fname_genie_pm1.at(i)] = new TH2D(histo_name.c_str(), ";Muon Angle (Truth); Muon Angle", 12, this_bins_mucostheta, 12, this_bins_mucostheta);

        histo_name = "bs_genie_pm1_true_reco_pangle_" + fname_genie_pm1.at(i);
        _event_histo_1d->bs_genie_pm1_true_reco_pangle[fname_genie_pm1.at(i)] = new TH2D(histo_name.c_str(), ";Proton Angle (Truth); Proton Angle", 9, this_bins_pcostheta, 9, this_bins_pcostheta);

        histo_name = "bs_genie_pm1_true_reco_thetamup_" + fname_genie_pm1.at(i);
        _event_histo_1d->bs_genie_pm1_true_reco_thetamup[fname_genie_pm1.at(i)] = new TH2D(histo_name.c_str(), ";Track Angle (Truth); Track Angle", 6, this_bins_muptheta, 6, this_bins_muptheta);


      }*/
      

      bs_genie_pm1_eff_mumom_num.SetWeightNames(fname_genie_pm1);
      bs_genie_pm1_eff_mumom_den.SetWeightNames(fname_genie_pm1);
      
      bs_genie_pm1_eff_pmom_num.SetWeightNames(fname_genie_pm1);
      bs_genie_pm1_eff_pmom_den.SetWeightNames(fname_genie_pm1);

      bs_genie_pm1_eff_muangle_num.SetWeightNames(fname_genie_pm1);
      bs_genie_pm1_eff_muangle_den.SetWeightNames(fname_genie_pm1);

      bs_genie_pm1_eff_pangle_num.SetWeightNames(fname_genie_pm1);
      bs_genie_pm1_eff_pangle_den.SetWeightNames(fname_genie_pm1);

      bs_genie_pm1_eff_thetamup_num.SetWeightNames(fname_genie_pm1);
      bs_genie_pm1_eff_thetamup_den.SetWeightNames(fname_genie_pm1);







    }
    //std::cout<<"libo test 6"<<std::endl;
    // Prepare the vector of weights to be used for bootstraps
    std::vector<double> wgts_genie_pm1;
    //ofstream outfile;
    //outfile.open("rewght_genie_rwfac.txt");
    if (!isdata && _fill_bootstrap_genie) {
      for (size_t i = 0; i < t->evtwgt_genie_pm1_weight.size(); i++) {
        wgts_genie_pm1.push_back(t->evtwgt_genie_pm1_weight.at(i).at(0));
        wgts_genie_pm1.push_back(t->evtwgt_genie_pm1_weight.at(i).at(1));
        //outfile<<setw(20)<<std::setprecision(12)<<t->evtwgt_genie_pm1_weight.at(i).at(0)<<setw(20)<<std::setprecision(12)<<t->evtwgt_genie_pm1_weight.at(i).at(1)<<std::endl;
      }
    }
    //outfile.close();

    // ************************
    //
    // Set weight names, prepare bootstraps -- GENIE MULTISIM
    //
    // ************************


    if (i == _initial_entry && !isdata && _fill_bootstrap_genie) {

      if (t->evtwgt_genie_multisim_nfunc == 1) {

        if (t->evtwgt_genie_multisim_funcname.at(0) != "genie_all_Genie") {
          std::cout << "GENIE Multisim: func name is " << t->evtwgt_genie_multisim_funcname.at(0) 
                    << " which is different than genie_all" << std::endl;
        }

        fname_genie_multisim.clear();
        fname_genie_multisim.resize(t->evtwgt_genie_multisim_nweight.at(0));

        std::ostringstream oss;
        for (size_t i_wgt = 0; i_wgt < fname_genie_multisim.size(); i_wgt++) {
          oss.str("");
          oss << "universe" << i_wgt;
          fname_genie_multisim.at(i_wgt) = oss.str();
        }

        LOG_NORMAL() << "GENIE Multisim Number of universes: " << fname_genie_multisim.size() << std::endl;


        // Number of events
        for (auto & iter : _event_histo_1d->hmap_trkmom_genie_multisim_bs /*map_bs_trkmom_genie_multisim*/) {


          std::string this_name = iter.first;

          // Now emplace the histograms for the variations
          for (size_t i = 0; i < fname_genie_multisim.size(); i++) {

            // Single - Muon Momentum
            std::string histo_name = "h_genie_multisim_trkmom_" + this_name + "_" + fname_genie_multisim.at(i);
            _event_histo_1d->hmap_trkmom_genie_multisim_bs[this_name][fname_genie_multisim.at(i)] = new TH1D(histo_name.c_str(), "; Track momentum;", n_bins_mumom, bins_mumom);

            // Signle - Muon Angle
            histo_name = "h_genie_multisim_trkangle_" + this_name + "_" + fname_genie_multisim.at(i); 
            _event_histo_1d->hmap_trkangle_genie_multisim_bs[this_name][fname_genie_multisim.at(i)] = new TH1D(histo_name.c_str(), "; Track angle;", n_bins_mucostheta, bins_mucostheta); 

            // Total
            histo_name = "h_genie_multisim_onebin_" + this_name + "_" + fname_genie_multisim.at(i); 
            _event_histo_1d->hmap_onebin_genie_multisim_bs[this_name][fname_genie_multisim.at(i)] = new TH1D(histo_name.c_str(), "; Track energy;", 1, 0, 1);

            // Double Diff
            histo_name = "h_genie_multisim_trkmom_trkangle_" + this_name + "_" + fname_genie_multisim.at(i); 
            _event_histo->hmap_trktheta_trkmom_genie_multisim_bs[this_name][fname_genie_multisim.at(i)] = new TH2D(histo_name.c_str(), "; Track angle;", n_bins_double_mucostheta, bins_double_mucostheta, n_bins_double_mumom, bins_double_mumom);

            // Double Diff (polybins)
            histo_name = "h_poly_genie_multisim_trkmom_trkangle_" + this_name + "_" + fname_genie_multisim.at(i); 
            _event_histo->hmap_trktheta_trkmom_poly_genie_multisim_bs[this_name][fname_genie_multisim.at(i)] = new UBTH2Poly(histo_name.c_str(), "; Track angle;", -1.0, 1.0, 0.0, 2.5);
          }

        }


        for (auto & iter : _event_histo_1d->hmap_trkpmom_genie_multisim_bs ){
          std::string this_name = iter.first;
          for(size_t i = 0; i < fname_genie_multisim.size(); i++){
           // Single - Proton Momentum
            std::string histo_name = "h_genie_multisim_trkpmom_" + this_name + "_"+fname_genie_multisim.at(i);
            _event_histo_1d->hmap_trkpmom_genie_multisim_bs[this_name][fname_genie_multisim.at(i)] = new TH1D(histo_name.c_str(), "; Track momentum;", n_bins_pmom, bins_pmom);

            // Signle - Proton Angle
            histo_name = "h_genie_multisim_trkpangle_" + this_name + "_" + fname_genie_multisim.at(i); 
            _event_histo_1d->hmap_trkpangle_genie_multisim_bs[this_name][fname_genie_multisim.at(i)] = new TH1D(histo_name.c_str(), "; Track angle;", n_bins_pcostheta, bins_pcostheta); 

            
            histo_name = "h_genie_multisim_thetamup_" + this_name + "_"+fname_genie_multisim.at(i);
            _event_histo_1d->hmap_thetamup_genie_multisim_bs[this_name][fname_genie_multisim.at(i)] = new TH1D(histo_name.c_str(), "; Thetamup;", n_bins_muptheta, bins_muptheta);
           }
        }

        
        for (size_t i = 0; i < fname_genie_multisim.size(); i++) {

          // Normal Bins
          std::string histo_name;
          histo_name = "bs_genie_multisim_reco_per_true_" + fname_genie_multisim.at(i);
          _event_histo->bs_genie_multisim_reco_per_true[fname_genie_multisim.at(i)].resize(n_bins_double_mucostheta, std::vector<TH2D*>(n_bins_double_mumom));
          
          for (int m = 0; m < n_bins_double_mucostheta; m++) {
            for (int n = 0; n < n_bins_double_mumom; n++) { 
              std::stringstream sstm;
              sstm << histo_name << "_" << m << "_" << n;
              _event_histo->bs_genie_multisim_reco_per_true[fname_genie_multisim.at(i)][m][n] = new TH2D(sstm.str().c_str(), "reco_per_true", n_bins_double_mucostheta, bins_double_mucostheta, n_bins_double_mumom, bins_double_mumom);
            }
          }


          // Poly bins
          histo_name = "bs_genie_multisim_poly_reco_per_true_" + fname_genie_multisim.at(i);
          _event_histo->bs_genie_multisim_poly_reco_per_true[fname_genie_multisim.at(i)].resize(_n_poly_bins);
          
          for (int m = 0; m < _n_poly_bins; m++) {
            _event_histo->bs_genie_multisim_poly_reco_per_true[fname_genie_multisim.at(i)].at(m).resize(_n_poly_bins, 0.);
          }

        }

      }

      _event_histo_1d->bs_genie_multisim_eff_onebin_num->SetWeightNames(fname_genie_multisim);
      _event_histo_1d->bs_genie_multisim_eff_onebin_den->SetWeightNames(fname_genie_multisim);

      _event_histo_1d->bs_genie_multisim_eff_mumom_num->SetWeightNames(fname_genie_multisim);  //set weight names for muon momentum num
      _event_histo_1d->bs_genie_multisim_eff_mumom_den->SetWeightNames(fname_genie_multisim);  //set weight names for muon momentum den

      _event_histo_1d->bs_genie_multisim_eff_muangle_num->SetWeightNames(fname_genie_multisim); //set weight names for muon angle num
      _event_histo_1d->bs_genie_multisim_eff_muangle_den->SetWeightNames(fname_genie_multisim); //set weight names for muon angle num

      _event_histo_1d->bs_genie_multisim_true_reco_mumom->SetWeightNames(fname_genie_multisim);  //set weight names for true-reco muon momentum 
      _event_histo_1d->bs_genie_multisim_true_reco_muangle->SetWeightNames(fname_genie_multisim); //set weight names for true-reco muon angle

      _event_histo->bs_genie_multisim_eff_muangle_mumom_num->SetWeightNames(fname_genie_multisim);
      _event_histo->bs_genie_multisim_eff_muangle_mumom_den->SetWeightNames(fname_genie_multisim);

      _event_histo->bs_genie_multisim_eff_poly_muangle_mumom_num->SetWeightNames(fname_genie_multisim);
      _event_histo->bs_genie_multisim_eff_poly_muangle_mumom_den->SetWeightNames(fname_genie_multisim);

      _event_histo_1d->bs_genie_multisim_eff_pmom_num->SetWeightNames(fname_genie_multisim);   //set weight names for proton momentum num
      _event_histo_1d->bs_genie_multisim_eff_pmom_den->SetWeightNames(fname_genie_multisim);   //set weight names for proton momentum den

      _event_histo_1d->bs_genie_multisim_eff_pangle_num->SetWeightNames(fname_genie_multisim); //set weight names for proton angle num
      _event_histo_1d->bs_genie_multisim_eff_pangle_den->SetWeightNames(fname_genie_multisim); //set weight names for proton angle den

      _event_histo_1d->bs_genie_multisim_eff_thetamup_num->SetWeightNames(fname_genie_multisim); //set weight names for thetamup num
      _event_histo_1d->bs_genie_multisim_eff_thetamup_den->SetWeightNames(fname_genie_multisim); //set weight names for thetamup den

      _event_histo_1d->bs_genie_multisim_true_reco_pmom->SetWeightNames(fname_genie_multisim);  //set weight names for true-reco proton momentum
      _event_histo_1d->bs_genie_multisim_true_reco_pangle->SetWeightNames(fname_genie_multisim);  //set weight names for true-reco proton angle
      _event_histo_1d->bs_genie_multisim_true_reco_thetamup->SetWeightNames(fname_genie_multisim);  //set weight names for true-reco thetamup 




    }

    // Prepare the vector of weights to be used for bootstraps
    std::vector<double> wgts_genie_multisim;
    if (!isdata && _fill_bootstrap_genie) {
      for (size_t i_wgt = 0; i_wgt < fname_genie_multisim.size(); i_wgt++) {
        double wgt = t->evtwgt_genie_multisim_weight.at(0).at(i_wgt);
        if (wgt>100 || wgt < 0){
           wgt = 1.;
        }
        wgts_genie_multisim.push_back(wgt);
      }
    }





    // ************************
    //
    // Set weight names, prepare bootstraps -- EXTRA SYSTS
    //
    // ************************

    if (i == _initial_entry && !isdata && _fill_bootstrap_extra_syst) {


        fname_extra_syst.clear();
        fname_extra_syst.resize(100/*t->evtwgt_extra_syst_multisim_nweight.at(i_func)*/);

        std::ostringstream oss;
        for (size_t i_wgt = 0; i_wgt < fname_extra_syst.size(); i_wgt++) {
          oss.str("");
          oss << "universe" << i_wgt;
          fname_extra_syst.at(i_wgt) = oss.str();
        }

        LOG_NORMAL() << "EXTRA SYST Number of universes: " << fname_extra_syst.size() << std::endl;

        // Number of events
        for (auto & iter : _event_histo_1d->hmap_trkmom_extra_syst_multisim_bs /*map_bs_trkmom_extra_syst*/) {


          std::string this_name = iter.first;

          // Now emplace the histograms for the variations
          for (size_t i = 0; i < fname_extra_syst.size(); i++) {

            // Single - Muon Momentum
            std::string histo_name = "h_extra_multisim_trkmom_" + this_name + "_" + fname_extra_syst.at(i);
            _event_histo_1d->hmap_trkmom_extra_syst_multisim_bs[this_name][fname_extra_syst.at(i)] = new TH1D(histo_name.c_str(), "; Track length;", n_bins_mumom, bins_mumom);

            // Signle - Muon Angle
            histo_name = "h_extra_syst_trkangle_" + this_name + "_" + fname_extra_syst.at(i); 
            _event_histo_1d->hmap_trkangle_extra_syst_multisim_bs[this_name][fname_extra_syst.at(i)] = new TH1D(histo_name.c_str(), "; Track angle;", n_bins_mucostheta, bins_mucostheta); 

            // Total
            histo_name = "h_extra_syst_onebin_" + this_name + "_" + fname_extra_syst.at(i); 
            _event_histo_1d->hmap_onebin_extra_syst_multisim_bs[this_name][fname_extra_syst.at(i)] = new TH1D(histo_name.c_str(), "; Track angle;", 1, 0, 1);

            // Double Diff
            histo_name = "h_extra_syst_trkmom_trkangle_" + this_name + "_" + fname_extra_syst.at(i); 
            _event_histo->hmap_trktheta_trkmom_extra_syst_multisim_bs[this_name][fname_extra_syst.at(i)] = new TH2D(histo_name.c_str(), "; Track angle;", n_bins_double_mucostheta, bins_double_mucostheta, n_bins_double_mumom, bins_double_mumom);

            histo_name = "h_poly_extra_syst_multisim_trkmom_trkangle_" + this_name + "_" + fname_extra_syst.at(i); 
            _event_histo->hmap_trktheta_trkmom_poly_extra_syst_multisim_bs[this_name][fname_extra_syst.at(i)] = new UBTH2Poly(histo_name.c_str(), "; Track angle;", -1.0, 1.0, 0.0, 2.5);
          }

          for (size_t i = 0; i < fname_extra_syst.size(); i++) {

            std::string histo_name;
            histo_name = "bs_extra_syst_multisim_reco_per_true_" + fname_extra_syst.at(i);
            _event_histo->bs_extra_syst_multisim_reco_per_true[fname_extra_syst.at(i)].resize(n_bins_double_mucostheta, std::vector<TH2D*>(n_bins_double_mumom));
          
            for (int m = 0; m < n_bins_double_mucostheta; m++) {
              for (int n = 0; n < n_bins_double_mumom; n++) { 
                std::stringstream sstm;
                sstm << histo_name << "_" << m << "_" << n;
                _event_histo->bs_extra_syst_multisim_reco_per_true[fname_extra_syst.at(i)][m][n] = new TH2D(sstm.str().c_str(), "reco_per_true", n_bins_double_mucostheta, bins_double_mucostheta, n_bins_double_mumom, bins_double_mumom);
              }
            }

            // Poly bins
            histo_name = "bs_extra_syst_multisim_poly_reco_per_true_" + fname_extra_syst.at(i);
            _event_histo->bs_extra_syst_multisim_poly_reco_per_true[fname_extra_syst.at(i)].resize(_n_poly_bins);
          
            for (int m = 0; m < _n_poly_bins; m++) {
              _event_histo->bs_extra_syst_multisim_poly_reco_per_true[fname_extra_syst.at(i)].at(m).resize(_n_poly_bins, 0.);
            }

         }

        }
        //std::cout<<"Start loop over the event_histo_1d->hmap_trkpmom_extra_syst_multisim_bs "<<std::endl;
        for (auto & iter : _event_histo_1d->hmap_trkpmom_extra_syst_multisim_bs /*map_bs_trkmom_extra_syst*/) {


          std::string this_name = iter.first;

          // Now emplace the histograms for the variations
          for (size_t i = 0; i < fname_extra_syst.size(); i++) {

            // Single - Proton Momentum
            std::string histo_name = "h_extra_multisim_trkpmom_" + this_name + "_" + fname_extra_syst.at(i);
            _event_histo_1d->hmap_trkpmom_extra_syst_multisim_bs[this_name][fname_extra_syst.at(i)] = new TH1D(histo_name.c_str(), "; Track length;", n_bins_pmom, bins_pmom);

            // Signle - Proton Angle
            histo_name = "h_extra_syst_trkpangle_" + this_name + "_" + fname_extra_syst.at(i); 
            _event_histo_1d->hmap_trkpangle_extra_syst_multisim_bs[this_name][fname_extra_syst.at(i)] = new TH1D(histo_name.c_str(), "; Track angle;", n_bins_pcostheta, bins_pcostheta); 

            // Signle - Proton Muon Angle
            histo_name = "h_extra_syst_thetamup_" + this_name + "_" + fname_extra_syst.at(i); 
            _event_histo_1d->hmap_thetamup_extra_syst_multisim_bs[this_name][fname_extra_syst.at(i)] = new TH1D(histo_name.c_str(), "; Track angle;", n_bins_muptheta, bins_muptheta); 

            

          }
        }
      //std::cout<<"End of loop over the event_histo_1d->hmap_trkpmom_extra_syst_multisim_bs "<<std::endl;
      _event_histo_1d->bs_extra_syst_multisim_eff_onebin_num->SetWeightNames(fname_extra_syst);
      _event_histo_1d->bs_extra_syst_multisim_eff_onebin_den->SetWeightNames(fname_extra_syst);

      _event_histo_1d->bs_extra_syst_multisim_eff_mumom_num->SetWeightNames(fname_extra_syst);
      _event_histo_1d->bs_extra_syst_multisim_eff_mumom_den->SetWeightNames(fname_extra_syst);

      _event_histo_1d->bs_extra_syst_multisim_eff_muangle_num->SetWeightNames(fname_extra_syst);
      _event_histo_1d->bs_extra_syst_multisim_eff_muangle_den->SetWeightNames(fname_extra_syst);

      _event_histo_1d->bs_extra_syst_true_reco_mumom->SetWeightNames(fname_extra_syst);
      _event_histo_1d->bs_extra_syst_true_reco_muangle->SetWeightNames(fname_extra_syst);

      _event_histo->bs_extra_syst_multisim_eff_muangle_mumom_num->SetWeightNames(fname_extra_syst);
      _event_histo->bs_extra_syst_multisim_eff_muangle_mumom_den->SetWeightNames(fname_extra_syst);

      _event_histo->bs_extra_syst_multisim_eff_poly_muangle_mumom_num->SetWeightNames(fname_extra_syst);
      _event_histo->bs_extra_syst_multisim_eff_poly_muangle_mumom_den->SetWeightNames(fname_extra_syst);

      _event_histo_1d->bs_extra_syst_multisim_eff_pmom_num->SetWeightNames(fname_extra_syst);
      _event_histo_1d->bs_extra_syst_multisim_eff_pmom_den->SetWeightNames(fname_extra_syst);
      
      _event_histo_1d->bs_extra_syst_multisim_eff_pangle_num->SetWeightNames(fname_extra_syst);
      _event_histo_1d->bs_extra_syst_multisim_eff_pangle_den->SetWeightNames(fname_extra_syst);

      _event_histo_1d->bs_extra_syst_multisim_eff_thetamup_num->SetWeightNames(fname_extra_syst);
      _event_histo_1d->bs_extra_syst_multisim_eff_thetamup_den->SetWeightNames(fname_extra_syst);

      _event_histo_1d->bs_extra_syst_true_reco_pmom->SetWeightNames(fname_extra_syst);
      _event_histo_1d->bs_extra_syst_true_reco_pangle->SetWeightNames(fname_extra_syst);
      _event_histo_1d->bs_extra_syst_true_reco_thetamup->SetWeightNames(fname_extra_syst);
    }

    // Prepare the vector of weights to be used for bootstraps
    std::vector<double> wgts_extra_syst;
    wgts_extra_syst.clear();
    wgts_extra_syst.resize(fname_extra_syst.size(), 1.);

    if (!isdata && _fill_bootstrap_extra_syst) {

      bool keep_all = false;
      if (_extra_syst_target_syst == "total") {
        keep_all = true;
      }

      // Loop over all the flux reweighting function names and find the one we want unlsee "total" was requested
      for (size_t i_func = 0; i_func < t->evtwgt_extra_syst_multisim_funcname.size(); i_func++) {

        std::string func_name = t->evtwgt_extra_syst_multisim_funcname.at(i_func);

        size_t found = std::string::npos;

        if (keep_all) {
          found = 0;
        } else {
          found = func_name.find(_extra_syst_target_syst);
        }

        if (found == std::string::npos) {
          continue;
        }

        if (i == _initial_entry) LOG_NORMAL() << "Filling bootstraps for extra systematic " << func_name << std::endl;

        // Always exclude the bnbcorrection weight, this is not a systematic, though should be applied to every event
        if (func_name == "bnbcorrection_FluxHist") {
          continue;
        }

        for (size_t i_wgt = 0; i_wgt < fname_extra_syst.size(); i_wgt++) {

          double wgt = t->evtwgt_extra_syst_multisim_weight.at(i_func).at(i_wgt);
          if(wgt > 100 || wgt < 0){
             wgt = 1.;
          }
          wgts_extra_syst.at(i_wgt) *= wgt;
        }
      }
    }

/*
*        0 *        0 *                                                                               bnbcorrection_FluxHist *
*        0 *        1 *                                                                     model_q0q3_ccmec_HistogramWeight *
*        0 *        2 *                                                                      model_q0q3_ccqe_HistogramWeight *
*/




    // ************************
    //
    // Set weight names, prepare bootstraps -- FLUX MULTISIM
    //
    // ************************

    if (i == _initial_entry && !isdata && _fill_bootstrap_flux) {

      fname_flux_multisim.clear();
      fname_flux_multisim.resize(t->evtwgt_flux_multisim_nweight.at(1));

      std::ostringstream oss;
      for (size_t i_wgt = 0; i_wgt < fname_flux_multisim.size(); i_wgt++) {
       oss.str("");
        oss << "universe" << i_wgt;
        fname_flux_multisim.at(i_wgt) = oss.str();
      }

      LOG_NORMAL()  << "FLUX Multisim Number of universes: " << fname_flux_multisim.size() << std::endl;

      // Number of events
      for (auto iter : _event_histo_1d->hmap_trkmom_flux_multisim_bs) {

        std::string this_name = iter.first;
        std::map<std::string, TH1D*> bs_map = iter.second;

        // Now emplace the histograms for the variations
        for (size_t i = 0; i < fname_flux_multisim.size(); i++) {

          // Single diff - Momentum
          std::string histo_name = "h_flux_multisim_trkmom_" + this_name + "_" + fname_flux_multisim.at(i);
          _event_histo_1d->hmap_trkmom_flux_multisim_bs[this_name][fname_flux_multisim.at(i)] = new TH1D(histo_name.c_str(), "; Track length;", n_bins_mumom, bins_mumom); 

          // Single diff - Angle
          histo_name = "h_flux_multisim_trkangle_" + this_name + "_" + fname_flux_multisim.at(i);
          _event_histo_1d->hmap_trkangle_flux_multisim_bs[this_name][fname_flux_multisim.at(i)] = new TH1D(histo_name.c_str(), "; Track angle;", n_bins_mucostheta, bins_mucostheta); 

          // Total
          histo_name = "h_flux_multisim_onebin_" + this_name + "_" + fname_flux_multisim.at(i);
          _event_histo_1d->hmap_onebin_flux_multisim_bs[this_name][fname_flux_multisim.at(i)] = new TH1D(histo_name.c_str(), "; Track angle;", 1, 0, 1);

          // Double Diff
          histo_name = "h_flux_multisim_trkmom_trkangle_" + this_name + "_" + fname_flux_multisim.at(i); 
          _event_histo->hmap_trktheta_trkmom_flux_multisim_bs[this_name][fname_flux_multisim.at(i)] = new TH2D(histo_name.c_str(), "; Track angle;", n_bins_double_mucostheta, bins_double_mucostheta, n_bins_double_mumom, bins_double_mumom);
 
          // Double Diff PolyBin
          histo_name = "h_poly_flux_multisim_trkmom_trkangle_" + this_name + "_" + fname_flux_multisim.at(i); 
          _event_histo->hmap_trktheta_trkmom_poly_flux_multisim_bs[this_name][fname_flux_multisim.at(i)] = new UBTH2Poly(histo_name.c_str(), "; Track angle;", -1.0, 1.0, 0.0, 2.5);

        }

      }

      for (auto iter : _event_histo_1d->hmap_trkpmom_flux_multisim_bs){
         std::string this_name = iter.first;
         std::map<std::string, TH1D*> bs_map = iter.second;
         for (size_t i = 0; i < fname_flux_multisim.size(); i++){
           //Single diff - Momentum
           std::string histo_name = "h_flux_multisim_trkpmom_" + this_name + "_" + fname_flux_multisim.at(i);
           _event_histo_1d->hmap_trkpmom_flux_multisim_bs[this_name][fname_flux_multisim.at(i)] = new TH1D(histo_name.c_str(), "; Track length;", n_bins_pmom, bins_pmom);    

          // Single diff - Angle
          histo_name = "h_flux_multisim_trkpangle_" + this_name + "_" + fname_flux_multisim.at(i);
          _event_histo_1d->hmap_trkpangle_flux_multisim_bs[this_name][fname_flux_multisim.at(i)] = new TH1D(histo_name.c_str(), "; Track angle;", n_bins_pcostheta, bins_pcostheta); 

          // Single diff - Proton Muon Angle
          histo_name = "h_flux_multisim_thetamup_" + this_name + "_" + fname_flux_multisim.at(i);
          _event_histo_1d->hmap_thetamup_flux_multisim_bs[this_name][fname_flux_multisim.at(i)] = new TH1D(histo_name.c_str(), "; Track angle;", n_bins_muptheta, bins_muptheta); 

         }

      }
      // Reco per true
      for (size_t i = 0; i < fname_flux_multisim.size(); i++) {

        // Normal Bins
        std::string histo_name;
        histo_name = "bs_flux_multisim_reco_per_true_" + fname_flux_multisim.at(i);
        _event_histo->bs_flux_multisim_reco_per_true[fname_flux_multisim.at(i)].resize(n_bins_double_mucostheta, std::vector<TH2D*>(n_bins_double_mumom));

        for (int m = 0; m < n_bins_double_mucostheta; m++) {
          for (int n = 0; n < n_bins_double_mumom; n++) { 
            std::stringstream sstm;
            sstm << histo_name << "_" << m << "_" << n;
            _event_histo->bs_flux_multisim_reco_per_true[fname_flux_multisim.at(i)][m][n] = new TH2D(sstm.str().c_str(), "reco_per_true", n_bins_double_mucostheta, bins_double_mucostheta, n_bins_double_mumom, bins_double_mumom);
          }
        }

        // Poly bins
        histo_name = "bs_flux_multisim_poly_reco_per_true_" + fname_flux_multisim.at(i);
        _event_histo->bs_flux_multisim_poly_reco_per_true[fname_flux_multisim.at(i)].resize(_n_poly_bins);

        for (int m = 0; m < _n_poly_bins; m++) {
          std::stringstream sstm;
          sstm << histo_name << "_" << m;
          _event_histo->bs_flux_multisim_poly_reco_per_true[fname_flux_multisim.at(i)].at(m).resize(_n_poly_bins, 0.);
        }

      }

      // Efficiency
      
      _event_histo_1d->bs_flux_multisim_eff_onebin_num->SetWeightNames(fname_flux_multisim);
      _event_histo_1d->bs_flux_multisim_eff_onebin_den->SetWeightNames(fname_flux_multisim);

      _event_histo_1d->bs_flux_multisim_eff_mumom_num->SetWeightNames(fname_flux_multisim);
      _event_histo_1d->bs_flux_multisim_eff_mumom_den->SetWeightNames(fname_flux_multisim);

      _event_histo_1d->bs_flux_multisim_eff_muangle_num->SetWeightNames(fname_flux_multisim);
      _event_histo_1d->bs_flux_multisim_eff_muangle_den->SetWeightNames(fname_flux_multisim);

      _event_histo_1d->bs_flux_multisim_true_reco_mumom->SetWeightNames(fname_flux_multisim);
      _event_histo_1d->bs_flux_multisim_true_reco_muangle->SetWeightNames(fname_flux_multisim);

      _event_histo->bs_flux_multisim_eff_muangle_mumom_num->SetWeightNames(fname_flux_multisim);
      _event_histo->bs_flux_multisim_eff_muangle_mumom_den->SetWeightNames(fname_flux_multisim);

      _event_histo->bs_flux_multisim_eff_poly_muangle_mumom_num->SetWeightNames(fname_flux_multisim);
      _event_histo->bs_flux_multisim_eff_poly_muangle_mumom_den->SetWeightNames(fname_flux_multisim);
      
      _event_histo_1d->bs_flux_multisim_eff_pmom_num->SetWeightNames(fname_flux_multisim);
      _event_histo_1d->bs_flux_multisim_eff_pmom_den->SetWeightNames(fname_flux_multisim);

      _event_histo_1d->bs_flux_multisim_eff_pangle_num->SetWeightNames(fname_flux_multisim);
      _event_histo_1d->bs_flux_multisim_eff_pangle_den->SetWeightNames(fname_flux_multisim);

      _event_histo_1d->bs_flux_multisim_eff_thetamup_num->SetWeightNames(fname_flux_multisim);
      _event_histo_1d->bs_flux_multisim_eff_thetamup_den->SetWeightNames(fname_flux_multisim);

      _event_histo_1d->bs_flux_multisim_true_reco_pmom->SetWeightNames(fname_flux_multisim);
      _event_histo_1d->bs_flux_multisim_true_reco_pangle->SetWeightNames(fname_flux_multisim);
      _event_histo_1d->bs_flux_multisim_true_reco_thetamup->SetWeightNames(fname_flux_multisim);

    }

    // Prepare the vector of weights to be used for bootstraps
    std::vector<double> wgts_flux_multisim;
    wgts_flux_multisim.clear();
    wgts_flux_multisim.resize(fname_flux_multisim.size(), 1.);

    if (!isdata && _fill_bootstrap_flux) {

      bool keep_all = false;
      if (_target_flux_syst == "total") {
        keep_all = true;
      }

      // Loop over all the flux reweighting function names and find the one we want unless "total" was requested
      for (size_t i_func = 0; i_func < t->evtwgt_flux_multisim_funcname.size(); i_func++) {

        std::string func_name = t->evtwgt_flux_multisim_funcname.at(i_func);

        size_t found = std::string::npos;

        if (keep_all) {
          found = 0;
        } else {
          found = func_name.find(_target_flux_syst);
        }

        if (found == std::string::npos) {
          continue;
        }

        // Always exclude the bnbcorrection weight, this is not a systematic, though should be applied to every event
        if (func_name == "bnbcorrection_FluxHist") {
          continue;
        }

        if (i == _initial_entry) LOG_NORMAL() << "Filling bootstraps for flux systematic " << func_name << std::endl;
 
        for (size_t i_wgt = 0; i_wgt < fname_flux_multisim.size(); i_wgt++) {
          
          double wgt = t->evtwgt_flux_multisim_weight.at(i_func).at(i_wgt);
          if (wgt>100 || wgt < 0){
            wgt = 1.;
          }
          wgts_flux_multisim.at(i_wgt) *= wgt;
          if (_reweigh_kaons
             && (t->evtwgt_flux_multisim_funcname.at(i_func) == "kminus_PrimaryHadronNormalization" 
             || t->evtwgt_flux_multisim_funcname.at(i_func) == "kplus_PrimaryHadronFeynmanScaling" 
             || t->evtwgt_flux_multisim_funcname.at(i_func) == "kzero_PrimaryHadronSanfordWang")) {
            // std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
            if (t->evtwgt_flux_multisim_weight.at(i_func).at(i_wgt) != 1) {
              is_from_kaon = true;
            }
          }
        }
      }
    }

    if (is_from_kaon && _reweigh_kaons) {
          event_weight *= _kaon_reweigh_factor;
    }


/*
*        0 *        0 *         bnbcorrection_FluxHist *                              1 *
*        0 *        1 *             expskin_FluxUnisim *                            100 *
*        0 *        2 *         horncurrent_FluxUnisim *                            100 *
*        0 *        3 * kminus_PrimaryHadronNormalizat *                            100 *
*        0 *        4 * kplus_PrimaryHadronFeynmanScal *                            100 *
*        0 *        5 * kzero_PrimaryHadronSanfordWang *                            100 *
*        0 *        6 *      nucleoninexsec_FluxUnisim *                            100 *
*        0 *        7 *       nucleonqexsec_FluxUnisim *                            100 *
*        0 *        8 *      nucleontotxsec_FluxUnisim *                            100 *
*        0 *        9 * piminus_PrimaryHadronSWCentral *                            100 *
*        0 *       10 *         pioninexsec_FluxUnisim *                            100 *
*        0 *       11 *          pionqexsec_FluxUnisim *                            100 *
*        0 *       12 *         piontotxsec_FluxUnisim *                            100 *
*        0 *       13 * piplus_PrimaryHadronSWCentralS *                            100 *
*/


    

    // ************************
    //
    // Set weight names, prepare bootstraps -- MC STAT MULTISIM
    //
    // ************************

    if (i == _initial_entry && !isdata && _fill_bootstrap_mc_stat) {

      if (_mc_stat_n_events > 0) {

        fname_mc_stat_multisim.clear();
        fname_mc_stat_multisim.resize(_mc_stat_n_events);

        std::ostringstream oss;
        for (size_t i_wgt = 0; i_wgt < fname_mc_stat_multisim.size(); i_wgt++) {
          oss.str("");
          oss << "universe" << i_wgt;
          fname_mc_stat_multisim.at(i_wgt) = oss.str();
        }

        LOG_NORMAL()  << "MC STAT Multisim Number of universes: " << fname_mc_stat_multisim.size() << std::endl;

        // Number of events
        for (auto & iter : _event_histo_1d->hmap_trkmom_mc_stat_multisim_bs) {


          std::string this_name = iter.first;

          // Now emplace the histograms for the variations
          for (size_t i = 0; i < fname_mc_stat_multisim.size(); i++) {

            // Single - Muon Momentum
            std::string histo_name = "h_mc_stat_multisim_trkmom_" + this_name + "_" + fname_mc_stat_multisim.at(i);
            _event_histo_1d->hmap_trkmom_mc_stat_multisim_bs[this_name][fname_mc_stat_multisim.at(i)] = new TH1D(histo_name.c_str(), "; Track length;", n_bins_mumom, bins_mumom);

            // Signle - Muon Angle
            histo_name = "h_mc_stat_multisim_trkangle_" + this_name + "_" + fname_mc_stat_multisim.at(i); 
            _event_histo_1d->hmap_trkangle_mc_stat_multisim_bs[this_name][fname_mc_stat_multisim.at(i)] = new TH1D(histo_name.c_str(), "; Track angle;", n_bins_mucostheta, bins_mucostheta); 

            // Total
            histo_name = "h_mc_stat_multisim_onebin_" + this_name + "_" + fname_mc_stat_multisim.at(i); 
            _event_histo_1d->hmap_onebin_mc_stat_multisim_bs[this_name][fname_mc_stat_multisim.at(i)] = new TH1D(histo_name.c_str(), "; Track angle;", 1, 0, 1);

            // Double Diff
            histo_name = "h_mc_stat_multisim_trkmom_trkangle_" + this_name + "_" + fname_mc_stat_multisim.at(i); 
            _event_histo->hmap_trktheta_trkmom_mc_stat_multisim_bs[this_name][fname_mc_stat_multisim.at(i)] = new TH2D(histo_name.c_str(), "; Track angle;", n_bins_double_mucostheta, bins_double_mucostheta, n_bins_double_mumom, bins_double_mumom);

            // Double Diff PolyBin
            histo_name = "h_poly_flux_multisim_trkmom_trkangle_" + this_name + "_" + fname_flux_multisim.at(i); 
            _event_histo->hmap_trktheta_trkmom_poly_mc_stat_multisim_bs[this_name][fname_flux_multisim.at(i)] = new UBTH2Poly(histo_name.c_str(), "; Track angle;", -1.0, 1.0, 0.0, 2.5);
          }

        }
        for (auto & iter : _event_histo_1d->hmap_trkpmom_mc_stat_multisim_bs) {
          std::string this_name = iter.first;
          // Now emplace the histograms for the variations
          for (size_t i = 0; i < fname_mc_stat_multisim.size(); i++) {

            // Single - Proton Momentum
            std::string histo_name = "h_mc_stat_multisim_trkpmom_" + this_name + "_" + fname_mc_stat_multisim.at(i);
            _event_histo_1d->hmap_trkpmom_mc_stat_multisim_bs[this_name][fname_mc_stat_multisim.at(i)] = new TH1D(histo_name.c_str(), "; Track length;", n_bins_mumom, bins_mumom);

            // Signle - Proton Angle
            histo_name = "h_mc_stat_multisim_trkpangle_" + this_name + "_" + fname_mc_stat_multisim.at(i); 
            _event_histo_1d->hmap_trkpangle_mc_stat_multisim_bs[this_name][fname_mc_stat_multisim.at(i)] = new TH1D(histo_name.c_str(), "; Track angle;", n_bins_pcostheta, bins_pcostheta); 

            // Signle - Proton Muon Angle
            histo_name = "h_mc_stat_multisim_thetamup_" + this_name + "_" + fname_mc_stat_multisim.at(i); 
            _event_histo_1d->hmap_thetamup_mc_stat_multisim_bs[this_name][fname_mc_stat_multisim.at(i)] = new TH1D(histo_name.c_str(), "; Track angle;", n_bins_muptheta, bins_muptheta); 

          }  

        }

        // Efficiency 

        for (size_t i = 0; i < fname_mc_stat_multisim.size(); i++) {

          // Normal bins
          std::string histo_name;
          histo_name = "bs_mc_stat_multisim_reco_per_true_" + fname_mc_stat_multisim.at(i);
          _event_histo->bs_mc_stat_multisim_reco_per_true[fname_mc_stat_multisim.at(i)].resize(n_bins_double_mucostheta, std::vector<TH2D*>(n_bins_double_mumom));
          
          for (int m = 0; m < n_bins_double_mucostheta; m++) {
            for (int n = 0; n < n_bins_double_mumom; n++) { 
              std::stringstream sstm;
              sstm << histo_name << "_" << m << "_" << n;
              _event_histo->bs_mc_stat_multisim_reco_per_true[fname_mc_stat_multisim.at(i)][m][n] = new TH2D(sstm.str().c_str(), "reco_per_true", n_bins_double_mucostheta, bins_double_mucostheta, n_bins_double_mumom, bins_double_mumom);
            }
          }

          // Poly bins
          histo_name = "bs_mc_stat_multisim_poly_reco_per_true_" + fname_flux_multisim.at(i);
          _event_histo->bs_mc_stat_multisim_poly_reco_per_true[fname_flux_multisim.at(i)].resize(_n_poly_bins);

          for (int m = 0; m < _n_poly_bins; m++) {
            std::stringstream sstm;
            sstm << histo_name << "_" << m;
            _event_histo->bs_mc_stat_multisim_poly_reco_per_true[fname_flux_multisim.at(i)].at(m).resize(_n_poly_bins, 0.);
          }

        }

      }

      _event_histo_1d->bs_mc_stat_multisim_eff_onebin_num->SetWeightNames(fname_mc_stat_multisim);
      _event_histo_1d->bs_mc_stat_multisim_eff_onebin_den->SetWeightNames(fname_mc_stat_multisim);

      _event_histo_1d->bs_mc_stat_multisim_eff_mumom_num->SetWeightNames(fname_mc_stat_multisim);
      _event_histo_1d->bs_mc_stat_multisim_eff_mumom_den->SetWeightNames(fname_mc_stat_multisim);

      _event_histo_1d->bs_mc_stat_multisim_eff_muangle_num->SetWeightNames(fname_mc_stat_multisim);
      _event_histo_1d->bs_mc_stat_multisim_eff_muangle_den->SetWeightNames(fname_mc_stat_multisim);

      _event_histo_1d->bs_mc_stat_multisim_true_reco_mumom->SetWeightNames(fname_mc_stat_multisim);
      _event_histo_1d->bs_mc_stat_multisim_true_reco_muangle->SetWeightNames(fname_mc_stat_multisim);

      _event_histo->bs_mc_stat_multisim_eff_muangle_mumom_num->SetWeightNames(fname_mc_stat_multisim);
      _event_histo->bs_mc_stat_multisim_eff_muangle_mumom_den->SetWeightNames(fname_mc_stat_multisim);

      _event_histo->bs_mc_stat_multisim_eff_poly_muangle_mumom_num->SetWeightNames(fname_mc_stat_multisim);
      _event_histo->bs_mc_stat_multisim_eff_poly_muangle_mumom_den->SetWeightNames(fname_mc_stat_multisim);
      
      _event_histo_1d->bs_mc_stat_multisim_eff_pmom_num->SetWeightNames(fname_mc_stat_multisim);
      _event_histo_1d->bs_mc_stat_multisim_eff_pmom_den->SetWeightNames(fname_mc_stat_multisim);

      _event_histo_1d->bs_mc_stat_multisim_eff_pangle_num->SetWeightNames(fname_mc_stat_multisim);
      _event_histo_1d->bs_mc_stat_multisim_eff_pangle_den->SetWeightNames(fname_mc_stat_multisim);

      _event_histo_1d->bs_mc_stat_multisim_eff_thetamup_num->SetWeightNames(fname_mc_stat_multisim);
      _event_histo_1d->bs_mc_stat_multisim_eff_thetamup_den->SetWeightNames(fname_mc_stat_multisim);

      _event_histo_1d->bs_mc_stat_multisim_true_reco_pmom->SetWeightNames(fname_mc_stat_multisim);
      _event_histo_1d->bs_mc_stat_multisim_true_reco_pangle->SetWeightNames(fname_mc_stat_multisim);
      _event_histo_1d->bs_mc_stat_multisim_true_reco_thetamup->SetWeightNames(fname_mc_stat_multisim);

    }

    // Prepare the vector of weights to be used for bootstraps
    std::vector<double> wgts_mc_stat_multisim;
    wgts_mc_stat_multisim.resize(fname_mc_stat_multisim.size());
    if (!isdata && _fill_bootstrap_mc_stat) {
      for (size_t i = 0; i < fname_mc_stat_multisim.size(); i++) {
        wgts_mc_stat_multisim.at(i) = _random_engine.PoissonD(1);
      }
    }
    //********************************************



    if (i == _initial_entry) {
      _event_histo->AddPolyBins();
      LOG_NORMAL() << "Polybins set." << std::endl;
    }


    // ************************
    //
    // Preliminary - Truth
    //
    // ************************
    
    // This variable will store if this is a signal event or not
    bool isSignal = false;
    if (_ana_int_type == "ccinclusive_analysis"){
    if (t->nupdg == 14 && t->ccnc == 0 && t->fv == 1 /*&& (t->tvtx_z[0] < 675 || t->tvtx_z[0] > 775)*/){

      nsignal += event_weight;
      isSignal = true;

      if (t->mode == 0) nsignal_qe += event_weight;
      if (t->mode == 1) nsignal_res += event_weight;
      if (t->mode == 2) nsignal_dis += event_weight;
      if (t->mode == 3) nsignal_coh += event_weight;
      if (t->mode == 10) nsignal_mec += event_weight;
      
    }
    }//end of if this is ccinclusive analysis
    if (_ana_int_type == "cc1unp_analysis"){
    if (t->nupdg == 14 && t->ccnc == 0 && t->fv == 1 /*&& (t->tvtx_z[0] < 675 || t->tvtx_z[0] > 775)*/){
      if(t->ngenie_muons>0 && t->ngenie_protons_300>0 &&(t->ngenie_electrons<1 && t->ngenie_pipms<1 && t->ngenie_pion0s<1)) { 
      //find out the leading proton momentum
      float ppindex=-999.0;
      float epindex=-999.0;
      for(size_t hh=0; hh<t->genie_mcpar_pdgcode.size(); hh++){
          if(abs(t->genie_mcpar_pdgcode[hh]==2212) && t->genie_mcpar_energy[hh]>epindex) {
             epindex=t->genie_mcpar_energy[hh];
             ppindex=TMath::Sqrt(t->genie_mcpar_px[hh]*t->genie_mcpar_px[hh]+
                                 t->genie_mcpar_py[hh]*t->genie_mcpar_py[hh]+
                                 t->genie_mcpar_pz[hh]*t->genie_mcpar_pz[hh]);
          }
      }
      if(t->true_muon_mom_matched >0.1 && ppindex<1.2){
      nsignal += event_weight;
      isSignal = true;

      if (t->mode == 0) nsignal_qe += event_weight;
      if (t->mode == 1) nsignal_res += event_weight;
      if (t->mode == 2) nsignal_dis += event_weight;
      if (t->mode == 3) nsignal_coh += event_weight;
      if (t->mode == 10) nsignal_mec += event_weight;
      }
      }
    }
    }//end of if this is cc1unp analysis
     
    if(isSignal && _ana_int_type == "cc1unp_analysis"){
      // loop over all the genie particles and get the denominator of the proton momentum efficiency
      temp_pmom=-999.0;
      temp_pangle=-999.0;
      temp_pphi=-999.0;
      temp_thetamup=-999.0;
      float temp_penergy=-999.0;
      for(size_t mpar=0; mpar<t->genie_mcpar_pdgcode.size(); mpar++){
           //std::cout<<"PDGCODE of GENIE Particle is "<<t->genie_mcpar_pdgcode[mpar]<<std::endl;
           if(abs(t->genie_mcpar_pdgcode[mpar])==2212 && t->genie_mcpar_energy[mpar]>temp_penergy) {
               temp_penergy=t->genie_mcpar_energy[mpar];
               temp_pmom=TMath::Sqrt(t->genie_mcpar_px[mpar]*t->genie_mcpar_px[mpar]+
                              t->genie_mcpar_py[mpar]*t->genie_mcpar_py[mpar]+
                              t->genie_mcpar_pz[mpar]*t->genie_mcpar_pz[mpar]);
               temp_pangle=t->genie_mcpar_pz[mpar]/temp_pmom;
               temp_pphi=TMath::ACos(t->genie_mcpar_px[mpar]/TMath::Sqrt(t->genie_mcpar_px[mpar]*t->genie_mcpar_px[mpar]+t->genie_mcpar_py[mpar]*t->genie_mcpar_py[mpar]));
               temp_thetamup=getAngle(temp_pmom, TMath::ACos(temp_pangle), temp_pphi,  t->true_muon_mom_matched, TMath::ACos(t->lep_costheta), t->lep_phi); 
           }
      }      
      _event_histo_1d->h_eff_pmom_den->Fill(temp_pmom, event_weight);
      _event_histo_1d->h_eff_pangle_den->Fill(temp_pangle, event_weight);
      //get the den of thetamup
      _event_histo_1d->h_eff_thetamup_den->Fill(temp_thetamup, event_weight);
      
      h_eff_allprotons_den->Fill(temp_pmom, event_weight);
      //loop over g4 stage and select the contained cc1unp and reinteracted or non reinteracted protons
      //vector<double>  geant_mcpar_pdgcode;
      temp_geantpmom=-999.0;
      temp_geantpenergy=-999.0;
      bool geant_pcontained=true;
      for(size_t g4par=0; g4par<t->geant_mcpar_pdgcode.size(); g4par++){
         if(t->geant_mcpar_pdgcode[g4par] != 2212) continue;
         if(!inCV(t->geant_mcpar_startx[g4par], t->geant_mcpar_starty[g4par], t->geant_mcpar_startz[g4par]) || !inCV(t->geant_mcpar_endx[g4par], t->geant_mcpar_endy[g4par], t->geant_mcpar_endz[g4par])) geant_pcontained=false;   
         //if(t->geant_mcpar_end_process[g4par] != "protonInelastic") continue;
         if(t->geant_mcpar_energy[g4par]>temp_geantpenergy) {
           temp_geantpenergy=t->geant_mcpar_energy[g4par];
           temp_geantpmom=TMath::Sqrt(t->geant_mcpar_px[g4par]*t->geant_mcpar_px[g4par]+
                                      t->geant_mcpar_py[g4par]*t->geant_mcpar_py[g4par]+
                                      t->geant_mcpar_pz[g4par]*t->geant_mcpar_pz[g4par]);
         }
      }
      h_eff_reintprotons_den->Fill(temp_geantpmom, event_weight);
      if(geant_pcontained==true){
         h_eff_containedprotons_den->Fill(temp_geantpmom, event_weight);
      }
      
    } //end of if this is signal and cc1unp analysis
    //std::cout<<"Calculated the proton momentum and angles at simulation level"<<std::endl;
    // Check if it's a nue event
    bool isNueCCFV = false;
    if (t->nupdg == 12 && t->ccnc == 0 && t->fv == 1) {
      isNueCCFV = true;
      nue_cc_fv+=t->bnb_weight;
    }
    bool isNue = false;
    if (t->nupdg == 12) {
      isNue = true;
    }
    //bool isOtherCC =false;
    //if (t->nupdg ==14 && t->ccnc == 0 && t->fv==1 && (t->ngenie_protons==0 || t->ngenie_pipms+t->ngenie_pion0s > 0)){ 
    //   isOtherCC = true; }
    //
    // Construct the denominator for the efficiency plots
    //
    // for (int j = 0; j < t->truth_nu_e.size(); j++) {
    //   if (t->truth_nupdg.at(j) == 14 && t->truth_ccnc.at(j) == 0 && t->truth_fv.at(j) == 1) {
    //     nsignal_all += event_weight;
    //   }
    // }
    if (isSignal) {
      
      _event_histo_1d->h_eff_onebin_den->Fill(0.5, event_weight);
      h_eff_den->Fill(t->nu_e, event_weight);
      
      if (!isdata && _fill_bootstrap_genie) bs_genie_pm1_eff_mumom_den.Fill(t->true_muon_mom, event_weight, wgts_genie_pm1);
      if (!isdata && _fill_bootstrap_genie) bs_genie_pm1_eff_pmom_den.Fill(temp_pmom, event_weight, wgts_genie_pm1);
      if (!isdata && _fill_bootstrap_genie) bs_genie_pm1_eff_muangle_den.Fill(t->lep_costheta, event_weight, wgts_genie_pm1);
      if (!isdata && _fill_bootstrap_genie) bs_genie_pm1_eff_pangle_den.Fill(temp_pangle, event_weight, wgts_genie_pm1);
      if (!isdata && _fill_bootstrap_genie) bs_genie_pm1_eff_thetamup_den.Fill(temp_thetamup, event_weight, wgts_genie_pm1);


      if (!isdata && _fill_bootstrap_genie) _event_histo_1d->bs_genie_multisim_eff_onebin_den->Fill(0.5, event_weight, wgts_genie_multisim);
      if (!isdata && _fill_bootstrap_genie) _event_histo_1d->bs_genie_multisim_eff_mumom_den->Fill(t->true_muon_mom, event_weight, wgts_genie_multisim);
      if (!isdata && _fill_bootstrap_genie) _event_histo_1d->bs_genie_multisim_eff_pmom_den->Fill(temp_pmom, event_weight, wgts_genie_multisim);
      if (!isdata && _fill_bootstrap_genie) _event_histo_1d->bs_genie_multisim_eff_muangle_den->Fill(t->lep_costheta, event_weight, wgts_genie_multisim);
      if (!isdata && _fill_bootstrap_genie) _event_histo_1d->bs_genie_multisim_eff_pangle_den->Fill(temp_pangle, event_weight, wgts_genie_multisim);
      if (!isdata && _fill_bootstrap_genie) _event_histo_1d->bs_genie_multisim_eff_thetamup_den->Fill(temp_thetamup, event_weight, wgts_genie_multisim);

      if (!isdata && _fill_bootstrap_genie) _event_histo->bs_genie_multisim_eff_muangle_mumom_den->Fill(t->lep_costheta, t->true_muon_mom, event_weight, wgts_genie_multisim);
      if (!isdata && _fill_bootstrap_genie) _event_histo->bs_genie_multisim_eff_poly_muangle_mumom_den->Fill(t->lep_costheta, t->true_muon_mom, event_weight, wgts_genie_multisim);

      if (!isdata && _fill_bootstrap_extra_syst) _event_histo_1d->bs_extra_syst_multisim_eff_onebin_den->Fill(0.5, event_weight, wgts_extra_syst);
      if (!isdata && _fill_bootstrap_extra_syst) _event_histo_1d->bs_extra_syst_multisim_eff_mumom_den->Fill(t->true_muon_mom, event_weight, wgts_extra_syst);
      if (!isdata && _fill_bootstrap_extra_syst) _event_histo_1d->bs_extra_syst_multisim_eff_pmom_den->Fill(temp_pmom, event_weight, wgts_extra_syst);
      if (!isdata && _fill_bootstrap_extra_syst) _event_histo_1d->bs_extra_syst_multisim_eff_muangle_den->Fill(t->lep_costheta, event_weight, wgts_extra_syst);
      if (!isdata && _fill_bootstrap_extra_syst) _event_histo_1d->bs_extra_syst_multisim_eff_pangle_den->Fill(temp_pangle, event_weight, wgts_extra_syst);
      if (!isdata && _fill_bootstrap_extra_syst) _event_histo_1d->bs_extra_syst_multisim_eff_thetamup_den->Fill(temp_thetamup, event_weight, wgts_extra_syst);

      if (!isdata && _fill_bootstrap_extra_syst) _event_histo->bs_extra_syst_multisim_eff_muangle_mumom_den->Fill(t->lep_costheta, t->true_muon_mom, event_weight, wgts_extra_syst);
      if (!isdata && _fill_bootstrap_extra_syst) _event_histo->bs_extra_syst_multisim_eff_poly_muangle_mumom_den->Fill(t->lep_costheta, t->true_muon_mom, event_weight, wgts_extra_syst);

      if (!isdata && _fill_bootstrap_flux) _event_histo_1d->bs_flux_multisim_eff_onebin_den->Fill(0.5, event_weight, wgts_flux_multisim);
      if (!isdata && _fill_bootstrap_flux) _event_histo_1d->bs_flux_multisim_eff_mumom_den->Fill(t->true_muon_mom, event_weight, wgts_flux_multisim);
      if (!isdata && _fill_bootstrap_flux) _event_histo_1d->bs_flux_multisim_eff_pmom_den->Fill(temp_pmom, event_weight, wgts_flux_multisim);
      if (!isdata && _fill_bootstrap_flux) _event_histo_1d->bs_flux_multisim_eff_muangle_den->Fill(t->lep_costheta, event_weight, wgts_flux_multisim);
      if (!isdata && _fill_bootstrap_flux) _event_histo_1d->bs_flux_multisim_eff_pangle_den->Fill(temp_pangle, event_weight, wgts_flux_multisim);
      if (!isdata && _fill_bootstrap_flux) _event_histo_1d->bs_flux_multisim_eff_thetamup_den->Fill(temp_thetamup, event_weight, wgts_flux_multisim);

      if (!isdata && _fill_bootstrap_flux) _event_histo->bs_flux_multisim_eff_muangle_mumom_den->Fill(t->lep_costheta, t->true_muon_mom, event_weight, wgts_flux_multisim);
      if (!isdata && _fill_bootstrap_flux) _event_histo->bs_flux_multisim_eff_poly_muangle_mumom_den->Fill(t->lep_costheta, t->true_muon_mom, event_weight, wgts_flux_multisim);

      if (!isdata && _fill_bootstrap_mc_stat) _event_histo_1d->bs_mc_stat_multisim_eff_onebin_den->Fill(0.5, event_weight, wgts_mc_stat_multisim);
      if (!isdata && _fill_bootstrap_mc_stat) _event_histo_1d->bs_mc_stat_multisim_eff_mumom_den->Fill(t->true_muon_mom, event_weight, wgts_mc_stat_multisim);
      if (!isdata && _fill_bootstrap_mc_stat) _event_histo_1d->bs_mc_stat_multisim_eff_pmom_den->Fill(temp_pmom, event_weight, wgts_mc_stat_multisim);
      if (!isdata && _fill_bootstrap_mc_stat) _event_histo_1d->bs_mc_stat_multisim_eff_muangle_den->Fill(t->lep_costheta, event_weight, wgts_mc_stat_multisim);
      if (!isdata && _fill_bootstrap_mc_stat) _event_histo_1d->bs_mc_stat_multisim_eff_pangle_den->Fill(temp_pangle, event_weight, wgts_mc_stat_multisim);
      if (!isdata && _fill_bootstrap_mc_stat) _event_histo_1d->bs_mc_stat_multisim_eff_thetamup_den->Fill(temp_thetamup, event_weight, wgts_mc_stat_multisim);

      if (!isdata && _fill_bootstrap_mc_stat) _event_histo->bs_mc_stat_multisim_eff_muangle_mumom_den->Fill(t->lep_costheta, t->true_muon_mom, event_weight, wgts_mc_stat_multisim);
      if (!isdata && _fill_bootstrap_mc_stat) _event_histo->bs_mc_stat_multisim_eff_poly_muangle_mumom_den->Fill(t->lep_costheta, t->true_muon_mom, event_weight, wgts_mc_stat_multisim);

      _event_histo_1d->h_eff_mumom_den->Fill(t->true_muon_mom, event_weight);
      _event_histo_1d->h_eff_muangle_den->Fill(t->lep_costheta, event_weight);
     
      _event_histo->h_eff_muangle_mumom_den->Fill(t->lep_costheta, t->true_muon_mom, event_weight);
      _event_histo->h_eff_muangle_mumom_poly_den->Fill(t->lep_costheta, t->true_muon_mom, event_weight);
      /*
      */



      //*******fill the den of the pmom and pangle locates in the loop to get the leading proton true momentum

      _event_histo_1d->h_eff_muphi_den->Fill(t->lep_phi, event_weight);

      h_eff_mult_den->Fill(t->genie_mult, event_weight);
      h_eff_mult_ch_den->Fill(t->genie_mult_ch, event_weight);

      if (t->mode == 0) h_eff_qe_den->Fill(t->nu_e, event_weight);
      if (t->mode == 1) h_eff_res_den->Fill(t->nu_e, event_weight);
      if (t->mode == 2) h_eff_dis_den->Fill(t->nu_e, event_weight);
      if (t->mode == 3) h_eff_coh_den->Fill(t->nu_e, event_weight);
      if (t->mode == 10) h_eff_mec_den->Fill(t->nu_e, event_weight);

      h_truth_xsec_mumom->Fill(t->true_muon_mom, event_weight);
      h_truth_xsec_muangle->Fill(t->lep_costheta, event_weight);

      h_truth_xsec_pmom->Fill(temp_pmom, event_weight);
      h_truth_xsec_pangle->Fill(temp_pangle, event_weight);

      h_truth_xsec_thetamup->Fill(temp_thetamup, event_weight);


      h_mueff_den->Fill(t->true_muon_mom, event_weight);
      h_mueff_angle_den->Fill(t->lep_costheta, event_weight);

      h_true_nu_eng_beforesel->Fill(t->nu_e, event_weight);
      
      if (t->muon_is_reco){
        h_mumom_nue->Fill(t->nu_e, t->true_muon_mom, event_weight);
        nSignalWMuonReco++;
        h_mueff_num->Fill(t->true_muon_mom, event_weight);
        h_mueff_angle_num->Fill(t->lep_costheta, event_weight);
        for (auto origin : t->slc_origin){
          if (origin == 0 || origin == 2) {
            h_mueff_2_num->Fill(t->true_muon_mom, event_weight);
            break;
          }
        }
        if (t->vtx_resolution > -1 && t->vtx_resolution < 10) nSignalMuonRecoVtxOk++;
        
        h_muon_track_eff->Fill(t->muon_reco_eff, event_weight);
        h_muon_track_pur->Fill(t->muon_reco_pur, event_weight);
        
        h_mu_eff_mom->Fill(t->true_muon_mom, t->muon_reco_eff, event_weight);
        h_mu_pur_mom->Fill(t->true_muon_mom, t->muon_reco_pur, event_weight);
      }
      else{
        //std::cout << "This is a signal event but the muon was not reconstructed. Event: " << event << std::endl;
      }


      // Also save the mc truth histogram per interaction type
      hmap_mctruth_nuenergy_gen["total"]->Fill(t->nu_e, event_weight);
      hmap_mctruth_mumom_gen["total"]->Fill(t->true_muon_mom, event_weight);
      hmap_mctruth_mucostheta_gen["total"]->Fill(t->lep_costheta, event_weight);
      hmap_mctruth_muphi_gen["total"]->Fill(t->lep_phi, event_weight);
      hmap_mctruth_chargedmult_gen["total"]->Fill(t->genie_mult_ch, event_weight);
      hmap_mctruth_mucostheta_mumom_gen["total"]->Fill(t->lep_costheta, t->true_muon_mom, event_weight);
      if (t->mode == 0) {
        hmap_mctruth_nuenergy_gen["qe"]->Fill(t->nu_e, event_weight);
        hmap_mctruth_mumom_gen["qe"]->Fill(t->true_muon_mom, event_weight);
        hmap_mctruth_mucostheta_gen["qe"]->Fill(t->lep_costheta, event_weight);
        hmap_mctruth_muphi_gen["qe"]->Fill(t->lep_phi, event_weight);
        hmap_mctruth_chargedmult_gen["qe"]->Fill(t->genie_mult_ch, event_weight);
        hmap_mctruth_mucostheta_mumom_gen["qe"]->Fill(t->lep_costheta, t->true_muon_mom, event_weight);
      }
      if (t->mode == 1) {
        hmap_mctruth_nuenergy_gen["res"]->Fill(t->nu_e, event_weight);
        hmap_mctruth_mumom_gen["res"]->Fill(t->true_muon_mom, event_weight);
        hmap_mctruth_mucostheta_gen["res"]->Fill(t->lep_costheta, event_weight);
        hmap_mctruth_muphi_gen["res"]->Fill(t->lep_phi, event_weight);
        hmap_mctruth_chargedmult_gen["res"]->Fill(t->genie_mult_ch, event_weight);
        hmap_mctruth_mucostheta_mumom_gen["res"]->Fill(t->lep_costheta, t->true_muon_mom, event_weight);
      }
      if (t->mode == 2) {
        hmap_mctruth_nuenergy_gen["dis"]->Fill(t->nu_e, event_weight);
        hmap_mctruth_mumom_gen["dis"]->Fill(t->true_muon_mom, event_weight);
        hmap_mctruth_mucostheta_gen["dis"]->Fill(t->lep_costheta, event_weight);
        hmap_mctruth_muphi_gen["dis"]->Fill(t->lep_phi, event_weight);
        hmap_mctruth_chargedmult_gen["dis"]->Fill(t->genie_mult_ch, event_weight);
        hmap_mctruth_mucostheta_mumom_gen["dis"]->Fill(t->lep_costheta, t->true_muon_mom, event_weight);
      }
      if (t->mode == 3) {
        hmap_mctruth_nuenergy_gen["coh"]->Fill(t->nu_e, event_weight);
        hmap_mctruth_mumom_gen["coh"]->Fill(t->true_muon_mom, event_weight);
        hmap_mctruth_mucostheta_gen["coh"]->Fill(t->lep_costheta, event_weight);
        hmap_mctruth_muphi_gen["coh"]->Fill(t->lep_phi, event_weight);
        hmap_mctruth_chargedmult_gen["coh"]->Fill(t->genie_mult_ch, event_weight);
        hmap_mctruth_mucostheta_mumom_gen["coh"]->Fill(t->lep_costheta, t->true_muon_mom, event_weight);
      }
      if (t->mode == 10) {
        hmap_mctruth_nuenergy_gen["mec"]->Fill(t->nu_e, event_weight);
        hmap_mctruth_mumom_gen["mec"]->Fill(t->true_muon_mom, event_weight);
        hmap_mctruth_mucostheta_gen["mec"]->Fill(t->lep_costheta, event_weight);
        hmap_mctruth_muphi_gen["mec"]->Fill(t->lep_phi, event_weight);
        hmap_mctruth_chargedmult_gen["mec"]->Fill(t->genie_mult_ch, event_weight);
        hmap_mctruth_mucostheta_mumom_gen["mec"]->Fill(t->lep_costheta, t->true_muon_mom, event_weight);
      }
    } // if is signal

    if(t->nupdg == 14 && t->ccnc == 0){
      nNumuCC++;
    }

    
    //std::cout<<"filled the histograms for signal"<<std::endl;

      
    
    // VTX before selection
    if(t->slc_nuvtx_x.size()>0) {
      hmap_vtxx_b["total"]->Fill(t->slc_nuvtx_x.at(0), event_weight);
      hmap_vtxy_b["total"]->Fill(t->slc_nuvtx_y.at(0), event_weight);
      hmap_vtxz_b["total"]->Fill(t->slc_nuvtx_z.at(0), event_weight);
      h_vtx_xz->Fill(t->slc_nuvtx_x.at(0), t->slc_nuvtx_z.at(0), event_weight);
      h_vtx_xy->Fill(t->slc_nuvtx_x.at(0), t->slc_nuvtx_y.at(0), event_weight);
    }
    
    // Number of TPCObjects
    hmap_ntpcobj["total"]->Fill(t->nslices, event_weight);
    if (isSignal) hmap_ntpcobj["signal"]->Fill(t->nslices, event_weight);
    else hmap_ntpcobj["background"]->Fill(t->nslices, event_weight);
    
    // Vertex Resolution plot
    if (isSignal) h_vtx_resolution->Fill(t->vtx_resolution, event_weight);
    
    
    for (int slc = 0; slc < t->nslices; slc ++) {
      
      if (t->slc_origin.at(slc) == 0 || t->slc_origin.at(slc) == 2) {
        h_slice_npfp->Fill(t->slc_npfp.at(slc), event_weight);
        h_slice_ntrack->Fill(t->slc_ntrack.at(slc), event_weight);
      } else {
        h_slice_npfp_others->Fill(t->slc_npfp.at(slc), event_weight);
        h_slice_ntrack_others->Fill(t->slc_ntrack.at(slc), event_weight);
      }
      

      
      if (t->slc_origin.at(slc) == 0) h_slice_origin->Fill(2., event_weight);
      if (t->slc_origin.at(slc) == 1) h_slice_origin->Fill(0., event_weight);
      if (t->slc_origin.at(slc) == 2) h_slice_origin->Fill(1., event_weight);
      
      
      if ((t->slc_origin.at(slc) == 0 || t->slc_origin.at(slc) == 2) && t->fv == 1) {
        n_slc_nu_origin ++;
        
      }
    }
    
    
    
    
    //std::cout<<"Start to perform the CCinclusive as precuts"<<std::endl; 
    
    

    
    
    
    // ************************
    //
    //  Selection
    //
    // ***********************


    if (isSignal) selected_signal_events_percut["initial"]+=event_weight;
    selected_events_percut["initial"]+=event_weight;
    
    
    //
    // Optical
    //
    //#1 number of spills greater than 0
    if (t->nbeamfls == 0) continue;
    

    int flashInBeamSpill = -1;
    double old_pe = -1;
    
    for (int fls = 0; fls < t->nbeamfls; fls ++){

      h_flsTime->Fill(t->beamfls_time.at(fls) - _flashShift, event_weight);
      if(t->beamfls_pe.at(fls) > _pe_cut) {
        h_flsTime_wcut->Fill(t->beamfls_time.at(fls) - _flashShift, event_weight);
      }
      if (t->beamfls_time.at(fls) > _beamSpillStarts && t->beamfls_time.at(fls) < _beamSpillEnds) {
        
        //flashInBeamSpill = fls;
        if (t->beamfls_pe.at(fls) >= _pe_cut) {
          if (t->beamfls_pe.at(fls) > old_pe) {
            flashInBeamSpill = fls;
            old_pe = t->beamfls_pe.at(fls);
          }
          nEvtsWFlashInBeamSpill++;
        }
      }
    }
    
    
    if (flashInBeamSpill == -1) continue;
    h_flsPe_wcut->Fill(t->beamfls_pe.at(flashInBeamSpill), event_weight);
    
    h_flsTime_flsPe_wcut->Fill(t->beamfls_time.at(flashInBeamSpill) - _flashShift, t->beamfls_pe.at(flashInBeamSpill), event_weight);

    h_flsTime_wcut_2->Fill(t->beamfls_time.at(flashInBeamSpill) - _flashShift, event_weight);
    

    
    
    if (isSignal) h_true_nu_eng_afterflash->Fill(t->nu_e, event_weight);

    if (isSignal) selected_signal_events_percut["beamflash"]+=event_weight;
    selected_events_percut["beamflash"]+=event_weight;
    
    
    if (t->nslices > 0) h_trklen_first->Fill(t->slc_longesttrack_length.at(0));
    if (t->nslices > 1) h_trklen_second->Fill(t->slc_longesttrack_length.at(1));

    //
    // Loop over TPC Events - Preliminary plots with only the flash cut applied
    //

    for (int slc = 0; slc < t->nslices; slc ++) {
      
      bool nu_origin = (t->slc_origin.at(slc) == 0 || t->slc_origin.at(slc) == 2);
      
      // PMTs
      if (flashInBeamSpill > -1 && t->slc_flsmatch_score.at(slc) > -1
          && t->slc_flsmatch_qllx.at(slc)!= -9999 && t->slc_flsmatch_tpcx.at(slc)!=-9999) {
        hmap_xdiff_b["total"]->Fill(t->slc_flsmatch_qllx.at(slc) - t->slc_flsmatch_tpcx.at(slc), event_weight);
        hmap_zdiff_b["total"]->Fill(t->slc_flsmatch_hypoz.at(slc) - t->beamfls_z.at(flashInBeamSpill), event_weight);
        if ( isSignal && nu_origin) {
          hmap_xdiff_b["signal"]->Fill(t->slc_flsmatch_qllx.at(slc) - t->slc_flsmatch_tpcx.at(slc), event_weight);
          hmap_zdiff_b["signal"]->Fill(t->slc_flsmatch_hypoz.at(slc) - t->beamfls_z.at(flashInBeamSpill), event_weight);
        } else {
          hmap_xdiff_b["background"]->Fill(t->slc_flsmatch_qllx.at(slc) - t->slc_flsmatch_tpcx.at(slc), event_weight);
          hmap_zdiff_b["background"]->Fill(t->slc_flsmatch_hypoz.at(slc) - t->beamfls_z.at(flashInBeamSpill), event_weight);
        }
      }
      
      //   nu
      if ( isSignal && nu_origin && flashInBeamSpill > -1 && t->slc_flsmatch_score.at(slc) > -1){
        for (int pmt = 0; pmt < 32; pmt++) {
          if ((t->slc_flshypo_spec.at(slc))[pmt] < 5 || (t->beamfls_spec.at(flashInBeamSpill))[pmt] < 5) continue;
          double mean = ((t->slc_flshypo_spec.at(slc))[pmt] + (t->beamfls_spec.at(flashInBeamSpill))[pmt]) / 2.;
          h_frac_diff->Fill(pmt, ( (t->slc_flshypo_spec.at(slc))[pmt] - (t->beamfls_spec.at(flashInBeamSpill))[pmt] ) / (mean) , event_weight);
        }
        if (t->slc_flsmatch_qllx.at(slc)!= -9999 && t->slc_flsmatch_tpcx.at(slc)!=-9999){
          h_xdiff->Fill(t->slc_flsmatch_qllx.at(slc) - t->slc_flsmatch_tpcx.at(slc), event_weight);
          h_zdiff->Fill(t->slc_flsmatch_hypoz.at(slc) - t->beamfls_z.at(flashInBeamSpill), event_weight);
        }
      }
      //   others
      if (t->slc_origin.at(slc) == 1 && flashInBeamSpill > -1 && t->slc_flsmatch_score.at(slc) > -1){
        for (int pmt = 0; pmt < 32; pmt++) {
          if ((t->slc_flshypo_spec.at(slc))[pmt] < 5 || (t->beamfls_spec.at(flashInBeamSpill))[pmt] < 5) continue;
          double mean = ((t->slc_flshypo_spec.at(slc))[pmt] + (t->beamfls_spec.at(flashInBeamSpill))[pmt]) / 2.;
          h_frac_diff_others->Fill(pmt, ( (t->slc_flshypo_spec.at(slc))[pmt] - (t->beamfls_spec.at(flashInBeamSpill))[pmt] ) / (mean) , event_weight);
        }
        if (t->slc_flsmatch_qllx.at(slc)!= -9999 && t->slc_flsmatch_tpcx.at(slc)!=-9999){
          h_xdiff_others->Fill(t->slc_flsmatch_qllx.at(slc) - t->slc_flsmatch_tpcx.at(slc), event_weight);
          h_zdiff_others->Fill(t->slc_flsmatch_hypoz.at(slc) - t->beamfls_z.at(flashInBeamSpill), event_weight);
        }
      }
      //  spec
      if (/*t->event==900 && t->run==5326*/t->event==-1 /*150801 4990051 3099969*/) {
        if (flashInBeamSpill > -1 && t->slc_flsmatch_score.at(slc) > -1){
          int instance = 0;
          std::cout << "tpcx " << t->slc_flsmatch_tpcx.at(instance) << std::endl;
          std::cout << "qllx " << t->slc_flsmatch_qllx.at(instance) << std::endl;
          std::cout << "SCORE IS " << t->slc_flsmatch_score.at(instance) << std::endl;
          for (int pmt = 0; pmt < 32; pmt++) {
            hypo_spec_x[pmt] = pmt;
            hypo_spec_y[pmt] = (t->slc_flshypo_spec.at(instance))[pmt]; ;//(t->slc_flshypo_spec.at(3))[pmt];
            meas_spec_x[pmt] = pmt;
            meas_spec_y[pmt] = (t->beamfls_spec.at(flashInBeamSpill))[pmt];
            numc_spec_x[pmt] = pmt;
            //numc_spec_y[pmt] = t->numc_flash_spec.at(pmt);
          }
        }
      }
      
      // CheckVertex
      hmap_vtxcheck_angle["total"]->Fill(t->slc_vtxcheck_angle.at(slc), event_weight);
      if (isSignal && nu_origin) {
        h_vtxcheck_angle_good->Fill(t->slc_vtxcheck_angle.at(slc), event_weight);
        hmap_vtxcheck_angle["signal"]->Fill(t->slc_vtxcheck_angle.at(slc), event_weight);
      } else {
        h_vtxcheck_angle_bad->Fill(t->slc_vtxcheck_angle.at(slc), event_weight);
        hmap_vtxcheck_angle["background"]->Fill(t->slc_vtxcheck_angle.at(slc), event_weight);
      }

      // Residuals std, mean
      hmap_residuals_std["total"]->Fill(t->slc_muoncandidate_residuals_std.at(slc), event_weight);
      hmap_residuals_mean["total"]->Fill(t->slc_muoncandidate_residuals_mean.at(slc), event_weight);
      hmap_perc_used_hits["total"]->Fill(t->slc_muoncandidate_perc_used_hits_in_cluster.at(slc), event_weight);
      if(t->slc_muoncandidate_contained.at(slc)) {
        hmap_mom_mcs_length["total"]->Fill(t->slc_muoncandidate_mom_mcs.at(slc) - t->slc_muoncandidate_mom_range.at(slc), event_weight);
      }
      if (isSignal && nu_origin) {
        hmap_residuals_std["signal"]->Fill(t->slc_muoncandidate_residuals_std.at(slc), event_weight);
        hmap_residuals_mean["signal"]->Fill(t->slc_muoncandidate_residuals_mean.at(slc), event_weight);
        hmap_perc_used_hits["signal"]->Fill(t->slc_muoncandidate_perc_used_hits_in_cluster.at(slc), event_weight);
        if(t->slc_muoncandidate_contained.at(slc)) {
          hmap_mom_mcs_length["signal"]->Fill(t->slc_muoncandidate_mom_mcs.at(slc) - t->slc_muoncandidate_mom_range.at(slc), event_weight);
        }
      } else {
        hmap_residuals_std["background"]->Fill(t->slc_muoncandidate_residuals_std.at(slc), event_weight);
        hmap_residuals_mean["background"]->Fill(t->slc_muoncandidate_residuals_mean.at(slc), event_weight);
        hmap_perc_used_hits["background"]->Fill(t->slc_muoncandidate_perc_used_hits_in_cluster.at(slc), event_weight);
        if(t->slc_muoncandidate_contained.at(slc))  {
          hmap_mom_mcs_length["background"]->Fill(t->slc_muoncandidate_mom_mcs.at(slc) - t->slc_muoncandidate_mom_range.at(slc), event_weight);
        }
      }
      
      // Track chi2
      h_chi2->Fill(t->slc_kalman_chi2.at(slc)/(double)t->slc_kalman_ndof.at(slc), event_weight);
      
      // Number of TPCObjects per event
      h_nslices->Fill(t->nslices, event_weight);
      
      
      
      
    } // slice loop
    

      
    int n_slc_flsmatch = 0;
    
    //
    // Find slice with maximum score
    //

    double score_max = -1;
    int scl_ll_max = -1;
    std::vector<double> temp_score; temp_score.clear();
    for (int slc = 0; slc < t->nslices; slc ++){

      temp_score.emplace_back(t->slc_flsmatch_score.at(slc));
      
      if (t->slc_flsmatch_score.at(slc) > 0.00000001) {
        n_slc_flsmatch++;
      }
      
      if (t->slc_flsmatch_score.at(slc) > score_max){
        scl_ll_max = slc;
        score_max = t->slc_flsmatch_score.at(slc);
      }
    }
    
    h_n_slc_flsmatch->Fill(n_slc_flsmatch, event_weight);

    std::sort(temp_score.begin(), temp_score.end(), std::greater<double>());


    //*******************************************
    //*******************************************
    //*******************************************
    // if (!t->is_selected) continue;
    //*******************************************
    //*******************************************
    //*******************************************
    //# 3 flash matched cut
    // In no flash-matched object, continue
    if (scl_ll_max == -1) continue;

    
    bool nu_origin = (t->slc_origin.at(scl_ll_max) == 0 || t->slc_origin.at(scl_ll_max) == 2);
    
    double dqdx_calib = t->slc_muoncandidate_dqdx_trunc.at(scl_ll_max) * _gainCalib;

    
    
    
    
    
    h_flsTime_wcut_3->Fill(t->beamfls_time.at(flashInBeamSpill) - _flashShift, event_weight);
    
    if (isSignal && nu_origin) nSignalFlashMatched ++;
    //# 4 flash matched cut
    // A score < 3e-4 or inf means no flash-matched object, continue
    if (score_max <= 3e-4) continue;
    if (std::isinf(score_max)) continue;


    if (isSignal && nu_origin) selected_signal_events_percut["flash_match"]+=event_weight;
    selected_events_percut["flash_match"]+=event_weight;
    
    
    h_flsTime_wcut_4->Fill(t->beamfls_time.at(flashInBeamSpill) - _flashShift, event_weight);
    h_deltax->Fill(t->slc_flsmatch_qllx.at(scl_ll_max) - t->slc_flsmatch_tpcx.at(scl_ll_max), event_weight);
    h_deltax_2d->Fill(t->slc_flsmatch_qllx.at(scl_ll_max), t->slc_flsmatch_tpcx.at(scl_ll_max), event_weight);
    h_deltaz_4->Fill(t->slc_flsmatch_hypoz.at(scl_ll_max) - t->beamfls_z.at(flashInBeamSpill), event_weight);
    
    // flash match deltax cut
    // If it doens't pass the flash-match deltaX cut, continue
    if(t->slc_flsmatch_qllx.at(scl_ll_max) - t->slc_flsmatch_tpcx.at(scl_ll_max) > 50) continue;
    
    h_flsTime_wcut_5->Fill(t->beamfls_time.at(flashInBeamSpill) - _flashShift, event_weight);

    if(t->slc_flsmatch_qllx.at(scl_ll_max) - t->slc_flsmatch_tpcx.at(scl_ll_max) < -100) continue;

    if (isSignal && nu_origin) selected_signal_events_percut["flash_match_deltax"]+=event_weight;
    selected_events_percut["flash_match_deltax"]+=event_weight;
    
    h_flsTime_wcut_6->Fill(t->beamfls_time.at(flashInBeamSpill) - _flashShift, event_weight);
    h_deltaz_6->Fill(t->slc_flsmatch_hypoz.at(scl_ll_max) - t->beamfls_z.at(flashInBeamSpill), event_weight);
    //flash match delta z cut
    // If it doens't pass the flash-match deltaZ cut, continue
    if(t->slc_flsmatch_hypoz.at(scl_ll_max) - t->beamfls_z.at(flashInBeamSpill) > 75) continue;
    
    h_flsTime_wcut_7->Fill(t->beamfls_time.at(flashInBeamSpill) - _flashShift, event_weight);

    if(t->slc_flsmatch_hypoz.at(scl_ll_max) - t->beamfls_z.at(flashInBeamSpill) < -75) continue;

    if (isSignal && nu_origin) selected_signal_events_percut["flash_match_deltaz"]+=event_weight;
    selected_events_percut["flash_match_deltaz"]+=event_weight;
    
    h_flsTime_wcut_8->Fill(t->beamfls_time.at(flashInBeamSpill) - _flashShift, event_weight);

    
      

    // if(t->slc_vtxcheck_angle.at(scl_ll_max) > 2.9) continue;
    
    //if(t->slc_vtxcheck_angle.at(scl_ll_max) < 0.05 && t->slc_vtxcheck_angle.at(scl_ll_max) !=-9999 ) continue;
    //#  there is at least one tracks in tpcobject   
    // If zero tracks in this tpcobject, continue
    if(t->slc_ntrack.at(scl_ll_max) == 0) continue;
    //# residual ans fraction of used hits in cluster
    // Cut on residuala ans fraction of used hits in cluster
    if (t->slc_muoncandidate_residuals_std.at(scl_ll_max) > 2.5) continue;
    // if (std::abs(t->slc_muoncandidate_residuals_mean.at(scl_ll_max)) > 0.7) continue;
    if (t->slc_muoncandidate_perc_used_hits_in_cluster.at(scl_ll_max) < 0.7) continue;
    
    //if(!t->slc_passed_min_track_quality.at(scl_ll_max)) continue;
    
    //if(!t->slc_passed_min_vertex_quality.at(scl_ll_max)) continue;

    if (isSignal && nu_origin) selected_signal_events_percut["quality"]+=event_weight;
    selected_events_percut["quality"]+=event_weight;
    
    //if (!t->slc_muoncandidate_contained.at(scl_ll_max)) continue;
    
    if(t->slc_muoncandidate_contained.at(scl_ll_max) && (t->slc_muoncandidate_mom_mcs.at(scl_ll_max) - t->slc_muoncandidate_mom_range.at(scl_ll_max) > 0.2)) continue;
    
    if (isSignal && nu_origin) selected_signal_events_percut["mcs_length_quality"]+=event_weight;
    selected_events_percut["mcs_length_quality"]+=event_weight;


    // DqDx cut
    std::vector<double> svm_x = {86300, 86050, 85850, 85600, 85400, 85150, 84950, 84700, 84500, 84300, 84100, 83850, 83650, 83450, 83250, 83050, 82850, 82650, 82450, 82250, 82050, 81900, 81700, 81500, 81300, 81150, 80950, 80750, 80600, 80400, 80250, 80050, 79900, 79750, 79550, 79400, 79250, 79050, 78900, 78750, 78600, 78400, 78250, 78100, 77950, 77800, 77650, 77500, 77350, 77200, 77050, 76900, 76750, 76650, 76500, 76350, 76200, 76050, 75950, 75800, 75650, 75550, 75400, 75300, 75150, 75000, 74900, 74750, 74650, 74500, 74400, 74250, 74150, 74050, 73900, 73800, 73700, 73550, 73450, 73350, 73250, 73100, 73000, 72900, 72800, 72700, 72550, 72450, 72350, 72250, 72150, 72050, 71950, 71850, 71750, 71650, 71550, 71450, 71350, 71250, 71150, 71100, 71000, 70900, 70800, 70700, 70600, 70550, 70450, 70350, 70250, 70200, 70100, 70000, 69950, 69850, 69750, 69700, 69600, 69550, 69450, 69350, 69300, 69200, 69150, 69050, 69000, 68900, 68850, 68750, 68700, 68600, 68550, 68500, 68400, 68350, 68250, 68200, 68150, 68050, 68000, 67950, 67850, 67800, 67750, 67700, 67600, 67550, 67500, 67450, 67350, 67300, 67250, 67200, 67150, 67100, 67000, 66950, 66900, 66850, 66800, 66750, 66700, 66650, 66600, 66550, 66500, 66450, 66400, 66350, 66300, 66250, 66200, 66150, 66100, 66050, 66000, 65950, 65900, 65850, 65800, 65750, 65750, 65700, 65650, 65600, 65550, 65500, 65450, 65450, 65400, 65350, 65300, 65250, 65250, 65200, 65150, 65100, 65100, 65050, 65000, 65000, 64950, 64900, 64850, 64850, 64800, 64750, 64750, 64700, 64700, 64650, 64600, 64600, 64550, 64500, 64500, 64450, 64450, 64400, 64400, 64350, 64350, 64300, 64250, 64250, 64200, 64200, 64150, 64150, 64100, 64100, 64100, 64050, 64050, 64000, 64000, 63950, 63950, 63900, 63900, 63900, 63850, 63850, 63800, 63800, 63800, 63750, 63750, 63750, 63700, 63700, 63700, 63650, 63650, 63650, 63600, 63600, 63600, 63550, 63550, 63550, 63500, 63500, 63500, 63500, 63450, 63450, 63450, 63450, 63400, 63400, 63400, 63400, 63400, 63350, 63350, 63350, 63350, 63350, 63350, 63300, 63300, 63300, 63300, 63300, 63300, 63250, 63250, 63250, 63250, 63250, 63250, 63250, 63250, 63250, 63250, 63200, 63200, 63200, 63200, 63200, 63200, 63200, 63200, 63200, 63200, 63200, 63200, 63200, 63200, 63200, 63200, 63200, 63200, 63200, 63200, 63200, 63200, 63200, 63200, 63200, 63200, 63200, 63200, 63200, 63200, 63200, 63200, 63200, 63250, 63250, 63250, 63250, 63250, 63250, 63250, 63250, 63250, 63250, 63300, 63300, 63300, 63300, 63300, 63300, 63300, 63350, 63350, 63350, 63350, 63350, 63350, 63400, 63400, 63400, 63400, 63400, 63400, 63450, 63450, 63450, 63450, 63500, 63500, 63500, 63500, 63500, 63550, 63550, 63550, 63550, 63600, 63600, 63600, 63600, 63650, 63650, 63650, 63650, 63700, 63700, 63700, 63750, 63750, 63750, 63800, 63800, 63800, 63800, 63850, 63850, 63850, 63900, 63900, 63900, 63950, 63950, 63950, 64000, 64000, 64000, 64050, 64050, 64100, 64100, 64100, 64150, 64150, 64150, 64200, 64200, 64250, 64250, 64250, 64300, 64300, 64350, 64350, 64350, 64400, 64400, 64450, 64450, 64450, 64500, 64500, 64550, 64550, 64600, 64600, 64650, 64650, 64650, 64700, 64700, 64750, 64750, 64800, 64800, 64850, 64850, 64900, 64900, 64950, 64950, 65000, 65000, 65050, 65050, 65100, 65100, 65150, 65150, 65200, 65200, 65250, 65250, 65300, 65300, 65350, 65350, 65400, 65400, 65450, 65500, 65500, 65550, 65550, 65600, 65600, 65650, 65650, 65700, 65750, 65750, 65800, 65800, 65850, 65850, 65900, 65950, 65950, 66000, 66000, 66050, 66100, 66100, 66150, 66150, 66200, 66250, 66250, 66300, 66350, 66350, 66400, 66400, 66450, 66500, 66500, 66550, 66600, 66600, 66650, 66700, 66700, 66750, 66750, 66800, 66850, 66850, 66900, 66950, 66950, 67000, 67050, 67050, 67100, 67150, 67150, 67200, 67250, 67300, 67300, 67350, 67400, 67400, 67450, 67500, 67500, 67550, 67600, 67650, 67650, 67700, 67750, 67750, 67800, 67850, 67900, 67900, 67950, 68000, 68050, 68050, 68100, 68150, 68150, 68200, 68250, 68300, 68300, 68350, 68400, 68450, 68450, 68500, 68550, 68600, 68650, 68650, 68700, 68750, 68800, 68800, 68850, 68900, 68950, 69000, 69000, 69050, 69100, 69150, 69150, 69200, 69250, 69300, 69350, 69350, 69400, 69450, 69500, 69550, 69600, 69600, 69650, 69700, 69750, 69800, 69800, 69850, 69900, 69950, 70000, 70050, 70050, 70100, 70150, 70200, 70250, 70300, 70300, 70350, 70400, 70450, 70500, 70550, 70600, 70600, 70650, 70700, 70750, 70800, 70850, 70900, 70900, 70950, 71000, 71050, 71100, 71150, 71200, 71250, 71300, 71300, 71350, 71400, 71450, 71500, 71550, 71600, 71650, 71700, 71700, 71750, 71800, 71850, 71900, 71950, 72000, 72050, 72100, 72150, 72200, 72200, 72250, 72300, 72350, 72400, 72450, 72500, 72550, 72600, 72650, 72700, 72750, 72800, 72850, 72850, 72900, 72950, 73000, 73050, 73100, 73150, 73200, 73250, 73300, 73350, 73400, 73450, 73500, 73550, 73600, 73650, 73700, 73750, 73800, 73850, 73900, 73900, 73950, 74000, 74050, 74100, 74150, 74200, 74250, 74300, 74350, 74400, 74450, 74500, 74550, 74600, 74650, 74700, 74750, 74800, 74850, 74900, 74950, 75000, 75050, 75100, 75150, 75200, 75250, 75300, 75350, 75400, 75450, 75500, 75550, 75600, 75650, 75700, 75750, 75800, 75850, 75900, 76000, 76050, 76100, 76150, 76200, 76250, 76300, 76350, 76400, 76450, 76500, 76550, 76600, 76650, 76700, 76750, 76800, 76850, 76900, 76950, 77000, 77050, 77100, 77200, 77250, 77300, 77350, 77400, 77450, 77500, 77550, 77600, 77650, 77700, 77750, 77800, 77850, 77900, 78000, 78050, 78100, 78150, 78200, 78250, 78300, 78350, 78400, 78450, 78500, 78550, 78600, 78700, 78750, 78800, 78850, 78900, 78950, 79000, 79050, 79100, 79150, 79250, 79300, 79350, 79400, 79450, 79500, 79550, 79600, 79650, 79750, 79800, 79850, 79900, 79950, 80000, 80050, 80100, 80150, 80250, 80300, 80350, 80400, 80450, 80500, 80550, 80600, 80700, 80750, 80800, 80850, 80900, 80950, 81000, 81100, 81150, 81200, 81250, 81300, 81350, 81400, 81450, 81550, 81600, 81650, 81700, 81750, 81800, 81900, 81950, 82000, 82050, 82100, 82150, 82200, 82300, 82350, 82400, 82450, 82500, 82550, 82650, 82700, 82750, 82800, 82850, 82900, 83000, 83050, 83100, 83150, 83200, 83250, 83350, 83400, 83450, 83500, 83550, 83650, 83700, 83750, 83800, 83850, 83900, 84000, 84050, 84100, 84150, 84200, 84300, 84350, 84400, 84450, 84500, 84600, 84650, 84700, 84750, 84800, 84900, 84950, 85000, 85050, 85100, 85200, 85250, 85300, 85350, 85400, 85500, 85550, 85600, 85650, 85700, 85800, 85850, 85900, 85950, 86000, 86100, 86150, 86200, 86250, 86350, 86400, 86450, 86500, 86550, 86650, 86700, 86750, 86800, 86900, 86950, 87000, 87050, 87150, 87200, 87250, 87300, 87350, 87450, 87500, 87550, 87600, 87700, 87750, 87800, 87850, 87950, 88000, 88050, 88100, 88200, 88250, 88300, 88350, 88450, 88500, 88550, 88600, 88700, 88750, 88800, 88850, 88950, 89000, 89050, 89100, 89200, 89250, 89300, 89350, 89450, 89500, 89550, 89600, 89700, 89750, 89800, 89900, 89950, 90000, 90050, 90150, 90200, 90250, 90300, 90400, 90450, 90500, 90600, 90650, 90700, 90750, 90850, 90900, 90950, 91000, 91100, 91150, 91200, 91300, 91350, 91400, 91450, 91550, 91600, 91650, 91750, 91800, 91850, 91950, 92000, 92050, 92100, 92200, 92250};

    double l = std::round(t->slc_muoncandidate_length.at(scl_ll_max));
    double dqdx_cut = 200000;
    if (l >= 0 && l < 1000) {
      dqdx_cut = svm_x.at(l);
    }
    //# truncated mean dqdx cut      
    if (dqdx_calib > dqdx_cut) continue;
    // if(t->slc_nuvtx_z.at(scl_ll_max) <= 500) continue;

    if (isSignal && nu_origin) selected_signal_events_percut["mip_consistency"]+=event_weight;
    selected_events_percut["mip_consistency"]+=event_weight;



    // Just before the FV cut, make distribution of vtxz
    if (t->slc_nuvtx_y.at(scl_ll_max) > 82) {
      hmap_vtxz_upborder["total"]->Fill(t->slc_nuvtx_z.at(scl_ll_max));
      hmap_vtxx_upborder["total"]->Fill(t->slc_nuvtx_z.at(scl_ll_max));
      if (isSignal && nu_origin) {
        hmap_vtxz_upborder["signal"]->Fill(t->slc_nuvtx_z.at(scl_ll_max));
        hmap_vtxx_upborder["signal"]->Fill(t->slc_nuvtx_z.at(scl_ll_max));
      } else {
        hmap_vtxz_upborder["background"]->Fill(t->slc_nuvtx_z.at(scl_ll_max));
        hmap_vtxx_upborder["background"]->Fill(t->slc_nuvtx_z.at(scl_ll_max));
      }
    }


    // FV cut
    if(t->slc_nuvtx_fv.at(scl_ll_max) == 0) continue;
    if(t->slc_nuvtx_z.at(scl_ll_max) > 675 && t->slc_nuvtx_z.at(scl_ll_max) < 775) continue;
    // if(t->slc_nuvtx_z.at(scl_ll_max) < 200) continue;

    if (isSignal && nu_origin) selected_signal_events_percut["fiducial_volume"]+=event_weight;
    selected_events_percut["fiducial_volume"]+=event_weight;

    // Select the bump only
    // if(t->slc_iscontained.at(scl_ll_max)) continue; // Uncontained
    // if(t->slc_ntrack.at(scl_ll_max) != 1) continue; // Multiplicity == 1
    // if(t->slc_nuvtx_y.at(scl_ll_max) > 0) continue; // Vertex in the bottom half of the detector
    // if(t->slc_nuvtx_x.at(scl_ll_max) > 128.175) continue; // Vertex in the anode half of the detector
    // if(t->slc_longesttrack_theta.at(scl_ll_max) > -0.6) continue; // cos(theta) < -0.6


    // if(t->slc_nuvtx_x.at(scl_ll_max) > 128.175) continue; // anode_vtx
    // if(t->slc_nuvtx_x.at(scl_ll_max) <= 128.175) continue; // cathode_vtx
    
    // if (abs(t->slc_longesttrack_phi.at(scl_ll_max)) < TMath::Pi()/2.) continue; // right (towards the anode)
    // if (abs(t->slc_longesttrack_phi.at(scl_ll_max)) >= TMath::Pi()/2.) continue; // left (towards the cathode)

    // if (!(t->slc_longesttrack_phi.at(scl_ll_max) > -TMath::Pi()/4 && t->slc_longesttrack_phi.at(scl_ll_max) < TMath::Pi()/4)) continue; // cathode
    // if (!(t->slc_longesttrack_phi.at(scl_ll_max) > TMath::Pi()/4 && t->slc_longesttrack_phi.at(scl_ll_max) < (3./4.)*TMath::Pi())) continue; // up
    // if (!(t->slc_longesttrack_phi.at(scl_ll_max) > -(3./4.)*TMath::Pi() && t->slc_longesttrack_phi.at(scl_ll_max) < -TMath::Pi()/4)) continue; // down
    // if (!((t->slc_longesttrack_phi.at(scl_ll_max) > (3./4.)*TMath::Pi() && t->slc_longesttrack_phi.at(scl_ll_max) < TMath::Pi())
    //         || (t->slc_longesttrack_phi.at(scl_ll_max) > -TMath::Pi() && t->slc_longesttrack_phi.at(scl_ll_max) < -(3./4.)*TMath::Pi()))) continue; // anode


    //if(t->slc_mult_track_tolerance.at(scl_ll_max) <= 1) continue;


    //if (t->slc_nuvtx_z.at(scl_ll_max) < 300) continue;

    //if(t->slc_longesttrack_length.at(scl_ll_max) < 25.) continue;
    
    // if(t->slc_iscontained.at(scl_ll_max)) continue;
    
    //if(t->slc_crosses_top_boundary.at(scl_ll_max) == 1) continue;

    // if (t->slc_muoncandidate_mom_mcs.at(scl_ll_max) > 2.5) continue;
    //===========================================================================================================

    // Remove flipped tracks (truth cut)
    // if (isSignal && (t->lep_costheta * t->slc_longesttrack_theta.at(scl_ll_max) < 0) ) continue;


    //std::cout<<"End of CCinclusive Selection and Start the CC1uNP Selection"<<std::endl;


    // CC1uNP Selection
    bool trackfromneutrino=true;
    if (t->pfp_truth_origin.size() != t->pfp_reco_ismuoncandidate.size()){
      std::cout << "Truth and Reco vectors not the same length!" << std::endl;
    }
    if (t->pfp_truth_origin.size() != t->pfp_reco_istrack.size()){
      std::cout << "Truth and Reco vectors not the same length!" << std::endl;
    }
     
    for(size_t npfp=0; npfp<t->pfp_truth_origin.size(); npfp++){
        if(t->pfp_reco_istrack[npfp]==0 &&t->pfp_reco_isshower[npfp]==0 ) continue;
        if(_showerastrack){
              if(t->pfp_truth_origin[npfp] !=1) {trackfromneutrino=false;}
        }
        else if(t->pfp_reco_istrack[npfp] && t->pfp_truth_origin[npfp]!=1) {trackfromneutrino=false; }
    }

    if(_ana_int_type=="ccinclusive_analysis"){trackfromneutrino=true;}
    muon_contained = false; 
    // AF - Move this upwards, before any cuts get applied
    for(unsigned int jj=0; jj<t->pfp_reco_ismuoncandidate.size(); jj++){
       if(t->pfp_reco_ismuoncandidate[jj]==0) continue;
       if(inCV(t->pfp_reco_endx[jj], t->pfp_reco_endy[jj], t->pfp_reco_endz[jj])) {muon_contained = true;}
       muind=jj;  
    }
   
    
    
    //#1 number of tracks>=2
    if(_showerastrack){
       if(_ana_int_type=="cc1unp_analysis" && t->num_pfp<2)  continue;
    } else {
       if(_ana_int_type=="cc1unp_analysis" && t->num_pfp_tracks<2)  continue;
    }
    if (isSignal && nu_origin && trackfromneutrino) {selected_signal_events_percut["ntrk2"] +=event_weight;}
    selected_events_percut["ntrk2"] +=event_weight;
    
    
    
    //#2 all the protons contained in FV
    int uncontained_proton=0;
    
    // loop over all the proton candidates and select the events with all the proton candidates 
    for(size_t n_trk_pfp=0; n_trk_pfp<t->pfp_reco_startx.size(); n_trk_pfp++){
       if(t->pfp_reco_ismuoncandidate[n_trk_pfp] == 1) continue;
       //if(t->pfp_reco_istrack[n_trk_pfp] !=1 && t->pfp_reco_isshower[n_trk_pfp] !=1 ) continue;
       if(_showerastrack){
         if(t->pfp_reco_numtracks[n_trk_pfp] !=1) continue;
       } else{
         if(t->pfp_reco_istrack[n_trk_pfp] !=1) continue;  
       }
       if(!inCV(t->pfp_reco_startx[n_trk_pfp], t->pfp_reco_starty[n_trk_pfp], t->pfp_reco_startz[n_trk_pfp]) ||
           !inCV(t->pfp_reco_endx[n_trk_pfp],   t->pfp_reco_endy[n_trk_pfp],   t->pfp_reco_endz[n_trk_pfp]))   {
           uncontained_proton +=1;
       }
       
    }
    
   
    if(_ana_int_type=="cc1unp_analysis" && uncontained_proton>0)  continue;
    /*
    bool pcontained_flag=true;

    for(size_t n_trk_pfp=0; n_trk_pfp<t->pfp_reco_startx.size(); n_trk_pfp++){
       if(t->pfp_reco_ismuoncandidate[n_trk_pfp] == 1) continue;
       //if(t->pfp_reco_istrack[n_trk_pfp] ==0 && t->pfp_reco_numtracks[n_trk_pfp] ==0) continue;
       //if(t->pfp_reco_istrack[n_trk_pfp] == 0) continue;
       if(t->pfp_reco_numtracks[n_trk_pfp] ==0) continue; 
        if(!inCV(t->pfp_reco_startx[n_trk_pfp], t->pfp_reco_starty[n_trk_pfp], t->pfp_reco_startz[n_trk_pfp]) ||
           !inCV(t->pfp_reco_endx[n_trk_pfp],   t->pfp_reco_endy[n_trk_pfp],   t->pfp_reco_endz[n_trk_pfp]))  pcontained_flag=false;
        
       
    }
 
    if(_ana_int_type=="cc1unp_analysis" && pcontained_flag==0) continue;
    */
    if (isSignal && nu_origin && trackfromneutrino) selected_signal_events_percut["pinCV"] +=event_weight;
    selected_events_percut["pinCV"] +=event_weight;
     
    //#3 find the leading proton candidate and perform the minCol cut
    pind=-999;
    float temp_length=-999.0;
    
    for(size_t np=0; np<t->pfp_reco_length.size(); np++){
        if(t->pfp_reco_ismuoncandidate[np]==1)  continue;
        
        if(_showerastrack){
           if(t->pfp_reco_istrack[np]==0 && t->pfp_reco_isshower[np] ==0) continue;      
           if(t->pfp_reco_numtracks[np] !=1) continue;
        } else {
           if(t->pfp_reco_istrack[np]==0 && t->pfp_reco_isshower[np] ==0) continue;      
           if(t->pfp_reco_istrack[np]==0) continue; 
        }
        if(t->pfp_reco_length[np]> temp_length) {
           temp_length=t->pfp_reco_length[np];
           pind=np;

        }
    }
     
    if(_ana_int_type=="cc1unp_analysis" && pind<0) continue;        
    if(_ana_int_type=="cc1unp_analysis" && pind>-999 && t->pfp_reco_dEdx[pind].size()<5)  continue; 
    

    if (isSignal && nu_origin && trackfromneutrino) selected_signal_events_percut["minCol"] +=event_weight;
    selected_events_percut["minCol"] +=event_weight;

    for(size_t nchi2=0; nchi2<t->pfp_reco_ismuoncandidate.size(); nchi2++){
      if(t->pfp_reco_ismuoncandidate[nchi2] ==1) continue;
      chi2_proton_hypothesis_pcand->Fill(t->pfp_reco_chi2_proton[nchi2], event_weight);
    }

    if(!isdata && _ana_int_type =="cc1unp_analysis"){
    //#4 check the chi2 distribution and perform the chi2 cut
    for(size_t npfp=0; npfp<t->pfp_truth_pdg.size(); npfp++){
        
                if(t->pfp_reco_ismuoncandidate[npfp]==1) continue;
                if(abs(t->pfp_truth_pdg[npfp])==2212) {
                  chi2_proton_hypothesis["proton"]->Fill(t->pfp_reco_chi2_proton[npfp],event_weight); 
                  chi2_muon_hypothesis["proton"]->Fill(t->pfp_reco_chi2_muon[npfp],event_weight); 
                  chi2_pion_hypothesis["proton"]->Fill(t->pfp_reco_chi2_pion[npfp],event_weight); 
                  chi2_kaon_hypothesis["proton"]->Fill(t->pfp_reco_chi2_kaon[npfp],event_weight); 
                  //std::cout<<"Found a Proton !!"<<std::endl;
                }
                else if(abs(t->pfp_truth_pdg[npfp])==13  ) {
                  chi2_proton_hypothesis["muon"]->Fill(t->pfp_reco_chi2_proton[npfp],event_weight); 
                  chi2_muon_hypothesis["muon"]->Fill(t->pfp_reco_chi2_muon[npfp],event_weight); 
                  chi2_pion_hypothesis["muon"]->Fill(t->pfp_reco_chi2_pion[npfp],event_weight); 
                  chi2_kaon_hypothesis["muon"]->Fill(t->pfp_reco_chi2_kaon[npfp],event_weight); 
                  //std::cout<<"Found a Muon !!"<<std::endl;
                }
                else if(abs(t->pfp_truth_pdg[npfp])==211 || abs(t->pfp_truth_pdg[npfp])==111) {
                  chi2_proton_hypothesis["pion"]->Fill(t->pfp_reco_chi2_proton[npfp],event_weight); 
                  chi2_muon_hypothesis["pion"]->Fill(t->pfp_reco_chi2_muon[npfp],event_weight); 
                  chi2_pion_hypothesis["pion"]->Fill(t->pfp_reco_chi2_pion[npfp],event_weight); 
                  chi2_kaon_hypothesis["pion"]->Fill(t->pfp_reco_chi2_kaon[npfp],event_weight); 
                  //std::cout<<"Found a pion !!"<<std::endl;
                }
                else if(abs(t->pfp_truth_pdg[npfp])==321 || abs(t->pfp_truth_pdg[npfp])==311) {
                  chi2_proton_hypothesis["kaon"]->Fill(t->pfp_reco_chi2_proton[npfp],event_weight); 
                  chi2_muon_hypothesis["kaon"]->Fill(t->pfp_reco_chi2_muon[npfp],event_weight); 
                  chi2_pion_hypothesis["kaon"]->Fill(t->pfp_reco_chi2_pion[npfp],event_weight); 
                  chi2_kaon_hypothesis["kaon"]->Fill(t->pfp_reco_chi2_kaon[npfp],event_weight); 
                  //std::cout<<"Found a kaon !!"<<std::endl;
                } else {
                  chi2_proton_hypothesis["other"]->Fill(t->pfp_reco_chi2_proton[npfp],event_weight); 
                  chi2_muon_hypothesis["other"]->Fill(t->pfp_reco_chi2_muon[npfp],event_weight); 
                  chi2_pion_hypothesis["other"]->Fill(t->pfp_reco_chi2_pion[npfp],event_weight); 
                  chi2_kaon_hypothesis["other"]->Fill(t->pfp_reco_chi2_kaon[npfp],event_weight); 
                  //std::cout<<"Found a PF particle not muon, proton, pion or kaon !!"<<std::endl;
                }
    }//end if loop over all the proton candidate
    } //end of it is not data file
 

    //try to perform nhits cut to all the proton candidates

    //loop over the hits of all the tracks and apply the number of hits cut
    /*bool minColflag=true; 
    bool chi2flag=true;
    if(t->pfp_reco_chi2_proton[muind]<88.000) {chi2flag=false;}
    if(t->pfp_reco_dEdx[muind].size()<5) {minColflag=false;}

    for(unsigned int ind_trk=0; ind_trk< t->pfp_reco_dEdx.size(); ind_trk++){ 
     if(t->pfp_reco_ismuoncandidate[ind_trk]==1) continue;
     if(t->pfp_reco_istrack[ind_trk]==0) continue;
     if(t->pfp_reco_dEdx[ind_trk].size()<5) {minColflag=false;}
     if(t->pfp_reco_chi2_proton[ind_trk]>88.000) {chi2flag=false;}   
    }

    if(_ana_int_type=="cc1unp_analysis" && minColflag==false) continue;
    if(_ana_int_type=="cc1unp_analysis" && chi2flag==false) continue;
    */


     
    Int_t npcand_fail_chi2=0;
    for(size_t ntrk=0; ntrk<t->pfp_reco_chi2_proton.size(); ntrk++){
        if(t->pfp_reco_ismuoncandidate[ntrk]==1) continue;
        if(t->pfp_reco_istrack[ntrk]==0 && t->pfp_reco_isshower[ntrk]==0) continue;
        if(_showerastrack){
          if(t->pfp_reco_numtracks[ntrk] !=1) continue;
        } else {
          if(t->pfp_reco_istrack[ntrk]==0) continue;
        }
        if(t->pfp_reco_dEdx[ntrk].size()>=5 && t->pfp_reco_chi2_proton[ntrk]>88) {npcand_fail_chi2++;}
    }
    

    if(_ana_int_type=="cc1unp_analysis" && npcand_fail_chi2>0) continue;
    

    





    if (isSignal && nu_origin && trackfromneutrino) selected_signal_events_percut["chi2"] +=event_weight;
    selected_events_percut["chi2"] +=event_weight;

    if(_ana_int_type=="cc1unp_analysis" && t->pfp_reco_Mom_proton[pind]<0.3) continue;
    if(_ana_int_type=="cc1unp_analysis" && t->pfp_reco_Mom_proton[pind]>1.2) continue;
    if(_ana_int_type=="cc1unp_analysis" && t->pfp_reco_Mom_muon[muind]<0.1) continue;

    
    

    if(_ana_int_type=="cc1unp_analysis"){
      thetamup=getAngle(t->pfp_reco_Mom_MCS[muind], t->pfp_reco_theta[muind], t->pfp_reco_phi[muind], t->pfp_reco_Mom_proton[pind], t->pfp_reco_theta[pind], t->pfp_reco_phi[pind]);
      //loop over all the particles and calculate PTmis
      float px_total=0;
      float py_total=0;
      float pxl=0;
      float pyl=0;
      float Esum=0;
      
       
      float muonmass=0.1069;
      float protonmass=0.938;
      ptmis=0.; alphat=0.;  phit=0.; enucal=0.;
      if(muon_contained){
      pxl     =t->pfp_reco_Mom_muon[muind]*TMath::Sin(t->pfp_reco_theta[muind])*TMath::Cos(t->pfp_reco_phi[muind]);
      pyl     =t->pfp_reco_Mom_muon[muind]*TMath::Sin(t->pfp_reco_theta[muind])*TMath::Sin(t->pfp_reco_phi[muind]);
      } else {
      pxl     =t->pfp_reco_Mom_MCS[muind]*TMath::Sin(t->pfp_reco_theta[muind])*TMath::Cos(t->pfp_reco_phi[muind]);
      pyl     =t->pfp_reco_Mom_MCS[muind]*TMath::Sin(t->pfp_reco_theta[muind])*TMath::Sin(t->pfp_reco_phi[muind]);
      }
      if(muon_contained) {
           px_total=t->pfp_reco_Mom_muon[muind]*TMath::Sin(t->pfp_reco_theta[muind])*TMath::Cos(t->pfp_reco_phi[muind]);
           py_total=t->pfp_reco_Mom_muon[muind]*TMath::Sin(t->pfp_reco_theta[muind])*TMath::Sin(t->pfp_reco_phi[muind]); 
      } else{ 
      px_total=t->pfp_reco_Mom_MCS[muind]*TMath::Sin(t->pfp_reco_theta[muind])*TMath::Cos(t->pfp_reco_phi[muind]);
      py_total=t->pfp_reco_Mom_MCS[muind]*TMath::Sin(t->pfp_reco_theta[muind])*TMath::Sin(t->pfp_reco_phi[muind]);
      }
      Esum=sqrt(t->pfp_reco_Mom_MCS[muind]*t->pfp_reco_Mom_MCS[muind]+muonmass*muonmass);

      pion_reco=false;
      pion_reint=false;
      for(size_t ii=0; ii<t->pfp_reco_ismuoncandidate.size(); ii++){
        if(t->pfp_reco_ismuoncandidate[ii]==1) continue;  //proton candidate
        if(t->pfp_reco_istrack[ii]==0 &&t->pfp_reco_isshower[ii]==0) continue;   
        if(_showerastrack){
          if(t->pfp_reco_numtracks[ii] !=1) continue;
        }else{ 
            if(t->pfp_reco_istrack[ii]==0) continue;
        }   
 
         if(t->pfp_reco_dEdx[ii].size()<5){ h_nhits_lt5->Fill(t->pfp_reco_dEdx[ii].size(), event_weight); }
         Esum=Esum+sqrt(t->pfp_reco_Mom_proton[ii]*t->pfp_reco_Mom_proton[ii]+protonmass*protonmass)-protonmass; 
         px_total=px_total+t->pfp_reco_Mom_proton[ii]*TMath::Sin(t->pfp_reco_theta[ii])*TMath::Cos(t->pfp_reco_phi[ii]);
         py_total=py_total+t->pfp_reco_Mom_proton[ii]*TMath::Sin(t->pfp_reco_theta[ii])*TMath::Sin(t->pfp_reco_phi[ii]);
         if(abs(t->pfp_truth_pdg[ii])==211){
             pion_reco=true;
             h_pimom_total->Fill(t->pfp_truth_mom[ii], event_weight);
             if(t->pfp_truth_endProcess[ii]=="pi+Inelastic" ||t->pfp_truth_endProcess[ii]=="pi-Inelastic" ||t->pfp_truth_endProcess[ii]=="pi0Inelastic") {
                pion_reint=true;
                h_pimom_reint->Fill(t->pfp_truth_mom[ii], event_weight);  
             }
         }
         if(abs(t->pfp_truth_pdg[ii])==2212){
           h_pmom_total->Fill(t->pfp_truth_mom[ii], event_weight);
           if(t->pfp_truth_endProcess[ii]=="protonInelastic"){
              h_pmom_reint->Fill(t->pfp_truth_mom[ii], event_weight);
           } 
           if(t->pfp_reco_dEdx[ii].size()<5){
             h_nhits_lt5_proton->Fill(t->pfp_reco_dEdx[ii].size(), event_weight);
           }
         }
      }
      
      ptmis=TMath::Sqrt(px_total*px_total+py_total*py_total); 
      alphat=TMath::ACos((-pxl*px_total-pyl*py_total)/
                         (TMath::Sqrt(pxl*pxl+pyl*pyl)*TMath::Sqrt(px_total*px_total+py_total*py_total)));
      phit = TMath::ACos((-pxl*(px_total-pxl)-pyl*(py_total-pyl))/
                         (TMath::Sqrt(pxl*pxl+pyl*pyl)*TMath::Sqrt((px_total-pxl)*(px_total-pxl)+(py_total-pyl)*(py_total-pyl))));

      enucal=Ecalomiss(Esum, ptmis, t->pfp_reco_ismuoncandidate.size()-1);
      if(t->pfp_reco_length[muind]>5 && t->pfp_reco_length[pind]>5){
         etatest=getEta(t->pfp_reco_dQdx, t->pfp_reco_RR, t->pfp_reco_length,  muind, pind);
      }
      //=============================================================================================================
      



    }
    
    //std::cout<<"End of CC1uNP Selection"<<" is data or not "<<isdata<<std::endl;
    if(!isdata){
    //get the momentum of the leading proton momentum 
    if(isSignal && trackfromneutrino && _ana_int_type == "cc1unp_analysis"){
      //check if the isprimary variable filled correctly
      // loop over all the genie particles and get the denominator of the proton momentum efficiency
      temp_pmom=-999.0;
      temp_pangle=-999.0;
      temp_pphi=-999.0;
      temp_thetamup=-999.0;
      float temp_penergy=-999.0;
      for(size_t mpar=0; mpar<t->genie_mcpar_pdgcode.size(); mpar++){
           
           if(abs(t->genie_mcpar_pdgcode[mpar])==2212 && t->genie_mcpar_energy[mpar]>temp_penergy) {
               temp_penergy=t->genie_mcpar_energy[mpar];
               temp_pmom=TMath::Sqrt(t->genie_mcpar_px[mpar]*t->genie_mcpar_px[mpar]+
                              t->genie_mcpar_py[mpar]*t->genie_mcpar_py[mpar]+
                              t->genie_mcpar_pz[mpar]*t->genie_mcpar_pz[mpar]);
               temp_pangle=t->genie_mcpar_pz[mpar]/temp_pmom;
               temp_pphi=TMath::ACos(t->genie_mcpar_px[mpar]/TMath::Sqrt(t->genie_mcpar_px[mpar]*t->genie_mcpar_px[mpar]+t->genie_mcpar_py[mpar]*t->genie_mcpar_py[mpar]));
               temp_thetamup=getAngle(temp_pmom, TMath::ACos(temp_pangle), temp_pphi,  t->true_muon_mom_matched, TMath::ACos(t->lep_costheta), t->lep_phi); 
            }
      }     
      //std::cout<<"TRUE MOMENTUM OF PROTON IS "<<temp_pmom<<std::endl; 
      _event_histo_1d->h_eff_pmom_num->Fill(temp_pmom, event_weight);
      _event_histo_1d->h_eff_pangle_num->Fill(temp_pangle, event_weight);
      _event_histo_1d->h_eff_thetamup_num->Fill(temp_thetamup, event_weight);
      
      h_eff_allprotons_num->Fill(temp_pmom, event_weight);
      //loop over g4 stage and select the contained cc1unp and reinteracted or non reinteracted protons
      temp_geantpmom=-999.0;
      temp_geantpenergy=-999.0;
      bool geant_pcontained = true;
      for(size_t g4par=0; g4par<t->geant_mcpar_pdgcode.size(); g4par++){
         if(t->geant_mcpar_pdgcode[g4par] != 2212) continue;
         if(!inCV(t->geant_mcpar_startx[g4par], t->geant_mcpar_starty[g4par], t->geant_mcpar_startz[g4par]) || !inCV(t->geant_mcpar_endx[g4par], t->geant_mcpar_endy[g4par], t->geant_mcpar_endz[g4par])) {geant_pcontained=false;}
         //if(t->geant_mcpar_end_process[g4par] != "protonInelastic") continue;
         if(t->geant_mcpar_energy[g4par]>temp_geantpenergy) {
           temp_geantpenergy=t->geant_mcpar_energy[g4par];
           temp_geantpmom=TMath::Sqrt(t->geant_mcpar_px[g4par]*t->geant_mcpar_px[g4par]+
                                      t->geant_mcpar_py[g4par]*t->geant_mcpar_py[g4par]+
                                      t->geant_mcpar_pz[g4par]*t->geant_mcpar_pz[g4par]);
         }
      }
      h_eff_reintprotons_num->Fill(temp_geantpmom, event_weight);
      if(geant_pcontained==true){
         h_eff_containedprotons_num->Fill(temp_geantpmom, event_weight);
      } 
      //fill dEdx_vs_residual_range
      
      /*for(size_t ii =0; ii< t->pfp_reco_ismuoncandidate.size(); ++ii){
         if(t->pfp_reco_ismuoncandidate[ii]==1) continue;
         if(t->pfp_reco_dEdx[ii].size()<5) continue;
         for(unsigned int jj=0; jj<t->pfp_reco_dEdx[ii].size(); ++jj){
        	 h_dEdx_vs_rr_pcand->Fill(t->pfp_reco_RR[ii][jj], t->pfp_reco_dEdx[ii][jj]);
         }
      } 
      */

    } //end of if this is signal and cc1unp analysis
    } //end of is not data
    std::cout<<"Get the truth information of the CC1uNP Selected Events"<<"pind= "<<pind<<"  "<<t->pfp_reco_istrack[pind]<<std::endl;
    unsigned int vec_length = t->pfp_reco_dEdx[pind].size();
    for(size_t i =0; i<vec_length; ++i){   
	 h_dEdx_vs_rr_pcand->Fill(t->pfp_reco_RR[pind][i],t->pfp_reco_dEdx[pind][i]);
	 h_dEdx_vs_rr_pcand_2->Fill(t->pfp_reco_RR[pind][i],t->pfp_reco_dEdx[pind][vec_length - i]);
    }
    h_tmdqdx_vs_rr_pcand->Fill(t->pfp_reco_length[pind],t->pfp_reco_trunmeandqdx[pind]);

    for(size_t ii=0; ii<t->pfp_reco_ismuoncandidate.size(); ii++){
        if(t->pfp_reco_ismuoncandidate[ii]==1) continue;
        if(t->pfp_reco_isshower[ii]==0) continue;
        if(t->pfp_reco_numtracks[ii]!=1) continue;
        if(t->pfp_reco_dEdx[ii].size()<5) continue;
        for(size_t jj=0; jj<t->pfp_reco_dEdx[ii].size(); jj++){
                 h_dEdx_vs_rr_pshower->Fill(t->pfp_reco_RR[ii][jj], t->pfp_reco_dEdx[ii][jj]);
        }  
    }

    std::cout<<"checking the momentum smearing in the first bin"<<std::endl;
    if(isSignal && _ana_int_type=="cc1unp_analysis"){
    //total number of hits, angle py to p
    if((t->pfp_reco_Mom_muon[muind]<0.3) &&t->pfp_reco_Mom_muon[muind]>0.18 && t->pfp_truth_mom[muind]<0.18){
    h_nhits_lowmu->Fill(t->pfp_reco_dEdx[muind].size());
    h_trklen_truevsreco_lowmu->Fill(TMath::Sqrt((t->pfp_truth_startx[muind]-t->pfp_truth_endx[muind])*(t->pfp_truth_startx[muind]-t->pfp_truth_endx[muind])+(t->pfp_truth_starty[muind]-t->pfp_truth_endy[muind])*(t->pfp_truth_starty[muind]-t->pfp_truth_endy[muind])+(t->pfp_truth_startz[muind]-t->pfp_truth_endz[muind])*(t->pfp_truth_startz[muind]-t->pfp_truth_endz[muind])), t->pfp_reco_length[muind]);

    h_costhetay_lowmu->Fill(TMath::Sin(t->pfp_truth_theta[muind])*TMath::Sin(t->pfp_truth_phi[muind]));
    h_vtx_resolution_lowmu->Fill(t->vtx_resolution);
    h_nhits_lowmu_uplane->Fill(t->pfp_reco_nhits[muind]);
    h_nhits_lowmu_vplane->Fill(t->pfp_reco_nhits[muind]);
    }
    
    if(t->pfp_reco_Mom_muon[muind]<0.18 && t->pfp_truth_mom[muind]<0.18){
    h_nhits_firstbin->Fill(t->pfp_reco_dEdx[muind].size());
    h_trklen_truevsreco_firstbin->Fill(TMath::Sqrt((t->pfp_truth_startx[muind]-t->pfp_truth_endx[muind])*(t->pfp_truth_startx[muind]-t->pfp_truth_endx[muind])+(t->pfp_truth_starty[muind]-t->pfp_truth_endy[muind])*(t->pfp_truth_starty[muind]-t->pfp_truth_endy[muind])+(t->pfp_truth_startz[muind]-t->pfp_truth_endz[muind])*(t->pfp_truth_startz[muind]-t->pfp_truth_endz[muind])), t->pfp_reco_length[muind]);

    h_costhetay_firstbin->Fill(TMath::Sin(t->pfp_truth_theta[muind])*TMath::Sin(t->pfp_truth_phi[muind]));
    h_vtx_resolution_firstbin->Fill(t->vtx_resolution);
    h_nhits_firstbin_uplane->Fill(t->pfp_reco_nhits[muind]);
    h_nhits_firstbin_vplane->Fill(t->pfp_reco_nhits[muind]);
    }
        
    }
    std::cout<<"Started fill histograms after finish the calculations"<<std::endl;
    //============================================================================================================
    //
    // EVENT IS SELECTED
    //
    
    // if (isSignal && nu_origin)std::cout << ">>>>>>>>>>>>>>>>> Event is selected, " << t->run << ", " << t->subrun << ", " << t->event << ", slice " << scl_ll_max << std::endl;
    
    _event_histo_1d->hmap_onebin["total"]->Fill(0.5, event_weight);
    hmap_trklen["total"]->Fill(t->slc_longesttrack_length.at(scl_ll_max), event_weight);
    if(muon_contained) {
    hmap_trkmom_classic["total"]->Fill(t->pfp_reco_Mom_muon[muind], event_weight);
    }else {
    hmap_trkmom_classic["total"]->Fill(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight);
    }
    hmap_trkphi["total"]->Fill(t->slc_longesttrack_phi.at(scl_ll_max), event_weight);
    hmap_trktheta_classic["total"]->Fill(t->slc_longesttrack_theta.at(scl_ll_max), event_weight);
    hmap_multpfp["total"]->Fill(t->slc_mult_pfp.at(scl_ll_max), event_weight);
    hmap_multtracktol["total"]->Fill(t->slc_mult_track_tolerance.at(scl_ll_max), event_weight);

    if(muon_contained){
    if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->pfp_reco_Mom_muon[muind], event_weight, _event_histo_1d->hmap_trkmom_genie_pm1_bs, "total", fname_genie_pm1, wgts_genie_pm1); 
    } else {
    if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkmom_genie_pm1_bs, "total", fname_genie_pm1, wgts_genie_pm1);
    }
    if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->pfp_reco_Mom_proton[pind], event_weight, _event_histo_1d->hmap_trkpmom_genie_pm1_bs, "total", fname_genie_pm1, wgts_genie_pm1);
    if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trktheta_genie_pm1_bs, "total", fname_genie_pm1, wgts_genie_pm1);
    if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->pfp_reco_costheta[pind], event_weight, _event_histo_1d->hmap_trkptheta_genie_pm1_bs, "total", fname_genie_pm1, wgts_genie_pm1);

    if (!isdata && _fill_bootstrap_genie) FillBootstrap(thetamup, event_weight, _event_histo_1d->hmap_thetamup_genie_pm1_bs, "total", fname_genie_pm1, wgts_genie_pm1);

    std::cout<<"Filled histogram for GENIE unisim"<<std::endl;
    if (!isdata && _fill_bootstrap_genie) FillBootstrap(0.5, event_weight, _event_histo_1d->hmap_onebin_genie_multisim_bs, "total", fname_genie_multisim, wgts_genie_multisim);
    if(muon_contained){
    if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->pfp_reco_Mom_muon[muind], event_weight, _event_histo_1d->hmap_trkmom_genie_multisim_bs, "total", fname_genie_multisim, wgts_genie_multisim);
    } else{
    if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkmom_genie_multisim_bs, "total", fname_genie_multisim, wgts_genie_multisim);
    }
    if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkangle_genie_multisim_bs, "total", fname_genie_multisim, wgts_genie_multisim);
    if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo->hmap_trktheta_trkmom_genie_multisim_bs, "total", fname_genie_multisim, wgts_genie_multisim);
    if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo->hmap_trktheta_trkmom_poly_genie_multisim_bs, "total", fname_genie_multisim, wgts_genie_multisim);

    if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(0.5, event_weight, _event_histo_1d->hmap_onebin_extra_syst_multisim_bs, "total", fname_extra_syst, wgts_extra_syst);
    if(muon_contained){
    if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(t->pfp_reco_Mom_muon[muind], event_weight, _event_histo_1d->hmap_trkmom_extra_syst_multisim_bs, "total", fname_extra_syst, wgts_extra_syst);
    } else{
    if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkmom_extra_syst_multisim_bs, "total", fname_extra_syst, wgts_extra_syst);
    }
    if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkangle_extra_syst_multisim_bs, "total", fname_extra_syst, wgts_extra_syst);
    if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo->hmap_trktheta_trkmom_extra_syst_multisim_bs, "total", fname_extra_syst, wgts_extra_syst);
    if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo->hmap_trktheta_trkmom_poly_extra_syst_multisim_bs, "total", fname_extra_syst, wgts_extra_syst);

    if (!isdata && _fill_bootstrap_flux) FillBootstrap(0.5, event_weight, _event_histo_1d->hmap_onebin_flux_multisim_bs, "total", fname_flux_multisim, wgts_flux_multisim);
    if(muon_contained){
    if (!isdata && _fill_bootstrap_flux) FillBootstrap(t->pfp_reco_Mom_muon[muind], event_weight, _event_histo_1d->hmap_trkmom_flux_multisim_bs, "total", fname_flux_multisim, wgts_flux_multisim);
    } else {
    if (!isdata && _fill_bootstrap_flux) FillBootstrap(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkmom_flux_multisim_bs, "total", fname_flux_multisim, wgts_flux_multisim);
    }
    if (!isdata && _fill_bootstrap_flux) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkangle_flux_multisim_bs, "total", fname_flux_multisim, wgts_flux_multisim);
    if (!isdata && _fill_bootstrap_flux) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo->hmap_trktheta_trkmom_flux_multisim_bs, "total", fname_flux_multisim, wgts_flux_multisim);
    if (!isdata && _fill_bootstrap_flux) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo->hmap_trktheta_trkmom_poly_flux_multisim_bs, "total", fname_flux_multisim, wgts_flux_multisim);

    if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(0.5, event_weight, _event_histo_1d->hmap_onebin_mc_stat_multisim_bs, "total", fname_mc_stat_multisim, wgts_mc_stat_multisim);
    if(muon_contained){
    if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(t->pfp_reco_Mom_muon[muind], event_weight, _event_histo_1d->hmap_trkmom_mc_stat_multisim_bs, "total", fname_mc_stat_multisim, wgts_mc_stat_multisim);
    } else {
    if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkmom_mc_stat_multisim_bs, "total", fname_mc_stat_multisim, wgts_mc_stat_multisim);
    }
    if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkangle_mc_stat_multisim_bs, "total", fname_mc_stat_multisim, wgts_mc_stat_multisim);
    if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo->hmap_trktheta_trkmom_mc_stat_multisim_bs, "total", fname_mc_stat_multisim, wgts_mc_stat_multisim);
    if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo->hmap_trktheta_trkmom_poly_mc_stat_multisim_bs, "total", fname_mc_stat_multisim, wgts_mc_stat_multisim);

    hmap_trkplen["total"]->Fill(t->pfp_reco_length[pind], event_weight);
    hmap_trkpmom_classic["total"]->Fill(t->pfp_reco_Mom_proton[pind], event_weight);
    hmap_trkptheta_classic["total"]->Fill(t->pfp_reco_costheta[pind], event_weight);
    hmap_trkpphi["total"]->Fill(t->pfp_reco_phi[pind], event_weight);
    hmap_thetamup["total"]->Fill(thetamup, event_weight);
    hmap_ptmis["total"]->Fill(ptmis, event_weight);
    hmap_etatest["total"]->Fill(etatest, event_weight);
    hmap_alphat["total"]->Fill(alphat, event_weight);
    hmap_phit["total"]->Fill(phit, event_weight);
    hmap_enucal["total"]->Fill(enucal, event_weight);
    hmap_nhits_leadingp["total"]->Fill(t->pfp_reco_dEdx[pind].size(), event_weight);


    if(_showerastrack) {
      hmap_pmult["total"]->Fill(t->num_pfp-1, event_weight);
    } else{
      hmap_pmult["total"]->Fill(t->num_pfp_tracks-1, event_weight);
    }
    std::cout<<"Filled histograms of the hmap pmult"<<std::endl;

    if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->pfp_reco_Mom_proton[pind], event_weight, _event_histo_1d->hmap_trkpmom_genie_multisim_bs, "total", fname_genie_multisim, wgts_genie_multisim);
    if (!isdata && _fill_bootstrap_flux) FillBootstrap(t->pfp_reco_Mom_proton[pind], event_weight, _event_histo_1d->hmap_trkpmom_flux_multisim_bs, "total", fname_flux_multisim, wgts_flux_multisim);
    if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(t->pfp_reco_Mom_proton[pind], event_weight, _event_histo_1d->hmap_trkpmom_extra_syst_multisim_bs, "total", fname_extra_syst, wgts_extra_syst);
    if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(t->pfp_reco_Mom_proton[pind], event_weight, _event_histo_1d->hmap_trkpmom_mc_stat_multisim_bs, "total", fname_mc_stat_multisim, wgts_mc_stat_multisim);

    if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->pfp_reco_costheta[pind], event_weight, _event_histo_1d->hmap_trkpangle_genie_multisim_bs, "total", fname_genie_multisim, wgts_genie_multisim);
    if (!isdata && _fill_bootstrap_flux) FillBootstrap(t->pfp_reco_costheta[pind], event_weight, _event_histo_1d->hmap_trkpangle_flux_multisim_bs, "total", fname_flux_multisim, wgts_flux_multisim);
    if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(t->pfp_reco_costheta[pind], event_weight, _event_histo_1d->hmap_trkpangle_extra_syst_multisim_bs, "total", fname_extra_syst, wgts_extra_syst);
    if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(t->pfp_reco_costheta[pind], event_weight, _event_histo_1d->hmap_trkpangle_mc_stat_multisim_bs, "total", fname_mc_stat_multisim, wgts_mc_stat_multisim);

    if (!isdata && _fill_bootstrap_genie) FillBootstrap(thetamup, event_weight, _event_histo_1d->hmap_thetamup_genie_multisim_bs, "total", fname_genie_multisim, wgts_genie_multisim);
    if (!isdata && _fill_bootstrap_flux) FillBootstrap(thetamup, event_weight, _event_histo_1d->hmap_thetamup_flux_multisim_bs, "total", fname_flux_multisim, wgts_flux_multisim);
    if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(thetamup, event_weight, _event_histo_1d->hmap_thetamup_extra_syst_multisim_bs, "total", fname_extra_syst, wgts_extra_syst);
    if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(thetamup, event_weight, _event_histo_1d->hmap_thetamup_mc_stat_multisim_bs, "total", fname_mc_stat_multisim, wgts_mc_stat_multisim);

    if(muon_contained){
    _event_histo_1d->hmap_trkmom["total"]->Fill(t->pfp_reco_Mom_muon[muind], event_weight);
    } else {
    _event_histo_1d->hmap_trkmom["total"]->Fill(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight);
    }
    _event_histo_1d->hmap_trktheta["total"]->Fill(t->slc_longesttrack_theta.at(scl_ll_max), event_weight);

    _event_histo->hmap_trktheta_trkmom["total"]->Fill(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight);
    _event_histo->hmap_trktheta_trkmom_poly["total"]->Fill(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight);

    _event_histo_1d->hmap_trkpmom["total"]->Fill(t->pfp_reco_Mom_proton[pind], event_weight);
    _event_histo_1d->hmap_trkptheta["total"]->Fill(t->pfp_reco_costheta[pind], event_weight);

    _event_histo_1d->hmap_thetamup["total"]->Fill(thetamup, event_weight); 

    if (isSignal&&trackfromneutrino) {
      h_true_nu_eng_aftersel->Fill(t->nu_e, event_weight);
    }
    //====================================================================================================================================================================

    std::cout<<"Filled histograms and tree leaves"<<"  isSignal and trackfromneutrino "<<isSignal<<" "<<trackfromneutrino<<std::endl;
    if (isSignal&&trackfromneutrino) {

      // Fill the true-reco TTree for the nominal case
      _mom_true = t->true_muon_mom;
      if(muon_contained){
      _mom_mcs = t->pfp_reco_Mom_muon[muind];
      } else {
      _mom_mcs = t->slc_muoncandidate_mom_mcs.at(scl_ll_max);
      }
      _pmom_true = temp_pmom;
      _pmom_reco = t->pfp_reco_Mom_proton[pind];
      _contained = t->slc_muoncandidate_contained.at(scl_ll_max);
      _selected = true;

      _angle_true = t->lep_costheta;
      _angle_reco = t->slc_longesttrack_theta.at(scl_ll_max); //t->slc_muoncandidate_theta.at(scl_ll_max);

      _event_weight_fortree = event_weight;
      
      _wgtsnames_genie_multisim = fname_genie_multisim;
      _wgts_genie_multisim = wgts_genie_multisim;

      _wgtsnames_extra_syst = fname_extra_syst;
      _wgts_extra_syst = wgts_extra_syst;

      _wgtsnames_flux_multisim = fname_flux_multisim;
      _wgts_flux_multisim = wgts_flux_multisim;

      _true_reco_tree->Fill();

      // *** Migr mat addition
      int m = _event_histo->h_reco_per_true[0][0]->GetXaxis()->FindBin(_angle_true) - 1; // true bin
      int n = _event_histo->h_reco_per_true[0][0]->GetYaxis()->FindBin(_mom_true) - 1;  // true bin
      if (m >= 0 && n >= 0 
          && m < _event_histo->h_reco_per_true[0][0]->GetNbinsX()    // Avoid overflows
          && n < _event_histo->h_reco_per_true[0][0]->GetNbinsY()) { // Avoid overflows
        _event_histo->h_reco_per_true[m][n]->Fill(_angle_reco, _mom_mcs, event_weight);
        if(!isdata && _fill_bootstrap_genie) FillBootstrap(_angle_reco, _mom_mcs, m, n, event_weight, _event_histo->bs_genie_multisim_reco_per_true, fname_genie_multisim, wgts_genie_multisim);
        if(!isdata && _fill_bootstrap_extra_syst) FillBootstrap(_angle_reco, _mom_mcs, m, n, event_weight, _event_histo->bs_extra_syst_multisim_reco_per_true, fname_extra_syst, wgts_extra_syst);
        if(!isdata && _fill_bootstrap_flux) FillBootstrap(_angle_reco, _mom_mcs, m, n, event_weight, _event_histo->bs_flux_multisim_reco_per_true, fname_flux_multisim, wgts_flux_multisim);
        if(!isdata && _fill_bootstrap_mc_stat) FillBootstrap(_angle_reco, _mom_mcs, m, n, event_weight, _event_histo->bs_mc_stat_multisim_reco_per_true, fname_mc_stat_multisim, wgts_mc_stat_multisim);
      }

      // For the migration matrix (poly)
      m = _event_histo->h_poly_reco_per_true[0]->FindBin(_angle_true, _mom_true); // true bin
      int i = _event_histo->h_poly_reco_per_true[0]->FindBin(_angle_reco, _mom_mcs); // reco bin
      if (m < 0) m = 0; // Negative bins are overflows, and are all added in entry 0 of the vector
      if (m < _event_histo->h_poly_reco_per_true[0]->GetNumberOfBins()+1) {
        _event_histo->h_poly_reco_per_true[m]->Fill(_angle_reco, _mom_mcs, event_weight);
        if(!isdata && _fill_bootstrap_genie) FillBootstrap(m, i, event_weight, _event_histo->bs_genie_multisim_poly_reco_per_true, fname_genie_multisim, wgts_genie_multisim);
        if(!isdata && _fill_bootstrap_flux) FillBootstrap(m, i, event_weight, _event_histo->bs_flux_multisim_poly_reco_per_true, fname_flux_multisim, wgts_flux_multisim);
        if(!isdata && _fill_bootstrap_extra_syst) FillBootstrap(m, i, event_weight, _event_histo->bs_extra_syst_multisim_poly_reco_per_true, fname_extra_syst, wgts_extra_syst);
        if(!isdata && _fill_bootstrap_mc_stat) FillBootstrap(m, i, event_weight, _event_histo->bs_mc_stat_multisim_poly_reco_per_true, fname_mc_stat_multisim, wgts_mc_stat_multisim);
      }

      // // Also fill the same tree for all te universes
      // FillTrueRecoTree(tmap_mom_tree_gene_multisim_bs, _mom_true, _mom_mcs, _angle_true, _angle_reco, fname_genie_multisim, wgts_genie_multisim);
      _event_histo_1d->h_true_reco_pmom->Fill(t->pfp_truth_mom[pind], t->pfp_reco_Mom_proton[pind], event_weight);
      _event_histo_1d->h_true_reco_pcostheta->Fill(temp_pangle, t->pfp_reco_costheta[pind], event_weight);
      _event_histo_1d->h_true_reco_thetamup->Fill(temp_thetamup, thetamup, event_weight);
      

      _event_histo_1d->h_true_reco_mom->Fill(_mom_true, _mom_mcs, event_weight);
      _event_histo_1d->h_true_reco_costheta->Fill(_angle_true, _angle_reco, event_weight);
      //===============================================================================================================================================
      std::cout<<"libo test 0"<<std::endl;
      /*if(!isdata && _fill_bootstrap_genie) FillBootstrap_test(_mom_true, _mom_mcs, event_weight, _event_histo_1d->bs_genie_pm1_true_reco_mom, fname_genie_pm1, wgts_genie_pm1);
      if(!isdata && _fill_bootstrap_genie) FillBootstrap_test(temp_pmom, t->pfp_reco_Mom_proton[pind], event_weight, _event_histo_1d->bs_genie_pm1_true_reco_pmom, fname_genie_pm1, wgts_genie_pm1);
      if(!isdata && _fill_bootstrap_genie) FillBootstrap_test(_angle_true, _angle_reco, event_weight, _event_histo_1d->bs_genie_pm1_true_reco_muangle, fname_genie_pm1, wgts_genie_pm1);
      if(!isdata && _fill_bootstrap_genie) FillBootstrap_test(temp_pangle, t->pfp_reco_costheta[pind], event_weight, _event_histo_1d->bs_genie_pm1_true_reco_pangle, fname_genie_pm1, wgts_genie_pm1);
      if(!isdata && _fill_bootstrap_genie) FillBootstrap_test(temp_thetamup, thetamup, event_weight, _event_histo_1d->bs_genie_pm1_true_reco_thetamup, fname_genie_pm1, wgts_genie_pm1);
      std::cout<<"libo test 1"<<std::endl;
      */
      if(!isdata && _fill_bootstrap_genie) FillBootstrap(_mom_true, _mom_mcs, event_weight, _event_histo_1d->bs_genie_pm1_true_reco_mom, fname_genie_pm1, wgts_genie_pm1);
      if(!isdata && _fill_bootstrap_genie) FillBootstrap(temp_pmom, t->pfp_reco_Mom_proton[pind], event_weight, _event_histo_1d->bs_genie_pm1_true_reco_pmom, fname_genie_pm1, wgts_genie_pm1);
      if(!isdata && _fill_bootstrap_genie) FillBootstrap(_angle_true, _angle_reco, event_weight, _event_histo_1d->bs_genie_pm1_true_reco_muangle, fname_genie_pm1, wgts_genie_pm1);
      if(!isdata && _fill_bootstrap_genie) FillBootstrap(temp_pangle, t->pfp_reco_costheta[pind], event_weight, _event_histo_1d->bs_genie_pm1_true_reco_pangle, fname_genie_pm1, wgts_genie_pm1);
      if(!isdata && _fill_bootstrap_genie) FillBootstrap(temp_thetamup, thetamup, event_weight, _event_histo_1d->bs_genie_pm1_true_reco_thetamup, fname_genie_pm1, wgts_genie_pm1);
       
      /*
     */

      if(!isdata && _fill_bootstrap_genie) _event_histo_1d->bs_genie_multisim_true_reco_mumom->Fill(_mom_true, _mom_mcs, event_weight, wgts_genie_multisim);
      if(!isdata && _fill_bootstrap_genie) _event_histo_1d->bs_genie_multisim_true_reco_muangle->Fill(_angle_true, _angle_reco, event_weight, wgts_genie_multisim);
      if(!isdata && _fill_bootstrap_extra_syst) _event_histo_1d->bs_extra_syst_true_reco_mumom->Fill(_mom_true, _mom_mcs, event_weight, wgts_extra_syst);
      if(!isdata && _fill_bootstrap_extra_syst) _event_histo_1d->bs_extra_syst_true_reco_muangle->Fill(_angle_true, _angle_reco, event_weight, wgts_extra_syst);
      if(!isdata && _fill_bootstrap_flux) _event_histo_1d->bs_flux_multisim_true_reco_mumom->Fill(_mom_true, _mom_mcs, event_weight, wgts_flux_multisim);
      if(!isdata && _fill_bootstrap_flux) _event_histo_1d->bs_flux_multisim_true_reco_muangle->Fill(_angle_true, _angle_reco, event_weight, wgts_flux_multisim);
      if(!isdata && _fill_bootstrap_mc_stat) _event_histo_1d->bs_mc_stat_multisim_true_reco_mumom->Fill(_mom_true, _mom_mcs, event_weight, wgts_mc_stat_multisim);
      if(!isdata && _fill_bootstrap_mc_stat) _event_histo_1d->bs_mc_stat_multisim_true_reco_muangle->Fill(_angle_true, _angle_reco, event_weight, wgts_mc_stat_multisim);

      if(!isdata && _fill_bootstrap_genie) _event_histo_1d->bs_genie_multisim_true_reco_pmom->Fill(temp_pmom, t->pfp_reco_Mom_proton[pind], event_weight, wgts_genie_multisim);
      if(!isdata && _fill_bootstrap_flux) _event_histo_1d->bs_flux_multisim_true_reco_pmom->Fill(temp_pmom, t->pfp_reco_Mom_proton[pind], event_weight, wgts_flux_multisim);
      if(!isdata && _fill_bootstrap_extra_syst) _event_histo_1d->bs_extra_syst_true_reco_pmom->Fill(temp_pmom, t->pfp_reco_Mom_proton[pind], event_weight, wgts_extra_syst);
      if(!isdata && _fill_bootstrap_mc_stat) _event_histo_1d->bs_mc_stat_multisim_true_reco_pmom->Fill(temp_pmom, t->pfp_reco_Mom_proton[pind], event_weight, wgts_mc_stat_multisim);

      if(!isdata && _fill_bootstrap_genie) _event_histo_1d->bs_genie_multisim_true_reco_pangle->Fill(temp_pangle, t->pfp_reco_costheta[pind], event_weight, wgts_genie_multisim);
      if(!isdata && _fill_bootstrap_flux) _event_histo_1d->bs_flux_multisim_true_reco_pangle->Fill(temp_pangle, t->pfp_reco_costheta[pind], event_weight, wgts_flux_multisim);
      if(!isdata && _fill_bootstrap_extra_syst) _event_histo_1d->bs_extra_syst_true_reco_pangle->Fill(temp_pangle, t->pfp_reco_costheta[pind], event_weight, wgts_extra_syst);
      if(!isdata && _fill_bootstrap_mc_stat) _event_histo_1d->bs_mc_stat_multisim_true_reco_pangle->Fill(temp_pangle, t->pfp_reco_costheta[pind], event_weight, wgts_mc_stat_multisim);

      if(!isdata && _fill_bootstrap_genie) _event_histo_1d->bs_genie_multisim_true_reco_thetamup->Fill(temp_thetamup, thetamup, event_weight, wgts_genie_multisim);
      if(!isdata && _fill_bootstrap_flux) _event_histo_1d->bs_flux_multisim_true_reco_thetamup->Fill(temp_thetamup, thetamup, event_weight, wgts_flux_multisim);
      if(!isdata && _fill_bootstrap_extra_syst) _event_histo_1d->bs_extra_syst_true_reco_thetamup->Fill(temp_thetamup, thetamup, event_weight, wgts_extra_syst);
      if(!isdata && _fill_bootstrap_mc_stat) _event_histo_1d->bs_mc_stat_multisim_true_reco_thetamup->Fill(temp_thetamup, thetamup, event_weight, wgts_mc_stat_multisim);


    }
    std::cout<<"Start to fill histograms for the CCinclusive variables"<<std::endl;
    int true_pdg = t->slc_muoncandidate_truepdg.at(scl_ll_max);
    //std::cout<<"Start to fill histograms of dQds"<<std::endl; 
    //
    // Fill dQ/ds histograms
    //
    hmap_dqdx_trunc["total"]->Fill(dqdx_calib, event_weight);
    h_dqdx_trunc_length->Fill(dqdx_calib, t->slc_muoncandidate_length.at(scl_ll_max), event_weight);
    if (true_pdg == 13 || true_pdg == -13) {
      hmap_dqdx_trunc["muon"]->Fill(dqdx_calib, event_weight);
      h_dqdx_trunc_length_muon->Fill(dqdx_calib, t->slc_muoncandidate_length.at(scl_ll_max), event_weight);
      if (dqdx_calib >= 0. && dqdx_calib <=200000. && t->slc_muoncandidate_length.at(scl_ll_max) < 2000) {
        _csvfile << dqdx_calib << ","
                 << t->slc_muoncandidate_length.at(scl_ll_max) << "," << "1" << std::endl;
      }
    } else if (true_pdg == 211 || true_pdg == -211) {
      hmap_dqdx_trunc["pion"]->Fill(dqdx_calib, event_weight);
    } else if (true_pdg == 2212) {
      hmap_dqdx_trunc["proton"]->Fill(dqdx_calib, event_weight);
      h_dqdx_trunc_length_proton->Fill(dqdx_calib, t->slc_muoncandidate_length.at(scl_ll_max), event_weight);
      if (dqdx_calib >= 0. && dqdx_calib <=200000. && t->slc_muoncandidate_length.at(scl_ll_max) < 2000) {
        _csvfile << dqdx_calib << ","
                 << t->slc_muoncandidate_length.at(scl_ll_max) << "," << "0" << std::endl;
      }
    } else if (true_pdg == 22) {
      hmap_dqdx_trunc["photon"]->Fill(dqdx_calib, event_weight);
    } else if (true_pdg == 11 || true_pdg == -11) {
      hmap_dqdx_trunc["electron"]->Fill(dqdx_calib, event_weight);
    } else {
      hmap_dqdx_trunc["else"]->Fill(dqdx_calib, event_weight);
    }
    
    double hypo_pe = 0;
    for (int pmt = 0; pmt < 32; pmt++) {
      hypo_pe += (t->slc_flshypo_spec.at(scl_ll_max))[pmt];
    }
    hmap_xdiff["total"]->Fill(t->slc_flsmatch_qllx.at(scl_ll_max) - t->slc_flsmatch_tpcx.at(scl_ll_max), event_weight);
    // reintro hmap_zdiff["total"]->Fill(t->slc_flsmatch_hypoz.at(scl_ll_max) - t->beamfls_z.at(flashInBeamSpill), event_weight);
    // reintro hmap_pediff["total"]->Fill(hypo_pe - t->beamfls_pe.at(flashInBeamSpill), event_weight);
        
    hmap_vtxx["total"]->Fill(t->slc_nuvtx_x.at(scl_ll_max), event_weight);
    hmap_vtxy["total"]->Fill(t->slc_nuvtx_y.at(scl_ll_max), event_weight);
    hmap_vtxz["total"]->Fill(t->slc_nuvtx_z.at(scl_ll_max), event_weight);
    
    double second_score = -9999, score_difference = -9999;
    if (temp_score.size() > 1) {
      second_score = temp_score.at(1);
      score_difference = temp_score.at(0) - temp_score.at(1);
    }

    hmap_flsmatch_score["total"]->Fill(t->slc_flsmatch_score.at(scl_ll_max), event_weight);
    hmap_flsmatch_score_second["total"]->Fill(second_score, event_weight);
    hmap_flsmatch_score_difference["total"]->Fill(score_difference, event_weight);


    // SIGNAL
    if ( isSignal && nu_origin && trackfromneutrino) {
      n_signal ++;
      hmap_xdiff["signal"]->Fill(t->slc_flsmatch_qllx.at(scl_ll_max) - t->slc_flsmatch_tpcx.at(scl_ll_max), event_weight);
      // reintro hmap_zdiff["signal"]->Fill(t->slc_flsmatch_hypoz.at(scl_ll_max) - t->beamfls_z.at(flashInBeamSpill), event_weight);
      // reintro hmap_pediff["signal"]->Fill(hypo_pe - t->beamfls_pe.at(flashInBeamSpill), event_weight);
      
      hmap_vtxx["signal"]->Fill(t->slc_nuvtx_x.at(scl_ll_max), event_weight);
      hmap_vtxy["signal"]->Fill(t->slc_nuvtx_y.at(scl_ll_max), event_weight);
      hmap_vtxz["signal"]->Fill(t->slc_nuvtx_z.at(scl_ll_max), event_weight);
      
      hmap_flsmatch_score["signal"]->Fill(t->slc_flsmatch_score.at(scl_ll_max), event_weight);
      hmap_flsmatch_score_second["signal"]->Fill(second_score, event_weight);
      hmap_flsmatch_score_difference["signal"]->Fill(score_difference, event_weight);
    }
    // BACKGROUND
    else {
      hmap_xdiff["background"]->Fill(t->slc_flsmatch_qllx.at(scl_ll_max) - t->slc_flsmatch_tpcx.at(scl_ll_max), event_weight);
      // reintro hmap_zdiff["background"]->Fill(t->slc_flsmatch_hypoz.at(scl_ll_max) - t->beamfls_z.at(flashInBeamSpill), event_weight);
      // reintro hmap_pediff["background"]->Fill(hypo_pe - t->beamfls_pe.at(flashInBeamSpill), event_weight);

      hmap_vtxx["background"]->Fill(t->slc_nuvtx_x.at(scl_ll_max), event_weight);
      hmap_vtxy["background"]->Fill(t->slc_nuvtx_y.at(scl_ll_max), event_weight);
      hmap_vtxz["background"]->Fill(t->slc_nuvtx_z.at(scl_ll_max), event_weight);
      
      hmap_flsmatch_score["background"]->Fill(t->slc_flsmatch_score.at(scl_ll_max), event_weight);
      hmap_flsmatch_score_second["background"]->Fill(second_score, event_weight);
      hmap_flsmatch_score_difference["background"]->Fill(score_difference, event_weight);
    }
    
    if (isNueCCFV && trackfromneutrino) {
      nue_cc_selected_total+=t->bnb_weight;
      if (t->nu_e >= 0.05 && t->nu_e <= 1.5){
        nue_cc_selected_total_energy_range+=t->bnb_weight;
      }
    }
    if (isNue && trackfromneutrino) {
      if (t->nu_e >= 0.05 && t->nu_e <= 1.5 && t->ccnc==0 
        && t->tvtx_x[0] > 0. && t->tvtx_x[0] < 256.35
        && t->tvtx_y[0] > -116.5 && t->tvtx_y[0] < 116.5
        && t->tvtx_z[0] > 0. && t->tvtx_y[0] < 1036.8){
        nue_selected_total_energy_range+=t->bnb_weight;

        h_nue_selected_energy->Fill(t->nu_e, event_weight); 
        if (std::abs(true_pdg) == 11) n_nue_electron+=t->bnb_weight;
        if (std::abs(true_pdg) == 2212) n_nue_proton+=t->bnb_weight;
        if (std::abs(true_pdg) == 211) n_nue_pion+=t->bnb_weight;
      }
    }

    std::cout<<"Start to fill histograms nuorigin"<<nu_origin<<isSignal<<trackfromneutrino<<std::endl;
    //
    // SIGNAL
    //
    if(nu_origin && isSignal && trackfromneutrino /*&& true_pdg==13 && true_origin == 0*/) {
      
      std::cout << "Is signal and is selected. event: " << t->event << std::endl;

      signal_sel += event_weight;
      _event_histo_1d->hmap_onebin["signal"]->Fill(0.5, event_weight);
      _event_histo_1d->h_eff_onebin_num->Fill(0.5, event_weight);
      h_eff_num->Fill(t->nu_e, event_weight);

      hmap_trklen["signal"]->Fill(t->slc_longesttrack_length.at(scl_ll_max), event_weight);
      if(muon_contained) {
      hmap_trkmom_classic["signal"]->Fill(t->pfp_reco_Mom_muon[muind], event_weight);
      } else {
      hmap_trkmom_classic["signal"]->Fill(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight);
      }
      hmap_trkphi["signal"]->Fill(t->slc_longesttrack_phi.at(scl_ll_max), event_weight);
      hmap_trktheta_classic["signal"]->Fill(t->slc_longesttrack_theta.at(scl_ll_max), event_weight);
      hmap_multpfp["signal"]->Fill(t->slc_mult_pfp.at(scl_ll_max), event_weight);
      hmap_multtracktol["signal"]->Fill(t->slc_mult_track_tolerance.at(scl_ll_max), event_weight);
       //hmap_trkmom_genie_pm1_bs["signal"]->Fill(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), 1., wgts_genie_pm1);
      //=======================================================================================================================================================

      if (!isdata && _fill_bootstrap_genie) bs_genie_pm1_eff_mumom_num.Fill(t->true_muon_mom, event_weight, wgts_genie_pm1);
      if (!isdata && _fill_bootstrap_genie) bs_genie_pm1_eff_pmom_num.Fill(temp_pmom, event_weight, wgts_genie_pm1);
      if (!isdata && _fill_bootstrap_genie) bs_genie_pm1_eff_muangle_num.Fill(t->lep_costheta, event_weight, wgts_genie_pm1);
      if (!isdata && _fill_bootstrap_genie) bs_genie_pm1_eff_pangle_num.Fill(temp_pangle, event_weight, wgts_genie_pm1);
      if (!isdata && _fill_bootstrap_genie) bs_genie_pm1_eff_thetamup_num.Fill(temp_thetamup, event_weight, wgts_genie_pm1);


      if (!isdata && _fill_bootstrap_genie) _event_histo_1d->bs_genie_multisim_eff_onebin_num->Fill(0.5, event_weight, wgts_genie_multisim);

      if (!isdata && _fill_bootstrap_genie) _event_histo_1d->bs_genie_multisim_eff_mumom_num->Fill(t->true_muon_mom, event_weight, wgts_genie_multisim);
      if (!isdata && _fill_bootstrap_genie) _event_histo_1d->bs_genie_multisim_eff_pmom_num->Fill(temp_pmom, event_weight, wgts_genie_multisim);
      if (!isdata && _fill_bootstrap_genie) _event_histo_1d->bs_genie_multisim_eff_muangle_num->Fill(t->lep_costheta, event_weight, wgts_genie_multisim);
      if (!isdata && _fill_bootstrap_genie) _event_histo_1d->bs_genie_multisim_eff_pangle_num->Fill(temp_pangle, event_weight, wgts_genie_multisim);
      if (!isdata && _fill_bootstrap_genie) _event_histo_1d->bs_genie_multisim_eff_thetamup_num->Fill(temp_thetamup, event_weight, wgts_genie_multisim);

      if (!isdata && _fill_bootstrap_genie) _event_histo->bs_genie_multisim_eff_muangle_mumom_num->Fill(t->lep_costheta, t->true_muon_mom, event_weight, wgts_genie_multisim);
      if (!isdata && _fill_bootstrap_genie) _event_histo->bs_genie_multisim_eff_poly_muangle_mumom_num->Fill(t->lep_costheta, t->true_muon_mom, event_weight, wgts_genie_multisim);

      if (!isdata && _fill_bootstrap_extra_syst) _event_histo_1d->bs_extra_syst_multisim_eff_onebin_num->Fill(0.5, event_weight, wgts_extra_syst);

      if (!isdata && _fill_bootstrap_extra_syst) _event_histo_1d->bs_extra_syst_multisim_eff_mumom_num->Fill(t->true_muon_mom, event_weight, wgts_extra_syst);
      if (!isdata && _fill_bootstrap_extra_syst) _event_histo_1d->bs_extra_syst_multisim_eff_pmom_num->Fill(temp_pmom, event_weight, wgts_extra_syst);
      if (!isdata && _fill_bootstrap_extra_syst) _event_histo_1d->bs_extra_syst_multisim_eff_muangle_num->Fill(t->lep_costheta, event_weight, wgts_extra_syst);
      if (!isdata && _fill_bootstrap_extra_syst) _event_histo_1d->bs_extra_syst_multisim_eff_pangle_num->Fill(temp_pangle, event_weight, wgts_extra_syst);
      if (!isdata && _fill_bootstrap_extra_syst) _event_histo_1d->bs_extra_syst_multisim_eff_thetamup_num->Fill(temp_thetamup, event_weight, wgts_extra_syst);

      if (!isdata && _fill_bootstrap_extra_syst) _event_histo->bs_extra_syst_multisim_eff_muangle_mumom_num->Fill(t->lep_costheta, t->true_muon_mom, event_weight, wgts_extra_syst);
      if (!isdata && _fill_bootstrap_extra_syst) _event_histo->bs_extra_syst_multisim_eff_poly_muangle_mumom_num->Fill(t->lep_costheta, t->true_muon_mom, event_weight, wgts_extra_syst);

      if (!isdata && _fill_bootstrap_flux) _event_histo_1d->bs_flux_multisim_eff_onebin_num->Fill(0.5, event_weight, wgts_flux_multisim);

      if (!isdata && _fill_bootstrap_flux) _event_histo_1d->bs_flux_multisim_eff_mumom_num->Fill(t->true_muon_mom, event_weight, wgts_flux_multisim);
      if (!isdata && _fill_bootstrap_flux) _event_histo_1d->bs_flux_multisim_eff_pmom_num->Fill(temp_pmom, event_weight, wgts_flux_multisim);
      if (!isdata && _fill_bootstrap_flux) _event_histo_1d->bs_flux_multisim_eff_muangle_num->Fill(t->lep_costheta, event_weight, wgts_flux_multisim);
      if (!isdata && _fill_bootstrap_flux) _event_histo_1d->bs_flux_multisim_eff_pangle_num->Fill(temp_pangle, event_weight, wgts_flux_multisim);
      if (!isdata && _fill_bootstrap_flux) _event_histo_1d->bs_flux_multisim_eff_thetamup_num->Fill(temp_thetamup, event_weight, wgts_flux_multisim);

      if (!isdata && _fill_bootstrap_flux) _event_histo->bs_flux_multisim_eff_muangle_mumom_num->Fill(t->lep_costheta, t->true_muon_mom, event_weight, wgts_flux_multisim);
      if (!isdata && _fill_bootstrap_flux) _event_histo->bs_flux_multisim_eff_poly_muangle_mumom_num->Fill(t->lep_costheta, t->true_muon_mom, event_weight, wgts_flux_multisim);

      if (!isdata && _fill_bootstrap_mc_stat) _event_histo_1d->bs_mc_stat_multisim_eff_onebin_num->Fill(0.5, event_weight, wgts_mc_stat_multisim);

      if (!isdata && _fill_bootstrap_mc_stat) _event_histo_1d->bs_mc_stat_multisim_eff_mumom_num->Fill(t->true_muon_mom, event_weight, wgts_mc_stat_multisim);
      if (!isdata && _fill_bootstrap_mc_stat) _event_histo_1d->bs_mc_stat_multisim_eff_pmom_num->Fill(temp_pmom, event_weight, wgts_mc_stat_multisim);
      if (!isdata && _fill_bootstrap_mc_stat) _event_histo_1d->bs_mc_stat_multisim_eff_muangle_num->Fill(t->lep_costheta, event_weight, wgts_mc_stat_multisim);
      if (!isdata && _fill_bootstrap_mc_stat) _event_histo_1d->bs_mc_stat_multisim_eff_pangle_num->Fill(temp_pangle, event_weight, wgts_mc_stat_multisim);
      if (!isdata && _fill_bootstrap_mc_stat) _event_histo_1d->bs_mc_stat_multisim_eff_thetamup_num->Fill(temp_thetamup, event_weight, wgts_mc_stat_multisim);

      if (!isdata && _fill_bootstrap_mc_stat) _event_histo->bs_mc_stat_multisim_eff_muangle_mumom_num->Fill(t->lep_costheta, t->true_muon_mom, event_weight, wgts_mc_stat_multisim);
      if (!isdata && _fill_bootstrap_mc_stat) _event_histo->bs_mc_stat_multisim_eff_poly_muangle_mumom_num->Fill(t->lep_costheta, t->true_muon_mom, event_weight, wgts_mc_stat_multisim);

      _event_histo_1d->h_eff_mumom_num->Fill(t->true_muon_mom, event_weight);
      _event_histo_1d->h_eff_muangle_num->Fill(t->lep_costheta, event_weight);
      _event_histo->h_eff_muangle_mumom_num->Fill(t->lep_costheta, t->true_muon_mom, event_weight);
      _event_histo->h_eff_muangle_mumom_poly_num->Fill(t->lep_costheta, t->true_muon_mom, event_weight);
      _event_histo_1d->h_eff_muphi_num->Fill(t->lep_phi, event_weight);
      /*
      */

      //*** pmom angle den already filled in the loop of getting true leading proton momentum

      h_eff_mult_num->Fill(t->genie_mult, event_weight);
      h_eff_mult_ch_num->Fill(t->genie_mult_ch, event_weight);
      h_mu_eff_mom_sel->Fill(t->true_muon_mom, t->muon_reco_eff, event_weight);

      if (t->mode == 0) {
        h_eff_qe_num->Fill(t->nu_e, event_weight);
        signal_sel_qe += event_weight;
      }
      if (t->mode == 1) {
        h_eff_res_num->Fill(t->nu_e, event_weight);
        signal_sel_res += event_weight;
      }
      if (t->mode == 2) {
        h_eff_dis_num->Fill(t->nu_e, event_weight);
        signal_sel_dis += event_weight;
      }
      if (t->mode == 3) {
        h_eff_coh_num->Fill(t->nu_e, event_weight);
        signal_sel_coh += event_weight;
      }
      if (t->mode == 10) {
        h_eff_mec_num->Fill(t->nu_e, event_weight);
        signal_sel_mec += event_weight;
      }
      std::cout<<"Filled histograms for efficiency calculation"<<std::endl;
      // Also save themc truth histogram per interaction type
      hmap_mctruth_nuenergy["total"]->Fill(t->nu_e, event_weight);
      hmap_mctruth_mumom["total"]->Fill(t->true_muon_mom, event_weight);
      hmap_mctruth_mucostheta["total"]->Fill(t->lep_costheta, event_weight);
      hmap_mctruth_muphi["total"]->Fill(t->lep_phi, event_weight);
      hmap_mctruth_chargedmult["total"]->Fill(t->genie_mult_ch, event_weight);
      hmap_mctruth_mucostheta_mumom["total"]->Fill(t->lep_costheta, t->true_muon_mom, event_weight);

      hmap_mctruth_pmom["total"]->Fill(temp_pmom, event_weight);
      hmap_mctruth_pcostheta["total"]->Fill(temp_pangle, event_weight);
      hmap_mctruth_pphi["total"]->Fill(temp_pphi, event_weight);
      hmap_mctruth_thetamup["total"]->Fill(temp_thetamup, event_weight);

      if (t->mode == 0) {
        hmap_mctruth_nuenergy["qe"]->Fill(t->nu_e, event_weight);
        hmap_mctruth_mumom["qe"]->Fill(t->true_muon_mom, event_weight);
        hmap_mctruth_mucostheta["qe"]->Fill(t->lep_costheta, event_weight);
        hmap_mctruth_muphi["qe"]->Fill(t->lep_phi, event_weight);
        hmap_mctruth_chargedmult["qe"]->Fill(t->genie_mult_ch, event_weight);
        hmap_mctruth_mucostheta_mumom["qe"]->Fill(t->lep_costheta, t->true_muon_mom, event_weight);

        hmap_mctruth_pmom["qe"]->Fill(temp_pmom, event_weight);
        hmap_mctruth_pcostheta["qe"]->Fill(temp_pangle, event_weight);
        hmap_mctruth_pphi["qe"]->Fill(temp_pphi, event_weight);
        hmap_mctruth_thetamup["qe"]->Fill(temp_thetamup, event_weight);
      }
      if (t->mode == 1) {
        hmap_mctruth_nuenergy["res"]->Fill(t->nu_e, event_weight);
        hmap_mctruth_mumom["res"]->Fill(t->true_muon_mom, event_weight);
        hmap_mctruth_mucostheta["res"]->Fill(t->lep_costheta, event_weight);
        hmap_mctruth_muphi["res"]->Fill(t->lep_phi, event_weight);
        hmap_mctruth_chargedmult["res"]->Fill(t->genie_mult_ch, event_weight);
        hmap_mctruth_mucostheta_mumom["res"]->Fill(t->lep_costheta, t->true_muon_mom, event_weight);

        hmap_mctruth_pmom["res"]->Fill(temp_pmom, event_weight);
        hmap_mctruth_pcostheta["res"]->Fill(temp_pangle, event_weight);
        hmap_mctruth_pphi["res"]->Fill(temp_pphi, event_weight);
        hmap_mctruth_thetamup["res"]->Fill(temp_thetamup, event_weight);
      }
      if (t->mode == 2) {
        hmap_mctruth_nuenergy["dis"]->Fill(t->nu_e, event_weight);
        hmap_mctruth_mumom["dis"]->Fill(t->true_muon_mom, event_weight);
        hmap_mctruth_mucostheta["dis"]->Fill(t->lep_costheta, event_weight);
        hmap_mctruth_muphi["dis"]->Fill(t->lep_phi, event_weight);
        hmap_mctruth_chargedmult["dis"]->Fill(t->genie_mult_ch, event_weight);
        hmap_mctruth_mucostheta_mumom["dis"]->Fill(t->lep_costheta, t->true_muon_mom, event_weight);

        hmap_mctruth_pmom["dis"]->Fill(temp_pmom, event_weight);
        hmap_mctruth_pcostheta["dis"]->Fill(temp_pangle, event_weight);
        hmap_mctruth_pphi["dis"]->Fill(temp_pphi, event_weight);
        hmap_mctruth_thetamup["dis"]->Fill(temp_thetamup, event_weight);
      }
      if (t->mode == 3) {
        hmap_mctruth_nuenergy["coh"]->Fill(t->nu_e, event_weight);
        hmap_mctruth_mumom["coh"]->Fill(t->true_muon_mom, event_weight);
        hmap_mctruth_mucostheta["coh"]->Fill(t->lep_costheta, event_weight);
        hmap_mctruth_muphi["coh"]->Fill(t->lep_phi, event_weight);
        hmap_mctruth_chargedmult["coh"]->Fill(t->genie_mult_ch, event_weight);
        hmap_mctruth_mucostheta_mumom["coh"]->Fill(t->lep_costheta, t->true_muon_mom, event_weight);

        hmap_mctruth_pmom["coh"]->Fill(temp_pmom, event_weight);
        hmap_mctruth_pcostheta["coh"]->Fill(temp_pangle, event_weight);
        hmap_mctruth_pphi["coh"]->Fill(temp_pphi, event_weight);
        hmap_mctruth_thetamup["coh"]->Fill(temp_thetamup, event_weight);
      }
      if (t->mode == 10) {
        hmap_mctruth_nuenergy["mec"]->Fill(t->nu_e, event_weight);
        hmap_mctruth_mumom["mec"]->Fill(t->true_muon_mom, event_weight);
        hmap_mctruth_mucostheta["mec"]->Fill(t->lep_costheta, event_weight);
        hmap_mctruth_muphi["mec"]->Fill(t->lep_phi, event_weight);
        hmap_mctruth_chargedmult["mec"]->Fill(t->genie_mult_ch, event_weight);
        hmap_mctruth_mucostheta_mumom["mec"]->Fill(t->lep_costheta, t->true_muon_mom, event_weight);

        hmap_mctruth_pmom["mec"]->Fill(temp_pmom, event_weight);
        hmap_mctruth_pcostheta["mec"]->Fill(temp_pangle, event_weight);
        hmap_mctruth_pphi["mec"]->Fill(temp_pphi, event_weight);
        hmap_mctruth_thetamup["mec"]->Fill(temp_thetamup, event_weight);
      }

      pEff->Fill(true, t->nu_e);
      if(muon_contained){
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->pfp_reco_Mom_muon[muind], event_weight, _event_histo_1d->hmap_trkmom_genie_pm1_bs, "signal", fname_genie_pm1, wgts_genie_pm1);
      }else{
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkmom_genie_pm1_bs, "signal", fname_genie_pm1, wgts_genie_pm1);
      }
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->lep_costheta, event_weight, _event_histo_1d->hmap_trktheta_genie_pm1_bs, "signal", fname_genie_pm1, wgts_genie_pm1);
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(temp_pmom, event_weight, _event_histo_1d->hmap_trkpmom_genie_pm1_bs, "signal", fname_genie_pm1, wgts_genie_pm1);
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(temp_pangle, event_weight, _event_histo_1d->hmap_trkptheta_genie_pm1_bs, "signal", fname_genie_pm1, wgts_genie_pm1);

      if (!isdata && _fill_bootstrap_genie) FillBootstrap(thetamup, event_weight, _event_histo_1d->hmap_thetamup_genie_pm1_bs, "signal", fname_genie_pm1, wgts_genie_pm1);



      if (!isdata && _fill_bootstrap_genie) FillBootstrap(0.5, event_weight, _event_histo_1d->hmap_onebin_genie_multisim_bs, "signal", fname_genie_multisim, wgts_genie_multisim);
      if(muon_contained){
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->pfp_reco_Mom_muon[muind], event_weight, _event_histo_1d->hmap_trkmom_genie_multisim_bs, "signal", fname_genie_multisim, wgts_genie_multisim);
      } else {
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkmom_genie_multisim_bs, "signal", fname_genie_multisim, wgts_genie_multisim);
      }
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkangle_genie_multisim_bs, "signal", fname_genie_multisim, wgts_genie_multisim);
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo->hmap_trktheta_trkmom_genie_multisim_bs, "signal", fname_genie_multisim, wgts_genie_multisim);      
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo->hmap_trktheta_trkmom_poly_genie_multisim_bs, "signal", fname_genie_multisim, wgts_genie_multisim);      

      if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(0.5, event_weight, _event_histo_1d->hmap_onebin_extra_syst_multisim_bs, "signal", fname_extra_syst, wgts_extra_syst);
      if(muon_contained){
      if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(t->pfp_reco_Mom_muon[muind], event_weight, _event_histo_1d->hmap_trkmom_extra_syst_multisim_bs, "signal", fname_extra_syst, wgts_extra_syst);
      } else {
      if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkmom_extra_syst_multisim_bs, "signal", fname_extra_syst, wgts_extra_syst);
      }
      if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkangle_extra_syst_multisim_bs, "signal", fname_extra_syst, wgts_extra_syst);
      if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo->hmap_trktheta_trkmom_extra_syst_multisim_bs, "signal", fname_extra_syst, wgts_extra_syst);      
      if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo->hmap_trktheta_trkmom_poly_extra_syst_multisim_bs, "signal", fname_extra_syst, wgts_extra_syst);
      
      if (!isdata && _fill_bootstrap_flux) FillBootstrap(0.5, event_weight, _event_histo_1d->hmap_onebin_flux_multisim_bs, "signal", fname_flux_multisim, wgts_flux_multisim);
      if(muon_contained){
      if (!isdata && _fill_bootstrap_flux) FillBootstrap(t->pfp_reco_Mom_muon[muind], event_weight, _event_histo_1d->hmap_trkmom_flux_multisim_bs, "signal", fname_flux_multisim, wgts_flux_multisim);
      }else {
      if (!isdata && _fill_bootstrap_flux) FillBootstrap(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkmom_flux_multisim_bs, "signal", fname_flux_multisim, wgts_flux_multisim);
      }
      if (!isdata && _fill_bootstrap_flux) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkangle_flux_multisim_bs, "signal", fname_flux_multisim, wgts_flux_multisim);
      if (!isdata && _fill_bootstrap_flux) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo->hmap_trktheta_trkmom_flux_multisim_bs, "signal", fname_flux_multisim, wgts_flux_multisim);
      if (!isdata && _fill_bootstrap_flux) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo->hmap_trktheta_trkmom_poly_flux_multisim_bs, "signal", fname_flux_multisim, wgts_flux_multisim);
 
      if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(0.5, event_weight, _event_histo_1d->hmap_onebin_mc_stat_multisim_bs, "signal", fname_mc_stat_multisim, wgts_mc_stat_multisim);
      if(muon_contained){
      if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(t->pfp_reco_Mom_muon[muind], event_weight, _event_histo_1d->hmap_trkmom_mc_stat_multisim_bs, "signal", fname_mc_stat_multisim, wgts_mc_stat_multisim);
      } else {
      if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkmom_mc_stat_multisim_bs, "signal", fname_mc_stat_multisim, wgts_mc_stat_multisim);
      }
      if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkangle_mc_stat_multisim_bs, "signal", fname_mc_stat_multisim, wgts_mc_stat_multisim);
      if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo->hmap_trktheta_trkmom_mc_stat_multisim_bs, "signal", fname_mc_stat_multisim, wgts_mc_stat_multisim);      
      if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo->hmap_trktheta_trkmom_poly_mc_stat_multisim_bs, "signal", fname_mc_stat_multisim, wgts_mc_stat_multisim);      

      hmap_trkplen["signal"]->Fill(t->pfp_reco_length[pind], event_weight);
      hmap_trkpmom_classic["signal"]->Fill(t->pfp_reco_Mom_proton[pind], event_weight);
      hmap_trkptheta_classic["signal"]->Fill(t->pfp_reco_costheta[pind], event_weight);
      hmap_trkpphi["signal"]->Fill(t->pfp_reco_phi[pind], event_weight);
      hmap_thetamup["signal"]->Fill(thetamup, event_weight);
      hmap_ptmis["signal"]->Fill(ptmis, event_weight);
      if(!isdata){hmap_etatest["signal"]->Fill(etatest, event_weight);}
      hmap_alphat["signal"]->Fill(alphat, event_weight);
      hmap_phit["signal"]->Fill(phit, event_weight);
      hmap_enucal["signal"]->Fill(enucal, event_weight);
      hmap_nhits_leadingp["signal"]->Fill(t->pfp_reco_dEdx[pind].size(), event_weight);
      if(_showerastrack){
        hmap_pmult["signal"]->Fill(t->num_pfp-1, event_weight);
      } else {
        hmap_pmult["signal"]->Fill(t->num_pfp_tracks-1, event_weight);
      }

      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->pfp_reco_Mom_proton[pind], event_weight, _event_histo_1d->hmap_trkpmom_genie_multisim_bs, "signal", fname_genie_multisim, wgts_genie_multisim);      
      if (!isdata && _fill_bootstrap_flux) FillBootstrap(t->pfp_reco_Mom_proton[pind], event_weight, _event_histo_1d->hmap_trkpmom_flux_multisim_bs, "signal", fname_flux_multisim, wgts_flux_multisim);
      if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(t->pfp_reco_Mom_proton[pind], event_weight, _event_histo_1d->hmap_trkpmom_extra_syst_multisim_bs, "signal", fname_extra_syst, wgts_extra_syst);
      if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(t->pfp_reco_Mom_proton[pind], event_weight, _event_histo_1d->hmap_trkpmom_mc_stat_multisim_bs, "signal", fname_mc_stat_multisim, wgts_mc_stat_multisim);
      

      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->pfp_reco_costheta[pind], event_weight, _event_histo_1d->hmap_trkpangle_genie_multisim_bs, "signal", fname_genie_multisim, wgts_genie_multisim);      
      if (!isdata && _fill_bootstrap_flux) FillBootstrap(t->pfp_reco_costheta[pind], event_weight, _event_histo_1d->hmap_trkpangle_flux_multisim_bs, "signal", fname_flux_multisim, wgts_flux_multisim);
      if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(t->pfp_reco_costheta[pind], event_weight, _event_histo_1d->hmap_trkpangle_extra_syst_multisim_bs, "signal", fname_extra_syst, wgts_extra_syst);
      if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(t->pfp_reco_costheta[pind], event_weight, _event_histo_1d->hmap_trkpangle_mc_stat_multisim_bs, "signal", fname_mc_stat_multisim, wgts_mc_stat_multisim);

      if (!isdata && _fill_bootstrap_genie) FillBootstrap(thetamup, event_weight, _event_histo_1d->hmap_thetamup_genie_multisim_bs, "signal", fname_genie_multisim, wgts_genie_multisim);      
      if (!isdata && _fill_bootstrap_flux) FillBootstrap(thetamup, event_weight, _event_histo_1d->hmap_thetamup_flux_multisim_bs, "signal", fname_flux_multisim, wgts_flux_multisim);
      if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(thetamup, event_weight, _event_histo_1d->hmap_thetamup_extra_syst_multisim_bs, "signal", fname_extra_syst, wgts_extra_syst);
      if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(thetamup, event_weight, _event_histo_1d->hmap_thetamup_mc_stat_multisim_bs, "signal", fname_mc_stat_multisim, wgts_mc_stat_multisim);
      if(muon_contained) {
      _event_histo_1d->hmap_trkmom["signal"]->Fill(t->pfp_reco_Mom_muon[muind], event_weight);
      }else {
      _event_histo_1d->hmap_trkmom["signal"]->Fill(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight);
      }
      _event_histo_1d->hmap_trkpmom["signal"]->Fill(t->pfp_reco_Mom_proton[pind], event_weight); 
      _event_histo_1d->hmap_trktheta["signal"]->Fill(t->slc_longesttrack_theta.at(scl_ll_max), event_weight);
      _event_histo_1d->hmap_trkptheta["signal"]->Fill(t->pfp_reco_costheta[pind], event_weight);
      _event_histo_1d->hmap_thetamup["signal"]->Fill(thetamup, event_weight);

      _event_histo->hmap_trktheta_trkmom["signal"]->Fill(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight);
      _event_histo->hmap_trktheta_trkmom_poly["signal"]->Fill(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight);
      //=============================================================================================================================== 
      if (t->slc_origin_extra.at(scl_ll_max) == 0) {
        _event_histo_1d->hmap_onebin["signal_stopmu"]->Fill(0.5, event_weight);
        hmap_trklen["signal_stopmu"]->Fill(t->slc_longesttrack_length.at(scl_ll_max), event_weight);
        _event_histo_1d->hmap_trkmom["signal_stopmu"]->Fill(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight);
        hmap_trkmom_classic["signal_stopmu"]->Fill(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight);
        //hmap_trkmom_genie_pm1_bs["signal_stopmu"]->Fill(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), 1., wgts_genie_pm1);

        hmap_trkphi["signal_stopmu"]->Fill(t->slc_longesttrack_phi.at(scl_ll_max), event_weight);
        _event_histo_1d->hmap_trktheta["signal_stopmu"]->Fill(t->slc_longesttrack_theta.at(scl_ll_max), event_weight);
        hmap_trktheta_classic["signal_stopmu"]->Fill(t->slc_longesttrack_theta.at(scl_ll_max), event_weight);
        hmap_multpfp["signal_stopmu"]->Fill(t->slc_mult_pfp.at(scl_ll_max), event_weight);
        hmap_multtracktol["signal_stopmu"]->Fill(t->slc_mult_track_tolerance.at(scl_ll_max), event_weight);
        _event_histo->hmap_trktheta_trkmom["signal_stopmu"]->Fill(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight);
        _event_histo->hmap_trktheta_trkmom_poly["signal_stopmu"]->Fill(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight);

      } else {
        _event_histo_1d->hmap_onebin["signal_nostopmu"]->Fill(0.5, event_weight);
        hmap_trklen["signal_nostopmu"]->Fill(t->slc_longesttrack_length.at(scl_ll_max), event_weight);
        _event_histo_1d->hmap_trkmom["signal_nostopmu"]->Fill(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight);
        hmap_trkmom_classic["signal_nostopmu"]->Fill(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight);
        //hmap_trkmom_genie_pm1_bs["signal_nostopmu"]->Fill(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), 1., wgts_genie_pm1);

        hmap_trkphi["signal_nostopmu"]->Fill(t->slc_longesttrack_phi.at(scl_ll_max), event_weight);
        _event_histo_1d->hmap_trktheta["signal_nostopmu"]->Fill(t->slc_longesttrack_theta.at(scl_ll_max), event_weight);
        hmap_trktheta_classic["signal_nostopmu"]->Fill(t->slc_longesttrack_theta.at(scl_ll_max), event_weight);
        hmap_multpfp["signal_nostopmu"]->Fill(t->slc_mult_pfp.at(scl_ll_max), event_weight);
        hmap_multtracktol["signal_nostopmu"]->Fill(t->slc_mult_track_tolerance.at(scl_ll_max), event_weight);
        _event_histo->hmap_trktheta_trkmom["signal_nostopmu"]->Fill(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight);
        _event_histo->hmap_trktheta_trkmom_poly["signal_nostopmu"]->Fill(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight);
      }
    }
    //
    // ANUMU
    //
    else if(nu_origin && t->ccnc==0 && t->nupdg==-14 && t->fv==1 && trackfromneutrino){
      bkg_anumu_sel += event_weight;  std::cout<<"ANUMU"<<std::endl;
      pEff->Fill(false, t->nu_e);
      _event_histo_1d->hmap_onebin["anumu"]->Fill(0.5, event_weight);

      hmap_trklen["anumu"]->Fill(t->slc_longesttrack_length.at(scl_ll_max), event_weight);
      if(muon_contained){
      hmap_trkmom_classic["anumu"]->Fill(t->pfp_reco_Mom_muon[muind], event_weight);
      } else {
      hmap_trkmom_classic["anumu"]->Fill(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight);
      }
      hmap_trkphi["anumu"]->Fill(t->slc_longesttrack_phi.at(scl_ll_max), event_weight);
      hmap_trktheta_classic["anumu"]->Fill(t->slc_longesttrack_theta.at(scl_ll_max), event_weight);
      hmap_multpfp["anumu"]->Fill(t->slc_mult_pfp.at(scl_ll_max), event_weight);
      hmap_multtracktol["anumu"]->Fill(t->slc_mult_track_tolerance.at(scl_ll_max), event_weight);

      //hmap_trkmom_genie_pm1_bs["anumu"]->Fill(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), 1., wgts_genie_pm1);
      if(muon_contained){
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->pfp_reco_Mom_muon[muind], event_weight, _event_histo_1d->hmap_trkmom_genie_pm1_bs, "anumu", fname_genie_pm1, wgts_genie_pm1);
      }else{
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkmom_genie_pm1_bs, "anumu", fname_genie_pm1, wgts_genie_pm1);
      }
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(temp_pmom, event_weight, _event_histo_1d->hmap_trkpmom_genie_pm1_bs, "anumu", fname_genie_pm1, wgts_genie_pm1);
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->lep_costheta, event_weight, _event_histo_1d->hmap_trktheta_genie_pm1_bs, "anumu", fname_genie_pm1, wgts_genie_pm1);
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(temp_pangle, event_weight, _event_histo_1d->hmap_trkptheta_genie_pm1_bs, "anumu", fname_genie_pm1, wgts_genie_pm1);

      if (!isdata && _fill_bootstrap_genie) FillBootstrap(thetamup, event_weight, _event_histo_1d->hmap_thetamup_genie_pm1_bs, "anumu", fname_genie_pm1, wgts_genie_pm1);




      if (!isdata && _fill_bootstrap_genie) FillBootstrap(0.5, event_weight, _event_histo_1d->hmap_onebin_genie_multisim_bs, "anumu", fname_genie_multisim, wgts_genie_multisim);
      if(muon_contained){
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->pfp_reco_Mom_muon[muind], event_weight, _event_histo_1d->hmap_trkmom_genie_multisim_bs, "anumu", fname_genie_multisim, wgts_genie_multisim);
      } else {
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkmom_genie_multisim_bs, "anumu", fname_genie_multisim, wgts_genie_multisim);
      }
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkangle_genie_multisim_bs, "anumu", fname_genie_multisim, wgts_genie_multisim);
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo->hmap_trktheta_trkmom_genie_multisim_bs, "anumu", fname_genie_multisim, wgts_genie_multisim);
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo->hmap_trktheta_trkmom_poly_genie_multisim_bs, "anumu", fname_genie_multisim, wgts_genie_multisim);

      if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(0.5, event_weight, _event_histo_1d->hmap_onebin_extra_syst_multisim_bs, "anumu", fname_extra_syst, wgts_extra_syst);
      if(muon_contained){
      if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(t->pfp_reco_Mom_muon[muind], event_weight, _event_histo_1d->hmap_trkmom_extra_syst_multisim_bs, "anumu", fname_extra_syst, wgts_extra_syst);
      } else {
      if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkmom_extra_syst_multisim_bs, "anumu", fname_extra_syst, wgts_extra_syst);
      }
      if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkangle_extra_syst_multisim_bs, "anumu", fname_extra_syst, wgts_extra_syst);
      if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo->hmap_trktheta_trkmom_extra_syst_multisim_bs, "anumu", fname_extra_syst, wgts_extra_syst);
      if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo->hmap_trktheta_trkmom_poly_extra_syst_multisim_bs, "anumu", fname_extra_syst, wgts_extra_syst);
      
      if (!isdata && _fill_bootstrap_flux) FillBootstrap(0.5, event_weight, _event_histo_1d->hmap_onebin_flux_multisim_bs, "anumu", fname_flux_multisim, wgts_flux_multisim);
      if(muon_contained){
      if (!isdata && _fill_bootstrap_flux) FillBootstrap(t->pfp_reco_Mom_muon[muind], event_weight, _event_histo_1d->hmap_trkmom_flux_multisim_bs, "anumu", fname_flux_multisim, wgts_flux_multisim);
      } else{
      if (!isdata && _fill_bootstrap_flux) FillBootstrap(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkmom_flux_multisim_bs, "anumu", fname_flux_multisim, wgts_flux_multisim);
      }
      if (!isdata && _fill_bootstrap_flux) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkangle_flux_multisim_bs, "anumu", fname_flux_multisim, wgts_flux_multisim);
      if (!isdata && _fill_bootstrap_flux) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo->hmap_trktheta_trkmom_flux_multisim_bs, "anumu", fname_flux_multisim, wgts_flux_multisim);
      if (!isdata && _fill_bootstrap_flux) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo->hmap_trktheta_trkmom_poly_flux_multisim_bs, "anumu", fname_flux_multisim, wgts_flux_multisim);

      if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(0.5, event_weight, _event_histo_1d->hmap_onebin_mc_stat_multisim_bs, "anumu", fname_mc_stat_multisim, wgts_mc_stat_multisim);
      if(muon_contained){
      if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(t->pfp_reco_Mom_muon[muind], event_weight, _event_histo_1d->hmap_trkmom_mc_stat_multisim_bs, "anumu", fname_mc_stat_multisim, wgts_mc_stat_multisim);
      } else {
      if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkmom_mc_stat_multisim_bs, "anumu", fname_mc_stat_multisim, wgts_mc_stat_multisim);
      }
      if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkangle_mc_stat_multisim_bs, "anumu", fname_mc_stat_multisim, wgts_mc_stat_multisim);
      if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo->hmap_trktheta_trkmom_mc_stat_multisim_bs, "anumu", fname_mc_stat_multisim, wgts_mc_stat_multisim);
      if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo->hmap_trktheta_trkmom_mc_stat_multisim_bs, "anumu", fname_mc_stat_multisim, wgts_mc_stat_multisim);
      hmap_trkplen["anumu"]->Fill(t->pfp_reco_length[pind], event_weight);
      hmap_trkpmom_classic["anumu"]->Fill(t->pfp_reco_Mom_proton[pind], event_weight);
      hmap_trkptheta_classic["anumu"]->Fill(t->pfp_reco_costheta[pind], event_weight);
      hmap_trkpphi["anumu"]->Fill(t->pfp_reco_phi[pind], event_weight);
      hmap_thetamup["anumu"]->Fill(thetamup, event_weight);
      hmap_ptmis["anumu"]->Fill(ptmis, event_weight);
      if(!isdata){hmap_etatest["anumu"]->Fill(etatest, event_weight);}
      hmap_alphat["anumu"]->Fill(alphat, event_weight);
      hmap_phit["anumu"]->Fill(phit, event_weight);
      hmap_enucal["anumu"]->Fill(enucal, event_weight);
      hmap_nhits_leadingp["anumu"]->Fill(t->pfp_reco_dEdx[pind].size(), event_weight);
      if(_showerastrack){
        hmap_pmult["anumu"]->Fill(t->num_pfp-1, event_weight);
      } else {
        hmap_pmult["anumu"]->Fill(t->num_pfp_tracks-1, event_weight);
      }

      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->pfp_reco_Mom_proton[pind], event_weight, _event_histo_1d->hmap_trkpmom_genie_multisim_bs, "anumu", fname_genie_multisim, wgts_genie_multisim);      
      if (!isdata && _fill_bootstrap_flux) FillBootstrap(t->pfp_reco_Mom_proton[pind], event_weight, _event_histo_1d->hmap_trkpmom_flux_multisim_bs, "anumu", fname_flux_multisim, wgts_flux_multisim);
      if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(t->pfp_reco_Mom_proton[pind], event_weight, _event_histo_1d->hmap_trkpmom_extra_syst_multisim_bs, "anumu", fname_extra_syst, wgts_extra_syst);
      if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(t->pfp_reco_Mom_proton[pind], event_weight, _event_histo_1d->hmap_trkpmom_mc_stat_multisim_bs, "anumu", fname_mc_stat_multisim, wgts_mc_stat_multisim);

      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->pfp_reco_costheta[pind], event_weight, _event_histo_1d->hmap_trkpangle_genie_multisim_bs, "anumu", fname_genie_multisim, wgts_genie_multisim);      
      if (!isdata && _fill_bootstrap_flux) FillBootstrap(t->pfp_reco_costheta[pind], event_weight, _event_histo_1d->hmap_trkpangle_flux_multisim_bs, "anumu", fname_flux_multisim, wgts_flux_multisim);
      if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(t->pfp_reco_costheta[pind], event_weight, _event_histo_1d->hmap_trkpangle_extra_syst_multisim_bs, "anumu", fname_extra_syst, wgts_extra_syst);
      if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(t->pfp_reco_costheta[pind], event_weight, _event_histo_1d->hmap_trkpangle_mc_stat_multisim_bs, "anumu", fname_mc_stat_multisim, wgts_mc_stat_multisim);

      if (!isdata && _fill_bootstrap_genie) FillBootstrap(thetamup, event_weight, _event_histo_1d->hmap_thetamup_genie_multisim_bs, "anumu", fname_genie_multisim, wgts_genie_multisim);      
      if (!isdata && _fill_bootstrap_flux) FillBootstrap(thetamup, event_weight, _event_histo_1d->hmap_thetamup_flux_multisim_bs, "anumu", fname_flux_multisim, wgts_flux_multisim);
      if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(thetamup, event_weight, _event_histo_1d->hmap_thetamup_extra_syst_multisim_bs, "anumu", fname_extra_syst, wgts_extra_syst);
      if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(thetamup, event_weight, _event_histo_1d->hmap_thetamup_mc_stat_multisim_bs, "anumu", fname_mc_stat_multisim, wgts_mc_stat_multisim);

      if(muon_contained){
      _event_histo_1d->hmap_trkmom["anumu"]->Fill(t->pfp_reco_Mom_muon[muind], event_weight);

      }else{
      _event_histo_1d->hmap_trkmom["anumu"]->Fill(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight);
      }
      _event_histo_1d->hmap_trkpmom["anumu"]->Fill(t->pfp_reco_Mom_proton[pind], event_weight);
      _event_histo_1d->hmap_trktheta["anumu"]->Fill(t->slc_longesttrack_theta.at(scl_ll_max), event_weight);
      _event_histo_1d->hmap_trkptheta["anumu"]->Fill(t->pfp_reco_costheta[pind], event_weight);
      _event_histo_1d->hmap_thetamup["anumu"]->Fill(thetamup, event_weight);

      _event_histo->hmap_trktheta_trkmom["anumu"]->Fill(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight);
      _event_histo->hmap_trktheta_trkmom_poly["anumu"]->Fill(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight);
    }
    //======================================================================================================
    //
    // CC0P or CCNpi
    //
    else if(nu_origin && t->ccnc==0 && t->nupdg==14 && t->fv==1 && trackfromneutrino && (t->ngenie_protons_300==0 ||(t->ngenie_pipms+t->ngenie_pion0s)>0 || t->ngenie_electrons>0)){
      bkg_ccother_sel += event_weight; std::cout<<"CCOther"<<std::endl;
      _event_histo_1d->hmap_onebin["cc_other"]->Fill(0.5, event_weight);
      hmap_trklen["cc_other"]->Fill(t->slc_longesttrack_length.at(scl_ll_max), event_weight);
      if(muon_contained){
      hmap_trkmom_classic["cc_other"]->Fill(t->pfp_reco_Mom_muon[muind], event_weight);
      } else {
      hmap_trkmom_classic["cc_other"]->Fill(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight);
      }
      hmap_trkphi["cc_other"]->Fill(t->slc_longesttrack_phi.at(scl_ll_max), event_weight);
      hmap_trktheta_classic["cc_other"]->Fill(t->slc_longesttrack_theta.at(scl_ll_max), event_weight);
      hmap_multpfp["cc_other"]->Fill(t->slc_mult_pfp.at(scl_ll_max), event_weight);
      hmap_multtracktol["cc_other"]->Fill(t->slc_mult_track_tolerance.at(scl_ll_max), event_weight);
       //hmap_trkmom_genie_pm1_bs["ccother"]->Fill(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), 1., wgts_genie_pm1);
      if(muon_contained){
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->pfp_reco_Mom_muon[muind], event_weight, _event_histo_1d->hmap_trkmom_genie_pm1_bs, "cc_other", fname_genie_pm1, wgts_genie_pm1);
      }else{   
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkmom_genie_pm1_bs, "cc_other", fname_genie_pm1, wgts_genie_pm1);
      }
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(temp_pmom, event_weight, _event_histo_1d->hmap_trkpmom_genie_pm1_bs, "cc_other", fname_genie_pm1, wgts_genie_pm1);
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->lep_costheta, event_weight, _event_histo_1d->hmap_trktheta_genie_pm1_bs, "cc_other", fname_genie_pm1, wgts_genie_pm1);
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(temp_pangle, event_weight, _event_histo_1d->hmap_trkptheta_genie_pm1_bs, "cc_other", fname_genie_pm1, wgts_genie_pm1);

      if (!isdata && _fill_bootstrap_genie) FillBootstrap(thetamup, event_weight, _event_histo_1d->hmap_thetamup_genie_pm1_bs, "cc_other", fname_genie_pm1, wgts_genie_pm1);



      if (!isdata && _fill_bootstrap_genie) FillBootstrap(0.5, event_weight, _event_histo_1d->hmap_onebin_genie_multisim_bs, "cc_other", fname_genie_multisim, wgts_genie_multisim);
      if(muon_contained){
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->pfp_reco_Mom_muon[muind], event_weight, _event_histo_1d->hmap_trkmom_genie_multisim_bs, "cc_other", fname_genie_multisim, wgts_genie_multisim);
      } else {
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkmom_genie_multisim_bs, "cc_other", fname_genie_multisim, wgts_genie_multisim);
      }  
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkangle_genie_multisim_bs, "cc_other", fname_genie_multisim, wgts_genie_multisim);
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo->hmap_trktheta_trkmom_genie_multisim_bs, "cc_other", fname_genie_multisim, wgts_genie_multisim);
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo->hmap_trktheta_trkmom_poly_genie_multisim_bs, "cc_other", fname_genie_multisim, wgts_genie_multisim);

      if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(0.5, event_weight, _event_histo_1d->hmap_onebin_extra_syst_multisim_bs, "cc_other", fname_extra_syst, wgts_extra_syst);
      if(muon_contained){
      if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(t->pfp_reco_Mom_muon[muind], event_weight, _event_histo_1d->hmap_trkmom_extra_syst_multisim_bs, "cc_other", fname_extra_syst, wgts_extra_syst);
      } else{
      if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkmom_extra_syst_multisim_bs, "cc_other", fname_extra_syst, wgts_extra_syst);
      }
      if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkangle_extra_syst_multisim_bs, "cc_other", fname_extra_syst, wgts_extra_syst);
      if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo->hmap_trktheta_trkmom_extra_syst_multisim_bs, "cc_other", fname_extra_syst, wgts_extra_syst);
      if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo->hmap_trktheta_trkmom_poly_extra_syst_multisim_bs, "cc_other", fname_extra_syst, wgts_extra_syst);
      
      if (!isdata && _fill_bootstrap_flux) FillBootstrap(0.5, event_weight, _event_histo_1d->hmap_onebin_flux_multisim_bs, "cc_other", fname_flux_multisim, wgts_flux_multisim);
      if(muon_contained){
      if (!isdata && _fill_bootstrap_flux) FillBootstrap(t->pfp_reco_Mom_muon[muind], event_weight, _event_histo_1d->hmap_trkmom_flux_multisim_bs, "cc_other", fname_flux_multisim, wgts_flux_multisim);
      }else {  
      if (!isdata && _fill_bootstrap_flux) FillBootstrap(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkmom_flux_multisim_bs, "cc_other", fname_flux_multisim, wgts_flux_multisim);
      }
      if (!isdata && _fill_bootstrap_flux) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkangle_flux_multisim_bs, "cc_other", fname_flux_multisim, wgts_flux_multisim);
      if (!isdata && _fill_bootstrap_flux) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo->hmap_trktheta_trkmom_flux_multisim_bs, "cc_other", fname_flux_multisim, wgts_flux_multisim);
      if (!isdata && _fill_bootstrap_flux) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo->hmap_trktheta_trkmom_poly_flux_multisim_bs, "cc_other", fname_flux_multisim, wgts_flux_multisim);

      if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(0.5, event_weight, _event_histo_1d->hmap_onebin_mc_stat_multisim_bs, "cc_other", fname_mc_stat_multisim, wgts_mc_stat_multisim);
      if(muon_contained){
      if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(t->pfp_reco_Mom_muon[muind], event_weight, _event_histo_1d->hmap_trkmom_mc_stat_multisim_bs, "cc_other", fname_mc_stat_multisim, wgts_mc_stat_multisim);
      } else {
      if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkmom_mc_stat_multisim_bs, "cc_other", fname_mc_stat_multisim, wgts_mc_stat_multisim);
      }
      if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkangle_mc_stat_multisim_bs, "cc_other", fname_mc_stat_multisim, wgts_mc_stat_multisim);
      if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo->hmap_trktheta_trkmom_mc_stat_multisim_bs, "cc_other", fname_mc_stat_multisim, wgts_mc_stat_multisim);
      if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo->hmap_trktheta_trkmom_poly_mc_stat_multisim_bs, "cc_other", fname_mc_stat_multisim, wgts_mc_stat_multisim);
      hmap_trkplen["cc_other"]->Fill(t->pfp_reco_length[pind], event_weight);
      hmap_trkpmom_classic["cc_other"]->Fill(t->pfp_reco_Mom_proton[pind], event_weight);
      hmap_trkptheta_classic["cc_other"]->Fill(t->pfp_reco_costheta[pind], event_weight);
      hmap_trkpphi["cc_other"]->Fill(t->pfp_reco_phi[pind], event_weight);
      hmap_thetamup["cc_other"]->Fill(thetamup, event_weight);
      hmap_ptmis["cc_other"]->Fill(ptmis, event_weight);
      if(!isdata){hmap_etatest["cc_other"]->Fill(etatest, event_weight);}
      hmap_alphat["cc_other"]->Fill(alphat, event_weight);
      hmap_phit["cc_other"]->Fill(phit, event_weight);
      hmap_enucal["cc_other"]->Fill(enucal, event_weight);
      hmap_nhits_leadingp["cc_other"]->Fill(t->pfp_reco_dEdx[pind].size(), event_weight);
      if(_showerastrack){ 
        hmap_pmult["cc_other"]->Fill(t->num_pfp-1, event_weight);
      } else{  
        hmap_pmult["cc_other"]->Fill(t->num_pfp_tracks-1, event_weight);
      }
      //std::cout<<"libo test ggg"<<std::endl;
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->pfp_reco_Mom_proton[pind], event_weight, _event_histo_1d->hmap_trkpmom_genie_multisim_bs, "cc_other", fname_genie_multisim, wgts_genie_multisim);      
      if (!isdata && _fill_bootstrap_flux) FillBootstrap(t->pfp_reco_Mom_proton[pind], event_weight, _event_histo_1d->hmap_trkpmom_flux_multisim_bs, "cc_other", fname_flux_multisim, wgts_flux_multisim);
      if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(t->pfp_reco_Mom_proton[pind], event_weight, _event_histo_1d->hmap_trkpmom_extra_syst_multisim_bs, "cc_other", fname_extra_syst, wgts_extra_syst);
      if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(t->pfp_reco_Mom_proton[pind], event_weight, _event_histo_1d->hmap_trkpmom_mc_stat_multisim_bs, "cc_other", fname_mc_stat_multisim, wgts_mc_stat_multisim);

      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->pfp_reco_costheta[pind], event_weight, _event_histo_1d->hmap_trkpangle_genie_multisim_bs, "cc_other", fname_genie_multisim, wgts_genie_multisim);      
      if (!isdata && _fill_bootstrap_flux) FillBootstrap(t->pfp_reco_costheta[pind], event_weight, _event_histo_1d->hmap_trkpangle_flux_multisim_bs, "cc_other", fname_flux_multisim, wgts_flux_multisim);
      if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(t->pfp_reco_costheta[pind], event_weight, _event_histo_1d->hmap_trkpangle_extra_syst_multisim_bs, "cc_other", fname_extra_syst, wgts_extra_syst);
      if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(t->pfp_reco_costheta[pind], event_weight, _event_histo_1d->hmap_trkpangle_mc_stat_multisim_bs, "cc_other", fname_mc_stat_multisim, wgts_mc_stat_multisim);

      if (!isdata && _fill_bootstrap_genie) FillBootstrap(thetamup, event_weight, _event_histo_1d->hmap_thetamup_genie_multisim_bs, "cc_other", fname_genie_multisim, wgts_genie_multisim);      
      if (!isdata && _fill_bootstrap_flux) FillBootstrap(thetamup, event_weight, _event_histo_1d->hmap_thetamup_flux_multisim_bs, "cc_other", fname_flux_multisim, wgts_flux_multisim);
      if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(thetamup, event_weight, _event_histo_1d->hmap_thetamup_extra_syst_multisim_bs, "cc_other", fname_extra_syst, wgts_extra_syst);
      if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(thetamup, event_weight, _event_histo_1d->hmap_thetamup_mc_stat_multisim_bs, "cc_other", fname_mc_stat_multisim, wgts_mc_stat_multisim);
      //std::cout<<"libo test jjj"<<std::endl;
      if(muon_contained){
      _event_histo_1d->hmap_trkmom["cc_other"]->Fill(t->pfp_reco_Mom_muon[muind], event_weight);
      } else{
      _event_histo_1d->hmap_trkmom["cc_other"]->Fill(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight);
      }
      _event_histo_1d->hmap_trkpmom["cc_other"]->Fill(t->pfp_reco_Mom_proton[pind], event_weight); 
      _event_histo_1d->hmap_trktheta["cc_other"]->Fill(t->slc_longesttrack_theta.at(scl_ll_max), event_weight);
      _event_histo_1d->hmap_trkptheta["cc_other"]->Fill(t->slc_longesttrack_theta.at(scl_ll_max), event_weight);
      _event_histo_1d->hmap_thetamup["cc_other"]->Fill(thetamup, event_weight);
      //std::cout<<"libo test kkk"<<std::endl;
      _event_histo->hmap_trktheta_trkmom["cc_other"]->Fill(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight);
      _event_histo->hmap_trktheta_trkmom_poly["cc_other"]->Fill(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight);
      //std::cout<<"libo test 0"<<std::endl;//characterize the cc background
      if(t->ngenie_protons_300==0 && (t->ngenie_pipms+t->ngenie_pion0s)==0&& t->ngenie_electrons==0) {  //CC0P0Pi
        hmap_trklen["cc_0proton"]->Fill(t->pfp_reco_length[muind], event_weight);
        hmap_trkmom_classic["cc_0proton"]->Fill(t->pfp_reco_Mom_MCS[muind], event_weight);
        hmap_trktheta_classic["cc_0proton"]->Fill(t->pfp_reco_costheta[pind], event_weight);
        hmap_trkphi["cc_0proton"]->Fill(t->pfp_reco_phi[pind], event_weight);
        hmap_trkplen["cc_0proton"]->Fill(t->pfp_reco_length[pind], event_weight);
        hmap_trkpmom_classic["cc_0proton"]->Fill(t->pfp_reco_Mom_proton[pind], event_weight);
        hmap_trkptheta_classic["cc_0proton"]->Fill(t->pfp_reco_costheta[pind], event_weight);
        hmap_trkpphi["cc_0proton"]->Fill(t->pfp_reco_phi[pind], event_weight);
        
        h_ngenie_proton->Fill(t->ngenie_protons, event_weight);
        
        int munum=0;
        for(unsigned int ii=0; ii<t->pfp_reco_ismuoncandidate.size(); ii++){
                if(t->pfp_truth_pdg[ii]==13) munum++;
        } 
        int pnum_03_new=0;
        for(unsigned int jj=0; jj<t->genie_mcpar_pdgcode.size(); jj++){
           if(abs(t->genie_mcpar_pdgcode[jj]) !=2212) continue;
           h_true_pmom_cc0p->Fill(TMath::Sqrt(t->genie_mcpar_px[jj]*t->genie_mcpar_px[jj]+t->genie_mcpar_py[jj]*t->genie_mcpar_py[jj]+t->genie_mcpar_pz[jj]*t->genie_mcpar_pz[jj]), event_weight);
           if(TMath::Sqrt(t->genie_mcpar_energy[jj]*t->genie_mcpar_energy[jj]-0.938*0.938)>0.3){
             pnum_03_new++; 
           } 
        }
        if(pnum_03_new>0) {
            std::cout<<"run number is "<<t->run<<std::endl;
            std::cout<<"subrun number is "<<t->subrun<<std::endl;
            std::cout<<"event number is "<<t->event<<std::endl;
            for(unsigned int ll=0; ll<t->genie_mcpar_pdgcode.size(); ll++){
               if(abs(t->genie_mcpar_pdgcode[ll]) !=2212) continue;
               std::cout<<"ll = "<<ll<<" PDGcode is "<<t->genie_mcpar_pdgcode[ll]<<" Momentum is "<<TMath::Sqrt(t->genie_mcpar_energy[ll]*t->genie_mcpar_energy[ll]-0.938*0.938)<<std::endl;
            }  
        }
        if(pnum_03_new==0 && munum>=2){
          if(t->pfp_truth_pdg[muind]==13 && t->pfp_truth_pdg[pind]==13){
            h_true_thetamup_cc0p->Fill(getAngle(t->pfp_truth_mom[muind], t->pfp_truth_theta[muind], t->pfp_truth_phi[muind], t->pfp_truth_mom[pind], t->pfp_truth_theta[pind], t->pfp_truth_phi[pind]), event_weight);
            h_reco_thetamup_cc0p->Fill(thetamup, event_weight);
          } 
          //else if(t->pfp_truth_pdg[muind]==13 && t->pfp_truth_pdg[pind]==22){
            //no this kind of events
            //h_true_photonmom_cc0p->Fill(t->pfp_truth_mom[pind], event_weight);
            //h_reco_photonmom_cc0p->Fill(t->pfp_reco_Mom_proton[pind], event_weight);
          //}
        } else if(pnum_03_new==0 && t->pfp_truth_pdg[muind]==2212 && t->pfp_truth_pdg[pind]==22){
            h_true_photonmom_cc0p->Fill(t->pfp_truth_mom[pind], event_weight);
            h_reco_photonmom_cc0p->Fill(t->pfp_reco_Mom_proton[pind], event_weight);
         
        } else if(pnum_03_new==0 && t->pfp_truth_pdg[muind]==13 && t->pfp_truth_pdg[pind]==2212){
            h_true_reco_protonmom_cc0p->Fill(t->pfp_truth_mom[pind], t->pfp_reco_Mom_proton[pind], event_weight);
            h_true_reco_protonlen_cc0p->Fill(TMath::Sqrt((t->pfp_truth_startx[pind]-t->pfp_truth_endx[pind])*(t->pfp_truth_startx[pind]-t->pfp_truth_endx[pind])
                                                         +(t->pfp_truth_starty[pind]-t->pfp_truth_endy[pind])*(t->pfp_truth_starty[pind]-t->pfp_truth_endy[pind])
                                                         +(t->pfp_truth_startz[pind]-t->pfp_truth_endz[pind])*(t->pfp_truth_startz[pind]-t->pfp_truth_endz[pind])),t->pfp_reco_length[pind], event_weight);
        }



      }
      if(t->ngenie_pipms+t->ngenie_pion0s>0 || t->ngenie_electrons>0){
        hmap_trklen["cc_pion"]->Fill(t->pfp_reco_length[muind], event_weight);
        hmap_trkmom_classic["cc_pion"]->Fill(t->pfp_reco_Mom_MCS[muind], event_weight);
        hmap_trktheta_classic["cc_pion"]->Fill(t->pfp_reco_costheta[pind], event_weight);
        hmap_trkphi["cc_pion"]->Fill(t->pfp_reco_phi[pind], event_weight);
        hmap_trkplen["cc_pion"]->Fill(t->pfp_reco_length[pind], event_weight);
        hmap_trkpmom_classic["cc_pion"]->Fill(t->pfp_reco_Mom_proton[pind], event_weight);
        hmap_trkptheta_classic["cc_pion"]->Fill(t->pfp_reco_costheta[pind], event_weight);
        hmap_trkpphi["cc_pion"]->Fill(t->pfp_reco_phi[pind], event_weight);
        if(t->ngenie_pipms>0 && t->ngenie_electrons==0){
           h_pion_reco->Fill(pion_reco, event_weight);   
           h_pion_reint->Fill(pion_reint, event_weight);
           //check the momentum of the pions that are not reconstructed
           if(pion_reco==false){ //check the pions not reconstructed interacted or not
              for(unsigned int kk=0; kk<t->genie_mcpar_pdgcode.size(); kk++){
                  if(abs(t->genie_mcpar_pdgcode[kk])==211) {
                      h_true_pimom_noreco->Fill(TMath::Sqrt(t->genie_mcpar_px[kk]*t->genie_mcpar_px[kk]+t->genie_mcpar_py[kk]*t->genie_mcpar_py[kk]+t->genie_mcpar_pz[kk]*t->genie_mcpar_pz[kk]), event_weight);
                      h_true_pilen_noreco->Fill(TMath::Sqrt((t->pfp_truth_startx[kk]-t->pfp_truth_endx[kk])*(t->pfp_truth_startx[kk]-t->pfp_truth_endx[kk])
                                                           +(t->pfp_truth_starty[kk]-t->pfp_truth_endy[kk])*(t->pfp_truth_starty[kk]-t->pfp_truth_endy[kk])
                                                           +(t->pfp_truth_startz[kk]-t->pfp_truth_endz[kk])*(t->pfp_truth_startz[kk]-t->pfp_truth_endz[kk])), event_weight);
                      h_true_pizlen_noreco->Fill(TMath::Sqrt(t->pfp_truth_startz[kk]-t->pfp_truth_endz[kk])*(t->pfp_truth_startz[kk]-t->pfp_truth_endz[kk]), event_weight);
                  }
              }
           }  
        }//end of there are charged pions produced
        //if(t->ngenie_pion0s>0){
           //std::cout<<"background with pi0 produced, the pdg code of leading proton candidate is : "<<t->pfp_truth_pdg[pind]<<" is it primary? "<<t->pfp_reco_isprimary[pind]<<std::endl;
        //}
      }
      //std::cout<<"libo test at the end of other cc background"<<std::endl;
    }
    //=======================================================================================================
    //
    // NUE
    //
    else if(nu_origin && t->ccnc==0 && (t->nupdg==-12 || t->nupdg==12) && t->fv==1 && trackfromneutrino){
      bkg_nue_sel += event_weight;  std::cout<<"NUE"<<std::endl;
      pEff->Fill(false, t->nu_e);
      _event_histo_1d->hmap_onebin["nue"]->Fill(0.5, event_weight);
      hmap_trklen["nue"]->Fill(t->slc_longesttrack_length.at(scl_ll_max), event_weight);
      if(muon_contained){
      hmap_trkmom_classic["nue"]->Fill(t->pfp_reco_Mom_muon[muind], event_weight);
      } else {
      hmap_trkmom_classic["nue"]->Fill(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight);
      }
      hmap_trkphi["nue"]->Fill(t->slc_longesttrack_phi.at(scl_ll_max), event_weight);
      hmap_trktheta_classic["nue"]->Fill(t->slc_longesttrack_theta.at(scl_ll_max), event_weight);
      hmap_multpfp["nue"]->Fill(t->slc_mult_pfp.at(scl_ll_max), event_weight);
      hmap_multtracktol["nue"]->Fill(t->slc_mult_track_tolerance.at(scl_ll_max), event_weight);
       //hmap_trkmom_genie_pm1_bs["nue"]->Fill(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), 1., wgts_genie_pm1);
      if(muon_contained){
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->pfp_reco_Mom_muon[muind], event_weight, _event_histo_1d->hmap_trkmom_genie_pm1_bs, "nue", fname_genie_pm1, wgts_genie_pm1);
      }else{
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkmom_genie_pm1_bs, "nue", fname_genie_pm1, wgts_genie_pm1);
      }
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(temp_pmom, event_weight, _event_histo_1d->hmap_trkpmom_genie_pm1_bs, "nue", fname_genie_pm1, wgts_genie_pm1);
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->lep_costheta, event_weight, _event_histo_1d->hmap_trktheta_genie_pm1_bs, "nue", fname_genie_pm1, wgts_genie_pm1);
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(temp_pangle, event_weight, _event_histo_1d->hmap_trkptheta_genie_pm1_bs, "nue", fname_genie_pm1, wgts_genie_pm1);

      if (!isdata && _fill_bootstrap_genie) FillBootstrap(thetamup, event_weight, _event_histo_1d->hmap_thetamup_genie_pm1_bs, "nue", fname_genie_pm1, wgts_genie_pm1);




      if (!isdata && _fill_bootstrap_genie) FillBootstrap(0.5, event_weight, _event_histo_1d->hmap_onebin_genie_multisim_bs, "nue", fname_genie_multisim, wgts_genie_multisim);
      if(muon_contained){
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->pfp_reco_Mom_muon[muind], event_weight, _event_histo_1d->hmap_trkmom_genie_multisim_bs, "nue", fname_genie_multisim, wgts_genie_multisim);
      } else{
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkmom_genie_multisim_bs, "nue", fname_genie_multisim, wgts_genie_multisim);
      } 
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkangle_genie_multisim_bs, "nue", fname_genie_multisim, wgts_genie_multisim);
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo->hmap_trktheta_trkmom_genie_multisim_bs, "nue", fname_genie_multisim, wgts_genie_multisim);
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo->hmap_trktheta_trkmom_poly_genie_multisim_bs, "nue", fname_genie_multisim, wgts_genie_multisim);

      if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(0.5, event_weight, _event_histo_1d->hmap_onebin_extra_syst_multisim_bs, "nue", fname_extra_syst, wgts_extra_syst);
      if(muon_contained){
      if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(t->pfp_reco_Mom_muon[muind], event_weight, _event_histo_1d->hmap_trkmom_extra_syst_multisim_bs, "nue", fname_extra_syst, wgts_extra_syst);
      } else{
      if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkmom_extra_syst_multisim_bs, "nue", fname_extra_syst, wgts_extra_syst);
      }
      if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkangle_extra_syst_multisim_bs, "nue", fname_extra_syst, wgts_extra_syst);
      if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo->hmap_trktheta_trkmom_extra_syst_multisim_bs, "nue", fname_extra_syst, wgts_extra_syst);
      if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo->hmap_trktheta_trkmom_poly_extra_syst_multisim_bs, "nue", fname_extra_syst, wgts_extra_syst);
      
      if (!isdata && _fill_bootstrap_flux) FillBootstrap(0.5, event_weight, _event_histo_1d->hmap_onebin_flux_multisim_bs, "nue", fname_flux_multisim, wgts_flux_multisim);
      if(muon_contained){
      if (!isdata && _fill_bootstrap_flux) FillBootstrap(t->pfp_reco_Mom_muon[muind], event_weight, _event_histo_1d->hmap_trkmom_flux_multisim_bs, "nue", fname_flux_multisim, wgts_flux_multisim);
      } else{
      if (!isdata && _fill_bootstrap_flux) FillBootstrap(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkmom_flux_multisim_bs, "nue", fname_flux_multisim, wgts_flux_multisim);
      }
      if (!isdata && _fill_bootstrap_flux) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkangle_flux_multisim_bs, "nue", fname_flux_multisim, wgts_flux_multisim);
      if (!isdata && _fill_bootstrap_flux) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo->hmap_trktheta_trkmom_flux_multisim_bs, "nue", fname_flux_multisim, wgts_flux_multisim);
      if (!isdata && _fill_bootstrap_flux) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo->hmap_trktheta_trkmom_poly_flux_multisim_bs, "nue", fname_flux_multisim, wgts_flux_multisim);

      if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(0.5, event_weight, _event_histo_1d->hmap_onebin_mc_stat_multisim_bs, "nue", fname_mc_stat_multisim, wgts_mc_stat_multisim);
      if(muon_contained){
      if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(t->pfp_reco_Mom_muon[muind], event_weight, _event_histo_1d->hmap_trkmom_mc_stat_multisim_bs, "nue", fname_mc_stat_multisim, wgts_mc_stat_multisim);
      } else{
      if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkmom_mc_stat_multisim_bs, "nue", fname_mc_stat_multisim, wgts_mc_stat_multisim);
      }
      if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkangle_mc_stat_multisim_bs, "nue", fname_mc_stat_multisim, wgts_mc_stat_multisim);
      if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo->hmap_trktheta_trkmom_mc_stat_multisim_bs, "nue", fname_mc_stat_multisim, wgts_mc_stat_multisim);
      if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo->hmap_trktheta_trkmom_poly_mc_stat_multisim_bs, "nue", fname_mc_stat_multisim, wgts_mc_stat_multisim);
      hmap_trkplen["nue"]->Fill(t->pfp_reco_length[pind], event_weight);
      hmap_trkpmom_classic["nue"]->Fill(t->pfp_reco_Mom_proton[pind], event_weight);
      hmap_trkptheta_classic["nue"]->Fill(t->pfp_reco_costheta[pind], event_weight);
      hmap_trkpphi["nue"]->Fill(t->pfp_reco_phi[pind], event_weight);
      hmap_thetamup["nue"]->Fill(thetamup, event_weight);
      hmap_ptmis["nue"]->Fill(ptmis, event_weight);
      if(!isdata){hmap_etatest["nue"]->Fill(etatest, event_weight);}
      hmap_alphat["nue"]->Fill(alphat, event_weight);
      hmap_phit["nue"]->Fill(phit, event_weight);
      hmap_enucal["nue"]->Fill(enucal, event_weight); 
      hmap_nhits_leadingp["nue"]->Fill(t->pfp_reco_dEdx[pind].size(), event_weight);
      if(_showerastrack){
        hmap_pmult["nue"]->Fill(t->num_pfp-1, event_weight);
      } else {
        hmap_pmult["nue"]->Fill(t->num_pfp_tracks-1, event_weight);
      }

      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->pfp_reco_Mom_proton[pind], event_weight, _event_histo_1d->hmap_trkpmom_genie_multisim_bs, "nue", fname_genie_multisim, wgts_genie_multisim);      
      if (!isdata && _fill_bootstrap_flux) FillBootstrap(t->pfp_reco_Mom_proton[pind], event_weight, _event_histo_1d->hmap_trkpmom_flux_multisim_bs, "nue", fname_flux_multisim, wgts_flux_multisim);
      if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(t->pfp_reco_Mom_proton[pind], event_weight, _event_histo_1d->hmap_trkpmom_extra_syst_multisim_bs, "nue", fname_extra_syst, wgts_extra_syst);
      if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(t->pfp_reco_Mom_proton[pind], event_weight, _event_histo_1d->hmap_trkpmom_mc_stat_multisim_bs, "nue", fname_mc_stat_multisim, wgts_mc_stat_multisim);

      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->pfp_reco_costheta[pind], event_weight, _event_histo_1d->hmap_trkpangle_genie_multisim_bs, "nue", fname_genie_multisim, wgts_genie_multisim);      
      if (!isdata && _fill_bootstrap_flux) FillBootstrap(t->pfp_reco_costheta[pind], event_weight, _event_histo_1d->hmap_trkpangle_flux_multisim_bs, "nue", fname_flux_multisim, wgts_flux_multisim);
      if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(t->pfp_reco_costheta[pind], event_weight, _event_histo_1d->hmap_trkpangle_extra_syst_multisim_bs, "nue", fname_extra_syst, wgts_extra_syst);
      if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(t->pfp_reco_costheta[pind], event_weight, _event_histo_1d->hmap_trkpangle_mc_stat_multisim_bs, "nue", fname_mc_stat_multisim, wgts_mc_stat_multisim);

      if (!isdata && _fill_bootstrap_genie) FillBootstrap(thetamup, event_weight, _event_histo_1d->hmap_thetamup_genie_multisim_bs, "nue", fname_genie_multisim, wgts_genie_multisim);      
      if (!isdata && _fill_bootstrap_flux) FillBootstrap(thetamup, event_weight, _event_histo_1d->hmap_thetamup_flux_multisim_bs, "nue", fname_flux_multisim, wgts_flux_multisim);
      if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(thetamup, event_weight, _event_histo_1d->hmap_thetamup_extra_syst_multisim_bs, "nue", fname_extra_syst, wgts_extra_syst);
      if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(thetamup, event_weight, _event_histo_1d->hmap_thetamup_mc_stat_multisim_bs, "nue", fname_mc_stat_multisim, wgts_mc_stat_multisim);

      if(muon_contained){
      _event_histo_1d->hmap_trkmom["nue"]->Fill(t->pfp_reco_Mom_muon[muind], event_weight);
      } else {
      _event_histo_1d->hmap_trkmom["nue"]->Fill(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight);
      }
      _event_histo_1d->hmap_trkpmom["nue"]->Fill(t->pfp_reco_Mom_proton[pind], event_weight);
      _event_histo_1d->hmap_trktheta["nue"]->Fill(t->slc_longesttrack_theta.at(scl_ll_max), event_weight);
      _event_histo_1d->hmap_trkptheta["nue"]->Fill(t->pfp_reco_costheta[pind], event_weight);
      _event_histo_1d->hmap_thetamup["nue"]->Fill(thetamup, event_weight); 

      _event_histo->hmap_trktheta_trkmom["nue"]->Fill(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight);
      _event_histo->hmap_trktheta_trkmom_poly["nue"]->Fill(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight);
      if (t->nupdg == 12)
        nue_cc_selected+=t->bnb_weight;
    }
    //
    // NC
    //
    else if(nu_origin && t->ccnc==1 && t->fv==1 && trackfromneutrino){
      bkg_nc_sel += event_weight;  std::cout<<"NC"<<std::endl;
      pEff->Fill(false, t->nu_e);
      _event_histo_1d->hmap_onebin["nc"]->Fill(0.5, event_weight);
      hmap_trklen["nc"]->Fill(t->slc_longesttrack_length.at(scl_ll_max), event_weight);
      if(muon_contained){
      hmap_trkmom_classic["nc"]->Fill(t->pfp_reco_Mom_muon[muind], event_weight);
      }else {
      hmap_trkmom_classic["nc"]->Fill(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight);
      }
      hmap_trkphi["nc"]->Fill(t->slc_longesttrack_phi.at(scl_ll_max), event_weight);
      hmap_trktheta_classic["nc"]->Fill(t->slc_longesttrack_theta.at(scl_ll_max), event_weight);
      hmap_multpfp["nc"]->Fill(t->slc_mult_pfp.at(scl_ll_max), event_weight);
      hmap_multtracktol["nc"]->Fill(t->slc_mult_track_tolerance.at(scl_ll_max), event_weight);
       //hmap_trkmom_genie_pm1_bs["nc"]->Fill(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), 1., wgts_genie_pm1);
      if(muon_contained){
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->pfp_reco_Mom_muon[muind], event_weight, _event_histo_1d->hmap_trkmom_genie_pm1_bs, "nc", fname_genie_pm1, wgts_genie_pm1);
      } else{
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkmom_genie_pm1_bs, "nc", fname_genie_pm1, wgts_genie_pm1);
      }
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(temp_pmom, event_weight, _event_histo_1d->hmap_trkpmom_genie_pm1_bs, "nc", fname_genie_pm1, wgts_genie_pm1);
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->lep_costheta, event_weight, _event_histo_1d->hmap_trktheta_genie_pm1_bs, "nc", fname_genie_pm1, wgts_genie_pm1);
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(temp_pangle, event_weight, _event_histo_1d->hmap_trkptheta_genie_pm1_bs, "nc", fname_genie_pm1, wgts_genie_pm1);

      if (!isdata && _fill_bootstrap_genie) FillBootstrap(thetamup, event_weight, _event_histo_1d->hmap_thetamup_genie_pm1_bs, "nc", fname_genie_pm1, wgts_genie_pm1);




      if (!isdata && _fill_bootstrap_genie) FillBootstrap(0.5, event_weight, _event_histo_1d->hmap_onebin_genie_multisim_bs, "nc", fname_genie_multisim, wgts_genie_multisim);
      if(muon_contained){
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->pfp_reco_Mom_muon[muind], event_weight, _event_histo_1d->hmap_trkmom_genie_multisim_bs, "nc", fname_genie_multisim, wgts_genie_multisim);
      } else {
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkmom_genie_multisim_bs, "nc", fname_genie_multisim, wgts_genie_multisim);
      }
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkangle_genie_multisim_bs, "nc", fname_genie_multisim, wgts_genie_multisim);
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo->hmap_trktheta_trkmom_genie_multisim_bs, "nc", fname_genie_multisim, wgts_genie_multisim);
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo->hmap_trktheta_trkmom_poly_genie_multisim_bs, "nc", fname_genie_multisim, wgts_genie_multisim);

      if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(0.5, event_weight, _event_histo_1d->hmap_onebin_extra_syst_multisim_bs, "nc", fname_extra_syst, wgts_extra_syst);
      if(muon_contained){
      if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(t->pfp_reco_Mom_muon[muind], event_weight, _event_histo_1d->hmap_trkmom_extra_syst_multisim_bs, "nc", fname_extra_syst, wgts_extra_syst);
      }else{
      if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkmom_extra_syst_multisim_bs, "nc", fname_extra_syst, wgts_extra_syst);
      }
      if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkangle_extra_syst_multisim_bs, "nc", fname_extra_syst, wgts_extra_syst);
      if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo->hmap_trktheta_trkmom_extra_syst_multisim_bs, "nc", fname_extra_syst, wgts_extra_syst);
      if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo->hmap_trktheta_trkmom_poly_extra_syst_multisim_bs, "nc", fname_extra_syst, wgts_extra_syst);
      
      if (!isdata && _fill_bootstrap_flux) FillBootstrap(0.5, event_weight, _event_histo_1d->hmap_onebin_flux_multisim_bs, "nc", fname_flux_multisim, wgts_flux_multisim);
      if(muon_contained){
      if (!isdata && _fill_bootstrap_flux) FillBootstrap(t->pfp_reco_Mom_muon[muind], event_weight, _event_histo_1d->hmap_trkmom_flux_multisim_bs, "nc", fname_flux_multisim, wgts_flux_multisim);
      }else{
      if (!isdata && _fill_bootstrap_flux) FillBootstrap(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkmom_flux_multisim_bs, "nc", fname_flux_multisim, wgts_flux_multisim);
      }
      if (!isdata && _fill_bootstrap_flux) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkangle_flux_multisim_bs, "nc", fname_flux_multisim, wgts_flux_multisim);
      if (!isdata && _fill_bootstrap_flux) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo->hmap_trktheta_trkmom_flux_multisim_bs, "nc", fname_flux_multisim, wgts_flux_multisim);
      if (!isdata && _fill_bootstrap_flux) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo->hmap_trktheta_trkmom_poly_flux_multisim_bs, "nc", fname_flux_multisim, wgts_flux_multisim);

      if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(0.5, event_weight, _event_histo_1d->hmap_onebin_mc_stat_multisim_bs, "nc", fname_mc_stat_multisim, wgts_mc_stat_multisim);
      if(muon_contained){
      if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(t->pfp_reco_Mom_muon[muind], event_weight, _event_histo_1d->hmap_trkmom_mc_stat_multisim_bs, "nc", fname_mc_stat_multisim, wgts_mc_stat_multisim);
      }else{
      if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkmom_mc_stat_multisim_bs, "nc", fname_mc_stat_multisim, wgts_mc_stat_multisim);
      }
      if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkangle_mc_stat_multisim_bs, "nc", fname_mc_stat_multisim, wgts_mc_stat_multisim);
      if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo->hmap_trktheta_trkmom_mc_stat_multisim_bs, "nc", fname_mc_stat_multisim, wgts_mc_stat_multisim);
      if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo->hmap_trktheta_trkmom_poly_mc_stat_multisim_bs, "nc", fname_mc_stat_multisim, wgts_mc_stat_multisim);
      hmap_trkplen["nc"]->Fill(t->pfp_reco_length[pind], event_weight);
      hmap_trkpmom_classic["nc"]->Fill(t->pfp_reco_Mom_proton[pind], event_weight);
      hmap_trkptheta_classic["nc"]->Fill(t->pfp_reco_costheta[pind], event_weight);
      hmap_trkpphi["nc"]->Fill(t->pfp_reco_phi[pind], event_weight);
      hmap_thetamup["nc"]->Fill(thetamup, event_weight);
      hmap_ptmis["nc"]->Fill(ptmis, event_weight);
      if(!isdata){hmap_etatest["nc"]->Fill(etatest, event_weight);}
      hmap_alphat["nc"]->Fill(alphat, event_weight);
      hmap_phit["nc"]->Fill(phit, event_weight);
      hmap_enucal["nc"]->Fill(enucal, event_weight);
      hmap_nhits_leadingp["nc"]->Fill(t->pfp_reco_dEdx[pind].size(), event_weight);
      if(_showerastrack){
        hmap_pmult["nc"]->Fill(t->num_pfp-1, event_weight);
      } else {
        hmap_pmult["nc"]->Fill(t->num_pfp_tracks-1, event_weight);
      }

      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->pfp_reco_Mom_proton[pind], event_weight, _event_histo_1d->hmap_trkpmom_genie_multisim_bs, "nc", fname_genie_multisim, wgts_genie_multisim);      
      if (!isdata && _fill_bootstrap_flux) FillBootstrap(t->pfp_reco_Mom_proton[pind], event_weight, _event_histo_1d->hmap_trkpmom_flux_multisim_bs, "nc", fname_flux_multisim, wgts_flux_multisim);
      if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(t->pfp_reco_Mom_proton[pind], event_weight, _event_histo_1d->hmap_trkpmom_extra_syst_multisim_bs, "nc", fname_extra_syst, wgts_extra_syst);
      if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(t->pfp_reco_Mom_proton[pind], event_weight, _event_histo_1d->hmap_trkpmom_mc_stat_multisim_bs, "nc", fname_mc_stat_multisim, wgts_mc_stat_multisim);

      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->pfp_reco_costheta[pind], event_weight, _event_histo_1d->hmap_trkpangle_genie_multisim_bs, "nc", fname_genie_multisim, wgts_genie_multisim);      
      if (!isdata && _fill_bootstrap_flux) FillBootstrap(t->pfp_reco_costheta[pind], event_weight, _event_histo_1d->hmap_trkpangle_flux_multisim_bs, "nc", fname_flux_multisim, wgts_flux_multisim);
      if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(t->pfp_reco_costheta[pind], event_weight, _event_histo_1d->hmap_trkpangle_extra_syst_multisim_bs, "nc", fname_extra_syst, wgts_extra_syst);
      if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(t->pfp_reco_costheta[pind], event_weight, _event_histo_1d->hmap_trkpangle_mc_stat_multisim_bs, "nc", fname_mc_stat_multisim, wgts_mc_stat_multisim);

      if (!isdata && _fill_bootstrap_genie) FillBootstrap(thetamup, event_weight, _event_histo_1d->hmap_thetamup_genie_multisim_bs, "nc", fname_genie_multisim, wgts_genie_multisim);      
      if (!isdata && _fill_bootstrap_flux) FillBootstrap(thetamup, event_weight, _event_histo_1d->hmap_thetamup_flux_multisim_bs, "nc", fname_flux_multisim, wgts_flux_multisim);
      if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(thetamup, event_weight, _event_histo_1d->hmap_thetamup_extra_syst_multisim_bs, "nc", fname_extra_syst, wgts_extra_syst);
      if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(thetamup, event_weight, _event_histo_1d->hmap_thetamup_mc_stat_multisim_bs, "nc", fname_mc_stat_multisim, wgts_mc_stat_multisim);
      
      if(muon_contained){
      _event_histo_1d->hmap_trkmom["nc"]->Fill(t->pfp_reco_Mom_muon[muind], event_weight);
      } else {
      _event_histo_1d->hmap_trkmom["nc"]->Fill(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight);
      }
      _event_histo_1d->hmap_trkpmom["nc"]->Fill(t->pfp_reco_Mom_proton[pind], event_weight);
      _event_histo_1d->hmap_trktheta["nc"]->Fill(t->slc_longesttrack_theta.at(scl_ll_max), event_weight);
      _event_histo_1d->hmap_trkptheta["nc"]->Fill(t->pfp_reco_costheta[pind], event_weight);
      _event_histo_1d->hmap_thetamup["nc"]->Fill(thetamup, event_weight);

      _event_histo->hmap_trktheta_trkmom["nc"]->Fill(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight);
      _event_histo->hmap_trktheta_trkmom_poly["nc"]->Fill(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight);
      // proton
      if (t->slc_origin_extra.at(scl_ll_max) == 3) {
        _event_histo_1d->hmap_onebin["nc_proton"]->Fill(0.5, event_weight);
        hmap_trklen["nc_proton"]->Fill(t->slc_longesttrack_length.at(scl_ll_max), event_weight);
        _event_histo_1d->hmap_trkmom["nc_proton"]->Fill(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight);
        hmap_trkmom_classic["nc_proton"]->Fill(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight);
        //hmap_trkmom_genie_pm1_bs["nc_proton"]->Fill(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), 1., wgts_genie_pm1);

        hmap_trkphi["nc_proton"]->Fill(t->slc_longesttrack_phi.at(scl_ll_max), event_weight);
        _event_histo_1d->hmap_trktheta["nc_proton"]->Fill(t->slc_longesttrack_theta.at(scl_ll_max), event_weight);
        hmap_trktheta_classic["nc_proton"]->Fill(t->slc_longesttrack_theta.at(scl_ll_max), event_weight);
        hmap_multpfp["nc_proton"]->Fill(t->slc_mult_pfp.at(scl_ll_max), event_weight);
        hmap_multtracktol["nc_proton"]->Fill(t->slc_mult_track_tolerance.at(scl_ll_max), event_weight);
        _event_histo->hmap_trktheta_trkmom["nc_proton"]->Fill(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight);
        _event_histo->hmap_trktheta_trkmom_poly["nc_proton"]->Fill(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight);
      }
      // pion
      else if (t->slc_origin_extra.at(scl_ll_max) == 2) {
        _event_histo_1d->hmap_onebin["nc_pion"]->Fill(0.5, event_weight);
        hmap_trklen["nc_pion"]->Fill(t->slc_longesttrack_length.at(scl_ll_max), event_weight);
        _event_histo_1d->hmap_trkmom["nc_pion"]->Fill(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight);
        hmap_trkmom_classic["nc_pion"]->Fill(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight);
        //hmap_trkmom_genie_pm1_bs["nc_pion"]->Fill(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), 1., wgts_genie_pm1);

        hmap_trkphi["nc_pion"]->Fill(t->slc_longesttrack_phi.at(scl_ll_max), event_weight);
        _event_histo_1d->hmap_trktheta["nc_pion"]->Fill(t->slc_longesttrack_theta.at(scl_ll_max), event_weight);
        hmap_trktheta_classic["nc_pion"]->Fill(t->slc_longesttrack_theta.at(scl_ll_max), event_weight);
        hmap_multpfp["nc_pion"]->Fill(t->slc_mult_pfp.at(scl_ll_max), event_weight);
        hmap_multtracktol["nc_pion"]->Fill(t->slc_mult_track_tolerance.at(scl_ll_max), event_weight);
        _event_histo->hmap_trktheta_trkmom["nc_pion"]->Fill(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight);
        _event_histo->hmap_trktheta_trkmom_poly["nc_pion"]->Fill(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight);
      }
      // other
      else {
        _event_histo_1d->hmap_onebin["nc_other"]->Fill(0.5, event_weight);
        hmap_trklen["nc_other"]->Fill(t->slc_longesttrack_length.at(scl_ll_max), event_weight);
        _event_histo_1d->hmap_trkmom["nc_other"]->Fill(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight);
        hmap_trkmom_classic["nc_other"]->Fill(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight);
        //hmap_trkmom_genie_pm1_bs["nc_other"]->Fill(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), 1., wgts_genie_pm1);

        hmap_trkphi["nc_other"]->Fill(t->slc_longesttrack_phi.at(scl_ll_max), event_weight);
        _event_histo_1d->hmap_trktheta["nc_other"]->Fill(t->slc_longesttrack_theta.at(scl_ll_max), event_weight);
        hmap_trktheta_classic["nc_other"]->Fill(t->slc_longesttrack_theta.at(scl_ll_max), event_weight);
        hmap_multpfp["nc_other"]->Fill(t->slc_mult_pfp.at(scl_ll_max), event_weight);
        hmap_multtracktol["nc_other"]->Fill(t->slc_mult_track_tolerance.at(scl_ll_max), event_weight);
        _event_histo->hmap_trktheta_trkmom["nc_other"]->Fill(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight);
        _event_histo->hmap_trktheta_trkmom_poly["nc_other"]->Fill(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight);
      }
    }
    //
    // OUTFV
    //
    else if(nu_origin && t->fv==0 && trackfromneutrino){
      bkg_outfv_sel += event_weight;  std::cout<<"OUTFV"<<std::endl;
      pEff->Fill(false, t->nu_e);
      _event_histo_1d->hmap_onebin["outfv"]->Fill(0.5, event_weight);
      hmap_trklen["outfv"]->Fill(t->slc_longesttrack_length.at(scl_ll_max), event_weight);
      if(muon_contained){
      hmap_trkmom_classic["outfv"]->Fill(t->pfp_reco_Mom_muon[muind], event_weight);
      } else {
      hmap_trkmom_classic["outfv"]->Fill(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight);
      }
      hmap_trkphi["outfv"]->Fill(t->slc_longesttrack_phi.at(scl_ll_max), event_weight);
      hmap_trktheta_classic["outfv"]->Fill(t->slc_longesttrack_theta.at(scl_ll_max), event_weight);
      hmap_multpfp["outfv"]->Fill(t->slc_mult_pfp.at(scl_ll_max), event_weight);
      hmap_multtracktol["outfv"]->Fill(t->slc_mult_track_tolerance.at(scl_ll_max), event_weight);
       //hmap_trkmom_genie_pm1_bs["outfv"]->Fill(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), 1., wgts_genie_pm1);
      if(muon_contained){
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->pfp_reco_Mom_muon[muind], event_weight, _event_histo_1d->hmap_trkmom_genie_pm1_bs, "outfv", fname_genie_pm1, wgts_genie_pm1);
      } else{
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkmom_genie_pm1_bs, "outfv", fname_genie_pm1, wgts_genie_pm1);
      }
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(temp_pmom, event_weight, _event_histo_1d->hmap_trkpmom_genie_pm1_bs, "outfv", fname_genie_pm1, wgts_genie_pm1);
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->lep_costheta, event_weight, _event_histo_1d->hmap_trktheta_genie_pm1_bs, "outfv", fname_genie_pm1, wgts_genie_pm1);
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(temp_pangle, event_weight, _event_histo_1d->hmap_trkptheta_genie_pm1_bs, "outfv", fname_genie_pm1, wgts_genie_pm1);

      if (!isdata && _fill_bootstrap_genie) FillBootstrap(thetamup, event_weight, _event_histo_1d->hmap_thetamup_genie_pm1_bs, "outfv", fname_genie_pm1, wgts_genie_pm1);




      if (!isdata && _fill_bootstrap_genie) FillBootstrap(0.5, event_weight, _event_histo_1d->hmap_onebin_genie_multisim_bs, "outfv", fname_genie_multisim, wgts_genie_multisim);
      if(muon_contained){
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->pfp_reco_Mom_muon[muind], event_weight, _event_histo_1d->hmap_trkmom_genie_multisim_bs, "outfv", fname_genie_multisim, wgts_genie_multisim);
      } else {
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkmom_genie_multisim_bs, "outfv", fname_genie_multisim, wgts_genie_multisim);
      }
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkangle_genie_multisim_bs, "outfv", fname_genie_multisim, wgts_genie_multisim);
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo->hmap_trktheta_trkmom_genie_multisim_bs, "outfv", fname_genie_multisim, wgts_genie_multisim);
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo->hmap_trktheta_trkmom_poly_genie_multisim_bs, "outfv", fname_genie_multisim, wgts_genie_multisim);

      if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(0.5, event_weight, _event_histo_1d->hmap_onebin_extra_syst_multisim_bs, "outfv", fname_extra_syst, wgts_extra_syst);
      if(muon_contained){
      if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(t->pfp_reco_Mom_muon[muind], event_weight, _event_histo_1d->hmap_trkmom_extra_syst_multisim_bs, "outfv", fname_extra_syst, wgts_extra_syst);
      } else{
      if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkmom_extra_syst_multisim_bs, "outfv", fname_extra_syst, wgts_extra_syst);
      }
      if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkangle_extra_syst_multisim_bs, "outfv", fname_extra_syst, wgts_extra_syst);
      if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo->hmap_trktheta_trkmom_extra_syst_multisim_bs, "outfv", fname_extra_syst, wgts_extra_syst);
      if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo->hmap_trktheta_trkmom_poly_extra_syst_multisim_bs, "outfv", fname_extra_syst, wgts_extra_syst);
      
      if (!isdata && _fill_bootstrap_flux) FillBootstrap(0.5, event_weight, _event_histo_1d->hmap_onebin_flux_multisim_bs, "outfv", fname_flux_multisim, wgts_flux_multisim);
      if(muon_contained){
      if (!isdata && _fill_bootstrap_flux) FillBootstrap(t->pfp_reco_Mom_muon[muind], event_weight, _event_histo_1d->hmap_trkmom_flux_multisim_bs, "outfv", fname_flux_multisim, wgts_flux_multisim);
      }else{
      if (!isdata && _fill_bootstrap_flux) FillBootstrap(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkmom_flux_multisim_bs, "outfv", fname_flux_multisim, wgts_flux_multisim);
      }
      if (!isdata && _fill_bootstrap_flux) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkangle_flux_multisim_bs, "outfv", fname_flux_multisim, wgts_flux_multisim);
      if (!isdata && _fill_bootstrap_flux) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo->hmap_trktheta_trkmom_flux_multisim_bs, "outfv", fname_flux_multisim, wgts_flux_multisim);
      if (!isdata && _fill_bootstrap_flux) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo->hmap_trktheta_trkmom_poly_flux_multisim_bs, "outfv", fname_flux_multisim, wgts_flux_multisim);

      if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(0.5, event_weight, _event_histo_1d->hmap_onebin_mc_stat_multisim_bs, "outfv", fname_mc_stat_multisim, wgts_mc_stat_multisim);
      if(muon_contained){
      if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(t->pfp_reco_Mom_muon[muind], event_weight, _event_histo_1d->hmap_trkmom_mc_stat_multisim_bs, "outfv", fname_mc_stat_multisim, wgts_mc_stat_multisim);
      } else{
      if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkmom_mc_stat_multisim_bs, "outfv", fname_mc_stat_multisim, wgts_mc_stat_multisim);
      }
      if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkangle_mc_stat_multisim_bs, "outfv", fname_mc_stat_multisim, wgts_mc_stat_multisim);
      if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo->hmap_trktheta_trkmom_mc_stat_multisim_bs, "outfv", fname_mc_stat_multisim, wgts_mc_stat_multisim);
      if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo->hmap_trktheta_trkmom_poly_mc_stat_multisim_bs, "outfv", fname_mc_stat_multisim, wgts_mc_stat_multisim);
      hmap_trkplen["outfv"]->Fill(t->pfp_reco_length[pind], event_weight);
      hmap_trkpmom_classic["outfv"]->Fill(t->pfp_reco_Mom_proton[pind], event_weight);
      hmap_trkptheta_classic["outfv"]->Fill(t->pfp_reco_costheta[pind], event_weight);
      hmap_trkpphi["outfv"]->Fill(t->pfp_reco_phi[pind], event_weight);
      hmap_thetamup["outfv"]->Fill(thetamup, event_weight);
      hmap_ptmis["outfv"]->Fill(ptmis, event_weight);
      if(!isdata){hmap_etatest["outfv"]->Fill(etatest, event_weight);}
      hmap_alphat["outfv"]->Fill(alphat, event_weight);
      hmap_phit["outfv"]->Fill(phit, event_weight);
      hmap_enucal["outfv"]->Fill(enucal, event_weight);
      hmap_nhits_leadingp["outfv"]->Fill(t->pfp_reco_dEdx[pind].size(), event_weight);
      if(_showerastrack){
        hmap_pmult["outfv"]->Fill(t->num_pfp-1, event_weight);
      } else {
        hmap_pmult["outfv"]->Fill(t->num_pfp_tracks-1, event_weight);
      }

      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->pfp_reco_Mom_proton[pind], event_weight, _event_histo_1d->hmap_trkpmom_genie_multisim_bs, "outfv", fname_genie_multisim, wgts_genie_multisim);      
      if (!isdata && _fill_bootstrap_flux) FillBootstrap(t->pfp_reco_Mom_proton[pind], event_weight, _event_histo_1d->hmap_trkpmom_flux_multisim_bs, "outfv", fname_flux_multisim, wgts_flux_multisim);
      if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(t->pfp_reco_Mom_proton[pind], event_weight, _event_histo_1d->hmap_trkpmom_extra_syst_multisim_bs, "outfv", fname_extra_syst, wgts_extra_syst);
      if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(t->pfp_reco_Mom_proton[pind], event_weight, _event_histo_1d->hmap_trkpmom_mc_stat_multisim_bs, "outfv", fname_mc_stat_multisim, wgts_mc_stat_multisim);

      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->pfp_reco_costheta[pind], event_weight, _event_histo_1d->hmap_trkpangle_genie_multisim_bs, "outfv", fname_genie_multisim, wgts_genie_multisim);      
      if (!isdata && _fill_bootstrap_flux) FillBootstrap(t->pfp_reco_costheta[pind], event_weight, _event_histo_1d->hmap_trkpangle_flux_multisim_bs, "outfv", fname_flux_multisim, wgts_flux_multisim);
      if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(t->pfp_reco_costheta[pind], event_weight, _event_histo_1d->hmap_trkpangle_extra_syst_multisim_bs, "outfv", fname_extra_syst, wgts_extra_syst);
      if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(t->pfp_reco_costheta[pind], event_weight, _event_histo_1d->hmap_trkpangle_mc_stat_multisim_bs, "outfv", fname_mc_stat_multisim, wgts_mc_stat_multisim);

      if (!isdata && _fill_bootstrap_genie) FillBootstrap(thetamup, event_weight, _event_histo_1d->hmap_thetamup_genie_multisim_bs, "outfv", fname_genie_multisim, wgts_genie_multisim);      
      if (!isdata && _fill_bootstrap_flux) FillBootstrap(thetamup, event_weight, _event_histo_1d->hmap_thetamup_flux_multisim_bs, "outfv", fname_flux_multisim, wgts_flux_multisim);
      if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(thetamup, event_weight, _event_histo_1d->hmap_thetamup_extra_syst_multisim_bs, "outfv", fname_extra_syst, wgts_extra_syst);
      if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(thetamup, event_weight, _event_histo_1d->hmap_thetamup_mc_stat_multisim_bs, "outfv", fname_mc_stat_multisim, wgts_mc_stat_multisim);
       
      if(muon_contained){
      _event_histo_1d->hmap_trkmom["outfv"]->Fill(t->pfp_reco_Mom_muon[muind], event_weight);
      } else{
      _event_histo_1d->hmap_trkmom["outfv"]->Fill(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight);
      }
      _event_histo_1d->hmap_trkpmom["outfv"]->Fill(t->pfp_reco_Mom_proton[pind], event_weight);
      _event_histo_1d->hmap_trktheta["outfv"]->Fill(t->slc_longesttrack_theta.at(scl_ll_max), event_weight);
      _event_histo_1d->hmap_trkptheta["outfv"]->Fill(t->pfp_reco_costheta[pind], event_weight);
      _event_histo_1d->hmap_thetamup["outfv"]->Fill(thetamup, event_weight);

      _event_histo->hmap_trktheta_trkmom["outfv"]->Fill(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight);
      _event_histo->hmap_trktheta_trkmom_poly["outfv"]->Fill(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight);
      
      if (t->slc_origin_extra.at(scl_ll_max) == 0) {
        _event_histo_1d->hmap_onebin["outfv_stopmu"]->Fill(0.5, event_weight);
        hmap_trklen["outfv_stopmu"]->Fill(t->slc_longesttrack_length.at(scl_ll_max), event_weight);
        _event_histo_1d->hmap_trkmom["outfv_stopmu"]->Fill(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight);
        hmap_trkmom_classic["outfv_stopmu"]->Fill(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight);
        //hmap_trkmom_genie_pm1_bs["outfv_stopmu"]->Fill(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), 1., wgts_genie_pm1);

        hmap_trkphi["outfv_stopmu"]->Fill(t->slc_longesttrack_phi.at(scl_ll_max), event_weight);
        _event_histo_1d->hmap_trktheta["outfv_stopmu"]->Fill(t->slc_longesttrack_theta.at(scl_ll_max), event_weight);
        hmap_trktheta_classic["outfv_stopmu"]->Fill(t->slc_longesttrack_theta.at(scl_ll_max), event_weight);
        hmap_multpfp["outfv_stopmu"]->Fill(t->slc_mult_pfp.at(scl_ll_max), event_weight);
        hmap_multtracktol["outfv_stopmu"]->Fill(t->slc_mult_track_tolerance.at(scl_ll_max), event_weight);
        _event_histo->hmap_trktheta_trkmom["outfv_stopmu"]->Fill(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight);
        _event_histo->hmap_trktheta_trkmom_poly["outfv_stopmu"]->Fill(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight);
      } else {
        _event_histo_1d->hmap_onebin["outfv_nostopmu"]->Fill(0.5, event_weight);
        hmap_trklen["outfv_nostopmu"]->Fill(t->slc_longesttrack_length.at(scl_ll_max), event_weight);
        _event_histo_1d->hmap_trkmom["outfv_nostopmu"]->Fill(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight);
        hmap_trkmom_classic["outfv_nostopmu"]->Fill(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight);

        hmap_trkphi["outfv_nostopmu"]->Fill(t->slc_longesttrack_phi.at(scl_ll_max), event_weight);
        _event_histo_1d->hmap_trktheta["outfv_nostopmu"]->Fill(t->slc_longesttrack_theta.at(scl_ll_max), event_weight);
        hmap_trktheta_classic["outfv_nostopmu"]->Fill(t->slc_longesttrack_theta.at(scl_ll_max), event_weight);
        hmap_multpfp["outfv_nostopmu"]->Fill(t->slc_mult_pfp.at(scl_ll_max), event_weight);
        hmap_multtracktol["outfv_nostopmu"]->Fill(t->slc_mult_track_tolerance.at(scl_ll_max), event_weight);
        _event_histo->hmap_trktheta_trkmom["outfv_nostopmu"]->Fill(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight);
        _event_histo->hmap_trktheta_trkmom_poly["outfv_nostopmu"]->Fill(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight);
      }
    }
    //
    // COSMIC
    //
    else {
      bkg_cosmic_sel += event_weight; std::cout<<"COSMIC"<<std::endl;   
      // Add extra weight to event_weight if we are scaling the cosmic background (for example from overlays)
      if (_scale_cosmics) event_weight *= _scale_factor_cosmic;
      if (t->slc_crosses_top_boundary.at(scl_ll_max) == 1 ) bkg_cosmic_top_sel++;
      pEff->Fill(false, t->nu_e);
      _event_histo_1d->hmap_onebin["cosmic"]->Fill(0.5, event_weight);
      hmap_trklen["cosmic"]->Fill(t->slc_longesttrack_length.at(scl_ll_max), event_weight);
      if(muon_contained){
      hmap_trkmom_classic["cosmic"]->Fill(t->pfp_reco_Mom_muon[muind], event_weight);
      } else {
      hmap_trkmom_classic["cosmic"]->Fill(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight);
      }
      hmap_trkphi["cosmic"]->Fill(t->slc_longesttrack_phi.at(scl_ll_max), event_weight);
      hmap_trktheta_classic["cosmic"]->Fill(t->slc_longesttrack_theta.at(scl_ll_max), event_weight);
      hmap_multpfp["cosmic"]->Fill(t->slc_mult_pfp.at(scl_ll_max), event_weight);
      hmap_multtracktol["cosmic"]->Fill(t->slc_mult_track_tolerance.at(scl_ll_max), event_weight);
       //hmap_trkmom_genie_pm1_bs["cosmic"]->Fill(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), 1., wgts_genie_pm1);
      if(muon_contained){
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->pfp_reco_Mom_muon[muind], event_weight, _event_histo_1d->hmap_trkmom_genie_pm1_bs, "cosmic", fname_genie_pm1, wgts_genie_pm1);
      }else{
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkmom_genie_pm1_bs, "cosmic", fname_genie_pm1, wgts_genie_pm1);
      }
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(temp_pmom, event_weight, _event_histo_1d->hmap_trkpmom_genie_pm1_bs, "cosmic", fname_genie_pm1, wgts_genie_pm1);
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->lep_costheta, event_weight, _event_histo_1d->hmap_trktheta_genie_pm1_bs, "cosmic", fname_genie_pm1, wgts_genie_pm1);
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(temp_pangle, event_weight, _event_histo_1d->hmap_trkptheta_genie_pm1_bs, "cosmic", fname_genie_pm1, wgts_genie_pm1);

      if (!isdata && _fill_bootstrap_genie) FillBootstrap(thetamup, event_weight, _event_histo_1d->hmap_thetamup_genie_pm1_bs, "cosmic", fname_genie_pm1, wgts_genie_pm1);




      if (!isdata && _fill_bootstrap_genie) FillBootstrap(0.5, event_weight, _event_histo_1d->hmap_onebin_genie_multisim_bs, "cosmic", fname_genie_multisim, wgts_genie_multisim);
      if(muon_contained){
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->pfp_reco_Mom_muon[muind], event_weight, _event_histo_1d->hmap_trkmom_genie_multisim_bs, "cosmic", fname_genie_multisim, wgts_genie_multisim);
      } else {
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkmom_genie_multisim_bs, "cosmic", fname_genie_multisim, wgts_genie_multisim);
      }
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkangle_genie_multisim_bs, "cosmic", fname_genie_multisim, wgts_genie_multisim);
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo->hmap_trktheta_trkmom_genie_multisim_bs, "cosmic", fname_genie_multisim, wgts_genie_multisim);
      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo->hmap_trktheta_trkmom_poly_genie_multisim_bs, "cosmic", fname_genie_multisim, wgts_genie_multisim);

      if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(0.5, event_weight, _event_histo_1d->hmap_onebin_extra_syst_multisim_bs, "cosmic", fname_extra_syst, wgts_extra_syst);
      if(muon_contained){
      if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(t->pfp_reco_Mom_muon[muind], event_weight, _event_histo_1d->hmap_trkmom_extra_syst_multisim_bs, "cosmic", fname_extra_syst, wgts_extra_syst);
      } else{
      if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkmom_extra_syst_multisim_bs, "cosmic", fname_extra_syst, wgts_extra_syst);
      }
      if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkangle_extra_syst_multisim_bs, "cosmic", fname_extra_syst, wgts_extra_syst);
      if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo->hmap_trktheta_trkmom_extra_syst_multisim_bs, "cosmic", fname_extra_syst, wgts_extra_syst);
      if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo->hmap_trktheta_trkmom_poly_extra_syst_multisim_bs, "cosmic", fname_extra_syst, wgts_extra_syst);
      
      if (!isdata && _fill_bootstrap_flux) FillBootstrap(0.5, event_weight, _event_histo_1d->hmap_onebin_flux_multisim_bs, "cosmic", fname_flux_multisim, wgts_flux_multisim);
      if(muon_contained){
      if (!isdata && _fill_bootstrap_flux) FillBootstrap(t->pfp_reco_Mom_muon[muind], event_weight, _event_histo_1d->hmap_trkmom_flux_multisim_bs, "cosmic", fname_flux_multisim, wgts_flux_multisim);
      }else{
      if (!isdata && _fill_bootstrap_flux) FillBootstrap(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkmom_flux_multisim_bs, "cosmic", fname_flux_multisim, wgts_flux_multisim);
      }
      if (!isdata && _fill_bootstrap_flux) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkangle_flux_multisim_bs, "cosmic", fname_flux_multisim, wgts_flux_multisim);
      if (!isdata && _fill_bootstrap_flux) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo->hmap_trktheta_trkmom_flux_multisim_bs, "cosmic", fname_flux_multisim, wgts_flux_multisim);
      if (!isdata && _fill_bootstrap_flux) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo->hmap_trktheta_trkmom_poly_flux_multisim_bs, "cosmic", fname_flux_multisim, wgts_flux_multisim);

      if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(0.5, event_weight, _event_histo_1d->hmap_onebin_mc_stat_multisim_bs, "cosmic", fname_mc_stat_multisim, wgts_mc_stat_multisim);
      if(muon_contained){
      if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(t->pfp_reco_Mom_muon[muind], event_weight, _event_histo_1d->hmap_trkmom_mc_stat_multisim_bs, "cosmic", fname_mc_stat_multisim, wgts_mc_stat_multisim);
      } else{
      if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkmom_mc_stat_multisim_bs, "cosmic", fname_mc_stat_multisim, wgts_mc_stat_multisim);
      }
      if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), event_weight, _event_histo_1d->hmap_trkangle_mc_stat_multisim_bs, "cosmic", fname_mc_stat_multisim, wgts_mc_stat_multisim);
      if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo->hmap_trktheta_trkmom_mc_stat_multisim_bs, "cosmic", fname_mc_stat_multisim, wgts_mc_stat_multisim);
      if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight, _event_histo->hmap_trktheta_trkmom_poly_mc_stat_multisim_bs, "cosmic", fname_mc_stat_multisim, wgts_mc_stat_multisim);
      hmap_trkplen["cosmic"]->Fill(t->pfp_reco_length[pind], event_weight);
      hmap_trkpmom_classic["cosmic"]->Fill(t->pfp_reco_Mom_proton[pind], event_weight);
      hmap_trkptheta_classic["cosmic"]->Fill(t->pfp_reco_costheta[pind], event_weight);
      hmap_trkpphi["cosmic"]->Fill(t->pfp_reco_phi[pind], event_weight);
      hmap_thetamup["cosmic"]->Fill(thetamup, event_weight);
      hmap_ptmis["cosmic"]->Fill(ptmis, event_weight);
      if(!isdata){hmap_etatest["cosmic"]->Fill(etatest, event_weight);}
      hmap_alphat["cosmic"]->Fill(alphat, event_weight);
      hmap_phit["cosmic"]->Fill(phit, event_weight);
      hmap_enucal["cosmic"]->Fill(enucal, event_weight);
      hmap_nhits_leadingp["cosmic"]->Fill(t->pfp_reco_dEdx[pind].size(), event_weight);
      if(_showerastrack){
        hmap_pmult["cosmic"]->Fill(t->num_pfp-1, event_weight);
      } else {
        hmap_pmult["cosmic"]->Fill(t->num_pfp_tracks-1, event_weight);
      }

      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->pfp_reco_Mom_proton[pind], event_weight, _event_histo_1d->hmap_trkpmom_genie_multisim_bs, "cosmic", fname_genie_multisim, wgts_genie_multisim);      
      if (!isdata && _fill_bootstrap_flux) FillBootstrap(t->pfp_reco_Mom_proton[pind], event_weight, _event_histo_1d->hmap_trkpmom_flux_multisim_bs, "cosmic", fname_flux_multisim, wgts_flux_multisim);
      if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(t->pfp_reco_Mom_proton[pind], event_weight, _event_histo_1d->hmap_trkpmom_extra_syst_multisim_bs, "cosmic", fname_extra_syst, wgts_extra_syst);
      if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(t->pfp_reco_Mom_proton[pind], event_weight, _event_histo_1d->hmap_trkpmom_mc_stat_multisim_bs, "cosmic", fname_mc_stat_multisim, wgts_mc_stat_multisim);

      if (!isdata && _fill_bootstrap_genie) FillBootstrap(t->pfp_reco_costheta[pind], event_weight, _event_histo_1d->hmap_trkpangle_genie_multisim_bs, "cosmic", fname_genie_multisim, wgts_genie_multisim);      
      if (!isdata && _fill_bootstrap_flux) FillBootstrap(t->pfp_reco_costheta[pind], event_weight, _event_histo_1d->hmap_trkpangle_flux_multisim_bs, "cosmic", fname_flux_multisim, wgts_flux_multisim);
      if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(t->pfp_reco_costheta[pind], event_weight, _event_histo_1d->hmap_trkpangle_extra_syst_multisim_bs, "cosmic", fname_extra_syst, wgts_extra_syst);
      if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(t->pfp_reco_costheta[pind], event_weight, _event_histo_1d->hmap_trkpangle_mc_stat_multisim_bs, "cosmic", fname_mc_stat_multisim, wgts_mc_stat_multisim);

      if (!isdata && _fill_bootstrap_genie) FillBootstrap(thetamup, event_weight, _event_histo_1d->hmap_thetamup_genie_multisim_bs, "cosmic", fname_genie_multisim, wgts_genie_multisim);      
      if (!isdata && _fill_bootstrap_flux) FillBootstrap(thetamup, event_weight, _event_histo_1d->hmap_thetamup_flux_multisim_bs, "cosmic", fname_flux_multisim, wgts_flux_multisim);
      if (!isdata && _fill_bootstrap_extra_syst) FillBootstrap(thetamup, event_weight, _event_histo_1d->hmap_thetamup_extra_syst_multisim_bs, "cosmic", fname_extra_syst, wgts_extra_syst);
      if (!isdata && _fill_bootstrap_mc_stat) FillBootstrap(thetamup, event_weight, _event_histo_1d->hmap_thetamup_mc_stat_multisim_bs, "cosmic", fname_mc_stat_multisim, wgts_mc_stat_multisim);
     
      if(muon_contained){
      _event_histo_1d->hmap_trkmom["cosmic"]->Fill(t->pfp_reco_Mom_muon[muind], event_weight);
      } else{
      _event_histo_1d->hmap_trkmom["cosmic"]->Fill(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight);
      }
      _event_histo_1d->hmap_trkpmom["cosmic"]->Fill(t->pfp_reco_Mom_proton[pind], event_weight);
      _event_histo_1d->hmap_trktheta["cosmic"]->Fill(t->slc_longesttrack_theta.at(scl_ll_max), event_weight);
      _event_histo_1d->hmap_trkptheta["cosmic"]->Fill(t->pfp_reco_costheta[pind], event_weight);
      _event_histo_1d->hmap_thetamup["cosmic"]->Fill(thetamup, event_weight);

      _event_histo->hmap_trktheta_trkmom["cosmic"]->Fill(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight);
      _event_histo->hmap_trktheta_trkmom_poly["cosmic"]->Fill(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight);
      //std::cout << "Is a cosmic but is selected. event: " << t->event << std::endl;
      
      if (t->slc_origin_extra.at(scl_ll_max) == 0) {
        _event_histo_1d->hmap_onebin["cosmic_stopmu"]->Fill(0.5, event_weight);
        hmap_trklen["cosmic_stopmu"]->Fill(t->slc_longesttrack_length.at(scl_ll_max), event_weight);
        _event_histo_1d->hmap_trkmom["cosmic_stopmu"]->Fill(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight);
        hmap_trkmom_classic["cosmic_stopmu"]->Fill(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight);
        //hmap_trkmom_genie_pm1_bs["cosmic_stopmu"]->Fill(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), 1., wgts_genie_pm1);

        hmap_trkphi["cosmic_stopmu"]->Fill(t->slc_longesttrack_phi.at(scl_ll_max), event_weight);
        _event_histo_1d->hmap_trktheta["cosmic_stopmu"]->Fill(t->slc_longesttrack_theta.at(scl_ll_max), event_weight);
        hmap_trktheta_classic["cosmic_stopmu"]->Fill(t->slc_longesttrack_theta.at(scl_ll_max), event_weight);
        hmap_multpfp["cosmic_stopmu"]->Fill(t->slc_mult_pfp.at(scl_ll_max), event_weight);
        hmap_multtracktol["cosmic_stopmu"]->Fill(t->slc_mult_track_tolerance.at(scl_ll_max), event_weight);
        _event_histo->hmap_trktheta_trkmom["cosmic_stopmu"]->Fill(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight);
        _event_histo->hmap_trktheta_trkmom_poly["cosmic_stopmu"]->Fill(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight);
      } else {
        _event_histo_1d->hmap_onebin["cosmic_nostopmu"]->Fill(0.5, event_weight);
        hmap_trklen["cosmic_nostopmu"]->Fill(t->slc_longesttrack_length.at(scl_ll_max), event_weight);
        _event_histo_1d->hmap_trkmom["cosmic_nostopmu"]->Fill(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight);
        hmap_trkmom_classic["cosmic_nostopmu"]->Fill(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight);
        //hmap_trkmom_genie_pm1_bs["cosmic_nostopmu"]->Fill(t->slc_muoncandidate_mom_mcs.at(scl_ll_max), 1., wgts_genie_pm1);

        hmap_trkphi["cosmic_nostopmu"]->Fill(t->slc_longesttrack_phi.at(scl_ll_max), event_weight);
        _event_histo_1d->hmap_trktheta["cosmic_nostopmu"]->Fill(t->slc_longesttrack_theta.at(scl_ll_max), event_weight);
        hmap_trktheta_classic["cosmic_nostopmu"]->Fill(t->slc_longesttrack_theta.at(scl_ll_max), event_weight);
        hmap_multpfp["cosmic_nostopmu"]->Fill(t->slc_mult_pfp.at(scl_ll_max), event_weight);
        hmap_multtracktol["cosmic_nostopmu"]->Fill(t->slc_mult_track_tolerance.at(scl_ll_max), event_weight);
        _event_histo->hmap_trktheta_trkmom["cosmic_nostopmu"]->Fill(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight);
        _event_histo->hmap_trktheta_trkmom_poly["cosmic_nostopmu"]->Fill(t->slc_longesttrack_theta.at(scl_ll_max), t->slc_muoncandidate_mom_mcs.at(scl_ll_max), event_weight);
      }
      std::cout<<"end of background check"<<std::endl;
      // Restore the event weight
      if (_scale_cosmics) event_weight /= _scale_factor_cosmic;
    } // end of else to fill the histogram with background events

    
  } // end of event loop
  
  std::cout<<"Start fill histograms of POT and print out section information"<<std::endl; 
  // Save POT and number of events 
  h_pot->SetBinContent(1, totalPOT);
  h_nevts->SetBinContent(1, total_events);


  

  
  
  
  
  // ************************
  //
  //  Printing
  //
  // ************************
  std::cout << std::endl << std::endl;
  LOG_NORMAL() << "Number of simulated signal events is " << nsignal << std::endl;
  // LOG_NORMAL() << "Number of ALL simulated signal events is " << nsignal_all << std::endl;
  int sel_tot = signal_sel + bkg_anumu_sel + bkg_nue_sel + bkg_nc_sel + bkg_outfv_sel + bkg_cosmic_sel + bkg_ccother_sel;
  LOG_NORMAL() << "Selected signal is " << signal_sel     << ", " << (double)signal_sel/(double)sel_tot * 100. << std::endl;
  LOG_NORMAL() << "Selected anumu is  " << bkg_anumu_sel  << ", " << (double)bkg_anumu_sel/(double)sel_tot * 100. << std::endl;
  LOG_NORMAL() << "Selected nue is    " << bkg_nue_sel    << ", " << (double)bkg_nue_sel/(double)sel_tot * 100. << std::endl;
  LOG_NORMAL() << "Selected nc is     " << bkg_nc_sel     << ", " << (double)bkg_nc_sel/(double)sel_tot * 100. << std::endl;
  LOG_NORMAL() << "Selected outfv is  " << bkg_outfv_sel  << ", " << (double)bkg_outfv_sel/(double)sel_tot * 100. << std::endl;
  LOG_NORMAL() << "Selected cosmic is " << bkg_cosmic_sel << ", " << (double)bkg_cosmic_sel/(double)sel_tot * 100. << std::endl << std::endl;
  LOG_NORMAL() << "Selected ccotheris " << bkg_ccother_sel<< ", " << (double)bkg_ccother_sel/(double)sel_tot * 100. << std::endl; 
  LOG_NORMAL() << "Efficiency: " << signal_sel/(double)nsignal << std::endl;
  LOG_NORMAL() << "Purity (does not include off-beam and dirt): " << signal_sel/(double)(signal_sel+bkg_anumu_sel+bkg_nue_sel+bkg_nc_sel+bkg_outfv_sel+bkg_cosmic_sel) << std::endl;
  LOG_NORMAL() << "Cosmic contamination: " << bkg_cosmic_sel/(double)(bkg_anumu_sel+bkg_nue_sel+bkg_nc_sel+bkg_outfv_sel+bkg_cosmic_sel) << std::endl;
  LOG_NORMAL() << "  of which crossing top: " << bkg_cosmic_top_sel/(double)bkg_cosmic_sel << std::endl;
  LOG_NORMAL() << "NC contamination: " << bkg_nc_sel/(double)(bkg_anumu_sel+bkg_nue_sel+bkg_nc_sel+bkg_outfv_sel+bkg_cosmic_sel) << std::endl;
  LOG_NORMAL() << "OUTFV contamination: " << bkg_outfv_sel/(double)(bkg_anumu_sel+bkg_nue_sel+bkg_nc_sel+bkg_outfv_sel+bkg_cosmic_sel) << std::endl << std::endl;
  
  LOG_NORMAL() << "Efficiency QE:  " << signal_sel_qe/(double)nsignal_qe   << " +- " << eff_uncertainty(signal_sel_qe, nsignal_qe) << std::endl;
  LOG_NORMAL() << "Efficiency RES: " << signal_sel_res/(double)nsignal_res << " +- " << eff_uncertainty(signal_sel_res, nsignal_res) << std::endl;
  LOG_NORMAL() << "Efficiency COH: " << signal_sel_coh/(double)nsignal_coh << " +- " << eff_uncertainty(signal_sel_coh, nsignal_coh) << std::endl;
  LOG_NORMAL() << "Efficiency DIS: " << signal_sel_dis/(double)nsignal_dis << " +- " << eff_uncertainty(signal_sel_dis, nsignal_dis) << std::endl;
  LOG_NORMAL() << "Efficiency MEC: " << signal_sel_mec/(double)nsignal_mec << " +- " << eff_uncertainty(signal_sel_mec, nsignal_mec) << std::endl << std::endl;

  LOG_NORMAL() << "Number of events with a flash in the beam spill: " << nEvtsWFlashInBeamSpill << std::endl;
  LOG_NORMAL() << "Number of events numu CC (all voulumes): " << nNumuCC << std::endl;
  LOG_NORMAL() << " Signal events that have a recon muon: " << nSignalWMuonReco << std::endl;
  LOG_NORMAL() << " Signal events that have a recon muon and a recon vertex 10 cm close in YZ plane: " << nSignalMuonRecoVtxOk << std::endl << std::endl;
  
  LOG_NORMAL() << "Number of signal events that were correctly flash-matched: " << nSignalFlashMatched << std::endl << std::endl;
  
  LOG_NORMAL() << "Number of neutrino origin slices in total: " << n_slc_nu_origin << std::endl;
    
  LOG_NORMAL() << "Number of simulated nue CC in FV (scaled to 6.6e20):                            " << nue_cc_fv                          * 6.6e20/totalPOT << std::endl;
  LOG_NORMAL() << "Number of selected nue CC in FV (as such) (scaled to 6.6e20):                   " << nue_cc_selected                    * 6.6e20/totalPOT << std::endl;
  LOG_NORMAL() << "Number of selected nue CC in FV (total) (scaled to 6.6e20):                     " << nue_cc_selected_total              * 6.6e20/totalPOT << std::endl;
  LOG_NORMAL() << "Number of selected nue CC in FV in [0.05, 1.5] GeV (total) (scaled to 6.6e20):  " << nue_cc_selected_total_energy_range * 6.6e20/totalPOT << std::endl;
  LOG_NORMAL() << "Number of selected nue in [0.05, 1.5] GeV (total) (scaled to 6.6e20):           " << nue_selected_total_energy_range    * 6.6e20/totalPOT << std::endl << std::endl << std::endl;
  
  LOG_NORMAL() << "Number of selected nue where an electron is selected (scaled to 6.6e20):        " << n_nue_electron                     * 6.6e20/totalPOT << std::endl;
  LOG_NORMAL() << "Number of selected nue where a proton is selected (scaled to 6.6e20):           " << n_nue_proton                       * 6.6e20/totalPOT << std::endl;
  LOG_NORMAL() << "Number of selected nue where a pion is selected (scaled to 6.6e20):             " << n_nue_pion                         * 6.6e20/totalPOT << std::endl;

  std::cout << std::endl;

  std::cout << std::endl;
  std::sort(run_numbers.begin(), run_numbers.end());
  LOG_NORMAL() << "First analysed run: " << run_numbers.at(0) << std::endl;
  LOG_NORMAL() << "Last analysed run: " << run_numbers.at(run_numbers.size()-1) << std::endl;


  // ************************
  //
  //  Plotting
  //
  // ************************
  
  TString temp2;
  
  TCanvas * canvas_efficiency = new TCanvas();
  TEfficiency* pEff2 = new TEfficiency(*h_eff_num,*h_eff_den);
  pEff2->SetTitle(";True Neutrino Energy [GeV];Efficiency");
  pEff2->SetLineColor(kGreen+3);
  pEff2->SetMarkerColor(kGreen+3);
  pEff2->SetMarkerStyle(20);
  pEff2->SetMarkerSize(0.5);
  pEff2->Draw("AP");
  gPad->Update();
  auto g = pEff2->GetPaintedGraph();
  g->SetMinimum(0);
  g->SetMaximum(1);
  gPad->Update();
  PlottingTools::DrawSimulationXSec();
  
  temp2 = "./output/efficiency";
  canvas_efficiency->SaveAs(temp2 + ".pdf");
  canvas_efficiency->SaveAs(temp2 + ".C","C");
  
  TCanvas * canvas_muon_reco_efficiency = new TCanvas();
  TEfficiency* pEff3 = new TEfficiency(*h_mueff_num,*h_mueff_den);
  pEff3->SetTitle(";True Muon Momentum [GeV];Reconstruction Efficiency");
  pEff3->Draw("AP");
  gPad->Update();
  g = pEff3->GetPaintedGraph();
  g->SetMinimum(0);
  g->SetMaximum(1);
  gPad->Update();
  PlottingTools::DrawSimulationXSec();
  
  temp2 = "./output/muon_reco_efficiency";
  canvas_muon_reco_efficiency->SaveAs(temp2 + ".pdf");
  canvas_muon_reco_efficiency->SaveAs(temp2 + ".C","C");

  TCanvas * canvas_muon_reco_efficiency_angle = new TCanvas();
  TEfficiency* pEff3_2 = new TEfficiency(*h_mueff_angle_num,*h_mueff_angle_den);
  pEff3_2->SetTitle(";True Muon cos(#theta);Reconstruction Efficiency");
  pEff3_2->Draw("AP");
  gPad->Update();
  g = pEff3_2->GetPaintedGraph();
  g->SetMinimum(0);
  g->SetMaximum(1);
  gPad->Update();
  PlottingTools::DrawSimulationXSec();
  
  temp2 = "./output/muon_reco_efficiency_angle";
  canvas_muon_reco_efficiency_angle->SaveAs(temp2 + ".pdf");
  canvas_muon_reco_efficiency_angle->SaveAs(temp2 + ".C","C");
  
  
  TCanvas * canvas_efficiency_mumom = new TCanvas();
  TEfficiency* pEff4 = new TEfficiency(*_event_histo_1d->h_eff_mumom_num,*_event_histo_1d->h_eff_mumom_den);
  pEff4->SetTitle(";True Muon Momentum [GeV];Efficiency");
  pEff4->SetLineColor(kGreen+3);
  pEff4->SetMarkerColor(kGreen+3);
  pEff4->SetMarkerStyle(20);
  pEff4->SetMarkerSize(0.5);
  pEff4->Draw("AP");
  gPad->Update();
  g = pEff4->GetPaintedGraph();
  g->SetMinimum(0);
  g->SetMaximum(1);
  gPad->Update();
  PlottingTools::DrawSimulationXSec();

  temp2 = "./output/efficiency_mumom";
  canvas_efficiency_mumom->SaveAs(temp2 + ".pdf");
  canvas_efficiency_mumom->SaveAs(temp2 + ".C","C");
  
  
  TCanvas * canvas_efficiency_muangle = new TCanvas();
  TEfficiency* pEff5 = new TEfficiency(*_event_histo_1d->h_eff_muangle_num,*_event_histo_1d->h_eff_muangle_den);
  pEff5->SetTitle(";True Muon cos(#theta);Efficiency");
  pEff5->SetLineColor(kGreen+3);
  pEff5->SetMarkerColor(kGreen+3);
  pEff5->SetMarkerStyle(20);
  pEff5->SetMarkerSize(0.5);
  pEff5->Draw("AP");
  gPad->Update();
  g = pEff5->GetPaintedGraph();
  g->SetMinimum(0);
  g->SetMaximum(1);
  gPad->Update();
  PlottingTools::DrawSimulationXSec();
  
  temp2 = "./output/efficiency_muangle";
  canvas_efficiency_muangle->SaveAs(temp2 + ".pdf");
  canvas_efficiency_muangle->SaveAs(temp2 + ".C","C");


  TCanvas * canvas_efficiency_muangle_mumom = new TCanvas();
  TEfficiency* pEff5_3 = new TEfficiency(*_event_histo->h_eff_muangle_mumom_num,*_event_histo->h_eff_muangle_mumom_den);
  pEff5_3->SetTitle("Efficiency;True Muon cos(#theta);True Muon Momentum [GeV]");
  pEff5_3->SetLineColor(kGreen+3);
  pEff5_3->SetMarkerColor(kGreen+3);
  pEff5_3->SetMarkerStyle(20);
  pEff5_3->SetMarkerSize(0.5);
  pEff5_3->Draw("colz");
  PlottingTools::DrawSimulationXSec();

  
  temp2 = "./output/efficiency_muangle_mumom";
  canvas_efficiency_muangle_mumom->SaveAs(temp2 + ".pdf");
  canvas_efficiency_muangle_mumom->SaveAs(temp2 + ".C","C");
  
  
  TCanvas * canvas_efficiency_muphi = new TCanvas();
  TEfficiency* pEff5_2 = new TEfficiency(*_event_histo_1d->h_eff_muphi_num,*_event_histo_1d->h_eff_muphi_den);
  pEff5_2->SetTitle(";True Muon #phi angle;Efficiency");
  pEff5_2->SetLineColor(kGreen+3);
  pEff5_2->SetMarkerColor(kGreen+3);
  pEff5_2->SetMarkerStyle(20);
  pEff5_2->SetMarkerSize(0.5);
  pEff5_2->Draw("AP");
  gPad->Update();
  g = pEff5_2->GetPaintedGraph();
  g->SetMinimum(0);
  g->SetMaximum(1);
  gPad->Update();
  PlottingTools::DrawSimulationXSec();
  
  temp2 = "./output/efficiency_muphi";
  canvas_efficiency_muphi->SaveAs(temp2 + ".pdf");
  canvas_efficiency_muphi->SaveAs(temp2 + ".C","C");
  
  
  TCanvas * canvas_efficiency_mult = new TCanvas();
  TEfficiency* pEff6 = new TEfficiency(*h_eff_mult_num,*h_eff_mult_den);
  pEff6->SetTitle(";True GENIE Particle Multiplicity;Efficiency");
  pEff6->SetLineColor(kGreen+3);
  pEff6->SetMarkerColor(kGreen+3);
  pEff6->SetMarkerStyle(20);
  pEff6->SetMarkerSize(0.5);
  pEff6->Draw("AP");
  gPad->Update();
  g = pEff6->GetPaintedGraph();
  g->SetMinimum(0);
  g->SetMaximum(1);
  gPad->Update();
  PlottingTools::DrawSimulationXSec();
  
  temp2 = "./output/efficiency_mult";
  canvas_efficiency_mult->SaveAs(temp2 + ".pdf");
  canvas_efficiency_mult->SaveAs(temp2 + ".C","C");
  
  
  TCanvas * canvas_efficiency_mult_ch = new TCanvas();
  TEfficiency* pEff7 = new TEfficiency(*h_eff_mult_ch_num,*h_eff_mult_ch_den);
  pEff7->SetTitle(";True GENIE Charged Particle Multiplicity;Efficiency");
  pEff7->SetLineColor(kGreen+3);
  pEff7->SetMarkerColor(kGreen+3);
  pEff7->SetMarkerStyle(20);
  pEff7->SetMarkerSize(0.5);
  pEff7->Draw("AP");
  gPad->Update();
  g = pEff7->GetPaintedGraph();
  g->SetMinimum(0);
  g->SetMaximum(1);
  gPad->Update();
  PlottingTools::DrawSimulationXSec();
  
  temp2 = "./output/efficiency_mult_ch";
  canvas_efficiency_mult_ch->SaveAs(temp2 + ".pdf");
  canvas_efficiency_mult_ch->SaveAs(temp2 + ".C","C");

  TCanvas * canvas_efficiency_pmom = new TCanvas();
  TEfficiency* pEff4_proton = new TEfficiency(*_event_histo_1d->h_eff_pmom_num,*_event_histo_1d->h_eff_pmom_den);
  pEff4_proton->SetTitle(";True Proton Momentum [GeV];Efficiency");
  pEff4_proton->SetLineColor(kGreen+3);
  pEff4_proton->SetMarkerColor(kGreen+3);
  pEff4_proton->SetMarkerStyle(20);
  pEff4_proton->SetMarkerSize(0.5);
  pEff4_proton->Draw("AP");
  gPad->Update();
  g = pEff4_proton->GetPaintedGraph();
  g->SetMinimum(0);
  g->SetMaximum(1);
  gPad->Update();
  PlottingTools::DrawSimulationXSec();

  temp2 = "./output/efficiency_pmom";
  canvas_efficiency_pmom->SaveAs(temp2 + ".pdf");
  canvas_efficiency_pmom->SaveAs(temp2 + ".C","C");
 
  TCanvas * canvas_efficiency_pangle = new TCanvas();
  TEfficiency* pEff5_proton = new TEfficiency(*_event_histo_1d->h_eff_pangle_num,*_event_histo_1d->h_eff_pangle_den);
  pEff5_proton->SetTitle(";True Proton CosTheta;Efficiency");
  pEff5_proton->SetLineColor(kGreen+3);
  pEff5_proton->SetMarkerColor(kGreen+3);
  pEff5_proton->SetMarkerStyle(20);
  pEff5_proton->SetMarkerSize(0.5);
  pEff5_proton->Draw("AP");
  gPad->Update();
  g = pEff5_proton->GetPaintedGraph();
  g->SetMinimum(0);
  g->SetMaximum(1);
  gPad->Update();
  PlottingTools::DrawSimulationXSec();

  temp2 = "./output/efficiency_pangle";
  canvas_efficiency_pangle->SaveAs(temp2 + ".pdf");
  canvas_efficiency_pangle->SaveAs(temp2 + ".C","C");

  TCanvas * canvas_efficiency_thetamup = new TCanvas();
  //TEfficiency* pEff_thetamup = new TEfficiency(*h_eff_thetamup_num,*h_eff_thetamup_den);

  TEfficiency* pEff_thetamup = new TEfficiency(*_event_histo_1d->h_eff_thetamup_num,*_event_histo_1d->h_eff_thetamup_den);
  pEff_thetamup->SetTitle(";Theta_{#mu p};Efficiency");
  pEff_thetamup->SetLineColor(kGreen+3);
  pEff_thetamup->SetMarkerColor(kGreen+3);
  pEff_thetamup->SetMarkerStyle(20);
  pEff_thetamup->SetMarkerSize(0.5);
  pEff_thetamup->Draw("AP");
  gPad->Update();
  g = pEff_thetamup->GetPaintedGraph();
  g->SetMinimum(0);
  g->SetMaximum(1);
  gPad->Update();
  PlottingTools::DrawSimulationXSec();

  temp2 = "./output/efficiency_thetamup";
  canvas_efficiency_thetamup->SaveAs(temp2 + ".pdf");
  canvas_efficiency_thetamup->SaveAs(temp2 + ".C","C");

 

  TCanvas * canvas_efficiency_mode = new TCanvas();
  TEfficiency* pEff_qe = new TEfficiency(*h_eff_qe_num,*h_eff_qe_den);
  pEff_qe->SetTitle(";True Neutrino Energy [GeV];Efficiency");
  pEff_qe->SetLineColor(kGreen+2); 
  pEff_qe->SetLineWidth(2);
  pEff_qe->SetMarkerColor(kGreen+2);
  pEff_qe->SetMarkerStyle(20);
  pEff_qe->SetMarkerSize(0.5);
  pEff_qe->Draw("ALP");
  gPad->Update();
  auto g_qe = pEff_qe->GetPaintedGraph();
  g_qe->SetMinimum(0);
  g_qe->SetMaximum(1);
  gPad->Update();

  TEfficiency* pEff_res = new TEfficiency(*h_eff_res_num,*h_eff_res_den);
  pEff_res->SetLineColor(kRed+1);
  pEff_res->SetMarkerColor(kRed+1);
  pEff_res->SetLineWidth(2);
  pEff_res->SetMarkerStyle(20);
  pEff_res->SetMarkerSize(0.5);
  pEff_res->Draw("LP same");

  TEfficiency* pEff_dis = new TEfficiency(*h_eff_dis_num,*h_eff_dis_den);
  pEff_dis->SetLineColor(kBlue+1);
  pEff_dis->SetMarkerColor(kBlue+1);
  pEff_dis->SetLineWidth(2);
  pEff_dis->SetMarkerStyle(20);
  pEff_dis->SetMarkerSize(0.5);
  pEff_dis->Draw("LP same");

  TEfficiency* pEff_coh = new TEfficiency(*h_eff_coh_num,*h_eff_coh_den);
  pEff_coh->SetLineColor(kOrange-3);
  pEff_coh->SetMarkerColor(kOrange-3);
  pEff_coh->SetLineWidth(2);
  pEff_coh->SetMarkerStyle(20);
  pEff_coh->SetMarkerSize(0.5);
  //pEff_coh->Draw("LP same");

  TEfficiency* pEff_mec = new TEfficiency(*h_eff_mec_num,*h_eff_mec_den);
  pEff_mec->SetLineColor(kMagenta+1); 
  pEff_mec->SetMarkerColor(kMagenta+1);
  pEff_mec->SetLineWidth(2);
  pEff_mec->SetMarkerStyle(20);
  pEff_mec->SetMarkerSize(0.5);
  pEff_mec->Draw("LP same");

  TLegend* leg_mode = new TLegend(0.6475645,0.1368421,0.8968481,0.3368421,NULL,"brNDC");
  leg_mode->AddEntry(pEff_qe,"GENIE QE","lep");
  leg_mode->AddEntry(pEff_res,"GENIE RES","lep");  
  leg_mode->AddEntry(pEff_dis,"GENIE DIS","lep");  
  // leg_mode->AddEntry(pEff_coh,"GENIE COH","lep");  
  leg_mode->AddEntry(pEff_mec,"GENIE MEC","lep");  
  leg_mode->Draw();

  PlottingTools::DrawSimulation();

  temp2 = "./output/efficiency_mode";
  canvas_efficiency_mode->SaveAs(temp2 + ".pdf");
  canvas_efficiency_mode->SaveAs(temp2 + ".C","C");

  
  
  TCanvas * canvas_chi2 = new TCanvas();
  h_chi2->Draw("histo");
  
  temp2 = "./output/chi2_mult";
  canvas_chi2->SaveAs(temp2 + ".pdf");
  canvas_chi2->SaveAs(temp2 + ".C","C");
  
  
  TCanvas * canvas_flsTime = new TCanvas();
  h_flsTime->Draw("histo");
  
  temp2 = "./output/flsTime";
  canvas_flsTime->SaveAs(temp2 + ".pdf");
  canvas_flsTime->SaveAs(temp2 + ".C","C");
  
  new TCanvas();
  h_nslices->Draw("histo");
  
  TCanvas * canvas_vtx_resolution = new TCanvas();
  h_vtx_resolution->Draw("histo");

  temp2 = "./output/vtx_resolution";
  canvas_vtx_resolution->SaveAs(temp2 + ".pdf");
  canvas_vtx_resolution->SaveAs(temp2 + ".C","C");
  
  new TCanvas();
  h_frac_diff->Draw("colz");
  
  new TCanvas();
  h_frac_diff_others->Draw("colz");
  
  // PE spec
  TCanvas * canvas_fm_pe_comparison = new TCanvas();
  TGraph* gr = new TGraph(32,hypo_spec_x,hypo_spec_y);
  TGraph* gr2 = new TGraph(32,meas_spec_x,meas_spec_y);
  TGraph* gr3 = new TGraph(32,numc_spec_x,numc_spec_y);
  gr->SetLineColor(kGreen+2);
  gr->SetLineWidth(2);
  gr->SetMarkerColor(kGreen+2);
  gr->SetMarkerSize(1.2);
  gr->SetMarkerStyle(20);
  gr->SetTitle("");
  gr->GetXaxis()->SetTitle("PMT ID");
  gr->GetYaxis()->SetTitle("PE Count");
  gr->Draw("ALP");
  gr2->SetLineColor(kBlue+2);
  gr2->SetLineWidth(2);
  gr2->SetMarkerColor(kBlue+2);
  gr2->SetMarkerSize(1.2);
  gr2->SetMarkerStyle(20);
  gr2->SetTitle("");
  gr2->GetXaxis()->SetTitle("PMT ID");
  gr2->GetYaxis()->SetTitle("PE Count");
  gr2->Draw("LP");
  gr3->SetLineColor(kRed+2);
  gr3->SetLineWidth(2);
  gr3->SetMarkerColor(kRed+2);
  gr3->SetMarkerSize(1.2);
  gr3->SetMarkerStyle(20);
  gr3->SetTitle("");
  gr3->GetXaxis()->SetTitle("PMT ID");
  gr3->GetYaxis()->SetTitle("PE Count");
  //gr3->Draw("LP");
  TLegend* leg = new TLegend(0.1,0.7,0.48,0.9);
  leg->AddEntry(gr,"Hypo flash","l");
  leg->AddEntry(gr2,"Reco flash","l");
  //leg->AddEntry(gr3,"Neutrino MCFlash","l");
  leg->Draw();
  temp2 = "./output/fm_pe_comparison";
  canvas_fm_pe_comparison->SaveAs(temp2 + ".pdf");
  canvas_fm_pe_comparison->SaveAs(temp2 + ".C","C");

  
  TCanvas * canvas_vtxcheck = new TCanvas();
  h_vtxcheck_angle_good->Scale(1./h_vtxcheck_angle_good->Integral());
  h_vtxcheck_angle_bad->Scale(1./h_vtxcheck_angle_bad->Integral());
  h_vtxcheck_angle_good->Draw("histo");
  h_vtxcheck_angle_bad->Draw("histo same");
  h_vtxcheck_angle_bad->SetLineColor(kRed);

  TLegend* vtx_leg = new TLegend(0.1,0.7,0.48,0.9);
  vtx_leg->AddEntry(h_vtxcheck_angle_good,"Signal","l");
  vtx_leg->AddEntry(h_vtxcheck_angle_bad,"Background","l");
  vtx_leg->Draw();
  
  temp2 = "./output/vtxcheck";
  canvas_vtxcheck->SaveAs(temp2 + ".pdf");
  canvas_vtxcheck->SaveAs(temp2 + ".C","C");
  
  
  TCanvas * canvas_mu_eff_mom = new TCanvas();
  h_mu_eff_mom->Draw("colz");
  
  temp2 = "./output/mu_eff_mom";
  canvas_mu_eff_mom->SaveAs(temp2 + ".pdf");
  canvas_mu_eff_mom->SaveAs(temp2 + ".C","C");
  
  
  TCanvas * canvas_mu_eff_mom_sel = new TCanvas();
  h_mu_eff_mom_sel->Draw("colz");
  
  temp2 = "./output/mu_eff_mom_sel";
  canvas_mu_eff_mom_sel->SaveAs(temp2 + ".pdf");
  canvas_mu_eff_mom_sel->SaveAs(temp2 + ".C","C");

  
  new TCanvas();
  //h_muon_track_pur->Draw();
  h_mu_pur_mom->Draw("colz");
  
  new TCanvas();
  h_mumom_nue->Draw("colz");
  
  new TCanvas();
  h_acpt_tagged->Draw("histo");
  
  new TCanvas();
  h_xdiff->Draw("histo");
  h_xdiff_others->Draw("histo same");
  h_xdiff_others->SetLineColor(kRed);
  //TCanvas *c16 = new TCanvas();
  //h_xdiff_others->Draw("histo");
  
  new TCanvas();
  h_zdiff->Draw("histo");
  h_zdiff_others->Draw("histo same");
  h_zdiff_others->SetLineColor(kRed);
  
  
  new TCanvas();
  h_slice_origin->Draw("histo");
  
  new TCanvas();
  h_slice_npfp->DrawNormalized("histo");
  h_slice_npfp_others->DrawNormalized("histo same");
  h_slice_npfp_others->SetLineColor(kRed);
  
  new TCanvas();
  h_slice_ntrack->DrawNormalized("histo");
  h_slice_ntrack_others->DrawNormalized("histo same");
  h_slice_ntrack_others->SetLineColor(kRed);
  
  new TCanvas();
  h_fm_score->Draw("histo");
  h_fm_score_others->Draw("histo same");
  h_fm_score_others->SetLineColor(kRed);
  
  new TCanvas();
  h_n_slc_flsmatch->Draw("histo");
  
  new TCanvas();
  h_fm_score_pe->Draw("colz");
  
  
  TCanvas * final1 = new TCanvas();
  THStack *hs_trklen = new THStack("hs_trklen",";Candidate Track Length [cm]; Selected Events");
  if (_makePlots) PlottingTools::DrawTHStack(hs_trklen, pot_scaling, _breakdownPlots, hmap_trklen);
  
  //
  // Construct legend
  // used basically for all plots
  //
  TLegend* leg2;
  if (_breakdownPlots){
    leg2 = new TLegend(0.56,0.37,0.82,0.82,NULL,"brNDC");
  } else {
    leg2 = new TLegend(0.56,0.54,0.82,0.82,NULL,"brNDC");
  }
  std::stringstream sstm;
  // numu
  if (_breakdownPlots) {
  leg2->AddEntry(hmap_trklen["signal_stopmu"],"#nu_{#mu} CC (stopping #mu)","f");
  leg2->AddEntry(hmap_trklen["signal_nostopmu"],"#nu_{#mu} CC (other)","f");
  } else {
    sstm << "#nu_{#mu} CC (signal), " << std::setprecision(2)  << hmap_trklen["signal"]->Integral() / hmap_trklen["total"]->Integral()*100. << "%";
    leg2->AddEntry(hmap_trklen["signal"],sstm.str().c_str(),"f");
    sstm.str("");
  }
  
  // nue
  sstm << "#nu_{e}, #bar{#nu}_{e} CC, " << std::setprecision(2)  << hmap_trklen["nue"]->Integral() / hmap_trklen["total"]->Integral()*100. << "%";
  leg2->AddEntry(hmap_trklen["nue"],sstm.str().c_str(),"f");
  sstm.str("");
  
  // anumu
  sstm << "#bar{#nu}_{#mu} CC, " << std::setprecision(2)  << hmap_trklen["anumu"]->Integral() / hmap_trklen["total"]->Integral()*100. << "%";
  leg2->AddEntry(hmap_trklen["anumu"],sstm.str().c_str(),"f");
  sstm.str("");
  
  // nc, outfv, cosmic
  if (_breakdownPlots) {
  leg2->AddEntry(hmap_trklen["nc_other"],"NC (other)","f");
  leg2->AddEntry(hmap_trklen["nc_pion"],"NC (pion)","f");
  leg2->AddEntry(hmap_trklen["nc_proton"],"NC (proton)","f");
  leg2->AddEntry(hmap_trklen["outfv_stopmu"],"OUTFV (stopping #mu)","f");
  leg2->AddEntry(hmap_trklen["outfv_nostopmu"],"OUTFV (other)","f");
  leg2->AddEntry(hmap_trklen["cosmic_stopmu"],"Cosmic (stopping #mu)","f");
  leg2->AddEntry(hmap_trklen["cosmic_nostopmu"],"Cosmic (other)","f");
  } else {
    sstm << "NC, " << std::setprecision(2)  << hmap_trklen["nc"]->Integral() / hmap_trklen["total"]->Integral()*100. << "%";
    leg2->AddEntry(hmap_trklen["nc"],sstm.str().c_str(),"f");
    sstm.str("");
    
    sstm << "OUTFV, " << std::setprecision(2)  << hmap_trklen["outfv"]->Integral() / hmap_trklen["total"]->Integral()*100. << "%";
    leg2->AddEntry(hmap_trklen["outfv"],sstm.str().c_str(),"f");
    sstm.str("");
    
    sstm << "Cosmic, " << std::setprecision(2)  << hmap_trklen["cosmic"]->Integral() / hmap_trklen["total"]->Integral()*100. << "%";
    leg2->AddEntry(hmap_trklen["cosmic"],sstm.str().c_str(),"f");
    sstm.str("");
  }
  leg2->AddEntry(hmap_trklen["total"],"MC Stat Unc.","f");
  leg2->Draw();
  DrawPOT2(totalPOT);
  
  temp2 = "./output/trklen";
  final1->SaveAs(temp2 + ".pdf");
  final1->SaveAs(temp2 + ".C","C");
  
  
  TCanvas * final1_1 = new TCanvas();
  THStack *hs_trkmom = new THStack("hs_trkmom",";Reconstructed Momentum [GeV]; Selected Events");
  if (_makePlots) PlottingTools::DrawTHStack(hs_trkmom, pot_scaling, _breakdownPlots, _event_histo_1d->hmap_trkmom);
  leg2->Draw();
  DrawPOT2(totalPOT);
  
  temp2 = "./output/trkmom";
  final1_1->SaveAs(temp2 + ".pdf");
  final1_1->SaveAs(temp2 + ".C","C");
  
  
  TCanvas * final2 = new TCanvas();
  THStack *hs_trkphi = new THStack("hs_trkphi",";Candidate Track #phi; Selected Events");
  if (_makePlots) PlottingTools::DrawTHStack(hs_trkphi, pot_scaling, _breakdownPlots, hmap_trkphi);
  leg2->Draw();
  DrawPOT2(totalPOT);
  
  temp2 = "./output/trkphi";
  final2->SaveAs(temp2 + ".pdf");
  final2->SaveAs(temp2 + ".C","C");
  
  
  
  TCanvas * final3 = new TCanvas();
  THStack *hs_trktheta = new THStack("hs_trktheta",";Candidate Track cos(#theta); Selected Events");
  if (_makePlots) PlottingTools::DrawTHStack(hs_trktheta, pot_scaling, _breakdownPlots, _event_histo_1d->hmap_trktheta);
  leg2->Draw();
  DrawPOT2(totalPOT);
  
  temp2 = "./output/trktheta";
  final3->SaveAs(temp2 + ".pdf");
  final3->SaveAs(temp2 + ".C","C");
  
  
  TCanvas * final4 = new TCanvas();
  THStack *hs_multpfp = new THStack("hs_multpfp",";PFP Multiplicity; Selected Events");
  if (_makePlots) PlottingTools::DrawTHStack(hs_multpfp, pot_scaling, _breakdownPlots, hmap_multpfp);
  leg2->Draw();
  DrawPOT2(totalPOT);
  
  temp2 = "./output/multpfp";
  final4->SaveAs(temp2 + ".pdf");
  final4->SaveAs(temp2 + ".C","C");
  
  
  TCanvas * final5 = new TCanvas();
  THStack *hs_multtracktol = new THStack("hs_multtracktol",";Track Multiplicity (5 cm); Selected Events");
  if (_makePlots) PlottingTools::DrawTHStack(hs_multtracktol, pot_scaling, _breakdownPlots, hmap_multtracktol);
  leg2->Draw();
  DrawPOT2(totalPOT);
  
  temp2 = "./output/multtracktol";
  final5->SaveAs(temp2 + ".pdf");
  final5->SaveAs(temp2 + ".C","C");
  
  
  TCanvas * canvas_dqdx = new TCanvas();
  THStack *hs_dqdx_trunc = new THStack("hs_dqdx_trunc",";Candidate Track <dQ/dx>_{trunc};Selected Events");
  if (_makePlots) PlottingTools::DrawTHStack3(hs_dqdx_trunc, pot_scaling, _breakdownPlots, hmap_dqdx_trunc);
  
  temp2 = "./output/dqdx_trunc";
  canvas_dqdx->SaveAs(temp2 + ".pdf");
  canvas_dqdx->SaveAs(temp2 + ".C","C");
  
  
  TCanvas * canvas_dqdx_length = new TCanvas();
  h_dqdx_trunc_length_muon->Draw("BOX");
  h_dqdx_trunc_length_muon->SetLineColor(kRed+2);
  h_dqdx_trunc_length_proton->Draw("BOX same");
  h_dqdx_trunc_length_proton->SetLineColor(kBlue+2);

  temp2 = "./output/dqdx_trunc_length";
  canvas_dqdx_length->SaveAs(temp2 + ".pdf");
  canvas_dqdx_length->SaveAs(temp2 + ".C","C");
  
  
  TCanvas * canvas_dqdx_length_muon = new TCanvas();
  h_dqdx_trunc_length_muon->Draw("colz");
  
  temp2 = "./output/dqdx_trunc_length_muon";
  canvas_dqdx_length_muon->SaveAs(temp2 + ".pdf");
  canvas_dqdx_length_muon->SaveAs(temp2 + ".C","C");

  
  TCanvas * canvas_dqdx_length_proton = new TCanvas();
  h_dqdx_trunc_length_proton->Draw("colz");
  
  temp2 = "./output/dqdx_trunc_length_proton";
  canvas_dqdx_length_proton->SaveAs(temp2 + ".pdf");
  canvas_dqdx_length_proton->SaveAs(temp2 + ".C","C");

  
  
  
  TCanvas * canvas_deltall_cosmic_stop = new TCanvas();
  h_deltall_cosmic_stop->Draw();
  
  temp2 = "./output/deltall_cosmic_stop";
  canvas_deltall_cosmic_stop->SaveAs(temp2 + ".pdf");
  canvas_deltall_cosmic_stop->SaveAs(temp2 + ".C","C");
  
  TCanvas * canvas_deltall_cosmic_nostop = new TCanvas();
  h_deltall_cosmic_nostop->Draw();
  
  temp2 = "./output/deltall_cosmic_nostop";
  canvas_deltall_cosmic_nostop->SaveAs(temp2 + ".pdf");
  canvas_deltall_cosmic_nostop->SaveAs(temp2 + ".C","C");

  TCanvas * canvas_deltall_nu = new TCanvas();
  h_deltall_nu->Draw();
  
  temp2 = "./output/deltall_nu";
  canvas_deltall_nu->SaveAs(temp2 + ".pdf");
  canvas_deltall_nu->SaveAs(temp2 + ".C","C");

  
  
  
  TCanvas * canvas_deltall_length_cosmic_stop = new TCanvas();
  h_deltall_length_cosmic_stop->Draw("colz");
  
  temp2 = "./output/deltall_cosmic_stop_length";
  canvas_deltall_length_cosmic_stop->SaveAs(temp2 + ".pdf");
  canvas_deltall_length_cosmic_stop->SaveAs(temp2 + ".C","C");
  
  TCanvas * canvas_deltall_length_cosmic_nostop = new TCanvas();
  h_deltall_length_cosmic_nostop->Draw("colz");
  
  temp2 = "./output/deltall_cosmic_nostop_length";
  canvas_deltall_length_cosmic_nostop->SaveAs(temp2 + ".pdf");
  canvas_deltall_length_cosmic_nostop->SaveAs(temp2 + ".C","C");
  
  TCanvas * canvas_deltall_length_nu = new TCanvas();
  h_deltall_length_nu->Draw("colz");
  
  temp2 = "./output/deltall_nu_length";
  canvas_deltall_length_nu->SaveAs(temp2 + ".pdf");
  canvas_deltall_length_nu->SaveAs(temp2 + ".C","C");
  
  
  
  
  
  TCanvas * canvas_nue_flash = new TCanvas();
  h_true_nu_eng_beforesel->SetLineColor(kBlue+2);
  h_true_nu_eng_beforesel->Draw();
  h_true_nu_eng_afterflash->SetLineColor(kRed+2);
  h_true_nu_eng_afterflash->Draw("same");
  h_true_nu_eng_aftersel->SetLineColor(kGreen+2);
  h_true_nu_eng_aftersel->Draw("same");
  TLegend* l01 = new TLegend(0.1,0.7,0.48,0.9);
  l01->AddEntry(h_true_nu_eng_beforesel,"Generated CC #nu_{#mu} events in FV","l");
  l01->AddEntry(h_true_nu_eng_afterflash,"CC #nu_{#mu} Passing Flash Finding","l");
  l01->AddEntry(h_true_nu_eng_aftersel,"Selected CC #nu_{#mu} events","l");
  l01->Draw();
  temp2 = "./output/nue_flash";
  canvas_nue_flash->SaveAs(temp2 + ".pdf");
  canvas_nue_flash->SaveAs(temp2 + ".C","C");


  TCanvas * canvas_nue_selected = new TCanvas();
  h_nue_selected_energy->Scale(6.6e20/totalPOT);
  h_nue_selected_energy->Draw("histo");
  DrawPOT2(totalPOT, 6.6e20);
  temp2 = "./output/nue_selected_contamination";
  canvas_nue_selected->SaveAs(temp2 + ".pdf");
  canvas_nue_selected->SaveAs(temp2 + ".C","C");

  std::cout<<"Start to print out cut table and Make Efficiency & Purity Plots"<<std::endl;

  // Efficiency for every cut

  TH1D * selected_percut = new TH1D("selected_percut", "selected_percut", 13, 0, 13);
  TH1D * selected_signal_percut = new TH1D("selected_signal_percut", "selected_percut", 13, 0, 13);
  // TH1D * generated_percut = new TH1D("generated_percut", "generated_percut", 8, 0, 7);
  TH1D * generated_signal_percut = new TH1D("generated_signal_percut", "generated_percut", 13, 0, 13);

  //double pot_scale = 35388924/72299264;
  selected_percut->SetBinContent(1, selected_events_percut["initial"]); // + 1280310);
  selected_percut->SetBinContent(2, selected_events_percut["beamflash"]); // + 821708);
  selected_percut->SetBinContent(3, selected_events_percut["flash_match"]); // + 194732);
  selected_percut->SetBinContent(4, selected_events_percut["flash_match_deltax"]); // + 154544);
  selected_percut->SetBinContent(5, selected_events_percut["flash_match_deltaz"]); // + 106802);
  selected_percut->SetBinContent(6, selected_events_percut["quality"]); // + 76023);
  selected_percut->SetBinContent(7, selected_events_percut["mcs_length_quality"]); // + 72577);
  selected_percut->SetBinContent(8, selected_events_percut["mip_consistency"]); // + 69692);
  selected_percut->SetBinContent(9, selected_events_percut["fiducial_volume"]); // + 22657);
  selected_percut->SetBinContent(10,selected_events_percut["ntrk2"]);
  selected_percut->SetBinContent(11,selected_events_percut["pinCV"]);
  selected_percut->SetBinContent(12,selected_events_percut["minCol"]);  
  selected_percut->SetBinContent(13,selected_events_percut["chi2"]);


  //std::cout<<"libo check"<<selected_events_percut["initial"]<<std::endl;


  selected_signal_percut->SetBinContent(1, selected_signal_events_percut["initial"]);
  selected_signal_percut->SetBinContent(2, selected_signal_events_percut["beamflash"]);
  selected_signal_percut->SetBinContent(3, selected_signal_events_percut["flash_match"]);
  selected_signal_percut->SetBinContent(4, selected_signal_events_percut["flash_match_deltax"]);
  selected_signal_percut->SetBinContent(5, selected_signal_events_percut["flash_match_deltaz"]);
  selected_signal_percut->SetBinContent(6, selected_signal_events_percut["quality"]);
  selected_signal_percut->SetBinContent(7, selected_signal_events_percut["mcs_length_quality"]);
  selected_signal_percut->SetBinContent(8, selected_signal_events_percut["mip_consistency"]);
  selected_signal_percut->SetBinContent(9, selected_signal_events_percut["fiducial_volume"]);
  selected_signal_percut->SetBinContent(10,selected_signal_events_percut["ntrk2"]);
  selected_signal_percut->SetBinContent(11,selected_signal_events_percut["pinCV"]);
  selected_signal_percut->SetBinContent(12,selected_signal_events_percut["minCol"]);  
  selected_signal_percut->SetBinContent(13,selected_signal_events_percut["chi2"]);


  std::vector<std::string> cut_names = {"initial", "beamflash", "flashmatch", "flashmatchdeltax", "flashmatchdeltaz", "quality", "mcslengthquality", "mipconsistency", "fiducialvolume", "ntrk2", "pinCV", "minCol", "chi2"};

  for (int i = 0; i < 13; i++) {
    std::cout << cut_names.at(i) << " & " << selected_signal_percut->GetBinContent(i+1) 
       << " => " << selected_percut->GetBinContent(i+1) 
              << " & " << selected_signal_percut->GetBinContent(i+1)/selected_signal_percut->GetBinContent(1) * 100 
              << " & " << selected_signal_percut->GetBinContent(i+1)/selected_signal_percut->GetBinContent(i) * 100  << "\\\\" << std::endl;
    generated_signal_percut->SetBinContent(i+1, (double)nsignal);
    //generated_percut->SetBinContent(i+1, (double)sel_tot);
  }

  TCanvas * canvas_eff_pur_graph_percut = new TCanvas();

  canvas_eff_pur_graph_percut->SetLeftMargin(0.05157593);
  canvas_eff_pur_graph_percut->SetRightMargin(0.1475645);
  canvas_eff_pur_graph_percut->SetTopMargin(0.04210526);
  canvas_eff_pur_graph_percut->SetBottomMargin(0.1578947);

  TH1F *h = new TH1F("h","",13, 0, 13);
  h->SetMaximum(1);
  h->GetXaxis()->SetBinLabel(1,"Initial");
  h->GetXaxis()->SetBinLabel(2,"Beam Flash");
  h->GetXaxis()->SetBinLabel(3,"Flash Match");
  h->GetXaxis()->SetBinLabel(4,"Flash Match #Deltax");
  h->GetXaxis()->SetBinLabel(5,"Flash Match #Deltaz");
  h->GetXaxis()->SetBinLabel(6,"Track Quality");
  h->GetXaxis()->SetBinLabel(7,"MCS-Length Quality");
  h->GetXaxis()->SetBinLabel(8,"MIP Consistency");
  h->GetXaxis()->SetBinLabel(9,"Fiducial Volume");
  h->GetXaxis()->SetBinLabel(10,"ntrk");
  h->GetXaxis()->SetBinLabel(11,"p-contained");
  h->GetXaxis()->SetBinLabel(12,"minCol");
  h->GetXaxis()->SetBinLabel(13,"Chi2");

  h->GetXaxis()->SetLabelOffset(0.009);
  h->GetXaxis()->SetLabelSize(0.06);

  h->Draw();
  std::cout<<"Calculating efficiency percut"<<std::endl;
  TEfficiency* pEff_percut = new TEfficiency(*selected_signal_percut,*generated_signal_percut);
  pEff_percut->SetTitle("EfficiencyPerCut;Cut index;Efficiency");
  pEff_percut->SetLineColor(kGreen+3);
  pEff_percut->SetMarkerColor(kGreen+3);
  pEff_percut->SetMarkerStyle(20);
  pEff_percut->SetMarkerSize(0.6);
  TGraphAsymmErrors * pEff_percut_graph = pEff_percut->CreateGraph();
  for (int i = 0; i < 13; i++) {
    pEff_percut_graph->SetPointEXhigh(i, 0.);
    pEff_percut_graph->SetPointEXlow(i, 0.);
  }
  auto axis = pEff_percut_graph->GetYaxis();
  axis->SetLimits(0.,1.); 

  pEff_percut_graph->GetHistogram()->GetXaxis()->SetBinLabel(1,"Initial");
  pEff_percut_graph->GetHistogram()->GetXaxis()->SetBinLabel(2,"BeamFlash");
  pEff_percut_graph->GetHistogram()->GetXaxis()->SetBinLabel(3,"FlashMatch");
  pEff_percut_graph->GetHistogram()->GetXaxis()->SetBinLabel(4,"FlashMatch#Deltax");
  pEff_percut_graph->GetHistogram()->GetXaxis()->SetBinLabel(5,"FlashMatch#Deltaz");
  pEff_percut_graph->GetHistogram()->GetXaxis()->SetBinLabel(6,"TrackQuality");
  pEff_percut_graph->GetHistogram()->GetXaxis()->SetBinLabel(7,"MCS-LengthQuality");
  pEff_percut_graph->GetHistogram()->GetXaxis()->SetBinLabel(8,"MIPConsistency");
  pEff_percut_graph->GetHistogram()->GetXaxis()->SetBinLabel(9,"FiducialVolume");
  pEff_percut_graph->GetHistogram()->GetXaxis()->SetBinLabel(10,"ntrk");
  pEff_percut_graph->GetHistogram()->GetXaxis()->SetBinLabel(11,"p-contained");
  pEff_percut_graph->GetHistogram()->GetXaxis()->SetBinLabel(12,"minCol");
  pEff_percut_graph->GetHistogram()->GetXaxis()->SetBinLabel(13,"Chi2");



  // pEff_percut_graph->GetXaxis()->SetLabelSize(0.07);
  pEff_percut_graph->Draw("PL");

  std::cout<<"Calculating purity percut"<<std::endl;  

  TEfficiency* pPur_percut = new TEfficiency(*selected_signal_percut,*selected_percut);
  pPur_percut->SetTitle("PurityPerCut;Cut index;Purity");
  pPur_percut->SetLineColor(kRed+3);
  pPur_percut->SetMarkerColor(kRed+3);
  pPur_percut->SetMarkerStyle(20);
  pPur_percut->SetMarkerSize(0.6);
  TGraphAsymmErrors * pPur_percut_graph = pPur_percut->CreateGraph();
  for (int i = 0; i < 13; i++) {
    pPur_percut_graph->SetPointEXhigh(i, 0.);
    pPur_percut_graph->SetPointEXlow(i, 0.);
  }

  pPur_percut_graph->Draw("PL");

  TLegend* l = new TLegend(0.4842407,0.8168421,0.777937,0.9221053,NULL,"brNDC");
  l->AddEntry(pEff_percut_graph,"Efficiency");
  l->AddEntry(pPur_percut_graph,"Purity");
  //leg->AddEntry(gr3,"Neutrino MCFlash","l");
  
  l->Draw();

  // PlottingTools::DrawSimulationXSec();

  TLatex* prelim = new TLatex(0.8524355,0.9810526, "MicroBooNE Simulation, Preliminary");
  prelim->SetTextFont(62);
  prelim->SetTextColor(kGray+2);
  prelim->SetNDC();
  prelim->SetTextSize(1/30.);
  prelim->SetTextAlign(32);
  // prelim->SetTextSize(0.04631579);
  prelim->Draw();

  temp2 = "./output/eff_pur_graph_percut";
  canvas_eff_pur_graph_percut->SaveAs(temp2 + ".pdf");
  canvas_eff_pur_graph_percut->SaveAs(temp2 + ".C","C");

  
  
  //
  // Save on file
  //

  LOG_NORMAL() << "Saving to file." << std::endl;
  
  file_out->cd();
  for (auto iter : hmap_trklen) {
    iter.second->Write();
  }
  h_pot->Write();
  h_nevts->Write();

  pEff4->Write();
  pEff5->Write();
  pEff4_proton->Write();
  pEff5_proton->Write();
  pEff5_3->Write();
  h_truth_xsec_mumom->Write();
  h_truth_xsec_muangle->Write();
 
  h_truth_xsec_pmom->Write(); 
  h_truth_xsec_pangle->Write();
 
  pEff_percut->Write();
  pPur_percut->Write();

  file_out->WriteObject(&hmap_trklen, "hmap_trklen");
  file_out->WriteObject(&hmap_trkplen, "hmap_trkplen");
  file_out->WriteObject(&hmap_trkmom_classic, "hmap_trkmom_classic");
  file_out->WriteObject(&hmap_trkpmom_classic, "hmap_trkpmom_classic");
  file_out->WriteObject(&hmap_trktheta_classic, "hmap_trktheta_classic");
  file_out->WriteObject(&hmap_trkptheta_classic, "hmap_trkptheta_classic");
  file_out->WriteObject(&hmap_trkphi, "hmap_trkphi");
  file_out->WriteObject(&hmap_trkpphi, "hmap_trkpphi");

  file_out->WriteObject(&chi2_proton_hypothesis, "chi2_proton_hypothesis");
  file_out->WriteObject(&chi2_muon_hypothesis, "chi2_muon_hypothesis");
  file_out->WriteObject(&chi2_pion_hypothesis, "chi2_pion_hypothesis");
  file_out->WriteObject(&chi2_kaon_hypothesis, "chi2_kaon_hypothesis");


  file_out->WriteObject(&hmap_thetamup, "hmap_thetamup");
  file_out->WriteObject(&hmap_ptmis, "hmap_ptmis");
  file_out->WriteObject(&hmap_etatest, "hmap_etatest");
  file_out->WriteObject(&hmap_alphat, "hmap_alphat");
  file_out->WriteObject(&hmap_phit, "hmap_phit");
  file_out->WriteObject(&hmap_enucal, "hmap_enucal");
  file_out->WriteObject(&hmap_nhits_leadingp, "hmap_nhits_leadingp");
  file_out->WriteObject(&hmap_pmult, "hmap_pmult");
  file_out->WriteObject(&hmap_multpfp, "hmap_multpfp");
  file_out->WriteObject(&hmap_multtracktol, "hmap_multtracktol");

  // Mc truth stacked in interaction type
  file_out->WriteObject(&hmap_mctruth_nuenergy, "hmap_mctruth_nuenergy");
  file_out->WriteObject(&hmap_mctruth_mumom, "hmap_mctruth_mumom");
  file_out->WriteObject(&hmap_mctruth_mucostheta, "hmap_mctruth_mucostheta");
  file_out->WriteObject(&hmap_mctruth_muphi, "hmap_mctruth_muphi");
  file_out->WriteObject(&hmap_mctruth_chargedmult, "hmap_mctruth_chargedmult");
  file_out->WriteObject(&hmap_mctruth_mucostheta_mumom, "hmap_mctruth_mucostheta_mumom");
  file_out->WriteObject(&hmap_mctruth_nuenergy_gen, "hmap_mctruth_nuenergy_gen");
  file_out->WriteObject(&hmap_mctruth_mumom_gen, "hmap_mctruth_mumom_gen");
  file_out->WriteObject(&hmap_mctruth_mucostheta_gen, "hmap_mctruth_mucostheta_gen");
  file_out->WriteObject(&hmap_mctruth_muphi_gen, "hmap_mctruth_muphi_gen");
  file_out->WriteObject(&hmap_mctruth_chargedmult_gen, "hmap_mctruth_chargedmult_gen");
  file_out->WriteObject(&hmap_mctruth_mucostheta_mumom_gen, "hmap_mctruth_mucostheta_mumom_gen");

  file_out->WriteObject(&hmap_mctruth_pmom, "hmap_mctruth_pmom");
  file_out->WriteObject(&hmap_mctruth_pcostheta, "hmap_mctruth_pcostheta");
  file_out->WriteObject(&hmap_mctruth_pphi, "hmap_mctruth_pphi");
  file_out->WriteObject(&hmap_mctruth_thetamup, "hmap_mctruth_thetamup"); 


  // Efficiency - GENIE pm1sigma
  file_out->WriteObject(&bs_genie_pm1_eff_mumom_num, "bs_genie_pm1_eff_mumom_num");
  file_out->WriteObject(&bs_genie_pm1_eff_mumom_den, "bs_genie_pm1_eff_mumom_den");

  file_out->WriteObject(&bs_genie_pm1_eff_pmom_num, "bs_genie_pm1_eff_pmom_num");
  file_out->WriteObject(&bs_genie_pm1_eff_pmom_den, "bs_genie_pm1_eff_pmom_den");

  file_out->WriteObject(&bs_genie_pm1_eff_muangle_num, "bs_genie_pm1_eff_muangle_num");
  file_out->WriteObject(&bs_genie_pm1_eff_muangle_den, "bs_genie_pm1_eff_muangle_den");

  file_out->WriteObject(&bs_genie_pm1_eff_pangle_num, "bs_genie_pm1_eff_pangle_num");
  file_out->WriteObject(&bs_genie_pm1_eff_pangle_den, "bs_genie_pm1_eff_pangle_den");

  file_out->WriteObject(&bs_genie_pm1_eff_thetamup_num, "bs_genie_pm1_eff_thetamup_num");
  file_out->WriteObject(&bs_genie_pm1_eff_thetamup_den, "bs_genie_pm1_eff_thetamup_den");


  // All MC Histo - GENIE pm1sigma
  /*file_out->WriteObject(&hmap_trkmom_genie_pm1_bs, "hmap_trkmom_genie_pm1_bs");
  file_out->WriteObject(&hmap_trktheta_genie_pm1_bs, "hmap_trktheta_genie_pm1_bs");
  file_out->WriteObject(&hmap_trkpmom_genie_pm1_bs, "hmap_trkpmom_genie_pm1_bs");
  file_out->WriteObject(&hmap_trkptheta_genie_pm1_bs, "hmap_trkptheta_genie_pm1_bs");
  file_out->WriteObject(&hmap_thetamup_genie_pm1_bs, "hmap_thetamup_genie_pm1_bs");

  // Reco-True - GENIE pm1sigma
  file_out->WriteObject(&bs_genie_pm1_true_reco_mom, "bs_genie_pm1_true_reco_mom");
  file_out->WriteObject(&bs_genie_pm1_true_reco_pmom, "bs_genie_pm1_true_reco_pmom");
  file_out->WriteObject(&bs_genie_pm1_true_reco_muangle, "bs_genie_pm1_true_reco_muangle");
  file_out->WriteObject(&bs_genie_pm1_true_reco_pangle, "bs_genie_pm1_true_reco_pangle");
  file_out->WriteObject(&bs_genie_pm1_true_reco_thetamup, "bs_genie_pm1_true_reco_thetamup");
  */


 



  file_out->WriteObject(&hmap_vtxcheck_angle, "hmap_vtxcheck_angle");
  file_out->WriteObject(&hmap_residuals_std, "hmap_residuals_std");
  file_out->WriteObject(&hmap_residuals_mean, "hmap_residuals_mean");
  file_out->WriteObject(&hmap_perc_used_hits, "hmap_perc_used_hits");
  file_out->WriteObject(&hmap_mom_mcs_length, "hmap_mom_mcs_length");

  file_out->WriteObject(&hmap_xdiff_b, "hmap_xdiff_b");
  file_out->WriteObject(&hmap_zdiff_b, "hmap_zdiff_b");
  file_out->WriteObject(&hmap_xdiff, "hmap_xdiff");
  file_out->WriteObject(&hmap_zdiff, "hmap_zdiff");
  file_out->WriteObject(&hmap_pediff, "hmap_pediff");
  
  file_out->WriteObject(&hmap_vtxx_b, "hmap_vtxx_b");
  file_out->WriteObject(&hmap_vtxx, "hmap_vtxx");
  file_out->WriteObject(&hmap_vtxy, "hmap_vtxy");
  file_out->WriteObject(&hmap_vtxz, "hmap_vtxz");
  h_vtx_xz->Write();
  h_vtx_xy->Write();
  file_out->WriteObject(&hmap_vtxz_upborder, "hmap_vtxz_upborder");
  file_out->WriteObject(&hmap_vtxx_upborder, "hmap_vtxx_upborder");
  
  file_out->WriteObject(&hmap_dqdx_trunc, "hmap_dqdx_trunc");
  h_dqdx_trunc_length->Write();
  
  file_out->WriteObject(&hmap_ntpcobj, "hmap_ntpcobj");

  file_out->WriteObject(&hmap_flsmatch_score, "hmap_flsmatch_score");
  file_out->WriteObject(&hmap_flsmatch_score_second, "hmap_flsmatch_score_second");
  file_out->WriteObject(&hmap_flsmatch_score_difference, "hmap_flsmatch_score_difference");

  h_trklen_first->Write();
  h_trklen_second->Write();

  h_flsTime->Write();
  h_flsTime_wcut->Write();
  h_flsTime_wcut_2->Write();
  h_flsTime_wcut_3->Write();
  h_flsTime_wcut_4->Write();
  h_flsTime_wcut_5->Write();
  h_flsTime_wcut_6->Write();
  h_flsTime_wcut_7->Write();
  h_flsTime_wcut_8->Write();
  
  h_flsPe_wcut->Write();
  h_flsTime_flsPe_wcut->Write();
  
  h_deltax->Write();
  h_deltax_2d->Write();
  h_deltaz_4->Write();
  h_deltaz_6->Write();

  _true_reco_tree->Write();




  LOG_NORMAL() << "Saving 1D Event Histo." << std::endl;
  
  file_out->WriteObject(_event_histo_1d, "UBXSecEventHisto1D");

  LOG_NORMAL() << "1D Event Histo saved." << std::endl;

 


  LOG_NORMAL() << "Saving 2D Event Histo." << std::endl;
 
  file_out->WriteObject(_event_histo, "UBXSecEventHisto");

  LOG_NORMAL() << "2D Event Histo saved." << std::endl;
  




  file_out->Write();

  LOG_NORMAL() << "All saved." << std::endl;
  
  file_out->Close();

  LOG_NORMAL() << "Output file closed." << std::endl;





  //outfile.close();

  
  
  // Computing time
  clock_t end = clock();
  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
  std::cout << std::endl << std::endl;
  LOG_NORMAL() << "Computing time: " << elapsed_secs << " seconds = " << elapsed_secs/60. << " minutes." << endl << endl;
  
  //rootapp->Run();
  //rootapp->Terminate(0);
  
  return;
}

#endif
