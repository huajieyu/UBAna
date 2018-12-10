//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Dec  6 11:40:57 2018 by ROOT version 6.06/08
// from TTree tree/
// found on file: /build/kirby/cc1muNp_ubxsec_pid_integration_test_larana/test_ntuples_Dec05/ubxsec_output_mc_bnbcosmic.root
//////////////////////////////////////////////////////////

#ifndef UBXSecEvent_h
#define UBXSecEvent_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

using namespace std;

// Header file for the classes stored in the TTree if any.

class UBXSecEvent {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
 //UBXSecEvent     *ubxsec_event_split;
   Int_t           run;
   Int_t           subrun;
   Int_t           event;
   TString         file_type;
   Bool_t          muon_is_reco;
   Double_t        muon_reco_pur;
   Double_t        muon_reco_eff;
   Double_t        true_muon_mom;
   Double_t        true_muon_mom_matched;
   Int_t           n_pfp;
   Int_t           n_pfp_primary;
   Int_t           n_primary_cosmic_pfp;
   Int_t           nPFPtagged;
   Int_t           muon_is_flash_tagged;
   Double_t        muon_tag_score;
   Double_t        fm_score;
   Int_t           fv;
   Int_t           fv_sce;
   Int_t           ccnc;
   Int_t           mode;
   Int_t           nupdg;
   Bool_t          is_signal;
   Double_t        nu_e;
   Double_t        lep_costheta;
   Double_t        lep_phi;
   Int_t           genie_mult;
   Int_t           genie_mult_ch;
   Int_t           ngenie_muons;
   Int_t           ngenie_protons;
   Int_t           ngenie_electrons;
   Int_t           ngenie_pipms;
   Int_t           ngenie_pion0s;
   Int_t           ngenie_protons_200;
   Int_t           ngenie_protons_300;
   Int_t           ngenie_protons_400;
   vector<double>  genie_mcpar_pdgcode;
   vector<double>  genie_mcpar_energy;
   vector<double>  genie_mcpar_px;
   vector<double>  genie_mcpar_py;
   vector<double>  genie_mcpar_pz;
   vector<double>  genie_mcpar_startx;
   vector<double>  genie_mcpar_starty;
   vector<double>  genie_mcpar_startz;
   Double_t        genie_mcpar_W;
   Double_t        genie_mcpar_QSqr;
   vector<double>  geant_mcpar_pdgcode;
   vector<double>  geant_mcpar_energy;
   vector<double>  geant_mcpar_px;
   vector<double>  geant_mcpar_py;
   vector<double>  geant_mcpar_pz;
   vector<double>  geant_mcpar_startx;
   vector<double>  geant_mcpar_starty;
   vector<double>  geant_mcpar_startz;
   vector<string>  geant_mcpar_end_process;
   Double_t        bnb_weight;
   Bool_t          is_selected;
   Int_t           selected_slice;
   Double_t        sce_corr_x;
   Double_t        sce_corr_y;
   Double_t        sce_corr_z;
   Int_t           mc_muon_contained;
   Int_t           is_swtriggered;
   Double_t        vtx_resolution;
   Int_t           n_tpcobj_nu_origin;
   Int_t           n_tpcobj_cosmic_origin;
   Int_t           nslices;
   vector<double>  slc_flsmatch_score;
   vector<double>  slc_flsmatch_qllx;
   vector<double>  slc_flsmatch_tpcx;
   vector<double>  slc_flsmatch_t0;
   vector<double>  slc_flsmatch_hypoz;
   vector<double>  slc_flsmatch_xfixed_chi2;
   vector<double>  slc_flsmatch_xfixed_ll;
   vector<double>  slc_flsmatch_cosmic_score;
   vector<double>  slc_flsmatch_cosmic_t0;
   vector<double>  slc_nuvtx_x;
   vector<double>  slc_nuvtx_y;
   vector<double>  slc_nuvtx_z;
   vector<int>     slc_nuvtx_fv;
   vector<double>  slc_vtxcheck_angle;
   vector<int>     slc_origin;
   vector<int>     slc_origin_extra;
   vector<int>     slc_nhits_u;
   vector<int>     slc_nhits_v;
   vector<int>     slc_nhits_w;
   vector<double>  slc_longesttrack_length;
   vector<double>  slc_longesttrack_phi;
   vector<double>  slc_longesttrack_theta;
   vector<bool>    slc_longesttrack_iscontained;
   vector<double>  slc_longestshower_length;
   vector<double>  slc_longestshower_phi;
   vector<double>  slc_longestshower_theta;
   vector<double>  slc_longestshower_openangle;
   vector<double>  slc_longestshower_startx;
   vector<double>  slc_longestshower_starty;
   vector<double>  slc_longestshower_startz;
   vector<int>     slc_acpt_outoftime;
   vector<int>     slc_crosses_top_boundary;
   vector<int>     slc_nuvtx_closetodeadregion_u;
   vector<int>     slc_nuvtx_closetodeadregion_v;
   vector<int>     slc_nuvtx_closetodeadregion_w;
   vector<double>  slc_kalman_chi2;
   vector<int>     slc_kalman_ndof;
   vector<bool>    slc_passed_min_track_quality;
   vector<bool>    slc_passed_min_vertex_quality;
   vector<double>  slc_n_intime_pe_closestpmt;
   vector<double>  slc_maxdistance_vtxtrack;
   vector<bool>    slc_geocosmictag;
   vector<bool>    slc_consistency;
   vector<double>  slc_consistency_score;
   vector<int>     slc_npfp;
   vector<int>     slc_ntrack;
   vector<int>     slc_nshower;
   vector<bool>    slc_iscontained;
   vector<int>     slc_mult_pfp;
   vector<int>     slc_mult_track;
   vector<int>     slc_mult_shower;
   vector<int>     slc_mult_track_tolerance;
   vector<bool>    slc_muoncandidate_exists;
   vector<double>  slc_muoncandidate_length;
   vector<double>  slc_muoncandidate_phi;
   vector<double>  slc_muoncandidate_theta;
   vector<double>  slc_muoncandidate_mom_range;
   vector<double>  slc_muoncandidate_mom_mcs;
   vector<double>  slc_muoncandidate_mom_mcs_pi;
   vector<double>  slc_muoncandidate_mcs_ll;
   vector<bool>    slc_muoncandidate_contained;
   vector<double>  slc_muoncandidate_dqdx_trunc;
   vector<double>  slc_muoncandidate_dqdx_u_trunc;
   vector<double>  slc_muoncandidate_dqdx_v_trunc;
   vector<vector<double> > slc_muoncandidate_dqdx_v;
   vector<bool>    slc_muoncandidate_mip_consistency;
   vector<bool>    slc_muoncandidate_mip_consistency2;
   vector<int>     slc_muoncandidate_truepdg;
   vector<int>     slc_muoncandidate_trueorigin;
   vector<double>  slc_muoncandidate_mcs_delta_ll;
   vector<double>  slc_muoncandidate_residuals_mean;
   vector<double>  slc_muoncandidate_residuals_std;
   vector<int>     slc_muoncandidate_wiregap;
   vector<int>     slc_muoncandidate_wiregap_dead;
   vector<double>  slc_muoncandidate_linearity;
   vector<double>  slc_muoncandidate_perc_used_hits_in_cluster;
   vector<double>  slc_muoncandidate_maxscatteringangle;
   vector<double>  slc_muoncandidate_truth_origin;
   vector<double>  slc_muoncandidate_truth_pdg;
   vector<double>  slc_muoncandidate_truth_time;
   vector<double>  slc_muoncandidate_truth_startx;
   vector<double>  slc_muoncandidate_truth_starty;
   vector<double>  slc_muoncandidate_truth_startz;
   vector<double>  slc_muoncandidate_truth_endx;
   vector<double>  slc_muoncandidate_truth_endy;
   vector<double>  slc_muoncandidate_truth_endz;
   vector<double>  slc_muoncandidate_truth_px;
   vector<double>  slc_muoncandidate_truth_py;
   vector<double>  slc_muoncandidate_truth_pz;
   Int_t           nbeamfls;
   vector<double>  beamfls_time;
   vector<double>  beamfls_pe;
   vector<double>  beamfls_z;
   Int_t           candidate_flash_time;
   Double_t        candidate_flash_z;
   Bool_t          no_mcflash_but_op_activity;
   vector<vector<double> > beamfls_spec;
   vector<double>  numc_flash_spec;
   vector<vector<double> > slc_flshypo_xfixed_spec;
   vector<vector<double> > slc_flshypo_spec;
   Int_t           nsignal;
   vector<double>  tvtx_x;
   vector<double>  tvtx_y;
   vector<double>  tvtx_z;
   Double_t        pot;
   Int_t           evtwgt_genie_pm1_nfunc;
   vector<string>  evtwgt_genie_pm1_funcname;
   vector<int>     evtwgt_genie_pm1_nweight;
   vector<vector<double> > evtwgt_genie_pm1_weight;
   Int_t           evtwgt_genie_multisim_nfunc;
   vector<string>  evtwgt_genie_multisim_funcname;
   vector<int>     evtwgt_genie_multisim_nweight;
   vector<vector<double> > evtwgt_genie_multisim_weight;
   Int_t           evtwgt_extra_syst_multisim_nfunc;
   vector<string>  evtwgt_extra_syst_multisim_funcname;
   vector<int>     evtwgt_extra_syst_multisim_nweight;
   vector<vector<double> > evtwgt_extra_syst_multisim_weight;
   Int_t           evtwgt_flux_multisim_nfunc;
   vector<string>  evtwgt_flux_multisim_funcname;
   vector<int>     evtwgt_flux_multisim_nweight;
   vector<vector<double> > evtwgt_flux_multisim_weight;
   Int_t           num_pfp;
   Int_t           num_pfp_tracks;
   Int_t           num_pfp_showers;
   vector<int>     pfp_truth_pdg;
   vector<int>     pfp_truth_origin;
   vector<int>     pfp_truth_status;
   vector<int>     pfp_truth_parId;
   vector<float>   pfp_truth_theta;
   vector<float>   pfp_truth_costheta;
   vector<float>   pfp_truth_phi;
   vector<float>   pfp_truth_mom;
   vector<float>   pfp_truth_startx;
   vector<float>   pfp_truth_starty;
   vector<float>   pfp_truth_startz;
   vector<float>   pfp_truth_endx;
   vector<float>   pfp_truth_endy;
   vector<float>   pfp_truth_endz;
   vector<float>   pfp_truth_endE;
   vector<string>  pfp_truth_endProcess;
   vector<float>   pfp_truth_KE;
   vector<float>   pfp_truth_Mass;
   vector<int>     pfp_reco_isprimary;
   vector<int>     pfp_reco_ndaughters;
   vector<int>     pfp_reco_ismuoncandidate;
   vector<int>     pfp_reco_istrack;
   vector<int>     pfp_reco_numtracks;
   vector<int>     pfp_reco_isshower;
   vector<int>     pfp_reco_numshowers;
   vector<int>     pfp_reco_upflag;
   vector<int>     pfp_reco_Id;
   vector<float>   pfp_reco_length;
   vector<float>   pfp_reco_theta;
   vector<float>   pfp_reco_costheta;
   vector<float>   pfp_reco_phi;
   vector<float>   pfp_reco_startx;
   vector<float>   pfp_reco_starty;
   vector<float>   pfp_reco_startz;
   vector<float>   pfp_reco_endx;
   vector<float>   pfp_reco_endy;
   vector<float>   pfp_reco_endz;
   vector<float>   pfp_reco_Mom;
   vector<float>   pfp_reco_Mom_proton;
   vector<float>   pfp_reco_Mom_muon;
   vector<float>   pfp_reco_Mom_MCS;
   vector<float>   pfp_reco_trunmeandqdx;
   vector<float>   pfp_reco_trunmeandqdx_U;
   vector<float>   pfp_reco_trunmeandqdx_V;
   vector<float>   pfp_reco_newpid_pida;
   vector<float>   pfp_reco_chi2_proton;
   vector<float>   pfp_reco_chi2_kaon;
   vector<float>   pfp_reco_chi2_muon;
   vector<float>   pfp_reco_chi2_pion;
   vector<float>   pfp_reco_bragg_ratio;
   vector<float>   pfp_reco_bragg_proton;
   vector<float>   pfp_reco_bragg_fwd_proton;
   vector<float>   pfp_reco_bragg_bwd_proton;
   vector<float>   pfp_reco_bragg_fwd_mip;
   vector<int>     pfp_reco_nhits;
   vector<vector<double> > pfp_reco_dEdx;
   vector<vector<double> > pfp_reco_dQdx;
   vector<vector<double> > pfp_reco_RR;
   Int_t           _default_value;

   // List of branches
   TBranch        *b_ubxsec_event_split_run;   //!
   TBranch        *b_ubxsec_event_split_subrun;   //!
   TBranch        *b_ubxsec_event_split_event;   //!
   TBranch        *b_ubxsec_event_split_file_type;   //!
   TBranch        *b_ubxsec_event_split_muon_is_reco;   //!
   TBranch        *b_ubxsec_event_split_muon_reco_pur;   //!
   TBranch        *b_ubxsec_event_split_muon_reco_eff;   //!
   TBranch        *b_ubxsec_event_split_true_muon_mom;   //!
   TBranch        *b_ubxsec_event_split_true_muon_mom_matched;   //!
   TBranch        *b_ubxsec_event_split_n_pfp;   //!
   TBranch        *b_ubxsec_event_split_n_pfp_primary;   //!
   TBranch        *b_ubxsec_event_split_n_primary_cosmic_pfp;   //!
   TBranch        *b_ubxsec_event_split_nPFPtagged;   //!
   TBranch        *b_ubxsec_event_split_muon_is_flash_tagged;   //!
   TBranch        *b_ubxsec_event_split_muon_tag_score;   //!
   TBranch        *b_ubxsec_event_split_fm_score;   //!
   TBranch        *b_ubxsec_event_split_fv;   //!
   TBranch        *b_ubxsec_event_split_fv_sce;   //!
   TBranch        *b_ubxsec_event_split_ccnc;   //!
   TBranch        *b_ubxsec_event_split_mode;   //!
   TBranch        *b_ubxsec_event_split_nupdg;   //!
   TBranch        *b_ubxsec_event_split_is_signal;   //!
   TBranch        *b_ubxsec_event_split_nu_e;   //!
   TBranch        *b_ubxsec_event_split_lep_costheta;   //!
   TBranch        *b_ubxsec_event_split_lep_phi;   //!
   TBranch        *b_ubxsec_event_split_genie_mult;   //!
   TBranch        *b_ubxsec_event_split_genie_mult_ch;   //!
   TBranch        *b_ubxsec_event_split_ngenie_muons;   //!
   TBranch        *b_ubxsec_event_split_ngenie_protons;   //!
   TBranch        *b_ubxsec_event_split_ngenie_electrons;   //!
   TBranch        *b_ubxsec_event_split_ngenie_pipms;   //!
   TBranch        *b_ubxsec_event_split_ngenie_pion0s;   //!
   TBranch        *b_ubxsec_event_split_ngenie_protons_200;   //!
   TBranch        *b_ubxsec_event_split_ngenie_protons_300;   //!
   TBranch        *b_ubxsec_event_split_ngenie_protons_400;   //!
   TBranch        *b_ubxsec_event_split_genie_mcpar_pdgcode;   //!
   TBranch        *b_ubxsec_event_split_genie_mcpar_energy;   //!
   TBranch        *b_ubxsec_event_split_genie_mcpar_px;   //!
   TBranch        *b_ubxsec_event_split_genie_mcpar_py;   //!
   TBranch        *b_ubxsec_event_split_genie_mcpar_pz;   //!
   TBranch        *b_ubxsec_event_split_genie_mcpar_startx;   //!
   TBranch        *b_ubxsec_event_split_genie_mcpar_starty;   //!
   TBranch        *b_ubxsec_event_split_genie_mcpar_startz;   //!
   TBranch        *b_ubxsec_event_split_genie_mcpar_W;   //!
   TBranch        *b_ubxsec_event_split_genie_mcpar_QSqr;   //!
   TBranch        *b_ubxsec_event_split_geant_mcpar_pdgcode;   //!
   TBranch        *b_ubxsec_event_split_geant_mcpar_energy;   //!
   TBranch        *b_ubxsec_event_split_geant_mcpar_px;   //!
   TBranch        *b_ubxsec_event_split_geant_mcpar_py;   //!
   TBranch        *b_ubxsec_event_split_geant_mcpar_pz;   //!
   TBranch        *b_ubxsec_event_split_geant_mcpar_startx;   //!
   TBranch        *b_ubxsec_event_split_geant_mcpar_starty;   //!
   TBranch        *b_ubxsec_event_split_geant_mcpar_startz;   //!
   TBranch        *b_ubxsec_event_split_geant_mcpar_end_process;   //!
   TBranch        *b_ubxsec_event_split_bnb_weight;   //!
   TBranch        *b_ubxsec_event_split_is_selected;   //!
   TBranch        *b_ubxsec_event_split_selected_slice;   //!
   TBranch        *b_ubxsec_event_split_sce_corr_x;   //!
   TBranch        *b_ubxsec_event_split_sce_corr_y;   //!
   TBranch        *b_ubxsec_event_split_sce_corr_z;   //!
   TBranch        *b_ubxsec_event_split_mc_muon_contained;   //!
   TBranch        *b_ubxsec_event_split_is_swtriggered;   //!
   TBranch        *b_ubxsec_event_split_vtx_resolution;   //!
   TBranch        *b_ubxsec_event_split_n_tpcobj_nu_origin;   //!
   TBranch        *b_ubxsec_event_split_n_tpcobj_cosmic_origin;   //!
   TBranch        *b_ubxsec_event_split_nslices;   //!
   TBranch        *b_ubxsec_event_split_slc_flsmatch_score;   //!
   TBranch        *b_ubxsec_event_split_slc_flsmatch_qllx;   //!
   TBranch        *b_ubxsec_event_split_slc_flsmatch_tpcx;   //!
   TBranch        *b_ubxsec_event_split_slc_flsmatch_t0;   //!
   TBranch        *b_ubxsec_event_split_slc_flsmatch_hypoz;   //!
   TBranch        *b_ubxsec_event_split_slc_flsmatch_xfixed_chi2;   //!
   TBranch        *b_ubxsec_event_split_slc_flsmatch_xfixed_ll;   //!
   TBranch        *b_ubxsec_event_split_slc_flsmatch_cosmic_score;   //!
   TBranch        *b_ubxsec_event_split_slc_flsmatch_cosmic_t0;   //!
   TBranch        *b_ubxsec_event_split_slc_nuvtx_x;   //!
   TBranch        *b_ubxsec_event_split_slc_nuvtx_y;   //!
   TBranch        *b_ubxsec_event_split_slc_nuvtx_z;   //!
   TBranch        *b_ubxsec_event_split_slc_nuvtx_fv;   //!
   TBranch        *b_ubxsec_event_split_slc_vtxcheck_angle;   //!
   TBranch        *b_ubxsec_event_split_slc_origin;   //!
   TBranch        *b_ubxsec_event_split_slc_origin_extra;   //!
   TBranch        *b_ubxsec_event_split_slc_nhits_u;   //!
   TBranch        *b_ubxsec_event_split_slc_nhits_v;   //!
   TBranch        *b_ubxsec_event_split_slc_nhits_w;   //!
   TBranch        *b_ubxsec_event_split_slc_longesttrack_length;   //!
   TBranch        *b_ubxsec_event_split_slc_longesttrack_phi;   //!
   TBranch        *b_ubxsec_event_split_slc_longesttrack_theta;   //!
   TBranch        *b_ubxsec_event_split_slc_longesttrack_iscontained;   //!
   TBranch        *b_ubxsec_event_split_slc_longestshower_length;   //!
   TBranch        *b_ubxsec_event_split_slc_longestshower_phi;   //!
   TBranch        *b_ubxsec_event_split_slc_longestshower_theta;   //!
   TBranch        *b_ubxsec_event_split_slc_longestshower_openangle;   //!
   TBranch        *b_ubxsec_event_split_slc_longestshower_startx;   //!
   TBranch        *b_ubxsec_event_split_slc_longestshower_starty;   //!
   TBranch        *b_ubxsec_event_split_slc_longestshower_startz;   //!
   TBranch        *b_ubxsec_event_split_slc_acpt_outoftime;   //!
   TBranch        *b_ubxsec_event_split_slc_crosses_top_boundary;   //!
   TBranch        *b_ubxsec_event_split_slc_nuvtx_closetodeadregion_u;   //!
   TBranch        *b_ubxsec_event_split_slc_nuvtx_closetodeadregion_v;   //!
   TBranch        *b_ubxsec_event_split_slc_nuvtx_closetodeadregion_w;   //!
   TBranch        *b_ubxsec_event_split_slc_kalman_chi2;   //!
   TBranch        *b_ubxsec_event_split_slc_kalman_ndof;   //!
   TBranch        *b_ubxsec_event_split_slc_passed_min_track_quality;   //!
   TBranch        *b_ubxsec_event_split_slc_passed_min_vertex_quality;   //!
   TBranch        *b_ubxsec_event_split_slc_n_intime_pe_closestpmt;   //!
   TBranch        *b_ubxsec_event_split_slc_maxdistance_vtxtrack;   //!
   TBranch        *b_ubxsec_event_split_slc_geocosmictag;   //!
   TBranch        *b_ubxsec_event_split_slc_consistency;   //!
   TBranch        *b_ubxsec_event_split_slc_consistency_score;   //!
   TBranch        *b_ubxsec_event_split_slc_npfp;   //!
   TBranch        *b_ubxsec_event_split_slc_ntrack;   //!
   TBranch        *b_ubxsec_event_split_slc_nshower;   //!
   TBranch        *b_ubxsec_event_split_slc_iscontained;   //!
   TBranch        *b_ubxsec_event_split_slc_mult_pfp;   //!
   TBranch        *b_ubxsec_event_split_slc_mult_track;   //!
   TBranch        *b_ubxsec_event_split_slc_mult_shower;   //!
   TBranch        *b_ubxsec_event_split_slc_mult_track_tolerance;   //!
   TBranch        *b_ubxsec_event_split_slc_muoncandidate_exists;   //!
   TBranch        *b_ubxsec_event_split_slc_muoncandidate_length;   //!
   TBranch        *b_ubxsec_event_split_slc_muoncandidate_phi;   //!
   TBranch        *b_ubxsec_event_split_slc_muoncandidate_theta;   //!
   TBranch        *b_ubxsec_event_split_slc_muoncandidate_mom_range;   //!
   TBranch        *b_ubxsec_event_split_slc_muoncandidate_mom_mcs;   //!
   TBranch        *b_ubxsec_event_split_slc_muoncandidate_mom_mcs_pi;   //!
   TBranch        *b_ubxsec_event_split_slc_muoncandidate_mcs_ll;   //!
   TBranch        *b_ubxsec_event_split_slc_muoncandidate_contained;   //!
   TBranch        *b_ubxsec_event_split_slc_muoncandidate_dqdx_trunc;   //!
   TBranch        *b_ubxsec_event_split_slc_muoncandidate_dqdx_u_trunc;   //!
   TBranch        *b_ubxsec_event_split_slc_muoncandidate_dqdx_v_trunc;   //!
   TBranch        *b_ubxsec_event_split_slc_muoncandidate_dqdx_v;   //!
   TBranch        *b_ubxsec_event_split_slc_muoncandidate_mip_consistency;   //!
   TBranch        *b_ubxsec_event_split_slc_muoncandidate_mip_consistency2;   //!
   TBranch        *b_ubxsec_event_split_slc_muoncandidate_truepdg;   //!
   TBranch        *b_ubxsec_event_split_slc_muoncandidate_trueorigin;   //!
   TBranch        *b_ubxsec_event_split_slc_muoncandidate_mcs_delta_ll;   //!
   TBranch        *b_ubxsec_event_split_slc_muoncandidate_residuals_mean;   //!
   TBranch        *b_ubxsec_event_split_slc_muoncandidate_residuals_std;   //!
   TBranch        *b_ubxsec_event_split_slc_muoncandidate_wiregap;   //!
   TBranch        *b_ubxsec_event_split_slc_muoncandidate_wiregap_dead;   //!
   TBranch        *b_ubxsec_event_split_slc_muoncandidate_linearity;   //!
   TBranch        *b_ubxsec_event_split_slc_muoncandidate_perc_used_hits_in_cluster;   //!
   TBranch        *b_ubxsec_event_split_slc_muoncandidate_maxscatteringangle;   //!
   TBranch        *b_ubxsec_event_split_slc_muoncandidate_truth_origin;   //!
   TBranch        *b_ubxsec_event_split_slc_muoncandidate_truth_pdg;   //!
   TBranch        *b_ubxsec_event_split_slc_muoncandidate_truth_time;   //!
   TBranch        *b_ubxsec_event_split_slc_muoncandidate_truth_startx;   //!
   TBranch        *b_ubxsec_event_split_slc_muoncandidate_truth_starty;   //!
   TBranch        *b_ubxsec_event_split_slc_muoncandidate_truth_startz;   //!
   TBranch        *b_ubxsec_event_split_slc_muoncandidate_truth_endx;   //!
   TBranch        *b_ubxsec_event_split_slc_muoncandidate_truth_endy;   //!
   TBranch        *b_ubxsec_event_split_slc_muoncandidate_truth_endz;   //!
   TBranch        *b_ubxsec_event_split_slc_muoncandidate_truth_px;   //!
   TBranch        *b_ubxsec_event_split_slc_muoncandidate_truth_py;   //!
   TBranch        *b_ubxsec_event_split_slc_muoncandidate_truth_pz;   //!
   TBranch        *b_ubxsec_event_split_nbeamfls;   //!
   TBranch        *b_ubxsec_event_split_beamfls_time;   //!
   TBranch        *b_ubxsec_event_split_beamfls_pe;   //!
   TBranch        *b_ubxsec_event_split_beamfls_z;   //!
   TBranch        *b_ubxsec_event_split_candidate_flash_time;   //!
   TBranch        *b_ubxsec_event_split_candidate_flash_z;   //!
   TBranch        *b_ubxsec_event_split_no_mcflash_but_op_activity;   //!
   TBranch        *b_ubxsec_event_split_beamfls_spec;   //!
   TBranch        *b_ubxsec_event_split_numc_flash_spec;   //!
   TBranch        *b_ubxsec_event_split_slc_flshypo_xfixed_spec;   //!
   TBranch        *b_ubxsec_event_split_slc_flshypo_spec;   //!
   TBranch        *b_ubxsec_event_split_nsignal;   //!
   TBranch        *b_ubxsec_event_split_tvtx_x;   //!
   TBranch        *b_ubxsec_event_split_tvtx_y;   //!
   TBranch        *b_ubxsec_event_split_tvtx_z;   //!
   TBranch        *b_ubxsec_event_split_pot;   //!
   TBranch        *b_ubxsec_event_split_evtwgt_genie_pm1_nfunc;   //!
   TBranch        *b_ubxsec_event_split_evtwgt_genie_pm1_funcname;   //!
   TBranch        *b_ubxsec_event_split_evtwgt_genie_pm1_nweight;   //!
   TBranch        *b_ubxsec_event_split_evtwgt_genie_pm1_weight;   //!
   TBranch        *b_ubxsec_event_split_evtwgt_genie_multisim_nfunc;   //!
   TBranch        *b_ubxsec_event_split_evtwgt_genie_multisim_funcname;   //!
   TBranch        *b_ubxsec_event_split_evtwgt_genie_multisim_nweight;   //!
   TBranch        *b_ubxsec_event_split_evtwgt_genie_multisim_weight;   //!
   TBranch        *b_ubxsec_event_split_evtwgt_extra_syst_multisim_nfunc;   //!
   TBranch        *b_ubxsec_event_split_evtwgt_extra_syst_multisim_funcname;   //!
   TBranch        *b_ubxsec_event_split_evtwgt_extra_syst_multisim_nweight;   //!
   TBranch        *b_ubxsec_event_split_evtwgt_extra_syst_multisim_weight;   //!
   TBranch        *b_ubxsec_event_split_evtwgt_flux_multisim_nfunc;   //!
   TBranch        *b_ubxsec_event_split_evtwgt_flux_multisim_funcname;   //!
   TBranch        *b_ubxsec_event_split_evtwgt_flux_multisim_nweight;   //!
   TBranch        *b_ubxsec_event_split_evtwgt_flux_multisim_weight;   //!
   TBranch        *b_ubxsec_event_split_num_pfp;   //!
   TBranch        *b_ubxsec_event_split_num_pfp_tracks;   //!
   TBranch        *b_ubxsec_event_split_num_pfp_showers;   //!
   TBranch        *b_ubxsec_event_split_pfp_truth_pdg;   //!
   TBranch        *b_ubxsec_event_split_pfp_truth_origin;   //!
   TBranch        *b_ubxsec_event_split_pfp_truth_status;   //!
   TBranch        *b_ubxsec_event_split_pfp_truth_parId;   //!
   TBranch        *b_ubxsec_event_split_pfp_truth_theta;   //!
   TBranch        *b_ubxsec_event_split_pfp_truth_costheta;   //!
   TBranch        *b_ubxsec_event_split_pfp_truth_phi;   //!
   TBranch        *b_ubxsec_event_split_pfp_truth_mom;   //!
   TBranch        *b_ubxsec_event_split_pfp_truth_startx;   //!
   TBranch        *b_ubxsec_event_split_pfp_truth_starty;   //!
   TBranch        *b_ubxsec_event_split_pfp_truth_startz;   //!
   TBranch        *b_ubxsec_event_split_pfp_truth_endx;   //!
   TBranch        *b_ubxsec_event_split_pfp_truth_endy;   //!
   TBranch        *b_ubxsec_event_split_pfp_truth_endz;   //!
   TBranch        *b_ubxsec_event_split_pfp_truth_endE;   //!
   TBranch        *b_ubxsec_event_split_pfp_truth_endProcess;   //!
   TBranch        *b_ubxsec_event_split_pfp_truth_KE;   //!
   TBranch        *b_ubxsec_event_split_pfp_truth_Mass;   //!
   TBranch        *b_ubxsec_event_split_pfp_reco_isprimary;   //!
   TBranch        *b_ubxsec_event_split_pfp_reco_ndaughters;   //!
   TBranch        *b_ubxsec_event_split_pfp_reco_ismuoncandidate;   //!
   TBranch        *b_ubxsec_event_split_pfp_reco_istrack;   //!
   TBranch        *b_ubxsec_event_split_pfp_reco_numtracks;   //!
   TBranch        *b_ubxsec_event_split_pfp_reco_isshower;   //!
   TBranch        *b_ubxsec_event_split_pfp_reco_numshowers;   //!
   TBranch        *b_ubxsec_event_split_pfp_reco_upflag;   //!
   TBranch        *b_ubxsec_event_split_pfp_reco_Id;   //!
   TBranch        *b_ubxsec_event_split_pfp_reco_length;   //!
   TBranch        *b_ubxsec_event_split_pfp_reco_theta;   //!
   TBranch        *b_ubxsec_event_split_pfp_reco_costheta;   //!
   TBranch        *b_ubxsec_event_split_pfp_reco_phi;   //!
   TBranch        *b_ubxsec_event_split_pfp_reco_startx;   //!
   TBranch        *b_ubxsec_event_split_pfp_reco_starty;   //!
   TBranch        *b_ubxsec_event_split_pfp_reco_startz;   //!
   TBranch        *b_ubxsec_event_split_pfp_reco_endx;   //!
   TBranch        *b_ubxsec_event_split_pfp_reco_endy;   //!
   TBranch        *b_ubxsec_event_split_pfp_reco_endz;   //!
   TBranch        *b_ubxsec_event_split_pfp_reco_Mom;   //!
   TBranch        *b_ubxsec_event_split_pfp_reco_Mom_proton;   //!
   TBranch        *b_ubxsec_event_split_pfp_reco_Mom_muon;   //!
   TBranch        *b_ubxsec_event_split_pfp_reco_Mom_MCS;   //!
   TBranch        *b_ubxsec_event_split_pfp_reco_trunmeandqdx;   //!
   TBranch        *b_ubxsec_event_split_pfp_reco_trunmeandqdx_U;   //!
   TBranch        *b_ubxsec_event_split_pfp_reco_trunmeandqdx_V;   //!
   TBranch        *b_ubxsec_event_split_pfp_reco_newpid_pida;   //!
   TBranch        *b_ubxsec_event_split_pfp_reco_chi2_proton;   //!
   TBranch        *b_ubxsec_event_split_pfp_reco_chi2_kaon;   //!
   TBranch        *b_ubxsec_event_split_pfp_reco_chi2_muon;   //!
   TBranch        *b_ubxsec_event_split_pfp_reco_chi2_pion;   //!
   TBranch        *b_ubxsec_event_split_pfp_reco_bragg_ratio;   //!
   TBranch        *b_ubxsec_event_split_pfp_reco_bragg_proton;   //!
   TBranch        *b_ubxsec_event_split_pfp_reco_bragg_fwd_proton;   //!
   TBranch        *b_ubxsec_event_split_pfp_reco_bragg_bwd_proton;   //!
   TBranch        *b_ubxsec_event_split_pfp_reco_bragg_fwd_mip;   //!
   TBranch        *b_ubxsec_event_split_pfp_reco_nhits;   //!
   TBranch        *b_ubxsec_event_split_pfp_reco_dEdx;   //!
   TBranch        *b_ubxsec_event_split_pfp_reco_dQdx;   //!
   TBranch        *b_ubxsec_event_split_pfp_reco_RR;   //!
   TBranch        *b_ubxsec_event_split__default_value;   //!

   UBXSecEvent(TTree *tree=0);
   virtual ~UBXSecEvent();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef UBXSecEvent_cxx
UBXSecEvent::UBXSecEvent(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/build/kirby/cc1muNp_ubxsec_pid_integration_test_larana/test_ntuples_Dec05/ubxsec_output_mc_bnbcosmic.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/build/kirby/cc1muNp_ubxsec_pid_integration_test_larana/test_ntuples_Dec05/ubxsec_output_mc_bnbcosmic.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("/build/kirby/cc1muNp_ubxsec_pid_integration_test_larana/test_ntuples_Dec05/ubxsec_output_mc_bnbcosmic.root:/UBXSec");
      dir->GetObject("tree",tree);

   }
   Init(tree);
}

UBXSecEvent::~UBXSecEvent()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t UBXSecEvent::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t UBXSecEvent::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void UBXSecEvent::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_ubxsec_event_split_run);
   fChain->SetBranchAddress("subrun", &subrun, &b_ubxsec_event_split_subrun);
   fChain->SetBranchAddress("event", &event, &b_ubxsec_event_split_event);
   fChain->SetBranchAddress("file_type", &file_type, &b_ubxsec_event_split_file_type);
   fChain->SetBranchAddress("muon_is_reco", &muon_is_reco, &b_ubxsec_event_split_muon_is_reco);
   fChain->SetBranchAddress("muon_reco_pur", &muon_reco_pur, &b_ubxsec_event_split_muon_reco_pur);
   fChain->SetBranchAddress("muon_reco_eff", &muon_reco_eff, &b_ubxsec_event_split_muon_reco_eff);
   fChain->SetBranchAddress("true_muon_mom", &true_muon_mom, &b_ubxsec_event_split_true_muon_mom);
   fChain->SetBranchAddress("true_muon_mom_matched", &true_muon_mom_matched, &b_ubxsec_event_split_true_muon_mom_matched);
   fChain->SetBranchAddress("n_pfp", &n_pfp, &b_ubxsec_event_split_n_pfp);
   fChain->SetBranchAddress("n_pfp_primary", &n_pfp_primary, &b_ubxsec_event_split_n_pfp_primary);
   fChain->SetBranchAddress("n_primary_cosmic_pfp", &n_primary_cosmic_pfp, &b_ubxsec_event_split_n_primary_cosmic_pfp);
   fChain->SetBranchAddress("nPFPtagged", &nPFPtagged, &b_ubxsec_event_split_nPFPtagged);
   fChain->SetBranchAddress("muon_is_flash_tagged", &muon_is_flash_tagged, &b_ubxsec_event_split_muon_is_flash_tagged);
   fChain->SetBranchAddress("muon_tag_score", &muon_tag_score, &b_ubxsec_event_split_muon_tag_score);
   fChain->SetBranchAddress("fm_score", &fm_score, &b_ubxsec_event_split_fm_score);
   fChain->SetBranchAddress("fv", &fv, &b_ubxsec_event_split_fv);
   fChain->SetBranchAddress("fv_sce", &fv_sce, &b_ubxsec_event_split_fv_sce);
   fChain->SetBranchAddress("ccnc", &ccnc, &b_ubxsec_event_split_ccnc);
   fChain->SetBranchAddress("mode", &mode, &b_ubxsec_event_split_mode);
   fChain->SetBranchAddress("nupdg", &nupdg, &b_ubxsec_event_split_nupdg);
   fChain->SetBranchAddress("is_signal", &is_signal, &b_ubxsec_event_split_is_signal);
   fChain->SetBranchAddress("nu_e", &nu_e, &b_ubxsec_event_split_nu_e);
   fChain->SetBranchAddress("lep_costheta", &lep_costheta, &b_ubxsec_event_split_lep_costheta);
   fChain->SetBranchAddress("lep_phi", &lep_phi, &b_ubxsec_event_split_lep_phi);
   fChain->SetBranchAddress("genie_mult", &genie_mult, &b_ubxsec_event_split_genie_mult);
   fChain->SetBranchAddress("genie_mult_ch", &genie_mult_ch, &b_ubxsec_event_split_genie_mult_ch);
   fChain->SetBranchAddress("ngenie_muons", &ngenie_muons, &b_ubxsec_event_split_ngenie_muons);
   fChain->SetBranchAddress("ngenie_protons", &ngenie_protons, &b_ubxsec_event_split_ngenie_protons);
   fChain->SetBranchAddress("ngenie_electrons", &ngenie_electrons, &b_ubxsec_event_split_ngenie_electrons);
   fChain->SetBranchAddress("ngenie_pipms", &ngenie_pipms, &b_ubxsec_event_split_ngenie_pipms);
   fChain->SetBranchAddress("ngenie_pion0s", &ngenie_pion0s, &b_ubxsec_event_split_ngenie_pion0s);
   fChain->SetBranchAddress("ngenie_protons_200", &ngenie_protons_200, &b_ubxsec_event_split_ngenie_protons_200);
   fChain->SetBranchAddress("ngenie_protons_300", &ngenie_protons_300, &b_ubxsec_event_split_ngenie_protons_300);
   fChain->SetBranchAddress("ngenie_protons_400", &ngenie_protons_400, &b_ubxsec_event_split_ngenie_protons_400);
   fChain->SetBranchAddress("genie_mcpar_pdgcode", &genie_mcpar_pdgcode, &b_ubxsec_event_split_genie_mcpar_pdgcode);
   fChain->SetBranchAddress("genie_mcpar_energy", &genie_mcpar_energy, &b_ubxsec_event_split_genie_mcpar_energy);
   fChain->SetBranchAddress("genie_mcpar_px", &genie_mcpar_px, &b_ubxsec_event_split_genie_mcpar_px);
   fChain->SetBranchAddress("genie_mcpar_py", &genie_mcpar_py, &b_ubxsec_event_split_genie_mcpar_py);
   fChain->SetBranchAddress("genie_mcpar_pz", &genie_mcpar_pz, &b_ubxsec_event_split_genie_mcpar_pz);
   fChain->SetBranchAddress("genie_mcpar_startx", &genie_mcpar_startx, &b_ubxsec_event_split_genie_mcpar_startx);
   fChain->SetBranchAddress("genie_mcpar_starty", &genie_mcpar_starty, &b_ubxsec_event_split_genie_mcpar_starty);
   fChain->SetBranchAddress("genie_mcpar_startz", &genie_mcpar_startz, &b_ubxsec_event_split_genie_mcpar_startz);
   fChain->SetBranchAddress("genie_mcpar_W", &genie_mcpar_W, &b_ubxsec_event_split_genie_mcpar_W);
   fChain->SetBranchAddress("genie_mcpar_QSqr", &genie_mcpar_QSqr, &b_ubxsec_event_split_genie_mcpar_QSqr);
   fChain->SetBranchAddress("geant_mcpar_pdgcode", &geant_mcpar_pdgcode, &b_ubxsec_event_split_geant_mcpar_pdgcode);
   fChain->SetBranchAddress("geant_mcpar_energy", &geant_mcpar_energy, &b_ubxsec_event_split_geant_mcpar_energy);
   fChain->SetBranchAddress("geant_mcpar_px", &geant_mcpar_px, &b_ubxsec_event_split_geant_mcpar_px);
   fChain->SetBranchAddress("geant_mcpar_py", &geant_mcpar_py, &b_ubxsec_event_split_geant_mcpar_py);
   fChain->SetBranchAddress("geant_mcpar_pz", &geant_mcpar_pz, &b_ubxsec_event_split_geant_mcpar_pz);
   fChain->SetBranchAddress("geant_mcpar_startx", &geant_mcpar_startx, &b_ubxsec_event_split_geant_mcpar_startx);
   fChain->SetBranchAddress("geant_mcpar_starty", &geant_mcpar_starty, &b_ubxsec_event_split_geant_mcpar_starty);
   fChain->SetBranchAddress("geant_mcpar_startz", &geant_mcpar_startz, &b_ubxsec_event_split_geant_mcpar_startz);
   fChain->SetBranchAddress("geant_mcpar_end_process", &geant_mcpar_end_process, &b_ubxsec_event_split_geant_mcpar_end_process);
   fChain->SetBranchAddress("bnb_weight", &bnb_weight, &b_ubxsec_event_split_bnb_weight);
   fChain->SetBranchAddress("is_selected", &is_selected, &b_ubxsec_event_split_is_selected);
   fChain->SetBranchAddress("selected_slice", &selected_slice, &b_ubxsec_event_split_selected_slice);
   fChain->SetBranchAddress("sce_corr_x", &sce_corr_x, &b_ubxsec_event_split_sce_corr_x);
   fChain->SetBranchAddress("sce_corr_y", &sce_corr_y, &b_ubxsec_event_split_sce_corr_y);
   fChain->SetBranchAddress("sce_corr_z", &sce_corr_z, &b_ubxsec_event_split_sce_corr_z);
   fChain->SetBranchAddress("mc_muon_contained", &mc_muon_contained, &b_ubxsec_event_split_mc_muon_contained);
   fChain->SetBranchAddress("is_swtriggered", &is_swtriggered, &b_ubxsec_event_split_is_swtriggered);
   fChain->SetBranchAddress("vtx_resolution", &vtx_resolution, &b_ubxsec_event_split_vtx_resolution);
   fChain->SetBranchAddress("n_tpcobj_nu_origin", &n_tpcobj_nu_origin, &b_ubxsec_event_split_n_tpcobj_nu_origin);
   fChain->SetBranchAddress("n_tpcobj_cosmic_origin", &n_tpcobj_cosmic_origin, &b_ubxsec_event_split_n_tpcobj_cosmic_origin);
   fChain->SetBranchAddress("nslices", &nslices, &b_ubxsec_event_split_nslices);
   fChain->SetBranchAddress("slc_flsmatch_score", &slc_flsmatch_score, &b_ubxsec_event_split_slc_flsmatch_score);
   fChain->SetBranchAddress("slc_flsmatch_qllx", &slc_flsmatch_qllx, &b_ubxsec_event_split_slc_flsmatch_qllx);
   fChain->SetBranchAddress("slc_flsmatch_tpcx", &slc_flsmatch_tpcx, &b_ubxsec_event_split_slc_flsmatch_tpcx);
   fChain->SetBranchAddress("slc_flsmatch_t0", &slc_flsmatch_t0, &b_ubxsec_event_split_slc_flsmatch_t0);
   fChain->SetBranchAddress("slc_flsmatch_hypoz", &slc_flsmatch_hypoz, &b_ubxsec_event_split_slc_flsmatch_hypoz);
   fChain->SetBranchAddress("slc_flsmatch_xfixed_chi2", &slc_flsmatch_xfixed_chi2, &b_ubxsec_event_split_slc_flsmatch_xfixed_chi2);
   fChain->SetBranchAddress("slc_flsmatch_xfixed_ll", &slc_flsmatch_xfixed_ll, &b_ubxsec_event_split_slc_flsmatch_xfixed_ll);
   fChain->SetBranchAddress("slc_flsmatch_cosmic_score", &slc_flsmatch_cosmic_score, &b_ubxsec_event_split_slc_flsmatch_cosmic_score);
   fChain->SetBranchAddress("slc_flsmatch_cosmic_t0", &slc_flsmatch_cosmic_t0, &b_ubxsec_event_split_slc_flsmatch_cosmic_t0);
   fChain->SetBranchAddress("slc_nuvtx_x", &slc_nuvtx_x, &b_ubxsec_event_split_slc_nuvtx_x);
   fChain->SetBranchAddress("slc_nuvtx_y", &slc_nuvtx_y, &b_ubxsec_event_split_slc_nuvtx_y);
   fChain->SetBranchAddress("slc_nuvtx_z", &slc_nuvtx_z, &b_ubxsec_event_split_slc_nuvtx_z);
   fChain->SetBranchAddress("slc_nuvtx_fv", &slc_nuvtx_fv, &b_ubxsec_event_split_slc_nuvtx_fv);
   fChain->SetBranchAddress("slc_vtxcheck_angle", &slc_vtxcheck_angle, &b_ubxsec_event_split_slc_vtxcheck_angle);
   fChain->SetBranchAddress("slc_origin", &slc_origin, &b_ubxsec_event_split_slc_origin);
   fChain->SetBranchAddress("slc_origin_extra", &slc_origin_extra, &b_ubxsec_event_split_slc_origin_extra);
   fChain->SetBranchAddress("slc_nhits_u", &slc_nhits_u, &b_ubxsec_event_split_slc_nhits_u);
   fChain->SetBranchAddress("slc_nhits_v", &slc_nhits_v, &b_ubxsec_event_split_slc_nhits_v);
   fChain->SetBranchAddress("slc_nhits_w", &slc_nhits_w, &b_ubxsec_event_split_slc_nhits_w);
   fChain->SetBranchAddress("slc_longesttrack_length", &slc_longesttrack_length, &b_ubxsec_event_split_slc_longesttrack_length);
   fChain->SetBranchAddress("slc_longesttrack_phi", &slc_longesttrack_phi, &b_ubxsec_event_split_slc_longesttrack_phi);
   fChain->SetBranchAddress("slc_longesttrack_theta", &slc_longesttrack_theta, &b_ubxsec_event_split_slc_longesttrack_theta);
   fChain->SetBranchAddress("slc_longesttrack_iscontained", &slc_longesttrack_iscontained, &b_ubxsec_event_split_slc_longesttrack_iscontained);
   fChain->SetBranchAddress("slc_longestshower_length", &slc_longestshower_length, &b_ubxsec_event_split_slc_longestshower_length);
   fChain->SetBranchAddress("slc_longestshower_phi", &slc_longestshower_phi, &b_ubxsec_event_split_slc_longestshower_phi);
   fChain->SetBranchAddress("slc_longestshower_theta", &slc_longestshower_theta, &b_ubxsec_event_split_slc_longestshower_theta);
   fChain->SetBranchAddress("slc_longestshower_openangle", &slc_longestshower_openangle, &b_ubxsec_event_split_slc_longestshower_openangle);
   fChain->SetBranchAddress("slc_longestshower_startx", &slc_longestshower_startx, &b_ubxsec_event_split_slc_longestshower_startx);
   fChain->SetBranchAddress("slc_longestshower_starty", &slc_longestshower_starty, &b_ubxsec_event_split_slc_longestshower_starty);
   fChain->SetBranchAddress("slc_longestshower_startz", &slc_longestshower_startz, &b_ubxsec_event_split_slc_longestshower_startz);
   fChain->SetBranchAddress("slc_acpt_outoftime", &slc_acpt_outoftime, &b_ubxsec_event_split_slc_acpt_outoftime);
   fChain->SetBranchAddress("slc_crosses_top_boundary", &slc_crosses_top_boundary, &b_ubxsec_event_split_slc_crosses_top_boundary);
   fChain->SetBranchAddress("slc_nuvtx_closetodeadregion_u", &slc_nuvtx_closetodeadregion_u, &b_ubxsec_event_split_slc_nuvtx_closetodeadregion_u);
   fChain->SetBranchAddress("slc_nuvtx_closetodeadregion_v", &slc_nuvtx_closetodeadregion_v, &b_ubxsec_event_split_slc_nuvtx_closetodeadregion_v);
   fChain->SetBranchAddress("slc_nuvtx_closetodeadregion_w", &slc_nuvtx_closetodeadregion_w, &b_ubxsec_event_split_slc_nuvtx_closetodeadregion_w);
   fChain->SetBranchAddress("slc_kalman_chi2", &slc_kalman_chi2, &b_ubxsec_event_split_slc_kalman_chi2);
   fChain->SetBranchAddress("slc_kalman_ndof", &slc_kalman_ndof, &b_ubxsec_event_split_slc_kalman_ndof);
   fChain->SetBranchAddress("slc_passed_min_track_quality", &slc_passed_min_track_quality, &b_ubxsec_event_split_slc_passed_min_track_quality);
   fChain->SetBranchAddress("slc_passed_min_vertex_quality", &slc_passed_min_vertex_quality, &b_ubxsec_event_split_slc_passed_min_vertex_quality);
   fChain->SetBranchAddress("slc_n_intime_pe_closestpmt", &slc_n_intime_pe_closestpmt, &b_ubxsec_event_split_slc_n_intime_pe_closestpmt);
   fChain->SetBranchAddress("slc_maxdistance_vtxtrack", &slc_maxdistance_vtxtrack, &b_ubxsec_event_split_slc_maxdistance_vtxtrack);
   fChain->SetBranchAddress("slc_geocosmictag", &slc_geocosmictag, &b_ubxsec_event_split_slc_geocosmictag);
   fChain->SetBranchAddress("slc_consistency", &slc_consistency, &b_ubxsec_event_split_slc_consistency);
   fChain->SetBranchAddress("slc_consistency_score", &slc_consistency_score, &b_ubxsec_event_split_slc_consistency_score);
   fChain->SetBranchAddress("slc_npfp", &slc_npfp, &b_ubxsec_event_split_slc_npfp);
   fChain->SetBranchAddress("slc_ntrack", &slc_ntrack, &b_ubxsec_event_split_slc_ntrack);
   fChain->SetBranchAddress("slc_nshower", &slc_nshower, &b_ubxsec_event_split_slc_nshower);
   fChain->SetBranchAddress("slc_iscontained", &slc_iscontained, &b_ubxsec_event_split_slc_iscontained);
   fChain->SetBranchAddress("slc_mult_pfp", &slc_mult_pfp, &b_ubxsec_event_split_slc_mult_pfp);
   fChain->SetBranchAddress("slc_mult_track", &slc_mult_track, &b_ubxsec_event_split_slc_mult_track);
   fChain->SetBranchAddress("slc_mult_shower", &slc_mult_shower, &b_ubxsec_event_split_slc_mult_shower);
   fChain->SetBranchAddress("slc_mult_track_tolerance", &slc_mult_track_tolerance, &b_ubxsec_event_split_slc_mult_track_tolerance);
   fChain->SetBranchAddress("slc_muoncandidate_exists", &slc_muoncandidate_exists, &b_ubxsec_event_split_slc_muoncandidate_exists);
   fChain->SetBranchAddress("slc_muoncandidate_length", &slc_muoncandidate_length, &b_ubxsec_event_split_slc_muoncandidate_length);
   fChain->SetBranchAddress("slc_muoncandidate_phi", &slc_muoncandidate_phi, &b_ubxsec_event_split_slc_muoncandidate_phi);
   fChain->SetBranchAddress("slc_muoncandidate_theta", &slc_muoncandidate_theta, &b_ubxsec_event_split_slc_muoncandidate_theta);
   fChain->SetBranchAddress("slc_muoncandidate_mom_range", &slc_muoncandidate_mom_range, &b_ubxsec_event_split_slc_muoncandidate_mom_range);
   fChain->SetBranchAddress("slc_muoncandidate_mom_mcs", &slc_muoncandidate_mom_mcs, &b_ubxsec_event_split_slc_muoncandidate_mom_mcs);
   fChain->SetBranchAddress("slc_muoncandidate_mom_mcs_pi", &slc_muoncandidate_mom_mcs_pi, &b_ubxsec_event_split_slc_muoncandidate_mom_mcs_pi);
   fChain->SetBranchAddress("slc_muoncandidate_mcs_ll", &slc_muoncandidate_mcs_ll, &b_ubxsec_event_split_slc_muoncandidate_mcs_ll);
   fChain->SetBranchAddress("slc_muoncandidate_contained", &slc_muoncandidate_contained, &b_ubxsec_event_split_slc_muoncandidate_contained);
   fChain->SetBranchAddress("slc_muoncandidate_dqdx_trunc", &slc_muoncandidate_dqdx_trunc, &b_ubxsec_event_split_slc_muoncandidate_dqdx_trunc);
   fChain->SetBranchAddress("slc_muoncandidate_dqdx_u_trunc", &slc_muoncandidate_dqdx_u_trunc, &b_ubxsec_event_split_slc_muoncandidate_dqdx_u_trunc);
   fChain->SetBranchAddress("slc_muoncandidate_dqdx_v_trunc", &slc_muoncandidate_dqdx_v_trunc, &b_ubxsec_event_split_slc_muoncandidate_dqdx_v_trunc);
   fChain->SetBranchAddress("slc_muoncandidate_dqdx_v", &slc_muoncandidate_dqdx_v, &b_ubxsec_event_split_slc_muoncandidate_dqdx_v);
   fChain->SetBranchAddress("slc_muoncandidate_mip_consistency", &slc_muoncandidate_mip_consistency, &b_ubxsec_event_split_slc_muoncandidate_mip_consistency);
   fChain->SetBranchAddress("slc_muoncandidate_mip_consistency2", &slc_muoncandidate_mip_consistency2, &b_ubxsec_event_split_slc_muoncandidate_mip_consistency2);
   fChain->SetBranchAddress("slc_muoncandidate_truepdg", &slc_muoncandidate_truepdg, &b_ubxsec_event_split_slc_muoncandidate_truepdg);
   fChain->SetBranchAddress("slc_muoncandidate_trueorigin", &slc_muoncandidate_trueorigin, &b_ubxsec_event_split_slc_muoncandidate_trueorigin);
   fChain->SetBranchAddress("slc_muoncandidate_mcs_delta_ll", &slc_muoncandidate_mcs_delta_ll, &b_ubxsec_event_split_slc_muoncandidate_mcs_delta_ll);
   fChain->SetBranchAddress("slc_muoncandidate_residuals_mean", &slc_muoncandidate_residuals_mean, &b_ubxsec_event_split_slc_muoncandidate_residuals_mean);
   fChain->SetBranchAddress("slc_muoncandidate_residuals_std", &slc_muoncandidate_residuals_std, &b_ubxsec_event_split_slc_muoncandidate_residuals_std);
   fChain->SetBranchAddress("slc_muoncandidate_wiregap", &slc_muoncandidate_wiregap, &b_ubxsec_event_split_slc_muoncandidate_wiregap);
   fChain->SetBranchAddress("slc_muoncandidate_wiregap_dead", &slc_muoncandidate_wiregap_dead, &b_ubxsec_event_split_slc_muoncandidate_wiregap_dead);
   fChain->SetBranchAddress("slc_muoncandidate_linearity", &slc_muoncandidate_linearity, &b_ubxsec_event_split_slc_muoncandidate_linearity);
   fChain->SetBranchAddress("slc_muoncandidate_perc_used_hits_in_cluster", &slc_muoncandidate_perc_used_hits_in_cluster, &b_ubxsec_event_split_slc_muoncandidate_perc_used_hits_in_cluster);
   fChain->SetBranchAddress("slc_muoncandidate_maxscatteringangle", &slc_muoncandidate_maxscatteringangle, &b_ubxsec_event_split_slc_muoncandidate_maxscatteringangle);
   fChain->SetBranchAddress("slc_muoncandidate_truth_origin", &slc_muoncandidate_truth_origin, &b_ubxsec_event_split_slc_muoncandidate_truth_origin);
   fChain->SetBranchAddress("slc_muoncandidate_truth_pdg", &slc_muoncandidate_truth_pdg, &b_ubxsec_event_split_slc_muoncandidate_truth_pdg);
   fChain->SetBranchAddress("slc_muoncandidate_truth_time", &slc_muoncandidate_truth_time, &b_ubxsec_event_split_slc_muoncandidate_truth_time);
   fChain->SetBranchAddress("slc_muoncandidate_truth_startx", &slc_muoncandidate_truth_startx, &b_ubxsec_event_split_slc_muoncandidate_truth_startx);
   fChain->SetBranchAddress("slc_muoncandidate_truth_starty", &slc_muoncandidate_truth_starty, &b_ubxsec_event_split_slc_muoncandidate_truth_starty);
   fChain->SetBranchAddress("slc_muoncandidate_truth_startz", &slc_muoncandidate_truth_startz, &b_ubxsec_event_split_slc_muoncandidate_truth_startz);
   fChain->SetBranchAddress("slc_muoncandidate_truth_endx", &slc_muoncandidate_truth_endx, &b_ubxsec_event_split_slc_muoncandidate_truth_endx);
   fChain->SetBranchAddress("slc_muoncandidate_truth_endy", &slc_muoncandidate_truth_endy, &b_ubxsec_event_split_slc_muoncandidate_truth_endy);
   fChain->SetBranchAddress("slc_muoncandidate_truth_endz", &slc_muoncandidate_truth_endz, &b_ubxsec_event_split_slc_muoncandidate_truth_endz);
   fChain->SetBranchAddress("slc_muoncandidate_truth_px", &slc_muoncandidate_truth_px, &b_ubxsec_event_split_slc_muoncandidate_truth_px);
   fChain->SetBranchAddress("slc_muoncandidate_truth_py", &slc_muoncandidate_truth_py, &b_ubxsec_event_split_slc_muoncandidate_truth_py);
   fChain->SetBranchAddress("slc_muoncandidate_truth_pz", &slc_muoncandidate_truth_pz, &b_ubxsec_event_split_slc_muoncandidate_truth_pz);
   fChain->SetBranchAddress("nbeamfls", &nbeamfls, &b_ubxsec_event_split_nbeamfls);
   fChain->SetBranchAddress("beamfls_time", &beamfls_time, &b_ubxsec_event_split_beamfls_time);
   fChain->SetBranchAddress("beamfls_pe", &beamfls_pe, &b_ubxsec_event_split_beamfls_pe);
   fChain->SetBranchAddress("beamfls_z", &beamfls_z, &b_ubxsec_event_split_beamfls_z);
   fChain->SetBranchAddress("candidate_flash_time", &candidate_flash_time, &b_ubxsec_event_split_candidate_flash_time);
   fChain->SetBranchAddress("candidate_flash_z", &candidate_flash_z, &b_ubxsec_event_split_candidate_flash_z);
   fChain->SetBranchAddress("no_mcflash_but_op_activity", &no_mcflash_but_op_activity, &b_ubxsec_event_split_no_mcflash_but_op_activity);
   fChain->SetBranchAddress("beamfls_spec", &beamfls_spec, &b_ubxsec_event_split_beamfls_spec);
   fChain->SetBranchAddress("numc_flash_spec", &numc_flash_spec, &b_ubxsec_event_split_numc_flash_spec);
   fChain->SetBranchAddress("slc_flshypo_xfixed_spec", &slc_flshypo_xfixed_spec, &b_ubxsec_event_split_slc_flshypo_xfixed_spec);
   fChain->SetBranchAddress("slc_flshypo_spec", &slc_flshypo_spec, &b_ubxsec_event_split_slc_flshypo_spec);
   fChain->SetBranchAddress("nsignal", &nsignal, &b_ubxsec_event_split_nsignal);
   fChain->SetBranchAddress("tvtx_x", &tvtx_x, &b_ubxsec_event_split_tvtx_x);
   fChain->SetBranchAddress("tvtx_y", &tvtx_y, &b_ubxsec_event_split_tvtx_y);
   fChain->SetBranchAddress("tvtx_z", &tvtx_z, &b_ubxsec_event_split_tvtx_z);
   fChain->SetBranchAddress("pot", &pot, &b_ubxsec_event_split_pot);
   fChain->SetBranchAddress("evtwgt_genie_pm1_nfunc", &evtwgt_genie_pm1_nfunc, &b_ubxsec_event_split_evtwgt_genie_pm1_nfunc);
   fChain->SetBranchAddress("evtwgt_genie_pm1_funcname", &evtwgt_genie_pm1_funcname, &b_ubxsec_event_split_evtwgt_genie_pm1_funcname);
   fChain->SetBranchAddress("evtwgt_genie_pm1_nweight", &evtwgt_genie_pm1_nweight, &b_ubxsec_event_split_evtwgt_genie_pm1_nweight);
   fChain->SetBranchAddress("evtwgt_genie_pm1_weight", &evtwgt_genie_pm1_weight, &b_ubxsec_event_split_evtwgt_genie_pm1_weight);
   fChain->SetBranchAddress("evtwgt_genie_multisim_nfunc", &evtwgt_genie_multisim_nfunc, &b_ubxsec_event_split_evtwgt_genie_multisim_nfunc);
   fChain->SetBranchAddress("evtwgt_genie_multisim_funcname", &evtwgt_genie_multisim_funcname, &b_ubxsec_event_split_evtwgt_genie_multisim_funcname);
   fChain->SetBranchAddress("evtwgt_genie_multisim_nweight", &evtwgt_genie_multisim_nweight, &b_ubxsec_event_split_evtwgt_genie_multisim_nweight);
   fChain->SetBranchAddress("evtwgt_genie_multisim_weight", &evtwgt_genie_multisim_weight, &b_ubxsec_event_split_evtwgt_genie_multisim_weight);
   fChain->SetBranchAddress("evtwgt_extra_syst_multisim_nfunc", &evtwgt_extra_syst_multisim_nfunc, &b_ubxsec_event_split_evtwgt_extra_syst_multisim_nfunc);
   fChain->SetBranchAddress("evtwgt_extra_syst_multisim_funcname", &evtwgt_extra_syst_multisim_funcname, &b_ubxsec_event_split_evtwgt_extra_syst_multisim_funcname);
   fChain->SetBranchAddress("evtwgt_extra_syst_multisim_nweight", &evtwgt_extra_syst_multisim_nweight, &b_ubxsec_event_split_evtwgt_extra_syst_multisim_nweight);
   fChain->SetBranchAddress("evtwgt_extra_syst_multisim_weight", &evtwgt_extra_syst_multisim_weight, &b_ubxsec_event_split_evtwgt_extra_syst_multisim_weight);
   fChain->SetBranchAddress("evtwgt_flux_multisim_nfunc", &evtwgt_flux_multisim_nfunc, &b_ubxsec_event_split_evtwgt_flux_multisim_nfunc);
   fChain->SetBranchAddress("evtwgt_flux_multisim_funcname", &evtwgt_flux_multisim_funcname, &b_ubxsec_event_split_evtwgt_flux_multisim_funcname);
   fChain->SetBranchAddress("evtwgt_flux_multisim_nweight", &evtwgt_flux_multisim_nweight, &b_ubxsec_event_split_evtwgt_flux_multisim_nweight);
   fChain->SetBranchAddress("evtwgt_flux_multisim_weight", &evtwgt_flux_multisim_weight, &b_ubxsec_event_split_evtwgt_flux_multisim_weight);
   fChain->SetBranchAddress("num_pfp", &num_pfp, &b_ubxsec_event_split_num_pfp);
   fChain->SetBranchAddress("num_pfp_tracks", &num_pfp_tracks, &b_ubxsec_event_split_num_pfp_tracks);
   fChain->SetBranchAddress("num_pfp_showers", &num_pfp_showers, &b_ubxsec_event_split_num_pfp_showers);
   fChain->SetBranchAddress("pfp_truth_pdg", &pfp_truth_pdg, &b_ubxsec_event_split_pfp_truth_pdg);
   fChain->SetBranchAddress("pfp_truth_origin", &pfp_truth_origin, &b_ubxsec_event_split_pfp_truth_origin);
   fChain->SetBranchAddress("pfp_truth_status", &pfp_truth_status, &b_ubxsec_event_split_pfp_truth_status);
   fChain->SetBranchAddress("pfp_truth_parId", &pfp_truth_parId, &b_ubxsec_event_split_pfp_truth_parId);
   fChain->SetBranchAddress("pfp_truth_theta", &pfp_truth_theta, &b_ubxsec_event_split_pfp_truth_theta);
   fChain->SetBranchAddress("pfp_truth_costheta", &pfp_truth_costheta, &b_ubxsec_event_split_pfp_truth_costheta);
   fChain->SetBranchAddress("pfp_truth_phi", &pfp_truth_phi, &b_ubxsec_event_split_pfp_truth_phi);
   fChain->SetBranchAddress("pfp_truth_mom", &pfp_truth_mom, &b_ubxsec_event_split_pfp_truth_mom);
   fChain->SetBranchAddress("pfp_truth_startx", &pfp_truth_startx, &b_ubxsec_event_split_pfp_truth_startx);
   fChain->SetBranchAddress("pfp_truth_starty", &pfp_truth_starty, &b_ubxsec_event_split_pfp_truth_starty);
   fChain->SetBranchAddress("pfp_truth_startz", &pfp_truth_startz, &b_ubxsec_event_split_pfp_truth_startz);
   fChain->SetBranchAddress("pfp_truth_endx", &pfp_truth_endx, &b_ubxsec_event_split_pfp_truth_endx);
   fChain->SetBranchAddress("pfp_truth_endy", &pfp_truth_endy, &b_ubxsec_event_split_pfp_truth_endy);
   fChain->SetBranchAddress("pfp_truth_endz", &pfp_truth_endz, &b_ubxsec_event_split_pfp_truth_endz);
   fChain->SetBranchAddress("pfp_truth_endE", &pfp_truth_endE, &b_ubxsec_event_split_pfp_truth_endE);
   fChain->SetBranchAddress("pfp_truth_endProcess", &pfp_truth_endProcess, &b_ubxsec_event_split_pfp_truth_endProcess);
   fChain->SetBranchAddress("pfp_truth_KE", &pfp_truth_KE, &b_ubxsec_event_split_pfp_truth_KE);
   fChain->SetBranchAddress("pfp_truth_Mass", &pfp_truth_Mass, &b_ubxsec_event_split_pfp_truth_Mass);
   fChain->SetBranchAddress("pfp_reco_isprimary", &pfp_reco_isprimary, &b_ubxsec_event_split_pfp_reco_isprimary);
   fChain->SetBranchAddress("pfp_reco_ndaughters", &pfp_reco_ndaughters, &b_ubxsec_event_split_pfp_reco_ndaughters);
   fChain->SetBranchAddress("pfp_reco_ismuoncandidate", &pfp_reco_ismuoncandidate, &b_ubxsec_event_split_pfp_reco_ismuoncandidate);
   fChain->SetBranchAddress("pfp_reco_istrack", &pfp_reco_istrack, &b_ubxsec_event_split_pfp_reco_istrack);
   fChain->SetBranchAddress("pfp_reco_numtracks", &pfp_reco_numtracks, &b_ubxsec_event_split_pfp_reco_numtracks);
   fChain->SetBranchAddress("pfp_reco_isshower", &pfp_reco_isshower, &b_ubxsec_event_split_pfp_reco_isshower);
   fChain->SetBranchAddress("pfp_reco_numshowers", &pfp_reco_numshowers, &b_ubxsec_event_split_pfp_reco_numshowers);
   fChain->SetBranchAddress("pfp_reco_upflag", &pfp_reco_upflag, &b_ubxsec_event_split_pfp_reco_upflag);
   fChain->SetBranchAddress("pfp_reco_Id", &pfp_reco_Id, &b_ubxsec_event_split_pfp_reco_Id);
   fChain->SetBranchAddress("pfp_reco_length", &pfp_reco_length, &b_ubxsec_event_split_pfp_reco_length);
   fChain->SetBranchAddress("pfp_reco_theta", &pfp_reco_theta, &b_ubxsec_event_split_pfp_reco_theta);
   fChain->SetBranchAddress("pfp_reco_costheta", &pfp_reco_costheta, &b_ubxsec_event_split_pfp_reco_costheta);
   fChain->SetBranchAddress("pfp_reco_phi", &pfp_reco_phi, &b_ubxsec_event_split_pfp_reco_phi);
   fChain->SetBranchAddress("pfp_reco_startx", &pfp_reco_startx, &b_ubxsec_event_split_pfp_reco_startx);
   fChain->SetBranchAddress("pfp_reco_starty", &pfp_reco_starty, &b_ubxsec_event_split_pfp_reco_starty);
   fChain->SetBranchAddress("pfp_reco_startz", &pfp_reco_startz, &b_ubxsec_event_split_pfp_reco_startz);
   fChain->SetBranchAddress("pfp_reco_endx", &pfp_reco_endx, &b_ubxsec_event_split_pfp_reco_endx);
   fChain->SetBranchAddress("pfp_reco_endy", &pfp_reco_endy, &b_ubxsec_event_split_pfp_reco_endy);
   fChain->SetBranchAddress("pfp_reco_endz", &pfp_reco_endz, &b_ubxsec_event_split_pfp_reco_endz);
   fChain->SetBranchAddress("pfp_reco_Mom", &pfp_reco_Mom, &b_ubxsec_event_split_pfp_reco_Mom);
   fChain->SetBranchAddress("pfp_reco_Mom_proton", &pfp_reco_Mom_proton, &b_ubxsec_event_split_pfp_reco_Mom_proton);
   fChain->SetBranchAddress("pfp_reco_Mom_muon", &pfp_reco_Mom_muon, &b_ubxsec_event_split_pfp_reco_Mom_muon);
   fChain->SetBranchAddress("pfp_reco_Mom_MCS", &pfp_reco_Mom_MCS, &b_ubxsec_event_split_pfp_reco_Mom_MCS);
   fChain->SetBranchAddress("pfp_reco_trunmeandqdx", &pfp_reco_trunmeandqdx, &b_ubxsec_event_split_pfp_reco_trunmeandqdx);
   fChain->SetBranchAddress("pfp_reco_trunmeandqdx_U", &pfp_reco_trunmeandqdx_U, &b_ubxsec_event_split_pfp_reco_trunmeandqdx_U);
   fChain->SetBranchAddress("pfp_reco_trunmeandqdx_V", &pfp_reco_trunmeandqdx_V, &b_ubxsec_event_split_pfp_reco_trunmeandqdx_V);
   fChain->SetBranchAddress("pfp_reco_newpid_pida", &pfp_reco_newpid_pida, &b_ubxsec_event_split_pfp_reco_newpid_pida);
   fChain->SetBranchAddress("pfp_reco_chi2_proton", &pfp_reco_chi2_proton, &b_ubxsec_event_split_pfp_reco_chi2_proton);
   fChain->SetBranchAddress("pfp_reco_chi2_kaon", &pfp_reco_chi2_kaon, &b_ubxsec_event_split_pfp_reco_chi2_kaon);
   fChain->SetBranchAddress("pfp_reco_chi2_muon", &pfp_reco_chi2_muon, &b_ubxsec_event_split_pfp_reco_chi2_muon);
   fChain->SetBranchAddress("pfp_reco_chi2_pion", &pfp_reco_chi2_pion, &b_ubxsec_event_split_pfp_reco_chi2_pion);
   fChain->SetBranchAddress("pfp_reco_bragg_ratio", &pfp_reco_bragg_ratio, &b_ubxsec_event_split_pfp_reco_bragg_ratio);
   fChain->SetBranchAddress("pfp_reco_bragg_proton", &pfp_reco_bragg_proton, &b_ubxsec_event_split_pfp_reco_bragg_proton);
   fChain->SetBranchAddress("pfp_reco_bragg_fwd_proton", &pfp_reco_bragg_fwd_proton, &b_ubxsec_event_split_pfp_reco_bragg_fwd_proton);
   fChain->SetBranchAddress("pfp_reco_bragg_bwd_proton", &pfp_reco_bragg_bwd_proton, &b_ubxsec_event_split_pfp_reco_bragg_bwd_proton);
   fChain->SetBranchAddress("pfp_reco_bragg_fwd_mip", &pfp_reco_bragg_fwd_mip, &b_ubxsec_event_split_pfp_reco_bragg_fwd_mip);
   fChain->SetBranchAddress("pfp_reco_nhits", &pfp_reco_nhits, &b_ubxsec_event_split_pfp_reco_nhits);
   fChain->SetBranchAddress("pfp_reco_dEdx", &pfp_reco_dEdx, &b_ubxsec_event_split_pfp_reco_dEdx);
   fChain->SetBranchAddress("pfp_reco_dQdx", &pfp_reco_dQdx, &b_ubxsec_event_split_pfp_reco_dQdx);
   fChain->SetBranchAddress("pfp_reco_RR", &pfp_reco_RR, &b_ubxsec_event_split_pfp_reco_RR);
   fChain->SetBranchAddress("_default_value", &_default_value, &b_ubxsec_event_split__default_value);
   Notify();
}

Bool_t UBXSecEvent::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void UBXSecEvent::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t UBXSecEvent::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef UBXSecEvent_cxx
