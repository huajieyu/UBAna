#ifndef __DATATYPES_UBXSECEVENTHISTO1D_CXX__
#define __DATATYPES_UBXSECEVENTHISTO1D_CXX__

#include "UBXSecEventHisto1D.h"


namespace DataTypes {


    void UBXSecEventHisto1D::InitializeBootstraps()
    {

      //
      // Total Cross Section
      //

      // Efficiency - Total Cross Section
      h_eff_onebin_num = new TH1D("h_eff_onebin_num", "h_eff_onebin_num", 1, 0, 1);
      h_eff_onebin_den = new TH1D("h_eff_onebin_den", "h_eff_onebin_den", 1, 0, 1);

      // Number of events per channel - Total Cross Section
      hmap_onebin["total"] = new TH1D("h_onebin_total", "; Track length;", 1, 0, 1);
      hmap_onebin["signal"] = new TH1D("h_onebin_signal", "; Track length;", 1, 0, 1);
      hmap_onebin["cosmic"] = new TH1D("h_onebin_cosmic", "; Track length;", 1, 0, 1);
      hmap_onebin["cosmic_stopmu"] = new TH1D("h_onebin_cosmic_stopmu", "; Track length;", 1, 0, 1);
      hmap_onebin["cosmic_nostopmu"] = new TH1D("h_onebin_cosmic_nostopmu", "; Track length;", 1, 0, 1);
      hmap_onebin["outfv"] = new TH1D("h_onebin_outfv", "; Track length;", 1, 0, 1);
      hmap_onebin["outfv_stopmu"] = new TH1D("h_onebin_outfv_stopmu", "; Track length;", 1, 0, 1);
      hmap_onebin["outfv_nostopmu"] = new TH1D("h_onebin_outfv_nostopmu", "; Track length;", 1, 0, 1);
      hmap_onebin["nc"] = new TH1D("h_onebin_nc", "; Track length;", 1, 0, 1);
      hmap_onebin["nc_proton"] = new TH1D("h_onebin_nc_proton", "; Track length;", 1, 0, 1);
      hmap_onebin["nc_pion"] = new TH1D("h_onebin_nc_pion", "; Track length;", 1, 0, 1);
      hmap_onebin["nc_other"] = new TH1D("h_onebin_nc_other", "; Track length;", 1, 0, 1);
      hmap_onebin["anumu"] = new TH1D("h_onebin_anumu", "; Track length;", 1, 0, 1);
      hmap_onebin["nue"] = new TH1D("h_onebin_nue", "; Track length;", 1, 0, 1);
      hmap_onebin["signal_stopmu"] = new TH1D("h_onebin_signal_stopmu", "; Track length;", 1, 0, 1);
      hmap_onebin["signal_nostopmu"] = new TH1D("h_onebin_signal_nostopmu", "; Track length;", 1, 0, 1);
      hmap_onebin["cc_other"] = new TH1D("h_onebin_ccother", "; Track length;", 1, 0, 1);
      // Efficiency - Total Cross Section - Multisim
      bs_genie_multisim_eff_onebin_num = new BootstrapTH1D("bs_genie_multisim_eff_onebin_num", "bs_genie_multisim_eff_onebin_num_title", 1, 0, 1);
      bs_genie_multisim_eff_onebin_den = new BootstrapTH1D("bs_genie_multisim_eff_onebin_den", "bs_genie_multisim_eff_onebin_den_title", 1, 0, 1);
      bs_flux_multisim_eff_onebin_num = new BootstrapTH1D("bs_flux_multisim_eff_onebin_num", "bs_flux_multisim_eff_onebin_num_title", 1, 0, 1);
      bs_flux_multisim_eff_onebin_den = new BootstrapTH1D("bs_flux_multisim_eff_onebin_den", "bs_flux_multisim_eff_onebin_den_title", 1, 0, 1);
      bs_extra_syst_multisim_eff_onebin_num = new BootstrapTH1D("bs_extra_syst_multisim_eff_onebin_num", "bs_extra_syst_eff_onebin_num_title", 1, 0, 1);
      bs_extra_syst_multisim_eff_onebin_den = new BootstrapTH1D("bs_extra_syst_multisim_eff_onebin_den", "bs_extra_syst_eff_onebin_den_title", 1, 0, 1);
      bs_mc_stat_multisim_eff_onebin_num = new BootstrapTH1D("bs_mc_stat_multisim_eff_onebin_num", "bs_mc_stat_eff_onebin_num_title", 1, 0, 1);
      bs_mc_stat_multisim_eff_onebin_den = new BootstrapTH1D("bs_mc_stat_multisim_eff_onebin_den", "bs_mc_stat_eff_onebin_den_title", 1, 0, 1);

      // Number of events per channel and universe - Total Cross Section - Genie Multisim
      hmap_onebin_genie_multisim_bs["total"]["nominal"] = new TH1D("h_onebin_total_genie_mulinominal", "; Track length;", 1, 0, 1);
      hmap_onebin_genie_multisim_bs["signal"]["nominal"] = new TH1D("h_onebin_signal_genie_mulinominal", "; Track length;", 1, 0, 1);
      hmap_onebin_genie_multisim_bs["cosmic"]["nominal"] = new TH1D("h_onebin_cosmic_genie_mulinominal", "; Track length;", 1, 0, 1);
      hmap_onebin_genie_multisim_bs["outfv"]["nominal"] = new TH1D("h_onebin_outfv_genie_mulinominal", "; Track length;", 1, 0, 1);
      hmap_onebin_genie_multisim_bs["nc"]["nominal"] = new TH1D("h_onebin_nc_genie_mulinominal", "; Track length;", 1, 0, 1);
      hmap_onebin_genie_multisim_bs["anumu"]["nominal"] = new TH1D("h_onebin_anumu_genie_mulinominal", "; Track length;", 1, 0, 1);
      hmap_onebin_genie_multisim_bs["nue"]["nominal"] = new TH1D("h_onebin_nue_genie_mulinominal", "; Track length;", 1, 0, 1);
      hmap_onebin_genie_multisim_bs["cc_other"]["nominal"] = new TH1D("h_onebin_ccother_genie_multinominal", "; Track length;", 1, 0, 1);
      // Number of events per channel and universe - Total Cross Section - Flux Multisim
      hmap_onebin_flux_multisim_bs["total"]["nominal"] = new TH1D("h_onebin_total_flux_mulinominal", "; Track length;", 1, 0, 1);
      hmap_onebin_flux_multisim_bs["signal"]["nominal"] = new TH1D("h_onebin_signal_flux_mulinominal", "; Track length;", 1, 0, 1);
      hmap_onebin_flux_multisim_bs["cosmic"]["nominal"] = new TH1D("h_onebin_cosmic_flux_mulinominal", "; Track length;", 1, 0, 1);
      hmap_onebin_flux_multisim_bs["outfv"]["nominal"] = new TH1D("h_onebin_outfv_flux_mulinominal", "; Track length;", 1, 0, 1);
      hmap_onebin_flux_multisim_bs["nc"]["nominal"] = new TH1D("h_onebin_nc_flux_mulinominal", "; Track length;", 1, 0, 1);
      hmap_onebin_flux_multisim_bs["anumu"]["nominal"] = new TH1D("h_onebin_anumu_flux_mulinominal", "; Track length;", 1, 0, 1);
      hmap_onebin_flux_multisim_bs["nue"]["nominal"] = new TH1D("h_onebin_nue_flux_mulinominal", "; Track length;", 1, 0, 1);
      hmap_onebin_flux_multisim_bs["cc_other"]["nominal"] = new TH1D("h_onebin_ccother_flux_multinominal", "; Track length;", 1, 0, 1);
      // Number of events per channel and universe - Total Cross Section - Extra Syst
      hmap_onebin_extra_syst_multisim_bs["total"]["nominal"] = new TH1D("h_onebin_total_extra_syst_multisim_nominal", "; Track length;", 1, 0, 1);
      hmap_onebin_extra_syst_multisim_bs["signal"]["nominal"] = new TH1D("h_onebin_signal_extra_syst_multisim_nominal", "; Track length;", 1, 0, 1);
      hmap_onebin_extra_syst_multisim_bs["cosmic"]["nominal"] = new TH1D("h_onebin_cosmic_extra_syst_multisim_nominal", "; Track length;", 1, 0, 1);
      hmap_onebin_extra_syst_multisim_bs["outfv"]["nominal"] = new TH1D("h_onebin_outfv_extra_syst_multisim_nominal", "; Track length;", 1, 0, 1);
      hmap_onebin_extra_syst_multisim_bs["nc"]["nominal"] = new TH1D("h_onebin_nc_extra_syst_multisim_nominal", "; Track length;", 1, 0, 1);
      hmap_onebin_extra_syst_multisim_bs["anumu"]["nominal"] = new TH1D("h_onebin_anumu_extra_syst_multisim_nominal", "; Track length;", 1, 0, 1);
      hmap_onebin_extra_syst_multisim_bs["nue"]["nominal"] = new TH1D("h_onebin_nue_extra_syst_multisim_nominal", "; Track length;", 1, 0, 1);
      hmap_onebin_extra_syst_multisim_bs["cc_other"]["nominal"] = new TH1D("h_onebin_ccother_extra_syst_multisim_nominal","; Track length;", 1, 0, 1);
      // Number of events per channel and universe - Total Cross Section - MC Stat
      hmap_onebin_mc_stat_multisim_bs["total"]["nominal"] = new TH1D("h_onebin_total_mc_stat_multisim_nominal", "; Track length;", 1, 0, 1);
      hmap_onebin_mc_stat_multisim_bs["signal"]["nominal"] = new TH1D("h_onebin_signal_mc_stat_multisim_nominal", "; Track length;", 1, 0, 1);
      hmap_onebin_mc_stat_multisim_bs["cosmic"]["nominal"] = new TH1D("h_onebin_cosmic_mc_stat_multisim_nominal", "; Track length;", 1, 0, 1);
      hmap_onebin_mc_stat_multisim_bs["outfv"]["nominal"] = new TH1D("h_onebin_outfv_mc_stat_multisim_nominal", "; Track length;", 1, 0, 1);
      hmap_onebin_mc_stat_multisim_bs["nc"]["nominal"] = new TH1D("h_onebin_nc_mc_stat_multisim_nominal", "; Track length;", 1, 0, 1);
      hmap_onebin_mc_stat_multisim_bs["anumu"]["nominal"] = new TH1D("h_onebin_anumu_mc_stat_multisim_nominal", "; Track length;", 1, 0, 1);
      hmap_onebin_mc_stat_multisim_bs["nue"]["nominal"] = new TH1D("h_onebin_nue_mc_stat_multisim_nominal", "; Track length;", 1, 0, 1);
      hmap_onebin_mc_stat_multisim_bs["cc_other"]["nominal"]= new TH1D("h_onebin_nue_mc_stat_multisim_nominal", "; Track length", 1, 0, 1);
      

      //==================================================================================================================================
      //
      // Single differential (mumom, pmom, thetamup)
      //

      // Efficiency - Single Differential (mumom)
      h_eff_mumom_num = new TH1D("h_eff_mumom_num", "h_eff_mumom_num", n_bins_mumom, bins_mumom);
      h_eff_mumom_den = new TH1D("h_eff_mumom_den", "h_eff_mumom_den", n_bins_mumom, bins_mumom);

      h_eff_pmom_num = new TH1D("h_eff_pmom_num", "h_eff_pmmom_num", n_bins_pmom, bins_pmom);
      h_eff_pmom_den = new TH1D("h_eff_pmom_den", "h_eff_pmmom_den", n_bins_pmom, bins_pmom);

      h_eff_thetamup_num = new TH1D("h_eff_thetamup_num", "h_eff_thetamup_num", n_bins_muptheta, bins_muptheta);
      h_eff_thetamup_den = new TH1D("h_eff_thetamup_den", "h_eff_thetamup_den", n_bins_muptheta, bins_muptheta);

      // Reco to true histograms - Single Differential (mumom)
      h_true_reco_mom= new TH2D("h_true_reco_mom", ";Muon Momentum (Truth) [GeV]; Muon Momentum (MCS) [GeV]", n_bins_mumom, bins_mumom, n_bins_mumom, bins_mumom);
      h_true_reco_pmom= new TH2D("h_true_reco_pmom", ";Proton Momentum (Truth) [GeV]; Proton Momentum (Reco) [GeV]", n_bins_pmom, bins_pmom, n_bins_pmom, bins_pmom);
      h_true_reco_thetamup = new TH2D("h_true_reco_thetamup", "#theta_{#mu,p}[Rad]", n_bins_muptheta, bins_muptheta, n_bins_muptheta, bins_muptheta);

      // Number of events per channel - Single Differential (mumom)
      hmap_trkmom["total"] = new TH1D("h_trkmom_total", "; Track momentum;", n_bins_mumom, bins_mumom);
      hmap_trkmom["signal"] = new TH1D("h_trkmom_signal", "; Track momentum;", n_bins_mumom, bins_mumom);
      hmap_trkmom["cosmic"] = new TH1D("h_trkmom_cosmic", "; Track momentum;", n_bins_mumom, bins_mumom);
      hmap_trkmom["cosmic_stopmu"] = new TH1D("h_trkmom_cosmic_stopmu", "; Track momentum;", n_bins_mumom, bins_mumom);
      hmap_trkmom["cosmic_nostopmu"] = new TH1D("h_trkmom_cosmic_nostopmu", "; Track momentum;", n_bins_mumom, bins_mumom);
      hmap_trkmom["outfv"] = new TH1D("h_trkmom_outfv", "; Track momentum;", n_bins_mumom, bins_mumom);
      hmap_trkmom["outfv_stopmu"] = new TH1D("h_trkmom_outfv_stopmu", "; Track momentum;", n_bins_mumom, bins_mumom);
      hmap_trkmom["outfv_nostopmu"] = new TH1D("h_trkmom_outfv_nostopmu", "; Track momentum;", n_bins_mumom, bins_mumom);
      hmap_trkmom["nc"] = new TH1D("h_trkmom_nc", "; Track momentum;", n_bins_mumom, bins_mumom);
      hmap_trkmom["nc_proton"] = new TH1D("h_trkmom_nc_proton", "; Track momentum;", n_bins_mumom, bins_mumom);
      hmap_trkmom["nc_pion"] = new TH1D("h_trkmom_nc_pion", "; Track momentum;", n_bins_mumom, bins_mumom);
      hmap_trkmom["nc_other"] = new TH1D("h_trkmom_nc_other", "; Track momentum;", n_bins_mumom, bins_mumom);
      hmap_trkmom["anumu"] = new TH1D("h_trkmom_anumu", "; Track momentum;", n_bins_mumom, bins_mumom);
      hmap_trkmom["nue"] = new TH1D("h_trkmom_nue", "; Track momentum;", n_bins_mumom, bins_mumom);
      hmap_trkmom["signal_stopmu"] = new TH1D("h_trkmom_signal_stopmu", "; Track momentum;", n_bins_mumom, bins_mumom);
      hmap_trkmom["signal_nostopmu"] = new TH1D("h_trkmom_signal_nostopmu", "; Track momentum;", n_bins_mumom, bins_mumom);
      hmap_trkmom["cc_other"] = new TH1D("h_trkmom_ccother", "; Track momentum;", n_bins_mumom, bins_mumom); 
      //==============================================================================================================
      hmap_trkpmom["total"] = new TH1D("h_trkpmom_total", "; Track-pcand momentum;", n_bins_pmom, bins_pmom);
      hmap_trkpmom["signal"] = new TH1D("h_trkpmom_signal", "; Track-pcand momentum;", n_bins_pmom, bins_pmom);
      hmap_trkpmom["cosmic"] = new TH1D("h_trkpmom_cosmic", "; Track-pcand momentum;", n_bins_pmom, bins_pmom);
      hmap_trkpmom["cosmic_stopmu"] = new TH1D("h_trkpmom_cosmic_stopmu", "; Track-pcand momentum;", n_bins_pmom, bins_pmom);
      hmap_trkpmom["cosmic_nostopmu"] = new TH1D("h_trkpmom_cosmic_nostopmu", "; Track-pcand momentum;", n_bins_pmom, bins_pmom);
      hmap_trkpmom["outfv"] = new TH1D("h_trkpmom_outfv", "; Track-pcand momentum;", n_bins_pmom, bins_pmom);
      hmap_trkpmom["outfv_stopmu"] = new TH1D("h_trkpmom_outfv_stopmu", "; Track-pcand momentum;", n_bins_pmom, bins_pmom);
      hmap_trkpmom["outfv_nostopmu"] = new TH1D("h_trkpmom_outfv_nostopmu", "; Track-pcand momentum;", n_bins_pmom, bins_pmom);
      hmap_trkpmom["nc"] = new TH1D("h_trkpmom_nc", "; Track-pcand momentum;", n_bins_pmom, bins_pmom);
      hmap_trkpmom["nc_proton"] = new TH1D("h_trkpmom_nc_proton", "; Track-pcand momentum;", n_bins_pmom, bins_pmom);
      hmap_trkpmom["nc_pion"] = new TH1D("h_trkpmom_nc_pion", "; Track-pcand momentum;", n_bins_pmom, bins_pmom);
      hmap_trkpmom["nc_other"] = new TH1D("h_trkpmom_nc_other", "; Track-pcand momentum;", n_bins_pmom, bins_pmom);
      hmap_trkpmom["anumu"] = new TH1D("h_trkpmom_anumu", "; Track-pcand momentum;", n_bins_pmom, bins_pmom);
      hmap_trkpmom["nue"] = new TH1D("h_trkpmom_nue", "; Track-pcand momentum;", n_bins_pmom, bins_pmom);
      hmap_trkpmom["signal_stopmu"] = new TH1D("h_trkpmom_signal_stopmu", "; Track-pcand momentum;", n_bins_pmom, bins_pmom);
      hmap_trkpmom["signal_nostopmu"] = new TH1D("h_trkpmom_signal_nostopmu", "; Track-pcand momentum;", n_bins_pmom, bins_pmom);
      hmap_trkpmom["cc_other"] = new TH1D("h_trkpmom_ccother", "; Track-pcand momentum;", n_bins_pmom, bins_pmom);
      //==============================================================================================================
      hmap_thetamup["total"] = new TH1D("h_thetamup_total", "; Track-pcand momentum;", n_bins_muptheta, bins_muptheta);
      hmap_thetamup["signal"] = new TH1D("h_thetamup_signal", "; Track-pcand momentum;", n_bins_muptheta, bins_muptheta);
      hmap_thetamup["cosmic"] = new TH1D("h_thetamup_cosmic", "; Track-pcand momentum;", n_bins_muptheta, bins_muptheta);
      hmap_thetamup["cosmic_stopmu"] = new TH1D("h_thetamup_cosmic_stopmu", "; Track-pcand momentum;", n_bins_muptheta, bins_muptheta);
      hmap_thetamup["cosmic_nostopmu"] = new TH1D("h_thetamup_cosmic_nostopmu", "; Track-pcand momentum;", n_bins_muptheta, bins_muptheta);
      hmap_thetamup["outfv"] = new TH1D("h_thetamup_outfv", "; Track-pcand momentum;", n_bins_muptheta, bins_muptheta);
      hmap_thetamup["outfv_stopmu"] = new TH1D("h_thetamup_outfv_stopmu", "; Track-pcand momentum;", n_bins_muptheta, bins_muptheta);
      hmap_thetamup["outfv_nostopmu"] = new TH1D("h_thetamup_outfv_nostopmu", "; Track-pcand momentum;", n_bins_muptheta, bins_muptheta);
      hmap_thetamup["nc"] = new TH1D("h_thetamup_nc", "; Track-pcand momentum;", n_bins_muptheta, bins_muptheta);
      hmap_thetamup["nc_proton"] = new TH1D("h_thetamup_nc_proton", "; Track-pcand momentum;", n_bins_muptheta, bins_muptheta);
      hmap_thetamup["nc_pion"] = new TH1D("h_thetamup_nc_pion", "; Track-pcand momentum;", n_bins_muptheta, bins_muptheta);
      hmap_thetamup["nc_other"] = new TH1D("h_thetamup_nc_other", "; Track-pcand momentum;", n_bins_muptheta, bins_muptheta);
      hmap_thetamup["anumu"] = new TH1D("h_thetamup_anumu", "; Track-pcand momentum;", n_bins_muptheta, bins_muptheta);
      hmap_thetamup["nue"] = new TH1D("h_thetamup_nue", "; Track-pcand momentum;", n_bins_muptheta, bins_muptheta);
      hmap_thetamup["signal_stopmu"] = new TH1D("h_thetamup_signal_stopmu", "; Track-pcand momentum;", n_bins_muptheta, bins_muptheta);
      hmap_thetamup["signal_nostopmu"] = new TH1D("h_thetamup_signal_nostopmu", "; Track-pcand momentum;", n_bins_muptheta, bins_muptheta);
      hmap_thetamup["cc_other"] = new TH1D("h_thetamup_ccother", "; Track-pcand momentum;", n_bins_muptheta, bins_muptheta);


     // Efficiency - Single Differential (mumom) - Multisim
      bs_genie_multisim_eff_mumom_num = new BootstrapTH1D("bs_genie_multisim_eff_mumom_num", "bs_genie_multisim_eff_mumom_num_title", n_bins_mumom, bins_mumom);
      bs_genie_multisim_eff_mumom_den = new BootstrapTH1D("bs_genie_multisim_eff_mumom_den", "bs_genie_multisim_eff_mumom_den_title", n_bins_mumom, bins_mumom);
      bs_flux_multisim_eff_mumom_num = new BootstrapTH1D("bs_flux_multisim_eff_mumom_num", "bs_flux_multisim_eff_mumom_num_title", n_bins_mumom, bins_mumom);
      bs_flux_multisim_eff_mumom_den = new BootstrapTH1D("bs_flux_multisim_eff_mumom_den", "bs_flux_multisim_eff_mumom_den_title", n_bins_mumom, bins_mumom);

      bs_extra_syst_multisim_eff_mumom_num = new BootstrapTH1D("bs_extra_syst_multisim_eff_mumom_num", "bs_extra_syst_eff_mumom_num_title", n_bins_mumom, bins_mumom);
      bs_extra_syst_multisim_eff_mumom_den = new BootstrapTH1D("bs_extra_syst_multisim_eff_mumom_den", "bs_extra_syst_eff_mumom_den_title", n_bins_mumom, bins_mumom);
      bs_mc_stat_multisim_eff_mumom_num = new BootstrapTH1D("bs_mc_stat_multisim_eff_mumom_num", "bs_mc_stat_eff_mumom_num_title", n_bins_mumom, bins_mumom);
      bs_mc_stat_multisim_eff_mumom_den = new BootstrapTH1D("bs_mc_stat_multisim_eff_mumom_den", "bs_mc_stat_eff_mumom_den_title", n_bins_mumom, bins_mumom);

      bs_genie_multisim_eff_pmom_num = new BootstrapTH1D("bs_genie_multisim_eff_pmom_num", "bs_genie_multisim_eff_pmom_num_title", n_bins_pmom, bins_pmom);
      bs_genie_multisim_eff_pmom_den = new BootstrapTH1D("bs_genie_multisim_eff_pmom_den", "bs_genie_multisim_eff_pmom_den_title", n_bins_pmom, bins_pmom);
      bs_flux_multisim_eff_pmom_num = new BootstrapTH1D("bs_flux_multisim_eff_pmom_num", "bs_flux_multisim_eff_pmom_num_title", n_bins_pmom, bins_pmom);
      bs_flux_multisim_eff_pmom_den = new BootstrapTH1D("bs_flux_multisim_eff_pmom_den", "bs_flux_multisim_eff_pmom_den_title", n_bins_pmom, bins_pmom);

      bs_extra_syst_multisim_eff_pmom_num = new BootstrapTH1D("bs_extra_syst_multisim_eff_pmom_num", "bs_extra_syst_eff_pmom_num_title", n_bins_pmom, bins_pmom);
      bs_extra_syst_multisim_eff_pmom_den = new BootstrapTH1D("bs_extra_syst_multisim_eff_pmom_den", "bs_extra_syst_eff_pmom_den_title", n_bins_pmom, bins_pmom);
      bs_mc_stat_multisim_eff_pmom_num = new BootstrapTH1D("bs_mc_stat_multisim_eff_pmom_num", "bs_mc_stat_eff_pmom_num_title", n_bins_pmom, bins_pmom);
      bs_mc_stat_multisim_eff_pmom_den = new BootstrapTH1D("bs_mc_stat_multisim_eff_pmom_den", "bs_mc_stat_eff_pmom_den_title", n_bins_pmom, bins_pmom);

  
      bs_genie_multisim_eff_thetamup_num = new BootstrapTH1D("bs_genie_multisim_eff_thetamup_num", "bs_genie_multisim_eff_thetamup_num_title", n_bins_muptheta, bins_muptheta);
      bs_genie_multisim_eff_thetamup_den = new BootstrapTH1D("bs_genie_multisim_eff_thetamup_den", "bs_genie_multisim_eff_thetamup_den_title", n_bins_muptheta, bins_muptheta);
      bs_flux_multisim_eff_thetamup_num = new BootstrapTH1D("bs_flux_multisim_eff_thetamup_num", "bs_flux_multisim_eff_thetamup_num_title", n_bins_muptheta, bins_muptheta);
      bs_flux_multisim_eff_thetamup_den = new BootstrapTH1D("bs_flux_multisim_eff_thetamup_den", "bs_flux_multisim_eff_thetamup_den_title", n_bins_muptheta, bins_muptheta);

      bs_extra_syst_multisim_eff_thetamup_num = new BootstrapTH1D("bs_extra_syst_multisim_eff_thetamup_num", "bs_extra_syst_eff_thetamup_num_title", n_bins_muptheta, bins_muptheta);
      bs_extra_syst_multisim_eff_thetamup_den = new BootstrapTH1D("bs_extra_syst_multisim_eff_thetamup_den", "bs_extra_syst_eff_thetamup_den_title", n_bins_muptheta, bins_muptheta);
      bs_mc_stat_multisim_eff_thetamup_num = new BootstrapTH1D("bs_mc_stat_multisim_eff_thetamup_num", "bs_mc_stat_eff_thetamup_num_title", n_bins_muptheta, bins_muptheta);
      bs_mc_stat_multisim_eff_thetamup_den = new BootstrapTH1D("bs_mc_stat_multisim_eff_thetamup_den", "bs_mc_stat_eff_thetamup_den_title", n_bins_muptheta, bins_muptheta);
  
      // Number of events per channel and universe - Single Differential (mumom) - Genie Multisim
      hmap_trkmom_genie_multisim_bs["total"]["nominal"] = new TH1D("h_trkmom_total_genie_mulinominal", "; Track length;", n_bins_mumom, bins_mumom); // 20, 0, 2.5
      hmap_trkmom_genie_multisim_bs["signal"]["nominal"] = new TH1D("h_trkmom_signal_genie_mulinominal", "; Track length;", n_bins_mumom, bins_mumom);
      hmap_trkmom_genie_multisim_bs["cosmic"]["nominal"] = new TH1D("h_trkmom_cosmic_genie_mulinominal", "; Track length;", n_bins_mumom, bins_mumom);
      hmap_trkmom_genie_multisim_bs["outfv"]["nominal"] = new TH1D("h_trkmom_outfv_genie_mulinominal", "; Track length;", n_bins_mumom, bins_mumom);
      hmap_trkmom_genie_multisim_bs["nc"]["nominal"] = new TH1D("h_trkmom_nc_genie_mulinominal", "; Track length;", n_bins_mumom, bins_mumom);
      hmap_trkmom_genie_multisim_bs["anumu"]["nominal"] = new TH1D("h_trkmom_anumu_genie_mulinominal", "; Track length;", n_bins_mumom, bins_mumom);
      hmap_trkmom_genie_multisim_bs["nue"]["nominal"] = new TH1D("h_trkmom_nue_genie_mulinominal", "; Track length;", n_bins_mumom, bins_mumom);
      hmap_trkmom_genie_multisim_bs["cc_other"]["nominal"] = new TH1D("h_trkmom_ccother_genie_mulinominal", "; Track length;", n_bins_mumom, bins_mumom);

      // Number of events per channel and universe - Single Differential (pmom) - Genie Multisim
      hmap_trkpmom_genie_multisim_bs["total"]["nominal"] = new TH1D("h_trkpmom_total_genie_mulinominal", "; Track length;", n_bins_pmom, bins_pmom); // 20, 0, 2.5
      hmap_trkpmom_genie_multisim_bs["signal"]["nominal"] = new TH1D("h_trkpmom_signal_genie_mulinominal", "; Track length;", n_bins_pmom, bins_pmom);
      hmap_trkpmom_genie_multisim_bs["cosmic"]["nominal"] = new TH1D("h_trkpmom_cosmic_genie_mulinominal", "; Track length;", n_bins_pmom, bins_pmom);
      hmap_trkpmom_genie_multisim_bs["outfv"]["nominal"] = new TH1D("h_trkpmom_outfv_genie_mulinominal", "; Track length;", n_bins_pmom, bins_pmom);
      hmap_trkpmom_genie_multisim_bs["nc"]["nominal"] = new TH1D("h_trkpmom_nc_genie_mulinominal", "; Track length;", n_bins_pmom, bins_pmom);
      hmap_trkpmom_genie_multisim_bs["anumu"]["nominal"] = new TH1D("h_trkpmom_anumu_genie_mulinominal", "; Track length;", n_bins_pmom, bins_pmom);
      hmap_trkpmom_genie_multisim_bs["nue"]["nominal"] = new TH1D("h_trkpmom_nue_genie_mulinominal", "; Track length;", n_bins_pmom, bins_pmom);
      hmap_trkpmom_genie_multisim_bs["cc_other"]["nominal"] = new TH1D("h_trkpmom_ccother_genie_mulinominal", "; Track length;", n_bins_pmom, bins_pmom);

      // Number of events per channel and universe - Single Differential (pmom) - Genie Multisim
      hmap_thetamup_genie_multisim_bs["total"]["nominal"] = new TH1D("h_thetamup_total_genie_mulinominal", "; MuP angle;", n_bins_muptheta, bins_muptheta); // 20, 0, 2.5
      hmap_thetamup_genie_multisim_bs["signal"]["nominal"] = new TH1D("h_thetamup_signal_genie_mulinominal", "; MuP angle;", n_bins_muptheta, bins_muptheta);
      hmap_thetamup_genie_multisim_bs["cosmic"]["nominal"] = new TH1D("h_thetamup_cosmic_genie_mulinominal", "; MuP angle;", n_bins_muptheta, bins_muptheta);
      hmap_thetamup_genie_multisim_bs["outfv"]["nominal"] = new TH1D("h_thetamup_outfv_genie_mulinominal", "; MuP angle;", n_bins_muptheta, bins_muptheta);
      hmap_thetamup_genie_multisim_bs["nc"]["nominal"] = new TH1D("h_thetamup_nc_genie_mulinominal", "; MuP angle;", n_bins_muptheta, bins_muptheta);
      hmap_thetamup_genie_multisim_bs["anumu"]["nominal"] = new TH1D("h_thetamup_anumu_genie_mulinominal", "; MuP angle;", n_bins_muptheta, bins_muptheta);
      hmap_thetamup_genie_multisim_bs["nue"]["nominal"] = new TH1D("h_thetamup_nue_genie_mulinominal", "; MuP angle;", n_bins_muptheta, bins_muptheta);
      hmap_thetamup_genie_multisim_bs["cc_other"]["nominal"] = new TH1D("h_thetamup_ccother_genie_mulinominal", "; MuP angle;", n_bins_muptheta, bins_muptheta);



      // Number of events per channel and universe - Single Differential (mumom) - Flux Multisim
      hmap_trkmom_flux_multisim_bs["total"]["nominal"] = new TH1D("h_trkmom_total_flux_mulinominal", "; Track length;", n_bins_mumom, bins_mumom); // 20, 0, 2.5
      hmap_trkmom_flux_multisim_bs["signal"]["nominal"] = new TH1D("h_trkmom_signal_flux_mulinominal", "; Track length;", n_bins_mumom, bins_mumom);
      hmap_trkmom_flux_multisim_bs["cosmic"]["nominal"] = new TH1D("h_trkmom_cosmic_flux_mulinominal", "; Track length;", n_bins_mumom, bins_mumom);
      hmap_trkmom_flux_multisim_bs["outfv"]["nominal"] = new TH1D("h_trkmom_outfv_flux_mulinominal", "; Track length;", n_bins_mumom, bins_mumom);
      hmap_trkmom_flux_multisim_bs["nc"]["nominal"] = new TH1D("h_trkmom_nc_flux_mulinominal", "; Track length;", n_bins_mumom, bins_mumom);
      hmap_trkmom_flux_multisim_bs["anumu"]["nominal"] = new TH1D("h_trkmom_anumu_flux_mulinominal", "; Track length;", n_bins_mumom, bins_mumom);
      hmap_trkmom_flux_multisim_bs["nue"]["nominal"] = new TH1D("h_trkmom_nue_flux_mulinominal", "; Track length;", n_bins_mumom, bins_mumom);
      hmap_trkmom_flux_multisim_bs["cc_other"]["nominal"] = new TH1D("h_trkmom_ccother_flux_mulinominal", "; Track length;", n_bins_mumom, bins_mumom);

      // Number of events per channel and universe - Single Differential (mumom) - Flux Multisim
      hmap_trkpmom_flux_multisim_bs["total"]["nominal"] = new TH1D("h_trkpmom_total_flux_mulinominal", "; Track length;", n_bins_pmom, bins_pmom); // 20, 0, 2.5
      hmap_trkpmom_flux_multisim_bs["signal"]["nominal"] = new TH1D("h_trkpmom_signal_flux_mulinominal", "; Track length;", n_bins_pmom, bins_pmom);
      hmap_trkpmom_flux_multisim_bs["cosmic"]["nominal"] = new TH1D("h_trkpmom_cosmic_flux_mulinominal", "; Track length;", n_bins_pmom, bins_pmom);
      hmap_trkpmom_flux_multisim_bs["outfv"]["nominal"] = new TH1D("h_trkpmom_outfv_flux_mulinominal", "; Track length;", n_bins_pmom, bins_pmom);
      hmap_trkpmom_flux_multisim_bs["nc"]["nominal"] = new TH1D("h_trkpmom_nc_flux_mulinominal", "; Track length;", n_bins_pmom, bins_pmom);
      hmap_trkpmom_flux_multisim_bs["anumu"]["nominal"] = new TH1D("h_trkpmom_anumu_flux_mulinominal", "; Track length;", n_bins_pmom, bins_pmom);
      hmap_trkpmom_flux_multisim_bs["nue"]["nominal"] = new TH1D("h_trkpmom_nue_flux_mulinominal", "; Track length;", n_bins_pmom, bins_pmom);
      hmap_trkpmom_flux_multisim_bs["cc_other"]["nominal"] = new TH1D("h_trkpmom_ccother_flux_mulinominal", "; Track length;", n_bins_pmom, bins_pmom);

      // Number of events per channel and universe - Single Differential (mumom) - Flux Multisim
      hmap_thetamup_flux_multisim_bs["total"]["nominal"] = new TH1D("h_thetamup_total_flux_mulinominal", "; Track length;", n_bins_muptheta, bins_muptheta); // 20, 0, 2.5
      hmap_thetamup_flux_multisim_bs["signal"]["nominal"] = new TH1D("h_thetamup_signal_flux_mulinominal", "; Track length;", n_bins_muptheta, bins_muptheta);
      hmap_thetamup_flux_multisim_bs["cosmic"]["nominal"] = new TH1D("h_thetamup_cosmic_flux_mulinominal", "; Track length;", n_bins_muptheta, bins_muptheta);
      hmap_thetamup_flux_multisim_bs["outfv"]["nominal"] = new TH1D("h_thetamup_outfv_flux_mulinominal", "; Track length;", n_bins_muptheta, bins_muptheta);
      hmap_thetamup_flux_multisim_bs["nc"]["nominal"] = new TH1D("h_thetamup_nc_flux_mulinominal", "; Track length;", n_bins_muptheta, bins_muptheta);
      hmap_thetamup_flux_multisim_bs["anumu"]["nominal"] = new TH1D("h_thetamup_anumu_flux_mulinominal", "; Track length;", n_bins_muptheta, bins_muptheta);
      hmap_thetamup_flux_multisim_bs["nue"]["nominal"] = new TH1D("h_thetamup_nue_flux_mulinominal", "; Track length;", n_bins_muptheta, bins_muptheta);
      hmap_thetamup_flux_multisim_bs["cc_other"]["nominal"] = new TH1D("h_thetamup_ccother_flux_mulinominal", "; Track length;", n_bins_muptheta, bins_muptheta);

       // Number of events per channel and universe - Single Differential (mumom) - Extra Syst
      hmap_trkmom_extra_syst_multisim_bs["total"]["nominal"] = new TH1D("h_trkmom_total_extra_syst_multisim_nominal", "; Track length;", n_bins_mumom, bins_mumom); // 20, 0, 2.5
      hmap_trkmom_extra_syst_multisim_bs["signal"]["nominal"] = new TH1D("h_trkmom_signal_extra_syst_multisim_nominal", "; Track length;", n_bins_mumom, bins_mumom);
      hmap_trkmom_extra_syst_multisim_bs["cosmic"]["nominal"] = new TH1D("h_trkmom_cosmic_extra_syst_multisim_nominal", "; Track length;", n_bins_mumom, bins_mumom);
      hmap_trkmom_extra_syst_multisim_bs["outfv"]["nominal"] = new TH1D("h_trkmom_outfv_extra_syst_multisim_nominal", "; Track length;", n_bins_mumom, bins_mumom);
      hmap_trkmom_extra_syst_multisim_bs["nc"]["nominal"] = new TH1D("h_trkmom_nc_extra_syst_multisim_nominal", "; Track length;", n_bins_mumom, bins_mumom);
      hmap_trkmom_extra_syst_multisim_bs["anumu"]["nominal"] = new TH1D("h_trkmom_anumu_extra_syst_multisim_nominal", "; Track length;", n_bins_mumom, bins_mumom);
      hmap_trkmom_extra_syst_multisim_bs["nue"]["nominal"] = new TH1D("h_trkmom_nue_extra_syst_multisim_nominal", "; Track length;", n_bins_mumom, bins_mumom);
      hmap_trkmom_extra_syst_multisim_bs["cc_other"]["nominal"] = new TH1D("h_trkmom_ccother_extra_syst_multisim_nominal", "; Track length;", n_bins_mumom, bins_mumom);

       // Number of events per channel and universe - Single Differential (pmom) - Extra Syst
      hmap_trkpmom_extra_syst_multisim_bs["total"]["nominal"] = new TH1D("h_trkpmom_total_extra_syst_multisim_nominal", "; Track length;", n_bins_pmom, bins_pmom); // 20, 0, 2.5
      hmap_trkpmom_extra_syst_multisim_bs["signal"]["nominal"] = new TH1D("h_trkpmom_signal_extra_syst_multisim_nominal", "; Track length;", n_bins_pmom, bins_pmom);
      hmap_trkpmom_extra_syst_multisim_bs["cosmic"]["nominal"] = new TH1D("h_trkpmom_cosmic_extra_syst_multisim_nominal", "; Track length;", n_bins_pmom, bins_pmom);
      hmap_trkpmom_extra_syst_multisim_bs["outfv"]["nominal"] = new TH1D("h_trkpmom_outfv_extra_syst_multisim_nominal", "; Track length;", n_bins_pmom, bins_pmom);
      hmap_trkpmom_extra_syst_multisim_bs["nc"]["nominal"] = new TH1D("h_trkpmom_nc_extra_syst_multisim_nominal", "; Track length;", n_bins_pmom, bins_pmom);
      hmap_trkpmom_extra_syst_multisim_bs["anumu"]["nominal"] = new TH1D("h_trkpmom_anumu_extra_syst_multisim_nominal", "; Track length;", n_bins_pmom, bins_pmom);
      hmap_trkpmom_extra_syst_multisim_bs["nue"]["nominal"] = new TH1D("h_trkpmom_nue_extra_syst_multisim_nominal", "; Track length;", n_bins_pmom, bins_pmom);
      hmap_trkpmom_extra_syst_multisim_bs["cc_other"]["nominal"] = new TH1D("h_trkpmom_ccother_extra_syst_multisim_nominal", "; Track length;", n_bins_pmom, bins_pmom);

       // Number of events per channel and universe - Single Differential (thetamup) - Extra Syst
      hmap_thetamup_extra_syst_multisim_bs["total"]["nominal"] = new TH1D("h_thetamup_total_extra_syst_multisim_nominal", "; Track length;", n_bins_muptheta, bins_muptheta); // 20, 0, 2.5
      hmap_thetamup_extra_syst_multisim_bs["signal"]["nominal"] = new TH1D("h_thetamup_signal_extra_syst_multisim_nominal", "; Track length;", n_bins_muptheta, bins_muptheta);
      hmap_thetamup_extra_syst_multisim_bs["cosmic"]["nominal"] = new TH1D("h_thetamup_cosmic_extra_syst_multisim_nominal", "; Track length;", n_bins_muptheta, bins_muptheta);
      hmap_thetamup_extra_syst_multisim_bs["outfv"]["nominal"] = new TH1D("h_thetamup_outfv_extra_syst_multisim_nominal", "; Track length;", n_bins_muptheta, bins_muptheta);
      hmap_thetamup_extra_syst_multisim_bs["nc"]["nominal"] = new TH1D("h_thetamup_nc_extra_syst_multisim_nominal", "; Track length;", n_bins_muptheta, bins_muptheta);
      hmap_thetamup_extra_syst_multisim_bs["anumu"]["nominal"] = new TH1D("h_thetamup_anumu_extra_syst_multisim_nominal", "; Track length;", n_bins_muptheta, bins_muptheta);
      hmap_thetamup_extra_syst_multisim_bs["nue"]["nominal"] = new TH1D("h_thetamup_nue_extra_syst_multisim_nominal", "; Track length;", n_bins_muptheta, bins_muptheta);
      hmap_thetamup_extra_syst_multisim_bs["cc_other"]["nominal"] = new TH1D("h_thetamup_ccother_extra_syst_multisim_nominal", "; Track length;", n_bins_muptheta, bins_muptheta);


      // Number of events per channel and universe - Single Differential (mumom) - MC Stat
      hmap_trkmom_mc_stat_multisim_bs["total"]["nominal"] = new TH1D("h_trkmom_total_mc_stat_multisim_nominal", "; Track length;", n_bins_mumom, bins_mumom); // 20, 0, 2.5
      hmap_trkmom_mc_stat_multisim_bs["signal"]["nominal"] = new TH1D("h_trkmom_signal_mc_stat_multisim_nominal", "; Track length;", n_bins_mumom, bins_mumom);
      hmap_trkmom_mc_stat_multisim_bs["cosmic"]["nominal"] = new TH1D("h_trkmom_cosmic_mc_stat_multisim_nominal", "; Track length;", n_bins_mumom, bins_mumom);
      hmap_trkmom_mc_stat_multisim_bs["outfv"]["nominal"] = new TH1D("h_trkmom_outfv_mc_stat_multisim_nominal", "; Track length;", n_bins_mumom, bins_mumom);
      hmap_trkmom_mc_stat_multisim_bs["nc"]["nominal"] = new TH1D("h_trkmom_nc_mc_stat_multisim_nominal", "; Track length;", n_bins_mumom, bins_mumom);
      hmap_trkmom_mc_stat_multisim_bs["anumu"]["nominal"] = new TH1D("h_trkmom_anumu_mc_stat_multisim_nominal", "; Track length;", n_bins_mumom, bins_mumom);
      hmap_trkmom_mc_stat_multisim_bs["nue"]["nominal"] = new TH1D("h_trkmom_nue_mc_stat_multisim_nominal", "; Track length;", n_bins_mumom, bins_mumom);
      hmap_trkmom_mc_stat_multisim_bs["cc_other"]["nominal"] = new TH1D("h_trkmom_ccother_mc_stat_multisim_nominal", "; Track length;", n_bins_mumom, bins_mumom);

      // Number of events per channel and universe - Single Differential (pmom) - MC Stat
      hmap_trkpmom_mc_stat_multisim_bs["total"]["nominal"] = new TH1D("h_trkpmom_total_mc_stat_multisim_nominal", "; Track length;", n_bins_pmom, bins_pmom); // 20, 0, 2.5
      hmap_trkpmom_mc_stat_multisim_bs["signal"]["nominal"] = new TH1D("h_trkpmom_signal_mc_stat_multisim_nominal", "; Track length;", n_bins_pmom, bins_pmom);
      hmap_trkpmom_mc_stat_multisim_bs["cosmic"]["nominal"] = new TH1D("h_trkpmom_cosmic_mc_stat_multisim_nominal", "; Track length;", n_bins_pmom, bins_pmom);
      hmap_trkpmom_mc_stat_multisim_bs["outfv"]["nominal"] = new TH1D("h_trkpmom_outfv_mc_stat_multisim_nominal", "; Track length;", n_bins_pmom, bins_pmom);
      hmap_trkpmom_mc_stat_multisim_bs["nc"]["nominal"] = new TH1D("h_trkpmom_nc_mc_stat_multisim_nominal", "; Track length;", n_bins_pmom, bins_pmom);
      hmap_trkpmom_mc_stat_multisim_bs["anumu"]["nominal"] = new TH1D("h_trkpmom_anumu_mc_stat_multisim_nominal", "; Track length;", n_bins_pmom, bins_pmom);
      hmap_trkpmom_mc_stat_multisim_bs["nue"]["nominal"] = new TH1D("h_trkpmom_nue_mc_stat_multisim_nominal", "; Track length;", n_bins_pmom, bins_pmom);
      hmap_trkpmom_mc_stat_multisim_bs["cc_other"]["nominal"] = new TH1D("h_trkpmom_ccother_mc_stat_multisim_nominal", "; Track length;", n_bins_pmom, bins_pmom);

      // Number of events per channel and universe - Single Differential (thetamup) - MC Stat
      hmap_thetamup_mc_stat_multisim_bs["total"]["nominal"] = new TH1D("h_thetamup_total_mc_stat_multisim_nominal", "; Track length;", n_bins_muptheta, bins_muptheta); // 20, 0, 2.5
      hmap_thetamup_mc_stat_multisim_bs["signal"]["nominal"] = new TH1D("h_thetamup_signal_mc_stat_multisim_nominal", "; Track length;", n_bins_muptheta, bins_muptheta);
      hmap_thetamup_mc_stat_multisim_bs["cosmic"]["nominal"] = new TH1D("h_thetamup_cosmic_mc_stat_multisim_nominal", "; Track length;", n_bins_muptheta, bins_muptheta);
      hmap_thetamup_mc_stat_multisim_bs["outfv"]["nominal"] = new TH1D("h_thetamup_outfv_mc_stat_multisim_nominal", "; Track length;", n_bins_muptheta, bins_muptheta);
      hmap_thetamup_mc_stat_multisim_bs["nc"]["nominal"] = new TH1D("h_thetamup_nc_mc_stat_multisim_nominal", "; Track length;", n_bins_muptheta, bins_muptheta);
      hmap_thetamup_mc_stat_multisim_bs["anumu"]["nominal"] = new TH1D("h_thetamup_anumu_mc_stat_multisim_nominal", "; Track length;", n_bins_muptheta, bins_muptheta);
      hmap_thetamup_mc_stat_multisim_bs["nue"]["nominal"] = new TH1D("h_thetamup_nue_mc_stat_multisim_nominal", "; Track length;", n_bins_muptheta, bins_muptheta);
      hmap_thetamup_mc_stat_multisim_bs["cc_other"]["nominal"] = new TH1D("h_thetamup_ccother_mc_stat_multisim_nominal", "; Track length;", n_bins_muptheta, bins_muptheta);




      // Reco to true histograms for every universe - Single Differential (mumom) - Genie Multisim
      bs_genie_multisim_true_reco_mumom = new BootstrapTH2D("bs_genie_multisim_true_reco_mumom", ";Muon Momentum (Truth) [GeV]; Muon Momentum (MCS) [GeV]", n_bins_mumom, bins_mumom, n_bins_mumom, bins_mumom);
      bs_flux_multisim_true_reco_mumom = new BootstrapTH2D("bs_flux_multisim_true_reco_mumom", ";Muon Momentum (Truth) [GeV]; Muon Momentum (MCS) [GeV]", n_bins_mumom, bins_mumom, n_bins_mumom, bins_mumom);
      bs_extra_syst_true_reco_mumom = new BootstrapTH2D("bs_extra_syst_true_reco_mumom", ";Muon Momentum (Truth) [GeV]; Muon Momentum (MCS) [GeV]", n_bins_mumom, bins_mumom, n_bins_mumom, bins_mumom);
      bs_mc_stat_multisim_true_reco_mumom = new BootstrapTH2D("bs_mc_stat_multisim_true_reco_mumom", ";Muon Momentum (Truth) [GeV]; Muon Momentum (MCS) [GeV]", n_bins_mumom, bins_mumom, n_bins_mumom, bins_mumom);

      bs_genie_multisim_true_reco_pmom = new BootstrapTH2D("bs_genie_multisim_true_reco_pmom", ";Proton Momentum (Truth) [GeV]; Proton Momentum (MCS) [GeV]", n_bins_pmom, bins_pmom, n_bins_pmom, bins_pmom);
      bs_flux_multisim_true_reco_pmom = new BootstrapTH2D("bs_flux_multisim_true_reco_pmom", ";Proton Momentum (Truth) [GeV]; Proton Momentum (MCS) [GeV]", n_bins_pmom, bins_pmom, n_bins_pmom, bins_pmom);
      bs_extra_syst_true_reco_pmom = new BootstrapTH2D("bs_extra_syst_true_reco_pmom", ";Proton Momentum (Truth) [GeV]; Proton Momentum (MCS) [GeV]", n_bins_pmom, bins_pmom, n_bins_pmom, bins_pmom);
      bs_mc_stat_multisim_true_reco_pmom = new BootstrapTH2D("bs_mc_stat_multisim_true_reco_pmom", ";Proton Momentum (Truth) [GeV]; Proton Momentum (MCS) [GeV]", n_bins_pmom, bins_pmom, n_bins_pmom, bins_pmom);

      bs_genie_multisim_true_reco_thetamup = new BootstrapTH2D("bs_genie_multisim_true_reco_thetamup", ";Proton Momentum (Truth) [GeV]; Proton Momentum (MCS) [GeV]", n_bins_muptheta, bins_muptheta, n_bins_muptheta, bins_muptheta);
      bs_flux_multisim_true_reco_thetamup = new BootstrapTH2D("bs_flux_multisim_true_reco_thetamup", ";Proton Momentum (Truth) [GeV]; Proton Momentum (MCS) [GeV]", n_bins_muptheta, bins_muptheta, n_bins_muptheta, bins_muptheta);
      bs_extra_syst_true_reco_thetamup = new BootstrapTH2D("bs_extra_syst_true_reco_thetamup", ";Proton Momentum (Truth) [GeV]; Proton Momentum (MCS) [GeV]", n_bins_muptheta, bins_muptheta, n_bins_muptheta, bins_muptheta);
      bs_mc_stat_multisim_true_reco_thetamup = new BootstrapTH2D("bs_mc_stat_multisim_true_reco_thetamup", ";Proton Momentum (Truth) [GeV]; Proton Momentum (MCS) [GeV]", n_bins_muptheta, bins_muptheta, n_bins_muptheta, bins_muptheta);



      //
      // Single differential (muangle)
      //

      // Efficiency - Single Differential (muangle)
      h_eff_muangle_num = new TH1D("h_eff_muangle_num", "h_eff_muangle_num", n_bins_mucostheta, bins_mucostheta);
      h_eff_muangle_den = new TH1D("h_eff_muangle_den", "h_eff_muangle_den", n_bins_mucostheta, bins_mucostheta);

      // Efficiency - Single Differential (pangle)
      h_eff_pangle_num = new TH1D("h_eff_pangle_num", "h_eff_pangle_num", n_bins_pcostheta, bins_pcostheta);
      h_eff_pangle_den = new TH1D("h_eff_pangle_den", "h_eff_pangle_den", n_bins_pcostheta, bins_pcostheta);


      // Reco to true histograms - Single Differential (muangle)
      h_true_reco_costheta= new TH2D("h_true_reco_costheta", ";Muon cos(#theta) (Truth) [GeV]; Muon cos(#theta) (MCS) [GeV]", n_bins_mucostheta, bins_mucostheta, n_bins_mucostheta, bins_mucostheta);

      // Reco to true histograms - Single Differential (pangle)
      h_true_reco_pcostheta= new TH2D("h_true_reco_pcostheta", ";Muon cos(#theta) (Truth) [GeV]; Proton cos(#theta) (MCS) [GeV]", n_bins_pcostheta, bins_pcostheta, n_bins_pcostheta, bins_pcostheta);


      // Number of events per channel - Single Differential (muangle)
      hmap_trktheta["total"] = new TH1D("h_trktheta_total", "; Track cos(#theta);", n_bins_mucostheta, bins_mucostheta);
      hmap_trktheta["signal"] = new TH1D("h_trktheta_signal", "; Track cos(#theta);", n_bins_mucostheta, bins_mucostheta);
      hmap_trktheta["cosmic"] = new TH1D("h_trktheta_cosmic", "; Track cos(#theta);", n_bins_mucostheta, bins_mucostheta);
      hmap_trktheta["outfv"] = new TH1D("h_trktheta_outfv", "; Track cos(#theta);", n_bins_mucostheta, bins_mucostheta);
      hmap_trktheta["nc"] = new TH1D("h_trktheta_nc", "; Track cos(#theta);", n_bins_mucostheta, bins_mucostheta);
      hmap_trktheta["anumu"] = new TH1D("h_trktheta_anumu", "; Track cos(#theta);", n_bins_mucostheta, bins_mucostheta);
      hmap_trktheta["nue"] = new TH1D("h_trktheta_nue", "; Track cos(#theta);", n_bins_mucostheta, bins_mucostheta);
      hmap_trktheta["cosmic_stopmu"] = new TH1D("h_trktheta_cosmic_stopmu", "; Track cos(#theta);", n_bins_mucostheta, bins_mucostheta);
      hmap_trktheta["cosmic_nostopmu"] = new TH1D("h_trktheta_cosmic_nostopmu", "; Track cos(#theta);", n_bins_mucostheta, bins_mucostheta);
      hmap_trktheta["outfv_stopmu"] = new TH1D("h_trktheta_outfv_stopmu", "; Track cos(#theta);", n_bins_mucostheta, bins_mucostheta);
      hmap_trktheta["outfv_nostopmu"] = new TH1D("h_trktheta_outfv_nostopmu", "; Track cos(#theta);", n_bins_mucostheta, bins_mucostheta);
      hmap_trktheta["nc_proton"] = new TH1D("h_trktheta_nc_proton", "; Track cos(#theta);", n_bins_mucostheta, bins_mucostheta);
      hmap_trktheta["nc_pion"] = new TH1D("h_trktheta_nc_pion", "; Track cos(#theta);", n_bins_mucostheta, bins_mucostheta);
      hmap_trktheta["nc_other"] = new TH1D("h_trktheta_nc_other", "; Track cos(#theta);", n_bins_mucostheta, bins_mucostheta);
      hmap_trktheta["signal_stopmu"] = new TH1D("h_trktheta_signal_stopmu", "; Track cos(#theta);", n_bins_mucostheta, bins_mucostheta);
      hmap_trktheta["signal_nostopmu"] = new TH1D("h_trktheta_signal_nostopmu", "; Track cos(#theta);", n_bins_mucostheta, bins_mucostheta);
      hmap_trktheta["cc_other"] = new TH1D("h_trktheta_ccother", "; Track cos(#theta);", n_bins_mucostheta, bins_mucostheta);

      // Number of events per channel - Single Differential (muangle)
      hmap_trkptheta["total"] = new TH1D("h_trkptheta_total", "; Track cos(#theta);", n_bins_pcostheta, bins_pcostheta);
      hmap_trkptheta["signal"] = new TH1D("h_trkptheta_signal", "; Track cos(#theta);", n_bins_pcostheta, bins_pcostheta);
      hmap_trkptheta["cosmic"] = new TH1D("h_trkptheta_cosmic", "; Track cos(#theta);", n_bins_pcostheta, bins_pcostheta);
      hmap_trkptheta["outfv"] = new TH1D("h_trkptheta_outfv", "; Track cos(#theta);", n_bins_pcostheta, bins_pcostheta);
      hmap_trkptheta["nc"] = new TH1D("h_trkptheta_nc", "; Track cos(#theta);", n_bins_pcostheta, bins_pcostheta);
      hmap_trkptheta["anumu"] = new TH1D("h_trkptheta_anumu", "; Track cos(#theta);", n_bins_pcostheta, bins_pcostheta);
      hmap_trkptheta["nue"] = new TH1D("h_trkptheta_nue", "; Track cos(#theta);", n_bins_pcostheta, bins_pcostheta);
      hmap_trkptheta["cosmic_stopmu"] = new TH1D("h_trkptheta_cosmic_stopmu", "; Track cos(#theta);", n_bins_pcostheta, bins_pcostheta);
      hmap_trkptheta["cosmic_nostopmu"] = new TH1D("h_trkptheta_cosmic_nostopmu", "; Track cos(#theta);", n_bins_pcostheta, bins_pcostheta);
      hmap_trkptheta["outfv_stopmu"] = new TH1D("h_trkptheta_outfv_stopmu", "; Track cos(#theta);", n_bins_pcostheta, bins_pcostheta);
      hmap_trkptheta["outfv_nostopmu"] = new TH1D("h_trkptheta_outfv_nostopmu", "; Track cos(#theta);", n_bins_pcostheta, bins_pcostheta);
      hmap_trkptheta["nc_proton"] = new TH1D("h_trkptheta_nc_proton", "; Track cos(#theta);", n_bins_pcostheta, bins_pcostheta);
      hmap_trkptheta["nc_pion"] = new TH1D("h_trkptheta_nc_pion", "; Track cos(#theta);", n_bins_pcostheta, bins_pcostheta);
      hmap_trkptheta["nc_other"] = new TH1D("h_trkptheta_nc_other", "; Track cos(#theta);", n_bins_pcostheta, bins_pcostheta);
      hmap_trkptheta["signal_stopmu"] = new TH1D("h_trkptheta_signal_stopmu", "; Track cos(#theta);", n_bins_pcostheta, bins_pcostheta);
      hmap_trkptheta["signal_nostopmu"] = new TH1D("h_trkptheta_signal_nostopmu", "; Track cos(#theta);", n_bins_pcostheta, bins_pcostheta);
      hmap_trkptheta["cc_other"] = new TH1D("h_trkptheta_ccother", "; Track cos(#theta);", n_bins_pcostheta, bins_pcostheta);

       // Efficiency - Single Differential (muangle) - Multisim
      bs_genie_multisim_eff_muangle_num = new BootstrapTH1D("bs_genie_multisim_eff_muangle_num", "bs_genie_multisim_eff_muangle_num_title", n_bins_mucostheta, bins_mucostheta);
      bs_genie_multisim_eff_muangle_den = new BootstrapTH1D("bs_genie_multisim_eff_muangle_den", "bs_genie_multisim_eff_muangle_den_title", n_bins_mucostheta, bins_mucostheta);
      bs_flux_multisim_eff_muangle_num = new BootstrapTH1D("bs_flux_multisim_eff_muangle_num", "bs_flux_multisim_eff_muangle_num_title", n_bins_mucostheta, bins_mucostheta);
      bs_flux_multisim_eff_muangle_den = new BootstrapTH1D("bs_flux_multisim_eff_muangle_den", "bs_flux_multisim_eff_muangle_den_title", n_bins_mucostheta, bins_mucostheta);
      bs_extra_syst_multisim_eff_muangle_num = new BootstrapTH1D("bs_extra_syst_multisim_eff_muangle_num", "bs_extra_syst_eff_muangle_num_title", n_bins_mucostheta, bins_mucostheta);
      bs_extra_syst_multisim_eff_muangle_den = new BootstrapTH1D("bs_extra_syst_multisim_eff_muangle_den", "bs_extra_syst_eff_muangle_den_title", n_bins_mucostheta, bins_mucostheta);
      bs_mc_stat_multisim_eff_muangle_num = new BootstrapTH1D("bs_mc_stat_multisim_eff_muangle_num", "bs_mc_stat_eff_muangle_num_title", n_bins_mucostheta, bins_mucostheta);
      bs_mc_stat_multisim_eff_muangle_den = new BootstrapTH1D("bs_mc_stat_multisim_eff_muangle_den", "bs_mc_stat_eff_muangle_den_title", n_bins_mucostheta, bins_mucostheta);

       // Efficiency - Single Differential (muangle) - Multisim
      bs_genie_multisim_eff_pangle_num = new BootstrapTH1D("bs_genie_multisim_eff_pangle_num", "bs_genie_multisim_eff_pangle_num_title", n_bins_pcostheta, bins_pcostheta);
      bs_genie_multisim_eff_pangle_den = new BootstrapTH1D("bs_genie_multisim_eff_pangle_den", "bs_genie_multisim_eff_pangle_den_title", n_bins_pcostheta, bins_pcostheta);
      bs_flux_multisim_eff_pangle_num = new BootstrapTH1D("bs_flux_multisim_eff_pangle_num", "bs_flux_multisim_eff_pangle_num_title", n_bins_pcostheta, bins_pcostheta);
      bs_flux_multisim_eff_pangle_den = new BootstrapTH1D("bs_flux_multisim_eff_pangle_den", "bs_flux_multisim_eff_pangle_den_title", n_bins_pcostheta, bins_pcostheta);
      bs_extra_syst_multisim_eff_pangle_num = new BootstrapTH1D("bs_extra_syst_multisim_eff_pangle_num", "bs_extra_syst_eff_pangle_num_title", n_bins_pcostheta, bins_pcostheta);
      bs_extra_syst_multisim_eff_pangle_den = new BootstrapTH1D("bs_extra_syst_multisim_eff_pangle_den", "bs_extra_syst_eff_pangle_den_title", n_bins_pcostheta, bins_pcostheta);
      bs_mc_stat_multisim_eff_pangle_num = new BootstrapTH1D("bs_mc_stat_multisim_eff_pangle_num", "bs_mc_stat_eff_pangle_num_title", n_bins_pcostheta, bins_pcostheta);
      bs_mc_stat_multisim_eff_pangle_den = new BootstrapTH1D("bs_mc_stat_multisim_eff_pangle_den", "bs_mc_stat_eff_pangle_den_title", n_bins_pcostheta, bins_pcostheta);
      //===============================================================================================================================================================
      // Number of events per channel and universe - Single Differential (muangle) - Genie Multisim
      hmap_trkangle_genie_multisim_bs["total"]["nominal"] = new TH1D("h_trkangle_total_genie_mulinominal", "; Track length;", n_bins_mucostheta, bins_mucostheta); // 20, 0, 2.5
      hmap_trkangle_genie_multisim_bs["signal"]["nominal"] = new TH1D("h_trkangle_signal_genie_mulinominal", "; Track length;", n_bins_mucostheta, bins_mucostheta);
      hmap_trkangle_genie_multisim_bs["cosmic"]["nominal"] = new TH1D("h_trkangle_cosmic_genie_mulinominal", "; Track length;", n_bins_mucostheta, bins_mucostheta);
      hmap_trkangle_genie_multisim_bs["outfv"]["nominal"] = new TH1D("h_trkangle_outfv_genie_mulinominal", "; Track length;", n_bins_mucostheta, bins_mucostheta);
      hmap_trkangle_genie_multisim_bs["nc"]["nominal"] = new TH1D("h_trkangle_nc_genie_mulinominal", "; Track length;", n_bins_mucostheta, bins_mucostheta);
      hmap_trkangle_genie_multisim_bs["anumu"]["nominal"] = new TH1D("h_trkangle_anumu_genie_mulinominal", "; Track length;", n_bins_mucostheta, bins_mucostheta);
      hmap_trkangle_genie_multisim_bs["nue"]["nominal"] = new TH1D("h_trkangle_nue_genie_mulinominal", "; Track length;", n_bins_mucostheta, bins_mucostheta);
      hmap_trkangle_genie_multisim_bs["cc_other"]["nominal"] = new TH1D("h_trkangle_ccother_genie_mulinominal", "; Track length;", n_bins_mucostheta, bins_mucostheta);
      // Number of events per channel and universe - Single Differential (muangle) - Flux Multisim
      hmap_trkangle_flux_multisim_bs["total"]["nominal"] = new TH1D("h_trkangle_total_flux_mulinominal", "; Track length;", n_bins_mucostheta, bins_mucostheta); // 20, 0, 2.5
      hmap_trkangle_flux_multisim_bs["signal"]["nominal"] = new TH1D("h_trkangle_signal_flux_mulinominal", "; Track length;", n_bins_mucostheta, bins_mucostheta);
      hmap_trkangle_flux_multisim_bs["cosmic"]["nominal"] = new TH1D("h_trkangle_cosmic_flux_mulinominal", "; Track length;", n_bins_mucostheta, bins_mucostheta);
      hmap_trkangle_flux_multisim_bs["outfv"]["nominal"] = new TH1D("h_trkangle_outfv_flux_mulinominal", "; Track length;", n_bins_mucostheta, bins_mucostheta);
      hmap_trkangle_flux_multisim_bs["nc"]["nominal"] = new TH1D("h_trkangle_nc_flux_mulinominal", "; Track length;", n_bins_mucostheta, bins_mucostheta);
      hmap_trkangle_flux_multisim_bs["anumu"]["nominal"] = new TH1D("h_trkangle_anumu_flux_mulinominal", "; Track length;", n_bins_mucostheta, bins_mucostheta);
      hmap_trkangle_flux_multisim_bs["nue"]["nominal"] = new TH1D("h_trkangle_nue_flux_mulinominal", "; Track length;", n_bins_mucostheta, bins_mucostheta);
      hmap_trkangle_flux_multisim_bs["cc_other"]["nominal"] = new TH1D("h_trkangle_ccother_flux_mulinominal", "; Track length;", n_bins_mucostheta, bins_mucostheta);
      // Number of events per channel and universe - Single Differential (muangle) - Extra Syst
      hmap_trkangle_extra_syst_multisim_bs["total"]["nominal"] = new TH1D("h_trkangle_total_extra_syst_multisim_nominal", "; Track length;", n_bins_mucostheta, bins_mucostheta); // 20, 0, 2.5
      hmap_trkangle_extra_syst_multisim_bs["signal"]["nominal"] = new TH1D("h_trkangle_signal_extra_syst_multisim_nominal", "; Track length;", n_bins_mucostheta, bins_mucostheta);
      hmap_trkangle_extra_syst_multisim_bs["cosmic"]["nominal"] = new TH1D("h_trkangle_cosmic_extra_syst_multisim_nominal", "; Track length;", n_bins_mucostheta, bins_mucostheta);
      hmap_trkangle_extra_syst_multisim_bs["outfv"]["nominal"] = new TH1D("h_trkangle_outfv_extra_syst_multisim_nominal", "; Track length;", n_bins_mucostheta, bins_mucostheta);
      hmap_trkangle_extra_syst_multisim_bs["nc"]["nominal"] = new TH1D("h_trkangle_nc_extra_syst_multisim_nominal", "; Track length;", n_bins_mucostheta, bins_mucostheta);
      hmap_trkangle_extra_syst_multisim_bs["anumu"]["nominal"] = new TH1D("h_trkangle_anumu_extra_syst_multisim_nominal", "; Track length;", n_bins_mucostheta, bins_mucostheta);
      hmap_trkangle_extra_syst_multisim_bs["nue"]["nominal"] = new TH1D("h_trkangle_nue_extra_syst_multisim_nominal", "; Track length;", n_bins_mucostheta, bins_mucostheta);
      hmap_trkangle_extra_syst_multisim_bs["cc_other"]["nominal"] = new TH1D("h_trkangle_ccother_extra_syst_multisim_nominal", "; Track length;", n_bins_mucostheta, bins_mucostheta);
      // Number of events per channel and universe - Single Differential (muangle) - MC Stat
      hmap_trkangle_mc_stat_multisim_bs["total"]["nominal"] = new TH1D("h_trkangle_total_mc_stat_multisim_nominal", "; Track length;", n_bins_mucostheta, bins_mucostheta); // 20, 0, 2.5
      hmap_trkangle_mc_stat_multisim_bs["signal"]["nominal"] = new TH1D("h_trkangle_signal_mc_stat_multisim_nominal", "; Track length;", n_bins_mucostheta, bins_mucostheta);
      hmap_trkangle_mc_stat_multisim_bs["cosmic"]["nominal"] = new TH1D("h_trkangle_cosmic_mc_stat_multisim_nominal", "; Track length;", n_bins_mucostheta, bins_mucostheta);
      hmap_trkangle_mc_stat_multisim_bs["outfv"]["nominal"] = new TH1D("h_trkangle_outfv_mc_stat_multisim_nominal", "; Track length;", n_bins_mucostheta, bins_mucostheta);
      hmap_trkangle_mc_stat_multisim_bs["nc"]["nominal"] = new TH1D("h_trkangle_nc_mc_stat_multisim_nominal", "; Track length;", n_bins_mucostheta, bins_mucostheta);
      hmap_trkangle_mc_stat_multisim_bs["anumu"]["nominal"] = new TH1D("h_trkangle_anumu_mc_stat_multisim_nominal", "; Track length;", n_bins_mucostheta, bins_mucostheta);
      hmap_trkangle_mc_stat_multisim_bs["nue"]["nominal"] = new TH1D("h_trkangle_nue_mc_stat_multisim_nominal", "; Track length;", n_bins_mucostheta, bins_mucostheta);
      hmap_trkangle_mc_stat_multisim_bs["cc_other"]["nominal"] = new TH1D("h_trkangle_ccother_mc_stat_multisim_nominal", "; Track length;", n_bins_mucostheta, bins_mucostheta);
      //=============================================================================================================================================================================
      // Number of events per channel and universe - Single Differential (muangle) - Genie Multisim
      hmap_trkpangle_genie_multisim_bs["total"]["nominal"] = new TH1D("h_trkpangle_total_genie_mulinominal", "; Track length;", n_bins_pcostheta, bins_pcostheta); // 20, 0, 2.5
      hmap_trkpangle_genie_multisim_bs["signal"]["nominal"] = new TH1D("h_trkpangle_signal_genie_mulinominal", "; Track length;", n_bins_pcostheta, bins_pcostheta);
      hmap_trkpangle_genie_multisim_bs["cosmic"]["nominal"] = new TH1D("h_trkpangle_cosmic_genie_mulinominal", "; Track length;", n_bins_pcostheta, bins_pcostheta);
      hmap_trkpangle_genie_multisim_bs["outfv"]["nominal"] = new TH1D("h_trkpangle_outfv_genie_mulinominal", "; Track length;", n_bins_pcostheta, bins_pcostheta);
      hmap_trkpangle_genie_multisim_bs["nc"]["nominal"] = new TH1D("h_trkpangle_nc_genie_mulinominal", "; Track length;", n_bins_pcostheta, bins_pcostheta);
      hmap_trkpangle_genie_multisim_bs["anumu"]["nominal"] = new TH1D("h_trkpangle_anumu_genie_mulinominal", "; Track length;", n_bins_pcostheta, bins_pcostheta);
      hmap_trkpangle_genie_multisim_bs["nue"]["nominal"] = new TH1D("h_trkpangle_nue_genie_mulinominal", "; Track length;", n_bins_pcostheta, bins_pcostheta);
      hmap_trkpangle_genie_multisim_bs["cc_other"]["nominal"] = new TH1D("h_trkpangle_ccother_genie_mulinominal", "; Track length;", n_bins_pcostheta, bins_pcostheta);
      // Number of events per channel and universe - Single Differential (muangle) - Flux Multisim
      hmap_trkpangle_flux_multisim_bs["total"]["nominal"] = new TH1D("h_trkpangle_total_flux_mulinominal", "; Track length;", n_bins_pcostheta, bins_pcostheta); // 20, 0, 2.5
      hmap_trkpangle_flux_multisim_bs["signal"]["nominal"] = new TH1D("h_trkpangle_signal_flux_mulinominal", "; Track length;", n_bins_pcostheta, bins_pcostheta);
      hmap_trkpangle_flux_multisim_bs["cosmic"]["nominal"] = new TH1D("h_trkpangle_cosmic_flux_mulinominal", "; Track length;", n_bins_pcostheta, bins_pcostheta);
      hmap_trkpangle_flux_multisim_bs["outfv"]["nominal"] = new TH1D("h_trkpangle_outfv_flux_mulinominal", "; Track length;", n_bins_pcostheta, bins_pcostheta);
      hmap_trkpangle_flux_multisim_bs["nc"]["nominal"] = new TH1D("h_trkpangle_nc_flux_mulinominal", "; Track length;", n_bins_pcostheta, bins_pcostheta);
      hmap_trkpangle_flux_multisim_bs["anumu"]["nominal"] = new TH1D("h_trkpangle_anumu_flux_mulinominal", "; Track length;", n_bins_pcostheta, bins_pcostheta);
      hmap_trkpangle_flux_multisim_bs["nue"]["nominal"] = new TH1D("h_trkpangle_nue_flux_mulinominal", "; Track length;", n_bins_pcostheta, bins_pcostheta);
      hmap_trkpangle_flux_multisim_bs["cc_other"]["nominal"] = new TH1D("h_trkpangle_ccother_flux_mulinominal", "; Track length;", n_bins_pcostheta, bins_pcostheta);
      // Number of events per channel and universe - Single Differential (muangle) - Extra Syst
      hmap_trkpangle_extra_syst_multisim_bs["total"]["nominal"] = new TH1D("h_trkpangle_total_extra_syst_multisim_nominal", "; Track length;", n_bins_pcostheta, bins_pcostheta); // 20, 0, 2.5
      hmap_trkpangle_extra_syst_multisim_bs["signal"]["nominal"] = new TH1D("h_trkpangle_signal_extra_syst_multisim_nominal", "; Track length;", n_bins_pcostheta, bins_pcostheta);
      hmap_trkpangle_extra_syst_multisim_bs["cosmic"]["nominal"] = new TH1D("h_trkpangle_cosmic_extra_syst_multisim_nominal", "; Track length;", n_bins_pcostheta, bins_pcostheta);
      hmap_trkpangle_extra_syst_multisim_bs["outfv"]["nominal"] = new TH1D("h_trkpangle_outfv_extra_syst_multisim_nominal", "; Track length;", n_bins_pcostheta, bins_pcostheta);
      hmap_trkpangle_extra_syst_multisim_bs["nc"]["nominal"] = new TH1D("h_trkpangle_nc_extra_syst_multisim_nominal", "; Track length;", n_bins_pcostheta, bins_pcostheta);
      hmap_trkpangle_extra_syst_multisim_bs["anumu"]["nominal"] = new TH1D("h_trkpangle_anumu_extra_syst_multisim_nominal", "; Track length;", n_bins_pcostheta, bins_pcostheta);
      hmap_trkpangle_extra_syst_multisim_bs["nue"]["nominal"] = new TH1D("h_trkpangle_nue_extra_syst_multisim_nominal", "; Track length;", n_bins_pcostheta, bins_pcostheta);
      hmap_trkpangle_extra_syst_multisim_bs["cc_other"]["nominal"] = new TH1D("h_trkpangle_ccother_extra_syst_multisim_nominal", "; Track length;", n_bins_pcostheta, bins_pcostheta);
      // Number of events per channel and universe - Single Differential (muangle) - MC Stat
      hmap_trkpangle_mc_stat_multisim_bs["total"]["nominal"] = new TH1D("h_trkpangle_total_mc_stat_multisim_nominal", "; Track length;", n_bins_pcostheta, bins_pcostheta); // 20, 0, 2.5
      hmap_trkpangle_mc_stat_multisim_bs["signal"]["nominal"] = new TH1D("h_trkpangle_signal_mc_stat_multisim_nominal", "; Track length;", n_bins_pcostheta, bins_pcostheta);
      hmap_trkpangle_mc_stat_multisim_bs["cosmic"]["nominal"] = new TH1D("h_trkpangle_cosmic_mc_stat_multisim_nominal", "; Track length;", n_bins_pcostheta, bins_pcostheta);
      hmap_trkpangle_mc_stat_multisim_bs["outfv"]["nominal"] = new TH1D("h_trkpangle_outfv_mc_stat_multisim_nominal", "; Track length;", n_bins_pcostheta, bins_pcostheta);
      hmap_trkpangle_mc_stat_multisim_bs["nc"]["nominal"] = new TH1D("h_trkpangle_nc_mc_stat_multisim_nominal", "; Track length;", n_bins_pcostheta, bins_pcostheta);
      hmap_trkpangle_mc_stat_multisim_bs["anumu"]["nominal"] = new TH1D("h_trkpangle_anumu_mc_stat_multisim_nominal", "; Track length;", n_bins_pcostheta, bins_pcostheta);
      hmap_trkpangle_mc_stat_multisim_bs["nue"]["nominal"] = new TH1D("h_trkpangle_nue_mc_stat_multisim_nominal", "; Track length;", n_bins_pcostheta, bins_pcostheta);
      hmap_trkpangle_mc_stat_multisim_bs["cc_other"]["nominal"] = new TH1D("h_trkpangle_ccother_mc_stat_multisim_nominal", "; Track length;", n_bins_pcostheta, bins_pcostheta);
 
      // Reco to true histograms for every universe - Single Differential (muangle) - Genie Multisim
      bs_genie_multisim_true_reco_muangle = new BootstrapTH2D("bs_genie_multisim_true_reco_muangle", ";Muon CosTheta (Truth); Muon CosTheta (Reco)", n_bins_mucostheta, bins_mucostheta, n_bins_mucostheta, bins_mucostheta);
      bs_flux_multisim_true_reco_muangle = new BootstrapTH2D("bs_flux_multisim_true_reco_muangle", ";Muon CosTheta (Truth); Muon CosTheta (Reco)", n_bins_mucostheta, bins_mucostheta, n_bins_mucostheta, bins_mucostheta);
      bs_extra_syst_true_reco_muangle = new BootstrapTH2D("bs_extra_syst_true_reco_muangle", ";Muon CosTheta (Truth); Muon CosTheta (Reco)", n_bins_mucostheta, bins_mucostheta, n_bins_mucostheta, bins_mucostheta);
      bs_mc_stat_multisim_true_reco_muangle = new BootstrapTH2D("bs_mc_stat_multisim_true_reco_muangle", ";Muon CosTheta (Truth); Muon CosTheta (Reco)", n_bins_mucostheta, bins_mucostheta, n_bins_mucostheta, bins_mucostheta);

      // Reco to true histograms for every universe - Single Differential (muangle) - Genie Multisim
      bs_genie_multisim_true_reco_pangle = new BootstrapTH2D("bs_genie_multisim_true_reco_pangle", ";Proton CosTheta (Truth); Proton CosTheta (Reco)", n_bins_pcostheta, bins_pcostheta, n_bins_pcostheta, bins_pcostheta);
      bs_flux_multisim_true_reco_pangle = new BootstrapTH2D("bs_flux_multisim_true_reco_pangle", ";Proton CosTheta (Truth); Proton CosTheta (Reco)", n_bins_pcostheta, bins_pcostheta, n_bins_pcostheta, bins_pcostheta);
      bs_extra_syst_true_reco_pangle = new BootstrapTH2D("bs_extra_syst_true_reco_pangle", ";Proton CosTheta (Truth); Proton CosTheta (Reco)", n_bins_pcostheta, bins_pcostheta, n_bins_pcostheta, bins_pcostheta);
      bs_mc_stat_multisim_true_reco_pangle = new BootstrapTH2D("bs_mc_stat_multisim_true_reco_pangle", ";Proton CosTheta (Truth); Proton CosTheta (Reco)", n_bins_pcostheta, bins_pcostheta, n_bins_pcostheta, bins_pcostheta);


      GENIE_par.push_back("genie_AGKYpT_Genie_p1");
      GENIE_par.push_back("genie_AGKYpT_Genie_m1");
      GENIE_par.push_back("genie_AGKYxF_Genie_p1");
      GENIE_par.push_back("genie_AGKYxF_Genie_m1");
      GENIE_par.push_back("genie_DISAth_Genie_p1");
      GENIE_par.push_back("genie_DISAth_Genie_m1");
      GENIE_par.push_back("genie_DISBth_Genie_p1");
      GENIE_par.push_back("genie_DISBth_Genie_m1");
      GENIE_par.push_back("genie_DISCv1u_Genie_p1");
      GENIE_par.push_back("genie_DISCv1u_Genie_m1");
      GENIE_par.push_back("genie_DISCv2u_Genie_p1");
      GENIE_par.push_back("genie_DISCv2u_Genie_m1");
      GENIE_par.push_back("genie_FermiGasModelKf_Genie_p1");
      GENIE_par.push_back("genie_FermiGasModelKf_Genie_m1");
      GENIE_par.push_back("genie_FermiGasModelSf_Genie_p1");
      GENIE_par.push_back("genie_FermiGasModelSf_Genie_m1");
      GENIE_par.push_back("genie_FormZone_Genie_p1");
      GENIE_par.push_back("genie_FormZone_Genie_m1");
      GENIE_par.push_back("genie_IntraNukeNabs_Genie_p1");
      GENIE_par.push_back("genie_IntraNukeNabs_Genie_m1");
      GENIE_par.push_back("genie_IntraNukeNcex_Genie_p1");
      GENIE_par.push_back("genie_IntraNukeNcex_Genie_m1");
      GENIE_par.push_back("genie_IntraNukeNel_Genie_p1");
      GENIE_par.push_back("genie_IntraNukeNel_Genie_m1");
      GENIE_par.push_back("genie_IntraNukeNinel_Genie_p1");
      GENIE_par.push_back("genie_IntraNukeNinel_Genie_m1");
      GENIE_par.push_back("genie_IntraNukeNmfp_Genie_p1");
      GENIE_par.push_back("genie_IntraNukeNmfp_Genie_m1");
      GENIE_par.push_back("genie_IntraNukeNpi_Genie_p1");
      GENIE_par.push_back("genie_IntraNukeNpi_Genie_m1");
      GENIE_par.push_back("genie_IntraNukePIabs_Genie_p1");
      GENIE_par.push_back("genie_IntraNukePIabs_Genie_m1");
      GENIE_par.push_back("genie_IntraNukePIcex_Genie_p1");
      GENIE_par.push_back("genie_IntraNukePIcex_Genie_m1");
      GENIE_par.push_back("genie_IntraNukePIel_Genie_p1");
      GENIE_par.push_back("genie_IntraNukePIel_Genie_m1");
      GENIE_par.push_back("genie_IntraNukePIinel_Genie_p1");
      GENIE_par.push_back("genie_IntraNukePIinel_Genie_m1");
      GENIE_par.push_back("genie_IntraNukePImfp_Genie_p1");
      GENIE_par.push_back("genie_IntraNukePImfp_Genie_m1");
      GENIE_par.push_back("genie_IntraNukePIpi_Genie_p1");
      GENIE_par.push_back("genie_IntraNukePIpi_Genie_m1");
      GENIE_par.push_back("genie_NC_Genie_p1");
      GENIE_par.push_back("genie_NC_Genie_m1");
      GENIE_par.push_back("genie_NonResRvbarp1pi_Genie_p1");
      GENIE_par.push_back("genie_NonResRvbarp1pi_Genie_m1");
      GENIE_par.push_back("genie_NonResRvbarp2pi_Genie_p1");
      GENIE_par.push_back("genie_NonResRvbarp2pi_Genie_m1");
      GENIE_par.push_back("genie_NonResRvp1pi_Genie_p1");
      GENIE_par.push_back("genie_NonResRvp1pi_Genie_m1");
      GENIE_par.push_back("genie_NonResRvp2pi_Genie_p1");
      GENIE_par.push_back("genie_NonResRvp2pi_Genie_m1");
      GENIE_par.push_back("genie_ResDecayEta_Genie_p1");
      GENIE_par.push_back("genie_ResDecayEta_Genie_m1");
      GENIE_par.push_back("genie_ResDecayGamma_Genie_p1");
      GENIE_par.push_back("genie_ResDecayGamma_Genie_m1");
      GENIE_par.push_back("genie_ResDecayTheta_Genie_p1");
      GENIE_par.push_back("genie_ResDecayTheta_Genie_m1");
      GENIE_par.push_back("genie_ccresAxial_Genie_p1");
      GENIE_par.push_back("genie_ccresAxial_Genie_m1");
      GENIE_par.push_back("genie_ccresVector_Genie_p1");
      GENIE_par.push_back("genie_ccresVector_Genie_m1");
      GENIE_par.push_back("genie_cohMA_Genie_p1");
      GENIE_par.push_back("genie_cohMA_Genie_m1");
      GENIE_par.push_back("genie_cohR0_Genie_p1");
      GENIE_par.push_back("genie_cohR0_Genie_m1");
      GENIE_par.push_back("genie_ncelAxial_Genie_p1");
      GENIE_par.push_back("genie_ncelAxial_Genie_m1");
      GENIE_par.push_back("genie_ncelEta_Genie_p1");
      GENIE_par.push_back("genie_ncelEta_Genie_m1");
      GENIE_par.push_back("genie_ncresAxial_Genie_p1");
      GENIE_par.push_back("genie_ncresAxial_Genie_m1");
      GENIE_par.push_back("genie_ncresVector_Genie_p1");
      GENIE_par.push_back("genie_ncresVector_Genie_m1");
      GENIE_par.push_back("genie_qema_Genie_p1");
      GENIE_par.push_back("genie_qema_Genie_m1");
      GENIE_par.push_back("genie_qevec_Genie_p1");
      GENIE_par.push_back("genie_qevec_Genie_m1");
      for(unsigned int i = 0; i<GENIE_par.size(); i++){
           bs_genie_pm1_true_reco_mom[GENIE_par.at(i)]= new TH2D("bs_genie_pm1_true_reco_mumom", ";Muon Momentum (Truth); Muon Momentum (Reco)", n_bins_mumom, bins_mumom, n_bins_mumom, bins_mumom);
           bs_genie_pm1_true_reco_pmom[GENIE_par.at(i)] = new TH2D("bs_genie_pm1_true_reco_pmom", ";Proton Momentum (Truth); Proton Momentum (Reco)", n_bins_pmom, bins_pmom, n_bins_pmom, bins_pmom);
           bs_genie_pm1_true_reco_muangle[GENIE_par.at(i)] = new TH2D("bs_genie_pm1_true_reco_mucostheta", ";Muon CosTheta (Truth); Muon CosTheta (Reco)", n_bins_mucostheta, bins_mucostheta, n_bins_mucostheta, bins_mucostheta);
           bs_genie_pm1_true_reco_pangle[GENIE_par.at(i)] = new TH2D("bs_genie_pm1_true_reco_pcostheta", ";Proton CosTheta (Truth); Proton CosTheta (Reco)", n_bins_pcostheta, bins_pcostheta, n_bins_pcostheta, bins_pcostheta);
           bs_genie_pm1_true_reco_thetamup[GENIE_par.at(i)] = new TH2D("bs_genie_pm1_true_reco_thetamup", ";Muon Momentum (Truth); Muon Momentum (Reco)", n_bins_muptheta, bins_muptheta, n_bins_muptheta, bins_muptheta);
       }


      //Reco to true histograms for every GENIE unisim
      bs_genie_pm1_true_reco_mom["nominal"] = new TH2D("bs_genie_pm1_true_reco_mumom", ";Muon Momentum (Truth); Muon Momentum (Reco)", n_bins_mumom, bins_mumom, n_bins_mumom, bins_mumom);
      bs_genie_pm1_true_reco_pmom["nominal"] = new TH2D("bs_genie_pm1_true_reco_pmom", ";Proton Momentum (Truth); Proton Momentum (Reco)", n_bins_pmom, bins_pmom, n_bins_pmom, bins_pmom);
      bs_genie_pm1_true_reco_muangle["nominal"] = new TH2D("bs_genie_pm1_true_reco_mucostheta", ";Muon CosTheta (Truth); Muon CosTheta (Reco)", n_bins_mucostheta, bins_mucostheta, n_bins_mucostheta, bins_mucostheta);
      bs_genie_pm1_true_reco_pangle["nominal"] = new TH2D("bs_genie_pm1_true_reco_pcostheta", ";Proton CosTheta (Truth); Proton CosTheta (Reco)", n_bins_pcostheta, bins_pcostheta, n_bins_pcostheta, bins_pcostheta);
      bs_genie_pm1_true_reco_thetamup["nominal"] = new TH2D("bs_genie_pm1_true_reco_thetamup", ";Muon Momentum (Truth); Muon Momentum (Reco)", n_bins_muptheta, bins_muptheta, n_bins_muptheta, bins_muptheta);
      
      //=============================================================================================================
      hmap_trkmom_genie_pm1_bs["total"]["nominal"] = new TH1D("h_trkmom_total_genie_pm1_nominal", "; Track momentum;", n_bins_mumom, bins_mumom);
      hmap_trkmom_genie_pm1_bs["signal"]["nominal"] = new TH1D("h_trkmom_signal_genie_pm1_nominal", "; Track momentum;", n_bins_mumom, bins_mumom);
      hmap_trkmom_genie_pm1_bs["cosmic"]["nominal"] = new TH1D("h_trkmom_cosmic_genie_pm1_nominal", "; Track momentum;", n_bins_mumom, bins_mumom);
      hmap_trkmom_genie_pm1_bs["cosmic_stopmu"]["nominal"] = new TH1D("h_trkmom_cosmic_stopmu_genie_pm1_nominal", "; Track momentum;", n_bins_mumom, bins_mumom);
      hmap_trkmom_genie_pm1_bs["cosmic_nostopmu"]["nominal"] = new TH1D("h_trkmom_cosmic_nostopmu_genie_pm1_nominal", "; Track momentum;", n_bins_mumom, bins_mumom);
      hmap_trkmom_genie_pm1_bs["outfv"]["nominal"] = new TH1D("h_trkmom_outfv_genie_pm1_nominal", "; Track momentum;", n_bins_mumom, bins_mumom);
      hmap_trkmom_genie_pm1_bs["outfv_stopmu"]["nominal"] = new TH1D("h_trkmom_outfv_stopmu_genie_pm1_nominal", "; Track momentum;", n_bins_mumom, bins_mumom);
      hmap_trkmom_genie_pm1_bs["outfv_nostopmu"]["nominal"] = new TH1D("h_trkmom_outfv_nostopmu_genie_pm1_nominal", "; Track momentum;", n_bins_mumom, bins_mumom);
      hmap_trkmom_genie_pm1_bs["nc"]["nominal"] = new TH1D("h_trkmom_nc_genie_pm1_nominal", "; Track momentum;", n_bins_mumom, bins_mumom);
      hmap_trkmom_genie_pm1_bs["nc_proton"]["nominal"] = new TH1D("h_trkmom_nc_proton_genie_pm1_nominal", "; Track momentum;", n_bins_mumom, bins_mumom);
      hmap_trkmom_genie_pm1_bs["nc_pion"]["nominal"] = new TH1D("h_trkmom_nc_pion_genie_pm1_nominal", "; Track momentum;", n_bins_mumom, bins_mumom);
      hmap_trkmom_genie_pm1_bs["nc_other"]["nominal"] = new TH1D("h_trkmom_nc_other_genie_pm1_nominal", "; Track momentum;", n_bins_mumom, bins_mumom);
      hmap_trkmom_genie_pm1_bs["anumu"]["nominal"] = new TH1D("h_trkmom_anumu_genie_pm1_nominal", "; Track momentum;", n_bins_mumom, bins_mumom);
      hmap_trkmom_genie_pm1_bs["nue"]["nominal"] = new TH1D("h_trkmom_nue_genie_pm1_nominal", "; Track momentum;", n_bins_mumom, bins_mumom);
      hmap_trkmom_genie_pm1_bs["signal_stopmu"]["nominal"] = new TH1D("h_trkmom_signal_stopmu_genie_pm1_nominal", "; Track momentum;", n_bins_mumom, bins_mumom);
      hmap_trkmom_genie_pm1_bs["signal_nostopmu"]["nominal"] = new TH1D("h_trkmom_signal_nostopmu_genie_pm1_nominal", "; Track momentum;", n_bins_mumom, bins_mumom);
      hmap_trkmom_genie_pm1_bs["cc_other"]["nominal"] = new TH1D("h_trkmom_ccother_genie_pm1_nominal", "; Track momentum;", n_bins_mumom, bins_mumom);


      hmap_trktheta_genie_pm1_bs["total"]["nominal"] = new TH1D("h_trktheta_total_genie_pm1_nominal", "; Track Angle;", n_bins_mucostheta, bins_mucostheta);
      hmap_trktheta_genie_pm1_bs["signal"]["nominal"] = new TH1D("h_trktheta_signal_genie_pm1_nominal", "; Track Angle;", n_bins_mucostheta, bins_mucostheta);
      hmap_trktheta_genie_pm1_bs["cosmic"]["nominal"] = new TH1D("h_trktheta_cosmic_genie_pm1_nominal", "; Track Angle;", n_bins_mucostheta, bins_mucostheta);
      hmap_trktheta_genie_pm1_bs["cosmic_stopmu"]["nominal"] = new TH1D("h_trktheta_cosmic_stopmu_genie_pm1_nominal", "; Track Angle;", n_bins_mucostheta, bins_mucostheta);
      hmap_trktheta_genie_pm1_bs["cosmic_nostopmu"]["nominal"] = new TH1D("h_trktheta_cosmic_nostopmu_genie_pm1_nominal", "; Track Angle;", n_bins_mucostheta, bins_mucostheta);
      hmap_trktheta_genie_pm1_bs["outfv"]["nominal"] = new TH1D("h_trktheta_outfv_genie_pm1_nominal", "; Track Angle;", n_bins_mucostheta, bins_mucostheta);
      hmap_trktheta_genie_pm1_bs["outfv_stopmu"]["nominal"] = new TH1D("h_trktheta_outfv_stopmu_genie_pm1_nominal", "; Track Angle;", n_bins_mucostheta, bins_mucostheta);
      hmap_trktheta_genie_pm1_bs["outfv_nostopmu"]["nominal"] = new TH1D("h_trktheta_outfv_nostopmu_genie_pm1_nominal", "; Track Angle;", n_bins_mucostheta, bins_mucostheta);
      hmap_trktheta_genie_pm1_bs["nc"]["nominal"] = new TH1D("h_trktheta_nc_genie_pm1_nominal", "; Track Angle;", n_bins_mucostheta, bins_mucostheta);
      hmap_trktheta_genie_pm1_bs["nc_proton"]["nominal"] = new TH1D("h_trktheta_nc_proton_genie_pm1_nominal", "; Track Angle;", n_bins_mucostheta, bins_mucostheta);
      hmap_trktheta_genie_pm1_bs["nc_pion"]["nominal"] = new TH1D("h_trktheta_nc_pion_genie_pm1_nominal", "; Track Angle;", n_bins_mucostheta, bins_mucostheta);
      hmap_trktheta_genie_pm1_bs["nc_other"]["nominal"] = new TH1D("h_trktheta_nc_other_genie_pm1_nominal", "; Track Angle;", n_bins_mucostheta, bins_mucostheta);
      hmap_trktheta_genie_pm1_bs["anumu"]["nominal"] = new TH1D("h_trktheta_anumu_genie_pm1_nominal", "; Track Angle;", n_bins_mucostheta, bins_mucostheta);
      hmap_trktheta_genie_pm1_bs["nue"]["nominal"] = new TH1D("h_trktheta_nue_genie_pm1_nominal", "; Track Angle;", n_bins_mucostheta, bins_mucostheta);
      hmap_trktheta_genie_pm1_bs["signal_stopmu"]["nominal"] = new TH1D("h_trktheta_signal_stopmu_genie_pm1_nominal", "; Track Angle;", n_bins_mucostheta, bins_mucostheta);
      hmap_trktheta_genie_pm1_bs["signal_nostopmu"]["nominal"] = new TH1D("h_trktheta_signal_nostopmu_genie_pm1_nominal", "; Track Angle;", n_bins_mucostheta, bins_mucostheta);
      hmap_trktheta_genie_pm1_bs["cc_other"]["nominal"] = new TH1D("h_trktheta_ccother_genie_pm1_nominal", "; Track Angle;", n_bins_mucostheta, bins_mucostheta);

      hmap_trkpmom_genie_pm1_bs["total"]["nominal"] = new TH1D("h_trkpmom_total_genie_pm1_nominal", "; Track momentum;", n_bins_pmom, bins_pmom);
      hmap_trkpmom_genie_pm1_bs["signal"]["nominal"] = new TH1D("h_trkpmom_signal_genie_pm1_nominal", "; Track momentum;", n_bins_pmom, bins_pmom);
      hmap_trkpmom_genie_pm1_bs["cosmic"]["nominal"] = new TH1D("h_trkpmom_cosmic_genie_pm1_nominal", "; Track momentum;", n_bins_pmom, bins_pmom);
      hmap_trkpmom_genie_pm1_bs["cosmic_stopmu"]["nominal"] = new TH1D("h_trkpmom_cosmic_stopmu_genie_pm1_nominal", "; Track momentum;", n_bins_pmom, bins_pmom);
      hmap_trkpmom_genie_pm1_bs["cosmic_nostopmu"]["nominal"] = new TH1D("h_trkpmom_cosmic_nostopmu_genie_pm1_nominal", "; Track momentum;", n_bins_pmom, bins_pmom);
      hmap_trkpmom_genie_pm1_bs["outfv"]["nominal"] = new TH1D("h_trkpmom_outfv_genie_pm1_nominal", "; Track momentum;", n_bins_pmom, bins_pmom);
      hmap_trkpmom_genie_pm1_bs["outfv_stopmu"]["nominal"] = new TH1D("h_trkpmom_outfv_stopmu_genie_pm1_nominal", "; Track momentum;", n_bins_pmom, bins_pmom);
      hmap_trkpmom_genie_pm1_bs["outfv_nostopmu"]["nominal"] = new TH1D("h_trkpmom_outfv_nostopmu_genie_pm1_nominal", "; Track momentum;", n_bins_pmom, bins_pmom);
      hmap_trkpmom_genie_pm1_bs["nc"]["nominal"] = new TH1D("h_trkpmom_nc_genie_pm1_nominal", "; Track momentum;", n_bins_pmom, bins_pmom);
      hmap_trkpmom_genie_pm1_bs["nc_proton"]["nominal"] = new TH1D("h_trkpmom_nc_proton_genie_pm1_nominal", "; Track momentum;", n_bins_pmom, bins_pmom);
      hmap_trkpmom_genie_pm1_bs["nc_pion"]["nominal"] = new TH1D("h_trkpmom_nc_pion_genie_pm1_nominal", "; Track momentum;", n_bins_pmom, bins_pmom);
      hmap_trkpmom_genie_pm1_bs["nc_other"]["nominal"] = new TH1D("h_trkpmom_nc_other_genie_pm1_nominal", "; Track momentum;", n_bins_pmom, bins_pmom);
      hmap_trkpmom_genie_pm1_bs["anumu"]["nominal"] = new TH1D("h_trkpmom_anumu_genie_pm1_nominal", "; Track momentum;", n_bins_pmom, bins_pmom);
      hmap_trkpmom_genie_pm1_bs["nue"]["nominal"] = new TH1D("h_trkpmom_nue_genie_pm1_nominal", "; Track momentum;", n_bins_pmom, bins_pmom);
      hmap_trkpmom_genie_pm1_bs["signal_stopmu"]["nominal"] = new TH1D("h_trkpmom_signal_stopmu_genie_pm1_nominal", "; Track momentum;", n_bins_pmom, bins_pmom);
      hmap_trkpmom_genie_pm1_bs["signal_nostopmu"]["nominal"] = new TH1D("h_trkpmom_signal_nostopmu_genie_pm1_nominal", "; Track momentum;", n_bins_pmom, bins_pmom);
      hmap_trkpmom_genie_pm1_bs["cc_other"]["nominal"] = new TH1D("h_trkpmom_ccother_genie_pm1_nominal", "; Track momentum;", n_bins_pmom, bins_pmom);
     
      hmap_trkptheta_genie_pm1_bs["total"]["nominal"] = new TH1D("h_trkptheta_total_genie_pm1_nominal", "; Track Angle;", n_bins_pcostheta, bins_pcostheta);
      hmap_trkptheta_genie_pm1_bs["signal"]["nominal"] = new TH1D("h_trkptheta_signal_genie_pm1_nominal", "; Track Angle;", n_bins_pcostheta, bins_pcostheta);
      hmap_trkptheta_genie_pm1_bs["cosmic"]["nominal"] = new TH1D("h_trkptheta_cosmic_genie_pm1_nominal", "; Track Angle;", n_bins_pcostheta, bins_pcostheta);
      hmap_trkptheta_genie_pm1_bs["cosmic_stopmu"]["nominal"] = new TH1D("h_trkptheta_cosmic_stopmu_genie_pm1_nominal", "; Track Angle;", n_bins_pcostheta, bins_pcostheta);
      hmap_trkptheta_genie_pm1_bs["cosmic_nostopmu"]["nominal"] = new TH1D("h_trkptheta_cosmic_nostopmu_genie_pm1_nominal", "; Track Angle;", n_bins_pcostheta, bins_pcostheta);
      hmap_trkptheta_genie_pm1_bs["outfv"]["nominal"] = new TH1D("h_trkptheta_outfv_genie_pm1_nominal", "; Track Angle;", n_bins_pcostheta, bins_pcostheta);
      hmap_trkptheta_genie_pm1_bs["outfv_stopmu"]["nominal"] = new TH1D("h_trkptheta_outfv_stopmu_genie_pm1_nominal", "; Track Angle;", n_bins_pcostheta, bins_pcostheta);
      hmap_trkptheta_genie_pm1_bs["outfv_nostopmu"]["nominal"] = new TH1D("h_trkptheta_outfv_nostopmu_genie_pm1_nominal", "; Track Angle;", n_bins_pcostheta, bins_pcostheta);
      hmap_trkptheta_genie_pm1_bs["nc"]["nominal"] = new TH1D("h_trkptheta_nc_genie_pm1_nominal", "; Track Angle;", n_bins_pcostheta, bins_pcostheta);
      hmap_trkptheta_genie_pm1_bs["nc_proton"]["nominal"] = new TH1D("h_trkptheta_nc_proton_genie_pm1_nominal", "; Track Angle;", n_bins_pcostheta, bins_pcostheta);
      hmap_trkptheta_genie_pm1_bs["nc_pion"]["nominal"] = new TH1D("h_trkptheta_nc_pion_genie_pm1_nominal", "; Track Angle;", n_bins_pcostheta, bins_pcostheta);
      hmap_trkptheta_genie_pm1_bs["nc_other"]["nominal"] = new TH1D("h_trkptheta_nc_other_genie_pm1_nominal", "; Track Angle;", n_bins_pcostheta, bins_pcostheta);
      hmap_trkptheta_genie_pm1_bs["anumu"]["nominal"] = new TH1D("h_trkptheta_anumu_genie_pm1_nominal", "; Track Angle;", n_bins_pcostheta, bins_pcostheta);
      hmap_trkptheta_genie_pm1_bs["nue"]["nominal"] = new TH1D("h_trkptheta_nue_genie_pm1_nominal", "; Track Angle;", n_bins_pcostheta, bins_pcostheta);
      hmap_trkptheta_genie_pm1_bs["signal_stopmu"]["nominal"] = new TH1D("h_trkptheta_signal_stopmu_genie_pm1_nominal", "; Track Angle;", n_bins_pcostheta, bins_pcostheta);
      hmap_trkptheta_genie_pm1_bs["signal_nostopmu"]["nominal"] = new TH1D("h_trkptheta_signal_nostopmu_genie_pm1_nominal", "; Track Angle;", n_bins_pcostheta, bins_pcostheta);
      hmap_trkptheta_genie_pm1_bs["cc_other"]["nominal"] = new TH1D("h_trkptheta_ccother_genie_pm1_nominal", "; Track Angle;", n_bins_pcostheta, bins_pcostheta);
    
      hmap_thetamup_genie_pm1_bs["total"]["nominal"] = new TH1D("h_thetamup_total_genie_pm1_nominal", "; Track Angle;", n_bins_muptheta, bins_muptheta);
      hmap_thetamup_genie_pm1_bs["signal"]["nominal"] = new TH1D("h_thetamup_signal_genie_pm1_nominal", "; Track Angle;", n_bins_muptheta, bins_muptheta);
      hmap_thetamup_genie_pm1_bs["cosmic"]["nominal"] = new TH1D("h_thetamup_cosmic_genie_pm1_nominal", "; Track Angle;", n_bins_muptheta, bins_muptheta);
      hmap_thetamup_genie_pm1_bs["cosmic_stopmu"]["nominal"] = new TH1D("h_thetamup_cosmic_stopmu_genie_pm1_nominal", "; Track Angle;", n_bins_muptheta, bins_muptheta);
      hmap_thetamup_genie_pm1_bs["cosmic_nostopmu"]["nominal"] = new TH1D("h_thetamup_cosmic_nostopmu_genie_pm1_nominal", "; Track Angle;", n_bins_muptheta, bins_muptheta);
      hmap_thetamup_genie_pm1_bs["outfv"]["nominal"] = new TH1D("h_thetamup_outfv_genie_pm1_nominal", "; Track Angle;", n_bins_muptheta, bins_muptheta);
      hmap_thetamup_genie_pm1_bs["outfv_stopmu"]["nominal"] = new TH1D("h_thetamup_outfv_stopmu_genie_pm1_nominal", "; Track Angle;", n_bins_muptheta, bins_muptheta);
      hmap_thetamup_genie_pm1_bs["outfv_nostopmu"]["nominal"] = new TH1D("h_thetamup_outfv_nostopmu_genie_pm1_nominal", "; Track Angle;", n_bins_muptheta, bins_muptheta);
      hmap_thetamup_genie_pm1_bs["nc"]["nominal"] = new TH1D("h_thetamup_nc_genie_pm1_nominal", "; Track Angle;", n_bins_muptheta, bins_muptheta);
      hmap_thetamup_genie_pm1_bs["nc_proton"]["nominal"] = new TH1D("h_thetamup_nc_proton_genie_pm1_nominal", "; Track Angle;", n_bins_muptheta, bins_muptheta);
      hmap_thetamup_genie_pm1_bs["nc_pion"]["nominal"] = new TH1D("h_thetamup_nc_pion_genie_pm1_nominal", "; Track Angle;", n_bins_muptheta, bins_muptheta);
      hmap_thetamup_genie_pm1_bs["nc_other"]["nominal"] = new TH1D("h_thetamup_nc_other_genie_pm1_nominal", "; Track Angle;", n_bins_muptheta, bins_muptheta);
      hmap_thetamup_genie_pm1_bs["anumu"]["nominal"] = new TH1D("h_thetamup_anumu_genie_pm1_nominal", "; Track Angle;", n_bins_muptheta, bins_muptheta);
      hmap_thetamup_genie_pm1_bs["nue"]["nominal"] = new TH1D("h_thetamup_nue_genie_pm1_nominal", "; Track Angle;", n_bins_muptheta, bins_muptheta);
      hmap_thetamup_genie_pm1_bs["signal_stopmu"]["nominal"] = new TH1D("h_thetamup_signal_stopmu_genie_pm1_nominal", "; Track Angle;", n_bins_muptheta, bins_muptheta);
      hmap_thetamup_genie_pm1_bs["signal_nostopmu"]["nominal"] = new TH1D("h_thetamup_signal_nostopmu_genie_pm1_nominal", "; Track Angle;", n_bins_muptheta, bins_muptheta);
      hmap_thetamup_genie_pm1_bs["cc_other"]["nominal"] = new TH1D("h_thetamup_ccother_genie_pm1_nominal", "; Track Angle;", n_bins_muptheta, bins_muptheta);
 







      //==============================================================================================================
 
      h_eff_muphi_num = new TH1D("h_eff_muphi_num", "h_eff_muphi_num", 15, -3.1415, 3.1415);
      h_eff_muphi_den = new TH1D("h_eff_muphi_den", "h_eff_muphi_den", 15, -3.1415, 3.1415);




    }

  }

#endif
