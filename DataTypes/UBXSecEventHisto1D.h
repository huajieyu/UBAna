/**
 * \file UBXSecEventHisto1D.h
 *
 * \ingroup DataTypes
 * 
 * \brief Class def header for a class UBXSecEventHisto1D
 *
 * @author deltutto
 */

/** \addtogroup DataTypes

    @{*/
#ifndef __DATATYPES_UBXSECEVENTHISTO1D_H__
#define __DATATYPES_UBXSECEVENTHISTO1D_H__

#include <iostream>
#include <sstream>
#include <string>
#include <ctime>
#include <vector>
#include <map>
#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <fstream>
#include <math.h>
#include <algorithm>

#include <TChain.h>
#include "THStack.h"
#include "TLegend.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include "TGraph.h"
#include <TString.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TFile.h>
#include "TH2Poly.h"

#include "UBTH2Poly.h"
#include "BootstrapTH1D.h"
#include "BootstrapTH2D.h"
#include "BootstrapTH2DPoly.h"

namespace DataTypes {

  /**
     \class UBXSecEventHisto1D
     User defined class UBXSecEventHisto1D ... these comments are used to generate
     doxygen documentation!
  */
  class UBXSecEventHisto1D{
    
  public:
    
    /// Default constructor
    UBXSecEventHisto1D() {}
    
    /// Default destructor
    ~UBXSecEventHisto1D(){}

    ///
    void InitializeBootstraps();



    double bins_mumom[7] = {0.1, 0.18, 0.30, 0.48, 0.75, 1.14, 2.50};
    double bins_mucostheta[13] = {-1.00, -0.82, -0.66, -0.39, -0.16, 0.05, 0.25, 0.43, 0.59, 0.73, 0.83, 0.91, 1.00};
    int n_bins_mumom = 6;
    int n_bins_mucostheta = 12;

    double bins_pmom[11] = {0.30, 0.41, 0.495, 0.56, 0.62, 0.68, 0.74, 0.80, 0.87, 0.93, 1.20};
    double bins_pcostheta[10]={-1.00, -0.50, 0.00, 0.27, 0.45, 0.62, 0.76, 0.86, 0.94, 1.00};
    int n_bins_pmom = 10;
    int n_bins_pcostheta = 9;
    
    double bins_muptheta[7] = {0.00, 0.8, 1.2, 1.57, 1.94, 2.34, 3.14};
    int n_bins_muptheta = 6;
     
  

    //
    // Total Cross Section
    //

    std::map<std::string,TH1D*> hmap_onebin; ///< Number of events per channel - Total Cross Section
    TH1D * h_eff_onebin_num = 0; ///< Efficiency Numerator - Total Cross Section
    TH1D * h_eff_onebin_den = 0; ///< Efficiency Denominator - Total Cross Section

    BootstrapTH1D * bs_genie_multisim_eff_onebin_num = 0; ///< Efficiency Numerator - Total Cross Section - Genie Multisim
    BootstrapTH1D * bs_genie_multisim_eff_onebin_den = 0; ///< Efficiency Denominator - Total Cross Section - Genie Multisim
    BootstrapTH1D * bs_flux_multisim_eff_onebin_num = 0; ///< Efficiency Numerator - Total Cross Section - Flux Multisim
    BootstrapTH1D * bs_flux_multisim_eff_onebin_den = 0; ///< Efficiency Denominator - Total Cross Section - Flux Multisim
    BootstrapTH1D * bs_extra_syst_multisim_eff_onebin_num = 0; ///< Efficiency Numerator - Total Cross Section - Extra Syst
    BootstrapTH1D * bs_extra_syst_multisim_eff_onebin_den = 0; ///< Efficiency Denominator - Total Cross Section - Extra Syst
    BootstrapTH1D * bs_mc_stat_multisim_eff_onebin_num = 0; ///< Efficiency Numerator - Total Cross Section - MC Stat
    BootstrapTH1D * bs_mc_stat_multisim_eff_onebin_den = 0; ///< Efficiency Denominator - Total Cross Section - MC Stat

    std::map<std::string,std::map<std::string,TH1D*>> hmap_onebin_genie_multisim_bs; ///< Number of events per channel and universe - Total Cross Section - Genie Multisim
    std::map<std::string,std::map<std::string,TH1D*>> hmap_onebin_flux_multisim_bs; ///< Number of events per channel and universe - Total Cross Section - Flux Multisim
    std::map<std::string,std::map<std::string,TH1D*>> hmap_onebin_extra_syst_multisim_bs; ///< Number of events per channel and universe - Total Cross Section - Extra Syst
    std::map<std::string,std::map<std::string,TH1D*>> hmap_onebin_mc_stat_multisim_bs; ///< Number of events per channel and universe - Total Cross Section - MC Stat

    //=================================================================================================================================================
    //
    // Single differential (mumom, pmom, thetamup)
    //
    std::map<std::string,TH1D*> hmap_trkmom;  ///< Number of events per channel - Single Differential (mumom) 
    std::map<std::string,TH1D*> hmap_trkpmom; ///< Number of events per channel - Single Differential (pmom)
    std::map<std::string,TH1D*> hmap_thetamup;
 
    TH1D* h_eff_mumom_num = 0; ///< Efficiency Numerator - Single Differential (mumom)
    TH1D* h_eff_mumom_den = 0; ///< Efficiency Denominator - Single Differential (mumom)

    TH1D* h_eff_pmom_num=0;
    TH1D* h_eff_pmom_den=0;

    TH1D* h_eff_thetamup_num=0;
    TH1D* h_eff_thetamup_den=0;

    TH2D * h_true_reco_mom = 0; ///< Reco to true histogram - Single Differential (mumom)
    TH2D * h_true_reco_pmom = 0;
    TH2D * h_true_reco_thetamup=0;

    BootstrapTH1D * bs_genie_multisim_eff_mumom_num = 0; ///< Efficiency Numerator - Single Differential (mumom) - Genie Multisim
    BootstrapTH1D * bs_genie_multisim_eff_mumom_den = 0; ///< Efficiency Denominator - Single Differential (mumom) - Genie Multisim
    BootstrapTH1D * bs_genie_multisim_eff_pmom_num = 0; ///< Efficiency Numerator - Single Differential (pmom) - Genie Multisim
    BootstrapTH1D * bs_genie_multisim_eff_pmom_den = 0; ///< Efficiency Denominator - Single Differential (pmom) - Genie Multisim
    BootstrapTH1D * bs_genie_multisim_eff_thetamup_num = 0; ///< Efficiency Numerator - Single Differential (thetamup) - Genie Multisim
    BootstrapTH1D * bs_genie_multisim_eff_thetamup_den = 0; ///< Efficiency Denominator - Single Differential (thetamup) - Genie Multisim

    BootstrapTH1D * bs_flux_multisim_eff_mumom_num = 0; ///< Efficiency Numerator - Single Differential (mumom) - Flux Multisim
    BootstrapTH1D * bs_flux_multisim_eff_mumom_den = 0; ///< Efficiency Denominator - Single Differential (mumom) - Flux Multisim
    BootstrapTH1D * bs_flux_multisim_eff_pmom_num = 0; ///< Efficiency Numerator - Single Differential (pmom) - Flux Multisim
    BootstrapTH1D * bs_flux_multisim_eff_pmom_den = 0; ///< Efficiency Denominator - Single Differential (pmom) - Flux Multisim
    BootstrapTH1D * bs_flux_multisim_eff_thetamup_num = 0; ///< Efficiency Numerator - Single Differential (thetamup) - Flux Multisim
    BootstrapTH1D * bs_flux_multisim_eff_thetamup_den = 0; ///< Efficiency Denominator - Single Differential (thetamup) - Flux Multisim

    BootstrapTH1D * bs_extra_syst_multisim_eff_mumom_num = 0; ///< Efficiency Numerator - Single Differential (mumom) - Extra Syst
    BootstrapTH1D * bs_extra_syst_multisim_eff_mumom_den = 0; ///< Efficiency Denominator - Single Differential (mumom) - Extra Syst
    BootstrapTH1D * bs_extra_syst_multisim_eff_pmom_num = 0; ///< Efficiency Numerator - Single Differential (pmom) - Extra Syst
    BootstrapTH1D * bs_extra_syst_multisim_eff_pmom_den = 0; ///< Efficiency Denominator - Single Differential (pmom) - Extra Syst
    BootstrapTH1D * bs_extra_syst_multisim_eff_thetamup_num = 0; ///< Efficiency Numerator - Single Differential (thetamup) - Extra Syst
    BootstrapTH1D * bs_extra_syst_multisim_eff_thetamup_den = 0; ///< Efficiency Denominator - Single Differential (thetamup) - Extra Syst

    BootstrapTH1D * bs_mc_stat_multisim_eff_mumom_num = 0; ///< Efficiency Numerator - Single Differential (mumom) - MC Stat
    BootstrapTH1D * bs_mc_stat_multisim_eff_mumom_den = 0; ///< Efficiency Denominator - Single Differential (mumom) - MC Stat
    BootstrapTH1D * bs_mc_stat_multisim_eff_pmom_num = 0; ///< Efficiency Numerator - Single Differential (mumom) - MC Stat
    BootstrapTH1D * bs_mc_stat_multisim_eff_pmom_den = 0; ///< Efficiency Denominator - Single Differential (mumom) - MC Stat
    BootstrapTH1D * bs_mc_stat_multisim_eff_thetamup_num = 0; ///< Efficiency Numerator - Single Differential (mumom) - MC Stat
    BootstrapTH1D * bs_mc_stat_multisim_eff_thetamup_den = 0; ///< Efficiency Denominator - Single Differential (mumom) - MC Stat
 
    std::map<std::string,std::map<std::string,TH1D*>> hmap_trkmom_genie_multisim_bs; ///< Number of events per channel and universe - Single Differential (mumom) - Genie Multisim
    std::map<std::string,std::map<std::string,TH1D*>> hmap_trkpmom_genie_multisim_bs; ///< Number of events per channel and universe - Single Differential (pmom) - Genie Multisim
    std::map<std::string,std::map<std::string,TH1D*>> hmap_thetamup_genie_multisim_bs; ///< Number of events per channel and universe - Single Differential (thetamup) - Genie Multisim

    std::map<std::string,std::map<std::string,TH1D*>> hmap_trkmom_flux_multisim_bs; ///< Number of events per channel and universe - Single Differential (mumom) - Flux Multisim
    std::map<std::string,std::map<std::string,TH1D*>> hmap_trkpmom_flux_multisim_bs; ///< Number of events per channel and universe - Single Differential (pmom) - Flux Multisim
    std::map<std::string,std::map<std::string,TH1D*>> hmap_thetamup_flux_multisim_bs; ///< Number of events per channel and universe - Single Differential (thetamup) - Flux Multisim

    std::map<std::string,std::map<std::string,TH1D*>> hmap_trkmom_extra_syst_multisim_bs; ///< Number of events per channel and universe - Single Differential (mumom) - Extra Syst
    std::map<std::string,std::map<std::string,TH1D*>> hmap_trkpmom_extra_syst_multisim_bs; ///< Number of events per channel and universe - Single Differential (pmom) - Extra Syst
    std::map<std::string,std::map<std::string,TH1D*>> hmap_thetamup_extra_syst_multisim_bs; ///< Number of events per channel and universe - Single Differential (thetamup) - Extra Syst


    std::map<std::string,std::map<std::string,TH1D*>> hmap_trkmom_mc_stat_multisim_bs; ///< Number of events per channel and universe - Single Differential (mumom) - MC Stat
    std::map<std::string,std::map<std::string,TH1D*>> hmap_trkpmom_mc_stat_multisim_bs; ///< Number of events per channel and universe - Single Differential (pmom) - MC Stat
    std::map<std::string,std::map<std::string,TH1D*>> hmap_thetamup_mc_stat_multisim_bs; ///< Number of events per channel and universe - Single Differential (thetamup) - MC Stat



    BootstrapTH2D * bs_genie_multisim_true_reco_mumom = 0; ///< Reco to true histograms for every universe - Single Differential (mumom) - Genie Multisim
    BootstrapTH2D * bs_flux_multisim_true_reco_mumom = 0; ///< Reco to true histograms for every universe - Single Differential (mumom) - Flux Multisim
    BootstrapTH2D * bs_extra_syst_true_reco_mumom = 0; ///< Reco to true histograms for every universe - Single Differential (mumom) - Extra Syst
    BootstrapTH2D * bs_mc_stat_multisim_true_reco_mumom = 0; ///< Reco to true histograms for every universe - Single Differential (mumom) - MC Stat

    BootstrapTH2D * bs_genie_multisim_true_reco_pmom = 0; ///< Reco to true histograms for every universe - Single Differential (pmom) - Genie Multisim
    BootstrapTH2D * bs_flux_multisim_true_reco_pmom = 0; ///< Reco to true histograms for every universe - Single Differential (pmom) - Flux Multisim
    BootstrapTH2D * bs_extra_syst_true_reco_pmom = 0; ///< Reco to true histograms for every universe - Single Differential (pmom) - Extra Syst
    BootstrapTH2D * bs_mc_stat_multisim_true_reco_pmom = 0; ///< Reco to true histograms for every universe - Single Differential (pmom) - MC Stat

    BootstrapTH2D * bs_genie_multisim_true_reco_thetamup = 0; ///< Reco to true histograms for every universe - Single Differential (thetamup) - Genie Multisim
    BootstrapTH2D * bs_flux_multisim_true_reco_thetamup = 0; ///< Reco to true histograms for every universe - Single Differential (thetamup) - Flux Multisim
    BootstrapTH2D * bs_extra_syst_true_reco_thetamup = 0; ///< Reco to true histograms for every universe - Single Differential (thetamup) - Extra Syst
    BootstrapTH2D * bs_mc_stat_multisim_true_reco_thetamup = 0; ///< Reco to true histograms for every universe - Single Differential (thetamup) - MC Stat



    //
    // Single differential (muangle)
    //

    std::map<std::string,TH1D*> hmap_trktheta; ///< Number of events per channel - Single Differential (mumom)
    std::map<std::string,TH1D*> hmap_trkptheta;

    TH1D* h_eff_muangle_num = 0; ///< Efficiency Numerator - Single Differential (mumom)
    TH1D* h_eff_muangle_den = 0; ///< Efficiency Denominator - Single Differential (mumom)

    TH1D* h_eff_pangle_num=0;
    TH1D* h_eff_pangle_den=0;

    TH2D * h_true_reco_costheta = 0; ///< Reco to true histogram - Single Differential (mumom)
    TH2D * h_true_reco_pcostheta = 0;

    BootstrapTH1D * bs_genie_multisim_eff_muangle_num = 0; ///< Efficiency Numerator - Single Differential (muangle) - Genie Multisim
    BootstrapTH1D * bs_genie_multisim_eff_muangle_den = 0; ///< Efficiency Denominator - Single Differential (muangle) - Genie Multisim
    BootstrapTH1D * bs_flux_multisim_eff_muangle_num = 0; ///< Efficiency Numerator - Single Differential (muangle) - Flux Multisim
    BootstrapTH1D * bs_flux_multisim_eff_muangle_den = 0; ///< Efficiency Denominator - Single Differential (muangle) - Flux Multisim
    BootstrapTH1D * bs_extra_syst_multisim_eff_muangle_num = 0; ///< Efficiency Numerator - Single Differential (muangle) - Extra Syst
    BootstrapTH1D * bs_extra_syst_multisim_eff_muangle_den = 0; ///< Efficiency Denominator - Single Differential (muangle) - Extra Syst
    BootstrapTH1D * bs_mc_stat_multisim_eff_muangle_num = 0; ///< Efficiency Numerator - Single Differential (mumangle - MC Stat
    BootstrapTH1D * bs_mc_stat_multisim_eff_muangle_den = 0; ///< Efficiency Denominator - Single Differential (muangle) - MC Stat

    BootstrapTH1D * bs_genie_multisim_eff_pangle_num = 0; ///< Efficiency Numerator - Single Differential (pangle) - Genie Multisim
    BootstrapTH1D * bs_genie_multisim_eff_pangle_den = 0; ///< Efficiency Denominator - Single Differential (pangle) - Genie Multisim
    BootstrapTH1D * bs_flux_multisim_eff_pangle_num = 0; ///< Efficiency Numerator - Single Differential (pangle) - Flux Multisim
    BootstrapTH1D * bs_flux_multisim_eff_pangle_den = 0; ///< Efficiency Denominator - Single Differential (pangle) - Flux Multisim
    BootstrapTH1D * bs_extra_syst_multisim_eff_pangle_num = 0; ///< Efficiency Numerator - Single Differential (pangle) - Extra Syst
    BootstrapTH1D * bs_extra_syst_multisim_eff_pangle_den = 0; ///< Efficiency Denominator - Single Differential (pangle) - Extra Syst
    BootstrapTH1D * bs_mc_stat_multisim_eff_pangle_num = 0; ///< Efficiency Numerator - Single Differential (mumangle - MC Stat
    BootstrapTH1D * bs_mc_stat_multisim_eff_pangle_den = 0; ///< Efficiency Denominator - Single Differential (pangle) - MC Stat

    std::map<std::string,std::map<std::string,TH1D*>> hmap_trkangle_genie_multisim_bs; ///< Number of events per channel and universe - Single Differential (muangle) - Genie Multisim
    std::map<std::string,std::map<std::string,TH1D*>> hmap_trkangle_flux_multisim_bs; ///< Number of events per channel and universe - Single Differential (muangle) - Flux Multisim
    std::map<std::string,std::map<std::string,TH1D*>> hmap_trkangle_extra_syst_multisim_bs; ///< Number of events per channel and universe - Single Differential (muangle) - Extra Syst
    std::map<std::string,std::map<std::string,TH1D*>> hmap_trkangle_mc_stat_multisim_bs; ///< Number of events per channel and universe - Single Differential (muangle) - MC Stat

    std::map<std::string,std::map<std::string,TH1D*>> hmap_trkpangle_genie_multisim_bs; ///< Number of events per channel and universe - Single Differential (muangle) - Genie Multisim
    std::map<std::string,std::map<std::string,TH1D*>> hmap_trkpangle_flux_multisim_bs; ///< Number of events per channel and universe - Single Differential (muangle) - Flux Multisim
    std::map<std::string,std::map<std::string,TH1D*>> hmap_trkpangle_extra_syst_multisim_bs; ///< Number of events per channel and universe - Single Differential (muangle) - Extra Syst
    std::map<std::string,std::map<std::string,TH1D*>> hmap_trkpangle_mc_stat_multisim_bs; ///< Number of events per channel and universe - Single Differential (muangle) - MC Stat

    BootstrapTH2D * bs_genie_multisim_true_reco_muangle = 0; ///< Reco to true histograms for every universe - Single Differential (muangle) - Genie Multisim
    BootstrapTH2D * bs_flux_multisim_true_reco_muangle = 0; ///< Reco to true histograms for every universe - Single Differential (muangle) - Flux Multisim
    BootstrapTH2D * bs_extra_syst_true_reco_muangle = 0; ///< Reco to true histograms for every universe - Single Differential (muangle) - Extra Syst
    BootstrapTH2D * bs_mc_stat_multisim_true_reco_muangle = 0; ///< Reco to true histograms for every universe - Single Differential (muangle) - MC Stat

    BootstrapTH2D * bs_genie_multisim_true_reco_pangle = 0; ///< Reco to true histograms for every universe - Single Differential (pangle) - Genie Multisim
    BootstrapTH2D * bs_flux_multisim_true_reco_pangle = 0; ///< Reco to true histograms for every universe - Single Differential (pangle) - Flux Multisim
    BootstrapTH2D * bs_extra_syst_true_reco_pangle = 0; ///< Reco to true histograms for every universe - Single Differential (pangle) - Extra Syst
    BootstrapTH2D * bs_mc_stat_multisim_true_reco_pangle = 0; ///< Reco to true histograms for every universe - Single Differential (pangle) - MC Stat
    //==================================================================================================================================================

    /*BootstrapTH2D * bs_genie_pm1_true_reco_mom = 0;
    BootstrapTH2D * bs_genie_pm1_true_reco_pmom = 0;
    BootstrapTH2D * bs_genie_pm1_true_reco_muangle = 0;
    BootstrapTH2D * bs_genie_pm1_true_reco_pangle = 0;
    BootstrapTH2D * bs_genie_pm1_true_reco_thetamup = 0;
    */
    //==================================================================================================================================================== 
    // Number of events histograms - Cross Section Muon Momentum - GENIE pm1sigma
    std::vector<std::string> GENIE_par;
    std::map<std::string,TH2D*> bs_genie_pm1_true_reco_mom;
    std::map<std::string,TH2D*> bs_genie_pm1_true_reco_pmom;
    std::map<std::string,TH2D*> bs_genie_pm1_true_reco_muangle;
    std::map<std::string,TH2D*> bs_genie_pm1_true_reco_pangle;
    std::map<std::string,TH2D*> bs_genie_pm1_true_reco_thetamup;




    std::map<std::string,std::map<std::string,TH1D*>> hmap_trkmom_genie_pm1_bs;
    std::map<std::string,std::map<std::string,TH1D*>> hmap_trktheta_genie_pm1_bs;
    std::map<std::string,std::map<std::string,TH1D*>> hmap_trkpmom_genie_pm1_bs;
    std::map<std::string,std::map<std::string,TH1D*>> hmap_trkptheta_genie_pm1_bs;
    std::map<std::string,std::map<std::string,TH1D*>> hmap_thetamup_genie_pm1_bs;



 

    TH1D* h_eff_muphi_num = 0;
    TH1D* h_eff_muphi_den = 0;
     

  protected:


    
  };
}

#endif
/** @} */ // end of doxygen group 

