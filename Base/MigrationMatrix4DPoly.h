/**
 * \file MigrationMatrix4DPoly.h
 *
 * \ingroup Base
 * 
 * \brief Class def header for a class MigrationMatrix4DPoly
 *
 * @author deltutto
 */

/** \addtogroup Base

    @{*/
#ifndef __BASE_MIGRATIONMATRIX4DPOLY_H__
#define __BASE_MIGRATIONMATRIX4DPOLY_H__

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <ctime>
#include <vector>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <iomanip>
#include <map>
#include <time.h>

#include <TSystem.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TRandom1.h>
#include <TSystem.h>
#include <TBrowser.h>
#include <TApplication.h>
#include <TChain.h>
#include "TThread.h"
#include "THStack.h"
#include "TLegend.h"
#include "TEfficiency.h"
#include "TGraph.h"
#include <TString.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <TMath.h>
#include "Math/SMatrix.h"
#include "TMatrix.h"
#include "TLine.h"
#include "TGaxis.h"

#include "Types.h"
#include "ubana/DataTypes/UBTH2Poly.h"
#include "LoggerFeature.h"


using namespace DataTypes;


namespace Base {

  /**
     \class MigrationMatrix4DPoly
     User defined class MigrationMatrix4DPoly ... these comments are used to generate
     doxygen documentation!
  */
  class MigrationMatrix4DPoly : public LoggerFeature {
    
  public:
    
    /// Default constructor
    MigrationMatrix4DPoly(std::string name = "MigrationMatrix4DPoly") 
    : LoggerFeature(name) {}
    
    /// Default destructor
    ~MigrationMatrix4DPoly(){}

    ///
    TMatrix CalculateMigrationMatrix(); 

    /// Sets the TTree for all signal selected events, to construct S (no longer used)
    void SetTTree(TTree*);

    /// Set 2D reco histogram for every true bin m
    void SetRecoPerTrueHistos(std::vector<UBTH2Poly*>);

    /// Set 1D unrolled vector for every true bin m (can be used instead of SetRecoPerTrueHistos)
    void SetRecoPerTrueVectors(std::vector<std::vector<double>>);

    ///
    void SetBins(int);

    /// Sets a template histo, just to keep track of bin sizes etc.
    void SetTemplateHisto(UBTH2Poly *h) { _th2poly_template = (UBTH2Poly *) h->Clone("_th2poly_template"); };

    /// Openes a text outpu files (then you can call PrintSmearingMatrixLatex)
    void SetLaTeXOutputFileName(std::string name);

    /// Sets the output directory
    void SetOutDir(std::string dir = "migration_matrix_poly_4D_plots");

    /// Prints the smeating matrix in latex format on a file opend by SetOutputFileName()
    void PrintSmearingMatrixLatex();

    ///
    void PlotMatrix();

    ///
    void DoMakePlots(bool option = true) { _do_make_plots = option; }

    /// If called uses the weights with name specified
    void UseWeights(std::string weight_name = "universe0", std::string weight_type = "genie_multisim");


  private:

    /// Checks the number of entries for a true bin.
    void CheckEntries(UBTH2Poly* h, int bin_number);

    std::string _prefix = "[MigrationMatrix4DPoly] ";

    TTree *_tree;
    std::vector<UBTH2Poly*> _h_reco_per_true; ///< 2D reco histogram for every true bin
    std::vector<std::vector<double>> _v_reco_per_true; ///< 1D unrolled vector for every true bin

    UBTH2Poly *_th2poly_template; ///< A template just to keep track of bin sizes etc.

    // std::vector<std::pair<double, double>> _var1_bins;
    // std::vector<std::pair<double, double>> _var2_bins;
    int _n_bins;

    UBTH2Poly *_reco_per_true = 0; ///< A UBTH2Poly, just used to retrive the bins

    TMatrix _S; ///< The smearing matrix

    std::ofstream _f_out; ///< The output file

    bool _do_make_plots = false;
    
    std::string _outdir;
    std::string _folder = "MigrationMatrix4DPolyPlots/";

    bool _use_weights = false; ///< If true uses additional wights (usually for multisim)
    std::string _weight_name = "not_set";
    std::string _weight_type = "not_set";
    
  };
}

#endif
/** @} */ // end of doxygen group 

