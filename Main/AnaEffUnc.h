/**
 * \file AnaEffUnc.h
 *
 * \ingroup Main
 * 
 * \brief Class def header for a class Analyse
 *
 * @author Libo Jiang
 */

/** \addtogroup Main

    @{*/
#ifndef __MAIN_ANAEFFUNC_H__
#define __MAIN_ANAEFFUNC_H__

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

#include <TRandom1.h>
#include <TSystem.h>
#include <TBrowser.h>
#include <TApplication.h>
#include <TChain.h>
#include "TThread.h"
#include "THStack.h"
#include "TLegend.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include "TGraph.h"
#include <TString.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TH2Poly.h>

#include "ubana/DataTypes/UBTH2Poly.h"
#include "ubana/DataTypes/BootstrapTH2DPoly.h"

#include "ubana/DataTypes/UBXSecEventHisto1D.h"
#include "ubana/DataTypes/UBXSecEventHisto.h"

#include "UBXSecEvent.h"
#include "ubana/DataTypes/BootstrapTH1D.h"
#include "ubana/DataTypes/BootstrapTH2D.h"
#include "ubana/Base/PlottingTools.h"
#include "ubana/Base/CrossSectionCalculator1D.h"
#include "ubana/Base/MigrationMatrix2D.h"
#include "ubana/Base/MigrationMatrix4D.h"
#include "ubana/Base/MigrationMatrix4DPoly.h"
#include "ubana/Base/CrossSectionCalculator2D.h"
#include "ubana/Base/ReweightingPlotter.h"
#include "ubana/Base/CovarianceCalculator2D.h"
#include "ubana/Base/CrossSectionBootstrapCalculator1D.h"
#include "ubana/Base/CrossSectionBootstrapCalculator2D.h"
#include "ubana/Base/CrossSectionBootstrapCalculator2DPoly.h"
#include "ubana/Base/CrossSectionCalculator2DPoly.h"
#include "ubana/Base/UncertaintyPlotter.h"

#include "ubana/Base/LoggerFeature.h"


using namespace DataTypes;
using namespace Base;


namespace Main {

  class AnaEffUnc : public LoggerFeature {

  public:
    
   /// Default constructor
   AnaEffUnc(std::string name = "AnaEffUnc") 
   : LoggerFeature(name) {}
    
   /// Default destructor
   ~AnaEffUnc(){}
   void SetBNBCosmicFile(std::string f);
   void SetOutputFile(std::string f);
   void SetInputBNBCosmicFile(std::string f);
    
   void SetPrefix(std::string);
    
   void DoAnaEffUnc();

   void GetRMS(std::map<std::string, TH1D*>* bnbcosmic);  
   private:
   
    std::string mc_bnbcosmic_file_name     = "ubxsecana_output.root";
    std::string mc_bnbcosmic_output_file_name   = "";
    std::string mc_bnbcosmic_input_file_name = "";


    


    bool _do_pm1sigma_plots = true;
    std::string _prefix = "";



  };

}
#endif















