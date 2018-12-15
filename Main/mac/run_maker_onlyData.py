import sys, os

from ROOT import gROOT
gROOT.ProcessLine(".x " + os.environ['MYSW_DIR'] + "/Utils/rootlogon.C")

from ROOT import Main

maker = Main.Maker()






# BNBON

maker.SetInputFile("/build/kirby/cc1muNp_ubxsec_pid_integration_test_larana/test_ntuples_Dec10/ubxsec_output_data_onbeam.root") # Run1 After Neutrino
maker.SetOutputFile("/uboone/data/users/jiangl/ubxsec_static/v06_26_01_22_Dec/ubxsecana_output_data_onbeam_ubcodev06_26_01_22.root") # Run1 After Neutrino

#maker.SetAnalysisType("cc1unp_analysis");
maker.SetAnalysisType("ccinclusive_analysis");
maker.SetEntries(-1)
maker.SetBeamSpillStart(3.2)    
maker.SetBeamSpillEnd(5.0)    
maker.SetFlashShift(0.)    
maker.SetGainCalibration(243)    
maker.SetCalculatePOT(False)    
maker.SetIsData(True)

maker.PrintConfig()

maker.MakeFile()


# EXTBNB
maker.SetInputFile("/build/kirby/cc1muNp_ubxsec_pid_integration_test_larana/test_ntuples_Dec10/ubxsec_output_data_offbeam.root");
maker.SetOutputFile("/uboone/data/users/jiangl/ubxsec_static/v06_26_01_22_Dec/ubxsecana_output_data_offbeam_ubcodev06_26_01_22.root") # Run1 After Neutrino
 
#maker.SetAnalysisType("cc1unp_analysis");
maker.SetAnalysisType("ccinclusive_analysis");
maker.SetEntries(-1)
maker.SetBeamSpillStart(3.6)    
maker.SetBeamSpillEnd(5.4)    
maker.SetFlashShift(0.406)    
maker.SetGainCalibration(243)    
maker.SetCalculatePOT(False)    
maker.SetIsData(True)

maker.PrintConfig()

maker.MakeFile()


