import sys, os

from ROOT import gROOT
gROOT.ProcessLine(".x " + os.environ['MYSW_DIR'] + "/Utils/rootlogon.C")

from ROOT import Main

maker = Main.Maker()






# BNBON

maker.SetInputFile("/build/kirby/cc1muNp_ubxsec_pid_integration_test_larana/UBXSec_ntuples_16Nov2018_v1/ubxsec_output_data_bnbon.root") # Run1 After Neutrino
# maker.SetOutputFile("/Users/deltutto/CCInclusiveFiles/Output/ubxsecana_output_data_bnbon_run1_ubcodev06_26_01_22__v1.root") # Run1 After Neutrino
#maker.SetInputFile("/Users/deltutto/CCInclusiveFiles/Input/ubxsec_output_data_bnbon_run1_ubcodev06_26_01_22__v4.root") # Run1 After Neutrino
#maker.SetOutputFile("/Users/deltutto/CCInclusiveFiles/Output/ubxsecana_output_data_bnbon_run1_ubcodev06_26_01_22__v4.root") # Run1 After Neutrino

maker.SetEntries(-1)
maker.SetBeamSpillStart(3.2)    
maker.SetBeamSpillEnd(5.0)    
maker.SetFlashShift(0.)    
maker.SetGainCalibration(243)    
maker.SetCalculatePOT(False)    
maker.SetIsData(True)

maker.PrintConfig()

maker.MakeFile()



