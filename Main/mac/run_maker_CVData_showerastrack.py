import sys, os

from ROOT import gROOT
gROOT.ProcessLine(".x " + os.environ['MYSW_DIR'] + "/Utils/rootlogon.C")

from ROOT import Main

maker = Main.Maker()






# BNBON


maker.SetInputFile("/uboone/data/users/kirby/cc1muNp_ubxsec_ntuples/ntuples_Mar20_merge/ubxsec_output_data_onbeam.root");
maker.SetOutputFile("/uboone/data/users/jiangl/ubxsec_static/v06_26_01_22_Apr25/ubxsecana_output_data_onbeam_ubcodev06_26_01_22.root") # Run1 After Neutrino


maker.SetAnalysisType("cc1unp_analysis");
#maker.SetAnalysisType("ccinclusive_analysis");
maker.RecoShowerAsTrack(True);
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
#maker.SetInputFile("/uboone/data/users/kirby/cc1muNp_ubxsec_ntuples/ntuples_Jan27_merge/ubxsec_output_data_extbnb.root");
#maker.SetOutputFile("/uboone/data/users/jiangl/ubxsec_static/v06_26_01_22_Feb/ubxsecana_output_data_offbeam_ubcodev06_26_01_22.root") # Run1 After Neutrino
 
#maker.SetInputFile("/uboone/data/users/kirby/cc1muNp_ubxsec_ntuples/ntuples_Mar14_merge/ubxsec_output_data_offbeam.root");
maker.SetInputFile("/uboone/data/users/kirby/cc1muNp_ubxsec_ntuples/ntuples_Mar20_merge/ubxsec_output_data_offbeam.root");
maker.SetOutputFile("/uboone/data/users/jiangl/ubxsec_static/v06_26_01_22_Apr25/ubxsecana_output_data_offbeam_ubcodev06_26_01_22.root") # Run1 After Neutrino
 
 
maker.SetAnalysisType("cc1unp_analysis");
#maker.SetAnalysisType("ccinclusive_analysis");
maker.RecoShowerAsTrack(True);
maker.SetEntries(-1)
maker.SetBeamSpillStart(3.6)    
maker.SetBeamSpillEnd(5.4)    
maker.SetFlashShift(0.406)    
maker.SetGainCalibration(243)    
maker.SetCalculatePOT(False)    
maker.SetIsData(True)

maker.PrintConfig()

maker.MakeFile()


