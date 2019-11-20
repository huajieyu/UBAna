import sys, os

from ROOT import gROOT
gROOT.ProcessLine(".x " + os.environ['MYSW_DIR'] + "/Utils/rootlogon.C")

from ROOT import Main

maker = Main.Maker()


samples = ["aup_bup","adown_bup","aup_bdown","adown_bdown","birks"]

for sample in samples:
  # BNBON
  maker.SetInputFile("/uboone/data/users/kirby/cc1muNp_ubxsec_ntuples/ntuples_Jul1_merge/ubxsec_output_data_onbeam_recomb_"+sample+".root") # Run1 After Neutrino
  maker.SetOutputFile("/uboone/data/users/afurmans/recombinationUncertainty_UBAna_files/ubxsecana_output_data_onbeam_recomb_"+sample+".root")
  
  #maker.SetAnalysisType("ccinclusive_analysis");
  maker.SetAnalysisType("cc1unp_analysis");
  maker.RecoShowerAsTrack(True)
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
  maker.SetInputFile("/uboone/data/users/kirby/cc1muNp_ubxsec_ntuples/ntuples_Jul1_merge/ubxsec_output_data_offbeam_recomb_"+sample+".root") # Run1 After Neutrino
  maker.SetOutputFile("/uboone/data/users/afurmans/recombinationUncertainty_UBAna_files/ubxsecana_output_data_offbeam_recomb_"+sample+".root")
   
  #maker.SetAnalysisType("ccinclusive_analysis");
  maker.SetAnalysisType("cc1unp_analysis");
  maker.RecoShowerAsTrack(True)
  maker.SetEntries(-1)
  maker.SetBeamSpillStart(3.6)    
  maker.SetBeamSpillEnd(5.4)    
  maker.SetFlashShift(0.406)    
  maker.SetGainCalibration(243)    
  maker.SetCalculatePOT(False)    
  maker.SetIsData(True)
  
  maker.PrintConfig()
  
  maker.MakeFile()


