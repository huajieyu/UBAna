import sys

from ROOT import Main

maker = Main.Maker()



# BNBComisc

maker.SetInputFile("../../ubxsecana/Files/ubxsec_output_mc_bnbcosmic_mcc8.7_test6.root")

maker.SetOutputFile("../../Files/Output/ubxsecana_output_bnbcosmic_mcc8.7_test6.root");
  
maker.SetEntries(-1)
maker.SetBeamSpillStart(3.2)    
maker.SetBeamSpillEnd(4.8)    
maker.SetFlashShift(0.)    
maker.SetGainCalibration(198)    
maker.SetCalculatePOT(True)    
maker.SetIsData(False)

maker.PrintConfig()

maker.MakeFile()



# BNBON

maker.SetInputFile("../../ubxsecana/Files/ubxsec_output_data_bnbon_mcc8.7_test6.root")

maker.SetOutputFile("../../Files/Output/ubxsecana_output_bnbon_mcc8.7_test6.root");
  
maker.SetEntries(-1)
maker.SetBeamSpillStart(3.3)    
maker.SetBeamSpillEnd(4.9)    
maker.SetFlashShift(0.)    
maker.SetGainCalibration(243)    
maker.SetCalculatePOT(False)    
maker.SetIsData(True)

maker.PrintConfig()

maker.MakeFile()



# EXTBNB

maker.SetInputFile("../../ubxsecana/Files/ubxsec_output_data_extbnb_mcc8.7_test6.root")

maker.SetOutputFile("../../Files/Output/ubxsecana_output_extbnb_mcc8.7_test6.root");
  
maker.SetEntries(100)
maker.SetBeamSpillStart(3.65)    
maker.SetBeamSpillEnd(5.25)    
maker.SetFlashShift(0.406)    
maker.SetGainCalibration(243)    
maker.SetCalculatePOT(False)    
maker.SetIsData(True)

maker.PrintConfig()

maker.MakeFile()





