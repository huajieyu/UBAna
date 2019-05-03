import sys, os

from ROOT import gROOT
gROOT.ProcessLine(".x " + os.environ['MYSW_DIR'] + "/Utils/rootlogon.C")

from ROOT import Main

maker = Main.Maker()




det_syst_list = ["CV", "DLdown", "DLup", "DTdown", "DTup", "LArG4BugFix", "altDeadChannels", "dataSCE", "downPEnoise", "enhancedexttpcvis","lifetime10ms", "noiseAmpDown", "noiseAmpUp", "squeezeResp", "stretchResp", "upPEnoise", "withDIC"]
#det_syst_list = ["CV", "DLdown"]
#det_syst_list = ["CV", "dataSCE", "withDIC", "squeezeResp", "DLdown", "DLup", "DTdown", "DTup", "LArG4BugFix", "downPEnoise", "upPEnoise", "noiseAmpDown", "noiseAmpUp"]



for systname in det_syst_list:

  
  maker.SetInputFile("/uboone/data/users/kirby/cc1muNp_ubxsec_ntuples/ntuples_Apr05_merge/ubxsec_output_mc_bnbcosmic_detsyst_" + systname + ".root")
  maker.SetOutputFile("/uboone/data/users/kirby/ubxsec_static/v06_26_01_22_Apr/detsyst/ubxsecana_output_mc_bnbcosmic_detsyst_" + systname + ".root");

  maker.SetEntries(-1)
  maker.SetInitialEntry(0)
  maker.SetBeamSpillStart(3.1)    
  maker.SetBeamSpillEnd(4.9)    
  maker.SetFlashShift(0.)    
  maker.SetGainCalibration(198)    
  maker.SetCalculatePOT(True)    
  maker.SetIsData(False)
  maker.SetExtraWeight(1.028); # Flux correction
  maker.RecoShowerAsTrack(True);
  maker.SetAnalysisType("cc1unp_analysis");
#maker.SetAnalysisType("ccinclusive_analysis");
#  maker.SetTargetFluxSystematic("total"); # not needed here

  maker.PrintConfig()

  maker.MakeFile()
