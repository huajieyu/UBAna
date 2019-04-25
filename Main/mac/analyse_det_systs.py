import sys, os, math

from ROOT import gROOT
gROOT.ProcessLine(".x " + os.environ['MYSW_DIR'] + "/Utils/rootlogon.C")

from ROOT import Main

analyser = Main.Analyse()


analyser.SetBNBCosmicFile     ("/uboone/data/users/jiangl/ubxsec_static/v06_26_01_22_Feb/ubxsecana_output_mc_bnbcosmic_ubcodev06_26_01_22_detsyst_CV.root");
#analyser.SetInTimeCosmicFile  ("/Users/deltutto/CCInclusiveFiles/Output/ubxsecana_output_bnbcosmic_mcc8.7_test6.root") # Just a placeholder
analyser.SetDirtFile          ("/uboone/data/users/jiangl/ubxsec_static/v06_26_01_22_Feb/ubxsecana_output_mc_bnbdirt_ubcodev06_26_01_22.root")


analyser.SetBNBONFile("/uboone/data/users/jiangl/ubxsec_static/v06_26_01_22_Feb/ubxsecana_output_data_onbeam_ubcodev06_26_01_22.root")    
analyser.SetEXTBNBFile("/uboone/data/users/jiangl/ubxsec_static/v06_26_01_22_Feb/ubxsecana_output_data_offbeam_ubcodev06_26_01_22.root")
analyser.SetBNBPOT(1.592e+20)    
analyser.SetBNBONTriggers(35388924.0)   
analyser.SetEXTBNBTriggers(72299264.0)
analyser.SetPrefix("cv")


analyser.SetFluxCorrectionWeight(1.028)

extra_unc = math.sqrt(0.02*0.02 + 0.06*0.06) # POT counting, beam window
# analyser.SetExtraUncertainty(extra_unc)

analyser.SetTargetFluxSystematic("total")
analyser.DoFluxSystematics(False)
analyser.ImportFluxSystematics(False)
analyser.DoGenieSystematics(False)
analyser.DoExtraSystematics(False)
analyser.ImportExtraSystematics(False)
analyser.ImportFluxSystematics(False)
analyser.SetExtraFluxUncertainty(0.)


det_syst_list = ["CV", "DLdown", "DLup", "DTup", "DTdown", "noiseAmpDown", "noiseAmpDown", "squeezeResp", "squeezeResp", "downPEnoise", "upPEnoise","LArG4BugFix","altDeadChannels", "dataSCE", "deadSaturatedChannels", "enhancedexttpcvis","lifetime10ms", "withDIC"]


for systname in det_syst_list:
  file_name = "/uboone/data/users/jiangl/ubxsec_static/v06_26_01_22_Feb/ubxsecana_output_mc_bnbcosmic_ubcodev06_26_01_22_detsyst_" + systname + ".root"
  print file_name
  analyser.SetBNBCosmicFile(file_name)
  analyser.SetPrefix(systname)
  analyser.DoAnalise()




raw_input("Please press enter to exit.")
