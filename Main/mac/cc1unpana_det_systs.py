import sys, os, math

from ROOT import gROOT
gROOT.ProcessLine(".x " + os.environ['MYSW_DIR'] + "/Utils/rootlogon.C")

from ROOT import Main
#gROOT.SetBatch(kTRUE)

analyser = Main.CC1uNPAna()

#analyser.SetBNBCosmicFile     ("/uboone/data/users/jiangl/ubxsec_static/v06_26_01_22_Feb/ubxsecana_output_mc_bnbcosmic_ubcodev06_26_01_22_detsyst_CV.root");
analyser.SetBNBCosmicFile("/uboone/data/users/kirby/ubxsec_static/v06_26_01_22_Apr/detsyst/ubxsecana_output_mc_bnbcosmic_ubcodev06_26_01_22.root")

#analyser.SetInTimeCosmicFile  ("/Users/deltutto/CCInclusiveFiles/Output/ubxsecana_output_bnbcosmic_mcc8.7_test6.root") # Just a placeholder
analyser.SetDirtFile          ("/uboone/data/users/jiangl/ubxsec_static/v06_26_01_22_Apr19/ubxsecana_output_mc_bnbdirt_ubcodev06_26_01_22.root")

analyser.SetBNBONFile("/uboone/data/users/kirby/ubxsec_static/v06_26_01_22_Apr/ubxsecana_output_data_onbeam_ubcodev06_26_01_22.root")
analyser.SetEXTBNBFile("/uboone/data/users/kirby/ubxsec_static/v06_26_01_22_Apr/ubxsecana_output_data_offbeam_ubcodev06_26_01_22.root")
analyser.SetBNBPOT(1.582e+20)    
analyser.SetBNBONTriggers(31865413.0)   
analyser.SetEXTBNBTriggers(71687142.0)
analyser.SetPrefix("CV")


analyser.SetFluxCorrectionWeight(1.028)

extra_unc = math.sqrt(0.02*0.02 + 0.06*0.06) # POT counting, beam window
# analyser.SetExtraUncertainty(extra_unc)

#analyser.SetTargetFluxSystematic("total")
#analyser.DoFluxSystematics(True)
#analyser.ImportFluxSystematics(True)
#analyser.DoGenieSystematics(True)
#analyser.DoExtraSystematics(True)
#analyser.ImportExtraSystematics(True)
#analyser.ImportFluxSystematics(True)
#analyser.SetExtraFluxUncertainty(0.)


#det_syst_list = ["CV", "DLdown", "DLup", "DTup", "DTdown", "noiseAmpUp", "noiseAmpDown", "squeezeResp", "stretchResp", "downPEnoise", "upPEnoise","LArG4BugFix","altDeadChannels", "dataSCE", "deadSaturatedChannels", "enhancedexttpcvis","lifetime10ms", "withDIC"] ###don't use this one because it has two extras that haven't been processed successfully yet
#det_syst_list = ["CV", "DLdown", "DLup", "DTdown", "DTup", "LArG4BugFix", "altDeadChannels", "dataSCE", "downPEnoise", "enhancedexttpcvis","lifetime10ms", "noiseAmpDown", "noiseAmpUp", "squeezeResp", "stretchResp", "upPEnoise", "withDIC"]
#det_syst_list = ["CV", "stretchResp"]
det_syst_list = ["enhancedexttpcvis","lifetime10ms", "withDIC"]

#det_syst_list = ["upPEnoise","LArG4BugFix","altDeadChannels", "dataSCE", "deadSaturatedChannels", "enhancedexttpcvis","lifetime10ms", "withDIC"]
#det_syst_list = ["deadSaturatedChannels"]

for systname in det_syst_list:
  file_name = "/uboone/data/users/kirby/ubxsec_static/v06_26_01_22_Apr/detsyst/ubxsecana_output_mc_bnbcosmic_detsyst_" + systname + ".root"
  print file_name
  analyser.SetBNBCosmicFile(file_name)
  analyser.SetPrefix(systname)
  analyser.DoAnalise()

raw_input("Please press enter to exit.")
