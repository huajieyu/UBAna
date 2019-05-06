import sys, os, math

from ROOT import gROOT
gROOT.ProcessLine(".x " + os.environ['MYSW_DIR'] + "/Utils/rootlogon.C")

from ROOT import Main

analyser = Main.CC1uNPAna()

#analyser.SetBNBCosmicFile("/uboone/data/users/jiangl/ubxsec_static/v06_26_01_22_Mar/ubxsecana_output_mc_bnbcosmic_ubcodev06_26_01_22_detsyst_CV.root") # tune 1
analyser.SetDirtFile ("/uboone/data/users/jiangl/ubxsec_static/v06_26_01_22_Apr16/ubxsecana_output_mc_bnbdirt_ubcodev06_26_01_22.root")
#analyser.SetInTimeCosmicFile  ("/Users/deltutto/RealWork/CCInclusiveEventSelection/Files/Output/ubxsecana_output_bnbcosmic_mcc8.7_test6.root") # Just a placeholder
analyser.SetBNBONFile("/uboone/data/users/jiangl/ubxsec_static/v06_26_01_22_Apr16/ubxsecana_output_data_onbeam_ubcodev06_26_01_22.root")    
analyser.SetEXTBNBFile("/uboone/data/users/jiangl/ubxsec_static/v06_26_01_22_Apr16/ubxsecana_output_data_offbeam_ubcodev06_26_01_22.root")




analyser.SetBNBPOT(1.584e20);    
analyser.SetBNBONTriggers(35222964.0)
analyser.SetEXTBNBTriggers(71687142.0)
analyser.DoFluxSystematics(True)
analyser.ImportFluxSystematics(True)
#extra_unc = math.sqrt(0.02*0.02 + 0.06*0.06 + 0.0699*0.0699) # POT counting, beam window, cosmics (overlay)
extra_unc = 0
analyser.SetExtraFluxUncertainty(extra_unc)
analyser.DoGenieSystematics(False)
analyser.SetFluxCorrectionWeight(1.028)

#extra_unc = math.sqrt(0.02*0.02 + 0.06*0.06) # POT counting, beam window, cosmics (overlay)
extra_unc = 0
analyser.SetExtraUncertainty(extra_unc)



#flux_syst_list = ["FluxUnisim", "kminus", "kplus", "kzero", "piminus", "piplus"]
#flux_syst_list = ["FluxUnisim"]
#flux_syst_list = ["kminus"]
#flux_syst_list = ["kplus"]
#flux_syst_list = ["kzero"]
#flux_syst_list = ["piminus"]
flux_syst_list = ["piplus"]






for systname in flux_syst_list:
  file_name = "/uboone/data/users/jiangl/ubxsec_static/v06_26_01_22_Apr16/ubxsecana_output_mc_bnbcosmic_ubcodev06_26_01_22_detsyst_flux_" + systname + ".root"
  analyser.SetBNBCosmicFile(file_name)
  analyser.SetPrefix(systname);
  analyser.SetTargetFluxSystematic(systname)
  analyser.DoAnalise()


raw_input("Please press enter to exit.")
