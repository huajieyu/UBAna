import sys, os, math

from ROOT import gROOT
gROOT.ProcessLine(".x " + os.environ['MYSW_DIR'] + "/Utils/rootlogon.C")

from ROOT import Main

analyser = Main.Analyse()


analyser.SetInTimeCosmicFile  ("/Users/deltutto/RealWork/CCInclusiveEventSelection/Files/Output/ubxsecana_output_bnbcosmic_mcc8.7_test6.root") # Just a placeholder
analyser.SetBNBONFile         ("/Users/deltutto/RealWork/CCInclusiveEventSelection/Files/Output/ubxsecana_output_bnbon_a_mcc8.9_test4.root")    
analyser.SetEXTBNBFile        ("/Users/deltutto/RealWork/CCInclusiveEventSelection/Files/Output/ubxsecana_output_extbnb_a_mcc8.9_test4.root")
analyser.SetBNBPOT(1.627e+20)   
analyser.SetBNBONTriggers(36177265)    
analyser.SetEXTBNBTriggers(33320382)
analyser.DoFluxSystematics(True)
analyser.ImportFluxSystematics(False)
extra_unc = math.sqrt(0.02*0.02 + 0.06*0.06 + 0.0699*0.0699) # POT counting, beam window, cosmics (overlay)
analyser.SetExtraFluxUncertainty(extra_unc)
analyser.DoGenieSystematics(False)
analyser.SetFluxCorrectionWeight(1.028)

extra_unc = math.sqrt(0.02*0.02 + 0.06*0.06) # POT counting, beam window, cosmics (overlay)
analyser.SetExtraUncertainty(extra_unc)



flux_syst_list = ["FluxUnisim", "kminus", "kplus", "kzero", "piminus", "piplus"]
#flux_syst_list = ["piplus"]

for systname in flux_syst_list:
  file_name = "/Users/deltutto/RealWork/CCInclusiveEventSelection/Files/Output/ubxsecana_output_bnbcosmic_mcc8.7_test6_" + systname + ".root"
  analyser.SetBNBCosmicFile(file_name)
  analyser.SetPrefix(systname);
  analyser.SetTargetFluxSystematic(systname)
  analyser.DoAnalise()


raw_input("Please press enter to exit.")
