import sys, os, math

from ROOT import gROOT
gROOT.ProcessLine(".x " + os.environ['MYSW_DIR'] + "/Utils/rootlogon.C")

from ROOT import Main

analyser = Main.Analyse()


# analyser.SetBNBCosmicFile     ("/Users/deltutto/CCInclusiveFiles/Output/ubxsecana_output_mc_bnbcosmic_ubcodev06_26_01_22__v2_nosyst.root") # Tune 1
# analyser.SetBNBCosmicFile     ("/Users/deltutto/CCInclusiveFiles/Output/ubxsecana_output_mc_bnbcosmic_ubcodev06_26_01_22__v4_all.root") # Tune 1 - full stat - full syst
analyser.SetBNBCosmicFile     ("/Users/deltutto/CCInclusiveFiles/Output/ubxsecana_output_mc_bnbcosmic_ubcodev06_26_01_22__v4_new.root") # Tune 1 - full stat
# analyser.SetBNBCosmicFile     ("/Users/deltutto/CCInclusiveFiles/Output/ubxsecana_output_mc_bnbcosmic_ubcodev06_26_01_22__v4_nodqdxcut.root") # Tune 1 - full stat
# analyser.SetBNBCosmicFile     ("/Users/deltutto/CCInclusiveFiles/Output/ubxsecana_output_mc_bnbcosmic_ubcodev06_26_01_22__v5.root") # Tune 1 - dev - full syst
# analyser.SetBNBCosmicFile     ("/Users/deltutto/CCInclusiveFiles/Output/ubxsecana_output_mc_bnbcosmic_tune3_full_ubcodev06_26_01_22__v3.root") # Tune 3 (full)

# analyser.SetBNBCosmicFile     ("/Users/deltutto/CCInclusiveFiles/Output/DetSyst/ubxsecana_output_mc_bnbcosmic_detsyst_" + "withDIC" + "_ubcodev06_26_01_22_nodqdxcut.root") # withDIC det syst sample w/o dqdx cuts

analyser.SetDirtFile          ("/Users/deltutto/CCInclusiveFiles/Output/ubxsecana_output_mc_dirt_ubcodev06_26_01_22__v3.root")

analyser.SetInTimeCosmicFile  ("/Users/deltutto/CCInclusiveFiles/Output/ubxsecana_output_bnbcosmic_mcc8.7_test6.root") # Just a placeholder

analyser.SetBNBONFile         ("/Users/deltutto/CCInclusiveFiles/Output/ubxsecana_output_data_bnbon_ubcodev06_26_01_22__v4.root")    
analyser.SetEXTBNBFile        ("/Users/deltutto/CCInclusiveFiles/Output/ubxsecana_output_data_extbnb_ubcodev06_26_01_22__v4.root")
analyser.SetBNBPOT(3.446e+19)    
analyser.SetBNBONTriggers(8098837.0)   
analyser.SetEXTBNBTriggers(58018807.0)
analyser.SetPrefix("cv")
analyser.SetFluxCorrectionWeight(1.028)

# analyser.ImportAlternativeMC("xsec_file_cv_tune3.root")

analyser.SetBeamOffSubtraction(False)
analyser.SetBreakdownPlots(True)

extra_unc = math.sqrt(0.02*0.02) # POT counting
# analyser.SetExtraUncertainty(extra_unc)

analyser.ImportDetectorSystematics(False)

analyser.ImportCosmicSystematics(False)

analyser.ImportDirtSystematics(False)

analyser.DoGenieSystematics(False)
analyser.ImportGenieSystematics(False)

analyser.DoExtraSystematics(False)
analyser.ImportExtraSystematics(False)

analyser.DoMCStatSystematics(False)
analyser.ImportMCStatSystematics(False)


analyser.DoFluxSystematics(False)
analyser.ImportFluxSystematics(False)
analyser.SetExtraFluxUncertainty(0.)
analyser.SetTargetFluxSystematic("total"); # Other options: "FluxUnisim", "kminus", "kplus", "kzero", "piminus", "piplus"

analyser.DoAnalise();

raw_input("Please press enter to exit.")