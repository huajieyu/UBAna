import sys, os, math

from ROOT import gROOT
gROOT.ProcessLine(".x " + os.environ['MYSW_DIR'] + "/Utils/rootlogon.C")

from ROOT import Main

analyser = Main.CC1uNPAna()

analyser.SetBNBCosmicFile("/uboone/data/users/jiangl/ubxsec_static/v06_26_01_22_Apr25/ubxsecana_output_mc_bnbcosmic_ubcodev06_26_01_22_Tune1_tot.root") # tune 1
analyser.SetDirtFile ("/uboone/data/users/jiangl/ubxsec_static/v06_26_01_22_Apr25/ubxsecana_output_mc_bnbdirt_ubcodev06_26_01_22.root")
#analyser.SetInTimeCosmicFile  ("/Users/deltutto/RealWork/CCInclusiveEventSelection/Files/Output/ubxsecana_output_bnbcosmic_mcc8.7_test6.root") # Just a placeholder

analyser.SetBNBONFile("/uboone/data/users/jiangl/ubxsec_static/v06_26_01_22_Apr25/ubxsecana_output_data_onbeam_ubcodev06_26_01_22.root")    
analyser.SetEXTBNBFile("/uboone/data/users/jiangl/ubxsec_static/v06_26_01_22_Apr25/ubxsecana_output_data_offbeam_ubcodev06_26_01_22.root")



analyser.SetBNBPOT(1.584e20);    
analyser.SetBNBONTriggers(35222964.0)
analyser.SetEXTBNBTriggers(71687142.0)
analyser.SetPrefix("CV");
analyser.SetFluxCorrectionWeight(1.028)

#extra_unc = math.sqrt(0.02*0.02 + 0.06*0.06) # POT counting, beam window, cosmics (overlay)
extra_unc = 0
analyser.SetExtraUncertainty(extra_unc)

analyser.DoGenieSystematics(False)
analyser.ImportGenieSystematics(False)

analyser.DoFluxSystematics(False)
analyser.ImportFluxSystematics(False)

analyser.DoExtraSystematics(False);
analyser.ImportExtraSystematics(False);

analyser.SetExtraFluxUncertainty(0.)
#analyser.SetTargetFluxSystematic("FluxUnisim");
#analyser.SetTargetFluxSystematic("kminus_PrimaryHadronNormalizat");
#analyser.SetTargetFluxSystematic("kplus_PrimaryHadronFeynmanScal");
#analyser.SetTargetFluxSystematic("kzero_PrimaryHadronSanfordWang");
#analyser.SetTargetFluxSystematic("piminus_PrimaryHadronSWCentral");
#analyser.SetTargetFluxSystematic("piplus_PrimaryHadronSWCentralS");
#analyser.SetTargetFluxSystematic("total");

analyser.DoAnalise();

raw_input("Please press enter to exit.")
