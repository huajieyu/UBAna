import sys, os, math

from ROOT import gROOT
gROOT.ProcessLine(".x " + os.environ['MYSW_DIR'] + "/Utils/rootlogon.C")

from ROOT import Main

analyser = Main.AnaEffUnc()

analyser.SetBNBCosmicFile("/uboone/data/users/jiangl/ubxsec_static/v06_26_01_22_Mar/ubxsecana_output_mc_bnbcosmic_ubcodev06_26_01_22_detsyst_CV_genie.root") # tune mic_output_file_name

analyser.SetInputBNBCosmicFile("/uboone/data/users/kirby/cc1muNp_ubxsec_ntuples/ntuples_Mar14_merge/ubxsec_output_mc_bnbcosmic.root") # tune 1


analyser.SetOutputFile("/uboone/data/users/jiangl/ubxsec_static/v06_26_01_22_Mar/ubxsecana_effunc_output_mc_bnbcosmic_ubcodev06_26_01_22.root")


#analyser.SetDirtFile ("/uboone/data/users/jiangl/ubxsec_static/v06_26_01_22_Dec/ubxsecana_output_mc_bnbdirt_ubcodev06_26_01_22.root")
#analyser.SetInTimeCosmicFile  ("/Users/deltutto/RealWork/CCInclusiveEventSelection/Files/Output/ubxsecana_output_bnbcosmic_mcc8.7_test6.root") # Just a placeholder
#analyser.SetBNBONFile("/uboone/data/users/jiangl/ubxsec_static/v06_26_01_22_Dec/ubxsecana_output_data_onbeam_ubcodev06_26_01_22.root")    
#analyser.SetEXTBNBFile("/uboone/data/users/jiangl/ubxsec_static/v06_26_01_22_Dec/ubxsecana_output_data_offbeam_ubcodev06_26_01_22.root")

analyser.SetPrefix("CV_test");



analyser.DoAnaEffUnc();

raw_input("Please press enter to exit.")
