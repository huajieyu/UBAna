import sys, os

from ROOT import gROOT
gROOT.ProcessLine(".x " + os.environ['MYSW_DIR'] + "/Utils/rootlogon.C")

from ROOT import Main

maker = Main.Maker()



# BNBComisc

maker.SetInputFile("/Users/deltutto/RealWork/CCInclusiveEventSelection/Files/Input/ubxsec_output_mc_bnbcosmic_mcc8.9_500k.root") # tune 1
maker.SetOutputFile("/Users/deltutto/RealWork/CCInclusiveEventSelection/Files/Output/ubxsecana_output_bnbcosmic_mcc8.9_500k.root"); # tune 1
# maker.SetInputFile("/Users/deltutto/RealWork/CCInclusiveEventSelection/Files/Input/ubxsec_output_mc_bnbcosmic_tune3_mcc8.9_test5.root") # tune 3
# maker.SetOutputFile("/Users/deltutto/RealWork/CCInclusiveEventSelection/Files/Output/ubxsecana_output_mc_bnbcosmic_tune3_mcc8.9_test5.root"); # tune 3
# maker.SetInputFile("/Users/deltutto/RealWork/CCInclusiveEventSelection/Files/Input/ubxsec_output_mc_bnbcosmic_mcc8.7_test6.root"); # tune 1 - genie flux
# maker.SetOutputFile("/Users/deltutto/RealWork/CCInclusiveEventSelection/Files/Output/ubxsecana_output_bnbcosmic_mcc8.7_test6.root"); # tune 1 - genie flux
# maker.SetInputFile("/uboone/data/users/mdeltutt/ubxsec_static/mcc8.9/ubxsec_output_mc_bnbcosmic_mcc8.9_test6.root") # tune 1 high stat weights
# maker.SetOutputFile("/uboone/data/users/mdeltutt/ubxsec_static/mcc8.9/ubxsecana_output_mc_bnbcosmic_mcc8.9_test6_maup_mecoff.root"); # tune 1 high stat weights
#maker.SetInputFile("/uboone/data/users/mdeltutt/ubxsec_static/mcc8.9/ubxsec_output_mc_bnbcosmic_mcc8.9_test6_tune3.root"); # tune 3 with flux syst
#maker.SetOutputFile("/uboone/data/users/mdeltutt/ubxsec_static/mcc8.9/ubxsecana_output_mc_bnbcosmic_mcc8.9_test6_tune3.root"); # tune 3 with flux syst
#maker.SetInputFile("/uboone/data/users/wketchum/det_var_files_Apr19/ubxsec_output_mc_bnbcosmic.root") # wes test
# maker.SetOutputFile("./ubxsecanatemp_output_mc_bnbcosmic.root")
# maker.SetInputFile("/Users/deltutto/RealWork/CCInclusiveEventSelection/Files/Input/ubxsec_output_mc_bnbcosmic_mcc8.9_test7.root") # tune 1 mcs symm
# maker.SetOutputFile("/Users/deltutto/RealWork/CCInclusiveEventSelection/Files/Output/ubxsecana_output_bnbcosmic_mcc8.9_test7.root"); # tune 1 mcs symm
# maker.SetInputFile("/Users/deltutto/RealWork/CCInclusiveEventSelection/Files/Input/ubxsec_output_data_overlay_mcc8.9_test10.root") # Overlay
# maker.SetOutputFile("/Users/deltutto/RealWork/CCInclusiveEventSelection/Files/Output/ubxsecana_output_data_overlay_mcc8.9_test10.root"); # Overlay


maker.SetEntries(-1)
maker.SetInitialEntry(0)
maker.SetBeamSpillStart(3.2)    
maker.SetBeamSpillEnd(4.8)
maker.SetFlashShift(0.)    
maker.SetGainCalibration(198)    
maker.SetCalculatePOT(True)    
maker.SetIsData(False)
maker.SetExtraWeight(1.028); # Flux correction

maker.SetMaUpMECOff(False)

#maker.SetTargetFluxSystematic("FluxUnisim");
#maker.SetTargetFluxSystematic("kminus_PrimaryHadronNormalizat");
#maker.SetTargetFluxSystematic("kplus_PrimaryHadronFeynmanScal");
#maker.SetTargetFluxSystematic("kzero_PrimaryHadronSanfordWang");
#maker.SetTargetFluxSystematic("piminus_PrimaryHadronSWCentral");
#maker.SetTargetFluxSystematic("piplus_PrimaryHadronSWCentralS");
maker.SetTargetFluxSystematic("total")



maker.PrintConfig()

maker.MakeFile()



# BNBON

# maker.SetInputFile("/Users/deltutto/RealWork/CCInclusiveEventSelection/Files/Input/ubxsec_output_data_bnbon_mcc8.7_test6.root") # 5e19
# maker.SetOutputFile("/Users/deltutto/RealWork/CCInclusiveEventSelection/Files/Output/ubxsecana_output_bnbon_mcc8.7_test6.root") # 5e19
# maker.SetInputFile("/Users/deltutto/RealWork/CCInclusiveEventSelection/Files/Input/ubxsec_output_data_bnbon_a_mcc8.9_test4.root") # 1.6e20
# maker.SetOutputFile("/Users/deltutto/RealWork/CCInclusiveEventSelection/Files/Output/ubxsecana_output_bnbon_a_mcc8.9_test4.root") # 1.6e20
# maker.SetInputFile("/Users/deltutto/RealWork/CCInclusiveEventSelection/Files/Input/ubxsec_output_data_bnbon_mcc8.9_test8.root") # 1.6e20, 2
# maker.SetOutputFile("/Users/deltutto/RealWork/CCInclusiveEventSelection/Files/Output/ubxsecana_output_bnbon_mcc8.9_test8.root") # 1.6e20, 2
# maker.SetInputFile("/Users/deltutto/RealWork/CCInclusiveEventSelection/Files/Input/ubxsec_output_data_bnbon_mcc8.9_test7.root")
# maker.SetOutputFile("/Users/deltutto/RealWork/CCInclusiveEventSelection/Files/Output/ubxsecana_output_bnbon_mcc8.9_test7.root")
maker.SetInputFile("/Users/deltutto/RealWork/CCInclusiveEventSelection/Files/Input/ubxsec_output_data_bnbon_mcc8.9_test11.root") # 5e19, ccv corr
maker.SetOutputFile("/Users/deltutto/RealWork/CCInclusiveEventSelection/Files/Output/ubxsecana_output_data_bnbon_mcc8.9_test11.root") # 5e19, ccv corr

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

# maker.SetInputFile("/Users/deltutto/RealWork/CCInclusiveEventSelection/Files/Input/ubxsec_output_data_extbnb_mcc8.7_test6.root") # 5e19
# maker.SetOutputFile("/Users/deltutto/RealWork/CCInclusiveEventSelection/Files/Output/ubxsecana_output_extbnb_mcc8.7_test6.root") # 5e19
# maker.SetInputFile("/Users/deltutto/RealWork/CCInclusiveEventSelection/Files/Input/ubxsec_output_data_extbnb_a_mcc8.9_test4.root") # 1.6e20
# maker.SetOutputFile("/Users/deltutto/RealWork/CCInclusiveEventSelection/Files/Output/ubxsecana_output_extbnb_a_mcc8.9_test4.root") # 1.6e20
# maker.SetInputFile("/Users/deltutto/RealWork/CCInclusiveEventSelection/Files/Input/ubxsec_output_data_extbnb_mcc8.9_test8.root") # 1.6e20, 2
# maker.SetOutputFile("/Users/deltutto/RealWork/CCInclusiveEventSelection/Files/Output/ubxsecana_output_extbnb_mcc8.9_test8.root") # 1.6e20, 2
# maker.SetInputFile("/Users/deltutto/RealWork/CCInclusiveEventSelection/Files/Input/ubxsec_output_data_extbnb_mcc8.9_test7.root")
# maker.SetOutputFile("/Users/deltutto/RealWork/CCInclusiveEventSelection/Files/Output/ubxsecana_output_extbnb_mcc8.9_test7.root")
maker.SetInputFile("/Users/deltutto/RealWork/CCInclusiveEventSelection/Files/Input/ubxsec_output_data_extbnb_mcc8.9_test11.root") # 5e19, ccv corr
maker.SetOutputFile("/Users/deltutto/RealWork/CCInclusiveEventSelection/Files/Output/ubxsecana_output_data_extbnb_mcc8.9_test11.root") # 5e19, ccv corr

maker.SetEntries(-1)
maker.SetBeamSpillStart(3.65)    
maker.SetBeamSpillEnd(5.25)    
maker.SetFlashShift(0.406)    
maker.SetGainCalibration(243)    
maker.SetCalculatePOT(False)    
maker.SetIsData(True)

maker.PrintConfig()

maker.MakeFile()


