import sys, os

from ROOT import gROOT
gROOT.ProcessLine(".x " + os.environ['MYSW_DIR'] + "/Utils/rootlogon.C")

from ROOT import Main

maker = Main.Maker()



# BNBComisc

maker.SetInputFile("/build/kirby/cc1muNp_ubxsec_pid_integration_test_larana/test_ntuples_Dec10/ubxsec_output_mc_bnbcosmic.root");
maker.SetOutputFile("/uboone/data/users/jiangl/ubxsec_static/v06_26_01_22_Dec/ubxsecana_output_mc_bnbcosmic_ubcodev06_26_01_22.root") # Run1 After Neutrino

#maker.SetInputFile("/build/kirby/cc1muNp_ubxsec_pid_integration_test_larana/test_ntuples_Dec10/ubxsec_output_mc_bnbdirt.root");
#maker.SetOutputFile("/uboone/data/users/jiangl/ubxsec_static/v06_26_01_22_Dec/ubxsecana_output_mc_bnbdirt_ubcodev06_26_01_22.root");


maker.SetEntries(-1)
maker.SetInitialEntry(0)
maker.SetBeamSpillStart(3.1)    
maker.SetBeamSpillEnd(4.9)
maker.SetFlashShift(0.)    
maker.SetGainCalibration(198)    
maker.SetCalculatePOT(True)    
maker.SetIsData(False)
maker.SetExtraWeight(1.028); # Flux correction
#maker.SetAnalysisType("cc1unp_analysis");
maker.SetAnalysisType("ccinclusive_analysis");
maker.SetMaUpMECOff(False)
# maker.ScaleCosmics(0.54548) # For overlay systematics

maker.FillBootstrapGenie(True)

maker.FillBootstrapExtraSyst(True)
#maker.SetTargetExtraSystematic("model_q0q3_ccmec_HistogramWeight")
# maker.SetTargetExtraSystematic("model_q0q3_ccqe_HistogramWeight")
# maker.SetTargetExtraSystematic("model_q0q3")
# maker.SetTargetExtraSystematic("reinteractions_proton")
# maker.SetTargetExtraSystematic("reinteractions_piplus")
# maker.SetTargetExtraSystematic("reinteractions_piminus")
maker.SetTargetExtraSystematic("total")

maker.FillBootstrapFlux(True)
#maker.SetTargetFluxSystematic("FluxUnisim");
#maker.SetTargetFluxSystematic("kminus_PrimaryHadronNormalizat");
#maker.SetTargetFluxSystematic("kplus_PrimaryHadronFeynmanScal");
#maker.SetTargetFluxSystematic("kzero_PrimaryHadronSanfordWang");
#maker.SetTargetFluxSystematic("piminus_PrimaryHadronSWCentral");
#maker.SetTargetFluxSystematic("piplus_PrimaryHadronSWCentralS");
maker.SetTargetFluxSystematic("total")

maker.FillBootstrapMCStat(False)
maker.SetNUniversesMCStat(1000)

maker.PrintConfig()

maker.MakeFile()




