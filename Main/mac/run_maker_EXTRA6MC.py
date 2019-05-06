import sys, os

from ROOT import gROOT
gROOT.ProcessLine(".x " + os.environ['MYSW_DIR'] + "/Utils/rootlogon.C")

from ROOT import Main

maker = Main.Maker()



# BNBComisc

#maker.SetInputFile("/uboone/data/users/kirby/cc1muNp_ubxsec_ntuples/ntuples_Feb11_merge/ubxsec_output_mc_bnbcosmic_detsyst_CV.root");
#maker.SetOutputFile("/uboone/data/users/jiangl/ubxsec_static/v06_26_01_22_Feb/ubxsecana_output_mc_bnbcosmic_ubcodev06_26_01_22_detsyst_CV.root") # Run1 After Neutrino

maker.SetInputFile("/uboone/data/users/kirby/cc1muNp_ubxsec_ntuples/ntuples_Mar14_merge/ubxsec_output_mc_bnbcosmic.root")
maker.SetOutputFile("/uboone/data/users/jiangl/ubxsec_static/v06_26_01_22_Mar/ubxsecana_output_mc_bnbcosmic_ubcodev06_26_01_22_detsyst_extra_total.root") # Run1 After Neutrino


#maker.SetInputFile("/uboone/data/users/kirby/cc1muNp_ubxsec_ntuples/ntuples_Feb11_merge/ubxsec_output_mc_bnbcosmic_detsyst_DLdown.root");
#maker.SetOutputFile("/uboone/data/users/jiangl/ubxsec_static/v06_26_01_22_Feb/ubxsecana_output_mc_bnbcosmic_ubcodev06_26_01_22_detsyst_DLdown.root");
# BNB DIRT Sample

#maker.SetInputFile("/uboone/data/users/kirby/cc1muNp_ubxsec_ntuples/ntuples_Feb11_merge/ubxsec_output_mc_bnbcosmic_detsyst_DLup.root");
#maker.SetOutputFile("/uboone/data/users/jiangl/ubxsec_static/v06_26_01_22_Feb/ubxsecana_output_mc_bnbcosmic_ubcodev06_26_01_22_detsyst_DLup.root");

#maker.SetInputFile("/uboone/data/users/kirby/cc1muNp_ubxsec_ntuples/ntuples_Feb11_merge/ubxsec_output_mc_bnbcosmic_detsyst_noiseAmpUp.root");
#maker.SetOutputFile("/uboone/data/users/jiangl/ubxsec_static/v06_26_01_22_Feb/ubxsecana_output_mc_bnbcosmic_ubcodev06_26_01_22_detsyst_noiseAmpUp.root");

#maker.SetInputFile("/uboone/data/users/kirby/cc1muNp_ubxsec_ntuples/ntuples_Feb11_merge/ubxsec_output_mc_bnbcosmic_detsyst_noiseAmpDown.root");
#maker.SetOutputFile("/uboone/data/users/jiangl/ubxsec_static/v06_26_01_22_Feb/ubxsecana_output_mc_bnbcosmic_ubcodev06_26_01_22_detsyst_noiseAmpDown.root");

#maker.SetInputFile("/uboone/data/users/kirby/cc1muNp_ubxsec_ntuples/ntuples_Feb11_merge/ubxsec_output_mc_bnbcosmic_detsyst_DTup.root");
#maker.SetOutputFile("/uboone/data/users/jiangl/ubxsec_static/v06_26_01_22_Feb/ubxsecana_output_mc_bnbcosmic_ubcodev06_26_01_22_detsyst_DTup.root");

#maker.SetInputFile("/uboone/data/users/kirby/cc1muNp_ubxsec_ntuples/ntuples_Feb11_merge/ubxsec_output_mc_bnbcosmic_detsyst_DTdown.root");
#maker.SetOutputFile("/uboone/data/users/jiangl/ubxsec_static/v06_26_01_22_Feb/ubxsecana_output_mc_bnbcosmic_ubcodev06_26_01_22_detsyst_DTdown.root");

#maker.SetInputFile("/uboone/data/users/kirby/cc1muNp_ubxsec_ntuples/ntuples_Feb11_merge/ubxsec_output_mc_bnbcosmic_detsyst_squeezeResp.root");
#maker.SetOutputFile("/uboone/data/users/jiangl/ubxsec_static/v06_26_01_22_Feb/ubxsecana_output_mc_bnbcosmic_ubcodev06_26_01_22_detsyst_squeezeResp.root");

#maker.SetInputFile("/uboone/data/users/kirby/cc1muNp_ubxsec_ntuples/ntuples_Feb11_merge/ubxsec_output_mc_bnbcosmic_detsyst_stretchResp.root");
#maker.SetOutputFile("/uboone/data/users/jiangl/ubxsec_static/v06_26_01_22_Feb/ubxsecana_output_mc_bnbcosmic_ubcodev06_26_01_22_detsyst_stretchResp.root");

#maker.SetInputFile("/uboone/data/users/kirby/cc1muNp_ubxsec_ntuples/ntuples_Feb11_merge/ubxsec_output_mc_bnbcosmic_detsyst_upPEnoise.root");
#maker.SetOutputFile("/uboone/data/users/jiangl/ubxsec_static/v06_26_01_22_Feb/ubxsecana_output_mc_bnbcosmic_ubcodev06_26_01_22_detsyst_upPEnoise.root");

#maker.SetInputFile("/uboone/data/users/kirby/cc1muNp_ubxsec_ntuples/ntuples_Feb11_merge/ubxsec_output_mc_bnbcosmic_detsyst_downPEnoise.root");
#maker.SetOutputFile("/uboone/data/users/jiangl/ubxsec_static/v06_26_01_22_Feb/ubxsecana_output_mc_bnbcosmic_ubcodev06_26_01_22_detsyst_downPEnoise.root");

#maker.SetInputFile("/uboone/data/users/kirby/cc1muNp_ubxsec_ntuples/ntuples_Feb11_merge/ubxsec_output_mc_bnbcosmic_detsyst_LArG4BugFix.root");
#maker.SetOutputFile("/uboone/data/users/jiangl/ubxsec_static/v06_26_01_22_Feb/ubxsecana_output_mc_bnbcosmic_ubcodev06_26_01_22_detsyst_LArG4BugFix.root");

#maker.SetInputFile("/uboone/data/users/kirby/cc1muNp_ubxsec_ntuples/ntuples_Feb11_merge/ubxsec_output_mc_bnbcosmic_detsyst_altDeadChannels.root");
#maker.SetOutputFile("/uboone/data/users/jiangl/ubxsec_static/v06_26_01_22_Feb/ubxsecana_output_mc_bnbcosmic_ubcodev06_26_01_22_detsyst_altDeadChannels.root");

#maker.SetInputFile("/uboone/data/users/kirby/cc1muNp_ubxsec_ntuples/ntuples_Feb11_merge/ubxsec_output_mc_bnbcosmic_detsyst_dataSCE.root");
#maker.SetOutputFile("/uboone/data/users/jiangl/ubxsec_static/v06_26_01_22_Feb/ubxsecana_output_mc_bnbcosmic_ubcodev06_26_01_22_detsyst_dataSCE.root");

#maker.SetInputFile("/uboone/data/users/kirby/cc1muNp_ubxsec_ntuples/ntuples_Feb11_merge/ubxsec_output_mc_bnbcosmic_detsyst_deadSaturatedChannels.root");
#maker.SetOutputFile("/uboone/data/users/jiangl/ubxsec_static/v06_26_01_22_Feb/ubxsecana_output_mc_bnbcosmic_ubcodev06_26_01_22_detsyst_deadSaturatedChannels.root");

#maker.SetInputFile("/uboone/data/users/kirby/cc1muNp_ubxsec_ntuples/ntuples_Feb11_merge/ubxsec_output_mc_bnbcosmic_detsyst_enhancedexttpcvis.root");
#maker.SetOutputFile("/uboone/data/users/jiangl/ubxsec_static/v06_26_01_22_Feb/ubxsecana_output_mc_bnbcosmic_ubcodev06_26_01_22_detsyst_enhancedexttpcvis.root");

#maker.SetInputFile("/uboone/data/users/kirby/cc1muNp_ubxsec_ntuples/ntuples_Feb11_merge/ubxsec_output_mc_bnbcosmic_detsyst_withDIC.root");
#maker.SetOutputFile("/uboone/data/users/jiangl/ubxsec_static/v06_26_01_22_Feb/ubxsecana_output_mc_bnbcosmic_ubcodev06_26_01_22_detsyst_withDIC.root");

#maker.SetInputFile("/uboone/data/users/kirby/cc1muNp_ubxsec_ntuples/ntuples_Feb11_merge/ubxsec_output_mc_bnbcosmic_detsyst_lifetime10ms.root");
#maker.SetOutputFile("/uboone/data/users/jiangl/ubxsec_static/v06_26_01_22_Feb/ubxsecana_output_mc_bnbcosmic_ubcodev06_26_01_22_detsyst_lifetime10ms.root");



# BNB DIRT Sample


#maker.SetInputFile("/uboone/data/users/kirby/cc1muNp_ubxsec_ntuples/ntuples_Jan27_merge/ubxsec_output_mc_dirt.root");
#maker.SetOutputFile("/uboone/data/users/jiangl/ubxsec_static/v06_26_01_22_Feb/ubxsecana_output_mc_bnbdirt_ubcodev06_26_01_22.root");

# Samples for detector systematics

#maker.SetInputFile();
#maker.SetOutputFile();



maker.SetEntries(-1)
maker.SetInitialEntry(0)
maker.SetBeamSpillStart(3.1)    
maker.SetBeamSpillEnd(4.9)
maker.SetFlashShift(0.)    
maker.SetGainCalibration(198)    
maker.SetCalculatePOT(True)    
maker.SetIsData(False)
maker.SetExtraWeight(1.028); # Flux correction
maker.SetAnalysisType("cc1unp_analysis");
#maker.SetAnalysisType("ccinclusive_analysis");
maker.RecoShowerAsTrack(False);

maker.SetMaUpMECOff(False)
# maker.ScaleCosmics(0.54548) # For overlay systematics

maker.FillBootstrapGenie(False)

maker.FillBootstrapExtraSyst(True)
#maker.SetTargetExtraSystematic("model_q0q3_ccmec_HistogramWeight")
#maker.SetTargetExtraSystematic("model_q0q3_ccqe_HistogramWeight")
#maker.SetTargetExtraSystematic("model_q0q3")
#maker.SetTargetExtraSystematic("reinteractions_proton")
#maker.SetTargetExtraSystematic("reinteractions_piplus")
#maker.SetTargetExtraSystematic("reinteractions_piminus")
maker.SetTargetExtraSystematic("total")

maker.FillBootstrapFlux(False)
#maker.SetTargetFluxSystematic("FluxUnisim");
#maker.SetTargetFluxSystematic("kminus_PrimaryHadronNormalizat");
#maker.SetTargetFluxSystematic("kplus_PrimaryHadronFeynmanScal");
#maker.SetTargetFluxSystematic("kzero_PrimaryHadronSanfordWang");
#maker.SetTargetFluxSystematic("piminus_PrimaryHadronSWCentral");
#maker.SetTargetFluxSystematic("piplus_PrimaryHadronSWCentralS");
#maker.SetTargetFluxSystematic("total")

maker.FillBootstrapMCStat(False)
maker.SetNUniversesMCStat(1000)

maker.PrintConfig()

maker.MakeFile()




