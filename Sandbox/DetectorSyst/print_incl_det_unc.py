import sys, os
import math

from ROOT import *
gROOT.ProcessLine(".x " + os.environ['MYSW_DIR'] + "/Utils/rootlogon.C")

from array import array

print "start to open the root file and get the histograms"

# Get the central value xsec
file_cv = TFile("xsec_file_CV.root");
xsec_onebin_cv = file_cv.Get("xsec_onebin_CV")
xsec_mumom_cv = file_cv.Get("xsec_mumom_CV")
xsec_muangle_cv = file_cv.Get("xsec_muangle_CV")
xsec_pmom_cv = file_cv.Get("xsec_pmom_CV")
xsec_pangle_cv = file_cv.Get("xsec_pangle_CV")
xsec_thetamup_cv = file_cv.Get("xsec_thetamup_CV")

print "Number of muon momentum bins", xsec_mumom_cv.GetNbinsX()
print "Number of muon angle bins", xsec_muangle_cv.GetNbinsX()
print "Number of proton momentum bins", xsec_pmom_cv.GetNbinsX()
print "Number of proton angle bins", xsec_pangle_cv.GetNbinsX()
print "Number of thetamup angle bins", xsec_thetamup_cv.GetNbinsX()

n_bins_mumom = xsec_mumom_cv.GetNbinsX()
n_bins_muangle = xsec_muangle_cv.GetNbinsX()
n_bins_pmom = xsec_pmom_cv.GetNbinsX()
n_bins_pangle = xsec_pangle_cv.GetNbinsX()
n_bins_thetamup = xsec_thetamup_cv.GetNbinsX()

print "start to build histograms for the covariance matrix and fractional covariance matrix"

# The cov matrix
cov_matrix_mumom = TH2D("cov_matrix_mumom", "", n_bins_mumom, 0, n_bins_mumom, n_bins_mumom, 0, n_bins_mumom)
cov_matrix_mumom_frac = TH2D("cov_matrix_mumom_frac", "", n_bins_mumom, 0, n_bins_mumom, n_bins_mumom, 0, n_bins_mumom)
cov_matrix_muangle = TH2D("cov_matrix_muangle", "", n_bins_muangle, 0, n_bins_muangle, n_bins_muangle, 0, n_bins_muangle)
cov_matrix_muangle_frac = TH2D("cov_matrix_muangle_frac", "", n_bins_muangle, 0, n_bins_muangle, n_bins_muangle, 0, n_bins_muangle)

cov_matrix_pmom = TH2D("cov_matrix_pmom", "", n_bins_pmom, 0, n_bins_pmom, n_bins_pmom, 0, n_bins_pmom)
cov_matrix_pmom_frac = TH2D("cov_matrix_pmom_frac", "", n_bins_pmom, 0, n_bins_pmom, n_bins_pmom, 0, n_bins_pmom)
cov_matrix_pangle = TH2D("cov_matrix_pangle", "", n_bins_pangle, 0, n_bins_pangle, n_bins_pangle, 0, n_bins_pangle)
cov_matrix_pangle_frac = TH2D("cov_matrix_pangle_frac", "", n_bins_pangle, 0, n_bins_pangle, n_bins_pangle, 0, n_bins_pangle)

cov_matrix_thetamup = TH2D("cov_matrix_thetamup", "", n_bins_thetamup, 0, n_bins_thetamup, n_bins_thetamup, 0, n_bins_thetamup)
cov_matrix_thetamup_frac = TH2D("cov_matrix_thetamup_frac", "", n_bins_thetamup, 0, n_bins_thetamup, n_bins_thetamup, 0, n_bins_thetamup)

print "CV cross section", xsec_onebin_cv.GetBinContent(1)

print "Systematic parameter & Cross section & Total xsec perc difference \\\\"

print "CV & ", xsec_onebin_cv.GetBinContent(1), " &  0.0  \\\\"

print "Start to reset the matrices to zero"

#
# Reset matrices to zero
#

for i in xrange(0, xsec_mumom_cv.GetNbinsX()):
	for j in xrange(0, xsec_mumom_cv.GetNbinsX()):
		cov_matrix_mumom.SetBinContent(i+1, j+1, 0) 
		cov_matrix_mumom_frac.SetBinContent(i+1, j+1, 0)

for i in xrange(0, xsec_muangle_cv.GetNbinsX()):
	for j in xrange(0, xsec_muangle_cv.GetNbinsX()):
		cov_matrix_muangle.SetBinContent(i+1, j+1, 0) 
		cov_matrix_muangle_frac.SetBinContent(i+1, j+1, 0)

for i in xrange(0, xsec_pmom_cv.GetNbinsX()):
	for j in xrange(0, xsec_pmom_cv.GetNbinsX()):
		cov_matrix_pmom.SetBinContent(i+1, j+1, 0) 
		cov_matrix_pmom_frac.SetBinContent(i+1, j+1, 0)

for i in xrange(0, xsec_pangle_cv.GetNbinsX()):
	for j in xrange(0, xsec_pangle_cv.GetNbinsX()):
		cov_matrix_pangle.SetBinContent(i+1, j+1, 0) 
		cov_matrix_pangle_frac.SetBinContent(i+1, j+1, 0)

for i in xrange(0, xsec_thetamup_cv.GetNbinsX()):
	for j in xrange(0, xsec_thetamup_cv.GetNbinsX()):
		cov_matrix_thetamup.SetBinContent(i+1, j+1, 0) 
		cov_matrix_thetamup_frac.SetBinContent(i+1, j+1, 0)

print "start to set the detector systematic list and statistic bin values for analysis"

det_syst_list = ["CV", "DLdown", "DLup", "DTdown", "DTup", "LArG4BugFix", "altDeadChannels", "dataSCE", "downPEnoise", "enhancedexttpcvis","lifetime10ms", "noiseAmpDown", "noiseAmpUp", "squeezeResp", "stretchResp", "upPEnoise", "withDIC"]
#det_syst_list = ["CV", "DLdown", "DLup", "DTup", "DTdown", "noiseAmpUp", "noiseAmpDown", "squeezeResp", "stretchResp", "downPEnoise", "upPEnoise","LArG4BugFix", "altDeadChannels", "dataSCE", "deadSaturatedChannels", "enhancedexttpcvis","lifetime10ms", "withDIC"]

######YOU HAVE TO RUN THESE SYSTEMATICS ONE AT A TIME RIGHT NOW!!!!!
######DO NOT RUN MORE THAN ONE AT A TIME.

#det_syst_list = [ "CV"]
#det_syst_list = [ "DLdown"]
#det_syst_list = ["DLup"]
#det_syst_list = ["DTup"]
#det_syst_list = ["DTdown"]
#det_syst_list = ["noiseAmpUp"]
#det_syst_list = ["noiseAmpDown"]
#det_syst_list = ["squeezeResp"]
#det_syst_list = ["stretchResp"]
#det_syst_list = ["downPEnoise"]
#det_syst_list = ["upPEnoise"]
#det_syst_list = ["LArG4BugFix"]
#det_syst_list = ["altDeadChannels"]
#det_syst_list = ["dataSCE"]
#det_syst_list = ["enhancedexttpcvis"]
#det_syst_list = ["lifetime10ms"]
#det_syst_list = ["withDIC"]

#stat_err_perbin_mumom = [0.0041, 0.013, 0.010, 0.0048, 0.0028, 0.00067]
#stat_err_perbin_muangle = [0.0034, 0.0026, 0.0060, 0.0038, 0.0055, 0.0087, 0.0080, 0.011, 0.018]
#stat_err_perbin_pmom = []
#stat_err_perbin_pangle = []
#stat_err_perbin_onebin = 0.0091

stat_err_perbin_mumom = [0, 0, 0, 0, 0, 0]
stat_err_perbin_muangle = [0, 0, 0, 0, 0, 0, 0, 0, 0]
stat_err_perbin_pmom = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
stat_err_perbin_pangle = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
stat_err_perbin_thetamup = [0, 0, 0, 0, 0, 0]
stat_err_perbin_onebin = 0

syst_total_xsec = 0
syst_mumom = [0, 0, 0, 0, 0, 0]

for syst_name in det_syst_list:


	file_name = "xsec_file_" + syst_name + ".root"
	file = TFile(file_name);
        if (file.IsOpen()):
	 	print "File with name", file_name, "is opened"
        else:
	 	print "Cannot find file", file_name
	 	exit(0)

 	xsec_onebin             = file.Get("xsec_onebin_" + syst_name)
 	xsec_mumom              = file.Get("xsec_mumom_" + syst_name)
 	xsec_muangle            = file.Get("xsec_muangle_" + syst_name)
 	xsec_pmom               = file.Get("xsec_pmom_" + syst_name)
        xsec_pangle             = file.Get("xsec_pangle_" + syst_name)
 	xsec_thetamup              = file.Get("xsec_thetamup_" + syst_name)
	
        print file_name 


	print " syst name & syst variation xsec & percent diff "
	perc_diff = (xsec_onebin_cv.GetBinContent(1) - xsec_onebin.GetBinContent(1)) / xsec_onebin_cv.GetBinContent(1)
	print syst_name, " & ", xsec_onebin.GetBinContent(1), " & ", perc_diff*100, "  \\\\"
	syst_total_xsec = syst_total_xsec + abs(xsec_onebin_cv.GetBinContent(1) - xsec_onebin.GetBinContent(1))**2
	print " absolute xsec diff is ", abs(xsec_onebin_cv.GetBinContent(1) - xsec_onebin.GetBinContent(1))

	# Loop over bins to calculate D0

