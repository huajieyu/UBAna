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

	for i in xrange(0, xsec_mumom.GetNbinsX()):
		for j in xrange(0, xsec_mumom.GetNbinsX()):
			d0 = (xsec_mumom.GetBinContent(i+1) - xsec_mumom_cv.GetBinContent(i+1))*(xsec_mumom.GetBinContent(j+1) - xsec_mumom_cv.GetBinContent(j+1))
			d0_frac = d0 / ((xsec_mumom_cv.GetBinContent(i+1))*(xsec_mumom_cv.GetBinContent(j+1)))
			cov_matrix_mumom.SetBinContent(i+1, j+1, cov_matrix_mumom.GetBinContent(i+1, j+1) + d0)
			cov_matrix_mumom_frac.SetBinContent(i+1, j+1, cov_matrix_mumom_frac.GetBinContent(i+1, j+1) + d0_frac)
			# syst_mumom[i] = syst_mumom[i] + d0_frac
  
	for i in xrange(0, xsec_muangle.GetNbinsX()):
		for j in xrange(0, xsec_muangle.GetNbinsX()):
			d0 = (xsec_muangle.GetBinContent(i+1) - xsec_muangle_cv.GetBinContent(i+1))*(xsec_muangle.GetBinContent(j+1) - xsec_muangle_cv.GetBinContent(j+1))
			d0_frac = d0 / ((xsec_muangle_cv.GetBinContent(i+1))*xsec_muangle_cv.GetBinContent(j+1))
			cov_matrix_muangle.SetBinContent(i+1, j+1, cov_matrix_muangle.GetBinContent(i+1, j+1) + d0)
			cov_matrix_muangle_frac.SetBinContent(i+1, j+1, cov_matrix_muangle_frac.GetBinContent(i+1, j+1) + d0_frac)


	for i in xrange(0, xsec_pmom.GetNbinsX()):
#		for j in xrange(0, xsec_pmom.GetNbinsX()):
		for j in xrange(0, xsec_pmom.GetNbinsX()):
			d0 = (xsec_pmom.GetBinContent(i+1) - xsec_pmom_cv.GetBinContent(i+1))*(xsec_pmom.GetBinContent(j+1) - xsec_pmom_cv.GetBinContent(j+1))
			if (((xsec_pmom_cv.GetBinContent(i+1))*(xsec_pmom_cv.GetBinContent(j+1))) == 0):
				print "The bin content was 0 so skipping bins i= ", i, " and j= ", j
			else:
				d0_frac = d0 / ((xsec_pmom_cv.GetBinContent(i+1))*(xsec_pmom_cv.GetBinContent(j+1)))
			cov_matrix_pmom.SetBinContent(i+1, j+1, cov_matrix_pmom.GetBinContent(i+1, j+1) + d0)
			cov_matrix_pmom_frac.SetBinContent(i+1, j+1, cov_matrix_pmom_frac.GetBinContent(i+1, j+1) + d0_frac)
			if (i==j):
				print "The CV xsec is: " , xsec_pmom_cv.GetBinContent(i+1)
				print "The variation xsec is: ", xsec_pmom.GetBinContent(i+1)
				print "The Cov matrix value: ", cov_matrix_pmom.GetBinContent(i+1,j+1)
				print "The Fractional Cov matrix value: ", cov_matrix_pmom_frac.GetBinContent(i+1,j+1)
			# syst_pmom[i] = syst_pmom[i] + d0_frac


  
	for i in xrange(0, xsec_pangle.GetNbinsX()):
		for j in xrange(0, xsec_pangle.GetNbinsX()):
			d0 = (xsec_pangle.GetBinContent(i+1) - xsec_pangle_cv.GetBinContent(i+1))*(xsec_pangle.GetBinContent(j+1) - xsec_pangle_cv.GetBinContent(j+1))
			d0_frac = d0 / ((xsec_pangle_cv.GetBinContent(i+1))*xsec_pangle_cv.GetBinContent(j+1))
			cov_matrix_pangle.SetBinContent(i+1, j+1, cov_matrix_pangle.GetBinContent(i+1, j+1) + d0)
			cov_matrix_pangle_frac.SetBinContent(i+1, j+1, cov_matrix_pangle_frac.GetBinContent(i+1, j+1) + d0_frac)


	for i in xrange(0, xsec_thetamup.GetNbinsX()):
		for j in xrange(0, xsec_thetamup.GetNbinsX()):
			d0 = (xsec_thetamup.GetBinContent(i+1) - xsec_thetamup_cv.GetBinContent(i+1))*(xsec_thetamup.GetBinContent(j+1) - xsec_thetamup_cv.GetBinContent(j+1))
			d0_frac = d0 / ((xsec_thetamup_cv.GetBinContent(i+1))*(xsec_thetamup_cv.GetBinContent(j+1)))
			cov_matrix_thetamup.SetBinContent(i+1, j+1, cov_matrix_thetamup.GetBinContent(i+1, j+1) + d0)
			cov_matrix_thetamup_frac.SetBinContent(i+1, j+1, cov_matrix_thetamup_frac.GetBinContent(i+1, j+1) + d0_frac)
			# syst_mumom[i] = syst_mumom[i] + d0_frac

	# gStyle.SetPalette(kDeepSea);
# alpha = 1
# # stops = [    0.0000,    0.1250,    0.2500,    0.3750,    0.5000,    0.6250,    0.7500,    0.8750,    1.0000];
# # red   = [   6./255.,   8./255.,  36./255.,  91./255., 169./255., 235./255., 246./255., 240./255., 233./255.];
# # green = [   0./255.,  46./255.,  99./255., 149./255., 194./255., 220./255., 183./255., 166./255., 147./255.];
# # blue  = [ 243./255., 243./255., 240./255., 240./255., 241./255., 239./255., 186./255., 151./255., 129./255.];
# stops = [    0.0000,    0.1250,    0.2500,    0.3750,    0.5000,    0.6250,    0.7500,    0.8750,    1.0000];
# red   = [   6./255.,   8./255.,  36./255.,  91./255., 255./255., 235./255., 246./255., 240./255., 233./255.];
# green = [   0./255.,  46./255.,  99./255., 149./255., 255./255., 220./255., 183./255., 166./255., 147./255.];
# blue  = [ 243./255., 243./255., 240./255., 240./255., 255./255., 239./255., 186./255., 151./255., 129./255.];
# stops_a = array('d', stops)
# red_a = array('d', red)
# green_a = array('d', green)
# blue_a = array('d', blue)
# Idx = TColor.CreateGradientColorTable(9, stops_a, red_a, green_a, blue_a, 255, alpha);

	NRGBs = 5
	NCont = 100
	mainColour = [ 1.00, 1.00, 1.00, 1.00, 1.00 ]
	otherColour = [ 0.99,0.80, 0.60, 0.40, 0.20 ]
        #otherOtherColour = [ 0.9,0.80, 0.80, 0.80, 0.80 ]
	stops = [ 0.00, 0.05, 0.1, 0.4, 1.00 ]
	mainColourArray = array('d', mainColour)
	otherColourArray = array('d', otherColour)
        # otherOtherColourArray = array('d', otherOtherColour)
	stopsArray = array('d', stops)
	TColor.CreateGradientColorTable(NRGBs, stopsArray, mainColourArray, otherColourArray, otherColourArray, NCont)
	gStyle.SetNumberContours(NCont)


# TH2F h = cov_matrix_muangle_mumom_frac.Clone("h");

# i_label_number = 0;
# j_label_number = 0;

# for i in xrange (0, n_bins_muangle+1):

#     for (int i = 0; i <  cov_matrix_histo->GetNbinsX()+1; i++) {
#       std::ostringstream oss;
#       oss << i_label_number << "," << j_label_number;
#       if (j_label_number % _bs.GetNbinsY() == 0) {
#         i_label_number ++;
#         j_label_number = 0;
#       }
#       j_label_number++;
#       std::string label = oss.str();
#       h->GetXaxis()->SetBinLabel(i,label.c_str());
#       h->GetYaxis()->SetBinLabel(i,label.c_str());
#     }



# cov_matrix_mumom.SetMarkerColor(kWhite);
# cov_matrix_mumom.SetMarkerSize(1.6);
# cov_matrix_mumom.GetXaxis().SetTitle("Bin i");
# cov_matrix_mumom.GetXaxis().SetTitle("Bin j");
# cov_matrix_mumom.GetXaxis().CenterTitle();
# cov_matrix_mumom.GetYaxis().CenterTitle();
# cov_matrix_mumom
# cov_matrix_mumom
# cov_matrix_mumom.Draw("colz TEXT")

	cov_file_name = "covariance_detector_" + syst_name +".root"
#	cov_file = TFile("covariance_detector.root", "RECREATE");
	cov_file = TFile(cov_file_name, "RECREATE");
	cov_file.cd();
	cov_matrix_mumom.Write("covariance_matrix_detector_mumom");
	cov_matrix_muangle.Write("covariance_matrix_detector_muangle");
	cov_matrix_pmom.Write("covariance_matrix_detector_pmom");
	cov_matrix_pangle.Write("covariance_matrix_detector_pangle");
	cov_matrix_thetamup.Write("covariance_matrix_detector_thetamup");
	cov_matrix_mumom_frac.Write("frac_covariance_matrix_detector_mumom");
	cov_matrix_muangle_frac.Write("frac_covariance_matrix_detector_muangle");
	cov_matrix_pmom_frac.Write("frac_covariance_matrix_detector_pmom");
	cov_matrix_pangle_frac.Write("frac_covariance_matrix_detector_pangle");
	cov_matrix_thetamup_frac.Write("frac_covariance_matrix_detector_thetamup");
	cov_file.Close();

	gStyle.SetPaintTextFormat("4.3f");
	c_mumom = TCanvas()
# cov_matrix_mumom_frac.SetMarkerColor(kWhite);
	cov_matrix_mumom_frac.SetMarkerSize(1.6);
	cov_matrix_mumom_frac.GetXaxis().SetTitle("Bin i");
	cov_matrix_mumom_frac.GetYaxis().SetTitle("Bin j");
	cov_matrix_mumom_frac.GetXaxis().CenterTitle();
	cov_matrix_mumom_frac.GetYaxis().CenterTitle();
	cov_matrix_mumom_frac.Draw("colz TEXT")

	gStyle.SetPaintTextFormat("4.3f");
	c_muangle = TCanvas()
# cov_matrix_muangle_frac.SetMarkerColor(kWhite);
	cov_matrix_muangle_frac.SetMarkerSize(1.6);
	cov_matrix_muangle_frac.GetXaxis().SetTitle("Bin i");
	cov_matrix_muangle_frac.GetYaxis().SetTitle("Bin j");
	cov_matrix_muangle_frac.GetXaxis().CenterTitle();
	cov_matrix_muangle_frac.GetYaxis().CenterTitle();
	cov_matrix_muangle_frac.Draw("colz TEXT")

	gStyle.SetPaintTextFormat("4.3f");
	c_pmom = TCanvas()
# cov_matrix_pmom_frac.SetMarkerColor(kWhite);
	cov_matrix_pmom_frac.SetMarkerSize(1.6);
	cov_matrix_pmom_frac.GetXaxis().SetTitle("Bin i");
	cov_matrix_pmom_frac.GetYaxis().SetTitle("Bin j");
	cov_matrix_pmom_frac.GetXaxis().CenterTitle();
	cov_matrix_pmom_frac.GetYaxis().CenterTitle();
	cov_matrix_pmom_frac.Draw("colz TEXT")


	gStyle.SetPaintTextFormat("4.3f");
	c_pangle = TCanvas()
# cov_matrix_pangle_frac.SetMarkerColor(kWhite);
	cov_matrix_pangle_frac.SetMarkerSize(1.6);
	cov_matrix_pangle_frac.GetXaxis().SetTitle("Bin i");
	cov_matrix_pangle_frac.GetYaxis().SetTitle("Bin j");
	cov_matrix_pangle_frac.GetXaxis().CenterTitle();
	cov_matrix_pangle_frac.GetYaxis().CenterTitle();
	cov_matrix_pangle_frac.Draw("colz TEXT")

	gStyle.SetPaintTextFormat("4.3f");
	c_thetamup = TCanvas()
# cov_matrix_thetamup_frac.SetMarkerColor(kWhite);
	cov_matrix_thetamup_frac.SetMarkerSize(1.6);
	cov_matrix_thetamup_frac.GetXaxis().SetTitle("Bin i");
	cov_matrix_thetamup_frac.GetYaxis().SetTitle("Bin j");
	cov_matrix_thetamup_frac.GetXaxis().CenterTitle();
	cov_matrix_thetamup_frac.GetYaxis().CenterTitle();
	cov_matrix_thetamup_frac.Draw("colz TEXT")

#gStyle.SetPaintTextFormat("4.2f");
#c_muangle_mumom = TCanvas()
# cov_matrix_muangle_mumom_frac.SetMarkerColor(kWhite);
#cov_matrix_muangle_mumom_frac.SetMarkerSize(1.1);
#cov_matrix_muangle_mumom_frac.GetXaxis().SetTitle("Bin ij");
#cov_matrix_muangle_mumom_frac.GetYaxis().SetTitle("Bin mn");
#cov_matrix_muangle_mumom_frac.GetXaxis().CenterTitle();
#cov_matrix_muangle_mumom_frac.GetYaxis().CenterTitle();
#cov_matrix_muangle_mumom_frac.Draw("colz TEXT")

#gStyle.SetPaintTextFormat("4.2f");
#c_muangle_mumom_poly = TCanvas()
# cov_matrix_muangle_mumom_frac.SetMarkerColor(kWhite);
#cov_matrix_poly_muangle_mumom_frac.SetMarkerSize(1.1);
#cov_matrix_poly_muangle_mumom_frac.GetXaxis().SetTitle("Bin i");
#cov_matrix_poly_muangle_mumom_frac.GetYaxis().SetTitle("Bin j");
#cov_matrix_poly_muangle_mumom_frac.GetXaxis().CenterTitle();
#cov_matrix_poly_muangle_mumom_frac.GetYaxis().CenterTitle();
#cov_matrix_poly_muangle_mumom_frac.SetMaximum(1.5)
#cov_matrix_poly_muangle_mumom_frac.SetMinimum(-1.5)
#cov_matrix_poly_muangle_mumom_frac.Draw("colz TEXT")

# gStyle.SetPaintTextFormat("4.3f");
# c_mumom = TCanvas()
# cov_matrix_mumom.SetMarkerColor(kWhite);
# cov_matrix_mumom.SetMarkerSize(1.6);
# cov_matrix_mumom.GetXaxis().SetTitle("Bin i");
# cov_matrix_mumom.GetYaxis().SetTitle("Bin j");
# cov_matrix_mumom.GetXaxis().CenterTitle();
# cov_matrix_mumom.GetYaxis().CenterTitle();
# cov_matrix_mumom.Draw("colz TEXT")

# gStyle.SetPaintTextFormat("4.3f");
# c_muangle = TCanvas()
# cov_matrix_muangle.SetMarkerColor(kWhite);
# cov_matrix_muangle.SetMarkerSize(1.6);
# cov_matrix_muangle.GetXaxis().SetTitle("Bin i");
# cov_matrix_muangle.GetYaxis().SetTitle("Bin j");
# cov_matrix_muangle.GetXaxis().CenterTitle();
# cov_matrix_muangle.GetYaxis().CenterTitle();
# cov_matrix_muangle.Draw("colz TEXT")

# gStyle.SetPaintTextFormat("4.2f");
# c_muangle_mumom = TCanvas()
# cov_matrix_muangle_mumom.SetMarkerColor(kWhite);
# cov_matrix_muangle_mumom.SetMarkerSize(1.1);
# cov_matrix_muangle_mumom.GetXaxis().SetTitle("Bin ij");
# cov_matrix_muangle_mumom.GetYaxis().SetTitle("Bin mn");
# cov_matrix_muangle_mumom.GetXaxis().CenterTitle();
# cov_matrix_muangle_mumom.GetYaxis().CenterTitle();
# cov_matrix_muangle_mumom.Draw("colz")

raw_input("Please press enter to exit.")
