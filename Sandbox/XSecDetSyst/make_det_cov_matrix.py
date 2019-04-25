import sys, os
import math

from ROOT import *
gROOT.ProcessLine(".x " + os.environ['MYSW_DIR'] + "/Utils/rootlogon.C")

from array import array

# Get the central value xsec
file_cv = TFile("xsec_file_CV.root");
xsec_onebin_cv = file_cv.Get("xsec_onebin_CV")
xsec_mumom_cv = file_cv.Get("xsec_mumom_CV")
xsec_muangle_cv = file_cv.Get("xsec_muangle_CV")
xsec_pmom_cv = file_cv.Get("xsec_pmom_CV")
xsec_pangle_cv = file_cv.Get("xsec_pangle_CV")


print "Number of muon momentum bins", xsec_mumom_cv.GetNbinsX()
print "Number of muon angle bins", xsec_muangle_cv.GetNbinsX()
print "Number of proton momentum bins", xsec_pmom_cv.GetNbinsX()
print "Number of proton angle bins", xsec_pangle_cv.GetNbinsX()


n_bins_mumom = xsec_mumom_cv.GetNbinsX()
n_bins_muangle = xsec_muangle_cv.GetNbinsX()
n_bins_pmom = xsec_pmom_cv.GetNbinsX()
n_bins_pangle = xsec_pangle_cv.GetNbinsX()

# The cov matrix
cov_matrix_mumom = TH2D("cov_matrix_mumom", "", n_bins_mumom, 0, n_bins_mumom, n_bins_mumom, 0, n_bins_mumom)
cov_matrix_mumom_frac = TH2D("cov_matrix_mumom_frac", "", n_bins_mumom, 0, n_bins_mumom, n_bins_mumom, 0, n_bins_mumom)
cov_matrix_muangle = TH2D("cov_matrix_muangle", "", n_bins_muangle, 0, n_bins_muangle, n_bins_muangle, 0, n_bins_muangle)
cov_matrix_muangle_frac = TH2D("cov_matrix_muangle_frac", "", n_bins_muangle, 0, n_bins_muangle, n_bins_muangle, 0, n_bins_muangle)

cov_matrix_pmom = TH2D("cov_matrix_pmom", "", n_bins_pmom, 0, n_bins_pmom, n_bins_pmom, 0, n_bins_pmom)
cov_matrix_pmom_frac = TH2D("cov_matrix_pmom_frac", "", n_bins_pmom, 0, n_bins_pmom, n_bins_pmom, 0, n_bins_pmom)
cov_matrix_pangle = TH2D("cov_matrix_pangle", "", n_bins_pangle, 0, n_bins_pangle, n_bins_pangle, 0, n_bins_pangle)
cov_matrix_pangle_frac = TH2D("cov_matrix_pangle_frac", "", n_bins_pangle, 0, n_bins_pangle, n_bins_pangle, 0, n_bins_pangle)

#print "CV cross section", xsec_onebin_cv.GetBinContent(1)

#print "Systematic parameter & Cross section & Total xsec perc difference \\\\"

print "CV & ", xsec_onebin_cv.GetBinContent(1), " &  0.0  \\\\"



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






det_syst_list = ["lifetime10ms", "DLdown", "DLup", "DTdown", "DTup"]



stat_err_perbin_mumom = [0.0041, 0.013, 0.010, 0.0048, 0.0028, 0.00067]
stat_err_perbin_muangle = [0.0034, 0.0026, 0.0060, 0.0038, 0.0055, 0.0087, 0.0080, 0.011, 0.018]
stat_err_perbin_onebin = 0.0091

stat_err_perbin_mumom = [0, 0, 0, 0, 0, 0]
stat_err_perbin_muangle = [0, 0, 0, 0, 0, 0, 0, 0, 0]
stat_err_perbin_onebin = 0

syst_total_xsec = 0
syst_mumom = [0, 0, 0, 0, 0, 0]

for syst_name in det_syst_list:

	# if syst_name == "nodeltaray" or syst_name == "nospacecharge" or syst_name == "stretchRes" or syst_name == "stretchResp" or syst_name == "noShortedResp":
	# 	continue

	file_name = "xsec_file_" + syst_name + ".root"
	file = TFile(file_name);
	# if (file.IsOpen()):
	# 	print "File with name", file_name, "is opened"
	# else:
	# 	print "Cannot find file", file_name
	# 	exit(0)

 	xsec_onebin        = file.Get("xsec_onebin_" + syst_name)
 	xsec_mumom         = file.Get("xsec_mumom_" + syst_name)
 	xsec_muangle       = file.Get("xsec_muangle_" + syst_name)
 	xsec_pmom         = file.Get("xsec_pmom_" + syst_name)
 	xsec_pangle       = file.Get("xsec_pangle_" + syst_name)


	# for i in xrange(0, n_bins_muangle):
	# 	for j in xrange(0, n_bins_mumom):
	# 		print syst_name, "Cross Section at ", i, j , xsec_muangle_mumom.GetBinContent(i+1, j+1), xsec_muangle_mumom_cv.GetBinContent(i+1, j+1)


	perc_diff = (xsec_onebin_cv.GetBinContent(1) - xsec_onebin.GetBinContent(1)) / xsec_onebin_cv.GetBinContent(1)
	print syst_name, " & ", xsec_onebin.GetBinContent(1), " & ", perc_diff*100, "  \\\\"
	syst_total_xsec = syst_total_xsec + abs(xsec_onebin_cv.GetBinContent(1) - xsec_onebin.GetBinContent(1))**2
	print "diff is ", abs(xsec_onebin_cv.GetBinContent(1) - xsec_onebin.GetBinContent(1))

	

	

	# Loop over bins to calculate D0

	for i in xrange(0, xsec_mumom.GetNbinsX()):
		for j in xrange(0, xsec_mumom.GetNbinsX()):
			d0 = (xsec_mumom.GetBinContent(i+1) - xsec_mumom_cv.GetBinContent(i+1))*(xsec_mumom.GetBinContent(j+1) - xsec_mumom_cv.GetBinContent(j+1))
			d0_frac = d0 / ((xsec_mumom_cv.GetBinContent(i+1))*(xsec_mumom_cv.GetBinContent(j+1)))
			cov_matrix_mumom.SetBinContent(i+1, j+1, cov_matrix_mumom.GetBinContent(i+1, j+1) + d0)
			cov_matrix_mumom_frac.SetBinContent(i+1, j+1, cov_matrix_mumom_frac.GetBinContent(i+1, j+1) + d0_frac)
			# syst_mumom[i] = syst_mumom[i] + d0_frac

        print "finish setting bin content for the muom momentum covariance matrix"
  
	for i in xrange(0, xsec_muangle.GetNbinsX()):
		for j in xrange(0, xsec_muangle.GetNbinsX()):
			d0 = (xsec_muangle.GetBinContent(i+1) - xsec_muangle_cv.GetBinContent(i+1))*(xsec_muangle.GetBinContent(j+1) - xsec_muangle_cv.GetBinContent(j+1))
			d0_frac = d0 / ((xsec_muangle_cv.GetBinContent(i+1))*xsec_muangle_cv.GetBinContent(j+1))
			cov_matrix_muangle.SetBinContent(i+1, j+1, cov_matrix_muangle.GetBinContent(i+1, j+1) + d0)
			cov_matrix_muangle_frac.SetBinContent(i+1, j+1, cov_matrix_muangle_frac.GetBinContent(i+1, j+1) + d0_frac)

        print "finish setting bin content for the muon angle covariance matrix"

	for i in xrange(0, xsec_pmom.GetNbinsX()):
		for j in xrange(0, xsec_pmom.GetNbinsX()):
			d0 = (xsec_pmom.GetBinContent(i+1) - xsec_pmom_cv.GetBinContent(i+1))*(xsec_pmom.GetBinContent(j+1) - xsec_pmom_cv.GetBinContent(j+1))
			d0_frac = d0 / ((xsec_pmom_cv.GetBinContent(i+1))*(xsec_pmom_cv.GetBinContent(j+1)))
			cov_matrix_pmom.SetBinContent(i+1, j+1, cov_matrix_pmom.GetBinContent(i+1, j+1) + d0)
			cov_matrix_pmom_frac.SetBinContent(i+1, j+1, cov_matrix_pmom_frac.GetBinContent(i+1, j+1) + d0_frac)
			# syst_pmom[i] = syst_pmom[i] + d0_frac

        print "finish setting bin content for the proton momentum covariance matrix"
  
	for i in xrange(0, xsec_pangle.GetNbinsX()):
		for j in xrange(0, xsec_pangle.GetNbinsX()):
			d0 = (xsec_pangle.GetBinContent(i+1) - xsec_pangle_cv.GetBinContent(i+1))*(xsec_pangle.GetBinContent(j+1) - xsec_pangle_cv.GetBinContent(j+1))
			d0_frac = d0 / ((xsec_pangle_cv.GetBinContent(i+1))*xsec_pangle_cv.GetBinContent(j+1))
			cov_matrix_pangle.SetBinContent(i+1, j+1, cov_matrix_pangle.GetBinContent(i+1, j+1) + d0)
			cov_matrix_pangle_frac.SetBinContent(i+1, j+1, cov_matrix_pangle_frac.GetBinContent(i+1, j+1) + d0_frac)

        print "finish setting bin content for the proton angle covariance matrix"



print "Total onebin syst err:", math.sqrt(syst_total_xsec)
print "Total onebin syst err (relative):", math.sqrt(syst_total_xsec) / xsec_onebin_cv.GetBinContent(1)

# for i in xrange(0, xsec_mumom_cv.GetNbinsX()):
# 	print "bin", i+1, "sqrt of d0_frac", math.sqrt(syst_mumom[i]) 




# gStyle.SetPalette(kDeepSea);

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

cov_file = TFile("covariance_detector_firstround.root", "RECREATE");
cov_file.cd();
cov_matrix_mumom.Write("covariance_matrix_detector_mumom");
cov_matrix_muangle.Write("covariance_matrix_detector_muangle");
cov_matrix_pmom.Write("covariance_matrix_detector_pmom");
cov_matrix_pangle.Write("covariance_matrix_detector_pangle");
cov_matrix_mumom_frac.Write("frac_covariance_matrix_detector_mumom");
cov_matrix_muangle_frac.Write("frac_covariance_matrix_detector_muangle");
cov_matrix_pmom_frac.Write("frac_covariance_matrix_detector_pmom");
cov_matrix_pangle_frac.Write("frac_covariance_matrix_detector_pangle");
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
# cov_matrix_pmom.SetMarkerColor(kWhite);
cov_matrix_pmom.SetMarkerSize(1.6);
cov_matrix_pmom.GetXaxis().SetTitle("Bin i");
cov_matrix_pmom.GetYaxis().SetTitle("Bin j");
cov_matrix_pmom.GetXaxis().CenterTitle();
cov_matrix_pmom.GetYaxis().CenterTitle();
cov_matrix_pmom.Draw("colz TEXT")

gStyle.SetPaintTextFormat("4.3f");
c_pangle = TCanvas()
# cov_matrix_pangle.SetMarkerColor(kWhite);
cov_matrix_pangle.SetMarkerSize(1.6);
cov_matrix_pangle.GetXaxis().SetTitle("Bin i");
cov_matrix_pangle.GetYaxis().SetTitle("Bin j");
cov_matrix_pangle.GetXaxis().CenterTitle();
cov_matrix_pangle.GetYaxis().CenterTitle();
cov_matrix_pangle.Draw("colz TEXT")

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










