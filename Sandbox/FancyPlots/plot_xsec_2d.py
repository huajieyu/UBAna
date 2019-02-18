import sys, os, math
import ROOT
from ROOT import *
# gROOT.ProcessLine(".x " + os.environ['MYSW_DIR'] + "/Utils/rootlogon.C")

f_ubana = ROOT.TFile("./xsec_file_cv.root")

xsec_data = []
xsec_data.append(f_ubana.Get("xsec_poly_muangle_mumom_cv_bin_0"))
xsec_data.append(f_ubana.Get("xsec_poly_muangle_mumom_cv_bin_1"))
xsec_data.append(f_ubana.Get("xsec_poly_muangle_mumom_cv_bin_2"))
xsec_data.append(f_ubana.Get("xsec_poly_muangle_mumom_cv_bin_3"))
xsec_data.append(f_ubana.Get("xsec_poly_muangle_mumom_cv_bin_4"))
xsec_data.append(f_ubana.Get("xsec_poly_muangle_mumom_cv_bin_5"))
xsec_data.append(f_ubana.Get("xsec_poly_muangle_mumom_cv_bin_6"))
xsec_data.append(f_ubana.Get("xsec_poly_muangle_mumom_cv_bin_7"))
xsec_data.append(f_ubana.Get("xsec_poly_muangle_mumom_cv_bin_8"))

xsec_data_unc = []
xsec_data_unc.append(f_ubana.Get("xsec_unc_poly_muangle_mumom_cv_bin_0"))
xsec_data_unc.append(f_ubana.Get("xsec_unc_poly_muangle_mumom_cv_bin_1"))
xsec_data_unc.append(f_ubana.Get("xsec_unc_poly_muangle_mumom_cv_bin_2"))
xsec_data_unc.append(f_ubana.Get("xsec_unc_poly_muangle_mumom_cv_bin_3"))
xsec_data_unc.append(f_ubana.Get("xsec_unc_poly_muangle_mumom_cv_bin_4"))
xsec_data_unc.append(f_ubana.Get("xsec_unc_poly_muangle_mumom_cv_bin_5"))
xsec_data_unc.append(f_ubana.Get("xsec_unc_poly_muangle_mumom_cv_bin_6"))
xsec_data_unc.append(f_ubana.Get("xsec_unc_poly_muangle_mumom_cv_bin_7"))
xsec_data_unc.append(f_ubana.Get("xsec_unc_poly_muangle_mumom_cv_bin_8"))

xsec_tune1 = []
xsec_tune1.append(f_ubana.Get("xsec_mc_poly_muangle_mumom_cv_bin_0"))
xsec_tune1.append(f_ubana.Get("xsec_mc_poly_muangle_mumom_cv_bin_1"))
xsec_tune1.append(f_ubana.Get("xsec_mc_poly_muangle_mumom_cv_bin_2"))
xsec_tune1.append(f_ubana.Get("xsec_mc_poly_muangle_mumom_cv_bin_3"))
xsec_tune1.append(f_ubana.Get("xsec_mc_poly_muangle_mumom_cv_bin_4"))
xsec_tune1.append(f_ubana.Get("xsec_mc_poly_muangle_mumom_cv_bin_5"))
xsec_tune1.append(f_ubana.Get("xsec_mc_poly_muangle_mumom_cv_bin_6"))
xsec_tune1.append(f_ubana.Get("xsec_mc_poly_muangle_mumom_cv_bin_7"))
xsec_tune1.append(f_ubana.Get("xsec_mc_poly_muangle_mumom_cv_bin_8"))

xsec_tune3 = []
xsec_tune3.append(f_ubana.Get("xsec_mc_alt_poly_muangle_mumom_cv_bin_0"))
xsec_tune3.append(f_ubana.Get("xsec_mc_alt_poly_muangle_mumom_cv_bin_1"))
xsec_tune3.append(f_ubana.Get("xsec_mc_alt_poly_muangle_mumom_cv_bin_2"))
xsec_tune3.append(f_ubana.Get("xsec_mc_alt_poly_muangle_mumom_cv_bin_3"))
xsec_tune3.append(f_ubana.Get("xsec_mc_alt_poly_muangle_mumom_cv_bin_4"))
xsec_tune3.append(f_ubana.Get("xsec_mc_alt_poly_muangle_mumom_cv_bin_5"))
xsec_tune3.append(f_ubana.Get("xsec_mc_alt_poly_muangle_mumom_cv_bin_6"))
xsec_tune3.append(f_ubana.Get("xsec_mc_alt_poly_muangle_mumom_cv_bin_7"))
xsec_tune3.append(f_ubana.Get("xsec_mc_alt_poly_muangle_mumom_cv_bin_8"))


# c_xsec_split = TCanvas("c_xsec_split", "multipads", 0, 0, 2800, 2000)
# c_xsec_split = TCanvas("c_xsec_split", "multipads", 0, 0, 1144,746)
c_xsec_split = TCanvas()
# c_xsec_split.Divide(3, 3, 0.01, 0.01)

x_bins = 9
costhetamu_ranges = ["-1.00 #leq cos(#theta_{#mu}^{reco}) < -0.50", "-0.50 #leq cos(#theta_{#mu}^{reco}) < 0.00", "0.00 #leq cos(#theta_{#mu}^{reco}) < 0.27",  "0.27 #leq cos(#theta_{#mu}^{reco}) < 0.45", "0.45 #leq cos(#theta_{#mu}^{reco}) < 0.62", "0.62 #leq cos(#theta_{#mu}^{reco}) < 0.76", "0.76 #leq cos(#theta_{#mu}^{reco}) < 0.86", "0.86 #leq cos(#theta_{#mu}^{reco}) < 0.94", "0.94 #leq cos(#theta_{#mu}^{reco}) #leq 1.00"]

labels = []

for i in xrange(1, x_bins+1):
# for i in xrange(1, 2):
  # c_xsec_split.cd(i)
  gPad.SetBottomMargin(0.16);
  gPad.SetLeftMargin(0.18);
  gPad.SetRightMargin(0.03);
  gPad.SetTopMargin(0.06);

  if i == 1:
    xsec_tune1[i-1].SetMinimum(-0.01)
    xsec_tune1[i-1].SetMaximum(0.45)
  if i == 2:
    xsec_tune1[i-1].SetMinimum(-0.03)
    xsec_tune1[i-1].SetMaximum(0.45)
  if i == 3:
    xsec_tune1[i-1].SetMinimum(-0.15)
    # xsec_tune1[i-1].SetMaximum(0.5)
  if i == 6:
    xsec_tune1[i-1].SetMinimum(-0.05)
    # xsec_tune1[i-1].SetMaximum(0.5)
  if i == 9:
    xsec_tune1[i-1].SetMinimum(0)
    # xsec_tune1[i-1].SetMaximum(0.5)

  xsec_tune1[i-1].SetTitle("")
  xsec_tune1[i-1].GetXaxis().SetTitle("p_{#mu}^{reco} (GeV)")
  xsec_tune1[i-1].GetYaxis().SetTitle("#frac{d^{2}#sigma}{dp_{#mu}^{reco}dcos(#theta_{#mu}^{reco})} [10^{-38} #frac{cm^{2}}{GeV n}]")
  xsec_tune1[i-1].GetXaxis().CenterTitle()
  xsec_tune1[i-1].GetYaxis().CenterTitle()
  xsec_tune1[i-1].GetXaxis().SetTitleOffset(1)
  xsec_tune1[i-1].GetXaxis().SetTitleSize(0.07)
  xsec_tune1[i-1].GetXaxis().SetLabelSize(0.06)
  xsec_tune1[i-1].GetYaxis().SetTitleOffset(1.0)
  xsec_tune1[i-1].GetYaxis().SetTitleSize(0.07)
  xsec_tune1[i-1].GetYaxis().SetLabelSize(0.06)

  gStyle.SetTitleFontSize(0.07)
  gStyle.SetTitleStyle(0)
  gStyle.SetEndErrorSize(5)

  xsec_tune1[i-1].SetLineColor(ROOT.kGreen+2)
  xsec_tune1[i-1].SetFillColor(0)
  xsec_tune3[i-1].SetLineColor(ROOT.kBlue+1)
  xsec_tune3[i-1].SetFillColor(0)
  xsec_tune3[i-1].SetLineStyle(9)

  xsec_data[i-1].SetMarkerStyle(20)
  # xsec_data[i-1].SetMarkerSize(1.3)
  xsec_data[i-1].SetMarkerSize(1.2)
  xsec_data[i-1].SetLineWidth(2)
  xsec_data_unc[i-1].SetMarkerStyle(20)
  xsec_data_unc[i-1].SetMarkerSize(0.1)
  xsec_data_unc[i-1].SetLineWidth(2)

  # genie_xsec[i-1].Draw("E2");
  # genie_xsec_main =  genie_xsec[i-1].Clone("genie_xsec_main");
  # genie_xsec_main.SetLineColor(kGreen+2);
  # genie_xsec_main.SetFillColor(0);
  # genie_xsec_main.Draw("histo same");

  xsec_tune1[i-1].Draw("histo")
  xsec_tune3[i-1].Draw("histo same")
#   nuwro_xsec[i-1].Draw("histo same")
  xsec_data[i-1].Draw("E1 X0 same")
  xsec_data_unc[i-1].Draw("E1 X0 same")

#   the_min = min(xsec_tune1[i-1].GetMinimum(), xsec_data[i-1].GetMinimum())
#   the_max = max(xsec_tune1[i-1].GetMaximum(), xsec_data[i-1].GetMaximum())

#   xsec_tune1[i-1].SetMinimum(the_min - abs(the_min) * 0.4);
#   xsec_tune1[i-1].SetMaximum(the_max * 1.5);


  costheta_label = TLatex(0.9500818,0.8657913, costhetamu_ranges[i-1]);
  costheta_label.SetNDC();
  # costheta_label.SetTextColor(kGray+2);
  costheta_label.SetTextSize(0.07);
  costheta_label.SetTextAlign(32);
  costheta_label.SetLineWidth(2);
  costheta_label.SetTextFont(42);
  costheta_label.Draw();
  labels.append(costheta_label)

  if i == 1:
    leg = TLegend(0.3705146,0.3828902,0.733819,0.6617582);
    leg.SetBorderSize(0);
    leg.SetTextSize(0.06389368);
    leg.SetLineColor(1);
    leg.SetLineStyle(1);
    leg.SetLineWidth(1);
    leg.SetFillColor(0);
    leg.SetFillStyle(0);
    leg.AddEntry(xsec_tune1[0], "GENIE Default + Emp. MEC")
    leg.AddEntry(xsec_tune3[0], "GENIE Alternative")
    leg.AddEntry(xsec_data[0],  "MicroBooNE (Stat. #oplus Syst.)", "ep")
    leg.Draw()

    ub_label = TLatex(0.92,0.71, "MicroBooNE #bf{1.6e20 POT}");
    ub_label.SetNDC();
    ub_label.SetTextColor(kGray+2);
    ub_label.SetTextSize(0.07);
    ub_label.SetTextAlign(32);
    ub_label.SetLineWidth(2);
    ub_label.SetTextFont(62);
    ub_label.Draw();

  c_xsec_split.Modified();
  c_xsec_split.cd();
  c_xsec_split.Draw()
  c_xsec_split.SaveAs('double_diff_prl' + str(i) +'.pdf')
  
  # l = TLatex(0.1,0.4,"C(x) = d #sqrt{#frac{2}{#lambdaD}}  #int^{x}_{0}cos(#frac{#pi}{2}t^{2})dt");
  # l.SetTextSize(0.15);
  # l.Draw()
  # labels.append(l)
  # pt = TPaveText(0.3739779,0.6860892,0.8386171,0.8016614,"blNDC");
  # pt.SetName("title");
  # pt.SetBorderSize(0);
  # pt.SetFillColor(0);
  # pt.SetFillStyle(0);
  # pt.SetTextFont(42);
  # pt.AddText(costhetamu_ranges[i-1]);
  # pt.Draw();

c_xsec_split.Modified();
c_xsec_split.cd();


c_xsec_split.Draw()
# ROOT.disableJSVis()
c_xsec_split.SaveAs('double_diff_prl.pdf')
c_xsec_split.SaveAs('double_diff_prl.png')

raw_input("Please press enter to exit.")
