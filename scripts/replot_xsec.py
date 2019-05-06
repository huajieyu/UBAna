import ROOT
import math
import array

#
#def GetError(hist_name)

#  if hist_name=="h_trkmom"
#     stat_err=[]
#     syst_err=[]
#  elif hist_name=="h_muangle"
#     stat_err=[]
#     syst_err=[]
#  elif hist_name=="h_trkpmom"
#     stat_err=[]
#     syst_err=[]
#  elif hist_name=="h_pangle"
#     stat_err=[]
#     syst_err=[]
#  else hist_name=="h_thetamup"
#     stat_err=[]
#     syst_err=[]
#  return

#get the histograms first
hist_names = []
hist_names.append("h_trkmom")
hist_names.append("h_muangle")
hist_names.append("h_pmom")
hist_names.append("h_pangle")
hist_names.append("h_thetamup")
#loop over the hist_names and set the syst and stat error bars
print hist_names   
	
f = ROOT.TFile("../xsec_file_cv.root")
ROOT.gROOT.SetBatch(0)
canvas = ROOT.TCanvas('canvas', '', 500, 500)
#loop over all the histograms

i=0
while i<len(hist_names):
	print'start to fill histogram',(hist_names[i])
	if hist_names=="h_trkmom":
		xsec_hist=f.Get("xsec_mumom_cv"); canvas.cd(); xsec_hist.Draw()
	i += 1


