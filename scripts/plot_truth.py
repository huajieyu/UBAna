import ROOT
import math
import array
from ROOT import THStack
from ROOT import gROOT
class HistogramFile(object):
	def __init__(self, filename):
		self.filename = filename
	def __enter__(self):
		self.file = ROOT.TFile.Open(self.filename, 'read')
		return self
	def __exit__(self, exception_type, exception_value, traceback):
		self.file.Close()
	def get_histogram(self, name):
		"""Return the histogram identified by name from the file.
		"""
		# The TFile::Get() method returns a pointer to an object stored in a ROOT file.
		hist = self.file.Get(name)
		if hist:
			return hist
		else:
			raise RuntimeError('Unable to retrieve histogram named {0} from {1}'.format(name, self.filename))
interaction_list = ['total','qe','res','dis','coh','mec','other']
hist_list = ['mumom', 'pmom', 'mucostheta', 'pcostheta', 'thetamup']

with HistogramFile('/uboone/data/users/jiangl/ubxsec_static/v06_26_01_22_Apr19/ubxsecana_output_mc_bnbcosmic_ubcodev06_26_01_22_detsyst_CV.root') as f:
	
	canvas = ROOT.TCanvas('canvas', '', 500, 500)
	hist = {}
	histnew = {}
	for k in range(len(hist_list)):
		for i in range(len(interaction_list)):		
			hist[i] = f.get_histogram('h_mctruth_'+hist_list[k]+'_'+interaction_list[i])
			hist[i].SetFillColor(i+1)
			hist[i].SetTitle("")
			hist[i].GetXaxis().SetTitle(hist_list[k])
	                if i==0:
				print ('total number of events histogram')
				continue
			else:
				hist[i].Divide(hist[0])

		ths1 = THStack(hist_list[k], hist_list[k])
		legend = ROOT.TLegend(0.1, 0.1, 1., 1.)
		for j in range(len(hist)):
			gROOT.cd()
			histnew[j] = hist[j].Clone()
			histnew[j].SetFillColor(j+1)
			
			if j>0:
				legend.AddEntry(histnew[j], interaction_list[j])
				ths1.Add(histnew[j])
		ths1.Draw("hist")
		canvas.SaveAs('h_mctruth_'+hist_list[k]+'.png')
		legend.Draw()
		canvas.SaveAs('legend.png')

