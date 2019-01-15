import ROOT
import math
import array

def setBins(bins):
  bins["protonMom"]=[0.30, 0.36, 0.41, 0.44, 0.49, 0.53, 0.56, 0.59, 0.63, 0.73, 0.81,1.27, 1.50]
  bins["muonMom"]=[0.00, 0.18, 0.30, 0.45, 0.77, 1.28, 2.50]
  bins["protonAngle"]=[-1.00, -0.50, 0.00, 0.27, 0.45, 0.62, 0.76, 0.86, 0.94, 1.00]
  bins["muonAngle"]=[-1.00, -0.50, 0.00, 0.27, 0.45, 0.62, 0.76, 0.86, 0.94, 1.00]
  bins["thetamup"]=[0.00, 0.18, 0.30, 0.45, 0.77, 1.28, 3.14]
  return

def getRecoValue(e, var):
  val = -9999
  if var=="protonMom":
    val = getLeadingRecoMom(e)
  elif var == "muonMom":
    for i in xrange(len(e.pfp_reco_theta)):
      if e.pfp_reco_ismuoncandidate[i]:
        if inCV(e.pfp_reco_endx, e.pfp_reco_endy, e.pfp_reco_endz):
          val = e.pfp_reco_Mom_muon[i]
        else:
          val = e.pfp_reco_Mom_MCS[i]
  elif var == "protonAngle":
    val = getLeadingRecoAngle(e)
  elif var == "muonAngle":
    for i in xrange(len(e.pfp_reco_theta)):
      if e.pfp_reco_ismuoncandidate[i]:
        val = math.cos(e.pfp_reco_theta[i])
  else:
    print "unknown variable",var,"returning default value of -9999"
  return val

def getTrueValue(e,var):
  val = -9999
  if var=="protonMom":
    val=getLeadingProtonMom(e)
  elif var == "muonMom":
    val = e.true_muon_mom
  elif var == "protonAngle":
    val=getLeadingProtonAngle(e)
  elif var == "muonAngle":
    val = e.lep_costheta
  else:
    print "unknown variable",var,"returning default value of -9999"
  return val

def inCV(x,y,z):
#  print x,y,z
  if x < 10:
    return False
  if x > 245:
    return False
  if abs(y)>100:
    return False
  if z<10:
    return False
  if z > 1030:
    return False
  return True

def isSignal(e,thr,z):
#  return 1 # DATA only
  if not e.is_signal:
    return 0
  if e.ngenie_pipms>0:
    return 0
  if e.ngenie_pion0s>0:
    return 0
  if z==0:
    if (getProtonsAbove(e,thr)>0):
      return 1
    return 0
  if z==1:
    if (getProtonsAboveZ(e,thr)>0):
      return 1
    return 0
  return 0



def getLeadingProtonMom(e):
  pmax=-1
  for i in xrange(len(e.genie_mcpar_pdgcode)):
    if e.genie_mcpar_pdgcode[i] == 2212:
      px = e.genie_mcpar_px[i]
      py = e.genie_mcpar_py[i]
      pz = e.genie_mcpar_pz[i]
      p = math.sqrt(px*px + py*py + pz*pz)
      if p > pmax:
        pmax = p
  return pmax

def getLeadingProtonAngle(e):
  pmax=-1
  pind=-1
  for i in xrange(len(e.genie_mcpar_pdgcode)):
    if e.genie_mcpar_pdgcode[i] == 2212:
      px = e.genie_mcpar_px[i]
      py = e.genie_mcpar_py[i]
      pz = e.genie_mcpar_pz[i]
      p = math.sqrt(px*px + py*py + pz*pz)
      if p > pmax:
        pmax = p
        pind=i
  if pind<0:
    print "can't find leading proton true angle"
    return -9999
  if pmax<0.0000000000000001:
    print "can't find leading proton true angle - momentum==0"
    return -9999
#  if e.genie_mcpar_pz[pind]<0:
  return e.genie_mcpar_pz[pind]/pmax

def getLeadingRecoMom(e):
  pmax=-1
  for i in xrange(len(e.pfp_reco_Mom_proton)):
    if e.pfp_reco_ismuoncandidate[i]:
      continue
    p = e.pfp_reco_Mom_proton[i]
    if p > pmax:
      pmax = p
  return pmax

def getLeadingRecoAngle(e):
  maxp=-1
  pind=-1
  for i in xrange(len(e.pfp_reco_Mom_proton)):
    if e.pfp_reco_ismuoncandidate[i]:
      continue
    p = e.pfp_reco_Mom_proton[i]
    if p > maxp:
      maxp=p
      pind=i
  if pind<0:
    print "can't find leading proton reco angle"
    return -9999
  return math.cos(e.pfp_reco_theta[pind])

def getProtonsAbove(e,thr):
  N = 0
  for i in xrange(len(e.genie_mcpar_pdgcode)):
    if e.genie_mcpar_pdgcode[i] == 2212:
      px = e.genie_mcpar_px[i]
      py = e.genie_mcpar_py[i]
      pz = e.genie_mcpar_pz[i]
      p = math.sqrt(px*px + py*py + pz*pz)
      if p > thr/1000.:
        N+=1
      #if p > 0.2 and thr==200:
      #  print "counting protons, this one has a momentum of", p, "N = ",
#  if N==0:
#    print "no protons above 0 MeV"
  return N

def getProtonsAboveZ(e,thr):
  N = 0
  for i in xrange(len(e.genie_mcpar_pdgcode)):
    if e.genie_mcpar_pdgcode[i] == 2212:
      pz = e.genie_mcpar_pz[i]
      if pz > thr/1000.:
        N+=1
  return N

def getRecoAbove(e,thr):
  N = 0
  for i in xrange(len(e.pfp_reco_chi2_proton)):
    if e.pfp_reco_ismuoncandidate[i]:
      continue
    if e.pfp_reco_istrack[i]==0:
      continue
    p = e.pfp_reco_Mom_proton[i]
    if p > thr/1000.:
      N+=1
  #print "returning value of ",N
  return N


def getRecoAboveZ(e,thr):
  N = 0
  for i in xrange(len(e.pfp_reco_chi2_proton)):
    if e.pfp_reco_ismuoncandidate[i]:
      continue
    if e.pfp_reco_istrack[i]==0:
      continue
    pz = e.pfp_reco_Mom_proton[i] * math.cos(e.pfp_reco_theta[i])
    if math.fabs(pz) > thr/1000.:
      N+=1
  return N


def passChi2Cut(e):
#  ntracks = len(e.pfp_reco_chi2_proton)
  ntracks = sum(1 for x in xrange(len(e.pfp_reco_istrack)) if bool(e.pfp_reco_istrack[x]))
  if ntracks < 2:
    return False
  for i in xrange(len(e.pfp_reco_chi2_proton)):
    nhits_CP = len(e.pfp_reco_dEdx[i])
    if nhits_CP<5:
      continue
    if not (e.pfp_reco_chi2_proton[i] < 100. or e.pfp_reco_ismuoncandidate[i]):
      return False
  return True

def passMinHitsCut(e):
  maxp=-1
  pind=-1
  for i in xrange(len(e.pfp_reco_Mom_proton)):
    p = e.pfp_reco_Mom_proton[i]
    if p > maxp:
      maxp=p
      pind=i
  if pind<0:
    return False
  if len(e.pfp_reco_dEdx[pind])<5:
    return False
  return True

def passContainmentCut(e):
  for i in xrange(len(e.pfp_reco_chi2_proton)):
    if e.pfp_reco_ismuoncandidate[i]:
      continue
    if not e.pfp_reco_isshower[i] and not e.pfp_reco_istrack[i]:
      continue
    if not inCV(e.pfp_reco_startx[i], e.pfp_reco_starty[i], e.pfp_reco_startz[i]):
      return False
    if not inCV(e.pfp_reco_endx[i], e.pfp_reco_endy[i], e.pfp_reco_endz[i]):
      return False
    return True

#f = ROOT.TFile("/build/kirby/cc1muNp_ubxsec_integration/ubxsec_test_ntuples/ubxsec_output_mc_bnbcosmic.root")
#f = ROOT.TFile("/build/kirby/cc1muNp_ubxsec_pid_integration_test_larana/testing_pid_larana/ubxsec_output_mc_bnbcosmic_Nov19_test_v1.root")
#f = ROOT.TFile("/build/kirby/cc1muNp_ubxsec_pid_integration_test_larana/testing_pid_larana/ubxsec_output_mc_bnbcosmic.root")
#f = ROOT.TFile("/build/kirby/cc1muNp_ubxsec_pid_integration_test_larana/test_ntuples_Nov26/ubxsec_output_mc_bnbcosmic.root")
#f = ROOT.TFile("/build/kirby/cc1muNp_ubxsec_pid_integration_test_larana/test_ntuples_Nov26/ubxsec_output_data_onbeam.root")
f = ROOT.TFile("/build/kirby/cc1muNp_ubxsec_pid_integration_test_larana/test_ntuples_Dec20/ubxsec_output_mc_bnbcosmic.root")
#f = ROOT.TFile("/build/kirby/cc1muNp_ubxsec_pid_integration_test_larana/test_ntuples_Dec20/ubxsec_output_mc_bnbdirt.root")

ROOT.gROOT.SetBatch(1)

#t = f.Get("SimpleAna/cc1unptree")
tree = f.Get("UBXSec/tree")

#tree.MakeClass("UBXSecEvent")
#ROOT.gROOT.ProcessLine(".L /uboone/app/users/afurmans/CCNproton_ccinc/software_builds/larsoft_06_26_01_22/srcs/uboonecode/uboone/UBXSec/DataTypes/UBXSecEvent.h")
ROOT.gROOT.ProcessLine(".L /build/kirby/cc1muNp_ubxsec_pid_integration_test_larana/srcs/uboonecode/uboone/UBXSec/DataTypes/UBXSecEvent.h")
#ROOT.gROOT.LoadMacro("UBXSec.C+")

h_reso = dict()
h_trueDist = dict()
variables = ["protonMom","muonMom","protonAngle","muonAngle"]#,"thetamup"] ## worry about thetamup later...
bins = dict()
setBins(bins)
for var in variables:
  for ibin in xrange(len(bins[var])-1):
    h_trueDist[var,ibin] = ROOT.TH1D("",";;",50,0,-1)
    h_reso[var,ibin] = ROOT.TH1D("",";;",50,0,-1)
    h_trueDist[var,ibin].SetDirectory(0)
    h_reso[var,ibin].SetDirectory(0)


N_entries = tree.GetEntries()
perc = int(N_entries/100.)
dec = int(N_entries/10.)
print "running over",N_entries,"events"
i=0
for event in tree:
  i+=1
  if i%perc == 0 and i>perc-10 and i < dec:
    print int(100*float(i)/N_entries),"%"
  if i%dec == 0 and i>dec-10:
    print int(100*float(i)/N_entries),"%"
    
  e = event.ubxsec_event_split
  proton_true_mom = getLeadingProtonMom(e)
  proton_reco_mom = getLeadingRecoMom(e)
  if e.is_selected:
    if not passContainmentCut(e):
      continue
    if passChi2Cut(e) and passMinHitsCut(e):
      if not isSignal(e,300,0):
        continue
      for var in variables:
        reco = getRecoValue(e, var)
        #if var == "protonAngle" and reco<0:
        #  print "========================================================="
        #  print "proton angle (reco)", reco
        #  print "proton mom (reco)", getRecoValue(e, "protonMom")
        true = getTrueValue(e, var)
        #if var == "protonAngle" and reco<0:
        #  print "proton angle (true)", true
        #  print "proton mom (true)", getTrueValue(e, "protonMom")
        #  print "========================================================="
        for ibin in xrange(len(bins[var])-1):
          if reco > bins[var][ibin] and reco <= bins[var][ibin+1]:
            h_trueDist[var,ibin].Fill(true)
            h_reso[var,ibin].Fill(true-reco)
#            if var == "muonAngle":
#              print "muon angle (t,r)", true, reco
            break

canvSplit = ROOT.TCanvas("","",3200,800)

leg=ROOT.TLegend(0.7, 0.7, 0.89, 0.89)
leg.SetFillColor(0)
leg.SetLineColor(0)
lowLine = dict()
highLine = dict()


for var in variables:
  nbins = len(bins[var])-1
  canvSplit.Clear()
  canvSplit.Divide(nbins,2)
  for ibin in xrange(nbins):
    ######## True distribution ##############
    lowLine[var,ibin,"true"] = ROOT.TLine(bins[var][ibin], 0.0, bins[var][ibin],h_trueDist[var,ibin].GetMaximum())
    highLine[var,ibin,"true"] = ROOT.TLine(bins[var][ibin+1], 0.0, bins[var][ibin+1],h_trueDist[var,ibin].GetMaximum())
    lowLine[var,ibin,"true"].SetLineColor(2)
    highLine[var,ibin,"true"].SetLineColor(2)
    lowLine[var,ibin,"true"].SetLineStyle(2)
    highLine[var,ibin,"true"].SetLineStyle(2)
#    lowLine[var,ibin,"true"].SetDirectory(0)
#    highLine[var,ibin,"true"].SetDirectory(0)
    #canv.cd()
    h_trueDist[var,ibin].SetTitle("reco bin "+str(bins[var][ibin])+" - "+str(bins[var][ibin+1])+";True "+var+";Events")
    h_trueDist[var,ibin].SetLineColor(1)
    #h_trueDist[var,ibin].Draw("hist")
    #lowLine[var,ibin,"true"].Draw()
    #highLine[var,ibin,"true"].Draw()
    #canv.SaveAs("figures/resolutions/"+var+"_trueDistribution_bin_"+str(ibin)+".eps")
    #canv.SaveAs("figures/resolutions/"+var+"_trueDistribution_bin_"+str(ibin)+".png")
    canvSplit.cd(ibin+1)
    h_trueDist[var,ibin].Draw("hist")
    lowLine[var,ibin,"true"].Draw()
    highLine[var,ibin,"true"].Draw()

    ######## True - reco ##############
    width = bins[var][ibin+1]-bins[var][ibin]
    lowLine[var,ibin,"reso"] = ROOT.TLine(-width/2., 0.0, -width/2.,h_reso[var,ibin].GetMaximum())
    highLine[var,ibin,"reso"] = ROOT.TLine(width/2., 0.0, width/2., h_reso[var,ibin].GetMaximum())
    lowLine[var,ibin,"reso"].SetLineColor(2)
    highLine[var,ibin,"reso"].SetLineColor(2)
    lowLine[var,ibin,"reso"].SetLineStyle(2)
    highLine[var,ibin,"reso"].SetLineStyle(2)
#    lowLine[var,ibin,"reso"].SetDirectory(0)
#    highLine[var,ibin,"reso"].SetDirectory(0)
    #canv.cd()
    h_reso[var,ibin].SetTitle("reco bin "+str(bins[var][ibin])+" - "+str(bins[var][ibin+1])+";(True - Reco) "+var+";Events")
    h_reso[var,ibin].SetLineColor(1)
    #h_reso[var,ibin].Draw("hist")
    #lowLine[var,ibin,"reso"].Draw()
    #highLine[var,ibin,"reso"].Draw()
    #canv.SaveAs("figures/resolutions/"+var+"_reso_bin_"+str(ibin)+".eps")
    #canv.SaveAs("figures/resolutions/"+var+"_reso_bin_"+str(ibin)+".png")
    canvSplit.cd(nbins+ibin+1)
    h_reso[var,ibin].Draw("hist")
    lowLine[var,ibin,"reso"].Draw()
    highLine[var,ibin,"reso"].Draw()
  canvSplit.SaveAs("figures/resolutions/"+var+"allPlots.eps")
  canvSplit.SaveAs("figures/resolutions/"+var+"allPlots.png")

canv = ROOT.TCanvas()
canv.cd()
for var in variables:
  nbins = len(bins[var])-1
  for ibin in xrange(nbins):
    h_trueDist[var,ibin].Draw("hist")
    lowLine[var,ibin,"true"].Draw()
    highLine[var,ibin,"true"].Draw()
    canv.SaveAs("figures/resolutions/"+var+"_trueDistribution_bin_"+str(ibin)+".eps")
    canv.SaveAs("figures/resolutions/"+var+"_trueDistribution_bin_"+str(ibin)+".png")
    
    h_reso[var,ibin].Draw("hist")
    lowLine[var,ibin,"reso"].Draw()
    highLine[var,ibin,"reso"].Draw()
    canv.SaveAs("figures/resolutions/"+var+"_reso_bin_"+str(ibin)+".eps")
    canv.SaveAs("figures/resolutions/"+var+"_reso_bin_"+str(ibin)+".png")
