######## Andy's python script for making efficiency plots

import ROOT
import math
import array

def dumpEvent(e):
  print "======== DUMPING EVENT =============="
  for i in xrange(len(e.genie_mcpar_pdgcode)):
    print "PDG", e.genie_mcpar_pdgcode[i]
    print "start position (", e.genie_mcpar_startx[i],",",e.genie_mcpar_starty[i],",", e.genie_mcpar_startz[i],")"
    p = math.sqrt(e.genie_mcpar_px[i]*e.genie_mcpar_px[i] + e.genie_mcpar_py[i]*e.genie_mcpar_py[i] +e.genie_mcpar_pz[i]*e.genie_mcpar_pz[i])
    print "momentum vector (", e.genie_mcpar_px[i],",",e.genie_mcpar_py[i],",",e.genie_mcpar_pz[i],") - total =",p
  print " == Reco info =="
  for i in xrange(len(e.pfp_reco_chi2_proton)):
    if not (e.pfp_reco_isshower[i] or e.pfp_reco_istrack[i]):
      continue
    if (e.pfp_reco_ismuoncandidate[i]):
      print "muon candidate PDG",e.pfp_truth_pdg[i],", status",e.pfp_truth_status[i]
    else:
      print "proton candidate",i,"PDG",e.pfp_truth_pdg[i],", status",e.pfp_truth_status[i]

  return

def modeToText(i):
  if i==0:
    return "QE"
  if i==1:
    return "RES"
  if i==2:
    return "DIS"
  if i==3:
    return "COH"
  if i==10:
    return "MEC"
  else:
    return "bkg"

def inCV(x,y,z):
#  print x,y,z
  if x <= 10.: 
    return False
  if x >= 246.35:
    return False
  if math.fabs(y)>=96.5:
    return False
  if z<=10.:
    return False
  if z >= 1026.8:
    return False
  return True
## These were my original numbers
#  if x < 10: 
#    return False
#  if x > 245:
#    return False
#  if abs(y)>100:
#    return False
#  if z<10:
#    return False
#  if z > 1030:
#    return False
#  return True

def isSignal(e,thr,z):
#  return 1 # DATA only
  if not e.is_signal:
    return 0
  if e.ngenie_pipms>0:
    return 0
  if e.ngenie_pion0s>0:
    return 0
  if e.ngenie_electrons>0:
    return 0
  if e.ngenie_muons<1:
    return 0
  if getTrueValue(e,"muonMom") < 0.1:
    return 0
  if getLeadingProtonMom(e)>1.2:
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

def isSignalOld(e,thr,z):
#  return 1 # DATA only
  if not e.is_signal:
    return 0
  if e.ngenie_pipms>0:
    return 0
  if e.ngenie_pion0s>0:
    return 0
#  if getTrueValue(e,"muonMom") < 0.1:
#    return 0
#  if getLeadingProtonMom(e)>1.2:
#    return 0
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


def getLeadingProtonEnergy(e):
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

def getProtonsAbove(e,thr):
  N = 0
  for i in xrange(len(e.genie_mcpar_pdgcode)):
    #if e.genie_mcpar_pdgcode[i] == 2212 and e.genie_mcpar_status[i]==1:
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
#  if e.is_signal:
#    if thr==300:
#      if e.ngenie_protons_300 != N:
#        print "Number of protons comparison - Libo, Me", e.ngenie_protons_300, N
  return N

def getProtonsAboveZ(e,thr):
  N = 0
  for i in xrange(len(e.genie_mcpar_pdgcode)):
    if e.genie_mcpar_pdgcode[i] == 2212:
      pz = e.genie_mcpar_pz[i]
      if pz > thr/1000.:
        N+=1
  return N

def getRecoAbove(e,thr, considerShowerTracks):
  N = 0
  for i in xrange(len(e.pfp_reco_chi2_proton)):
    if not e.pfp_reco_isshower[i] and not e.pfp_reco_istrack[i]:
      continue
    if e.pfp_reco_ismuoncandidate[i]:
      continue
    if not considerShowerTracks and e.pfp_reco_istrack[i]==0:
      continue
    p = e.pfp_reco_Mom_proton[i]
    if p > thr/1000.:
      N+=1
  #print "returning value of ",N
  return N


def getRecoAboveZ(e,thr, considerShowerTracks):
  N = 0
  for i in xrange(len(e.pfp_reco_chi2_proton)):
    if not e.pfp_reco_isshower[i] and not e.pfp_reco_istrack[i]:
      continue
    if e.pfp_reco_ismuoncandidate[i]:
      continue
    if not considerShowerTracks and e.pfp_reco_istrack[i]==0:
      continue
    pz = e.pfp_reco_Mom_proton[i] * math.cos(e.pfp_reco_theta[i])
    if math.fabs(pz) > thr/1000.:
      N+=1
  return N

def getRecoProtonDir(e):
  maxp=-1
  pind=-1
  vec = ROOT.TVector3(0,0,0)
  for i in xrange(len(e.pfp_reco_Mom_proton)):
    if e.pfp_reco_ismuoncandidate[i]:
      continue
    p = e.pfp_reco_Mom_proton[i]
    if p > maxp:
      maxp=p
      pind=i
  if pind<0:
    #print "can't find leading proton reco angle"
    return vec
  vec.SetMagThetaPhi(1, e.pfp_reco_theta[pind], e.pfp_reco_phi[pind])
  return vec # math.cos(e.pfp_reco_theta[pind])

def getTrueValue(e,var):
  val = -9999
  if var=="protonMom":
    val=getLeadingProtonMom(e)
  elif var == "muonMom":
    val = e.true_muon_mom
  elif var == "protonAngle":
    val=math.cos(getTrueProtonDir(e).Theta())
  elif var == "muonAngle":
    val = e.lep_costheta
  elif var == "thetamup":
    val = getTrueMuonDir(e).Angle(getTrueProtonDir(e))
  elif var == "costhetamup":
    val = getTrueMuonDir(e).Dot(getTrueProtonDir(e))
  elif var == "Enu":
    val = e.nu_e
  else:
    print "unknown variable",var,"returning default value of -9999"
  return val


def passPhaseSpaceCuts(e, considerShowerTracks):
  maxp=-1
  pind=-1
  for i in xrange(len(e.pfp_reco_Mom_proton)):
    if not e.pfp_reco_isshower[i] and not e.pfp_reco_istrack[i]:
      continue # skip the neutrino...
    if not considerShowerTracks and e.pfp_reco_isshower[i]:
      continue
    if e.pfp_reco_ismuoncandidate[i]:
      if inCV(e.pfp_reco_endx[i], e.pfp_reco_endy[i], e.pfp_reco_endz[i]):
        muonMom = e.pfp_reco_Mom_muon[i]
      else:
        muonMom = e.pfp_reco_Mom_MCS[i]
      if muonMom<0.1:
        return False
    else:
      p = e.pfp_reco_Mom_proton[i]
      if p > maxp:
        maxp=p
        pind=i
  if pind<0:
    return False
  if e.pfp_reco_Mom_proton[pind]<0.3:
    return False
  if e.pfp_reco_Mom_proton[pind]>1.2:
    return False
  return True

def passChi2Cut(e, considerShowerTracks):
#  ntracks = len(e.pfp_reco_chi2_proton)
  #ntracks = sum(1 for x in xrange(len(e.pfp_reco_istrack)) if bool(e.pfp_reco_istrack[x]))
  #if ntracks < 2:
  #  return False
  for i in xrange(len(e.pfp_reco_chi2_proton)):
    if e.pfp_reco_ismuoncandidate[i]:
      continue
    if not e.pfp_reco_isshower[i] and not e.pfp_reco_istrack[i]:
      continue
    if not considerShowerTracks and e.pfp_reco_isshower[i]:
      continue
    nhits_CP = len(e.pfp_reco_dEdx[i])
    if nhits_CP<5:
      continue
    if not (e.pfp_reco_chi2_proton[i] < 88.):
      return False
  return True

def passMinHitsCut(e, considerShowerTracks):
  maxp=-1
  pind=-1
  for i in xrange(len(e.pfp_reco_Mom_proton)):
    if e.pfp_reco_ismuoncandidate[i]:
      continue
    if not considerShowerTracks and e.pfp_reco_isshower[i]:
      continue
    p = e.pfp_reco_Mom_proton[i]
    if p > maxp:
      maxp=p
      pind=i
  if pind<0:
    return False
  if len(e.pfp_reco_dEdx[pind])<5:
    return False
  return True

def passNtracksCut(e, considerShowerTracks):
  for i in xrange(len(e.pfp_reco_chi2_proton)):
    if e.pfp_reco_ismuoncandidate[i]:
      continue
    if not e.pfp_reco_isshower[i] and not e.pfp_reco_istrack[i]:
      continue
    if not considerShowerTracks and e.pfp_reco_isshower[i]:
      continue
    return True # if you got to this point, you found something that is not the muon and not the neutrino
  return False # if you got to this point, the event only has a muon candidate and a neutrino
    
  

def passContainmentCut(e, considerShowerTracks):
#  return True # HERE - Temporary check...
  for i in xrange(len(e.pfp_reco_chi2_proton)):
    if e.pfp_reco_ismuoncandidate[i]:
      continue
    if e.pfp_reco_isshower[i]==0 and e.pfp_reco_istrack[i]==0:
      continue
    if e.pfp_reco_numtracks[i]!=1:
      continue
    if not considerShowerTracks and e.pfp_reco_isshower[i]:
      continue
    if not inCV(e.pfp_reco_startx[i], e.pfp_reco_starty[i], e.pfp_reco_startz[i]):
      #print "Event fails inCV cut. Start point x,y,z", e.pfp_reco_startx[i], e.pfp_reco_starty[i], e.pfp_reco_startz[i]
      #print "istrack, isshower",e.pfp_reco_isshower[i], e.pfp_reco_istrack[i]
      #print type(e.pfp_reco_istrack[i])
      #print "ismuondcandidate",e.pfp_reco_ismuoncandidate[i]
      #print "numtracks",e.pfp_reco_numtracks[i]
      #print "end x,y,z", e.pfp_reco_endx[i], e.pfp_reco_endy[i], e.pfp_reco_endz[i]
      #print "=========== event, subrun = ",e.event,e.subrun,"==================="
      return False # HERE - temporary - check if start containment does anything
    if not inCV(e.pfp_reco_endx[i], e.pfp_reco_endy[i], e.pfp_reco_endz[i]):
      return False # HERE - temporary - check if end containment does anything
  return True

def IsCosmic(e, considerShowerTracks):
  for i in xrange(len(e.pfp_reco_chi2_proton)):
    if not (e.pfp_reco_isshower[i] or e.pfp_reco_istrack[i]):
      continue
    if not considerShowerTracks and e.pfp_reco_isshower[i]:
      continue
    if e.pfp_truth_origin[i]!=1:
      return True
  return False


def getRecoProtonIndex(e, considerShowerTracks):
  maxp=-1
  pind=-1
  for i in xrange(len(e.pfp_reco_Mom_proton)):
    if (not considerShowerTracks) and e.pfp_reco_isshower[i]:
      continue
    if e.pfp_reco_ismuoncandidate[i]:
      continue
    p = e.pfp_reco_Mom_proton[i]
    if p > maxp:
      maxp=p
      pind=i
  if pind<0:
    print "can't find leading proton index"
  return pind

def getTrueProtonIndex(e):
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
#  if pind<0:
#    print "can't find true leading proton index"
  return pind

def getTrueProtonDir(e):
  pind = getTrueProtonIndex(e)
  direction = ROOT.TVector3(0,0,0)
  if pind<0:
    return direction
  px = e.genie_mcpar_px[pind]
  py = e.genie_mcpar_py[pind]
  pz = e.genie_mcpar_pz[pind]
  p = math.sqrt(px*px + py*py + pz*pz)
  if p>0:
    direction.SetXYZ(px/p, py/p, pz/p)
  return direction


def muCandpCandTruth(e, considerShowerTracks):
  for i in xrange(len(e.pfp_reco_chi2_proton)):
    if not (e.pfp_reco_isshower[i] or e.pfp_reco_istrack[i]):
      continue
    if not considerShowerTracks and e.pfp_reco_isshower[i]:
      continue
    if (e.pfp_reco_ismuoncandidate[i]):
      if(e.pfp_truth_pdg[i]!=13 or e.pfp_truth_status[i]!=1):
        return False
    else:
      if(e.pfp_truth_pdg[i]!=2212 or e.pfp_truth_status[i]!=1):
        return False
  return True

def protonIsLeading(e, considerShowerTracks):
  pmass = 0.938270329696
  pind_truth = getTrueProtonIndex(e)
  true_E = e.genie_mcpar_energy[pind_truth] - pmass
  pind = getRecoProtonIndex(e, considerShowerTracks)
  reco_true_E = e.pfp_truth_KE[pind]
  if math.fabs(true_E - reco_true_E) < 0.01:
    return True
  else:
    #print "not found leading proton. Selected proton true E =",reco_true_E,"but true leading E =",true_E
    #print "difference = ", math.fabs(true_E - reco_true_E)
    return False


def getLeadingPionIndexGENIE(e):
  imax=-1
  pmax=-1
  for i in xrange((len(e.genie_mcpar_pdgcode))):
    if abs(e.genie_mcpar_pdgcode[i])==211:
      px = e.genie_mcpar_px[i]
      py = e.genie_mcpar_py[i]
      pz = e.genie_mcpar_pz[i]
      p = math.sqrt(px*px+py*py+pz*pz)
      if p > pmax:
        pmax=p
        imax=i
  return imax

def getLeadingPionIndexGEANT(e):
  p_genie = getLeadingPionMom(e)
  for i in xrange(len(e.geant_mcpar_px)):
    px = e.geant_mcpar_px[i]
    py = e.geant_mcpar_py[i]
    pz = e.geant_mcpar_pz[i]
    p = math.sqrt(px*px + py*py + pz*pz)
    if math.fabs(p-p_genie)<1e-3:
      i_geant = i
      return i_geant


def getLeadingPionMom(e):
  i = getLeadingPionIndexGENIE(e)
  if i<0:
    print "can't find leading pion"
    return -1
  px = e.genie_mcpar_px[i]
  py = e.genie_mcpar_py[i]
  pz = e.genie_mcpar_pz[i]
  p = math.sqrt(px*px+py*py+pz*pz)
  return p

def getLeadingPionAngle(e):
  i = getLeadingPionIndexGENIE(e)
  if i<0:
    print "can't find leading pion"
    return -1
  px = e.genie_mcpar_px[i]
  py = e.genie_mcpar_py[i]
  pz = e.genie_mcpar_pz[i]
  p = math.sqrt(px*px+py*py+pz*pz)
  costheta = pz/p
  return costheta

def simpleParticleIndex(pdg):
  if pdg==-999:
    return -1
  elif abs(pdg)==13:
    return 0
  elif abs(pdg)==211:
    return 1
  elif abs(pdg)==2212:
    return 2
  elif abs(pdg)==11:
    return 3
  else:
    return 4

#f = ROOT.TFile("/build/kirby/cc1muNp_ubxsec_integration/ubxsec_test_ntuples/ubxsec_output_mc_bnbcosmic.root")
#f = ROOT.TFile("/build/kirby/cc1muNp_ubxsec_pid_integration_test_larana/testing_pid_larana/ubxsec_output_mc_bnbcosmic_Nov19_test_v1.root")
#f = ROOT.TFile("/build/kirby/cc1muNp_ubxsec_pid_integration_test_larana/testing_pid_larana/ubxsec_output_mc_bnbcosmic.root")
#f = ROOT.TFile("/build/kirby/cc1muNp_ubxsec_pid_integration_test_larana/test_ntuples_Nov26/ubxsec_output_mc_bnbcosmic.root")
#f = ROOT.TFile("/build/kirby/cc1muNp_ubxsec_pid_integration_test_larana/test_ntuples_Nov26/ubxsec_output_data_onbeam.root")
#f = ROOT.TFile("/build/kirby/cc1muNp_ubxsec_pid_integration_test_larana/test_ntuples_Dec20/ubxsec_output_mc_bnbcosmic.root")
f = ROOT.TFile("/uboone/data/users/kirby/cc1muNp_ubxsec_ntuples/ntuples_Mar20_merge/ubxsec_output_mc_bnbcosmic.root")
#fT3 = ROOT.TFile("/uboone/data/users/kirby/cc1muNp_ubxsec_ntuples/ntuples_Mar20_merge/ubxsec_output_mc_bnbcosmic.root")
fT3 = ROOT.TFile("/uboone/data/users/kirby/cc1muNp_ubxsec_ntuples/ntuples_Apr05_merge/ubxsec_output_mc_bnbcosmic_Tune3_merged.root")

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptStat(0)
considerShowerTracks = True

#t = f.Get("SimpleAna/cc1unptree")
tree = dict()
tree["MC"] = f.Get("UBXSec/tree")
tree["T3"] = fT3.Get("UBXSec/tree")

#tree.MakeClass("UBXSecEvent")
#ROOT.gROOT.ProcessLine(".L /uboone/app/users/afurmans/CCNproton_ccinc/software_builds/larsoft_06_26_01_22/srcs/uboonecode/uboone/UBXSec/DataTypes/UBXSecEvent.h")
#ROOT.gROOT.ProcessLine(".L /build/kirby/cc1muNp_ubxsec_pid_integration_test_larana/srcs/uboonecode/uboone/UBXSec/DataTypes/UBXSecEvent.h")
ROOT.gROOT.ProcessLine(".L /build/kirby/cc1muNp_pandora_update_v27_test/srcs/uboonecode/uboone/UBXSec/DataTypes/UBXSecEvent.h")
#ROOT.gROOT.LoadMacro("UBXSec.C+")

thresholds = [0,200,250,300,350,400]
N_above=dict()
N_above_z=dict()
N_reco=dict()
N_reco_z=dict()

h_protonSmearing=dict()
h_signal=dict()
for thr in thresholds:
  h_protonSmearing[thr,"Z"] = ROOT.TH2D("h_protonSmearing"+str(thr)+"Z","h_protonSmearing"+str(thr)+"Z",4,0,4,4,0,4)
  h_protonSmearing[thr,"tot"] = ROOT.TH2D("h_protonSmearing"+str(thr)+"tot","h_protonSmearing"+str(thr)+"tot",4,0,4,4,0,4)

for thr in thresholds:
  N_above[thr]=0
  N_above_z[thr]=0

variables = ["pz"]

xbins = array.array("d",[0.30, 0.36, 0.41, 0.44, 0.49, 0.53, 0.56, 0.59, 0.63, 0.73, 0.81,1.27, 1.50])#[0.2, 0.3, 0.38, 0.42, 0.45, 0.5, 0.54, 0.57, 0.6, 0.65, 0.75, 0.85, 1.25, 1.5])
nbins = len(xbins)-1

h_all_true_mom = ROOT.TH1D("h_all_true_mom",";p_{p} / GeV c^{-1};",nbins,xbins)
h_sel_true_mom = h_all_true_mom.Clone("h_sel_true_mom")

h_all_true_mom_thr = ROOT.TH1D("h_all_true_mom_thr",";p_{p} / GeV c^{-1};",nbins,xbins)
h_sel_true_mom_thr = h_all_true_mom.Clone("h_sel_true_mom_thr")

#tree.Draw("abs(ubxsec_event_split.geant_mcpar_pz)","ubxsec_event_split.geant_mcpar_pdgcode==2212")
#raw_input()
#tree.Draw("ubxsec_event_split.pfp_reco_Mom_proton","ubxsec_event_split.pfp_reco_chi2_proton<100 && ubxsec_event_split.pfp_reco_Mom_proton>0")
#raw_input()
N_entries = tree["MC"].GetEntries()
perc = int(N_entries/100.)
dec = int(N_entries/10.)
print "running over",N_entries,"events"
i=0
N_signal=dict()
N_signal_z=dict()
N_selected=dict()
N_selected_z=dict()
h_protonMom = dict()
h_protonAngle = dict()
h_muonMom = dict()
h_muonAngle = dict()
h_protonMom_old = dict()
h_protonMom_fullPS = dict()
h_protonAngle_fullPS = dict()
h_muonMom_fullPS = dict()
h_muonAngle_fullPS = dict()
bins=dict()
bins["protonMom"]=array.array("d",[0.30, 0.41, 0.495, 0.56, 0.62, 0.68, 0.74, 0.80, 0.87, 0.93, 1.2])
bins["muonMom"]=array.array("d",[0.10, 0.18, 0.30, 0.48, 0.75, 1.14, 2.50])
bins["protonAngle"]=array.array("d",[-1.00, -0.50, 0.00, 0.27, 0.45, 0.62, 0.76, 0.86, 0.94, 1.00])
bins["muonAngle"]=array.array("d",[-1.00, -0.82, -0.66, -0.39, -0.16, 0.05, 0.25, 0.43, 0.59, 0.73, 0.83, 0.91, 1.00])
bins["thetamup"]=array.array("d",[0.00, 0.8, 1.2, 1.57, 1.94, 2.34, 3.14])
bins["costhetamup"]=array.array("d",[-1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0])

for cut in xrange(5):
  for cut2 in xrange(2):
    h_protonMom[cut,cut2] = ROOT.TH1D("h_protonMom_cut"+str(cut)+str(cut2), "h_protonMom cut "+str(cut)+str(cut2)+";proton momentum [GeV/c];Events",len(bins["protonMom"])-1,bins["protonMom"])
    h_protonAngle[cut,cut2] = ROOT.TH1D("h_protonAngle_cut"+str(cut)+str(cut2), "h_protonAngle cut "+str(cut)+str(cut2)+";proton cos(#theta);Events",len(bins["protonAngle"])-1,bins["protonAngle"])
    h_muonMom[cut,cut2] = ROOT.TH1D("h_muonMom_cut"+str(cut)+str(cut2), "h_muonMom cut "+str(cut)+str(cut2)+";muon momentum [GeV/c];Events",len(bins["muonMom"])-1,bins["muonMom"])
    h_muonAngle[cut,cut2] = ROOT.TH1D("h_muonAngle_cut"+str(cut)+str(cut2), "h_muonAngle cut "+str(cut)+str(cut2)+";muon cos(#theta);Events",len(bins["muonAngle"])-1,bins["muonAngle"])
    h_protonMom_old[cut,cut2] = ROOT.TH1D("h_protonMom_cut_old"+str(cut)+str(cut2), "h_protonMom cut "+str(cut)+str(cut2)+";momentum [GeV/c];Events",len(bins["protonMom"])-1,bins["protonMom"])
    h_protonMom_fullPS[cut,cut2] = ROOT.TH1D("h_protonMom_fullPS_cut"+str(cut)+str(cut2), "h_protonMom_fullPS cut "+str(cut)+str(cut2)+";proton momentum [GeV/c];Events",len(bins["protonMom"])-1,bins["protonMom"])
    h_protonAngle_fullPS[cut,cut2] = ROOT.TH1D("h_protonAngle_fullPS_cut"+str(cut)+str(cut2), "h_protonAngle_fullPS cut "+str(cut)+str(cut2)+";proton cos(#theta);Events",len(bins["protonAngle"])-1,bins["protonAngle"])
    h_muonMom_fullPS[cut,cut2] = ROOT.TH1D("h_muonMom_fullPS_cut"+str(cut)+str(cut2), "h_muonMom_fullPS cut "+str(cut)+str(cut2)+";muon momentum [GeV/c];Events",len(bins["muonMom"])-1,bins["muonMom"])
    h_muonAngle_fullPS[cut,cut2] = ROOT.TH1D("h_muonAngle_fullPS_cut"+str(cut)+str(cut2), "h_muonAngle_fullPS cut "+str(cut)+str(cut2)+";muon cos(#theta);Events",len(bins["muonAngle"])-1,bins["muonAngle"])
    ## 0 - all, 1, CC-inc, 2 - contained, 3 - PID
    ## cut2: 0 - event passed cuts, 1, correct particles selected

for mode in ["QE","MEC","RES","DIS","COH"]:
  h_protonMom[mode] = h_protonMom[0,0].Clone("h_protonMom_"+mode+"_precuts")
  h_protonMom[mode,"post"] = h_protonMom[0,0].Clone("h_protonMom_"+mode+"_postcuts")
  h_protonAngle[mode] = h_protonAngle[0,0].Clone("h_protonAngle_"+mode+"_precuts")
  h_protonAngle[mode,"post"] = h_protonAngle[0,0].Clone("h_protonAngle_"+mode+"_postcuts")
  h_muonMom[mode] = h_muonMom[0,0].Clone("h_muonMom_"+mode+"_precuts")
  h_muonMom[mode,"post"] = h_muonMom[0,0].Clone("h_muonMom_"+mode+"_postcuts")
  h_muonAngle[mode] = h_muonAngle[0,0].Clone("h_muonAngle_"+mode+"_precuts")
  h_muonAngle[mode,"post"] = h_muonAngle[0,0].Clone("h_muonAngle_"+mode+"_postcuts")

h_protonMom_leadingRight = ROOT.TH1D("h_protonMom_leadingRight", "h_protonMom_leadingRight;momentum [GeV/c];Events",20,0,1.5)
h_protonMom_leadingWrong = ROOT.TH1D("h_protonMom_leadingWrong", "h_protonMom_leadingWrong;momentum [GeV/c];Events",20,0,1.5)
h_protonAngle_leadingRight = ROOT.TH1D("h_protonAngle_leadingRight", "h_protonAngle_leadingRight;momentum [GeV/c];Events",20,-1,1)
h_protonAngle_leadingWrong = ROOT.TH1D("h_protonAngle_leadingWrong", "h_protonAngle_leadingWrong;momentum [GeV/c];Events",20,-1,1)
h_protonAngleMom2D_leadingRight = ROOT.TH2D("h_protonAngleMom2D_leadingRight", "h_protonAngleMom2D_leadingRight;cos(#theta);Momentum [GeV/c]",20,-1,1, 20, 0, 1.5)
h_protonAngleMom2D_leadingWrong = ROOT.TH2D("h_protonAngleMom2D_leadingWrong", "h_protonAngleMom2D_leadingWrong;cos(#theta);Momentum [GeV/c]",20,-1,1, 20, 0, 1.5)

h_selprotonMom_leadingRight = ROOT.TH1D("h_selprotonMom_leadingRight", "h_protonMom_leadingRight;momentum [GeV/c];Events",20,0,1.5)
h_selprotonMom_leadingWrong = ROOT.TH1D("h_selprotonMom_leadingWrong", "h_protonMom_leadingWrong;momentum [GeV/c];Events",20,0,1.5)
h_selprotonAngle_leadingRight = ROOT.TH1D("h_selprotonAngle_leadingRight", "h_protonAngle_leadingRight;momentum [GeV/c];Events",20,-1,1)
h_selprotonAngle_leadingWrong = ROOT.TH1D("h_selprotonAngle_leadingWrong", "h_protonAngle_leadingWrong;momentum [GeV/c];Events",20,-1,1)
h_selprotonAngleMom2D_leadingRight = ROOT.TH2D("h_selprotonAngleMom2D_leadingRight", "h_protonAngleMom2D_leadingRight;cos(#theta);Momentum [GeV/c]",20,-1,1, 20, 0, 1.5)
h_selprotonAngleMom2D_leadingWrong = ROOT.TH2D("h_selprotonAngleMom2D_leadingWrong", "h_protonAngleMom2D_leadingWrong;cos(#theta);Momentum [GeV/c]",20,-1,1, 20, 0, 1.5)

h_rightVsWrong_angle = ROOT.TH2D("h_rightVsWrong_angle","h_rightVsWrong_angle;leading cos(#theta);sub-leading cos(#theta)",20,-1,1, 20, -1, 1)
h_rightVsWrong_mom = ROOT.TH2D("h_rightVsWrong_mom","h_rightVsWrong_mom;leading momentum;sub-leading momentum",40,0.2,1.5, 40, 0.2, 1.5)



h_pionMom = ROOT.TH1D("h_pionMom","h_pionMom",20,0,0.5)
h_pionAngle =  ROOT.TH1D("h_pionAngle","h_pionAngle",20,-1,1)
h_pionMomAngle = ROOT.TH2D("h_pionMomAngle","h_pionMomAngle",20,0,0.5,20,-1,1)
h_pionLength = ROOT.TH1D("h_pionLength","h_pionLength",40,0,100)
h_pionLengthZ = ROOT.TH1D("h_pionLengthZ","h_pionLengthZ",40,0,100)
h_pionLengthMom= ROOT.TH2D("h_pionLengthMom","h_pionLengthMom",40,0,100,20,0,0.5)


h_muonMomAngle_pionIsMuCand = ROOT.TH2D("h_muonMomAngle_pionIsMuCand","h_muonMomAngle_pionIsMuCand",20,0,2, 20,-1,1)
h_pionLength_muonIsMuCand = ROOT.TH1D("h_pionLength_muonIsMuCand","h_pionLength_muonIsMuCand",40,0,100)
h_pionLengthZ_muonIsMuCand = ROOT.TH1D("h_pionLengthZ_muonIsMuCand","h_pionLengthZ_muonIsMuCand",40,0,100)
h_pionLengthMom_muonIsMuCand = ROOT.TH2D("h_pionLengthMom_muonIsMuCand","h_pionLengthMom_muonIsMuCand",40,0,100,20,0,0.5)
h_pionLengthZ_pionIsMuCand = ROOT.TH1D("h_pionLengthZ_pionIsMuCand","h_pionLengthZ_pionIsMuCand",40,0,100)


h_mupdg = ROOT.TH1D("h_mupdg","h_mupdg",6,-1,5)
h_mupdg.GetXaxis().SetBinLabel(1,"no match")
h_mupdg.GetXaxis().SetBinLabel(2,"muon")
h_mupdg.GetXaxis().SetBinLabel(3,"pion")
h_mupdg.GetXaxis().SetBinLabel(4,"proton")
h_mupdg.GetXaxis().SetBinLabel(5,"electron")
h_mupdg.GetXaxis().SetBinLabel(6,"other")

h_picount_ch_neut = ROOT.TH1D("h_picount_ch_neut","h_picount_ch_neut",4,0,4)
h_picount_ch_neut.GetXaxis().SetBinLabel(1,"total CC background")
h_picount_ch_neut.GetXaxis().SetBinLabel(2,"neutral pion")
h_picount_ch_neut.GetXaxis().SetBinLabel(3,"charged pion")
h_picount_ch_neut.GetXaxis().SetBinLabel(4,"both")

h_process_pionIsProtonCand = ROOT.TH1D("h_process_pionIsProtonCand","h_process_pionIsProtonCand",2,0,2)
h_process_pionIsProtonCand.GetXaxis().SetBinLabel(2,"Inelastic scatter")
h_process_pionIsProtonCand.GetXaxis().SetBinLabel(1,"stopped")

N_signal_cut=dict()
N_signal_cut_pluscos=dict()
N_cut=dict()
for thr in thresholds:
  for cut in xrange(10):
    N_signal_cut[thr,cut]=0
    N_signal_cut_pluscos[thr,cut]=0
    N_cut[thr,cut]=0

for thr in thresholds:
  N_signal[thr]=0
  N_signal_z[thr]=0
  N_selected[thr]=0
  N_selected_z[thr]=0
for event in tree["MC"]:
  i+=1
  if i%perc == 0 and i>perc-10 and i < dec:
    print int(100*float(i)/N_entries),"%"
  if i%dec == 0 and i>dec-10:
    print int(100*float(i)/N_entries),"%"
  if i==10:
    print "processed 10 events"
  if i==100:
    print "processed 100 events"
  if i==1000:
    print "processed 1000 events"
  if i==10000:
    print "processed 10000 events"
  if i==100000:
    print "processed 100000 events"
    
#  if i > dec:
#    break
    
  e = event.ubxsec_event_split
  #print dir(e)
  #exit()
#  if not e.is_signal:
#    continue
#  if not e.ngenie_protons:
#    continue
#  if (e.ngenie_pion0s>0 or e.ngenie_pipms>0):
#    continue
  proton_true_mom = getLeadingProtonEnergy(e)
  #proton_true_costheta = math.cos(getRecoProtonDir(e).Theta())
  proton_true_costheta = math.cos(getTrueProtonDir(e).Theta())
  muon_true_mom = getTrueValue(e, "muonMom")
  muon_true_costheta = getTrueValue(e, "muonAngle")
  h_all_true_mom.Fill(proton_true_mom)
  is_cosmic = IsCosmic(e, considerShowerTracks)
  is_cosmic_old = IsCosmic(e, False)
  particlesCorrect = muCandpCandTruth(e, considerShowerTracks)
  particlesCorrect_old = muCandpCandTruth(e, False)
  mode = modeToText(e.mode)
  for thr in thresholds:
    N_above[thr]=0
    N_above_z[thr]=0
    N_reco[thr]=0
    N_reco_z[thr]=0
    N_above[thr]=getProtonsAbove(e,thr)
    N_above_z[thr]=getProtonsAboveZ(e,thr)
    N_signal[thr]+=isSignal(e,thr,0)
    N_signal_z[thr]+=isSignal(e,thr,1)
    N_signal_cut[thr,0]+=isSignal(e,thr,0)*(not is_cosmic)
    N_signal_cut_pluscos[thr,0]+=isSignal(e,thr,0)
    N_cut[thr,0]+=1 # isSignal(e,thr,0)
  #print mode, isSignal(e,300,0)
  if isSignal(e,300,0):
    h_all_true_mom_thr.Fill(proton_true_mom)
    h_protonMom[mode].Fill(proton_true_mom,isSignal(e,300,0))
    h_protonAngle[mode].Fill(proton_true_costheta,isSignal(e,300,0))
    h_muonMom[mode].Fill(muon_true_mom,isSignal(e,300,0))
    h_muonAngle[mode].Fill(muon_true_costheta,isSignal(e,300,0))
  h_protonMom[0,0].Fill(proton_true_mom,isSignal(e,300,0))# and particlesCorrect and not is_cosmic)
  h_protonMom[0,1].Fill(proton_true_mom,isSignal(e,300,0) and particlesCorrect and not is_cosmic)
  h_protonAngle[0,0].Fill(proton_true_costheta,isSignal(e,300,0))# and particlesCorrect and not is_cosmic)
  h_protonAngle[0,1].Fill(proton_true_costheta,isSignal(e,300,0) and particlesCorrect and not is_cosmic)
  h_protonMom_old[0,0].Fill(proton_true_mom,isSignal(e,300,0))# and particlesCorrect and not is_cosmic)
  h_protonMom_old[0,1].Fill(proton_true_mom,isSignal(e,300,0) and particlesCorrect_old and not is_cosmic_old)
  h_muonMom[0,0].Fill(muon_true_mom,isSignal(e,300,0))# and particlesCorrect and not is_cosmic)
  h_muonMom[0,1].Fill(muon_true_mom,isSignal(e,300,0) and particlesCorrect and not is_cosmic)
  h_muonAngle[0,0].Fill(muon_true_costheta,isSignal(e,300,0))# and particlesCorrect and not is_cosmic)
  h_muonAngle[0,1].Fill(muon_true_costheta,isSignal(e,300,0) and particlesCorrect and not is_cosmic)
  h_protonMom_fullPS[0,0].Fill(proton_true_mom,isSignalOld(e,300,0))# and particlesCorrect and not is_cosmic)
  h_protonMom_fullPS[0,1].Fill(proton_true_mom,isSignalOld(e,300,0) and particlesCorrect and not is_cosmic)
  h_protonAngle_fullPS[0,0].Fill(proton_true_costheta,isSignalOld(e,300,0))# and particlesCorrect and not is_cosmic)
  h_protonAngle_fullPS[0,1].Fill(proton_true_costheta,isSignalOld(e,300,0) and particlesCorrect and not is_cosmic)
  h_muonMom_fullPS[0,0].Fill(muon_true_mom,isSignalOld(e,300,0))# and particlesCorrect and not is_cosmic)
  h_muonMom_fullPS[0,1].Fill(muon_true_mom,isSignalOld(e,300,0) and particlesCorrect and not is_cosmic)
  h_muonAngle_fullPS[0,0].Fill(muon_true_costheta,isSignalOld(e,300,0))# and particlesCorrect and not is_cosmic)
  h_muonAngle_fullPS[0,1].Fill(muon_true_costheta,isSignalOld(e,300,0) and particlesCorrect and not is_cosmic)
  #print "next event"
  if e.is_selected:
    pind = getTrueProtonIndex(e)
    g4pind=-1
    if pind>-1:
      for ig4 in xrange(len(e.geant_mcpar_px)): ## Get G4 proton index to check if it rescattered
        if e.geant_mcpar_pdgcode[ig4] != 2212:
          continue
        if math.fabs(e.geant_mcpar_energy[ig4] - e.genie_mcpar_energy[pind])<0.0001:
          g4pind = ig4
    #if g4pind>-1:
    #  if e.geant_mcpar_end_process[g4pind] == "protonInelastic":
    #    particlesCorrect=False  # These three lines are temporary
    N_signal_cut[300,1]+=isSignal(e,300,0)*(not is_cosmic) # count selected signal at CC-inclusive
    N_signal_cut_pluscos[300,1]+=isSignal(e,300,0) # count selected signal at CC-inclusive
    N_cut[300,1]+=1 # count all events at CC-inclusive
    h_protonMom[1,0].Fill(proton_true_mom,isSignal(e,300,0))
    h_protonMom[1,1].Fill(proton_true_mom,isSignal(e,300,0) and particlesCorrect and not is_cosmic)
    h_protonAngle[1,0].Fill(proton_true_costheta,isSignal(e,300,0))
    h_protonAngle[1,1].Fill(proton_true_costheta,isSignal(e,300,0) and particlesCorrect and not is_cosmic)
    h_protonMom_old[1,0].Fill(proton_true_mom,isSignal(e,300,0))
    h_protonMom_old[1,1].Fill(proton_true_mom,isSignal(e,300,0) and particlesCorrect_old and not is_cosmic_old)
    h_muonMom[1,0].Fill(muon_true_mom,isSignal(e,300,0))
    h_muonMom[1,1].Fill(muon_true_mom,isSignal(e,300,0) and particlesCorrect and not is_cosmic)
    h_muonAngle[1,0].Fill(muon_true_costheta,isSignal(e,300,0))
    h_muonAngle[1,1].Fill(muon_true_costheta,isSignal(e,300,0) and particlesCorrect and not is_cosmic)
    h_protonMom_fullPS[1,0].Fill(proton_true_mom,isSignalOld(e,300,0))
    h_protonMom_fullPS[1,1].Fill(proton_true_mom,isSignalOld(e,300,0) and particlesCorrect and not is_cosmic)
    h_protonAngle_fullPS[1,0].Fill(proton_true_costheta,isSignalOld(e,300,0))
    h_protonAngle_fullPS[1,1].Fill(proton_true_costheta,isSignalOld(e,300,0) and particlesCorrect and not is_cosmic)
    h_muonMom_fullPS[1,0].Fill(muon_true_mom,isSignalOld(e,300,0))
    h_muonMom_fullPS[1,1].Fill(muon_true_mom,isSignalOld(e,300,0) and particlesCorrect and not is_cosmic)
    h_muonAngle_fullPS[1,0].Fill(muon_true_costheta,isSignalOld(e,300,0))
    h_muonAngle_fullPS[1,1].Fill(muon_true_costheta,isSignalOld(e,300,0) and particlesCorrect and not is_cosmic)
    if not passNtracksCut(e, considerShowerTracks):
      continue
    N_signal_cut[300,2]+=isSignal(e,300,0)*(not is_cosmic) # count selected signal at 2-tracks
    N_signal_cut_pluscos[300,2]+=isSignal(e,300,0) # count selected signal at 2-tracks
    N_cut[300,2]+=1 # count all events at 2-tracks
    if not passContainmentCut(e, considerShowerTracks):
      continue
    N_signal_cut[300,3]+=isSignal(e,300,0)*(not is_cosmic) # count selected signal at protons contained
    N_signal_cut_pluscos[300,3]+=isSignal(e,300,0) # count selected signal at protons contained
    N_cut[300,3]+=1 # count all events at protons contained
    h_protonMom[2,0].Fill(proton_true_mom,isSignal(e,300,0))# and particlesCorrect and not is_cosmic)
    h_protonMom[2,1].Fill(proton_true_mom,isSignal(e,300,0) and particlesCorrect and not is_cosmic)
    h_protonAngle[2,0].Fill(proton_true_costheta,isSignal(e,300,0))# and particlesCorrect and not is_cosmic)
    h_protonAngle[2,1].Fill(proton_true_costheta,isSignal(e,300,0) and particlesCorrect and not is_cosmic)
    h_muonMom[2,0].Fill(muon_true_mom,isSignal(e,300,0))# and particlesCorrect and not is_cosmic)
    h_muonMom[2,1].Fill(muon_true_mom,isSignal(e,300,0) and particlesCorrect and not is_cosmic)
    h_muonAngle[2,0].Fill(muon_true_costheta,isSignal(e,300,0))# and particlesCorrect and not is_cosmic)
    h_muonAngle[2,1].Fill(muon_true_costheta,isSignal(e,300,0) and particlesCorrect and not is_cosmic)
    h_protonMom_fullPS[2,0].Fill(proton_true_mom,isSignalOld(e,300,0))# and particlesCorrect and not is_cosmic)
    h_protonMom_fullPS[2,1].Fill(proton_true_mom,isSignalOld(e,300,0) and particlesCorrect and not is_cosmic)
    h_protonAngle_fullPS[2,0].Fill(proton_true_costheta,isSignalOld(e,300,0))# and particlesCorrect and not is_cosmic)
    h_protonAngle_fullPS[2,1].Fill(proton_true_costheta,isSignalOld(e,300,0) and particlesCorrect and not is_cosmic)
    h_muonMom_fullPS[2,0].Fill(muon_true_mom,isSignalOld(e,300,0))# and particlesCorrect and not is_cosmic)
    h_muonMom_fullPS[2,1].Fill(muon_true_mom,isSignalOld(e,300,0) and particlesCorrect and not is_cosmic)
    h_muonAngle_fullPS[2,0].Fill(muon_true_costheta,isSignalOld(e,300,0))# and particlesCorrect and not is_cosmic)
    h_muonAngle_fullPS[2,1].Fill(muon_true_costheta,isSignalOld(e,300,0) and particlesCorrect and not is_cosmic)
    if passMinHitsCut(e, considerShowerTracks):
      N_signal_cut[300,4]+=isSignal(e,300,0)*(not is_cosmic) # count selected signal at 5 hit cut
      N_signal_cut_pluscos[300,4]+=isSignal(e,300,0) # count selected signal at 5 hit cut
      N_cut[300,4]+=1 # count selected event at 5 hit cut
      if passChi2Cut(e, considerShowerTracks):
        N_signal_cut[300,5]+=isSignal(e,300,0)*(not is_cosmic) # count selected signal at chi2 cut
        N_signal_cut_pluscos[300,5]+=isSignal(e,300,0) # count selected signal at chi2 cut
        N_cut[300,5]+=1 # count all selected events at chi2 cut
    if passChi2Cut(e, considerShowerTracks) and passMinHitsCut(e, considerShowerTracks) and passPhaseSpaceCuts(e, considerShowerTracks):
      N_signal_cut[300,6]+=isSignal(e,300,0)*(not is_cosmic) # count selected signal at phase space cuts
      N_signal_cut_pluscos[300,6]+=isSignal(e,300,0) # count selected signal at phase space cuts
      N_cut[300,6]+=1 # count all selected events at phase space cuts
      h_protonMom[3,0].Fill(proton_true_mom,isSignal(e,300,0))# and particlesCorrect and not is_cosmic)
      h_protonMom[3,1].Fill(proton_true_mom,isSignal(e,300,0) and particlesCorrect and not is_cosmic)
      h_protonAngle[3,0].Fill(proton_true_costheta,isSignal(e,300,0))# and particlesCorrect and not is_cosmic)
      h_protonAngle[3,1].Fill(proton_true_costheta,isSignal(e,300,0) and particlesCorrect and not is_cosmic)
      h_muonMom[3,0].Fill(muon_true_mom,isSignal(e,300,0))# and particlesCorrect and not is_cosmic)
      h_muonMom[3,1].Fill(muon_true_mom,isSignal(e,300,0) and particlesCorrect and not is_cosmic)
      h_muonAngle[3,0].Fill(muon_true_costheta,isSignal(e,300,0))# and particlesCorrect and not is_cosmic)
      h_muonAngle[3,1].Fill(muon_true_costheta,isSignal(e,300,0) and particlesCorrect and not is_cosmic)
      h_protonMom_fullPS[3,0].Fill(proton_true_mom,isSignalOld(e,300,0))# and particlesCorrect and not is_cosmic)
      h_protonMom_fullPS[3,1].Fill(proton_true_mom,isSignalOld(e,300,0) and particlesCorrect and not is_cosmic)
      h_protonAngle_fullPS[3,0].Fill(proton_true_costheta,isSignalOld(e,300,0))# and particlesCorrect and not is_cosmic)
      h_protonAngle_fullPS[3,1].Fill(proton_true_costheta,isSignalOld(e,300,0) and particlesCorrect and not is_cosmic)
      h_muonMom_fullPS[3,0].Fill(muon_true_mom,isSignalOld(e,300,0))# and particlesCorrect and not is_cosmic)
      h_muonMom_fullPS[3,1].Fill(muon_true_mom,isSignalOld(e,300,0) and particlesCorrect and not is_cosmic)
      h_muonAngle_fullPS[3,0].Fill(muon_true_costheta,isSignalOld(e,300,0))# and particlesCorrect and not is_cosmic)
      h_muonAngle_fullPS[3,1].Fill(muon_true_costheta,isSignalOld(e,300,0) and particlesCorrect and not is_cosmic)
      # Mode variables
      if isSignal(e,300,0):
        #h_protonMom[mode,"post"].Fill(proton_true_mom,isSignal(e,300,0) and particlesCorrect and not is_cosmic)
        #h_protonAngle[mode,"post"].Fill(proton_true_costheta,isSignal(e,300,0) and particlesCorrect and not is_cosmic)
        #h_muonMom[mode,"post"].Fill(muon_true_mom,isSignal(e,300,0) and particlesCorrect and not is_cosmic)
        #h_muonAngle[mode,"post"].Fill(muon_true_costheta,isSignal(e,300,0) and particlesCorrect and not is_cosmic)
        h_protonMom[mode,"post"].Fill(proton_true_mom,isSignal(e,300,0) and particlesCorrect and not is_cosmic)
        h_protonAngle[mode,"post"].Fill(proton_true_costheta,isSignal(e,300,0) and particlesCorrect and not is_cosmic)
        h_muonMom[mode,"post"].Fill(muon_true_mom,isSignal(e,300,0) and particlesCorrect and not is_cosmic)
        h_muonAngle[mode,"post"].Fill(muon_true_costheta,isSignal(e,300,0) and particlesCorrect and not is_cosmic)
      #print "next selected event"
      h_sel_true_mom.Fill(proton_true_mom)
      
      pind_truth = getRecoProtonIndex(e, considerShowerTracks)
      true_E = e.genie_mcpar_energy[pind] - 0.938270329696
      reco_true_E = e.pfp_truth_KE[pind_truth]
      if math.fabs(true_E - reco_true_E) < 0.01:
        protonCorrect=True
      else:
        protonCorrect=False
      
      #protonCorrect = protonIsLeading(e, considerShowerTracks)
      
      reco_proton_true_costheta = e.pfp_truth_costheta[pind_truth]
      reco_proton_true_mom = e.pfp_truth_mom[pind_truth]
#      reco_proton_true_px = e.genie_mcpar_px[pind_truth]
#      reco_proton_true_py = e.genie_mcpar_py[pind_truth]
#      reco_proton_true_pz = e.genie_mcpar_pz[pind_truth]
#      reco_proton_true_mom = math.sqrt(reco_proton_true_px*reco_proton_true_px + reco_proton_true_py*reco_proton_true_py + reco_proton_true_pz*reco_proton_true_pz)
#      if (reco_proton_true_mom):
#        reco_proton_true_costheta = reco_proton_true_pz/reco_proton_true_mom
#      else:
#        reco_proton_true_costheta=-999

      h_protonMom_leadingRight.Fill(proton_true_mom, isSignal(e,300,0) and protonCorrect and particlesCorrect)
      h_protonMom_leadingWrong.Fill(proton_true_mom, isSignal(e,300,0) and particlesCorrect and not protonCorrect)
      h_protonAngle_leadingRight.Fill(proton_true_costheta, isSignal(e,300,0) and protonCorrect and particlesCorrect)
      h_protonAngle_leadingWrong.Fill(proton_true_costheta, isSignal(e,300,0) and particlesCorrect and not protonCorrect)
      h_protonAngleMom2D_leadingRight.Fill(proton_true_costheta, proton_true_mom, isSignal(e,300,0) and protonCorrect and particlesCorrect)
      h_protonAngleMom2D_leadingWrong.Fill(proton_true_costheta, proton_true_mom, isSignal(e,300,0) and particlesCorrect and not protonCorrect)

      h_selprotonMom_leadingRight.Fill(reco_proton_true_mom, isSignal(e,300,0) and protonCorrect and particlesCorrect)
      h_selprotonMom_leadingWrong.Fill(reco_proton_true_mom, isSignal(e,300,0) and particlesCorrect and not protonCorrect)
      h_selprotonAngle_leadingRight.Fill(reco_proton_true_costheta, isSignal(e,300,0) and protonCorrect and particlesCorrect)
      h_selprotonAngle_leadingWrong.Fill(reco_proton_true_costheta, isSignal(e,300,0) and particlesCorrect and not protonCorrect)
      h_selprotonAngleMom2D_leadingRight.Fill(reco_proton_true_costheta, proton_true_mom, isSignal(e,300,0) and protonCorrect and particlesCorrect)
      h_selprotonAngleMom2D_leadingWrong.Fill(reco_proton_true_costheta, proton_true_mom, isSignal(e,300,0) and particlesCorrect and not protonCorrect)

      h_rightVsWrong_angle.Fill(proton_true_costheta, reco_proton_true_costheta, isSignal(e,300,0) and particlesCorrect and not protonCorrect)
      h_rightVsWrong_mom.Fill(proton_true_mom, reco_proton_true_mom, isSignal(e,300,0) and particlesCorrect and not protonCorrect)

      if isSignal(e,300,0):
        h_sel_true_mom_thr.Fill(proton_true_mom)
      for thr in thresholds:
        N_reco[thr] = getRecoAbove(e,thr, considerShowerTracks)   
        N_reco_z[thr] = getRecoAboveZ(e,thr, considerShowerTracks)   
        N_selected[thr] += isSignal(e,thr,0)*(N_reco[thr]>0)
        N_selected_z[thr] += isSignal(e,thr,1)*(N_reco_z[thr]>0)
        #if thr==200:
          #print "thresh = ",thr
          #print "N_reco =", N_reco[thr]
          #print "N_true =", N_above[thr]
        h_protonSmearing[thr,"Z"].Fill(N_above_z[thr],N_reco_z[thr])
        h_protonSmearing[thr,"tot"].Fill(N_above[thr],N_reco[thr])

      if e.is_signal and not isSignal(e,300,0) and not IsCosmic(e,True): # CC, in FV, but not 0pNp
        #if e.ngenie_pion0s==0 and e.ngenie_pipms==0:
        #  print "Whaaaat HERE"
        #  if particlesCorrect:
        #    print "My code says we DID reco the right particles"
        #  else:
        #    print "My code says we DID NOT reco the right particles"
          #dumpEvent(e)
        h_picount_ch_neut.Fill(0)
        if e.ngenie_pion0s>0 and e.ngenie_pipms>0:
          h_picount_ch_neut.Fill(3)
        if e.ngenie_pion0s>0:
          h_picount_ch_neut.Fill(1)
        if e.ngenie_pipms>0:
          h_picount_ch_neut.Fill(2)
          p_pi = getLeadingPionMom(e)
          costheta_pi = getLeadingPionAngle(e)
          h_pionMom.Fill(p_pi)
          h_pionAngle.Fill(costheta_pi)
          h_pionMomAngle.Fill(p_pi, costheta_pi)
          geant_index = getLeadingPionIndexGEANT(e)
          lengthx = math.fabs(e.geant_mcpar_startx[geant_index] - e.geant_mcpar_endx[geant_index])
          lengthy = math.fabs(e.geant_mcpar_starty[geant_index] - e.geant_mcpar_endy[geant_index])
          lengthz = math.fabs(e.geant_mcpar_startz[geant_index] - e.geant_mcpar_endz[geant_index])
          process = e.geant_mcpar_end_process[geant_index]
          pi_length = math.sqrt(lengthx*lengthx + lengthy*lengthy + lengthz*lengthz )
          h_pionLength.Fill(pi_length)
          h_pionLengthZ.Fill(lengthz)
          h_pionLengthMom.Fill(pi_length, p_pi)
          i_mu=-1
          protonCandIsPi=False
          for pfp_i in xrange(len(e.pfp_reco_Mom_proton)):
            if e.pfp_reco_ismuoncandidate[pfp_i]:
              i_mu = pfp_i
            elif abs(e.pfp_truth_pdg[i_mu])==211:
              protonCandIsPi=True
          if i_mu>=0:
            muonpdg = e.pfp_truth_pdg[i_mu]
          else:
            muonpdg=-999
          h_mupdg.Fill(simpleParticleIndex(muonpdg))
          if math.fabs(muonpdg)==211:
            h_muonMomAngle_pionIsMuCand.Fill(muon_true_mom, muon_true_costheta)
            h_pionLengthZ_pionIsMuCand.Fill(lengthz)
          if muonpdg==13:
            h_pionLength_muonIsMuCand.Fill(pi_length)
            h_pionLengthZ_muonIsMuCand.Fill(lengthz)
            h_pionLengthMom_muonIsMuCand.Fill(pi_length, p_pi)
          if protonCandIsPi:
            if process=="pi+Inelastic" or process=="pi-Inelastic":
              h_process_pionIsProtonCand.Fill(1)
            else:
              h_process_pionIsProtonCand.Fill(0)
            


  #print event.ubxsec_event_split.event
  #print event.ubxsec_event_split.ccnc
  #print event.ubxsec_event_split.genie_mult
  #print event.ubxsec_event_split.genie_mcpar_pdgcode[0]
  #print event.ubxsec_event_split.nslices
  #e = event.ubxsec_event_split
  #print e.ngenie_muons
  #for i in xrange(e.nslices):
  #  print i,  e.slc_muoncandidate_theta[i]
# # e.Show(0)
  #print e._default_value
  #if e.is_selected:
  #  print e.n_pfp
  #  num_pfp = getattr(e,"num_pfp")
  #  print e.pfp_truth_pdg[0]
  #print event.ubxsec_event_split.num_pfp
  #print pfp_pdg[num_pfp]
  #print event.ubxsec_event_split.track_pfp_bragg_ratio[0]
canv = ROOT.TCanvas()


h_pionMom.Draw("e")
canv.SaveAs("figures/h_pionMom.eps")
h_pionAngle.Draw("e")
canv.SaveAs("figures/h_pionAngle.eps")
h_pionMomAngle.Draw("colz")
canv.SaveAs("figures/h_pionMomAngle.eps")
h_pionLength.Draw("e")
canv.SaveAs("figures/h_pionLength.eps")
h_pionLengthZ.Draw("e")
canv.SaveAs("figures/h_pionLengthZ.eps")
h_pionLengthMom.Draw("colz")
canv.SaveAs("figures/h_pionLengthMom.eps")

h_mupdg.Draw("hist")
canv.SaveAs("figures/h_muonPdg_pionInEvent.eps")

h_muonMomAngle_pionIsMuCand.Draw("col")
canv.SaveAs("figures/h_muonMomAngle_pionIsMuCand.eps")
h_pionLength_muonIsMuCand.Draw("e")
canv.SaveAs("figures/h_pionLength_muonIsMuCand.eps")
h_pionLengthZ_muonIsMuCand.Draw("e")
canv.SaveAs("figures/h_pionLengthZ_muonIsMuCand.eps")
h_pionLengthMom_muonIsMuCand.Draw("col")
canv.SaveAs("figures/h_pionLengthMom_muonIsMuCand.eps")


h_pionLengthZ.SetLineColor(1)
h_pionLengthZ_muonIsMuCand.SetLineColor(2)
h_pionLengthZ_pionIsMuCand.SetLineColor(3)
h_pionLengthZ.Draw("hist")
h_pionLengthZ_muonIsMuCand.Draw("hist same")
h_pionLengthZ_pionIsMuCand.Draw("hist same")
canv.SaveAs("figures/h_pionLengthZ_allVsmissed.eps")


h_picount_ch_neut.SetLineColor(1)
h_picount_ch_neut.Draw("hist")
canv.SaveAs("figures/h_picount_ch_neut.eps")

h_process_pionIsProtonCand.SetLineColor(1)
h_process_pionIsProtonCand.SetMinimum(0)
h_process_pionIsProtonCand.Draw("hist")
canv.SaveAs("figures/h_process_pionIsProtonCand")


print N_signal[0], N_selected[0]
print "300 MeV threshold - signal ",N_signal[300], " of which selected ", N_selected[300]

print "========================================================="
print "===================== Cut table ========================="
print "========================================================="
print "No cuts -",N_signal_cut[300,0], "of total", N_cut[300,0], "and including cosmics",N_signal_cut_pluscos[300,0]
print "CC-inclusive -",N_signal_cut[300,1], "of total", N_cut[300,1], "and including cosmics",N_signal_cut_pluscos[300,1]
print "Ntracks -",N_signal_cut[300,2], "of total", N_cut[300,2], "and including cosmics",N_signal_cut_pluscos[300,2]
print "Contained -",N_signal_cut[300,3], "of total", N_cut[300,3], "and including cosmics",N_signal_cut_pluscos[300,3]
print "5 hits -",N_signal_cut[300,4], "of total", N_cut[300,4], "and including cosmics",N_signal_cut_pluscos[300,4]
print "Chi2 -",N_signal_cut[300,5], "of total", N_cut[300,5], "and including cosmics",N_signal_cut_pluscos[300,5]
print "Phase space -",N_signal_cut[300,6], "of total", N_cut[300,6], "and including cosmics",N_signal_cut_pluscos[300,6]
print "========================================================="

for thr in thresholds:
  eff = N_selected[thr]/(float(N_signal[thr]))
  print "threshold",thr,"eff = ",eff
for thr in thresholds:
  eff = N_selected_z[thr]/(float(N_signal_z[thr]))
  print "threshold_z",thr,"eff = ",eff

h_all_true_mom.Draw("hist")
#raw_input()

h_sel_true_mom.Divide(h_all_true_mom)
h_sel_true_mom.Draw("hist")
canv.SaveAs("figures/efficiency_protonmom.eps")
#raw_input()

h_sel_true_mom_thr.Divide(h_all_true_mom_thr)
h_sel_true_mom_thr.Draw("hist")
canv.SaveAs("figures/efficiency_protonmom_thr.eps")
#raw_input()

###############################
## Proton momentum eff curves
effpmom=dict()
effpangle=dict()
effmumom=dict()
effmuangle=dict()
effpmom_fullPS=dict()
effpangle_fullPS=dict()
effmumom_fullPS=dict()
effmuangle_fullPS=dict()

for i in xrange(4):
  for j in xrange(2):
    effpmom[i,j]=ROOT.TEfficiency(h_protonMom[i,j], h_protonMom[0,0])
    effpmom[i,j].SetLineColor(i+1)
    effpangle[i,j]=ROOT.TEfficiency(h_protonAngle[i,j], h_protonAngle[0,0])
    effpangle[i,j].SetLineColor(i+1)
    effmumom[i,j]=ROOT.TEfficiency(h_muonMom[i,j], h_muonMom[0,0])
    effmumom[i,j].SetLineColor(i+1)
    effmuangle[i,j]=ROOT.TEfficiency(h_muonAngle[i,j], h_muonAngle[0,0])
    effmuangle[i,j].SetLineColor(i+1)
    effpmom_fullPS[i,j]=ROOT.TEfficiency(h_protonMom_fullPS[i,j], h_protonMom_fullPS[0,0])
    effpmom_fullPS[i,j].SetLineColor(i+1)
    effpangle_fullPS[i,j]=ROOT.TEfficiency(h_protonAngle_fullPS[i,j], h_protonAngle_fullPS[0,0])
    effpangle_fullPS[i,j].SetLineColor(i+1)
    effmumom_fullPS[i,j]=ROOT.TEfficiency(h_muonMom_fullPS[i,j], h_muonMom_fullPS[0,0])
    effmumom_fullPS[i,j].SetLineColor(i+1)
    effmuangle_fullPS[i,j]=ROOT.TEfficiency(h_muonAngle_fullPS[i,j], h_muonAngle_fullPS[0,0])
    effmuangle_fullPS[i,j].SetLineColor(i+1)
  effpmom[i,1].SetLineStyle(2)
  effmumom[i,1].SetLineStyle(2)
  effpangle[i,1].SetLineStyle(2)
  effmuangle[i,1].SetLineStyle(2)
  effpmom_fullPS[i,1].SetLineStyle(2)
  effmumom_fullPS[i,1].SetLineStyle(2)
  effpangle_fullPS[i,1].SetLineStyle(2)
  effmuangle_fullPS[i,1].SetLineStyle(2)

for mode in ["QE","MEC","RES","DIS"]:
  effpmom[mode] = ROOT.TEfficiency(h_protonMom[mode,"post"], h_protonMom[mode])
  effmumom[mode] = ROOT.TEfficiency(h_muonMom[mode,"post"], h_muonMom[mode])
  effpangle[mode] = ROOT.TEfficiency(h_protonAngle[mode,"post"], h_protonAngle[mode])
  effmuangle[mode] = ROOT.TEfficiency(h_muonAngle[mode,"post"], h_muonAngle[mode])

legvscut=ROOT.TLegend(0.55, 0.65, 0.89, 0.89)
legvscut.AddEntry(effpmom[0,0],"No cuts","l")
legvscut.AddEntry(effpmom[1,0],"CC-inclusive","l")
legvscut.AddEntry(effpmom[2,0],"Containment","l")
legvscut.AddEntry(effpmom[3,0],"PID","l")
legvscut.AddEntry(effpmom[0,1],"correct particles","l")


f_Kirby = ROOT.TFile("Tune1_Tune3_eff.root")
fout_Andy = ROOT.TFile("AndyEfficiency.root","recreate")

effpmom[0,0].SetTitle(";Proton momentum [GeV/c];Efficiency")
effpmom[0,0].Draw("A ")
ROOT.gPad.Update()
graph = effpmom[0,0].GetPaintedGraph()
graph.SetMinimum(0)
graph.SetMaximum(1.4)

for i in xrange(4):
  for j in xrange(2):
    if i==0 and j==0:
      continue
    effpmom[i,j].Draw("e same")
legvscut.Draw()
canv.SaveAs("figures/efficiency/efficiency_protonMom_vsCut.eps")
canv.SaveAs("figures/efficiency/efficiency_protonMom_vsCut.png")

graph.SetMaximum(0.1)
ROOT.gPad.Update()
canv.SaveAs("figures/efficiency/efficiency_protonMom_vsCut_zoom.eps")
canv.SaveAs("figures/efficiency/efficiency_protonMom_vsCut_zoom.png")

canv.Clear()
effpmom[3,1].SetTitle(";Proton cos(#theta);Efficiency")
effpmom[3,1].Draw("A e")
canv.SaveAs("figures/efficiency/efficiency_protonMom_final.eps")
canv.SaveAs("figures/efficiency/efficiency_protonMom_final.png")
h_tmp = f_Kirby.Get("h_eff_trkpmom_den_smear_Tune1")
g = effpmom[3,1].GetPaintedGraph()
g.SetMaximum(0.5)
ROOT.gPad.Update()
h_tmp.Draw("same")
canv.SaveAs("figures/efficiency/efficiency_protonMom_final_vsKirby.eps")
canv.SaveAs("figures/efficiency/efficiency_protonMom_final_vsKirby.png")
fout_Andy.cd()
effpmom[3,0].Write()
effpmom[3,1].Write()


## Proton cos(theta)
effpangle[0,0].SetTitle(";Proton cos(#theta);Efficiency")
effpangle[0,0].Draw("A")
ROOT.gPad.Update()
graph = effpangle[0,0].GetPaintedGraph()
graph.SetMinimum(0)
graph.SetMaximum(1.4)

for i in xrange(4):
  for j in xrange(2):
    if i==0 and j==0:
      continue
    effpangle[i,j].Draw("e same")
legvscut.Draw()
canv.SaveAs("figures/efficiency/efficiency_protonAngle_vsCut.eps")
canv.SaveAs("figures/efficiency/efficiency_protonAngle_vsCut.png")

canv.Clear()
effpangle[3,1].SetTitle(";Proton cos(#theta);Efficiency")
effpangle[3,1].Draw("A e")
canv.SaveAs("figures/efficiency/efficiency_protonAngle_final.eps")
canv.SaveAs("figures/efficiency/efficiency_protonAngle_final.png")
h_tmp = f_Kirby.Get("h_eff_trkpcostheta_den_smear_Tune1")
h_tmp.Draw("same")
g = effpangle[3,1].GetPaintedGraph()
g.SetMaximum(0.5)
ROOT.gPad.Update()
canv.SaveAs("figures/efficiency/efficiency_protonAngle_final_vsKirby.eps")
canv.SaveAs("figures/efficiency/efficiency_protonAngle_final_vsKirby.png")
fout_Andy.cd()
effpangle[3,0].Write()
effpangle[3,1].Write()

## Muon momentum
effmumom[0,0].SetTitle(";Muon momentum [GeV/c];Efficiency")
effmumom[0,0].Draw("A")
ROOT.gPad.Update()
graph = effmumom[0,0].GetPaintedGraph()
graph.SetMinimum(0)
graph.SetMaximum(1.4)

for i in xrange(4):
  for j in xrange(2):
    if i==0 and j==0:
      continue
    effmumom[i,j].Draw("e same")
legvscut.Draw()
canv.SaveAs("figures/efficiency/efficiency_muonMom_vsCut.eps")
canv.SaveAs("figures/efficiency/efficiency_muonMom_vsCut.png")

canv.Clear()
effmumom[3,1].SetTitle(";Muon momentum [GeV/c];Efficiency")
effmumom[3,1].Draw("A e")
canv.SaveAs("figures/efficiency/efficiency_muonMom_final.eps")
canv.SaveAs("figures/efficiency/efficiency_muonMom_final.png")
h_tmp = f_Kirby.Get("h_eff_trkmom_den_smear_Tune1")
g = effmumom[3,1].GetPaintedGraph()
g.SetMaximum(0.5)
ROOT.gPad.Update()
h_tmp.Draw("same")
canv.SaveAs("figures/efficiency/efficiency_muonMom_final_vsKirby.png")
canv.SaveAs("figures/efficiency/efficiency_muonMom_final_vsKirby.eps")
fout_Andy.cd()
effmumom[3,0].Write()
effmumom[3,1].Write()

## Muon cos(theta)
effmuangle[0,0].SetTitle(";Muon cos(#theta);Efficiency")
effmuangle[0,0].Draw("A")
ROOT.gPad.Update()
graph = effmuangle[0,0].GetPaintedGraph()
graph.SetMinimum(0)
graph.SetMaximum(1.4)

for i in xrange(4):
  for j in xrange(2):
    if i==0 and j==0:
      continue
    effmuangle[i,j].Draw("e same")
legvscut.Draw()
canv.SaveAs("figures/efficiency/efficiency_muonAngle_vsCut.eps")
canv.SaveAs("figures/efficiency/efficiency_muonAngle_vsCut.png")

canv.Clear()
effmuangle[3,1].SetTitle(";Muon cos(#theta);Efficiency")
effmuangle[3,1].Draw("A e")
canv.SaveAs("figures/efficiency/efficiency_muonAngle_final.eps")
canv.SaveAs("figures/efficiency/efficiency_muonAngle_final.png")
h_tmp = f_Kirby.Get("h_eff_trkcostheta_den_smear_Tune1")
g = effmuangle[3,1].GetPaintedGraph()
g.SetMaximum(0.5)
ROOT.gPad.Update()
h_tmp.Draw("same")
canv.SaveAs("figures/efficiency/efficiency_muonAngle_final_vsKirby.eps")
canv.SaveAs("figures/efficiency/efficiency_muonAngle_final_vsKirby.png")
fout_Andy.cd()
effmuangle[3,0].Write()
effmuangle[3,1].Write()

# Full-phase-space comparisons?
effpmom_fullPS[0,0].SetTitle(";Proton momentum [GeV/c];Efficiency")
effpmom_fullPS[0,0].Draw("A ")
ROOT.gPad.Update()
graph = effpmom_fullPS[0,0].GetPaintedGraph()
graph.SetMinimum(0)
graph.SetMaximum(1.4)

for i in xrange(4):
  for j in xrange(2):
    if i==0 and j==0:
      continue
    effpmom_fullPS[i,j].Draw("e same")
legvscut.Draw()
canv.SaveAs("figures/efficiency_protonMom_vsCut_fullPS.eps")
canv.SaveAs("figures/efficiency_protonMom_vsCut_fullPS.png")

graph.SetMaximum(0.1)
ROOT.gPad.Update()
canv.SaveAs("figures/efficiency_protonMom_vsCut_zoom_fullPS.eps")
canv.SaveAs("figures/efficiency_protonMom_vsCut_zoom_fullPS.png")

canv.Clear()
effpmom_fullPS[3,1].SetTitle(";Proton momentum [GeV/c];Efficiency")
effpmom_fullPS[3,1].Draw("A e")
canv.SaveAs("figures/efficiency_protonMom_final_fullPS.eps")
canv.SaveAs("figures/efficiency_protonMom_final_fullPS.png")

## Proton cos(theta)
effpangle_fullPS[0,0].SetTitle(";Proton cos(#theta);Efficiency")
effpangle_fullPS[0,0].Draw("A")
ROOT.gPad.Update()
graph = effpangle_fullPS[0,0].GetPaintedGraph()
graph.SetMinimum(0)
graph.SetMaximum(1.4)

for i in xrange(4):
  for j in xrange(2):
    if i==0 and j==0:
      continue
    effpangle_fullPS[i,j].Draw("e same")
legvscut.Draw()
canv.SaveAs("figures/efficiency_protonAngle_vsCut_fullPS.eps")
canv.SaveAs("figures/efficiency_protonAngle_vsCut_fullPS.png")

canv.Clear()
effpangle_fullPS[3,1].SetTitle(";Proton cos(#theta);Efficiency")
effpangle_fullPS[3,1].Draw("A e")
canv.SaveAs("figures/efficiency_protonAngle_final_fullPS.eps")
canv.SaveAs("figures/efficiency_protonAngle_final_fullPS.png")

## Muon momentum
effmumom_fullPS[0,0].SetTitle(";Muon momentum [GeV/c];Efficiency")
effmumom_fullPS[0,0].Draw("A")
ROOT.gPad.Update()
graph = effmumom_fullPS[0,0].GetPaintedGraph()
graph.SetMinimum(0)
graph.SetMaximum(1.4)

for i in xrange(4):
  for j in xrange(2):
    if i==0 and j==0:
      continue
    effmumom_fullPS[i,j].Draw("e same")
legvscut.Draw()
canv.SaveAs("figures/efficiency_muonMom_vsCut_fullPS.eps")
canv.SaveAs("figures/efficiency_muonMom_vsCut_fullPS.png")

canv.Clear()
effmumom_fullPS[3,1].SetTitle(";Muon momentum [GeV/c];Efficiency")
effmumom_fullPS[3,1].Draw("A e")
canv.SaveAs("figures/efficiency_muonMom_final_fullPS.eps")
canv.SaveAs("figures/efficiency_muonMom_final_fullPS.png")

## Muon cos(theta)
effmuangle_fullPS[0,0].SetTitle(";Muon cos(#theta);Efficiency")
effmuangle_fullPS[0,0].Draw("A")
ROOT.gPad.Update()
graph = effmuangle_fullPS[0,0].GetPaintedGraph()
graph.SetMinimum(0)
graph.SetMaximum(1.4)

for i in xrange(4):
  for j in xrange(2):
    if i==0 and j==0:
      continue
    effmuangle_fullPS[i,j].Draw("e same")
legvscut.Draw()
canv.SaveAs("figures/efficiency_muonAngle_vsCut_fullPS.eps")
canv.SaveAs("figures/efficiency_muonAngle_vsCut_fullPS.png")

canv.Clear()
effmuangle_fullPS[3,1].SetTitle(";Muon cos(#theta);Efficiency")
effmuangle_fullPS[3,1].Draw("A e")
canv.SaveAs("figures/efficiency_muonAngle_final_fullPS.eps")
canv.SaveAs("figures/efficiency_muonAngle_final_fullPS.png")


# Comparisons of full and restricted phase space efficiencies
print "Got to full vs restricted phase space plots..."
effpmom[3,1].Draw("A e")
effpmom_fullPS[3,1].SetLineColor(1)
effpmom_fullPS[3,1].Draw("same")
ROOT.gPad.Update()
graph = effpmom[3,1].GetPaintedGraph()
graph.SetMinimum(0)
canv.SaveAs("figures/efficiency/efficiency_protonMom_fullVsRestrictedPS.eps")
canv.SaveAs("figures/efficiency/efficiency_protonMom_fullVsRestrictedPS.png")

effpangle[3,1].Draw("A e")
effpangle_fullPS[3,1].SetLineColor(1)
effpangle_fullPS[3,1].Draw("same")
ROOT.gPad.Update()
graph = effpangle[3,1].GetPaintedGraph()
graph.SetMinimum(0)
canv.SaveAs("figures/efficiency/efficiency_protonAngle_fullVsRestrictedPS.eps")
canv.SaveAs("figures/efficiency/efficiency_protonAngle_fullVsRestrictedPS.png")

effmumom[3,1].Draw("A e")
effmumom_fullPS[3,1].SetLineColor(1)
effmumom_fullPS[3,1].Draw("same")
ROOT.gPad.Update()
graph = effmumom[3,1].GetPaintedGraph()
graph.SetMinimum(0)
canv.SaveAs("figures/efficiency/efficiency_muonMom_fullVsRestrictedPS.eps")
canv.SaveAs("figures/efficiency/efficiency_muonMom_fullVsRestrictedPS.png")

effmuangle[3,1].Draw("A e")
effmuangle_fullPS[3,1].SetLineColor(1)
effmuangle_fullPS[3,1].Draw("same")
ROOT.gPad.Update()
graph = effmuangle[3,1].GetPaintedGraph()
graph.SetMinimum(0)
canv.SaveAs("figures/efficiency/efficiency_muonAngle_fullVsRestrictedPS.eps")
canv.SaveAs("figures/efficiency/efficiency_muonAngle_fullVsRestrictedPS.png")


#  PID efficiency only
PIDeff_good = ROOT.TEfficiency(h_protonMom[3,1],h_protonMom[2,1])
PIDeff_all = ROOT.TEfficiency(h_protonMom[3,0],h_protonMom[2,0])
tmp3 = h_protonMom[3,0].Clone()
tmp2 = h_protonMom[2,0].Clone()
tmp3.Add(h_protonMom[3,1], -1)
tmp2.Add(h_protonMom[2,1], -1)
PIDeff_bad = ROOT.TEfficiency(tmp3, tmp2)
PIDeff_all.SetLineColor(4)
PIDeff_good.SetLineColor(4)
PIDeff_good.SetLineStyle(2)
PIDeff_bad.SetLineColor(1)
PIDeff_all.Draw("A")
ROOT.gPad.Update()
g=PIDeff_all.GetPaintedGraph()
g.SetMinimum(0)
g.SetMaximum(1.2)
PIDeff_good.Draw("e same")
PIDeff_bad.Draw("e same")
canv.SaveAs("figures/efficiency_PIDcut_goodbad.eps")
canv.SaveAs("figures/efficiency_PIDcut_goodbad.png")

###########################
## Proton eff, tracks vs showers, after CCinc
###########################
canv.Clear()
eff_showerTracks = ROOT.TEfficiency(h_protonMom[0,1], h_protonMom[0,0])
eff_onlyTracks = ROOT.TEfficiency(h_protonMom_old[0,1], h_protonMom_old[0,0])

eff_showerTracks.SetLineColor(2)
eff_onlyTracks.SetLineColor(1)
eff_showerTracks.Draw("A")
eff_onlyTracks.Draw("e same")
ROOT.gPad.Update()
g=eff_showerTracks.GetPaintedGraph()
g.SetMinimum(0)
g.SetMaximum(1.2)
legTrSh = ROOT.TLegend(0.7, 0.7, 0.89, 0.89)
legTrSh.SetLineColor(0)
legTrSh.SetFillColor(0)
legTrSh.AddEntry(eff_onlyTracks,"only tracks","l")
legTrSh.AddEntry(eff_showerTracks,"showers and tracks","l")
legTrSh.Draw()
canv.SaveAs("figures/efficiency_trackshower_comparison.eps")
canv.SaveAs("figures/efficiency_trackshower_comparison.png")

###########################
## Curves broken up by mode
###########################
legModes=ROOT.TLegend(0.7, 0.7, 0.89, 0.89)
legModes.SetLineColor(0)
legModes.SetFillColor(0)
i=1
for mode in ["QE","MEC","RES","DIS"]:
  i+=1
  effpmom[mode].SetLineColor(i)
  effpangle[mode].SetLineColor(i)
  effmumom[mode].SetLineColor(i)
  effmuangle[mode].SetLineColor(i)
  effpmom[mode].SetLineWidth(2)
  effpangle[mode].SetLineWidth(2)
  effmumom[mode].SetLineWidth(2)
  effmuangle[mode].SetLineWidth(2)
  legModes.AddEntry(effpmom[mode],mode,"l")

canv.Clear()
i=-1
for mode in ["QE","MEC","RES","DIS"]:
  i+=1
  if not i:
    effpmom[mode].Draw("A")
    ROOT.gPad.Update()
    effpmom[mode].GetPaintedGraph().SetMaximum(0.6)
    effpmom[mode].GetPaintedGraph().SetMinimum(0.0)
  else:
    effpmom[mode].Draw("e same")
legModes.Draw()
canv.SaveAs("figures/efficiency/efficiency_protonMom_modes.eps")
canv.SaveAs("figures/efficiency/efficiency_protonMom_modes.png")

canv.Clear()
i=-1
for mode in ["QE","MEC","RES","DIS"]:
  i+=1
  if not i:
    effpangle[mode].Draw("A")
    ROOT.gPad.Update()
    effpangle[mode].GetPaintedGraph().SetMaximum(0.6)
    effpangle[mode].GetPaintedGraph().SetMinimum(0.0)
  else:
    effpangle[mode].Draw("e same")
legModes.Draw()
canv.SaveAs("figures/efficiency/efficiency_protonAngle_modes.eps")
canv.SaveAs("figures/efficiency/efficiency_protonAngle_modes.png")

canv.Clear()
i=-1
for mode in ["QE","MEC","RES","DIS"]:
  i+=1
  if not i:
    effmumom[mode].Draw("A")
    ROOT.gPad.Update()
    effmumom[mode].GetPaintedGraph().SetMaximum(0.6)
    effmumom[mode].GetPaintedGraph().SetMinimum(0.0)
  else:
    effmumom[mode].Draw("e same")
legModes.Draw()
canv.SaveAs("figures/efficiency/efficiency_muonMom_modes.eps")
canv.SaveAs("figures/efficiency/efficiency_muonMom_modes.png")

canv.Clear()
i=-1
for mode in ["QE","MEC","RES","DIS"]:
  i+=1
  if not i:
    effmuangle[mode].Draw("A")
    ROOT.gPad.Update()
    effmuangle[mode].GetPaintedGraph().SetMaximum(0.6)
    effmuangle[mode].GetPaintedGraph().SetMinimum(0.0)
  else:
    effmuangle[mode].Draw("e same")
legModes.Draw()
canv.SaveAs("figures/efficiency/efficiency_muonAngle_modes.eps")
canv.SaveAs("figures/efficiency/efficiency_muonAngle_modes.png")

###########################
## Proton smearing
###########################
canv.Clear()
canv.Divide(2)
for thr in thresholds:
  canv.cd(1)
  h_protonSmearing[thr,"Z"].Draw("colz")
  h_protonSmearing[thr,"Z"].Draw("text same")
  canv.cd(2)
  h_protonSmearing[thr,"tot"].Draw("colz")
  h_protonSmearing[thr,"tot"].Draw("text same")
  canv.SaveAs("figures/protonSmearing/"+str(thr)+"thresh_smearing.eps")
  canv.SaveAs("figures/protonSmearing/"+str(thr)+"thresh_smearing.png")
  canv.SaveAs("figures/protonSmearing/"+str(thr)+"thresh_smearing.pdf")


###########################
## Proton swapping?
###########################
print h_protonMom_leadingWrong.Integral(), "protons wrong, and", h_protonMom_leadingRight.Integral(), "protons correct"
print "fraction = ",h_protonMom_leadingWrong.Integral()/(h_protonMom_leadingRight.Integral() + h_protonMom_leadingWrong.Integral())

canv.Clear()
h_protonMom_leadingRight.Sumw2()
h_protonMom_leadingWrong.Sumw2()
h_protonMom_leadingRight.SetLineColor(1)
h_protonMom_leadingWrong.SetLineColor(2)
h_protonMom_leadingRight.Draw("e")
h_protonMom_leadingWrong.Draw("e same")
canv.SaveAs("figures/protonRightOrWrong_vs_momentum.eps")
canv.SaveAs("figures/protonRightOrWrong_vs_momentum.pdf")

canv.Clear()
tmp = h_protonMom_leadingRight.Clone()
tmp.Add(h_protonMom_leadingWrong)
tmp2 = h_protonMom_leadingWrong.Clone()
tmp2.Divide(tmp)
tmp2.SetMinimum(0)
tmp2.SetMaximum(0.2)
tmp2.SetLineColor(1)
tmp2.Draw("e")
canv.SaveAs("figures/protonConfusion_vs_momentum.eps")
canv.SaveAs("figures/protonConfusion_vs_momentum.pdf")

###########################
## Proton swapping - angle
###########################
canv.Clear()
h_protonAngle_leadingRight.Sumw2()
h_protonAngle_leadingWrong.Sumw2()
h_protonAngle_leadingRight.SetLineColor(1)
h_protonAngle_leadingWrong.SetLineColor(2)
h_protonAngle_leadingRight.Draw("e")
h_protonAngle_leadingWrong.Draw("e same")
canv.SaveAs("figures/protonRightOrWrong_vs_costheta.eps")
canv.SaveAs("figures/protonRightOrWrong_vs_costheta.pdf")

canv.Clear()
tmp = h_protonAngle_leadingRight.Clone()
tmp.Add(h_protonAngle_leadingWrong)
tmp2 = h_protonAngle_leadingWrong.Clone()
tmp2.Divide(tmp)
tmp2.SetMinimum(0)
tmp2.SetMaximum(0.2)
tmp2.SetLineColor(1)
tmp2.Draw("e")
canv.SaveAs("figures/protonConfusion_vs_costheta.eps")
canv.SaveAs("figures/protonConfusion_vs_costheta.pdf")


###########################
## Proton swapping - 2D
###########################
#canv.Clear()
#h_protonAngleMom2D_leadingRight.Sumw2()
#h_protonAngleMom2D_leadingWrong.Sumw2()
#h_protonAngleMom2D_leadingRight.SetLineColor(1)
#h_protonAngleMom2D_leadingWrong.SetLineColor(2)
#h_protonAngleMom2D_leadingRight.Draw("e")
#h_protonAngleMom2D_leadingWrong.Draw("e same")
#canv.SaveAs("figures/protonRightOrWrong_vs_costheta.eps")
#canv.SaveAs("figures/protonRightOrWrong_vs_costheta.pdf")

canv.Clear()
tmp = h_protonAngleMom2D_leadingRight.Clone()
tmp.Add(h_protonAngleMom2D_leadingWrong)
tmp2 = h_protonAngleMom2D_leadingWrong.Clone()
tmp2.Divide(tmp)
tmp2.SetMinimum(0)
tmp2.SetMaximum(0.5)
tmp2.SetTitle("proton confusion rate")
tmp2.SetLineColor(1)
tmp2.Draw("colz")
canv.SaveAs("figures/protonConfusion_2D_costheta_mom.eps")
canv.SaveAs("figures/protonConfusion_2D_costheta_mom.pdf")



#################################################################################
## Proton swapping - this time, what was the selected proton?
#################################################################################
canv.Clear()
h_selprotonMom_leadingRight.Sumw2()
h_selprotonMom_leadingWrong.Sumw2()
print "selprotonMom correct underflow = ", h_selprotonMom_leadingRight.GetBinContent(0)
print "selprotonMom correct overflow = ", h_selprotonMom_leadingRight.GetBinContent(h_selprotonMom_leadingRight.GetNbinsX()+1)

canv.Clear()
tmp = h_selprotonMom_leadingRight.Clone()
tmp.Add(h_selprotonMom_leadingWrong)
tmp2 = h_selprotonMom_leadingWrong.Clone()
tmp2.Divide(tmp)
tmp2.SetMinimum(0)
tmp2.SetMaximum(0.2)
tmp2.SetLineColor(1)
tmp2.Draw("e")
canv.SaveAs("figures/protonConfusion_vs_selmomentum.eps")
canv.SaveAs("figures/protonConfusion_vs_selmomentum.pdf")

###########################
## Proton swapping - angle
###########################
canv.Clear()
h_selprotonAngle_leadingRight.Sumw2()
h_selprotonAngle_leadingWrong.Sumw2()
print "selprotonAngle correct underflow = ", h_selprotonAngle_leadingRight.GetBinContent(0)
print "selprotonAngle correct overflow = ", h_selprotonAngle_leadingRight.GetBinContent(h_selprotonAngle_leadingRight.GetNbinsX()+1)

canv.Clear()
tmp = h_selprotonAngle_leadingRight.Clone()
tmp.Add(h_selprotonAngle_leadingWrong)
tmp2 = h_selprotonAngle_leadingWrong.Clone()
tmp2.Divide(tmp)
tmp2.SetMinimum(0)
tmp2.SetMaximum(0.2)
tmp2.SetLineColor(1)
tmp2.Draw("e")
canv.SaveAs("figures/protonConfusion_vs_selcostheta.eps")
canv.SaveAs("figures/protonConfusion_vs_selcostheta.pdf")


###########################
## Proton swapping - 2D
###########################
#canv.Clear()
h_selprotonAngleMom2D_leadingRight.Sumw2()
h_selprotonAngleMom2D_leadingWrong.Sumw2()
#h_protonAngleMom2D_leadingRight.SetLineColor(1)
#h_protonAngleMom2D_leadingWrong.SetLineColor(2)
#h_protonAngleMom2D_leadingRight.Draw("e")
#h_protonAngleMom2D_leadingWrong.Draw("e same")
#canv.SaveAs("figures/protonRightOrWrong_vs_costheta.eps")
#canv.SaveAs("figures/protonRightOrWrong_vs_costheta.pdf")

canv.Clear()
tmp = h_selprotonAngleMom2D_leadingRight.Clone()
tmp.Add(h_selprotonAngleMom2D_leadingWrong)
tmp2 = h_selprotonAngleMom2D_leadingWrong.Clone()
tmp2.Divide(tmp)
tmp2.SetMinimum(0)
tmp2.SetMaximum(0.5)
tmp2.SetTitle("proton confusion rate")
tmp2.SetLineColor(1)
tmp2.Draw("colz")
canv.SaveAs("figures/protonConfusion_2D_selcostheta_selmom.eps")
canv.SaveAs("figures/protonConfusion_2D_selcostheta_selmom.pdf")



canv.Clear()
h_rightVsWrong_angle.Draw("colz")
canv.SaveAs("figures/protonConfusion_trueThetaVsSelTheta.eps")
canv.SaveAs("figures/protonConfusion_trueThetaVsSelTheta.pdf")


canv.Clear()
h_rightVsWrong_mom.Draw("colz")
canv.SaveAs("figures/protonConfusion_trueMomVsSelMom.eps")
canv.SaveAs("figures/protonConfusion_trueMomVsSelMom.pdf")



