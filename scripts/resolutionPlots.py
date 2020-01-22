import ROOT
import math
import array

#TPaveStats *s = (TPaveStats*) gPad->GetPrimitive("stats");
#s->SetTextSize(0.025);
#s->SetX1NDC(0.55);
#s-SetX2NDC(0.75); 
#s->SetY1NDC(0.701271);
#s->SetY2NDC(0.800847);



def setBins(bins):
#  bins["protonMom"]=[0.30, 0.36, 0.41, 0.44, 0.49, 0.53, 0.56, 0.59, 0.63, 0.73, 0.81,1.27, 1.50] # original
  bins["protonMom"]=[0.30, 0.41, 0.495, 0.56, 0.62, 0.68, 0.74, 0.80, 0.87, 0.93, 1.50]
#  bins["muonMom"]=[0.00, 0.18, 0.30, 0.45, 0.77, 1.28, 2.50] # original
  bins["muonMom"]=[0.00, 0.18, 0.30, 0.48, 0.75, 1.14, 2.50]

  bins["protonAngle"]=[-1.00, -0.50, 0.00, 0.27, 0.45, 0.62, 0.76, 0.86, 0.94, 1.00]
  bins["muonAngle"]=[-1.00, -0.82, -0.66, -0.39, -0.16, 0.05, 0.25, 0.43, 0.59, 0.73, 0.83, 0.91, 1.00]
#  -1.00, -0.82, -0.66, -0.39, -0.16, 0.05, 0.25, 0.43, 0.59, 0.73, 0.83, 0.91, 1.00 # from Maker?
  bins["thetamup"]=[0.00, 0.8, 1.2, 1.57, 1.94, 2.34, 3.14]
  bins["costhetamup"]=[-1.0, -0.7, -0.35, 0.0, 0.35, 0.7, 1.0]
  return

def getRecoValue(e, var):
  val = -9999
  if var=="protonMom":
    val = getLeadingRecoMom(e)
  elif var == "muonMom":
    for i in xrange(len(e.pfp_reco_theta)):
      if e.pfp_reco_ismuoncandidate[i]:
        if inCV(e.pfp_reco_endx[i], e.pfp_reco_endy[i], e.pfp_reco_endz[i]):
          val = e.pfp_reco_Mom_muon[i]
        else:
          val = e.pfp_reco_Mom_MCS[i]
  elif var == "protonAngle":
#    val = getLeadingRecoAngle(e)
    val = math.cos(getRecoProtonDir(e).Theta())
  elif var == "muonAngle":
    for i in xrange(len(e.pfp_reco_theta)):
      if e.pfp_reco_ismuoncandidate[i]:
        val = math.cos(e.pfp_reco_theta[i])
  elif var == "thetamup":
    val = getRecoMuonDir(e).Angle(getRecoProtonDir(e))
  elif var == "costhetamup":
    val = getRecoMuonDir(e).Dot(getRecoProtonDir(e))
  elif var == "protonStart":
    pind = getRecoProtonIndex(e)
    vec = ROOT.TVector3(-999,-999,-999)
    vec.SetXYZ(e.pfp_reco_startx[pind], e.pfp_reco_starty[pind], e.pfp_reco_startz[pind])
    return vec
  elif var == "protonEnd":
    pind = getRecoProtonIndex(e)
    vec = ROOT.TVector3(-999,-999,-999)
    vec.SetXYZ(e.pfp_reco_endx[pind], e.pfp_reco_endy[pind], e.pfp_reco_endz[pind])
    return vec
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
    val=math.cos(getTrueProtonDir(e).Theta())
  elif var == "muonAngle":
    val = e.lep_costheta
  elif var == "thetamup":
    val = getTrueMuonDir(e).Angle(getTrueProtonDir(e))
  elif var == "costhetamup":
    val = getTrueMuonDir(e).Dot(getTrueProtonDir(e))
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

def getTrueProtonDir(e):
  pmax=-1
  pind=-1
  vec = ROOT.TVector3(0,0,0)
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
    return vec
  if pmax<0.0000000000000001:
    print "can't find leading proton true angle - momentum==0"
    return vec
  vec.SetXYZ(e.genie_mcpar_px[pind], e.genie_mcpar_py[pind], e.genie_mcpar_pz[pind])
  vec*=1./vec.Mag()
  return vec

def leadingProtonStopped(e, g):
  #print "finding G4 leading proton"
  pmax=-1
  pind=-1
  g4pind=-1
  for i in xrange(len(e.genie_mcpar_pdgcode)):
    if e.genie_mcpar_pdgcode[i] == 2212:
      px = e.genie_mcpar_px[i]
      py = e.genie_mcpar_py[i]
      pz = e.genie_mcpar_pz[i]
      p = math.sqrt(px*px + py*py + pz*pz)
      if p > pmax:
        pmax = p
        pind=i
  #print "genie leading proton index", pind
  #print "genie leading proton energy", e.genie_mcpar_energy[pind]
  if pind<0:
    print "can't find leading proton"
    return False
  if pmax<0.0000000000000001:
    print "can't find leading proton"
    return False
  for i in xrange(len(e.geant_mcpar_px)):
    #print "i", i, "energy = ",e.geant_mcpar_energy[i]
    if e.geant_mcpar_pdgcode[i] != 2212:
      continue
    if math.fabs(e.geant_mcpar_energy[i] - e.genie_mcpar_energy[pind])<0.0001:
      g4pind = i
    #  print "=================="
    #  print "G4:GENIE..."
    #  print "px", e.geant_mcpar_px[i], e.genie_mcpar_px[pind] 
    #  print "py", e.geant_mcpar_py[i], e.genie_mcpar_py[pind] 
    #  print "pz", e.geant_mcpar_pz[i], e.genie_mcpar_pz[pind] 
    #  print "=================="
      break
  if g4pind<0:
    print "can't find leading proton in g4"
    return False
  if e.geant_mcpar_end_process[g4pind] == "protonInelastic":
    return False
  startx = e.geant_mcpar_startx[g4pind]
  starty = e.geant_mcpar_starty[g4pind]
  startz = e.geant_mcpar_startz[g4pind]
  px = e.geant_mcpar_px[g4pind]
  py = e.geant_mcpar_py[g4pind]
  pz = e.geant_mcpar_pz[g4pind]
  endPos = calcEndPoint(startx, starty, startz, px, py, pz, g)
  if not inCV(endPos.X(), endPos.Y(), endPos.Z()):
    return False
  return True

def allProtonsStopped(e, g):  # Check that all protons stopped!
  for i in xrange(len(e.genie_mcpar_pdgcode)):
    if e.genie_mcpar_pdgcode[i] != 2212:
      continue
    gpx = e.genie_mcpar_px[i]
    gpy = e.genie_mcpar_py[i]
    gpz = e.genie_mcpar_pz[i]
    genie_p = math.sqrt(gpx*gpx + gpy*gpy + gpz*gpz)
    if genie_p==0: # GENIE proton has no momentum, so it can't be not stopping...
      continue
    g4pind=-1
    for j in xrange(len(e.geant_mcpar_px)): # Find the geant ID
      #print "i", i, "energy = ",e.geant_mcpar_energy[i]
      if e.geant_mcpar_pdgcode[j] != 2212:
        continue
      if math.fabs(e.geant_mcpar_energy[j] - e.genie_mcpar_energy[i])<0.0001:
        g4pind = j # G4 ID is this one
        break
    if g4pind<0:
      print "Can't find the G4 track for that GENIE proton"
      continue # Not sure, for now let's skip the GENIE protons that we can't find the G4 one for
    if e.geant_mcpar_end_process[g4pind] == "protonInelastic":
      return False
    startx = e.geant_mcpar_startx[g4pind]
    starty = e.geant_mcpar_starty[g4pind]
    startz = e.geant_mcpar_startz[g4pind]
    px = e.geant_mcpar_px[g4pind]
    py = e.geant_mcpar_py[g4pind]
    pz = e.geant_mcpar_pz[g4pind]
    endPos = calcEndPoint(startx, starty, startz, px, py, pz, g)
    if not inCV(endPos.X(), endPos.Y(), endPos.Z()):
      return False
  return True



def calcEndPoint(startx, starty, startz, px, py, pz, g):
  p = math.sqrt(px*px + py*py + pz*pz)
  if p==0:
    vec = ROOT.TVector3(startx, starty, startz)
    print "proton g4 momentum=0!"
    print "start point", startx, starty, startz
    return vec
  length = getLength(p, g)
  xlength = length*(px/p)
  ylength = length*(py/p)
  zlength = length*(pz/p)
  endX = startx + xlength
  endY = starty + ylength
  endZ = startz + zlength
  vec = ROOT.TVector3(endX, endY, endZ)
  return vec

def getLength(p, g):
  length = g.Eval(p)
  #print "momentum, length =", p, length
  return length

def getTrueMuonDir(e):
  vec = ROOT.TVector3(0,0,0)
  if math.fabs(e.lep_costheta)>1:
    print e.lep_costheta
    return vec
  if (e.lep_phi)<-7: # can't be less than 2pi
    print e.lep_phi
    return vec
  vec.SetMagThetaPhi(1, math.acos(e.lep_costheta), e.lep_phi)
  return vec

def getRecoProtonIndex(e):
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
    print "can't find leading proton index"
  return pind

def getMuonIndex(e):
  for i in xrange(len(e.pfp_reco_Mom_proton)):
    if e.pfp_reco_ismuoncandidate[i]:
      return i
  return -1
  

def getLeadingRecoMom(e):
  pmax=-1
  for i in xrange(len(e.pfp_reco_Mom_proton)):
    if e.pfp_reco_ismuoncandidate[i]:
      continue
    p = e.pfp_reco_Mom_proton[i]
    if p > pmax:
      pmax = p
  return pmax

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
    print "can't find leading proton reco angle"
    return vec
  vec.SetMagThetaPhi(1, e.pfp_reco_theta[pind], e.pfp_reco_phi[pind])
  return vec # math.cos(e.pfp_reco_theta[pind])

def getRecoMuonDir(e):
  vec = ROOT.TVector3(0,0,0)
  for i in xrange(len(e.pfp_reco_ismuoncandidate)):
    if e.pfp_reco_ismuoncandidate[i]:
      vec.SetMagThetaPhi(1, e.pfp_reco_theta[i],e.pfp_reco_phi[i])
      return vec
  print "No Muon Candidate Found!!!"
  return vec

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
  ntracks = sum(1 for x in xrange(len(e.pfp_reco_istrack)) if bool(e.pfp_reco_istrack[x] or e.pfp_reco_isshower[x]))
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

def passThreshCut(e):
  pind = getRecoProtonIndex(e)
  if e.pfp_reco_Mom_proton[pind]<0.3:
    return False
  return True

def IsCosmic(e):
  for i in xrange(len(e.pfp_reco_chi2_proton)):
    if not (e.pfp_reco_isshower[i] or e.pfp_reco_istrack[i]):
      continue
    if e.pfp_truth_origin[i]!=1:
      return True
  return False


ROOT.gStyle.SetLabelSize(0.045,"X")
ROOT.gStyle.SetLabelFont(62,"X")
ROOT.gStyle.SetTitleSize(0.045,"X")
ROOT.gStyle.SetTitleFont(62,"X")
ROOT.gStyle.SetTitleOffset(0.85,"X")
ROOT.gStyle.SetLabelSize(0.045,"Y")
ROOT.gStyle.SetLabelFont(62,"Y")
ROOT.gStyle.SetTitleSize(0.045,"Y")
ROOT.gStyle.SetTitleFont(62,"Y")
ROOT.gStyle.SetTitleOffset(1.0,"Y")
#ROOT.gStyle.SetTitleX(0.22)
#ROOT.gStyle.SetTitleY(0.98)
ROOT.gStyle.SetTitleSize(0.07,"t")
ROOT.gStyle.SetTitleBorderSize(0)
ROOT.gStyle.SetCanvasBorderSize(0)


#f = ROOT.TFile("/build/kirby/cc1muNp_ubxsec_integration/ubxsec_test_ntuples/ubxsec_output_mc_bnbcosmic.root")
#f = ROOT.TFile("/build/kirby/cc1muNp_ubxsec_pid_integration_test_larana/testing_pid_larana/ubxsec_output_mc_bnbcosmic_Nov19_test_v1.root")
#f = ROOT.TFile("/build/kirby/cc1muNp_ubxsec_pid_integration_test_larana/testing_pid_larana/ubxsec_output_mc_bnbcosmic.root")
#f = ROOT.TFile("/build/kirby/cc1muNp_ubxsec_pid_integration_test_larana/test_ntuples_Nov26/ubxsec_output_mc_bnbcosmic.root")
#f = ROOT.TFile("/build/kirby/cc1muNp_ubxsec_pid_integration_test_larana/test_ntuples_Nov26/ubxsec_output_data_onbeam.root")
#f = ROOT.TFile("/build/kirby/cc1muNp_ubxsec_pid_integration_test_larana/test_ntuples_Dec20/ubxsec_output_mc_bnbcosmic.root")
f = ROOT.TFile("/uboone/data/users/kirby/cc1muNp_ubxsec_ntuples/ntuples_Mar20_merge/ubxsec_output_mc_bnbcosmic.root")
f_T3 = ROOT.TFile("/uboone/data/users/kirby/cc1muNp_ubxsec_ntuples/ntuples_Apr05_merge/ubxsec_output_mc_bnbcosmic_Tune3_merged.root")
#f = ROOT.TFile("/build/kirby/cc1muNp_ubxsec_pid_integration_test_larana/test_ntuples_Dec20/ubxsec_output_mc_bnbdirt.root")

ROOT.gROOT.SetBatch(1)
#ROOT.gStyle.SetOptStat(0)

#t = f.Get("SimpleAna/cc1unptree")
tree = dict()
tree{"MC"] = f.Get("UBXSec/tree")
tree{"T3"] = f_T3.Get("UBXSec/tree")

#tree.MakeClass("UBXSecEvent")
#ROOT.gROOT.ProcessLine(".L /uboone/app/users/afurmans/CCNproton_ccinc/software_builds/larsoft_06_26_01_22/srcs/uboonecode/uboone/UBXSec/DataTypes/UBXSecEvent.h")
#ROOT.gROOT.ProcessLine(".L /build/kirby/cc1muNp_ubxsec_pid_integration_test_larana/srcs/uboonecode/uboone/UBXSec/DataTypes/UBXSecEvent.h")
ROOT.gROOT.ProcessLine(".L /build/kirby/cc1muNp_pandora_update_v27_test/srcs/uboonecode/uboone/UBXSec/DataTypes/UBXSecEvent.h")
#ROOT.gROOT.LoadMacro("UBXSec.C+")

f_conversionGraph = ROOT.TFile("protonLengthFromMom.root")
g_lengthFromMom = f_conversionGraph.Get("Graph")

h_reso = dict()
h_trueDist = dict()
h_reso400 = dict()
h_trueDist400 = dict()
h_reso500 = dict()
h_trueDist500 = dict()
h_trueDistBadProtons = dict()
h_trueRecoBadProtons = dict()
h_trueDistGoodProtons = dict()
h_trueRecoGoodProtons = dict()

h_reso_openCut = dict()
h_trueDist_openCut = dict()
h_reso_showers = dict()
h_trueDist_showers = dict()

variables = ["protonMom","muonMom","protonAngle","muonAngle","thetamup", "costhetamup"] ## worry about thetamup later...
bins = dict()
setBins(bins)
b_arr = dict()
for var in variables:
  b_arr[var] = array.array("d",bins[var])
#b_pmom = array.array("d",bins["protonMom"])
#b_mumom = array.array("d",bins["muonMom"])
#b_muangle = array.array("d",bins["muonAngle"])
#b_muangle = array.array("d",bins["muonAngle"])
h_2D_smearing = dict()
h_2D_smearing_showers = dict()

for var in variables:
  h_trueDistBadProtons[var] = ROOT.TH1D("",";"+var+";",50,0,-1)
  h_trueRecoBadProtons[var] = ROOT.TH2D("",var+";True;Reco",50,0,-1, 50,0,-1)
  h_trueDistGoodProtons[var] = ROOT.TH1D("",";"+var+";",50,0,-1)
  h_trueRecoGoodProtons[var] = ROOT.TH2D("",var+";True;Reco",50,0,-1, 50,0,-1)
  h_2D_smearing[var] = ROOT.TH2D("SmearingMatrix"+var,var+";True;reco",len(b_arr[var])-1, b_arr[var], len(b_arr[var])-1, b_arr[var])
  h_2D_smearing_showers[var] = ROOT.TH2D("SmearingMatrix_showers_"+var,var+";True;reco",len(b_arr[var])-1, b_arr[var], len(b_arr[var])-1, b_arr[var])
  for ibin in xrange(len(bins[var])-1):
    h_trueDist[var,ibin] = ROOT.TH1D("",";;",50,0,-1)
    if var == "protonMom":
      h_reso[var,ibin] = ROOT.TH1D("",";;",100,-0.3,1)
    else:
      h_reso[var,ibin] = ROOT.TH1D("",";;",50,0,-1)
    h_trueDist[var,ibin].SetDirectory(0)
    h_reso[var,ibin].SetDirectory(0)
    h_trueDist400[var,ibin] = ROOT.TH1D("",";;",50,0,-1)
    h_reso400[var,ibin] = ROOT.TH1D("",";;",50,0,-1)
    h_trueDist400[var,ibin].SetDirectory(0)
    h_reso400[var,ibin].SetDirectory(0)
    h_trueDist500[var,ibin] = ROOT.TH1D("",";;",50,0,-1)
    h_reso500[var,ibin] = ROOT.TH1D("",";;",50,0,-1)
    h_trueDist500[var,ibin].SetDirectory(0)
    h_reso500[var,ibin].SetDirectory(0)
    h_trueDist_openCut[var,ibin] = ROOT.TH1D("",";;",50,0,-1)
    h_reso_openCut[var,ibin] = ROOT.TH1D("",";;",50,0,-1)
    h_trueDist_openCut[var,ibin].SetDirectory(0)
    h_reso_openCut[var,ibin].SetDirectory(0)
    h_trueDist_showers[var,ibin] = ROOT.TH1D("",";;",50,0,-1)
    h_reso_showers[var,ibin] = ROOT.TH1D("",";;",50,0,-1)
    h_trueDist_showers[var,ibin].SetDirectory(0)
    h_reso_showers[var,ibin].SetDirectory(0)

h_thetamup = ROOT.TH1D("h_thetamup","#theta_{#mu-p}",50,0,3.142)
h_thetamup_cos = ROOT.TH1D("h_thetamup_cos",";cos(#theta_{#mu-p});Events",50,-1,1)

h_protonAngle_2D_good = ROOT.TH2D("h_protonAngle_2D_good", "", 2, -1, 1, 2, -1,1)
h_protonAngle_2D_bad = h_protonAngle_2D_good.Clone("h_protonAngle_2D_bad")

h_protonMom_good = ROOT.TH1D("h_protonMom_good",";Momentum;Events",len(b_arr["protonMom"])-1,b_arr["protonMom"])
h_protonMom_bad = h_protonMom_good.Clone("h_protonMom_bad")

h_thetaMuP_trueReco_all = ROOT.TH2D("h_thetaMuP_trueReco_all","",50,-1,1,50,-1,1)
h_thetaMuP_trueReco_good = ROOT.TH2D("h_thetaMuP_trueReco_good","",50,-1,1,50,-1,1)
h_thetaMuP_trueReco_bad = ROOT.TH2D("h_thetaMuP_trueReco_bad","",50,-1,1,50,-1,1)

h_LLR_bwdfwd_good = ROOT.TH1D("h_LLR_bwdfwd_good","h_LLR_bwdfwd_good;LLR_{bwd} - LLR_{fwd};",50,-3,3)
h_LLR_bwdfwd_bad = h_LLR_bwdfwd_good.Clone("h_LLR_bwdfwd_bad")
h_LLR_bwdfwd_good_bin0 = h_LLR_bwdfwd_good.Clone("h_LLR_bwdfwd_good_bin0")
h_LLR_bwdfwd_bad_bin0 = h_LLR_bwdfwd_good.Clone("h_LLR_bwdfwd_bad_bin0")

h_trueRecoMom_badOA = ROOT.TH2D("h_trueRecoMom_badOA","h_trueRecoMom_badOA;True momentum;Reco momentum",30,0.2,1.4,30,0.2,1.4)
h_trueRecoMom_goodOA = ROOT.TH2D("h_trueRecoMom_goodOA","h_trueRecoMom_goodOA;True momentum;Reco momentum",30,0.2,1.4,30,0.2,1.4)

h_stopped_trueMom = ROOT.TH1D("h_stopped_trueMom","",20,0,1.4)
h_int_trueMom = h_stopped_trueMom.Clone("h_int_trueMom")
h_stopped_trueMom_sel = h_stopped_trueMom.Clone("h_stopped_trueMom_sel")
h_int_trueMom_sel = h_stopped_trueMom.Clone("h_int_trueMom_sel")

h_trueRecoMom_stopped = h_trueRecoMom_badOA.Clone("h_trueRecoMom_stopped")
h_trueRecoMom_int = h_trueRecoMom_badOA.Clone("h_trueRecoMom_int")

h_reso_protonMom_bwdBin = ROOT.TH1D("",";;",50,0,-1)
h_reso_protonMom_bwdBinAndBadReco = ROOT.TH1D("",";;",50,0,-1)

N_entries = tree["MC"].GetEntries()
perc = int(N_entries/100.)
dec = int(N_entries/10.)
print "running over",N_entries,"events"
i=0

h_protonMom_correctProton = h_stopped_trueMom.Clone("h_protonMom_correctProton")
h_protonMom_incorrectProton = h_stopped_trueMom.Clone("h_protonMom_incorrectProton")
h_protonMom_incorrectProton_mu = h_stopped_trueMom.Clone("h_protonMom_incorrectProton_mu")
h_protonMom_incorrectProton_cos = h_stopped_trueMom.Clone("h_protonMom_incorrectProton_cos")
h_protonMom_incorrectMuon_p = h_stopped_trueMom.Clone("h_protonMom_incorrectMuons_p")
h_protonMom_incorrectMuon_cos = h_stopped_trueMom.Clone("h_protonMom_incorrectMuons_cos")

h_muonMom = ROOT.TH1D("h_muonMom","h_muonMom",len(b_arr["muonMom"])-1, b_arr["muonMom"])
h_muonMom_correctProton = h_muonMom.Clone("h_muonMom_correctProton")
h_muonMom_incorrectProton = h_muonMom.Clone("h_muonMom_incorrectProton")
h_muonMom_incorrectProton_mu = h_muonMom.Clone("h_muonMom_incorrectProton_mu")
h_muonMom_incorrectProton_cos = h_muonMom.Clone("h_muonMom_incorrectProton_cos")
h_muonMom_incorrectMuon_p = h_muonMom.Clone("h_muonMom_incorrectMuons_p")
h_muonMom_incorrectMuon_cos = h_muonMom.Clone("h_muonMom_incorrectMuons_cos")

#b_muangle = array.array("d",bins["muonAngle"])
h_muonAngle = ROOT.TH1D("h_muonAngle","h_muonAngle",len(b_arr["muonAngle"])-1, b_arr["muonAngle"])
h_muonAngle_correctProton = h_muonAngle.Clone("h_muonAngle_correctProton")
h_muonAngle_incorrectProton = h_muonAngle.Clone("h_muonAngle_incorrectProton")
h_muonAngle_incorrectProton_mu = h_muonAngle.Clone("h_muonAngle_incorrectProton_mu")
h_muonAngle_incorrectProton_cos = h_muonAngle.Clone("h_muonAngle_incorrectProton_cos")
h_muonAngle_incorrectMuon_p = h_muonAngle.Clone("h_muonAngle_incorrectMuons_p")
h_muonAngle_incorrectMuon_cos = h_muonAngle.Clone("h_muonAngle_incorrectMuons_cos")

h_CCinc_truedotreco = ROOT.TH1D("h_CCinc_truedotreco","h_CCinc_truedotreco;true dot reco;event fraction",50,-1,1)
h_CCinc_npfp_goodreco = ROOT.TH1D("h_CCinc_npfp_goodreco","h_CCinc_npfp_goodreco;N_{pfp};event fraction",5,0,5)
h_CCinc_npfp_badreco = ROOT.TH1D("h_CCinc_npfp_badreco","h_CCinc_npfp_badreco;N_{pfp};event fraction",5,0,5)

h_CCinc_process_badreco = ROOT.TH1D("h_CCinc_process_badreco","h_CCinc_process_badreco;mode;event fraction",11,0,11)
h_CCinc_process_goodreco = ROOT.TH1D("h_CCinc_process_goodreco","h_CCinc_process_goodreco;mode;event fraction",11,0,11)

N_scatter = 0
N_noscatter = 0

for event in tree["MC"]:
  i+=1
  if i%perc == 0 and i>perc-10 and i < dec:
    print int(100*float(i)/N_entries),"%"
  if i%dec == 0 and i>dec-10:
    #break
    print int(100*float(i)/N_entries),"%"
    
  e = event.ubxsec_event_split
  proton_true_mom = getLeadingProtonMom(e)
  proton_reco_mom = getLeadingRecoMom(e)
  muon_true_mom = getTrueValue(e, "muonMom")
  muon_true_costheta = getTrueValue(e, "muonAngle")
  if isSignal(e,300,0):
    #if leadingProtonStopped(e, g_lengthFromMom):
    if allProtonsStopped(e, g_lengthFromMom): # check ALL protons :)
      h_stopped_trueMom.Fill(proton_true_mom)
    else:
      h_int_trueMom.Fill(proton_true_mom)
  if e.is_selected:
    npfp = sum(1 for x in xrange(len(e.pfp_reco_istrack)) if bool(e.pfp_reco_istrack[x] or e.pfp_reco_isshower[x]))
    if (e.is_signal and not IsCosmic(e)):
      ## CCinc plots
      muonTrueDir = getTrueMuonDir(e)
      muonRecoDir = getRecoMuonDir(e)
      #if muonTrueDir.Mag()<0.5:
      #  print "true length", muonTrueDir.Mag()
      #if muonRecoDir.Mag()<0.5:
      #  print "reco length", muonRecoDir.Mag()
      dotProduct = muonTrueDir.Dot(muonRecoDir)
      h_CCinc_truedotreco.Fill(dotProduct)
      if math.fabs(dotProduct)>1:
        print "uh-oh, dot products can't be greater than 1!"
      if dotProduct<0:
        h_CCinc_npfp_badreco.Fill(npfp)
        h_CCinc_process_badreco.Fill(e.mode)
      else:
        h_CCinc_npfp_goodreco.Fill(npfp)
        h_CCinc_process_goodreco.Fill(e.mode)
    ## CCNp plots
    if npfp>1 and isSignal(e,300,0):
      pind = getRecoProtonIndex(e)
      muind = getMuonIndex(e)
      if pind>=0:
        if e.pfp_truth_pdg[pind]==2212 and e.pfp_truth_origin[pind]==1:
          h_protonMom_correctProton.Fill(proton_true_mom)
          h_muonMom_correctProton.Fill(muon_true_mom)
          h_muonAngle_correctProton.Fill(muon_true_costheta)
        else:
          h_protonMom_incorrectProton.Fill(proton_true_mom)
          h_muonMom_incorrectProton.Fill(muon_true_mom)
          h_muonAngle_incorrectProton.Fill(muon_true_costheta)
        if abs(e.pfp_truth_pdg[pind])==13 and e.pfp_truth_origin[pind]==1:
          h_protonMom_incorrectProton_mu.Fill(proton_true_mom)
          h_muonMom_incorrectProton_mu.Fill(muon_true_mom)
          h_muonAngle_incorrectProton_mu.Fill(muon_true_costheta)
        if e.pfp_truth_origin[pind]!=1:
          h_protonMom_incorrectProton_cos.Fill(proton_true_mom)
          h_muonMom_incorrectProton_cos.Fill(muon_true_mom)
          h_muonAngle_incorrectProton_cos.Fill(muon_true_costheta)
        if e.pfp_truth_origin[muind]!=1:
          h_protonMom_incorrectMuon_cos.Fill(proton_true_mom)
          h_muonMom_incorrectMuon_cos.Fill(muon_true_mom)
          h_muonAngle_incorrectMuon_cos.Fill(muon_true_costheta)
        if e.pfp_truth_pdg[muind]==2212:
          h_protonMom_incorrectMuon_p.Fill(proton_true_mom)
          h_muonMom_incorrectMuon_p.Fill(muon_true_mom)
          h_muonAngle_incorrectMuon_p.Fill(muon_true_costheta)
    if not passContainmentCut(e):
      continue
    if passChi2Cut(e) and passMinHitsCut(e):
      if not isSignal(e,300,0):
        continue
      if not passThreshCut(e):
        continue
      if not muCandpCandTruth(e, True):
        continue
      tmp = getRecoValue(e,"thetamup")
      tmp_cos = getRecoValue(e,"costhetamup")
      h_thetamup.Fill(tmp)
      h_thetamup_cos.Fill(tmp_cos)
      ## Backwards proton check...
      muPos = ROOT.TVector3(-999,-999,-999)
      pPos = getRecoValue(e,"protonStart")
      pEndPos = getRecoValue(e,"protonEnd")
      for i_pfp in e.pfp_reco_ismuoncandidate:
        if e.pfp_reco_ismuoncandidate[i_pfp]:
          muPos.SetXYZ(e.pfp_reco_startx[i_pfp], e.pfp_reco_starty[i_pfp], e.pfp_reco_startz[i_pfp])
      startDist = muPos - pPos
      endDist = muPos - pEndPos
      if startDist.Mag() < endDist.Mag():
        h_protonAngle_2D_good.Fill(getTrueValue(e,"protonAngle"), getRecoValue(e,"protonAngle"))
      if startDist.Mag() > endDist.Mag():
        h_protonAngle_2D_bad.Fill(getTrueValue(e,"protonAngle"), getRecoValue(e,"protonAngle"))
      proton_recoTrueDot = getRecoProtonDir(e).Dot(getTrueProtonDir(e))
      h_thetaMuP_trueReco_all.Fill(getTrueValue(e,"costhetamup"), getRecoValue(e,"costhetamup"))
      if proton_recoTrueDot>0:
        h_thetaMuP_trueReco_good.Fill(getTrueValue(e,"costhetamup"), getRecoValue(e,"costhetamup"))
        h_protonMom_good.Fill(getTrueValue(e,"protonMom"))
      else:
        h_thetaMuP_trueReco_bad.Fill(getTrueValue(e,"costhetamup"), getRecoValue(e,"costhetamup"))
        h_protonMom_bad.Fill(getTrueValue(e,"protonMom"))
      if proton_recoTrueDot>0:
        if e.pfp_reco_bragg_fwd_proton[getRecoProtonIndex(e)]>0 and e.pfp_reco_bragg_bwd_proton[getRecoProtonIndex(e)]>0:
          h_LLR_bwdfwd_good.Fill(math.log(e.pfp_reco_bragg_bwd_proton[getRecoProtonIndex(e)] / e.pfp_reco_bragg_fwd_proton[getRecoProtonIndex(e)]))
          if getRecoValue(e,"protonAngle")<-0.5:
            h_LLR_bwdfwd_good_bin0.Fill(math.log(e.pfp_reco_bragg_bwd_proton[getRecoProtonIndex(e)] / e.pfp_reco_bragg_fwd_proton[getRecoProtonIndex(e)]))
      else:
        if e.pfp_reco_bragg_fwd_proton[getRecoProtonIndex(e)]>0 and e.pfp_reco_bragg_bwd_proton[getRecoProtonIndex(e)]>0:
          h_LLR_bwdfwd_bad.Fill(math.log(e.pfp_reco_bragg_bwd_proton[getRecoProtonIndex(e)] / e.pfp_reco_bragg_fwd_proton[getRecoProtonIndex(e)]))
          if getRecoValue(e,"protonAngle")<-0.5:
            h_LLR_bwdfwd_bad_bin0.Fill(math.log(e.pfp_reco_bragg_bwd_proton[getRecoProtonIndex(e)] / e.pfp_reco_bragg_fwd_proton[getRecoProtonIndex(e)]))
      if getRecoValue(e,"costhetamup")<-0.8 and getTrueValue(e,"costhetamup")>-0.8:
        h_trueRecoMom_badOA.Fill(getTrueValue(e,"protonMom"),getRecoValue(e,"protonMom"))
      if getRecoValue(e,"costhetamup")<-0.8 and getTrueValue(e,"costhetamup")<-0.8:
        h_trueRecoMom_goodOA.Fill(getTrueValue(e,"protonMom"),getRecoValue(e,"protonMom"))
      #if leadingProtonStopped(e, g_lengthFromMom):
      if allProtonsStopped(e, g_lengthFromMom): # checking ALL protons :)
        h_stopped_trueMom_sel.Fill(proton_true_mom)
        h_trueRecoMom_stopped.Fill(proton_true_mom, getRecoValue(e,"protonMom"))
      else:
        h_int_trueMom_sel.Fill(proton_true_mom)
        h_trueRecoMom_int.Fill(proton_true_mom, getRecoValue(e,"protonMom"))
      for var in variables:
        reco = getRecoValue(e, var)
        true = getTrueValue(e, var)
        ptrue = getTrueValue(e, "protonAngle")
        preco = getRecoValue(e, "protonAngle")
        if preco > bins["protonAngle"][0] and preco <= bins["protonAngle"][1]: ## quick test to look at bad proton reco events
          if var=="protonMom":
            h_reso_protonMom_bwdBin.Fill(true - reco)
            if ptrue > bins["protonAngle"][2]:
              h_reso_protonMom_bwdBinAndBadReco.Fill(true - reco)
          if ptrue > bins["protonAngle"][1]:
            h_trueDistBadProtons[var].Fill(true)
            h_trueRecoBadProtons[var].Fill(true, reco)
            N_scatter += not leadingProtonStopped(e,g_lengthFromMom)
            N_noscatter += leadingProtonStopped(e,g_lengthFromMom)
          else:
            h_trueDistGoodProtons[var].Fill(true)
            h_trueRecoGoodProtons[var].Fill(true, reco)
        h_2D_smearing[var].Fill(true, reco)
        if e.pfp_reco_isshower[pind]:
          h_2D_smearing_showers[var].Fill(true, reco)
        for ibin in xrange(len(bins[var])-1):
          if reco > bins[var][ibin] and reco <= bins[var][ibin+1]:
            h_trueDist[var,ibin].Fill(true)
            h_reso[var,ibin].Fill(true-reco)
            if isSignal(e,400,0): ## Checking what happens if we have a higher threshold
              h_trueDist400[var,ibin].Fill(true)
              h_reso400[var,ibin].Fill(true-reco)
            if isSignal(e,500,0): ## Checking what happens if we have a higher threshold
              h_trueDist500[var,ibin].Fill(true)
              h_reso500[var,ibin].Fill(true-reco)
#            if var == "muonAngle":
#              print "muon angle (t,r)", true, reco
            if getRecoValue(e, "thetamup") < 2.7:
              h_trueDist_openCut[var,ibin].Fill(true)
              h_reso_openCut[var,ibin].Fill(true - reco)
            if e.pfp_reco_isshower[pind]:
              h_trueDist_showers[var,ibin].Fill(true)
              h_reso_showers[var,ibin].Fill(true - reco)
            break

print N_scatter,"bad protons scattered,", N_noscatter,"didn't"

canv = ROOT.TCanvas("","",800,400)

h_CCinc_truedotreco.Scale(1./h_CCinc_truedotreco.Integral())
h_CCinc_truedotreco.Draw("hist")
canv.SaveAs("figures/CCinclusive/h_CCinc_truedotreco.eps")
canv.SaveAs("figures/CCinclusive/h_CCinc_truedotreco.png")

canv.Clear()
h_CCinc_npfp_goodreco.SetLineColor(1)
h_CCinc_npfp_badreco.SetLineColor(2)
legGoodBad = ROOT.TLegend(0.6, 0.6, 0.89, 0.89)
legGoodBad.AddEntry(h_CCinc_npfp_goodreco,"good reco dir","l")
legGoodBad.AddEntry(h_CCinc_npfp_badreco,"bad reco dir","l")
h_CCinc_npfp_goodreco.Scale(1./h_CCinc_npfp_goodreco.Integral())
h_CCinc_npfp_badreco.Scale(1./h_CCinc_npfp_badreco.Integral())
h_CCinc_npfp_badreco.Draw("hist")
h_CCinc_npfp_goodreco.Draw("hist same")
legGoodBad.Draw()
canv.SaveAs("figures/CCinclusive/h_CCinc_npfp_badreco.eps")
canv.SaveAs("figures/CCinclusive/h_CCinc_npfp_badreco.png")

canv.Clear()
h_CCinc_process_goodreco.SetLineColor(1)
h_CCinc_process_badreco.SetLineColor(2)
h_CCinc_process_goodreco.Scale(1./h_CCinc_process_goodreco.Integral())
h_CCinc_process_badreco.Scale(1./h_CCinc_process_badreco.Integral())
h_CCinc_process_badreco.GetXaxis().SetBinLabel(1,"QE")
h_CCinc_process_badreco.GetXaxis().SetBinLabel(2,"RES")
h_CCinc_process_badreco.GetXaxis().SetBinLabel(3,"DIS")
h_CCinc_process_badreco.GetXaxis().SetBinLabel(4,"COH")
h_CCinc_process_badreco.GetXaxis().SetBinLabel(11,"MEC")
h_CCinc_process_badreco.Draw("hist")
h_CCinc_process_goodreco.Draw("hist same")
legGoodBad.Draw()
canv.SaveAs("figures/CCinclusive/h_CCinc_process_badreco.eps")
canv.SaveAs("figures/CCinclusive/h_CCinc_process_badreco.png")


canv.Clear()
canv.Divide(2)
h_thetamup.SetMinimum(0)
h_thetamup_cos.SetMinimum(0)
canv.cd(1)
h_thetamup.Draw("e")
canv.cd(2)
h_thetamup_cos.Draw("e")
canv.SaveAs("h_thetamup_vs_cos.eps")
h_thetamup_bincont = ROOT.TH1D("h_thetamup_bincont",";Bin Content / Events;Bins",50,0,-1)
h_thetamup_bincont_cos = ROOT.TH1D("h_thetamup_bincont_cos",";Bin Content / Events;Bins",50,0,-1)
for i in xrange(h_thetamup.GetNbinsX()+1):
  h_thetamup_bincont.Fill(h_thetamup.GetBinContent(i))
for i in xrange(h_thetamup_cos.GetNbinsX()+1):
  h_thetamup_bincont_cos.Fill(h_thetamup_cos.GetBinContent(i))
print "thetamup, mean and rms bin content = ", h_thetamup_bincont.GetMean(), h_thetamup_bincont.GetRMS()
print "cos(thetamup), mean and rms bin content = ", h_thetamup_bincont_cos.GetMean(), h_thetamup_bincont_cos.GetRMS()
canv.cd(1)
h_protonAngle_2D_good.Draw("col")
h_protonAngle_2D_good.Draw("text same")
canv.cd(2)
h_protonAngle_2D_bad.Draw("col")
h_protonAngle_2D_bad.Draw("text same")
canv.SaveAs("protonAngleSmearing_vtxReso.eps")


canv.Clear()
canv.Divide(3)
canv.cd(1)
h_thetaMuP_trueReco_all.Draw("colz")
h_thetaMuP_trueReco_all.Draw("text same")
canv.cd(2)
h_thetaMuP_trueReco_good.Draw("colz")
h_thetaMuP_trueReco_good.Draw("text same")
canv.cd(3)
h_thetaMuP_trueReco_bad.Draw("colz")
h_thetaMuP_trueReco_bad.Draw("text same")
canv.SaveAs("h_thetamup_2D_goodProtonDir_badProtonDir.eps")


canv.Clear()
h_protonMom_all = h_protonMom_good.Clone("h_protonMom_all")
h_protonMom_all.Sumw2()
h_protonMom_good.Sumw2()
h_protonMom_bad.Sumw2()
h_protonMom_all.Add(h_protonMom_bad)
h_protonMom_bad.Divide(h_protonMom_all)
h_protonMom_bad.SetLineColor(1)
h_protonMom_bad.SetLineWidth(2)
h_protonMom_bad.SetTitle("Fraction of protons reco'd backwards;Momentum;Fraction")
h_protonMom_bad.Draw("e")
canv.SaveAs("h_protonMom_goodBadDir.eps")

canv.Clear()
h_LLR_bwdfwd_good.SetLineColor(1)
h_LLR_bwdfwd_bad.SetLineColor(2)
h_LLR_bwdfwd_good.Sumw2()
h_LLR_bwdfwd_bad.Sumw2()
#h_LLR_bwdfwd_good.Scale(1./h_LLR_bwdfwd_good.Integral())
#h_LLR_bwdfwd_bad.Scale(1./h_LLR_bwdfwd_bad.Integral())
h_LLR_bwdfwd_good.Draw("hist")
h_LLR_bwdfwd_bad.Draw("hist same")
canv.SetLogy()
canv.SaveAs("h_LLR_bwdfwd_comparison_log.eps")
canv.SetLogy(0)
canv.SaveAs("h_LLR_bwdfwd_comparison.eps")

canv.Clear()
h_LLR_bwdfwd_good_bin0.SetLineColor(1)
h_LLR_bwdfwd_bad_bin0.SetLineColor(2)
h_LLR_bwdfwd_good_bin0.Sumw2()
h_LLR_bwdfwd_bad_bin0.Sumw2()
#h_LLR_bwdfwd_good.Scale(1./h_LLR_bwdfwd_good.Integral())
#h_LLR_bwdfwd_bad.Scale(1./h_LLR_bwdfwd_bad.Integral())
h_LLR_bwdfwd_good_bin0.Draw("hist")
h_LLR_bwdfwd_bad_bin0.Draw("hist same")
canv.SetLogy()
canv.SaveAs("h_LLR_bwdfwd_comparison_bin0_log.eps")
canv.SetLogy(0)
canv.SaveAs("h_LLR_bwdfwd_comparison_bin0.eps")

canv.Clear()
canv.Divide(2)
canv.cd(1)
h_trueRecoMom_goodOA.Draw("colz")
canv.cd(2)
h_trueRecoMom_badOA.Draw("colz")
canv.SaveAs("figures/resolutions/h_trueRecoMom_goodbadOA.eps")
canv.SaveAs("figures/resolutions/h_trueRecoMom_goodbadOA.png")


canv.Clear()
h_eff_stopped_trueMom = h_stopped_trueMom_sel.Clone()
h_eff_stopped_trueMom.Sumw2()
h_stopped_trueMom.Sumw2()
h_eff_stopped_trueMom.Divide(h_stopped_trueMom)
h_eff_int_trueMom = h_int_trueMom_sel.Clone()
h_eff_int_trueMom.Sumw2()
h_int_trueMom.Sumw2()
h_eff_int_trueMom.Divide(h_int_trueMom)
h_eff_stopped_trueMom.SetLineColor(1)
h_eff_int_trueMom.SetLineColor(2)
h_eff_stopped_trueMom.Draw("e")
h_eff_int_trueMom.Draw("e same")
leg = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
leg.AddEntry(h_eff_stopped_trueMom,"stopping","l")
leg.AddEntry(h_eff_int_trueMom,"interacting","l")
leg.Draw()
canv.SaveAs("figures/eff_protonMom_stoppedVsNot.eps")



h_protonMom_misIDrate = h_protonMom_incorrectProton.Clone("h_protonMom_misIDrate;Leading proton momemtum / GeV;miss-ID rate")
h_protonMom_muIDrate = h_protonMom_incorrectProton_mu.Clone("h_protonMom_muIDrate;Leading proton momemtum / GeV;muon-proton confusion rate")
h_protonMom_cosIDrate = h_protonMom_incorrectProton_cos.Clone("h_protonMom_muIDrate;Leading proton momemtum / GeV;muon-proton confusion rate")
h_protonMom_MuonpIDrate = h_protonMom_incorrectMuon_p.Clone("h_protonMom_MuonpIDrate;Leading proton momemtum / GeV;muon-proton confusion rate")
h_protonMom_MuoncosIDrate = h_protonMom_incorrectMuon_cos.Clone("h_protonMom_MuoncosIDrate;Leading proton momemtum / GeV;muon-proton confusion rate")
h_protonMom_all = h_protonMom_incorrectProton.Clone("h_protonMom_all")
h_protonMom_all.Add(h_protonMom_correctProton)
h_protonMom_misIDrate.Sumw2()
h_protonMom_muIDrate.Sumw2()
h_protonMom_cosIDrate.Sumw2()
h_protonMom_MuonpIDrate.Sumw2()
h_protonMom_MuoncosIDrate.Sumw2()
h_protonMom_all.Sumw2()
h_protonMom_misIDrate.Divide(h_protonMom_all)
h_protonMom_misIDrate.SetLineColor(1)
h_protonMom_misIDrate.Draw("e")
h_protonMom_muIDrate.Divide(h_protonMom_all)
h_protonMom_muIDrate.SetLineColor(2)
h_protonMom_muIDrate.Draw("e same")
h_protonMom_cosIDrate.Divide(h_protonMom_all)
h_protonMom_cosIDrate.SetLineColor(3)
h_protonMom_cosIDrate.Draw("e same")
legmissID = ROOT.TLegend(0.1, 0.7, 0.3, 0.89)
legmissID.AddEntry(h_protonMom_misIDrate,"all mis-ID","l")
legmissID.AddEntry(h_protonMom_muIDrate,"mis-ID true mu","l")
legmissID.AddEntry(h_protonMom_cosIDrate,"mis-ID true cos","l")
legmissID.Draw()
canv.SaveAs("figures/h_protonMom_misIDrate.eps")

## Muon mis ID as a function of proton momentum...
canv.Clear()
h_protonMom_MuonpIDrate.Divide(h_protonMom_all)
h_protonMom_MuonpIDrate.SetLineColor(2)
h_protonMom_MuonpIDrate.Draw("e")
h_protonMom_MuoncosIDrate.Divide(h_protonMom_all)
h_protonMom_MuoncosIDrate.SetLineColor(3)
h_protonMom_MuoncosIDrate.Draw("e same")
legMuID = ROOT.TLegend(0.1, 0.7, 0.3, 0.89)
legMuID.AddEntry(h_protonMom_MuonpIDrate,"muon cand is proton","l")
legMuID.AddEntry(h_protonMom_MuoncosIDrate,"muon cand is cosmic","l")
legMuID.Draw()
canv.SaveAs("figures/h_protonMom_MuonmisIDrate.eps")

### Muon mis-ID as a function of muon momentum
canv.Clear()
h_muonMom_MuonpIDrate = h_muonMom_incorrectMuon_p.Clone("h_muonMom_MuonpIDrate;Muon momemtum / GeV;muon-proton confusion rate")
h_muonMom_MuonpIDrate.Sumw2()
h_muonMom_all = h_muonMom_incorrectProton.Clone("h_muonMom_all")
h_muonMom_all.Sumw2()
h_muonMom_correctProton.Sumw2()
h_muonMom_all.Add(h_muonMom_correctProton)
h_muonMom_MuonpIDrate.Divide(h_muonMom_all)
h_muonMom_MuonpIDrate.SetLineColor(2)
h_muonMom_MuonpIDrate.Draw("e")
canv.SaveAs("figures/h_muonMom_MuonmisIDrate.eps")

### Muon mis-ID as a function of muon cos theta
canv.Clear()
h_muonAngle_MuonpIDrate = h_muonAngle_incorrectMuon_p.Clone("h_muonAngle_MuonpIDrate;Muon cos(#theta);muon-proton confusion rate")
h_muonAngle_MuonpIDrate.Sumw2()
h_muonAngle_all = h_muonAngle_incorrectProton.Clone("h_muonAngle_all")
h_muonAngle_all.Sumw2()
h_muonAngle_correctProton.Sumw2()
h_muonAngle_all.Add(h_muonAngle_correctProton)
h_muonAngle_MuonpIDrate.Divide(h_muonAngle_all)
h_muonAngle_MuonpIDrate.SetLineColor(2)
h_muonAngle_MuonpIDrate.Draw("e")
canv.SaveAs("figures/h_muonAngle_MuonmisIDrate.eps")

####################
canv.Clear()
canv.Divide(2)
canv.cd(1)
h_trueRecoMom_stopped.SetTitle("Stopping")
h_trueRecoMom_stopped.Draw("colz")
canv.cd(2)
h_trueRecoMom_int.SetTitle("Interacting")
h_trueRecoMom_int.Draw("colz")
canv.SaveAs("figures/protonMomSmearing_intVsStopping.eps")


canv.Clear()
h_thetamup_bincont.SetLineColor(1)
h_thetamup_bincont.Draw("hist")
h_thetamup_bincont_cos.SetLineColor(2)
h_thetamup_bincont_cos.Draw("hist same")
legtmp = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
legtmp.AddEntry(h_thetamup_bincont,"thetamup","l")
legtmp.AddEntry(h_thetamup_bincont_cos,"cos(thetamup)","l")
legtmp.Draw()
canv.SaveAs("binContents_thetaVsCos.eps")

canvSplit = ROOT.TCanvas("","",3200,800)

leg=ROOT.TLegend(0.7, 0.7, 0.89, 0.89)
leg.SetFillColor(0)
leg.SetLineColor(0)
lowLine = dict()
highLine = dict()
lowLinesh = dict()
highLinesh = dict()



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


for var in variables:
  nbins = len(bins[var])-1
  canvSplit.Clear()
  canvSplit.Divide(nbins,2)
  for ibin in xrange(nbins):
    ######## True distribution ##############
    lowLinesh[var,ibin,"true"] = ROOT.TLine(bins[var][ibin], 0.0, bins[var][ibin],h_trueDist_showers[var,ibin].GetMaximum())
    highLinesh[var,ibin,"true"] = ROOT.TLine(bins[var][ibin+1], 0.0, bins[var][ibin+1],h_trueDist_showers[var,ibin].GetMaximum())
    lowLinesh[var,ibin,"true"].SetLineColor(2)
    highLinesh[var,ibin,"true"].SetLineColor(2)
    lowLinesh[var,ibin,"true"].SetLineStyle(2)
    highLinesh[var,ibin,"true"].SetLineStyle(2)
#    lowLine[var,ibin,"true"].SetDirectory(0)
#    highLine[var,ibin,"true"].SetDirectory(0)
    #canv.cd()
    h_trueDist_showers[var,ibin].SetTitle("reco bin "+str(bins[var][ibin])+" - "+str(bins[var][ibin+1])+";True "+var+";Events")
    h_trueDist_showers[var,ibin].SetLineColor(1)
    #h_trueDist[var,ibin].Draw("hist")
    #lowLine[var,ibin,"true"].Draw()
    #highLine[var,ibin,"true"].Draw()
    #canv.SaveAs("figures/resolutions/"+var+"_trueDistribution_bin_"+str(ibin)+".eps")
    #canv.SaveAs("figures/resolutions/"+var+"_trueDistribution_bin_"+str(ibin)+".png")
    canvSplit.cd(ibin+1)
    h_trueDist_showers[var,ibin].Draw("hist")
    lowLinesh[var,ibin,"true"].Draw()
    highLinesh[var,ibin,"true"].Draw()

    ######## True - reco ##############
    width = bins[var][ibin+1]-bins[var][ibin]
    lowLinesh[var,ibin,"reso"] = ROOT.TLine(-width/2., 0.0, -width/2.,h_reso_showers[var,ibin].GetMaximum())
    highLinesh[var,ibin,"reso"] = ROOT.TLine(width/2., 0.0, width/2., h_reso_showers[var,ibin].GetMaximum())
    lowLinesh[var,ibin,"reso"].SetLineColor(2)
    highLinesh[var,ibin,"reso"].SetLineColor(2)
    lowLinesh[var,ibin,"reso"].SetLineStyle(2)
    highLinesh[var,ibin,"reso"].SetLineStyle(2)
#    lowLine[var,ibin,"reso"].SetDirectory(0)
#    highLine[var,ibin,"reso"].SetDirectory(0)
    #canv.cd()
    h_reso_showers[var,ibin].SetTitle("reco bin "+str(bins[var][ibin])+" - "+str(bins[var][ibin+1])+";(True - Reco) "+var+";Events")
    h_reso_showers[var,ibin].SetLineColor(1)
    #h_reso[var,ibin].Draw("hist")
    #lowLine[var,ibin,"reso"].Draw()
    #highLine[var,ibin,"reso"].Draw()
    #canv.SaveAs("figures/resolutions/"+var+"_reso_bin_"+str(ibin)+".eps")
    #canv.SaveAs("figures/resolutions/"+var+"_reso_bin_"+str(ibin)+".png")
    canvSplit.cd(nbins+ibin+1)
    h_reso_showers[var,ibin].Draw("hist")
    lowLinesh[var,ibin,"reso"].Draw()
    highLinesh[var,ibin,"reso"].Draw()
  canvSplit.SaveAs("figures/resolutions/"+var+"_showers_allPlots.eps")
  canvSplit.SaveAs("figures/resolutions/"+var+"_showers_allPlots.png")

canv = ROOT.TCanvas()
canv.cd()
for var in variables:
  nbins = len(bins[var])-1
  for ibin in xrange(nbins):
    h_trueDist[var,ibin].SetMaximum(h_trueDist[var,ibin].GetMaximum()*1.3)
    h_trueDist[var,ibin].Draw("hist")
    s = ROOT.gPad.GetPrimitive("stats");
    #s.SetTextSize(0.025);
    s.SetX1NDC(0.55);
    s.SetX2NDC(0.89); 
    s.SetY1NDC(0.7);
    s.SetY2NDC(0.89);

    lowLine[var,ibin,"true"].Draw()
    highLine[var,ibin,"true"].Draw()
    canv.SaveAs("figures/resolutions/"+var+"_trueDistribution_bin_"+str(ibin)+".eps")
    canv.SaveAs("figures/resolutions/"+var+"_trueDistribution_bin_"+str(ibin)+".png")
    
    h_reso[var,ibin].SetMaximum(h_reso[var,ibin].GetMaximum()*1.3)
    h_reso[var,ibin].Draw("hist")
    s = ROOT.gPad.GetPrimitive("stats");
    #s.SetTextSize(0.025);
    s.SetX1NDC(0.55);
    s.SetX2NDC(0.89); 
    s.SetY1NDC(0.7);
    s.SetY2NDC(0.89);
    if var=="protonAngle" and ibin>2:
      s.SetX1NDC(0.11)
      s.SetX2NDC(0.45)
    lowLine[var,ibin,"reso"].Draw()
    highLine[var,ibin,"reso"].Draw()
    canv.SaveAs("figures/resolutions/"+var+"_reso_bin_"+str(ibin)+".eps")
    canv.SaveAs("figures/resolutions/"+var+"_reso_bin_"+str(ibin)+".png")

canv.cd()
for var in variables:
  nbins = len(bins[var])-1
  for ibin in xrange(nbins):
    h_trueDist400[var,ibin].Draw("hist")
    lowLine[var,ibin,"true"].Draw()
    highLine[var,ibin,"true"].Draw()
    canv.SaveAs("figures/resolutions/"+var+"_trueDistribution_400thresh_bin_"+str(ibin)+".eps")
    canv.SaveAs("figures/resolutions/"+var+"_trueDistribution_400thresh_bin_"+str(ibin)+".png")
    
    h_reso400[var,ibin].Draw("hist")
    lowLine[var,ibin,"reso"].Draw()
    highLine[var,ibin,"reso"].Draw()
    canv.SaveAs("figures/resolutions/"+var+"_reso_400thresh_bin_"+str(ibin)+".eps")
    canv.SaveAs("figures/resolutions/"+var+"_reso_400thresh_bin_"+str(ibin)+".png")


canv.cd()
for var in variables:
  nbins = len(bins[var])-1
  for ibin in xrange(nbins):
    h_trueDist500[var,ibin].Draw("hist")
    lowLine[var,ibin,"true"].Draw()
    highLine[var,ibin,"true"].Draw()
    canv.SaveAs("figures/resolutions/"+var+"_trueDistribution_500thresh_bin_"+str(ibin)+".eps")
    canv.SaveAs("figures/resolutions/"+var+"_trueDistribution_500thresh_bin_"+str(ibin)+".png")
    
    h_reso500[var,ibin].Draw("hist")
    lowLine[var,ibin,"reso"].Draw()
    highLine[var,ibin,"reso"].Draw()
    canv.SaveAs("figures/resolutions/"+var+"_reso_500thresh_bin_"+str(ibin)+".eps")
    canv.SaveAs("figures/resolutions/"+var+"_reso_500thresh_bin_"+str(ibin)+".png")


canv.Clear()
for var in variables:
  h_trueDistBadProtons[var].SetMinimum(0)
  h_trueDistBadProtons[var].Draw("hist")
  canv.SaveAs("figures/resolutions/"+var+"_protonRecoBackwardsTrueForwards.eps")
  canv.SaveAs("figures/resolutions/"+var+"_protonRecoBackwardsTrueForwards.png")
canv.Clear()
# add in the line of good protons
for var in variables:
  h_trueDistBadProtons[var].SetMinimum(0)
  h_trueDistBadProtons[var].Draw("hist")
  h_trueDistGoodProtons[var].Scale(h_trueDistBadProtons[var].GetMaximum() / h_trueDistGoodProtons[var].GetMaximum() )
  h_trueDistGoodProtons[var].SetLineColor(2)
  h_trueDistGoodProtons[var].Draw("hist same")
  canv.SaveAs("figures/resolutions/"+var+"_goodVsBadProtonReco.eps")
  canv.SaveAs("figures/resolutions/"+var+"_goodVsBadProtonReco.png")


canv.Clear()
for var in variables:
  h_trueRecoBadProtons[var].SetMinimum(0)
  h_trueRecoBadProtons[var].Draw("colz")
  canv.SaveAs("figures/resolutions/"+var+"_TrueVsReco_protonRecoBackwardsTrueForwards.eps")
  canv.SaveAs("figures/resolutions/"+var+"_TrueVsReco_protonRecoBackwardsTrueForwards.png")


for var in variables:
  nbins = len(bins[var])-1
  for ibin in xrange(nbins):
    h_trueDist_openCut[var,ibin].Draw("hist")
    lowLine[var,ibin,"true"].Draw()
    highLine[var,ibin,"true"].Draw()
    canv.SaveAs("figures/resolutions/"+var+"_trueDistribution_openCutthresh_bin_"+str(ibin)+".eps")
    canv.SaveAs("figures/resolutions/"+var+"_trueDistribution_openCutthresh_bin_"+str(ibin)+".png")
    
    h_reso_openCut[var,ibin].Draw("hist")
    lowLine[var,ibin,"reso"].Draw()
    highLine[var,ibin,"reso"].Draw()
    canv.SaveAs("figures/resolutions/"+var+"_reso_openCutthresh_bin_"+str(ibin)+".eps")
    canv.SaveAs("figures/resolutions/"+var+"_reso_openCutthresh_bin_"+str(ibin)+".png")

for var in variables:
  nbins = len(bins[var])-1
  for ibin in xrange(nbins):
    h_trueDist_showers[var,ibin].Draw("hist")
    lowLine[var,ibin,"true"].Draw()
    highLine[var,ibin,"true"].Draw()
    canv.SaveAs("figures/resolutions/"+var+"_trueDistribution_showers_bin_"+str(ibin)+".eps")
    canv.SaveAs("figures/resolutions/"+var+"_trueDistribution_showers_bin_"+str(ibin)+".png")
    
    h_reso_showers[var,ibin].Draw("hist")
    lowLine[var,ibin,"reso"].Draw()
    highLine[var,ibin,"reso"].Draw()
    canv.SaveAs("figures/resolutions/"+var+"_reso_showers_bin_"+str(ibin)+".eps")
    canv.SaveAs("figures/resolutions/"+var+"_reso_showers_bin_"+str(ibin)+".png")


canv.Clear()
for var in variables:
  h_2D_smearing[var].SetMinimum(0)
  h_2D_smearing[var].Draw("colz")
  h_2D_smearing[var].Draw("text same")
  canv.SaveAs("figures/resolutions/"+var+"_2D_smearing.eps")
  canv.SaveAs("figures/resolutions/"+var+"_2D_smearing.png")


canv.Clear()
for var in variables:
  h_2D_smearing_showers[var].SetMinimum(0)
  h_2D_smearing_showers[var].Draw("colz")
  h_2D_smearing_showers[var].Draw("text same")
  canv.SaveAs("figures/resolutions/"+var+"_2D_smearing_showers.eps")
  canv.SaveAs("figures/resolutions/"+var+"_2D_smearing_showers.png")


########################################
## resolution of proton mom for bwd reco'd protons
########################################


h_reso_protonMom_bwdBin.SetTitle("Reco proton angle < -0.5;True - Reco momentum [GeV/c];Events")
h_reso_protonMom_bwdBin.SetLineColor(1)
h_reso_protonMom_bwdBin.Draw("hist")
#lowLine["protonMom",0,"reso"].Draw()
#highLine["protonMom",0,"reso"].Draw()
canv.SaveAs("figures/resolutions/protonMom_resolution_backwardsBin.eps")
canv.SaveAs("figures/resolutions/protonMom_resolution_backwardsBin.png")

h_reso_protonMom_bwdBinAndBadReco.SetTitle("Reco proton angle < -0.5 and True proton angle > 0;True - Reco momentum [GeV/c];Events")
h_reso_protonMom_bwdBinAndBadReco.SetLineColor(1)
h_reso_protonMom_bwdBinAndBadReco.Draw("hist")
#lowLine["protonMom",0,"reso"].Draw()
#highLine["protonMom",0,"reso"].Draw()
canv.SaveAs("figures/resolutions/protonMom_resolution_backwardsBinBadReco.eps")
canv.SaveAs("figures/resolutions/protonMom_resolution_backwardsBinBadReco.png")

h_tmp = h_reso["protonMom",0].Clone()
for i in range(1,len(bins["protonMom"])-1):
  h_tmp.Add(h_reso["protonMom",i])

h_tmp.Scale(1./h_tmp.Integral())
h_reso_protonMom_bwdBin.Scale(1./h_reso_protonMom_bwdBin.Integral())
h_reso_protonMom_bwdBinAndBadReco.Scale(1./h_reso_protonMom_bwdBinAndBadReco.Integral())

h_tmp.SetLineWidth(2)
h_reso_protonMom_bwdBin.SetLineWidth(2)
h_reso_protonMom_bwdBinAndBadReco.SetLineWidth(2)
h_tmp.Draw("hist")
h_reso_protonMom_bwdBin.SetLineColor(2)
h_reso_protonMom_bwdBin.Draw("hist same")
h_reso_protonMom_bwdBinAndBadReco.SetLineColor(4)
h_reso_protonMom_bwdBinAndBadReco.Draw("hist same")

leg_protonReso = ROOT.TLegend(0.55, 0.55, 0.89, 0.89)
leg_protonReso.AddEntry(h_tmp,"all protons","l")
leg_protonReso.AddEntry(h_reso_protonMom_bwdBin,"backwards-most bin","l")
leg_protonReso.AddEntry(h_reso_protonMom_bwdBinAndBadReco,"backwards-most bin, and bad reco","l")
leg_protonReso.Draw()
canv.SaveAs("figures/resolutions/protonMom_resolution_comparison.eps")
canv.SaveAs("figures/resolutions/protonMom_resolution_comparison.png")
