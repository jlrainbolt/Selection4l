from __future__ import print_function
from __future__ import division

import sys

from ROOT import TFile, TTree, TH1D, TNamed, TParameter, TLorentzVector
from array import array

from LHEUtils import *
from SelectionTools import *
from Cuts import *


##
##  INFO
##

n_max = int(sys.argv[1])

selection = ["4l", "4m", "2m2e", "4e"]
channels = {"4l":5, "4m":6, "2m2e":7, "4e":9}



##
##  READ IN LHE
##

infile = LHEfile("unweighted_events.lhe")
infile.setMax(n_max)
events = infile.readEvents()


# Parameters
infile.readInfo()

model = infile.Model
if "U" in model:
    isBSM = True
else:
    isBSM = False

model_name = TNamed("model", model)
print("Model:", model)

xsec = infile.xsec
x_sec = TParameter("Float_t")("xsec", infile.xsec)
print("Cross section:", xsec, "pb")

if isBSM:
    MUb = infile.MU
    m_U = TParameter("Float_t")("mU", MUb)
    print("U mass:", MUb, "GeV")

    WUb = infile.WU
    w_U = TParameter("Float_t")("wU", WUb)
    print("U width:", WUb, "GeV")

    gUbe = infile.gUbe
    g_e = TParameter("Float_t")("ge", gUbe)
    print("Electron coupling:", gUbe)

    gUbmu = infile.gUbmu
    g_mu = TParameter("Float_t")("gmu", gUbmu)
    print("Muon coupling:", gUbmu)



##
##  OUTPUT
##

if isBSM:
    outname = model + "_M-" + "{:g}".format(MUb) + "_ge-" + "{:g}".format(gUbe).replace("0.", "")
    outname = outname + "_gmu-" + "{:g}".format(gUbmu).replace("0.", "") + ".root"
else:
    outname = model + ".root"

outfile = TFile(outname, "RECREATE")
outfile.cd()


# Trees and branches

runNum, evtNum, lumiSec = array("i", [0]), array("i", [0]), array("i", [0])
weight, channel, isFiducial = array("f", [1]), array("I", [0]), array("H", [0])

zzp4, z1p4, z2p4 = TLorentzVector(), TLorentzVector(), TLorentzVector()
z1pdg, z2pdg = array("h", [0]), array("h", [0])

l1p4, l2p4, l3p4, l4p4 = TLorentzVector(), TLorentzVector(), TLorentzVector(), TLorentzVector()
l1pdg, l2pdg, l3pdg, l4pdg = array("h", [0]), array("h", [0]), array("h", [0]), array("h", [0])

tree = {}
for sel in selection:
    tree[sel] = TTree(sel + "_" + model, "")

    tree[sel].Branch("runNum", runNum, "runNum/I")
    tree[sel].Branch("evtNum", evtNum, "evtNum/I")
    tree[sel].Branch("lumiSec", lumiSec, "lumiSec/I")
    tree[sel].Branch("channel", channel, "channel/i")
    tree[sel].Branch("weight", weight, "weight/F")
    tree[sel].Branch("isFiducial", isFiducial, "isFiducial/O")

    tree[sel].Branch("zzp4", zzp4)
    tree[sel].Branch("z1p4", z1p4)
    tree[sel].Branch("z1pdg", z1pdg, "z1pdg/S")
    tree[sel].Branch("z2p4", z2p4)
    tree[sel].Branch("z2pdg", z2pdg, "z2pdg/S")

    tree[sel].Branch("l1p4", l1p4)
    tree[sel].Branch("l1pdg", l1pdg, "l1pdg/S")
    tree[sel].Branch("l2p4", l2p4)
    tree[sel].Branch("l2pdg", l2pdg, "l2pdg/S")
    tree[sel].Branch("l3p4", l3p4)
    tree[sel].Branch("l3pdg", l3pdg, "l3pdg/S")
    tree[sel].Branch("l4p4", l4p4)
    tree[sel].Branch("l4pdg", l4pdg, "l4pdg/S")


# Histograms

hPhaseSpaceEvents = TH1D("PhaseSpaceEvents_" + model, "", 10, 0.5, 10.5)
hFiducialEvents = TH1D("FiducialEvents_" + model, "", 10, 0.5, 10.5)





####
####
####    EVENT LOOP
####
####


print("\nLooping over", n_max, "events...", end='')
sys.stdout.flush()

evt_num = 0
for this_event in events:
    event = LHEevent()
    event.fillEvent(this_event)


    # Get leptons

    particles = event.Particles
    elecs, muons = [], []
    lep_q, lep_id = {}, {}

    for p in particles:
        if abs(p['ID']) == 11:
            elecs.append(TLorentzVector(p['Px'], p['Py'], p['Pz'], p['E']))
            lep_q[elecs[-1]] = -1 * sign(p['ID'])
            lep_id[elecs[-1]] = p['ID']

        elif abs(p['ID']) == 13:
            muons.append(TLorentzVector(p['Px'], p['Py'], p['Pz'], p['E']))
            lep_q[muons[-1]] = -1 * sign(p['ID'])
            lep_id[muons[-1]] = p['ID']
 

    # Make pairs

    leptons = muons + elecs
    leptons.sort(key=lambda fourvec: fourvec.P(), reverse=True)

    if len(muons) == 4:                         # 4mu
        z1pdg[0], z2pdg[0] = 13, 13
        pair1, pair2, p1_leps, p2_leps = GetMassPairs(leptons, lep_q)
        sel = "4m"

    elif len(elecs) == 4:                       # 4e
        z1pdg[0], z2pdg[0] = 11, 11
        pair1, pair2, p1_leps, p2_leps = GetMassPairs(leptons, lep_q)
        sel = "4e"

    elif len(muons) == 2 and len(elecs) == 2:   # 2m2e
        muon_pair, elec_pair = GetSum(muons), GetSum(elecs)
        sel = "2m2e"

        if muon_pair.M() > elec_pair.M():
            z1pdg[0], z2pdg[0] = 13, 11
            pair1, pair2 = muon_pair, elec_pair
            p1_leps, p2_leps = muons, elecs

        else:
            z1pdg[0], z2pdg[0] = 11, 13
            pair1, pair2 = elec_pair, muon_pair
            p1_leps, p2_leps = elecs, muons

    else:   # we have a problem...
        print("ERROR: wrong number of leptons")
        continue

    system = GetSum([pair1, pair2])

    hPhaseSpaceEvents.Fill(channels[sel])
    hPhaseSpaceEvents.Fill(channels["4l"])



    ##
    ##  FIDUCIAL CUTS
    ##

    is_fiducial = True

    if leptons[0].Pt() < FID_PT1_MIN:
        is_fiducial = False

    if leptons[1].Pt() < FID_PT2_MIN:
        is_fiducial = False

    if leptons[2].Pt() < FID_PT_MIN:
        is_fiducial = False

    if leptons[3].Pt() < FID_PT_MIN:
        is_fiducial = False

    for lep in leptons:
        if abs(lep.Eta()) > FID_ETA_MAX:
            is_fiducal = False

    if is_fiducial:
        hFiducialEvents.Fill(channels[sel])
        hFiducialEvents.Fill(channels["4l"])



    ##
    ##  FILL TREE
    ##

    # Annoying arrays
    evtNum[0] = evt_num
    channel[0] = channels[sel]
    isFiducial[0] = is_fiducial


    # Momenta
    # (simply writing "p4 = _p4[0]" seems to produce a seg fault...)

    zzp4.SetPxPyPzE(system.Px(), system.Py(), system.Pz(), system.E())
    z1p4.SetPxPyPzE(pair1.Px(), pair1.Py(), pair1.Pz(), pair1.E())
    z2p4.SetPxPyPzE(pair2.Px(), pair2.Py(), pair2.Pz(), pair2.E())

    l1p4.SetPxPyPzE(leptons[0].Px(), leptons[0].Py(), leptons[0].Pz(), leptons[0].E())
    l2p4.SetPxPyPzE(leptons[1].Px(), leptons[1].Py(), leptons[1].Pz(), leptons[1].E())
    l3p4.SetPxPyPzE(leptons[2].Px(), leptons[2].Py(), leptons[2].Pz(), leptons[2].E())
    l4p4.SetPxPyPzE(leptons[3].Px(), leptons[3].Py(), leptons[3].Pz(), leptons[3].E())


    # PDG codes
    # (remember: charge and ID are stored as dictionaries)

    l1pdg[0] = lep_id[leptons[0]]
    l2pdg[0] = lep_id[leptons[1]]
    l3pdg[0] = lep_id[leptons[2]]
    l4pdg[0] = lep_id[leptons[3]]


    tree[sel].Fill()
    tree["4l"].Fill()
    evt_num += 1

    del event, this_event

print("done!")

print("\nAcceptance:")
for sel in selection:
    if tree[sel].GetEntries() > 0:
        print("\t" + sel + ":", tree[sel].GetEntries("isFiducial") / tree[sel].GetEntries())
    else:
        print("\t" + sel + ": 0")



##
##  WRITE
##

model_name.Write()
x_sec.Write()
if isBSM:
    m_U.Write()
    w_U.Write()
    g_mu.Write()
    g_e.Write()

for sel in selection:
    tree[sel].Write()
hPhaseSpaceEvents.Write()
hFiducialEvents.Write()

outfile.Close()

print("\nWrote trees to", outname)
