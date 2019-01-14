from __future__ import print_function
from __future__ import division

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import rc

from ROOT import TFile, TH1



##
## CONSTANTS (FIXME)
##

MU_SUFF = "muon_2017"
EL_SUFF = "electron_2017"
MC_SUFF = [ "zz_4l", "zjets_m-50", "ggH_zz_4l", "vbfH_zz_4l", "ttbar", "ww_2l2nu", "wz_2l2q",
            "wz_3lnu", "zz_2l2q"
            ]
MUON_TRIG_LUMI = 36.735
ELEC_TRIG_LUMI = 41.529
ELEC_TRIG_SF = 0.991
NGEN = [    6.89477e+06, 1.85093e+07, 990084, 234457, 5.64074e+07, 1.99229e+06, 1.65162e+07, 
        6.79248e+06, 1.77959e+07
        ]
XSEC = [1.212, 5765.4, 0.01212, 0.001034, 831.76, 12.178, 5.595, 4.42965, 3.22]

# Real transparency
lAlpha = 0.75
lBlue = (0, 0.4470, 0.7410, lAlpha)
lOrange = (0.8500, 0.3250, 0.0980, lAlpha)
lYellow = (0.9290, 0.6940, 0.1250, lAlpha)
lPurple = (0.4940, 0.1840, 0.5560, lAlpha)
lGreen = (0.4660, 0.6740, 0.1880, lAlpha)
lLightBlue = (0.3010, 0.7450, 0.9330, lAlpha)
lRed = (0.6350, 0.0780, 0.1840, lAlpha)

# Fake transparency
#lBlue = (0.25, 0.58525, 0.80575)
#lOrange = (0.8875, 0.49375, 0.3235)
#lYellow = (0.94675, 0.7705, 0.34375)
#lPurple = (0.6205, 0.388, 0.667)
#lGreen = (0.5995, 0.7555, 0.391)
#lLightBlue = (0.47575, 0.80875, 0.94975)
#lRed = (0.72625, 0.3085, 0.388)

COLOR = [lLightBlue, lYellow, lPurple, lPurple, lGreen, lOrange, lOrange, lOrange, lOrange]


##
##  SAMPLE INFO
##

selection = ["4m", "2m2e", "2e2m", "4e"]



##
##  LOAD DATA
##

prefix = "../SMPV/unscaled4l"

# Muon file
muName = prefix + "_" + MU_SUFF + ".root"
muFile = TFile(muName, "READ")
print("Opened", muName)

# Electron file
elName = prefix + "_" + EL_SUFF + ".root"
elFile = TFile(elName, "READ")
print("Opened", elName)

# Get histograms
hname = "zzm"

hdata = []
for sel in selection:
    if sel == "4m" or sel == "2m2e":
        hist = muFile.Get(sel + "/" + hname + "_" + MU_SUFF)
    elif sel == "2e2m" or sel == "4e":
        hist = elFile.Get(sel + "/" + hname + "_" + EL_SUFF)
    hist.SetDirectory(0)
    hdata.append(hist)

muFile.Close()
elFile.Close()
print("Got data histograms")
print("")

# Add channels #FIXME
data = hdata[0]
for i in range(1, len(hdata)):
    data.Add(hdata[i])



##
##  LOAD MC
##

mc = []

# Loop over all samples
for suff, xsec, ngen in zip(MC_SUFF, XSEC, NGEN):
    inName = prefix + "_" + suff + ".root"
    inFile = TFile(inName, "READ")
    print("Opened", inName)

    # Get histograms
    hmc = []
    for sel in selection:
        hist = inFile.Get(sel + "/" + hname + "_" + suff)
        hist.SetDirectory(0)

        if sel == "4m" or sel == "2m2e":
            lumi = MUON_TRIG_LUMI
        elif sel == "2e2m" or sel == "4e":
            lumi = ELEC_TRIG_LUMI * ELEC_TRIG_SF
        sf = lumi * 1000 * xsec / ngen

        hist.Scale(sf)
        hmc.append(hist)
    inFile.Close()

    # Add channels
    h4l = hmc[0]
    for i in range(1, len(hmc)):
        h4l.Add(hmc[i])
    mc.append(h4l)

print("Got MC histograms")
print("")



##
##  RATIO PLOT
##

total = mc[0].Clone()
for i in range(1, len(mc)):
    total.Add(mc[i])
ratio = data.Clone()
ratio.Divide(total)



##
##  GET BIN CONTENT
##

# Data
x_data, y_data, yerr_data = [], [], []
for i in range(1, data.GetNbinsX()+1):
    x_data.append(data.GetBinCenter(i))
    y_data.append(data.GetBinContent(i))
    yerr_data.append(data.GetBinError(i))

# MC
x_mc, y_mc = [], []
for i in range(1, total.GetNbinsX()+1):
    x_mc.append(total.GetBinLowEdge(i))
for hist in mc:
    y_mc_ = []
    for i in range(1, hist.GetNbinsX()+1):
        binc = hist.GetBinContent(i)
        if binc < 0:
            binc = 0
        y_mc_.append(binc)
    y_mc.append(y_mc_)

# Bottoms
# FIXME use numpy
bot_mc = [[0] * len(x_mc)]
for j in range(0, len(y_mc)):
    bot_mc_ = []
    for i in range(0, len(x_mc)):
        bot_mc_.append(bot_mc[j][i] + y_mc[j][i])
    bot_mc.append(bot_mc_)

# Ratio
y_ratio, yerr_ratio, xerr_ratio = [], [], []
for i in range(1, ratio.GetNbinsX()+1):
    y_ratio.append(ratio.GetBinContent(i))
    yerr_ratio.append(ratio.GetBinError(i))
    xerr_ratio.append(ratio.GetBinWidth(i)/2)



####
####
####    MAKE PLOTS
####
####


##
##  SETUP
##

# Font/TeX setup
# (https://3diagramsperpage.wordpress.com/2015/04/11/)
rc('text', usetex=True)
mpl.rcParams['text.latex.preamble'] = [ r'\usepackage{helvet}', r'\usepackage{sansmath}', 
                                        r'\sansmath'] 
# Figure size, aspect ratio
mpl.rcParams["figure.figsize"] = [6, 6]

# Subplot
fig, (ax_top, ax_bot) = plt.subplots(2, sharex=True, gridspec_kw={'height_ratios':[3, 1]})
fig.subplots_adjust(hspace=0)



##
##  TOP PLOTS
##

# Data
p_data = ax_top.errorbar(x_data, y_data, yerr=yerr_data, 
        marker='o', markersize=10, markeredgecolor='black', markerfacecolor='black',
        linewidth=0, ecolor='black', elinewidth=2, capsize=0) 

# MC
p_mc = []
for i in range(0, len(mc)):
    p_mc.append(ax_top.bar(x_mc, y_mc[i], 1.0, align='edge', color=COLOR[i], linewidth=0,
        bottom=bot_mc[i]))



##
##  BOTTOM PLOTS
##

# Ratio plot
ax_bot.errorbar(x_data, y_ratio, yerr=yerr_ratio, xerr=xerr_ratio, 
        marker='o', markersize=10, markeredgecolor='black', markerfacecolor='black',
        linewidth=0, ecolor='black', elinewidth=2, capsize=0) 

# Horizontal lines
ax_bot.axhline(1, color='black')
ax_bot.axhline(0.8, color='black', linestyle=':')
ax_bot.axhline(1.2, color='black', linestyle=':')

# Vertical range
ax_bot.set_ylim(0.5, 1.5)



##
##  LABELS
##

# Titles
ax_top.set_title(r'\textbf{CMS} \textit{Work in Progress}', loc='left')
ax_top.set_title(r'41.5\,fb$^{-1}$ (13\,TeV)', loc='right')

# Top y axis
ax_top.set_ylabel(r'Events $/$ GeV', horizontalalignment='right')
ax_top.yaxis.set_label_coords(-0.075, 1)

# Bottom y axis
ax_bot.set_ylabel(r'Data $/$ MC')
ax_bot.yaxis.set_label_coords(-0.075, 0.5)

# Shared x axis
plt.xlabel(r'$m_{4\ell}$ (GeV)', horizontalalignment='right', x=1)
plt.xticks(np.arange(x_mc[0], x_mc[-1]+2, step=2))



##
##  LEGEND
##

ax_top.legend((p_data, p_mc[0], p_mc[1], p_mc[2], p_mc[4], p_mc[5]),
        (r'Data', r'$\mbox{ZZ} \to 4\ell$', r'$\mbox{Z} \to \ell^{+} \ell^{-}$', 
        r'H', r'$\mbox{t}\bar{\mbox{t}}$', r'VV'),
        loc=2, numpoints=1, frameon=False)



##
##  WRITE TO FILE
##

fig.savefig('mass.pdf')
