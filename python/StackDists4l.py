from __future__ import print_function
from __future__ import division

import numpy as np

from ROOT import TFile, TH1, TKey

from PlotUtils import *
from Cuts2017 import *



##
##  SAMPLE INFO
##

L4, M4, ME, EM, E4  = 0, 1, 2, 3, 4
selection = ["4l", "4m", "2m2e", "2e2m", "4e"]
N = len(selection)



##
##  LOAD DATA
##

prefix = "~/nobackup/Selection2017/SMPV/unscaled4l"

# Muon file
muName = prefix + "_" + MU_SUFF + ".root"
muFile = TFile(muName, "READ")
print("Opened", muName)

# Electron file
elName = prefix + "_" + EL_SUFF + ".root"
elFile = TFile(elName, "READ")
print("Opened", elName)



##
##  GET KEYS
##

keyDir = muFile.GetDirectory("/" + selection[M4], True, "GetDirectory")

hnames = []
for key in keyDir.GetListOfKeys():
    hname = key.GetName()
    hnames.append(hname.replace("_" + MU_SUFF, "")

H = len(hnames)



hdata = []
for sel in selection:
    hists = []
    for hname in hnames:
        if sel == "4m" or sel == "2m2e":
            hist = muFile.Get(sel + "/" + hname + "_" + MU_SUFF)
        elif sel == "2e2m" or sel == "4e":
            hist = elFile.Get(sel + "/" + hname + "_" + EL_SUFF)
        hist.SetDirectory(0)
        hists.append(hist)
    hdata.append(hists)

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
        hists = []
        for hname in hnames:
            hist = inFile.Get(sel + "/" + hname + "_" + suff)
            hist.SetDirectory(0)

            if sel == "4m" or sel == "2m2e":
                lumi = MUON_TRIG_LUMI
            elif sel == "2e2m" or sel == "4e":
                lumi = ELEC_TRIG_LUMI * ELEC_TRIG_SF
            sf = lumi * 1000 * xsec / ngen

            hist.Scale(sf)
            hists.append(hist)
        hmc.append(hists)
    inFile.Close()

    # Add channels
    h4l = hmc[0]
    for i in range(1, N_MC):
        h4l.Add(hmc[i])
    mc.append(h4l)

print("Got MC histograms")
print("")





####
####
####    DISTRIBUTION LOOP
####
####



##
##  RATIO PLOT
##

ratios = []

for i in range(N):
    hists = []
    for h in range(H):
        total = mc[0].Clone()

        for j in range(1, N_MC):
            total.Add(mc[i])

        ratio = hdata[i][h].Clone()
        ratio.Divide(total)
        hists.append(ratio)
    ratios.append(hists)



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


fig, (ax_top, ax_bot) = plt.subplots(2, sharex = True, gridspec_kw = lRatioGridSpec)

fig.subplots_adjust(    left = lLeftMargin, right = lRightMargin,   bottom = lBottomMargin,
                        top = lTopMargin,   hspace = lHorizSpace
                        )



##
##  TOP PLOTS
##

# Data
p_data = ax_top.errorbar(   x_data, y_data, yerr = yerr_data, 
                            linewidth = 0,  ecolor = lMarkerColor,  elinewidth = lErrorLineWidth,
                            marker = 'o',   capsize = lCapSize,     markersize = lMarkerSize,
                            markeredgecolor = lMarkerColor,     markerfacecolor = lMarkerColor
                            )

# MC
p_mc = []
for i in range(0, len(mc)):
    p_mc.append(
                    ax_top.bar( x_mc,   y_mc[i],    1.0,            bottom=bot_mc[i],
                                align = 'edge',     linewidth=0,    color = COLOR[i]    )
                    )



##
##  BOTTOM PLOTS
##

# Ratio plot
ax_bot.errorbar(    x_data, y_ratio,    yerr = yerr_ratio,  xerr = xerr_ratio, 
                    linewidth = 0,  ecolor = lMarkerColor,  elinewidth = lErrorLineWidth,
                    marker = 'o',   capsize = lCapSize,     markersize = lMarkerSize,
                    markeredgecolor = lMarkerColor,     markerfacecolor = lMarkerColor
                    )

# Horizontal lines
ax_bot.axhline(lRatioUpper, color = lRatioLineColor, linestyle = ':')
ax_bot.axhline(lRatioMid,   color = lRatioLineColor)
ax_bot.axhline(lRatioLower, color = lRatioLineColor, linestyle = ':')

# Vertical range
ax_bot.set_ylim(lRatioMin, lRatioMax)



##
##  LABELS
##

# Titles
ax_top.set_title(r'\textbf{CMS} \Large{\textit{Work in Progress}}', loc='left')
ax_top.set_title(r'\Large{41.5\,fb$^{-1}$ (13\,TeV)}', loc='right')

# Top y axis
#ytitle = '$' + hdata[0].GetYaxis().GetTitle() + '$'
ytitle = "Events"
ax_top.set_ylabel(ytitle, horizontalalignment='right')
ax_top.yaxis.set_label_coords(-0.08, 1)

# Bottom y axis
ax_bot.set_ylabel(r'Data$/$MC')
ax_bot.yaxis.set_label_coords(-0.08, 0.5)

# Shared x axis
xtitle = '$' + hdata[0].GetXaxis().GetTitle() + '$'
ax_bot.set_xlabel(xtitle, horizontalalignment='right')
ax_bot.xaxis.set_label_coords(1, -0.3)
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
