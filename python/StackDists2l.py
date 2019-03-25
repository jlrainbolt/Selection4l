from __future__ import print_function
from __future__ import division
import sys

import numpy as np

from ROOT import TFile, TH1, TKey

from PlotUtils import *
from Cuts2017 import *
#from Cuts2016 import *



##
##  SAMPLE INFO
##

selection = ["mumu", "ee"]

T = np.dtype([(sel, object) for sel in selection])
V = np.dtype([("x", 'f4'), ("y", 'f4'), ("ex", 'f4'), ("ey", 'f4'), ("b", 'f4')])

year = sys.argv[1]
if year != YEAR_STR:
    print("Wrong year in header file")

#tag = "noQt"


##
##  DATA
##

prefix = "2l"

# Muon file
muName = prefix + "_" + year + "_" + MU_SUFF + ".root"
muFile = TFile(muName, "READ")
print("Opened", muName)

# Electron file
elName = prefix + "_" + year + "_" + EL_SUFF + ".root"
elFile = TFile(elName, "READ")
print("Opened", elName)


# Get keys
#keyDir = muFile.GetDirectory("/mumu", True, "GetDirectory")
#hnames = []
#for key in keyDir.GetListOfKeys():
#    hname = key.GetName()
#    hnames.append(hname.replace("_" + MU_SUFF, ""))
hnames = ["z1m", "z1pt", "l1pt", "l2pt", "l1eta", "l2eta"]#, "dphi", "z1y", "nPV"]
H = len(hnames)


# get histograms
data = np.empty(H, dtype=T)
h = 0

for hname in hnames:
    data[h]['mumu'] = muFile.Get("mumu/" + hname + "_" + MU_SUFF)
    data[h]['ee']   = elFile.Get("ee/" + hname + "_" + EL_SUFF)

    for sel in selection:
        data[h][sel].SetDirectory(0)

    h = h + 1

muFile.Close()
elFile.Close()
print("Got data histograms")
print("")



##
##  MONTE CARLO
##

mc_arr = np.empty((N_MC, H), dtype=T)
mc = {}
h, j = 0, 0

# Loop over all samples
for suff in MC_SUFF_2L:
    inName = prefix + "_" + year + "_" + suff + ".root"
    inFile = TFile.Open(inName)
    print("Opened", inName)

    # Get histograms
    for sel in selection:
        if sel == "mumu":
            lumi = MUON_TRIG_LUMI
        elif sel == "ee":
            lumi = ELEC_TRIG_LUMI * ELEC_TRIG_SF

        sf = lumi * 1000 * XSEC[suff] / NGEN[suff]

        for hname in hnames:
            mc_arr[j][h][sel] = inFile.Get(sel + "/" + hname + "_" + suff)
            mc_arr[j][h][sel].SetDirectory(0)
            mc_arr[j][h][sel].Scale(sf)
            h = h + 1
        h = 0

    mc[suff] = mc_arr[j]
    j = j + 1
    inFile.Close()

print("Got MC histograms")
print("")


# Get total
total, ratio = np.empty(H, dtype=T), np.empty(H, dtype=T)

for sel in selection:
    for h in range(H):
        for suff in MC_SUFF_2L:
            if suff == "zjets_m-50":
                total[h][sel] = mc[suff][h][sel].Clone()
            else:
                total[h][sel].Add(mc[suff][h][sel])

        ratio[h][sel] = data[h][sel].Clone()
        ratio[h][sel].Divide(total[h][sel])



####
####
####    LOOP OVER STACKS
####
####


for sel in selection:
    print("Drawing", sel, "plots...")

    if sel == "mumu":
        lumi = MUON_TRIG_LUMI
    elif sel == "ee":
        lumi = ELEC_TRIG_LUMI

    for h in range(H):

        ##
        ##  GET BIN CONTENT
        ##

        # Data
        v_data = np.zeros(data[h][sel].GetNbinsX(), dtype=V)
        for i in range(len(v_data)):
            v_data[i]['x']  = data[h][sel].GetBinCenter(i+1)
            v_data[i]['y']  = data[h][sel].GetBinContent(i+1)
            v_data[i]['ey'] = data[h][sel].GetBinError(i+1)

        # MC
        v_mc_arr = np.zeros([N_MC, total[h][sel].GetNbinsX()], dtype=V)
        v_mc = {}
        j = 0
        for suff in MC_SUFF_2L:
            for i in range(len(v_mc_arr[0])):
                v_mc_arr[j][i]['x'] = total[h][sel].GetBinLowEdge(i+1)
                v_mc_arr[j][i]['y'] = mc[suff][h][sel].GetBinContent(i+1)

            v_mc[suff] = v_mc_arr[j]
            j = j + 1

        # "Bottoms"
        for j in range(N_MC - 1):
            v_mc_arr[j]['b'] = np.sum(v_mc_arr[j+1:]['y'], axis=0)

        # Ratio
        v_ratio = np.zeros(ratio[h][sel].GetNbinsX(), dtype=V)
        for i in range(len(v_ratio)):
            v_ratio[i]['x']     = ratio[h][sel].GetBinCenter(i+1)
            v_ratio[i]['ex']    = ratio[h][sel].GetBinWidth(i+1) / 2
            v_ratio[i]['y']     = ratio[h][sel].GetBinContent(i+1)
            v_ratio[i]['ey']    = ratio[h][sel].GetBinError(i+1)



        ##
        ##  MAKE PLOTS
        ##

        width = total[h][sel].GetBinWidth(1)

        fig, (ax_top, ax_bot) = plt.subplots(2, sharex = True, gridspec_kw = lRatioGridSpec)
        fig.subplots_adjust(left = lLeftMargin, right = lRightMargin,   bottom = lBottomMargin,
                            top = lTopMargin,   hspace = lHorizSpace
                            )

        # Top plots
        p_data = ax_top.errorbar(   v_data['x'],    v_data['y'],    yerr = v_data['ey'], 
                            linewidth = 0,  ecolor = lMarkerColor,  elinewidth = lErrorLineWidth2l,
                            marker = 'o',   capsize = lCapSize,     markersize = lMarkerSize2l,
                            markeredgecolor = lMarkerColor,         markerfacecolor = lMarkerColor
                            )

        p_mc = {}
        for suff in MC_SUFF_2L:
            p_mc[suff] = ax_top.bar(    v_mc[suff]['x'],    v_mc[suff]['y'],    width,
                                bottom = v_mc[suff]['b'],   align = 'edge',     linewidth=0,
                                color = COLOR[suff]
                                )

        top_min, top_max = ax_top.get_ylim()

        if "eta" in hnames[h]:
            top_max = top_max * 1.7
#       elif "y" in hnames[h]:
#           top_max = top_max * 1.8
        elif "pt" in hnames[h]:
            top_max = top_max * 1.3
#       elif "phi" in hnames[h]:
#           top_max = top_max * 1.2
        elif "z1m" in hnames[h]:
            top_max = top_max * 1.2

        ax_top.set_ylim(0, top_max)


        # Ratio plot
        ax_bot.errorbar(v_ratio['x'],   v_ratio['y'],   xerr = v_ratio['ex'],   yerr = v_ratio['ey'], 
                    linewidth = 0,  ecolor = lMarkerColor,  elinewidth = lErrorLineWidth2l,
                    marker = 'o',   capsize = lCapSize,     markersize = lMarkerSize2l,
                    markeredgecolor = lMarkerColor,         markerfacecolor = lMarkerColor
                    )

        ax_bot.axhline(lRatioMid,   color = lRatioLineColor, linestyle = ':')
        ax_bot.set_ylim(lRatioMin2l, lRatioMax2l)



        ##
        ##  LABELS
        ##

        # Titles
        ax_top.text(    0.025,  0.95,
                r'\LARGE{\textbf{CMS}}' + '\n' + r'\Large{\textit{Work in Progress}}',
                verticalalignment = 'top', transform = ax_top.transAxes)
        ax_top.set_title(r'\Large{' + '%.1f' % lumi + '\,fb$^{-1}$ (13\,TeV, ' + YEAR_STR + ')}',
                loc='right')

        # Top y axis
        unit = '$' + mc['zjets_m-50'][h][sel].GetYaxis().GetTitle() + '$'
        ytitle = r"Events$/$" + '%g' % width + " " + unit
        ax_top.set_ylabel(ytitle, horizontalalignment='right')
        ax_top.yaxis.set_label_coords(-0.08, 1)
        ax_top.minorticks_on()

        # Bottom y axis
        ax_bot.set_ylabel(r'Data$/$MC')
        ax_bot.yaxis.set_label_coords(-0.08, 0.5)

        # Shared x axis
        xtitle = '$' + mc['zjets_m-50'][h][sel].GetXaxis().GetTitle() + '$'
        if "Delta" in xtitle:
            xtitle = xtitle.replace("Delta", "triangle")
        ax_bot.set_xlabel(xtitle, horizontalalignment='right')
        ax_bot.xaxis.set_label_coords(1, -0.3)

        

        ##
        ##  TICKS
        ##

        # x axes
        plt.xlim(v_mc['zjets_m-50']['x'][0], v_mc['zjets_m-50']['x'][-1] + width)

        major_step, minor_step = 5 * width, width
        if not (major_step.is_integer()):
            major_step = 4 * width

        for ax in [ax_bot.xaxis, ax_top.xaxis]:
            ax.set_ticks( np.arange(
                            v_mc['zjets_m-50']['x'][0],
                            v_mc['zjets_m-50']['x'][-1] + major_step,
                            step = major_step)
                            )
            ax.set_ticks( np.arange(
                            v_mc['zjets_m-50']['x'][0],
                            v_mc['zjets_m-50']['x'][-1] + minor_step,
                            step = minor_step),
                        minor = True)

        # Top y axis
        ax_top.ticklabel_format(axis = 'y', style = 'sci')
        ax_top.yaxis.get_major_formatter().set_powerlimits((0, 1))

        # Bottom y axis
        ax_bot.yaxis.set_ticks( np.arange(lRatioMin2l+0.1, lRatioMax2l, step = 0.1) )
        ax_bot.yaxis.set_ticks( np.arange(lRatioMin2l+0.05, lRatioMax2l, step = 0.05),
                        minor = True  )



        ##
        ##  LEGEND
        ##

        if hnames[h] in ["l2pt", "z1m"]:
            leg_loc = 'center left'
        else:
            leg_loc = 'upper right'

        ax_top.legend(
                (   p_data,                         p_mc['zz_4l'],
                    p_mc['zjets_m-50'],             p_mc['ttbar'],
                    p_mc['ww_2l2nu'],               p_mc['zzz_4l2nu'],
                    p_mc['ggH_zz_4l']
                    ),
                (   r'Data',                        r'$\mbox{ZZ}\to4\ell$',
                    r'$\mbox{Z}\to\ell^+\ell^-$',   r'$\mbox{t}\bar{\mbox{t}}$(V)',
                    r'VV',                          r'VVV',
                    r'H'
                    ),
                loc = leg_loc, numpoints = 1, frameon = False)

        fig.savefig(year + "_" + hnames[h] + "_" + sel + ".pdf")
#       fig.savefig(year + "_" + tag + "_" + hnames[h] + "_" + sel + ".pdf")
        plt.clf()
