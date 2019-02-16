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

selection = ["4l", "4m", "2m2e", "2e2m", "4e"]

T = np.dtype([(sel, object) for sel in selection])
V = np.dtype([("x", 'f4'), ("y", 'f4'), ("ex", 'f4'), ("ey", 'f4'), ("b", 'f4')])

year = sys.argv[1]
if year != YEAR_STR:
    print("Wrong year in header file")



##
##  DATA
##

prefix = "4l"

# Muon file
muName = prefix + "_" + year + "_" + MU_SUFF + ".root"
muFile = TFile(muName, "READ")
print("Opened", muName)

# Electron file
elName = prefix + "_" + year + "_" + EL_SUFF + ".root"
elFile = TFile(elName, "READ")
print("Opened", elName)


# Get keys
#keyDir = muFile.GetDirectory("/4m", True, "GetDirectory")

#hnames = []
#for key in keyDir.GetListOfKeys():
#    hname = key.GetName()
#    hnames.append(hname.replace("_" + MU_SUFF, ""))
#hnames = ["zzm", "zzpt", "z1m", "z2m", "sin_phi", "sin_phi_2"]
hnames = ["sin_phi", "sin_phi_2"]

H = len(hnames)


# Get histograms
data = np.empty(H, dtype=T)
h = 0

for sel in selection:
    if sel == "4l":
        continue

    for hname in hnames:
        if sel in ["4m", "2m2e"]:
            data[h][sel] = muFile.Get(sel + "/" + hname + "_" + MU_SUFF)
        elif sel in ["4e", "2e2m"]:
            data[h][sel] = elFile.Get(sel + "/" + hname + "_" + EL_SUFF)

        data[h][sel].SetDirectory(0)
        h = h + 1
    h = 0

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
for suff in MC_SUFF_4L:
    inName = prefix + "_" +year + "_" + suff + ".root"
    inFile = TFile.Open(inName)
    print("Opened", inName)

    # Get histograms
    for sel in selection:
        if sel == "4l":
            continue
        elif sel in ["4m", "2m2e"]:
            lumi = MUON_TRIG_LUMI
        elif sel in ["4e", "2e2m"]:
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



##
##  ADD CHANNELS
##

# Get 4l and 2m2e, rebin 4e
for h in range(H):
    data[h]['2m2e'].Add(data[h]['2e2m'])
    data[h]['4l'] = data[h]['2m2e'].Clone()
    data[h]['4l'].Add(data[h]['4m'])
    data[h]['4l'].Add(data[h]['4e'])

    if hnames[h] != "sin_phi_2":
        data[h]['4e'].Rebin(2)

    for suff in MC_SUFF_4L:
        mc[suff][h]['2m2e'].Add(mc[suff][h]['2e2m'])
        mc[suff][h]['4l'] = mc[suff][h]['2m2e'].Clone()
        mc[suff][h]['4l'].Add(mc[suff][h]['4m'])
        mc[suff][h]['4l'].Add(mc[suff][h]['4e'])

        if hnames[h] != "sin_phi_2":
            mc[suff][h]['4e'].Rebin(2)

# Get total
total, ratio = np.empty(H, dtype=T), np.empty(H, dtype=T)

for sel in selection:
    if sel == "2e2m":
        continue

    for h in range(H):
        for suff in MC_SUFF_4L:
            if suff == MC_SUFF_4L[0]:
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
#for sel in ["4e"]:
    if sel == "2e2m":
        continue

    if (MUON_TRIG_LUMI == ELEC_TRIG_LUMI):
        lumi = '%.1f' % MUON_TRIG_LUMI
    elif sel == "4m":
        lumi = '%.1f' % MUON_TRIG_LUMI
    elif sel == "4e":
        lumi = '%.1f' % ELEC_TRIG_LUMI
    elif sel in ["4l", "2m2e"]:
        lumi = '%.1f' % MUON_TRIG_LUMI + " + " + '%.1f' % ELEC_TRIG_LUMI


    print("Drawing", sel, "plots...")

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
        for suff in MC_SUFF_4L:
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
                            linewidth = 0,  ecolor = lMarkerColor,  elinewidth = lErrorLineWidth4l,
                            marker = 'o',   capsize = lCapSize,     markersize = lMarkerSize4l,
                            markeredgecolor = lMarkerColor,         markerfacecolor = lMarkerColor
                            )

        p_mc = {}
        for suff in MC_SUFF_4L:
            p_mc[suff] = ax_top.bar(    v_mc[suff]['x'],    v_mc[suff]['y'],    width,
                                bottom = v_mc[suff]['b'],   align = 'edge',     linewidth=0,
                                color = COLOR[suff]
                                )

        top_min, top_max = ax_top.get_ylim()

        if hnames[h] == "zzpt":
            top_max = top_max * 1.2
        elif hnames[h] == "z2m":
            top_max = top_max * 1.3
        elif hnames[h] == "sin_phi_2":
            top_max = top_max * 2
        ax_top.set_ylim(0, top_max)


        # Ratio plot
        ax_bot.errorbar(v_ratio['x'],   v_ratio['y'],   xerr = v_ratio['ex'],   yerr = v_ratio['ey'], 
                    linewidth = 0,  ecolor = lMarkerColor,  elinewidth = lErrorLineWidth4l,
                    marker = 'o',   capsize = lCapSize,     markersize = lMarkerSize4l,
                    markeredgecolor = lMarkerColor,         markerfacecolor = lMarkerColor
                    )

        ax_bot.axhline(lRatioMid,   color = lRatioLineColor, linestyle = ':')
        if sel == "4e":
            ax_bot.set_ylim(0, 3)
        else:
            ax_bot.set_ylim(lRatioMin4l, lRatioMax4l)



        ##
        ##  LABELS
        ##

        # Titles
        ax_top.text(    0.025,  0.95,
                r'\LARGE{\textbf{CMS}}' + '\n' + r'\Large{\textit{Work in Progress}}',
                verticalalignment = 'top', transform = ax_top.transAxes)
        ax_top.set_title(r'\Large{' + lumi + '\,fb$^{-1}$ (13\,TeV, ' + YEAR_STR + ')}',
                loc='right')

        # Top y axis
        if sel == "4e":
            ytitle = r'\mbox{Events / }' + '%.0f' % width + r'\mbox{ }\mbox{GeV}'
        else:
            ytitle = '$' + mc['zz_4l'][h][sel].GetYaxis().GetTitle() + '$'
        ax_top.set_ylabel(ytitle, horizontalalignment='right')
        ax_top.yaxis.set_label_coords(-0.08, 1)
        ax_top.minorticks_on()

        # Bottom y axis
        ax_bot.set_ylabel(r'Data$/$MC')
        ax_bot.yaxis.set_label_coords(-0.08, 0.5)

        # Shared x axis
        xtitle = '$' + mc['zz_4l'][h][sel].GetXaxis().GetTitle() + '$'
        if "Delta" in xtitle:
            xtitle = xtitle.replace("Delta", "bigtriangleup")
        ax_bot.set_xlabel(xtitle, horizontalalignment='right')
        ax_bot.xaxis.set_label_coords(1, -0.3)

        

        ##
        ##  TICKS
        ##

        # x axes
        plt.xlim(v_mc['zz_4l']['x'][0], v_mc['zz_4l']['x'][-1] + width)

        major_step, minor_step = 2 * width, width
        if sel == "4e":
            major_step = width

        for ax in [ax_bot.xaxis, ax_top.xaxis]:
            ax.set_ticks( np.arange(
                            v_mc['zz_4l']['x'][0],
                            v_mc['zz_4l']['x'][-1] + major_step,
                            step = major_step)
                            )
            ax.set_ticks( np.arange(
                            v_mc['zz_4l']['x'][0],
                            v_mc['zz_4l']['x'][-1] + minor_step,
                            step = minor_step),
                        minor = True)

        # Top y axis
#       ax_top.ticklabel_format(axis = 'y', style = 'sci')
#       ax_top.yaxis.get_major_formatter().set_powerlimits((0, 1))

        # Bottom y axis
        if sel == "4e":
            ax_bot.yaxis.set_ticks( np.arange(lRatioMid, lRatioMid+2, step = 1) )
        else:
            ax_bot.yaxis.set_ticks( np.arange(lRatioMin4l+0.5, lRatioMax4l, step = 0.5) )
#       ax_bot.yaxis.set_ticks( np.arange(lRatioMin2l+0.05, lRatioMax4l, step = 0.05),
#                       minor = True  )



        ##
        ##  LEGEND
        ##

        if hnames[h] in ["zzm", "z1m"]:
            leg_loc = 'center left'
        elif hnames[h] in ["sin_phi", "sin_phi_2"]:
            leg_loc = 'upper center'
        else:
            leg_loc = 'upper right'

        if year == "2017" and hnames[h] == "zzm" and sel == "4e":
            leg_loc = 'upper right'

        ax_top.legend(
                (   p_data,                         p_mc['zz_4l'],
                    p_mc['zjets_m-50'],             p_mc['ttbar'],
                    p_mc['ww_2l2nu'],               p_mc['ggH_zz_4l']
                    ),
                (   r'Data',                        r'$\mbox{ZZ}\to4\ell$',
                    r'$\mbox{Z}\to\ell^+\ell^-$',   r'$\mbox{t}\bar{\mbox{t}}$', 
                    r'VV',                          r'H'
                    ),
                loc = leg_loc, numpoints = 1, frameon = False)

        fig.savefig(year + "_" + hnames[h] + "_" + sel + ".pdf")
        plt.clf()
