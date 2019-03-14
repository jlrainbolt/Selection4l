from __future__ import print_function
from __future__ import division
import sys

import numpy as np

from ROOT import TFile, TH1, TKey

from PlotUtils import *
#from Cuts2017 import *
from Cuts2016 import *



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
##  UNFOLDED DATA
##

prefix = "unfolding"

ufName = prefix + "_" + year + ".root"
ufFile = TFile(ufName, "READ")
print("Opened", ufName)

hnames = ["b_ttm", "b_l1p", "cos_theta_z1", "cos_theta_z2", 
            "angle_z1leps", "angle_z2leps", "angle_z1l2_z2"]
H = len(hnames)


# Get histograms
data, pred = np.empty(H, dtype=T), np.empty(H, dtype=T)
h = 0

for hname in hnames:
    data[h]['4l'] = ufFile.Get(hname + "/" + hname + "_result")
    data[h]['4l'].SetDirectory(0)

    pred[h]['4l'] = ufFile.Get(hname + "/" + hname + "_reco")
    pred[h]['4l'].SetDirectory(0)

    h = h + 1
h = 0

ufFile.Close()
print("Got data histograms")
print("")



##
##  ACC * EFF
##

prefix = "4l"

# Unscaled signal events
zzName = prefix + "_" + year + "_zz_4l.root"
zzFile = TFile(zzName, "READ")
print("Opened", zzName)

axe = np.empty(H, dtype=T)
h = 0

for sel in selection:
    if sel == "4l":
        continue

    for hname in hnames:
        axe[h][sel] = zzFile.Get(sel + "/" + hname + "_zz_4l")
        axe[h][sel].SetDirectory(0)
        axe[h][sel].SetName(hname + "_acc_x_eff")

        h = h + 1
    h = 0

zzFile.Close()


# Phase space events
ps = np.empty(H, dtype=T)
h = 0

psName = prefix + "_" + year + "_phase_space.root"
psFile = TFile(psName, "READ")
print("Opened", psName)

for sel in selection:
    if sel == "4l":
        continue

    for hname in hnames:
        ps[h][sel] = psFile.Get(sel + "/" + hname + "_phase_space")
        ps[h][sel].SetDirectory(0)

        h = h + 1
    h = 0

psFile.Close()

print("Got acc * eff histograms")
print("")



##
##  ADD CHANNELS
##

# Get 4l and 2m2e
for h in range(H):
    for sample in [ps, axe]:
        sample[h]['2m2e'].Add(sample[h]['2e2m'])
        sample[h]['4l'] = sample[h]['2m2e'].Clone()
        sample[h]['4l'].Add(sample[h]['4m'])
        sample[h]['4l'].Add(sample[h]['4e'])



##
##  SCALING
##

for sel in ["4l"]:
    for h in range(H):
        axe[h][sel].Divide(ps[h][sel])
        scale = CAP_K * GAMMA_Z

        for sample in [data, pred]:
            sample[h][sel].Divide(axe[h][sel])
            sample[h][sel].Scale(scale)


# Get ratio
ratio = np.empty(H, dtype=T)

for sel in ["4l"]:
    for h in range(H):
        ratio[h][sel] = data[h][sel].Clone()
        ratio[h][sel].Divide(pred[h][sel])





####
####
####    LOOP OVER DISTS
####
####


for sel in ["4l"]:
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
        v_pred = np.zeros(pred[h][sel].GetNbinsX(), dtype=V)
        for i in range(len(v_pred)):
            v_pred[i]['x'] = pred[h][sel].GetBinLowEdge(i+1)
            v_pred[i]['y'] = pred[h][sel].GetBinContent(i+1)

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

        width = data[h][sel].GetBinWidth(1)

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

        p_pred = ax_top.bar(        v_pred['x'],        v_pred['y'],        width,
                                    align = 'edge',     linewidth=0,        color = COLOR['zz_4l']
                                )

        top_min, top_max = ax_top.get_ylim()

        if hnames[h] in ["angle_z1l2_z2", "angle_z2leps", "b_ttm"]:
            top_max = top_max * 1.2
        elif hnames[h] in ["cos_theta_z1", "cos_theta_z2"]:
            top_max = top_max * 1.5
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
                verticalalignment = 'top', transform = ax_top.transAxes, fontname = "Helvetica")
        ax_top.set_title(r'\Large{' + lumi + '\,fb$^{-1}$ (13\,TeV, ' + YEAR_STR + ')}',
                loc='right', fontname = "Helvetica")

        # Shared x axis
        xtitle = '$' + data[h][sel].GetXaxis().GetTitle() + '$'
        if "Delta" in xtitle:
            xtitle = xtitle.replace("Delta", "bigtriangleup")
        ax_bot.set_xlabel(xtitle, horizontalalignment='right')
        ax_bot.xaxis.set_label_coords(1, -0.3)

        # Top y axis
        ytitle = r'$d\Gamma / d$' + xtitle
        ytitle = ytitle.replace("$$", "")
        ytitle = ytitle.replace("(", "$ (keV$/" + '%g' % width + r"\ ")
        if hnames[h] in ["cos_theta_z1", "cos_theta_z2"]:
            ytitle = ytitle + r" (keV$/$ " + '%g' % width + " units)"
        ax_top.set_ylabel(ytitle, horizontalalignment='right')
        ax_top.yaxis.set_label_coords(-0.08, 1)
        ax_top.minorticks_on()

        # Bottom y axis
        ax_bot.set_ylabel(r'Data$/$MC')
        ax_bot.yaxis.set_label_coords(-0.08, 0.5)

        

        ##
        ##  TICKS
        ##

        # x axes
        plt.xlim(v_pred['x'][0], v_pred['x'][-1] + width)

        major_step, minor_step = 2 * width, width
        if sel == "4e":
            major_step = width

        for ax in [ax_bot.xaxis, ax_top.xaxis]:
            ax.set_ticks( np.arange(
                            v_pred['x'][0],
                            v_pred['x'][-1] + major_step,
                            step = major_step)
                            )
            ax.set_ticks( np.arange(
                            v_pred['x'][0],
                            v_pred['x'][-1] + minor_step,
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

        if hnames[h] in ["angle_z1leps", "b_l1p"]:
            leg_loc = 'center left'
#       elif hnames[h] in ["cos_theta_z1", "cos_theta_z2"]:
#           leg_loc = 'upper center'
        else:
            leg_loc = 'upper right'

        if year == "2017" and hnames[h] == "zzm" and sel == "4e":
            leg_loc = 'upper right'

        ax_top.legend(
                (   p_data,     p_pred, ),
                (   'Measured', 'Predicted',
                    ),
                loc = leg_loc, numpoints = 1, frameon = False)

        fig.savefig(year + "_" + hnames[h] + "_" + sel + ".pdf")
        plt.clf()
