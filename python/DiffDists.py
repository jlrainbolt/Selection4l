from __future__ import print_function
from __future__ import division
import sys

import numpy as np

from ROOT import TFile, TH1, TKey

from PlotUtils import *
#from Cuts2018 import *
from Cuts2017 import *
#from Cuts2016 import *
#from Cuts2012 import *


mpl.rcParams["legend.fontsize"] = "x-large"

##
##  SAMPLE INFO
##

selection = ["4l"]

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

#hnames = ["cos_theta_z2"]
hnames = ["b_ttm", "b_l1p", "cos_theta_z1", "cos_theta_z2", 
            "angle_z1leps", "angle_z2leps", "angle_z1l2_z2"]
H = len(hnames)

# Single-parameter result
infile = "combination_1.npz"
npzfile = np.load(infile)

alpha, bf_pred = npzfile['alpha_total'], np.sum(npzfile['bf_pred'])
delta_syst = npzfile['delta_syst']


# Get histograms
data, pred, stat = np.empty(H, dtype=T), np.empty(H, dtype=T), np.empty(H, dtype=T)
h = 0

for hname in hnames:
    data[h]['4l'] = ufFile.Get(hname + "/" + hname + "_result")
    data[h]['4l'].SetDirectory(0)

    stat[h]['4l'] = ufFile.Get(hname + "/" + hname + "_stat")
    stat[h]['4l'].SetDirectory(0)

    pred[h]['4l'] = ufFile.Get(hname + "/" + hname + "_gen")
    pred[h]['4l'].SetDirectory(0)

    h = h + 1
h = 0

ufFile.Close()
print("Got data histograms")
print("")



##
##  ACC * EFF
##

prefix = "migration"

# Unscaled signal events
zzName = prefix + "_" + year + "_zz_4l.root"
zzFile = TFile(zzName, "READ")
print("Opened", zzName)

axe = np.empty(H, dtype=T)
h = 0

for sel in selection:
    for hname in hnames:
        axe[h][sel] = zzFile.Get(sel + "/" + hname + "_gen")
        axe[h][sel].SetDirectory(0)
        axe[h][sel].SetName(hname + "_acc_x_eff")

        h = h + 1
    h = 0

zzFile.Close()


# Phase space events
ps = np.empty(H, dtype=T)
h = 0

prefix = "4l"

psName = prefix + "_" + year + "_phase_space.root"
psFile = TFile(psName, "READ")
print("Opened", psName)

sf = INT_LUMI * 1000 * XSEC['zz_4l'] / NGEN['zz_4l']

for sel in selection:
    for hname in hnames:
        ps[h][sel] = psFile.Get(sel + "/" + hname + "_phase_space")
        ps[h][sel].SetDirectory(0)
        ps[h][sel].Scale(sf)

        h = h + 1
    h = 0

psFile.Close()

print("Got acc * eff histograms")
print("")



##
##  SCALING
##

scale = alpha * bf_pred * GAMMA_Z

for sel in ["4l"]:
    for h in range(H):
        axe[h][sel].Divide(ps[h][sel])

        for sample in [data, pred, stat]:
            sample[h][sel].Divide(axe[h][sel])
            sample[h][sel].Scale(scale / sample[h][sel].Integral())

        ps[h][sel].Scale(scale / ps[h][sel].Integral())
#        pred[h][sel].Scale(scale / pred[h][sel].Integral())


# Systemtatic uncertainty
for sel in ["4l"]:
    for h in range(H):

        # Add total systematic
        for i in range(data[h][sel].GetNbinsX()):
            data[h][sel].SetBinError(i + 1,
                    np.sqrt((data[h][sel].GetBinContent(i + 1) * delta_syst) ** 2
                        + data[h][sel].GetBinError(i + 1) ** 2))

        # Get rid of the prediction uncertainty
        for i in range(pred[h][sel].GetNbinsX()):
            pred[h][sel].SetBinError(i+1, 0)
            ps[h][sel].SetBinError(i+1, 0)


# Get ratio
ratio, ratio_stat = np.empty(H, dtype=T), np.empty(H, dtype=T)

for sel in ["4l"]:
    for h in range(H):
        ratio_stat[h][sel] = stat[h][sel].Clone()
#       ratio_stat[h][sel].Divide(pred[h][sel])
        ratio_stat[h][sel].Divide(ps[h][sel])

        ratio[h][sel] = data[h][sel].Clone()
#       ratio[h][sel].Divide(pred[h][sel])
        ratio[h][sel].Divide(ps[h][sel])





####
####
####    LOOP OVER DISTS
####
####


for sel in ["4l"]:
    lumi = '%.1f' % INT_LUMI
    sqrt_s = '%i' % SQRT_S


    print("Drawing", sel, "plots...")

    for h in range(H):

        ##
        ##  GET BIN CONTENT
        ##

        # Data
        v_data = np.zeros(data[h][sel].GetNbinsX(), dtype=V)
        v_stat = np.zeros(data[h][sel].GetNbinsX(), dtype=V)
        for i in range(len(v_data)):
            v_data[i]['x']  = data[h][sel].GetBinCenter(i+1)
            v_data[i]['y']  = data[h][sel].GetBinContent(i+1)
            v_data[i]['ey'] = data[h][sel].GetBinError(i+1)
            v_stat[i]['ey'] = stat[h][sel].GetBinError(i+1)

        # MC
        v_pred = np.zeros(ps[h][sel].GetNbinsX(), dtype=V)
        for i in range(len(v_pred)):
#           v_pred[i]['x']  = pred[h][sel].GetBinLowEdge(i+1)
#           v_pred[i]['y']  = pred[h][sel].GetBinContent(i+1)
#           v_pred[i]['ey'] = pred[h][sel].GetBinContent(i+1)
            v_pred[i]['x']  = ps[h][sel].GetBinLowEdge(i+1)
            v_pred[i]['y']  = ps[h][sel].GetBinContent(i+1)
            v_pred[i]['ey'] = ps[h][sel].GetBinContent(i+1)

        # Ratio
        v_ratio = np.zeros(ratio[h][sel].GetNbinsX(), dtype=V)
        v_ratio_stat = np.zeros(ratio[h][sel].GetNbinsX(), dtype=V)
        for i in range(len(v_ratio)):
            v_ratio[i]['x']     = ratio[h][sel].GetBinCenter(i+1)
            v_ratio[i]['ex']    = ratio[h][sel].GetBinWidth(i+1) / 2
            v_ratio[i]['y']     = ratio[h][sel].GetBinContent(i+1)
            v_ratio[i]['ey']    = ratio[h][sel].GetBinError(i+1)
            v_ratio_stat[i]['ey'] = ratio_stat[h][sel].GetBinError(i+1)



        ##
        ##  MAKE PLOTS
        ##

        width = data[h][sel].GetBinWidth(1)

        fig, (ax_top, ax_bot) = plt.subplots(2, sharex = True, gridspec_kw = lRatioGridSpec)
        fig.subplots_adjust(left = lLeftMargin, right = lRightMargin,   bottom = lBottomMargin,
                            top = lTopMargin,   hspace = lHorizSpace
                            )

        # Top plots
        p_pred = ax_top.errorbar(   v_data['x'],    v_pred['y'],    xerr = v_ratio['ex'], 
                            linewidth = 0,  ecolor = lBlue,
                            fmt = 'None',   capsize = lCapSize,
                            elinewidth = 2 * lErrorLineWidth4l
                            )
        ax_top.errorbar(    v_data['x'],    v_data['y'],            yerr = v_data['ey'],
                            linewidth = 0,  ecolor = '#C0C0C0',     elinewidth = 4 * lErrorLineWidth4l,
                            marker = None,   capsize = 0
                            )
        p_data = ax_top.errorbar(   v_data['x'],    v_data['y'],    yerr = v_stat['ey'], 
                            linewidth = 0,  ecolor = lMarkerColor,  elinewidth = lErrorLineWidth4l,
                            marker = 'o',   capsize = lCapSize,     markersize = lMarkerSize2l,
                            markeredgecolor = lMarkerColor,         markerfacecolor = lMarkerColor
                            )

        top_min, top_max = ax_top.get_ylim()

        if hnames[h] in ["b_ttm", "angle_z1l2_z2"]:
            top_max = 3.5
        elif hnames[h] == "b_l1p":
            top_max = 5.5
        elif hnames[h] == "angle_z1leps":
            top_max = 8
        elif hnames[h] in ["angle_z2leps", "cos_theta_z1"]:
            top_max = 4
        elif hnames[h] == "cos_theta_z2":
            top_max = 3

        ax_top.set_ylim(0, top_max)


        # Ratio plot

        ax_bot.set_ylim(lRatioMin4l, lRatioMax4l)

        ax_bot.axhline(lRatioMid,   color = lBlue,     linewidth = 2 * lErrorLineWidth4l)
        ax_bot.errorbar(v_ratio['x'],   v_ratio['y'],       yerr = v_ratio['ey'],
                    linewidth = 0,  ecolor = '#C0C0C0',     elinewidth = 4 * lErrorLineWidth4l,
                    marker = None,   capsize = 0
                    )
        ax_bot.errorbar(v_ratio['x'],   v_ratio['y'],       yerr = v_ratio_stat['ey'],
                    linewidth = 0,  ecolor = lMarkerColor,  elinewidth = lErrorLineWidth4l,
                    marker = 'o',   capsize = lCapSize,     markersize = lMarkerSize2l,
                    markeredgecolor = lMarkerColor,         markerfacecolor = lMarkerColor
                    )



        ##
        ##  LABELS
        ##

        # Titles
        ax_top.text(0.025,  0.95,   "CMS",
                size = "xx-large",  weight = "bold",
                verticalalignment = 'top', transform = ax_top.transAxes, usetex = False)
        ax_top.text(0.025,  0.875,  "Work in Progress",
                size = "x-large",   style = "italic",
                verticalalignment = 'top', transform = ax_top.transAxes, usetex = False)
        ax_top.set_title(r'\Large{' + lumi + r'\,fb$^{-1}$ (' + sqrt_s + r'\,TeV, ' + YEAR_STR + ')}',
                loc='right')

        # Shared x axis
        if hnames[h] == "b_l1p":
            xtitle = r"$p_{\ell_{1}}$ (GeV)"
            ytitle = r"$d\Gamma/dp_{\ell_{1}}$ (keV$/$GeV)"
        elif hnames[h] == "b_ttm":
            xtitle = r"$m_{\ell_{2,3,4}}$ (GeV)"
            ytitle = r"$d\Gamma/dm_{\ell_{2,3,4}}$ (keV$/$GeV)"
        elif hnames[h] == "angle_z1l2_z2":
            xtitle = r"$\beta$ ($\pi$ rad)"
            ytitle = r"$d\Gamma/d\beta$ (keV$/\pi$ rad)"
        elif hnames[h] == "angle_z1leps":
            xtitle = r"$\alpha_{\mathrm{Z}_{1}}$ ($\pi$ rad)"
            ytitle = r"$d\Gamma/d\alpha_{\mathrm{Z}_{1}}$ (keV$/\pi$ rad)"
        elif hnames[h] == "angle_z2leps":
            xtitle = r"$\alpha_{\mathrm{Z}_{2}}$ ($\pi$ rad)"
            ytitle = r"$d\Gamma/d\alpha_{\mathrm{Z}_{2}}$ (keV$/\pi$ rad)"
        elif hnames[h] == "cos_theta_z1":
            xtitle = r"$\cos\theta_{\mathrm{Z}_{1}}$"
            ytitle = r"$d\Gamma/d\cos\theta_{\mathrm{Z}_{1}}$ (keV$/$unit)"
        elif hnames[h] == "cos_theta_z2":
            xtitle = r"$\cos\theta_{\mathrm{Z}_{2}}$"
            ytitle = r"$d\Gamma/d\cos\theta_{\mathrm{Z}_{2}}$ (keV$/$unit)"
        ax_bot.set_xlabel(xtitle, horizontalalignment='right')
        ax_bot.xaxis.set_label_coords(1, -0.3)

        # Top y axis
        ax_top.set_ylabel(ytitle, horizontalalignment='right')
        if hnames[h] in ["b_ttm", "b_l1p", "angle_z1leps", "angle_z2leps", "cos_theta_z1", "cos_theta_z2"]:    
            ax_top.yaxis.set_label_coords(-0.065, 1)
        else:
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
        ax_bot.yaxis.set_ticks( np.arange(lRatioMin4l+0.5, lRatioMax4l, step = 0.5) )



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
                (   'Measured', 'POWHEG',
                    ),
                loc = leg_loc, numpoints = 1, frameon = False)

        fig.savefig(year + "_" + hnames[h] + "_ddr.pdf")
        plt.clf()
