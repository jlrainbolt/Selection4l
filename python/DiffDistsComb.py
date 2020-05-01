from __future__ import print_function
from __future__ import division
import sys

import numpy as np

from ROOT import TFile, TH1

from PlotUtils import *
from Cuts2018 import *

mpl.rcParams["legend.fontsize"] = "x-large"


##
##  SAMPLE INFO
##

selection = ["4l"]

nameTeX = {"b_z1m":r"m_{\mathrm{Z}_{1}}", "b_z2m":r"m_{\mathrm{Z}_{2}}",
        "b_l1p":r"p_{\ell_{1}}", "b_ttm":r"m_{\ell_{2,3,4}}",
        "angle_z1l2_z2":r"\beta", "angle_z1leps":r"\alpha_{\mathrm{Z}_{1}}",
        "angle_z2leps":r"\alpha_{\mathrm{Z}_{2}}",
        "cos_theta_z1":r"\cos\theta_{\mathrm{Z}_{1}}",
        "cos_theta_z2":r"\cos\theta_{\mathrm{Z}_{2}}",
        "sin_phi":r"\sin\phi", "sin_phi_10":r"\sin\phi"}

T = np.dtype([(sel, object) for sel in selection])
V = np.dtype([("x", 'f4'), ("y", 'f4'), ("ex", 'f4'), ("eyu", 'f4'), ("eyd", 'f4'), ("b", 'f4')])

hnames  = ["b_z1m", "b_z2m", "b_l1p", "b_ttm", "cos_theta_z1", "cos_theta_z2",
                "angle_z1leps", "angle_z2leps", "angle_z1l2_z2", "sin_phi"]
H = len(hnames)

# Single-parameter result
infile = "combination_1.npz"
npzfile = np.load(infile)

alpha, bf_pred = npzfile['alpha_total'], np.sum(npzfile['bf_pred'])
delta_syst = npzfile['delta_syst']


##
##  UNFOLDED DATA
##

prefix = "unfolding"

ufName = prefix + "_comb.root"
ufFile = TFile(ufName, "READ")
print("Opened", ufName)

data, axe, stat = np.empty(H, dtype=T), np.empty(H, dtype=T), np.empty(H, dtype=T)

h = 0
for hname in hnames:
    data[h]['4l'] = ufFile.Get(hname + "/" + hname + "_result")
    data[h]['4l'].SetDirectory(0)

    stat[h]['4l'] = ufFile.Get(hname + "/" + hname + "_stat")
    stat[h]['4l'].SetDirectory(0)

    h = h + 1
ufFile.Close()

print("Got data histograms")
print("")



##
##  ACC * EFF
##

prefix = "4l"

# Unscaled signal events
zzName = prefix + "_" + YEAR_STR + "_zz_4l.root"
zzFile = TFile(zzName, "READ")
print("Opened", zzName)

axe = np.empty(H, dtype=T)

h = 0
for hname in hnames:
    axe[h]['4l'] = zzFile.Get("4l/" + hname + "_zz_4l")
    axe[h]['4l'].SetDirectory(0)
    axe[h]['4l'].SetName(hname + "_acc_x_eff")

    h = h + 1
h = 0
zzFile.Close()


for year in ["2017", "2016", "2012"]:
    zzName = prefix + "_" + year + "_zz_4l.root"
    zzFile = TFile(zzName, "READ")
    print("Opened", zzName)

    h = 0
    for hname in hnames:
        hist = zzFile.Get("4l/" + hname + "_zz_4l")
        hist.SetDirectory(0)
        axe[h]['4l'].Add(hist)

        h = h + 1
    h = 0
    zzFile.Close()


# Phase space events
ps, amc = np.empty(H, dtype=T), np.empty(H, dtype=T)
h = 0

psName = prefix + "_" + YEAR_STR + "_phase_space.root"
psFile = TFile(psName, "READ")
print("Opened", psName)

amcName = prefix + "_2016_phase_space_aMC.root"
amcFile = TFile(amcName, "READ")
print("Opened", amcName)

for hname in hnames:
    ps[h]['4l'] = psFile.Get("4l/" + hname + "_phase_space")
    ps[h]['4l'].SetDirectory(0)

    amc[h]['4l'] = amcFile.Get("4l/" + hname + "_phase_space_aMC")
    amc[h]['4l'].SetDirectory(0)

    h = h + 1
h = 0
psFile.Close()


for year in ["2017", "2016", "2012"]:
    psName = prefix + "_" + year + "_phase_space.root"
    psFile = TFile(psName, "READ")
    print("Opened", psName)

    h = 0
    for sel in selection:
        for hname in hnames:
            hist = psFile.Get(sel + "/" + hname + "_phase_space")
            hist.SetDirectory(0)
            ps[h][sel].Add(hist)

            h = h + 1
        h = 0
    psFile.Close()

print("Got acc * eff histograms")
print("")



##
##  SCALING
##

scale = alpha * bf_pred * GAMMA_Z
scl = []

for sel in ["4l"]:
    for h in range(H):

        # Divide by acceptance and efficiency
        axe[h][sel].Divide(ps[h][sel])
        for sample in [data, stat]:
            sample[h][sel].Divide(axe[h][sel])

        nbins = data[h][sel].GetNbinsX() + 2
        scl.append(np.ones(nbins))

        for i in range(nbins):
            bc = axe[h][sel].GetBinContent(i)
            if bc != 0:
                scl[-1][i] /= bc
            else:
                scl[-1][i] = 0

        # Scale to partial width
        scl[-1] *= scale / data[h][sel].Integral()

        for sample in [data, stat]:
            sample[h][sel].Scale(scale / sample[h][sel].Integral())
        
        ps[h][sel].Scale(scale / ps[h][sel].Integral())
        amc[h][sel].Scale(scale / amc[h][sel].Integral())



# Systematic uncertainty
for sel in ["4l"]:
    for h in range(H):

        # Add total systematic
        for i in range(data[h][sel].GetNbinsX()):
            data[h][sel].SetBinError(i + 1,
                    np.sqrt((data[h][sel].GetBinContent(i + 1) * delta_syst) ** 2
                        + data[h][sel].GetBinError(i + 1) ** 2))


# Get ratio
ratio, ratio_stat = np.empty(H, dtype=T), np.empty(H, dtype=T)

for sel in ["4l"]:
    for h in range(H):
        ratio_stat[h][sel] = stat[h][sel].Clone()
        ratio_stat[h][sel].Divide(ps[h][sel])

        ratio[h][sel] = data[h][sel].Clone()
        ratio[h][sel].Divide(ps[h][sel])

        for i in range(ratio[h][sel].GetNbinsX()):
            ratio_stat[h][sel].SetBinContent(i, ratio[h][sel].GetBinContent(i))





####
####
####    LOOP OVER DISTS
####
####


for sel in ["4l"]:

    print("Drawing", sel, "plots...")

    for h in range(H):

        ##
        ##  GET BIN CONTENT
        ##

        # Data
        v_data = np.zeros(data[h][sel].GetNbinsX(), dtype=V)
        v_stat = np.zeros(data[h][sel].GetNbinsX(), dtype=V)
        for i in range(len(v_data)):
            v_data[i]['x']      = data[h][sel].GetBinCenter(i+1)
            v_data[i]['y']      = data[h][sel].GetBinContent(i+1)
            v_data[i]['eyu']    = data[h][sel].GetBinErrorUp(i+1)
            v_data[i]['eyd']    = data[h][sel].GetBinErrorLow(i+1)
            v_stat[i]['y']      = stat[h][sel].GetBinContent(i+1)
            v_stat[i]['eyu']    = stat[h][sel].GetBinErrorUp(i+1)
            v_stat[i]['eyd']    = stat[h][sel].GetBinErrorLow(i+1)

#       print(v_stat['y'] - v_data['y'])

        # MC
        v_pred = np.zeros(ps[h][sel].GetNbinsX(), dtype=V)
        for i in range(len(v_pred)):
            v_pred[i]['x']      = ps[h][sel].GetBinLowEdge(i+1)
            v_pred[i]['y']      = ps[h][sel].GetBinContent(i+1)
            v_pred[i]['eyu']    = ps[h][sel].GetBinError(i+1)
            v_pred[i]['eyd']    = ps[h][sel].GetBinError(i+1)

        v_amc = np.zeros(amc[h][sel].GetNbinsX(), dtype=V)
        for i in range(len(v_amc)):
            v_amc[i]['x']       = amc[h][sel].GetBinLowEdge(i+1)
            v_amc[i]['y']       = amc[h][sel].GetBinContent(i+1)
            v_amc[i]['eyu']     = amc[h][sel].GetBinError(i+1)
            v_amc[i]['eyd']     = amc[h][sel].GetBinError(i+1)

        # Ratio
        v_ratio = np.zeros(ratio[h][sel].GetNbinsX(), dtype=V)
        v_ratio_stat = np.zeros(ratio[h][sel].GetNbinsX(), dtype=V)
        for i in range(len(v_ratio)):
            v_ratio[i]['x']         = ratio[h][sel].GetBinCenter(i+1)
            v_ratio[i]['ex']        = ratio[h][sel].GetBinWidth(i+1) / 2
            v_ratio[i]['y']         = ratio[h][sel].GetBinContent(i+1)
            v_ratio[i]['eyu']       = ratio[h][sel].GetBinErrorUp(i+1)
            v_ratio[i]['eyd']       = ratio[h][sel].GetBinErrorLow(i+1)
            v_ratio_stat[i]['eyu']  = ratio_stat[h][sel].GetBinErrorUp(i+1)
            v_ratio_stat[i]['eyd']  = ratio_stat[h][sel].GetBinErrorLow(i+1)



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
                            fmt = 'None',   capsize = lCapSize,     elinewidth = 2 * lErrorLineWidth4l
                            )
        p_amc = ax_top.errorbar(    v_data['x'],    v_amc['y'],     xerr = v_ratio['ex'], 
                            linewidth = 0,  ecolor = lRed,
                            fmt = 'None',   capsize = lCapSize,     elinewidth = 2 * lErrorLineWidth4l
                            )
        ax_top.errorbar(    v_data['x'],    v_data['y'],            yerr = (v_data['eyd'], v_data['eyu']), 
                            linewidth = 0,  ecolor = '#C0C0C0',     elinewidth = 4 * lErrorLineWidth4l,
                            marker = None,   capsize = 0
                            )
        p_data = ax_top.errorbar(   v_data['x'],    v_data['y'],    yerr = (v_stat['eyd'], v_stat['eyu']), 
                            linewidth = 0,  ecolor = lMarkerColor,  elinewidth = lErrorLineWidth4l,
                            marker = 'o',   capsize = 0,            markersize = lMarkerSize2l,
                            markeredgecolor = lMarkerColor,         markerfacecolor = lMarkerColor
                            )

        top_min, top_max = ax_top.get_ylim()

        if hnames[h] in ["b_ttm", "angle_z1l2_z2"]:
            top_max = 3.5
        elif hnames[h] == "b_l1p":
            top_max = 5.5
        elif hnames[h] in ["sin_phi_10"]:
            top_max = 5
        elif hnames[h] == "b_z2m":
            top_max = 10
        elif hnames[h] == "angle_z1leps":
            top_max = 8
        elif hnames[h] in ["angle_z2leps", "cos_theta_z1", "cos_theta_z2", "sin_phi", "b_z1m"]:
            top_max = 4

        ax_top.set_ylim(0, top_max)


        # Ratio plot

        ax_bot.set_ylim(lRatioMin4l, lRatioMax4l)

        ax_bot.axhline(lRatioMid,   color = lBlue,          linewidth=2 * lErrorLineWidth4l)
        ax_bot.errorbar(v_ratio['x'],   v_ratio['y'],       yerr = (v_ratio['eyu'], v_ratio['eyd']), 
                    linewidth = 0,  ecolor = '#C0C0C0',     elinewidth = 4 * lErrorLineWidth4l,
                    marker = None,   capsize = 0
                    )
        ax_bot.errorbar(v_ratio['x'],   v_ratio['y'],       yerr = (v_ratio_stat['eyu'], v_ratio_stat['eyd']), 
                    linewidth = 0,  ecolor = lMarkerColor,  elinewidth = lErrorLineWidth4l,
                    marker = 'o',   capsize = lCapSize,     markersize = lMarkerSize2l,
                    markeredgecolor = lMarkerColor,         markerfacecolor = lMarkerColor
                    )



        ##
        ##  LABELS
        ##

        # Titles
#       ax_top.text(    0.025,  0.95,
#               r'\LARGE{\textbf{CMS}}' + '\n' + r'\Large{\textit{Work in Progress}}',
#               verticalalignment = 'top', transform = ax_top.transAxes)
        ax_top.text(0.025,  0.95,   "CMS",
                size = "xx-large",  weight = "bold",
#               fontproperties = helvet_bold,
                verticalalignment = 'top', transform = ax_top.transAxes, usetex = False)
        ax_top.text(0.025,  0.875,  "Work in Progress",
                size = "x-large",   style = "italic",
#               fontproperties = helvet_bold,
                verticalalignment = 'top', transform = ax_top.transAxes, usetex = False)
        ax_top.set_title(r'\Large{19.7\,fb$^{-1}$ (8\,TeV) $+$ 137\,fb$^{-1}$ (13\,TeV)}', loc='right')

        # Shared x axis
        if hnames[h] in ["b_l1p", "b_ttm", "b_z1m", "b_z2m"]:
            xtitle = "$" + nameTeX[hnames[h]] + "$ (GeV)"
            ytitle = r"$d\Gamma/d" + nameTeX[hnames[h]] + r"$ (keV$/$GeV)"
        elif hnames[h] in ["angle_z1leps", "angle_z2leps", "angle_z1l2_z2"]:
            xtitle = "$" + nameTeX[hnames[h]] + "$ ($\pi$ rad)"
            ytitle = r"$d\Gamma/d" + nameTeX[hnames[h]] + r"$ (keV$/\pi$ rad)"
        elif hnames[h] in ["cos_theta_z1", "cos_theta_z2", "sin_phi", "sin_phi_10"]:
            xtitle = "$" + nameTeX[hnames[h]] + "$"
            ytitle = r"$d\Gamma/d" + nameTeX[hnames[h]] + r"$ (keV$/$unit)"

        ax_bot.set_xlabel(xtitle, horizontalalignment='right')
        ax_bot.xaxis.set_label_coords(1, -0.3)

        # Top y axis
        ax_top.set_ylabel(ytitle, horizontalalignment='right')
        if hnames[h] in ["b_ttm", "b_l1p", "angle_z1leps", "angle_z2leps", "cos_theta_z1",
                "cos_theta_z2", "b_z1m", "b_z2m"]:
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
#       ax_bot.yaxis.set_ticks( np.arange(lRatioMin2l+0.05, lRatioMax4l, step = 0.05),
#                       minor = True  )



        ##
        ##  LEGEND
        ##

        if hnames[h] in ["angle_z1leps", "b_l1p", "b_z1m"]:
            leg_loc = 'center left'
#       elif hnames[h] in ["cos_theta_z1", "cos_theta_z2"]:
#           leg_loc = 'upper center'
        else:
            leg_loc = 'upper right'

        ax_top.legend(
                (   p_data,     p_pred,     p_amc),
                (   'Measured', 'POWHEG',   'aMC@NLO',
                    ),
                loc = leg_loc, numpoints = 1, frameon = False)

        fig.savefig("comb_" + hnames[h] + "_ddr.pdf")
        plt.clf()

print(scl)


# Save arrays
outfile = "scaling.npz"
np.savez(outfile, scl_vec_b_z1m=scl[0], scl_vec_b_z2m=scl[1], scl_vec_b_l1p=scl[2],
        scl_vec_b_ttm=scl[3], scl_vec_cos_theta_z1=scl[4],
        scl_vec_cos_theta_z2=scl[5], scl_vec_angle_z1leps=scl[6],
        scl_vec_angle_z2leps=scl[7], scl_vec_angle_z1l2_z2=scl[8],
        scl_vec_sin_phi=scl[9])

print("Wrote arrays to", outfile)
