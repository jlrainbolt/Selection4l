from __future__ import print_function
from __future__ import division
import sys

import numpy as np

from ROOT import TFile, TH1

from PlotUtils import *
from Cuts2018 import *
#from Cuts2017 import *
#from Cuts2016 import *
#from Cuts2012 import *



##
##  SAMPLE INFO
##

selection = ["4l"]

T = np.dtype([(sel, object) for sel in selection])
V = np.dtype([("x", 'f4'), ("y", 'f4'), ("ex", 'f4'), ("ey", 'f4'), ("b", 'f4')])

#tightOnly = True
tightOnly = False



##
##  DATA
##

if tightOnly:
    prefix = "bkg_tight"
else:
    prefix = "bkg_loose"

# Muon file
muName = prefix + "_" + YEAR_STR + "_muon_" + YEAR_STR + ".root"
muFile = TFile(muName, "READ")
print("Opened", muName)

# Electron file
elName = prefix + "_" + YEAR_STR + "_electron_" + YEAR_STR + ".root"
elFile = TFile(elName, "READ")
print("Opened", elName)


# Get histograms
hnames = ["zzm", "zzpt", "z1m", "z2m", "z1pt", "z2pt"]

H = len(hnames)

data = np.empty(H, dtype=T)
h = 0

for sel in selection:
    for hname in hnames:
        data[h][sel] = muFile.Get(sel + "/" + hname + "_muon_" + YEAR_STR)
        data[h][sel].Add(elFile.Get(sel + "/" + hname + "_electron_" + YEAR_STR))

        data[h][sel].SetDirectory(0)
        h = h + 1
    h = 0

muFile.Close()
elFile.Close()



##
##  MONTE CARLO
##

mc_arr = np.empty((N_MC, H), dtype=T)
mc = {}
h, j = 0, 0

# Loop over all samples
for suff in MC_SUFF:
    inName = prefix + "_" + YEAR_STR + "_" + suff + ".root"
    inFile = TFile.Open(inName)
    print("Opened", inName)

    # Get histograms
    for sel in selection:
        sf = INT_LUMI * 1000 * XSEC[suff] / NGEN[suff]

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
##  REBIN
##
'''
if not tightOnly:
    for h in range(H):
        if (YEAR_STR == "2012"):
            for sel in selection:
                data[h][sel].Rebin(2)

        for suff in MC_SUFF:
            if (YEAR_STR == "2012"):
                for sel in selection:
                    mc[suff][h][sel].Rebin(2)
'''

# Get total
total, ratio = np.empty(H, dtype=T), np.empty(H, dtype=T)

for sel in selection:
    for h in range(H):
        for suff in MC_SUFF:
            if suff == "zz_4l":
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


for sel in ["4l"]:
#for sel in selection:


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
        for suff in MC_SUFF:
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
        for suff in MC_SUFF:
            p_mc[suff] = ax_top.bar(    v_mc[suff]['x'],    v_mc[suff]['y'],    width,
                                bottom = v_mc[suff]['b'],   align = 'edge',     linewidth=0,
                                color = COLOR[suff]
                                )

        top_min, top_max = ax_top.get_ylim()

        if hnames[h] == "zzm":
            top_max = 1.2 * top_max
        else:
            top_max = 1.3 * top_max

        ax_top.set_ylim(0, top_max)


        # Ratio plot
        ax_bot.errorbar(v_ratio['x'],   v_ratio['y'],   xerr = v_ratio['ex'],   yerr = v_ratio['ey'], 
                    linewidth = 0,  ecolor = lMarkerColor,  elinewidth = lErrorLineWidth4l,
                    marker = 'o',   capsize = lCapSize,     markersize = lMarkerSize4l,
                    markeredgecolor = lMarkerColor,         markerfacecolor = lMarkerColor
                    )

        ax_bot.axhline(lRatioMid,   color = lRatioLineColor, linestyle = ':')
#       if sel == "4e":
#           ax_bot.set_ylim(0, 3)
#       else:
        ax_bot.set_ylim(lRatioMinBkg, lRatioMaxBkg)



        ##
        ##  LABELS
        ##

        # Titles
        ax_top.text(0.025,  0.95,   "CMS",
                size = "xx-large",  weight = "bold",    family = "Liberation Sans",
                verticalalignment = 'top', transform = ax_top.transAxes, usetex = False)
        ax_top.text(0.025,  0.875,  "Work in Progress",
                size = "x-large",   style = "italic",    family = "Liberation Sans",
                verticalalignment = 'top', transform = ax_top.transAxes, usetex = False)
        ax_top.set_title(r'\Large{' + '%.1f' % INT_LUMI + r'\,fb$^{-1}$ (' + '%i' % SQRT_S
                + r'\,TeV, ' + YEAR_STR + ')}', loc='right')

        # Top y axis
        unit = '$' + mc['zz_4l'][h][sel].GetYaxis().GetTitle() + '$'
        ytitle = r"Events$/$" + '%g' % width + " " + unit
        if hnames[h] == "zzm":
            if sel == "4e":
                ytitle = r'Events$/$2 GeV'
            else:
                ytitle = r'Events$/$GeV'
        if hnames[h] == "sin_phi":
            ytitle = r'Events$/$0.1 units'
        ax_top.set_ylabel(ytitle, horizontalalignment='right')
        ax_top.yaxis.set_label_coords(-0.08, 1)
        ax_top.minorticks_on()

        # Bottom y axis
        ax_bot.set_ylabel(r'Data$/$MC')
        ax_bot.yaxis.set_label_coords(-0.08, 0.5)

        # Shared x axis
        xtitle = '$' + mc['zz_4l'][h][sel].GetXaxis().GetTitle() + '$'
        if hnames[h] == "zzm":
            if sel == "4l":
                xtitle = r'$m_{4\ell}$ (GeV)'
            elif sel == "4m":
                xtitle = r'$m_{4\mu}$ (GeV)'
            elif sel == "2m2e":
                xtitle = r'$m_{2\mu 2\mathrm{e}}$ (GeV)'
            elif sel == "4e":
                xtitle = r'$m_{4\mathrm{e}}$ (GeV)'
#       if "Delta" in xtitle:
#           xtitle = xtitle.replace("Delta", "bigtriangleup")
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
        ax_bot.yaxis.set_ticks( np.arange(lRatioMinBkg+1, lRatioMaxBkg, step = 1) )



        ##
        ##  LEGEND
        ##

        if hnames[h] in ["zzm", "z1m", "angle_z1leps", "b_l1p"]:#, "cos_theta_z2"]:
            leg_loc = 'center left'
        elif hnames[h] in ["sin_phi_2"]:
            leg_loc = 'upper center'
        else:
            leg_loc = 'upper right'

        if hnames[h] in ["sin_phi", "cos_theta_z1", "cos_theta_z2"]:
            leg_ncol = 2
        else:
            leg_ncol = 1

        handles = [ p_data,     p_mc['zjets_m-50'],                 p_mc['zz_4l'],
                                p_mc['ww_2l2nu'],
                    ]
        labels = [  r'Data',    r'$\mbox{Z}\to\ell^+\ell^-$',       r'$\mbox{ZZ}\to4\ell$',
                                r'VV',
                    ]

        if YEAR_STR != "2012":
            handles.append(p_mc['zzz_4l2nu'])
            labels.append(r'VVV')
            handles.append(p_mc['ggH_zz_4l'])
            labels.append(r'H')
        handles.append(p_mc['ttz_2l2nu'])
        labels.append(r'$\mbox{t}\bar{\mbox{t}}\mbox{Z}$')

        ax_top.legend(handles, labels, loc = leg_loc, numpoints = 1, frameon = False, ncol = leg_ncol)

        fig_name = YEAR_STR + "_" + hnames[h] + "_" + prefix + ".pdf"
        fig.savefig(fig_name)
        plt.clf()
