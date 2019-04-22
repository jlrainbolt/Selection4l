from __future__ import print_function
from __future__ import division
import sys

import numpy as np

from ROOT import TFile, TH1

from PlotUtils import *
from CutsComb import *



##
##  SAMPLE INFO
##

selection = ["4l", "4m", "2m2e", "2e2m", "4e"]

T = np.dtype([(sel, object) for sel in selection])
V = np.dtype([("x", 'f4'), ("y", 'f4'), ("ex", 'f4'), ("ey", 'f4'), ("b", 'f4')])



##
##  DATA
##

prefix = "4l"

# Muon file
muName = prefix + "_" + YEAR_STR + "_muon_" + YEAR_STR + ".root"
muFile = TFile(muName, "READ")
print("Opened", muName)

# Electron file
elName = prefix + "_" + YEAR_STR + "_electron_" + YEAR_STR + ".root"
elFile = TFile(elName, "READ")
print("Opened", elName)


# Get histograms for 2017
#hnames = ["zzm"]
hnames = ["sin_phi"]
#hnames = ["b_ttm", "b_l1p", "cos_theta_z1", "cos_theta_z2",
#            "angle_z1leps", "angle_z2leps", "angle_z1l2_z2"]

H = len(hnames)

data = np.empty(H, dtype=T)
h = 0

for sel in selection:
    if sel == "4l":
        continue

    for hname in hnames:
        if sel in ["4m", "2m2e"]:
            data[h][sel] = muFile.Get(sel + "/" + hname + "_muon_" + YEAR_STR)
        elif sel in ["4e", "2e2m"]:
            data[h][sel] = elFile.Get(sel + "/" + hname + "_electron_" + YEAR_STR)

        data[h][sel].SetDirectory(0)
        h = h + 1
    h = 0

muFile.Close()
elFile.Close()


# Add 2016 & 2012
for year in ["2016", "2012"]:
    # Muon file
    muName = prefix + "_" + year + "_muon_" + year + ".root"
    muFile = TFile(muName, "READ")
    print("Opened", muName)

    # Electron file
    elName = prefix + "_" + year + "_electron_" + year + ".root"
    elFile = TFile(elName, "READ")
    print("Opened", elName)

    # Get histograms
    h = 0

    for sel in selection:
        if sel == "4l":
            continue

        for hname in hnames:
            if sel in ["4m", "2m2e"]:
                hist = muFile.Get(sel + "/" + hname + "_muon_" + year)
            elif sel in ["4e", "2e2m"]:
                hist = elFile.Get(sel + "/" + hname + "_electron_" + year)

            hist.SetDirectory(0)
            data[h][sel].Add(hist)
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
for suff in MC_SUFF:
    if suff == "zjets_m-50":
        mc[suff] = mc_arr[j]
        j = j + 1
        continue

    inName = prefix + "_" + YEAR_STR + "_" + suff + ".root"
    inFile = TFile.Open(inName)
    print("Opened", inName)

    # Get histograms
    for sel in selection:
        if sel == "4l":
            continue
        elif sel in ["4m", "2m2e"]:
            lumi = MUON_TRIG_LUMI_2017
        elif sel in ["4e", "2e2m"]:
            lumi = ELEC_TRIG_LUMI_2017 * ELEC_TRIG_SF_2017

        sf = lumi * 1000 * XSEC[suff] / NGEN[suff]

        for hname in hnames:
            mc_arr[j][h][sel] = inFile.Get(sel + "/" + hname + "_" + suff)
            mc_arr[j][h][sel].SetDirectory(0)
            mc_arr[j][h][sel].Scale(sf)

            h = h + 1
        h = 0

    if suff in MC_SUFF_2012:
        year = "2012"
        inName = prefix + "_" + year + "_" + suff + ".root"
        inFile = TFile.Open(inName)
        print("Opened", inName)
 
        # Get histograms
        for sel in selection:
            if sel == "4l":
                continue
            elif sel in ["4m", "2m2e"]:
                lumi = MUON_TRIG_LUMI_2012
            elif sel in ["4e", "2e2m"]:
                lumi = ELEC_TRIG_LUMI_2012 * ELEC_TRIG_SF_2012
 
            sf = lumi * 1000 * XSEC_2012[suff] / NGEN_2012[suff]
 
            for hname in hnames:
                hist = inFile.Get(sel + "/" + hname + "_" + suff)
                hist.SetDirectory(0)
                hist.Scale(sf)

                mc_arr[j][h][sel].Add(hist)
 
                h = h + 1
            h = 0

    if suff in MC_SUFF_2016:
        year = "2016"
        inName = prefix + "_" + year + "_" + suff + ".root"
        inFile = TFile.Open(inName)
        print("Opened", inName)
 
        # Get histograms
        for sel in selection:
            if sel == "4l":
                continue
            elif sel in ["4m", "2m2e"]:
                lumi = MUON_TRIG_LUMI_2016
            elif sel in ["4e", "2e2m"]:
                lumi = ELEC_TRIG_LUMI_2016 * ELEC_TRIG_SF_2016
 
            sf = lumi * 1000 * XSEC_2016[suff] / NGEN_2016[suff]
 
            for hname in hnames:
                hist = inFile.Get(sel + "/" + hname + "_" + suff)
                hist.SetDirectory(0)
                hist.Scale(sf)

                mc_arr[j][h][sel].Add(hist)
 
                h = h + 1
            h = 0


    mc[suff] = mc_arr[j]
    j = j + 1
    inFile.Close()

print("Got MC histograms")
print("")



##
##  BACKGROUND
##

prefix = "bkg_all"

# Muon file
muName = prefix + "_" + YEAR_STR + "_muon_" + YEAR_STR + ".root"
muFile = TFile(muName, "READ")
print("Opened", muName)

# Electron file
elName = prefix + "_" + YEAR_STR + "_electron_" + YEAR_STR + ".root"
elFile = TFile(elName, "READ")
print("Opened", elName)


# Get histograms for 2017
h = 0
for sel in selection:
    if sel == "4l":
        continue

    for hname in hnames:
        if sel in ["4m", "2m2e"]:
            mc['zjets_m-50'][h][sel] = muFile.Get(sel + "/" + hname + "_muon_" + YEAR_STR)
        elif sel in ["4e", "2e2m"]:
            mc['zjets_m-50'][h][sel] = elFile.Get(sel + "/" + hname + "_electron_" + YEAR_STR)

        mc['zjets_m-50'][h][sel].SetDirectory(0)
        sf = npt_2017[sel] / mc['zjets_m-50'][h][sel].Integral()
        mc['zjets_m-50'][h][sel].Scale(sf)

        h = h + 1
    h = 0

muFile.Close()
elFile.Close()


# Add 2016 & 2012
for year in ["2016", "2012"]:
    # Muon file
    muName = prefix + "_" + year + "_muon_" + year + ".root"
    muFile = TFile(muName, "READ")
    print("Opened", muName)

    # Electron file
    elName = prefix + "_" + year + "_electron_" + year + ".root"
    elFile = TFile(elName, "READ")
    print("Opened", elName)

    # Get histograms
    h = 0

    for sel in selection:
        if sel == "4l":
            continue

        for hname in hnames:
            if sel in ["4m", "2m2e"]:
                hist = muFile.Get(sel + "/" + hname + "_muon_" + year)
            elif sel in ["4e", "2e2m"]:
                hist = elFile.Get(sel + "/" + hname + "_electron_" + year)

            hist.SetDirectory(0)
            if year == "2012":
                sf = npt_2012[sel] / hist.Integral()
            elif year == "2016":
                sf = npt_2016[sel] / hist.Integral()

            hist.Scale(sf)
            mc['zjets_m-50'][h][sel].Add(hist)
            h = h + 1
        h = 0

    muFile.Close()
    elFile.Close()
print("Got data histograms")
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

    data[h]['4e'].Rebin(2)

    for suff in MC_SUFF:
        mc[suff][h]['2m2e'].Add(mc[suff][h]['2e2m'])
        mc[suff][h]['4l'] = mc[suff][h]['2m2e'].Clone()
        mc[suff][h]['4l'].Add(mc[suff][h]['4m'])
        mc[suff][h]['4l'].Add(mc[suff][h]['4e'])

        mc[suff][h]['4e'].Rebin(2)


# Get total
total, ratio = np.empty(H, dtype=T), np.empty(H, dtype=T)

for sel in selection:
    if sel == "2e2m":
        continue

    for h in range(H):
        for suff in MC_SUFF:
            if suff == "zz_4l":
                total[h][sel] = mc[suff][h][sel].Clone()
            else:
                total[h][sel].Add(mc[suff][h][sel])

        ratio[h][sel] = data[h][sel].Clone()
        ratio[h][sel].Divide(total[h][sel])

        print(mc['zjets_m-50'][h][sel].Integral())





####
####
####    LOOP OVER STACKS
####
####


#for sel in ["4e"]:
for sel in selection:
    if sel == "2e2m":
        continue


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
            if sel == "4l":
                top_max = 200
            elif sel == "4m":
                top_max = 120
            elif sel == "2m2e":
                top_max = 70
            elif sel == "4e":
                top_max = 40

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
        ax_bot.set_ylim(lRatioMin4l, lRatioMax4l)



        ##
        ##  LABELS
        ##

        # Titles
        ax_top.text(    0.025,  0.95,
                r'\LARGE{\textbf{CMS}}' + '\n' + r'\Large{\textit{Work in Progress}}',
                verticalalignment = 'top', transform = ax_top.transAxes)
#       ax_top.set_title(r'\Large{19.7\,fb$^{-1}$ (8\,TeV) $+$ 77.8\,fb$^{-1}$ (13\,TeV, 2016--17)', loc='right')
        ax_top.set_title(r'\Large{19.7\,fb$^{-1}$ (8\,TeV) $+$ 77.8\,fb$^{-1}$ (13\,TeV)', loc='right')

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
                            v_mc['zz_4l']['x'][-1] + 2 * width,
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
        ax_bot.yaxis.set_ticks( np.arange(lRatioMin4l+0.5, lRatioMax4l, step = 0.5) )
#       ax_bot.yaxis.set_ticks( np.arange(lRatioMin2l+0.05, lRatioMax4l, step = 0.05),
#                       minor = True  )



        ##
        ##  LEGEND
        ##

        if hnames[h] == "sin_phi":
            leg_loc = 'upper right'
            bbox = (0.86, 1)
        else:
            leg_loc = 'center left'

        ax_top.legend(
                (   p_data,
                    p_mc['zz_4l'],          p_mc['zjets_m-50'], p_mc['ww_2l2nu'],
                    p_mc['zzz_4l2nu'],      p_mc['ggH_zz_4l'],  p_mc['ttz_2l2nu']
                    ),
                (   r'Data',
                    r'$\mbox{Z}\to4\ell$',  r'Nonprompt',       r'VV',
                    r'VVV',                 r'H',               r'$\mbox{t}\bar{\mbox{t}}\mbox{Z}$' 
                    ),
                loc = leg_loc, numpoints = 1, frameon = False#, bbox_to_anchor = bbox
            )

        fig_name = "comb_" + hnames[h] + "_" + sel + ".pdf"
        fig.savefig(fig_name)
        plt.clf()
