from __future__ import print_function
from __future__ import division
import sys

import numpy as np

from ROOT import TFile, TH1

from PlotUtils import *



##
##  SAMPLE INFO
##

selection = ["4l", "4m", "2m2e", "4e"]
period = ["2016", "2017", "2018"]

T = np.dtype([(sel, object) for sel in selection])
V = np.dtype([("x", 'f4'), ("y", 'f4'), ("ex", 'f4'), ("ey", 'f4'), ("b", 'f4')])


# Get header info
lumi, xsec, ngen, mc_suff = {}, {}, {}, {}

from Cuts2016 import *
lumi[YEAR_STR] = INT_LUMI
xsec[YEAR_STR] = XSEC
ngen[YEAR_STR] = NGEN
mc_suff[YEAR_STR] = MC_SUFF_AMC

from Cuts2017 import *
lumi[YEAR_STR] = INT_LUMI
xsec[YEAR_STR] = XSEC
ngen[YEAR_STR] = NGEN
mc_suff[YEAR_STR] = MC_SUFF_AMC

from Cuts2018 import *
lumi[YEAR_STR] = INT_LUMI
xsec[YEAR_STR] = XSEC
ngen[YEAR_STR] = NGEN
mc_suff[YEAR_STR] = MC_SUFF_AMC



##
##  DATA
##

prefix = "4l"
year = "2018"

# Muon file
muName = prefix + "_" + year + "_muon_" + year + ".root"
muFile = TFile(muName, "READ")
print("Opened", muName)

# Electron file
elName = prefix + "_" + year + "_electron_" + year + ".root"
elFile = TFile(elName, "READ")
print("Opened", elName)


# Get histograms for 2018
#hnames = ["zzm"]
#hnames = ["sin_phi"]
hnames = ["zzm", "zzpt", "z1m", "z2m", "z1pt", "z2pt", "l1pt", "l2pt", "l3pt", "l4pt", "l1eta", "l2eta", "l3eta", "l4eta"]
#hnames = ["b_ttm", "b_l1p", "cos_theta_z1", "cos_theta_z2", "angle_z1leps", "angle_z2leps", "angle_z1l2_z2"]
#hnames = ["b_z1m"]

H = len(hnames)

data = np.empty(H, dtype=T)
h = 0

for sel in selection:
    for hname in hnames:
        data[h][sel] = muFile.Get(sel + "/" + hname + "_muon_" + year)
        data[h][sel].Add(elFile.Get(sel + "/" + hname + "_electron_" + year))

        data[h][sel].SetDirectory(0)
#       data[h][sel].SetBinErrorOpt(kPoisson)

        h = h + 1
    h = 0

muFile.Close()
elFile.Close()


# Add 2016 & 2017
for year in ["2017", "2016"]:
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
        for hname in hnames:
            hist = muFile.Get(sel + "/" + hname + "_muon_" + year)
            hist.SetDirectory(0)
            data[h][sel].Add(hist)

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
year = "2018"

# Loop over all samples
for suff in mc_suff[year]:
    if suff == "zjets_m-50":
        mc[suff] = mc_arr[j]
        j = j + 1
        continue

    # Get 2018 histograms
    year = "2018"
    inName = prefix + "_" + year + "_" + suff + ".root"
    inFile = TFile.Open(inName)
    print("Opened", inName)

    sf = lumi[year] * 1000 * xsec[year][suff] / ngen[year][suff]

    for sel in selection:
        for hname in hnames:
            mc_arr[j][h][sel] = inFile.Get(sel + "/" + hname + "_" + suff)
            mc_arr[j][h][sel].SetDirectory(0)
            mc_arr[j][h][sel].Scale(sf)

            h = h + 1
        h = 0

    # Add 2016 & 2017
    for year in ["2017", "2016"]:
        if suff not in mc_suff[year]:
            continue

        inName = prefix + "_" + year + "_" + suff + ".root"
        inFile = TFile.Open(inName)
        print("Opened", inName)

        sf = lumi[year] * 1000 * xsec[year][suff] / ngen[year][suff]

        for sel in selection:
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

# Get number of nonprompt events
npt = {}

for year in period:
    infile = "nonprompt" + year + ".npz"
    npzfile = np.load(infile)
    npt[year] = npzfile['npt']


prefix = "bkg_all"

# Muon file
muName = prefix + "_" + YEAR_STR + "_muon_" + YEAR_STR + ".root"
muFile = TFile(muName, "READ")
print("Opened", muName)

# Electron file
elName = prefix + "_" + YEAR_STR + "_electron_" + YEAR_STR + ".root"
elFile = TFile(elName, "READ")
print("Opened", elName)


# Get histograms for 2018
h = 0
year = "2018"
for sel in selection:
    for hname in hnames:
        mc['zjets_m-50'][h][sel] = muFile.Get(sel + "/" + hname + "_muon_" + YEAR_STR)
        mc['zjets_m-50'][h][sel].Add(elFile.Get(sel + "/" + hname + "_electron_" + YEAR_STR))

        mc['zjets_m-50'][h][sel].SetDirectory(0)
        sf = npt[year][sel] / mc['zjets_m-50'][h][sel].Integral()
        mc['zjets_m-50'][h][sel].Scale(sf)

        h = h + 1
    h = 0

muFile.Close()
elFile.Close()


# Add 2016 & 2017
for year in ["2017", "2016"]:
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
        for hname in hnames:
            hist = muFile.Get(sel + "/" + hname + "_muon_" + year)
            hist.Add(elFile.Get(sel + "/" + hname + "_electron_" + year))

            hist.SetDirectory(0)
            sf = npt[year][sel] / hist.Integral()
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

year = "2018"

# Rebin 4e
for h in range(H):
    data[h]['4e'].Rebin(2)
    for suff in mc_suff[year]:
        mc[suff][h]['4e'].Rebin(2)


# Get total
total, ratio = np.empty(H, dtype=T), np.empty(H, dtype=T)

for sel in selection:
    for h in range(H):
        for suff in mc_suff[year]:
            if suff == "zz_4l_aMC":
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


#for sel in ["4e"]:
for sel in selection:
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
        for suff in MC_SUFF_AMC:
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
        for suff in MC_SUFF_AMC:
            p_mc[suff] = ax_top.bar(    v_mc[suff]['x'],    v_mc[suff]['y'],    width,
                                bottom = v_mc[suff]['b'],   align = 'edge',     linewidth=0,
                                color = COLOR[suff]
                                )

        top_min, top_max = ax_top.get_ylim()

        if hnames[h] == "zzm":
            if sel == "4l":
                top_max = 350
            elif sel == "4m":
                top_max = 200
            elif sel == "2m2e":
                top_max = 140
            elif sel == "4e":
                top_max = 70
        elif hnames[h] == "sin_phi":
            if sel == "4l":
                top_max = 550
            elif sel == "4m":
                top_max = 275
            elif sel == "2m2e":
                top_max = 225
            elif sel == "4e":
                top_max = 85
        elif hnames[h] in ["b_ttm", "cos_theta_z2"]:
            top_max = 1.5 * top_max
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
        ax_bot.set_ylim(lRatioMin4l, lRatioMax4l)



        ##
        ##  LABELS
        ##

        # Titles
        ax_top.text(0.025,  0.95,   "CMS",
                size = "xx-large",  weight = "bold",    family = "Liberation Sans",
                verticalalignment = 'top', transform = ax_top.transAxes, usetex = False)
        ax_top.text(0.025,  0.875,  "Work in Progress",
                size = "x-large",   style = "italic",   family = "Liberation Sans",
                verticalalignment = 'top', transform = ax_top.transAxes, usetex = False)
        ax_top.set_title(r'\Large{137\,fb$^{-1}$ (13\,TeV)', loc='right')

        # Top y axis
        unit = '$' + mc['zz_4l_aMC'][h][sel].GetYaxis().GetTitle() + '$'
        ytitle = r"Events$/$" + '%g' % width + " " + unit
        if hnames[h] == "zzm":
            if sel == "4e":
                ytitle = r'Events$/$2 GeV'
            else:
                ytitle = r'Events$/$GeV'
        if hnames[h] == "sin_phi":
            if sel == "4e":
                ytitle = r'Events$/$0.2 units'
            else:
                ytitle = r'Events$/$0.1 units'
        ax_top.set_ylabel(ytitle, horizontalalignment='right')
        ax_top.yaxis.set_label_coords(-0.08, 1)
        ax_top.minorticks_on()

        # Bottom y axis
        ax_bot.set_ylabel(r'Data$/$MC')
        ax_bot.yaxis.set_label_coords(-0.08, 0.5)

        # Shared x axis
        xtitle = '$' + mc['zz_4l_aMC'][h][sel].GetXaxis().GetTitle() + '$'
        if hnames[h] == "zzm":
            if sel == "4l":
                xtitle = r'$m_{4\ell}$ (GeV)'
            elif sel == "4m":
                xtitle = r'$m_{4\mu}$ (GeV)'
            elif sel == "2m2e":
                xtitle = r'$m_{2\mu 2\mathrm{e}}$ (GeV)'
            elif sel == "4e":
                xtitle = r'$m_{4\mathrm{e}}$ (GeV)'
        xtitle = xtitle.replace("mbox", "mathrm")
        ax_bot.set_xlabel(xtitle, horizontalalignment='right')
        ax_bot.xaxis.set_label_coords(1, -0.3)

        

        ##
        ##  TICKS
        ##

        # x axes
        plt.xlim(v_mc['zz_4l_aMC']['x'][0], v_mc['zz_4l_aMC']['x'][-1] + width)

        major_step, minor_step = 2 * width, width
        if sel == "4e" and hnames[h] == "zzm":
            major_step = width

        for ax in [ax_bot.xaxis, ax_top.xaxis]:
            ax.set_ticks( np.arange(
                            v_mc['zz_4l_aMC']['x'][0],
                            v_mc['zz_4l_aMC']['x'][-1] + major_step,
                            step = major_step)
                            )
            ax.set_ticks( np.arange(
                            v_mc['zz_4l_aMC']['x'][0],
                            v_mc['zz_4l_aMC']['x'][-1] + minor_step,
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

        if hnames[h] in ["zzm", "z1m", "angle_z1leps", "b_l1p", "b_z1m"]:
            leg_loc = 'center left'
        else:
            leg_loc = 'upper right'

        if hnames[h] in ["sin_phi", "cos_theta_z1", "cos_theta_z2", "b_ttm", "l1eta", "l2eta", "l3eta", "l4eta"]:
            leg_ncol = 2
        else:
            leg_ncol = 1

        ax_top.legend(
                (   p_data,
                    p_mc['zz_4l_aMC'],      p_mc['zjets_m-50'], p_mc['ww_2l2nu'],
                    p_mc['zzz_4l2nu'],      p_mc['ggH_zz_4l'],  p_mc['ttz_2l2nu']
                    ),
                (   r'Data',
                    r'$\mbox{Z}\to4\ell$ (a\textsc{mc@nlo)',    r'Nonprompt',       r'VV',
                    r'VVV',                 r'H',               r'$\mbox{t}\bar{\mbox{t}}\mbox{Z}$' 
                    ),
                loc = leg_loc, numpoints = 1, frameon = False, ncol = leg_ncol
            )

        fig_name = "comb_aMC_" + hnames[h] + "_" + sel + ".pdf"
        fig.savefig(fig_name)
        plt.clf()
