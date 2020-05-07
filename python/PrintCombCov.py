from __future__ import print_function
from __future__ import division

import numpy as np

from ROOT import TFile, TH2, TH2D

from PlotUtils import *
import itertools

#from Cuts2018 import *


np.set_printoptions(precision=6, suppress=True, linewidth=200)

##
##  SAMPLE INFO
##

selTeX  = ["4\\mu", "2\\mu2\\mbox{e}", "4\\mbox{e}"]
period  = ["2012", "2016", "2017", "2018"]

period  = list(itertools.chain.from_iterable(itertools.repeat(x, 3) for x in period))
selTeX  = selTeX * 4



##
##  INPUT FILE
##

infile = "combination_1.npz"
npzfile = np.load(infile)

cov_total, cov_stat, cov_syst = npzfile["cov_total"], npzfile["cov_stat"], npzfile["cov_syst"]



##
##  DRAW
##

nbins = 12

# Total covariance matrix
h_cov_total = TH2D("cov_total", "", nbins, 0, nbins, nbins, 0, nbins)
for i in range(nbins):
    h_cov_total.GetXaxis().SetBinLabel(i + 1, selTeX[i] + "\\mbox{ }" + period[i])
    h_cov_total.GetYaxis().SetBinLabel(i + 1, selTeX[i] + "\\mbox{ }" + period[i])
    for j in range(nbins):
        h_cov_total.SetBinContent(i + 1, j + 1, cov_total[i][j])


# Systematic covariance matrix
h_cov_syst = TH2D("cov_syst", "", nbins, 0, nbins, nbins, 0, nbins)
for i in range(nbins):
    h_cov_syst.GetXaxis().SetBinLabel(i + 1, selTeX[i] + "\\mbox{ }" + period[i])
    h_cov_syst.GetYaxis().SetBinLabel(i + 1, selTeX[i] + "\\mbox{ }" + period[i])
    for j in range(nbins):
        h_cov_syst.SetBinContent(i + 1, j + 1, cov_syst[i][j])


# Statistical covariance matrix
h_cov_stat = TH2D("cov_stat", "", nbins, 0, nbins, nbins, 0, nbins)
for i in range(nbins):
    h_cov_stat.GetXaxis().SetBinLabel(i + 1, selTeX[i] + "\\mbox{ }" + period[i])
    h_cov_stat.GetYaxis().SetBinLabel(i + 1, selTeX[i] + "\\mbox{ }" + period[i])
    for j in range(nbins):
        h_cov_stat.SetBinContent(i + 1, j + 1, cov_stat[i][j])



##
##  OUTPUT FILE
##

# Create file
outName = "fit_covariance.root"
outFile = TFile(outName, "RECREATE")

h_cov_total.Write()
h_cov_syst.Write()
h_cov_stat.Write()

outFile.Close()
print("Wrote output to", outName)



##
##  PRINT TEXT
##

filePref = "fit_cov"

name = ["total", "syst", "stat"]
mat = [cov_total, cov_syst, cov_stat]

for i in range(3):
    fileName = filePref + "_" + name[i] + ".tex"
    np.savetxt(fileName, mat[i], fmt='% .2e', comments='', header=r'\begin{verbatim}',
            footer=r'\end{verbatim}')

print("Wrote arrays to", filePref + "_*.tex")
