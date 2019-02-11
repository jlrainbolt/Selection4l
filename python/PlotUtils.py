##
##  COLORS
##

lAlpha = 0.75
lBlue = (0, 0.4470, 0.7410, lAlpha)
lOrange = (0.8500, 0.3250, 0.0980, lAlpha)
lYellow = (0.9290, 0.6940, 0.1250, lAlpha)
lPurple = (0.4940, 0.1840, 0.5560, lAlpha)
lGreen = (0.4660, 0.6740, 0.1880, lAlpha)
lLightBlue = (0.3010, 0.7450, 0.9330, lAlpha)
lRed = (0.6350, 0.0780, 0.1840, lAlpha)



##
##  MPL SETUP
##

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import rc

# Font/TeX setup
# (https://3diagramsperpage.wordpress.com/2015/04/11/)
rc('text', usetex=True)
mpl.rcParams['text.latex.preamble'] = [ r'\usepackage{helvet}', r'\usepackage{sansmath}', 
                                        r'\sansmath'] 
# Figure size, aspect ratio
mpl.rcParams["figure.figsize"] = [6, 6]
mpl.rcParams["axes.labelsize"], mpl.rcParams["axes.titlesize"] = "xx-large", "xx-large"
mpl.rcParams["xtick.labelsize"], mpl.rcParams["ytick.labelsize"] = "large", "large"



##
##  CONSTANTS
##

lLeftMargin, lRightMargin, lBottomMargin, lTopMargin    =   0.13,   0.93,   0.11,   0.91
lHorizSpace = 0.05

lRatioGridSpec = {'height_ratios':[3, 1]}

lMarkerSize = 10
lMarkerColor = 'black'
lErrorLineWidth = 2
lCapSize = 0

lRatioLineColor = 'black'
lRatioMin,      lRatioMax                       =   0.5,    1.5
lRatioUpper,    lRatioMid,      lRatioLower     =   1.2,    1,      0.8
