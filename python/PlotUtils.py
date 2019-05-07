##
##  COLORS
##

'''
lAlpha = 0.75
lBlue = (0, 0.4470, 0.7410, lAlpha)
lOrange = (0.8500, 0.3250, 0.0980, lAlpha)
lYellow = (0.9290, 0.6940, 0.1250, lAlpha)
lPurple = (0.4940, 0.1840, 0.5560, lAlpha)
lGreen = (0.4660, 0.6740, 0.1880, lAlpha)
lLightBlue = (0.3010, 0.7450, 0.9330, lAlpha)
lRed = (0.6350, 0.0780, 0.1840, lAlpha)
'''
lBlue       = (0.25,    0.58525, 0.80575)
lOrange     = (0.8875,  0.49375, 0.3235)
lYellow     = (0.94675, 0.7705,  0.34375) 
lPurple     = (0.6205,  0.388,   0.667)
lGreen      = (0.5995,  0.7555,  0.391) 
lLightBlue  = (0.47575, 0.80875, 0.94975) 
lRed        = (0.72625, 0.3085,  0.388)


##
##  MPL SETUP
##

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import rc

# Font/TeX setup
# (https://3diagramsperpage.wordpress.com/2015/04/11/)
rc('text', usetex=True)
mpl.rcParams["font.family"] = "serif"
mpl.rcParams["font.serif"] = "Palatino"
mpl.rcParams["font.sans-serif"] = "Helvetica"
#mpl.rcParams['text.latex.preamble'] = [ r'\usepackage{palatino}']
#mpl.rcParams['text.latex.preamble'] = [ r'\usepackage{helvet}',]
#                                       r'\usepackage{sansmath}', r'\sansmath'] 

# Tick marks
mpl.rcParams["xtick.top"], mpl.rcParams["ytick.right"] = True, True
mpl.rcParams["xtick.direction"], mpl.rcParams["ytick.direction"] = "in", "in"
mpl.rcParams["xtick.labelsize"], mpl.rcParams["ytick.labelsize"] = "large", "large"

# Figure size, aspect ratio
mpl.rcParams["figure.figsize"] = [6, 6]
mpl.rcParams["axes.labelsize"], mpl.rcParams["axes.titlesize"] = "xx-large", "xx-large"

# Legend
mpl.rcParams["legend.fontsize"] = "large"



##
##  CONSTANTS
##

lLeftMargin, lRightMargin, lBottomMargin, lTopMargin    =   0.13,   0.93,   0.11,   0.91
lHorizSpace = 0.05

lRatioGridSpec = {'height_ratios':[3, 1]}

lMarkerSize2l, lMarkerSize4l = 5, 10
lMarkerColor = 'black'
lErrorLineWidth2l, lErrorLineWidth4l = 2, 1
lCapSize = 0

lRatioLineColor = 'black'
lRatioMin4l,    lRatioMax4l                     =   0,      2
lRatioMin2l,    lRatioMax2l                     =   0.8,    1.2
lRatioUpper,    lRatioMid,      lRatioLower     =   1.5,    1,      0.5
