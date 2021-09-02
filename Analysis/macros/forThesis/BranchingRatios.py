from tdrstyle_all import *
from array import array

import tdrstyle_all as TDR
TDR.writeExtraText = True
ForThesis(TDR)

def BranchingRatios():
    rt.gROOT.SetBatch(rt.kTRUE)
    rt.gStyle.SetOptStat(0)
    rt.gStyle.SetOptFit(0)

    BRs = {
        "H": {
            "bb":         58.1,
            "WW":         21.5,
            "gg":         8.2,
            "tautau":     6.3,
            "cc":         2.9,
            "ZZ":         2.6,
            "gammagamma": 0.2,
            "mumu":       0.02,
            "Zgamma":     0.01,
        },
        "Z": {
            "qq":         69.9,
            "nunu":       20.0,
            "ll":         6.74,
            "ee":         3.37,
            "mumu":       3.37,
            "tautau":     3.37,
        },
        "W" : {
            "qq":         67.41,
            "lnu":        21.72,
            "enu":        10.86,
            "munu":       10.86,
            "taunu":      10.86,
        },
        "tau":{
            "lnunu":      35.21,
            "enunu":      17.82,
            "mununu":     17.39,
            "qq":         64.79,
        },
    }

    print "ll+jets"
    print "\t", "ZZ",   "\t", 2*BRs["Z"]["ll"]*BRs["Z"]["qq"]/100.
    print "\t", "ZW",   "\t",   BRs["Z"]["ll"]*BRs["W"]["qq"]/100.
    print "\t", "ZHbb", "\t",   BRs["Z"]["ll"]*BRs["H"]["bb"]/100.
    print "\t", "ZH4q", "\t",   BRs["Z"]["ll"]*BRs["H"]["WW"]*(BRs["W"]["qq"]/100)**2/100.
    print "\t", "ZHcc", "\t",   BRs["Z"]["ll"]*BRs["H"]["cc"]/100.

    print "lnu+jets"
    print "\t", "WW",   "\t", 2*BRs["W"]["lnu"]*BRs["W"]["qq"]/100.
    print "\t", "WZ",   "\t",   BRs["W"]["lnu"]*BRs["Z"]["qq"]/100.
    print "\t", "WHbb", "\t",   BRs["W"]["lnu"]*BRs["H"]["bb"]/100.
    print "\t", "WH4q", "\t",   BRs["W"]["lnu"]*BRs["H"]["WW"]*(BRs["W"]["qq"]/100)**2/100.
    print "\t", "WHcc", "\t",   BRs["W"]["lnu"]*BRs["H"]["cc"]/100.

    print "\t", "WW+tau",   "\t", 2*(BRs["W"]["lnu"]+BRs["W"]["taunu"]*(BRs["tau"]["lnunu"])/100)*BRs["W"]["qq"]/100.
    print "\t", "WZ+tau",   "\t",   (BRs["W"]["lnu"]+BRs["W"]["taunu"]*(BRs["tau"]["lnunu"])/100)*BRs["Z"]["qq"]/100.
    print "\t", "WHbb+tau", "\t",   (BRs["W"]["lnu"]+BRs["W"]["taunu"]*(BRs["tau"]["lnunu"])/100)*BRs["H"]["bb"]/100.
    print "\t", "WH4q+tau", "\t",   (BRs["W"]["lnu"]+BRs["W"]["taunu"]*(BRs["tau"]["lnunu"])/100)*BRs["H"]["WW"]*(BRs["W"]["qq"]/100)**2/100.
    print "\t", "WHcc+tau", "\t",   (BRs["W"]["lnu"]+BRs["W"]["taunu"]*(BRs["tau"]["lnunu"])/100)*BRs["H"]["cc"]/100.

    print "nunu+jets"
    print "\t", "ZZ",   "\t", 2*BRs["Z"]["nunu"]*BRs["W"]["qq"]/100.
    print "\t", "ZW",   "\t",   BRs["Z"]["nunu"]*BRs["W"]["qq"]/100.
    print "\t", "ZHbb", "\t",   BRs["Z"]["nunu"]*BRs["H"]["bb"]/100.
    print "\t", "ZH4q", "\t",   BRs["Z"]["nunu"]*BRs["H"]["WW"]*(BRs["W"]["qq"]/100)**2/100.
    print "\t", "ZHcc", "\t",   BRs["Z"]["nunu"]*BRs["H"]["cc"]/100.

    print "all hadronic"
    print "\t", "ZZ",   "\t",  BRs["Z"]["qq"]*BRs["Z"]["qq"]/100.
    print "\t", "WW",   "\t",  BRs["W"]["qq"]*BRs["W"]["qq"]/100.
    print "\t", "ZW",   "\t",  BRs["Z"]["qq"]*BRs["W"]["qq"]/100.
    print "\t", "ZHbb", "\t",  BRs["Z"]["qq"]*BRs["H"]["bb"]/100.
    print "\t", "ZH4q", "\t",  BRs["Z"]["qq"]*BRs["H"]["WW"]*(BRs["W"]["qq"]/100)**2/100.
    print "\t", "ZHcc", "\t",  BRs["Z"]["qq"]*BRs["H"]["cc"]/100.
    print "\t", "WHbb", "\t",  BRs["W"]["qq"]*BRs["H"]["bb"]/100.
    print "\t", "WH4q", "\t",  BRs["W"]["qq"]*BRs["H"]["WW"]*(BRs["W"]["qq"]/100)**2/100.
    print "\t", "WHcc", "\t",  BRs["W"]["qq"]*BRs["H"]["cc"]/100.



if __name__ == '__main__':
    BranchingRatios()
