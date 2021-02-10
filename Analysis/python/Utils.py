import time, sys, os, glob
import functools, argparse
from sklearn.externals import joblib

import math
import numpy as np
import pandas as pd
from array import array

sys.path.append(os.environ["CMSSW_BASE"]+"/src/UHH2/VHResonances/Analysis/macros/")
from ModuleRunnerBase import *

from tdrstyle_all import *

import ROOT
ROOT.gInterpreter.ProcessLine('#include "'+os.environ["CMSSW_BASE"]+'/src/UHH2/VHResonances/include/constants.hpp"')
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)

def prettydic(d, indent=8):
    space = max([0]+[len(str(x)) for x in d])+2
    for key, value in d.items():
        print(" "*indent + str(key)),
        if isinstance(value, dict):
            print ""
            prettydic(value, len(" "*indent + str(key)+" "*(space+1-len(str(key)))))
        else:
            print(" "*(space-len(str(key))) + str(value))

def timeit(method):
    @functools.wraps(method)
    def timed(*args, **kw):
        print "Start", method.__name__
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()
        if 'log_time' in kw:
            name = kw.get('log_name', method.__name__.upper())
            kw['log_time'][name] = int((te - ts))
        else:
            print '%r  %2.2f s' % \
                  (method.__name__, (te - ts))
        return result
    return timed


def DoControl(controls, control_, channel, sample):
    ''' This part is a bit arbitrary. The idea is to catch all the combinations for different years and channels'''
    check = False
    if all(not control in control_ for control in controls):
        check = True
    if "invisible" in channel and "MC_DY" in sample and not "MC_DY_inv" in sample:
        check = not ("MC_DY_201" in sample or "MC_DY_RunII" in sample)
    if "invisible" in channel and "MC_ZprimeToZH" in sample and not "_inv" in sample:
        check = True
    if not "invisible" in channel and "_inv" in sample:
        check = True
    if not ("invisible" in channel or "charm" in channel) and "WJets" in sample:
        check = True
    if not ("invisible" in channel or "charm" in channel) and "QCD" in sample:
        check = True
    if "electron"  in channel and "DATA" in sample and (not "SingleElectron" in sample and not "SinglePhoton" in sample): check = True
    if "muon"      in channel and "DATA" in sample and not "SingleMuon" in sample: check = True
    if "invisible" in channel and "DATA" in sample and not "MET" in sample: check = True
    if "charm"     in channel and "DATA" in sample and not "JetHT" in sample: check = True
    if not "muon"  in channel and "MuonScale" in control_: check = True
    if not "muon"  in channel and "isolation" in control_: check = True
    if not "muon"  in channel and "tracking"  in control_: check = True
    if not "muon"  in channel and "reco"      in control_: check = True
    return check



def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--years',          action='append', dest="years",          default=[])
    parser.add_argument('--Collections',    action='append', dest="Collections",    default=[])
    parser.add_argument('--Channels',       action='append', dest="Channels",       default=[])
    parser.add_argument('--histFolders',    action='append', dest="histFolders",    default=[])
    return parser.parse_args()

def PrintFormattedLine(listArgs=[], space=10):
    for x in listArgs: print x, " "*(space-len(str(x)) if space-len(str(x))>0 else 2*space-len(str(x))),
    print "\t"



def deltaPhi(p1, p2):
    deltaphi = math.fabs(p1.phi() - p2.phi())
    if deltaphi > math.pi:
        deltaphi = 2* math.pi - deltaphi
    return deltaphi




def deltaR(p1, p2):
    deltaeta = p1.eta() - p2.eta()
    dphi = deltaPhi(p1, p2)
    return math.sqrt(deltaeta * deltaeta + dphi * dphi)




# inline double deltaPhi(const T & p1, const U & p2){
#     double deltaphi = fabs(p1.phi() - p2.phi());
#     if(deltaphi > M_PI) deltaphi = 2* M_PI - deltaphi;
#     return deltaphi;
# }
#
# /// distance in eta-phi space. works for any types which have 'eta' and 'phi' routines
# // T and U have to have a 'phi()' and an 'eta()' method, e.g. Particle, LorentzVector, etc.
# template<typename T, typename U>
# inline double deltaR(const T & p1, const U & p2){
#     double deltaeta = p1.eta() - p2.eta();
#     double dphi = deltaPhi(p1, p2);
#     return sqrt(deltaeta * deltaeta + dphi * dphi);
# }
#
# }
