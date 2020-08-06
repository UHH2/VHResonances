import time, sys, os, glob
import functools, argparse

import numpy as np

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
        check = not "MC_DY_201" in sample
    if "invisible" in channel and "MC_ZprimeToZH" in sample and not "_inv" in sample:
        check = True
    if not "invisible" in channel and "_inv" in sample:
        check = True
    if not "invisible" in channel and "WJets" in sample:
        check = True
    if "electron"  in channel and "DATA" in sample and not "SingleElectron" in sample: check = True
    if "muon"      in channel and "DATA" in sample and not "SingleMuon" in sample: check = True
    if "invisible" in channel and "DATA" in sample and not "MET" in sample: check = True
    if not "muon"  in channel and "MuonScale" in control_: check = True
    if not "muon"  in channel and "tracking" in control_: check = True
    if not "muon"  in channel and "reco" in control_: check = True
    return check



def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--years',          action='append', dest="years",          default=[])
    parser.add_argument('--Collections',    action='append', dest="Collections",    default=[])
    parser.add_argument('--Channels',       action='append', dest="Channels",       default=[])
    parser.add_argument('--histFolders',    action='append', dest="histFolders",    default=[])
    return parser.parse_args()
