import time, sys, os, functools, ROOT

sys.path.append(os.environ["CMSSW_BASE"]+"/src/UHH2/VHResonances/Analysis/macros/")
from ModuleRunnerBase import *

ROOT.gInterpreter.ProcessLine('#include "'+os.environ["CMSSW_BASE"]+'/src/UHH2/VHResonances/include/constants.hpp"')
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.gStyle.SetOptStat(0)

from tdrstyle_all import *

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
