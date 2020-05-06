#!/usr/bin/env python
from ModuleRunner import *
import time

"""
This macro and guide describes and steers the complete progress of the VHResonances analysis.
Modify these items, everything else will work. The files given here must already be adapted to new settings in case there are any.
"""

# time.sleep(5*60*60)

Collections = ["Puppi"]
Channels = ["muon", "electron"]
# Channels = ["muon"]
# Channels = ["electron"]
# Systematics = ["nominal", "JER_up", "JER_down", "JEC_up", "JEC_down"]
# Systematics = ["JER_up","JER_down", "JEC_up", "JEC_down"]
Systematics = ["nominal"]
# Systematics = ["JER_up"]
# Systematics = ["JER_down"]
# Systematics = ["JEC_up"]
# Systematics = ["JEC_down"]

controls = [""] # default
# controls = ["2018"]
# controls = ["MC_ZprimeToZH"]
# controls = ["MC_ZprimeToZH", "MC_TTbar", "incl", "MC_WW" ] # To Run Local
# controls = ["MC_DY", "MC_WZ", "MC_ZZ", "DATA"] # To Submit

#################################################
#                                               #
#                  Selections                   #
#                                               #
#################################################

for year in ["2016","2017","2018"]:
# for year in ["2016","2017"]:
# for year in ["2018"]:

    isNice=True
    # isNice=False
    nProcess=20

    # Modules = ModuleRunner(year,controls)

    # Modules.SetModule("GenericCleaning", Collections, Channels, Systematics)
    # Modules.SetModule("Preselection", Collections, Channels, Systematics)
    # Modules.SetModule("Selection", Collections, Channels, Systematics)
    # Modules.SetModule("SignalRegion", Collections, Channels, Systematics)
    # Modules.SetModule("ProbeNN", Collections, Channels, Systematics)
    # Modules.SetModule("Test", Collections, Channels, Systematics)
    # Modules.SetModule("SF", Collections, Channels, Systematics)

    # Modules.DeleteWorkdirs()
    # Modules.CreateConfigFiles()
    # Modules.CondorControl("")
    # Modules.RunLocal("Check")
    # Modules.CondorControl("Submit")
    # Modules.CondorControl("")
    # Modules.RunLocal("Check")
    # Modules.RunLocal("Local")
    # Modules.RunLocal("Local", skip=False, nProcess=nProcess, isNice=isNice)
    # time.sleep(9*60*60)
    # Modules.CondorControl("Resubmit")
    # Modules.RunLocal("Check")
    # Modules.CondorControl("List")

    # Modules.StoreModuleOutput()
    # Modules.CreateXml()

    # Modules.SecureMerge(mergeCategory=False)
    # Modules.SecureMerge(mergeCategory=True) # For Preselection only
    # Modules.MakePlots()

    # Modules.HowToSpeedCondor()

#################################################
#                                               #
#                   Analysis                    #
#                                               #
#################################################

controls = [""]
controls = ["SignalRegion"]

Modules = ModuleRunner(controls=controls)

# Modules.MakeRunII(Collections, Channels, Systematics)

# Modules.RunCommand("TaggerCutStudy", isPython=True)
# Modules.RunCommand("CalculateSignalEfficiencies")
# Modules.RunCommand("CreateWorkspace", histFolders=["btag_DeepBoosted_H4qvsQCD"])
# Modules.RunCommand("CreateWorkspace", histFolders=["btag_DeepBoosted_H4qvsQCD", "btag_DeepBoosted_H4qvsQCDp2", "btag_DeepBoosted_H4qvsQCDp02"])
# Modules.RunCommand("CreateDataCards", isPython=True)
# Modules.RunCommand("CompareCombineInputs", isPython=True)
# Modules.RunCommand("PlotLimits")
# Modules.RunCommand("CreateWorkspace", histFolders=["btag_DeepBoosted_H4qvsQCD", "btag_DeepBoosted_H4qvsQCDp2", "btag_DeepBoosted_H4qvsQCDp02", "btag_DeepBoosted_H4qvsQCDpt1000", "btag_DeepBoosted_H4qvsQCDpt1000p2", "btag_DeepBoosted_H4qvsQCDpt1000p02"])
