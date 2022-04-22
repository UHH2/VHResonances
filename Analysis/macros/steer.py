#!/usr/bin/env python
from ModuleRunner import *
import time

"""
This macro and guide describes and steers the complete progress of the VHResonances analysis.
Modify these items, everything else will work. The files given here must already be adapted to new settings in case there are any.
"""

# time.sleep(3*60*60)

Collections = ["Puppi"]
Channels = ["muon", "electron", "invisible"]
# Channels = ["muon", "invisible"]
# Channels = ["muon", "electron"]
# Channels = ["muon"]
# Channels = ["electron"]
# Channels = ["invisible"]
# Channels = ["charm"]
# Systematics = ["nominal", "JER_up", "JER_down", "JEC_up", "JEC_down", "MuonScale_up", "MuonScale_down"]
# Systematics = ["JER_up", "JER_down", "JEC_up", "JEC_down", "MuonScale_up", "MuonScale_down"]
# Systematics = ["JER_up","JER_down", "JEC_up", "JEC_down"]
# Systematics = ["JER_up","JER_down"]
# Systematics = ["JEC_up", "JEC_down"]
# Systematics = ["MuonScale_up", "MuonScale_down"]
Systematics = ["nominal"]
# Systematics = ["JER_up"]
# Systematics = ["JER_down"]
# Systematics = ["JEC_up"]
# Systematics = ["JEC_down"]

controls = [""] # default
# controls = ["UL18"]
# controls = ["DATA"]
# controls = ["RunD"]
# controls = ["MC"]
# controls = ["MET"]
# controls = ["MC_ZprimeToZH"]
# controls = ["MC_DY"]
# controls = ["MC_DY", "MC_ZprimeToZH"]
# controls = ["MC_TT", "MC_WW", "MC_WZ", "MC_ZZ"]
# controls = ["DATA", "MC_TT", "MC_WW", "MC_WZ", "MC_ZZ"]
# controls = ["MC_WJets"]
# controls = ["MC_QCD"]
# controls = ["MC_TT"]
# controls = ["MC_TT", "MC_WZ", "MC_ZZ", "DATA"]


#################################################
#                                               #
#                  Selections                   #
#                                               #
#################################################
# for year in ["UL16preVFP", "UL16postVFP", "UL17","UL18"]:
# for year in ["UL16postVFP", "UL17","UL18"]:
for year in ["UL18"]:

    isNice=True
    # isNice=False
    nProcess=20
    Modules = ModuleRunner(year,controls)
    # Modules = ModuleRunner(year,controls, isAnalysis=False)

    # Modules.SetModule("GenericCleaning", Collections, Channels, Systematics)
    Modules.SetModule("HEMIssueStudy", Collections, Channels, Systematics)
    # Modules.SetModule("PDFReweight", Collections, Channels, Systematics)
    # Modules.SetModule("Preselection", Collections, Channels, Systematics)
    # Modules.SetModule("Selection", Collections, Channels, Systematics)
    # Modules.SetModule("SignalRegion", Collections, Channels, Systematics)
    # Modules.SetModule("ProbeNN", Collections, Channels, Systematics)
    # Modules.SetModule("Test", Collections, Channels, Systematics)
    # Modules.SetModule("SF", Collections, Channels, Systematics)
    # Modules.SetModule("LeptonIDStudies", Collections, Channels, Systematics)

    # Modules.SetModule("HccSFSelection", Collections, Channels, Systematics)

    Modules.DeleteWorkdirs()
    Modules.CreateConfigFiles()
    Modules.CondorControl("")
    Modules.RunLocal("Check")
    # time.sleep(6*60*60)
    # Modules.CondorControl("Submit")
    # Modules.CondorControl("")
    # Modules.RunLocal("Check")
    # Modules.RunLocal("Local")
    # Modules.RunLocal("Local", skip=False, nProcess=nProcess, isNice=isNice)
    # Modules.CondorControl("Resubmit")
    # Modules.RunLocal("Check")
    # Modules.CondorControl("List")

    # time.sleep(10*60)
    # Modules.StoreModuleOutput()
    # Modules.DoChecks()
    # Modules.DoChecks()
    # Modules.CreateXml()

    # Modules.SecureMerge(mergeCategory=False)
    # Modules.SecureMerge(mergeCategory=True) # For Preselection only
    # Modules.MakePlots(doPlots=False)
    # Modules.MakePlots()

    # Modules.HowToSpeedCondor()

#################################################
#                                               #
#                   Analysis                    #
#                                               #
#################################################

years = ["UL16","UL17","UL18", "RunII"]
controls = [""]
# controls = ["Preselection", "Selection"]
# controls = ["Preselection"]
# controls = ["Selection"]
# controls = ["SignalRegion"]
# controls = ["PDFReweight"]
controls = ["LeptonIDStudies"]

Modules = ModuleRunner(controls=controls)

# Modules.MakeRunII(Collections, Channels, Systematics, doPlots=False)
# Modules.MakeRunII(Collections, Channels, Systematics, doPlots=True)

# histFolders = ["btag_DeepBoosted_H4qvsQCDmassdep", "btag_DeepBoosted_H4qvsQCDmassdep_cc", "btag_DeepBoosted_H4qvsQCDmassdep_cc1", "btag_DeepBoosted_H4qvsQCDmassdep_cc2", "btag_DeepBoosted_H4qvsQCDmassdep_ccMD"]
# histFolders = ["btag_DeepBoosted_H4qvsQCD_cc", "btag_DeepBoosted_H4qvsQCD_ccMD", "btag_DeepBoosted_H4qvsQCD_ccMD2", "btag_DeepBoosted_H4qvsQCD_ccMD3", "btag_DeepBoosted_H4qvsQCDmassdep_cc"]
# histFolders = ["btag_DeepBoosted_H4qvsQCD_cc"]
# histFolders = ["btag_DeepBoosted_H4qvsQCD_cc","btag_DeepBoosted_H4qvsQCD_ccMD","btag_DeepBoosted_H4qvsQCD_ccMD2"]

# histFolders = ["DeepAk8_H4qvsQCD_massdep", "DeepAk8_ZHccvsQCD_MD", "DeepAk8_HccvsQCD_MD", "DeepAk8_H4qvsQCD_MD",
#                "DeepAk8_H4qvsQCD_massdep_HccvsQCD_MD", "DeepAk8_H4qvsQCD", "DeepAk8_HccvsQCD", "DeepAk8_ZHccvsQCD",
               # "DeepAk8_H4qvsQCD_massdep_HccvsQCD", "DeepAk8_H4qvsQCD_massdep_ZHccvsQCD", "DeepAk8_H4qvsQCD_massdep_ZHccvsQCD_MD"]

# histFolders = ["DeepAk8_HccvsQCD", "DeepAk8_HccvsQCD2", "DeepAk8_ZHccvsQCD_MD", "DeepAk8_ZHccvsQCD_MD2"]

# histFolders = ["DeepAk8_ZHccvsQCD_MD", "DeepAk8_ZHccvsQCD_MD2", "DeepAk8_HccvsQCD", "DeepAk8_H4qvsQCD", "DeepAk8_H4qvsQCD_massdep", "DeepAk8_H4qvsQCD_massdep_HccvsQCD", "tau42"]
histFolders = ["DeepAk8_ZHccvsQCD_MD", "DeepAk8_H4qvsQCD", "DeepAk8_H4qvsQCD_massdep"]
# histFolders = ["DeepAk8_ZHccvsQCD_MD"]

# Modules.RunCommand("PlotNLOCorrections", isPython=True)
# Modules.RunCommand("PlotBTagEfficiencies", isPython=True)
# Modules.RunCommand("PlotHccTagEfficiencies", isPython=True)
# Modules.RunCommand("PlotLeptonIDEfficiency", isPython=True, years=years, Channels=Channels)
# Modules.RunCommand("VerifySignalNormalization", isPython=True, histFolders=histFolders, Channels=["all"], years=["all"], Collections=["all"])
# Modules.RunCommand("ExtractTaggerInformation", isPython=True)
# Modules.RunCommand("ConfusionMatrix", isPython=True)
# Modules.RunCommand("TaggerCutStudy", isPython=True)
# Modules.RunCommand("PlotBkgCuts", isPython=True, Channels=["all"], years=["all"], Collections=["all"])
# Modules.RunCommand("CalculateSignalEfficiencies", histFolders=histFolders)
# Modules.RunCommand("CreateWorkspace", histFolders=histFolders, Channels=["all"], years=["all"], Collections=["all"])
# Modules.RunCommand("CreateWorkspace", histFolders=histFolders, Channels=Channels, years=["RunII"], Collections=Collections)
# Modules.RunCommand("PlotSystematics", isPython=True, years=years, Channels=Channels, histFolders=histFolders, Systematics=Systematics, Collections=Collections)
# Modules.RunCommand("CalculateSystematicEffects", isPython=True, years=years, Channels=Channels, histFolders=histFolders, Collections=Collections)
# Modules.RunCommand("CreateDataCards", isPython=True)
# Modules.RunCommand("CreateDataCards", isPython=True, RunCombine=False)
# Modules.RunCommand("CreateDataCards", isPython=True, RunCombine=True)
# Modules.RunCommand("CompareCombineInputs", isPython=True)
# Modules.RunCommand("PlotLimits")
# Modules.RunCommand("Plotter")
