import sys,os, time
import numpy as np

from fileManipulation import *
from parallelise import *
from tdrstyle_all import *

sys.path.append(os.environ["CMSSW_BASE"]+"/src/UHH2/VHResonances/Analysis/macros")
from ModuleRunner import *

ROOT.gInterpreter.ProcessLine('#include "'+os.environ["CMSSW_BASE"]+'/src/UHH2/VHResonances/include/constants.hpp"')
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.gStyle.SetOptStat(0)

class PlotSystematics(ModuleRunnerBase):
    def __init__(self,year, studies = "nominal", histFolders=[], module="SignalRegion"):
        ModuleRunnerBase.__init__(self,year)
        self.histFolders = histFolders
        self.module = module
        self.histos = {}
        self.rebin = 25
        self.x_name = {"pt_jet":"p_{T}^{jet} (GeV)", "Zprime_mass_rebin2":"M_{Z'} (GeV)"}
        self.y_name = "Events"
        self.color = {"nominal":    ROOT.kBlack,
                      "JER_up":     ROOT.kRed+1,
                      "JER_down":   ROOT.kRed,
                      "JEC_up":     ROOT.kGreen+2,
                      "JEC_down":   ROOT.kGreen+1,
                      "PU_up":      ROOT.kBlue+1,
                      "PU_down":    ROOT.kAzure+10,
                      }
        self.HistType = {"Preselection": "nTopJet", "Selection": "ZprimeCandidate", "SignalRegion": "ZprimeCandidate"}
        self.HistName = {"Preselection": "pt_jet",  "Selection": "Zprime_mass_rebin2",  "SignalRegion": "Zprime_mass_rebin2"}
        self.outdir = self.Path_ANALYSIS+"Analysis/OtherPlots/Systematics/"
        os.system("mkdir -p "+self.outdir)
    def LoadHistos(self):
        for collection in self.Collections:
            self.histos.setdefault(collection,{})
            for channel in self.Channels:
                self.histos[collection].setdefault(channel, {})
                commonpath = self.Path_STORAGE+self.year+"/"+self.module+"/"+collection+"/"+channel+"channel/"
                # print commonpath
                for syst in self.Systematics+self.Systematics_Scale:
                    if "PU" in syst and self.module!="SignalRegion": continue
                    isNominalFolder = syst in self.Systematics_Scale
                    isNominalSyst = syst=="nominal";
                    self.histos[collection][channel].setdefault(syst, {})
                    for sample in self.Samples_Category[self.year if self.year!="RunII" else "2016"]:
                        if sample =="MC_WW_2018": continue # TODO
                        if sample =="MC_WW_incl_2018": continue # TODO
                        sample = sample.replace("2016",self.year)
                        if "muon" in channel and "SingleElectron" in sample: continue
                        if "electron" in channel and "SingleMuon" in sample: continue
                        mode = "MC" if "MC" in sample else "DATA"
                        filename = commonpath+syst+"/"+self.PrefixrootFile+mode+"."+sample+"_noTree.root"
                        if isNominalFolder: filename = filename.replace(syst,"nominal")
                        file_ = ROOT.TFile(filename)
                        self.histos[collection][channel][syst].setdefault(sample, {})
                        for histFolder in self.histFolders:
                            hname = self.HistType[self.module]
                            if isNominalFolder: hname  += "_"+syst
                            if isNominalFolder: hname += histFolder+"/"+self.HistName[self.module]
                            else: hname += "_"+histFolder+"/"+self.HistName[self.module] # TODO += "_"+histFolder
                            # print filename, hname
                            h_ = file_.Get(hname)
                            h_.SetDirectory(0)
                            if self.HistName[self.module]=="Zprime_mass_rebin2":
                                h_.Rebin(self.rebin)
                            if self.signal in sample:
                                h_.Scale(ROOT.xsec_ref.at("default_value"))
                            # if syst!="nominal":
                            #     h_.Divide(self.histos[collection][channel]["nominal"][sample][histFolder])
                            #     h_.Scale(10000)
                            self.histos[collection][channel][syst][sample][histFolder] = h_
                        file_.Close()
    def PlotHistos(self):
        self.LoadHistos()
        for histFolder in self.histFolders:
            for collection in self.Collections:
                for channel in self.Channels:
                    for sample in self.Samples_Category[self.year if self.year!="RunII" else "2016"]:
                        if sample =="MC_WW_2018": continue # TODO
                        if sample =="MC_WW_incl_2018": continue # TODO
                        sample = sample.replace("2016",self.year)
                        if "muon" in channel and "SingleElectron" in sample: continue
                        if "electron" in channel and "SingleMuon" in sample: continue
                        mode = "MC" if "MC" in sample else "DATA"
                        if self.signal in sample:
                            mean = float(sample.replace("MC_ZprimeToZH_M", "").replace("_"+self.year,""))
                            sigma = 5*(mean*0.02+20.)
                            ymax = 2
                        else:
                            mean = 2000
                            sigma = 2000
                            ymax = 250
                        ymin = 0.01
                        ymax *= self.lumi_fb/round(float(self.lumi_map["RunII"]["lumi_fb"]),1)
                        if self.HistName[self.module]=="pt_jet":
                            # ymin = 5000
                            ymax = 5
                            mean /= 2
                            sigma *= 2
                        # mean = 500
                        # sigma = 500
                        # ymin = 5000
                        # ymax = 200000
                        # mean = 500
                        # sigma = 500
                        canv = tdrCanvas("canv_"+self.year+self.module+collection+channel+sample+histFolder, mean-sigma, mean+sigma, ymin, ymax, self.x_name[self.HistName[self.module]],self.y_name)
                        leg = tdrLeg(0.70,0.70,0.89,0.89, 0.030, 42, ROOT.kBlack);
                        # canv.SetLogy(True)
                        for syst in self.Systematics+self.Systematics_Scale:
                            if "PU" in syst and self.module!="SignalRegion": continue
                            # if syst == "JEC_down": continue
                            # if "PU" in syst: continue
                            # if "JEC" in syst: continue
                            # if "JER" in syst: continue
                            isNominalFolder = syst in self.Systematics_Scale
                            isNominalSyst = syst=="nominal"
                            tdrDraw(self.histos[collection][channel][syst][sample][histFolder], "", ROOT.kDot, self.color[syst], 1, self.color[syst], 0, self.color[syst])
                            leg.AddEntry(self.histos[collection][channel][syst][sample][histFolder], syst, "l")
                        leg.Draw("same")
                        canv.SaveAs(self.outdir+"Syst_"+self.year+"_"+self.module+"_"+collection+"_"+channel+"channel_"+sample+"_"+histFolder+".pdf")


        # # self.histFolders = ["btag_DeepBoosted_H4qvsQCD", "NN", "NN_1","NN_2", "CNN", "tau42"]
        # # self.histFolders = ["btag_DeepBoosted_H4qvsQCD"]
        # self.histFolders = ["btag_DeepBoosted_H4qvsQCDptdep"]
        # self.histFolders=["btag_DeepBoosted_H4qvsQCD", "btag_DeepBoosted_H4qvsQCDp02", "btag_DeepBoosted_H4qvsQCDptdep_x3", "btag_DeepBoosted_H4qvsQCDptdep_x2x3", "btag_DeepBoosted_H4qvsQCDptdep_x1x3", "btag_DeepBoosted_H4qvsQCDmassdep_x3", "btag_DeepBoosted_H4qvsQCDmassdep2_x3", "btag_DeepBoosted_H4qvsQCDmassdep_x2x3", "btag_DeepBoosted_H4qvsQCDmassdep_x1x3", "btag_DeepBoosted_H4qvsQCDmassdep_x1x2"]
        # self.histFolders=["btag_DeepBoosted_H4qvsQCDp02",]
        # # self.histFolders = ["btag_DeepBoosted_H4qvsQCD", "btag_DeepBoosted_H4qvsQCDp2", "btag_DeepBoosted_H4qvsQCDp02", "btag_DeepBoosted_H4qvsQCDpt1000", "btag_DeepBoosted_H4qvsQCDpt1000p2", "btag_DeepBoosted_H4qvsQCDpt1000p02"]
        # # self.histFolders = ["btag_DeepBoosted_H4qvsQCD", "btag_DeepBoosted_H4qvsQCDp2", "btag_DeepBoosted_H4qvsQCDp02"]
        # # self.histFolders = ["btag_DeepBoosted_H4qvsQCD", "btag_DeepBoosted_H4qvsQCDptdep", "btag_DeepBoosted_H4qvsQCDp02"]
        # self.years = ["2016", "2017", "2018", "RunII"]
        # self.channels = ["muonchannel", "electronchannel"]
        # self.collections = ["Puppi"]
        # # self.years = ["2016"]
        # # self.channels = ["muonchannel"]
        # # self.mode = "CB"
        # # self.mode = "Exp_3"
        # self.mode = "Exp_2"
        # self.ResetLists()





if __name__ == '__main__':
    studies = "nominal"
    module = "Selection"
    # module = "SignalRegion"
    years = ["2016","2017","2018", "RunII"]
    years = ["2016"]
    # histFolders=["btag_DeepBoosted_H4qvsQCDmassdep_x3_SR"]
    histFolders=["PTMassCut"]
    for year in years:
        PlotSyst = PlotSystematics(year=year, studies=studies, histFolders=histFolders, module=module)
        PlotSyst.PlotHistos()
