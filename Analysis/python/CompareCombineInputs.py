#!/usr/bin/env python
import os, ROOT, glob, subprocess, math
from tdrstyle_all import *
from Utils import *

ROOT.gInterpreter.ProcessLine('#include "'+os.environ["CMSSW_BASE"]+'/src/UHH2/VHResonances/include/constants.hpp"')
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.gStyle.SetOptStat(0)

import tdrstyle_all
tdrstyle_all.writeExtraText = True
tdrstyle_all.extraText  = "Work in progress"

colors = {  "h_bkg_inp"  : (ROOT.kFullSquare,   ROOT.kBlue+1),
            "h_sig_inp"  : (ROOT.kFullSquare,   ROOT.kGreen+3),
            "h_bkg_pre"  : (ROOT.kFullDotLarge, ROOT.kRed+1),
            "h_sig_pre"  : (ROOT.kFullDotLarge, ROOT.kGreen+1),
            "h_bkg_post" : (ROOT.kPlus,         ROOT.kViolet),
            "h_sig_post" : (ROOT.kPlus,         ROOT.kCyan+1),
            }

class CompareCombineInputs(ModuleRunnerBase):
    def __init__(self, year="2016", histFolder="btag_DeepBoosted_H4qvsQCDptdep_x3", nameFolders="SignalRegion/Puppi/muonchannel/nominal/"):
        ModuleRunnerBase.__init__(self,year)
        self.histoName      = "Zprime_mass_rebin30"
        self.fitFunction    = "Exp_2"
        self.histFolder     = histFolder
        self.nameFolders    = nameFolders
        self.channel        = "muon" if "muon" in nameFolders else "electron"
        self.histos         = {}
        self.min = 600 #TODO take it from file
        self.max = 4000 #TODO take it from file
        self.nEventsSR = 1235. #TODO take it from file
        self.xsec_ref = 0.001 #TODO take it from file
        self.outdir = self.Path_ANALYSIS+"Analysis/OtherPlots/CompareCombineInputs/"
        os.system("mkdir -p "+self.outdir)

    def RunAll(self):
        self.LoadHistos()
        self.CreateCanvas()
        self.Plot()
        self.Save()

    def LoadHistos(self):
        # file_ = ROOT.TFile(self.Path_STORAGE+self.year+"/"+self.nameFolders+"uhh2.AnalysisModuleRunner.DATA.DATA_Single"+self.channel.capitalize()+"_"+self.year+"_noTree.root")
        # self.histos["data"] = file_.Get("ZprimeCandidate_"+self.histFolder+"_SR/"+self.histoName).Clone("data")
        # self.histos["data"].SetDirectory(0)
        file_ = ROOT.TFile(self.Path_STORAGE+self.year+"/"+self.nameFolders+"uhh2.AnalysisModuleRunner.MC.MC_DY_"+self.year+"_noTree.root")
        self.histos["bkg_SR"] = file_.Get("ZprimeCandidate_"+self.histFolder+"_SR/"+self.histoName).Clone("bkg_SR")
        self.histos["bkg_SR"].SetDirectory(0)
        self.histos["DATA_CR"] = file_.Get("ZprimeCandidate_"+self.histFolder+"_CR/"+self.histoName).Clone("DATA_CR")
        self.histos["DATA_CR"].SetDirectory(0)
        self.histos["DATA_CR"].Scale(self.nEventsSR/(self.histos["DATA_CR"].Integral(self.histos["DATA_CR"].FindBin(self.min),self.histos["DATA_CR"].FindBin(self.max))))
        for massPoint in self.SignalSamples:
            file_ = ROOT.TFile(self.Path_STORAGE+self.year+"/"+self.nameFolders+"uhh2.AnalysisModuleRunner.MC."+massPoint+"_"+self.year+"_noTree.root")
            self.histos[massPoint] = file_.Get("ZprimeCandidate_"+self.histFolder+"_SR/"+self.histoName).Clone(massPoint)
            self.histos[massPoint].SetDirectory(0)
            self.histos[massPoint].Scale(self.xsec_ref)
            fileCombine = ROOT.TFile(self.Path_ANALYSIS+"Analysis/Limits/nominal/"+self.year+"/"+self.nameFolders.replace("SignalRegion","").replace("nominal","")+self.histFolder+"/datacards/fitDiagnostics"+massPoint.replace(self.signal+"_","")+"_"+self.fitFunction+".root")
            self.histos[massPoint+"sign_prefit"] = fileCombine.Get("shapes_prefit/"+self.channel+"_"+self.year+"/total_signal")
            self.histos[massPoint+"sign_prefit"].SetDirectory(0)
            self.histos[massPoint+"sign_prefit"].Scale(self.histos[massPoint+"sign_prefit"].GetBinWidth(1))
            # self.histos[massPoint+"shapes_fit_s"] = fileCombine.Get("shapes_fit_s/"+self.channel+"_"+self.year+"/total_signal")
            # self.histos[massPoint+"shapes_fit_s"].SetDirectory(0)
            # self.histos[massPoint+"shapes_fit_s"].Scale(self.histos[massPoint+"shapes_fit_s"].GetBinWidth(1))
            self.histos["bkg_prefit"] = fileCombine.Get("shapes_fit_s/"+self.channel+"_"+self.year+"/total_background")
            self.histos["bkg_prefit"].SetDirectory(0)
            self.histos["bkg_prefit"].Scale(self.histos["bkg_prefit"].GetBinWidth(1))



    def CreateCanvas(self):
        tdrstyle_all.lumi_13TeV  = str(self.lumi_fb)+" fb^{-1}"
        self.canv = tdrCanvas("canv", 300, 10000, 1e-3, 1e06, "M(Z')", "Events")
        self.canv.SetLogy(1)
        self.leg = tdrLeg(0.40,0.70,0.89,0.89, 0.030, 42, ROOT.kBlack)
        self.leg.SetNColumns(4)

    def Plot(self):
        for name,hist in self.histos.items():
            col = ROOT.kRed+1 if "Zprime" in name else ROOT.kBlue+1 if "DATA_CR" in name else ROOT.kBlack
            if "bkg_SR" in name : col = ROOT.kGreen+2
            if "bkg_prefit" in name : col = ROOT.kOrange+1
            line = ROOT.kDashed if not "fit" else ROOT.kSolid
            tdrDraw(hist,  "hist", ROOT.kFullDotLarge, col, line, col, 0, col)
            if not self.signal in name or "bkg" in name: self.leg.AddEntry(hist, name.replace(self.signal+"_",""), "l")

    def Save(self):
        self.leg.Draw("same")
        self.canv.SaveAs(self.outdir+"bkg_pred_signals_"+self.histFolder+"_"+self.year+"_"+self.channel+".pdf")


class CompareDistibutionsOverYear(VariablesBase):
    def __init__(self, nameFolders="SignalRegion/Puppi/muonchannel/nominal/", sampleName="MC_DY", histFolder="ZprimeCandidate_Selection", extraname="muon"):
        VariablesBase.__init__(self)
        self.defYear        = "2016"
        self.nameFolders    = nameFolders
        self.sampleName     = sampleName
        self.histFolder     = histFolder
        self.extraname      = extraname
        self.h_inp          = {}
        self.h_ratio        = {}
        self.files          = {}
        self.outdir = self.Path_ANALYSIS+"Analysis/OtherPlots/CompareCombineInputs/"
        os.system("mkdir -p "+self.outdir)

    def RunAll(self):
        self.LoadFiles()
        self.LoadHistos()
        self.CreateCanvas()
        self.Plot()
        self.Save()

    def LoadFiles(self):
        for year in self.years:
            fName = self.Path_STORAGE+year+"/"+self.nameFolders+"uhh2.AnalysisModuleRunner.MC."+self.sampleName+"_"+year+"_noTree.root"
            if "DATA" in self.sampleName: fName = fName.replace(".MC.", ".DATA.")
            self.files[year] = ROOT.TFile(fName)
            if not self.files[year]:
                raise RuntimeError("fName = "+fName+" not found.")

    def LoadHistos(self):
        if "ZprimeCandidate" in self.histFolder:
            self.hname = self.histFolder+"/Zprime_mass_rebin30"
        elif "nTopJet" in self.histFolder:
            self.hname = self.histFolder+"/SDmass_jet"
        else:
            raise ValueError("histFolder = "+self.histFolder+" not expected")
        for year in self.years:
            self.h_inp[year] = self.files[year].Get(self.hname).Clone("inp")
            self.h_ratio[year] = self.files[self.defYear].Get(self.hname).Clone("ratio")
            if not self.h_inp[year]:
                raise RuntimeError("hname = "+self.hname+" not found in "+self.files[year].GetName())
            self.h_inp[year].SetDirectory(0)
            self.h_ratio[year].SetDirectory(0)
            self.h_ratio[year].Divide(self.h_inp[year])
            self.h_ratio[year].Scale(math.sqrt(ModuleRunnerBase(year).lumi_fb/ModuleRunnerBase(self.defYear).lumi_fb))

    def CreateCanvas(self):
        tdrstyle_all.lumi_13TeV  = str(ModuleRunnerBase("RunII").lumi_fb)+" fb^{-1}"
        if "Zprime" in self.hname:
            self.canv = tdrCanvas(self.hname, 300, 4500, 1e-01, 1e05, "M(Z')", "Events")
        elif "jet" in self.hname:
            self.canv = tdrCanvas(self.hname, 0, 200, 1e-01, 1e05, "m(jet)", "Events")
        else:
            raise ValueError("Hist case not expected: "+self.hname)
        self.canv.SetLogy(1)
        self.leg = tdrLeg(0.40,0.70,0.89,0.89, 0.030, 42, ROOT.kBlack)
        self.leg.SetNColumns(3)

    def Plot(self):
        for year in self.years:
            col = ROOT.kRed+1 if year=="2016" else (ROOT.kGreen+1 if year=="2017" else (ROOT.kOrange+1 if year=="2018" else (ROOT.kBlue+1)))
            tdrDraw(self.h_ratio[year],  "hist", ROOT.kFullDotLarge,col, ROOT.kDashed, col, 0, col)
            tdrDraw(self.h_inp[year],  "hist", ROOT.kFullDotLarge,col, ROOT.kSolid, col, 0, col)
            self.leg.AddEntry(self.h_inp[year], year, "l")

    def Save(self):
        self.leg.Draw("same")
        self.canv.SaveAs(self.outdir+self.sampleName+"_"+self.histFolder+"_"+self.extraname+".pdf")



def main():
    for channel in ["muon","electron"]:
        for sample in ["MC_DY","DATA_Single"]:
            if "DATA" in sample: sample += channel.capitalize()
            CompareDistibutionsOverYear(nameFolders="SignalRegion/Puppi/"+channel+"channel/nominal/", sampleName=sample, histFolder="ZprimeCandidate_Selection", extraname=channel).RunAll()

    CompareCombineInputs().RunAll()


    # for channel in channels:
    #     for histFolder in histFolders:
    #             CompareCombineInputs(year,studies,histFolder,channel)

if __name__ == '__main__':
    main()
