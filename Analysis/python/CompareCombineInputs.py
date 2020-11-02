from Utils import *

import tdrstyle_all as TDR
TDR.writeExtraText = True
TDR.extraText = "Work in progress"

colors = {  "h_bkg_inp"  : (ROOT.kFullSquare,   ROOT.kBlue+1),
            "h_sig_inp"  : (ROOT.kFullSquare,   ROOT.kGreen+3),
            "h_bkg_pre"  : (ROOT.kFullDotLarge, ROOT.kRed+1),
            "h_sig_pre"  : (ROOT.kFullDotLarge, ROOT.kGreen+1),
            "h_bkg_post" : (ROOT.kPlus,         ROOT.kViolet),
            "h_sig_post" : (ROOT.kPlus,         ROOT.kCyan+1),
            }

class CompareCombineInputs(ModuleRunnerBase):
    def __init__(self, year="RunII", histFolder="btag_DeepBoosted_H4qvsQCDmassdep_x3", channel="muonchannel"):
        VariablesBase.__init__(self)
        self.year           = year
        self.histoName      = "Zprime_mass_rebin30"
        # self.histoName      = "Zprime_mass_rebin100"
        self.fitFunction    = "Exp_2"
        self.histFolder     = histFolder
        self.histoPath      = self.Path_STORAGE+self.year+"/SignalRegion/Puppi/"+channel+"/nominal/"
        self.channel        = "muon" if "muon" in channel else "electron"
        self.histos         = {}
        self.min = 600 #TODO take it from file
        self.max = 4000 #TODO take it from file
        self.nEventsSR = 1235. #TODO take it from file
        self.xsec_ref = 0.001 #TODO take it from file
        self.normFroPlot = 0.01 #Used for display purposed
        self.outdir = self.Path_ANALYSIS+"Analysis/OtherPlots/CompareCombineInputs/"
        os.system("mkdir -p "+self.outdir)

    def RunAll(self):
        self.LoadHistos()
        self.CreateCanvas()
        self.Plot()
        self.Save()

    def LoadHistos(self):
        file_ = ROOT.TFile(self.histoPath+self.PrefixrootFile+"MC.MC_DY_"+self.year+"_noTree.root")
        self.histos["bkg_SR"] = file_.Get("ZprimeCandidate_"+self.histFolder+"_SR/"+self.histoName).Clone("bkg_SR")
        self.histos["bkg_SR"].SetDirectory(0)
        self.histos["bkg_CR"] = file_.Get("ZprimeCandidate_"+self.histFolder+"_CR/"+self.histoName).Clone("bkg_CR")
        self.histos["bkg_CR"].SetDirectory(0)
        file_.Close()
        dataName = "DATA_SingleMuon" if "muon" in self.channel else "DATA_SingleElectron"
        file_ = ROOT.TFile(self.histoPath+self.PrefixrootFile+"DATA."+dataName+"_"+self.year+"_noTree.root")
        self.histos["DATA_CR"] = file_.Get("ZprimeCandidate_"+self.histFolder+"_CR/"+self.histoName).Clone("DATA_CR")
        self.histos["DATA_CR"].SetDirectory(0)
        self.histos["DATA_SR"] = file_.Get("ZprimeCandidate_"+self.histFolder+"_SR/"+self.histoName).Clone("DATA_SR")
        self.histos["DATA_SR"].SetDirectory(0)
        file_.Close()
        for massPoint in self.SignalSamples:
            if "inv" in massPoint: continue
            if massPoint=="MC_ZprimeToZH_M600": continue
            if massPoint=="MC_ZprimeToZH_M800": continue
            file_ = ROOT.TFile(self.histoPath+self.PrefixrootFile+"MC."+massPoint+"_"+self.year+"_noTree.root")
            self.histos[massPoint] = file_.Get("ZprimeCandidate_"+self.histFolder+"_SR/"+self.histoName).Clone(massPoint)
            self.histos[massPoint].SetDirectory(0)
            self.histos[massPoint].Scale(self.xsec_ref)
            self.histos[massPoint].Scale(self.normFroPlot)
            file_.Close()
            fileCombine = ROOT.TFile(self.histoPath.replace("SignalRegion","").replace("nominal","").replace(self.Path_STORAGE,self.Path_ANALYSIS+"Analysis/Limits/nominal/")+self.histFolder+"/datacards/fitDiagnostics"+massPoint.replace(self.Signal+"_","")+"_"+self.fitFunction+".root")
            self.histos[massPoint+"sign_prefit"] = fileCombine.Get("shapes_prefit/"+self.channel+"_"+self.year+"/total_signal")
            self.histos[massPoint+"sign_prefit"].SetDirectory(0)
            self.histos[massPoint+"sign_prefit"].Scale(self.histos[massPoint+"sign_prefit"].GetBinWidth(1))
            self.histos[massPoint+"sign_prefit"].Scale(self.normFroPlot)
            self.histos["bkg_prefit"] = fileCombine.Get("shapes_prefit/"+self.channel+"_"+self.year+"/total_background")
            self.histos["bkg_prefit"].SetDirectory(0)
            self.histos["bkg_prefit"].Scale(self.histos["bkg_prefit"].GetBinWidth(1))
            fileCombine.Close()



    def CreateCanvas(self):
        TDR.lumi_13TeV  = str(round(float(self.lumi_map["RunII"]["lumi_fb"]),1))+" fb^{-1}"
        self.canv = tdrCanvas("canv", 300, 10000, 1e-3, 1e06, "M(Z')", "Events")
        self.canv.SetLogy(1)
        self.leg = tdrLeg(0.40,0.70,0.89,0.89, 0.030, 42, ROOT.kBlack)
        self.leg.SetNColumns(4)

    def Plot(self):
        for name,hist in self.histos.items():
            col = ROOT.kRed+1 if "Zprime" in name else ROOT.kBlue+1 if "DATA_CR" in name else ROOT.kBlack
            if "bkg_SR" in name : col = ROOT.kGreen+2
            if "bkg_CR" in name : col = ROOT.kViolet+1
            if "bkg_prefit" in name : col = ROOT.kOrange+1
            line = ROOT.kDashed if not "fit" else ROOT.kSolid
            tdrDraw(hist,  "hist", ROOT.kFullDotLarge, col, line, col, 0, col)
            if not self.Signal in name or "bkg" in name: self.leg.AddEntry(hist, name.replace(self.Signal+"_",""), "l")

    def Save(self):
        if self.normFroPlot!=1 : self.leg.AddEntry(0, "Signal x "+str(self.normFroPlot), "")
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
            fName = self.Path_STORAGE+year+"/"+self.nameFolders+self.PrefixrootFile+"MC."+self.sampleName+"_"+year+"_noTree.root"
            if "Preselection" in self.nameFolders and not (year=="2016" and self.sampleName=="MC_TTbar") and not (year!="2016" and (self.sampleName=="MC_WW" or self.sampleName=="MC_WZ" or self.sampleName=="MC_ZZ")): fName = fName.replace(".root","_merge.root")
            if "DATA" in self.sampleName: fName = fName.replace(".MC.", ".DATA.")
            self.files[year] = ROOT.TFile(fName)
            if not self.files[year]:
                raise RuntimeError("fName = "+fName+" not found.")

    def LoadHistos(self):
        if "ZprimeCandidate" in self.histFolder:
            self.hname = self.histFolder+"/Zprime_mass_rebin30"
            if "Preselection" in self.nameFolders:
                self.hname = self.histFolder+"/sum_event_weights"
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
        TDR.lumi_13TeV  = str(round(float(self.lumi_map["RunII"]["lumi_fb"]),1))+" fb^{-1}"
        print self.hname,
        if "sum_event_weights" in self.hname:
            self.canv = tdrCanvas(self.hname, 0.5, 1.5, 1e-01, 1e07, "counting Experiment", "Events")
        elif "jet" in self.hname:
            self.canv = tdrCanvas(self.hname, 0, 200, 1e-01, 1e05, "m(jet)", "Events")
        elif "Zprime" in self.hname:
            self.canv = tdrCanvas(self.hname, 300, 4500, 1e-01, 1e05, "M(Z')", "Events")
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
    # for channel in ["muon","electron"]:
    #     # for sample in ["MC_DY","DATA_Single"]:
    #     for sample in ["MC_DY","DATA_Single", "MC_TTbar", "MC_WW", "MC_WZ", "MC_ZZ"]:
    #         if "DATA" in sample: sample += channel.capitalize()
    #         CompareDistibutionsOverYear(nameFolders="Preselection/Puppi/"+channel+"channel/nominal/", sampleName=sample, histFolder="ZprimeCandidate_JetDiLeptonPhiAngular", extraname=channel).RunAll()
    #         CompareDistibutionsOverYear(nameFolders="Preselection/Puppi/"+channel+"channel/nominal/", sampleName=sample, histFolder="nTopJet_JetDiLeptonPhiAngular", extraname=channel).RunAll()
    #         CompareDistibutionsOverYear(nameFolders="Selection/Puppi/"+channel+"channel/nominal/", sampleName=sample, histFolder="ZprimeCandidate_PTMassCut", extraname=channel).RunAll()
    #         CompareDistibutionsOverYear(nameFolders="SignalRegion/Puppi/"+channel+"channel/nominal/", sampleName=sample, histFolder="ZprimeCandidate_Selection", extraname=channel).RunAll()

    # CompareCombineInputs().RunAll()
    # CompareCombineInputs().RunAll()

    # for histFolder in ["btag_DeepBoosted_H4qvsQCDmassdep", "btag_DeepBoosted_H4qvsQCDmassdep_cc"]:
    #     CompareCombineInputs(histFolder=histFolder).RunAll()
    # CompareCombineInputs(histFolder="btag_DeepBoosted_H4qvsQCDmassdep_cc_2").RunAll()
    CompareCombineInputs(histFolder="btag_DeepBoosted_H4qvsQCDmassdep_cc").RunAll()


    # for channel in channels:
    #     for histFolder in histFolders:
    #             CompareCombineInputs(year,studies,histFolder,channel)

if __name__ == '__main__':
    main()
