from Utils import *

import tdrstyle_all as TDR
TDR.writeExtraText = True
TDR.extraText  = "Simulation"
TDR.extraText2 = "Work in progress"

'''
Module to plot Signal Efficiencies
- Need the Selectionas input
'''

class PlotSignalEfficiencies(VariablesBase):
    def __init__(self):
        VariablesBase.__init__(self)
        self.MassPoints = np.array(self.MassPoints)
        self.Samples = ["M"+str(m) for m in self.MassPoints]
        self.SampleMask = (self.MassPoints>=1000)*(self.MassPoints<=5000)
        # self.SampleMask = (self.MassPoints>=100)*(self.MassPoints<=50000)
        self.SamplesPlot = self.MassPoints[self.SampleMask]
        self.years = self.years+["RunII"]
        # self.Channels = ["invisible"]

        self.outdir = self.Path_ANALYSIS+"Analysis/OtherPlots/SignalEfficiencies/"
        os.system("mkdir -p "+self.outdir)

        self.colors = {
            "Inc"    : (ROOT.kRed+1,"H #rightarrow inclusive"),
            "WW"     : (ROOT.kOrange+1,"H #rightarrow 4q (merged)"),
            "bb"     : (ROOT.kOrange-2,"H #rightarrow b#bar{b}"),
            "cc"     : (ROOT.kBlue+1,"H #rightarrow c#bar{c}"),
            "gg"     : (ROOT.kGreen+3,"H #rightarrow gg"),
            "tautau" : (ROOT.kRed+1,"H #rightarrow #tau#tau"),
            "else"   : (ROOT.kAzure+1,"H #rightarrow 4q (semi-lep)"),
        }

        self.decays = []
        self.Hdecays = ["Inc", "WW", "bb", "cc", "gg", "tautau", "else"]
        self.Zdecays = ["ee", "mumu", "else"]
        for H in self.Hdecays:
            self.decays.append(("Hto"+H).replace("toelse", "else"))
        for Z in ["ee", "mumu", "else"]:
            self.decays.append("Z"+Z)
        for H in self.Hdecays:
            if "Inc" in H: continue
            for Z in ["ee", "mumu", "else"]:
                self.decays.append(("Z"+Z+"Hto"+H).replace("toelse", "else"))
                self.decays.append("Z"+Z)

        self.CutsList = ["Triggers", "Lepton selection", "Angular cuts", "b tag veto", "H4qvsQCD", "ZHccvsQCD"]
        self.Cuts = {
            "weights": {
                "module"   : "Preselection",
                "color"    : ROOT.kBlack,
                "folder"   : "weights",
                },
            "Triggers": {
                "module"   : "Preselection",
                "color"    : ROOT.kRed+1,
                "folder"   : "Trigger",
                },
            "Lepton selection": {
                "module"   : "Preselection",
                "color"    : ROOT.kOrange+1,
                "folder"   : "NLeptonSel",
                },
            "Angular cuts": {
                "module"   : "Selection",
                "color"    : ROOT.kOrange-2,
                "folder"   : "MuonScale",
                },
            "b tag veto": {
                "module"   : "Selection",
                "color"    : ROOT.kBlue+1,
                "folder"   : "ZprimeSelection",
                },
            "PTMassCut": {
                "module"   : "Selection",
                "color"    : ROOT.kBlack,
                "folder"   : "PTMassCut",
                },
            "ScaleFactors": {
                "module"   : "Selection",
                "color"    : ROOT.kBlack,
                "folder"   : "ScaleFactors",
                },
            "Selection": {
                "module"   : "SignalRegion",
                "color"    : ROOT.kBlack,
                "folder"   : "Selection",
                },
            "H4qvsQCD": {
                "module"   : "SignalRegion",
                "color"    : ROOT.kGreen+3,
                "folder"   : "DeepAk8_H4qvsQCD_massdep_SR",
                },
            "ZHccvsQCD": {
                "module"   : "SignalRegion",
                "color"    : ROOT.kAzure+1,
                "folder"   : "DeepAk8_ZHccvsQCD_MD_SR",
                },
        }

        self.values = {}
        self.LoadValues()
        # print self.values.keys()

    def LoadValues(self):
        for year in self.years:
            for collection in self.Collections:
                for channel in self.Channels:
                    for cut in self.Cuts:
                        module = self.Cuts[cut]["module"]
                        folder = self.Cuts[cut]["folder"]
                        for decay in self.decays:
                            unique_name = year+collection+channel+folder+decay
                            if not unique_name in self.values:
                                self.values[unique_name] = np.array([0.]*len(self.Samples))
                            unique_name = year+collection+"chargedlepton"+folder+decay
                            if not unique_name in self.values:
                                self.values[unique_name] = np.array([0.]*len(self.Samples))
                        for index, sample in enumerate(self.Samples):
                            if DoControl([""], year+channel, channel, sample): continue
                            fname = self.Path_STORAGE+year+"/"+module+"/"+collection+"/"+channel+"channel/nominal/uhh2.AnalysisModuleRunner.MC."+self.Signal+("_inv" if "inv" in channel else "")+"_"+sample+"_"+year+"_noTree.root"
                            f_ = ROOT.TFile(fname)
                            for decay in self.decays:
                                unique_name = year+collection+channel+folder+decay
                                h_ = f_.Get("ZprimeCandidate_"+folder+"/sum_event_weights"+("_"+decay if not "Inc" in decay else ""))
                                val = h_.GetBinContent(1)
                                self.values[unique_name][index] += val
                                if not "inv" in channel:
                                    if cut=="weights" and "muon" in channel: continue
                                    self.values[unique_name.replace(channel,"chargedlepton")][index] += val
                            f_.Close()

    def CreateCanvas(self, CanvName, isLow=False, isBR=False):
        self.canv = tdrCanvas(CanvName, 800, 5200, 0.001 if isLow else 0.01 , 1.2 if isBR else 10, "M(Z') (GeV)", "Selection efficiency")
        self.canv.SetLogy(not isBR)
        self.leg = tdrLeg(0.40,0.68,0.89,0.89, 0.030, 42, ROOT.kBlack)
        self.leg.SetNColumns(3)

    def PlotEfficiencyCuts(self):
        plotName = "Cuts_"
        for year in self.years:
            for collection in self.Collections:
                for channel in self.Channels+["chargedlepton"]:
                    self.CreateCanvas(year+collection+channel+"Cuts")
                    grs = []
                    for cut in self.CutsList:
                        folder  = self.Cuts[cut]["folder"]
                        color   = self.Cuts[cut]["color"]
                        eff = self.values[year+collection+channel+folder+"HtoInc"][self.SampleMask]
                        eff /= self.values[year+collection+channel+"weights"+"HtoInc"][self.SampleMask]
                        if "SignalRegion" == self.Cuts[cut]["module"]:
                            eff *= self.values[year+collection+channel+"PTMassCut"+"HtoInc"][self.SampleMask]
                            eff /= self.values[year+collection+channel+"ScaleFactors"+"HtoInc"][self.SampleMask]
                        gr = ROOT.TGraphErrors(len(self.SamplesPlot), array('d',self.SamplesPlot), array('d',eff))
                        gr.SetLineWidth(2)
                        tdrDraw(gr, "lp", ROOT.kFullDotLarge, color, ROOT.kSolid, color, 0, color);
                        self.leg.AddEntry(gr, cut,"lp")
                        grs.append(gr)
                    self.canv.SaveAs(self.outdir+"SignalEfficiencies_"+plotName+year+"_"+collection+"_"+channel+".pdf")

    def PlotEfficiencyHiggsDecays(self, isBR=False, isSR=False):
        plotName = "HiggsDecays_"+("BR_" if isBR else "")+("SR_" if isSR else "")
        cut  = "b tag veto"
        cut  = "Selection"
        cut  = "PTMassCut" if not isSR else "ZHccvsQCD"
        folder  = self.Cuts[cut]["folder"]
        for year in self.years:
            for collection in self.Collections:
                for channel in self.Channels:
                    self.CreateCanvas(year+collection+channel+plotName, isLow= True, isBR = isBR)
                    grs = []
                    for decay in self.Hdecays:
                        if "Inc" in decay and isBR: continue
                        color   = self.colors[decay][0]
                        legName = self.colors[decay][1]
                        decay = ("Hto" if not "else" in decay else "H")+decay
                        decaynorm = decay if isBR else "HtoInc"
                        if isBR and not "Inc" in decay:
                            decaynorm = ("Zee" if "ele" in channel else ("Zmumu" if "muon" in channel else "") ) +decaynorm
                        eff = self.values[year+collection+channel+folder+decay][self.SampleMask]
                        eff /= self.values[year+collection+channel+"weights"+decaynorm][self.SampleMask]
                        gr = ROOT.TGraphErrors(len(self.SamplesPlot), array('d',self.SamplesPlot), array('d',eff))
                        gr.SetLineWidth(2)
                        tdrDraw(gr, "lp", ROOT.kFullDotLarge, color, ROOT.kSolid, color, 0, color);
                        self.leg.AddEntry(gr, legName,"lp")
                        grs.append(gr)
                    self.canv.SaveAs(self.outdir+"SignalEfficiencies_"+plotName+year+"_"+collection+"_"+channel+".pdf")



if __name__ == '__main__':
    PlotBkg = PlotSignalEfficiencies()
    PlotBkg.PlotEfficiencyCuts()
    PlotBkg.PlotEfficiencyHiggsDecays()
    PlotBkg.PlotEfficiencyHiggsDecays(isSR=True)
    PlotBkg.PlotEfficiencyHiggsDecays(isBR = True)
    PlotBkg.PlotEfficiencyHiggsDecays(isBR = True, isSR=True)
