from Utils import *

import tdrstyle_all
tdrstyle_all.writeExtraText = True
tdrstyle_all.extraText  = "Work in progress"

ROOT.gInterpreter.ProcessLine('#include "'+os.environ["CMSSW_BASE"]+'/src/UHH2/common/include/JetHists.h"')

from array import array

colors = {"2016":       ROOT.kGreen+1,
          "2017":       ROOT.kRed+1,
          "2018":       ROOT.kOrange+1,
          "lepton":     ROOT.kFullCircle,
          "muon":       ROOT.kFullTriangleDown,
          "electron":   ROOT.kFullTriangleUp,

}

class PlotBTagEfficiencies(ModuleRunnerBase):
    def __init__(self):
        VariablesBase.__init__(self)
        tdrstyle_all.lumi_13TeV  = str(round(float(self.lumi_map["RunII"]["lumi_fb"]),1))+" fb^{-1}"
        self.outdir     = self.Path_ANALYSIS+"Analysis/OtherPlots/BTag/"
        self.nameXaxis  = "p_{T} [GeV]"
        self.nameYaxis  = "Efficiency"
        self.Xaxis_min  = 0
        self.Xaxis_max  = 1000
        self.Yaxis_min  = 0.3
        self.Yaxis_max  = 1
        self.lepton     = "lepton"
        self.Cuts       = ["NLeptonSel", "NBoostedJet","JetDiLeptonPhiAngular"]
        self.Modes      = ["BTagMCEffFlavBPassing", "BTagMCEffFlavBTotal"]
        self.histos     = {}

        os.system("mkdir -p "+self.outdir)

    def LoadHistos(self):
        for year in self.years:
            for channel in self.Channels:
                filename = self.Path_STORAGE+year+"/Preselection/Puppi/"+channel+"channel/nominal/uhh2.AnalysisModuleRunner.MC.MC_DY_"+year+"_noTree_merge.root"
                file_ = ROOT.TFile(filename)
                for cut in self.Cuts:
                    for mode in self.Modes:
                        hname = year+channel+cut+mode
                        self.histos[hname] = file_.Get("BTagEff_"+cut+"/"+mode).Clone(hname)
                        self.histos[hname].SetDirectory(0)
                file_.Close()
            for channel in self.Channels:
                for cut in self.Cuts:
                    for mode in self.Modes:
                        hname  = year+channel+cut+mode
                        hname_lep = hname.replace(channel,self.lepton)
                        if hname_lep in self.histos:
                            self.histos[hname_lep].Add(self.histos[hname])
                        else:
                            self.histos[hname_lep] = self.histos[hname].Clone(hname_lep)
                            self.histos[hname_lep].SetDirectory(0)
            for channel in self.Channels+[self.lepton]:
                for cut in self.Cuts:
                    hname = year+channel+cut
                    self.histos[hname] = self.histos[hname+"BTagMCEffFlavBPassing"].Clone(hname)
                    self.histos[hname].Divide(self.histos[hname+"BTagMCEffFlavBTotal"])
                    self.histos[hname].SetDirectory(0)


    def ResetCanvas(self, name="canv"):
        self.canv = tdrCanvas(name, self.Xaxis_min, self.Xaxis_max, self.Yaxis_min, self.Yaxis_max, self.nameXaxis,self.nameYaxis)
        self.leg = tdrLeg(0.40,0.70,0.95,0.89, 0.025, 42, ROOT.kBlack)
        self.leg.SetNColumns(3)

    def PlotHistos(self):
        # For Debug
        # for [hname,hist] in self.histos.items():
        #     self.canv = tdrCanvas("canv_"+hname, self.Xaxis_min, self.Xaxis_max, -2.4, 2.4, self.nameXaxis,"eta")
        #     self.canv.SetRightMargin(0.15)
        #     self.canv.SetLogz(True)
        #     hist.SetMinimum(self.Yaxis_min)
        #     hist.SetMaximum(self.Yaxis_max)
        #     hist.Draw("colz")
        #     self.canv.SaveAs(self.outdir+hname+".pdf")
        self.ResetCanvas()
        for year in self.years:
            for channel in self.Channels+[self.lepton]:
                hname = year+channel+"JetDiLeptonPhiAngular"
                self.histos[hname+"pt"] = ROOT.TH1D(hname+"pt",hname+"pt", ROOT.BTagMCEffBinsPt.size()-1,array('d',list(ROOT.BTagMCEffBinsPt)))
                for x in range(1,self.histos[hname+"pt"].GetNbinsX()+1):
                    self.histos[hname+"pt"].SetBinContent(x,self.histos[hname].GetBinContent(x,1))
                    self.histos[hname+"pt"].SetBinError(x,self.histos[hname].GetBinError(x,1))
                for col in colors:
                    if col in hname:
                        if col.isdigit(): color = colors[col]
                        else : point = colors[col]
                tdrDraw(self.histos[hname+"pt"], "P", point, color, 1, color, 0, color)
                self.leg.AddEntry(self.histos[hname+"pt"], hname.replace("JetDiLeptonPhiAngular",""), "lp")
        self.canv.SaveAs(self.outdir+"Years_lepton.pdf")





def main():
    PlotSyst = PlotBTagEfficiencies()
    PlotSyst.LoadHistos()
    PlotSyst.PlotHistos()


if __name__ == '__main__':
    main()
