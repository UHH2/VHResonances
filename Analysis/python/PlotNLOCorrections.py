from Utils import *

import tdrstyle_all as TDR
TDR.writeExtraText = False

'''
Module for visualization of NLO Corrections

- EWK corrections are available for W+jets, Z+jets, gamma+jets samples in the "merged_kfactors_*.root" files under the name "kfactor_monojet_ewk"
- QCD NLO corrections are available for W+jets, Z+jets, gamma+jets in the same files under the name "kfactor_monojet_qcd". Those are calculated for 2016 samples.
- 2017 version of QCD NLO corrections are available for Z+jets (ll + nunu cases) in the "kfac_*_filter" files.
- QCD NNLO corrections are in the "lindert_qcd_nnlo_sf" file with the following convention:
    - eej -> Z(ll) +jets
    - vvj -> Z(nunu) +jets
    - evj -> W +jets
    - aj -> gamma +jets
- QCD NNLO corrections need to be applied on top of EWK corrections for NLO samples and on top of EWK + QCD NLO corrections for LO samples.
- According to Andreas "I do not apply the NNLO corrections. I have not seen any evidence that they actually improve data/MC agreement. I do not trust them."
- For W+Jets @LO for 2017 and 2018: wjet_dress_monojet or wjet_dress_inclusive in "2017_gen_v_pt_qcd_sf.root"
- In the "merged_kfactors_*.root" file, for Z and W + jets, the qcd_ewk histograms are also present: qcd_ewk = QCD * EWK
- taken from https://github.com/bu-cms/bucoffea/tree/master/bucoffea/data/sf/theory
- relative to those studies https://arxiv.org/abs/1705.04664
'''

class PlotNLOCorrections(VariablesBase):
    def __init__(self):
        VariablesBase.__init__(self)
        self.histos = {}
        self.histFolder = os.environ["CMSSW_BASE"]+"/src/UHH2/VHResonances/Analysis/Theory/"
        self.nameXaxis = "p_{T}^{V} (GeV)"
        self.nameYaxis = "Theory correction factor"
        self.Xaxis_min = {"all": 30.0,   "qcd_nnlo": 30.0,   "qcd": 50.0,   "ewk": 50.0,   "nlo_simple": 100.0,  "nlo_pdfwgt": 50.0,   "qcd_ewk": 100.0}
        self.Xaxis_max = {"all": 2100.0, "qcd_nnlo": 6000.0, "qcd": 2100.0, "ewk": 1300.0, "nlo_simple": 1260.0, "nlo_pdfwgt": 2100.0, "qcd_ewk": 1300.0}
        self.Yaxis_min = {"all": 0.2,    "qcd_nnlo": 1.02,   "qcd": 0.2,    "ewk": 0.65,   "nlo_simple": 0.60,   "nlo_pdfwgt": 0.2,    "qcd_ewk": 0.2}
        self.Yaxis_max = {"all": 1.65,   "qcd_nnlo": 1.20,   "qcd": 1.65,   "ewk": 1.15,   "nlo_simple": 1.60,   "nlo_pdfwgt": 1.65,   "qcd_ewk": 1.65}
        # self.color = {"z_qcd_ewk":    (ROOT.kMagenta+1,   ROOT.kFullTriangleUp),
        #               "w_qcd_ewk":    (ROOT.kPink+1,      ROOT.kFullTriangleUp),
        #               "z_ewk":        (ROOT.kRed+2,       ROOT.kFullSquare),
        #               "w_ewk":        (ROOT.kOrange+1,    ROOT.kFullSquare),
        #               "g_ewk":        (ROOT.kOrange+2,    ROOT.kFullSquare),
        #               "z_qcd":        (ROOT.kAzure+1,     ROOT.kFullCircle),
        #               "w_qcd":        (ROOT.kAzure-7,     ROOT.kFullCircle),
        #               "g_qcd":        (ROOT.kViolet+1,    ROOT.kFullCircle),
        #               "dy_qcd_2017":  (ROOT.kAzure-2,     ROOT.kFullCircle),
        #               "w_qcd_2017":   (ROOT.kAzure+9,     ROOT.kFullCircle),
        #               "znn_qcd_2017": (ROOT.kBlue-9,      ROOT.kFullCircle),
        #               "eej_qcd_nnlo": (ROOT.kGreen+1,     ROOT.kFullTriangleDown),
        #               "evj_qcd_nnlo": (ROOT.kGreen+2,     ROOT.kFullTriangleDown),
        #               "vvj_qcd_nnlo": (ROOT.kSpring-8,    ROOT.kFullTriangleDown),
        #               "aj_qcd_nnlo":  (ROOT.kGreen+3,     ROOT.kFullTriangleDown),
        #               }
        # ["dy_qcd_2017","znn_qcd_2017","w_qcd_2017"] ["z_qcd","z_ewk","w_qcd", "w_ewk"]
        self.color = {"w_qcd":        (ROOT.kGreen+2,     ROOT.kFullSquare),
                      "z_qcd":        (ROOT.kAzure+2,     ROOT.kFullCircle),
                      "w_ewk":        (ROOT.kRed+1,       ROOT.kFullSquare),
                      "z_ewk":        (ROOT.kOrange+1,    ROOT.kFullCircle),
                      "w_qcd_2017":   (ROOT.kGreen+3,     ROOT.kFullSquare),
                      "dy_qcd_2017":  (ROOT.kAzure-3,     ROOT.kFullTriangleUp),
                      "znn_qcd_2017": (ROOT.kAzure+3,     ROOT.kFullTriangleDown),

                      "g_qcd":        (ROOT.kViolet+1,    ROOT.kFullCross),
                      "g_ewk":        (ROOT.kOrange+2,    ROOT.kFullCross),
                      "w_qcd_ewk":    (ROOT.kAzure+2,     ROOT.kFullSquare),
                      "z_qcd_ewk":    (ROOT.kMagenta+1,   ROOT.kFullCircle),
                      "w_qcd_ewk":    (ROOT.kAzure+2,     ROOT.kFullSquare),
                      "eej_qcd_nnlo": (ROOT.kGreen+1,     ROOT.kFullCross),
                      "evj_qcd_nnlo": (ROOT.kGreen+2,     ROOT.kFullCross),
                      "vvj_qcd_nnlo": (ROOT.kSpring-8,    ROOT.kFullCross),
                      "aj_qcd_nnlo":  (ROOT.kGreen+3,     ROOT.kFullCross),
                      }

    def LoadHistos(self):
        for proc in ["g", "w", "z"]:
            file_ = ROOT.TFile(self.histFolder+"merged_kfactors_"+proc+"jets.root")
            for corr in ["ewk", "qcd", "qcd_ewk"]:
                if corr == "qcd_ewk" and proc == "g":
                    continue
                self.histos[proc+"_"+corr] = file_.Get("kfactor_monojet_"+corr)
                self.histos[proc+"_"+corr].SetDirectory(0)
            file_.Close()
        for proc in ["dy", "znn"]:
            file_ = ROOT.TFile(self.histFolder+"kfac_"+proc+"_filter.root")
            self.histos[proc+"_qcd_2017"] = file_.Get("kfac_"+proc+"_filter")
            self.histos[proc+"_qcd_2017"].SetDirectory(0)
            file_.Close()
        file_ = ROOT.TFile(self.histFolder+"lindert_qcd_nnlo_sf.root")
        for proc in ["eej", "evj", "vvj", "aj"]:
            self.histos[proc+"_qcd_nnlo"] = file_.Get(proc)
            self.histos[proc+"_qcd_nnlo"].SetDirectory(0)
        file_.Close()
        file_ = ROOT.TFile(self.histFolder+"2017_gen_v_pt_qcd_sf.root")
        self.histos["w_qcd_2017"] = file_.Get("wjet_dress_inclusive")
        self.histos["w_qcd_2017"].SetDirectory(0)
        file_.Close()

    def ResetCanvas(self, name="all"):
        self.canv = tdrCanvas(name, self.Xaxis_min[name], self.Xaxis_max[name], self.Yaxis_min[name], self.Yaxis_max[name], self.nameXaxis, self.nameYaxis, iPos=0)
        self.leg = tdrLeg(0.40, 0.70, 0.95, 0.89)
        if "all" == name:
            self.leg.SetNColumns(2)
        if "nnlo" in name:
            self.canv.SetLogx(True)

    def PlotHistos(self):
        graphs = []
        self.ResetCanvas()
        for [hname, hist] in self.histos.items():
            graph = ROOT.TGraph(hist)
            graphs.append(graph)
            tdrDraw(graph, "CP", self.color[hname][1], self.color[hname][0], 1, self.color[hname][0], 0, self.color[hname][0])
            self.leg.AddEntry(graph, hname, "lp")
        self.leg.Draw("same")
        self.canv.SaveAs(self.histFolder+"all.pdf")
        for corr in ["ewk", "qcd", "qcd_nnlo", "qcd_ewk"]:
            self.ResetCanvas(corr)
            for [hname, hist] in self.histos.items():
                if hname[-len(corr)::] != corr and not ("qcd" == corr and "2017" in hname):
                    continue
                if "qcd" in hname and "ewk" in hname and "qcd_ewk" != corr:
                    continue
                graph = ROOT.TGraph(hist)
                graphs.append(graph)
                graph.SetLineWidth(2)
                tdrDraw(graph, "CP", self.color[hname][1], self.color[hname][0], 1, self.color[hname][0], 0, self.color[hname][0])
                hname = hname.replace("vvj", "Z_{#nu#nu}+Jets").replace("eej","DY+Jets").replace("aj", "#gamma+Jets").replace("evj", "W+Jets")
                self.leg.AddEntry(graph, hname, "lp")
            self.leg.Draw("same")
            self.canv.SaveAs(self.histFolder+corr+".pdf")

        corr = "nlo_simple"
        self.ResetCanvas(corr)
        self.leg  = tdrLeg(0.50, 0.70, 0.77, 0.89, 0.05)
        self.leg2 = tdrLeg(0.70, 0.70, 0.95, 0.89, 0.05)
        tdrHeader(self.leg, "QCD",  textSize = 0.05)
        tdrHeader(self.leg2, "EWK", textSize = 0.05)

        for hname in ["w_qcd", "w_ewk", "z_qcd","z_ewk"]:
            hist = self.histos[hname]
            graph = ROOT.TGraph(hist)
            graphs.append(graph)
            graph.SetLineWidth(5)
            graph.SetMarkerSize(1.7)
            color = self.color[hname][0]
            linestyle = ROOT.kDashed if "w_" in hname else ROOT.kSolid
            tdrDraw(graph, "CP", self.color[hname][1], color, linestyle, color)
            nleg = hname.replace("_", " ").upper().split()[0]+"+jets"
            if "qcd" in hname:
                self.leg.AddEntry(graph, nleg, "lp")
            else:
                self.leg2.AddEntry(graph, nleg, "lp")
        self.canv.SaveAs(self.histFolder+corr+".pdf")


        corr = "nlo_pdfwgt"
        self.ResetCanvas(corr)
        self.leg = tdrLeg(0.65, 0.65, 0.95, 0.89, 0.05)
        tdrHeader(self.leg, "QCD",  textSize = 0.05)
        for hname in ["w_qcd_2017","dy_qcd_2017","znn_qcd_2017"]:
            hist = self.histos[hname]
            graph = ROOT.TGraph(hist)
            graphs.append(graph)
            graph.SetLineWidth(5)
            graph.SetMarkerSize(1.7)
            color = self.color[hname][0]
            linestyle = ROOT.kDashed if "w_" in hname else ROOT.kSolid
            markerstyle = ROOT.kFullSquare if "w_" in hname else (ROOT.kFullTriangleUp if "dy_" in hname else ROOT.kFullTriangleDown)
            tdrDraw(graph, "CP", markerstyle, color, linestyle, color)
            nleg = hname.replace("_", " ").replace("dy", "Z(ll)").replace("w","W").replace("znn", "Z(#nu#nu)").split()[0]+"+jets"
            self.leg.AddEntry(graph, nleg, "lp")

        self.canv.SaveAs(self.histFolder+corr+".pdf")

if __name__ == '__main__':
    NLO = PlotNLOCorrections()
    NLO.LoadHistos()
    NLO.PlotHistos()
