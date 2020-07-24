from Utils import *

import tdrstyle_all as TDR
TDR.writeExtraText = False

'''
Module for visualization of NLO Corrections

- EWK corrections are available for W+jets, Z+jets, gamma+jets samples in the "merged_kfactors_*.root" files under the name “kfactor_monojet_ewk”
- QCD NLO corrections are available for W+jets, Z+jets, gamma+jets in the same files under the name “kfactor_monojet_qcd”. Those are calculated for 2016 samples.
- 2017 version of QCD NLO corrections are available for Z+jets (ll + nunu cases) in the "kfac_*_filter” files.
- QCD NNLO corrections are in the “lindert_qcd_nnlo_sf” file with the following convention:
    - eej -> Z(ll) +jets
    - vvj -> Z(nunu) +jets
    - evj -> W +jets
    - aj -> gamma +jets
- QCD NNLO corrections need to be applied on top of EWK corrections for NLO samples and on top of EWK + QCD NLO corrections for LO samples.
- According to Andreas "I do not apply the NNLO corrections. I have not seen any evidence that they actually improve data/MC agreement. I do not trust them."
- For W+Jets @LO for 2017 and 2018: wjet_dress_monojet or wjet_dress_inclusive in https://github.com/bu-cms/bucoffea/blob/master/bucoffea/data/sf/theory/2017_gen_v_pt_qcd_sf.root
- In the "merged_kfactors_*.root” file, for Z and W + jets, the qcd_ewk histograms are also present: qcd_ewk = QCD * EWK
- taken from https://github.com/bu-cms/bucoffea/tree/master/bucoffea/data/sf/theory
- relative to those studies https://arxiv.org/abs/1705.04664
'''

class PlotNLOCorrections(VariablesBase):
    def __init__(self):
        VariablesBase.__init__(self)
        self.histos = {}
        self.histFolder = os.environ["CMSSW_BASE"]+"/src/UHH2/VHResonances/Analysis/Theory/"
        self.nameXaxis = "p_{T}^{V} (GeV)"
        self.nameYaxis = "Theory Corrections"
        self.Xaxis_min = {"all": 30.0,   "qcd_nnlo": 30.0,
                          "qcd": 100.0,  "ewk": 100.0,  "qcd_ewk": 100.0}
        self.Xaxis_max = {"all": 2000.0, "qcd_nnlo": 6000.0,
                          "qcd": 2100.0, "ewk": 1300.0, "qcd_ewk": 1300.0}
        self.Yaxis_min = {"all": 0.3,    "qcd_nnlo": 1.02,
                          "qcd": 0.3,    "ewk": 0.7,    "qcd_ewk": 0.3}
        self.Yaxis_max = {"all": 1.6,    "qcd_nnlo": 1.20,
                          "qcd": 1.6,    "ewk": 1.1,    "qcd_ewk": 1.6}
        self.color = {"w_qcd_ewk":    (ROOT.kPink+1,      ROOT.kFullTriangleUp),
                      "z_qcd_ewk":    (ROOT.kMagenta+1,   ROOT.kFullTriangleUp),
                      "g_ewk":        (ROOT.kRed+1,       ROOT.kFullCross),
                      "w_ewk":        (ROOT.kOrange+2,    ROOT.kFullCross),
                      "z_ewk":        (ROOT.kOrange-1,    ROOT.kFullCross),
                      "g_qcd":        (ROOT.kViolet+1,    ROOT.kFullCircle),
                      "w_qcd":        (ROOT.kAzure+10,    ROOT.kFullCircle),
                      "z_qcd":        (ROOT.kBlue+1,      ROOT.kFullCircle),
                      "dy_qcd_2017":  (ROOT.kCyan-7,      ROOT.kFullCircle),
                      "znn_qcd_2017": (ROOT.kAzure-7,     ROOT.kFullCircle),
                      "w_qcd_2017":   (ROOT.kBlue-9,      ROOT.kFullCircle),
                      "eej_qcd_nnlo": (ROOT.kGreen+1,     ROOT.kFullTriangleDown),
                      "evj_qcd_nnlo": (ROOT.kGreen+2,     ROOT.kFullTriangleDown),
                      "vvj_qcd_nnlo": (ROOT.kSpring-8,    ROOT.kFullTriangleDown),
                      "aj_qcd_nnlo":  (ROOT.kGreen+3,     ROOT.kFullTriangleDown),
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
        self.canv = tdrCanvas(name, self.Xaxis_min[name], self.Xaxis_max[name], self.Yaxis_min[name], self.Yaxis_max[name], self.nameXaxis, self.nameYaxis, iPeriod=0)
        self.leg = tdrLeg(0.40, 0.70, 0.95, 0.89, 0.025, 42, ROOT.kBlack)
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
            tdrDraw(graph, "CP", self.color[hname][1], self.color[hname]
                    [0], 1, self.color[hname][0], 0, self.color[hname][0])
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
                tdrDraw(graph, "CP", self.color[hname][1], self.color[hname]
                        [0], 1, self.color[hname][0], 0, self.color[hname][0])
                hname = hname.replace("vvj", "Z_{#nu#nu}+Jets").replace("eej",
                                                                        "DY+Jets").replace("aj", "#gamma+Jets").replace("evj", "W+Jets")
                self.leg.AddEntry(graph, hname, "lp")
            self.leg.Draw("same")
            self.canv.SaveAs(self.histFolder+corr+".pdf")
         # if df['is_lo_w']:
         #        all_weights["theory"] = evaluator["qcd_nlo_w_2017"](gen_v_pt) * evaluator["qcd_nnlo_w"](gen_v_pt)
         #    elif df['is_lo_z']:
         #        all_weights["theory"] = evaluator["qcd_nlo_z_2017"](gen_v_pt) * evaluator["qcd_nnlo_z"](gen_v_pt)
         #    elif df['is_lo_g']:
         #        all_weights["theory"] = evaluator["ewk_nlo_g"](gen_v_pt) * evaluator["qcd_nlo_g"](gen_v_pt) * evaluator["qcd_nnlo_g"](gen_v_pt)
         #    else:
         #        all_weights["theory"] = np.ones(df.size)


if __name__ == '__main__':
    NLO = PlotNLOCorrections()
    NLO.LoadHistos()
    NLO.PlotHistos()
