from Utils import *

import tdrstyle_all as TDR
TDR.writeExtraText = True
TDR.extraText  = "Simulation"
TDR.extraText2 = "Work in progress"

def GetSRScore(jetpt):
    return 5.03e-03+1.7e07*(jetpt*2)**(-3)

class ConfusionMatrix(VariablesBase):
    def __init__(self):
        VariablesBase.__init__(self)
        self.isFast = True
        # self.isFast = False
        self.fName = "TaggerVariables"
        if not self.isFast: self.fName += "_all"
        self.Samples = filter(lambda x: self.MainBkg in x or self.Signal in x, self.Processes_Year_Dict["2016"]) # 2016 as default. They are all the same
        self.outdir = self.Path_ANALYSIS+"Analysis/OtherPlots/TaggerInfo/" #TODO where do we want that?
        os.system("mkdir -p "+self.outdir)

        self.TaggerScores = {
            "T"   : ["jet_Tbqq", "jet_Tbcq"],
            "T1"  : ["jet_Tbq", "jet_Tbc"],
            "W"   : ["jet_Wqq","jet_Wcq"],
            "Z"   : ["jet_Zcc", "jet_Zqq", "jet_Zbb"],
            "Zbb" : ["jet_Zbb"],
            "H"   : ["jet_Hbb", "jet_Hcc", "jet_Hqqqq"],
            "Hbb" : ["jet_Hbb"],
            "Hnb" : ["jet_Hcc", "jet_Hqqqq"],
            "H4q" : ["jet_Hqqqq"],
            "QCD" : ["jet_QCDb","jet_QCDbb", "jet_QCDc", "jet_QCDcc", "jet_QCDqq"],
        }

    def LoadVars(self):
        self.df = pd.read_pickle(self.outdir+self.fName+".pkl")
        # print "Loaded."
        # print self.df
        for score in self.TaggerScores:
            self.df["jet_"+score+"score"] = self.df.loc[:,self.TaggerScores[score]].sum(axis=1)
        self.df["jet_HccvsQCD"] = self.df["jet_Hcc"]/(self.df["jet_Hcc"]+self.df["jet_QCDscore"])
        self.df["jet_HnbvsQCD"] = self.df["jet_Hnbscore"]/(self.df["jet_Hnbscore"]+self.df["jet_QCDscore"])
        self.df["jet_HvsQCD"] = self.df["jet_Hscore"]/(self.df["jet_Hscore"]+self.df["jet_QCDscore"])
        for score in ["bb", "b", "lepb", "uds", "g", "c"]:
            self.df["jet_"+score] = (self.df["subjet_0_"+score] + self.df["subjet_0_"+score])/self.df["subjet_size"]

        for year in ["2016","2017","2018", "RunII"]:
            for sample in ["DY","ZprimeToZH_M1000","ZprimeToZH_M2000","ZprimeToZH_M4000"]:
                df = self.df[(self.df["channel"]=="muon")]
                if year=="RunII":
                    df = df[(df["sample"]=="MC_"+sample+"_2016") | (df["sample"]=="MC_"+sample+"_2017")| (df["sample"]=="MC_"+sample+"_2018")]
                else:
                    df = df[(df["year"]==year) & (df["sample"]=="MC_"+sample+"_"+year)]
                tot = df["weight_GLP"].sum()
                df = df[(df["jet_H4qvsQCD"] > GetSRScore(df["jet_pt"]))]
                SR = df["weight_GLP"].sum()
                df = df[(df["jet_pt"] > 400)]
                SR_pt = df["weight_GLP"].sum()
                print year, sample, tot
                print year, sample, SR/tot
                print year, sample, SR_pt/SR
                print year, sample, df[(df["jet_HccvsQCD"] > (-3.5e-04*df["jet_pt"]*2+0.7))]["weight_GLP"].sum()/SR_pt
                # print (-3.5e-04*df["jet_pt"]*2+0.7)
                for cut in [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]:
                    print year, sample, cut, df[(df["jet_HccvsQCD"] > cut)]["weight_GLP"].sum()/SR_pt
                    print year, sample, cut, df[((df["subjet_0_g"]<cut) | (df["subjet_1_g"]<cut))]["weight_GLP"].sum()/SR_pt


    def CheckDefinitions(self):
        #if you want to filter, otherwise it takes some time!!!
        # self.df = self.df[(self.df["year"]=="2016") & (self.df["channel"]=="muon") & (self.df["sample"]=="MC_ZprimeToZH_M1000_2016")]
        histos = {}
        for x in ["TOT","TvsQCD","WvsQCD","ZvsQCD","ZbbvsQCD","HbbvsQCD","H4qvsQCD"]:
            histos[x] = ROOT.TH1D(x,x,1000,0,2)
            histos[x].SetDirectory(0)
        for ind, df in self.df.iterrows():
            histos["TOT"].Fill(df["jet_Tscore"]+df["jet_T1score"]+df["jet_Wscore"]+df["jet_Zscore"]+df["jet_Hscore"]+df["jet_QCDscore"])
            histos["TvsQCD"].Fill(df["jet_TvsQCD"]/(df["jet_Tscore"]/(df["jet_Tscore"]+df["jet_QCDscore"])) if df["jet_Tscore"]!=0 else 1)
            histos["WvsQCD"].Fill(df["jet_WvsQCD"]/(df["jet_Wscore"]/(df["jet_Wscore"]+df["jet_QCDscore"])) if df["jet_Wscore"]!=0 else 1)
            histos["ZvsQCD"].Fill(df["jet_ZvsQCD"]/(df["jet_Zscore"]/(df["jet_Zscore"]+df["jet_QCDscore"])) if df["jet_Zscore"]!=0 else 1)
            histos["ZbbvsQCD"].Fill(df["jet_ZbbvsQCD"]/(df["jet_Zbb"]/(df["jet_Zbb"]+df["jet_QCDscore"])) if df["jet_Zbb"]!=0 else 1)
            histos["HbbvsQCD"].Fill(df["jet_HbbvsQCD"]/(df["jet_Hbb"]/(df["jet_Hbb"]+df["jet_QCDscore"])) if df["jet_Hbb"]!=0 else 1)
            histos["H4qvsQCD"].Fill(df["jet_H4qvsQCD"]/(df["jet_Hqqqq"]/(df["jet_Hqqqq"]+df["jet_QCDscore"])) if df["jet_Hqqqq"]!=0 else 1)

        for hn in histos:
            canv = tdrCanvas(hn, 0, 2, 1e-10, 1e10, hn,"Events")
            canv.SetLogy(1)
            tdrDraw(histos[hn], "P")
            canv.SaveAs(self.outdir+hn+".pdf")

    @timeit
    def CreatePlots(self):
        self.samples = ["DY", "ZprimeToZH_M1000", "ZprimeToZH_M2000", "ZprimeToZH_M3000"]
        self.scores = ["H4qvsQCD","HbbvsQCD", "HccvsQCD", "HnbvsQCD", "HvsQCD", "bb", "b", "lepb", "uds", "g", "c"]
        self.modes = ["","_SR", "_SR_pt", "_GR"]
        self.ptbins = ["inc", "low", "med", "high"]
        self.histos = {}
        Ybins = sorted([x for x in np.arange(0.1,1,0.005)]+[x for x in np.arange(0.001,0.1,0.0005)])
        for sample in self.samples:
            for score in self.scores:
                for mode in self.modes:
                    for ptbin in self.ptbins:
                        hn = sample+"_"+score+mode+"_"+ptbin
                        self.histos[hn] = ROOT.TH1D(hn,hn,20,0.00001,1)
                        self.histos[hn].SetDirectory(0)
                        hn = sample+"_"+score+mode+"_"+ptbin+"_H4q_2D"
                        self.histos[hn] = ROOT.TH2D(hn,hn,20,0.001,1, 20,0.001,1)
                        self.histos[hn].SetDirectory(0)
                        hn = sample+"_"+score+mode+"_"+ptbin+"_Mass_2D"
                        self.histos[hn] = ROOT.TH2D(hn,hn,92, 700, 9900,20,0.001,1)
                        self.histos[hn].SetDirectory(0)
                        if "vs" in score: continue
                        hn = sample+"_"+score+mode+"_"+ptbin+"_sj_2D"
                        self.histos[hn] = ROOT.TH2D(hn,hn,20,0.001,1,20,0.001,1)
                        self.histos[hn].SetDirectory(0)
                        for sj in ["0","1"]:
                            hn = sample+"_"+score+mode+"_"+ptbin+"_sj"+sj
                            self.histos[hn] = ROOT.TH1D(hn,hn,20,0.001,1)
                            self.histos[hn].SetDirectory(0)

    @timeit
    def PlotVariables(self):
        self.CreatePlots()
        colors = {"DY" : rt.kRed+1, "ZprimeToZH_M1000" : rt.kGreen+1, "ZprimeToZH_M2000" : rt.kBlue+1, "ZprimeToZH_M3000" : rt.kOrange-1}
        year = "2016"
        self.df = self.df[(self.df["year"]==year) & (self.df["channel"]=="muon")]
        for sample in self.samples:
            for ind, df in self.df[(self.df["sample"]=="MC_"+sample+"_"+year)].iterrows():
                # print (df["subjet_0_g"]<0.1 or df["subjet_1_g"]<0.1), df["subjet_0_g"], df["subjet_1_g"]
                weight = df["weight_GLP"]
                pt = df["jet_pt"]
                H4q = df["jet_H4qvsQCD"]
                mass = df["Zprime_mass"]
                for score in self.scores:
                    sj0, sj1, jet_score = (-1,-1,-1)
                    isJet = "vs" in score
                    if isJet:
                        jet_score = df["jet_"+score]
                    else:
                        jet_score = (df["subjet_0_"+score]+df["subjet_1_"+score])/2
                        sj0 = df["subjet_0_"+score]
                        sj1 = df["subjet_1_"+score]
                    for mode in self.modes:
                        if "SR" in mode and H4q < GetSRScore(pt): continue
                        if "pt" in mode and pt<350: continue
                        if "GR" in mode and (df["subjet_0_g"]>0.1 and df["subjet_1_g"]>0.1): continue
                        for ptbin in self.ptbins:
                            if "low" in ptbin and pt>600: continue
                            if "med" in ptbin and (pt<600 or pt>1000): continue
                            if "high" in ptbin and pt<1000: continue
                            self.histos[sample+"_"+score+mode+"_"+ptbin].Fill(jet_score, weight)
                            self.histos[sample+"_"+score+mode+"_"+ptbin+"_H4q_2D"].Fill(H4q,jet_score, weight)
                            self.histos[sample+"_"+score+mode+"_"+ptbin+"_Mass_2D"].Fill(mass,jet_score, weight)
                            if not isJet:
                                self.histos[sample+"_"+score+mode+"_"+ptbin+"_sj_2D"].Fill(sj0,sj1,weight)
                                self.histos[sample+"_"+score+mode+"_"+ptbin+"_sj0"].Fill(sj0,weight)
                                self.histos[sample+"_"+score+mode+"_"+ptbin+"_sj1"].Fill(sj1,weight)
    @timeit
    def SaveCanvas(self):
        canv1D = {}
        leg1D = {}
        for hn,hist in self.histos.items():
            hn = hn.replace("ZprimeToZH_","")
            if "2D" in hn:
                namex = "H4qvsQCD" if not "sj" in hn else hn.split("_")[1]+"_sj0"
                namey = hn.split("_")[1] if not "sj" in hn else hn.split("_")[1]+"_sj1"
                if "Mass" in hn: canv = tdrCanvas(hn, 0.001, 9000, 0.001, 1, "Mass", namey)
                else: canv = tdrCanvas(hn, 0.001, 1, 0.001, 1, namex, namey)
                canv.SetLogz(1)
                # canv.SetLogx(1)
                # canv.SetLogy(1)
                canv.SetRightMargin(0.05)
                if (hist.Integral()!=0): hist.Scale(1./hist.Integral())
                hist.SetMinimum(1e-05 if "SR" in hn else 1e-06)
                tdrDraw(hist, "colz")
                canv.SaveAs(self.outdir+hn+".pdf")
        #     else:
        #         cname = hn.replace(hn.split("_")[0]+"_","")
        #         canv1D[cname].setdefault(tdrCanvas(cname, 0, 1, 1e-3, 10, hn.split("_")[1],"Events")).cd()
        #         canv1D[cname].SetLogy(1)
        #         leg1D[cname].setdefault(tdrLeg(0.40,0.70,0.95,0.89, 0.025, 42, ROOT.kBlack))
        #         if "SR" in hn and not "SR" in mode: continue
        #         if not "SR" in hn and "SR" in mode: continue
        #         if (histos[hn].Integral()!=0): histos[hn].Scale(1./histos[hn].Integral())
        #         histos[hn].SetLineWidth(2)
        #         tdrDraw(histos[hn], "hist", rt.kFullCircle, colors[hn.replace("_SR","")], rt.kSolid if "SR" in hn else rt.kDashed, colors[hn.replace("_SR","")], 0)
        #         leg.AddEntry(histos[hn], hn, "lp")
        # canv.SaveAs(self.outdir+score+"_comparison"+mode+".pdf")
        #
        # for mode in ["","_SR"]:
        #     canv = tdrCanvas(score+mode, 0, 1, 1e-3, 10, score,"Events")
        #     canv.SetLogy(1)
        #     leg = tdrLeg(0.40,0.70,0.95,0.89, 0.025, 42, ROOT.kBlack)
        #     for hn in histos:
        #         if "2D" in hn: continue
        #         if "SR" in hn and not "SR" in mode: continue
        #         if not "SR" in hn and "SR" in mode: continue
        #         if (histos[hn].Integral()!=0): histos[hn].Scale(1./histos[hn].Integral())
        #         histos[hn].SetLineWidth(2)
        #         # tdrDraw(histos[hn], "P", rt.kFullTriangleUp if "SR" in hn else rt.kFullCircle, colors[hn.replace("_SR","")])
        #         tdrDraw(histos[hn], "hist", rt.kFullCircle, colors[hn.replace("_SR","")], rt.kSolid if "SR" in hn else rt.kDashed, colors[hn.replace("_SR","")], 0)
        #         leg.AddEntry(histos[hn], hn, "lp")
        #     canv.SaveAs(self.outdir+score+"_comparison"+mode+".pdf")
        #
        # for hn in histos:
        #     if not "2D" in hn: continue
        #     canv = tdrCanvas(score+hn, 0, 1, 0, 1, "H4qvsQCD", score)
        #     canv.SetLogz(1)
        #     canv.SetRightMargin(0.05)
        #     if (histos[hn].Integral()!=0): histos[hn].Scale(1./histos[hn].Integral())
        #     histos[hn].SetMinimum(1e-05 if "SR" in hn else 1e-06)
        #     tdrDraw(histos[hn], "colz")
        #     canv.SaveAs(self.outdir+score+"_comparison_"+hn+".pdf")

            # nBins = histos["DY"].GetNbinsX()+1
            # for x in range(0, nBins):
            #     print " "*6,
            #     for hn in histos:
            #         print hn.replace("ZprimeToZH_",""), "\t",
            #     print ""
            #     print x, round(histos[hn].GetBinCenter(x),2),
            #     for hn in histos:
            #         print round(histos[hn].Integral(x,nBins),2), "\t",
            #     print ""


def main():
    CM = ConfusionMatrix()
    CM.LoadVars()
    # # CM.CheckDefinitions()
    # CM.PlotVariables()
    # CM.SaveCanvas()


if __name__ == '__main__':
    main()
