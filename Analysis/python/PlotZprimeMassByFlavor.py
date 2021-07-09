from Utils import *

import tdrstyle_all as TDR
TDR.writeExtraText = True
TDR.extraText  = "Simulation"
TDR.extraText2 = "Work in progress"

f_xmin = 1000
f_xmax = 6000
fNorms = {}

'''
Module to plot Tagger shapes
- Need the Selectionas input
'''

class PlotZprimeMassByFlavor(VariablesBase):
    def __init__(self):
        VariablesBase.__init__(self)
        self.SignalSamples = [self.Signal+mode+"_M"+str(mass) for mass in [1400, 1600, 1800, 2000, 2500, 3000, 3500, 4000, 4500, 5000] for mode in ["","_inv"]]
        # self.SignalSamples = [self.Signal+mode+"_M"+str(mass) for mass in [3000] for mode in ["","_inv"]]
        # self.years = ["2018"]
        # self.Channels = ["muon"]
        # self.Channels = ["muon", "invisible"]
        # self.Channels = ["invisible"]
        # self.MatchesList = ["HggMatch", "HccMatch", "H4qMatch", "HbbMatch",  "HVVMatch","HtautauMatch", "noMatch"]
        self.MatchesList = ["noMatch", "HtautauMatch","HVVMatch", "HggMatch", "HccMatch", "H4qMatch", "HbbMatch"]
        self.LegendList = list(reversed(self.MatchesList))

        self.Cuts = ["","SR"]

        self.Matches = {"Inclusive"    : (ROOT.kBlack,"H #rightarrow inclusive"),
                        "noMatch"      : (ROOT.kAzure+1,"not matched"),
                        "H4qMatch"     : (ROOT.kOrange+1,"H #rightarrow 4q (merged)"),
                        "HVVMatch"     : (ROOT.kOrange+2,"H #rightarrow 4q (not merged)"),
                        "HbbMatch"     : (ROOT.kOrange-2,"H #rightarrow b#bar{b}"),
                        "HccMatch"     : (ROOT.kBlue+1,"H #rightarrow c#bar{c}"),
                        "HggMatch"     : (ROOT.kGreen+3,"H #rightarrow gg"),
                        "HtautauMatch" : (ROOT.kRed+1,"H #rightarrow #tau#tau"),
                        "HelseMatch"   : (ROOT.kAzure-7,"H #rightarrow else"),
                        }
        self.outdir = self.Path_ANALYSIS+"Analysis/OtherPlots/SignalShapes/"
        os.system("mkdir -p "+self.outdir)

        self.LoadHistos()
        # print self.hists.keys()

    def LoadHistos(self):
        self.hists = {}
        for year in self.years:
            for channel in self.Channels:
                for collection in self.Collections:
                    for sample in self.SignalSamples:
                        if DoControl([""], year+channel+sample, channel, sample): continue
                        print "Running on", year+channel+sample
                        for cut in self.Cuts:
                            for match in self.Matches.keys():
                                for y_ in [year,"RunII"]:
                                    unique_name = y_+"_"+channel+"_"+collection+"_"+sample+"_"+cut+"_"+match
                                    if unique_name in self.hists: continue
                                    self.hists[unique_name] = ROOT.TH1D(unique_name,unique_name, 100, 0,10000)
                                    self.hists[unique_name].SetDirectory(0)
                                    self.hists[unique_name].SetLineWidth(2)
                        for fname in glob.glob(self.Path_STORAGE+"/"+year+"/Selection/"+collection+"/"+channel+"channel/nominal/workdir_Selection_"+sample+"_"+year+"/*.root"):
                            f_ = ROOT.TFile(fname)
                            t_ = f_.Get("AnalysisTree")
                            for ev in t_:
                                Hdecay = rt.ZprimeDecayToString(int(ev.HDecay))
                                if ev.ZprimeCandidate.size()!=1 :
                                    continue
                                for zp in ev.ZprimeCandidate:
                                    weight = ev.weight_GLP
                                    mass = zp.Zprime_mass()
                                    tagger = zp.discriminator("btag_DeepBoosted_ZHccvsQCD_MD")
                                    match = rt.MatchingToString(rt.FloatToMatching(zp.discriminator("Match")))
                                    MatchingStatus = rt.MatchingStatusToString(rt.FloatToMatching(zp.discriminator("MatchingStatus")))
                                    if match == "HZZMatch" or match == "HWWMatch":
                                        if MatchingStatus == "Hadronic":
                                            match = "H4qMatch"
                                        else: match = "HVVMatch"
                                    if not match in self.MatchesList:
                                        match = "noMatch"
                                    for cut in self.Cuts:
                                        if cut=="SR" and tagger<0.8:continue
                                        unique_name = year+"_"+channel+"_"+collection+"_"+sample+"_"+cut+"_"+match
                                        self.hists[unique_name].Fill(mass,weight)
                                        self.hists[unique_name.replace(match,"Inclusive")].Fill(mass,weight)
                                        self.hists[unique_name.replace(year,"RunII")].Fill(mass,weight)
                                        self.hists[unique_name.replace(year,"RunII").replace(match,"Inclusive")].Fill(mass,weight)

                            f_.Close()


    def PlotHistos(self):
        for year in self.years+["RunII"]:
            for channel in self.Channels:
                isInv = "inv" in channel
                for collection in self.Collections:
                    for sample in self.SignalSamples:
                        for cut in self.Cuts:
                            unique_name = year+"_"+channel+"_"+collection+"_"+sample+"_"+cut
                            if DoControl([""], unique_name, channel, sample): continue
                            stack = ROOT.THStack(unique_name, "")
                            norm = self.hists[unique_name+"_Inclusive"].Integral()
                            TDR.lumi_13TeV  = str(round(float(self.lumi_map[year]["lumi_fb"]),1))+" fb^{-1}"
                            mass = float(sample.replace(self.Signal+("_inv" if isInv else "")+"_M",""))
                            canv = tdrCanvas(unique_name, int(mass/3.), int(mass*1.5), 0.0001, 1 if not isInv else 1, "M(Z') [GeV]" if not isInv else "M_{T}(Z') [GeV]" ,"A.U.", isExtraSpace=True)
                            canv.SetLogy(1)
                            for match in self.MatchesList:
                                color = self.Matches[match][0]
                                hist = self.hists[unique_name+"_"+match]
                                hist.Scale(1./norm)
                                hist_stack = hist.Clone(unique_name+"Stack")
                                tdrDraw(hist, "hist", 9, color, 1, color, 0, color)
                                hist_stack.SetFillColorAlpha(color, 0.99)
                                hist_stack.SetFillStyle(1001)
                                hist_stack.SetLineColor(color)
                                hist_stack.SetLineWidth(0)
                                stack.Add(hist_stack)
                            canv.RedrawAxis()
                            canv.SaveAs(self.outdir+"ZprimeMassByFlavor_"+unique_name+".pdf")
                            # canv = tdrCanvas(unique_name+"Stack", int(mass/3.), int(mass*1.5), 0.0001, 1 if not isInv else 1, "M(Z') [GeV]" if not isInv else "M_{T}(Z') [GeV]" ,"A.U.", isExtraSpace=True)
                            # canv.SetLogy(1)
                            canv = tdrCanvas(unique_name+"Stack", 1000 if not isInv else 800, 3500, 0.0001, 0.25 if not isInv else 0.14, "M(Z') [GeV]" if not isInv else "M_{T}(Z') [GeV]" ,"A.U.", isExtraSpace=True)
                            leg = tdrLeg(0.43, 0.6,0.8,0.89, 0.035, 42, ROOT.kBlack)
                            # leg = tdrLeg(0.19, 0.7 - 0.04*len(self.LegendList),0.45,0.7, 0.035, 42, ROOT.kBlack)
                            # leg = tdrLeg(0.63,0.60,0.92,0.89, 0.035, 42, ROOT.kBlack)
                            stack.Draw("hist same")
                            match = "Inclusive"
                            self.hists[unique_name+"_"+match].Scale(1./norm)
                            tdrDraw(self.hists[unique_name+"_"+match], "hist", 9, self.Matches[match][0], 1, self.Matches[match][0], 0, self.Matches[match][0])
                            leg.AddEntry(self.hists[unique_name+"_"+match], self.Matches[match][1], "l")
                            for match in self.LegendList:
                                hist = self.hists[unique_name+"_"+match]
                                hist.SetFillStyle(1001)
                                leg.AddEntry(hist, self.Matches[match][1], "f")
                            canv.RedrawAxis()
                            canv.SaveAs(self.outdir+"ZprimeMassByFlavor_Stack_"+unique_name+".pdf")



if __name__ == '__main__':
    PlotBkg = PlotZprimeMassByFlavor()
    PlotBkg.PlotHistos()
    # PlotBkg.PlotFuncions()
    # PlotBkg.PlotParameters()
