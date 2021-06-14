from Utils import *

import tdrstyle_all as TDR
TDR.writeExtraText = True
TDR.extraText  = "Simulation"
TDR.extraText2 = "Work in progress"

# ForThesis(TDR)

'''
Module to plot Tagger shapes
- Need the Selectionas input
'''

class PlotTaggerShapes(VariablesBase):
    def __init__(self):
        VariablesBase.__init__(self)
        self.Samples = ["DY", "WJets", "ZprimeToZH"]
        # self.Samples = ["ZprimeToZH_M2000"]
        # self.Samples = ["DY", "ZprimeToZH_M2000"]
        # self.Discriminators = ["H4qvsQCD","H4qvsQCD_MD","HccvsQCD","HccvsQCD_MD","ZHccvsQCD","ZHccvsQCD_MD"]
        self.Discriminators = ["H4qvsQCD","ZHccvsQCD_MD"]

        self.color  = {"Background": ROOT.kOrange-2,
                       "Signal":     ROOT.kAzure+2,
                       "DY":               ROOT.kOrange+1,
                       "WJets":            ROOT.kAzure+7,
                       "ZprimeToZH_M2000": ROOT.kAzure+2,
                       "ZprimeToZH_M3000": ROOT.kRed+1,
                       "ZprimeToZH_M4000": ROOT.kGreen+2,

                       "DeepAk8_H4qvsQCD_massdep":              ROOT.kAzure+1,
                       "DeepAk8_ZHccvsQCD_MD":                  ROOT.kGreen+3,
                       "DeepAk8_ZHccvsQCD_MD2":                 ROOT.kOrange+1,
                       "DeepAk8_HccvsQCD_MD":                   ROOT.kRed+1,
                       "DeepAk8_H4qvsQCD_MD":                   ROOT.kAzure-7,
                       "DeepAk8_H4qvsQCD_massdep_HccvsQCD_MD":  ROOT.kOrange+1,
                       "DeepAk8_H4qvsQCD":                      ROOT.kCyan+3,
                       "DeepAk8_HccvsQCD":                      ROOT.kViolet+1,
                       "DeepAk8_HccvsQCD2":                     ROOT.kMagenta+2,
                       "DeepAk8_ZHccvsQCD":                     ROOT.kCyan+1,
                       "DeepAk8_H4qvsQCD_massdep_HccvsQCD":     ROOT.kOrange-1,
                       "DeepAk8_H4qvsQCD_massdep_ZHccvsQCD":    ROOT.kMagenta+2,
                       "DeepAk8_H4qvsQCD_massdep_ZHccvsQCD_MD": ROOT.kGreen+1,
                       }
        self.outdir = self.Path_ANALYSIS+"Analysis/OtherPlots/TaggerShape/"
        os.system("mkdir -p "+self.outdir)

    def LoadHistos(self):
        self.hists = {}
        for sample in ["Background","Signal"]:
            for channel in self.Channels+["lepton"]:
                for disc in self.Discriminators:
                    name_ = sample+channel+disc
                    self.hists[name_] = ROOT.TH1D(name_,name_,30, -0.01, 1.01)
                    self.hists[name_].SetDirectory(0)
        for sample in self.Samples:
            for channel in self.Channels:
                for year in self.years:
                    for collection in self.Collections:
                        filename = self.Path_STORAGE+year+"/Selection/"+collection+"/"+channel+"channel/nominal/work*/uhh2.AnalysisModuleRunner.MC.MC_"+sample+"*.root"
                        if channel == "invisible":
                            filename = filename.replace("ZprimeToZH", "ZprimeToZH_inv")
                        elif sample=="WJets": continue
                        for fname in glob.glob(filename):
                            print fname
                            f_ = ROOT.TFile(fname)
                            t_ = f_.Get("AnalysisTree")
                            for ev in t_:
                                if ev.ZprimeCandidate.size()!=1 :
                                    continue
                                for zp in ev.ZprimeCandidate:
                                    for disc in self.Discriminators:
                                        name = ("Signal" if "ZprimeToZH" in sample else "Background")+channel+disc
                                        val = zp.discriminator("btag_DeepBoosted_"+disc)
                                        self.hists[name].Fill(val,ev.weight_GLP)
                                        if not "inv" in channel:
                                            self.hists[name.replace(channel,"lepton")].Fill(val,ev.weight_GLP)
                            f_.Close()

    def PlotHistos(self):
        self.LoadHistos()
        for disc in self.Discriminators:
            for channel in self.Channels+["lepton"]:
                TDR.lumi_13TeV  = str(round(float(self.lumi_map["RunII"]["lumi_fb"]),1))+" fb^{-1}"
                canv = tdrCanvas(disc+channel, -0.05, 1.05, 0.0001, 20, disc.replace("_MD",""),"A.U.")
                canv.SetLogy(1)
                leg = tdrLeg(0.70,0.70,0.89,0.89, 0.030, 42, ROOT.kBlack);
                for sample in ["Background", "Signal"]:
                    color = self.color[sample]
                    self.hists[sample+channel+disc].Scale(1./self.hists[sample+channel+disc].Integral())
                    tdrDraw(self.hists[sample+channel+disc], "hist", 0, color, 1, color, 1001, color)
                    self.hists[sample+channel+disc].SetFillColorAlpha(color, 0.4)
                    leg.AddEntry(self.hists[sample+channel+disc], sample, "f")
                canv.SaveAs(self.outdir+"TaggerShape_"+disc+"_"+channel+".pdf")



if __name__ == '__main__':
    PlotBkg = PlotTaggerShapes()
    PlotBkg.PlotHistos()
    # args = parse_arguments()
    # years       = args.years if len(args.years)!=0 else ["2016"]
    # histFolders = args.histFolders if len(args.histFolders)!=0 else []
    # Channels    = args.Channels if len(args.Channels)!=0 else ["muonchannel"]
    # Collections = args.Collections if len(args.Collections)!=0 else ["Puppi"]

    # for year in years:
    #     for channel in Channels:
    #         for collection in Collections:
    #             if "inv" in channel: continue
    #             PlotBkg = PlotTaggerShapes(year=year, channel=channel, collection=collection)
    #             PlotBkg.PlotHistos()
