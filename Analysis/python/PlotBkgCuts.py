from Utils import *

import tdrstyle_all as TDR
TDR.writeExtraText = True
TDR.extraText  = "Simulation"
TDR.extraText2 = "Work in progress"

'''
Module to plot Systematics

- Need the Selection or SignalRegion output as input
- Specify in the steer file which years, channels and histFolders you want to run over.

'''
# TODO invisible channel not fully implemented yet
# TODO fix yaxis


class PlotBkgCuts(VariablesBase):
    def __init__(self,year, channel="", collection=""):
        VariablesBase.__init__(self)
        self.year = year
        self.channel = channel
        self.collection = collection
        # self.histFolders = ["btag_DeepBoosted_H4qvsQCDmassdep_2",    "btag_DeepBoosted_H4qvsQCDmassdep_3",    "btag_DeepBoosted_H4qvsQCDmassdep",
        #                     # "btag_DeepBoosted_H4qvsQCDmassdep_bb_2", "btag_DeepBoosted_H4qvsQCDmassdep_bb_3", "btag_DeepBoosted_H4qvsQCDmassdep_bb",
        #                     "btag_DeepBoosted_H4qvsQCDmassdep_cc_2", "btag_DeepBoosted_H4qvsQCDmassdep_cc_3", "btag_DeepBoosted_H4qvsQCDmassdep_cc",
        #                     # "btag_DeepBoosted_H4qvsQCDmassdep_gg_2", "btag_DeepBoosted_H4qvsQCDmassdep_gg_3", "btag_DeepBoosted_H4qvsQCDmassdep_gg"
        #                     "btag_DeepBoosted_H4qvsQCDmassdep2_cc", "btag_DeepBoosted_H4qvsQCDmassdep3_cc",
                            # ]
        # self.histFolders = ["btag_DeepBoosted_H4qvsQCDmassdep", "btag_DeepBoosted_H4qvsQCDmassdep_2",
        #                     "btag_DeepBoosted_H4qvsQCDmassdep_3", "btag_DeepBoosted_H4qvsQCDmassdep_cc",
        #                     "btag_DeepBoosted_H4qvsQCDmassdep_cc_2", "btag_DeepBoosted_H4qvsQCDmassdep_cc_3",
        #                     "btag_DeepBoosted_H4qvsQCDmassdep_cc2", "btag_DeepBoosted_H4qvsQCDmassdep_cc2_2", "btag_DeepBoosted_H4qvsQCDmassdep_cc2_3",
        #                     ]

        self.histFolders = ["Selection","ExtraCleaning", "btag_DeepBoosted_H4qvsQCDmassdep",
                            "btag_DeepBoosted_H4qvsQCDmassdep_cc", "btag_DeepBoosted_H4qvsQCDmassdep_cc1",
                            "btag_DeepBoosted_H4qvsQCDmassdep_cc2", "btag_DeepBoosted_H4qvsQCDmassdep_ccMD"]

        self.color  = {"btag_DeepBoosted_H4qvsQCDmassdep":      ROOT.kRed+1,
                       "btag_DeepBoosted_H4qvsQCDmassdep_bb":   ROOT.kGreen+1,
                       "btag_DeepBoosted_H4qvsQCDmassdep_cc":   ROOT.kViolet+1,
                       "btag_DeepBoosted_H4qvsQCDmassdep_gg":   ROOT.kAzure+1,
                       "btag_DeepBoosted_H4qvsQCDmassdep_2":    ROOT.kBlack,
                       "btag_DeepBoosted_H4qvsQCDmassdep_bb_2": ROOT.kOrange+1,
                       "btag_DeepBoosted_H4qvsQCDmassdep_cc_2": ROOT.kCyan+1,
                       "btag_DeepBoosted_H4qvsQCDmassdep_gg_2": ROOT.kMagenta+1,
                       "btag_DeepBoosted_H4qvsQCDmassdep_3":    ROOT.kBlue,
                       "btag_DeepBoosted_H4qvsQCDmassdep_bb_3": ROOT.kOrange-1,
                       "btag_DeepBoosted_H4qvsQCDmassdep_cc_3": ROOT.kCyan+3,
                       "btag_DeepBoosted_H4qvsQCDmassdep_gg_3": ROOT.kMagenta+2,

                       "btag_DeepBoosted_H4qvsQCDmassdep2_cc":  ROOT.kGreen+1,
                       "btag_DeepBoosted_H4qvsQCDmassdep2":     ROOT.kAzure+1,
                       "btag_DeepBoosted_H4qvsQCDmassdep3_cc":  ROOT.kOrange+1,
                       "btag_DeepBoosted_H4qvsQCDmassdep3":     ROOT.kCyan+1,
                       "btag_DeepBoosted_HccvsQCDmassdep4_cc":  ROOT.kMagenta+1,
                       "btag_DeepBoosted_HccvsQCDmassdep4":     ROOT.kBlue+1,
                       "btag_DeepBoosted_H4qvsQCDmassdep_cc2":  ROOT.kGreen+1,
                       "btag_DeepBoosted_H4qvsQCDmassdep_cc2_2":ROOT.kOrange+1,
                       "btag_DeepBoosted_H4qvsQCDmassdep_cc2_3":ROOT.kOrange-1,
                       "Selection": ROOT.kBlack,
                       "ExtraCleaning": ROOT.kMagenta+1,
                       "btag_DeepBoosted_H4qvsQCDmassdep_cc1": ROOT.kBlue+1,
                       "btag_DeepBoosted_H4qvsQCDmassdep_ccMD": ROOT.kOrange+1,

                       }
        self.outdir = self.Path_ANALYSIS+"Analysis/SignalEfficiencies/Bkg/"
        os.system("mkdir -p "+self.outdir)

    def LoadHistos(self):
        filename = self.Path_STORAGE+self.year+"/SignalRegion/"+collection+"/"+channel+"/nominal/uhh2.AnalysisModuleRunner.MC.MC_DY_"+self.year+"_noTree.root"
        print filename
        file_ = ROOT.TFile(filename)
        self.histos = {}
        for histFolder in self.histFolders:
            hname = "ZprimeCandidate_"+histFolder+"_SR/Zprime_mass_rebin30"
            hname = "ZprimeCandidate_"+histFolder+"_SR/Zprime_mass_rebin100"
            if not "btag" in hname : hname=hname.replace("_SR","")
            print hname
            h_ = file_.Get(hname)
            h_.SetDirectory(0)
            self.histos[histFolder] = h_
        file_.Close()

    def PlotHistos(self):
        self.LoadHistos()
        TDR.lumi_13TeV  = str(round(float(self.lumi_map[self.year]["lumi_fb"]),1))+" fb^{-1}"
        canv = tdrCanvas("canv", 300, 4000, 1e-01, 3*1e3, "M_{Z'} (GeV)","Events")
        leg = tdrLeg(0.50,0.60,0.89,0.89, 0.030, 42, ROOT.kBlack);
        leg.SetNColumns(3)
        canv.SetLogy(1)
        for hname,hist in self.histos.items():
            hist.SetLineWidth(2)
            tdrDraw(hist, "hist", ROOT.kDot, self.color[hname], 1, self.color[hname], 0, self.color[hname])
            leg.AddEntry(hist, hname.replace("btag_DeepBoosted_", "").replace("vsQCDmassdep",""), "l")
        canv.SaveAs(self.outdir+"Bkg_"+self.year+"_"+self.collection+"_"+self.channel+".pdf")



if __name__ == '__main__':
    args = parse_arguments()
    years       = args.years if len(args.years)!=0 else ["2016"]
    histFolders = args.histFolders if len(args.histFolders)!=0 else []
    Channels    = args.Channels if len(args.Channels)!=0 else ["muonchannel"]
    Collections = args.Collections if len(args.Collections)!=0 else ["Puppi"]

    for year in years:
        for channel in Channels:
            for collection in Collections:
                if "inv" in channel: continue
                PlotBkg = PlotBkgCuts(year=year, channel=channel, collection=collection)
                PlotBkg.PlotHistos()
