from Utils import *

import tdrstyle_all as TDR
TDR.writeExtraText = True
TDR.extraText  = "Simulation"
TDR.extraText2 = "Work in progress"

ForThesis(TDR)

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

        # self.histFolders = ["DeepAk8_H4qvsQCD_massdep", "DeepAk8_ZHccvsQCD_MD",
        #                     "DeepAk8_HccvsQCD_MD", "DeepAk8_H4qvsQCD_MD", "DeepAk8_H4qvsQCD_massdep_HccvsQCD_MD"]
        self.histFolders = ["Selection", "DeepAk8_HccvsQCD", "DeepAk8_HccvsQCD2", "DeepAk8_ZHccvsQCD_MD", "DeepAk8_ZHccvsQCD_MD2"]

        self.color  = {"Selection":                             ROOT.kBlack,
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
        self.outdir = self.Path_ANALYSIS+"Analysis/SignalEfficiencies/Bkg/"
        os.system("mkdir -p "+self.outdir)

    def LoadHistos(self):
        filename = self.Path_STORAGE+self.year+"/SignalRegion/"+collection+"/"+channel+"/nominal/uhh2.AnalysisModuleRunner.MC.MC_DY_"+self.year+"_noTree.root"
        print filename
        file_ = ROOT.TFile(filename)
        self.histos = {}
        for histFolder in self.histFolders:
            # hname = "ZprimeCandidate_"+histFolder+"_SR/Zprime_mass_rebin30"
            hname = "ZprimeCandidate_"+histFolder+"_SR/Zprime_mass_rebin100"
            if not "DeepAk8" in hname : hname=hname.replace("_SR","")
            print hname
            h_ = file_.Get(hname)
            # h_.Rebin(4)
            h_.SetDirectory(0)
            self.histos[histFolder] = h_
        file_.Close()

    def PlotHistos(self):
        self.LoadHistos()
        TDR.lumi_13TeV  = str(round(float(self.lumi_map[self.year]["lumi_fb"]),1))+" fb^{-1}"
        canv = tdrCanvas("canv", 300, 4000, 1e-03, 6*1e4, "M_{Z'} (GeV)","Events")
        leg = tdrLeg(0.50,0.60,0.89,0.89, 0.030, 42, ROOT.kBlack);
        leg.SetNColumns(2)
        canv.SetLogy(1)
        for histFolder in self.histFolders:
            hist = self.histos[histFolder]
            hist.SetLineWidth(2)
            tdrDraw(hist, "hist", ROOT.kDot, self.color[histFolder], 1, self.color[histFolder], 0, self.color[histFolder])
            leg.AddEntry(hist, histFolder.replace("DeepAk8_", "").replace("vsQCD", "").replace("massdep", "mdc"), "l")
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
