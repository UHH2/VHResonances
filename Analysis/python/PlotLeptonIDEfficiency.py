from Utils import *

import tdrstyle_all as TDR
TDR.writeExtraText = True
TDR.extraText  = "Simulation"
TDR.extraText2 = "Work in progress"

ForThesis(TDR)

'''
Module to plot LeptonID Efficiency

- Need the LeptonIDStudiesModule output as input
- Specify in the main function which years and channels you want to run over.
- Default collection is Puppi.
- Efficiency as function of DeltaR or pt. Choose option
- ID_kincut as reference to calculate efficiency

'''
# TODO invisible channel not fully implemented yet

# map ID name to line color and marker style
colors = {"ele_ID_kincut":                  (ROOT.kBlack,ROOT.kDot),
          "ele_ID_MVA_loose_noIso":         (ROOT.kAzure+10,ROOT.kFullTriangleUp),
          "ele_ID_MVA_loose_iso":           (ROOT.kRed+1,ROOT.kFullCircle),
          "ele_ID_HEEP":                    (ROOT.kAzure+7,ROOT.kFullTriangleDown),
          "ele_ID_veto_noIso":              (ROOT.kRed,ROOT.kOpenCircle),
          "ele_ID_veto":                    (ROOT.kOrange+1,ROOT.kFullCross),
          "ele_ID_loose_noIso":             (ROOT.kBlue+1,ROOT.kOpenSquare),
          "ele_ID_loose":                   (ROOT.kGreen+2,ROOT.kFullSquare),
          "ele_ID_medium_noIso":            (ROOT.kGreen+1,ROOT.kDot),
          "ele_ID_medium":                  (ROOT.kGray,ROOT.kFullDiamond),
          "ele_ID_tight_noIso":             (ROOT.kMagenta,ROOT.kDot),
          "ele_ID_tight":                   (ROOT.kOrange,ROOT.kDot),

          "muon_ID_kincut":                 (ROOT.kBlack,ROOT.kDot),
          "muon_ID_CutBasedIdLoose":        (ROOT.kGreen+2,ROOT.kOpenCircle),
          "muon_ID_CutBasedIdMedium":       (ROOT.kAzure+10,ROOT.kDot),
          "muon_ID_CutBasedIdMediumPrompt": (ROOT.kGray,ROOT.kDot),
          "muon_ID_CutBasedIdTight":        (ROOT.kOrange-2,ROOT.kDot),
          "muon_ID_CutBasedIdGlobalHighPt": (ROOT.kRed+1,ROOT.kOpenSquare),
          "muon_ID_CutBasedIdTrkHighPt":    (ROOT.kOrange+1,ROOT.kFullTriangleUp),
          "muon_ID_MvaLoose":               (ROOT.kBlue+1,ROOT.kFullTriangleDown),
          "muon_ID_MvaMedium":              (ROOT.kGreen,ROOT.kDot),
          "muon_ID_MvaTight":               (ROOT.kMagenta,ROOT.kDot),
          "muon_ID_TkIsoLoose":             (ROOT.kBlue+1,ROOT.kFullTriangleDown),
          "muon_ID_TkIsoTight":             (ROOT.kRed,ROOT.kDot),

}

class PlotLeptonIDEfficiency(VariablesBase):
    def __init__(self,year="2016", collection="Puppi", channel="muonchannel", isDR=True):
        VariablesBase.__init__(self)
        self.year       = year
        TDR.lumi_13TeV  = str(round(float(self.lumi_map[self.year]["lumi_fb"]),1))+" fb^{-1}"
        self.outdir     = self.Path_ANALYSIS+"Analysis/OtherPlots/LeptonID/"
        self.RefID      = "ID_kincut"
        self.HistType   = "LeptonID"
        self.doubleID   = "_2IDs"
        self.collection = collection
        self.channel    = channel
        self.isDR       = isDR
        self.isMuon     = self.channel=="muonchannel"
        self.isEle      = self.channel=="electronchannel"
        self.isInv      = self.channel=="invisiblechannel"
        #Create lists of all IDs and select the ones relative to the channel
        self.IDs        = colors.keys() + [x+self.doubleID for x in colors.keys() if "muon" in x ]
        self.IDs        = list(filter(lambda x: (self.isMuon and "muon" in x) or (self.isEle and "ele" in x), self.IDs))
        self.IDs        = list(reversed(sorted(self.IDs)))
        self.MinSetIDs  = ["ID_MVA", "veto", "loose", "loose_noIso", "HEEP", "IdLoose", "HighPt", "MvaLoose"]
        # self.HistName   = "_DR12" if self.isDR else "_pt12"
        # self.HistName   = "_DR12_sel" if self.isDR else "_pt12_sel"
        self.HistName   = "_DR12_match" if self.isDR else "_pt12_match"
        self.nameXaxis  = ("#DeltaR(e_{1},e_{2})" if self.isEle else "#DeltaR(#mu_{1},#mu_{2})") if isDR else "p_{T} [GeV]"
        self.nameYaxis  = "Efficiency"
        self.Xaxis_min  = 0 if self.isDR else 0
        self.Xaxis_max  = 1 if self.isDR else 5000
        self.Yaxis_min  = 0.
        self.Yaxis_max  = 1.2
        self.histos     = {}

        self.SignalSamples.append(self.MainBkg)
        self.SignalSamples = list(filter(lambda x: (self.isInv and "inv" in x) or (not self.isInv and not "inv" in x), self.SignalSamples))
        self.SignalSamples = list(reversed(sorted(self.SignalSamples)))

        os.system("mkdir -p "+self.outdir)

    def LoadHistos(self):
        for sample in self.SignalSamples:
            filename = self.Path_STORAGE+self.year+"/LeptonIDStudies/"+self.collection+"/"+self.channel+"/nominal/uhh2.AnalysisModuleRunner.MC."+sample+"_"+self.year+"_noTree.root"
            if self.MainBkg in sample and self.year!="RunII": filename = filename.replace(".root","_merge.root")
            file_ = ROOT.TFile(filename)
            for ID in self.IDs:
                hname = self.HistType+"_NLeptonSel/"+ID+self.HistName
                if self.doubleID in ID: hname = hname.replace(self.doubleID, "")+self.doubleID
                self.histos[sample+ID] = file_.Get(hname).Clone(sample+ID)
                self.histos[sample+ID].SetDirectory(0)
            file_.Close()

    def NormHistos(self):
        # First add and then divide all the signal Samples
        for ID in self.IDs:
            for sample in self.SignalSamples:
                if self.MainBkg == sample: continue
                if ID in self.histos:
                    self.histos[ID].Add(self.histos[sample+ID])
                else:
                    self.histos[ID] = self.histos[sample+ID].Clone(ID)

        if self.isMuon: h_den = self.histos["muon_"+self.RefID].Clone("h_den")
        if self.isEle: h_den = self.histos["ele_"+self.RefID].Clone("h_den")
        for ID in self.IDs:
            # self.histos[ID].Divide(h_den)
            self.histos[ID].Divide(self.histos[ID],h_den,1,1,"B")

        # Calculate efficiency for each signal sample separately
        for sample in self.SignalSamples:
            if self.isMuon: h_den = self.histos[sample+"muon_"+self.RefID].Clone("h_den"+sample)
            if self.isEle: h_den = self.histos[sample+"ele_"+self.RefID].Clone("h_den"+sample)
            for ID in self.IDs:
                self.histos[sample+ID].Divide(self.histos[sample+ID],h_den,1,1,"B")

    def ResetCanvas(self, name):
        self.canv = tdrCanvas(name, self.Xaxis_min, self.Xaxis_max, self.Yaxis_min, self.Yaxis_max, self.nameXaxis,self.nameYaxis)
        self.leg = tdrLeg(0.35,0.30,0.95,0.6, 0.035, 42, ROOT.kBlack)
        self.leg.SetNColumns(2 if "simple" in name else 3)

    def PlotHistos(self):
        #Plot all histos
        self.ResetCanvas("canv_"+self.year+self.channel+self.HistName)
        for [ID,hist] in self.histos.items():
            if self.Signal in ID or self.MainBkg in ID: continue
            hist.SetMarkerSize(1.2)
            color = colors[ID.replace(self.doubleID,"")][0]
            marker = colors[ID.replace(self.doubleID,"")][1]
            if "Loose_2IDs" in ID: marker = ROOT.kFullCircle
            if "GlobalHighPt_2IDs" in ID: marker = ROOT.kFullSquare
            tdrDraw(hist, "hist", marker, color, 2 if self.doubleID in ID else 1, color, 0, color)
            self.leg.AddEntry(hist, ID.replace("ID_","").replace("muon_","").replace("ele_","").replace("CutBasedId",""), "l")
        self.canv.SaveAs(self.outdir+self.year+"_"+self.channel+self.HistName+".pdf")

        #Plot only short list of histos
        self.ResetCanvas("canv_simple_"+self.year+self.channel+self.HistName)
        for [ID,hist] in self.histos.items():
            if self.Signal in ID or self.MainBkg in ID: continue
            if all(not x in ID for x in self.MinSetIDs): continue
            if "noIso" in ID and not "loose" in ID : continue
            if "TrkHighPt" in ID and self.doubleID in ID : continue
            if "Tk" in ID and self.doubleID in ID : continue
            if "Mva" in ID and self.doubleID in ID : continue
            hist.SetMarkerSize(1.2)
            color = colors[ID.replace(self.doubleID,"")][0]
            marker = colors[ID.replace(self.doubleID,"")][1]
            if "Loose_2IDs" in ID: marker = ROOT.kFullCircle
            if "GlobalHighPt_2IDs" in ID: marker = ROOT.kFullSquare
            tdrDraw(hist, "P", marker, color, 2 if self.doubleID in ID else 1, color, 0, color)
            self.leg.AddEntry(hist, ID.replace("ID_","").replace("muon_","").replace("ele_","").replace("CutBasedId","").replace(self.doubleID,"+TrkHighPt"), "lp")
        self.canv.SaveAs(self.outdir+self.year+"_"+self.channel+self.HistName+"_simple.pdf")

        if True: # Used only as a cross-check for each sample
            for sample in self.SignalSamples:
                self.ResetCanvas("canv_simple_"+self.year+self.channel+self.HistName+sample)
                for ID in self.IDs:
                    if self.isMuon and not "muon" in ID: continue
                    if self.isEle and not "ele" in ID: continue
                    hist = self.histos[sample+ID]
                    hist.SetMarkerSize(1.2)
                    if all(not x in ID for x in self.MinSetIDs): continue
                    color = colors[ID.replace(self.doubleID,"")][0]
                    marker = colors[ID.replace(self.doubleID,"")][1]
                    if "Loose_2IDs" in ID: marker = ROOT.kFullCircle
                    if "GlobalHighPt_2IDs" in ID: marker = ROOT.kFullSquare
                    tdrDraw(hist, "P", marker, color, 2 if self.doubleID in ID else 1, color, 0, color)
                    self.leg.AddEntry(hist, ID.replace("ID_","").replace("muon_","").replace("ele_","").replace("CutBasedId",""), "lp")
                self.canv.SaveAs(self.outdir+self.year+"_"+self.channel+"_"+sample+self.HistName+"_simple.pdf")




def main():
    args = parse_arguments()
    years    = args.years
    Channels = args.Channels

    # years = ["2016","2017","2018", "RunII"]
    # Channels = ["muonchannel", "electronchannel"]
    # years = ["2016"]
    # Channels = ["muonchannel"]

    for year,channel,isDR in list(itertools.product(years, Channels, [True,False])):
    # for year,channel,isDR in list(itertools.product(years, Channels, [True])):
        PlotSyst = PlotLeptonIDEfficiency(year=year, channel=channel, isDR=isDR)
        PlotSyst.LoadHistos()
        PlotSyst.NormHistos()
        PlotSyst.PlotHistos()


if __name__ == '__main__':
    main()
