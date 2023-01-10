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


class PlotSystematics(VariablesBase):
    def __init__(self,year, studies = "nominal", histFolders=[], module="SignalRegion", Channels="", Collections=""):
        VariablesBase.__init__(self)
        self.year = year
        TDR.cms_lumi = self.lumi_map[self.year]['lumiPlot']+' fb^{-1}'
        self.module = module
        self.histFolders = histFolders
        self.SystVar = ["up","down"]
        if self.module=="SignalRegion": self.histFolders = [x+"_SR" for x in self.histFolders]
        if Channels!="": self.Channels = Channels
        if Collections!="": self.Collections = Collections
        self.histos = {}
        self.x_name = {"pt_jet":"p_{T}^{jet} (GeV)", "Zprime_mass_rebin30":"M_{Z'} (GeV)"}
        self.y_name = "Events"
        self.color  = {"nominal":           ROOT.kBlack,
                       "JER_up":            ROOT.kRed+1,
                       "JER_down":          ROOT.kRed,
                       "JEC_up":            ROOT.kGreen+2,
                       "JEC_down":          ROOT.kGreen+1,
                       "MuonScale_up":      ROOT.kBlue+1,
                       "MuonScale_down":    ROOT.kAzure+10,
                       "pu_up":             ROOT.kViolet-3,
                       "pu_down":           ROOT.kMagenta+1,
                       "btag_up":           ROOT.kRed+1,
                       "btag_down":         ROOT.kRed,
                       "prefiring_up":      ROOT.kCyan+2,
                       "prefiring_down":    ROOT.kCyan+1,
                       "id_up":             ROOT.kOrange+1,
                       "id_down":           ROOT.kOrange,
                       "isolation_up":      ROOT.kCyan+3,
                       "isolation_down":    ROOT.kCyan,
                       "tracking_up":       ROOT.kOrange+3,
                       "tracking_down":     ROOT.kOrange+4,
                       "trigger_up":        ROOT.kGreen+2,
                       "trigger_down":      ROOT.kGreen+1,
                       "reco_up":           ROOT.kBlue+1,
                       "reco_down":         ROOT.kAzure+10,
                       }

        self.HistType = {"Preselection": "nTopJet", "Selection": "ZprimeCandidate", "SignalRegion": "ZprimeCandidate"}
        self.HistName = {"Preselection": "pt_jet",  "Selection": "Zprime_mass_rebin30",  "SignalRegion": "Zprime_mass_rebin30"}
        self.outdir = self.Path_ANALYSIS+"Analysis/OtherPlots/Systematics/"
        os.system("mkdir -p "+self.outdir)

    def LoadHistos(self):
        for collection, channel, syst, var, sample in list(itertools.product(self.Collections, self.Channels, self.Systematics+self.Systematics_Scale, self.SystVar, self.Processes_Year_Dict[self.year if self.year!="RunII" else "2016"])):
            isNominalFolder = syst in self.Systematics_Scale
            if not "nominal" in syst: syst += "_"+var

            if DoControl([""], collection+channel+syst+sample, channel, sample): continue
            if "invisible" in channel: continue #TODO
            if isNominalFolder and self.module!="SignalRegion": continue

            commonpath = self.Path_STORAGE+self.year+"/"+self.module+"/"+collection+"/"+channel+"/"

            sample = sample.replace("2016",self.year)
            mode = "MC" if "MC" in sample else "DATA"
            filename = commonpath+syst+"/"+self.PrefixrootFile+mode+"."+sample+"_noTree.root"
            if isNominalFolder: filename = filename.replace(syst,"nominal")
            file_ = ROOT.TFile(filename)
            for histFolder in self.histFolders:
                hname = self.HistType[self.module]
                if isNominalFolder: hname  += "_"+syst.lower()
                hname += "_"+histFolder+"/"+self.HistName[self.module]
                h_ = file_.Get(hname)
                if not h_: continue #TODO
                print file_, h_, filename, hname
                h_.SetDirectory(0)
                if self.Signal in sample:
                    h_.Scale(ROOT.xsec_ref.at("default_value"))
                self.histos.setdefault(collection,{}).setdefault(channel, {}).setdefault(syst, {}).setdefault(sample, {})[histFolder] = h_
            file_.Close()

    def PlotHistos(self):
        self.LoadHistos()
        for histFolder, collection, channel, sample in list(itertools.product(self.histFolders, self.Collections, self.Channels, self.Processes_Year_Dict[self.year if self.year!="RunII" else "2016"])):
            sample = sample.replace("2016",self.year)
            if DoControl([""], histFolder+collection+channel+sample, channel, sample): continue
            if "invisible" in channel: continue #TODO

            mode = "MC" if "MC" in sample else "DATA"
            if self.Signal in sample:
                mean = float(sample.replace("MC_ZprimeToZH_M", "").replace("_"+self.year,""))
                sigma = 5*(mean*0.02+20.)
                ymax = 1.2
            else:
                mean = 2000
                sigma = 2000
                ymax = 15000 if self.module!="SignalRegion" else 500
            ymin = 0.01
            ymax *= round(float(self.lumi_map[self.year]["lumi_fb"])/float(self.lumi_map["RunII"]["lumi_fb"]),1)
            if self.HistName[self.module]=="pt_jet":
                # ymin = 5000
                ymax = 5
                mean /= 2
                sigma *= 2
            for listSyst in [self.Systematics,self.Systematics_Scale]:
                isNominalFolder = listSyst==self.Systematics_Scale
                extraText = "_scale" if isNominalFolder else "_shape"
                canv = tdrCanvas("canv_"+self.year+self.module+collection+channel+sample+histFolder+extraText, mean-sigma, mean+sigma, ymin, ymax, self.x_name[self.HistName[self.module]],self.y_name)
                leg = tdrLeg(0.70,0.40,0.89,0.89, 0.030, 42, ROOT.kBlack);
                # canv.SetLogy(True)
                isNominalFolder = listSyst==self.Systematics_Scale
                for syst_ in listSyst:
                    if "nominal" in syst_: continue
                    for var in self.SystVar:
                        syst = syst_+"_"+var
                        if DoControl([""], syst, channel, sample): continue
                        if isNominalFolder and self.module!="SignalRegion": continue
                        tdrDraw(self.histos[collection][channel][syst][sample][histFolder], "hist", ROOT.kDot, self.color[syst], 1, self.color[syst], 0, self.color[syst])
                        leg.AddEntry(self.histos[collection][channel][syst][sample][histFolder], syst, "l")
                syst = "nominal"
                tdrDraw(self.histos[collection][channel][syst][sample][histFolder], "hist", ROOT.kDot, self.color[syst], 1, self.color[syst], 0, self.color[syst])
                leg.AddEntry(self.histos[collection][channel][syst][sample][histFolder], syst, "l")
                leg.Draw("same")
                canv.SaveAs(self.outdir+"Syst_"+self.year+"_"+self.module+"_"+collection+"_"+channel+"_"+sample+"_"+histFolder+extraText+".pdf")





if __name__ == '__main__':
    args = parse_arguments()
    years       = args.years if len(args.years)!=0 else ["2016"]
    histFolders = args.histFolders if len(args.histFolders)!=0 else ["btag_DeepBoosted_H4qvsQCDmassdep_cc"]
    Channels    = args.Channels if len(args.Channels)!=0 else ["muonchannel"]
    Collections = args.Collections if len(args.Collections)!=0 else ["Puppi"]

    studies = "nominal"
    # module = "Selection"
    module = "SignalRegion"
    # histFolders=["PTMassCut"]
    # histFolders=["btag_DeepBoosted_H4qvsQCDmassdep_x3"]

    for year in years:
        PlotSyst = PlotSystematics(year=year, studies=studies, histFolders=histFolders, module=module, Channels=Channels, Collections=Collections)
        PlotSyst.PlotHistos()
