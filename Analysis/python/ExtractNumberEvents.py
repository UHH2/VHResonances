from Utils import *

import tdrstyle_all as TDR
TDR.writeExtraText = True
TDR.extraText  = "Simulation"
TDR.extraText2 = "Work in progress"

ForThesis(TDR)

'''
Module to extract Number of events from histograms.

'''
# TODO invisible channel not fully implemented yet

colors = {"HbbMatchall":       (rt.kOrange+1, rt.kFullTriangleUp,   "Hbb"),
          "HqqMatchall":       (rt.kRed+1,    rt.kFullTriangleDown, "Hqq"),
          "HWWMatchHadronic3": (rt.kBlue+1,   rt.kOpenSquare,       "HWW3q"),
          "HWWMatchHadronic2": (rt.kGray,     rt.kOpenCircle,       "HWW2q"),
          "4qHadronic":        (rt.kAzure+7,  rt.kFullCircle,       "H4q"),
          "Non4qHadronic":     (rt.kGreen+2,  rt.kFullSquare,       "HNon4q"),
          # "HtautauMatchall":   (rt.kAzure+7,  rt.kFullTriangleDown, "H#tau#tau"),
          # "HZZMatchHadronic":  (rt.kRed,      rt.kOpenCircle,       "HZZ4q"),
          # "HWWMatchHadronic":  (rt.kAzure+10, rt.kFullTriangleUp,   "HWW4q"),
          # "HZZMatchHadronic3": (rt.kOrange+1, rt.kFullCross,        "HZZ3q"),
          # "HZZMatchSemiLep":   (rt.kGreen+2,  rt.kFullSquare,       "HZZSemiLep"),
          # "HWWMatchSemiLep":   (rt.kGreen+1,  rt.kFullCircle,       "HWWSemiLep"),
          # "HZZMatchFullLep":   (rt.kGray,     rt.kFullDiamond,      "HZZFullLep"),
          # "HWWMatchFullLep":   (rt.kMagenta,  rt.kFullCircle,       "HWWFullLep"),
}

class PrintEventNumber(VariablesBase):
    def __init__(self, years=[], Channels=[], Collections=["Puppi"], histFolders=[]):
        VariablesBase.__init__(self)
        TDR.lumi_13TeV  = str(round(float(self.lumi_map["RunII"]["lumi_fb"]),1))+" fb^{-1}"
        self.years = years
        self.Channels = Channels
        self.Collections = Collections
        self.histFolders = histFolders
        self.Samples = list(filter(lambda x: not self.Signal in x, self.Processes_Year_Dict["2016"])) # 2016 as default. They are all the same
        self.mode = "Selection"
        # self.mode = "btag_DeepBoosted_H4qvsQCDmassdep_cc_SR"
        # self.mode = "btag_DeepBoosted_H4qvsQCD_cc_SR"
        self.outdir = self.Path_ANALYSIS+"Analysis/SignalEfficiencies/"
        os.system("mkdir -p "+self.outdir)

    def ExtractInfo(self):
        info = {}
        for year, collection, channel, in list(itertools.product(self.years, self.Collections, self.Channels)):
            if "invisible" in channel: continue #TODO
            info[year+collection+channel+"BKG_"+year] = 0
            info[year+collection+channel+"MC_VV_"+year] = 0
            for sample in self.Samples:
                sample = sample.replace("2016",year)
                if DoControl([""], year+channel+sample, channel, sample): continue
                fname = self.Path_STORAGE+year+"/SignalRegion/"+collection+"/"+channel+"/nominal/"+self.PrefixrootFile+"MC."+sample+"_noTree.root"
                if "DATA" in sample: fname = fname.replace("MC.", "DATA.")
                f_ = rt.TFile(fname)
                hist = f_.Get("ZprimeCandidate_"+self.mode+"/sum_event_weights")
                info[year+collection+channel+sample] = hist.GetBinContent(1)
                if not "DATA" in sample:
                    if any(x in sample for x in ["WW", "WZ", "ZZ"]):
                        info[year+collection+channel+"MC_VV_"+year] += hist.GetBinContent(1)
                    info[year+collection+channel+"BKG_"+year] += hist.GetBinContent(1)
        for sample_ in self.Samples+["MC_VV","BKG"]:
            if "MET" in sample_: continue #TODO
            sample_ = sample_.replace("_2016","").replace("_SingleMuon","").replace("_SingleElectron","")
            print sample_, (" "*(20-len(sample_))),
            for year, collection, channel, in list(itertools.product(self.years, self.Collections, self.Channels)):
                sample = sample_.replace("DATA", "DATA_SingleMuon" if "muon" in channel else ("DATA_SingleElectron" if "ele" in channel else ""))
                sample = sample+"_"+year
                if "invisible" in channel: continue #TODO
                if DoControl([""], year+channel+sample, channel, sample): continue
                # var = str(round(info[year+collection+channel+sample],2))+channel[0]+year[-1]
                var = str(round(info[year+collection+channel+sample],2))
                # var2 = str(round(info[year+collection+channel+sample]/info[year+collection+channel+"BKG_"+year]*100,2))
                print "&", var, (" "*(10-len(var))),
                # print "&", var, (" "*(10-len(var))), "&", var2, (" "*(7-len(var2))),
            print "\\\\"

    def LoadVars(self):
        # self.df = pd.read_pickle(self.Path_ANALYSIS+"Analysis/OtherPlots/TaggerInfo/TaggerVariables.pkl")
        # self.df = pd.read_pickle(self.Path_ANALYSIS+"Analysis/OtherPlots/TaggerInfo/TaggerVariables_all.pkl")
        self.df = pd.read_pickle(self.Path_ANALYSIS+"Analysis/OtherPlots/TaggerInfo/TaggerVariables_Preselection.pkl")
        self.df = self.df[self.df["Match"]!="noMatch"]
        self.df = self.df[self.df["Match"]!="gluonMatch"]
        self.df = self.df[self.df["Match"]!="qMatch"]
        self.df["jet_HccvsQCD"] = self.df["jet_Hcc"]/(self.df["jet_Hcc"]+self.df.loc[:,["jet_QCDb","jet_QCDbb", "jet_QCDc", "jet_QCDcc", "jet_QCDqq"]].sum(axis=1))
        print len(self.df),
        self.df = self.df[self.df["jet_HccvsQCD"]>0]
        print len(self.df),
        self.df = self.df[self.df["jet_H4qvsQCD"]>0]
        print len(self.df)
        self.DY = self.df[self.df["sample"].str.contains("DY")]
        self.df = self.df[~self.df["sample"].str.contains("DY")]
        print len(self.df)
        for match in ["noMatch","gluonMatch","qMatch","HMatch","HbbMatch","HWWMatch","HqqMatch","HtautauMatch","HZZMatch","topMatch","tWbMatch","WMatch","WqqMatch","WllMatch","ZMatch","ZqqMatch","ZllMatch"]:
            print match, len(self.df[self.df["Match"]==match])

    def LoadVars2(self):
        vars = {}
        for year, channel, mass in list(itertools.product(self.years, self.Channels, self.MassPoints)):
            sample = self.Signal+"_M"+str(mass)+"_"+year
            print year, channel, sample
            for filename in glob.glob(self.Path_STORAGE+year+"/Preselection/Puppi/"+channel+"/nominal/workdir_Preselection_"+sample+"/*.root"):
                f_ = ROOT.TFile(filename)
                t_ = f_.Get("AnalysisTree")
                for ev in t_:
                    for jet in ev.jetsAk8PuppiSubstructure_SoftDropPuppi:
                        vars.setdefault("sample",[]).append(sample)
                        vars.setdefault("Match",[]).append(rt.MatchingToString(rt.FloatToMatching(jet.get_tag(33))))
                        vars.setdefault("MatchingStatus",[]).append(rt.MatchingStatusToString(rt.FloatToMatching(jet.get_tag(34))))
                        vars.setdefault("jet_QCDb",      []).append(jet.btag_DeepBoosted_probQCDb())
                        vars.setdefault("jet_QCDbb",     []).append(jet.btag_DeepBoosted_probQCDbb())
                        vars.setdefault("jet_QCDc",      []).append(jet.btag_DeepBoosted_probQCDc())
                        vars.setdefault("jet_QCDcc",     []).append(jet.btag_DeepBoosted_probQCDcc())
                        vars.setdefault("jet_QCDqq",     []).append(jet.btag_DeepBoosted_probQCDothers())
                        vars.setdefault("jet_Hcc",       []).append(jet.btag_DeepBoosted_probHcc())
                        vars.setdefault("jet_HbbvsQCD",  []).append(jet.btag_DeepBoosted_HbbvsQCD())
                        vars.setdefault("jet_H4qvsQCD",  []).append(jet.btag_DeepBoosted_H4qvsQCD())
        df = pd.DataFrame(data=vars)
        df.to_pickle(self.Path_ANALYSIS+"Analysis/OtherPlots/TaggerInfo/TaggerVariables_Preselection.pkl")


    def ExtractComposition(self, DB="H4qvsQCD"):
        grs = []
        canv_frac = tdrCanvas("frac_"+DB, 600, 8500, 0, 1.5, "M[GeV]", "Fraction")
        leg_frac = tdrLeg(0.40, 0.70, 0.95, 0.89, 0.025, 42, rt.kBlack)
        canv_comp = tdrCanvas("comp_"+DB, 600, 8500, 0, 1.5, "M[GeV]", "H4qvsQCD")
        leg_comp = tdrLeg(0.40, 0.70, 0.95, 0.89, 0.025, 42, rt.kBlack)
        leg_frac.SetNColumns(2)
        leg_comp.SetNColumns(2)
        DY = np.mean(self.DY[(self.DY["sample"].str.contains("DY")) & (self.DY["jet_"+DB] > 0.001)]["jet_"+DB])
        line_ =  rt.TLine(600, DY, 8500, DY)
        line_.SetLineWidth(2)
        line_.SetLineColor(rt.kBlack)
        canv_comp.cd()
        leg_comp.AddEntry(line_, "DY", "l")
        line_.Draw("same")
        for match in ["all", "4q", "Non4q", "HWWMatch", "HbbMatch",  "HqqMatch", "HtautauMatch", "HZZMatch",]:
            for MS in ["all", "Hadronic", "Hadronic3", "SemiLep", "FullLep"]:
                name = match+MS
                if not name in colors: continue
                print name
                fracs = []
                comps = []
                comps_err = []
                for mass in self.MassPoints:
                    sample = self.Signal+"_M"+str(mass)
                    df = self.df[self.df["sample"].str.contains(sample+"_")]
                    NTOT = len(df)
                    # print name, sample, len(self.df), NTOT
                    if NTOT==0:
                        print "ERROR in ", sample
                        continue
                    if match!="all":
                        if "Non4q" in match:
                            df = df[(df["Match"]!="HWWMatch") & (df["Match"]!="HZZMatch")]
                        elif "4q" in match:
                            df = df[(df["Match"]=="HWWMatch") | (df["Match"]=="HZZMatch")]
                        else:
                            df = df[df["Match"]==match]
                    NMatch = len(df)
                    if NMatch==0:
                        print "ERROR match in ", sample
                        continue
                    if "HWW" in match or "HZZ" in match or ("4q" == match):
                        df = df[df["MatchingStatus"]==MS]
                    frac = NMatch*1./NTOT
                    fracs.append(frac)
                    df = df[df["jet_"+DB] > df["jet_"+DB].quantile(0.2)]
                    df = df[df["jet_"+DB] < df["jet_"+DB].quantile(0.8)]
                    # print name, sample, frac, np.mean(df["jet_"+DB]), np.std(df["jet_"+DB])
                    comps.append(np.mean(df["jet_"+DB]))
                    comps_err.append(np.std(df["jet_"+DB])/2.)
                if len (fracs)==0: continue
                fracs = np.nan_to_num(np.array(fracs))
                comps = np.nan_to_num(np.array(comps))
                comps_err = np.nan_to_num(np.array(comps_err))
                gr_frac = rt.TGraph(len(self.MassPoints), array('d',self.MassPoints), array('d',fracs))
                # gr_frac = rt.TGraphErrors(len(self.MassPoints), array('d',self.MassPoints), array('d',comps), array('d',np.zeros(len(self.MassPoints))), array('d',comps_err))
                # gr_comp = rt.TGraphErrors(len(self.MassPoints), array('d',self.MassPoints), array('d',comps), array('d',np.zeros(len(self.MassPoints))), array('d',comps_err))
                gr_comp = rt.TGraph(len(self.MassPoints), array('d',self.MassPoints), array('d',comps))
                canv_frac.cd()
                tdrDraw(gr_frac, "CP", colors[name][1], colors[name][0], 1, colors[name][0], 0, colors[name][0])
                leg_frac.AddEntry(gr_frac, colors[name][2], "lp")
                canv_comp.cd()
                tdrDraw(gr_comp, "CP", colors[name][1], colors[name][0], 1, colors[name][0], 0, colors[name][0])
                leg_comp.AddEntry(gr_comp, colors[name][2], "lp")
                grs.append(gr_frac)
                grs.append(gr_comp)
        canv_frac.SaveAs(self.outdir+"Fraction_"+DB+".pdf")
        canv_comp.SaveAs(self.outdir+"Composition_"+DB+".pdf")


if __name__ == '__main__':
    args = parse_arguments()
    years       = args.years if len(args.years)!=0 else ["2016","2017", "2018"]
    histFolders = args.histFolders if len(args.histFolders)!=0 else []
    Channels    = args.Channels if len(args.Channels)!=0 else ["muonchannel","electronchannel"]
    Collections = args.Collections if len(args.Collections)!=0 else ["Puppi"]

    PEN = PrintEventNumber(years=years, Channels=Channels, Collections=Collections, histFolders= histFolders)
    PEN.ExtractInfo()
    # PEN.LoadVars2()
    # PEN.LoadVars()
    # PEN.ExtractComposition(DB="H4qvsQCD")
    # PEN.ExtractComposition(DB="HccvsQCD")
