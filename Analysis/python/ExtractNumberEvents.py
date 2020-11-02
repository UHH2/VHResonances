from Utils import *

'''
Module to extract Number of events from histograms.

'''
# TODO invisible channel not fully implemented yet

class PrintEventNumber(VariablesBase):
    def __init__(self, years=[], Channels=[], Collections=["Puppi"], histFolders=[]):
        VariablesBase.__init__(self)
        self.years = years
        self.Channels = Channels
        self.Collections = Collections
        self.histFolders = histFolders
        self.Samples = filter(lambda x: not self.Signal in x, self.Processes_Year_Dict["2016"]) # 2016 as default. They are all the same
        self.mode = "Selection"
        self.mode = "btag_DeepBoosted_H4qvsQCDmassdep_cc_2_SR"

    def ExtractInfo(self):
        info = {}
        for year, collection, channel, in list(itertools.product(self.years, self.Collections, self.Channels)):
            if "invisible" in channel: continue #TODO
            info[year+collection+channel+"BKG_"+year] = 0
            for sample in self.Samples:
                sample = sample.replace("2016",year)
                if DoControl([""], year+channel+sample, channel, sample): continue
                fname = self.Path_STORAGE+year+"/SignalRegion/"+collection+"/"+channel+"/nominal/"+self.PrefixrootFile+"MC."+sample+"_noTree.root"
                if "DATA" in sample: fname = fname.replace("MC.", "DATA.")
                f_ = ROOT.TFile(fname)
                hist = f_.Get("ZprimeCandidate_"+self.mode+"/sum_event_weights")
                info[year+collection+channel+sample] = hist.GetBinContent(1)
                if not "DATA" in sample:
                    info[year+collection+channel+"BKG_"+year] += hist.GetBinContent(1)
        for sample_ in self.Samples+["BKG"]:
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


if __name__ == '__main__':
    args = parse_arguments()
    years       = args.years if len(args.years)!=0 else ["2016","2017", "2018"]
    histFolders = args.histFolders if len(args.histFolders)!=0 else []
    Channels    = args.Channels if len(args.Channels)!=0 else ["muonchannel","electronchannel"]
    Collections = args.Collections if len(args.Collections)!=0 else ["Puppi"]

    PEN = PrintEventNumber(years=years, Channels=Channels, Collections=Collections, histFolders= histFolders)
    PEN.ExtractInfo()
