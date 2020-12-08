from Utils import *

# import tdrstyle_all as TDR
# TDR.writeExtraText = True
# TDR.extraText  = "Simulation"
# TDR.extraText2 = "Work in progress"

'''
Quick module to verify the content for signal.
It can also be used to verify the efficiencies of the cxx module.

- Mind which module you run it on.
- For each decay, the NEvents, SigEff and SigEff/BR are calculated.
- For the "weights" step, the Eff/BR is expected to be 1.
- The output is print on screen. This is only a quick check.

'''


class VerifySignalNormalization(VariablesBase):
    def __init__(self, year = "RunII", channel = "muonchannel", collection="Puppi",histFolder=None):
        VariablesBase.__init__(self)
        self.year = year
        self.channel = channel.replace("channel","")
        self.collection = collection
        self.histFolder = histFolder

        self.BR = { ""      : 1.,
                    "Zee"   : 1./3,
                    "Zmumu" : 1./3,
                    "Htobb" : 58./100,
                    "HtoWW" : (21.5/100)*(68./100)*(68./100),
                    }
        self.BR["Zelse"] = 1-self.BR["Zee"]-self.BR["Zmumu"]
        self.BR["Helse"] = 1-self.BR["Htobb"]-self.BR["HtoWW"]
        self.Zmodes = ["", "Zee", "Zmumu", "Zelse"]
        self.Hmodes = ["","HtoWW", "Htobb", "Helse"]
        # self.Zmodes = [""]
        # self.Hmodes = [""]

        self.cuts = {
            # "Preselection" : ["nocuts", "weights", "cleaned", "Veto", "NLeptonSel", "NBoostedJet", "DeltaRDiLepton", "JetDiLeptonPhiAngular"],
            # "Selection"    : ["Preselection", "Trigger", "MuonScale", "ZprimeReco", "ZprimeSelection", "PTMassCut", "ScaleFactors"],
            # "SignalRegion" : ["Selection"]+ [histFolder+"_SR"] if histFolder else ["btag_DeepBoosted_H4qvsQCDmassdep_SR"]
            "Preselection" : ["weights"],
            "SignalRegion" : [histFolder+"_SR"] if histFolder else ["btag_DeepBoosted_H4qvsQCDmassdep_SR"]
            }
        self.Modules  = ["Preselection", "SignalRegion"] #Needed to mantain the order

    def ReadFile(self):
        self.Values = {}
        print year, channel, collection, histFolder
        PrintFormattedLine(["massPoint", "SelectionStep", "decay", "NEvents", "Efficiency", "Eff/BR"],15)
        for mp in self.MassPoints:
            massPoint = str(mp)
            self.Values[massPoint] = {}
            for module in self.Modules:
                fname = self.Path_STORAGE+self.year+"/"+module+"/"+self.collection+"/"+self.channel+"channel/nominal/"+self.PrefixrootFile+"MC."+self.Signal+("_inv" if "invisible" in channel else "")+"_M"+massPoint+"_"+self.year+"_noTree.root"
                self.Values[massPoint][module] = {}
                f_ = ROOT.TFile(str(fname))
                for cut in self.cuts[module]:
                    self.Values[massPoint][module][cut] = {}
                    for Zmode in self.Zmodes:
                        for Hmode in self.Hmodes:
                            mode = Zmode+Hmode
                            if Zmode!="" or Hmode!="": mode = "_"+mode
                            h_ = f_.Get("ZprimeCandidate_"+cut+"/sum_event_weights"+mode)
                            bin = h_.GetBinContent(1)
                            mode = mode.replace("_","") if mode!="" else "all"
                            self.Values[massPoint][module][cut][mode] = round(bin,2)
                            BR = self.BR[Zmode]*self.BR[Hmode] if Zmode!="" or Hmode!="" else self.BR[""]
                            eff = bin/self.Values[massPoint]["Preselection"]["weights"]["all"]
                            PrintFormattedLine([massPoint, cut.replace("btag_DeepBoosted_H4qvsQCDmassdep_x3_",""), mode, round(bin,2), round(eff,2), round(eff/BR,2)],15)
                f_.Close()


if __name__ == '__main__':

    years = ["2016","2017","2018"]
    Channels = ["muonchannel","electronchannel", "invisiblechannel"]
    Collections = ["Puppi"]
    histFolders = ["btag_DeepBoosted_H4qvsQCDmassdep_x3"]

    years = ["2016"]
    Channels = ["muonchannel"]

    for year in years:
        for channel in Channels:
            for collection in Collections:
                for histFolder in histFolders:
                    VSN = VerifySignalNormalization(year=year, channel=channel, collection=collection, histFolder=histFolder)
                    VSN.ReadFile()
