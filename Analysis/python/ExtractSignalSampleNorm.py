from Utils import *
import json

'''
Module to extract Number of events from histograms.
PDFReweight module required
'''

class ExtractSignalSampleNorm(VariablesBase):
    def __init__(self):
        VariablesBase.__init__(self)
        self.Samples = list(filter(lambda x: self.Signal in x, self.Processes_Year_Dict["2016"])) # 2016 as default. They are all the same
        self.module = "PDFReweight"
        self.PDF_variations = 100
        self.outdir = self.Path_ANALYSIS+"Analysis/OtherPlots/SignalNormalization/"
        os.system("mkdir -p "+self.outdir)

    @timeit
    def ExtractInfo(self):
        info = {}
        for year, col, ch, sample in list(itertools.product(self.years+["RunII"], self.Collections, self.Channels, self.Samples)):
            if DoControl([""], year+ch+sample, ch, sample): continue
            sample = sample.replace("2016",year)
            info.setdefault(year,{}).setdefault(col,{}).setdefault(ch,{}).setdefault(sample,{})
            fname = self.Path_STORAGE+year+"/"+self.module+"/"+col+"/"+ch+"channel/nominal/"+self.PrefixrootFile+"MC."+sample+"_noTree.root"
            f_ = rt.TFile(fname)
            for mode in ["nocuts", "weights"]:
                for syst in ["nominal", "PDF", "NNPDF", "NNPDF31_lo_as_0130", "PDF4LHC15_nnlo_100"]:
                    isNominal = syst=="nominal"
                    for var in range (1 if isNominal else self.PDF_variations):
                        if isNominal:
                            tag = "nTopJet_"+mode
                        else:
                            tag = "ZprimeCandidate_"+syst+"_"+str(var)+"_"+mode
                        hist = f_.Get(tag+"/sum_event_weights")
                        tag = tag.replace("nTopJet_","").replace("ZprimeCandidate_","")
                        info[year][col][ch][sample][tag] = hist.GetBinContent(1)
            if "inv" in ch: continue
            if info[year][col][ch][sample]["nocuts"]!=info[year][col][ch][sample]["weights"]:
                print "This is not supposed to happen:. Check", year, col, ch, sample, info[year][col][ch][sample]["nocuts"], info[year][col][ch][sample]["weights"]

        # prettydic(info)
        with open(self.outdir+'normalization.json', 'w') as f_:
            json.dump(info, f_)

    @timeit
    def PrintNorm(self, channel="invisible", collection="Puppi"):
        with open(self.outdir+'normalization.json', 'r') as f_:
            info = json.load(f_)
            for year, sample in list(itertools.product(self.years, self.Samples)):
                sample = sample.replace("2016",year)
                if DoControl([""], year+channel+sample, channel, sample): continue
                print year, sample, info[year][collection][channel][sample]["weights"]


if __name__ == '__main__':

    Extractor = ExtractSignalSampleNorm()
    Extractor.ExtractInfo()
    Extractor.PrintNorm()
