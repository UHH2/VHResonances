from Utils import *

import tdrstyle_all as TDR
import math

from ROOT import *

"""
Script to create ratio plots relative to the muonchannel for the DeepBoosted
variables.
Takes as an input the JetFlavour module output.
Loops over "match"  for the TopJets and different histoNames (name of the
DeepBoosted variable).
"""

TDR.writeExtraText = True
TDR.extraText = "Work in progress"

colors = {  "invisible"  : ROOT.kBlue+1,
            "muon"  : ROOT.kRed+1,
            "electron"  : ROOT.kOrange+1,
            }

class CompareBTagVariables(ModuleRunnerBase):
    def __init__(self, year="RunII", match="incl", histoName="btag_MassDecorrelatedDeepBoosted_ZHccvsQCD"):
        print "** CompareBTagVariables**"
        VariablesBase.__init__(self)
        self.year           = year
        self.match          = (match if match!="" else "incl")
        self.histFolders    = ["nocuts", "weights", "Trigger", "HEM", "cleaned", "Veto", "NLeptonSel", "NBoostedJet", "METCut", "DeltaRDiLepton", "JetDiLeptonPhiAngular", "QCDRejection"]
        self.histoPath      = self.Path_STORAGE+self.year+"/JetFlavour/Puppi/"
        self.histoName      = histoName
        self.histos         = {}
        self.canv_ratio     = {}
        self.legend         = {}
        self.outdir = self.Path_ANALYSIS+"Analysis/OtherPlots/CompareBTagVariables/"
        self.samples = ["MC_DY", "MC_ZprimeToZH_M2000", "MC_ZprimeToZH_inv_M2000"]
        os.system("mkdir -p "+self.outdir)
        print "Will store the output to:", self.outdir

    def RunAll(self):
        self.LoadHistos()
        self.normalizeAndRatio()
        self.plotHistos()


    def LoadHistos(self):
        print ""
        print "Loading histograms from", self.histoPath
        print "year:", self.year
        print "match:", self.match
        print "steps:", self.histFolders
        print "histogram:", self.histoName

        for channel in ["electron", "muon", "invisible"]:
            for sample in self.samples:
                for histFolder in self.histFolders:
                    # skip the samples for invisiblechannel in lepton channel
                    if "inv" in sample and not "invisible" in channel: continue
                    if not "inv" in sample and "invisible" in channel and "Zprime" in sample: continue
                    suffix = ("_merge" if "DY" in sample else "")
                    match_rootFileName = (self.match if self.match !="incl" else "")
                    file_ = ROOT.TFile(self.histoPath+channel+"channel/nominal/"+""+self.PrefixrootFile+"MC."+sample+"_"+self.year+"_noTree"+suffix+".root")
                    # print "Loading file", self.histoPath+channel+"channel/nominal/"+""+self.PrefixrootFile+"MC."+sample+"_"+self.year+"_noTree"+suffix+".root"
                    self.histos[channel+"_"+sample.replace("_inv","")+"_"+histFolder+"_"+self.year+"_"+self.match] = file_.Get("nTopJet"+match_rootFileName+"_"+histFolder+"/"+self.histoName+"_jet").Clone(channel+"_"+sample+"_"+histFolder+"_"+self.year+"_"+self.match)
                    self.histos[channel+"_"+sample.replace("_inv","")+"_"+histFolder+"_"+self.year+"_"+self.match].SetDirectory(0)
                    file_.Close()
        print "Loaded files."

        # Proof that the files are loaded:
        # for channel in ["electron", "muon", "invisible"]:
        #     for sample in self.samples:
        #         for histFolder in self.histFolders:
        #             if "inv" in sample and not "invisible" in channel: continue
        #             if not "inv" in sample and "invisible" in channel and "Zprime" in sample: continue
        #             print "Value of bin 0 of histogram", channel+"_"+sample+"_"+histFolder+"_"+self.year+"_"+self.match, "is", self.histos[channel+"_"+sample+"_"+histFolder+"_"+self.year+"_"+self.match].GetBinContent(0)


    def normalizeAndRatio(self):

    # Normalize area to 1
        for channel in ["electron", "muon", "invisible"]:
            for sample in self.samples:
                if "inv" in sample: continue # always skip the inv samples
                for histFolder in self.histFolders:
                    # print "Normalizing ", sample, histFolder, "Channel:",channel

                    # print "Value of bin 0 of histogram", channel+"_"+sample+"_"+histFolder+"_"+self.year+"_"+self.match, "is", self.histos[channel+"_"+sample+"_"+histFolder+"_"+self.year+"_"+self.match].GetBinContent(0)

                    # Compute the integral over the whole range
                    area = self.histos[channel+"_"+sample+"_"+histFolder+"_"+self.year+"_"+self.match].Integral()
                    if (area == 0): print "Area is 0 for ", channel+"_"+sample+"_"+histFolder+"_"+self.year+"_"+self.match, "skipping normalization"
                    if (area == 0): continue

                    # Scale according to documentation:
                    # "One can scale an histogram such that the bins integral is equal to the normalization parameter via TH1::Scale(Double_t norm), where norm is the desired normalization divided by the integral of the histogram."
                    self.histos[channel+"_"+sample+"_"+histFolder+"_"+self.year+"_"+self.match].Scale(1/area)
                    # print "Value of bin 0 of histogram", channel+"_"+sample+"_"+histFolder+"_"+self.year+"_"+self.match, "is", self.histos[channel+"_"+sample+"_"+histFolder+"_"+self.year+"_"+self.match].GetBinContent(0)

                    # TODO check next line
                    # self.histos["bkg_prefit"].Scale(self.histos["bkg_prefit"].GetBinWidth(1))

        # Create the ratio relative to muon histogram
        for channel in ["electron", "muon", "invisible"]:
            for sample in self.samples:
                if "inv" in sample: continue
                for histFolder in self.histFolders:
                    histTitle = channel+"_"+sample+"_"+histFolder+"_"+self.year+"_"+self.match

                    self.histos[histTitle+"_ratio"] = self.histos[histTitle].Clone(histTitle+"_ratio")

                    # divide by the muon histogram: replace then channel with "muon" and delete the "_inv"
                    self.histos[histTitle+"_ratio"].Divide(self.histos[histTitle.replace(channel, "muon").replace("_inv","")])
        print "Normalized histos and calculated ratios."


    def plotHistos(self):

        # Create the canvas, skip for invisiblechannel
        for sample in self.samples:
            if "inv" in sample: continue
            for histFolder in self.histFolders:
                TDR.lumi_13TeV  = str(round(float(self.lumi_map[self.year]["lumi_fb"]),1))+" fb^{-1}"
                self.canv_ratio[sample+"_"+self.year+"_"+histFolder+self.match] = tdrCanvas("canv_ratio_"+sample+"_"+self.year+"_"+histFolder+"_"+self.match, 0, 1, -1, 10, self.histoName, "ratio")
                self.legend[sample+"_"+self.year+"_"+self.histoName+"_"+histFolder] = tdrLeg(0.30,0.80,0.89,0.89, 0.030, 42, ROOT.kBlack)

        # Loop over all samples (skipping the invisiblechannel samples), create
        # the plots and store them
        for sample in self.samples:
            if "inv" in sample: continue

            for histFolder in self.histFolders:

                for channel in ["electron", "muon", "invisible"]:
                    histTitle = channel+"_"+sample+"_"+histFolder+"_"+self.year+"_"+self.match

                    # Change canvas
                    self.canv_ratio[sample+"_"+self.year+"_"+histFolder+self.match].cd()

                    colorHistogram = colors[channel]

                    tdrDraw(self.histos[histTitle+"_ratio"], "hist", ROOT.kFullDotLarge, colorHistogram, colorHistogram, colorHistogram, 0, colorHistogram)

                    self.legend[sample+"_"+self.year+"_"+self.histoName+"_"+histFolder].AddEntry(self.histos[histTitle+"_ratio"], "ratio " + histTitle,"l")

                self.canv_ratio[sample+"_"+self.year+"_"+histFolder+self.match].cd()
                self.legend[sample+"_"+self.year+"_"+self.histoName+"_"+histFolder].Draw()
                self.canv_ratio[sample+"_"+self.year+"_"+histFolder+self.match].SaveAs(self.outdir+"ratio_"+self.histoName+"_"+sample+"_"+histFolder+"_"+self.year+"_"+self.match+".pdf")
        print "Plotted and stored histos."


def main():
    # Settings - also change them in the definition of CompareBTagVariables
    # years = ["2018"]
    years = ["2016", "2017", "2018"]
    histoNames = ["btag_DeepBoosted_ZHccvsQCD"]
    # histoNames = ["btag_DeepBoosted_HccvsQCD"]
    # histoNames = ["btag_MassDecorrelatedDeepBoosted_ZHccvsQCD"]
    # histoNames = ["btag_MassDecorrelatedDeepBoosted_HccvsQCD"]
    # histoNames = ["btag_MassDecorrelatedDeepBoosted_HccvsQCD", "btag_MassDecorrelatedDeepBoosted_ZHccvsQCD", "btag_DeepBoosted_HccvsQCD", "btag_DeepBoosted_ZHccvsQCD"]
    matches = ["", "Hbb", "Hcc", "Hqq", "HWW"]

    for year in years:
        for histoName in histoNames:
            for match in matches:
                CompareBTagVariables(year=year, match=match, histoName=histoName).RunAll()

if __name__ == '__main__':
    main()
