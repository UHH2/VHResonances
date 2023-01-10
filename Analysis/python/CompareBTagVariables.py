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
            "muon"  :      ROOT.kRed+1,
            "electron"  :  ROOT.kOrange+1,
            }

class CompareBTagVariables(ModuleRunnerBase):
    def __init__(self, match="incl", histoName="btag_MassDecorrelatedDeepBoosted_ZHccvsQCD", debug=False):
        VariablesBase.__init__(self)
        self.years          = ["2016", "2017", "2018"] # include all years
        self.match          = (match if match!="" else "incl")
        # self.histFolders    = ["nocuts", "weights", "Trigger", "HEM", "cleaned", "Veto", "NLeptonSel", "NBoostedJet", "METCut", "DeltaRDiLepton", "JetDiLeptonPhiAngular", "QCDRejection"]
        self.histFolders    = ["nocuts", "weights", "QCDRejection"]
        self.histoPath      = self.Path_STORAGE+"YEAR/JetFlavour/Puppi/"
        self.histoName      = histoName
        self.histos         = {}
        self.canv_ratio     = {}
        self.legend         = {}
        self.outdir = self.Path_ANALYSIS+"Analysis/OtherPlots/CompareBTagVariables/"
        self.samples = ["MC_DY", "MC_ZprimeToZH_M2000", "MC_ZprimeToZH_inv_M2000"]
        self.debug          = debug

        os.system("mkdir -p "+self.outdir)
        if self.debug: self.samples, match, histoName

    def doAll(self):

        for year, channel, sample in list(itertools.product(["2016", "2017", "2018"], ["muon", "electron", "invisible"], self.samples)):
            prefixSample = "DATA" if "DATA" in sample else "MC"

            sample_names_data ={"electron": "DATA_SingleElectron", "muon": "DATA_SingleMuon", "invisible": "DATA_MET"}
            if "DATA" in sample: sample = sample_names_data[channel]

            # skip the samples for invisiblechannel in lepton channel
            if "inv" in sample and not "invisible" in channel: continue
            if not "inv" in sample and "invisible" in channel and "Zprime" in sample: continue
            suffix = ("_merge" if ("DY" in sample or "DATA" in sample or "WJets" in sample) else "")
            match_rootFileName = (self.match if self.match !="incl" else "")

            filepath = self.histoPath.replace("YEAR", year)+channel+"channel/nominal/"+""+self.PrefixrootFile+prefixSample+"."+sample+"_"+year+"_noTree"+suffix+".root"

            # WJets is only available for the invisible channel, do other denominators: Use the DY sample, either invisiblechannel or electronchannel as a denominator.
            # Do this by changing the filename of the file to load, in order to maintain all other strings (sample, channel, ...)
            if "WJets" in sample and "electron" in channel: filepath = self.histoPath.replace("YEAR", year)+channel+"channel/nominal/"+""+self.PrefixrootFile+prefixSample+".MC_DY_"+year+"_noTree_merge.root"
            if "WJets" in sample and "muon" in channel: filepath = self.histoPath.replace("YEAR", year)+"invisiblechannel/nominal/"+""+self.PrefixrootFile+prefixSample+".MC_DY_"+year+"_noTree_merge.root"

            file_ = ROOT.TFile(filepath)
            if self.debug: print "Loading file", filepath

            if "DATA" in sample: sample = "DATA"

            # normalise
            for histFolder in self.histFolders:
                self.histos[channel+sample.replace("_inv", "")+histFolder+year+self.match] = file_.Get("nTopJet"+match_rootFileName+"_"+histFolder+"/"+self.histoName+"_jet").Clone(channel+sample+histFolder+year+self.match)
                self.histos[channel+sample.replace("_inv", "")+histFolder+year+self.match].SetDirectory(0)

                if "inv" in sample: continue #  skip the inv samples (e.g. the signal)

                histoName = channel+sample+histFolder+year+self.match
                area = self.histos[histoName].Integral()
                if (area == 0): print "Area is 0 for ", histoName, "skipping normalization"
                if (area == 0): continue

                self.histos[histoName].Scale(1/area)
            file_.Close()

        # Create the ratio relative to muon histogram
        for year, channel, sample in list(itertools.product(["2016", "2017", "2018"], ["electron", "muon", "invisible"], self.samples)):
            if "WJets" in sample: continue # treat this special case later
            if "inv" in sample: continue
            for histFolder in self.histFolders:
                histTitle = channel+sample+histFolder+year+self.match
                self.histos[histTitle+"ratio"] = self.histos[histTitle].Clone(histTitle+"ratio")
                # divide by the muon histogram
                self.histos[histTitle+"ratio"].Divide(self.histos[histTitle.replace(channel, "muon").replace("inv","")])

        # WJets with different denominators
        channel = "invisible"
        for year in ["2016", "2017", "2018"]:
            if not "MC_Jets" in self.samples: continue
            sample = "MC_WJets"

            for histFolder in self.histFolders:
                histTitle = channel+sample+histFolder+year+self.match
                self.histos[histTitle+"D1ratio"] = self.histos[histTitle].Clone(histTitle+"D1ratio")
                self.histos[histTitle+"D2ratio"] = self.histos[histTitle].Clone(histTitle+"D2ratio")

                # Use D1 (DY of electronchannel) as a denominator and WJets of invisiblechannel as numerator, D1 is stored in WJets for the electronchannel.
                self.histos[histTitle+"D1ratio"].Divide(self.histos[histTitle.replace(channel, "electron")])

                # Use D2 (DY of invisiblechannel) as a denominator and WJets of invisiblechannel as numerator, D2 is stored in WJets for the muonchannel.
                self.histos[histTitle+"D2ratio"].Divide(self.histos[histTitle.replace(channel, "muon")])

        # Create the canvas, skip for samples of invisiblechannel (like signal sample)
        for year, sample in list(itertools.product(self.years, self.samples)):
            if "inv" in sample: continue
            for histFolder in self.histFolders:
                TDR.cms_lumi = self.lumi_map[year]['lumiPlot']+' fb^{-1}'
                self.canv_ratio[sample+year+histFolder+self.match] = tdrCanvas("canvratio"+sample+year+histFolder+self.match, 0, 1, -1, 10, self.histoName, "ratio")
                self.legend[sample+year+self.histoName+histFolder] = tdrLeg(0.30,0.80,0.89,0.89, 0.030, 42, ROOT.kBlack)

        # Loop over all samples (skipping the samples with "inv" in the name, such as the signal samples), create the plots and store them
        for year, sample, histFolder in list(itertools.product(self.years, self.samples, self.histFolders)):
            if "inv" in sample: continue

            for channel in ["electron", "muon", "invisible"]:
                histTitle = channel+sample+histFolder+year+self.match
                colorHistogram = colors[channel]

                # Plots with different nominators for WJets
                if "WJets" in sample:
                    if "invisible" in channel: continue
                    denominatorSuffix = ("D1" if "electron" in channel else "D2")
                    channel = "invisible"
                    histTitle = channel+sample+histFolder+year+self.match+denominatorSuffix

                # Change canvas
                self.canv_ratio[sample+year+histFolder+self.match].cd()

                tdrDraw(self.histos[histTitle+"ratio"], "p", ROOT.kFullDotLarge, colorHistogram, colorHistogram, colorHistogram, 0, colorHistogram)
                tdrDraw(self.histos[histTitle+"ratio"], "hist", ROOT.kFullDotLarge, colorHistogram, colorHistogram, colorHistogram, 0, colorHistogram)

                self.legend[sample+year+self.histoName+histFolder].AddEntry(self.histos[histTitle+"ratio"], "ratio " + histTitle,"l")

            self.canv_ratio[sample+year+histFolder+self.match].cd()
            self.legend[sample+year+self.histoName+histFolder].Draw()
            self.canv_ratio[sample+year+histFolder+self.match].SaveAs(self.outdir+"ratio_"+self.histoName.replace("MassDecorrelated", "MD").replace("DeepBoosted", "DB")+"_"+sample+"_"+histFolder+"_"+year+"_"+self.match+".pdf")

        # Take 2017 of a channel as denominator and 2016 or 2018 as nominator of same channel.
        # histograms are normalised
        TDR.cms_lumi = ""
        for channel, sample, histFolder in list(itertools.product(["electron","muon", "invisible"], self.samples, self.histFolders)):
            if "inv" in sample: continue
            histo_den = self.histos[channel+sample+histFolder+"2017"+self.match]
            for year in ["2016", "2018"]:
                canvas = tdrCanvas("canvas", 0, 1, 0.001, 10, self.histoName, "ratio")
                canvas.SetLogy()
                legend = tdrLeg(0.30,0.80,0.89,0.89, 0.030, 42, ROOT.kBlack)
                histo_nom = self.histos[channel+sample+histFolder+year+self.match]
                histo_ratio = histo_nom.Clone("ratio")
                histo_ratio.Divide(histo_den)
                colorHistogram = ROOT.kOrange; tdrDraw(histo_nom  , "p", ROOT.kFullDotLarge, colorHistogram, colorHistogram, colorHistogram, 0, colorHistogram)
                colorHistogram = ROOT.kRed;    tdrDraw(histo_den  , "p", ROOT.kFullDotLarge, colorHistogram, colorHistogram, colorHistogram, 0, colorHistogram)
                colorHistogram = ROOT.kBlack;  tdrDraw(histo_ratio, "p", ROOT.kFullDotLarge, colorHistogram, colorHistogram, colorHistogram, 0, colorHistogram)
                legend.AddEntry(histo_nom, channel+ " "+sample+ " "+ histFolder + " " + year,"l")
                legend.AddEntry(histo_den, channel+ " "+sample+ " "+ histFolder + " 2017","l")
                legend.AddEntry(histo_ratio, channel+ " "+sample+ " "+ histFolder + " ratio "+year+"/2017","l")
                canvas.SaveAs(self.outdir+"ratio_"+channel+"_"+self.histoName.replace("MassDecorrelated", "MD").replace("DeepBoosted", "DB")+"_"+sample+"_"+histFolder+"_ratio_"+year+"_to_2017.pdf")

def main():
    print("Using silent ROOT mode.")
    ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = " + str(ROOT.kError) + ";")

    histoNames = ["btag_MassDecorrelatedDeepBoosted_ZHccvsQCD"]
    histoNames = ["btag_MassDecorrelatedDeepBoosted_HccvsQCD", "btag_MassDecorrelatedDeepBoosted_ZHccvsQCD", "btag_DeepBoosted_HccvsQCD", "btag_DeepBoosted_ZHccvsQCD"]

    # matches = ["", "Hbb", "Hcc", "Hqq", "HWW"]
    matches = [""]

    debug = False
    debug = True


    for histoName in histoNames:
        for match in matches:
            CompareBTagVariables(match=match, histoName=histoName, debug=debug).doAll()

if __name__ == '__main__':
    main()
