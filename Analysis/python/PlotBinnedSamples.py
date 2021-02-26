from Utils import *

import tdrstyle_all as TDR
TDR.writeExtraText = True
TDR.extraText  = "Simulation"
TDR.extraText2 = "Work in progress"

'''
Script to plot binned samples, such as the DY sample, which is binned in HT.

'''

class PlotBinnedSamples(VariablesBase):
    def __init__(self,year, steps):
        VariablesBase.__init__(self)
        self.year = year
        self.steps = steps
        TDR.lumi_13TeV  = str(round(float(self.lumi_map[self.year]["lumi_fb"]),1))+" fb^{-1}"
        self.module = module
        self.histFolders = histFolders
        if Channels!="": self.Channels = Channels
        self.histos = {}
        self.color  = {0:       ROOT.kBlack,
                       1:       ROOT.kRed+1,
                       2:       ROOT.kRed,
                       3:       ROOT.kGreen+2,
                       4:       ROOT.kGreen+1,
                       5:       ROOT.kBlue+1,
                       6:       ROOT.kAzure+10,
                       7:       ROOT.kViolet-3,
                       8:       ROOT.kMagenta+1,
                       }
        # self.HistType = {"Preselection": "nTopJet", "Selection": "ZprimeCandidate"}
        # self.HistName = {"Preselection": "pt_jet",  "Selection": "Zprime_mass_rebin30"}
        self.outdir = self.Path_ANALYSIS+"Analysis/OtherPlots/BinnedSamples/"
        os.system("mkdir -p "+self.outdir)
        print("Files will be written to : " + self.outdir)

    def PlotHistos(self):

        self.samples = ["MC_DY", "MC_DY_inv_HT200to400", "MC_DY_inv_HT400to600", "MC_DY_inv_HT600to800", "MC_DY_inv_HT800to1200", "MC_DY_inv_HT1200to2500", "MC_DY_inv_HT2500toInf"]
        # self.samples = ["MC_DY", "MC_DY_HT200to400", "MC_DY_HT400to600", "MC_DY_HT600to800", "MC_DY_HT800to1200", "MC_DY_HT1200to2500", "MC_DY_HT2500toInf"]

        extraText = ""

        ymin = 1
        ymax = 10**8

        commonpath = self.Path_STORAGE+self.year+"/"+"Preselection"+"/"+"Puppi/invisiblechannel/"

        # sampleName = "pt_Z"
        sampleName = "HT"

        print "Will load from " + commonpath+"nominal"+"/"

        numberSample = 0
        for step in steps:
            canv = tdrCanvas("canv_"+"MC_DY_inv+"+step+"extraText", 0, 4000, ymin, ymax , sampleName, "Events")
            # canv = tdrCanvas("canv_"+"MC_DY_muon+"+step+"extraText", 0, 2000, ymin, ymax , "p_T Z", "Events")
            leg = tdrLeg(0.60,0.40,0.89,0.89, 0.030, 42, ROOT.kBlack);
            canv.SetLogy(True)
            self.histos={}
            for sample in (self.samples):

                mode = "MC" if "MC" in sample else "DATA"
                suffix = ("_noTree.root" if "HT" in sample else "_noTree_merge.root")
                filename = commonpath+"nominal"+"/"+self.PrefixrootFile+mode+"."+sample+"_"+year+suffix

                print "Loading: " + sample
                file_ = ROOT.TFile(filename)

                h_ = file_.Get("gen_"+step+"/"+sampleName)

                h_.SetDirectory(0) # .setdefault(sample, {})
                self.histos.setdefault("MC_DY_inv",{}).setdefault("gen_"+step, {})["MC_DY"] = h_
                tdrDraw(self.histos["MC_DY_inv"]["gen_"+step]["MC_DY"], "hist", ROOT.kDot, self.color[numberSample], 1, self.color[numberSample], 0, self.color[numberSample])
                leg.AddEntry(self.histos["MC_DY_inv"]["gen_"+step]["MC_DY"], sample, "l")
                numberSample += 1
                if numberSample == 9: numberSample = 0
                file_.Close()
            leg.Draw("same")

            canv.SaveAs(self.outdir+"MC_DY_inv_gen_"+step+"_"+sampleName+"_"+extraText+".pdf")
            print "Saved in " + self.outdir+"MC_DY_inv_gen_"+step+"_"+sampleName+"_"+extraText+".pdf"

            # canv.SaveAs(self.outdir+"MC_DY_muon_gen_"+step+"_"+sampleName+"_"+extraText+".pdf")
            # print "Saved in " + self.outdir+"MC_DY_muon_gen_"+step+"_"+sampleName+"_"+extraText+".pdf"


if __name__ == '__main__':
    args = parse_arguments()
    years       = args.years
    histFolders = args.histFolders
    Channels    = args.Channels

    studies = "nominal"
    module = "PlotBinnedSamples"
    years = ["2016"]
    # years = ["RunII"]

    steps = ["nocuts", "weights", "JetDiLeptonPhiAngular"]

    print "Years: " + str(years)
    print "Steps: " + str(steps)

    for year in years:
        PlotSyst = PlotBinnedSamples(year=year, steps=steps)
        PlotSyst.PlotHistos()
