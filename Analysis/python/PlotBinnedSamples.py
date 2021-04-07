from Utils import *

import tdrstyle_all as TDR
TDR.writeExtraText = True
TDR.extraText  = "Simulation"
TDR.extraText2 = "Work in progress"

'''
Script to plot binned samples, such as the DY sample, which is binned in HT.

'''

class PlotBinnedSamples(VariablesBase):
    def __init__(self,year, steps, samples, sampleName, channel, debug=False):
        VariablesBase.__init__(self)
        self.year = year
        self.steps = steps
        TDR.lumi_13TeV  = str(round(float(self.lumi_map[self.year]["lumi_fb"]),1))+" fb^{-1}"
        self.module = module
        self.histFolders = histFolders
        self.samples = samples
        self.sampleName = sampleName
        self.channel = channel
        self.debug = debug
        self.selectionModule = "JetFlavour"

        if self.debug: print "SelectionModule", self.selectionModule

        if not self.debug: ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = " + str(ROOT.kError) + ";")

        if Channels!="": self.Channels = Channels
        self.histos = {}
        self.color  = {0:       ROOT.kBlack,
                       1:       ROOT.kRed,
                       2:       ROOT.kGreen,
                       3:       ROOT.kBlue,
                       4:       ROOT.kYellow,
                       5:       ROOT.kRed   +2,
                       6:       ROOT.kGreen +2,
                       7:       ROOT.kBlue  +2,
                       8:       ROOT.kYellow+2,
                       }
        self.outdir = self.Path_ANALYSIS+"Analysis/OtherPlots/BinnedSamples/"
        os.system("mkdir -p "+self.outdir)
        if self.debug: print("Files will be written to : " + self.outdir)

    def PlotHistos(self):

        ymin = 1
        ymax = 10**8

        commonpath = self.Path_STORAGE+self.year+"/"+self.selectionModule+"/"+"Puppi/" + self.channel + "/"

        # sampleName = "pt_Z"
        sampleName = "HT"

        if self.debug: print "Will load from " + commonpath+"nominal"+"/"

        for step in steps:
            numberSample = 0
            canv = tdrCanvas("canv"+self.sampleName+step+self.year+self.channel, 0, 4000, ymin, ymax , sampleName + " | " + step, "Events")
            leg = tdrLeg(0.60,0.40,0.89,0.89, 0.030, 42, ROOT.kBlack);
            canv.SetLogy(True)
            self.histos={}
            for sample in (self.samples):

                mode = "MC" if "MC" in sample else "DATA"
                suffix = ("_noTree.root" if "HT" in sample else "_noTree_merge.root")
                filename = commonpath+"nominal"+"/"+self.PrefixrootFile+mode+"."+sample+"_"+year+suffix

                if self.debug: print "Loading: " + sample
                file_ = ROOT.TFile(filename)

                h_ = file_.Get("gen_"+step+"/"+sampleName)

                h_.SetDirectory(0) # .setdefault(sample, {})
                self.histos.setdefault(self.sampleName,{}).setdefault("gen_"+step, {})[self.sampleName] = h_
                tdrDraw(self.histos[self.sampleName]["gen_"+step][self.sampleName], "hist", ROOT.kDot, self.color[numberSample], 1,  self.color[numberSample], 0, self.color[numberSample])
                self.histos[self.sampleName]["gen_"+step][self.sampleName].SetLineWidth(3)
                leg.AddEntry(self.histos[self.sampleName]["gen_"+step][self.sampleName], sample, "l")
                numberSample += 1
                if numberSample == 9: numberSample = 0
                file_.Close()
            leg.Draw("same")

            filename = self.outdir+self.sampleName+"_"+sampleName+"_"+step+"_"+self.year+"_"+self.channel+".pdf"
            canv.SaveAs(filename)
            if self.debug: print "Saved in", filename


if __name__ == '__main__':
    args = parse_arguments()
    years       = args.years
    histFolders = args.histFolders
    Channels    = args.Channels

    studies = "nominal"
    module = "PlotBinnedSamples"
    years = ["2016", "2017", "2018"]
    # years = ["RunII"]

    steps = ["nocuts", "weights", "JetDiLeptonPhiAngular"]

    print "Years: " + str(years)
    print "Steps: " + str(steps)

    samples = ["MC_DY_inv_HT100to200", "MC_DY_inv_HT200to400", "MC_DY_inv_HT400to600", "MC_DY_inv_HT600to800", "MC_DY_inv_HT800to1200", "MC_DY_inv_HT1200to2500", "MC_DY_inv_HT2500toInf"]
    sampleName = "MC_DY"
    channel = "invisiblechannel"

    debug = False

    for year in years:
        PlotSyst = PlotBinnedSamples(year=year, steps=steps, samples=samples, sampleName=sampleName, channel=channel, debug=debug)
        PlotSyst.PlotHistos()
