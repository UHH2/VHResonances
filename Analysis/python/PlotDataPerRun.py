from Utils import *

import tdrstyle_all as TDR
TDR.writeExtraText = True
TDR.extraText  = "Simulation"
TDR.extraText2 = "Work in progress"

'''
Module to plot binned variable

'''

class PlotDataPerRun(VariablesBase):
    def __init__(self,year, steps, variables):
        VariablesBase.__init__(self)
        self.year = year
        self.steps = steps
        self.variables = variables
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
        self.outdir = self.Path_ANALYSIS+"Analysis/OtherPlots/PlotDataPerRun/"
        os.system("mkdir -p "+self.outdir)
        print("Files will be written to : " + self.outdir)

    def PlotHistos(self):

        self.samples = ["DATA_MET_RunA", "DATA_MET_RunB", "DATA_MET_RunC", "DATA_MET_RunD"]

        extraText = ""

        ymin = 1
        ymax = 10**10

        for variable in variables:
            for step in steps:
                canv = tdrCanvas("canv_"+"DATA+"+step+"extraText", -4, 4, ymin, ymax , variable+"_jet", "Events")
                leg = tdrLeg(0.70,0.9,0.89,0.6, 0.030, 42, ROOT.kBlack);
                canv.SetLogy(True)
                self.histos={}
                numberSample = 0
                for sample in (self.samples):

                    commonpath = self.Path_STORAGE+self.year+"/"+"Preselection"+"/"+"Puppi/invisiblechannel/"
                    mode = "MC" if "MC" in sample else "DATA"
                    suffix = "_noTree.root"
                    filename = commonpath+"nominal"+"/"+self.PrefixrootFile+mode+"."+sample+"_"+year+suffix

                    print "Loading: " + sample
                    file_ = ROOT.TFile(filename)

                    h_ = file_.Get("nTopJet_"+step+"/"+variable+"_jet")

                    h_.SetDirectory(0)
                    self.histos.setdefault("Data",{}).setdefault("nTopJet_"+step, {})[sample]= h_
                    tdrDraw(self.histos["Data"]["nTopJet_"+step][sample], "hist", ROOT.kDot, self.color[numberSample], 1, self.color[numberSample], 0, self.color[numberSample])
                    leg.AddEntry(self.histos["Data"]["nTopJet_"+step][sample], sample, "l")
                    numberSample += 1
                    if numberSample == 9: numberSample = 0
                    file_.Close()
                leg.Draw("same")

                canv.SaveAs(self.outdir+"Data_nTopJet_"+step+"_"+ variable +extraText+".pdf")
                print "Saved in " + self.outdir+"Data_nTopJet_"+step+"_"+variable+extraText+".pdf"


if __name__ == '__main__':
    args = parse_arguments()
    years       = args.years
    histFolders = args.histFolders
    Channels    = args.Channels

    studies = "nominal"
    module = "PlotDataPerRun"
    years = ["2018"]
    # years = ["RunII"]

    steps = ["weights", "JetDiLeptonPhiAngular"]
    variables = ["eta", "phi"]

    print "Years: " + str(years)
    print "Steps: " + str(steps)
    print "Variables: " + str(variables)

    for year in years:
        PlotSyst = PlotDataPerRun(year=year, steps=steps, variables=variables)
        PlotSyst.PlotHistos()
