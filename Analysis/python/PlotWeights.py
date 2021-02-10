from Utils import *

import tdrstyle_all as TDR
TDR.writeExtraText = True
TDR.extraText  = "Simulation"
TDR.extraText2 = "Work in progress"

'''
Script to plot weights (Preselection or TestModule output)
Change the steps, baseHistoStepName and the year as needed.
A ratio can be drawn, but only for a TestModule with correct input files,
this is not possible for the Preselection.
'''

class PlotWeights(VariablesBase):
    def __init__(self,year, steps, variables, samples, module, subfolder=""):
        VariablesBase.__init__(self)
        self.year = year
        self.steps = steps
        self.variables = variables
        self.samples = samples
        self.module = module
        TDR.lumi_13TeV  = str(round(float(self.lumi_map[self.year]["lumi_fb"]),1))+" fb^{-1}"
        self.module = module
        self.histFolders = histFolders
        if Channels!="": self.Channels = Channels
        self.histos = {}
        self.color  = {0:       ROOT.kBlack,
                       1:       ROOT.kGreen+2,
                       2:       ROOT.kBlue+1,
                       3:       ROOT.kRed+1,
                       4:       ROOT.kAzure+10,
                       5:       ROOT.kRed,
                       6:       ROOT.kGreen+1,
                       7:       ROOT.kViolet-3,
                       8:       ROOT.kMagenta+1,
                       }
        self.subfolder = subfolder
        self.outdir = os.path.join(self.Path_ANALYSIS, "Analysis/OtherPlots/PlotWeights/", self.year, self.module, self.subfolder)
        os.system("mkdir -p "+self.outdir)
        print("Files will be written to : " + self.outdir)

        self.typeAndVariable  =  {"pt_jet":     "nTopJet",
                                  "Weights":    "event",
                                  "WeightsLogBins":    "event"
                                    }

        self.xMin  =  {"pt_jet":           0,
                        "Weights":         0,
                        "WeightsLogBins": -1
                        }

        self.xMax  =  {"pt_jet":          5*10**3,
                       "Weights":               2,
                       "WeightsLogBins":      110
                        }

        self.yMin  =  {"pt_jet":          -20000,
                        "Weights":        -10**4,
                        "WeightsLogBins":      1
                        }

        self.yMax  =  {"pt_jet":               50000,
                       "Weights":          1.2*10**5,
                       "WeightsLogBins":  1.4* 10**5
                        }

        self.doLog  =  {"pt_jet":         False,
                        "Weights":        False,
                        "WeightsLogBins":  True
                        }

        self.ratioMinY  =  {"pt_jet":   -2,
                            "Weights":  -5
                        }

        self.ratioMaxY  = {"pt_jet":     5,
                           "Weights":   50
                        }
    def PlotHistos(self):

        for variable in self.variables:
            for step in self.steps:
                canv = tdrCanvas("canv_"+self.year+"_"+step+"_"+self.subfolder+"_"+variable, self.xMin[variable], self.xMax[variable], self.yMin[variable], self.yMax[variable], variable, "Events")
                # no extra entries leg = tdrLeg(0.85,0.9,0.95, 0.9 - len(self.samples) * 0.04, 0.030, 42, ROOT.kBlack);
                leg = tdrLeg(0.65,0.9,0.95, 0.9 - len(self.samples) * 0.04, 0.030, 42, ROOT.kBlack);
                canv.SetLogy(self.doLog[variable])
                self.histos = {}
                self.histos.setdefault(self.module, {})
                numberSample = 0
                for sample in (self.samples):

                    histFolder = self.typeAndVariable[variable] + "_" + step

                    commonpath = self.Path_STORAGE+self.year+"/"+self.module+"/"+"Puppi/invisiblechannel/"

                    mode = "MC" if "MC" in sample else "DATA"
                    suffix = "_noTree.root"
                    filename = commonpath+"nominal"+"/"+self.PrefixrootFile+mode+"."+sample+"_"+year+suffix

                    # print "Loading: " + sample + " step: " + step
                    file_ = ROOT.TFile(filename)

                    h_ = file_.Get(histFolder+"/"+variable)

                    h_.SetDirectory(0)
                    self.histos[self.module].setdefault(histFolder, {})[sample]= h_
                    tdrDraw(self.histos[self.module][histFolder][sample], "hist", ROOT.kDot, self.color[numberSample], 1, self.color[numberSample], 0, self.color[numberSample])
                    # leg.AddEntry(self.histos[self.module][histFolder][sample], sample.replace("MC_ZprimeToZH_inv_", ""), "l")
                    hist = self.histos[self.module][histFolder][sample]
                    leg.AddEntry (hist, sample.replace("MC_ZprimeToZH_inv_", "")+ " u: "+str(round(hist.GetBinContent(0),2))+" o: "+str(round(hist.GetBinContent(hist.GetNbinsX()+1),2)), "lp")
                    numberSample += 1
                    if numberSample == 9: numberSample = 0
                    file_.Close()
                leg.Draw("same")

                canv.SaveAs(os.path.join(self.outdir,histFolder+"_"+ variable +".pdf"))
                # print "Saved in " +os.path.join(self.outdir,histFolder+"_"+variable+".pdf")

    def PlotRatioHistos(self):

        for variable in self.variables:
            self.histos = {}

            self.histos.setdefault(self.module, {})

            for sample in (self.samples):
                numberSample = 0
                for step in self.steps:

                    histFolder = self.typeAndVariable[variable] + "_" + step

                    commonpath = self.Path_STORAGE+self.year+"/"+self.module+"/"+"Puppi/invisiblechannel/"
                    mode = "MC" if "MC" in sample else "DATA"
                    suffix = "_noTree.root"
                    filename = commonpath+"nominal"+"/"+self.PrefixrootFile+mode+"."+sample+"_"+year+suffix

                    file_ = ROOT.TFile(filename)

                    h_ = file_.Get(histFolder+"/"+variable)

                    h_.SetDirectory(0)
                    # Store the histogram
                    self.histos[self.module].setdefault(histFolder, {})[sample]= h_

                    file_.Close()


            for sample in (self.samples):
                numberSample = 0

                baseHistoStepName = "nocuts"
                baseHistoStepName = "weights"
                baseHisto = self.histos[self.module][self.typeAndVariable[variable] + "_" + baseHistoStepName][sample].Clone()
                # print "The baseHisto is called " + self.typeAndVariable[variable] + "_" + baseHistoStepName + " with sample " + sample
                # print "baseHisto should not change! Initial value in bin 10: "  + str(baseHisto.GetBinContent(10))

                canv = tdrCanvas("canv_"+sample+"ratio", self.xMin[variable], self.xMax[variable], self.ratioMinY[variable], 7000 if "80x00" in sample else self.ratioMaxY[variable], variable, "ratio")
                leg = tdrLeg(0.6,0.9,0.95, 0.9 - len(self.samples) * 0.04, 0.030, 42, ROOT.kBlack);
                canv.SetLogy(self.doLog[variable])
                for step in self.steps:
                    histFolder = self.typeAndVariable[variable] + "_" + step
                    if self.histos[self.module][histFolder][sample] == baseHisto:
                        print "Not dividing histogram by itself."
                        raise Exception("The histogram is divided by itself. This means it was not copied as a separate clone.")
                    else:
                        self.histos[self.module][histFolder][sample].Divide(baseHisto)
                        tdrDraw(self.histos[self.module][histFolder][sample], "hist", ROOT.kDot, self.color[numberSample], 1, self.color[numberSample], 0, self.color[numberSample])
                        hist = self.histos[self.module][histFolder][sample]
                        leg.AddEntry(hist, sample.replace("MC_ZprimeToZH_inv_", "") + ": " + step + "/" + baseHistoStepName , "lp")
                        numberSample += 1
                        if numberSample == 9: numberSample = 0

                leg.Draw("same")
                outFile = os.path.join(self.outdir,self.typeAndVariable[variable]+"_"+ variable +"_"+sample.replace("MC_ZprimeToZH_inv_", "")+ "_ratio.pdf")
                canv.SaveAs(outFile)
                canv.Close()
                # print "baseHisto should not change!  Final value in bin 10: "  + str(baseHisto.GetBinContent(10))
                # print "Saved in " + outFile

if __name__ == '__main__':
    args = parse_arguments()
    years       = args.years
    histFolders = args.histFolders
    Channels    = args.Channels

    studies = "nominal"
    module = "PlotWeights"
    years = ["2016", "2017", "2018"]
    # years = ["2017"]

    steps = ["weights", "nocuts"]

    variables = ["pt_jet", "Weights", "WeightsLogBins"]

    allSamples = [
    "MC_ZprimeToZH_inv_M600",  "MC_ZprimeToZH_inv_M800",  "MC_ZprimeToZH_inv_M1000",
    "MC_ZprimeToZH_inv_M1200", "MC_ZprimeToZH_inv_M1400", "MC_ZprimeToZH_inv_M1600",
    "MC_ZprimeToZH_inv_M1800", "MC_ZprimeToZH_inv_M2000", "MC_ZprimeToZH_inv_M2500",
    "MC_ZprimeToZH_inv_M3000", "MC_ZprimeToZH_inv_M3500", "MC_ZprimeToZH_inv_M4000",
    "MC_ZprimeToZH_inv_M4500", "MC_ZprimeToZH_inv_M5000", "MC_ZprimeToZH_inv_M5500",
    "MC_ZprimeToZH_inv_M6000", "MC_ZprimeToZH_inv_M7000", "MC_ZprimeToZH_inv_M8000"]
    allSubfolder = "allMassPoints"

    samples = ["MC_ZprimeToZH_inv_M600", "MC_ZprimeToZH_inv_M3000", "MC_ZprimeToZH_inv_M8000"]
    subfolder = "someMassPoints"

    ratioSamples = ["MC_ZprimeToZH_inv_M600", "MC_ZprimeToZH_inv_M1000", "MC_ZprimeToZH_inv_M2000", "MC_ZprimeToZH_inv_M3000", "MC_ZprimeToZH_inv_M8000"]
    ratioSubfolder = "someMassPoints-ratio"
    ratioVariables = ["Weights", "pt_jet"]

    modules = ["Test", "Preselection"]
    modules = ["Test"]
    print "Years: " + str(years)
    print "Steps: " + str(steps)
    print "Variables: " + str(variables)
    print "Module: " + str(module)

    print("Using silent ROOT mode.")
    ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = " + str(ROOT.kError) + ";")

    for year in years:
        for module in modules:
            Plotter = PlotWeights(year=year, steps=steps, variables=variables, samples=allSamples, module=module, subfolder=allSubfolder)
            Plotter.PlotHistos()

            Plotter = PlotWeights(year=year, steps=steps, variables=variables, samples=samples, module=module, subfolder=subfolder)
            Plotter.PlotHistos()

        module="Test" # only here is the ratio possible.
        Plotter = PlotWeights(year=year, steps=steps, variables=ratioVariables, samples=ratioSamples, module=module, subfolder=ratioSubfolder)
        # Plotter.PlotHistos()
        Plotter.PlotRatioHistos()

