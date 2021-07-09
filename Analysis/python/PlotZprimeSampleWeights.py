from Utils import *

import tdrstyle_all as TDR
import glob
import re
from datetime import datetime
import xml.etree.ElementTree as ET

TDR.writeExtraText = True
TDR.extraText  = "Simulation"
TDR.extraText2 = "Work in progress"

'''
Module to look at Zprime samples and their weights

'''

class PlotZprimeSampleWeights(VariablesBase):
    def __init__(self, year, beVerbose=False):
        VariablesBase.__init__(self)
        self.year = year
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
        self.outdir = self.Path_ANALYSIS+"Analysis/OtherPlots/PlotZprimeSampleWeights/"
        self.beVerbose = beVerbose
        os.system("mkdir -p "+self.outdir)
        print("Files will be written to : " + self.outdir)

    def AnalyzeYear(self):
        inputPath = {
        '2016': self.Path_ANALYSIS + '../common/UHH2-datasets/RunII_102X_v2/2016v3/ZprimeToZH/',
        '2017': self.Path_ANALYSIS + '../common/UHH2-datasets/RunII_102X_v2/2017/ZprimeToZH/',
        '2018': self.Path_ANALYSIS + '../common/UHH2-datasets/RunII_102X_v2/2018/ZprimeToZH/'}

        extraText = ""

        ymin = 1
        ymax = 10**8

        print "Path: " + inputPath[year]

        self.samples = glob.glob(inputPath[year] + "/*Hall*")#

        numberOfEventsYear = 0
        numberOfNtupleFilesYear = 0
        numberOfXMLfilesYear = 0
        numberOfEventsWithNegWeightYear = 0
        numberOfEventsWithZeroWeightYear = 0
        numberOfEventsWithWeightOneYear = 0
        sumOfNegWeightsYear = 0
        sumOfAllWeightsYear = 0

        # per year
        numberOfEvents = 0
        numberOfNtupleFiles = 0
        numberOfXMLfiles = 0
        numberOfEventsWithNegWeight = 0
        numberOfEventsWithZeroWeight = 0
        numberOfEventsWithWeightOne = 0
        sumOfNegWeights = 0
        sumOfAllWeights = 0

        # for progress
        currentXMLfile = 0

        numberSamples = len(self.samples)

        for sample in self.samples:

            # store them
            numberOfEventsYear += numberOfEvents
            numberOfNtupleFilesYear += numberOfNtupleFiles
            numberOfXMLfilesYear += numberOfXMLfiles
            numberOfEventsWithNegWeightYear += numberOfEventsWithNegWeight
            numberOfEventsWithWeightOneYear += numberOfEventsWithWeightOne
            numberOfEventsWithZeroWeightYear += numberOfEventsWithZeroWeight
            sumOfNegWeightsYear += sumOfNegWeights
            sumOfAllWeightsYear += sumOfAllWeights

            # have them per file
            numberOfEvents = 0
            numberOfNtupleFiles = 0
            numberOfXMLfiles = 0
            numberOfEventsWithNegWeight = 0
            numberOfEventsWithWeightOne = 0
            sumOfNegWeights = 0
            sumOfAllWeights = 0

            print "Opening XML file " + str(currentXMLfile) + " out of " + str(numberSamples)
            numberOfXMLfiles += 1
            currentXMLfile += 1

            canv = tdrCanvas("canv_"+"ZprimeToZH_"+sample, 0, 2000, ymin, ymax , "weights", "Events")
            leg = tdrLeg(0.60,0.40,0.89,0.89, 0.030, 42, ROOT.kBlack);
            canv.SetLogy(True)
            self.histos={}

            filename = os.path.join(inputPath[year], sample)

            file_list_ = open(filename)
            fl = file_list_.readlines()

            maxfiles = len(fl)
            currentfile = 0

            for currentLine in fl:
                currentfile+=1
                x = re.findall("FileName=\"(.*)\" Lumi", currentLine)
                if len(x)==0: continue
                # print x[0]

                # x contains all the ntuple file names
                file_ = ROOT.TFile(x[0])
                numberOfNtupleFiles += 1
                # print "    opened " + x[0] + " (File " + str(currentfile) + " out of " + str(maxfiles) + " files)"
                print "  ntuple file "  +str(currentfile) + " out of " + str(maxfiles) + " files"

                tree_ = file_.Get("AnalysisTree")
                for ev in tree_:
                    numberOfEvents += 1
                    if numberOfEvents % 1000 == 0:
                        print "    processed " + str(numberOfEvents) + " events."
                    weight = ev.genInfo.weights().at(0)
                    sumOfAllWeights += weight
                    if weight == -1:
                        numberOfEventsWithWeightOne +=1
                        if self.beVerbose:
                                print "      Weight: " + str(weight)
                    if weight < 0:
                        numberOfEventsWithNegWeight += 1
                        sumOfNegWeights += weight
                        if self.beVerbose:
                            print "      Weight: " + str(weight)
                            for el in ev.genInfo.weights():
                                print "      all weights: " + str(el)
                    elif weight == 0:
                        numberOfEventsWithZeroWeight += 1
                        if self.beVerbose:
                            print "      Weight: " + str(weight)
                            for el in ev.genInfo.weights():
                                print "      all weights: " + str(el)


                file_.Close()

            file_list_.close()

            sampleFileNameOnly =  sample.rsplit('/',1)[1]
            f = open(os.path.join(self.outdir, "stats"+str(year) + "_" + sampleFileNameOnly +".txt"), "w")
            f.write("Created: " + str(datetime.now()) + "\n")
            f.write("Year: " + str(year) + "\n")
            f.write("This sample: " + sampleFileNameOnly + "\n")
            f.write("Total number of events: " + str(numberOfEvents) + "\n")
            f.write("Number of ntuple files: " + str(numberOfNtupleFiles) + "\n")
            f.write("Number of XML files: " +str(numberOfXMLfiles) + "\n")
            f.write("Number of events with negative weights: " +str(numberOfEventsWithNegWeight) + "\n")
            f.write("Number of events with weight == - 1: " +str(numberOfEventsWithWeightOne) + "\n")
            f.write("Sum of negative weights: " +str(sumOfNegWeights) + "\n")
            f.write("Sum of all weights: " +str(sumOfAllWeights) + "\n")
            f.close()

        f = open(os.path.join(self.outdir, "stats"+str(year) + ".txt"), "w")
        f.write("Created: " + str(datetime.now()) + "\n")
        f.write("Year: " + str(year) + "\n")
        f.write("Samples: " + str(self.samples) + "\n")
        f.write("Total number of events: " + str(numberOfEventsYear) + "\n")
        f.write("Number of ntuple files: " + str(numberOfNtupleFilesYear) + "\n")
        f.write("Number of XML files: " +str(numberOfXMLfilesYear) + "\n")
        f.write("Number of events with negative weights: " +str(numberOfEventsWithNegWeightYear) + "\n")
        f.write("Number of events with weight == - 1: " +str(numberOfEventsWithWeightOneYear) + "\n")
        f.write("Sum of negative weights: " +str(sumOfNegWeightsYear) + "\n")
        f.write("Sum of all weights: " +str(sumOfAllWeightsYear) + "\n")

        f.close()
        print("Info written to " + os.path.join(self.outdir, "stats"+str(year) +".txt"))


if __name__ == '__main__':
    args = parse_arguments()
    years       = args.years
    histFolders = args.histFolders
    Channels    = args.Channels

    studies = "nominal"
    module = "PlotZprimeSampleWeights"
    years = ["2016", "2017", "2018"]
    # years = ["2016"]

    print "Years: " + str(years)

    for year in years:
        AnalyzeYear = PlotZprimeSampleWeights(year=year, beVerbose=False)
        AnalyzeYear.AnalyzeYear()
