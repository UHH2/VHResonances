from Utils import *

import tdrstyle_all as TDR
import datetime
TDR.writeExtraText = True
TDR.extraText  = "Simulation"
TDR.extraText2 = "Private work"

#ForThesis(TDR)

'''
Module to plot the effect of PDFReweight on the invisiblechannel

'''

class PlotHistograms(VariablesBase):
    def __init__(self,years):
        VariablesBase.__init__(self)
        self.outdir = self.Path_ANALYSIS+"Analysis/OtherPlots/PDFScaleVariations/"
        os.system("mkdir -p "+self.outdir)
        self.years = years
        self.module = "PDFReweight"
        self.Modes = ["nocuts", "weights"]
	self.Norm = True
        self.Variables = ["pt_jet","eta_jet"]
        self.MassPoints = [str(m) for m in self.MassPoints ]
        self.histos = {}
        self.color  = {"600":  rt.kBlack,
                       "800":  rt.kGreen+2,
                       "1000": rt.kViolet+1,
                       "1200": rt.kRed+1,
                       "1400": rt.kAzure+10,
                       "1600": rt.kMagenta+1,
                       "1800": rt.kGreen+1,
                       "2000": rt.kGreen+2,
                       "2500": rt.kOrange+1,
                       "3000": rt.kOrange+1,
                       "3500": rt.kOrange+1,
                       "4000": rt.kGray,
                       "4500": rt.kOrange+1,
                       "5000": rt.kRed+1,
                       "5500": rt.kOrange+1,
                       "6000": rt.kGray,
                       "7000": rt.kOrange+1,
                       "8000": rt.kAzure-2,
                       }

    def LoadHistos(self):
        for year in (self.years):
            for mass in self.MassPoints:
                filename = self.Path_STORAGE+year+"/"+self.module+"/"+"Puppi/invisiblechannel/"+"nominal"+"/"+self.PrefixrootFile+"MC."+self.Signal+"_inv_M"+mass+"_"+year+"_noTree.root"
                file_ = ROOT.TFile(filename)
                for mode in self.Modes:
                    for var in self.Variables:
                        h_ = file_.Get("nTopJet_"+mode+"/"+var)
                        h_.SetDirectory(0)
                        self.histos[year+mass+mode+var] = h_
                file_.Close()

    def CreateCanvas(self, nameVar):
        xName = "#eta^{jet}" if "eta" in nameVar else "p_{T}^{jet}"
        xMin = -3 if "eta" in nameVar else 50
        xMax =  3 if "eta" in nameVar else 4500

	if self.Norm:
	    yMin = -0.1 if "eta" in nameVar else -0.15
            yMax =  0.15 if "eta" in nameVar else 0.35
	else:
            yMin = -8000 if "eta" in nameVar else -10000
            yMax =  10000 if "eta" in nameVar else 30000

        self.canv = tdrCanvas("canv_"+nameVar, xMin,xMax,yMin,yMax, xName, "A.U." if self.Norm else "Events", isExtraSpace=True)
        self.leg = tdrLeg(0.77, 0.55, 0.95, 0.85, 0.045, 42, ROOT.kBlack)
        self.leg2 = tdrLeg(0.6, 0.75, 0.7, 0.85, 0.045, 42, ROOT.kBlack)

    def PlotHistos(self, isReduced=True):
        if isReduced:
            MassPoints = list(filter(lambda x: "2000" in x or "3000" in x or "5000" in x or "8000" in x, self.MassPoints))
        else: MassPoints = self.MassPoints
        for year in years:
            for var in self.Variables:
                TDR.lumi_13TeV  = str(round(float(self.lumi_map[year]["lumi_fb"]),1))+" fb^{-1}"
                self.CreateCanvas(year+var)
                for mode in self.Modes:
                    for mass in MassPoints:
                        lineStyle = rt.kSolid if "weights" in mode else rt.kDashed
                        h = self.histos[year+mass+mode+var]
                        h.SetLineWidth(2)
                        if self.Norm: h.Scale(1./self.histos[year+mass+"weights"+var].Integral())
                        tdrDraw(h, "hist" , rt.kFullCircle, self.color[mass], lineStyle, self.color[mass], 0, self.color[mass])
                        if "weights" in mode:
                            self.leg.AddEntry(h, "Z' "+mass.replace("000","").replace("00","0").replace("50",".5")+"TeV", "l")
                        if mass == MassPoints[0]:
                            h_ = h.Clone(h.GetName()+"_leg")
                            h_.SetDirectory(0)
                            h_.SetLineColor(rt.kBlack)
                            h_.SetLineStyle(lineStyle)
                            self.histos[year+mass+mode+var+"_leg"] = h_
                            self.leg2.AddEntry(h_, "Reweight" if "weights" in mode else "Original" , "l")

                self.leg.Draw("same")
                self.canv.SaveAs(self.outdir+"PDFReweight_"+year+"_"+var+".pdf")

if __name__ == '__main__':
    args = parse_arguments()
    years       = args.years if len(args.years)!=0 else ["2016","2017","2018"]

    PDF = PlotHistograms(years=years)
    PDF.LoadHistos()
    PDF.PlotHistos()
