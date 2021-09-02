from Utils import *

import tdrstyle_all as TDR
TDR.writeExtraText = True
TDR.extraText  = "Simulation"
TDR.extraText2 = "Work in progress"

'''
Module to plot Tagger shapes
- Need the Selectionas input
'''

class PlotTaggerPtDependence(VariablesBase):
    def __init__(self):
        VariablesBase.__init__(self)
        self.outdir = self.Path_ANALYSIS+"Analysis/OtherPlots/TaggerShape/"
        os.system("mkdir -p "+self.outdir)
        self.order = ["x3","x13","x23"]
        self.info = {"x3":  {"color": ROOT.kAzure+2,
                             "style": ROOT.kSolid,
                             "func":  "[0]+[1]*TMath::Power(x,-3)",
                             "name":  "a + b*m^{-3}",
                             },
                     "x13": {"color": ROOT.kGreen+2,
                             "style": ROOT.kDashed,
                             "func":  "[0]+[1]*TMath::Power(x,-1)+[2]*TMath::Power(x,-3)",
                             "name":  "a + b*m^{-1}+c*m^{-3}",
                             },
                     "x23": {"color": ROOT.kOrange+1,
                             "style": ROOT.kDashed,
                            "func":  "[0]+[1]*TMath::Power(x,-2)+[2]*TMath::Power(x,-3)",
                            "name":   "a + b*m^{-2}+c*m^{-3}",
                            },
                     }

    def Plot(self):
        x_values = np.array([800,1000,1200,1400,1600,1800,2000,2500,3000,3500,4000,4500,5000])
        y_values = np.array([0.0360,0.0193,0.0113,0.0075,0.0055,0.0044,0.0034,0.0035,0.0024,0.0034,0.0025,0.0034,0.0044])
        y_errors = np.array([0.0546,0.0242,0.0139,0.0089,0.0064,0.0065,0.0044,0.0044,0.0029,0.0044,0.0029,0.0039,0.0071])
        TDR.lumi_13TeV  = "MC RunII"
        canv = tdrCanvas("PlotPt", 600, 5200, 1e-03, 2*1e-01, "M(Z')[GeV]","H4qvsQCD threshold")
        canv.SetLogy(1)
        leg = tdrLeg(0.60,0.60,0.89,0.89, 0.045, 42, ROOT.kBlack)
        graph = rt.TGraphErrors(len(x_values), array('d',x_values), array('d',y_values), array('d',np.zeros(len(x_values))), array('d',np.abs(y_values-y_errors)))
        graph_fit = rt.TGraphErrors(len(x_values), array('d',x_values), array('d',y_values), array('d',np.zeros(len(x_values))), array('d',np.abs(y_values-y_errors)))
        leg.AddEntry(graph, "Simulation", "lp")
        funcs = {}
        rt.gStyle.SetOptFit(0)
        for fname in reversed(self.order):
            funcs[fname] = rt.TF1(fname,self.info[fname]["func"],700,5100)
            graph_fit.Fit(funcs[fname],"RQ")
            funcs[fname].SetLineColor(self.info[fname]["color"])
            funcs[fname].SetLineStyle(self.info[fname]["style"])
            funcs[fname].SetLineWidth(3)
            funcs[fname].Draw("same")
        for fname in self.order:
            leg.AddEntry(funcs[fname], "f = "+self.info[fname]["name"], "l")
        tdrDraw(graph, "P", ROOT.kFullCircle, ROOT.kBlack, ROOT.kDashed)
        canv.SaveAs(self.outdir+"TaggerPtDependence.pdf")



if __name__ == '__main__':
    PlotPt = PlotTaggerPtDependence()
    PlotPt.Plot()
