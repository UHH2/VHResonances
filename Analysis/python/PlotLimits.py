from Utils import *

import tdrstyle_all as TDR
TDR.writeExtraText = True
TDR.extraText  = "Work in progress"
'''
Module to plot Limits
It takes as an input a txt file
'''

class PlotLimits(VariablesBase):
    def __init__(self, isObs=False):
        VariablesBase.__init__(self)
        self.isObs = isObs
        self.MassPoints = array('d',[float(x) for x in self.MassPointsReduced])
        self.MassPoints = array('d',[float(x) for x in [1400, 1600, 1800, 2000, 2500, 3000, 3500, 4000, 4500, 5000]])
        self.nPoints = len(self.MassPoints)
        self.outdir = self.Path_ANALYSIS+"Analysis/Limits/"
        os.system("mkdir -p "+self.outdir)
        TDR.cms_lumi = self.lumi_map['RunII']['lumiPlot']+' fb^{-1}'
        self.CreateDefaultGraphs()

        self.colors = { "muon":            ROOT.kOrange+1,
                        "electron":        ROOT.kAzure+1,
                        "charged lepton":  ROOT.kGreen+2,
                        "invisible":       ROOT.kViolet-3,
                        "neutrino":        ROOT.kViolet-3,
                        "lepton":          ROOT.kRed+1,
        }

    def CreateDefaultGraphs(self):
        xval = {}
        yval = {}
        self.graphs = {}
        self.PlotDetails = {}

        self.graphs["Default_obs"] = ROOT.TGraphErrors()
        self.graphs["Default_obs"].SetLineColor(ROOT.kBlack)
        self.graphs["Default_obs"].SetLineStyle(ROOT.kSolid)
        self.graphs["Default_obs"].SetFillStyle(0)
        self.graphs["Default_obs"].SetFillColor(ROOT.kBlack)

        self.graphs["Default_exp"] = ROOT.TGraphErrors()
        self.graphs["Default_exp"].SetLineColor(ROOT.kBlack)
        self.graphs["Default_exp"].SetLineStyle(ROOT.kDashed)
        self.graphs["Default_exp"].SetFillStyle(0)
        self.graphs["Default_exp"].SetFillColor(ROOT.kBlack)

        self.graphs["Default_band1"] = ROOT.TGraphErrors()
        self.graphs["Default_band1"].SetLineColor(ROOT.kGreen+2)
        self.graphs["Default_band1"].SetLineStyle(ROOT.kSolid)
        self.graphs["Default_band1"].SetFillStyle(1000)
        self.graphs["Default_band1"].SetFillColor(ROOT.kGreen+2)

        self.graphs["Default_band2"] = ROOT.TGraphErrors()
        self.graphs["Default_band2"].SetLineColor(ROOT.kOrange)
        self.graphs["Default_band2"].SetLineStyle(ROOT.kSolid)
        self.graphs["Default_band2"].SetFillStyle(1000)
        self.graphs["Default_band2"].SetFillColor(ROOT.kOrange)

        xval["theo"] = [800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000, 3100, 3200, 3300, 3400, 3500, 3600, 3700, 3800, 3900, 4000, 4100, 4200, 4300, 4400, 4500, 4600, 4700, 4800, 4900, 5000]
        yval["HVT_A"] = [578.0999, 358.3873, 231.6927, 154.8033, 106.2204, 74.50991, 53.23298, 38.62051, 28.39542, 21.11701, 15.85439, 12.00795, 9.163542, 7.039392, 5.439621, 4.225485, 3.297668, 2.584068, 2.032477, 1.604015, 1.269687, 1.009061, 0.8013809, 0.6388036, 0.5101927, 0.4081523, 0.3269965, 0.2623107, 0.2106453, 0.1693069, 0.1361855, 0.109623, 0.08828075, 0.07113544, 0.05733185, 0.04620448, 0.03724663, 0.03002855, 0.02420639, 0.01950782, 0.01571446, 0.01265338, 0.01018433]
        yval["HVT_B"] = [485.1621, 367.2001, 263.8817, 188.6818, 135.7415, 98.58745, 72.32583, 53.57184, 40.04648, 30.18608, 22.91748, 17.51997, 13.47566, 10.42173, 8.100062, 6.323773, 4.956944, 3.899242, 3.077397, 2.435993, 1.93346, 1.540306, 1.225945, 0.9791569, 0.7834139, 0.627741, 0.5036664, 0.4045769, 0.3252906, 0.2617513, 0.2107641, 0.1698193, 0.1368799, 0.1103879, 0.08903561, 0.07180615, 0.0579232, 0.0467275, 0.03768939, 0.03039015, 0.02449322, 0.01973149, 0.01588843]
        yval["theo_xsec_inv"]     = [15703.333, 5512.333, 2544.333, 1340.000, 763.866, 456.833, 283.366, 179.533, 62.176, 23.075, 9.119, 3.824, 1.716, 0.803, 0.385, 0.185, 0.041, 0.007]
        yval["theo_xsec_err_inv"] = [26.24, 13.02, 4.9, 3.5, 2.2, 0.7, 0.8, 0.2, 0.036, 0.005, 0.013, 0.007, 0.004, 0.001, 0.0007, 0.0004, 5.55e-05, 1.35e-05]

        xval["ZeeH0b"] = [ 800, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 5000]
        yval["ZeeH0b"] = [  60,   60,   25,   11,  7.5,  5.5,  4.5,    4,    4]

        xval["Hbb"] = [ 800, 900, 1000, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000, 3100, 3200, 3300, 3400, 3500, 3600, 3700, 3800, 3900, 4000, 4100, 4200, 4300, 4400, 4500, 4600, 4700, 4800, 4900, 5000]
        yval["Hbb_exp"]   = [ 18.675, 10.9125, 7.3125, 4.2937, 3.4781, 2.8594, 2.4281, 2.0719, 1.7906, 1.5609, 1.3781, 1.2328, 1.1156, 1.0125, 0.9281, 0.8578, 0.8016, 0.7547, 0.7125, 0.6797, 0.6539, 0.6328, 0.6141, 0.6000, 0.5906, 0.5813, 0.5766, 0.5813, 0.5859, 0.5953, 0.6094, 0.6234, 0.6492, 0.6797, 0.7172, 0.7594, 0.8156, 0.8203, 0.8297, 0.8438, 0.8578, 0.8766]
        yval["Hbb0b_exp"] = [ 14.325, 15.8625, 17.850, 17.625, 14.625, 11.962, 9.7312, 7.9313, 6.3938, 5.2875, 4.4437, 3.7875, 3.2625, 2.8312, 2.4844, 2.1844, 1.9312, 1.7062, 1.5141, 1.3500, 1.2047, 1.0828, 0.9797, 0.8930, 0.8203, 0.7570, 0.7031, 0.6609, 0.6281, 0.5977, 0.5719, 0.5484, 0.5344, 0.5238, 0.5180, 0.5156, 0.5203, 0.4945, 0.4734, 0.4570, 0.4453, 0.4336]

        self.graphs["HVT_A"]  = ROOT.TGraphErrors(len(xval["theo"]),   array('d',xval["theo"]),   array('d',yval["HVT_A"]))
        self.graphs["HVT_B"]  = ROOT.TGraphErrors(len(xval["theo"]),   array('d',xval["theo"]),   array('d',yval["HVT_B"]))
        self.graphs["Hbb"]    = ROOT.TGraphErrors(len(xval["Hbb"]),    array('d',xval["Hbb"]),    array('d',yval["Hbb_exp"]))
        self.graphs["Hbb0b"]  = ROOT.TGraphErrors(len(xval["Hbb"]),    array('d',xval["Hbb"]),    array('d',yval["Hbb0b_exp"]))
        self.graphs["ZeeH0b"] = ROOT.TGraphErrors(len(xval["ZeeH0b"]), array('d',xval["ZeeH0b"]), array('d',yval["ZeeH0b"]))

        # self.PlotDetails["HVT_A"]  = {"color": ROOT.kViolet-9, "lstyle": ROOT.kSolid, "LegName": "HVT model A"}
        # self.PlotDetails["HVT_B"]  = {"color": ROOT.kViolet-1, "lstyle": ROOT.kSolid, "LegName": "HVT model B"}
        self.PlotDetails["HVT_A"]  = {"color": ROOT.kRed+2,    "lstyle": ROOT.kSolid,  "LegName": "HVT model A"}
        self.PlotDetails["HVT_B"]  = {"color": ROOT.kAzure+10, "lstyle": ROOT.kSolid,  "LegName": "HVT model B"}
        self.PlotDetails["Hbb"]    = {"color": ROOT.kAzure,    "lstyle": ROOT.kDotted, "LegName": "B2G-19-006 2b cat."}
        self.PlotDetails["Hbb0b"]  = {"color": ROOT.kAzure,    "lstyle": ROOT.kDashed, "LegName": "B2G-19-006 0b cat."}
        self.PlotDetails["ZeeH0b"] = {"color": ROOT.kAzure,    "lstyle": ROOT.kDashed, "LegName": "#splitline{B2G-19-006}{0b cat. (Z#rightarrow ee)}"}
# arXiv:2102.08198, CMS-B2G-19-006, CERN-EP-2021-009

    def CreateGraphDefault(self, GraphName = "chargedlepton", FileName = "Limits_ZHcc_MD"):
        vals = {}
        for line in open(self.outdir+FileName+"_"+GraphName+".txt").readlines():
            name = line.split()[0]
            vals[name] = array('d',[float(x.replace(",","").replace("[","").replace("]","")) for x in line.split()[2:]])
        if self.isObs and not "obs" in vals:
            vals["obs"] = array('d',np.array([-1.]*len(vals["exp"])))
        dummy = array('d',np.array([0]*len(vals["exp"])))
        self.CreateGraph(vals, GraphName)

    def CreateGraphIterative(self, channel = "muon", FileName = "PN_ZHccvsQCD_MD", year ="RunII", mode = "Exp_2_AsimovNoSys_AsymptoticLimits", extraName=""):
        vals = {}
        for x in ["obs", "exp", "exp_1sigmaNeg", "exp_2sigmaNeg", "exp_1sigmaPos", "exp_2sigmaPos"]:
            vals[x] = np.array([-1.]*self.nPoints)
        for index, mass in enumerate(self.MassPoints):
            fname = self.outdir+"nominal/"+year+"/Puppi/"+channel+"channel/"+FileName+"/datacards/out/DataCard_"+year+"_Puppi_"+channel+"channel_"+FileName+"_M"+str(int(mass))+"_"+mode+".out"
            for line in open(fname).readlines():
                if not "r <" in line: continue
                if "Observed Limit:" in line: name = "obs"
                if "Expected  2.5%:" in line: name = "exp_2sigmaNeg"
                if "Expected 16.0%:" in line: name = "exp_1sigmaNeg"
                if "Expected 50.0%:" in line: name = "exp"
                if "Expected 84.0%:" in line: name = "exp_1sigmaPos"
                if "Expected 97.5%:" in line: name = "exp_2sigmaPos"
                vals[name][index] = line.split()[4]
        # print channel, FileName
        for x in vals:
            if "sigmaNeg" in x: vals[x] = vals["exp"] - vals[x]
            if "sigmaPos" in x: vals[x] = vals[x] - vals["exp"]
            vals[x] = array('d',vals[x])
            # print x, "=", [round(el,3) for el in vals[x]]
        self.CreateGraph(vals, channel, extraName)

    def CreateGraph(self, vals, GraphName, extraName=""):
        nPoints = self.nPoints
        xval = self.MassPoints
        # if nPoints != len(vals["exp"]):
        #     nPoints = len(vals["exp"])
        #     xval = array('d',[float(x) for x in self.MassPointsReduced])
        dummy = array('d',np.array([0]*nPoints))
        self.graphs[extraName+GraphName+"band1"] = ROOT.TGraphAsymmErrors(nPoints, xval, vals["exp"], dummy, dummy, vals["exp_1sigmaNeg"], vals["exp_1sigmaPos"])
        self.graphs[extraName+GraphName+"band2"] = ROOT.TGraphAsymmErrors(nPoints, xval, vals["exp"], dummy, dummy, vals["exp_2sigmaNeg"], vals["exp_2sigmaPos"])
        self.graphs[extraName+GraphName+"exp"]   = ROOT.TGraphErrors(nPoints, xval, vals["exp"])
        if self.isObs:
            self.graphs[extraName+GraphName+"obs"]   = ROOT.TGraphErrors(nPoints, xval, vals["obs"])

        self.PlotDetails[extraName+GraphName+"band1"] = {"color": ROOT.kGreen+2, "lstyle": ROOT.kSolid, "LegName": "68% Expected"}
        self.PlotDetails[extraName+GraphName+"band2"] = {"color": ROOT.kOrange,  "lstyle": ROOT.kSolid, "LegName": "95% Expected"}
        self.PlotDetails[extraName+GraphName+"exp"]   = {"color": ROOT.kBlack,   "lstyle": ROOT.kDashed, "LegName": "Expected"}
        if self.isObs:
            self.PlotDetails[extraName+GraphName+"obs"]   = {"color": ROOT.kBlack,   "lstyle": ROOT.kSolid, "LegName": "Observed"}

    def SaveCanvas(self, FileName, isStored):
        self.canv.Update()
        self.canv.RedrawAxis()
        self.canv.SaveAs(self.outdir+"Limits_"+FileName+("_Final" if isStored else "")+("_obs" if self.isObs else "")+".pdf")

    def CreateCanvas(self, PlotName, headerName = "", isStored=False):
        self.canv = tdrCanvas(PlotName+("_Final" if isStored else ""), 1200, 5200, 0.5e-01, 400, "M(Z') [GeV]", "#sigma#left(pp#rightarrowZ'#right) #times Br#left(Z'#rightarrow ZH #right)#left[fb#right]")
        self.canv.SetLogy(1)
        isDeepAk8 = "PN" in headerName
        nElementsLegend = 0.07*3 if isDeepAk8 else 0.075*3 if "Expected" in headerName else 0.1*3
        nElementsTheory = 0.08*2 if isDeepAk8 else (0.05*4 if "muon" in PlotName or "ele" in PlotName else 0.047*4)
        if headerName!="" and isDeepAk8:
            self.leg  = tdrLeg(0.63, 0.85-nElementsLegend, 0.89, 0.85, 0.04, 42, ROOT.kBlack)
        else:
            self.leg  = tdrLeg(0.68 if not "Expected" in headerName else 0.70, 0.85-nElementsLegend, 0.89, 0.85, 0.04, 42, ROOT.kBlack)
        self.leg.SetHeader("95% CL upper limit" if headerName=="" else headerName, "L")
        # if not "Expected" in headerName or "Channels" o:
        self.leg_theo = tdrLeg(0.40, 0.85-nElementsTheory, 0.60, 0.85, 0.04, 42, ROOT.kBlack)
        self.leg_theo.SetHeader("Theory prediction", "L")
        if not "Expected" in headerName:
            if (self.isObs): self.leg.AddEntry(self.graphs["Default_obs"], "Observed", "L")
            self.leg.AddEntry(self.graphs["Default_exp"], "Expected", "L")
            self.leg.AddEntry(self.graphs["Default_band1"], "68% Expected", "CF")
            self.leg.AddEntry(self.graphs["Default_band2"], "95% Expected", "CF")

    def PlotLines(self, GraphName, isBand = False, LegName="", Color=None, isObs=False, LStyle=None):
        if not GraphName in self.graphs:
            print "ERROR: "+GraphName+" not found"
        else:
            self.canv.cd()
            self.graphs[GraphName].SetLineWidth(2)
            color = self.PlotDetails[GraphName]["color"]
            lstyle = self.PlotDetails[GraphName]["lstyle"]
            if LegName=="":
                LegName = self.PlotDetails[GraphName]["LegName"]
            else:
                lstyle = LStyle if LStyle else ROOT.kSolid
                color = Color if Color else self.colors[LegName]
            tdrDraw(self.graphs[GraphName], "e3" if isBand else ("" if isObs else "C"), ROOT.kFullCircle, color, lstyle, color, 1000 if isBand else 0, color)
            if "Expected" in self.PlotDetails[GraphName]["LegName"]:
                if not isBand and LegName !="Expected":
                    self.leg.AddEntry(self.graphs[GraphName], LegName, "CF" if isBand else "L")
            elif not "Observed" in LegName:
                self.leg_theo.AddEntry(self.graphs[GraphName], LegName, "l")

    def PlotBand(self, GraphName):
        self.PlotLines(GraphName+"band2", isBand=True)
        self.PlotLines(GraphName+"band1", isBand=True)
        self.PlotLines(GraphName+"exp")
        if self.isObs:
            self.PlotLines(GraphName+"obs",isObs=True)

    def PlotTheoryLines(self):
        self.PlotLines("HVT_A")
        self.PlotLines("HVT_B")

    def PlotReferenceLines(self, isLep=False):
        if not isLep:
            self.PlotLines("Hbb")
            self.PlotLines("Hbb0b")
        else:
            self.PlotLines("ZeeH0b")

    def PlotLimitsStandard(self, PlotName, isStored =True):
        self.CreateCanvas(PlotName, isStored=isStored)
        if isStored:
            self.CreateGraphDefault(PlotName)
        else:
            self.CreateGraphIterative(PlotName)
        self.PlotBand(PlotName)
        self.PlotTheoryLines()
        if "chargedlepton" in PlotName:
            self.PlotLines("muon"+"exp",   LegName= "Expected Z#rightarrow #mu#mu", Color = ROOT.kAzure, LStyle=9)
            self.PlotLines("electron"+"exp",   LegName= "Expected Z#rightarrow ee", Color = ROOT.kRed+1, LStyle=9)
        self.SaveCanvas(PlotName+"_NoRef", isStored)
        self.PlotReferenceLines(("muon"in PlotName or "ele" in PlotName))
        self.SaveCanvas(PlotName, isStored)

    def PlotLimitsCompare(self, PlotName, isStored = False):
        self.CreateCanvas(PlotName, headerName="PN score")
        # histFolders = ["DeepAk8_H4qvsQCD", "DeepAk8_H4qvsQCD_massdep", "DeepAk8_H4qvsQCD_massdep_HccvsQCD", "DeepAk8_HccvsQCD", "DeepAk8_ZHccvsQCD_MD2", "tau42"]
        # histFolders = ["DeepAk8_H4qvsQCD", "DeepAk8_H4qvsQCD_massdep", 'PN_ZHccvsQCD_MD']
        histFolders = ['PN_ZHccvsQCD_MD']
        infos = {"DeepAk8_H4qvsQCD":         {"LegName":  "H4qvsQCD   1-dim.",    "Color": ROOT.kOrange+1},
                 "DeepAk8_H4qvsQCD_massdep": {"LegName":  "H4qvsQCD   p_{T}-dep.", "Color": ROOT.kAzure+2},
                 "PN_ZHccvsQCD_MD": {"LegName":  "ZHccvsQCD PN", "Color": ROOT.kRed+1},
                  }
        for fname in histFolders:
            LegName = infos[fname]["LegName"]
            if isStored:
                self.CreateGraphDefault(PlotName)
            else:
                self.CreateGraphIterative("lepton", FileName = fname, extraName = LegName)
            self.PlotLines(LegName+"lepton"+"exp", LegName = LegName, Color = infos[fname]["Color"])
        self.PlotLines("lepton"+"exp", LegName= "ZHccvsQCD 1-dim.", Color=ROOT.kGreen+2)
        self.PlotTheoryLines()
        self.SaveCanvas(PlotName, isStored)

    def CompareChannels(self, extraText=""):
        self.CreateCanvas("Channels"+extraText, "Expected")
        self.PlotLines("muon"+"exp",   LegName= "muon")
        self.PlotLines("electron"+"exp",   LegName= "electron")
        self.PlotLines("chargedlepton"+"exp",   LegName= "charged lepton")
        self.PlotLines("invisible"+"exp",   LegName= "neutrino")
        self.PlotLines("lepton"+"exp",   LegName= "lepton")
        # self.PlotTheoryLines()
        self.PlotReferenceLines()
        self.SaveCanvas("Channels"+extraText,False)


if __name__ == '__main__':
    PlotBkg = PlotLimits(isObs = True)
    for channel in ["muon", "electron", "chargedlepton", "invisible", "lepton"]:
        # PlotBkg.PlotLimitsStandard(channel, isStored=True)
        PlotBkg.PlotLimitsStandard(channel, isStored=False)
    # PlotBkg.PlotLimitsCompare("DeepAk8scores", isStored=False)
    PlotBkg.CompareChannels()
