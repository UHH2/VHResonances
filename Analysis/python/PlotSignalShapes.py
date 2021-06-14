from Utils import *

import tdrstyle_all as TDR
TDR.writeExtraText = True
TDR.extraText  = "Simulation"
TDR.extraText2 = "Work in progress"

f_xmin = 1000
f_xmax = 6000
fNorms = {}

def CrystalBallIntegral(x, par):
    mean = par[0]
    sigma = par[1]
    alpha = par[2]
    n = par[3]
    std =(x[0]-mean)/sigma
    if (alpha < 0): std *= -1
    alpha = ROOT.TMath.Abs(alpha)
    result = 0
    if std >= -alpha:
        result = ROOT.TMath.Exp(-0.5*std*std)
    else:
        A = ROOT.TMath.Power(n/alpha, n)*ROOT.TMath.Exp(-0.5*alpha*alpha)
        B = n/alpha-alpha
        result = A/ROOT.TMath.Power(B-std, n)
    return result

def CrystalBall(x, par):
    norm = 1 if len(list(par))<=4 else par[4]
    mean = par[0]
    sigma = par[1]
    alpha = par[2]
    n = par[3]
    std =(x[0]-mean)/sigma
    if (alpha < 0): std *= -1
    alpha = ROOT.TMath.Abs(alpha)
    result = 0
    if std >= -alpha:
        result = ROOT.TMath.Exp(-0.5*std*std)
    else:
        A = ROOT.TMath.Power(n/alpha, n)*ROOT.TMath.Exp(-0.5*alpha*alpha)
        B = n/alpha-alpha
        result = A/ROOT.TMath.Power(B-std, n)
    # return result
    return norm*result*fNorms[str(par[0])+str(par[1])+str(par[2])+str(par[3])]

'''
Module to plot Tagger shapes
- Need the Selectionas input
'''

class PlotSignalShapes(VariablesBase):
    def __init__(self):
        VariablesBase.__init__(self)
        self.Samples = ["M"+str(m) for m in self.MassPointsReduced]
        self.Samples_reduced = ["M1600", "M2000", "M2500", "M3000", "M3500", "M4000", "M4500", "M5000"]
        self.Syst_list = self.Systematics[1:]+self.Systematics_Scale
        # self.Syst_list = ["trigger"]
        self.Syst = [""]+[x+v for x in self.Syst_list for v in ["Up", "Down"]]
        print self.Syst
        self.Folder = "DeepAk8_ZHccvsQCD_MD2"

        self.color  = {"M1600": ROOT.kAzure+2,
                       "M2000": ROOT.kAzure+7,
                       "M2500": ROOT.kGreen+3,
                       "M3000": ROOT.kGreen+2,
                       "M3500": ROOT.kOrange+1,
                       "M4000": ROOT.kOrange-2,
                       "M4500": ROOT.kRed+2,
                       "M5000": ROOT.kRed+1,
                       }
        self.style  = {"":           ROOT.kSolid,
                       "up":         ROOT.kDashed,
                       "down":       ROOT.kDotted,
                       "Nominal":    ROOT.kSolid,
                       "Syst. up":   ROOT.kDashed,
                       "Syst. down": ROOT.kDotted,
                       }
        self.outdir = self.Path_ANALYSIS+"Analysis/OtherPlots/SignalShapes/"
        os.system("mkdir -p "+self.outdir)

    def LoadFuncions(self):
        self.funcs = {}
        for year in ["RunII"]:
            for channel in self.Channels:
                for collection in self.Collections:
                    fname = self.Path_ANALYSIS+"/Analysis/Limits/nominal/"+year+"/"+collection+"/"+channel+"channel/"+self.Folder+"/datacards/SignalProperties_"+year+"_"+self.Folder+".txt"
                    with open(fname, 'r') as f_:
                        lines = f_.readlines()
                        for sample in self.Samples:
                            for sys in self.Syst:
                                if DoControl([""], year+channel+sys, channel, sample): continue
                                name = year+channel+collection+sample+sys
                                # print name
                                self.funcs[name] = ROOT.TF1(name, CrystalBall, f_xmin, f_xmax, 5)
                        for line in lines:
                            # if not any(x in line for x in self.Syst_list): continue
                            if "Up" in line and not any(x in line for x in self.Syst_list): continue
                            if "Down" in line and not any(x in line for x in self.Syst_list): continue
                            # if "Down" in line and not "btag" in line: continue
                            if "corrected" in line:
                                sys = line.split()[0].split("0")[-1]
                                sample = line.split()[0].replace(sys,"")
                                if not any(x in sample for x in self.Samples): continue
                                if not any(x in sys for x in self.Syst): continue
                                if DoControl([""], year+channel+sys, channel, sample): continue
                                self.funcs[year+channel+collection+sample+sys].SetParameter(4,float(line.split()[-1]))
                            if not "param" in line: continue
                            par = line.split()[0].replace("_"+channel+"_"+year,"")
                            var = float(line.split()[2])
                            sys = par.split("0")[-1]
                            sample = par.split("_")[2].replace(sys,"")
                            if not any(x in sample for x in self.Samples): continue
                            if not any(x in sys for x in self.Syst): continue
                            if DoControl([""], year+channel+sys, channel, sample): continue
                            # print line.split(), sample, sys
                            self.funcs[year+channel+collection+sample+sys].SetParameter(int(par.split("_")[1][1]),var)
                    for sample in self.Samples:
                        for sys in self.Syst:
                            if DoControl([""], year+channel+sys, channel, sample): continue
                            name = year+channel+collection+sample+sys
                            func = self.funcs[name]
                            fname = str(func.GetParameter(0))+str(func.GetParameter(1))+str(func.GetParameter(2))+str(func.GetParameter(3))
                            f_ = ROOT.TF1("", CrystalBallIntegral, f_xmin,f_xmax, 4)
                            f_.SetParameters(func.GetParameter(0),func.GetParameter(1),func.GetParameter(2),func.GetParameter(3))

                            fNorms[fname] = 1./f_.Integral(f_xmin,f_xmax)
                            # print name, fname, func.GetParameter(4), func.Integral(f_xmin,f_xmax), f_.Integral(f_xmin,f_xmax)

    def PlotFuncions(self):
        self.LoadFuncions()
        for channel in self.Channels:
            TDR.lumi_13TeV  = str(round(float(self.lumi_map["RunII"]["lumi_fb"]),1))+" fb^{-1}"
            isInv = "inv" in channel
            canv = tdrCanvas(channel, f_xmin, f_xmax, 0.0001, 0.036 if not isInv else 0.056, "M(Z') [GeV]" if not isInv else "M_{T}(Z') [GeV]" ,"A.U.", isExtraSpace=True)
            hframe = ROOT.gROOT.FindObject("hframe")
            hframe.GetXaxis().SetNdivisions(5)
            leg2 = tdrLeg(0.42,0.66,0.65,0.89, 0.035, 42, ROOT.kBlack)
            leg = tdrLeg(0.60,0.60,0.92,0.89, 0.035, 42, ROOT.kBlack)
            leg.SetNColumns(2)
            for sys in ["Nominal","Syst. up", "Syst. down"]:
                f_ = ROOT.TF1()
                f_.SetLineColor(ROOT.kBlack)
                f_.SetLineStyle(self.style[sys])
                leg2.AddEntry(f_, sys,"l")
            for year in ["RunII"]:
                for collection in self.Collections:
                    for sample in self.Samples_reduced:
                        for sys in self.Syst:
                            if DoControl([""], year+channel+sys, channel, sample): continue
                            func = self.funcs[year+channel+collection+sample+sys]
                            func.SetLineColor(self.color[sample])
                            func.SetLineStyle(self.style["down" if "Down" in sys else ("up" if "Up" in sys else "")])
                            func.SetNpx(1000)
                            func.Draw("same")
                            if sys=="": leg.AddEntry(func, "Z' "+str(round(float(sample.replace("M",""))/1000.,1))+" TeV", "l")
            canv.SaveAs(self.outdir+"SignalShape_"+channel+".pdf")



if __name__ == '__main__':
    PlotBkg = PlotSignalShapes()
    PlotBkg.PlotFuncions()
