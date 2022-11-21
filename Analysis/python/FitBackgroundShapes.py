from Utils import *

import tdrstyle_all as TDR
TDR.writeExtraText = True
TDR.extraText  = "Work in progress"

'''
Module to plot Tagger shapes
- Need the Selectionas input
'''


p_xmin = 700
p_xmax = 4000

def PolinomialExponent(x, par):
    xp = x[0]/1000.
    result = 0
    for i_ in range(len(list(par))):
        result += par[i_]*ROOT.TMath.Power(xp, i_)
    return ROOT.TMath.Exp(result)


def dijet(x, par):
    xp = x[0]/1000.
    return par[2]*(ROOT.TMath.Power(5.-xp,par[0])/ROOT.TMath.Power(xp,par[1]))

def dijet_p2(x, par):
    xp = x[0]/1000.
    return par[2]*ROOT.TMath.Exp(-(0.8026*par[0] + 0.5965*par[1])*xp + (0.5965*par[0] - 0.8026*par[1])*xp*xp)

def dijet_p3(x, par):
    xp = x[0]/1000.
    return par[3]*(ROOT.TMath.Power(5.-xp,par[0])/(ROOT.TMath.Power(xp,par[1]+par[2]*ROOT.TMath.Log(xp))))

def dijet_p4(x, par):
    xp = x[0]/1000.
    return par[4]*(ROOT.TMath.Power(5.-xp,par[0])/(ROOT.TMath.Power(xp,par[1]+par[2]*ROOT.TMath.Log(xp)+ par[3]*ROOT.TMath.Log(xp)*ROOT.TMath.Log(xp))))

def FTest(chi2_1, chi2_2, npar_1, npar_2, n_):
    return 1. - ROOT.TMath.FDistI(((chi2_1-chi2_2)/(npar_2-npar_1))/(chi2_2/(n_-npar_2-1)), npar_2-npar_1, n_-npar_2)


class FitBackgroundShapes(VariablesBase):
    def __init__(self, xmin=-1, xmax=-1):
        VariablesBase.__init__(self)
        self.xmin = xmin
        self.xmax = xmax
        self.Folder = "DeepAk8_ZHccvsQCD_MD"
        self.Folder = "PN_ZHccvsQCD_MD"
        self.degrees = ["1","2","3"]

        self.color  = {"1": ROOT.kAzure+2,
                       "2": ROOT.kGreen+2,
                       "3": ROOT.kRed+1,
                       "4": ROOT.kOrange+1,
                       "5": ROOT.kOrange-2,
                       }
        self.style  = {"":           ROOT.kSolid,
                       "up":         ROOT.kDashed,
                       "down":       ROOT.kDotted,
                       "Nominal":    ROOT.kSolid,
                       "Syst. up":   ROOT.kDashed,
                       "Syst. down": ROOT.kDotted,
                       }

        self.ranges  = {
            "muon": {
                "MC_SR":   (1000, 2700),
                "MC_CR":   (1100, 6000),
                "DATA_SR": (1000, 3000),
                "DATA_CR": (1000, 3000),
                },
            "electron": {
                "MC_SR":   (1000, 3000),
                "MC_CR":   (1100, 3600),
                "DATA_SR": (900, 3000),
                "DATA_CR": (1000, 3300),
                },
            "invisible": {
                "MC_SR":   (1000, 4000),
                "MC_CR":   (1100, 6000),
                "DATA_SR": (1000, 3000),
                "DATA_CR": (1200, 5000),
                },
            }
        self.outdir = self.Path_ANALYSIS+"Analysis/OtherPlots/BackgroundShapes/"
        os.system("mkdir -p "+self.outdir)

    def LoadHistos(self, year, channel, collection="Puppi"):
        self.hists = {}
        self.bins = {}
        fname = self.Path_STORAGE+"/"+year+"/SignalRegion/"+collection+"/"+channel+"channel/nominal/"+self.PrefixrootFile
        for fn in glob.glob(fname+"*.root"):
            if "Zprime" in fn: continue
            if "VV" in fn: continue
            if "Ele" in fn and not "ele" in channel: continue
            if "Muo" in fn and not "muon" in channel: continue
            if "MET" in fn and not "inv" in channel: continue
            if "WJets" in fn and not "inv" in channel: continue
            f_ = ROOT.TFile(fn)
            for sr in ["_SR", "_CR"]:
                h_ = f_.Get("ZprimeCandidate_"+self.Folder+sr+"/Zprime_mass"+("" if not "inv" in channel else "_transversal")+"_rebin100")
                hname = fn.replace(fname,"").replace("_"+year+"_noTree.root","").replace("DATA.","").replace("MC.","").split("_")[0]
                hname += sr
                for bin in range(h_.GetNbinsX()):
                    if h_.GetBinCenter(bin) < 500 or h_.GetBinCenter(bin) > self.ranges[channel][hname][1]:
                        h_.SetBinContent(bin,0)
                        h_.SetBinError(bin,0)
                self.bins[hname] = []
                if hname in self.hists:
                    self.hists[hname].Add(h_)
                else:
                    self.hists[hname] = h_
                    self.hists[hname].SetDirectory(0)
            f_.Close()
        for bin in range(self.hists["MC_CR"].GetNbinsX()):
            if "muon" in channel:
                self.hists["MC_CR"].SetBinError(bin,1.5*self.hists["MC_CR"].GetBinError(bin))
            if "inv" in channel:
                if self.hists["MC_CR"].GetBinCenter(bin) > 2800:
                    self.hists["MC_CR"].SetBinError(bin,2*self.hists["MC_CR"].GetBinError(bin))
                elif self.hists["MC_CR"].GetBinCenter(bin) > 2500:
                    self.hists["MC_CR"].SetBinError(bin,3*self.hists["MC_CR"].GetBinError(bin))
                elif self.hists["MC_CR"].GetBinCenter(bin) > 2000:
                    self.hists["MC_CR"].SetBinError(bin,5*self.hists["MC_CR"].GetBinError(bin))
                elif self.hists["MC_CR"].GetBinCenter(bin) > 1500:
                    self.hists["MC_CR"].SetBinError(bin,10*self.hists["MC_CR"].GetBinError(bin))
                else: self.hists["MC_CR"].SetBinError(bin,25*self.hists["MC_CR"].GetBinError(bin))

    def DoFits(self):
        ftest_tot = {}
        ftest_store = {}
        for d_ in self.degrees:
            ftest_tot[d_] = 0
        for year in ["RunII"]:
            for channel in self.Channels:
            # for channel in ["invisible"]:
                for collection in self.Collections:
                    self.LoadHistos(year,channel,collection)
                    for hname in self.hists:
                        unique_name = year+"_"+channel+"_"+collection+"_"+hname
                        hist = self.hists[hname]
                        self.funcs = {}
                        self.chi2 = {}
                        f_xmin, f_xmax = self.ranges[channel][hname]
                        if self.xmin!=-1: f_xmin = self.xmin
                        if self.xmax!=-1: f_xmax = self.xmax
                        for func_name in self.degrees:
                            self.funcs[func_name] = ROOT.TF1(unique_name+func_name, PolinomialExponent, f_xmin, f_xmax, int(func_name)+1)
                        self.errors = {}
                        for func_name in self.degrees if not "DATA_SR" in hname else reversed(self.degrees) if "ele" in channel else ["3","1","2"]:
                            func = self.funcs[func_name]
                            func.SetLineColor(self.color[func_name])
                            if "DATA_SR" in hname and "ele" in channel: func.SetLineColor(self.color["2"])
                            fitRes = hist.Fit(unique_name+func_name, "RMQS" if hname == "DATA_SR" else "RMQS+")
                            self.errors["band"+func_name], self.errors["band_pull"+func_name], self.errors["pull"+func_name] = ComputeHistWithCL(unique_name+func_name, func, fitRes, hist, cl=0.68)
                            self.errors["band2"+func_name], self.errors["band_pull2"+func_name], self.errors["pull2"+func_name] = ComputeHistWithCL(unique_name+func_name, func, fitRes, hist, cl=0.95)
                            self.chi2[func_name] = (func.GetChisquare(), func.GetNDF(), func.GetProb())
                        TDR.lumi_13TeV  = str(round(float(self.lumi_map[year]["lumi_fb"]),1))+" fb^{-1}"
                        isInv = "inv" in channel
                        ymax = 5*1e02
                        ymin = 0.00101
                        if "CR" in hname:
                            ymin = 0.101
                            ymax = 1e05
                        if "DATA" in hname:
                            ymin = 1
                            ymax = 1e04
                        if "DATA_SR" in hname:
                            ymin = 0.101
                            ymax = 1e03
                        if "MC_SR" in hname:
                            ymin = 0.0501
                        if isInv:
                            ymin = 0.101
                            ymax = 1e07
                            if "DATA_SR" in hname:
                                ymax = 1e05
                        # if not isInv and "SR" in hname: ymax = 1e02
                        canv = tdrDiCanvas(unique_name, p_xmin, f_xmax+150, ymin, ymax, -1, 3, "M(Z') [GeV]" if not isInv else "M_{T}(Z') [GeV]" ,"Events", "Data / Fit" if "DATA" in hname else "Sim. / Fit")
                        canv.cd(1).SetLogy(1)
                        leg = tdrLeg(0.43,0.89-0.08*(2 if "DATA_SR" in unique_name else 4),0.83,0.89, 0.040, 42, ROOT.kBlack)
			leg.AddEntry(hist, "Data" if "DATA" in unique_name else "Simulation" ,"lp")
                        for func_name in self.degrees:
                            if func_name!="2" and "DATA_SR" in unique_name: continue
                            canv.cd(1)
                            color = self.color[func_name]
                            func = self.funcs[func_name]
                            band = self.errors["band"+func_name]
                            band2 = self.errors["band2"+func_name]
                            pull = self.errors["pull"+func_name]
                            band_pull = self.errors["band_pull"+func_name]
                            band_pull2 = self.errors["band_pull2"+func_name]
                            if "DATA_SR" in unique_name and "ele" in unique_name:
                                func = self.funcs["1"]
                                band = self.errors["band"+"1"]
                                band2 = self.errors["band2"+"1"]
                                pull = self.errors["pull"+"1"]
                                band_pull = self.errors["band_pull"+"1"]
                                band_pull2 = self.errors["band_pull2"+"1"]
                            chi2_red = self.chi2[func_name][0]/self.chi2[func_name][1]
                            prob_new = self.chi2[func_name][2]
                            ftest  = FTest(self.chi2[str(int(func_name)-1)][0], self.chi2[func_name][0], int(func_name), int(func_name)+1, band.GetNbinsX()) if func_name!="1" else 0
                            ftest_store[unique_name+func_name+"vs"+str(int(func_name)-1)] = round(ftest*100,2)
                            func.Draw("same")
                            if func_name != "1" and hname!="DATA_SR":
                                d_ = int(func_name)
                                if ftest>0.05:
                                    ftest_tot[str(d_-1)] +=1
                                    print unique_name+func_name+"vs"+str(int(func_name)-1), ftest, str(d_-1)
                                else:
                                    prob_new = self.chi2[func_name][2]
                                    chi2_old = self.chi2[str(d_-1)][0]/self.chi2[str(d_-1)][1]
                                    if d_==3 and (prob_new>0.95 or prob_new<0.05 or (chi2_old>0.9 and chi2_old<1.1) ):
                                        ftest_tot[str(d_-1)] +=1
                                        print unique_name+func_name+"vs"+str(int(func_name)-1), ftest, str(d_-1)
                                    else:
                                        ftest_tot[str(d_)] +=1
                                        print unique_name+func_name+"vs"+str(int(func_name)-1), ftest, str(d_)
                            # leg.AddEntry(func, "N = "+func_name+" F-test = "+str(round(ftest,2))+" #chi^{2}"+"/n.d.f = "+str(round(chi2_red,2))+" p-value = "+str(round(self.chi2[func_name][2],2)), "l")
                            #if hname == "DATA_SR" and func_name=="1" and "ele" in channel:
                            #    leg.AddEntry(func, "N = 2,  #chi^{2}"+"/n.d.f = "+str(int(chi2_red*100)/100.)+", p-value = "+str(int(prob_new*100))+" %", "l")
                            #else:
                            #    leg.AddEntry(func, "N = "+func_name+",  #chi^{2}"+"/n.d.f = "+str(int(chi2_red*100)/100.)+", p-value = "+str(int(prob_new*100))+" %", "l")
                            leg.AddEntry(func, "f_{"+func_name+"} : #chi^{2}/n.d.f. = "+str(int(chi2_red*100)/100.)+", p-value = "+str(int(prob_new*100))+" %", "l")
                            # tdrDraw(band2, "e3", 9, ROOT.kWhite, 1, color, 1001, color)
                            band2.SetMarkerSize(0)
                            tdrDraw(band, "e3", 9, ROOT.kWhite, 1, color+1, 1001, color+1)
                            band.SetMarkerSize(0)
                            band.SetFillColorAlpha(color+1,0.35)
                            band2.SetFillColorAlpha(color,0.35)
                            canv.cd(2)
                            for n in range(0,pull.GetNbinsX()+1):
                                if pull.GetBinContent(n)==0 : pull.SetBinContent(n, -10)
                            if "DATA_SR" in unique_name:
                                tdrDraw(pull, "e0", ROOT.kFullCircle, ROOT.kBlack, 1, ROOT.kBlack, 0, ROOT.kBlack)
                            else:
                                tdrDraw(pull, "e0", ROOT.kFullCircle, color, 1, color, 0, color)
                            # tdrDraw(band_pull2, "e3", 9, ROOT.kWhite, 1, color, 1001, color)
                            band_pull2.SetMarkerSize(0)
                            tdrDraw(band_pull, "e3", 9, ROOT.kWhite, 1, color+1, 1001, color+1)
                            band_pull.SetMarkerSize(0)
                            band_pull.SetFillColorAlpha(color+1,0.35)
                            band_pull2.SetFillColorAlpha(color,0.35)
                            line_ =  rt.TLine(p_xmin, 1, f_xmax+150, 1)
                            line_.SetLineColor(ROOT.kBlack)
                            line_.SetLineStyle(ROOT.kSolid)
                            line_.Draw("same")
                        canv.cd(1)
                        for i in range(hist.GetNbinsX()):
                            if i in self.bins[hname]:
                                hist.SetBinContent(i,0)
                                hist.SetBinError(i,0)
                        tdrDraw(hist, "", ROOT.kFullCircle, ROOT.kBlack, 1, ROOT.kBlack, 0, ROOT.kBlack)
                        if self.xmin!=-1 or self.xmax!=-1:
                            canv.SaveAs(self.outdir+"BackgroundShape_"+unique_name+"_"+str(self.xmin)+"_"+str(self.xmax)+"_"+self.Folder+".pdf")
                        else:
                            canv.SaveAs(self.outdir+"BackgroundShape_"+unique_name+"_"+self.Folder+".pdf")
        print ftest_tot
        print ftest_store

if __name__ == '__main__':
    PlotBkg = FitBackgroundShapes()
    PlotBkg.DoFits()
    # for xmin in range(1000,1200+1,100):
    #     for xmax in range(5900,6000+1,100):
    #         if xmin >= xmax: continue
    #         PlotBkg = FitBackgroundShapes(xmin,xmax)
    #         PlotBkg.DoFits()
