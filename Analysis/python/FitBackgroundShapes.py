from Utils import *

import tdrstyle_all as TDR
TDR.writeExtraText = True
TDR.extraText  = "Simulation"
TDR.extraText2 = "Work in progress"

'''
Module to plot Tagger shapes
- Need the Selectionas input
'''


p_xmin = 350
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


# dijetfunction_altp3  f = "[3]*( ROOT.TMath.Power(5.-xp + [2]*x*x/1e6,[0])/(ROOT.TMath.Power(xp,[1]) ))", f_xmin, f_xmax);
# dijetfunction_p4     f = "[4]*( ROOT.TMath.Power(5.-xp,[0])/(ROOT.TMath.Power(xp,[1] + [2]*ROOT.TMath.Log(xp) + [3]*ROOT.TMath.Log(xp)*ROOT.TMath.Log(xp) ) ) )", f_xmin, f_xmax);
# dijetfunction_altp4  f = "[4]*( ROOT.TMath.Power(5.-xp + [2]*x*x/1e6,[0])/(ROOT.TMath.Power(xp,[1]+[3]*ROOT.TMath.Log(xp)) ))", f_xmin, f_xmax);
# dijetfunction_p5     f = "[5]*( ROOT.TMath.Power(5.-xp,[0])) / (ROOT.TMath.Power(xp,[1] + [2]*ROOT.TMath.Log(xp) + [3]*ROOT.TMath.Log(xp)*ROOT.TMath.Log(xp) + [4]*ROOT.TMath.Power(ROOT.TMath.Log(xp),3) ) )", f_xmin, f_xmax);

def CountBinsMinMax(hist, min, max):
    Nbins = 0
    xmin = 1e6
    xmax = 0
    for i in range(hist.GetNbinsX()+1):
        if hist.GetXaxis().GetBinLowEdge(i)>=min:
            if (xmin>1e5): xmin = hist.GetXaxis().GetBinLowEdge(i)
            if hist.GetXaxis().GetBinUpEdge(i)<=max:
              xmax = hist.GetXaxis().GetBinUpEdge(i)
              Nbins += 1
    return (Nbins,xmin,xmax)

def GetConfidenceIntervals(func, fitRes, Nbins, xcenters, ci, cl = 0.68):
    npar = func.GetNumberFreeParameters()
    npar_real = func.GetNpar()
    fixed = 0
    if npar_real != npar:
        print "ERROR. Check other function."

    covmatr = fitRes.GetCovarianceMatrix()
    rho = fitRes.GetCorrelationMatrix()
    tStudent = ROOT.TMath.StudentQuantile(0.5 + cl/2, func.GetNDF())
    chindf = ROOT.TMath.Sqrt(func.GetChisquare()/func.GetNDF())

    grad = array('d',[0.]*npar_real)
    sum_vector = array('d',[0.]*npar)
    for ipoint in range(Nbins):
        ci_=0
        func.GradientPar(array('d',[xcenters[ipoint]]), grad)
        # multiply the covariance matrix by gradient
        for irow in range(npar):
            sum_vector[irow]=0;
            for icol in range(npar):
                igrad=0
                ifree=0
                if fixed:
                    print "ERROR. Check other function."
                else:
                    igrad = icol;
                sum_vector[irow] += covmatr[irow][icol]*grad[igrad]
        igrad = 0;
        for i_ in range(npar):
            igrad=0
            ifree=0
            if fixed:
                print "ERROR. Check other function."
            else:
                igrad = i_
            ci_ += grad[igrad]*sum_vector[i_]
        ci_ = ROOT.TMath.Sqrt(ci_)
        ci[ipoint] = ci_*tStudent*chindf

def ComputeHistWithCL(name, func, fitRes, hist, cl=0.68):
  # create a histogram for the fit region only
  Nbins,xmin,xmax = CountBinsMinMax(hist, func.GetXmin(), func.GetXmax())
  band = ROOT.TH1F("band"+name+str(cl),name+str(cl), Nbins, xmin, xmax)
  band_pull  = ROOT.TH1F("band_pull"+name+str(cl),name+str(cl), Nbins, xmin, xmax)
  pull = ROOT.TH1F("pull"+name+str(cl),name+str(cl), Nbins, xmin, xmax)

  # now get an array with the bin centers and compute the CLs
  xcenters = array('d',[0.]*Nbins)
  ci = array('d',[0.]*Nbins)
  for i in range(1,Nbins+1):
      xcenters[i-1] = band.GetXaxis().GetBinCenter(i)

  GetConfidenceIntervals(func, fitRes, Nbins, xcenters, ci, cl);

  for i in range(1,Nbins+1):
      x_ = band.GetXaxis().GetBinCenter(i)
      y_ = func.Eval(x_)
      ibin = hist.GetXaxis().FindBin(x_)
      band.SetBinContent(i, y_)
      band.SetBinError(i, ci[i-1])
      band_pull.SetBinContent(i, 1)
      band_pull.SetBinError(i, ci[i-1]/y_)
      pull.SetBinContent(i, hist.GetBinContent(ibin)/y_)
      pull.SetBinError(i, hist.GetBinError(ibin)/y_)
      # print name, x_, pull.GetBinContent(i), band_pull.GetBinError(i)
  return band, band_pull,pull



class FitBackgroundShapes(VariablesBase):
    def __init__(self, xmin=-1, xmax=-1):
        VariablesBase.__init__(self)
        self.xmin = xmin
        self.xmax = xmax
        self.Folder = "DeepAk8_ZHccvsQCD_MD"
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
                "MC_CR":   (1000, 3800),
                "DATA_SR": (900,  2100),
                "DATA_CR": (900,  2800),
                },
            "electron": {
                "MC_SR":   (1000, 2900),
                "MC_CR":   (1000, 3600),
                "DATA_SR": (1000, 1800),
                "DATA_CR": (800,  3300),
                },
            "invisible": {
                "MC_SR":   (1000, 4100),
                "MC_CR":   (1200, 5700),
                "DATA_SR": (1000, 3100),
                "DATA_CR": (1000, 4900),
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
            if "Ele" in fn and not "ele" in channel: continue
            if "Muo" in fn and not "muon" in channel: continue
            if "MET" in fn and not "inv" in channel: continue
            if "WJets" in fn and not "inv" in channel: continue
            f_ = ROOT.TFile(fn)
            for sr in ["_SR", "_CR"]:
                h_ = f_.Get("ZprimeCandidate_"+self.Folder+sr+"/Zprime_mass"+("" if not "inv" in channel else "_transversal")+"_rebin100")
                for bin in range(h_.GetNbinsX()):
                    if h_.GetBinCenter(bin) < 500:
                        h_.SetBinContent(bin,0)
                        h_.SetBinError(bin,0)
                hname = fn.replace(fname,"").replace("_"+year+"_noTree.root","").replace("DATA.","").replace("MC.","").split("_")[0]
                hname += sr
                self.bins[hname] = []
                if hname in self.hists:
                    self.hists[hname].Add(h_)
                else:
                    self.hists[hname] = h_
                    self.hists[hname].SetDirectory(0)
            f_.Close()
        for bin in range(self.hists["MC_CR"].GetNbinsX()):
            self.hists["MC_CR"].SetBinError(bin,1.5*self.hists["MC_CR"].GetBinError(bin))
        if channel=="muon":
            h_ = self.hists["MC_CR"]
            bin = h_.FindBin(2650)
            h_.SetBinContent(bin,0.6)
            h_.SetBinError(bin,h_.GetBinError(bin+1))
            h_ = self.hists["MC_SR"]
            bin = h_.FindBin(3050)
            h_.SetBinContent(bin,0)
            h_.SetBinError(bin,0)
        # for i in range(self.hists["MC_CR"].GetNbinsX()):
        #     if self.hists["MC_CR"].GetBinCenter(i)<1000:continue
        #     if self.hists["MC_CR"].GetBinCenter(i)>2000:continue
        #     print self.hists["MC_CR"].GetBinError(i)/self.hists["DATA_CR"].GetBinError(i), self.hists["MC_CR"].GetBinError(i), self.hists["DATA_CR"].GetBinError(i)
            # h_ = self.hists["MC_CR"]
        #     if h_.GetBinContent(i)<2*1e-03:
        #         h_.SetBinContent(i,0)
        #         h_.SetBinError(i,0)
        #         self.bins[hname].append(i)
                # h_.SetBinContent(i,(h_.GetBinContent(i-1)+h_.GetBinContent(i+1))/2)
                # h_.SetBinError(i,(h_.GetBinError(i-1)+h_.GetBinError(i+1))/2)

    def DoFits(self):
        ftest_tot = {}
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
                        for func_name in self.degrees:
                            func = self.funcs[func_name]
                            func.SetLineColor(self.color[func_name])
                            fitRes = hist.Fit(unique_name+func_name, "RMQS+")
                            self.errors["band"+func_name], self.errors["band_pull"+func_name], self.errors["pull"+func_name] = ComputeHistWithCL(unique_name+func_name, func, fitRes, hist, cl=0.68)
                            self.errors["band2"+func_name], self.errors["band_pull2"+func_name], self.errors["pull2"+func_name] = ComputeHistWithCL(unique_name+func_name, func, fitRes, hist, cl=0.95)
                            self.chi2[func_name] = (func.GetChisquare(), func.GetNDF(), func.GetProb())
                        TDR.lumi_13TeV  = str(round(float(self.lumi_map[year]["lumi_fb"]),1))+" fb^{-1}"
                        isInv = "inv" in channel
                        ymax = 1e07
                        if not isInv and "SR" in hname: ymax = 1e05
                        canv = tdrDiCanvas(unique_name, p_xmin, f_xmax+150, 0.00101, ymax, -1, 3, "M(Z') [GeV]" if not isInv else "M_{T}(Z') [GeV]" ,"Events", "Hist/Fit")
                        canv.cd(1).SetLogy(1)
                        leg = tdrLeg(0.40,0.60,0.92,0.89, 0.035, 42, ROOT.kBlack)
                        for func_name in self.degrees:
                            canv.cd(1)
                            color = self.color[func_name]
                            func = self.funcs[func_name]
                            func.Draw("same")
                            band = self.errors["band"+func_name]
                            band2 = self.errors["band2"+func_name]
                            pull = self.errors["pull"+func_name]
                            band_pull = self.errors["band_pull"+func_name]
                            band_pull2 = self.errors["band_pull2"+func_name]
                            chi2_red = self.chi2[func_name][0]/self.chi2[func_name][1]
                            prob_new = self.chi2[func_name][2]
                            ftest  = FTest(self.chi2[str(int(func_name)-1)][0], self.chi2[func_name][0], int(func_name), int(func_name)+1, band.GetNbinsX()) if func_name!="1" else 0
                            if func_name != "1" and hname!="DATA_SR" and channel!="invisible":
                                d_ = int(func_name)
                                prob_old = self.chi2[str(d_-1)][2]
                                if ftest>0.05 and prob_old>0.06 and prob_new < 0.94 and prob_old>1.05*prob_new:
                                    if prob_old>0.95*prob_new and prob_old<1.05*prob_new and d_>=2 and d_<=3:
                                        ftest_tot["2"] +=1
                                    else:
                                        ftest_tot[str(d_-1)] +=1
                                else:
                                    if (prob_old>0.95*prob_new and prob_old<1.05*prob_new and d_>=2 and d_<=3) or (prob_old<0.50*prob_new or prob_old>1.50*prob_new) :
                                        ftest_tot["2"] +=1
                                    else:
                                        ftest_tot[str(d_-1)] +=1
                            # leg.AddEntry(func, "N = "+func_name+" F-test = "+str(round(ftest,2))+" #chi^{2}"+"/n.d.f = "+str(round(chi2_red,2))+" p-value = "+str(round(self.chi2[func_name][2],2)), "l")
                            leg.AddEntry(func, "N = "+func_name+"     #chi^{2}"+"/n.d.f = "+str(int(chi2_red*100)/100.)+"     p-value = "+str(int(prob_new*100))+" %", "l")
                            tdrDraw(band2, "e3", 9, ROOT.kWhite, 1, color, 1001, color)
                            band2.SetMarkerSize(0)
                            tdrDraw(band, "e3", 9, ROOT.kWhite, 1, color+1, 1001, color+1)
                            band.SetMarkerSize(0)
                            band.SetFillColorAlpha(color+1,0.35)
                            band2.SetFillColorAlpha(color,0.35)
                            canv.cd(2)
                            tdrDraw(pull, "", ROOT.kFullCircle, color, 1, color, 0, color)
                            tdrDraw(band_pull2, "e3", 9, ROOT.kWhite, 1, color, 1001, color)
                            band_pull2.SetMarkerSize(0)
                            tdrDraw(band_pull, "e3", 9, ROOT.kWhite, 1, color+1, 1001, color+1)
                            band_pull.SetMarkerSize(0)
                            band_pull.SetFillColorAlpha(color+1,0.35)
                            band_pull2.SetFillColorAlpha(color,0.35)
                            line_ =  rt.TLine(p_xmin, 1, p_xmax, 1)
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
                            canv.SaveAs(self.outdir+"BackgroundShape_"+unique_name+"_"+str(self.xmin)+"_"+str(self.xmax)+".pdf")
                        else:
                            canv.SaveAs(self.outdir+"BackgroundShape_"+unique_name+".pdf")
        print ftest_tot

if __name__ == '__main__':
    PlotBkg = FitBackgroundShapes()
    PlotBkg.DoFits()
    # for xmin in range(800,1100+1,100):
    #     for xmax in range(5100,6000+1,100):
    #         if xmin >= xmax: continue
    #         PlotBkg = FitBackgroundShapes(xmin,xmax)
    #         PlotBkg.DoFits()
