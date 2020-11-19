from Utils import *

import tdrstyle_all as TDR
TDR.writeExtraText = True
TDR.extraText  = "Simulation"
TDR.extraText2 = "Work in progress"

ForThesis(TDR)

'''
Module to study optimal pt-dependent cut

- Need the ExtractTaggerInformation module output as input
- Producing all histograms takes a lot of time (~1h) and memory (~35GB)

'''

import numpy as np
import pandas as pd
from math import sqrt, log, pow
from array import array

def PtToMass(pt):
    return pt*2

def PtDependentCut(pt):
    return 5.03*1e-03+1.70*1e07*pow(PtToMass(pt),-3)

def MassDependentCut(mass):
    return 5.08*1e-03+2.72*1e07*pow(mass,-3)


colors = {"2016":       rt.kGreen+1,
          "2017":       rt.kRed+1,
          "2018":       rt.kOrange+1,
          "RunII":      rt.kBlue+1,
          "lepton":     rt.kFullCircle,
          "muon":       rt.kFullTriangleDown,
          "electron":   rt.kFullTriangleUp,
          "invisible":  rt.kFullSquare,
}

class TaggerCutStudy(VariablesBase):
    def __init__(self):
        VariablesBase.__init__(self)
        TDR.lumi_13TeV  = str(round(float(self.lumi_map["RunII"]["lumi_fb"]),1))+" fb^{-1}"
        self.hfName = "Histos"
        self.isFast = True
        self.isFast = False
        self.fName = "TaggerVariables"
        if not self.isFast: self.fName += "_all"
        self.Samples = filter(lambda x: self.MainBkg in x or self.Signal in x, self.Processes_Year_Dict["2016"])
        self.Ybins = sorted([x for x in np.arange(0,1,0.001)])
        self.Xbins = sorted([self.MassPoints[i] - (self.MassPoints[i]-self.MassPoints[i-1])/2 for i in range(1,len(self.MassPoints))]+[500,8500])
        print self.Xbins
        self.Xvars = ["Mass"]
        self.Yvars = ["H4q", "Hcc", "H4q_SR", "Hcc_SR"]

        self.indir  = self.Path_ANALYSIS+"Analysis/OtherPlots/TaggerInfo/"
        self.outdir = self.Path_ANALYSIS+"Analysis/OtherPlots/TaggerCutStudy/"
        os.system("mkdir -p "+self.outdir)

    def LoadVars(self):
        self.df = joblib.load(self.indir+self.fName+".pkl")
        self.df["jet_HccvsQCD"] = self.df["jet_Hcc"]/(self.df["jet_Hcc"]+self.df.loc[:,["jet_QCDb","jet_QCDbb", "jet_QCDc", "jet_QCDcc", "jet_QCDqq"]].sum(axis=1))

    @timeit
    def CreateHistos(self):
        self.LoadVars()

        histos = {}
        for Xvar, Yvar, year, channel, sample in list(itertools.product(self.Xvars, self.Yvars, self.years+["RunII"], self.Channels+(["lepton"] if self.Channels!=["invisible"] else []), self.Samples+[self.Signal, self.Signal+"_inv"])):
            sample = sample.replace("2016",year)
            if DoControl([""], year+channel+sample, channel, sample): continue
            name = Xvar+Yvar+year+channel+sample
            histos[name] = rt.TH2D(name,name,len(self.Xbins)-1,array('d',self.Xbins), len(self.Ybins)-1, array('d',self.Ybins))
            histos[name+"bins"] = rt.TH2D(name+"bins", name+"bins", 1000, 0, 10000, 2000, 0, 1)
            histos[name].SetDirectory(0)
            histos[name+"bins"].SetDirectory(0)

        for year, channel, sample in list(itertools.product(self.years, self.Channels, self.Samples)):
            sample = sample.replace("2016",year)
            if DoControl([""], year+channel+sample, channel, sample): continue
            print year, channel, sample

            for ind, entry in self.df[(self.df["year"]==year) & (self.df["channel"]==channel) & (self.df["sample"]==sample)].iterrows():
                sample, pt, mass, DB_H4q, DB_Hcc, w = (entry["sample"], entry["jet_pt"],entry["Zprime_mass"],entry["jet_H4qvsQCD"],entry["jet_HccvsQCD"],entry["weight_GLP"])
                for Xvar in self.Xvars:
                    x_val = mass
                    for Yvar in self.Yvars:
                        y_val = DB_H4q if "H4q" in Yvar else DB_Hcc
                        if "SR" in Yvar and DB_H4q<PtDependentCut(pt): continue
                        for y_ in [year, "RunII"]:
                            for c_ in [channel, "lepton"]:
                                signal = self.Signal
                                # Change signal name for the invisible channel
                                if "invisible" in channel: signal =  signal + "_inv"
                                if "invisible" in channel and c_=="lepton":
                                    continue
                                for s_ in [sample, signal]:
                                    if not signal in sample and s_==signal: continue
                                    histos[Xvar+Yvar+y_+c_+s_.replace(year,y_)].Fill(x_val, y_val, w)
                                    histos[Xvar+Yvar+y_+c_+s_.replace(year,y_)+"bins"].Fill(x_val, y_val, w)

        file_ = rt.TFile(self.outdir+self.hfName+".root", "RECREATE")
        for name, histo in histos.items():
            canv = tdrCanvas("canv"+name, 500, 10000, 0, 1.1, "M (GeV)",   "DeepBoosted")
            canv.SetRightMargin(0.15)
            canv.SetLogy()
            canv.SetLogz()
            histo.Draw("colz")
            canv.SaveAs(self.outdir+name+".pdf")
            histo.Write(name)
        file_.Close()

    def LoadHistos(self, Xvar="Mass"):
        file_ = rt.TFile(self.outdir+self.hfName+".root", "OPEN")
        self.histos = {}
        for Yvar, year, channel, sample in list(itertools.product(self.Yvars, self.years+["RunII"], self.Channels+(["lepton"] if self.Channels!=["invisible"] else []), self.Samples+[self.Signal, self.Signal+"_inv"])):
            sample = sample.replace("2016",year)
            if DoControl([""], year+channel+sample, channel, sample): continue
            name = Xvar+Yvar+year+channel+sample
            self.histos[name] = file_.Get(name)
            # self.histos[name] = file_.Get(name+"bins")
            self.histos[name].SetDirectory(0)
            # self.histos[name].RebinX(100)
            if self.Signal in sample: self.histos[name].Scale(0.1) #TODO
        file_.Close()

    @timeit
    def SensitivityScan(self, Xvar="Mass", DB="H4q"):
        self.LoadHistos(Xvar)
        significance = {}
        pt_max       = {}
        DB_cut       = {}
        err          = {}
        err2         = {}
        histo = rt.TH2D("SensitivityScan","SensitivityScan",len(self.Xbins)-1,array('d',self.Xbins), len(self.Ybins)-1, array('d',self.Ybins))
        x_min = 0
        x_max = 6500
        canv = tdrCanvas("canv", x_min, x_max, 0.0001, 10, Xvar, DB+"vsQCD")
        canv.SetLogy()
        leg = tdrLeg(0.40,0.70,0.95,0.89, 0.025, 42, rt.kBlack)
        # leg.SetNColumns(3)
        dic_gr = {}
        for year, channel in list(itertools.product(self.years+["RunII"], self.Channels +(["lepton"] if self.Channels!=["invisible"] else []))):
        # for year, channel in list(itertools.product(["RunII"],["lepton"])):

            h_sig = self.histos[Xvar+DB+year+channel+self.Signal+("_inv" if channel=="invisible" else "")]
            h_bkg = self.histos[Xvar+DB+year+channel+self.MainBkg+"_"+year]
            n_yBins = h_bkg.GetNbinsY()

            for xbin in range(1,h_bkg.GetNbinsX()+1):
                mean = h_sig.GetXaxis().GetBinCenter(xbin)
                sigma = 1*(mean*0.02+20.) # 1-2-3 sigma is rather stable
                xbin_lo = h_sig.GetXaxis().FindBin(mean-sigma)
                xbin_hi = h_sig.GetXaxis().FindBin(mean+sigma)
                name = year+channel+str(mean)
                significance.setdefault(name,[])
                pt_max.setdefault(name,[])
                DB_cut.setdefault(name,[])
                err.setdefault(name,[])
                err2.setdefault(name,[])
                for ybin in range(1,n_yBins+1):
                    y = h_bkg.GetYaxis().GetBinCenter(ybin)
                    s = h_sig.Integral(xbin_lo,xbin_hi,ybin,n_yBins)
                    b = h_bkg.Integral(xbin_lo,xbin_hi,ybin,n_yBins)
                    sig = sqrt(2*((s+b)*log(1+s/b)-s)) if b>0 else sqrt(s)
                    significance[name].append(sig)
                    pt_max[name].append(mean)
                    DB_cut[name].append(y)
                    if year=="RunII" and (channel=="lepton" or channel=="invisible"):
                        histo.SetBinContent(histo.GetXaxis().FindBin(mean),histo.GetYaxis().FindBin(y),sig)
                significance[name] = np.array(significance[name])
                pt_max[name] = np.array(pt_max[name])
                DB_cut[name] = np.array(DB_cut[name])
                index = np.abs(significance[name] - significance[name].max()).argmin()
                index_var = np.abs(significance[name] - significance[name].max()+significance[name].std()).argmin()
                err[name] = abs(DB_cut[name][index]-DB_cut[name][index_var])/sqrt(12)
                significance[name] = significance[name][index]
                pt_max[name] = pt_max[name][index]
                DB_cut[name] = DB_cut[name][index]
                if year!="RunII" and channel!="lepton":
                    err2[name].append(DB_cut[name])

            filtered_pt_max = []
            filtered_DB_cut = []
            filtered_err = []
            all_pt_max = []
            all_DB_cut = []
            all_err = []
            for xbin in range(1,h_sig.GetNbinsX()+1):
                mean = h_sig.GetXaxis().GetBinCenter(xbin)
                name = year+channel+str(mean)
                all_err.append(err[name])
                all_pt_max.append(pt_max[name])
                all_DB_cut.append(DB_cut[name])
                if mean<700 or mean>6500: continue
                if channel=="invisible":
                    if mean==5500: continue
                    if mean==3500: continue
                filtered_pt_max.append(pt_max[name])
                filtered_DB_cut.append(DB_cut[name])
                # filtered_err.append(err[name])
                if year=="RunII" and (channel=="lepton" or channel=="invisible"):
                    values = [err2[x][0] for x in err2 if str(mean) in x and len(err2[x])==1]
                    if np.std(values)!=0:
                        values = list(filter(lambda x: x-np.mean(values)<2*np.std(values) and x>1e-04, values))
                        filtered_err.append((np.max(values)-np.min(values))/2)
                    else: filtered_err.append(values[0])
            gr = rt.TGraph(len(all_DB_cut), array('d',all_pt_max), array('d',all_DB_cut))
            # gr = rt.TGraphErrors(len(all_DB_cut), array('d',all_pt_max), array('d',all_DB_cut), array('d',np.zeros(len(all_DB_cut))), array('d',all_err));
            dic_gr[year+channel] = gr
            # tdrDraw(gr, "P",  colors[channel], colors[year], 2, colors[year], 1000, colors[year])

            if year=="RunII" and (channel=="lepton" or channel=="invisible"):
                # gr_forfit = rt.TGraph(len(filtered_DB_cut), array('d',filtered_pt_max), array('d',filtered_DB_cut))
                gr_forfit = rt.TGraphErrors(len(filtered_DB_cut), array('d',filtered_pt_max), array('d',filtered_DB_cut), array('d',np.zeros(len(filtered_DB_cut))), array('d',filtered_err))
                dic_gr[year+channel+"fit"] = gr_forfit
                func_ref = rt.TF1("1/x3","[0]+[1]*TMath::Power(x,-3)",x_min, x_max)
                func_ref2 = rt.TF1("1/x3_2","[0]+[1]*TMath::Power(x,-3)",x_min, x_max)
                func_ref3 = rt.TF1("1/x3_3","[0]+[1]*TMath::Power(x,-3)",x_min, x_max)
                funcs = {}
                funcs["1/x3"] = rt.TF1("1/x3","[0]+[1]*TMath::Power(x,-3)",x_min, x_max)
                funcs["1/x3"].SetLineColor(rt.kBlue+1)
                funcs["1/x2x3"] = rt.TF1("1/x2x3","[0]+[1]*TMath::Power(x,-2)+[2]*TMath::Power(x,-3)",x_min, x_max)
                funcs["1/x2x3"].SetLineColor(rt.kOrange+1)
                funcs["1/x1x3"] = rt.TF1("1/x1x3","[0]+[1]*TMath::Power(x,-1)+[2]*TMath::Power(x,-3)",x_min, x_max)
                funcs["1/x1x3"].SetLineColor(rt.kRed+1)
                if Xvar=="Mass":
                    func_ref = rt.TF1("1/x3","[0]+[1]*TMath::Power(x/2,-3)",600,8000)
                func_ref.SetParameters(5.08*1e-03,2.72*1e07)
                func_ref2.SetParameters(2.52*1e-03, 1.39*1e07)
                func_ref3.SetParameters(0.0032618254201,63499451.4601)
                rt.gStyle.SetOptFit(0)
                for func in funcs:
                    gr_forfit.Fit(funcs[func],"RQ")
                    funcs[func].Draw("same")
                    name = ""
                    if func=="1/x3":   name = "[0]+[1]*x^{-3}"
                    if func=="1/x2x3": name = "[0]+[1]*x^{-2}+[2]*x^{-3}"
                    if func=="1/x1x3": name = "[0]+[1]*x^{-1}+[2]*x^{-3}"
                    print "#"*50, "\n", func
                    for n_ in range(0,funcs[func].GetNpar()):
                        print n_, "{:.2e}".format(funcs[func].GetParameter(n_))
                        name = name.replace("["+str(n_)+"]", "{:.2e}".format(funcs[func].GetParameter(n_)))
                    print "#"*50
                    leg.AddEntry(funcs[func], name, "l")
                tdrDraw(gr_forfit, "P",  colors[channel], rt.kBlack, 2, rt.kBlack, 1000, rt.kBlack)
                # tdrDraw(gr_forfit, "P",  colors[channel], colors[year], 2, colors[year], 1000, colors[year])
                func_ref.SetLineColor(rt.kRed+1)
                func_ref2.SetLineColor(rt.kGreen+2)
                func_ref3.SetLineColor(rt.kOrange-1)
                # func_ref.Draw("same")
                # func_ref2.Draw("same")
                # func_ref3.Draw("same")
            # leg.AddEntry(gr, year+"_"+channel, "lp")
        canv.SaveAs(self.outdir+DB+"Vs"+Xvar+".pdf")
        canv = tdrCanvas("canv_Significance", 0, 1.1, x_min, x_max, "DeepBoosted","pT")
        canv.SetRightMargin(0.15)
        canv.SetLogy()
        histo.Draw("colz")
        canv.SaveAs(self.outdir+"Significance"+Xvar+".pdf")
        canv.SaveAs(self.outdir+"Significance"+Xvar+".root")


def main():
    TCS = TaggerCutStudy()
    # TCS.CreateHistos()
    TCS.SensitivityScan()
    TCS.SensitivityScan(DB="Hcc")

if __name__ == '__main__':
    main()
