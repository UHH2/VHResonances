from Utils import *

import tdrstyle_all as TDR
TDR.writeExtraText = True
TDR.extraText  = "Simulation"
TDR.extraText2 = "Work in progress"

'''
Module to study optimal pt-dependent cut

- Need the Selection output as input
- The StoreVars function is time consuming (~20 min for everything). Run it only once
- Default collection is Puppi.
- Efficiency as function of DeltaR or pt. Choose option
- ID_kincut as reference to calculate efficiency

'''
# TODO invisible channel not fully implemented yet

import numpy as np
import pandas as pd
from math import sqrt, log, pow
from array import array

def PtToMass(pt):
    return pt*2

def MassToPt(mass):
    return mass/2

def PtToMass2(pt):
    return (pt-22.5)/0.36

def MassToPt2(mass):
    return mass*0.36 +22.5

def PtDependentCut(pt):
    return 5.08*1e-03+2.72*1e07*pow(PtToMass(pt),-3)

def PtDependentCut2(pt):
    return 5.08*1e-03+2.72*1e07*pow(PtToMass2(pt),-3)

def MassDependentCut(mass):
    return 5.08*1e-03+2.72*1e07*pow(mass,-3)


colors = {"2016":       ROOT.kGreen+1,
          "2017":       ROOT.kRed+1,
          "2018":       ROOT.kOrange+1,
          "all":        ROOT.kBlue+1,
          "lepton":     ROOT.kFullCircle,
          "muon":       ROOT.kFullTriangleDown,
          "electron":   ROOT.kFullTriangleUp,
}

class TaggerCutStudy(VariablesBase):
    def __init__(self):
        VariablesBase.__init__(self)
        TDR.lumi_13TeV  = str(round(float(self.lumi_map["RunII"]["lumi_fb"]),1))+" fb^{-1}"
        self.outdir = self.Path_ANALYSIS+"Analysis/OtherPlots/TaggerCutStudy/"
        self.fName = "Variables"
        self.Samples = filter(lambda x: self.MainBkg in x or self.Signal in x, self.Processes_Year_Dict["2016"])
        self.Ybins = sorted([x for x in np.arange(0.1,1,0.005)]+[x for x in np.arange(0.001,0.1,0.0005)])
        self.Xbins = sorted([self.MassPoints[i] - (self.MassPoints[i]-self.MassPoints[i-1])/2 for i in range(1,len(self.MassPoints))]+[500,8500])

        os.system("mkdir -p "+self.outdir)

    @timeit # circa 40 minutes
    def StoreVars(self):
        vars = {}
        for year, channel, sample in list(itertools.product(self.years, self.Channels, self.Samples)):
            sample = sample.replace("2016",year)
            if DoControl([""], year+channel+sample, channel, sample): continue
            if "invisible" in channel: continue #TODO
            print year,channel,sample
            for filename in glob.glob(self.Path_STORAGE+year+"/Selection/Puppi/"+channel+"channel/nominal/workdir_Selection_"+sample+"/*.root"):
                f_ = ROOT.TFile(filename)
                t_ = f_.Get("AnalysisTree")
                for ev in t_:
                    if ev.ZprimeCandidate.size()!=1 :
                        print "Unexpected number of ZprimeCandidate."
                        continue
                    for zp in ev.ZprimeCandidate:
                        jet = zp.H()
                        vars.setdefault("year",[]).append(year)
                        vars.setdefault("channel",[]).append(channel)
                        vars.setdefault("sample",[]).append(sample)
                        vars.setdefault("weight_GLP",[]).append(ev.weight_GLP)
                        vars.setdefault("jet_pt",[]).append(jet.pt())
                        vars.setdefault("jet_DB",[]).append(jet.btag_DeepBoosted_H4qvsQCD())
                        vars.setdefault("Zprime_mass",[]).append(zp.Zprime_mass())
        df = pd.DataFrame(data=vars)
        df.to_pickle(self.outdir+self.fName+".pkl")


    def LoadVars(self):
        self.df = pd.read_pickle(self.outdir+self.fName+".pkl")

    def LoadHistos(self, var):
        file_ = ROOT.TFile(self.outdir+"Histos.root", "OPEN")
        self.histos = {}
        for year, channel in list(itertools.product(self.years+["all"], self.Channels+["lepton"])):
            if "invisible" in channel: continue #TODO
            histoname = var+year+channel+self.Signal
            histoname = histoname.replace("all","").replace("lepton","")
            self.histos[year+channel+self.Signal] = file_.Get(histoname)
            if year == "all": histoname = histoname.replace(self.Signal,self.MainBkg)
            else: histoname = histoname.replace(self.Signal,self.MainBkg+"_"+year)
            self.histos[year+channel+self.MainBkg] = file_.Get(histoname)
            if not self.histos[year+channel+self.Signal]: continue
            self.histos[year+channel+self.Signal].Scale(0.1)
            self.histos[year+channel+self.Signal].RebinX(30)
            self.histos[year+channel+self.MainBkg].RebinX(30)
            self.histos[year+channel+self.Signal].SetDirectory(0)
            self.histos[year+channel+self.MainBkg].SetDirectory(0)
        file_.Close()

    @timeit # circa 20 minutes
    def PlotVar(self):
        self.LoadVars()
        def CreateHisto(var,name):
            cut_min, cut_max, cut_bins = (0.,1.,2000)
            pt_min, pt_max, pt_bins = (0,5000,1000)
            mass_min, mass_max, mass_bins = (0,10000,1000)
            if var == "PT":   hist = ROOT.TH2D(name, name, pt_bins, pt_min, pt_max, cut_bins, cut_min, cut_max)
            if var == "Mass": hist = ROOT.TH2D(name, name, mass_bins, mass_min, mass_max, cut_bins, cut_min, cut_max)
            hist.SetDirectory(0)
            return hist

        histos = {}
        for var in ["PT", "Mass"]:
            for s_ in [self.MainBkg, self.Signal]:
                name = var+s_
                histos[name] = CreateHisto(var,name)
            for year, channel in list(itertools.product(self.years, self.Channels)):
                if "invisible" in channel: continue #TODO
                name = var+year+channel+self.Signal
                histos[name] = CreateHisto(var,name)

        for year, channel, sample in list(itertools.product(self.years, self.Channels, self.Samples)):
            sample = sample.replace("2016",year)
            if DoControl([""], year+channel+sample, channel, sample): continue
            if "invisible" in channel: continue #TODO
            print year, channel, sample

            for var in ["PT", "Mass"]:
                name = var+year+channel+sample
                histos[name] = CreateHisto(var,name)
                histos[name].SetDirectory(0)

            for ind, entry in self.df[(self.df["year"]==year) & (self.df["channel"]==channel) & (self.df["sample"]==sample)].iterrows():
                pt, mass, DB, w = (entry["jet_pt"],entry["Zprime_mass"],entry["jet_DB"],entry["weight_GLP"])
                for var in ["PT", "Mass"]:
                    x_val = pt if var=="PT" else mass
                    histos[var+year+channel+sample].Fill(x_val, DB, w)
                    for s_ in [self.MainBkg, self.Signal]:
                        if s_ in sample: histos[var+s_].Fill(x_val, DB, w)
                    if self.Signal in sample:
                        histos[var+year+channel+self.Signal].Fill(x_val, DB, w)

        file_ = ROOT.TFile(self.outdir+"Histos.root", "RECREATE")
        for name, histo in histos.items():
            if "PT" in name:   canv = tdrCanvas("canv"+name, 200, 5000,  0, 1.1, "p_T (GeV)", "DeepBoosted")
            if "Mass" in name: canv = tdrCanvas("canv"+name, 500, 10000, 0, 1.1, "M (GeV)",   "DeepBoosted")
            canv.SetRightMargin(0.15)
            canv.SetLogy()
            canv.SetLogz()
            histo.Draw("colz")
            canv.SaveAs(self.outdir+name+".pdf")
            histo.Write(name)
        file_.Close()

    @timeit
    def SensitivityScan(self, var="PT"):
        self.LoadVars()
        self.LoadHistos(var)

        significance = {}
        pt_max       = {}
        DB_cut       = {}
        err          = {}
        err2         = {}
        histo = ROOT.TH2D("SensitivityScan","SensitivityScan",len(self.Xbins)-1,array('d',self.Xbins), len(self.Ybins)-1, array('d',self.Ybins))
        x_min = 0
        x_max = 3000 if var=="PT" else 8500
        canv = tdrCanvas("canv", x_min, x_max, 0.0001, 10, var,"DeepBoosted")
        canv.SetLogy()
        canv.SetTickx(1)
        canv.SetTicky(1)
        leg = tdrLeg(0.40,0.70,0.95,0.89, 0.025, 42, ROOT.kBlack)
        leg.SetNColumns(3)
        dic_gr = {}
        for year, channel in list(itertools.product(self.years+["all"], self.Channels+["lepton"])):
            if "invisible" in channel: continue #TODO
            print year, channel
            h_sig = self.histos[year+channel+self.Signal]
            h_bkg = self.histos[year+channel+self.MainBkg]
            if not h_sig or not h_bkg: continue
            for xbin in range(1,h_sig.GetNbinsX()+1):
                mass = h_sig.GetXaxis().GetBinCenter(xbin)
                mean = MassToPt(mass) if var=="PT" else mass
                sigma = 1*(mean*0.02+20.) # TODO taken from fits
                xbin_lo = h_sig.GetXaxis().FindBin(mean-sigma)
                xbin_hi = h_sig.GetXaxis().FindBin(mean+sigma)
                y_var = "jet_pt" if var=="PT" else "Zprime_mass"
                name = year+channel+str(mass)
                significance.setdefault(name,[])
                pt_max.setdefault(name,[])
                DB_cut.setdefault(name,[])
                err.setdefault(name,[])
                err2.setdefault(name,[])
                for ybin in range(1,h_sig.GetNbinsY()+1):
                    y = h_sig.GetYaxis().GetBinCenter(ybin)
                    s = h_sig.Integral(xbin_lo,xbin_hi,ybin,h_sig.GetNbinsY())
                    b = h_bkg.Integral(xbin_lo,xbin_hi,ybin,h_sig.GetNbinsY())
                    sig = sqrt(2*((s+b)*log(1+s/b)-s)) if b>0 else sqrt(s)
                    significance[name].append(sig)
                    pt_max[name].append(mean)
                    DB_cut[name].append(y)
                    if year=="all" and channel=="lepton":
                        histo.SetBinContent(histo.GetXaxis().FindBin(mass),histo.GetYaxis().FindBin(y),sig)
                significance[name] = np.array(significance[name])
                pt_max[name] = np.array(pt_max[name])
                DB_cut[name] = np.array(DB_cut[name])
                index = np.abs(significance[name] - significance[name].max()).argmin()
                index_var = np.abs(significance[name] - significance[name].max()+significance[name].std()).argmin()
                err[name] = abs(DB_cut[name][index]-DB_cut[name][index_var])/sqrt(12)
                significance[name] = significance[name][index]
                pt_max[name] = pt_max[name][index]
                DB_cut[name] = DB_cut[name][index]
                if year!="all" and channel!="lepton":
                    err2[name].append(DB_cut[name])

            filtered_pt_max = []
            filtered_DB_cut = []
            filtered_err = []
            all_pt_max = []
            all_DB_cut = []
            all_err = []
            for xbin in range(1,h_sig.GetNbinsX()+1):
                mass = h_sig.GetXaxis().GetBinCenter(xbin)
                name = year+channel+str(mass)
                all_err.append(err[name])
                all_pt_max.append(pt_max[name])
                all_DB_cut.append(DB_cut[name])
                if mass<600 or mass>6000: continue
                filtered_pt_max.append(pt_max[name])
                filtered_DB_cut.append(DB_cut[name])
                if year=="all" and channel=="lepton":
                    values = [err2[x][0] for x in err2 if str(mass) in x and err2[x]]
                    values = list(filter(lambda x: x-np.mean(values)<2*np.std(values) and x>1e-04, values))
                    filtered_err.append((np.max(values)-np.min(values))/2)
            gr = ROOT.TGraph(len(all_DB_cut), array('d',all_pt_max), array('d',all_DB_cut))
            # gr = ROOT.TGraphErrors(len(all_DB_cut), array('d',all_pt_max), array('d',all_DB_cut), array('d',np.zeros(len(all_DB_cut))), array('d',all_err));
            dic_gr[year+channel] = gr
            # tdrDraw(gr, "P",  colors[channel], colors[year], 2, colors[year], 1000, colors[year])

            if year=="all" and channel=="lepton":
                # gr_forfit = ROOT.TGraph(len(filtered_DB_cut), array('d',filtered_pt_max), array('d',filtered_DB_cut))
                gr_forfit = ROOT.TGraphErrors(len(filtered_DB_cut), array('d',filtered_pt_max), array('d',filtered_DB_cut), array('d',np.zeros(len(filtered_DB_cut))), array('d',filtered_err))
                dic_gr[year+channel+"fit"] = gr_forfit
                func_ref = ROOT.TF1("1/x3","[0]+[1]*TMath::Power(x,-3)",x_min, x_max)
                funcs = {}
                funcs["1/x3"] = ROOT.TF1("1/x3","[0]+[1]*TMath::Power(x,-3)",x_min, x_max)
                funcs["1/x3"].SetLineColor(ROOT.kBlue+1)
                funcs["1/x2x3"] = ROOT.TF1("1/x2x3","[0]+[1]*TMath::Power(x,-2)+[2]*TMath::Power(x,-3)",x_min, x_max)
                funcs["1/x2x3"].SetLineColor(ROOT.kOrange+1)
                funcs["1/x1x3"] = ROOT.TF1("1/x1x3","[0]+[1]*TMath::Power(x,-1)+[2]*TMath::Power(x,-3)",x_min, x_max)
                funcs["1/x1x3"].SetLineColor(ROOT.kViolet+1)
                if var=="Mass":
                    func_ref = ROOT.TF1("1/x3","[0]+[1]*TMath::Power(x/2,-3)",600,8000)
                func_ref.SetParameters(5.08*1e-03,2.72*1e07)
                ROOT.gStyle.SetOptFit(0)
                for func in funcs:
                    gr_forfit.Fit(funcs[func],"RQ")
                    print "#"*50, "\n", func, "\n", funcs[func].GetParameter(0), "\n", funcs[func].GetParameter(1), "\n", funcs[func].GetParameter(2), "\n", "#"*50
                    funcs[func].Draw("same")
                tdrDraw(gr_forfit, "P",  colors[channel], ROOT.kBlack, 2, ROOT.kBlack, 1000, ROOT.kBlack)
                # tdrDraw(gr_forfit, "P",  colors[channel], colors[year], 2, colors[year], 1000, colors[year])
                func_ref.Draw("same")
            # leg.AddEntry(gr, year+"_"+channel, "lp")
        canv.SaveAs(self.outdir+"YVs"+var+".pdf")
        canv = tdrCanvas("canv_Significance", 0, 1.1, x_min, x_max, "DeepBoosted","pT")
        canv.SetRightMargin(0.15)
        canv.SetLogy()
        histo.Draw("colz")
        canv.SaveAs(self.outdir+"Significance"+var+".pdf")
        canv.SaveAs(self.outdir+"Significance"+var+".root")


def main():
    TCS = TaggerCutStudy()
    # TCS.StoreVars()
    # TCS.PlotVar()
    TCS.SensitivityScan("Mass")
    TCS.SensitivityScan("PT")

if __name__ == '__main__':
    main()
