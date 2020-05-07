import os, ROOT, glob
from math import sqrt, log
from numpy import arange
from array import array
from tdrstyle_all import *
from Utils import *

ROOT.gInterpreter.ProcessLine('#include "'+os.environ["CMSSW_BASE"]+'/src/UHH2/VHResonances/include/constants.hpp"')
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.gStyle.SetOptStat(0)

PlotFolder = os.environ["CMSSW_BASE"]+"/src/UHH2/VHResonances/Analysis/OtherPlots/"

color = {600:  400,
         800:  416,
         1000: 432,
         1200: 500,
         1400: 616,
         1600: 632,
         1800: 600,
         2000: 616,
         2500: 632,
         3000: 632,
         3500: 632,
         4000: 750,
         4500: 800,
         5000: 820,
         5500: 840,
         6000: 860,
         7000: 880,
         8000: 900,
         }

@timeit
def GetDeepBoostedvsPT():
    # imax = 100
    imax = 1000
    # imax = -1
    canv_PTSelection = tdrCanvas("canv_PTSelection", 200, 5000, 0.01, 7000, "pT","event")
    # canv_PTSelection.SetLogy()

    DeepBoostedvsPT_DY = ROOT.TH2F("DeepBoostedvsPT_DY", "", 101,0,1.01, 50,0,5000)
    DeepBoostedvsPT_DY.SetDirectory(0)
    for f_DY_name in glob.glob("/nfs/dust/cms/user/amalara//WorkingArea/File//Analysis/2016/Selection/Puppi/muonchannel/nominal/workdir_Selection_MC_DY_2016/uhh2.AnalysisModuleRunner.MC.MC_DY_2016_*"):
        print f_DY_name
        f_DY = ROOT.TFile(f_DY_name)
        t_DY = f_DY.Get("AnalysisTree")
        i_ = 0
        for ev in t_DY:
            i_+=1
            if i_>imax and imax>0: break
            for zp in ev.ZprimeCandidate:
                jet = zp.H()
                a = DeepBoostedvsPT_DY.Fill(jet.btag_DeepBoosted_H4qvsQCD(), jet.pt(),ev.weight_GLP)
    canv_DY = tdrCanvas("canv_DY", 0, 1.1, 200, 5000, "DeepBoosted","pT")
    canv_DY.SetRightMargin(0.15)
    canv_DY.SetLogz()
    DeepBoostedvsPT_DY.Draw("colz")
    canv_DY.SaveAs(PlotFolder+"DeepBoostedvsPT_DY.pdf")


    DeepBoostedvsPT = ROOT.TH2F("DeepBoostedvsPT", "", 101,0,1.01, 50,0,5000)
    DeepBoostedvsPT.SetDirectory(0)
    DeepBoosted = ROOT.TH1F("DeepBoosted", "", 101,0,1.01)
    PTSelection = ROOT.TH1F("PTSelection", "", 50,0,5000)
    h_pts = []
    DeepBoostedvsPT_all = {}
    for mass in ROOT.MassPoints:
        massName = str(int(mass))
        DeepBoostedvsPT_all[massName] = ROOT.TH2F("DeepBoostedvsPT"+massName, "", 101,0,1.01, 50,0,5000)
        DeepBoostedvsPT_all[massName].SetDirectory(0)
        print mass
        # if mass <3000: continue
        f = ROOT.TFile("/nfs/dust/cms/user/amalara//WorkingArea/File//Analysis/2016/Selection/Puppi/muonchannel/nominal/workdir_Selection_MC_ZprimeToZH_M"+massName+"_2016/uhh2.AnalysisModuleRunner.MC.MC_ZprimeToZH_M"+massName+"_2016_0.root")
        t = f.Get("AnalysisTree")
        h_pt = ROOT.TH1F("PTSelection"+massName, "; pT; events", 50, 0, 5000)
        h_pt.SetDirectory(0)
        h_pts.append(h_pt)
        i_ = 0
        for ev in t:
            i_+=1
            if i_>imax and imax>0: break
            for zp in ev.ZprimeCandidate:
                jet = zp.H()
                a = DeepBoostedvsPT.Fill(jet.btag_DeepBoosted_H4qvsQCD(), jet.pt(),ev.weight_GLP)
                a = DeepBoostedvsPT_all[massName].Fill(jet.btag_DeepBoosted_H4qvsQCD(), jet.pt(),ev.weight_GLP)
                a = DeepBoosted.Fill(jet.btag_DeepBoosted_H4qvsQCD(), ev.weight_GLP)
                a = PTSelection.Fill(jet.pt(), ev.weight_GLP)
                a = h_pt.Fill(jet.pt(), ev.weight_GLP)
        # h_pt.Draw("same")
        tdrDraw(h_pt, "p same", ROOT.kFullDotLarge, color[mass], 1, color[mass], 0, color[mass])


    PTSelection.Draw("same")
    # tdrDraw(PTSelection, "hist same", ROOT.kFullDotLarge, ROOT.kBlack, 1, ROOT.kBlack, 0, ROOT.kBlack)
    canv_PTSelection.SaveAs(PlotFolder+"PTSelection.pdf")

    canv = tdrCanvas("canv", 0, 1.1, 200, 5000, "DeepBoosted","pT")
    canv.SetRightMargin(0.15)
    canv.SetLogz()
    DeepBoostedvsPT.Draw("colz")
    canv.SaveAs(PlotFolder+"DeepBoostedvsPT.pdf")
    canv_DeepBoosted = tdrCanvas("canv_DeepBoosted", 0, 1.1, 0, 10000, "DeepBoosted","event")
    DeepBoosted.Draw("")
    canv_DeepBoosted.SaveAs(PlotFolder+"DeepBoosted.pdf")
    File_DeepBoosted = ROOT.TFile.Open(PlotFolder+"DeepBoosted.root" ,"RECREATE")
    File_DeepBoosted.cd()
    DeepBoostedvsPT_DY.Write()
    DeepBoostedvsPT.Write()
    for name in DeepBoostedvsPT_all:
        DeepBoostedvsPT_all[name].Write()
    File_DeepBoosted.Close()
    return (DeepBoostedvsPT_DY,DeepBoostedvsPT,DeepBoostedvsPT_all)

@timeit
def SensitivityScan():
    File_DeepBoosted = ROOT.TFile.Open(PlotFolder+"DeepBoosted.root" ,"READ")
    h_bkg = File_DeepBoosted.Get("DeepBoostedvsPT_DY")
    h_sig = File_DeepBoosted.Get("DeepBoostedvsPT")
    h_sig_all = {}
    for mass in ROOT.MassPoints:
        massName = str(int(mass))
        h_sig_all[massName] = File_DeepBoosted.Get("DeepBoostedvsPT"+massName)
    max = -1
    bin_max = -1
    bin_max_y = -1
    for x in range(1,h_sig.GetNbinsX()):
        s = h_sig.Integral(x,h_sig.GetNbinsX(),1,h_sig.GetNbinsY())
        b = h_bkg.Integral(x,h_bkg.GetNbinsX(),1,h_bkg.GetNbinsY())
        if b==0: continue
        # sig = s/sqrt(b);
        sig = sqrt(2*((s+b)*log(1+s/b)-s));
        # return isfinite(sig)? sig : 0;
        if sig>max:
            max = sig
            bin_max = x
        # print "x", x, h_sig.GetXaxis().GetBinCenter(x), "s", s, "b", b, "s%", s/h_sig.Integral(), "b%", b/h_bkg.Integral(), "sig", sig, "max", max
    print bin_max, bin_max_y, h_sig.GetXaxis().GetBinCenter(bin_max), h_sig.GetYaxis().GetBinCenter(bin_max_y), max, h_sig.Integral(x,h_sig.GetNbinsX(),1,h_sig.GetNbinsY())/h_sig.Integral(), h_bkg.Integral(x,h_bkg.GetNbinsX(),1,h_bkg.GetNbinsY())/h_bkg.Integral()
    for name in sorted(h_sig_all.keys()):
        hs = h_sig_all[name]
        max = -1
        bin_max = -1
        bin_max_y = -1
        max_all = -1
        bin_max_all = -1
        bin_max_y_all = -1
        max_1000 = -1
        bin_max_1000 = -1
        bin_max_y_1000 = -1
        for y in range(1,hs.GetNbinsY()):
            for x in range(1,hs.GetNbinsX()):
                s = hs.Integral(x,hs.GetNbinsX(),y,hs.GetNbinsY())
                b = h_bkg.Integral(x,h_bkg.GetNbinsX(),y,h_bkg.GetNbinsY())
                if b==0: continue
                # sig = s/sqrt(b);
                sig = sqrt(2*((s+b)*log(1+s/b)-s));
                # return isfinite(sig)? sig : 0;
                if sig>max:
                    max = sig
                    bin_max =  hs.GetXaxis().GetBinCenter(x)
                    bin_max_y = hs.GetYaxis().GetBinCenter(y)
        for x in range(1,hs.GetNbinsX()):
            s = hs.Integral(x,hs.GetNbinsX(),1,hs.GetNbinsY())
            b = h_bkg.Integral(x,h_bkg.GetNbinsX(),1,h_bkg.GetNbinsY())
            if b==0: continue
            # sig = s/sqrt(b);
            sig = sqrt(2*((s+b)*log(1+s/b)-s));
            # return isfinite(sig)? sig : 0;
            if sig>max_all:
                max_all = sig
                bin_max_all = hs.GetXaxis().GetBinCenter(x)
                bin_max_y_all = hs.GetYaxis().GetBinCenter(1)
        for x in range(1,hs.GetNbinsX()):
            s = hs.Integral(x,hs.GetNbinsX(),hs.GetYaxis().FindBin(1000),hs.GetNbinsY())
            b = h_bkg.Integral(x,h_bkg.GetNbinsX(),h_bkg.GetYaxis().FindBin(1000),h_bkg.GetNbinsY())
            if b==0: continue
            # sig = s/sqrt(b);
            sig = sqrt(2*((s+b)*log(1+s/b)-s));
            # return isfinite(sig)? sig : 0;
            if sig>max_1000:
                max_1000 = sig
                bin_max_1000 = hs.GetXaxis().GetBinCenter(x)
                bin_max_y_1000 = 1000
        print name, "\t", bin_max, "\t", bin_max_y, "\t", round(max,2), "\t", bin_max_1000, "\t", bin_max_y_1000, "\t", round(max_1000,2), "\t", bin_max_all, "\t", bin_max_y_all, "\t", round(max_all,2), "\t",
        s = hs.Integral(hs.GetXaxis().FindBin(0.2),hs.GetNbinsX(),hs.GetYaxis().FindBin(1000),hs.GetNbinsY())
        b = h_bkg.Integral(h_bkg.GetXaxis().FindBin(0.2),h_bkg.GetNbinsX(),h_bkg.GetYaxis().FindBin(1000),h_bkg.GetNbinsY())
        if b==0: continue
        sig_my = sqrt(2*((s+b)*log(1+s/b)-s));
        print sig_my, "\t", s, "\t", b,
        s = hs.Integral(hs.GetXaxis().FindBin(0.02),hs.GetNbinsX(),hs.GetYaxis().FindBin(1000),hs.GetNbinsY())
        b = h_bkg.Integral(h_bkg.GetXaxis().FindBin(0.02),h_bkg.GetNbinsX(),h_bkg.GetYaxis().FindBin(1000),h_bkg.GetNbinsY())
        if b==0: continue
        sig_my = sqrt(2*((s+b)*log(1+s/b)-s));
        print sig_my, "\t", s, "\t", b


@timeit
def Plot2DVar(folder,var):
    canv_Plot2DVar = tdrCanvas("canv_Plot2DVar", 200, 5000, 0.01, 7000, "pT","event")
    canv_Plot2DVar.SetLogz()
    canv_Plot2DVar.SetRightMargin(0.15)

    isFirst = True
    hs = []
    for mass in ROOT.MassPoints:
        massName = str(int(mass))
        print mass,
        # if mass >5000: continue
        f = ROOT.TFile("/nfs/dust/cms/user/amalara//WorkingArea/File//Analysis/2016/Preselection/Puppi/muonchannel/nominal/uhh2.AnalysisModuleRunner.MC.MC_ZprimeToZH_M"+massName+"_2016_noTree.root")
        h = f.Get(folder+"/"+var)
        print h.Integral()
        h.SetDirectory(0)
        hs.append(h)
        if isFirst:
            h2D = h
            isFirst = False
        else:
            h2D.Add(h)


    h2D.Draw("colz")
    canv_Plot2DVar.SaveAs(PlotFolder+"Plot2DVar_"+folder+"_"+var+".pdf")


@timeit
def SensitivityScan2D():
    doFast = True
    doFast = False
    imax = -1
    # imax = 20000
    path = "/nfs/dust/cms/user/amalara//WorkingArea/File//Analysis/2016/Selection/Puppi/muonchannel/nominal/"

    extraText = "_Total" if doFast else ""
    bins = arange(0,1,0.001)
    histo2D = ROOT.TH2D("btag_DeepBoosted_H4qvsQCD","btag_DeepBoosted_H4qvsQCD",ROOT.MassPoints.size(),array('d',list(ROOT.MassPoints)+[9000]), len(bins)-1, array('d',bins))
    # for x in range(histo2D.GetNbinsX()+1):
    #     print x, histo2D.GetXaxis().GetBinCenter(x)
    # for x in range(histo2D.GetNbinsY()):
    #     print x, histo2D.GetYaxis().GetBinCenter(x)


    file_bkg = ROOT.TFile(path+"uhh2.AnalysisModuleRunner.MC.MC_DY_2016_noTree.root")
    h_bkg = file_bkg.Get("ZprimeCandidate_PTMassCut/H_btag_DeepBoosted_H4qvsQCD")
    DeepBoosted_bkg = ROOT.TH2D("DeepBoosted_bkg","DeepBoosted_bkg", 500, 0, 5000, 1000, -0.01,1.01,)
    if not doFast:
        for n_bkg in glob.glob(path+"workdir_Selection_MC_DY_2016/uhh2.AnalysisModuleRunner.MC.MC_DY_2016_*"):
            print n_bkg
            f_bkg = ROOT.TFile(n_bkg)
            t_bkg = f_bkg.Get("AnalysisTree")
            i_ = 0
            for ev in t_bkg:
                i_+=1
                if i_>imax and imax>0: break
                for zp in ev.ZprimeCandidate:
                    a = DeepBoosted_bkg.Fill(zp.Zprime_mass(), zp.H().btag_DeepBoosted_H4qvsQCD(), ev.weight_GLP)

    if not doFast:
        canv2D_bkg = tdrCanvas("canv2D_bkg", 600, 10000, 0.01, 1, "M_{Z'}", "cut")
        canv2D_bkg.SetLogz(1)
        DeepBoosted_bkg.SetContour(200)
        DeepBoosted_bkg.GetZaxis().SetRangeUser(0.001, 300)
        canv2D_bkg.SetRightMargin(0.15)
        DeepBoosted_bkg.Draw("colz")
        canv2D_bkg.SaveAs(PlotFolder+"SensitivityScan2D_DY"+extraText+".pdf")
        canv2D_bkg.SetLogy(1)
        canv2D_bkg.SaveAs(PlotFolder+"SensitivityScan2D_DY"+extraText+"_log.pdf")


    DeepBoosted_sig_all = ROOT.TH2D("DeepBoosted_sig_All","DeepBoosted_sig_All", 1000, 0, 10000, 1000, -0.01,1.01)
    for mass in ROOT.MassPoints:
        massName = str(int(mass))
        sigma = mass*0.02+20.
        print mass, sigma
        f_sig = ROOT.TFile(path+"workdir_Selection_MC_ZprimeToZH_M"+massName+"_2016/uhh2.AnalysisModuleRunner.MC.MC_ZprimeToZH_M"+massName+"_2016_0.root")
        hs = f_sig.Get("ZprimeCandidate_PTMassCut/H_btag_DeepBoosted_H4qvsQCD")
        t_sig = f_sig.Get("AnalysisTree")
        DeepBoosted_sig = ROOT.TH2D("DeepBoosted_sig"+massName,"DeepBoosted_sig"+massName, 1000, 0, 10000, 1000, -0.01,1.01)
        if not doFast:
            i_ = 0
            for ev in t_sig:
                i_+=1
                if i_>imax and imax>0: break
                for zp in ev.ZprimeCandidate:
                    a = DeepBoosted_sig.Fill(zp.Zprime_mass(), zp.H().btag_DeepBoosted_H4qvsQCD(), ev.weight_GLP)
                    a = DeepBoosted_sig_all.Fill(zp.Zprime_mass(), zp.H().btag_DeepBoosted_H4qvsQCD(), ev.weight_GLP)
        if not doFast:
            canv2D_sig = tdrCanvas("canv2D_sig"+massName, 600, 10000, 0.01, 1, "M_{Z'}", "cut")
            canv2D_sig.SetLogz(1)
            DeepBoosted_sig.SetContour(200)
            DeepBoosted_sig.GetZaxis().SetRangeUser(0.001, 300)
            canv2D_sig.SetRightMargin(0.15)
            DeepBoosted_sig.Draw("colz")
            canv2D_sig.SaveAs(PlotFolder+"SensitivityScan2D_M"+massName+extraText+".pdf")
            canv2D_sig.SetLogy(1)
            canv2D_sig.SaveAs(PlotFolder+"SensitivityScan2D_M"+massName+extraText+"_log.pdf")
        for cut in bins:
            if doFast:
                s = hs.Integral(hs.GetXaxis().FindBin(cut),hs.GetNbinsX())
                b = h_bkg.Integral(h_bkg.GetXaxis().FindBin(cut),h_bkg.GetNbinsX())
            else:
                s = DeepBoosted_sig.Integral(DeepBoosted_sig.GetXaxis().FindBin(mass-3*sigma),DeepBoosted_sig.GetXaxis().FindBin(mass+3*sigma), DeepBoosted_sig.GetYaxis().FindBin(cut),DeepBoosted_sig.GetNbinsY())
                b = DeepBoosted_bkg.Integral(DeepBoosted_bkg.GetXaxis().FindBin(mass-3*sigma),DeepBoosted_bkg.GetXaxis().FindBin(mass+3*sigma), DeepBoosted_bkg.GetYaxis().FindBin(cut),DeepBoosted_bkg.GetNbinsY())
            if b==0:
                sig = 0
            else:
                # sig = s/sqrt(b)
                sig = sqrt(2*((s+b)*log(1+s/b)-s))
            # print mass, cut, sig, DeepBoosted_sig.Integral(), DeepBoosted_bkg.Integral()
            histo2D.Fill(mass,cut,round(sig,1))
    if not doFast:
        canv2D_sig_all = tdrCanvas("canv2D_sig_all", 600, 10000, 0.01, 1, "M_{Z'}", "cut")
        canv2D_sig_all.SetLogz(1)
        DeepBoosted_sig_all.SetContour(200)
        DeepBoosted_sig_all.GetZaxis().SetRangeUser(0.001, 300)
        canv2D_sig_all.SetRightMargin(0.15)
        DeepBoosted_sig_all.Draw("colz")
        canv2D_sig_all.SaveAs(PlotFolder+"SensitivityScan2D_Sig_all"+extraText+".pdf")
        canv2D_sig_all.SetLogy(1)
        canv2D_sig_all.SaveAs(PlotFolder+"SensitivityScan2D_Sig_all"+extraText+"_log.pdf")
    canv2D = tdrCanvas("btag_DeepBoosted_H4qvsQCD", 600, 10000, 0.01, 1, "M_{Z'}", "cut")
    histo2D.SetContour(200)
    if not doFast:
        histo2D.GetZaxis().SetRangeUser(50, 300)
    canv2D.SetRightMargin(0.15)
    histo2D.Draw("colz")
    canv2D.SaveAs(PlotFolder+"SensitivityScan2D"+extraText+".pdf")
    canv2D.SetLogy(1)
    canv2D.SaveAs(PlotFolder+"SensitivityScan2D"+extraText+"_log.pdf")
    canv_projx = tdrCanvas("btag_DeepBoosted_H4qvsQCD_projx", 0, 1, 0.01, 100, "M_{Z'}", "cut")
    histo2D.ProjectionX("_px",0,-1,"e").Draw()
    canv_projy = tdrCanvas("btag_DeepBoosted_H4qvsQCD_projy", 0, 1, 0.01, 100, "M_{Z'}", "cut")
    histo2D.ProjectionY("_py",0,-1,"e").Draw()
    canv_projx.SaveAs(PlotFolder+"SensitivityScan2D_projX"+extraText+".pdf")
    canv_projy.SaveAs(PlotFolder+"SensitivityScan2D_projY"+extraText+".pdf")

def main():
    GetDeepBoostedvsPT()
    # SensitivityScan()
    # SensitivityScan2D()
    # Plot2DVar("nTopJet_DeltaRDiLepton","jetptvsdeltaphi_dilep_jet")
    # Plot2DVar("nTopJet_JetDiLeptonPhiAngular","jetptvsdeltaphi_dilep_jet")
    # Plot2DVar("nTopJet_DeltaRDiLepton","jetptvsdeltaphi_dilep_1")
    # Plot2DVar("nTopJet_JetDiLeptonPhiAngular","jetptvsdeltaphi_dilep_1")

if __name__ == '__main__':
    main()
