from tdrstyle_all import *
from array import array

import tdrstyle_all as TDR
TDR.writeExtraText = True
ForThesis(TDR)

def FutureColliders():
    rt.gROOT.SetBatch(rt.kTRUE)
    rt.gStyle.SetOptStat(0)
    rt.gStyle.SetOptFit(0)


    canv = tdrCanvas("FutureColliders", 50, 5000, 0.5, 1000, "#sqrt{s} [GeV]", "Luminosity [10^{34} cm^{-2}s^{-1} ]", kSquare, iPeriod=-1)
    canv.SetLogx(True)
    canv.SetLogy(True)

    leg = tdrLeg(0.60, 0.70, 0.95, 0.89, 0.035, 42, rt.kBlack)
    #leg.SetNColumns(2)
    graphs = {}
    colors = {}
    graphs["CLIC"] = rt.TGraph(3, array('d',[380,1500,3000]), array('d',[1.5,3.7,5.9]))
    colors["CLIC"] = rt.kOrange+1
    graphs["ILC"]  = rt.TGraph(3, array('d',[250,350,500]), array('d',[0.75,1.0,1.8]))
    colors["ILC"]  = rt.kGreen+2
    graphs["CEPC (2IPs)"]  = rt.TGraph(3, array('d',[91,160,240]), array('d',[64, 20, 6]))
    colors["CEPC (2IPs)"]  = rt.kRed+1
    graphs["FCC-cc (2IPs)"]  = rt.TGraph(5, array('d',[91, 160, 240, 350, 365]), array('d',[450,56,17,3.8,3.1]))
    colors["FCC-cc (2IPs)"]  = rt.kAzure+7

    txts = {}
    txts["Z"] = rt.TLatex(95, 500, "Z (91 GeV)")
    txts["W"] = rt.TLatex(165, 60, "WW (160 GeV)")
    txts["H"] = rt.TLatex(260, 18, "ZH (240 GeV)")
    txts["t"] = rt.TLatex(370, 4, "t#bar{t} (350 GeV)")

    for n,t in txts.items():
        t.SetTextFont(42)
        t.SetTextSize(0.035)
        #t.SetTextColor(rt.kRed+1)
        t.SetTextColor(rt.kBlack)
        t.Draw("same")

    for n,h in graphs.items():
        h.SetLineWidth(2)
        tdrDraw(h, "LP", mcolor=colors[n], lcolor=colors[n])
        leg.AddEntry(h, n, "lp")
    leg.Draw("same")
    canv.SaveAs("FutureColliders.pdf")


if __name__ == '__main__':
    FutureColliders()
