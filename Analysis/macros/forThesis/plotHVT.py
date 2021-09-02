from tdrstyle_all import *
from array import array
import numpy as np
import json

import tdrstyle_all as TDR
TDR.lumi_13TeV = ""
TDR.lumi_sqrtS = ""
ForThesis(TDR)

def plotHVT():
    rt.gROOT.SetBatch(rt.kTRUE)
    rt.gStyle.SetOptStat(0)
    rt.gStyle.SetOptFit(0)

    with open('Gammas_HVT.txt', 'r') as infile:
        gammas = json.load(infile)

        graphs = {}
        colors = {"1": rt.kRed+1,
                  "2": rt.kOrange+1,
                  "3": rt.kGreen+2,
                  "5": rt.kAzure,
                  "8": rt.kViolet+2,
                  "ee":  rt.kOrange+1,
                  "bb":  rt.kRed+1,
                  "ZH":  rt.kGreen+3,
                  "WW":  rt.kAzure+7,
                  "eve": rt.kOrange+1,
                  "tb":  rt.kRed+1,
                  "HW":  rt.kGreen+3,
                  "WZ":  rt.kAzure+7,
                  "tot": rt.kBlack,
                  }

        for particle in list(gammas.keys())+["Vprime"]:
            canv = tdrCanvas("Gamma"+particle, 800, 8300, 0.5, 3000, "M_{"+particle.replace("prime","'")+"} [GeV]", "#Gamma("+particle.replace("prime","'")+"#rightarrow XY) [GeV]", False, iPeriod=-1)
            canv.SetLogy(True)
            # legA = tdrLeg(0.50, 0.18, 0.70, 0.38, 0.035, 42, rt.kBlack)
            # legB = tdrLeg(0.65, 0.18, 0.85, 0.38, 0.035, 42, rt.kBlack)
            # legC = tdrLeg(0.80, 0.18, 0.95, 0.38, 0.035, 42, rt.kBlack)
            # tdrHeader(legA, "Model A", textSize=0.035)
            # tdrHeader(legB, "Model B", textSize=0.035)
            # tdrHeader(legC, "Model C", textSize=0.035)

            leg_gv    = tdrLeg(0.60, 0.15, 0.80, 0.38, 0.045, 42, rt.kBlack)
            leg_model = tdrLeg(0.75, 0.225, 0.95, 0.385, 0.045, 42, rt.kBlack)
            tdrHeader(leg_model, "Scenario")
            tdrHeader(leg_gv, "Coupling")

            isV = particle=="Vprime"
            if isV: particle = "Zprime"
            for model in sorted(gammas[particle]):
                model = str(model)
		if "1" in model: continue
                isA = "A" in model
                isB = "B" in model
                isC = "C" in model
                for gv in sorted(gammas[particle][model]):
                    Masses = np.array([float(x) for x in gammas[particle][model][gv]["mass"]])
                    Gammas = np.array([float(x) for x in gammas[particle][model][gv]["tot"]])
                    gv = str(gv)
                    mask = Masses>1000
                    if gv=="5": mask = Masses>1400
                    if gv=="8": mask = Masses>=2200
                    print particle, model, gv, Gammas[Masses==4000]
                    graphs[particle+model+gv] = rt.TGraph(len(Masses[mask]), array('d',Masses[mask]), array('d',Gammas[mask]))
                    graphs[particle+model+gv].SetLineWidth(2)
                    tdrDraw(graphs[particle+model+gv], "L", mcolor=colors[gv], lcolor=colors[gv], lstyle=rt.kSolid if isA else (rt.kDashed if isB else rt.kDotted))
                    # if isA: legA.AddEntry(graphs[model+gv], "g_{V} = "+str(gv), "l")
                    # elif isB: legB.AddEntry(graphs[model+gv], "g_{V} = "+str(gv), "l")
                    # else: legC.AddEntry(graphs[model+gv], "g_{V} = "+str(gv), "l")
                    if isC:
                        graphs[particle+model+gv+"leg"] = rt.TGraph(len(Masses[mask]), array('d',Masses[mask]), array('d',Gammas[mask]))
                        graphs[particle+model+gv+"leg"].SetLineStyle(rt.kSolid)
                        graphs[particle+model+gv+"leg"].SetLineColor(colors[gv])
                        leg_gv.AddEntry(graphs[particle+model+gv+"leg"], "g_{V} = "+str(gv), "l")
                    if gv=="3":
                        graphs[particle+model+gv+"leg2"] = rt.TGraph(len(Masses[mask]), array('d',Masses[mask]), array('d',Gammas[mask]))
                        graphs[particle+model+gv+"leg2"].SetLineStyle(rt.kSolid if isA else (rt.kDashed if isB else rt.kDotted))
                        graphs[particle+model+gv+"leg2"].SetLineColor(rt.kBlack)
                        leg_model.AddEntry(graphs[particle+model+gv+"leg2"], "Model "+model, "l")


            # legA.Draw("same")
            # legB.Draw("same")
            if isV: canv.SaveAs("HVT_Gamma.pdf")
            else: canv.SaveAs("HVT_Gamma"+particle+".pdf")

        # for model in ["A","B", "C", "A1", "B1"]:
        for model in ["A","B", "C"]:
            isA = "A" in model
            gv = "1" if isA else "3"
            canv = tdrCanvas("BR"+model, 800, 8300, 0.005 if isA else 0.001 , 0.26 if isA else 1, "M_{V'} [GeV]", "BR(V'#rightarrow XY) [GeV]", False, iPeriod=-1)
            canv.SetLogy(not isA)
            leg_gv = tdrLeg(0.18, 0.726, 0.30, 0.80, 0.045, 42, rt.kBlack)
            legZ   = tdrLeg(0.30, 0.60, 0.60, 0.80, 0.045, 42, rt.kBlack)
            legW   = tdrLeg(0.60, 0.60, 0.95, 0.80, 0.045, 42, rt.kBlack)
            tdrHeader(leg_gv, "Model "+model)
            tdrHeader(leg_gv, "g_{V} = "+gv, isToRemove = False)
            tdrHeader(legZ, "Z'")
            tdrHeader(legW, "W'")
            for particle in gammas:
                isZ = particle=="Zprime"
                Masses = np.array([float(x) for x in gammas[particle][model][gv]["mass"]])
                Gammas = np.array([float(x) for x in gammas[particle][model][gv]["tot"]])
                mask = Masses>1000
                if isZ:
                    decays = ["ZH","WW","ee","bb"]
                else: decays = ["HW","WZ","eve","tb"]
                for decay in decays:
                    BR = np.array([float(x) for x in gammas[particle][model][gv][decay]])/Gammas
                    graphs[model+decay] = rt.TGraph(len(Masses[mask]), array('d',Masses[mask]), array('d',BR[mask]))
                    graphs[model+decay].SetLineWidth(2 if isZ else 4)
                    tdrDraw(graphs[model+decay], "L", mcolor=colors[decay], lcolor=colors[decay], lstyle=rt.kSolid if isZ else rt.kDashed)
                    if isZ: legZ.AddEntry(graphs[model+decay], decay.replace("ee","ll").replace("bb","qq"), "l")
                    else: legW.AddEntry(graphs[model+decay], decay.replace("eve","l#nu").replace("tb","qq'").replace("HW","WH"), "l")

            legZ.Draw("same")
            legW.Draw("same")
            canv.SaveAs("HVT_BR_Model_"+model+".pdf")


if __name__ == '__main__':
    plotHVT()
