from Utils import *
from array import array

import tdrstyle_all as TDR
TDR.writeExtraText = True
TDR.extraText  = "Simulation"
TDR.extraText2 = "Work in progress"

ForThesis(TDR)

class HEMIssueImpact(VariablesBase):
    def __init__(self):
        VariablesBase.__init__(self)
        TDR.lumi_13TeV  = str("2018 RunCD")
        self.outdir = self.Path_ANALYSIS+"Analysis/OtherPlots/HEMIssueImpact/"
        os.system("mkdir -p "+self.outdir)
        self.histfolder = "nTopJet_DeepAk8_ZHccvsQCD_MD_SR"
        self.hname = "etaphi_1"
        self.defaultYear = "2018"


    def CalculateImpact(self):
        for ch in self.Channels:
            for mass in self.MassPointsReduced:
                fname = self.Path_STORAGE+self.defaultYear+"/SignalRegion/Puppi/"+ch+"channel/nominal/"+self.PrefixrootFile+"MC."+self.Signal+("_inv" if "inv" in ch else "")+"_M"+str(mass)+"_"+self.defaultYear+"_"+"noTree.root"
                if not os.path.isfile(fname):
                    print "NOT FOUND,", fname
                    continue
                f_ = rt.TFile(fname)

                h2D = f_.Get(self.histfolder+"/"+self.hname)
                print ch, mass, round(100.*h2D.Integral(h2D.GetXaxis().FindBin(-3.2),h2D.GetXaxis().FindBin(-1.3),h2D.GetYaxis().FindBin(-1.57),h2D.GetYaxis().FindBin(-0.87))/h2D.Integral(),2), "%"

    def Plot2D(self):
        for ch in self.Channels:
            sampleList = list(filter(lambda x: "DATA" in x, self.SubSamples_Year_Dict[self.defaultYear]))
            sampleList = list(filter(lambda x: "RunC" in x or "RunD" in x, sampleList))
            if "muo" in ch: sampleList = list(filter(lambda x: "Muon" in x, sampleList))
            if "inv" in ch: sampleList = list(filter(lambda x: "MET" in x, sampleList))
            if "ele" in ch: sampleList = list(filter(lambda x: "Ele" in x or "Pho" in x, sampleList))
            h2D_sum = 0
            for sample in sampleList:
                fname = self.Path_STORAGE+self.defaultYear+"/Preselection/Puppi/"+ch+"channel/nominal/"+self.PrefixrootFile+"DATA."+sample+"_noTree.root"
                f_ = rt.TFile(fname)
                h2D = f_.Get("nTopJet_JetDiLeptonPhiAngular/"+self.hname)
                if h2D_sum==0:
                    h2D_sum = h2D.Clone(ch)
                    h2D_sum.SetDirectory(0)
                else:
                    h2D_sum.Add(h2D)
                f_.Close()
            canv = tdrCanvas(ch, -2.6, 2.6, -3.3, 3.3, "#eta^{jet}","#phi^{jet}", square=kSquare, is2D=True)
            SetAlternative2DColor(h2D_sum)
            tdrDraw(h2D_sum, "colz")
            canv.Update()

            canv.SaveAs(self.outdir+ch+"_RunCD_"+self.defaultYear+".pdf")


def main():
    HEM = HEMIssueImpact()
    HEM.CalculateImpact()
    HEM.Plot2D()

if __name__ == '__main__':
    main()
