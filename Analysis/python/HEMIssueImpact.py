from Utils import *
from array import array

import tdrstyle_all as TDR
TDR.writeExtraText = True
TDR.extraText = "Work in progress"

# ForThesis(TDR)

class HEMIssueImpact(VariablesBase):
    def __init__(self):
        VariablesBase.__init__(self)
        TDR.lumi_13TeV  = str("2018 RunCD")
        self.outdir = self.Path_ANALYSIS+"Analysis/OtherPlots/HEMIssueImpact/"
        os.system("mkdir -p "+self.outdir)
        # self.histfolders = ["nocuts", "weights", "HEM", "cleaned", "JetDiLeptonPhiAngular", "QCDRejection", "ScaleFactors", "ExtraCleaning", "DeepAk8_ZHccvsQCD_MD_SR"]
        self.histfolders = ["cleaned", "ScaleFactors", "ExtraCleaning"]
        self.hname = "etaphi_1"
        self.defaultYear = "2018"

    def Plot2D(self):
        for mode in ["HEM","noHEM"]:
            for ch in self.Channels:
                sampleList = list(filter(lambda x: "DATA" in x, self.Processes_Year_Dict[self.defaultYear]))
                if "muo" in ch: sampleList = list(filter(lambda x: "Muon" in x, sampleList))
                if "inv" in ch: sampleList = list(filter(lambda x: "MET" in x, sampleList))
                if "ele" in ch: sampleList = list(filter(lambda x: "Ele" in x, sampleList))
                if len(sampleList)>1: raise "ERROR"
                fname = self.Path_STORAGE+self.defaultYear+"/HEMIssueStudy_"+mode+"/Puppi/"+ch+"channel/nominal/"+self.PrefixrootFile+"DATA."+sampleList[0]+"_noTree_merge.root"
                f_ = rt.TFile(fname)
                for hfolder in self.histfolders:
                    h2D = f_.Get("nTopJet_"+hfolder+"/"+self.hname)
                    print mode, ch, hfolder, round(100.*h2D.Integral(h2D.GetXaxis().FindBin(-3.2),h2D.GetXaxis().FindBin(-1.3),h2D.GetYaxis().FindBin(-1.57),h2D.GetYaxis().FindBin(-0.87))/h2D.Integral(),2), "%"
                    canv = tdrCanvas(mode+ch+hfolder, -2.6, 2.6, -3.3, 3.3, "#eta^{jet}","#phi^{jet}", square=kSquare, is2D=True, iPos=0)
                    #canv.SetLogz(1)
                    SetAlternative2DColor(h2D)
                    h2D.GetZaxis().SetTitle("Number of jets")
                    h2D.GetZaxis().SetTitleOffset(1.35)
                    h2D.RebinX(10)
                    h2D.RebinY(15)
                    tdrDraw(h2D, "colz")
                    canv.Update()
                    canv.SaveAs(self.outdir+ch+"_RunCD_"+self.defaultYear+"_"+hfolder+"_"+mode+".pdf")
                f_.Close()

def main():
    HEM = HEMIssueImpact()
    HEM.Plot2D()

if __name__ == '__main__':
    main()
