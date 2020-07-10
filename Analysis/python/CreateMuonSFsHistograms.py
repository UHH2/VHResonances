from Utils import *
import math


'''
More details here https://twiki.cern.ch/twiki/bin/viewauth/CMS/MuonReferenceSelectionAndCalibrationsRun2#Special_systematic_uncertainties

Iso SFs:
    Conservative uncertainties: statistics (0.5%) + systematic (2%) as described here (slide 11).
'''

class CreateMuonSFsHistograms():
    def __init__(self):
        self.dir = os.environ["CMSSW_BASE"]+"/src/UHH2/VHResonances/Analysis/ScaleFactors/Muons/"
        self.fname = "Muon_mode_SF_2016.root"
        self.hname = "mode_SF"

    def TrackingSF(self, mode="Tracking"):
        fname = self.fname.replace("mode",mode)
        hname = self.hname.replace("mode",mode)
        hist = ROOT.TH2D(hname,hname,1,0,2.4,2,0,600)
        for eta in range(hist.GetNbinsX()+1):
            for pt in range(hist.GetNbinsY()+1):
                hist.SetBinContent(eta,pt,1)
                hist.SetBinError(eta,pt,0.005 if (hist.GetYaxis().GetBinUpEdge(pt)<=300) else 0.01 )
                # print eta, pt, hist.GetXaxis().GetBinCenter(eta), hist.GetYaxis().GetBinCenter(pt), hist.GetXaxis().GetBinUpEdge(eta), hist.GetYaxis().GetBinUpEdge(pt)
        for year in ["2016", "2017", "2018"]:
            file_ = ROOT.TFile(self.dir+fname.replace("2016",year), "RECREATE")
            hist.Write()
            file_.Close()

    def ReconstructionSF(self, mode="Reconstruction"):
        def SF(a_Data, b_Data, a_MC, b_MC, pt):
            if a_Data==1 and b_Data==0: return (a_MC + b_MC*pt)
            else: return (a_Data + b_Data*pt)/(a_MC + b_MC*pt)

        fname = self.fname.replace("mode",mode)
        hname = self.hname.replace("mode",mode)
        hist = ROOT.TH2D(hname,hname,3,0,2.4,200,0,1000) # arbitrary choice of pt bins of 5GeV
        for eta in range(hist.GetNbinsX()+1):
            x = hist.GetXaxis().GetBinCenter(eta)
            for pt in range(hist.GetNbinsY()+1):
                y = hist.GetYaxis().GetBinCenter(pt)
                binerror = 0.
                if x<=1.6:
                    if y>=100:
                        binerror1 = abs(1-SF(1, 0, 0.9936, -3.71e-6, pt))
                        binerror2 = abs(1-SF(0.9828, -1.947e-5, 0.989, -2.399e-6, pt))
                        binerror = math.sqrt(binerror1**2+binerror2**2)
                else:
                    if y>=200:
                        binerror = abs(1-SF(0.9784, -4.73e-5, 0.9908, -1.26e-5, pt))
                        if y>=275:
                            binerror2 = abs(1-SF(0.9893, -3.666e-5, 0.9974, -1.721e-5, pt))
                            binerror = math.sqrt(binerror**2+binerror2**2)
                hist.SetBinContent(eta,pt,1)
                hist.SetBinError(eta,pt,binerror)
                # print eta, pt, x, y, hist.GetXaxis().GetBinUpEdge(eta), hist.GetYaxis().GetBinUpEdge(pt)
        for year in ["2016", "2017", "2018"]:
            file_ = ROOT.TFile(self.dir+fname.replace("2016",year), "RECREATE")
            hist.Write()
            file_.Close()


def main():
    SFCreator = CreateMuonSFsHistograms()
    SFCreator.TrackingSF()
    SFCreator.ReconstructionSF()



if __name__ == '__main__':
    main()
