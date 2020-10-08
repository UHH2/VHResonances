from Utils import *
from array import array


# h1*Lumi1+h2*Lumi2/(Lumi1+Lumi2)

def main():
    dir = os.environ["CMSSW_BASE"]+"/src/UHH2/VHResonances/Analysis/ScaleFactors/Muons/"
    hname = "NUM_HighPtID_DEN_genTracks_eta_pair_newTuneP_probe_pt"
    file1 = ROOT.TFile(dir+"Muon_ID_SF_2016_RunBCDEF.root")
    file2 = ROOT.TFile(dir+"Muon_ID_SF_2016_RunGH.root")
    h1 = file1.Get(hname)
    h2 = file2.Get(hname)
    Lumi1 = 19.65606276
    Lumi2 = 16.226452636
    h1.Scale(Lumi1)
    h2.Scale(Lumi2)
    h1.Add(h2)
    h1.Scale(1./(Lumi1+Lumi2))
    file_ = ROOT.TFile(dir+"Muon_ID_SF_2016_RunBCDEFGH.root", "RECREATE")
    h1.Write(hname)
    file_.Close()
    file1.Close()
    file2.Close()

    hname = "NUM_TightRelIso_DEN_TightIDandIPCut_eta_pt"
    file1 = ROOT.TFile(dir+"Muon_Isolation_SF_2016_RunBCDEF.root")
    file2 = ROOT.TFile(dir+"Muon_Isolation_SF_2016_RunGH.root")
    h1 = file1.Get(hname)
    h2 = file2.Get(hname)
    Lumi1 = 19.65606276
    Lumi2 = 16.226452636
    h1.Scale(Lumi1)
    h2.Scale(Lumi2)
    h1.Add(h2)
    h1.Scale(1./(Lumi1+Lumi2))
    file_ = ROOT.TFile(dir+"Muon_Isolation_SF_2016_RunBCDEFGH.root", "RECREATE")
    h1.Write(hname)
    file_.Close()
    file1.Close()
    file2.Close()

    hdir  = "Mu50_OR_TkMu50_PtEtaBins/"
    hname = "abseta_pt_ratio"
    file1 = ROOT.TFile(dir+"Muon_Trigger_SF_2016_RunBCDEF.root")
    file2 = ROOT.TFile(dir+"Muon_Trigger_SF_2016_RunGH.root")
    h1 = file1.Get(hdir+hname)
    h2 = file2.Get(hdir+hname)
    h1.Scale(Lumi1)
    h2.Scale(Lumi2)
    h1.Add(h2)
    h1.Scale(1./(Lumi1+Lumi2))
    file_ = ROOT.TFile(dir+"Muon_Trigger_SF_2016_RunBCDEFGH.root", "RECREATE")
    ROOT.gDirectory.mkdir(hdir)
    file_.cd(hdir)
    h1.Write(hname)
    file_.Close()
    file1.Close()
    file2.Close()

    hdir  = "Mu50_OR_OldMu100_OR_TkMu100_PtEtaBins/"
    hname = "abseta_pt_ratio"
    file1 = ROOT.TFile(dir+"Muon_Trigger_SF_2018_BeforeMuonHLTUpdate.root")
    file2 = ROOT.TFile(dir+"Muon_Trigger_SF_2018_AfterMuonHLTUpdate.root")
    h1 = file1.Get(hdir+hname)
    h2 = file2.Get(hdir+hname)
    Lumi1 = 8.950818835
    Lumi2 = 59.832475339 - Lumi1
    h1.Scale(Lumi1)
    h2.Scale(Lumi2)
    h1.Add(h2)
    h1.Scale(1./(Lumi1+Lumi2))
    file_ = ROOT.TFile(dir+"Muon_Trigger_SF_2018.root", "RECREATE")
    ROOT.gDirectory.mkdir(hdir)
    file_.cd(hdir)
    h1.Write(hname)
    file_.Close()
    file1.Close()
    file2.Close()

    dir = os.environ["CMSSW_BASE"]+"/src/UHH2/VHResonances/Analysis/ScaleFactors/Electrons/"
    for year in ["2016", "2017", "2018"]:
        fname = "Electron_Trigger_SF_"+year+"_original.root"
        hname  = "SF_TH2F"
        file1 = ROOT.TFile(dir+fname)
        h1 = file1.Get(hname+"_Barrel")
        h2 = file1.Get(hname+"_EndCap")
        eta_bins = [-2.5,-1.44,0,1.44,2.5]
        pt_bins = [0]+[x for x in range(50,70,5)]+[x for x in range(70,250,10)]+[250,300,500,2000]
        hist = ROOT.TH2D(hname,hname,len(eta_bins)-1,array('d',eta_bins), len(pt_bins)-1, array('d',pt_bins))
        for eta in range(hist.GetNbinsX()+1):
            eta_bin = hist.GetXaxis().GetBinCenter(eta)
            for pt in range(hist.GetNbinsY()+1):
                pt_bin = hist.GetYaxis().GetBinCenter(pt)
                h1_cont = h1.GetBinContent(h1.GetXaxis().FindBin(eta_bin),h1.GetYaxis().FindBin(pt_bin))
                h2_cont = h2.GetBinContent(h2.GetXaxis().FindBin(eta_bin),h2.GetYaxis().FindBin(pt_bin))
                if (h1_cont>0 and h2_cont>0):
                    raise RuntimeError("Both histos share same bin. Not expected.")
                elif (h1_cont>0):
                    hist.SetBinContent(eta,pt,h1_cont)
                    hist.SetBinError(eta,pt,h1.GetBinError(h1.GetXaxis().FindBin(eta_bin),h1.GetYaxis().FindBin(pt_bin)))
                elif (h2_cont>0):
                    hist.SetBinContent(eta,pt,h2_cont)
                    hist.SetBinError(eta,pt,h2.GetBinError(h1.GetXaxis().FindBin(eta_bin),h1.GetYaxis().FindBin(pt_bin)))
                else:
                    hist.SetBinContent(eta,pt,1)
                    hist.SetBinError(eta,pt,0)
        file_ = ROOT.TFile(dir+fname.replace("_original",""), "RECREATE")
        hist.Write(hname)
        file_.Close()
        file1.Close()



if __name__ == '__main__':
    main()
