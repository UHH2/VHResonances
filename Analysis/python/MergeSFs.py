from Utils import *


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



if __name__ == '__main__':
    main()
