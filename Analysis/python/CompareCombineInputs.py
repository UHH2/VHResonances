import os, ROOT, glob, subprocess, math
from tdrstyle_all import *

ROOT.gInterpreter.ProcessLine('#include "'+os.environ["CMSSW_BASE"]+'/src/UHH2/VHResonances/include/constants.hpp"')
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.gStyle.SetOptStat(0)

colors = {  "h_bkg_inp"  : (ROOT.kFullSquare,   ROOT.kBlue+1),
            "h_sig_inp"  : (ROOT.kFullSquare,   ROOT.kGreen+3),
            "h_bkg_pre"  : (ROOT.kFullDotLarge, ROOT.kRed+1),
            "h_sig_pre"  : (ROOT.kFullDotLarge, ROOT.kGreen+1),
            "h_bkg_post" : (ROOT.kPlus,         ROOT.kViolet),
            "h_sig_post" : (ROOT.kPlus,         ROOT.kCyan+1),
            }

def CompareCombineInputs(year,studies,histFolder,channel,collection="Puppi", mass="M2000", mode="Exp_2"):
    Path_STORAGE = "/nfs/dust/cms/user/"+os.environ["USER"]+"/"+"WorkingArea/File/Analysis/"
    Path_ANALYSIS = os.environ["CMSSW_BASE"]+"/src/UHH2/VHResonances/Analysis/"
    AnalysisDir = Path_ANALYSIS+"Limits/"+studies+"/"
    workingDir = AnalysisDir+"/"+year+"/"+collection+"/"+channel+"/"+histFolder+"/datacards/"
    DataCard = workingDir+"DataCard"+year+"_"+mass+"_"+mode+".txt"
    # print DataCard
    # print "combine -M FitDiagnostics -d "+DataCard+ " -t -1 -n "+mass+" --saveWorkspace --saveShapes"
    # print "combine -d "+DataCard+ " -t -1 -n "+mass
    process = subprocess.Popen("combine -M FitDiagnostics -d "+DataCard+ " -t -1 -n "+mass+" --saveWorkspace --saveShapes", shell=True)
    process.wait()
    cwd = os.getcwd()
    fileAnalisis_bkg = ROOT.TFile(Path_STORAGE+year+"/SignalRegion/"+collection+"/"+channel+"/nominal/uhh2.AnalysisModuleRunner.MC.MC_DY_"+year+"_noTree.root")
    fileAnalisis_sig = ROOT.TFile(Path_STORAGE+year+"/SignalRegion/"+collection+"/"+channel+"/nominal/uhh2.AnalysisModuleRunner.MC.MC_ZprimeToZH_"+mass+"_"+year+"_noTree.root")
    h_bkg_inp = fileAnalisis_bkg.Get("ZprimeCandidate_"+histFolder+"_CR/Zprime_mass_rebin2")
    h_sig_inp = fileAnalisis_sig.Get("ZprimeCandidate_"+histFolder+"_SR/Zprime_mass_rebin2")
    h_bkg_inp.Rebin(30)
    h_sig_inp.Rebin(30)
    h_sig_inp.Scale(0.01)
    outHistFile = ROOT.TFile.Open(cwd+"/Theta_"+year+".root" ,"RECREATE")
    h_bkg_inp.Write("H__DY")
    h_sig_inp.Write("H__SIG"+mass)
    outHistFile.Close()

    fileCombine = ROOT.TFile(cwd+"/fitDiagnostics"+mass+".root")
    h_bkg_pre = fileCombine.Get("shapes_prefit/"+channel+"_"+year+"/total_background")
    h_sig_pre = fileCombine.Get("shapes_prefit/"+channel+"_"+year+"/total_signal")
    h_bkg_post = fileCombine.Get("shapes_fit_s/"+channel+"_"+year+"/total_background")
    h_sig_post = fileCombine.Get("shapes_fit_s/"+channel+"_"+year+"/total_signal")
    h_bkg_pre.Scale(h_bkg_pre.GetBinWidth(1))
    # h_bkg_pre.Scale(6156.6)
    h_sig_pre.Scale(h_sig_pre.GetBinWidth(1))
    h_bkg_post.Scale(h_bkg_post.GetBinWidth(1))
    h_sig_post.Scale(h_sig_post.GetBinWidth(1))
    # h_bkg_pre.Scale(1./0.439519)
    # print h_bkg_pre.Integral(), h_bkg_post.Integral()
    # h_bkg_post.Scale(1./0.439519)

    canv = tdrCanvas(workingDir, 300, 8200, 1e-10, 1e03, "M(Z')", "Events");
    canv.SetLogy(1)
    leg = tdrLeg(0.40,0.70,0.89,0.89, 0.030, 42, ROOT.kBlack);
    leg.SetNColumns(2);
    tdrDraw(h_bkg_pre,  "hist", colors["h_bkg_pre"][0], colors["h_bkg_pre"][1], 1, colors["h_bkg_pre"][1], 0, colors["h_bkg_pre"][1])
    leg.AddEntry(h_bkg_pre, "h_bkg_pre", "l")
    tdrDraw(h_sig_pre,  "hist", colors["h_sig_pre"][0], colors["h_sig_pre"][1], 1, colors["h_sig_pre"][1], 0, colors["h_sig_pre"][1])
    leg.AddEntry(h_sig_pre, "h_sig_pre", "l")
    tdrDraw(h_bkg_post, "hist", colors["h_bkg_post"][0], colors["h_bkg_post"][1], 1, colors["h_bkg_post"][1], 0, colors["h_bkg_post"][1])
    leg.AddEntry(h_bkg_post, "h_bkg_post", "l")
    tdrDraw(h_sig_post, "hist", colors["h_sig_post"][0], colors["h_sig_post"][1], 1, colors["h_sig_post"][1], 0, colors["h_sig_post"][1])
    leg.AddEntry(h_sig_post, "h_sig_post", "l")
    tdrDraw(h_bkg_inp, "hist", colors["h_bkg_inp"][0], colors["h_bkg_inp"][1], 1, colors["h_bkg_inp"][1], 0, colors["h_bkg_inp"][1])
    leg.AddEntry(h_bkg_inp, "h_bkg_inp", "l")
    tdrDraw(h_sig_inp, "hist", colors["h_sig_inp"][0], colors["h_sig_inp"][1], 1, colors["h_sig_inp"][1], 0, colors["h_sig_inp"][1])
    leg.AddEntry(h_sig_inp, "h_sig_inp", "l")
    leg.Draw("same")
    canv.SaveAs(cwd+"/OtherPlots/CompareCombineInputs_"+studies+"_"+year+"_"+channel+"_"+histFolder+"_"+mass+"_"+mode+".pdf")




def main():
    # histFolders = ["btag_DeepBoosted_H4qvsQCD", "btag_DeepBoosted_H4qvsQCDp2", "btag_DeepBoosted_H4qvsQCDp02", "btag_DeepBoosted_H4qvsQCDpt1000", "btag_DeepBoosted_H4qvsQCDpt1000p2", "btag_DeepBoosted_H4qvsQCDpt1000p02"]
    # channels = ["muonchannel", "electronchannel"]
    # years = ["2016", "2017", "2018", "RunII"]
    histFolders = ["btag_DeepBoosted_H4qvsQCD"]
    channels = ["muonchannel"]
    # years = ["2016","2017","2018"]
    years = ["2016"]
    studies="nominal"
    Path_STORAGE = "/nfs/dust/cms/user/"+os.environ["USER"]+"/"+"WorkingArea/File/Analysis/"

    canv = {}
    leg = {}
    h_bkg_inp = {}
    h_bkg_ratio = {}
    fileAnalisis_ratio = {}
    fileAnalisis_ratio["DY"] = ROOT.TFile(Path_STORAGE+"2016/SignalRegion/Puppi/muonchannel/nominal/uhh2.AnalysisModuleRunner.MC.MC_DY_2016_noTree.root")
    fileAnalisis_ratio["data"] = ROOT.TFile(Path_STORAGE+"2016/SignalRegion/Puppi/muonchannel/nominal/uhh2.AnalysisModuleRunner.DATA.DATA_SingleMuon_2016_noTree.root")
    fileAnalisis_ratio["DY_Presel"] = ROOT.TFile(Path_STORAGE+"2016/Selection/Puppi/muonchannel/nominal/uhh2.AnalysisModuleRunner.MC.MC_DY_2016_noTree.root")
    fileAnalisis_ratio["data_Presel"] = ROOT.TFile(Path_STORAGE+"2016/Selection/Puppi/muonchannel/nominal/uhh2.AnalysisModuleRunner.DATA.DATA_SingleMuon_2016_noTree.root")

    for mode in ["DY", "data", "DY_Presel","data_Presel"]:
        for reg in ["SR", "CR", "Sel", "Presel_Preselection", "Presel_TopJetcleaner", "Presel_NBoostedJet", "Presel_Trigger", "Presel_ZprimeReco", "Presel_ZprimeSelection", "Presel_PTMassCut"]:
            if "Presel" in mode and not "Presel_" in reg : continue
            if not "Presel" in mode and "Presel_" in reg : continue
            name = mode+"_"+reg
            if any(x in reg for x in ["Presel_Preselection", "Presel_TopJetcleaner", "Presel_NBoostedJet", "Presel_Trigger"]):
                canv[name] = tdrCanvas(name, 0, 200, 1e-01, 1e05, "m(jet)", "Events");
            else:
                canv[name] = tdrCanvas(name, 300, 4500, 1e-01, 1e05, "M(Z')", "Events");
            canv[name].SetLogy(1)
            leg[name] = tdrLeg(0.40,0.70,0.89,0.89, 0.030, 42, ROOT.kBlack);
            leg[name].SetNColumns(3);

    for year in years:
        col = ROOT.kRed+1 if year=="2016" else (ROOT.kGreen+1 if year=="2017" else (ROOT.kOrange+1 if year=="2018" else (ROOT.kBlue+1)))
        lumi = 36 if year=="2016" else (41 if year=="2017" else (60 if year=="2018" else (139)))
        line = ROOT.kSolid
        fileAnalisis = {}
        fileAnalisis["DY"] = ROOT.TFile(Path_STORAGE+year+"/SignalRegion/Puppi/muonchannel/nominal/uhh2.AnalysisModuleRunner.MC.MC_DY_"+year+"_noTree.root")
        fileAnalisis["data"] = ROOT.TFile(Path_STORAGE+year+"/SignalRegion/Puppi/muonchannel/nominal/uhh2.AnalysisModuleRunner.DATA.DATA_SingleMuon_"+year+"_noTree.root")
        fileAnalisis["DY_Presel"] = ROOT.TFile(Path_STORAGE+year+"/Selection/Puppi/muonchannel/nominal/uhh2.AnalysisModuleRunner.MC.MC_DY_"+year+"_noTree.root")
        fileAnalisis["data_Presel"] = ROOT.TFile(Path_STORAGE+year+"/Selection/Puppi/muonchannel/nominal/uhh2.AnalysisModuleRunner.DATA.DATA_SingleMuon_"+year+"_noTree.root")
        for mode in ["DY", "data", "DY_Presel","data_Presel"]:
            for reg in ["SR", "CR", "Sel", "Presel_Preselection", "Presel_TopJetcleaner", "Presel_NBoostedJet", "Presel_Trigger", "Presel_ZprimeReco", "Presel_ZprimeSelection", "Presel_PTMassCut"]:
                if ("Presel" in mode and not "Presel_" in reg) or (not "Presel" in mode and "Presel_" in reg) : continue
                name = mode+"_"+reg
                fullname = year+"_"+mode+"_"+reg
                doRescale = any(x in reg for x in ["Presel_Preselection", "Presel_TopJetcleaner", "Presel_NBoostedJet", "Presel_Trigger"])
                print fullname, name, mode, reg, doRescale
                if reg == "Sel":
                    h_bkg_inp[fullname] = fileAnalisis[mode].Get("ZprimeCandidate_Selection/Zprime_mass_rebin2")
                elif ("Presel" in mode or "Presel_" in reg):
                    if any(x in reg for x in ["Presel_ZprimeReco", "Presel_ZprimeSelection", "Presel_PTMassCut"]):
                        h_bkg_inp[fullname] = fileAnalisis[mode].Get("ZprimeCandidate_"+reg.replace("Presel_","")+"/Zprime_mass_rebin2")
                    else:
                        print "FIND", "nTopJet_"+reg.replace("Presel_","")+"/SDmass_jet"
                        h_bkg_inp[fullname] = fileAnalisis[mode].Get("nTopJet_"+reg.replace("Presel_","")+"/SDmass_jet")
                else:
                    h_bkg_inp[fullname] = fileAnalisis[mode].Get("ZprimeCandidate_btag_DeepBoosted_H4qvsQCD_"+reg+"/Zprime_mass_rebin2")
                h_bkg_inp[fullname].SetDirectory(0)
                if not doRescale:
                    h_bkg_inp[fullname].Rebin(30)
                canv[name].cd()
                tdrDraw(h_bkg_inp[fullname],  "hist", ROOT.kFullDotLarge,col, line, col, 0, col)
                leg[name].AddEntry(h_bkg_inp[fullname], year, "l")

                if reg == "Sel":
                    h_bkg_ratio[fullname] = fileAnalisis_ratio[mode].Get("ZprimeCandidate_Selection/Zprime_mass_rebin2")
                elif ("Presel" in mode or "Presel_" in reg):
                    if any(x in reg for x in ["Presel_ZprimeReco", "Presel_ZprimeSelection", "Presel_PTMassCut"]):
                        h_bkg_ratio[fullname] = fileAnalisis_ratio[mode].Get("ZprimeCandidate_"+reg.replace("Presel_","")+"/Zprime_mass_rebin2")
                    else:
                        print "FIND", "nTopJet_"+reg.replace("Presel_","")+"/SDmass_jet"
                        h_bkg_ratio[fullname] = fileAnalisis_ratio[mode].Get("nTopJet_"+reg.replace("Presel_","")+"/SDmass_jet")
                else:
                    h_bkg_ratio[fullname] = fileAnalisis_ratio[mode].Get("ZprimeCandidate_btag_DeepBoosted_H4qvsQCD_"+reg+"/Zprime_mass_rebin2")
                h_bkg_ratio[fullname].SetDirectory(0)
                if not doRescale:
                    h_bkg_ratio[fullname].Rebin(30)
                h_bkg_ratio[fullname].Divide(h_bkg_inp[fullname])
                h_bkg_ratio[fullname].Scale(math.sqrt(lumi*1./36))
                tdrDraw(h_bkg_ratio[fullname],  "hist", ROOT.kFullDotLarge,col, line, col, 0, col)

        for channel in channels:
            for histFolder in histFolders:
                    CompareCombineInputs(year,studies,histFolder,channel)
    cwd = os.getcwd()
    for mode in ["DY", "data", "DY_Presel","data_Presel"]:
        for reg in ["SR", "CR", "Sel", "Presel_Preselection", "Presel_TopJetcleaner", "Presel_NBoostedJet", "Presel_Trigger", "Presel_ZprimeReco", "Presel_ZprimeSelection", "Presel_PTMassCut"]:
            if "Presel" in mode and not "Presel_" in reg : continue
            if not "Presel" in mode and "Presel_" in reg : continue
            name = mode+"_"+reg
            canv[name].cd()
            leg[name].Draw("same")
            canv[name].SaveAs(cwd+"/OtherPlots/CompareCombineInputs_"+name+"_years.pdf")

if __name__ == '__main__':
    main()
