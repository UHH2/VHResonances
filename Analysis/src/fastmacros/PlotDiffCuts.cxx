#include "/nfs/dust/cms/user/amalara/WorkingArea/UHH2_102X_v2/CMSSW_10_2_16/src/UHH2/PersonalCode/tdrstyle_all.C"
#include "TAttLine.h"

void PlotDiffCuts() {

  writeExtraText = true;       // if extra text
  extraText  = "Work in progress" ;//"Preliminary";
  lumi_13TeV  = "";

  std::map<std::string, std::pair<int, int> > colors;
  colors["btag_DeepBoosted_H4qvsQCD"]           = std::pair<int, int>(kFullDotLarge,    kRed+1);
  colors["btag_DeepBoosted_H4qvsQCDp2"]         = std::pair<int, int>(kFullDotLarge,    kGreen+3);
  colors["btag_DeepBoosted_H4qvsQCDp02"]        = std::pair<int, int>(kFullDotLarge,    kBlue+1);
  colors["btag_DeepBoosted_H4qvsQCDpt1000"]     = std::pair<int, int>(kFullSquare,      kViolet);
  colors["btag_DeepBoosted_H4qvsQCDpt1000p2"]   = std::pair<int, int>(kFullSquare,      kGreen+1);
  colors["btag_DeepBoosted_H4qvsQCDpt1000p02"]  = std::pair<int, int>(kFullSquare,      kCyan+1);
  colors["btag_DeepBoosted_H4qvsQCDpt1500"]     = std::pair<int, int>(kFullTriangleUp,  kMagenta);
  colors["btag_DeepBoosted_H4qvsQCDpt1500p2"]   = std::pair<int, int>(kFullTriangleUp,  kTeal-7);
  colors["btag_DeepBoosted_H4qvsQCDpt1500p02"]  = std::pair<int, int>(kFullTriangleUp,  kAzure+7);

  TFile *f_DY = new TFile("/nfs/dust/cms/user/amalara/WorkingArea/File/Analysis/RunII/SignalRegion/Puppi/muonchannel/nominal/uhh2.AnalysisModuleRunner.MC.MC_DY_RunII_noTree.root");
  for (std::string hist : {"Zprime_mass_rebin_full", "H_SDmass", "H_pt"}) {
    for (std::string mode : {"CR", "SR"}) {
      TCanvas* canv = tdrCanvas((mode+hist+"DY").c_str(), 0, (hist=="H_SDmass")? 300:4000, 5*1e-01,(mode=="SR")? 1e05:1e06, (hist=="H_SDmass")?"M(H)":((hist=="H_pt")?"p_T(H)":"M(Z')"), "Events");
      canv->SetLogy();
      TLegend *leg = tdrLeg(0.50,0.60,0.78,0.9, 0.03);
      for (auto [name,color]: colors) {
        TH1F* h = (TH1F*)f_DY->Get(("ZprimeCandidate_"+name+"_"+mode+"/"+hist).c_str());
        tdrDraw(h, "", color.first, color.second, kSolid, color.second, 1, color.second);
        leg->AddEntry(h, name.c_str() ,"lp");
      }
      leg->Draw("same");
      canv->SaveAs(("PlotDiffCuts_DY_"+mode+"_"+hist+".pdf").c_str());
    }
  }

  TFile *f_data = new TFile("/nfs/dust/cms/user/amalara/WorkingArea/File/Analysis/RunII/SignalRegion/Puppi/muonchannel/nominal/uhh2.AnalysisModuleRunner.DATA.DATA_SingleMuon_RunII_noTree.root");
  for (std::string hist : {"Zprime_mass_rebin_full", "H_SDmass", "H_pt"}) {
    for (std::string mode : {"CR", "SR"}) {
      TCanvas* canv = tdrCanvas((mode+hist+"data").c_str(), 0, (hist=="H_SDmass")? 300:4000, 5*1e-01, (mode=="SR")? 1e05:1e06, (hist=="H_SDmass")?"M(H)":((hist=="H_pt")?"p_T(H)":"M(Z')"), "Events");
      canv->SetLogy();
      TLegend *leg = tdrLeg(0.50,0.60,0.78,0.9, 0.03);
      for (auto [name,color]: colors) {
        TH1F* h = (TH1F*)f_data->Get(("ZprimeCandidate_"+name+"_"+mode+"/"+hist).c_str());
        tdrDraw(h, "", colors[name].first, colors[name].second, kSolid, colors[name].second, 1, colors[name].second);
        leg->AddEntry(h, name.c_str() ,"lp");
      }
      leg->Draw("same");
      canv->SaveAs(("PlotDiffCuts_Data_"+mode+"_"+hist+".pdf").c_str());
    }
  }
}
