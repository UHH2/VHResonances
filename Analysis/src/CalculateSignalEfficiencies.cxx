#include "CalculateSignalEfficiencies.hpp"


double EfficiencyError(double count, double N)  { double eff = count/N; return TMath::Sqrt(eff*(1-eff)/N); }

bool isCSRegion(std::string tag) {return (tag.find("SR")!=std::string::npos || tag.find("CR")!=std::string::npos); }

void CalculateSignalEfficiencies(std::string histFolder) {

  std::string user            = std::getenv("USER");
  std::string Path_NFS        = "/nfs/dust/cms/user/"+user+"/";
  std::string Path_STORAGE    = Path_NFS+"WorkingArea/File/Analysis/";

  std::string outdir = "./SignalEfficiencies/";
  std::string prefix = "uhh2.AnalysisModuleRunner.MC.";
  std::string syst = "nominal";

  std::vector<std::string> years =  {"2016", "2017", "2018", "RunII"};
  // std::vector<std::string> years =  {"2016"};
  // std::vector<std::string> collections =  {"Puppi", "CHS", "HOTVR"};
  std::vector<std::string> collections =  {"Puppi"};
  std::vector<std::string> channels =  {"muonchannel", "electronchannel", "leptonchannel" };
  //std::vector<std::string> channels =  {"muonchannel"};
  std::vector<std::string> decaymodes =  {"bb", "WW", "Inc" };

  double x_min = 300;
  double y_min = 8500;
  double x_max = 0.001;
  double y_max = 30;
  TString x_name = "M_{Z'} (GeV)";
  TString y_name = "Selection efficiency";

  writeExtraText = true;       // if extra text
  extraText  = "Simulation";
  extraText2 = "Work in progress";
  lumi_13TeV  = "";

  std::map<std::string, mypair_I > Cuts;

  // Cuts.insert(std::pair<std::string, mypair_I>("0_OppositeLeptonVeto", mypair_I("Veto",kBlue+1)));
  Cuts.insert(std::pair<std::string, mypair_I>("1_N_{l} #geq 2",                        mypair_I("NLeptonSel",                    kRed+1)));
  Cuts.insert(std::pair<std::string, mypair_I>("2_#it{p}_{T}^{#it{jet}} #geq 200 GeV",  mypair_I("NBoostedJet",                   kOrange+1)));
  Cuts.insert(std::pair<std::string, mypair_I>("3_#Delta R(ll) < 1.0",                  mypair_I("DeltaRDiLepton",                kOrange-2)));
  Cuts.insert(std::pair<std::string, mypair_I>("4_#Delta#phi(ll,#it{jet}) #geq #pi/2",  mypair_I("JetDiLeptonPhiAngular",         kGreen+1)));
  // Cuts.insert(std::pair<std::string, mypair_I>("5_preselection",                        mypair_I("Preselection",                  kBlack)));
  Cuts.insert(std::pair<std::string, mypair_I>("6_Triggers",                            mypair_I("Trigger",                       kGreen+3)));
  // Cuts.insert(std::pair<std::string, mypair_I>("7_Z' Reconstruction",                   mypair_I("ZprimeReco",                    kGreen+2)));
  Cuts.insert(std::pair<std::string, mypair_I>("7_Z' Selection",                        mypair_I("ZprimeSelection",               kAzure+1)));
  Cuts.insert(std::pair<std::string, mypair_I>("8_#it{p}_{T}^{#it{ll}}/m_{Z'} #geq 0.2",mypair_I("PTMassCut",                     kBlue+1)));
  Cuts.insert(std::pair<std::string, mypair_I>("9_DeepBoosted",                         mypair_I(histFolder+"_SR",  kViolet+1)));
  Cuts.insert(std::pair<std::string, mypair_I>("00_H4qvsQCD_SR",                        mypair_I(histFolder+"_SR",  kGreen+3)));
  Cuts.insert(std::pair<std::string, mypair_I>("01_H4qvsQCD_CR",                        mypair_I(histFolder+"_CR",  kGreen+1)));
  Cuts.insert(std::pair<std::string, mypair_I>("10_HbbvsQCD_SR",                        mypair_I("btag_DeepBoosted_HbbvsQCD_SR",  kViolet+1)));
  Cuts.insert(std::pair<std::string, mypair_I>("11_HbbvsQCD_CR",                        mypair_I("btag_DeepBoosted_HbbvsQCD_CR",  kMagenta+1)));
  Cuts.insert(std::pair<std::string, mypair_I>("20_tau42_SR",                           mypair_I("tau42_SR",                      kOrange+1)));
  Cuts.insert(std::pair<std::string, mypair_I>("21_tau42_CR",                           mypair_I("tau42_CR",                      kOrange-2)));


  std::map<std::string, TGraph* > Plot_ComparisonFinal;
  std::map<std::string, std::vector<double> > SignalEfficiencies;
  std::map<std::string, std::vector<double> > SignalEfficiencies_err;
  std::vector<std::string> order_norm;
  std::vector<double> keepDen;
  std::vector<double> keepNum;
  std::vector<double> dummy(MassPoints.size(),0.001);

  for (std::string year: years) {
    for (std::string collection: collections) {
      for (std::string channel: channels) {
        for (std::string decaymode: decaymodes) {
          std::string namePlot = decaymode+"_"+collection+"_"+channel+"_"+year+"_"+histFolder;
          TString namePlotShort = "Hto"+namePlot; namePlotShort.ReplaceAll("Puppi_","").ReplaceAll("muon","#mu-").ReplaceAll("electron","e-").ReplaceAll("lepton","l-").ReplaceAll("_"," ");

          bool isInc = decaymode=="Inc";
          bool isLeptonChannel = channel=="leptonchannel";

          TString hname = (channel=="muonchannel") ? "sum_event_weights_ZmumuHto"+decaymode : "sum_event_weights_ZeeHto"+decaymode;
          if (isLeptonChannel) hname = "sum_event_weights_Hto"+decaymode;
          if (isInc) hname = "sum_event_weights";

          //TODO do we want these plots or the "sum_event_weights" inclusive?

          std::string PresectionStorePath = Path_STORAGE+year+"/Preselection/"+collection+"/"+channel+"/"+syst+"/";
          std::string SectionStorePath    = Path_STORAGE+year+"/Selection/"+collection+"/"+channel+"/"+syst+"/";
          std::string CSRStorePath        = Path_STORAGE+year+"/SignalRegion/"+collection+"/"+channel+"/"+syst+"/";

          // Reset values
          SignalEfficiencies.clear();
          SignalEfficiencies_err.clear();
          order_norm.clear();
          keepDen.clear(); keepNum.clear();
          keepDen = std::vector<double>(MassPoints.size(), 0.);
          keepNum = std::vector<double>(MassPoints.size(), 0.);

          for (std::pair<std::string, mypair_I> element : Cuts) {
            if (!isCSRegion(element.first)) order_norm.push_back(element.first);
            SignalEfficiencies[element.first] = std::vector<double>(MassPoints.size(), 0);
            // SignalEfficiencies_err[element.first] = std::vector<double>(MassPoints.size(), 0);
          }

          int loop = isLeptonChannel?2:1;
          for (int lep = 0; lep < loop ; lep++) {

            int Mass_index = 0;
            for (int MassPoint : MassPoints) {
              std::string MassName  = std::to_string((int)MassPoint);
              TString fn_presel = PresectionStorePath+prefix+"MC_ZprimeToZH_M"+MassName+"_"+year+"_noTree.root";
              TString fn_sel    = SectionStorePath+prefix+"MC_ZprimeToZH_M"+MassName+"_"+year+"_noTree.root";
              // TString fn_csr    = CSRStorePath+prefix+"MC_ZprimeToZHTo"+decaymode+"_M"+MassName+"_noTree.root";
              TString fn_csr    = CSRStorePath+prefix+"MC_ZprimeToZH_M"+MassName+"_"+year+"_noTree.root";

              if (isLeptonChannel) {
                fn_presel.ReplaceAll("leptonchannel","electronchannel");
                fn_sel.ReplaceAll("leptonchannel","electronchannel");
                fn_csr.ReplaceAll("leptonchannel","electronchannel");
                if (lep==1) {
                  fn_presel.ReplaceAll("Zee","Zmumu");
                  fn_sel.ReplaceAll("Zee","Zmumu");
                  fn_csr.ReplaceAll("Zee","Zmumu");
                  fn_presel.ReplaceAll("electronchannel","muonchannel");
                  fn_sel.ReplaceAll("electronchannel","muonchannel");
                  fn_csr.ReplaceAll("electronchannel","muonchannel");
                }
              }

              TFile *file_presel  = new TFile(fn_presel);
              TFile *file_sel     = new TFile(fn_sel);
              TFile *file_csr     = new TFile(fn_csr);
              // Load histograms
              for (std::pair<std::string, mypair_I> element : Cuts) {
                std::string tag = element.first;
                std::string cut = element.second.first;
                bool isCSR = isCSRegion(tag);
                TH1F* h_tot = (TH1F*)file_presel->Get("ZprimeCandidate_weights/"+hname);
                if (isCSR) h_tot = (TH1F*)file_csr->Get("ZprimeCandidate_Selection/"+hname);
                double tot_event = h_tot->GetBinContent(1);

                TH1F* h_preselection = (TH1F*)file_presel->Get("ZprimeCandidate_"+cut+"/"+hname);
                TH1F* h_selection = (TH1F*)file_sel->Get("ZprimeCandidate_"+cut+"/"+hname);
                if (isCSR) h_selection = (TH1F*)file_csr->Get("ZprimeCandidate_"+cut+"/"+hname);
                TH1F* h_sel;
                if (h_preselection) h_sel = h_preselection;
                else if (h_selection) h_sel = h_selection;
                else h_sel = (TH1F*)file_csr->Get("ZprimeCandidate_"+cut+"/"+hname);

                if (tag=="NBoostedJet") h_sel = h_preselection;

                double sel_event = h_sel->GetBinContent(1);
                if (isLeptonChannel && lep==1) SignalEfficiencies[tag][Mass_index] = (sel_event+keepNum[Mass_index])/(tot_event+keepDen[Mass_index]);
                else SignalEfficiencies[tag][Mass_index] = sel_event/tot_event;
                if (isLeptonChannel && lep==0) keepDen[Mass_index] = tot_event;
                if (isLeptonChannel && lep==0) keepNum[Mass_index] = sel_event;
              }

              file_presel->Close(); file_sel->Close(); file_csr->Close();
              Mass_index++;

            }
          }

          // Plot efficiencies
          lumi_13TeV  = TString::Format("%.1f fb^{-1}", lumi_map.at(year).at("lumi_fb"));
          y_max = isInc? 30: 1.5;

          TCanvas* canv_eff = tdrCanvas(("canv_eff"+namePlot).c_str(), x_min, y_min, x_max, y_max, x_name, y_name);
          if (isInc) canv_eff->SetLogy(1);

          TCanvas* canv_eff_rel = tdrCanvas(("canv_eff_rel"+namePlot).c_str(), x_min, y_min, x_max, y_max, x_name, y_name);
          if (isInc) canv_eff_rel->SetLogy(1);

          TLegend *leg_eff = tdrLeg(0.50,0.68,0.89,0.89, 0.030, 42, kBlack);
          leg_eff->SetNColumns(2);
          TLegend *leg_eff_rel = tdrLeg(0.40,0.70,0.89,0.89, 0.030, 42, kBlack);
          leg_eff_rel->SetNColumns(2);

          for (std::pair<std::string, mypair_I> element : Cuts) {
            std::string tag = element.first;
            if (tag.find("CR")!=std::string::npos) continue; // TODO plot also CR
            int color = element.second.second;
            bool isCSR = isCSRegion(tag);
            // TGraphErrors* gr_eff = new TGraphErrors(MassPoints.size(), &(MassPoints[0]), &(SignalEfficiencies[tag][0]), &(dummy[0]), &(SignalEfficiencies_err[tag][0]));
            TGraph* gr_eff = new TGraphErrors(MassPoints.size(), &(MassPoints[0]), &(SignalEfficiencies[tag][0]));
            gr_eff->SetLineWidth(2);
            if (isCSR) canv_eff_rel->cd();
            else canv_eff->cd();
            tdrDraw(gr_eff, "lp", kFullDotLarge, color, kSolid, color, 1000, color);
            if(tag=="9_DeepBoosted") {
              Plot_ComparisonFinal[namePlot] = gr_eff;
            }
            if (isCSR) leg_eff_rel->AddEntry(gr_eff, tag.substr(tag.find("_")+1).c_str(),"lp");
            else leg_eff->AddEntry(gr_eff, tag.substr(tag.find("_")+1).c_str(),"lp");
          }

          canv_eff->cd();
          leg_eff->Draw("same");
          canv_eff->SaveAs((outdir+"Eff_HTo"+namePlot+".pdf").c_str());
          leg_eff->Delete();

          canv_eff_rel->cd();
          leg_eff_rel->Draw("same");
          canv_eff_rel->SaveAs((outdir+"Eff_rel_HTo"+namePlot+".pdf").c_str());
          leg_eff_rel->Delete();


          TCanvas* canv_eff_norm = tdrCanvas(("canv_eff_norm"+namePlot).c_str(), x_min, y_min, x_max, y_max, x_name, y_name);
          if (isInc) canv_eff_norm->SetLogy(1);

          TLegend *leg_eff_norm = tdrLeg(0.50,0.68,0.89,0.89, 0.030, 42, kBlack);
          leg_eff_norm->SetNColumns(2);

          for (unsigned int i = 0; i < order_norm.size(); i++) {
            std::string tag = order_norm[i];
            std::string tag_norm = order_norm[std::max((int)(i-1),0)];
            int color = Cuts[tag].second;
            std::vector<double> eff_norm(MassPoints.size(), 0);
            for (unsigned int m = 0; m < MassPoints.size(); m++) eff_norm[m] = SignalEfficiencies[tag][m]/SignalEfficiencies[tag_norm][m];
            TGraph* gr_eff = new TGraphErrors(MassPoints.size(), &(MassPoints[0]), &(eff_norm[0]));
            gr_eff->SetLineWidth(2);
            tdrDraw(gr_eff, "lp", kFullDotLarge, color, kSolid, color, 1000, color);
            leg_eff_norm->AddEntry(gr_eff, tag.substr(tag.find("_")+1).c_str(),"lp");
          }

          leg_eff_norm->Draw("same");
          canv_eff_norm->SaveAs((outdir+"Eff_norm_HTo"+namePlot+".pdf").c_str());
          leg_eff_norm->Delete();

        }
      }
    }
  }

  // Now all the eff are stored. One can plot extra comparisons

  TCanvas* canv_ComparisonFinal;
  TLegend* leg_ComparisonFinal;
  int color, lineStyle, markerSyle=kFullDotLarge;
  std::string namePlot, nameLeg;

  if (FindInVector(years, "RunII")>0) {
    lumi_13TeV  = TString::Format("%.1f fb^{-1}", lumi_map.at("RunII").at("lumi_fb"));
    canv_ComparisonFinal = tdrCanvas("canv_ComparisonFinal_RunII", x_min, y_min, x_max, y_max, x_name, y_name);
    canv_ComparisonFinal->SetLogy(1);
    leg_ComparisonFinal = tdrLeg(0.40,0.68,0.89,0.89, 0.030, 42, kBlack);
    leg_ComparisonFinal->SetNColumns(2);

    color = kRed+1;     lineStyle = kSolid;  namePlot = "WW_Puppi_muonchannel_RunII_"+histFolder;     nameLeg = "HToWW #mu-channel"; tdrDraw(Plot_ComparisonFinal[namePlot], "lp", markerSyle, color, lineStyle, color, 1000, color); leg_ComparisonFinal->AddEntry(Plot_ComparisonFinal[namePlot], nameLeg.c_str(),"lp");
    color = kOrange+1;  lineStyle = kDashed; namePlot = "bb_Puppi_muonchannel_RunII_"+histFolder;     nameLeg = "HTobb #mu-channel"; tdrDraw(Plot_ComparisonFinal[namePlot], "lp", markerSyle, color, lineStyle, color, 1000, color); leg_ComparisonFinal->AddEntry(Plot_ComparisonFinal[namePlot], nameLeg.c_str(),"lp");
    color = kBlue+1;    lineStyle = kSolid;  namePlot = "WW_Puppi_electronchannel_RunII_"+histFolder; nameLeg = "HToWW e-channel";   tdrDraw(Plot_ComparisonFinal[namePlot], "lp", markerSyle, color, lineStyle, color, 1000, color); leg_ComparisonFinal->AddEntry(Plot_ComparisonFinal[namePlot], nameLeg.c_str(),"lp");
    color = kViolet+1;  lineStyle = kDashed; namePlot = "bb_Puppi_electronchannel_RunII_"+histFolder; nameLeg = "HTobb e-channel";   tdrDraw(Plot_ComparisonFinal[namePlot], "lp", markerSyle, color, lineStyle, color, 1000, color); leg_ComparisonFinal->AddEntry(Plot_ComparisonFinal[namePlot], nameLeg.c_str(),"lp");
    color = kGreen+3;   lineStyle = kSolid;  namePlot = "WW_Puppi_leptonchannel_RunII_"+histFolder;   nameLeg = "HToWW lep-channel"; tdrDraw(Plot_ComparisonFinal[namePlot], "lp", markerSyle, color, lineStyle, color, 1000, color); leg_ComparisonFinal->AddEntry(Plot_ComparisonFinal[namePlot], nameLeg.c_str(),"lp");
    color = kGreen+1;   lineStyle = kDashed; namePlot = "bb_Puppi_leptonchannel_RunII_"+histFolder;   nameLeg = "HTobb lep-channel"; tdrDraw(Plot_ComparisonFinal[namePlot], "lp", markerSyle, color, lineStyle, color, 1000, color); leg_ComparisonFinal->AddEntry(Plot_ComparisonFinal[namePlot], nameLeg.c_str(),"lp");
    leg_ComparisonFinal->Draw("same");
    canv_ComparisonFinal->SaveAs((outdir+"Eff_ComparisonFinal_RunII_"+histFolder+".pdf").c_str());
  }
  // color = kOrange-2;  lineStyle = kSolid; namePlot = "Inc_Puppi_muonchannel_RunII"; nameLeg = "HToInc #mu-channel"; tdrDraw(Plot_ComparisonFinal[namePlot], "lp", markerSyle, color, lineStyle, color, 1000, color); leg_ComparisonFinal->AddEntry(Plot_ComparisonFinal[namePlot], nameLeg.c_str(),"lp");

  lumi_13TeV  = TString::Format("%.1f fb^{-1}", lumi_map.at("RunII").at("lumi_fb"));
  canv_ComparisonFinal = tdrCanvas("canv_ComparisonFinal_Inc", x_min, y_min, x_max, y_max, x_name, y_name);
  canv_ComparisonFinal->SetLogy(1);
  leg_ComparisonFinal = tdrLeg(0.40,0.68,0.89,0.89, 0.025, 42, kBlack);
  leg_ComparisonFinal->SetNColumns(3);

  for (std::string year: years) {
    if (year=="2016")  color = kOrange+1;
    if (year=="2017")  color = kGreen+1;
    if (year=="2018")  color = kAzure+1;
    if (year=="RunII") color = kRed+1;
    for (std::string channel: channels) {
      namePlot = "Inc_Puppi_"+channel+"_"+year+"_"+histFolder;
      nameLeg = year+std::string(6-year.size(), ' ' );
      if (channel=="muonchannel") {     color = kDashed; nameLeg += "#mu-channel";}
      if (channel=="electronchannel") { color = kDotted; nameLeg += "e-channel";}
      if (channel=="leptonchannel") {   color = kSolid; nameLeg += "l-channel";}
      tdrDraw(Plot_ComparisonFinal[namePlot], "lp", markerSyle, color, lineStyle, color, 1000, color); leg_ComparisonFinal->AddEntry(Plot_ComparisonFinal[namePlot], nameLeg.c_str(),"lp");
    }
  }

  leg_ComparisonFinal->Draw("same");
  canv_ComparisonFinal->SaveAs((outdir+"Eff_ComparisonFinal_Inc_"+histFolder+".pdf").c_str());
}

//TODO import this func from Utils
int FindInVector(const std::vector<std::string>& vec, const std::string& el) {
  int index = -1;
  // Find given element in vector
  auto it = std::find(vec.begin(), vec.end(), el);
  if (it != vec.end()) index = distance(vec.begin(), it);
  return index;
}


int main(int argc, char** argv){
  gSystem->Exec("mkdir -p ./SignalEfficiencies");

  std::vector<std::string> histFolders = {"btag_DeepBoosted_H4qvsQCD", "btag_DeepBoosted_H4qvsQCDptdep", "btag_DeepBoosted_H4qvsQCDp2",
  "btag_DeepBoosted_H4qvsQCDp02", "btag_DeepBoosted_H4qvsQCDpt1000", "btag_DeepBoosted_H4qvsQCDpt1000p2", "btag_DeepBoosted_H4qvsQCDpt1000p02",
  "btag_DeepBoosted_H4qvsQCDptdep_x3", "btag_DeepBoosted_H4qvsQCDptdep_x2x3", "btag_DeepBoosted_H4qvsQCDptdep_x1x3", "btag_DeepBoosted_H4qvsQCDmassdep_x3",
  "btag_DeepBoosted_H4qvsQCDmassdep2_x3", "btag_DeepBoosted_H4qvsQCDmassdep_x2x3", "btag_DeepBoosted_H4qvsQCDmassdep_x1x3", "btag_DeepBoosted_H4qvsQCDmassdep_x1x2" };

  if (argc>1) {
    std::string histFolder;
    for (int i = 1; i < argc; i++) {
      if (std::find(histFolders.begin(), histFolders.end(), argv[i]) != histFolders.end() ) histFolder = argv[i];
    }
    std::cout << histFolder << '\n';
    CalculateSignalEfficiencies(histFolder);

  } else {
    std::cout << "more " << '\n';
    for (std::string histFolder: histFolders) CalculateSignalEfficiencies(histFolder);
  }
  return 0;
}
