#include "CalculateSignalEfficiencies.hpp"


double EfficiencyError(double count, double N)  { double eff = count/N; return TMath::Sqrt(eff*(1-eff)/N); }

bool isCSRegion(std::string tag) {return (tag.find("SR")!=std::string::npos || tag.find("CR")!=std::string::npos); }

void CalculateSignalEfficiencies() {

  std::string outdir = "./SignalEfficiencies/";
  std::string prefix = "uhh2.AnalysisModuleRunner.MC.";
  std::string syst = "nominal";

  std::vector<std::string> years =  {"2016", "2017", "2018", "RunII"};
  // std::vector<std::string> years =  {"2018"};
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
  extraText  = "Work in progress" ;//"Preliminary";
  lumi_13TeV  = "";

  std::map<std::string, TGraph* > Plot_ComparisonFinal;


  for (std::string year: years) {
    for (std::string collection: collections) {
      for (std::string channel: channels) {
        for (std::string decaymode: decaymodes) {
          std::string namePlot = decaymode+"_"+collection+"_"+channel+"_"+year;
          TString namePlotShort = "Hto"+namePlot; namePlotShort.ReplaceAll("Puppi_","").ReplaceAll("muon","#mu-").ReplaceAll("electron","e-").ReplaceAll("lepton","l-").ReplaceAll("_"," ");

          bool isInc = decaymode=="Inc";
          bool isLeptonChannel = channel=="leptonchannel";

          TString hname = (channel=="muonchannel") ? "sum_event_weights_ZmumuHto"+decaymode : "sum_event_weights_ZeeHto"+decaymode;
          if (isLeptonChannel) hname = "sum_event_weights_Hto"+decaymode;
          if (isInc) hname = "sum_event_weights";

          //TODO do we want these plots or the "sum_event_weights" inclusive?

          std::string user            = std::getenv("USER");
          std::string Path_NFS        = "/nfs/dust/cms/user/"+user+"/";
          std::string Path_STORAGE    = Path_NFS+"WorkingArea/File/Analysis/";
          std::string PresectionStorePath = Path_STORAGE+year+"/Preselection/"+collection+"/"+channel+"/"+syst+"/";
          std::string SectionStorePath    = Path_STORAGE+year+"/Selection/"+collection+"/"+channel+"/"+syst+"/";
          std::string CSRStorePath        = Path_STORAGE+year+"/SignalRegion/"+collection+"/"+channel+"/"+syst+"/";

          std::map<std::string, std::vector<double> > SignalEfficiencies;
          // std::map<std::string, std::vector<double> > SignalEfficiencies_err;
          std::map<std::string, mypair_I > Cuts;




          // kViolet+1
          // kMagenta+1

          // Cuts.insert(std::pair<std::string, mypair_I>("0_OppositeLeptonVeto", mypair_I("Veto",kBlue+1)));
          Cuts.insert(std::pair<std::string, mypair_I>("1_N_{l} #geq 2",                        mypair_I("NLeptonSel",                    kRed+1)));
          Cuts.insert(std::pair<std::string, mypair_I>("2_#it{p}_{T}^{#it{jet}} #geq 200 GeV",  mypair_I("NBoostedJet",                   kOrange+1)));
          Cuts.insert(std::pair<std::string, mypair_I>("3_#Delta R(ll) < 1.0",                  mypair_I("DeltaRDiLepton",                kOrange-2)));
          Cuts.insert(std::pair<std::string, mypair_I>("4_#Delta#phi(ll,#it{jet}) #geq #pi/2",  mypair_I("JetDiLeptonPhiAngular",         kGreen+1)));
          // Cuts.insert(std::pair<std::string, mypair_I>("5_preselection",                        mypair_I("Preselection",                  kBlack)));
          // Cuts.insert(std::pair<std::string, mypair_I>("6_NBoostedJet",                        mypair_I("NBoostedJet",                  kBlack)));
          Cuts.insert(std::pair<std::string, mypair_I>("6_Triggers",                            mypair_I("Trigger",                       kGreen+3)));
          // Cuts.insert(std::pair<std::string, mypair_I>("7_Z' Reconstruction",                   mypair_I("ZprimeReco",                    kGreen+2)));
          Cuts.insert(std::pair<std::string, mypair_I>("7_Z' Selection",                        mypair_I("ZprimeSelection",               kAzure+1)));
          Cuts.insert(std::pair<std::string, mypair_I>("8_#it{p}_{T}^{#it{ll}}/m_{Z'} #geq 0.2",mypair_I("PTMassCut",                     kBlue+1)));
          Cuts.insert(std::pair<std::string, mypair_I>("9_DeepBoosted",                         mypair_I("btag_DeepBoosted_H4qvsQCD_SR",  kViolet+1)));
          Cuts.insert(std::pair<std::string, mypair_I>("00_H4qvsQCD_SR",                        mypair_I("btag_DeepBoosted_H4qvsQCD_SR",  kGreen+3)));
          Cuts.insert(std::pair<std::string, mypair_I>("01_H4qvsQCD_CR",                        mypair_I("btag_DeepBoosted_H4qvsQCD_CR",  kGreen+1)));
          Cuts.insert(std::pair<std::string, mypair_I>("10_HbbvsQCD_SR",                        mypair_I("btag_DeepBoosted_HbbvsQCD_SR",  kViolet+1)));
          Cuts.insert(std::pair<std::string, mypair_I>("11_HbbvsQCD_CR",                        mypair_I("btag_DeepBoosted_HbbvsQCD_CR",  kMagenta+1)));
          Cuts.insert(std::pair<std::string, mypair_I>("20_tau42_SR",                           mypair_I("tau42_SR",                      kOrange+1)));
          Cuts.insert(std::pair<std::string, mypair_I>("21_tau42_CR",                           mypair_I("tau42_CR",                      kOrange-2)));
          // Cuts.insert(std::pair<std::string, mypair_I>("20_NN_SR",    mypair_I("NN_SR",                         kAzure+1)));
          // Cuts.insert(std::pair<std::string, mypair_I>("21_NN_CR",    mypair_I("NN_CR",                         kBlue+1)));
          // Cuts.insert(std::pair<std::string, mypair_I>("40_NN_SR_1",  mypair_I("NN_1_SR",                       kOrange+1)));
          // Cuts.insert(std::pair<std::string, mypair_I>("41_NN_CR_1",  mypair_I("NN_1_CR",                       kOrange)));
          // Cuts.insert(std::pair<std::string, mypair_I>("50_NN_SR_2",  mypair_I("NN_2_SR",                       kCyan+1)));
          // Cuts.insert(std::pair<std::string, mypair_I>("51_NN_CR_2",  mypair_I("NN_2_CR",                       kCyan)));
          // Cuts.insert(std::pair<std::string, mypair_I>("60_CNN_SR",   mypair_I("CNN_SR",                        kBlue+1)));
          // Cuts.insert(std::pair<std::string, mypair_I>("61_CNN_CR",   mypair_I("CNN_CR",                        kBlue)));


          std::vector<std::string> order_norm;
          for (std::pair<std::string, mypair_I> element : Cuts) {
            if (!isCSRegion(element.first)) order_norm.push_back(element.first);
            SignalEfficiencies[element.first] = std::vector<double>(MassPoints.size(), 0);
            // SignalEfficiencies_err[element.first] = std::vector<double>(MassPoints.size(), 0);
          }

          std::vector<double> keepDen(MassPoints.size(), 0);
          std::vector<double> keepNum(MassPoints.size(), 0);

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
              }

              if (isLeptonChannel && lep==1) {
                fn_presel.ReplaceAll("Zee","Zmumu");
                fn_sel.ReplaceAll("Zee","Zmumu");
                fn_csr.ReplaceAll("Zee","Zmumu");
                fn_presel.ReplaceAll("electronchannel","muonchannel");
                fn_sel.ReplaceAll("electronchannel","muonchannel");
                fn_csr.ReplaceAll("electronchannel","muonchannel");
              }

              TFile *file_presel    = new TFile(fn_presel);
              TFile *file_sel       = new TFile(fn_sel);
              TFile *file_csr       = new TFile(fn_csr);

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
                //TODO Check inclusive plots
                // std::cout << sel_event << " " << sel_event_inc << " " << tot_event << " " << inc_event << " " << sel_event/tot_event << " " << sel_event_inc/inc_event << '\n';
                // if (isCSR) {
                //   std::cout << Mass_index << " " << tag << " " << cut << " " << tot_event << " " << sel_event << " " << tot_event-sel_event << '\n';
                // }// TODO questa cosa nnon torna. Controlla CR vs SR

              }

              file_presel->Close(); file_sel->Close(); file_csr->Close();
              Mass_index++;

            }
          }


          lumi_13TeV  = TString::Format("%.1f fb^{-1}", lumi_map.at(year).at("lumi_fb"));

          y_max = isInc? 30: 1.4;

          TCanvas* canv_eff = tdrCanvas(("canv_eff"+namePlot).c_str(), x_min, y_min, x_max, y_max, x_name, y_name);
          if (isInc) canv_eff->SetLogy(1);

          TCanvas* canv_eff_rel = tdrCanvas(("canv_eff_rel"+namePlot).c_str(), x_min, y_min, x_max, y_max, x_name, y_name);
          if (isInc) canv_eff_rel->SetLogy(1);

          TLegend *leg_eff = tdrLeg(0.50,0.68,0.89,0.89, 0.030, 42, kBlack);
          leg_eff->SetNColumns(2);
          TLegend *leg_eff_rel = tdrLeg(0.40,0.70,0.89,0.89, 0.030, 42, kBlack);
          leg_eff_rel->SetNColumns(2);

          std::vector<double> dummy(MassPoints.size(),0.001);
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
            // 1_N_{l} #geq 2
            // 2_#it{p}_{T}^{#it{jet}} #geq 200 GeV
            // 3_#Delta R(ll) < 1.0
            // 4_#Delta#phi(ll,#it{jet}) #geq #pi/2
            // 6_Triggers
            // 7_Z' Selection
            // 8_#it{p}_{T}^{#it{ll}}/m_{Z'} #geq 0.2
            // 9_DeepBoosted
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

  // std::unordered_map<std::string, int> Colors;
  // Colors.insert(mypair_I("2016_muonchannel",      kRed+1));
  // Colors.insert(mypair_I("2016_electronchannel",  kOrange+1));
  // Colors.insert(mypair_I("2016_leptonchannel",    kOrange-2));
  // Colors.insert(mypair_I("2017_muonchannel",      kGreen+1));
  // Colors.insert(mypair_I("2017_electronchannel",  kGreen+3));
  // Colors.insert(mypair_I("2017_leptonchannel",    kGreen+2));
  // Colors.insert(mypair_I("2018_muonchannel",      kAzure+1));
  // Colors.insert(mypair_I("2018_electronchannel",  kBlue+1));
  // Colors.insert(mypair_I("2018_leptonchannel",    kViolet+1));
  // Colors.insert(mypair_I("RunII_muonchannel",     kRed+1));
  // Colors.insert(mypair_I("RunII_electronchannel", kOrange+1));
  // Colors.insert(mypair_I("RunII_leptonchannel",   kOrange-2));
  //
  //
  //
  // Colors.insert(mypair_I("bb_RunII_muonchannel",      kRed+1));
  // Colors.insert(mypair_I("bb_RunII_electronchannel",  kOrange+1));
  // Colors.insert(mypair_I("bb_RunII_leptonchannel",    kOrange-2));
  // Colors.insert(mypair_I("WW_RunII_muonchannel",      kAzure+1));
  // Colors.insert(mypair_I("WW_RunII_electronchannel",  kBlue+1));
  // Colors.insert(mypair_I("WW_RunII_leptonchannel",    kViolet+1));
  // Colors.insert(mypair_I("Inc_RunII_muonchannel",     kGreen+1));
  // Colors.insert(mypair_I("Inc_RunII_electronchannel", kGreen+3));
  // Colors.insert(mypair_I("Inc_RunII_leptonchannel",   kGreen+2));
  //
  // std::map<std::string, TLegend*> leg_ComparisonFinal;
  // std::map<std::string, TCanvas*> canv_ComparisonFinal;
  // canv_ComparisonFinal.insert(std::pair<std::string, TCanvas*>("RunII",     tdrCanvas("canv_ComparisonFinal_RunII",     x_min, y_min, x_max, y_max, x_name, y_name)));
  // canv_ComparisonFinal.insert(std::pair<std::string, TCanvas*>("muonchannel",     tdrCanvas("canv_ComparisonFinal_muonchannel",     x_min, y_min, x_max, y_max, x_name, y_name)));
  // canv_ComparisonFinal.insert(std::pair<std::string, TCanvas*>("electronchannel", tdrCanvas("canv_ComparisonFinal_electronchannel", x_min, y_min, x_max, y_max, x_name, y_name)));
  // canv_ComparisonFinal.insert(std::pair<std::string, TCanvas*>("leptonchannel",   tdrCanvas("canv_ComparisonFinal_leptonchannel",   x_min, y_min, x_max, y_max, x_name, y_name)));
  //
  // for (auto canv: canv_ComparisonFinal) {
  //   canv.second->SetLogy(1);
  //   leg_ComparisonFinal.insert(std::pair<std::string, TLegend*>(canv.first, tdrLeg(0.40,0.68,0.89,0.89, 0.025, 42, kBlack)));
  //   leg_ComparisonFinal[canv.first]->SetNColumns(2);
  // }

  // for (std::string name_ : {"RunII"} ) {
  //   if(std::string(namePlotShort).find(name_) == std::string::npos) continue;
  //   if(canv_ComparisonFinal.find(name_) != canv_ComparisonFinal.end()) {
  //     canv_ComparisonFinal[name_]->cd();
  //     TGraph* gr_eff_temp = gr_eff;
  //     tdrDraw(gr_eff_temp, "lp", kFullDotLarge, Colors[decaymode+"_"+year+"_"+channel], isInc? kSolid: (decaymode=="WW"? kDashed:kDotted), Colors[decaymode+"_"+year+"_"+channel], 1000, Colors[decaymode+"_"+year+"_"+channel]);
  //     leg_ComparisonFinal[name_]->AddEntry(gr_eff_temp, namePlotShort,"lp");
  //   }
  // }
  //
  // for (auto x: Plot_ComparisonFinal) {
  //   x.second->cd();
  //   leg_ComparisonFinal[canv.first]->Draw("same");
  //   x.second->SaveAs((outdir+"Eff_ComparisonFinal_"+canv.first+".pdf").c_str());
  //   leg_ComparisonFinal[canv.first]->Delete();
  // }
  TCanvas* canv_ComparisonFinal;
  TLegend* leg_ComparisonFinal;
  int color, lineStyle, markerSyle=kFullDotLarge;
  std::string namePlot, nameLeg;

  lumi_13TeV  = TString::Format("%.1f fb^{-1}", lumi_map.at("RunII").at("lumi_fb"));
  canv_ComparisonFinal = tdrCanvas("canv_ComparisonFinal_RunII", x_min, y_min, x_max, y_max, x_name, y_name);
  canv_ComparisonFinal->SetLogy(1);
  leg_ComparisonFinal = tdrLeg(0.40,0.68,0.89,0.89, 0.030, 42, kBlack);
  leg_ComparisonFinal->SetNColumns(2);
  color = kRed+1;     lineStyle = kSolid;  namePlot = "WW_Puppi_muonchannel_RunII";     nameLeg = "HToWW #mu-channel"; tdrDraw(Plot_ComparisonFinal[namePlot], "lp", markerSyle, color, lineStyle, color, 1000, color); leg_ComparisonFinal->AddEntry(Plot_ComparisonFinal[namePlot], nameLeg.c_str(),"lp");
  color = kOrange+1;  lineStyle = kDashed; namePlot = "bb_Puppi_muonchannel_RunII";     nameLeg = "HTobb #mu-channel"; tdrDraw(Plot_ComparisonFinal[namePlot], "lp", markerSyle, color, lineStyle, color, 1000, color); leg_ComparisonFinal->AddEntry(Plot_ComparisonFinal[namePlot], nameLeg.c_str(),"lp");
  color = kBlue+1;    lineStyle = kSolid;  namePlot = "WW_Puppi_electronchannel_RunII"; nameLeg = "HToWW e-channel";   tdrDraw(Plot_ComparisonFinal[namePlot], "lp", markerSyle, color, lineStyle, color, 1000, color); leg_ComparisonFinal->AddEntry(Plot_ComparisonFinal[namePlot], nameLeg.c_str(),"lp");
  color = kViolet+1;  lineStyle = kDashed; namePlot = "bb_Puppi_electronchannel_RunII"; nameLeg = "HTobb e-channel";   tdrDraw(Plot_ComparisonFinal[namePlot], "lp", markerSyle, color, lineStyle, color, 1000, color); leg_ComparisonFinal->AddEntry(Plot_ComparisonFinal[namePlot], nameLeg.c_str(),"lp");
  color = kGreen+3;   lineStyle = kSolid;  namePlot = "WW_Puppi_leptonchannel_RunII";   nameLeg = "HToWW lep-channel"; tdrDraw(Plot_ComparisonFinal[namePlot], "lp", markerSyle, color, lineStyle, color, 1000, color); leg_ComparisonFinal->AddEntry(Plot_ComparisonFinal[namePlot], nameLeg.c_str(),"lp");
  color = kGreen+1;   lineStyle = kDashed; namePlot = "bb_Puppi_leptonchannel_RunII";   nameLeg = "HTobb lep-channel"; tdrDraw(Plot_ComparisonFinal[namePlot], "lp", markerSyle, color, lineStyle, color, 1000, color); leg_ComparisonFinal->AddEntry(Plot_ComparisonFinal[namePlot], nameLeg.c_str(),"lp");
  leg_ComparisonFinal->Draw("same");
  canv_ComparisonFinal->SaveAs((outdir+"Eff_ComparisonFinal_RunII.pdf").c_str());

  // color = kOrange-2;  lineStyle = kSolid; namePlot = "Inc_Puppi_muonchannel_RunII"; nameLeg = "HToInc #mu-channel"; tdrDraw(Plot_ComparisonFinal[namePlot], "lp", markerSyle, color, lineStyle, color, 1000, color); leg_ComparisonFinal->AddEntry(Plot_ComparisonFinal[namePlot], nameLeg.c_str(),"lp");

  lumi_13TeV  = TString::Format("%.1f fb^{-1}", lumi_map.at("RunII").at("lumi_fb"));
  canv_ComparisonFinal = tdrCanvas("canv_ComparisonFinal_Inc", x_min, y_min, x_max, y_max, x_name, y_name);
  canv_ComparisonFinal->SetLogy(1);
  leg_ComparisonFinal = tdrLeg(0.40,0.68,0.89,0.89, 0.025, 42, kBlack);
  leg_ComparisonFinal->SetNColumns(3);
  color = kOrange+1;  lineStyle = kDashed;  namePlot = "Inc_Puppi_muonchannel_2016";      nameLeg = "2016  #mu-channel"; tdrDraw(Plot_ComparisonFinal[namePlot], "lp", markerSyle, color, lineStyle, color, 1000, color); leg_ComparisonFinal->AddEntry(Plot_ComparisonFinal[namePlot], nameLeg.c_str(),"lp");
  color = kOrange+1;  lineStyle = kDotted;  namePlot = "Inc_Puppi_electronchannel_2016";  nameLeg = "2016  e-channel";   tdrDraw(Plot_ComparisonFinal[namePlot], "lp", markerSyle, color, lineStyle, color, 1000, color); leg_ComparisonFinal->AddEntry(Plot_ComparisonFinal[namePlot], nameLeg.c_str(),"lp");
  color = kOrange+1;  lineStyle = kSolid;   namePlot = "Inc_Puppi_leptonchannel_2016";    nameLeg = "2016  l-channel";   tdrDraw(Plot_ComparisonFinal[namePlot], "lp", markerSyle, color, lineStyle, color, 1000, color); leg_ComparisonFinal->AddEntry(Plot_ComparisonFinal[namePlot], nameLeg.c_str(),"lp");
  color = kGreen+1;   lineStyle = kDashed;  namePlot = "Inc_Puppi_muonchannel_2017";      nameLeg = "2017  #mu-channel"; tdrDraw(Plot_ComparisonFinal[namePlot], "lp", markerSyle, color, lineStyle, color, 1000, color); leg_ComparisonFinal->AddEntry(Plot_ComparisonFinal[namePlot], nameLeg.c_str(),"lp");
  color = kGreen+1;   lineStyle = kDotted;  namePlot = "Inc_Puppi_electronchannel_2017";  nameLeg = "2017  e-channel";   tdrDraw(Plot_ComparisonFinal[namePlot], "lp", markerSyle, color, lineStyle, color, 1000, color); leg_ComparisonFinal->AddEntry(Plot_ComparisonFinal[namePlot], nameLeg.c_str(),"lp");
  color = kGreen+1;   lineStyle = kSolid;   namePlot = "Inc_Puppi_leptonchannel_2017";    nameLeg = "2017  l-channel";   tdrDraw(Plot_ComparisonFinal[namePlot], "lp", markerSyle, color, lineStyle, color, 1000, color); leg_ComparisonFinal->AddEntry(Plot_ComparisonFinal[namePlot], nameLeg.c_str(),"lp");
  color = kAzure+1;   lineStyle = kDashed;  namePlot = "Inc_Puppi_muonchannel_2018";      nameLeg = "2018  #mu-channel"; tdrDraw(Plot_ComparisonFinal[namePlot], "lp", markerSyle, color, lineStyle, color, 1000, color); leg_ComparisonFinal->AddEntry(Plot_ComparisonFinal[namePlot], nameLeg.c_str(),"lp");
  color = kAzure+1;   lineStyle = kDotted;  namePlot = "Inc_Puppi_electronchannel_2018";  nameLeg = "2018  e-channel";   tdrDraw(Plot_ComparisonFinal[namePlot], "lp", markerSyle, color, lineStyle, color, 1000, color); leg_ComparisonFinal->AddEntry(Plot_ComparisonFinal[namePlot], nameLeg.c_str(),"lp");
  color = kAzure+1;   lineStyle = kSolid;   namePlot = "Inc_Puppi_leptonchannel_2018";    nameLeg = "2018  l-channel";   tdrDraw(Plot_ComparisonFinal[namePlot], "lp", markerSyle, color, lineStyle, color, 1000, color); leg_ComparisonFinal->AddEntry(Plot_ComparisonFinal[namePlot], nameLeg.c_str(),"lp");
  color = kRed+1;     lineStyle = kDashed;  namePlot = "Inc_Puppi_muonchannel_RunII";     nameLeg = "RunII #mu-channel"; tdrDraw(Plot_ComparisonFinal[namePlot], "lp", markerSyle, color, lineStyle, color, 1000, color); leg_ComparisonFinal->AddEntry(Plot_ComparisonFinal[namePlot], nameLeg.c_str(),"lp");
  color = kRed+1;     lineStyle = kDotted;  namePlot = "Inc_Puppi_electronchannel_RunII"; nameLeg = "RunII e-channel";   tdrDraw(Plot_ComparisonFinal[namePlot], "lp", markerSyle, color, lineStyle, color, 1000, color); leg_ComparisonFinal->AddEntry(Plot_ComparisonFinal[namePlot], nameLeg.c_str(),"lp");
  color = kRed+1;     lineStyle = kSolid;   namePlot = "Inc_Puppi_leptonchannel_RunII";   nameLeg = "RunII l-channel";   tdrDraw(Plot_ComparisonFinal[namePlot], "lp", markerSyle, color, lineStyle, color, 1000, color); leg_ComparisonFinal->AddEntry(Plot_ComparisonFinal[namePlot], nameLeg.c_str(),"lp");
  leg_ComparisonFinal->Draw("same");
  canv_ComparisonFinal->SaveAs((outdir+"Eff_ComparisonFinal_Inc.pdf").c_str());
  //
  // for (int i = 0; i < 18; i++) {
  //   double x,y;
  //   Plot_ComparisonFinal["WW_Puppi_muonchannel_RunII"]->GetPoint(i,x,y); std::cout << i << " " << x << " " << y << '\n';
  //   Plot_ComparisonFinal["WW_Puppi_electronchannel_RunII"]->GetPoint(i,x,y); std::cout << i << " " << x << " " << y << '\n';
  //   Plot_ComparisonFinal["WW_Puppi_leptonchannel_RunII"]->GetPoint(i,x,y); std::cout << i << " " << x << " " << y << '\n';
  //   std::cout << '\n';
  // }

}


int main() {
  gSystem->Exec("mkdir -p ./SignalEfficiencies");
  CalculateSignalEfficiencies();
  return 0;
}
