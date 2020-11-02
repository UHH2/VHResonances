#include "PlotLimits.hpp"



void PlotLimits(bool doObs, bool isHbb) {

  // double plot_hi = 2800;
  const double plot_lo = 300;
  const double plot_hi = 8200;

  const double yaxis_lo = 1e-01;
  const double yaxis_hi = 5e03;

  std::map<std::string, int> colors = {
    {"2016",kViolet-9+1}, {"2017",kOrange+1}, {"2018",kCyan+1}, {"RunII",kGreen+2}, {"fullRunII", kGreen+1},
    {"NN",kOrange+1}, {"NN_1",kViolet-9+1}, {"NN_2",kCyan+1}, {"CNN",kBlue+1}, {"tau42",kBlack},
    {"btag_DeepBoosted_probHbb",kViolet-9+1}, {"tau32",kBlue+1}, {"tau21",kOrange+1},
    {"btag_DeepBoosted_H4qvsQCDp2", kOrange+1}, {"btag_DeepBoosted_H4qvsQCDp02", kYellow+1}, {"btag_DeepBoosted_H4qvsQCDpt1000", kCyan+1},
    {"btag_DeepBoosted_H4qvsQCDpt1000p2", kBlue+1}, {"btag_DeepBoosted_H4qvsQCDpt1000p02", kRed+1},
    {"btag_DeepBoosted_H4qvsQCDptdep", kRed+1},
    {"btag_DeepBoosted_H4qvsQCD",               kAzure+7},
    {"btag_DeepBoosted_H4qvsQCDmassdep",        kRed+1},
    {"btag_DeepBoosted_H4qvsQCDmassdep_2",      kViolet-9},
    {"btag_DeepBoosted_H4qvsQCDmassdep_3",      kOrange-1},
    {"btag_DeepBoosted_H4qvsQCDmassdep_bb",     kGreen+2},
    {"btag_DeepBoosted_H4qvsQCDmassdep_cc",     kViolet+1},
    {"btag_DeepBoosted_H4qvsQCDmassdep_gg",     kAzure+1},
    {"btag_DeepBoosted_H4qvsQCDmassdep_cc_2",   kViolet-7},
    {"btag_DeepBoosted_H4qvsQCDmassdep_gg_2",   kGreen+1},
    {"btag_DeepBoosted_H4qvsQCDmassdep_cc_3",   kCyan+1},
    {"btag_DeepBoosted_H4qvsQCDmassdep_gg_3",   kBlue+1},
    {"btag_DeepBoosted_H4qvsQCDmassdep_cc2_3",  kOrange+1},
    {"btag_DeepBoosted_H4qvsQCDmassdep_cc2",    kBlue+1},
    {"btag_DeepBoosted_H4qvsQCDmassdep_cc2_2",  kCyan+2},
    {"btag_DeepBoosted_H4qvsQCDmassdep_bb_2",   kViolet-1},
    {"btag_DeepBoosted_H4qvsQCDmassdep_cc1",    kViolet-1},
    {"btag_DeepBoosted_H4qvsQCDmassdep_ccMD",   kCyan-1},
  };

  TGraph* gr_theo    = new TGraphErrors(nPoints, &(MyMassPoints[0]), &theo_xsec[0]);
  TGraph* gr_HVT_A   = new TGraphErrors(theo_mass.size(), &(theo_mass[0]), &HVT_A_xsec[0]);
  TGraph* gr_HVT_B   = new TGraphErrors(theo_mass.size(), &(theo_mass[0]), &HVT_B_xsec[0]);
  TGraph* gr_Hbb     = new TGraphErrors(MassPoints_Hbb.size(), &(MassPoints_Hbb[0]), &expectedHbb[0]);
  TGraph* gr_Hbb0b   = new TGraphErrors(MassPoints_Hbb.size(), &(MassPoints_Hbb[0]), &expectedHbb0b[0]);
  TGraph* gr_ZeeH0b  = new TGraphErrors(ee_mass.size(), &(ee_mass[0]), &ee_xsec[0]);
  // TGraph* gr_H0lll2b = new TGraphErrors(MassPoints_Hbb.size(), &(MassPoints_Hbb[0]), &expectedH0lll2b[0]);
  // TGraph* gr_H0lll0b = new TGraphErrors(MassPoints_Hbb.size(), &(MassPoints_Hbb[0]), &expectedH0lll0b[0]);
  // TGraph* gr_Hcomb   = new TGraphErrors(MassPoints_Hbb.size(), &(MassPoints_Hbb[0]), &expectedHcomb[0]);

  gr_theo->SetLineWidth(2);
  gr_HVT_A->SetLineWidth(2);
  gr_HVT_B->SetLineWidth(2);
  gr_Hbb->SetLineWidth(2);
  gr_Hbb0b->SetLineWidth(2);
  gr_ZeeH0b->SetLineWidth(2);
  // gr_H0lll2b->SetLineWidth(2);
  // gr_H0lll0b->SetLineWidth(2);
  // gr_Hcomb->SetLineWidth(2);

  // std::string mode = "CB";
  std::string mode = "Exp_2";
  std::string Path_ANALYSIS = std::getenv("CMSSW_BASE"); Path_ANALYSIS += "/src/UHH2/VHResonances/Analysis/";
  std::string studies = "nominal";
  // std::string studies = "noSignalFlatUncertainty";

  std::string Method = "AsymptoticLimits";
  // std::string extraOptionsText = "Asimov";
  // std::string extraOptionsText = "AsimovNoSys";
  // std::string extraOptionsText = "Expected";
  // std::string extraOptionsText = "ExpectedSys0";
  // std::string extraOptionsText = "ExpectedSys1";
  // std::string extraOptionsText = "ExpectedSys3";
  // std::string extraOptionsText = "ExpectedSys5";
  std::string extraOptionsText = "ExpectedSys9";
  // std::string extraOptionsText = "ExpectedNoSys";
  // std::string extraOptionsText = "ExpectedNoSys0";
  std::string AnalysisDir = Path_ANALYSIS+"Limits/"+studies+"/";
  if (isHbb)  AnalysisDir = Path_ANALYSIS+"Limits/"+studies+"/Hbb/";
  // std::vector<std::string> collections = {"Puppi", kRed+1},
  std::vector<std::string> collections = {"Puppi"};
  // std::vector<std::string> channels = {"muonchannel", "electronchannel"};
  std::vector<std::string> channels = {"muonchannel", "electronchannel", "leptonchannel"};
  // std::vector<std::string> channels = {"leptonchannel"};
  // std::vector<std::string> years = {"2016", "2017", "2018", "RunII", "fullRunII"};
  std::vector<std::string> years = {"RunII"};
  // std::vector<std::string> years = {"2016"};
  // std::vector<std::string> channels = {"muonchannel"};
  // std::vector<std::string> histFolders = {"btag_DeepBoosted_H4qvsQCD", "NN", "tau21","tau31", "tau41", "tau32", "tau42", "tau43"};
  // std::vector<std::string> histFolders = {"btag_DeepBoosted_H4qvsQCD"};
  // std::vector<std::string> histFolders = {"btag_DeepBoosted_H4qvsQCD", "btag_DeepBoosted_H4qvsQCDp2", "btag_DeepBoosted_H4qvsQCDp02", "btag_DeepBoosted_H4qvsQCDpt1000", "btag_DeepBoosted_H4qvsQCDpt1000p2", "btag_DeepBoosted_H4qvsQCDpt1000p02"};
  // std::vector<std::string> histFolders = {"btag_DeepBoosted_H4qvsQCDptdep", "btag_DeepBoosted_H4qvsQCDp02",
  // "btag_DeepBoosted_H4qvsQCDptdep_x3", "btag_DeepBoosted_H4qvsQCDptdep_x2x3", "btag_DeepBoosted_H4qvsQCDptdep_x1x3", "btag_DeepBoosted_H4qvsQCDmassdep_x3",
  // "btag_DeepBoosted_H4qvsQCDmassdep2_x3", "btag_DeepBoosted_H4qvsQCDmassdep_x2x3", "btag_DeepBoosted_H4qvsQCDmassdep_x1x3", "btag_DeepBoosted_H4qvsQCDmassdep_x1x2" };
  // std::vector<std::string> histFolders = {"btag_DeepBoosted_H4qvsQCDmassdep_x3","btag_DeepBoosted_H4qvsQCDptdep",
  // "btag_DeepBoosted_H4qvsQCDptdep_x3", "btag_DeepBoosted_H4qvsQCDptdep_x2x3", "btag_DeepBoosted_H4qvsQCDptdep_x1x3", "btag_DeepBoosted_H4qvsQCDmassdep_x3",
  // "btag_DeepBoosted_H4qvsQCDmassdep2_x3", "btag_DeepBoosted_H4qvsQCDmassdep_x2x3", "btag_DeepBoosted_H4qvsQCDmassdep_x1x3", "btag_DeepBoosted_H4qvsQCDmassdep_x1x2" };
  // std::vector<std::string> histFolders = {"btag_DeepBoosted_H4qvsQCDmassdep_x3"};
  // std::vector<std::string> histFolders = {"btag_DeepBoosted_H4qvsQCD", "btag_DeepBoosted_H4qvsQCDmassdep", "btag_DeepBoosted_H4qvsQCDmassdep_2", "btag_DeepBoosted_H4qvsQCDmassdep_3", "btag_DeepBoosted_H4qvsQCDmassdep_bb", "btag_DeepBoosted_H4qvsQCDmassdep_cc", "btag_DeepBoosted_H4qvsQCDmassdep_gg", "btag_DeepBoosted_H4qvsQCDmassdep_cc_2", "btag_DeepBoosted_H4qvsQCDmassdep_gg_2", "btag_DeepBoosted_H4qvsQCDmassdep_cc_3", "btag_DeepBoosted_H4qvsQCDmassdep_gg_3", "tau42"};
  // std::vector<std::string> histFolders = {"btag_DeepBoosted_H4qvsQCD", "btag_DeepBoosted_H4qvsQCDmassdep", "btag_DeepBoosted_H4qvsQCDmassdep_2", "btag_DeepBoosted_H4qvsQCDmassdep_3", "btag_DeepBoosted_H4qvsQCDmassdep_cc", "btag_DeepBoosted_H4qvsQCDmassdep_cc_2", "btag_DeepBoosted_H4qvsQCDmassdep_gg_2", "btag_DeepBoosted_H4qvsQCDmassdep_cc_3", "tau42"};
  // std::vector<std::string> histFolders = {"btag_DeepBoosted_H4qvsQCD", "btag_DeepBoosted_H4qvsQCDmassdep", "btag_DeepBoosted_H4qvsQCDmassdep_cc", "tau42"};
  // std::vector<std::string> histFolders = {"btag_DeepBoosted_H4qvsQCDmassdep_3", "btag_DeepBoosted_H4qvsQCDmassdep_cc", "btag_DeepBoosted_H4qvsQCDmassdep_bb", "btag_DeepBoosted_H4qvsQCDmassdep_gg"};
  // std::vector<std::string> histFolders = {"btag_DeepBoosted_H4qvsQCDmassdep_3","btag_DeepBoosted_H4qvsQCDmassdep_cc2_3","btag_DeepBoosted_H4qvsQCDmassdep_cc_3","btag_DeepBoosted_H4qvsQCDmassdep","btag_DeepBoosted_H4qvsQCDmassdep_cc2","btag_DeepBoosted_H4qvsQCDmassdep_cc","btag_DeepBoosted_H4qvsQCDmassdep_2","btag_DeepBoosted_H4qvsQCDmassdep_cc2_2","btag_DeepBoosted_H4qvsQCDmassdep_cc_2",

  // std::vector<std::string> histFolders = {"btag_DeepBoosted_H4qvsQCDmassdep_cc", "btag_DeepBoosted_H4qvsQCDmassdep_cc2_2",
  // "btag_DeepBoosted_H4qvsQCDmassdep_cc_2", "btag_DeepBoosted_H4qvsQCDmassdep_2"};

  // std::vector<std::string> histFolders = { "btag_DeepBoosted_H4qvsQCDmassdep_cc_2", "btag_DeepBoosted_H4qvsQCDmassdep_3"};
  std::vector<std::string> histFolders = {"btag_DeepBoosted_H4qvsQCDmassdep", "btag_DeepBoosted_H4qvsQCDmassdep_cc", "btag_DeepBoosted_H4qvsQCDmassdep_cc1", "btag_DeepBoosted_H4qvsQCDmassdep_cc2", "btag_DeepBoosted_H4qvsQCDmassdep_ccMD"};
  // std::vector<std::string> histFolders = {"btag_DeepBoosted_H4qvsQCDmassdep", "btag_DeepBoosted_H4qvsQCDmassdep_cc"};
  // std::vector<std::string> histFolders = {"btag_DeepBoosted_H4qvsQCDmassdep_cc"};


  // "btag_DeepBoosted_H4qvsQCDmassdep_bb","btag_DeepBoosted_H4qvsQCDmassdep_bb_2","btag_DeepBoosted_H4qvsQCDmassdep_gg_3","btag_DeepBoosted_H4qvsQCDmassdep_gg","btag_DeepBoosted_H4qvsQCDmassdep_gg_2"};
  // std::vector<std::string> histFolders = {"btag_DeepBoosted_H4qvsQCD", "btag_DeepBoosted_H4qvsQCDmassdep", "btag_DeepBoosted_H4qvsQCDmassdep_2", "btag_DeepBoosted_H4qvsQCDmassdep_3"};
  // std::vector<std::string> histFolders = {"btag_DeepBoosted_H4qvsQCDmassdep_cc"};
  // std::vector<std::string> histFolders = {"btag_DeepBoosted_H4qvsQCDmassdep_x3_3","btag_DeepBoosted_H4qvsQCDmassdep_x3"};
  // std::vector<std::string> histFolders = {"btag_DeepBoosted_H4qvsQCDmassdep_old_x3", "btag_DeepBoosted_H4qvsQCD", "btag_DeepBoosted_H4qvsQCDp2", "btag_DeepBoosted_H4qvsQCDp02", "btag_DeepBoosted_H4qvsQCDptdep_old_x3", "btag_DeepBoosted_H4qvsQCDmassdep_x3", "btag_DeepBoosted_H4qvsQCDmassdep_x2x3", "btag_DeepBoosted_H4qvsQCDmassdep_x1x3"};
  // // if (isHbb) histFolders = {"btag_DeepBoosted_HbbvsQCD", "btag_DeepBoosted_probHbb", "tau42", "tau32", "tau21" };

  std::unordered_map<std::string, std::unordered_map<std::string, std::vector<double>>> Limits;

  for (std::string year: years) {
    for (std::string collection: collections) {
      for (std::string channel: channels) {
        double BR = 0.1;
        for (std::string histFolder : histFolders) {
          std::string workingDir = year+"/"+collection+"/"+channel+"/"+histFolder+"/";

          double pbTofb = 1000.;
          // double xsec_ref_ = 0.1; // 0.1 comes from normalising the signal strenght
          double xsec_ref_ = (xsec_ref.find(histFolder) != xsec_ref.end())? xsec_ref.at(histFolder): xsec_ref.at("default_value"); // default_value = 1
          // std::cout << BR << std::endl;
          double norm = xsec_ref_*pbTofb/BR;

          Limits[workingDir]["obs"] = std::vector<double>(nPoints,  norm);
          Limits[workingDir]["xsec"] = std::vector<double>(nPoints, norm);
          Limits[workingDir]["xsecNeg1"] = std::vector<double>(nPoints, norm);
          Limits[workingDir]["xsecNeg2"] = std::vector<double>(nPoints, norm);
          Limits[workingDir]["xsecPos1"] = std::vector<double>(nPoints, norm);
          Limits[workingDir]["xsecPos2"] = std::vector<double>(nPoints, norm);
          Limits[workingDir]["expNeg1"] = std::vector<double>(nPoints, 0);
          Limits[workingDir]["expNeg2"] = std::vector<double>(nPoints, 0);
          Limits[workingDir]["expPos1"] = std::vector<double>(nPoints, 0);
          Limits[workingDir]["expPos2"] = std::vector<double>(nPoints, 0);

          // Extract values for limits from txt file TODO Implement from root file
          for (unsigned int i=0; i<MyMassPoints.size(); ++i) {
            std::string filename = AnalysisDir+workingDir+"/datacards/out/DataCard_"+year+"_"+collection+"_"+channel+"_"+histFolder+"_M"+std::to_string((int)MyMassPoints[i])+"_"+mode+"_"+extraOptionsText+"_"+Method+".out";
            std::ifstream txtfile(filename.c_str(), std::ios::in);

            std::string line;
            bool found= false;
            int count = 0;

            while (count<1000 && !found && !txtfile.eof() && !gSystem->AccessPathName(filename.c_str())) {
              getline(txtfile, line);
              if (line.find("-- AsymptoticLimits ( CLs ) --")!=std::string::npos) found=true;
              count++;
            }

            if (count<1000) { // TODO
              getline(txtfile, line);
              if (line.find("Observed")!=std::string::npos) {
                Limits[workingDir]["obs"][i]    = atof(line.substr(line.find("<")+1).c_str())*norm; getline(txtfile, line);
              }
              Limits[workingDir]["xsecNeg2"][i] = atof(line.substr(line.find("<")+1).c_str())*norm; getline(txtfile, line);
              Limits[workingDir]["xsecNeg1"][i] = atof(line.substr(line.find("<")+1).c_str())*norm; getline(txtfile, line);
              Limits[workingDir]["xsec"][i]     = atof(line.substr(line.find("<")+1).c_str())*norm; getline(txtfile, line);
              Limits[workingDir]["xsecPos1"][i] = atof(line.substr(line.find("<")+1).c_str())*norm; getline(txtfile, line);
              Limits[workingDir]["xsecPos2"][i] = atof(line.substr(line.find("<")+1).c_str())*norm; getline(txtfile, line);


              Limits[workingDir]["expNeg2"][i]=Limits[workingDir]["xsec"][i]-Limits[workingDir]["xsecNeg2"][i];
              Limits[workingDir]["expNeg1"][i]=Limits[workingDir]["xsec"][i]-Limits[workingDir]["xsecNeg1"][i];
              Limits[workingDir]["expPos1"][i]=Limits[workingDir]["xsecPos1"][i]-Limits[workingDir]["xsec"][i];
              Limits[workingDir]["expPos2"][i]=Limits[workingDir]["xsecPos2"][i]-Limits[workingDir]["xsec"][i];

            }

            txtfile.close();
          }

          for (auto [name,lims]: Limits[workingDir] ) {
            if (name!="xsec") continue;
            std::cout << workingDir << " " << name << std::endl;
            for (size_t i = 0; i < lims.size(); i++) std::cout << lims[i] << ", ";
            std::cout << std::endl;
          }

        }
      }
    }
  }

  TString nameXaxix = "m(Z') [GeV]";
  // TString nameYaxix = "#sigma\(pp#rightarrowX\) #times Br\(Z'#rightarrowZ\(ll\) H\(WW\)\) \(fb\)";
  TString nameYaxix = "#sigma#left(pp#rightarrowX#right) #times Br#left(Z'#rightarrow ZH #right)#left(fb#right)";

  writeExtraText = true;       // if extra text
  extraText  = "Work in progress" ;//"Preliminary";
  lumi_13TeV  = TString::Format("%.1f fb^{-1}", lumi_map.at("RunII").at("lumi_fb"));

  std::unordered_map<std::string, TCanvas*> Canvas_Limits_Comparison;
  std::unordered_map<std::string, TLegend*> Legend_Limits_Comparison;

  Canvas_Limits_Comparison["all"] = tdrCanvas("c_limits_comparison_all", plot_lo, plot_hi, yaxis_lo, yaxis_hi, nameXaxix, nameYaxix);
  Legend_Limits_Comparison["all"] = tdrLeg(0.65,0.5,0.9,0.85, 0.025, 42, kBlack);

  for (auto [name,canvas]: Canvas_Limits_Comparison ) {
    canvas->SetLogy(1);
    canvas->SetGridx(1);
    canvas->SetGridy(1);
    Legend_Limits_Comparison[name]->SetFillColor(kWhite);
    Legend_Limits_Comparison[name]->SetFillStyle(1);
    Legend_Limits_Comparison[name]->SetFillStyle(0);
    tdrDraw(gr_theo,    "C", kFullDotLarge, kRed+1,    kSolid,  kRed+1);
    tdrDraw(gr_HVT_A,   "C", kFullDotLarge, kViolet-1,   kSolid,  kViolet-1);
    tdrDraw(gr_HVT_B,   "C", kFullDotLarge, kViolet-9,kSolid,  kViolet-9);
    tdrDraw(gr_Hbb,     "C", kFullDotLarge, kBlue+1,   kDotted, kBlue+1);
    tdrDraw(gr_Hbb0b,   "C", kFullDotLarge, kBlue+1,   kDashed, kBlue+1);
    tdrDraw(gr_ZeeH0b,  "C", kFullDotLarge, kBlue+1,   kDashDotted, kBlue+1);
    // tdrDraw(gr_H0lll2b, "C", kFullDotLarge, kOrange+1, kSolid, kOrange+1);
    // tdrDraw(gr_H0lll0b, "C", kFullDotLarge, kOrange+1, kSolid, kOrange+1);
    // tdrDraw(gr_Hcomb,   "C", kFullDotLarge, kGreen+1,  kSolid, kGreen+1);

    TLegend *leg_ExtraPlot = tdrLeg(0.45,0.7,0.60,0.85, 0.03, 42, kBlack);
    canvas->cd();
    leg_ExtraPlot->AddEntry(gr_theo,    "Prediction", "l");
    leg_ExtraPlot->AddEntry(gr_HVT_A,   "HVT model A", "l");
    leg_ExtraPlot->AddEntry(gr_HVT_B,   "HVT model B", "l");
    leg_ExtraPlot->AddEntry(gr_Hbb,     "Hbb", "l");
    leg_ExtraPlot->AddEntry(gr_Hbb0b,   "Hbb cat 0b", "l");
    leg_ExtraPlot->AddEntry(gr_ZeeH0b,  "Zee cat 0b", "l");
    // leg_ExtraPlot->AddEntry(gr_H0lll2b, "Hbb 0l2l2b", "l");
    // leg_ExtraPlot->AddEntry(gr_H0lll0b, "Hbb 0l2l0b", "l");
    // leg_ExtraPlot->AddEntry(gr_Hcomb,   "Hbb comb", "l");

    leg_ExtraPlot->Draw();
  }

  for (std::string year: years) {
    for (std::string collection: collections) {
      for (std::string channel: channels) {
        writeExtraText = true;       // if extra text
        extraText  = "Work in progress" ;//"Preliminary";
        lumi_13TeV  = TString::Format("%.1f fb^{-1}", lumi_map.at((year=="fullRunII")?"RunII":year).at("lumi_fb"));

        // TCanvas* c_xsec_modecomparison = tdrCanvas("c_xsec_modecomparison", plot_lo, plot_hi, yaxis_lo, yaxis_hi, nameXaxix, nameYaxix);
        // TCanvas* c_xsec_modecomparison = tdrCanvas("c_xsec_modecomparison", plot_lo, plot_hi, 0.01, 10.3, nameXaxix, nameYaxix);
        TCanvas* c_xsec_modecomparison = tdrCanvas("c_xsec_modecomparison", 1000, plot_hi, 0.1, 500, nameXaxix, nameYaxix);
        c_xsec_modecomparison->SetLogy(1);
        c_xsec_modecomparison->SetGridx(1);
        c_xsec_modecomparison->SetGridy(1);

        TLegend *leg_modecomparison = tdrLeg(0.45,0.65,0.9,0.85, 0.02, 42, kBlack);
        leg_modecomparison->SetFillStyle(1); leg_modecomparison->SetFillColor(kWhite);
        leg_modecomparison->SetFillStyle(0);

        tdrDraw(gr_theo,    "C", kFullDotLarge, kRed+1,    kSolid,  kRed+1);
        tdrDraw(gr_HVT_A,   "C", kFullDotLarge, kViolet-1,   kSolid,  kViolet-1);
        tdrDraw(gr_HVT_B,   "C", kFullDotLarge, kViolet-9,kSolid,  kViolet-9);
        tdrDraw(gr_Hbb,     "C", kFullDotLarge, kBlue+1,   kDotted, kBlue+1);
        tdrDraw(gr_Hbb0b,   "C", kFullDotLarge, kBlue+1,   kDashed, kBlue+1);
        tdrDraw(gr_ZeeH0b,  "C", kFullDotLarge, kBlue+1,   kDashDotted, kBlue+1);
        // tdrDraw(gr_H0lll2b, "C", kFullDotLarge, kOrange+1, kSolid, kOrange+1);
        // tdrDraw(gr_H0lll0b, "C", kFullDotLarge, kOrange+1, kSolid, kOrange+1);
        // tdrDraw(gr_Hcomb,   "C", kFullDotLarge, kGreen+1,  kSolid, kGreen+1);
        leg_modecomparison->AddEntry(gr_theo,    "Prediction", "l");
        leg_modecomparison->AddEntry(gr_HVT_A,   "HVT model A", "l");
        leg_modecomparison->AddEntry(gr_HVT_B,   "HVT model B", "l");
        leg_modecomparison->AddEntry(gr_Hbb,     "Hbb", "l");
        leg_modecomparison->AddEntry(gr_Hbb0b,   "Hbb cat 0b", "l");
        leg_modecomparison->AddEntry(gr_ZeeH0b,  "Zee cat 0b", "l");

        std::string dir_modecomparison = AnalysisDir+"/"+year+"/"+collection+"/"+channel+"/";

        TGraphAsymmErrors *g_xsec_mode_ref;

        for (std::string histFolder : histFolders) {
          std::string workingDir = year+"/"+collection+"/"+channel+"/"+histFolder+"/";

          TCanvas* c_xsec = tdrCanvas(("xsec"+workingDir).c_str(), plot_lo, plot_hi, yaxis_lo, yaxis_hi, nameXaxix, nameYaxix);
          c_xsec->SetLogy(1);
          c_xsec->SetGridx(1);
          c_xsec->SetGridy(1);
          c_xsec->cd();

          TGraph *g_obs = new TGraph(nPoints, &(MyMassPoints[0]), &(Limits[workingDir]["obs"][0]));
          std::vector<double> dummy(nPoints,0);
          TGraphAsymmErrors *g_xsec_1sigma = new TGraphAsymmErrors(nPoints, &(MyMassPoints[0]), &(Limits[workingDir]["xsec"][0]), &(dummy[0]), &(dummy[0]), &(Limits[workingDir]["expNeg1"][0]), &(Limits[workingDir]["expPos1"][0]));
          TGraphAsymmErrors *g_xsec_2sigma = new TGraphAsymmErrors(nPoints, &(MyMassPoints[0]), &(Limits[workingDir]["xsec"][0]), &(dummy[0]), &(dummy[0]), &(Limits[workingDir]["expNeg2"][0]), &(Limits[workingDir]["expPos2"][0]));
          TGraphAsymmErrors *g_xsec = new TGraphAsymmErrors(nPoints, &(MyMassPoints[0]), &(Limits[workingDir]["xsec"][0]));
          TGraphAsymmErrors *g_xsec_mode = new TGraphAsymmErrors(nPoints, &(MyMassPoints[0]), &(Limits[workingDir]["xsec"][0]));;
          g_xsec->SetLineWidth(2);
          g_xsec_mode->SetLineWidth(2);

          tdrDraw(g_xsec_2sigma,  "3",  20, kOrange, kSolid, kOrange, 1000, kOrange);
          tdrDraw(g_xsec_1sigma,  "3",  20, kGreen+1, kSolid, kGreen+1, 1000, kGreen+1);
          tdrDraw(g_xsec,  "L",  20, kBlack, 2, kBlack, 1000, kBlack);
          if (doObs) tdrDraw(g_obs,  "LP", kFullDotLarge, kBlack, kSolid, kBlack);
          tdrDraw(gr_theo,    "C", kFullDotLarge, kRed+1,    kDotted,  kRed+1);
          tdrDraw(gr_HVT_A,   "C", kFullDotLarge, kViolet-1, kDotted,  kViolet-1);
          tdrDraw(gr_HVT_B,   "C", kFullDotLarge, kViolet-9, kDotted,  kViolet-9);
          tdrDraw(gr_Hbb,     "C", kFullDotLarge, kBlue+1,   kDotted, kBlue+1);
          tdrDraw(gr_Hbb0b,   "C", kFullDotLarge, kBlue+1,   kDashed, kBlue+1);
          tdrDraw(gr_ZeeH0b,  "C", kFullDotLarge, kBlue+1,   kDashDotted, kBlue+1);
          // tdrDraw(gr_H0lll2b, "C", kFullDotLarge, kOrange+1, kSolid, kOrange+1);
          // tdrDraw(gr_H0lll0b, "C", kFullDotLarge, kOrange+1, kSolid, kOrange+1);
          // tdrDraw(gr_Hcomb,   "C", kFullDotLarge, kGreen+1,  kSolid, kGreen+1);

          // g_obs->SetLineWidth(2);
          // g_obs->SetMarkerStyle(20);
          //g_obs->Draw("LP SAME");

          // TLegend *leg = tdrLeg(0.65,0.45,0.9,0.7, 0.03, 42, kBlack);
          TLegend *leg = tdrLeg(0.65,0.7,0.9,0.85, 0.03, 42, kBlack);
          leg->SetFillStyle(1); leg->SetFillColor(kWhite);
          leg->SetFillStyle(0);
          tdrHeader(leg,"95\% CL Upper Limit", 12);
          if (doObs) leg->AddEntry(g_obs, "Observed", "LP");
          leg->AddEntry(g_xsec, "Expected", "L");
          leg->AddEntry(g_xsec_1sigma, "68\% Expected", "CF");
          leg->AddEntry(g_xsec_2sigma, "95\% Expected", "CF");
          leg->Draw();

          // TLegend *leg_ExtraPlot = tdrLeg(0.65,0.7,0.9,0.9, 0.03, 42, kBlack);
          TLegend *leg_ExtraPlot = tdrLeg(0.45,0.7,0.60,0.85, 0.03, 42, kBlack);
          leg_ExtraPlot->AddEntry(gr_theo,    "Prediction", "l");
          leg_ExtraPlot->AddEntry(gr_HVT_A,   "HVT model A", "l");
          leg_ExtraPlot->AddEntry(gr_HVT_B,   "HVT model B", "l");
          leg_ExtraPlot->AddEntry(gr_Hbb,     "Hbb", "l");
          leg_ExtraPlot->AddEntry(gr_Hbb0b,   "Hbb cat 0b", "l");
          leg_ExtraPlot->AddEntry(gr_ZeeH0b,  "Zee cat 0b", "l");
          // leg_ExtraPlot->AddEntry(gr_H0lll2b, "Hbb 0l2l2b", "l");
          // leg_ExtraPlot->AddEntry(gr_H0lll0b, "Hbb 0l2l0b", "l");
          // leg_ExtraPlot->AddEntry(gr_Hcomb,   "Hbb comb", "l");
          leg_ExtraPlot->Draw();

          c_xsec->SaveAs((dir_modecomparison+histFolder+"/"+"UpperLimit_"+histFolder+"_"+extraOptionsText+".pdf").c_str());
          c_xsec->SaveAs((dir_modecomparison+histFolder+"/"+"UpperLimit_"+histFolder+"_"+extraOptionsText+".root").c_str());
          c_xsec->SaveAs((dir_modecomparison+histFolder+"/"+"UpperLimit_"+histFolder+"_"+extraOptionsText+".C").c_str());

          // TFile *file=new TFile(Form("UpperLimits_xsec_%s.root", background.c_str()), "RECREATE");
          // g_obs->Write("g_obs");
          // g_xsec->Write("g_xsec");
          // g_xsec_1sigma->Write("g_xsec_1sigma");
          // g_xsec_2sigma->Write("g_xsec_2sigma");
          // file->Close();
          c_xsec_modecomparison->cd();
          // if (histFolder=="btag_DeepBoosted_H4qvsQCDptdep") g_xsec_mode_ref = new TGraphAsymmErrors(nPoints, &(MyMassPoints[0]), &(Limits[workingDir]["xsec"][0]));
          // if (histFolder=="btag_DeepBoosted_H4qvsQCDmassdep_x3") g_xsec_mode_ref = new TGraphAsymmErrors(nPoints, &(MyMassPoints[0]), &(Limits[workingDir]["xsec"][0]));
          // if (histFolder=="btag_DeepBoosted_H4qvsQCDmassdep_x3") g_xsec_mode_ref = new TGraphAsymmErrors(nPoints, &(MyMassPoints[0]), &(Limits[workingDir]["xsec"][0]));
          if (histFolder==histFolders[0]) g_xsec_mode_ref = new TGraphAsymmErrors(nPoints, &(MyMassPoints[0]), &(Limits[workingDir]["xsec"][0]));
          // if (histFolder=="btag_DeepBoosted_H4qvsQCDptdep") continue;
          TGraphAsymmErrors* g_xsec_mode_scale  = new TGraphAsymmErrors(nPoints, &(MyMassPoints[0]), &(Limits[workingDir]["xsec"][0]));
          for (size_t i = 0; i < nPoints; i++) {
            double x_num, y_num, x_den, y_den;
            g_xsec_mode_scale->GetPoint(i,x_num, y_num);
            g_xsec_mode_ref->GetPoint(i,x_den, y_den);
            // g_xsec_mode_scale->SetPoint(i,x_num, y_num/y_den);
            g_xsec_mode_scale->SetPoint(i,x_num, y_num);
            g_xsec_mode_scale->GetPoint(i,x_num, y_num);
          }
          g_xsec_mode_scale->SetLineWidth(2);
          tdrDraw(g_xsec_mode_scale,  "L",  20, colors[histFolder], kSolid, colors[histFolder], 1000, colors[histFolder]);
          leg_modecomparison->AddEntry(g_xsec_mode_scale, histFolder.c_str(), "L");
          // tdrDraw(g_xsec_mode,  "L",  20, colors[histFolder], 2, colors[histFolder], 1000, colors[histFolder]);
          // leg_modecomparison->AddEntry(g_xsec_mode, histFolder.c_str(), "L");

          // if (histFolder=="btag_DeepBoosted_H4qvsQCDmassdep_x3") {
          if (histFolder==histFolders[0]) {
            Canvas_Limits_Comparison["all"]->cd();
            tdrDraw(g_xsec,  "L",  20, colors[year], (channel=="muonchannel")? kSolid: ((channel=="electronchannel")? kDotted: kDashed), colors[year], 1000, colors[year]);
            // Legend_Limits_Comparison["all"]->AddEntry(g_xsec, TString(workingDir).ReplaceAll("Puppi","").ReplaceAll("btag_DeepBoosted_H4qvsQCDmassdep_x3","").ReplaceAll("channel","").ReplaceAll("/"," "), "L");
            Legend_Limits_Comparison["all"]->AddEntry(g_xsec, TString(workingDir).ReplaceAll("Puppi","").ReplaceAll(histFolders[0],"").ReplaceAll("channel","").ReplaceAll("/"," "), "L");
          }
        }
        c_xsec_modecomparison->cd();
        // tdrDraw(gr_theo,  "C", kFullDotLarge, kRed+1, kSolid, kRed+1);
        // tdrDraw(gr_Hbb,  "C", kFullDotLarge, kBlue+1, kDotted, kBlue+1);
        // tdrDraw(gr_Hbb0b,  "C", kFullDotLarge, kBlue+1, kDashed, kBlue+1);
        // leg_modecomparison->AddEntry(gr_theo, "Prediction", "l");
        // leg_modecomparison->AddEntry(gr_Hbb, "Hbb", "l");
        // leg_modecomparison->AddEntry(gr_Hbb0b, "Hbb cat 0b", "l");
        // leg_modecomparison->Draw("");
        c_xsec_modecomparison->SaveAs((dir_modecomparison+"UpperLimit_modecomparison_"+extraOptionsText+".pdf").c_str());
        c_xsec_modecomparison->SaveAs((dir_modecomparison+"UpperLimit_modecomparison_"+extraOptionsText+".root").c_str());
        c_xsec_modecomparison->SaveAs((dir_modecomparison+"UpperLimit_modecomparison_"+extraOptionsText+".C").c_str());
      }
    }
  }
  Canvas_Limits_Comparison["all"]->cd();
  tdrDraw(gr_theo,    "C", kFullDotLarge, kRed+1,    kSolid,  kRed+1);
  tdrDraw(gr_HVT_A,   "C", kFullDotLarge, kViolet-1,   kSolid,  kViolet-1);
  tdrDraw(gr_HVT_B,   "C", kFullDotLarge, kViolet-9,kSolid,  kViolet-9);
  tdrDraw(gr_Hbb,     "C", kFullDotLarge, kBlue+1,   kDotted, kBlue+1);
  tdrDraw(gr_Hbb0b,   "C", kFullDotLarge, kBlue+1,   kDashed, kBlue+1);
  tdrDraw(gr_ZeeH0b,  "C", kFullDotLarge, kBlue+1,   kDashDotted, kBlue+1);
  // tdrDraw(gr_H0lll2b, "C", kFullDotLarge, kOrange+1, kSolid, kOrange+1);
  // tdrDraw(gr_H0lll0b, "C", kFullDotLarge, kOrange+1, kSolid, kOrange+1);
  // tdrDraw(gr_Hcomb,   "C", kFullDotLarge, kGreen+1,  kSolid, kGreen+1);

  Legend_Limits_Comparison["all"]->AddEntry(gr_theo,    "Prediction", "l");
  Legend_Limits_Comparison["all"]->AddEntry(gr_HVT_A,   "HVT model A", "l");
  Legend_Limits_Comparison["all"]->AddEntry(gr_HVT_B,   "HVT model B", "l");
  Legend_Limits_Comparison["all"]->AddEntry(gr_Hbb,     "Hbb", "l");
  Legend_Limits_Comparison["all"]->AddEntry(gr_Hbb0b,   "Hbb cat 0b", "l");
  Legend_Limits_Comparison["all"]->AddEntry(gr_ZeeH0b,  "Zee cat 0b", "l");
  // Legend_Limits_Comparison["all"]->AddEntry(gr_H0lll2b, "Hbb 0l2l2b", "l");
  // Legend_Limits_Comparison["all"]->AddEntry(gr_H0lll0b, "Hbb 0l2l0b", "l");
  // Legend_Limits_Comparison["all"]->AddEntry(gr_Hcomb,   "Hbb comb", "l");
  Legend_Limits_Comparison["all"]->Draw("");
  Canvas_Limits_Comparison["all"]->SaveAs((AnalysisDir+"UpperLimit_comparison_"+extraOptionsText+".pdf").c_str());
  Canvas_Limits_Comparison["all"]->SaveAs((AnalysisDir+"UpperLimit_comparison_"+extraOptionsText+".root").c_str());
  Canvas_Limits_Comparison["all"]->SaveAs((AnalysisDir+"UpperLimit_comparison_"+extraOptionsText+".C").c_str());

  std::unordered_map<std::string, TCanvas*> Canvas_Limits_Ratio;
  std::unordered_map<std::string, TLegend*> Legend_Limits_Ratio;

  std::string collection = "Puppi";
  std::string channel = "muonchannel";
  std::string histFolder = histFolders[0];
  // std::string histFolder = "btag_DeepBoosted_H4qvsQCDmassdep_x3";
  // std::string histFolder = "btag_DeepBoosted_H4qvsQCDp2";
  // std::string histFolder = "btag_DeepBoosted_H4qvsQCDp02";
  std::string workingDir = collection+"/"+channel+"/"+histFolder+"/";
  Canvas_Limits_Ratio[workingDir] = tdrCanvas("c_limits_ratio", plot_lo, plot_hi, yaxis_lo, 5, nameXaxix, "Ratio Limits");
  Legend_Limits_Ratio[workingDir] = tdrLeg(0.65,0.5,0.9,0.85, 0.025, 42, kBlack);

  for (std::string year: years) {
    std::vector<double> ratio(nPoints,0);
    for (size_t i = 0; i < nPoints; i++) ratio[i] = Limits[year+"/"+workingDir]["xsec"][i]/Limits["RunII/"+workingDir]["xsec"][i];

    TGraph *g_ratio = new TGraph(nPoints, &(MyMassPoints[0]), &(ratio[0]));
    tdrDraw(g_ratio,  "L",  20, colors[year], (channel=="muonchannel")? kSolid: ((channel=="electronchannel")? kDotted: kDashed), colors[year], 1000, colors[year]);
    Legend_Limits_Ratio[workingDir]->AddEntry(g_ratio, year.c_str(), "L");
    double lumiFact = lumi_map.at("RunII").at("lumi_fb")/lumi_map.at((year=="fullRunII")?"RunII":year).at("lumi_fb");
    lumiFact = TMath::Sqrt(lumiFact);
    TLine *line=new TLine(plot_lo, lumiFact, plot_hi, lumiFact); line->SetLineStyle(kDashed); line->SetLineWidth(2); line->SetLineColor(colors[year]); line->Draw("same");
  }
  Legend_Limits_Ratio[workingDir]->Draw("");
  Canvas_Limits_Ratio[workingDir]->SaveAs((AnalysisDir+"UpperLimit_comparison_"+(TString(workingDir).ReplaceAll("/","_"))+"_"+extraOptionsText+".pdf"));
  Canvas_Limits_Ratio[workingDir]->SaveAs((AnalysisDir+"UpperLimit_comparison_"+(TString(workingDir).ReplaceAll("/","_"))+"_"+extraOptionsText+".root"));
  Canvas_Limits_Ratio[workingDir]->SaveAs((AnalysisDir+"UpperLimit_comparison_"+(TString(workingDir).ReplaceAll("/","_"))+"_"+extraOptionsText+".C"));


}



void PlotExpectedLines(TCanvas* canv, TLegend* leg, const std::vector<double> & x_val, const std::vector<double> & y_val, std::string name,  int lineColor=kBlack, int lineStyle=kSolid) {
  TGraph* gr   = new TGraphErrors(x_val.size(), &(x_val[0]), &y_val[0]);
  gr->SetLineWidth(2);
  canv->cd();
  tdrDraw(gr, "C", kFullDotLarge, lineColor, lineStyle, lineColor);
  leg->AddEntry(gr, name.c_str(), "l");
}


void PlotTheoryLines(TCanvas* canv, TLegend* leg) {
  PlotExpectedLines(canv, leg, theo_mass, HVT_A_xsec, "HVT model A", kViolet-8,    kSolid);
  PlotExpectedLines(canv, leg, theo_mass, HVT_B_xsec, "HVT model B", kViolet-1, kSolid);
}


void PlotRefLines(TCanvas* canv, TLegend* leg) {
  PlotExpectedLines(canv, leg, MassPoints_Hbb, expectedHbb,     "B2G-19-006 bb", kBlue+1,   kDotted);
  PlotExpectedLines(canv, leg, MassPoints_Hbb, expectedHbb0b,   "B2G-19-006 0b", kBlue+1,   kDashed);
  // PlotExpectedLines(canv, leg, MassPoints_Hbb, expectedH0lll2b, "Hbb 0l2l2b", kOrange+1, kSolid);
  // PlotExpectedLines(canv, leg, MassPoints_Hbb, expectedH0lll0b, "Hbb 0l2l0b", kOrange+1, kSolid);
  // PlotExpectedLines(canv, leg, MassPoints_Hbb, expectedHcomb,   "Hbb comb",   kGreen+1,  kSolid);

  PlotExpectedLines(canv, leg, ee_mass,        ee_xsec,   "B2G-19-006 0b (Zee)",   kBlue+1,  kDashDotted);

}

void PlotExpectedError(TCanvas* canv, TLegend* leg, const std::vector<double> & x_val, const std::vector<double> & y_val, const std::vector<double> & y_1sigmaPos, const std::vector<double> & y_1sigmaNeg, const std::vector<double> & y_2sigmaPos, const std::vector<double> & y_2sigmaNeg) {
  unsigned int nPoints = x_val.size();

  std::vector<double> dummy(nPoints,0);

  TGraphAsymmErrors *g_xsec_1sigma = new TGraphAsymmErrors(nPoints, &(x_val[0]), &(y_val[0]), &(dummy[0]), &(dummy[0]), &(y_1sigmaNeg[0]), &(y_1sigmaPos[0]));
  TGraphAsymmErrors *g_xsec_2sigma = new TGraphAsymmErrors(nPoints, &(x_val[0]), &(y_val[0]), &(dummy[0]), &(dummy[0]), &(y_2sigmaNeg[0]), &(y_2sigmaPos[0]));
  TGraphAsymmErrors *g_xsec        = new TGraphAsymmErrors(nPoints, &(x_val[0]), &(y_val[0]));

  canv->cd();
  tdrDraw(g_xsec_2sigma,  "3",  20, kOrange,   kSolid, kOrange, 1000, kOrange);
  tdrDraw(g_xsec_1sigma,  "3",  20, kGreen+1,  kSolid, kGreen+1, 1000, kGreen+1);
  tdrDraw(g_xsec,         "L",  20, kBlack, 2, kBlack, 1000, kBlack);
  g_xsec->SetLineWidth(2);

  tdrHeader(leg,"95\% CL Upper Limit", 12);
  leg->AddEntry(g_xsec, "Expected", "L");
  leg->AddEntry(g_xsec_1sigma, "68\% Expected", "CF");
  leg->AddEntry(g_xsec_2sigma, "95\% Expected", "CF");
}


void PlotLimitsFinal(){
  std::string Path_ANALYSIS = std::getenv("CMSSW_BASE"); Path_ANALYSIS += "/src/UHH2/VHResonances/Analysis/";
  std::string AnalysisDir = Path_ANALYSIS+"Limits/nominal/";

  lumi_13TeV  = TString::Format("%.1f fb^{-1}", lumi_map.at("RunII").at("lumi_fb"));

  TCanvas* canv = tdrCanvas("canv_final", 1000, 5000, 1e-01, 700, "m(Z') [GeV]", "#sigma#left(pp#rightarrowX#right) #times Br#left(Z'#rightarrow ZH #right)#left(fb#right)");
  TLegend* leg_theo = tdrLeg(0.45, 0.7, 0.60, 0.85, 0.025, 42, kBlack);
  TLegend* leg_comp = tdrLeg(0.45, 0.6, 0.60, 0.70, 0.025, 42, kBlack);
  TLegend* leg = tdrLeg(0.65, 0.7, 0.9, 0.85, 0.025, 42, kBlack);


  canv->SetLogy(1);

  const std::vector<double> & x_val         = {1000.0,   1200.0,  1400.0, 1600.0,   1800.0, 2000.0, 2500.0, 3000.0, 3500.0, 4000.0, 4500.0, 5000.0, 5500.0, 6000.0, 7000.0, 8000.0};

  const std::vector<double> & y_val_old    = {101.562,  28.828,  16.875,  11.523,  8.867,  7.09,   3.799,  2.109,  1.313,  0.991,  0.884,  0.884,  0.962,  1.025,  1.338,  1.68};

  //bin100
  const std::vector<double> & y_val_cc_1   = {86.875,23.359,14.57,10.078,8.125,107.812,26.484,3.057,2.295,1.846,1.67,1.543,1.572,1.65,1.992,4.297};
  // const std::vector<double> & y_val_cc1_1  = {86.25,25.703,16.016,10.977,8.711,7.188,4.551,3.301,2.422,1.953,1.685,1.504,1.499,1.533,2.168};
  // const std::vector<double> & y_val_cc2_1  = {80.625,21.484,13.672,9.648,7.5,7.93,3.867,2.715,2.031,1.611,41.094,1.406,1.46,1.65,30.078,4.863};
  // const std::vector<double> & y_val_ccMD_1 = {121.25,26.953,17.109,12.109,9.805,8.477,5.801,4.238,3.213,2.461,1.973,1.699,1.631,1.504,1.621,1.934};

  const std::vector<double> & y_val_qq_1   = {10,10,10.039,18.047,14.297,12.148,7.891,5.469,4.004,2.959,2.422,2.09,2.061,1.973,2.393,3.086};
  // const std::vector<double> & y_val_cc   = {86.875,23.359,14.57,10.078,8.125,107.812,26.484,3.057,2.295,1.846,1.67,1.543,1.572,1.65,1.992,4.297};
  const std::vector<double> & y_val_cc1  = {75.312, 20.859, 13.32, 9.336, 7.168, 5.781, 3.525, 2.432, 1.802, 1.47, 1.299, 1.23, 1.27, 1.406, 1.777, 2.1};
  const std::vector<double> & y_val_cc2  = {71.875, 19.922, 12.578, 8.789, 6.699, 5.391, 3.379, 2.441, 1.885, 1.567, 1.416, 1.382, 1.455, 1.65, 2.158, 2.637,};
  const std::vector<double> & y_val_ccMD = {71.875, 21.203, 14.57, 10.508, 8.281, 6.992, 4.473, 2.891, 1.924, 1.372, 1.108, 0.991, 1.006, 1.084, 1.353, 1.67, };

  //bin30
  // const std::vector<double> & y_val_cc2  = {76.562,21.016,13.594,9.883,7.812,6.445,4.238,3.115,2.344,1.802,1.514,1.416,1.465,1.65,2.158,2.637};
  // const std::vector<double> & y_val_ccMD = {140.625,31.875,19.766,14.336,11.367,9.805,6.816,5.059,3.906,3.066,2.568,2.236,2.148,2.129,2.197,2.412};

  const std::vector<double> & y_val_MD_3    = {70,       49.688,  30,      20.547,  14.609,11.328, 6.039, 3.291,  2.09,   1.455,  1.162,  1.006,  1.011,  0.981,  1.177, 1.533};
  const std::vector<double> & y_val_qq      = {35.781,23.828,17.031,13.477,11.523,7.695,0,4.043,3.145,2.617,2.314,2.256,2.363,2.822,12.461};
  // const std::vector<double> & y_val_cc_s1   = {87.812,24.531,15.781,11.445,8.984,7.441,4.805,3.447,2.627,2.148,1.88,1.758,1.777,1.934,2.285,2.578};
  // const std::vector<double> & y_val_cc_s3   = {96.562,26.875,17.266,12.5,9.805,8.086,5.195,3.711,2.822,2.305,2.002,1.875,1.895,2.061,2.441,2.734};
  // const std::vector<double> & y_val_cc_s5   = {109.688,30.625,19.609,14.18,11.094,9.141,5.82,4.141,3.125,2.539,2.207,2.07,2.09,2.266,2.676,2.998};


  const std::vector<double> & y_val_cc_ns   = {70.938, 19.766, 12.539, 8.867, 6.836, 5.566, 3.467, 2.461, 1.865, 1.543, 1.367, 1.299, 1.328, 1.475, 1.875, 2.227};
  const std::vector<double> & y_val_cc_s0   = {70.938, 19.766, 12.539, 8.867, 6.836, 5.566, 3.467, 2.461, 1.865, 1.543, 1.367, 1.299, 1.328, 1.475, 1.875, 2.227, };
  const std::vector<double> & y_val_cc_s1   = {71.875, 20.078, 12.734, 8.984, 6.934, 5.645, 3.516, 2.49, 1.885, 1.558, 1.377, 1.309, 1.343, 1.484, 1.885, 2.246};
  const std::vector<double> & y_val_cc_s3   = {79.062, 21.953, 13.906, 9.805, 7.559, 6.133, 3.799, 2.676, 2.021, 1.66, 1.47, 1.392, 1.426, 1.577, 1.992, 2.383, };
  const std::vector<double> & y_val_cc_s5   = {89.688, 25, 15.781, 11.094, 8.555, 6.914, 4.258, 2.979, 2.236, 1.826, 1.611, 1.523, 1.553, 1.719, 2.168, 2.588, };
  const std::vector<double> & y_val_cc_s9   = {118.125, 32.812, 20.703, 14.531, 11.133, 8.984, 5.469, 3.789, 2.812, 2.275, 1.992, 1.875, 1.904, 2.1, 2.646, 3.145};

  const std::vector<double> & y_val_cc_muo   = {98.438, 29.062, 19.609, 14.219, 11.68, 10, 6.836, 5.059, 3.877, 3.232, 2.871, 2.588, 2.422, 2.334, 2.305, 2.412};
  const std::vector<double> & y_val_cc_ele   = {122.812, 32.344, 19.531, 13.633, 10.195, 8.086, 5, 3.652, 2.91, 2.539, 2.363, 2.441, 2.871, 4.004, 10.664, 37.344};


  const std::vector<double> & y_2sigmaPos   = {87.108, 24.911, 16.186, 11.789, 9.45, 7.93, 5.536, 4.363, 3.606, 3.194, 2.837, 2.689, 2.759, 3.056, 3.87, 4.604};
  const std::vector<double> & y_1sigmaPos   = {36.435, 10.403, 6.751, 4.897, 3.896, 3.262, 2.217, 1.707, 1.416, 1.254, 1.169, 1.15, 1.223, 1.379, 1.807, 2.167};
  const std::vector<double> & y_val         = {73.125, 20.391, 12.93, 9.102, 7.031, 5.723, 3.564, 2.52, 1.909, 1.572, 1.396, 1.323, 1.357, 1.504, 1.904, 2.266};
  const std::vector<double> & y_1sigmaNeg   = {22.852, 6.503, 4.155, 2.97, 2.329, 1.948, 1.291, 0.969, 0.772, 0.665, 0.615, 0.597, 0.627, 0.703, 0.905, 1.083};
  const std::vector<double> & y_2sigmaNeg   = {36.563, 10.275, 6.566, 4.693, 3.68, 3.041, 1.991, 1.477, 1.163, 1.001, 0.916, 0.889, 0.922, 1.034, 1.331, 1.593};



  // PlotExpectedError(canv, leg, x_val, y_val, y_1sigmaPos, y_1sigmaNeg, y_2sigmaPos, y_2sigmaNeg);

  // PlotExpectedLines(canv, leg_comp, x_val, y_val_cc_muo,  "H4q+cc (#mu ch.)", kRed+1, kSolid);
  // PlotExpectedLines(canv, leg_comp, x_val, y_val_cc_ele,  "H4q+cc (e ch.)",   kAzure+1, kSolid);
  // PlotExpectedLines(canv, leg_comp, x_val, y_val_old,  "H4q+cc_old", kOrange-1, kSolid);
  //
  // PlotExpectedLines(canv, leg_comp, x_val, y_val_cc_1,  "H4q+cc_1", kRed+1, kDashed);
  //
  // PlotExpectedLines(canv, leg_comp, x_val, y_val_qq_1, "H4q_1", kBlack, kSolid);
  // PlotExpectedLines(canv, leg_comp, x_val, y_val_qq,   "H4q", kBlack, kDotted);

  // PlotExpectedLines(canv, leg_comp, x_val, y_val,      "H4q+cc", kRed+1, kSolid);
  // PlotExpectedLines(canv, leg_comp, x_val, y_val_cc1,  "H4q+cc_1",   kOrange+1, kSolid);
  // PlotExpectedLines(canv, leg_comp, x_val, y_val_cc2,  "H4q+cc_2",   kGreen+1, kSolid);
  // PlotExpectedLines(canv, leg_comp, x_val, y_val_ccMD, "H4q+cc_MD",  kAzure+1, kSolid);

  // PlotExpectedLines(canv, leg_comp, x_val, y_val_MD_3,    "H4q only", kRed+1, kDashed);

  PlotExpectedLines(canv, leg_comp, x_val, y_val_cc_ns, "H4q+cc_noSys",kBlack, kSolid);
  PlotExpectedLines(canv, leg_comp, x_val, y_val_cc_s0, "H4q+cc_s0",   kViolet+1, kSolid);
  PlotExpectedLines(canv, leg_comp, x_val, y_val,       "H4q+cc_s1.5", kRed+1, kSolid);
  PlotExpectedLines(canv, leg_comp, x_val, y_val_cc_s1, "H4q+cc_s1",   kOrange+1, kSolid);
  PlotExpectedLines(canv, leg_comp, x_val, y_val_cc_s3, "H4q+cc_s3",   kGreen+1, kSolid);
  PlotExpectedLines(canv, leg_comp, x_val, y_val_cc_s5, "H4q+cc_s5",   kAzure+1, kSolid);
  PlotExpectedLines(canv, leg_comp, x_val, y_val_cc_s9, "H4q+cc_s9",   kBlue+1, kSolid);

  // PlotTheoryLines(canv,leg_theo);
  // PlotRefLines(canv,leg_theo);

  // PlotExpectedLines(canv, leg, MassPoints_Hbb, expectedHbbTest,   "H4q_only", kBlue+1,   kDashed);

  canv->Update();
  canv->RedrawAxis();
  canv->SaveAs((AnalysisDir+"UpperLimit_final.pdf").c_str());
  canv->SaveAs((AnalysisDir+"UpperLimit_final.root").c_str());
  canv->SaveAs((AnalysisDir+"UpperLimit_final.C").c_str());

}


int main(){
  PlotLimits();
  PlotLimitsFinal();
  return 0;
}
