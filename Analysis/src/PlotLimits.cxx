#include "PlotLimits.hpp"

void PlotLimits(bool doObs, bool isHbb) {

  // double plot_hi = 2800;
  const double plot_lo = 300;
  const double plot_hi = 8200;

  const double yaxis_lo = 1e-01;
  const double yaxis_hi = 5e03;

  std::map<std::string, int> colors = {
    {"2016",kViolet-9+1}, {"2017",kOrange+1}, {"2018",kCyan+1}, {"RunII",kGreen+2}, {"fullRunII", kGreen+2},
    {"NN",kOrange+1}, {"NN_1",kViolet-9+1}, {"NN_2",kCyan+1}, {"CNN",kBlue+1}, {"tau42",kBlack},
    {"btag_DeepBoosted_probHbb",kViolet-9+1}, {"tau32",kBlue+1}, {"tau21",kOrange+1},
    {"btag_DeepBoosted_H4qvsQCDp2", kOrange+1}, {"btag_DeepBoosted_H4qvsQCDp02", kYellow+1}, {"btag_DeepBoosted_H4qvsQCDpt1000", kCyan+1},
    {"btag_DeepBoosted_H4qvsQCDpt1000p2", kBlue+1}, {"btag_DeepBoosted_H4qvsQCDpt1000p02", kRed+1},
    {"btag_DeepBoosted_H4qvsQCDptdep", kRed+1},
    {"btag_DeepBoosted_H4qvsQCD",               kAzure+7},
    {"btag_DeepBoosted_H4qvsQCD_cc",            kAzure-9},
    {"btag_DeepBoosted_H4qvsQCDmassdep",        kRed+1},
    {"btag_DeepBoosted_H4qvsQCDmassdep_2",      kViolet-9},
    {"btag_DeepBoosted_H4qvsQCDmassdep_3",      kOrange-1},
    {"btag_DeepBoosted_H4qvsQCDmassdep_bb",     kGreen+2},
    {"btag_DeepBoosted_H4qvsQCDmassdep_cc",     kViolet+1},
    {"btag_DeepBoosted_H4qvsQCDmassdep_gg",     kAzure+1},
    {"btag_DeepBoosted_H4qvsQCDmassdep_cc_2",   kViolet-7},
    {"btag_DeepBoosted_H4qvsQCDmassdep_gg_2",   kGreen+2},
    {"btag_DeepBoosted_H4qvsQCDmassdep_cc_3",   kCyan+1},
    {"btag_DeepBoosted_H4qvsQCDmassdep_gg_3",   kBlue+1},
    {"btag_DeepBoosted_H4qvsQCDmassdep_cc2_3",  kOrange+1},
    {"btag_DeepBoosted_H4qvsQCDmassdep_cc2",    kBlue+1},
    {"btag_DeepBoosted_H4qvsQCDmassdep_cc2_2",  kCyan+2},
    {"btag_DeepBoosted_H4qvsQCDmassdep_bb_2",   kViolet-1},
    {"btag_DeepBoosted_H4qvsQCDmassdep_cc1",    kViolet-1},
    {"btag_DeepBoosted_H4qvsQCDmassdep_ccMD",   kCyan-1},
    {"btag_DeepBoosted_H4qvsQCD_ccMD",          kCyan-1},
    {"btag_DeepBoosted_H4qvsQCD_ccMD2",         kCyan-1},
    {"btag_DeepBoosted_H4qvsQCD_ccMD3",         kCyan-1},
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
  // std::string extraOptionsText = "Observed";
  std::string extraOptionsText = "Expected";
  // std::string extraOptionsText = "ExpectedSys0";
  // std::string extraOptionsText = "ExpectedSys1";
  // std::string extraOptionsText = "ExpectedSys3";
  // std::string extraOptionsText = "ExpectedSys5";
  // std::string extraOptionsText = "ExpectedSys9";
  // std::string extraOptionsText = "ExpectedNoSys";
  // std::string extraOptionsText = "ExpectedNoSys0";
  std::string AnalysisDir = Path_ANALYSIS+"Limits/"+studies+"/";
  if (isHbb)  AnalysisDir = Path_ANALYSIS+"Limits/"+studies+"/Hbb/";
  // std::vector<std::string> collections = {"Puppi", kRed+1},
  std::vector<std::string> collections = {"Puppi"};
  // std::vector<std::string> channels = {"muonchannel", "electronchannel"};
  std::vector<std::string> channels = {"muonchannel", "electronchannel", "chargedleptonchannel", "invisiblechannel", "leptonchannel"};
  // std::vector<std::string> channels = {"leptonchannel"};
  // std::vector<std::string> years = {"2016", "2017", "2018", "RunII", "fullRunII"};
  std::vector<std::string> years = {"RunII"};
  // std::vector<std::string> years = {"2016"};
  // std::vector<std::string> channels = {"muonchannel"};


  // std::vector<std::string> histFolders = { "DeepAk8_H4qvsQCD_massdep", "DeepAk8_ZHccvsQCD_MD", "DeepAk8_HccvsQCD_MD", "DeepAk8_H4qvsQCD_MD",
  // "DeepAk8_H4qvsQCD_massdep_HccvsQCD_MD", "DeepAk8_H4qvsQCD", "DeepAk8_HccvsQCD", "DeepAk8_ZHccvsQCD", "DeepAk8_H4qvsQCD_massdep_HccvsQCD",
  // "DeepAk8_H4qvsQCD_massdep_ZHccvsQCD", "DeepAk8_H4qvsQCD_massdep_ZHccvsQCD_MD"};
  // std::vector<std::string> histFolders = {"DeepAk8_ZHccvsQCD_MD"};
  // std::vector<std::string> histFolders = {"DeepAk8_HccvsQCD", "DeepAk8_HccvsQCD2", "DeepAk8_ZHccvsQCD_MD", "DeepAk8_ZHccvsQCD_MD2"};
  // std::vector<std::string> histFolders = {"DeepAk8_HccvsQCD2", "DeepAk8_ZHccvsQCD_MD2"};
  std::vector<std::string> histFolders = {"DeepAk8_ZHccvsQCD_MD"};

  // // if (isHbb) histFolders = {"btag_DeepBoosted_HbbvsQCD", "btag_DeepBoosted_probHbb", "tau42", "tau32", "tau21" };

  std::unordered_map<std::string, std::unordered_map<std::string, std::vector<double>>> Limits;

  for (std::string year: years) {
    for (std::string collection: collections) {
      for (std::string channel: channels) {
        // double BR = 0.1;
        // if (channel=="invisiblechannel"){BR = 0.2;}
        for (std::string histFolder : histFolders) {
          std::string workingDir = year+"/"+collection+"/"+channel+"/"+histFolder+"/";

          double pbTofb = 1000.;
          // double xsec_ref_ = 0.1; // 0.1 comes from normalising the signal strenght
          double xsec_ref_ = (xsec_ref.find(histFolder) != xsec_ref.end())? xsec_ref.at(histFolder): xsec_ref.at("default_value"); // default_value = 1
          // std::cout << BR << std::endl;
          // double norm = xsec_ref_*pbTofb/BR;
          double norm = xsec_ref_*pbTofb;

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
            if (channel!="invisiblechannel") continue;
            // if (name!="xsec") continue;
            std::cout << workingDir << " " << name << std::endl;
            for (size_t i = 0; i < lims.size(); i++) std::cout << lims[i] << ", ";
            std::cout << std::endl;
          }

        }
      }
    }
  }

  TString nameXaxis = "M(Z') [GeV]";
  TString nameYaxis = "#sigma#left(pp#rightarrowZ'#right) #times Br#left(Z'#rightarrow ZH #right)#left(fb#right)";

  writeExtraText = true;       // if extra text
  extraText  = "Work in progress" ;//"Preliminary";
  lumi_13TeV  = TString::Format("%.1f fb^{-1}", lumi_map.at("RunII").at("lumi_fb"));

  std::unordered_map<std::string, TCanvas*> Canvas_Limits_Comparison;
  std::unordered_map<std::string, TLegend*> Legend_Limits_Comparison;

  Canvas_Limits_Comparison["all"] = tdrCanvas("c_limits_comparison_all", plot_lo, plot_hi, yaxis_lo, yaxis_hi, nameXaxis, nameYaxis);
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
    // tdrDraw(gr_Hcomb,   "C", kFullDotLarge, kGreen+2,  kSolid, kGreen+2);

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

        // TCanvas* c_xsec_modecomparison = tdrCanvas("c_xsec_modecomparison", plot_lo, plot_hi, yaxis_lo, yaxis_hi, nameXaxis, nameYaxis);
        // TCanvas* c_xsec_modecomparison = tdrCanvas("c_xsec_modecomparison", plot_lo, plot_hi, 0.01, 10.3, nameXaxis, nameYaxis);
        TCanvas* c_xsec_modecomparison = tdrCanvas("c_xsec_modecomparison", 1000, plot_hi, 0.1, 500, nameXaxis, nameYaxis);
        c_xsec_modecomparison->SetLogy(1);
        c_xsec_modecomparison->SetGridx(1);
        c_xsec_modecomparison->SetGridy(1);

        TLegend *leg_modecomparison = tdrLeg(0.45,0.65,0.9,0.85, 0.02, 42, kBlack);
        leg_modecomparison->SetFillStyle(1); leg_modecomparison->SetFillColor(kWhite);
        leg_modecomparison->SetFillStyle(0);

        tdrDraw(gr_theo,    "C", kFullDotLarge, kRed+1,    kSolid,      kRed+1);
        tdrDraw(gr_HVT_A,   "C", kFullDotLarge, kViolet-1, kSolid,      kViolet-1);
        tdrDraw(gr_HVT_B,   "C", kFullDotLarge, kViolet-9, kSolid,      kViolet-9);
        tdrDraw(gr_Hbb,     "C", kFullDotLarge, kBlue+1,   kDotted,     kBlue+1);
        tdrDraw(gr_Hbb0b,   "C", kFullDotLarge, kBlue+1,   kDashed,     kBlue+1);
        tdrDraw(gr_ZeeH0b,  "C", kFullDotLarge, kBlue+1,   kDashDotted, kBlue+1);
        // tdrDraw(gr_H0lll2b, "C", kFullDotLarge, kOrange+1, kSolid, kOrange+1);
        // tdrDraw(gr_H0lll0b, "C", kFullDotLarge, kOrange+1, kSolid, kOrange+1);
        // tdrDraw(gr_Hcomb,   "C", kFullDotLarge, kGreen+2,  kSolid, kGreen+2);
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

          TCanvas* c_xsec = tdrCanvas(("xsec"+workingDir).c_str(), plot_lo, plot_hi, yaxis_lo, yaxis_hi, nameXaxis, nameYaxis);
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
          tdrDraw(g_xsec_1sigma,  "3",  20, kGreen+2, kSolid, kGreen+2, 1000, kGreen+2);
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
          // tdrDraw(gr_Hcomb,   "C", kFullDotLarge, kGreen+2,  kSolid, kGreen+2);

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

          // Plot lepton channel as comparison
          if (channel=="invisiblechannel"){
            const std::vector<double> & x_val         = {1000.0,   1200.0,  1400.0, 1600.0,   1800.0, 2000.0, 2500.0, 3000.0, 3500.0, 4000.0, 4500.0, 5000.0, 5500.0, 6000.0, 7000.0, 8000.0};
            const std::vector<double> & y_val         = {18.828, 20.391, 12.93, 9.102, 7.031, 5.723, 3.564, 2.52, 1.909, 1.572, 1.396, 1.323, 1.357, 1.504, 1.904, 2.266};
            TGraph* gr   = new TGraphErrors(x_val.size(), &(x_val[0]), &y_val[0]);
            gr->SetLineWidth(2);
            c_xsec->cd();
            tdrDraw(gr, "C", kFullDotLarge, kRed+1, kDashed, kRed+1);
            leg_ExtraPlot->AddEntry(gr,  "lepton channel", "l");
          }

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
  // tdrDraw(gr_Hcomb,   "C", kFullDotLarge, kGreen+2,  kSolid, kGreen+2);

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
  std::string channel = "muonchannel";//TODO
  std::string histFolder = histFolders[0];
  // std::string histFolder = "btag_DeepBoosted_H4qvsQCDmassdep_x3";
  // std::string histFolder = "btag_DeepBoosted_H4qvsQCDp2";
  // std::string histFolder = "btag_DeepBoosted_H4qvsQCDp02";
  std::string workingDir = collection+"/"+channel+"/"+histFolder+"/";
  Canvas_Limits_Ratio[workingDir] = tdrCanvas("c_limits_ratio", plot_lo, plot_hi, yaxis_lo, 5, nameXaxis, "Ratio Limits");
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
  // PlotExpectedLines(canv, leg, theo_mass, HVT_A_xsec, "HVT model A", kRed+2,    kSolid);
  // PlotExpectedLines(canv, leg, theo_mass, HVT_B_xsec, "HVT model B", kAzure+10, kSolid);

  PlotExpectedLines(canv, leg, theo_mass, HVT_A_xsec, "HVT model A", kViolet-9, kSolid);
  PlotExpectedLines(canv, leg, theo_mass, HVT_B_xsec, "HVT model B", kViolet-1, kSolid);
}


void PlotRefLines(TCanvas* canv, TLegend* leg) {
  // arXiv:2102.08198
  // B2G-19-006
  // PlotExpectedLines(canv, leg, MassPoints_Hbb, expectedHbb,     "B2G-19-006 Hbb", kViolet-9, kDashed);
  // PlotExpectedLines(canv, leg, MassPoints_Hbb, expectedHbb0b,   "B2G-19-006 H0b", kViolet-1, kDashed);

  // PlotExpectedLines(canv, leg, MassPoints_Hbb, expectedHbb,     "B2G-19-006 Hbb", kAzure+10, kDashed);
  PlotExpectedLines(canv, leg, MassPoints_Hbb, expectedHbb,     "B2G-19-006 Hbb", kAzure, kDashed);
  PlotExpectedLines(canv, leg, MassPoints_Hbb, expectedHbb0b,   "B2G-19-006 H0b", kRed+2, kDashed);
  // PlotExpectedLines(canv, leg, MassPoints_Hbb, expectedH0lll2b, "Hbb 0l2l2b", kOrange+1, kSolid);
  // PlotExpectedLines(canv, leg, MassPoints_Hbb, expectedH0lll0b, "Hbb 0l2l0b", kOrange+1, kSolid);
  // PlotExpectedLines(canv, leg, MassPoints_Hbb, expectedHcomb,   "Hbb comb",   kGreen+2,  kSolid);

  // PlotExpectedLines(canv, leg, ee_mass,        ee_xsec,   "B2G-19-006 (ZeeH0b)",   kBlue+1,  kDashDotted);

}

void PlotExpectedError(TCanvas* canv, TLegend* leg, const std::vector<double> & x_val, const std::vector<double> & y_val, const std::vector<double> & y_1sigmaPos, const std::vector<double> & y_1sigmaNeg, const std::vector<double> & y_2sigmaPos, const std::vector<double> & y_2sigmaNeg) {
  unsigned int nPoints = x_val.size();

  std::vector<double> dummy(nPoints,0);

  TGraphAsymmErrors *g_xsec_1sigma = new TGraphAsymmErrors(nPoints, &(x_val[0]), &(y_val[0]), &(dummy[0]), &(dummy[0]), &(y_1sigmaNeg[0]), &(y_1sigmaPos[0]));
  TGraphAsymmErrors *g_xsec_2sigma = new TGraphAsymmErrors(nPoints, &(x_val[0]), &(y_val[0]), &(dummy[0]), &(dummy[0]), &(y_2sigmaNeg[0]), &(y_2sigmaPos[0]));
  TGraphAsymmErrors *g_xsec        = new TGraphAsymmErrors(nPoints, &(x_val[0]), &(y_val[0]));

  canv->cd();
  tdrDraw(g_xsec_2sigma,  "3",  20, kOrange,   kSolid, kOrange, 1000, kOrange);
  tdrDraw(g_xsec_1sigma,  "3",  20, kGreen+2,  kSolid, kGreen+2, 1000, kGreen+2);
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

  writeExtraText = true;       // if extra text
  extraText  = "Work in progress" ;//"Preliminary";
  lumi_13TeV  = TString::Format("%.1f fb^{-1}", lumi_map.at("RunII").at("lumi_fb"));

  TCanvas* canv = tdrCanvas("canv_final", 1000, 5200, 1e-01, 700, "M(Z') [GeV]", "#sigma#left(pp#rightarrowZ'#right) #times Br#left(Z'#rightarrow ZH #right)#left(fb#right)");
  TLegend* leg_theo = tdrLeg(0.40, 0.7, 0.60, 0.85, 0.035, 42, kBlack);
  TLegend* leg_comp = tdrLeg(0.40, 0.5, 0.60, 0.70, 0.035, 42, kBlack);
  TLegend* leg = tdrLeg(0.65, 0.7, 0.9, 0.85, 0.035, 42, kBlack);
  canv->SetLogy(1);

  // const std::vector<double> & x_val         = {1000.0,   1200.0,  1400.0, 1600.0,   1800.0, 2000.0, 2500.0, 3000.0, 3500.0, 4000.0, 4500.0, 5000.0, 5500.0, 6000.0, 7000.0, 8000.0};
  //
  // const std::vector<double> & y_val_old    = {101.562,  28.828,  16.875,  11.523,  8.867,  7.09,   3.799,  2.109,  1.313,  0.991,  0.884,  0.884,  0.962,  1.025,  1.338,  1.68};
  //
  // //bin100
  // const std::vector<double> & y_val_cc_1   = {86.875,23.359,14.57,10.078,8.125,107.812,26.484,3.057,2.295,1.846,1.67,1.543,1.572,1.65,1.992,4.297};
  // // const std::vector<double> & y_val_cc1_1  = {86.25,25.703,16.016,10.977,8.711,7.188,4.551,3.301,2.422,1.953,1.685,1.504,1.499,1.533,2.168};
  // // const std::vector<double> & y_val_cc2_1  = {80.625,21.484,13.672,9.648,7.5,7.93,3.867,2.715,2.031,1.611,41.094,1.406,1.46,1.65,30.078,4.863};
  // // const std::vector<double> & y_val_ccMD_1 = {121.25,26.953,17.109,12.109,9.805,8.477,5.801,4.238,3.213,2.461,1.973,1.699,1.631,1.504,1.621,1.934};
  //
  // const std::vector<double> & y_val_qq_1   = {10,10,10.039,18.047,14.297,12.148,7.891,5.469,4.004,2.959,2.422,2.09,2.061,1.973,2.393,3.086};
  // // const std::vector<double> & y_val_cc   = {86.875,23.359,14.57,10.078,8.125,107.812,26.484,3.057,2.295,1.846,1.67,1.543,1.572,1.65,1.992,4.297};
  // const std::vector<double> & y_val_cc1  = {75.312, 20.859, 13.32, 9.336, 7.168, 5.781, 3.525, 2.432, 1.802, 1.47, 1.299, 1.23, 1.27, 1.406, 1.777, 2.1};
  // const std::vector<double> & y_val_cc2  = {71.875, 19.922, 12.578, 8.789, 6.699, 5.391, 3.379, 2.441, 1.885, 1.567, 1.416, 1.382, 1.455, 1.65, 2.158, 2.637,};
  // const std::vector<double> & y_val_ccMD = {71.875, 21.203, 14.57, 10.508, 8.281, 6.992, 4.473, 2.891, 1.924, 1.372, 1.108, 0.991, 1.006, 1.084, 1.353, 1.67, };
  //
  // //bin30
  // // const std::vector<double> & y_val_cc2  = {76.562,21.016,13.594,9.883,7.812,6.445,4.238,3.115,2.344,1.802,1.514,1.416,1.465,1.65,2.158,2.637};
  // // const std::vector<double> & y_val_ccMD = {140.625,31.875,19.766,14.336,11.367,9.805,6.816,5.059,3.906,3.066,2.568,2.236,2.148,2.129,2.197,2.412};
  //
  // const std::vector<double> & y_val_MD_3    = {70,       49.688,  30,      20.547,  14.609,11.328, 6.039, 3.291,  2.09,   1.455,  1.162,  1.006,  1.011,  0.981,  1.177, 1.533};
  // const std::vector<double> & y_val_qq      = {35.781,23.828,17.031,13.477,11.523,7.695,0,4.043,3.145,2.617,2.314,2.256,2.363,2.822,12.461};
  // // const std::vector<double> & y_val_cc_s1   = {87.812,24.531,15.781,11.445,8.984,7.441,4.805,3.447,2.627,2.148,1.88,1.758,1.777,1.934,2.285,2.578};
  // // const std::vector<double> & y_val_cc_s3   = {96.562,26.875,17.266,12.5,9.805,8.086,5.195,3.711,2.822,2.305,2.002,1.875,1.895,2.061,2.441,2.734};
  // // const std::vector<double> & y_val_cc_s5   = {109.688,30.625,19.609,14.18,11.094,9.141,5.82,4.141,3.125,2.539,2.207,2.07,2.09,2.266,2.676,2.998};
  //
  //
  // const std::vector<double> & y_val_cc_ns   = {70.938, 19.766, 12.539, 8.867, 6.836, 5.566, 3.467, 2.461, 1.865, 1.543, 1.367, 1.299, 1.328, 1.475, 1.875, 2.227};
  // const std::vector<double> & y_val_cc_s0   = {70.938, 19.766, 12.539, 8.867, 6.836, 5.566, 3.467, 2.461, 1.865, 1.543, 1.367, 1.299, 1.328, 1.475, 1.875, 2.227, };
  // const std::vector<double> & y_val_cc_s1   = {71.875, 20.078, 12.734, 8.984, 6.934, 5.645, 3.516, 2.49, 1.885, 1.558, 1.377, 1.309, 1.343, 1.484, 1.885, 2.246};
  // const std::vector<double> & y_val_cc_s3   = {79.062, 21.953, 13.906, 9.805, 7.559, 6.133, 3.799, 2.676, 2.021, 1.66, 1.47, 1.392, 1.426, 1.577, 1.992, 2.383, };
  // const std::vector<double> & y_val_cc_s5   = {89.688, 25, 15.781, 11.094, 8.555, 6.914, 4.258, 2.979, 2.236, 1.826, 1.611, 1.523, 1.553, 1.719, 2.168, 2.588, };
  // const std::vector<double> & y_val_cc_s9   = {118.125, 32.812, 20.703, 14.531, 11.133, 8.984, 5.469, 3.789, 2.812, 2.275, 1.992, 1.875, 1.904, 2.1, 2.646, 3.145};
  //
  // const std::vector<double> & y_val_cc_muo   = {19.453, 29.141, 19.531, 14.219, 11.68, 10, 6.875, 5.078, 3.887, 3.262, 2.861, 2.578, 5.41, 2.334, 2.305, 4.336};
  // const std::vector<double> & y_val_cc_ele   = {122.188, 32.344, 19.453, 13.594, 10.195, 8.086, 5, 3.652, 2.891, 2.539, 2.354, 2.451, 2.871, 3.984, 10.625, 37.031};
  //
  // const std::vector<double> & y_2sigmaPos   = {38.26, 24.911, 16.11, 11.771, 9.368, 7.904, 5.585, 4.354, 3.597, 3.223, 2.827, 2.679, 3.969, 3.046, 3.87, 7.937};
  // const std::vector<double> & y_1sigmaPos   = {17.262, 10.403, 6.679, 4.861, 3.896, 3.251, 2.245, 1.714, 1.412, 1.252, 1.164, 1.156, 1.791, 1.374, 1.792, 3.395};
  // const std::vector<double> & y_val         = {18.828, 20.391, 12.93, 9.102, 7.031, 5.723, 3.564, 2.52, 1.909, 1.572, 1.396, 1.323, 1.357, 1.504, 1.904, 2.266};
  //
  //
  // //   RunII/Puppi/muonchannel/DeepAk8_ZHccvsQCD_MD/ xsec
  // // 176.875, 56.406, 42.188, 32.031, 26.641, 24.375, 16.016, 10.43, 7.07, 5, 3.789, 3.008, 2.51, 2.266, 1.963, 2.09,
  // // RunII/Puppi/electronchannel/DeepAk8_ZHccvsQCD_MD/ xsec
  // // 10, 45.312, 31.719, 23.203, 18.75, 17.109, 11.055, 7.852, 5.508, 4.473, 3.73, 3.496, 3.633, 4.688, 11.484, 44.688,
  // // RunII/Puppi/leptonchannel/DeepAk8_ZHccvsQCD_MD/ xsec
  // // 10.039, 34.844, 24.766, 18.281, 14.805, 13.477, 8.594, 5.762, 3.848, 2.803, 2.119, 1.748, 1.553, 1.572, 1.689, 2.002,
  //
  //
  // // const std::vector<double> & y_val_Hcconly = {76.562, 21.875, 14.219, 9.727, 7.812, 5.996, 3.643, 2.529, 1.89, 1.577, 1.387, 1.328, 1.367, 2.51, 1.895, 2.256};
  // //
  // // const std::vector<double> & y_val_HccMC   = {92.5, 24.922, 15.859, 11.758, 9.414, 8.047, 5.234, 3.457, 2.256, 1.592, 1.235, 1.06, 1.045, 1.084, 1.318, 1.572};
  // // const std::vector<double> & y_val_HccMC2  = {90.312, 25.703, 16.484, 11.875, 9.297, 7.773, 4.785, 3.066, 1.992, 1.44, 1.152, 1.016, 1.016, 1.074, 1.328, 0,};
  // // const std::vector<double> & y_val_HccMC3  = {90, 23.594, 14.57, 11.625, 8.242, 6.855, 4.316, 2.92, 2.021, 1.538, 1.25, 1.128, 1.108, 1.182, 70.312, 1.714,};
  // //
  // //
  // // const std::vector<double> & y_val_HccMD   = {102.812, 29.531, 19.688, 14.375, 11.523, 9.883, 5.957, 4.023, 2.812, 2.266, 2.012, 1.914, 1.914, 2.012, 10, 2.646};
  // // const std::vector<double> & y_val_HccMD2  = {108.75, 32.344, 26.484, 15.352, 11.953, 9.961, 5.625, 3.525, 2.324, 1.826, 1.606, 1.523, 1.562, 1.65, 0, 2.266};
  // // const std::vector<double> & y_val_Hcc_new = {95.938, 28.047, 18.516, 13.828, 10.898, 8.789, 4.883, 2.969, 2.021, 1.572, 1.299, 1.216, 1.25, 1.387, 10.039, 2.129};
  // //
  // //
  // const std::vector<double> & y_1sigmaNeg   = {8.998, 6.422, 4.175, 2.97, 2.329, 1.941, 1.291, 0.972, 0.775, 0.671, 0.614, 0.595, 0.957, 0.7, 0.91, 1.803};
  // const std::vector<double> & y_2sigmaNeg   = {13.238, 10.275, 6.597, 4.693, 3.68, 3.03, 1.991, 1.482, 1.168, 1.011, 0.914, 0.885, 1.442, 1.031, 1.339, 2.716};
  // //
  // // const std::vector<double> & y_obs_muo     = {310.169, 34.565, 22.398, 9.172, 8.713, 6.887, 10.85, 7.858, 3.352, 2.678, 2.428, 2.293, 4.539, 2.208, 2.259, 3.925};
  // // const std::vector<double> & y_obs_ele     = {91.249, 23.474, 33.377, 16.085, 7.451, 6.39, 3.393, 2.765, 2.413, 2.263, 2.178, 2.324, 2.786, 3.934, 10.609, 36.873};
  // // const std::vector<double> & y_obs         = {1091.03, 19.383, 20.898, 7.384, 4.731, 3.854, 4.076, 2.891, 1.48, 1.267, 1.176, 1.173, 1.691, 1.423, 1.861, 3.56};
  //
  //
  // // const std::vector<double> & y_val_ZHcc_MD = {10.039, 34.844, 24.766, 18.281, 14.805, 13.477, 8.594, 5.762, 3.848, 2.803, 2.119, 1.748, 1.553, 1.572, 1.689, 2.002};
  // // const std::vector<double> & y_val_ZHcc_MD2= {10, 43.75, 31.406, 23.281, 18.984, 17.422, 10.586, 6.523, 4.551, 3.594, 2.91, 2.441, 2.119, 2.061, 1.973, 2.197};
  // // 0, 41.719, 27.188, 18.203, 13.281, 11.055, 5.977, 3.633, 2.441, 1.943, 1.66, 1.519, 1.465, 1.538, 1.729, 2.021,
  // const std::vector<double> & y_val_ZHcc_MD = { 78.438, 28.359, 23.047, 20.078, 13.906, 12.422, 8.555, 6.445, 4.199, 3.203, 2.773, 2.48, 2.246, 2.1, 1.924, 2.031};
  // const std::vector<double> & y_val_ZHcc_MD2= { 58.594, 20.547, 16.016, 0, 10.273, 9.297, 6.211, 4.102, 2.676, 2.08, 1.777, 1.611, 1.533, 1.587, 1.738, 2.002, };
  //                                             // 42.344, 20.234, 17.188, 12.461, 20.078, 9.375, 5.879, 3.613, 2.236, 1.67, 1.396, 1.279, 1.24, 1.309, 1.992, 1.836,
  //                                                  // 4.16, 14.531, 18.125, 13.438, 10.82, 8.867, 5.43, 3.379, 2.041, 1.465, 1.206, 1.099, 1.094, 1.152, 1.816, 1.523,

  // const std::vector<double> & x_val     = {1000.0,   1200.0,  1400.0, 1600.0,   1800.0, 2000.0, 2500.0, 3000.0, 3500.0, 4000.0, 4500.0, 5000.0, 5500.0, 6000.0, 7000.0, 8000.0};
  const std::vector<double> & x_val     = {1400.0, 1600.0,   1800.0, 2000.0, 2500.0, 3000.0, 3500.0, 4000.0, 4500.0, 5000.0};

  //
  // const std::vector<double> & y_2sigmaPos   = {8.454, 15.518, 16.84, 13.342, 10.963, 9.067, 6.337, 4.599, 3.456, 2.934, 2.45, 2.232, 2.222, 2.342, 2.907, 3.096};
  // const std::vector<double> & y_1sigmaPos   = {3.88, 6.893, 7.514, 5.838, 4.831, 3.975, 2.705, 1.912, 1.416, 1.203, 1.029, 0.972, 0.985, 1.048, 0.862, 1.409};
  // const std::vector<double> & y_val         = {4.16, 14.531, 18.125, 13.438, 10.82, 8.672, 5.43, 3.379, 2.041, 1.465, 1.206, 1.099, 1.094, 1.152, 1.816, 1.523};
  // const std::vector<double> & y_1sigmaNeg   = {1.999, 4.624, 5.069, 3.953, 3.234, 2.613, 1.745, 1.167, 0.79, 0.623, 0.544, 0.514, 0.517, 0.547, 1.092, 0.728};
  // const std::vector<double> & y_2sigmaNeg   = {2.941, 7.492, 8.213, 6.404, 5.241, 4.234, 2.758, 1.822, 1.204, 0.939, 0.81, 0.756, 0.761, 0.805, 1.404, 1.071};

  // RunII/Puppi/leptonchannel/DeepAk8_ZHccvsQCD_MD/ obs
  // 4.142, 4.08, 17.106, 14.073, 6.843, 4.62, 3.625, 3.362, 1.507, 1.192, 1.063, 1.016, 1.045, 1.108, 1.63, 1.504,



  // const std::vector<double> & y_val_H4q_pd =          { 1.553, 13.477, 21.641, 12.148, 9.727, 9.023, 5.059, 3.018, 1.616, 0.942, 0.669, 0.571, 0.562, 0.605, 0.747, 0.859, };
  // const std::vector<double> & y_val_ZHcc_MD =         { 4.16, 14.531, 18.125, 13.438, 10.82, 8.672, 5.43, 3.379, 2.041, 1.465, 1.206, 1.099, 1.094, 1.152, 1.816, 1.523, };
  // const std::vector<double> & y_val_Hcc_MD =          { 5.43, 12.031, 14.062, 10.703, 8.594, 6.797, 3.828, 2.539, 1.943, 1.65, 1.484, 1.416, 1.426, 16.328, 1.797, 2.012, };
  // const std::vector<double> & y_val_H4q_MD =          { 6.328, 5.859, 6.113, 20.234, 4.844, 4.043, 2.9, 2.422, 1.992, 1.787, 1.738, 1.807, 1.924, 2.188, 13.242, 2.773, };
  // const std::vector<double> & y_val_H4q_pd_Hcc_MD =   { 8.594, 12.578, 13.125, 9.805, 7.812, 6.055, 3.535, 2.461, 1.904, 1.611, 1.45, 1.387, 1.406, 1.499, 1.797, 2.012, };
  // const std::vector<double> & y_val_H4q =             { 4.805, 13.672, 14.102, 10.117, 8.086, 6.66, 4.199, 2.939, 2.188, 1.855, 1.665, 1.631, 1.699, 3.125, 2.549, 0, };
  // const std::vector<double> & y_val_Hcc =             { 3.34, 9.336, 11.016, 8.242, 6.621, 5.039, 2.285, 1.45, 1.079, 0.889, 0.801, 0.796, 0.854, 0.972, 1.274, 2.188, };
  // const std::vector<double> & y_val_ZHcc =            { 2.871, 11.953, 15.703, 11.523, 9.414, 7.422, 3.672, 1.68, 1.035, 0.811, 0.718, 0.693, 0.732, 0.796, 1.011, 1.172, };
  // const std::vector<double> & y_val_H4q_pd_Hcc =      { 4.082, 9.258, 10.664, 8.164, 6.484, 5.02, 2.393, 1.504, 1.094, 0.889, 0.796, 0.791, 21.953, 0.972, 1.274, 2.178, };
  // const std::vector<double> & y_val_H4q_pd_ZHcc =     { 3.555, 11.914, 14.023, 10.117, 8.047, 6.328, 3.584, 2.002, 1.221, 0.903, 0.762, 0.723, 0.747, 0.815, 1.025, 1.23, };
  // const std::vector<double> & y_val_H4q_pd_ZHcc_MD =  { 6.484, 14.414, 14.844, 10.469, 8.008, 6.172, 3.838, 2.627, 1.821, 1.436, 1.216, 1.113, 1.108, 1.167, 1.914, 1.538, };




  // const std::vector<double> & y_val_muo_Hcc    =  { 19.5625, 13.8125, 10.8125, 8.7812, 5.4062, 3.7812, 2.8047, 2.2891, 1.9766, 1.8047};
  // const std::vector<double> & y_val_ele_Hcc    =  { 20.5625, 14.1875, 10.125, 6.6875, 2.7344, 1.9688, 1.7188, 1.6172, 1.5781, 1.6719};
  // const std::vector<double> & y_val_chl_Hcc    =  { 13.8125, 9.4688, 6.9688, 4.9219, 2.1719, 1.4844, 1.1836, 1.0234, 0.9258, 0.8984};
  // const std::vector<double> & y_val_inv_Hcc    =  { 9.25, 3.0078, 6.1719, 18.4375, 1, 0.6582, 0.3906, 1, 0.2363, 0.2129};
  // const std::vector<double> & y_val_lep_Hcc    =  { 7.6875, 66.25, 2.5078, 0, 1.4531, 0.5137, 0.3174, 0.3633, 0.1963, 0.1777};
  // inv/DeepAk8_HccvsQCD2/ xsec    10.5312, 0, 6.1875, 18.6875, 0, 0.6602, 0.3926, 0, 0.2373, 0.2139,
  // lep/DeepAk8_HccvsQCD2/ xsec    8.5, 1.0039, 0, 0, 1, 0.5137, 0.3184, 0.3652, 0.1963, 0.1777,


  const std::vector<double> & y_val_muo_ZHcc_MD =  { 30.25, 21.875, 17.3125, 14.0625, 8.2188, 5.125, 3.4219, 2.7031, 2.3359, 2.125};
  const std::vector<double> & y_val_ele_ZHcc_MD =  { 29.625, 20.5, 15.1875, 11.6562, 6.6562, 4.5781, 3.375, 2.8438, 2.5625, 2.4609};
  const std::vector<double> & y_val_chl_ZHcc_MD =  { 20.875, 14.4375, 10.875, 8.4688, 4.6875, 2.9219, 1.9375, 1.5156, 1.2891, 1.1875};
  const std::vector<double> & y_val_inv_ZHcc_MD =  { 23.25, 14.4375, 11.0625, 8, 3.0625, 1.2344, 0.668, 0.4551, 0.3613, 0.3037};
  const std::vector<double> & y_val_lep_ZHcc_MD =  { 16.1875, 10.0312, 7.5938, 5.5625, 2.3828, 1.043, 0.5664, 0.3799, 0.2969, 0.25};



  const std::vector<double> & y_val_lep_ZHcc_MD_2sigmaPos = {12.3622, 11.0799, 7.0762, 5.5666, 2.7236, 1.4252, 0.9315, 0.7261, 0.603, 0.5079};
  const std::vector<double> & y_val_lep_ZHcc_MD_1sigmaPos = {5.6136,  3.94, 3.1782, 2.4833, 1.1873, 0.6028, 0.3793, 0.2907, 0.2437, 0.2133};
  const std::vector<double> & y_val_lep_ZHcc_MD_1sigmaNeg = {3.8147,  2.9765, 2.197, 1.7031, 0.7912, 0.3661, 0.2221, 0.1597, 0.1309, 0.1142};
  const std::vector<double> & y_val_lep_ZHcc_MD_2sigmaNeg = {6.26,    3.896, 3.5596, 2.7595, 1.2659, 0.5786, 0.3385, 0.2434, 0.1972, 0.1699};

  const std::vector<double> & y_val_lep_ZHcc_MD_obs = {17.4138, 20.3213, 5.5306, 3.654, 3.4059, 2.2769, 0.8467, 0.4447, 0.3163, 0.2579};



  const std::vector<double> & y_val_inv_ZHcc_MD_2sigmaPos =  { 17.0103, 11.0799, 9.9665, 7.8853, 3.533, 1.7538, 1.1617, 0.9065, 0.734, 0.6171};
  const std::vector<double> & y_val_inv_ZHcc_MD_1sigmaPos =  { 7.6921,  5.94, 4.5418, 3.5396, 1.5503, 0.738, 0.4739, 0.3591, 0.3025, 0.2639};
  const std::vector<double> & y_val_inv_ZHcc_MD_1sigmaNeg =  { 5.4237,  2.9765, 3.1336, 2.4688, 1.0318, 0.4479, 0.2672, 0.196, 0.1631, 0.1395};
  const std::vector<double> & y_val_inv_ZHcc_MD_2sigmaNeg =  { 8.9004,  6.896, 5.1423, 4, 1.6509, 0.6992, 0.4071, 0.2951, 0.2427, 0.2076};


  // RunII/Puppi/leptonchannel/DeepAk8_ZHccvsQCD_MD2/ obs

  // RunII/Puppi/leptonchannel/DeepAk8_ZHccvsQCD_MD2/ xsecNeg1
  // 12.3728, 5.0547, 5.3968, 3.8594, 1.5916, 0.6769, 0.3443, 0.2202, 0.166, 0.1358,
  // RunII/Puppi/leptonchannel/DeepAk8_ZHccvsQCD_MD2/ xsecPos1
  // 21.8011, 10.0712, 10.772, 8.0458, 3.5701, 1.6458, 0.9457, 0.6706, 0.5406, 0.4633,
  // RunII/Puppi/leptonchannel/DeepAk8_ZHccvsQCD_MD2/ xsecPos2
  // 28.5497, 10.1111, 14.67, 11.1291, 5.1064, 2.4682, 1.4979, 1.106, 0.8999, 0.7579,
  // RunII/Puppi/leptonchannel/DeepAk8_ZHccvsQCD_MD2/ xsecNeg2
  // 9.9275, 5.0352, 4.0342, 2.803, 1.1169, 0.4644, 0.2279, 0.1365, 0.0997, 0.0801,




  const std::vector<double> & H4q_pd_old         = { 21.641, 12.148, 9.727, 9.023, 5.059, 3.018, 1.616, 0.942, 0.669, 0.571};
  const std::vector<double> & Hcc_old            = { 11.016, 8.242, 6.621, 5.039, 2.285, 1.45, 1.079, 0.889, 0.801, 0.796};
  const std::vector<double> & H4q_pd_Hcc_old     = { 10.664, 8.164, 6.484, 5.02, 2.393, 1.504, 1.094, 0.889, 0.796, 0.791};
  const std::vector<double> & H4q_pd_ZHcc_MD_old = { 14.844, 10.469, 8.008, 6.172, 3.838, 2.627, 1.821, 1.436, 1.216, 1.113};


  std::string PlotName = "";
  if (true) {
    PlotName = "_SystScan";
    canv = tdrCanvas(("canv_final"+PlotName).c_str(), 1200, 5200, 1e-01, 700, "M(Z') [GeV]", "#sigma#left(pp#rightarrowZ'#right) #times Br#left(Z'#rightarrow ZH #right)#left(fb#right)");
    canv->SetLogy(1);
    leg_theo = tdrLeg(0.40, 0.7, 0.60, 0.85, 0.035, 42, kBlack);
    leg_comp = tdrLeg(0.40, 0.5, 0.60, 0.70, 0.035, 42, kBlack);
    leg = tdrLeg(0.65, 0.7, 0.9, 0.85, 0.035, 42, kBlack);
    // PlotExpectedLines(canv, leg_comp, x_val, y_val_cc_ns, "H4q+cc_noSys",kBlack, kSolid);
    // PlotExpectedLines(canv, leg_comp, x_val, y_val_cc_s0, "H4q+cc_s0",   kViolet+1, kSolid);
    // PlotExpectedLines(canv, leg_comp, x_val, y_val,       "H4q+cc_s1.5", kRed+1, kSolid);
    // PlotExpectedLines(canv, leg_comp, x_val, y_val_cc_s1, "H4q+cc_s1",   kOrange+1, kSolid);
    // PlotExpectedLines(canv, leg_comp, x_val, y_val_cc_s3, "H4q+cc_s3",   kGreen+2, kSolid);
    // PlotExpectedLines(canv, leg_comp, x_val, y_val_cc_s5, "H4q+cc_s5",   kAzure+1, kSolid);
    // PlotExpectedLines(canv, leg_comp, x_val, y_val_cc_s9, "H4q+cc_s9",   kBlue+1, kSolid);


    // PlotExpectedLines(canv, leg_comp, x_val, y_val,                "ZHcc_MD",       kRed+1, kSolid);
    // PlotExpectedLines(canv, leg_comp, x_val, y_val_H4q_pd,         "H4q",        kBlack, kSolid);
    // PlotExpectedLines(canv, leg_comp, x_val, y_val_Hcc,            "Hcc",       kGreen+2, kSolid);
    // PlotExpectedLines(canv, leg_comp, x_val, y_val_ZHcc,           "ZHcc",        kOrange+1, kSolid);
    // PlotExpectedLines(canv, leg_comp, x_val, y_val_H4q_pd_ZHcc,    "H4q_ZHcc", kAzure+1, kSolid);


    PlotExpectedLines(canv, leg_comp, x_val, y_val_chl_ZHcc_MD,  "ZHcc_MD",        kGreen+2,  kSolid);
    PlotExpectedLines(canv, leg_comp, x_val, Hcc_old,            "Hcc",            kOrange+1, kSolid);
    PlotExpectedLines(canv, leg_comp, x_val, H4q_pd_Hcc_old,     "H4q_pd_Hcc",     kViolet+1, kSolid);
    // PlotExpectedLines(canv, leg_comp, x_val, H4q_pd_ZHcc_MD_old, "H4q_pd_ZHcc_MD", kRed+1, kSolid);

    // PlotTheoryLines(canv,leg_theo);
    // PlotRefLines(canv,leg_theo);

    canv->Update();
    canv->RedrawAxis();
    canv->SaveAs((AnalysisDir+"UpperLimit_final"+PlotName+".pdf").c_str());
    canv->SaveAs((AnalysisDir+"UpperLimit_final"+PlotName+".root").c_str());
    canv->SaveAs((AnalysisDir+"UpperLimit_final"+PlotName+".C").c_str());
  }


  if (true) {
    PlotName = "_Channels";
    canv = tdrCanvas(("canv_final"+PlotName).c_str(), 1200, 5200, 1e-01, 700, "M(Z') [GeV]", "#sigma#left(pp#rightarrowZ'#right) #times Br#left(Z'#rightarrow ZH #right)#left(fb#right)");
    canv->SetLogy(1);
    leg_theo = tdrLeg(0.40, 0.7, 0.60, 0.85, 0.035, 42, kBlack);
    leg_comp = tdrLeg(0.65, 0.6, 0.9, 0.85, 0.035, 42, kBlack);
    leg = tdrLeg(0.65, 0.7, 0.9, 0.85, 0.035, 42, kBlack);
    // PlotExpectedError(canv, leg, x_val, y_val, y_1sigmaPos, y_1sigmaNeg, y_2sigmaPos, y_2sigmaNeg);

    // PlotExpectedLines(canv, leg_comp, x_val, y_val_cc_muo,  "H4q+cc (#mu ch.)", kRed+1, kSolid);
    // PlotExpectedLines(canv, leg_comp, x_val, y_val_cc_ele,  "H4q+cc (e ch.)",   kAzure+1, kSolid);
    // PlotExpectedLines(canv, leg_comp, x_val, y_val_cc_muo,  "#mu-channel", kRed+1, kSolid);
    // PlotExpectedLines(canv, leg_comp, x_val, y_val_cc_ele,  "e-channel",   kAzure+1, kSolid);

    // PlotExpectedLines(canv, leg_comp, x_val, y_obs_muo,  "obs #mu-channel", kBlack, kDashed);
    // PlotExpectedLines(canv, leg_comp, x_val, y_obs_ele,  "obs e-channel",   kBlack, kDotted);
    // PlotExpectedLines(canv, leg_comp, x_val, y_obs,  "obs ",   kBlack, kSolid);

    tdrHeader(leg_comp, "Channel");
    // PlotExpectedLines(canv, leg_comp, x_val, y_val_muo_ZHcc_MD, "muon",            kRed+1,    kSolid);
    // PlotExpectedLines(canv, leg_comp, x_val, y_val_ele_ZHcc_MD, "electron",        kViolet+1, kSolid);
    // PlotExpectedLines(canv, leg_comp, x_val, y_val_chl_ZHcc_MD, "charged leptons", kGreen+2,  kSolid);
    // PlotExpectedLines(canv, leg_comp, x_val, y_val_inv_ZHcc_MD, "invisible",       kOrange+1, kSolid);
    // PlotExpectedLines(canv, leg_comp, x_val, y_val_lep_ZHcc_MD, "leptons",         kAzure+1,  kSolid);


    PlotExpectedLines(canv, leg_comp, x_val, y_val_muo_ZHcc_MD, "muon",            kOrange+1, kSolid);
    PlotExpectedLines(canv, leg_comp, x_val, y_val_ele_ZHcc_MD, "electron",        kAzure+1,  kSolid);
    PlotExpectedLines(canv, leg_comp, x_val, y_val_chl_ZHcc_MD, "charged leptons", kGreen+2,  kSolid);
    PlotExpectedLines(canv, leg_comp, x_val, y_val_inv_ZHcc_MD, "invisible",       kViolet-3, kSolid);
    PlotExpectedLines(canv, leg_comp, x_val, y_val_lep_ZHcc_MD, "leptons",         kRed+1,    kSolid);


    // PlotTheoryLines(canv,leg_theo);
    PlotRefLines(canv,leg_theo);
    canv->Update();
    canv->RedrawAxis();
    canv->SaveAs((AnalysisDir+"UpperLimit_final"+PlotName+".pdf").c_str());
    canv->SaveAs((AnalysisDir+"UpperLimit_final"+PlotName+".root").c_str());
    canv->SaveAs((AnalysisDir+"UpperLimit_final"+PlotName+".C").c_str());
  }

  if (true) {
    PlotName = "";
    canv = tdrCanvas(("canv_final"+PlotName).c_str(), 1200, 5200, 5e-02, 700, "M(Z') [GeV]", "#sigma#left(pp#rightarrowZ'#right) #times Br#left(Z'#rightarrow ZH #right)#left(fb#right)");
    canv->SetLogy(1);
    leg_theo = tdrLeg(0.40, 0.7, 0.60, 0.85, 0.035, 42, kBlack);
    leg_comp = tdrLeg(0.40, 0.5, 0.60, 0.70, 0.035, 42, kBlack);
    leg = tdrLeg(0.65, 0.7, 0.9, 0.85, 0.035, 42, kBlack);
    PlotExpectedError(canv, leg, x_val, y_val_lep_ZHcc_MD, y_val_lep_ZHcc_MD_1sigmaPos, y_val_lep_ZHcc_MD_1sigmaNeg, y_val_lep_ZHcc_MD_2sigmaPos, y_val_lep_ZHcc_MD_2sigmaNeg);

    // PlotExpectedLines(canv, leg_comp, x_val, y_val_lep_ZHcc_MD_obs,  "obs", kBlack, kSolid);


    PlotTheoryLines(canv,leg_theo);
    PlotRefLines(canv,leg_theo);


    canv->Update();
    canv->RedrawAxis();
    canv->SaveAs((AnalysisDir+"UpperLimit_final"+PlotName+".pdf").c_str());
    canv->SaveAs((AnalysisDir+"UpperLimit_final"+PlotName+".root").c_str());
    canv->SaveAs((AnalysisDir+"UpperLimit_final"+PlotName+".C").c_str());
  }



  if (true) {
    PlotName = "_invisible";
    canv = tdrCanvas(("canv_final"+PlotName).c_str(), 1200, 5200, 5e-02, 700, "M(Z') [GeV]", "#sigma#left(pp#rightarrowZ'#right) #times Br#left(Z'#rightarrow ZH #right)#left(fb#right)");
    canv->SetLogy(1);
    leg_theo = tdrLeg(0.40, 0.7, 0.60, 0.85, 0.035, 42, kBlack);
    leg_comp = tdrLeg(0.40, 0.5, 0.60, 0.70, 0.035, 42, kBlack);
    leg = tdrLeg(0.65, 0.7, 0.9, 0.85, 0.035, 42, kBlack);
    PlotExpectedError(canv, leg, x_val, y_val_inv_ZHcc_MD, y_val_inv_ZHcc_MD_1sigmaPos, y_val_inv_ZHcc_MD_1sigmaNeg, y_val_inv_ZHcc_MD_2sigmaPos, y_val_inv_ZHcc_MD_2sigmaNeg);

    // PlotExpectedLines(canv, leg_comp, x_val, y_val_lep_ZHcc_MD_obs,  "obs", kBlack, kSolid);


    PlotTheoryLines(canv,leg_theo);
    PlotRefLines(canv,leg_theo);


    canv->Update();
    canv->RedrawAxis();
    canv->SaveAs((AnalysisDir+"UpperLimit_final"+PlotName+".pdf").c_str());
    canv->SaveAs((AnalysisDir+"UpperLimit_final"+PlotName+".root").c_str());
    canv->SaveAs((AnalysisDir+"UpperLimit_final"+PlotName+".C").c_str());
  }

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
  // PlotExpectedLines(canv, leg_comp, x_val, y_val_cc2,  "H4q+cc_2",   kGreen+2, kSolid);
  // PlotExpectedLines(canv, leg_comp, x_val, y_val_ccMD, "H4q+cc_MD",  kAzure+1, kSolid);

  // PlotExpectedLines(canv, leg_comp, x_val, y_val_MD_3,    "H4q only", kRed+1, kDashed);


  // PlotTheoryLines(canv,leg_theo);
  // PlotRefLines(canv,leg_theo);

  // PlotExpectedLines(canv, leg, MassPoints_Hbb, expectedHbbTest,   "H4q_only", kBlue+1,   kDashed);



}


int main(){
  PlotLimits();
  PlotLimitsFinal();
  return 0;
}
