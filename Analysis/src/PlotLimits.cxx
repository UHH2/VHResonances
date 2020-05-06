#include "PlotLimits.hpp"



void PlotLimits(bool doObs, bool isHbb) {

  std::map<std::string, int> colors = {
    {"2016",kViolet+1}, {"2017",kOrange+1}, {"2018",kCyan+1}, {"RunII",kGreen+2}, {"fullRunII", kGreen+1},
    {"btag_DeepBoosted_H4qvsQCD",kGreen+1}, {"NN",kOrange+1}, {"NN_1",kViolet+1}, {"NN_2",kCyan+1}, {"CNN",kBlue+1}, {"tau42",kRed+1},
    {"btag_DeepBoosted_HbbvsQCD",kGreen+1}, {"btag_DeepBoosted_probHbb",kViolet+1}, {"tau32",kBlue+1}, {"tau21",kOrange+1},
    {"btag_DeepBoosted_H4qvsQCDp2", kOrange+1}, {"btag_DeepBoosted_H4qvsQCDp02", kViolet+1}, {"btag_DeepBoosted_H4qvsQCDpt1000", kCyan+1},
    {"btag_DeepBoosted_H4qvsQCDpt1000p2", kBlue+1}, {"btag_DeepBoosted_H4qvsQCDpt1000p02", kRed+1},
  };


  TGraph* gr_theo = new TGraphErrors(nPoints, &(MassPoints[0]), &theo_xsec[0]);
  TGraph* gr_Hbb = new TGraphErrors(MassPoints_Hbb.size(), &(MassPoints_Hbb[0]), &expectedHbb[0]);
  TGraph* gr_Hbb0b = new TGraphErrors(MassPoints_Hbb.size(), &(MassPoints_Hbb[0]), &expectedHbb0b[0]);
  gr_theo->SetLineWidth(2);
  gr_Hbb->SetLineWidth(2);
  gr_Hbb0b->SetLineWidth(2);


  // std::string mode = "CB";
  std::string mode = "Exp_2";
  std::string Path_ANALYSIS = std::getenv("CMSSW_BASE"); Path_ANALYSIS += "/src/UHH2/VHResonances/Analysis/";
  std::string studies = "nominal";
  // std::string studies = "noSignalFlatUncertainty";

  std::string Method = "AsymptoticLimits";
  std::string extraOptionsText = "Asimov";
  std::string AnalysisDir = Path_ANALYSIS+"Limits/"+studies+"/";
  if (isHbb)  AnalysisDir = Path_ANALYSIS+"Limits/"+studies+"/Hbb/";
  // std::vector<std::string> collections = {"Puppi", "CHS"};
  std::vector<std::string> collections = {"Puppi"};
  // std::vector<std::string> channels = {"muonchannel", "electronchannel"};
  std::vector<std::string> channels = {"muonchannel", "electronchannel", "leptonchannel"};
  std::vector<std::string> years = {"2016", "2017", "2018", "RunII", "fullRunII"};
  // std::vector<std::string> years = {"2016"};
  // std::vector<std::string> channels = {"muonchannel"};
  // std::vector<std::string> histFolders = {"btag_DeepBoosted_H4qvsQCD", "NN", "tau21","tau31", "tau41", "tau32", "tau42", "tau43"};
  // std::vector<std::string> histFolders = {"btag_DeepBoosted_H4qvsQCD"};
  std::vector<std::string> histFolders = {"btag_DeepBoosted_H4qvsQCD", "btag_DeepBoosted_H4qvsQCDp2", "btag_DeepBoosted_H4qvsQCDp02", "btag_DeepBoosted_H4qvsQCDpt1000", "btag_DeepBoosted_H4qvsQCDpt1000p2", "btag_DeepBoosted_H4qvsQCDpt1000p02"};
  if (isHbb) histFolders = {"btag_DeepBoosted_HbbvsQCD", "btag_DeepBoosted_probHbb", "tau42", "tau32", "tau21" };

  std::unordered_map<std::string, std::unordered_map<std::string, std::vector<double>>> Limits;

  for (std::string year: years) {
    for (std::string collection: collections) {
      for (std::string channel: channels) {
        for (std::string histFolder : histFolders) {
          std::string workingDir = year+"/"+collection+"/"+channel+"/"+histFolder+"/";

          double pbTofb = 1000.;
          // double xsec_ref_ = 0.1; // 0.1 comes from normalising the signal strenght
          double xsec_ref_ = (xsec_ref.find(histFolder) != xsec_ref.end())? xsec_ref.at(histFolder): xsec_ref.at("default_value"); // default_value = 1
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
          for (unsigned int i=0; i<MassPoints.size(); ++i) {
            std::string filename = AnalysisDir+workingDir+"/datacards/out/DataCard_"+year+"_"+collection+"_"+channel+"_"+histFolder+"_M"+std::to_string((int)MassPoints[i])+"_"+mode+"_"+Method+"_"+extraOptionsText+".out";
            std::ifstream txtfile(filename.c_str(), std::ios::in);

            std::string line;
            bool found= false;
            int count = 0;

            while (count<200 && !found && !txtfile.eof() && !gSystem->AccessPathName(filename.c_str())) {
              getline(txtfile, line);
              std::size_t pos = line.find("-- AsymptoticLimits ( CLs ) --");
              if (pos!=std::string::npos) found=true;
              count++;
            }

            if (count<100) { // TODO
              getline(txtfile, line); Limits[workingDir]["obs"][i] = atof(line.substr(line.find("<")+1).c_str())*norm;
              getline(txtfile, line); Limits[workingDir]["xsecNeg2"][i] = atof(line.substr(line.find("<")+1).c_str())*norm;
              getline(txtfile, line); Limits[workingDir]["xsecNeg1"][i] = atof(line.substr(line.find("<")+1).c_str())*norm;
              getline(txtfile, line); Limits[workingDir]["xsec"][i] = atof(line.substr(line.find("<")+1).c_str())*norm;
              getline(txtfile, line); Limits[workingDir]["xsecPos1"][i] = atof(line.substr(line.find("<")+1).c_str())*norm;
              getline(txtfile, line); Limits[workingDir]["xsecPos2"][i] = atof(line.substr(line.find("<")+1).c_str())*norm;

              Limits[workingDir]["expNeg2"][i]=Limits[workingDir]["xsec"][i]-Limits[workingDir]["xsecNeg2"][i];
              Limits[workingDir]["expNeg1"][i]=Limits[workingDir]["xsec"][i]-Limits[workingDir]["xsecNeg1"][i];
              Limits[workingDir]["expPos1"][i]=Limits[workingDir]["xsecPos1"][i]-Limits[workingDir]["xsec"][i];
              Limits[workingDir]["expPos2"][i]=Limits[workingDir]["xsecPos2"][i]-Limits[workingDir]["xsec"][i];
            }

            txtfile.close();
          }
        }
      }
    }
  }

  TString nameXaxix = "m(Z')";
  // TString nameYaxix = "#sigma\(pp#rightarrowX\) #times Br\(Z'#rightarrowZ\(ll\) H\(WW\)\) \(fb\)";
  TString nameYaxix = "#sigma#left(pp#rightarrowX#right) #times Br#left(Z'#rightarrow ZH #right)#left(fb#right)";

  writeExtraText = true;       // if extra text
  extraText  = "Work in progress" ;//"Preliminary";
  lumi_13TeV  = TString::Format("%.1f fb^{-1}", lumi_map.at("RunII").at("lumi_fb"));

  std::unordered_map<std::string, TCanvas*> Canvas_Limits_Comparison;
  std::unordered_map<std::string, TLegend*> Legend_Limits_Comparison;

  Canvas_Limits_Comparison["all"] = tdrCanvas("c_limits_comparison_all", plot_lo, plot_hi, 1e-01, 1e3, nameXaxix, nameYaxix);
  Legend_Limits_Comparison["all"] = tdrLeg(0.65,0.5,0.9,0.85, 0.025, 42, kBlack);

  for (auto [name,canvas]: Canvas_Limits_Comparison ) {
    canvas->SetLogy(1);
    canvas->SetGridx(1);
    canvas->SetGridy(1);
    Legend_Limits_Comparison[name]->SetFillColor(kWhite);
    Legend_Limits_Comparison[name]->SetFillStyle(1);
    Legend_Limits_Comparison[name]->SetFillStyle(0);
    tdrDraw(gr_theo,  "C", kFullDotLarge, kRed+1, kSolid, kRed+1);
    tdrDraw(gr_Hbb,  "C", kFullDotLarge, kBlue+1, kDotted, kBlue+1);
    tdrDraw(gr_Hbb0b,  "C", kFullDotLarge, kBlue+1, kDashed, kBlue+1);
    TLegend *leg_ExtraPlot = tdrLeg(0.45,0.7,0.60,0.85, 0.03, 42, kBlack);
    canvas->cd();
    leg_ExtraPlot->AddEntry(gr_theo, "Prediction", "l");
    leg_ExtraPlot->AddEntry(gr_Hbb, "Hbb", "l");
    leg_ExtraPlot->AddEntry(gr_Hbb0b, "Hbb cat 0b", "l");
    leg_ExtraPlot->Draw();
  }


  for (std::string year: years) {
    for (std::string collection: collections) {
      for (std::string channel: channels) {

        writeExtraText = true;       // if extra text
        extraText  = "Work in progress" ;//"Preliminary";
        lumi_13TeV  = TString::Format("%.1f fb^{-1}", lumi_map.at((year=="fullRunII")?"RunII":year).at("lumi_fb"));

        TCanvas* c_xsec_modecomparison = tdrCanvas("c_xsec_modecomparison", plot_lo, plot_hi, 1e-01, 1e5, nameXaxix, nameYaxix);
        c_xsec_modecomparison->SetLogy(1);
        c_xsec_modecomparison->SetGridx(1);
        c_xsec_modecomparison->SetGridy(1);

        TLegend *leg_modecomparison = tdrLeg(0.45,0.65,0.9,0.85, 0.03, 42, kBlack);
        leg_modecomparison->SetFillStyle(1); leg_modecomparison->SetFillColor(kWhite);
        leg_modecomparison->SetFillStyle(0);

        std::string dir_modecomparison = AnalysisDir+"/"+year+"/"+collection+"/"+channel+"/";

        for (std::string histFolder : histFolders) {
          std::string workingDir = year+"/"+collection+"/"+channel+"/"+histFolder+"/";

          TCanvas* c_xsec = tdrCanvas(("xsec"+workingDir).c_str(), plot_lo, plot_hi, 1e-01, 1e5, nameXaxix, nameYaxix);
          c_xsec->SetLogy(1);
          c_xsec->SetGridx(1);
          c_xsec->SetGridy(1);
          c_xsec->cd();

          TGraph *g_obs = new TGraph(nPoints, &(MassPoints[0]), &(Limits[workingDir]["obs"][0]));

          std::vector<double> dummy(nPoints,0);
          TGraphAsymmErrors *g_xsec_1sigma = new TGraphAsymmErrors(nPoints, &(MassPoints[0]), &(Limits[workingDir]["xsec"][0]), &(dummy[0]), &(dummy[0]), &(Limits[workingDir]["expNeg1"][0]), &(Limits[workingDir]["expPos1"][0]));
          TGraphAsymmErrors *g_xsec_2sigma = new TGraphAsymmErrors(nPoints, &(MassPoints[0]), &(Limits[workingDir]["xsec"][0]), &(dummy[0]), &(dummy[0]), &(Limits[workingDir]["expNeg2"][0]), &(Limits[workingDir]["expPos2"][0]));
          TGraphAsymmErrors *g_xsec = new TGraphAsymmErrors(nPoints, &(MassPoints[0]), &(Limits[workingDir]["xsec"][0]));
          TGraphAsymmErrors *g_xsec_mode = new TGraphAsymmErrors(nPoints, &(MassPoints[0]), &(Limits[workingDir]["xsec"][0]));;
          g_xsec->SetLineWidth(2);
          g_xsec_mode->SetLineWidth(2);


          tdrDraw(g_xsec_2sigma,  "3",  20, kOrange, kSolid, kOrange, 1000, kOrange);
          tdrDraw(g_xsec_1sigma,  "3",  20, kGreen+1, kSolid, kGreen+1, 1000, kGreen+1);
          tdrDraw(g_xsec,  "L",  20, kBlack, 2, kBlack, 1000, kBlack);
          if (doObs) tdrDraw(g_obs,  "LP", kFullDotLarge, kBlack, kSolid, kBlack);
          tdrDraw(gr_theo,  "C", kFullDotLarge, kRed+1, kSolid, kRed+1);
          tdrDraw(gr_Hbb,  "C", kFullDotLarge, kBlue+1, kDotted, kBlue+1);
          tdrDraw(gr_Hbb0b,  "C", kFullDotLarge, kBlue+1, kDashed, kBlue+1);


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
          leg_ExtraPlot->AddEntry(gr_theo, "Prediction", "l");
          leg_ExtraPlot->AddEntry(gr_Hbb, "Hbb", "l");
          leg_ExtraPlot->AddEntry(gr_Hbb0b, "Hbb cat 0b", "l");
          leg_ExtraPlot->Draw();


          c_xsec->SaveAs((dir_modecomparison+histFolder+"/"+"UpperLimit_"+histFolder+".pdf").c_str());
          // c_xsec->SaveAs(Form("UpperLimit_%s.root", background.c_str()));

          // TFile *file=new TFile(Form("UpperLimits_xsec_%s.root", background.c_str()), "RECREATE");
          // g_obs->Write("g_obs");
          // g_xsec->Write("g_xsec");
          // g_xsec_1sigma->Write("g_xsec_1sigma");
          // g_xsec_2sigma->Write("g_xsec_2sigma");
          // file->Close();
          c_xsec_modecomparison->cd();
          tdrDraw(g_xsec_mode,  "L",  20, colors[histFolder], 2, colors[histFolder], 1000, colors[histFolder]);
          leg_modecomparison->AddEntry(g_xsec_mode, histFolder.c_str(), "L");

          if (histFolder=="btag_DeepBoosted_H4qvsQCD") {
            Canvas_Limits_Comparison["all"]->cd();
            tdrDraw(g_xsec,  "L",  20, colors[year], (channel=="muonchannel")? kSolid: ((channel=="electronchannel")? kDotted: kDashed), colors[year], 1000, colors[year]);
            Legend_Limits_Comparison["all"]->AddEntry(g_xsec, TString(workingDir).ReplaceAll("Puppi","").ReplaceAll("btag_DeepBoosted_H4qvsQCD","").ReplaceAll("channel","").ReplaceAll("/"," "), "L");
          }
        }
        c_xsec_modecomparison->cd();
        leg_modecomparison->Draw("");
        c_xsec_modecomparison->SaveAs((dir_modecomparison+"UpperLimit_modecomparison.pdf").c_str());
      }
    }
  }
  Legend_Limits_Comparison["all"]->Draw("");
  Canvas_Limits_Comparison["all"]->SaveAs((AnalysisDir+"UpperLimit_comparison.pdf").c_str());



  std::unordered_map<std::string, TCanvas*> Canvas_Limits_Ratio;
  std::unordered_map<std::string, TLegend*> Legend_Limits_Ratio;

  std::string collection = "Puppi";
  std::string channel = "muonchannel";
  std::string histFolder = "btag_DeepBoosted_H4qvsQCD";
  // std::string histFolder = "btag_DeepBoosted_H4qvsQCDp2";
  // std::string histFolder = "btag_DeepBoosted_H4qvsQCDp02";
  std::string workingDir = collection+"/"+channel+"/"+histFolder+"/";
  Canvas_Limits_Ratio[workingDir] = tdrCanvas("c_limits_ratio", plot_lo, plot_hi, 1e-01, 5, nameXaxix, "Ratio Limits");
  Legend_Limits_Ratio[workingDir] = tdrLeg(0.65,0.5,0.9,0.85, 0.025, 42, kBlack);

  for (std::string year: years) {
    std::vector<double> ratio(nPoints,0);
    for (size_t i = 0; i < nPoints; i++) ratio[i] = Limits[year+"/"+workingDir]["xsec"][i]/Limits["RunII/"+workingDir]["xsec"][i];

    TGraph *g_ratio = new TGraph(nPoints, &(MassPoints[0]), &(ratio[0]));
    tdrDraw(g_ratio,  "L",  20, colors[year], (channel=="muonchannel")? kSolid: ((channel=="electronchannel")? kDotted: kDashed), colors[year], 1000, colors[year]);
    Legend_Limits_Ratio[workingDir]->AddEntry(g_ratio, year.c_str(), "L");
    double lumiFact = lumi_map.at("RunII").at("lumi_fb")/lumi_map.at((year=="fullRunII")?"RunII":year).at("lumi_fb");
    lumiFact = TMath::Sqrt(lumiFact);
    TLine *line=new TLine(plot_lo, lumiFact, plot_hi, lumiFact); line->SetLineStyle(kDashed); line->SetLineWidth(2); line->SetLineColor(colors[year]); line->Draw("same");
  }
  Legend_Limits_Ratio[workingDir]->Draw("");
  Canvas_Limits_Ratio[workingDir]->SaveAs((AnalysisDir+"UpperLimit_comparison_"+(TString(workingDir).ReplaceAll("/","_"))+".pdf"));



}


int main(){
  PlotLimits();
  return 0;
}
