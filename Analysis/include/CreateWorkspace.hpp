#pragma once

#include <stdio.h>
#include <string.h>
#include <cstring>
#include <cstdio>
#include <fstream>

#include <TH1F.h>
#include <TROOT.h>
#include <TFile.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TProfile.h>
#include <iostream>
#include <TStyle.h>
#include <TPaveText.h>
#include <THStack.h>
#include <TF1.h>
#include <TSystem.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TLine.h>
#include <TString.h>
#include <RooRealVar.h>
#include <RooFit.h>
#include <RooWorkspace.h>
#include <RooBinning.h>
#include <RooAbsBinning.h>
#include <RooDataHist.h>
#include <TGraphAsymmErrors.h>
#include <RooAbsPdf.h>
#include <RooPlot.h>
#include <RooHist.h>
#include <RooGaussian.h>
#include <RooLandau.h>
#include <RooCBShape.h>
#include <RooFitResult.h>
#include <RooHistPdf.h>
#include <TError.h>
#include <RooGlobalFunc.h>
#include <TRandom3.h>

#include "RooNovosibirsk.h"
#include "HiggsAnalysis/CombinedLimit/interface/RevCrystalBall.hpp"
#include "HiggsAnalysis/CombinedLimit/interface/ExpGaussExp.hpp"
#include "HiggsAnalysis/CombinedLimit/interface/PolinomialExponent_1p.hpp"
#include "HiggsAnalysis/CombinedLimit/interface/PolinomialExponent_2p.hpp"
#include "HiggsAnalysis/CombinedLimit/interface/PolinomialExponent_3p.hpp"
#include "HiggsAnalysis/CombinedLimit/interface/PolinomialExponent_4p.hpp"
#include "HiggsAnalysis/CombinedLimit/interface/PolinomialExponent_5p.hpp"
#include "HiggsAnalysis/CombinedLimit/interface/PolinomialExponent_6p.hpp"
#include "HiggsAnalysis/CombinedLimit/interface/BkgPdf4p.hpp"
#include "HiggsAnalysis/CombinedLimit/interface/BkgPdf3p.hpp"

#include "UHH2/core/include/Utils.h"
#include "UHH2/VHResonances/include/constants.hpp"
#include "UHH2/VHResonances/Analysis/include/tdrstyle_all.hpp"


/*
******************************************
/         Background Estimation           /
/-----------------------------------------/
/  fitMC: bkg_pred= MC_SR* Data_SR/MC_SR  /
/ !fitMC: bkg_pred= Data_SR               /
/  doObs: show Data_SR                    /
/  doObs: data_obs= Data_SR               /
/ !doObs: data_obs= bkg_pred*Norm(Data_SR)/
*******************************************
*/

/*
******************************************
/        Explaination of Integrals       /
/-----------------------------------------/
/ CalculateIntegral(histo_map[mode],x_lo,x_hi,doBinWidth) = Integral over the range where the histogram is defined, using TH1F
/ rooHist_map[mode]->sum(!doBinWidth) = Integral over the range where the histogram is defined, using RooDataHist
/ CalculateIntegral(histo_map[mode],fit_lo,fit_hi,doBinWidth) = Integral over the fit range, using TH1F
/ CalculateFractionArea(histo_map[mode],fit_lo,fit_hi, x_lo, x_hi,doBinWidth) = Fraction of fitting area over the full range where the histo is defined
/ CalculateIntegral(histo_map[mode],fit_lo,fit_hi,doBinWidth)/CalculateFractionArea(histo_map[mode],fit_lo,fit_hi, x_lo, x_hi,doBinWidth) = Should give the same value of the integral in the full range
/ Normalization in the full range is rooHist_map[mode]->sum(!doBinWidth). Needed in combine because the histo is defined here.
/ i_tot = Integral over the full range, using RooAbsPdf. Should be 1 because it's normalized
/ i_fit = Integral over the fit range, using RooAbsPdf
/ CalculateIntegral(histo_map[mode],fit_lo,fit_hi,doBinWidth)*i_tot/i_fit = Normalization of the RooAbsPdf, such that it has the correct norm in the fitting range.
/ nEventsSR*i_tot/i_fit = Normalization of the RooAbsPdf, such that it has the correct norm in the fitting range as the event in the SR.
*******************************************
*/


// const std::vector<double> MyMassPoints = {1000, 1200, 1400, 1600, 1800, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500, 6000, 7000, 8000};
const std::vector<double> MyMassPoints = {1400, 1600, 1800, 2000, 2500, 3000, 3500, 4000, 4500, 5000};
// const std::vector<double> MyMassPoints = {1000, 1200, 1400, 1600, 1800, 2000, 2500, 3000, 3500, 4000, 4500, 5000};

double GetRange(TH1F* h, double x);

double CalculateIntegral(TH1F* h, double min, double max, bool dorebin);
double CalculateFractionArea(TH1F* h, double min, double max, double x_min, double x_max, bool dorebin);
double CalculateFractionAreaPDF(RooAbsPdf* PDF, RooRealVar x_var, double fit_lo, double fit_hi);

double DoFTest(double chi2_1, double chi2_2, double npar_1, double npar_2, double n);
void CalculateChiSquare(double& chi2, int& nbins, RooHist* hpull, double xmin, double xmax);

std::string GetSgName(int mass, std::string syst="nominal");

std::unordered_map<std::string, int> Colors = {
  { "Landau",      kRed+1},
  { "Flatte",      kRed+1},
  { "LN",          kRed+1},
  { "LE",          kBlue+1},
  { "GE",          kGreen+1},
  { "BkgPdf3p",    kAzure+7},
  { "BkgPdf4p",    kMagenta+1},

  { "NO",          kAzure+7},
  { "CB",          kOrange+1},
  { "Exp_1",       kMagenta+1},
  { "Exp_2",       kRed+1},
  { "Exp_3",       kGreen+3},
  { "Exp_4",       kOrange+1},
  { "Exp_5",       kAzure+7},
  { "Exp_6",       kBlack},

  { "bkg_pred",    kBlue+1},
  { "data",        kRed+1},
  { "main_bkg_CR", kGreen+1},
  { "main_bkg_SR", kOrange+1},
  { "main_bkg",    kAzure+7},
  { "extra_bkg",   kOrange+1},
  { "DY_CR",       kGreen-2},
  { "DY_SR",       kOrange-2},
  { "DY",          kOrange-2},
  { "TTbar",       kOrange+10},
  { "VV",          kGreen+2},
  { "WW",          kGreen},
  { "WZ",          kGreen+3},
  { "ZZ",          kGreen-10},

  { "nominal",     kRed},
  { "JEC_up",      kOrange+1},
  { "JEC_down",    kOrange-2},
  { "JER_up",      kGreen+2},
  { "JER_down",    kGreen+1},
  { "PU_up",       kBlue+1},
  { "PU_down",     kAzure+10},
};



class CreateRooWorkspace {

public:
  // CreateRooWorkspace(std::string year, std::string collection, std::string channel, std::string histFolder);
  CreateRooWorkspace(std::string year, std::string collection, std::string channel, std::string histFolder, std::string min, std::string max);
  ~CreateRooWorkspace();

  void Process();
  void SetEnv();
  void LoadFiles();
  void LoadHistos();
  void PrepocessHistos();
  void NormaliseData();
  void DoRebin();
  void InitializePDFs();
  void CreateRooDataHist();
  void DoFits();
  void ImportToWorkspace();
  void DoPlots();
  void PlotBkgFit();
  void PlotSignals(std::string syst="nominal");
  void PlotSgPars();
  void PlotControl();
  void PlotTranferFunction();
  void InputDatacards();
  void CalculateSignalFittingRange(double mass, double& rangeLo, double& rangeHi, double& plotLo, double& plotHi, double& ymax, std::string& name);
  inline bool FindInVector(const std::vector<std::string>& vec, const std::string& str) {return (std::find(vec.begin(), vec.end(), str) != vec.end());};
  inline bool FindInString(const std::string& search, const std::string& str) {return str.find(search)!=std::string::npos ;}
  inline bool isNominalSyst(std::string syst) { return FindInString("nominal",syst);}; // TODO
  inline bool isNominalFolder(std::string syst) { bool var = isNominalSyst(syst); for (const auto& e: SystematicsScale) var+=FindInString(e,syst); return var;};


private:
  int myMin, myMax;

  std::string year, collection, channel, histFolder;
  std::string user, Path_ANALYSIS, Path_NFS, Path_STORAGE, PrefixrootFile;
  std::string workingDir, unique_name_complete, unique_name, filepath;
  std::string SRname, CRname;

  std::vector<std::string> SystNames = {"nominal", "all"};
  // std::vector<std::string> SystematicsScale = {"pu", "btag", "prefiring", "id", "isolation", "tracking", "trigger", "reco", "taggerSF", "mur", "muf", "murmuf" };
  std::vector<std::string> SystematicsScale = {"pu", "btag", "prefiring", "id", "isolation", "tracking", "trigger", "reco", "taggerSF", "murmuf", "NNPDF"};
  std::vector<std::string> SystematicsShape = {"JEC", "JER", "MuonScale"};
  std::vector<std::string> SystematicsAll;
  std::vector<std::string> Var_murmuf = {"upup", "upnone", "noneup", "nonedown", "downnone", "downdown"};
  const int PDF_variations = 100;


  // std::vector<std::string> BkgNames = {"DY", "WJets", "TTbar", "VV", "WZ", "WW", "ZZ"};
  std::vector<std::string> BkgNames = {"DY", "WJets", "TTbar"};
  std::vector<std::string> Modes = {"bkg_pred", "data", "data_CR", "MC_CR", "MC_fake_CR", "MC_SR", "DY_CR", "DY_SR", "WJets_CR", "WJets_SR"};

  std::string dataName, dataFileName;
  std::string BkgName = "DY";
  std::string FitSignal = "CB";

  std::unique_ptr<RooRealVar> x_var;
  std::unique_ptr<RooDataHist> data_obs;
  std::unique_ptr<RooPlot> plotter;

  std::unordered_map<std::string, std::unique_ptr<TH1F> > histo_map;
  std::unordered_map<std::string, std::unordered_map<std::string, bool> > doFits_map;
  std::unordered_map<std::string, std::unique_ptr<RooDataHist> > rooHist_map;
  std::unordered_map<std::string, std::unordered_map<std::string, std::unique_ptr<RooAbsPdf> > > Fits_map;
  std::unordered_map<std::string, std::unordered_map<std::string, std::unique_ptr<RooFitResult> > > FitRes_map;

  std::unordered_map<std::string, std::vector<std::unique_ptr<RooRealVar>> > fitPars;
  std::unordered_map<std::string, std::vector<double>> SgPars;
  std::vector<std::string> NameSgPars = {"Masses","nevents", "nevents_err", "mean","sigma","alpha","k","mean_err","sigma_err","alpha_err","k_err","chi2","pvalue","fit_min", "fit_max"};

  std::unordered_map<std::string, double> nEventsSignal;
  std::unordered_map<std::string, double> fit_min, fit_max, plot_min, plot_max, y_max;

  std::ofstream DataCard, output, SignalProperties;
  std::unique_ptr<RooWorkspace> ws;

  std::string studies="nominal";
  bool isHbb = false;
  bool fitMC = true;
  // bool doObs = false;
  bool doObs = true; //TODO


  std::string Module, Histtype, HistName;

  // bool dorebin = true;
  bool dorebin = false;
  bool doBinWidth = false;

  bool doCheckPlots = true;
  bool doBkgPlots = true;
  bool doFtest = false;
  // bool doPlotRatio = true;
  bool doPlotRatio = false;

  std::string plotting_mode = "pdf";
  // std::string plotting_mode_2 = "eps";
  std::string plotting_mode_2 = "";


  TString nameXaxis, nameYaxis, nameRatioaxis;

  TString extra_text = doFtest? "_Ftest_": "";

  // double x_lo     = 800;
  double x_lo     = 1000;
  double x_hi     = 6000;
  double x_lo_short = 1000;
  double x_hi_short = 6000;
  double plot_lo  = 1000;
  double plot_hi  = 6000;
  double plot_ylo = 0.101;
  double plot_yhi = 1e07;
  // double fit_lo   = 600;
  // double fit_lo   = 1230;
  // double fit_hi   = 3500;

  // double fit_lo   = 1360; DY_CR bin30 ok also for DR_SR 30
  // double fit_hi   = 5500; DY_CR bin30 ok also for DR_SR 30
  // double fit_lo   = 1200; data_CR bin30
  // double fit_hi   = 4000; data_CR bin30
  // double fit_lo   = 1300; DR_SR bin30
  // double fit_hi   = 4500; DR_SR bin30
  // double fit_SR   = 810;
  // double fit_min_CR = 1400;
  // double fit_max_CR = 5000;

  double fit_lo_SR;
  double fit_lo_CR;
  double fit_hi_SR;
  double fit_hi_CR;
  double show_lo_SR;
  double show_hi_SR;
  double show_lo_CR;
  double show_hi_CR;

  const std::unordered_map<std::string, std::unordered_map<std::string, std::unordered_map<std::string, double>>> ranges = {
    { "SR", {
      { "muonchannel", {
        { "fit_lo",  1000},
        { "fit_hi",  3000},
        { "show_lo", 1000},
        { "show_hi", 4500},
      }},
      { "electronchannel", {
        { "fit_lo",  900},
        { "fit_hi",  2800},
        { "show_lo", 1000},
        { "show_hi", 4500},
      }},
      { "invisiblechannel", {
        { "fit_lo",  1000},
        { "fit_hi",  3000},
        { "show_lo", 700},
        { "show_hi", 3150},
      }},
    }},
    { "CR", {
      { "muonchannel", {
        { "fit_lo",  1300},
        { "fit_hi",  4200},
        { "show_lo", 1200},
        { "show_hi", 6000},
        { "show_lo_data", 1200},
        { "show_hi_data", 3500},
      }},
      { "electronchannel", {
        { "fit_lo",  1300},
        { "fit_hi",  3300},
        { "show_lo", 1200},
        { "show_hi", 6000},
        { "show_lo_data", 1200},
        { "show_hi_data", 4500},
      }},
      { "invisiblechannel", {
        { "fit_lo",  1300},
        { "fit_hi",  4200},
        { "show_lo", 1200},
        { "show_hi", 6000},
        { "show_lo_data", 1200},
        { "show_hi_data", 3500},
      }},
    }},
  };

  // double fit_SR   = 1000;
  double fit_min_CR = 1400;
  double fit_max_CR = 6000;

  bool debug = false;
  // bool debug = true;

  int rebin;
  int bin;
  int bin2;
  std::vector<double> bins_Zprime_rebin;

  double nEventsSR, nEventsSR_fit, xsec_ref_;


};
