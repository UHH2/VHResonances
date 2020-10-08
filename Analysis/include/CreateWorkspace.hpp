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

#include "RooNovosibirsk.h"
#include "HiggsAnalysis/CombinedLimit/interface/RevCrystalBall.hpp"
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
/  fitCR: bkg_pred= Data_CR* MC_SR/MC_CR  /
/ !fitCR: bkg_pred= MC_SR* MC_SR/MC_SR    /
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
  { "BkgPdf3p",    kGreen+1},
  { "BkgPdf4p",    kRed+1},

  { "NO",          kBlue+1},
  { "CB",          kOrange+1},
  { "Exp_1",       kMagenta+1},
  { "Exp_2",       kRed+1},
  { "Exp_3",       kGreen+2},
  { "Exp_4",       kGreen+1},
  { "Exp_5",       kYellow+1},
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
  CreateRooWorkspace(std::string year, std::string collection, std::string channel, std::string histFolder);
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
  void CalculateSignalFittingRange(double mass, double& rangeLo, double& rangeHi, double& plotLo, double& plotHi, double& ymax);
  inline bool FindInVector(const std::vector<std::string>& vec, const std::string& str) {return (std::find(vec.begin(), vec.end(), str) != vec.end());};
  inline bool FindInString(const std::string& search, const std::string& str) {return str.find(search)!=std::string::npos ;}
  inline bool isNominalFolder(std::string syst) {return (isNominalSyst(syst) || FindInString("pu",syst) || FindInString("btag",syst) || FindInString("prefiring",syst) || FindInString("id",syst) || FindInString("isolation",syst) || FindInString("tracking",syst) || FindInString("trigger",syst) || FindInString("reco",syst));};
  inline bool isNominalSyst(std::string syst) { return FindInString("nominal",syst);}; // TODO

private:
  std::string year, collection, channel, histFolder;
  std::string user, Path_ANALYSIS, Path_NFS, Path_STORAGE, PrefixrootFile;
  std::string workingDir, unique_name_complete, unique_name, filepath;
  std::string SRname, CRname;

  std::vector<std::string> SystNames = {"nominal", "all"};
  // std::vector<std::string> SystNames = {"nominal"};
  // std::vector<std::string> BkgNames = {"DY", "TTbar", "WZ", "WW", "ZZ"}; //TODO
  std::vector<std::string> BkgNames = {"DY", "TTbar", "WZ","ZZ"};
  std::vector<std::string> Modes = {"bkg_pred", "data", "main_bkg_CR", "main_bkg_SR", "DY_CR", "DY_SR"};

  //TODO
  // const std::unordered_map<int, int> colors = {{600, kRed}, {800, kGreen}, {1000, kViolet}, {1200, kBlue}, {1400, kBlack}, {1600, kOrange}, {1800, kAzure}, {2000, kSpring}, {2500, kPink},
  // {3000, kRed}, {3500, kGreen}, {4000, kViolet}, {4500, kBlue}, {5000, kBlack}, {5500, kOrange}, {6000, kAzure}, {7000, kSpring}, {8000, kPink}};

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

  // TFile file_WS("WS.root","RECREATE");
  std::ofstream DataCard, output, SignalProperties;
  std::unique_ptr<RooWorkspace> ws;

  // TODO
  std::string studies="nominal";
  bool isHbb = false;
  bool fitCR = true;
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

  double x_lo     = 1000;
  double x_hi     = 10000;
  double plot_lo  = 1000;
  double plot_hi  = 4000;
  double plot_ylo = 1.1*1e-03;
  double plot_yhi = 1e07;
  // double fit_lo   = 600;
  double fit_lo   = 1100;
  double fit_hi   = 4000;
  double fit_SR   = 810;

  bool debug = false;
  // bool debug = true;



  int rebin;
  int bin;
  int bin2;
  std::vector<double> bins_Zprime_rebin;

  double nEventsSR, xsec_ref_;


};
