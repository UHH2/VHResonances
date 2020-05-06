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

void CreateWorkspace(std::string studies="nominal", std::string histFolder = "btag_DeepBoosted_H4qvsQCD", std::string channel = "muonchannel", std::string collection = "Puppi", std::string year = "2016", bool isHbb = false, bool fitCR = true, bool doObs = false);

RooRealVar* ReDoVar(RooRealVar& var, TString& name);

double GetRange(TH1F* h, double x);

double CalculateIntegral(TH1F* h, double min, double max, bool dorebin);
double CalculateFractionArea(TH1F* h, double min, double max, double x_min, double x_max, bool dorebin);

double DoFTest(double chi2_1, double chi2_2, double npar_1, double npar_2, double n);
void CalculateChiSquare(double& chi2, int& nbins, RooHist* hpull, double xmin, double xmax);
