#pragma once
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <limits>

#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TH1D.h"
#include "THStack.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TFile.h"
#include "TString.h"
#include "TObjString.h"
#include "TLatex.h"
#include "UHH2/core/include/Utils.h"
#include "UHH2/VHResonances/include/constants.hpp"
#include "UHH2/VHResonances/Analysis/include/tdrstyle_all.hpp"

enum Type{typeMC=0, typeSignal=1, typeData=2};

template<typename T>
void PrintGreen(const T& x) {std::cout << green << x << reset << std::endl;}

template<typename T>
void PrintGreen(std::vector<T> vec) {std::cout << green; for (auto& x: vec) {std::cout<< x << " ";}; std::cout<< reset << std::endl;}

template<typename T>
void PrintInfo(const std::string& name, const T& x, int tot = 25, int start = 4) {
  int nspace = tot-start-2-name.size();
  if (nspace<0) nspace=0;
  std::cout << green << std::string(start,' ')+name+std::string(nspace,' ')+": " << x << reset << std::endl;
}

inline bool FindInVector(const std::vector<std::string>& vec, const std::string& str) {return (std::find(vec.begin(), vec.end(), str) != vec.end());};
inline bool FindInString(const std::string& search, const std::string& str) {return str.find(search)!=std::string::npos ;}

class SampleInfo {
public:
  SampleInfo(std::string uniqueName, std::string fileName, double weight, Type type, std::string legName, int color, int linestyle);
  ~SampleInfo() {};
  void SetUniqueName(std::string x) {UniqueName = x;}
  void SetLegName(std::string x) {LegName = x;}
  void SetFileName(std::string x) {FileName = x;}
  void SetColor(int x) {Color = x;}
  void SetLinestyle(int x) {Linestyle = x;}
  void SetWeight(double x) {weight = x;}
  void SetIsMC()     {IsMC = true;  IsData = false; IsSignal = false; IsToStack = true;}
  void SetIsData()   {IsMC = false; IsData = true;  IsSignal = false; IsToStack = false;}
  void SetIsSignal() {IsMC = true;  IsData = false; IsSignal = true;  IsToStack = false;}

  std::string GetUniqueName() const {return UniqueName;}
  std::string GetLegName() const {return LegName;}
  std::string GetFileName() const {return FileName;}
  int GetColor() const {return Color;}
  int GetLinestyle() const {return Linestyle;}
  double GetWeight() const {return weight;};
  bool GetIsToStack() const {return IsToStack;}
  bool GetIsMC() const {return IsMC;}
  bool GetIsData() const {return IsData;}
  bool GetIsSignal() const {return IsSignal;}

private:
  std::string UniqueName;
  std::string LegName;
  std::string FileName;
  int Color;
  int Linestyle;
  double weight;
  bool IsMC;
  bool IsData;
  bool IsSignal;
  bool IsToStack;

};

class Plotter {

public:
  Plotter(std::string module, std::string hname, std::string year, std::string channel, std::string collections="Puppi");
  ~Plotter();

  void SetEnv();

  void AddSample(std::string uniqueName, std::string fileName, double weight, Type type, std::string legName, int color, int linestyle);
  void LoadHists();
  void MakeStackHist();
  void MakePlot();
  void FindRanges();
  void Process();

  void swapLegend(bool x) {swapLegend_ = x;};
  bool isToSwapLegend() {return swapLegend_;};

  void SetXmin(double x) {xmin = x; xmin_set=xmin>0;};
  void SetXmax(double x) {xmax = x; xmax_set=xmax>0;};
  void SetYmin(double x) {ymin = x; ymin_set=ymin>0;};
  void SetYmax(double x) {ymax = x; ymax_set=ymax>0;};
  void SetXRange(double min, double max) {SetXmin(min); SetXmax(max);};
  void SetYRange(double min, double max) {SetYmin(min); SetYmax(max);};
  void SetXTitle(std::string x) {xTitle = TString(x);};

  void SetRebin(int x) {rebin = x;};
  int  GetRebin() const {return rebin;};

  bool IsSetXmin() {return xmin_set;};
  bool IsSetXmax() {return xmax_set;};
  bool IsSetYmin() {return ymin_set;};
  bool IsSetYmax() {return ymax_set;};

  std::string GetChannel() const {return channel;};
  std::string GetModule() const {return module;};
  std::string GetYear() const {return year;};

  inline bool isNominalSyst(std::string syst) { return FindInString("nominal",syst);};
  inline bool isNominalFolder(std::string syst) { bool var = isNominalSyst(syst); for (const auto& e: SystematicsScale) var+=FindInString(e,syst); return var;};


private:
  std::string module;
  std::string hname;
  std::string year;
  std::string channel;
  std::string collections;
  std::string user;
  std::string Path_ANALYSIS;
  std::string Path_NFS;
  std::string Path_STORAGE;
  std::string PrefixrootFile;
  std::string inputdir;
  std::string outputdir;
  std::string relPath;
  double lumi_unc;

  std::unordered_map<std::string, std::unique_ptr<SampleInfo>> Samples;
  std::unordered_map<std::string, std::unique_ptr<TH1D> > histos;
  std::unordered_map<std::string, std::unordered_map<std::string, std::unique_ptr<TH1D>> > uncertainties;
  std::vector<std::string> histnames;
  std::vector<std::string> histstacknames;
  std::unique_ptr<THStack> stack;
  std::unique_ptr<THStack> stack_unc;
  std::unique_ptr<TH1F> h_err;
  std::unique_ptr<TH1F> h_err_syst;
  std::unordered_map<std::string, std::unique_ptr<THStack>> stacksts;

  TString xTitle;

  double xmin=std::numeric_limits<double>::max();
  double xmax=-std::numeric_limits<double>::max();
  double ymin=std::numeric_limits<double>::max();
  double ymax=-std::numeric_limits<double>::max();

  bool swapLegend_ = false;

  std::vector<std::string> systematics = {"lumi", "JEC", "JER","MuonScale"};
  // std::vector<std::string> systematics = {"lumi"};
  std::vector<std::string> SystematicsScale = {"pu", "btag", "prefiring", "id", "isolation", "tracking", "trigger", "reco", "taggerSF", "murmuf", "NNPDF"};
  // std::vector<std::string> SystematicsScale = {};
  std::vector<std::string> Var_murmuf = {"upup", "upnone", "noneup", "nonedown", "downnone", "downdown"};
  const int PDF_variations = 100;
  bool xmin_set = false;
  bool xmax_set = false;
  bool ymin_set = false;
  bool ymax_set = false;
  int rebin=1;
  int nbins=1;
  double binwidth=0;


  // bool debug=true;
  bool debug=false;
  bool plotSyst=false;
  // bool plotSyst=true;

};
