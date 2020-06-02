#include "CreateWorkspace.hpp"
#include "TSpline.h"

double GetRange(TH1F* h, double x) { return h->GetXaxis()->GetBinLowEdge(h->FindBin(x)); }

double CalculateIntegral(TH1F* h, double min, double max, bool doBinWidth) { return h->Integral(h->FindBin(min),h->FindBin(max),doBinWidth?"width":"");}
double CalculateFractionArea(TH1F* h, double min, double max, double x_min, double x_max, bool doBinWidth) { return CalculateIntegral(h, min, max, doBinWidth)/CalculateIntegral(h, x_min, x_max, doBinWidth);}

double DoFTest(double chi2_1, double chi2_2, double npar_1, double npar_2, double n) {return 1. - TMath::FDistI(((chi2_1-chi2_2)/(npar_2-npar_1))/(chi2_2/(n-npar_2-1)), npar_2-npar_1, n-npar_2); }

void CalculateChiSquare(double& chi2, int& nbins, RooHist* hpull, double xmin, double xmax) {
  // This method should take only the fitted points. plotter->chiSquare(npf) is on the full range.
  // diff in chi2 wrt to sum(hpull) are due to avarage(pull) vs interpolate(chiSquare) values in calculations
  chi2 = 0; nbins=0;
  for(int i=0 ; i<hpull->GetN() ; i++) {
    double x,pull;
    hpull->GetPoint(i,x,pull);
    if (x<xmin || x>xmax) continue;
    // std::cout << "FIND ME2 " << x << " " << pull << " " << pull*pull << " " << chi2 << '\n';
    chi2+=pull*pull;
    nbins++;
  }
}


void CreateRooWorkspace::CalculateSignalFittingRange(double mass, double& rangeLo, double& rangeHi, double& plotLo, double& plotHi, double& ymax) {


  // rangeLo = mass*(1-5./28.);
  // rangeHi = mass*(1+3.5/28.);
  rangeLo = 0.776*mass-56;
  rangeHi = 1.085*mass+46;

  if (mass==1200) rangeLo = 790;
  if (mass==1200) rangeHi = 1390;


  plotLo = mass*(1-10./28.);
  plotHi = mass*(1+10./28.);
  ymax = (50-mass*1./100)*1.5;
  if (mass<1200) ymax = (-20+mass*6./100)*1.4;
  if (mass>3500) ymax = (30-mass*3./1000)*1.3;

  // fit mass:600	390 -- 690
  // fit mass:800	540 -- 900
  // fit mass:1000	720 -- 1110
  // fit mass:1200	870 -- 1320
  // fit mass:1400	1020 -- 1560
  // fit mass:1600	1170 -- 1770
  // fit mass:1800	1320 -- 1980
  // fit mass:2000	1470 -- 2190
  // fit mass:2500	1860 -- 2730
  // fit mass:3000	2250 -- 3300
  // fit mass:3500	2640 -- 3840
  // fit mass:4000	3030 -- 4380
  // fit mass:4500	3420 -- 4920
  // fit mass:5000	3810 -- 5460
  // fit mass:5500	4200 -- 6000
  // fit mass:6000	4590 -- 6540
  // fit mass:7000	5370 -- 7620
  // fit mass:8000	6150 -- 8700

  ymax *= lumi_map.at(year).at("lumi_fb")/lumi_map.at("RunII").at("lumi_fb");
  if (year=="2017") ymax *= 1.1;

  ymax *= xsec_ref_/0.1;

  if ((histFolder.find("ptdep")!=std::string::npos || histFolder.find("massdep")!=std::string::npos)) ymax *= 6;


  rangeLo = GetRange(histo_map[GetSgName(mass)].get(), rangeLo);
  rangeHi = GetRange(histo_map[GetSgName(mass)].get(), rangeHi);
  plotLo  = GetRange(histo_map[GetSgName(mass)].get(), plotLo);
  plotHi  = GetRange(histo_map[GetSgName(mass)].get(), plotHi);

  // rangeLo = GetRange(histo_map[GetSgName(mass)].get(), 0.8*mass-88);
  // rangeHi = GetRange(histo_map[GetSgName(mass)].get(), 1.1*mass+18);

  // rangeLo = GetRange(histo_map[GetSgName(mass)].get(), 0.8*mass-60);
  // rangeHi = GetRange(histo_map[GetSgName(mass)].get(), 1.1*mass+44);

  // rangeLo = GetRange(histo_map[GetSgName(mass)].get(), 0.8*mass-103);
  // rangeHi = GetRange(histo_map[GetSgName(mass)].get(), 1.1*mass+18);

  // 2016 ele
  // if (year=="2016" && channel=="electronchannel") rangeLo = GetRange(histo_map[GetSgName(mass)].get(), 0.8*mass-88);
  // if (year=="2016" && channel=="electronchannel") rangeHi = GetRange(histo_map[GetSgName(mass)].get(), 1.1*mass+18);
  //   fit_max p0 = 43.636 +- 31.361
  // fit_max p1 = 1.086 +- 0.008
  // fit_min p0 = -60.489 +- 26.919
  // fit_min p1 = 0.777 +- 0.007
  //
  // fit_max p0 = 17.868 +- 31.564
  // fit_max p1 = 1.099 +- 0.008
  // fit_min p0 = -103.616 +- 27.194
  // fit_min p1 = 0.798 +- 0.007

}

std::string GetSgName(int mass) { return "M"+std::to_string(mass); }


CreateRooWorkspace::CreateRooWorkspace(std::string year_, std::string collection_, std::string channel_, std::string histFolder_) : year(year_), collection(collection_), channel(channel_), histFolder(histFolder_) {
  /*
  #  ######   #######  ##    ##  ######  ######## ########  ##     ##  ######  ########  #######  ########
  # ##    ## ##     ## ###   ## ##    ##    ##    ##     ## ##     ## ##    ##    ##    ##     ## ##     ##
  # ##       ##     ## ####  ## ##          ##    ##     ## ##     ## ##          ##    ##     ## ##     ##
  # ##       ##     ## ## ## ##  ######     ##    ########  ##     ## ##          ##    ##     ## ########
  # ##       ##     ## ##  ####       ##    ##    ##   ##   ##     ## ##          ##    ##     ## ##   ##
  # ##    ## ##     ## ##   ### ##    ##    ##    ##    ##  ##     ## ##    ##    ##    ##     ## ##    ##
  #  ######   #######  ##    ##  ######     ##    ##     ##  #######   ######     ##     #######  ##     ##
  */


  std::cout << "****************************************" << std::endl;
  std::cout << "             CreateWorkspace            " << std::endl;
  std::cout << "----------------------------------------" << std::endl;
  std::cout << "year"       << std::string(15-4, ' ') << year       << '\n';
  std::cout << "channel"    << std::string(15-7, ' ') << channel    << '\n';
  std::cout << "collection" << std::string(15-10,' ') << collection << '\n';
  std::cout << "histFolder" << std::string(15-10,' ') << histFolder << '\n';
  std::cout << "****************************************\n" << std::endl;

  SetEnv();

};

CreateRooWorkspace::~CreateRooWorkspace(){
  ws->writeToFile((workingDir+"/datacards/ws_"+histFolder+".root").c_str());
  DataCard.close();
  if (doFtest) output.close();
  SignalProperties.close();
  std::cout << "End of CreateRooWorkspace" << '\n';
}


void CreateRooWorkspace::SetEnv() {

  writeExtraText = true;       // if extra text
  extraText  = "Work in progress" ;//"Preliminary";
  lumi_13TeV  = TString::Format("%.1f fb^{-1}", lumi_map.at(year).at("lumi_fb"));

  user            = std::getenv("USER");
  Path_ANALYSIS   = std::getenv("CMSSW_BASE"); Path_ANALYSIS += "/src/UHH2/VHResonances/";
  Path_NFS        = "/nfs/dust/cms/user/"+user+"/";
  Path_STORAGE    = Path_NFS+"WorkingArea/File/Analysis/";
  PrefixrootFile  = "uhh2.AnalysisModuleRunner.";

  workingDir      = Path_ANALYSIS+"Analysis/"+(isHbb? "Limits/Hbb_": "Limits/")+studies+"/"+year+"/"+collection+"/"+channel+"/"+histFolder+"/";
  gSystem->Exec(("mkdir -p "+workingDir+"/datacards").c_str());

  Module    = "SignalRegion";
  Histtype  = "ZprimeCandidate";
  // HistName  = "Zprime_mass_rebin_full";
  // HistName  = "Zprime_mass";
  // HistName  = "Zprime_mass_rebin1";
  // HistName  = "Zprime_mass_rebin2";
  HistName  = "Zprime_mass_rebin30"; // TODO


  unique_name_complete = year+"_"+collection+"_"+channel+"_"+histFolder;
  unique_name = "_"+year+"_"+collection+"_"+channel;
  filepath    = Path_STORAGE+year+"/"+Module+"/"+collection+"/"+channel+"/nominal/";


  dataName  = channel.substr(0, channel.find("channel")); dataName[0] = toupper(dataName[0]) ; dataName = "DATA_Single"+dataName;
  dataFileName = PrefixrootFile+"DATA."+dataName+"_"+year+"_noTree.root";

  x_var.reset(new RooRealVar("x_var", "m_{Zprime} (GeV)", x_lo, x_hi));
  ws.reset(new RooWorkspace((channel+"_"+year).c_str()));
  DataCard.open (workingDir+"datacards/OutputFit_"+histFolder+".txt");
  DataCard << "=== RooFit data fit result to be entered in datacard === " << std::endl;
  if (doFtest) output.open(workingDir+"datacards/output_"+year+"_"+histFolder+".txt");
  SignalProperties.open(workingDir+"datacards/SignalProperties_"+year+"_"+histFolder+".txt");




  // rebin = 30;
  rebin = 0;
  bin = 50;
  bin2 = 100;
  for (double i = 0; i < 1500; i+=bin) bins_Zprime_rebin.push_back(i);
  for (double i = 1500+bin; i <= 4000; i+=bin2) bins_Zprime_rebin.push_back(i);
  bins_Zprime_rebin.push_back(4000);

  for (std::string name : NameSgPars) SgPars[name] = std::vector<double>(MassPoints.size(), 0);


  // xsec_ref_ = 0.1; // this is to mantain the signal strenght close to 1; (remember to multiply for this Normalization when plotting)
  xsec_ref_ = (xsec_ref.find(histFolder) != xsec_ref.end())? xsec_ref.at(histFolder): xsec_ref.at("default_value"); // default_value = 1
  DataCard << "xsec_ref " << " " << xsec_ref_ <<std::endl;


  for (auto x : Colors) { for (auto mode : Modes) doFits_map[mode][x.first]  = false; }

  // for (auto x : { "NO", "CB", "Exp_1", "Exp_2", "Exp_3", "Exp_4", "Exp_5", "Exp_6"}) { for (auto mode : Modes) doFits_map[mode][x]  = true; }
  // for (auto x : {"CB", "Exp_1", "Exp_2", "Exp_3" }) { for (auto mode : Modes) doFits_map[mode][x]  = true; }
  for (auto x : {"Exp_2", "Exp_3" }) { for (auto mode : Modes) doFits_map[mode][x]  = true; }
  // for (auto x : { "NO", "CB", "Exp_3"}) { for (auto mode : Modes) doFits_map[mode][x]  = true; }
  // BkgPdf4p, BkgPdf3p TODO
  // for (auto x : { "NO", "CB", "Exp_4"}) { for (auto mode : Modes) doFits_map[mode][x]  = true; }
  // for (auto x : { "CB"}) { for (auto mode : Modes) doFits_map[mode][x]  = true; }


};




void CreateRooWorkspace::LoadFiles() {
  /*
  # ##        #######     ###    ########     ##     ## ####  ######  ########  #######   ######   ########     ###    ##     ##  ######
  # ##       ##     ##   ## ##   ##     ##    ##     ##  ##  ##    ##    ##    ##     ## ##    ##  ##     ##   ## ##   ###   ### ##    ##
  # ##       ##     ##  ##   ##  ##     ##    ##     ##  ##  ##          ##    ##     ## ##        ##     ##  ##   ##  #### #### ##
  # ##       ##     ## ##     ## ##     ##    #########  ##   ######     ##    ##     ## ##   #### ########  ##     ## ## ### ##  ######
  # ##       ##     ## ######### ##     ##    ##     ##  ##        ##    ##    ##     ## ##    ##  ##   ##   ######### ##     ##       ##
  # ##       ##     ## ##     ## ##     ##    ##     ##  ##  ##    ##    ##    ##     ## ##    ##  ##    ##  ##     ## ##     ## ##    ##
  # ########  #######  ##     ## ########     ##     ## ####  ######     ##     #######   ######   ##     ## ##     ## ##     ##  ######
  */


  std::unordered_map<std::string, std::unique_ptr<TFile>> f_map;

  SRname = Histtype+"_"+histFolder+"_SR/"+HistName;
  CRname = Histtype+"_"+histFolder+"_CR/"+HistName;

  if (debug) {
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "Opening\t" << filepath+dataFileName << '\n';
  }

  f_map["data"].reset(new TFile((filepath+dataFileName).c_str()));

  for (auto bkg: BkgNames) {
    std::string fname = filepath+PrefixrootFile+"MC.MC_"+bkg+"_"+year+"_noTree.root";
    if (debug) std::cout << "Opening\t" << fname << '\n';
    f_map[bkg].reset(new TFile(fname.c_str()));
    histo_map[bkg+"_CR"].reset((TH1F*)((TH1F*)f_map[bkg]->Get((CRname).c_str()))->Clone((bkg+" CR").c_str()));
    histo_map[bkg+"_SR"].reset((TH1F*)((TH1F*)f_map[bkg]->Get((SRname).c_str()))->Clone((bkg+" SR").c_str()));
    histo_map[bkg+"_ratio"].reset((TH1F*)(histo_map[bkg+"_SR"])->Clone((bkg+"_ratio").c_str()));
    for (std::string reg: {"CR","SR"}) {
      if (bkg==BkgName) {
        histo_map["main_bkg_"+reg].reset((TH1F*)((TH1F*)f_map[BkgName]->Get((CRname).c_str()))->Clone(("main_bkg_"+reg+": "+BkgName).c_str()));
        continue;
      }
      histo_map["main_bkg_"+reg]->Add(histo_map[bkg+"_"+reg].get());
      histo_map["main_bkg_"+reg]->SetName((std::string(histo_map["main_bkg_"+reg]->GetName())+"+"+bkg).c_str());
      if (bkg=="TTbar") {
        histo_map["extra_bkg_"+reg].reset((TH1F*)((TH1F*)f_map["TTbar"]->Get((CRname).c_str()))->Clone(("extra_bkg_"+reg+": TTbar").c_str()));
        continue;
      }
      histo_map["extra_bkg_"+reg]->Add(histo_map[bkg+"_"+reg].get());
      histo_map["extra_bkg_"+reg]->SetName((std::string(histo_map["main_bkg_"+reg]->GetName())+"+"+bkg).c_str());
      if (bkg=="WZ") {
        histo_map["VV_"+reg].reset((TH1F*)((TH1F*)f_map["WZ"]->Get((CRname).c_str()))->Clone(("extra_bkg_"+reg+": VV").c_str()));
        continue;
      }
      histo_map["VV_"+reg]->Add(histo_map[bkg+"_"+reg].get());
    }
    // TODO
    // for (auto syst: SystNames) {
    //   TString fname = filepath+PrefixrootFile+"MC.MC_"+bkg+"_"+year+"_noTree.root";
    //   fname = fname.ReplaceAll("nominal",syst);
    //   if (debug) std::cout << "Opening\t" << fname << '\n';
    //   f_map[bkg+syst] = new TFile(fname);
    // }
  }

  for (const int & mass : MassPoints) {

    std::string SgName = GetSgName(mass);
    std::string fname = filepath+PrefixrootFile+"MC.MC_ZprimeToZH_"+SgName+"_"+year+"_noTree.root";
    if (debug) std::cout << "Opening\t" << fname << '\n';
    f_map[SgName].reset(new TFile(fname.c_str()));
    // histo_map[SgName] = (TH1F*)f_map[SgName]->Get(SRname.c_str());
    histo_map[SgName].reset((TH1F*)((TH1F*)f_map[SgName]->Get((SRname).c_str()))->Clone((SgName+" SR").c_str()));
    // TODO
    // for (auto syst: SystNames) {
    //   TString fname = filepath+PrefixrootFile+"MC.MC_ZprimeToZH_"+SgName+"_"+year+"_noTree.root";
    //   fname = fname.ReplaceAll("nominal",syst);
    // if (debug) std::cout << "Opening\t" << fname << '\n';
    //
    //   histo_map[SgName+syst]=(TH1F*)TFile(fname).Get(SRname.c_str());
    //   std::cout << "FIND ME " << histo_map[SgName+syst] << '\n';
    //   // histo_map[SgName+syst]->SetDirectory(0);
    // }


    // if (isNominalFolder) Histtype  += "_"+syst;
    //
    // SRname = Histtype+"_HWW"+histFolder+"_SR/"+HistName; // WE DON'T WANT THE BR INVOLVED
    // SRname = Histtype+"_"+histFolder+"_SR/"+HistName;
    // SRname = Histtype+"_"+histFolder+"_SR/"+HistName;
    // CRname = Histtype+"_"+histFolder+"_CR/"+HistName;
    //
    // if (isNominalFolder) SRname = Histtype+histFolder+"_SR/"+HistName;
    // if (isNominalFolder) SRname = Histtype+histFolder+"_SR/"+HistName;
    // if (isNominalFolder) CRname = Histtype+histFolder+"_CR/"+HistName;
  }

  histo_map["norm"].reset((TH1F*)((TH1F*)f_map["data"]->Get((SRname).c_str()))->Clone("data: data SR"));

  if (doObs) histo_map["data"].reset((TH1F*)((TH1F*)f_map["data"]->Get((SRname).c_str()))->Clone("data: data SR"));
  else histo_map["data"].reset((TH1F*)((TH1F*)f_map[BkgName]->Get((SRname).c_str()))->Clone(("h_data_fake: "+BkgName+" SR").c_str()));

  if (fitCR) histo_map["bkg_pred"].reset((TH1F*)((TH1F*)f_map["data"]->Get((CRname).c_str()))->Clone("bkg_pred: data CR"));
  else histo_map["bkg_pred"].reset((TH1F*)((TH1F*)f_map[BkgName]->Get((SRname).c_str()))->Clone(("bkg_pred: "+BkgName+" SR").c_str()));


  for (std::string name : {"main", "extra"}) histo_map[name+"_bkg_ratio"].reset((TH1F*)(histo_map[name+"_bkg_SR"])->Clone((name+"_bkg_ratio").c_str()));

  // Need to make the histogram available after closing the files.
  for (const auto& x : histo_map) x.second->SetDirectory(0);

  if (debug){
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "CR hist_name:\t" << CRname << '\n';
    std::cout << "SR hist_name:\t" << SRname << '\n';
    std::cout << "fitCR:\t"        << BoolToString(fitCR) << '\n';
    std::cout << "doObs:\t"        << BoolToString(doObs) << '\n';
    std::cout << "bkg_pred on: "   << histo_map["bkg_pred"]->GetName() << '\n';
    std::cout << "data on: "       << histo_map["data"]->GetName() << '\n';
    std::cout << "----------------------------------------" << std::endl;
  }

  if (debug) for (const auto& x : histo_map) std::cout << "histo " << x.first << " " << std::string(20-x.first.size(),' ') << x.second->Integral() << '\n';

}

void CreateRooWorkspace::PrepocessHistos() {
  /*
  # ########  ########  ######## ########   #######   ######  ########  ######   ######     ##     ## ####  ######  ########  #######   ######
  # ##     ## ##     ## ##       ##     ## ##     ## ##    ## ##       ##    ## ##    ##    ##     ##  ##  ##    ##    ##    ##     ## ##    ##
  # ##     ## ##     ## ##       ##     ## ##     ## ##       ##       ##       ##          ##     ##  ##  ##          ##    ##     ## ##
  # ########  ########  ######   ########  ##     ## ##       ######    ######   ######     #########  ##   ######     ##    ##     ##  ######
  # ##        ##   ##   ##       ##        ##     ## ##       ##             ##       ##    ##     ##  ##        ##    ##    ##     ##       ##
  # ##        ##    ##  ##       ##        ##     ## ##    ## ##       ##    ## ##    ##    ##     ##  ##  ##    ##    ##    ##     ## ##    ##
  # ##        ##     ## ######## ##         #######   ######  ########  ######   ######     ##     ## ####  ######     ##     #######   ######
  */

  for (const auto& x : histo_map) {
    std::string mode = x.first;
    if (debug) std::cout << "pre " << mode << "\t" << CalculateIntegral(x.second.get(),fit_lo,fit_hi,doBinWidth) << '\n';
    if (dorebin) {
      if (rebin) histo_map[mode]->Rebin(rebin);
      else histo_map[mode].reset(dynamic_cast<TH1F*>(histo_map[mode]->Rebin(bins_Zprime_rebin.size()-1, x.second->GetName(), &bins_Zprime_rebin[0])));

      histo_map[mode]->Scale(1,doBinWidth?"width":"");
    }
    if (debug) std::cout << "rebin " << mode << "\t" << CalculateIntegral(x.second.get(),fit_lo,fit_hi,doBinWidth) << '\n';

    // Removing bins with low stat TODO
    if (mode=="DY_SR") {
      for (int i = 0; i < x.second->GetNbinsX()+1; i++) {
        if (x.second->GetBinContent(i)<2*1e-02) { histo_map[mode]->SetBinContent(i,0); histo_map[mode]->SetBinError(i,0); } // TODO
        // if (h->GetBinCenter(i)>1280 && h->GetBinCenter(i)<1320) { h->SetBinContent(i,0); h->SetBinError(i,0); }
      }
    }
  }

  for (const auto& x : histo_map) {
    if (x.first.find("ratio") != std::string::npos) {
      histo_map[x.first]->Divide(histo_map[TString(histo_map[x.first]->GetName()).ReplaceAll("SR","CR").Data()].get());
    }
  }


  // Normalize signal to arbitraty xsec.
  for (const int & mass : MassPoints) histo_map[GetSgName(mass)]->Scale(xsec_ref_);

}

void CreateRooWorkspace::NormaliseData() {

  /*
  # ########     ###    ########    ###              #######  ########   ######
  # ##     ##   ## ##      ##      ## ##            ##     ## ##     ## ##    ##
  # ##     ##  ##   ##     ##     ##   ##           ##     ## ##     ## ##
  # ##     ## ##     ##    ##    ##     ##          ##     ## ########   ######
  # ##     ## #########    ##    #########          ##     ## ##     ##       ##
  # ##     ## ##     ##    ##    ##     ##          ##     ## ##     ## ##    ##
  # ########  ##     ##    ##    ##     ##           #######  ########   ######
  */

  // Normalize h_MC_SR to h_Data_SR pretending it's data but has shape of bkg_pred in SR TODO
  nEventsSR  = CalculateIntegral(histo_map["norm"].get(),fit_lo,fit_hi,doBinWidth);
  if (!doObs) {
    histo_map["data"]->Scale(nEventsSR/CalculateIntegral(histo_map["data"].get(),fit_lo,fit_hi,doBinWidth));
    // histo_map["data"]->Scale(1./CalculateFractionArea(histo_map["data"],fit_lo,fit_hi, x_lo, x_hi,doBinWidth));
    histo_map["data"]->SetName((std::string(histo_map["data"]->GetName())+"*nEventsSR").c_str());
  }

  data_obs.reset(new RooDataHist("data_obs", histo_map["data"]->GetName(), RooArgList(*x_var), histo_map["data"].get()));//TODO check it's imported correctly

  if (debug) {
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "data_obs\t"  << histo_map["data"]->GetName() << '\n';
    std::cout << "nEventsSR\t" << nEventsSR << '\n';
    std::cout << "norm_data\t" << CalculateIntegral(histo_map["data"].get(),fit_lo,fit_hi,doBinWidth) << '\n';
    std::cout << "----------------------------------------" << std::endl;
  }

  DataCard << BkgName+" Background number of events = " << nEventsSR << std::endl;

}


void CreateRooWorkspace::DoRebin() {
  /*
  # ########   #######     ########  ######## ########  #### ##    ##
  # ##     ## ##     ##    ##     ## ##       ##     ##  ##  ###   ##
  # ##     ## ##     ##    ##     ## ##       ##     ##  ##  ####  ##
  # ##     ## ##     ##    ########  ######   ########   ##  ## ## ##
  # ##     ## ##     ##    ##   ##   ##       ##     ##  ##  ##  ####
  # ##     ## ##     ##    ##    ##  ##       ##     ##  ##  ##   ###
  # ########   #######     ##     ## ######## ########  #### ##    ##
  */


  if (! dorebin || rebin==0) return;

  // RooBinning binning(bins_Zprime_rebin.size()-1, &bins_Zprime_rebin[0], "rebin");
  int Nbins = 0;
  double xmin = 1e6;
  double xmax = 0;
  for (int i=0; i<histo_map["bkg_pred"]->GetNbinsX()+1; ++i){
    if (histo_map["bkg_pred"]->GetXaxis()->GetBinLowEdge(i)>=fit_lo){
      xmin = std::min(xmin, histo_map["bkg_pred"]->GetXaxis()->GetBinLowEdge(i));
      if (histo_map["bkg_pred"]->GetXaxis()->GetBinUpEdge(i)<=fit_hi){
        xmax = std::max(xmax, histo_map["bkg_pred"]->GetXaxis()->GetBinUpEdge(i));
        ++Nbins;
      }
    }
  }
  RooBinning binningFitting(Nbins, xmin, xmax, "rebin");
  if (debug) std::cout << "binning fitting " << Nbins << " " << xmin << " " << xmax << binningFitting << '\n';
  x_var->setBinning(binningFitting);

  for (int i=0; i<histo_map["bkg_pred"]->GetNbinsX()+1; ++i){
    if (histo_map["bkg_pred"]->GetXaxis()->GetBinLowEdge(i)>=x_lo){
      xmin = std::min(xmin, histo_map["bkg_pred"]->GetXaxis()->GetBinLowEdge(i));
      if (histo_map["bkg_pred"]->GetXaxis()->GetBinUpEdge(i)<=x_hi){
        xmax = std::max(xmax, histo_map["bkg_pred"]->GetXaxis()->GetBinUpEdge(i));
        ++Nbins;
      }
    }
  }
  RooBinning binning(Nbins, xmin, xmax, "rebin");
  if (debug) std::cout << "binning " << Nbins << " " << xmin << " " << xmax << binning << '\n';
  x_var->setBinning(binning);


}



void CreateRooWorkspace::InitializePDFs() {
  /*
  # #### ##    ## #### ######## ####    ###    ##       #### ######## ########    ########  ########  ########  ######
  #  ##  ###   ##  ##       ##   ##    ## ##   ##        ##       ##  ##          ##     ## ##     ## ##       ##    ##
  #  ##  ####  ##  ##      ##    ##   ##   ##  ##        ##      ##   ##          ##     ## ##     ## ##       ##
  #  ##  ## ## ##  ##     ##     ##  ##     ## ##        ##     ##    ######      ########  ##     ## ######    ######
  #  ##  ##  ####  ##    ##      ##  ######### ##        ##    ##     ##          ##        ##     ## ##             ##
  #  ##  ##   ###  ##   ##       ##  ##     ## ##        ##   ##      ##          ##        ##     ## ##       ##    ##
  # #### ##    ## #### ######## #### ##     ## ######## #### ######## ########    ##        ########  ##        ######
  */

  std::string FitName, ParName;

  for (const int & mass : MassPoints) {
    FitName = GetSgName(mass);
    ParName = "sg"+unique_name+FitName;
    fitPars[FitName].emplace_back(new RooRealVar((ParName+"_p0").c_str(), (ParName+"_p0").c_str(), mass, mass*0.8, mass*1.2));
    fitPars[FitName].emplace_back(new RooRealVar((ParName+"_p1").c_str(), (ParName+"_p1").c_str(), (mass>3000)? 150: 30, (mass>3000)? 120: 10., (mass>3000)? 400:110.));
    fitPars[FitName].emplace_back(new RooRealVar((ParName+"_p2").c_str(), (ParName+"_p2").c_str(), 1, -10, 10));
    fitPars[FitName].emplace_back(new RooRealVar((ParName+"_p3").c_str(), (ParName+"_p3").c_str(), 0.1, -10, 10));
    Fits_map[FitName]["Gauss"].reset(new RooGaussian(FitName.c_str(), "Signal Prediction", *x_var, *fitPars[FitName][0], *fitPars[FitName][1]));
    Fits_map[FitName]["CB"].reset(new RooCBShape(FitName.c_str(), "Signal Prediction", *x_var, *fitPars[FitName][0], *fitPars[FitName][1], *fitPars[FitName][2], *fitPars[FitName][3]));
  }


  for (auto mode: Modes) {
    ParName = mode+unique_name;

    // Old functions
    // FitName = mode+"NO";
    // fitPars[FitName].emplace_back(new RooRealVar((ParName+"_Novo_p1").c_str(), (ParName+"_Novo_p1").c_str(),  500, 200,  1000.));
    // fitPars[FitName].emplace_back(new RooRealVar((ParName+"_Novo_p2").c_str(), (ParName+"_Novo_p2").c_str(),  100,  10,  200.));
    // fitPars[FitName].emplace_back(new RooRealVar((ParName+"_Novo_p3").c_str(), (ParName+"_Novo_p3").c_str(), -0.6,  -2,  2. ));
    // Fits_map[mode]["NO"].reset(new RooNovosibirsk((ParName+"_NO").c_str(),(ParName+"_NO").c_str(),*x_var, *fitPars[FitName][0], *fitPars[FitName][1], *fitPars[FitName][2]));

    //TODO To be put in the old functions
    FitName = mode+"CB";
    fitPars[FitName].emplace_back(new RooRealVar((ParName+"_CB_p1").c_str(), (ParName+"_CB_p1").c_str(), +1,   0,   10));
    fitPars[FitName].emplace_back(new RooRealVar((ParName+"_CB_p2").c_str(), (ParName+"_CB_p2").c_str(), +10,  0,   100));
    fitPars[FitName].emplace_back(new RooRealVar((ParName+"_CB_p3").c_str(), (ParName+"_CB_p3").c_str(), +500, 100, 1000));
    fitPars[FitName].emplace_back(new RooRealVar((ParName+"_CB_p4").c_str(), (ParName+"_CB_p4").c_str(), +100, 10,  500));
    Fits_map[mode]["CB"].reset(new RevCrystalBall((ParName+"_CB").c_str(), (ParName+"_CB").c_str(), *x_var, *fitPars[FitName][0], *fitPars[FitName][1], *fitPars[FitName][2], *fitPars[FitName][3]));
    FitName = mode+"Exp_1";
    fitPars[FitName].emplace_back(new RooRealVar((ParName+"_Exp_1_p1").c_str(), (ParName+"_Exp_1_p1").c_str(), -4., -100, 100));
    Fits_map[mode]["Exp_1"].reset(new PolinomialExponent_1p((ParName+"_Exp_1").c_str(),(ParName+"_Exp_1").c_str(),*x_var, *fitPars[FitName][0]));


    FitName = mode+"Exp_2";
    fitPars[FitName].emplace_back(new RooRealVar((ParName+"_Exp_2_p1").c_str(), (ParName+"_Exp_2_p1").c_str(), -4., -100, 100));
    fitPars[FitName].emplace_back(new RooRealVar((ParName+"_Exp_2_p2").c_str(), (ParName+"_Exp_2_p2").c_str(), 0.4, -100, 100));
    Fits_map[mode]["Exp_2"].reset(new PolinomialExponent_2p((ParName+"_Exp_2").c_str(),(ParName+"_Exp_2").c_str(),*x_var, *fitPars[FitName][0], *fitPars[FitName][1]));

    FitName = mode+"Exp_3";
    fitPars[FitName].emplace_back(new RooRealVar((ParName+"_Exp_3_p1").c_str(), (ParName+"_Exp_3_p1").c_str(), -4., -100, 100));
    fitPars[FitName].emplace_back(new RooRealVar((ParName+"_Exp_3_p2").c_str(), (ParName+"_Exp_3_p2").c_str(), 0.4, -100, 100));
    fitPars[FitName].emplace_back(new RooRealVar((ParName+"_Exp_3_p3").c_str(), (ParName+"_Exp_3_p3").c_str(), -0.1, -100, 100));
    Fits_map[mode]["Exp_3"].reset(new PolinomialExponent_3p((ParName+"_Exp_3").c_str(),(ParName+"_Exp_3").c_str(),*x_var, *fitPars[FitName][0], *fitPars[FitName][1], *fitPars[FitName][2]));

    FitName = mode+"Exp_4";
    fitPars[FitName].emplace_back(new RooRealVar((ParName+"_Exp_4_p1").c_str(), (ParName+"_Exp_4_p1").c_str(), 4., -100, 100));
    fitPars[FitName].emplace_back(new RooRealVar((ParName+"_Exp_4_p2").c_str(), (ParName+"_Exp_4_p2").c_str(), -10, -100,  10));
    fitPars[FitName].emplace_back(new RooRealVar((ParName+"_Exp_4_p3").c_str(), (ParName+"_Exp_4_p3").c_str(), 4, -100, 100));
    fitPars[FitName].emplace_back(new RooRealVar((ParName+"_Exp_4_p4").c_str(), (ParName+"_Exp_4_p4").c_str(), -0.4, -100, 100));
    Fits_map[mode]["Exp_4"].reset(new PolinomialExponent_4p((ParName+"_Exp_4").c_str(),(ParName+"_Exp_4").c_str(),*x_var, *fitPars[FitName][0], *fitPars[FitName][1], *fitPars[FitName][2], *fitPars[FitName][3]));

    FitName = mode+"Exp_5";
    fitPars[FitName].emplace_back(new RooRealVar((ParName+"_Exp_5_p1").c_str(), (ParName+"_Exp_5_p1").c_str(), -8.46, -9.00, -8.00));
    fitPars[FitName].emplace_back(new RooRealVar((ParName+"_Exp_5_p2").c_str(), (ParName+"_Exp_5_p2").c_str(),  7.75,  7.00,  9.00));
    fitPars[FitName].emplace_back(new RooRealVar((ParName+"_Exp_5_p3").c_str(), (ParName+"_Exp_5_p3").c_str(), -7.45, -8.00, -6.00));
    fitPars[FitName].emplace_back(new RooRealVar((ParName+"_Exp_5_p4").c_str(), (ParName+"_Exp_5_p4").c_str(),  3.35,  3.00,  4.00));
    fitPars[FitName].emplace_back(new RooRealVar((ParName+"_Exp_5_p5").c_str(), (ParName+"_Exp_5_p5").c_str(), -5.37, -6.00, -4.00));
    Fits_map[mode]["Exp_5"].reset(new PolinomialExponent_5p((ParName+"_Exp_5").c_str(),(ParName+"_Exp_5").c_str(),*x_var, *fitPars[FitName][0], *fitPars[FitName][1], *fitPars[FitName][2], *fitPars[FitName][3], *fitPars[FitName][4]));

    FitName = mode+"Exp_6";
    fitPars[FitName].emplace_back(new RooRealVar((ParName+"_Exp_6_p1").c_str(), (ParName+"_Exp_6_p1").c_str(), -3.68, -4.0, -3.0 ));
    fitPars[FitName].emplace_back(new RooRealVar((ParName+"_Exp_6_p2").c_str(), (ParName+"_Exp_6_p2").c_str(), -1.57, -2.0, -1.0 ));
    fitPars[FitName].emplace_back(new RooRealVar((ParName+"_Exp_6_p3").c_str(), (ParName+"_Exp_6_p3").c_str(), +1.56, +1.0, +2.0 ));
    fitPars[FitName].emplace_back(new RooRealVar((ParName+"_Exp_6_p4").c_str(), (ParName+"_Exp_6_p4").c_str(), -1.14, -1.5, -1.0 ));
    fitPars[FitName].emplace_back(new RooRealVar((ParName+"_Exp_6_p5").c_str(), (ParName+"_Exp_6_p5").c_str(), +5.29, +5.0, +5.5 ));
    fitPars[FitName].emplace_back(new RooRealVar((ParName+"_Exp_6_p6").c_str(), (ParName+"_Exp_6_p6").c_str(), -7.05, -9.5, -9.0 ));
    Fits_map[mode]["Exp_6"].reset(new PolinomialExponent_6p((ParName+"_Exp_6").c_str(),(ParName+"_Exp_6").c_str(),*x_var, *fitPars[FitName][0], *fitPars[FitName][1], *fitPars[FitName][2], *fitPars[FitName][3], *fitPars[FitName][4], *fitPars[FitName][5]));

    // Not fitting functions
    // FitName = mode+"BkgPdf4p";
    // fitPars[FitName].emplace_back(new RooRealVar((ParName+"_BkgPdf4p_p1").c_str(), (ParName+"_BkgPdf4p_p1").c_str(), +2.0, +2.0, +4.0));
    // fitPars[FitName].emplace_back(new RooRealVar((ParName+"_BkgPdf4p_p2").c_str(), (ParName+"_BkgPdf4p_p2").c_str(), +4.2, +3.0, +5.0));
    // fitPars[FitName].emplace_back(new RooRealVar((ParName+"_BkgPdf4p_p3").c_str(), (ParName+"_BkgPdf4p_p3").c_str(), +1.4, +0.5, +2.0));
    // fitPars[FitName].emplace_back(new RooRealVar((ParName+"_BkgPdf4p_p4").c_str(), (ParName+"_BkgPdf4p_p4").c_str(), -1.0, -2.0, -0.5));
    // Fits_map[mode]["BkgPdf4p"].reset(new BkgPdf4p((ParName+"_BkgPdf4p").c_str(),(ParName+"_BkgPdf4p").c_str(),*x_var, *fitPars[FitName][0], *fitPars[FitName][1], *fitPars[FitName][2], *fitPars[FitName][3]));
    //
    // FitName = mode+"BkgPdf3p";
    // fitPars[FitName].emplace_back(new RooRealVar((ParName+"_BkgPdf3p_p1").c_str(), (ParName+"_BkgPdf3p_p1").c_str(), -10., -11., -9.0));
    // fitPars[FitName].emplace_back(new RooRealVar((ParName+"_BkgPdf3p_p2").c_str(), (ParName+"_BkgPdf3p_p2").c_str(), +7.0, +6.0, +8.0));
    // fitPars[FitName].emplace_back(new RooRealVar((ParName+"_BkgPdf3p_p3").c_str(), (ParName+"_BkgPdf3p_p3").c_str(), +3.7, +3.0, +4.0));
    //
    // fitPars[FitName].emplace_back(new RooRealVar((ParName+"_BkgPdf3p_p1").c_str(), (ParName+"_BkgPdf3p_p1").c_str(), -13.2, -1000, 1000));
    // fitPars[FitName].emplace_back(new RooRealVar((ParName+"_BkgPdf3p_p2").c_str(), (ParName+"_BkgPdf3p_p2").c_str(), +9.1, -1000, 1000));
    // fitPars[FitName].emplace_back(new RooRealVar((ParName+"_BkgPdf3p_p3").c_str(), (ParName+"_BkgPdf3p_p3").c_str(), +2.5, -100, 100));
    //
    // fitPars[FitName].emplace_back(new RooRealVar((ParName+"_BkgPdf3p_p1").c_str(), (ParName+"_BkgPdf3p_p1").c_str(), -4.3, -5.3, -3.3));
    // fitPars[FitName].emplace_back(new RooRealVar((ParName+"_BkgPdf3p_p2").c_str(), (ParName+"_BkgPdf3p_p2").c_str(), +5.7, +4.7, +7.7));
    // fitPars[FitName].emplace_back(new RooRealVar((ParName+"_BkgPdf3p_p3").c_str(), (ParName+"_BkgPdf3p_p3").c_str(), +2.5, +1.5, +3.5));
    // Fits_map[mode]["BkgPdf3p"].reset(new BkgPdf3p((ParName+"_BkgPdf3p").c_str(),(ParName+"_BkgPdf3p").c_str(),*x_var, *fitPars[FitName][0], *fitPars[FitName][1], *fitPars[FitName][2]));


  }
}


void CreateRooWorkspace::CreateRooDataHist() {
  /*
  #  ######  ########  ########    ###    ######## ######## ########   #######   #######  ########     ###    ########    ###    ##     ## ####  ######  ########
  # ##    ## ##     ## ##         ## ##      ##    ##       ##     ## ##     ## ##     ## ##     ##   ## ##      ##      ## ##   ##     ##  ##  ##    ##    ##
  # ##       ##     ## ##        ##   ##     ##    ##       ##     ## ##     ## ##     ## ##     ##  ##   ##     ##     ##   ##  ##     ##  ##  ##          ##
  # ##       ########  ######   ##     ##    ##    ######   ########  ##     ## ##     ## ##     ## ##     ##    ##    ##     ## #########  ##   ######     ##
  # ##       ##   ##   ##       #########    ##    ##       ##   ##   ##     ## ##     ## ##     ## #########    ##    ######### ##     ##  ##        ##    ##
  # ##    ## ##    ##  ##       ##     ##    ##    ##       ##    ##  ##     ## ##     ## ##     ## ##     ##    ##    ##     ## ##     ##  ##  ##    ##    ##
  #  ######  ##     ## ######## ##     ##    ##    ######## ##     ##  #######   #######  ########  ##     ##    ##    ##     ## ##     ## ####  ######     ##
  */


  if (dorebin && rebin!=0) {
    for (auto mode: Modes) rooHist_map[mode].reset(new RooDataHist(mode.c_str(), mode.c_str(), RooArgList(*x_var), RooFit::Import(*histo_map[mode].get(),doBinWidth)));
    for (const int & mass : MassPoints) rooHist_map[GetSgName(mass)].reset(new RooDataHist(GetSgName(mass).c_str(), GetSgName(mass).c_str(), RooArgList(*x_var), RooFit::Import(*histo_map[GetSgName(mass)].get(),doBinWidth)));
  } else {
    for (auto mode: Modes) rooHist_map[mode].reset(new RooDataHist(mode.c_str(), mode.c_str(), RooArgList(*x_var), histo_map[mode].get()));
    for (const int & mass : MassPoints) rooHist_map[GetSgName(mass)].reset(new RooDataHist(GetSgName(mass).c_str(), GetSgName(mass).c_str(), RooArgList(*x_var), histo_map[GetSgName(mass)].get()));
  }

}

void CreateRooWorkspace::DoFits() {
  /*
  # ########   #######  ######## #### ########  ######
  # ##     ## ##     ## ##        ##     ##    ##    ##
  # ##     ## ##     ## ##        ##     ##    ##
  # ##     ## ##     ## ######    ##     ##     ######
  # ##     ## ##     ## ##        ##     ##          ##
  # ##     ## ##     ## ##        ##     ##    ##    ##
  # ########   #######  ##       ####    ##     ######
  */

  for (const auto& x : histo_map) {
    fit_min[x.first] = GetRange(x.second.get(),fit_lo);
    fit_max[x.first] = GetRange(x.second.get(),fit_hi);
    if ("DY_SR"==x.first) fit_min[x.first] = GetRange(x.second.get(),800);
    std::string mass = x.first;
    if (mass.compare(0,1,"M")==0) {
      CalculateSignalFittingRange(std::stoi(mass.substr(1, mass.size()-1)), fit_min[x.first], fit_max[x.first], plot_min[x.first], plot_max[x.first], y_max[x.first]); // TODO
      // nEventsSignal[x.first] = CalculateIntegral(histo_map[x.first].get(),x_lo,x_hi,doBinWidth); //TODO what range?
      nEventsSignal[x.first] = CalculateIntegral(histo_map[x.first].get(),fit_min[x.first]-100,fit_max[x.first]+100,doBinWidth);// TODO!!!!
    }
  }


  std::cout << "********************" << '\n';
  std::cout << "*  Background Fits *" << '\n';
  std::cout << "********************" << '\n';

  if (debug) {
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "Fitting Modes:\t\t";
    for (auto x : doFits_map) std::cout << x.first << '\t';
    std::cout << "\nFitting Functions:\t";
    for (auto x : doFits_map["bkg_pred"]) if (x.second) std::cout << x.first << '\t';
    std::cout << '\n';
    std::cout << "Fitting range:\t fit_min\fit_max" << '\n';
    // std::cout << "sub_fit" << std::string(10,' ') << fit_min_turnOn << "\t" << fit_max_tails << std::endl; // not used
    for (auto mode : Modes) std::cout << mode << std::string(17-mode.size(),' ') << fit_min[mode] << "\t" << fit_max[mode] << std::endl;
    std::cout << "----------------------------------------" << std::endl;
  }

  for (auto mode: Modes) {
    if (mode!="bkg_pred" && mode!="DY_CR" && mode!="DY_SR") continue;
    for (auto const& [model,dofit] : doFits_map[mode] ) {
      if (!dofit) continue;
      std::cout << "PDF FIT " << mode << "\t" << model << '\n';
      if (mode=="bkg_pred" || mode=="data") { //TODO check this
        FitRes_map[mode][model].reset(Fits_map[mode][model]->fitTo(*rooHist_map[mode], RooFit::Range(fit_min[mode], fit_max[mode]), RooFit::SumW2Error(kFALSE), RooFit::Minimizer("Minuit2"), RooFit::Save(), RooFit::Verbose(kFALSE), RooFit::PrintEvalErrors(-1)));
      } else {
        FitRes_map[mode][model].reset(Fits_map[mode][model]->fitTo(*rooHist_map[mode], RooFit::Range(fit_min[mode], fit_max[mode]), RooFit::SumW2Error(kTRUE), RooFit::Save(), RooFit::Verbose(kFALSE), RooFit::PrintEvalErrors(-1)));
      }
    }
  }

  std::cout << "********************" << '\n';
  std::cout << "*    Signal Fits   *" << '\n';
  std::cout << "********************" << '\n';

  for (const int & mass : MassPoints) {
    if (debug) std::cout << "fit mass:" << mass << "\t" << fit_min[GetSgName(mass)] << " -- " << fit_max[GetSgName(mass)] << '\n';
    FitRes_map[GetSgName(mass)][FitSignal].reset(Fits_map[GetSgName(mass)][FitSignal]->fitTo(*rooHist_map[GetSgName(mass)], RooFit::Range(fit_min[GetSgName(mass)], fit_max[GetSgName(mass)]), RooFit::SumW2Error(kTRUE), RooFit::Save(), RooFit::Verbose(kFALSE), RooFit::PrintEvalErrors(-1)));
  }

}


void CreateRooWorkspace::ImportToWorkspace() {
  /*
  # #### ##     ## ########   #######  ########  ######## ########  #######  ##      ##  #######  ########  ##    ##  ######  ########     ###     ######  ########
  #  ##  ###   ### ##     ## ##     ## ##     ##    ##       ##    ##     ## ##  ##  ## ##     ## ##     ## ##   ##  ##    ## ##     ##   ## ##   ##    ## ##
  #  ##  #### #### ##     ## ##     ## ##     ##    ##       ##    ##     ## ##  ##  ## ##     ## ##     ## ##  ##   ##       ##     ##  ##   ##  ##       ##
  #  ##  ## ### ## ########  ##     ## ########     ##       ##    ##     ## ##  ##  ## ##     ## ########  #####     ######  ########  ##     ## ##       ######
  #  ##  ##     ## ##        ##     ## ##   ##      ##       ##    ##     ## ##  ##  ## ##     ## ##   ##   ##  ##         ## ##        ######### ##       ##
  #  ##  ##     ## ##        ##     ## ##    ##     ##       ##    ##     ## ##  ##  ## ##     ## ##    ##  ##   ##  ##    ## ##        ##     ## ##    ## ##
  # #### ##     ## ##         #######  ##     ##    ##       ##     #######   ###  ###   #######  ##     ## ##    ##  ######  ##        ##     ##  ######  ########
  */

  for (const int & mass : MassPoints) ws->import(*Fits_map[GetSgName(mass)][FitSignal], RooFit::Silence());
  ws->import(*data_obs.get());

  //TODO decide weather to go for shaped analysis or not
  // ws->import(*(new RooHistPdf(SgName.c_str(),"",*x_var,*rooHist_map[GetSgName(mass)])));
  // for (auto [mode,h]: rooHist_map) {
  //   // Store RooDataHist format for bkg_pred and data et all
  //   ws->import(*h);
  //   // Store TH1F format for bkg_pred and data et all
  //   ws->import(*histo_map[mode],Form("h_%s",+h->GetName()));
  //   // Store RooHistPdf format for bkg_pred and data et all
  //   ws->import(*(new RooHistPdf( Form("hp_%s",+h->GetName()),"",RooArgSet(*x_var),*h)));
  // }

  if (debug) ws->Print();
}


void CreateRooWorkspace::DoPlots() {
  /*
  # ########   #######  ########  ##        #######  ########  ######
  # ##     ## ##     ## ##     ## ##       ##     ##    ##    ##    ##
  # ##     ## ##     ## ##     ## ##       ##     ##    ##    ##
  # ##     ## ##     ## ########  ##       ##     ##    ##     ######
  # ##     ## ##     ## ##        ##       ##     ##    ##          ##
  # ##     ## ##     ## ##        ##       ##     ##    ##    ##    ##
  # ########   #######  ##        ########  #######     ##     ######
  */

  PlotBkgFit();
  PlotSignals();
  PlotSgPars();
}

void CreateRooWorkspace::PlotSignals() {
  // c_sg_all.reset(tdrDiCanvas("signal", plot_lo, x_hi, 1*1e-03, 1, -6, 6, nameXaxis, nameYaxis, "Pull"));
  //TODO
  RooPlot* Pullplotter_all = x_var->frame(RooFit::Range(x_lo-10,x_hi+10), RooFit::Name("sgALL_pull"));
  std::unique_ptr<RooPlot> plotter_all(x_var->frame(RooFit::Range(x_lo,x_hi), RooFit::Name("signal_all")));
  std::unique_ptr<TCanvas> c_sg_all(tdrDiCanvas("signal", plot_lo, x_hi, 1*1e-03, 1, -6, 6, nameXaxis, nameYaxis, "Pull"));

  int i_mass = -1;
  for (const int & mass : MassPoints) {
    i_mass++;
    std::string SgName = GetSgName(mass);

    // std::cout << "FIND pre " << Pullplotter.get() << '\n';
    // Pullplotter.reset(x_var->frame(RooFit::Range(plot_min[SgName],plot_max[SgName]),RooFit::Name(SgName.c_str())));
    // std::cout << "FIND post " << Pullplotter.get() << '\n'; //TODO

    RooPlot* Pullplotter = x_var->frame(RooFit::Range(plot_min[SgName],plot_max[SgName]),RooFit::Name(SgName.c_str()));

    if (debug) std::cout << "plotting mass:" << mass << "\t" << plot_min[SgName] << " -- " << plot_max[SgName] << "\tymax " << y_max[SgName] << '\n';

    plotter.reset(x_var->frame(plot_min[SgName],plot_max[SgName]));
    TCanvas* c_sg = tdrDiCanvas(SgName.c_str(), plot_min[SgName],plot_max[SgName], 2*1e-03, y_max[SgName], -6, 6, nameXaxis, nameYaxis, "Pull");
    // c_sg->cd(1)->SetLogy(1);
    // rooHist_map[SgName]->plotOn(plotter.get());
    rooHist_map[SgName]->plotOn(plotter.get(), RooFit::LineColor(kBlack), RooFit::MarkerColor(kBlack));
    Fits_map[SgName][FitSignal]->plotOn(plotter.get(), RooFit::VisualizeError(*FitRes_map[SgName][FitSignal], 1), RooFit::FillColor(kRed-7), RooFit::FillStyle(3001));
    Fits_map[SgName][FitSignal]->plotOn(plotter.get(), RooFit::VisualizeError(*FitRes_map[SgName][FitSignal], 2), RooFit::FillColor(kRed-9), RooFit::FillStyle(3001));
    Fits_map[SgName][FitSignal]->plotOn(plotter.get(), RooFit::LineColor(kRed));
    rooHist_map[SgName]->plotOn(plotter.get(), RooFit::LineColor(kBlack), RooFit::MarkerColor(kBlack));
    RooHist* hpull = plotter->pullHist();
    int npf = FitRes_map[SgName][FitSignal]->floatParsFinal().getSize();
    // int ndf = hpull->GetN()-npf;
    int ndf;
    double chi2;
    CalculateChiSquare(chi2, ndf, hpull, fit_min[SgName],fit_max[SgName]); //TODO do the same for bkg plotter->chiSquare(npf); fit_min or plot_min?
    ndf -=npf;
    double pv = TMath::Prob(chi2, ndf)*100;
    chi2 /=ndf;
    std::cout << "FIND ME " << mass << " " << year << " " << channel << " " << chi2 << " " << chi2*ndf << " " << ndf << '\n';
    SignalProperties << SgName+" signal number of events = " << nEventsSignal[SgName] <<""<<std::endl;//TODO
    SignalProperties << " chi2-ndf-pvalue "  << chi2 << " " << ndf << " " << pv <<std::endl;
    SignalProperties << fitPars[SgName][0]->GetName() << "   param   " << fitPars[SgName][0]->getVal() << " " << fitPars[SgName][0]->getError() << std::endl;
    SignalProperties << fitPars[SgName][1]->GetName() << "   param   " << fitPars[SgName][1]->getVal() << " " << fitPars[SgName][1]->getError() << std::endl;
    SignalProperties << fitPars[SgName][2]->GetName() << "   param   " << fitPars[SgName][2]->getVal() << " " << fitPars[SgName][2]->getError() << std::endl;
    SignalProperties << fitPars[SgName][3]->GetName() << "   param   " << fitPars[SgName][3]->getVal() << " " << fitPars[SgName][3]->getError() << std::endl;

    std::unique_ptr<TPaveText> pave(new TPaveText(0.7,0.6,0.8,0.8,"NDC"));
    pave->SetBorderSize(0); pave->SetTextSize(0.03); pave->SetLineColor(1); pave->SetLineStyle(1);
    pave->SetLineWidth(2); pave->SetFillColor(0); pave->SetFillStyle(0);
    pave->AddText(TString::Format("Fit range = [%.0f,%.0f]", fit_min[SgName],fit_max[SgName]));
    pave->AddText(TString::Format("#chi^{2}/n.d.f. = %.1f", chi2));
    pave->AddText(TString::Format("p-value = %.1f", pv));
    pave->AddText(TString::Format("#mu = %2.3f +- %2.3f", fitPars[SgName][0]->getVal(),fitPars[SgName][0]->getError()));
    pave->AddText(TString::Format("#sigma = %2.3f +- %2.3f", fitPars[SgName][1]->getVal(),fitPars[SgName][1]->getError()));
    pave->AddText(TString::Format("#alpha = %2.3f +- %2.3f", fitPars[SgName][2]->getVal(),fitPars[SgName][2]->getError()));
    pave->AddText(TString::Format("k = %2.3f +- %2.3f", fitPars[SgName][3]->getVal(),fitPars[SgName][3]->getError()));
    SgPars["Masses"].at(i_mass)       = (double)mass;
    // SgPars["nevents"].at(i_mass)      = CalculateIntegral(histo_map[SgName],rangeLo-100,rangeHi+100,doBinWidth);//TODO
    SgPars["nevents"].at(i_mass)      = CalculateIntegral(histo_map[SgName].get(),x_lo,x_hi,doBinWidth);//TODO what range?
    SgPars["nevents_err"].at(i_mass)  = TMath::Sqrt(SgPars["nevents"].at(i_mass));
    SgPars["mean"].at(i_mass)         = fitPars[SgName][0]->getVal();
    SgPars["mean_err"].at(i_mass)     = fitPars[SgName][0]->getError();
    SgPars["sigma"].at(i_mass)        = fitPars[SgName][1]->getVal();
    SgPars["sigma_err"].at(i_mass)    = fitPars[SgName][1]->getError();
    SgPars["alpha"].at(i_mass)        = fitPars[SgName][2]->getVal();
    SgPars["alpha_err"].at(i_mass)    = fitPars[SgName][2]->getError();
    SgPars["k"].at(i_mass)            = fitPars[SgName][3]->getVal();
    SgPars["k_err"].at(i_mass)        = fitPars[SgName][3]->getError();
    SgPars["chi2"].at(i_mass)         = chi2;
    SgPars["pvalue"].at(i_mass)       = pv;
    SgPars["fit_min"].at(i_mass)      = fit_min[SgName];
    SgPars["fit_max"].at(i_mass)      = fit_max[SgName];

    Fits_map[SgName][FitSignal]->plotOn(plotter_all.get(), RooFit::LineColor(kRed));
    c_sg->cd(1);
    plotter->Draw("same");
    histo_map[SgName]->Draw("same");
    pave->Draw("same");
    c_sg->cd(2);
    Pullplotter->addPlotable(hpull,"P same");
    Pullplotter->SetNdivisions(505,"Y");
    Pullplotter->Draw("same");

    std::unique_ptr<TLine> line(new TLine(plot_min[SgName], 0, plot_max[SgName], 0));
    line->SetLineWidth(2);
    line->Draw("same");

    c_sg_all->cd(1);
    plotter_all->Draw("same");
    c_sg_all->cd(2);
    Pullplotter_all->addPlotable(hpull,"P same");
    Pullplotter_all->SetNdivisions(505,"Y");
    Pullplotter_all->Draw("same");
    line->Draw("same");

    c_sg->SaveAs(workingDir+"Fit_Sg"+SgName+"_"+histFolder+extra_text+"."+plotting_mode);
    if (plotting_mode_2!="") c_sg->SaveAs(workingDir+"Fit_Sg"+SgName+"_"+histFolder+extra_text+"."+plotting_mode_2);

  }

  c_sg_all->SaveAs((workingDir+"Fit_Sg_all_"+histFolder+"."+plotting_mode).c_str());
  if (plotting_mode_2!="") c_sg_all->SaveAs((workingDir+"Fit_Sg_all_"+histFolder+"."+plotting_mode_2).c_str());

}

void CreateRooWorkspace::PlotSgPars() {

  gStyle->SetOptFit(kFALSE);
  std::unique_ptr<TCanvas> c_sg_par;
  std::unique_ptr<TGraphErrors> gr_par;
  std::unique_ptr<TF1> fit_par;
  std::unique_ptr<TSpline3> fit_par_spline;
  std::unique_ptr<TPaveText> pave_par;
  std::vector<double> dummy(SgPars["Masses"].size(),50);
  std::vector<double> dummy2(SgPars["Masses"].size(),0);

  for (auto x: SgPars) {
    if (x.first=="Masses" || x.first.find("_err") != std::string::npos) continue;
    double y_max = (x.first=="mean")? 9000: ((x.first=="sigma")? 500:((x.first=="pvalue")? 100:((x.first=="nevents")? 40*lumi_map.at(year).at("lumi_fb")/lumi_map.at("RunII").at("lumi_fb"):5)));
    y_max = (x.first=="fit_min")? 9000: y_max;
    y_max = (x.first=="fit_max")? 9000: y_max;
    c_sg_par.reset(tdrCanvas(("c_sg_"+x.first).c_str(), plot_lo, x_hi, 2*1e-03, y_max, nameXaxis, (x.first).c_str()));
    if (x.first=="chi2" || x.first=="pvalue" || x.first=="fit_min" || x.first=="fit_max") gr_par.reset(new TGraphErrors(SgPars["Masses"].size(), &(SgPars["Masses"][0]), &(x.second[0]), &(dummy2[0]), &(dummy2[0])));
    else gr_par.reset(new TGraphErrors(SgPars["Masses"].size(), &(SgPars["Masses"][0]), &(x.second[0]), &(dummy[0]), &(SgPars[x.first+"_err"][0])));
    if (x.first!="chi2" && x.first!="pvalue") {
      fit_par.reset(new TF1((x.first).c_str(),(x.first=="nevents")? "pol3":((x.first=="sigma")?"pol2":"pol1"),500, 8100));
      fit_par_spline.reset(new TSpline3(("spline_"+x.first).c_str(), &(SgPars["Masses"][0]), &(x.second[0]), SgPars["Masses"].size(), "b2e2"));
      fit_par_spline->SetLineColor(kOrange+1);
      fit_par_spline->Draw("same");
      // if (x.first=="nevents") fit_par->SetParameters();
      fit_par->SetLineColor(kRed+1); fit_par->SetLineWidth(3);
      gr_par->Fit(fit_par.get(),"RSQ");
      pave_par.reset(new TPaveText(0.7,0.6,0.8,0.8,"NDC"));
      pave_par->SetBorderSize(0); pave_par->SetTextSize(0.03); pave_par->SetLineColor(1); pave_par->SetLineStyle(1);
      pave_par->SetLineWidth(2); pave_par->SetFillColor(0); pave_par->SetFillStyle(0);
      for (int i = 0; i < fit_par->GetNumberFreeParameters(); i++) {
        if (fabs(fit_par->GetParameter(i)<1e-03)) pave_par->AddText(TString::Format("p%d = (%2.3f +- %2.3f)*10^{-3}", i,fit_par->GetParameter(i)*1e03,fit_par->GetParError(i)*1e03));
        else pave_par->AddText(TString::Format("p%d = %2.3f +- %2.3f", i,fit_par->GetParameter(i),fit_par->GetParError(i)));
        if (x.first=="fit_min" || x.first=="fit_max") {
          std::cout << x.first << " " << TString::Format("p%d = %2.3f +- %2.3f", i,fit_par->GetParameter(i),fit_par->GetParError(i)) << '\n';
        }
      }
      pave_par->AddText(TString::Format("#chi^{2}/n.d.f. = %.1f", fit_par->GetChisquare()/fit_par->GetNDF()));
      pave_par->AddText(TString::Format("p-value = %.2f", fit_par->GetProb()));
      gStyle->SetOptFit(kFALSE);
      pave_par->Draw("same");
    }
    tdrDraw(gr_par.get(), "P", kFullDotLarge, kRed+1, kSolid, kRed+1, 3004, kRed+1);
    c_sg_par->SaveAs((workingDir+"Fit_Sg_"+x.first+"_"+histFolder+"."+plotting_mode).c_str());
    if (plotting_mode_2!="") c_sg_par->SaveAs((workingDir+"Fit_Sg_"+x.first+"_"+histFolder+"."+plotting_mode_2).c_str());
  }
}

void CreateRooWorkspace::PlotControl() {
  /*
  #  ######  ##     ## ########  ######  ##    ##    ########  ##        #######  ########  ######
  # ##    ## ##     ## ##       ##    ## ##   ##     ##     ## ##       ##     ##    ##    ##    ##
  # ##       ##     ## ##       ##       ##  ##      ##     ## ##       ##     ##    ##    ##
  # ##       ######### ######   ##       #####       ########  ##       ##     ##    ##     ######
  # ##       ##     ## ##       ##       ##  ##      ##        ##       ##     ##    ##          ##
  # ##    ## ##     ## ##       ##    ## ##   ##     ##        ##       ##     ##    ##    ##    ##
  #  ######  ##     ## ########  ######  ##    ##    ##        ########  #######     ##     ######
  */

  TCanvas* c_data_obs = tdrCanvas("data_obs", plot_lo, plot_hi, plot_ylo, 1e03, nameXaxis, nameYaxis);
  c_data_obs->SetLogy(1);
  // plotter = x_var->frame(plot_lo,plot_hi);
  plotter.reset(x_var->frame(plot_lo,plot_hi));
  data_obs->plotOn(plotter.get(), RooFit::LineColor(kBlack), RooFit::MarkerColor(kBlack));
  plotter->Draw("same");
  std::unique_ptr<TLegend> leg_data_obs(tdrLeg(0.50,0.70,0.9,0.9, 0.035));
  leg_data_obs->AddEntry(data_obs.get(),  Form("%s: %s", data_obs->GetName(), data_obs->GetTitle()) ,"lp");
  leg_data_obs->Draw("same");

  c_data_obs->SaveAs((workingDir+"Plot_data_obs_"+histFolder+"."+plotting_mode).c_str());
  if (plotting_mode_2!="") c_data_obs->SaveAs((workingDir+"Plot_data_obs_"+histFolder+"."+plotting_mode_2).c_str());

  // TODO Do Transfer functions
  std::unique_ptr<TCanvas> c_inputs(tdrCanvas("inputs DataCard", plot_lo, plot_hi, plot_ylo, plot_yhi, nameXaxis, nameYaxis));
  c_inputs->SetLogy(1);
  // histo_map["bkg_pred"]->SetLineColor(kBlue+1);
  // histo_map["data"]->SetLineColor(kRed+1);
  // histo_map["main_bkg_SR"]->SetLineColor(kGreen+1);
  histo_map["bkg_pred"]->Draw("same");
  histo_map["data"]->Draw("same");
  histo_map["main_bkg_SR"]->Draw("same");
  std::unique_ptr<TLegend> leg_inputs(tdrLeg(0.50,0.70,0.78,0.9, 0.03));
  leg_inputs->AddEntry(histo_map["bkg_pred"].get(), Form("dofit: %s", histo_map["bkg_pred"]->GetName()) ,"lp");
  leg_inputs->AddEntry(histo_map["data"].get(), Form("extract obs. Limits: %s", histo_map["data"]->GetName()) ,"lp");
  leg_inputs->AddEntry(histo_map["main_bkg_SR"].get(), Form("extract exp. Limits: %s", histo_map["main_bkg_SR"]->GetName()) ,"lp");
  leg_inputs->Draw("same");
  c_inputs->SaveAs((workingDir+"Plot_Inputs_"+histFolder+"."+plotting_mode).c_str());
  if (plotting_mode_2!="") c_inputs->SaveAs((workingDir+"Plot_Inputs_"+histFolder+"."+plotting_mode_2).c_str());


  TCanvas* c_CRvsSR = tdrCanvas(("bkg_pred vs data"+unique_name_complete).c_str(), plot_lo, plot_hi, plot_ylo, plot_yhi, nameXaxis, nameYaxis);
  c_CRvsSR->SetLogy(1);
  std::unique_ptr<TLegend> leg_CRvsSR(tdrLeg(0.50,0.70,0.78,0.9, 0.03));
  for (auto mode : Modes) {
    int color = Colors[mode];
    histo_map[mode]->SetMarkerSize(0.5);
    tdrDraw(histo_map[mode].get(), "", kFullDotLarge, color, kSolid, color, 3004, color);
    leg_CRvsSR->AddEntry(histo_map[mode].get(), Form("%s", histo_map[mode]->GetName()) ,"lp");
  }
  // histo_map["norm"]->SetMarkerSize(0.5);
  // tdrDraw(histo_map["norm"], "", kFullDotLarge, kBlack, kSolid, kBlack, 3004, kBlack);
  // leg_CRvsSR->AddEntry(histo_map["norm"], Form("%s", histo_map["norm"]->GetName()) ,"lp");
  leg_CRvsSR->Draw("same");
  c_CRvsSR->SaveAs((workingDir+"Plot_CRvsSR_"+histFolder+"."+plotting_mode).c_str());
  if (plotting_mode_2!="") c_CRvsSR->SaveAs((workingDir+"Plot_CRvsSR_"+histFolder+"."+plotting_mode_2).c_str());


  TCanvas* c_ratios = tdrDiCanvas(("ratios SR/CR"+unique_name_complete).c_str(), plot_lo, plot_hi, plot_ylo, 1e3,0,2, nameXaxis, nameYaxis, "MC/pol0");
  c_ratios->cd(1)->SetLogy(1);
  gStyle->SetOptFit(kFALSE);
  std::unique_ptr<TLegend> leg_ratios(tdrLeg(0.50,0.70,0.78,0.9, 0.03));
  for (const auto& x: histo_map) {
    std::string hname = x.first;
    if (hname.find("SR")==std::string::npos) continue;
    std::string name = TString(hname).ReplaceAll("_SR","").Data();
    if (name=="WW") continue;
    if (name=="WZ") continue;
    if (name=="ZZ") continue;
    int color = Colors[name];
    std::unique_ptr<TGraphAsymmErrors> ratio(new TGraphAsymmErrors());
    std::unique_ptr<TGraphAsymmErrors> ratio2(new TGraphAsymmErrors());
    std::unique_ptr<TH1F> h_SR((TH1F*)histo_map[x.first]->Clone());
    std::unique_ptr<TH1F> h_CR((TH1F*)histo_map[name+"_CR"]->Clone());
    h_SR->Rebin(5);
    h_CR->Rebin(5);
    ratio->Divide(h_SR.get(), h_CR.get(), "pois");
    ratio2->Divide(h_SR.get(), h_CR.get(), "pois");
    // TCanvas* c_ = tdrCanvas(hname, plot_lo, plot_hi, plot_ylo, 1e2, nameXaxis, nameYaxis);
    // c_->SetLogy(1);
    std::unique_ptr<TF1> f1(new TF1("f1","pol0",fit_min[hname],fit_max[hname]));
    // TF1 *f2 = new TF1("f2","[0]+TMath::Log10([1]+[2]*x)",500,plot_hi);
    // TF1 *f4 = new TF1("f3","[0]-TMath::Exp([1]+[2]*x)",500,plot_hi);
    // TF1 *f4 = new TF1("f4","[0]-TMath::Exp([1]+TMath::Exp([2]+[3]*x))",500,plot_hi);
    f1->SetParameter(0,0.01);
    f1->SetLineColor(color);
    f1->SetLineWidth(3);
    ratio->Fit(f1.get(),"RS0Q");
    // (TVirtualFitter::GetFitter())->GetConfidenceIntervals(f1,0.68);
    if (BkgName == name) {
      DataCard << name+" rateParam " <<  f1->GetParameter(0) << std::endl;
      // Next number should be ~ 1
      // DataCard << name+" rateParam " <<  f1->GetParameter(0)*CalculateIntegral(histo_map[name+"_CR"],fit_lo,fit_hi,doBinWidth)/CalculateIntegral(histo_map[name+"_SR"],fit_lo,fit_hi,doBinWidth) << std::endl;
    }
    c_ratios->cd(1);
    tdrDraw(ratio.get(), "P0", kFullDotLarge, color, kSolid, color, 3004, color);
    // f4->Draw("same");
    leg_ratios->AddEntry(ratio.get(), name+std::string(name=="TTbar"?8:(name.size()==2?13:3),' ')+TString::Format("#chi^{2}/n.d.f.=%.2f p-value=%.2f", f1->GetChisquare()/f1->GetNDF(), f1->GetProb()) ,"lp");

    c_ratios->cd(2);
    for (int i = 0; i < ratio->GetN(); i++) {
      double x,y;
      // double err;
      ratio->GetPoint(i,x,y);
      double xh = ratio->GetErrorXhigh(i);
      double xl = ratio->GetErrorXlow(i);
      double yh = ratio->GetErrorYhigh(i);
      double yl = ratio->GetErrorYlow(i);
      double f_ = f1->Eval(x);
      if (y<=0 || x<fit_min[hname] || x>fit_max[hname]) y =-5;
      ratio2->SetPoint(i,x,y/f_);
      ratio2->SetPointError(i,xh,xl,yh/f_,yl/f_);
    }
    tdrDraw(ratio2.get(), "P0", kFullDotLarge, color, kSolid, color, 3004, color);
  }
  c_ratios->cd(1);
  leg_ratios->Draw("same");

  c_ratios->SaveAs((workingDir+"Plot_Ratios_"+histFolder+"."+plotting_mode).c_str());
  if (plotting_mode_2!="") c_ratios->SaveAs((workingDir+"Plot_Ratios_"+histFolder+"."+plotting_mode_2).c_str());

}


void CreateRooWorkspace::PlotTranferFunction() {

  std::string mode_CR = BkgName+"_CR";
  std::string mode_SR = BkgName+"_SR";

  std::unique_ptr<TCanvas> c_TF(tdrCanvas("TF", plot_lo, plot_hi, 1e-03, 1e00, nameXaxis, "A.U"));
  c_TF->SetLogy(1);
  std::unique_ptr<TLegend> leg_TF(tdrLeg(0.40,0.60,0.7,0.85, 0.035));
  tdrHeader(leg_TF.get(), "Transfer Fuctions", 12, 0.04, 42, kBlack, true);
  gStyle->SetOptFit(kFALSE); //TODO put it in one place

  std::unique_ptr<TGraphAsymmErrors> ratio_hist(new TGraphAsymmErrors());
  ratio_hist->Divide(histo_map[mode_SR].get(), histo_map[mode_CR].get(), "pois");
  std::unique_ptr<TF1> f1(new TF1("TF","pol0",fit_min[mode_CR],fit_max[mode_CR]));
  ratio_hist->Fit(f1.get(),"RQ");
  f1->SetLineColor(kRed);
  f1->SetLineWidth(4);
  f1->Draw("same");

  // for (auto const& [model,dofit] : doFits_map[mode_CR] ) {
  //   if (!dofit || (model.Contains("Exp") && model!="Exp_2" && model!="Exp_3" && model!="Exp_4") ) continue;
  //   long int nEvents = 100000000000;
  //   TGraphAsymmErrors* ratio_func = new TGraphAsymmErrors();
  //   TH1F* h_func_CR = (TH1F*)Fits_map[mode_CR][model]->generateBinned(*x_var,nEvents,true)->createHistogram("func_"+model+"_CR",*x_var);
  //   TH1F* h_func_SR = (TH1F*)Fits_map[mode_SR][model]->generateBinned(*x_var,nEvents,true)->createHistogram("func_"+model+"_CR",*x_var);
  //   h_func_CR->Scale(histo_map[mode_CR]->Integral(doBinWidth?"width":"")/h_func_CR->Integral(doBinWidth?"width":""));
  //   h_func_SR->Scale(histo_map[mode_SR]->Integral(doBinWidth?"width":"")/h_func_SR->Integral(doBinWidth?"width":""));
  //   ratio_func->Divide(h_func_SR, h_func_CR);
  //   tdrDraw(ratio_func, "L", kFullDotLarge, Colors[model], kSolid, Colors[model], 3000, Colors[model]);
  //   leg_TF->AddEntry(ratio_func, model ,"l");
  // }
  tdrDraw(ratio_hist.get(), "P0", kFullDotLarge, kBlack, kSolid, kBlack, 3004, kBlack);
  leg_TF->AddEntry(ratio_hist.get(), Form("%s/%s", histo_map[mode_SR]->GetName(),histo_map[mode_CR]->GetName()) ,"lep");
  c_TF->SaveAs(workingDir+"TFs_"+histFolder+extra_text+"."+plotting_mode);
  if (plotting_mode_2!="") c_TF->SaveAs(workingDir+"TFs_"+histFolder+extra_text+"."+plotting_mode_2);

}


void CreateRooWorkspace::InputDatacards(){
  /*
  # #### ##    ## ########  ##     ## ######## ########     ###    ########    ###     ######     ###    ########  ########   ######
  #  ##  ###   ## ##     ## ##     ##    ##    ##     ##   ## ##      ##      ## ##   ##    ##   ## ##   ##     ## ##     ## ##    ##
  #  ##  ####  ## ##     ## ##     ##    ##    ##     ##  ##   ##     ##     ##   ##  ##        ##   ##  ##     ## ##     ## ##
  #  ##  ## ## ## ########  ##     ##    ##    ##     ## ##     ##    ##    ##     ## ##       ##     ## ########  ##     ##  ######
  #  ##  ##  #### ##        ##     ##    ##    ##     ## #########    ##    ######### ##       ######### ##   ##   ##     ##       ##
  #  ##  ##   ### ##        ##     ##    ##    ##     ## ##     ##    ##    ##     ## ##    ## ##     ## ##    ##  ##     ## ##    ##
  # #### ##    ## ##         #######     ##    ########  ##     ##    ##    ##     ##  ######  ##     ## ##     ## ########   ######
  */


  for (auto mode: Modes) {
    for (auto const& [model,dofit] : doFits_map[mode] ) {
      if (!dofit || (model.find("Exp")!= std::string::npos && model!="Exp_2" && model!="Exp_3" && model!="Exp_4") ) continue;

      x_var->setRange("fitting", fit_lo, fit_hi);
      double i_fit = ((RooAbsReal*)((RooAbsPdf*)Fits_map[mode][model].get())->createIntegral(*x_var, RooFit::NormSet(*x_var), RooFit::Range("fitting")))->getVal();
      x_var->setRange("total", x_lo, x_hi);
      double i_tot = ((RooAbsReal*)((RooAbsPdf*)Fits_map[mode][model].get())->createIntegral(*x_var, RooFit::NormSet(*x_var), RooFit::Range("total")))->getVal();
      // std::cout << mode << " " << model << " test integral " << i_tot << " " << i_fit << std::endl;

      // EXPLAINATION OF INTEGRALS:
      // CalculateIntegral(histo_map[mode],x_lo,x_hi,doBinWidth) = Integral over the range where the histogram is defined, using TH1F
      // rooHist_map[mode]->sum(!doBinWidth) = Integral over the range where the histogram is defined, using RooDataHist
      // CalculateIntegral(histo_map[mode],fit_lo,fit_hi,doBinWidth) = Integral over the fit range, using TH1F
      // CalculateFractionArea(histo_map[mode],fit_lo,fit_hi, x_lo, x_hi,doBinWidth) = Fraction of fitting area over the full range where the histo is defined
      // CalculateIntegral(histo_map[mode],fit_lo,fit_hi,doBinWidth)/CalculateFractionArea(histo_map[mode],fit_lo,fit_hi, x_lo, x_hi,doBinWidth) = Should give the same value of the integral in the full range
      // Normalization in the full range is rooHist_map[mode]->sum(!doBinWidth). Needed in combine because the histo is defined here.
      // i_tot = Integral over the full range, using RooAbsPdf. Should be 1 because it's normalized
      // i_fit = Integral over the fit range, using RooAbsPdf
      // CalculateIntegral(histo_map[mode],fit_lo,fit_hi,doBinWidth)*i_tot/i_fit = Normalization of the RooAbsPdf, such that it has the correct norm in the fitting range.
      // nEventsSR*i_tot/i_fit = Normalization of the RooAbsPdf, such that it has the correct norm in the fitting range as the event in the SR.

      if (debug) std::cout << mode+unique_name+"_"+model << " integral " << nEventsSR*i_tot/i_fit << " " << CalculateIntegral(histo_map[mode].get(),fit_lo,fit_hi,doBinWidth)*i_tot/i_fit << " " << CalculateIntegral(histo_map[mode].get(),fit_lo,fit_hi,doBinWidth)/CalculateFractionArea(histo_map[mode].get(),fit_lo,fit_hi, x_lo, x_hi,doBinWidth) << " " << rooHist_map[mode]->sum(!doBinWidth) << std::endl;
      DataCard  << mode+unique_name+"_"+model << " integral " << nEventsSR*i_tot/i_fit << " " << CalculateIntegral(histo_map[mode].get(),fit_lo,fit_hi,doBinWidth)*i_tot/i_fit << " " << CalculateIntegral(histo_map[mode].get(),fit_lo,fit_hi,doBinWidth)/CalculateFractionArea(histo_map[mode].get(),fit_lo,fit_hi, x_lo, x_hi,doBinWidth) << " " << rooHist_map[mode]->sum(!doBinWidth) << std::endl;


      std::unique_ptr<RooArgSet> model_params(Fits_map[mode][model]->getParameters(*x_var));
      TString name; int from = 0;
      while (TString(model_params->contentsString()).Tokenize(name, from, ",")) {
        RooRealVar* par = dynamic_cast<RooRealVar*>(model_params->find(name));
        DataCard << Form("%s%sparam %.3f %.3f", name.Data(), std::string(60-name.Length(),' ').c_str(), par->getVal(), par->getError()) <<std::endl;
      }
    }
  }


  for (const int & mass : MassPoints) {
    std::string SgName = GetSgName(mass);
    //TODO fix the norm of signal
    // DataCard << SgName+" signal number of events = " << histo_map[SgName]->GetSumOfWeights()<<""<<std::endl;
    // DataCard << SgName+" signal number of events = " << histo_map[SgName]->Integral(histo_map[SgName]->FindBin(rangeLo-100),histo_map[SgName]->FindBin(rangeHi+100), doBinWidth?"width":"")<<""<<std::endl;
    DataCard << SgName+" signal number of events = " << nEventsSignal[SgName] <<""<<std::endl;
    for (unsigned int i = 0; i < fitPars[GetSgName(mass)].size(); i++) {
      DataCard << SgName << " " << fitPars[GetSgName(mass)][i]->GetName() << "   param   " <<fitPars[GetSgName(mass)][i]->getVal() << " " << fitPars[GetSgName(mass)][i]->getError() <<  std::endl;
    }
  }

}





void CreateRooWorkspace::Process() {
  /*
  # ########  ########   #######   ######  ########  ######   ######
  # ##     ## ##     ## ##     ## ##    ## ##       ##    ## ##    ##
  # ##     ## ##     ## ##     ## ##       ##       ##       ##
  # ########  ########  ##     ## ##       ######    ######   ######
  # ##        ##   ##   ##     ## ##       ##             ##       ##
  # ##        ##    ##  ##     ## ##    ## ##       ##    ## ##    ##
  # ##        ##     ##  #######   ######  ########  ######   ######
  */


  LoadFiles();
  PrepocessHistos();
  NormaliseData();
  DoRebin(); //TODO put it in the correct place
  InitializePDFs();
  CreateRooDataHist();
  DoFits();
  ImportToWorkspace();
  InputDatacards();
  DoPlots();
};













void CreateRooWorkspace::PlotBkgFit() {
  plotter.reset(x_var->frame(plot_lo,plot_hi));
  double pull_min = GetRange(histo_map["bkg_pred"].get(),pull_lo);
  double pull_max = GetRange(histo_map["bkg_pred"].get(),pull_hi);

  for (auto mode: Modes) {
    if (mode!="bkg_pred" && mode!="DY_CR" && mode!="DY_SR") continue;
    // if (mode!="bkg_pred" && mode!="DY_CR") continue;
    // if (mode!="bkg_pred") continue;

    std::unordered_map<std::string, RooHist*> hpull;
    std::unordered_map<std::string, RooHist*> hratio;
    std::unordered_map<std::string, double> chi2_map;

    TCanvas* c_bg = tdrDiCanvas(("Events"+mode).c_str(), plot_lo, plot_hi, plot_ylo, plot_yhi, doPlotRatio?0.8:-6, doPlotRatio?1.2:6, nameXaxis, nameYaxis, nameRatioaxis);
    c_bg->cd(1)->SetLogy(1);
    // plotter = x_var->frame(plot_lo,plot_hi);
    plotter.reset(x_var->frame(plot_lo,plot_hi));

    if (mode=="bkg_pred" || mode=="data") rooHist_map[mode]->plotOn(plotter.get(),RooFit::DataError(RooAbsData::Poisson));
    else rooHist_map[mode]->plotOn(plotter.get(),RooFit::DataError(RooAbsData::SumW2));
    // rooHist_map[mode]->plotOn(plotter.get(),RooFit::DataError(RooAbsData::SumW2));
    // rooHist_map[mode]->plotOn(plotter.get());

    std::unique_ptr<RooCurve> curve;
    std::unique_ptr<RooHist> hist_;
    for (auto const& [model,dofit] : doFits_map[mode] ) {
      if (dofit) {
        int color = Colors[model];
        // if (!doFtest && model.Contains("Exp") && model!="Exp_2" && model!="Exp_3" && model!="Exp_4" ) continue;
        Fits_map[mode][model]->plotOn(plotter.get(), RooFit::LineColor(color), RooFit::Range(pull_min, pull_max, kFALSE));
        // Fits_map[mode][model]->plotOn(plotter, RooFit::LineColor(color), RooFit::Range(fit_min[mode], fit_max[mode], kFALSE));

        // CREATE RATIO BETWEEN CURVE AND HIST
        double xstart,xstop,y ;
        curve.reset((RooCurve*) ( (RooCurve*) plotter->findObject(0,RooCurve::Class()))->Clone(""));
        hist_.reset((RooHist*) ( (RooHist*) plotter->findObject(0,RooHist::Class()))->Clone(""));
        curve->GetPoint(0,xstart,y) ;
        curve->GetPoint(curve->GetN()-1,xstop,y) ;
        hratio[model] = new RooHist() ;
        for(int i=0 ; i<hist_->GetN() ; i++) {
          double x,point;
          hist_->GetPoint(i,x,point) ;
          if (x<xstart || x>xstop) continue ;
          double yy = point/ curve->interpolate(x) ;
          double dyl = 2*hist_->GetErrorYlow(i)/ curve->interpolate(x);
          double dyh = 2*hist_->GetErrorYhigh(i)/ curve->interpolate(x);
          hratio[model]->addBinWithError(x,yy,dyl,dyh);
        }

        // CREATE PULL BETWEEN CURVE AND HIST
        hpull[model] = plotter->pullHist();
        int npf = FitRes_map[mode][model]->floatParsFinal().getSize();
        int ndf = hpull[model]->GetN()-npf;
        chi2_map[model] = plotter->chiSquare(npf);
        double pv = TMath::Prob(chi2_map[model]*ndf, ndf)*100;
        TString name = model+std::string(7-model.size(), ' ' )+" #chi^{2}/n.d.f. = "+TString::Format("%.1f",chi2_map[model])+" p-value = "+TString::Format("%.2f",pv)+"\%";
        // TString name = model+std::string(5-model.size(), ' ' )+"np="+TString::Format("%d",npf);
        hpull[model]->SetName(name);
        hpull[model]->SetMarkerColor(color);
        hpull[model]->SetLineColor(color);
        hratio[model]->SetName(name);
        hratio[model]->SetMarkerColor(color);
        hratio[model]->SetLineColor(color);
        Fits_map[mode][model]->plotOn(plotter.get(), RooFit::Range(fit_min[mode], fit_max[mode], kFALSE), RooFit::FillColor(color), RooFit::VisualizeError(*FitRes_map[mode][model], 1), RooFit::FillStyle(3001), RooFit::VLines());
      }
    }

    plotter->Draw("same");

    std::unique_ptr<TLegend>leg_bg(tdrLeg(0.40,0.55,0.7,0.85, 0.032));
    tdrHeader(leg_bg.get(), histo_map[mode]->GetName(), 12, 0.04, 42, kBlack, true);

    std::unique_ptr<RooPlot> Pullplotter_bkg(x_var->frame(RooFit::Range(plot_lo,plot_hi), RooFit::Name("bkg")));
    for (auto x : doPlotRatio?hratio:hpull) {
      std::string model = x.first;
      if (model.find("Exp_")!=std::string::npos || model.find("BkgPdf")!=std::string::npos) {
        int deg = atoi(new char((model.find("Exp_")!=std::string::npos)? model[4] : model[6]));
        std::string old_model = TString(x.first).ReplaceAll(std::to_string(deg),std::to_string(deg-1)).Data();
        if (chi2_map.find(old_model)!= chi2_map.end()) {
          double FTest = DoFTest(chi2_map[old_model],chi2_map[model],deg-1,deg,hpull[model]->GetN());
          if (doFtest) output << mode << " " << model << " FTest " << FTest << " chi2: " << chi2_map[model] << '\n';
          x.second->SetName(x.second->GetName()+TString::Format(" F-test = %.2f",FTest));
        }
      }
      leg_bg->AddEntry(x.second, x.second->GetName() ,"l");
      Pullplotter_bkg->addPlotable(x.second,"P same");
    }
    leg_bg->Draw("same");

    c_bg->cd(2);
    Pullplotter_bkg->SetNdivisions(505,"Y");
    Pullplotter_bkg->Draw("same");
    std::unique_ptr<TLine> line(new TLine(plot_lo, 0, plot_hi, 0));
    line->SetLineWidth(2);
    line->Draw("same");

    c_bg->SaveAs(workingDir+"Fit_Bg_all_"+mode+"_"+histFolder+extra_text+"."+plotting_mode);
    if (plotting_mode_2!="") c_bg->SaveAs(workingDir+"Fit_Bg_all_"+mode+"_"+histFolder+extra_text+"."+plotting_mode_2);


    if (doCheckPlots) {
      for (auto [model,pull] : hpull) {
        TCanvas* c_gauss = tdrCanvas(("gaus_"+mode+"_"+model).c_str(), -6, 6, 0.00001,100, nameXaxis, nameYaxis);
        std::unique_ptr<TH1F> gauss(new TH1F(model.c_str(), model.c_str(), 13, -6, 6));
        std::unique_ptr<TF1> f1(new TF1("f1","gaus",-5,5));
        for (int i = 0; i < pull->GetN(); i++) {
          double x,y;
          pull->GetPoint(i,x,y);
          gauss->Fill(y);
        }
        gauss->Fit(f1.get(),"RMQ");
        tdrDraw(gauss.get(), "Hist", kSolid, kBlack, kSolid, kBlack, 3000, kBlack);
        std::unique_ptr<TLegend> leg_gauss(tdrLeg(0.70,0.60,0.9,0.85, 0.035));
        tdrHeader(leg_gauss.get(), mode+" "+model, 12, 0.04, 42, kBlack, true);
        c_gauss->SaveAs(workingDir+"Gauss_Bg_all_"+mode+"_"+model+"_"+histFolder+extra_text+"."+plotting_mode);
        if (plotting_mode_2!="") c_gauss->SaveAs(workingDir+"Gauss_Bg_all_"+mode+"_"+model+"_"+histFolder+extra_text+"."+plotting_mode_2);
      }
    }
  }
}


int main(int argc, char** argv){
  gErrorIgnoreLevel = kFatal;
  RooMsgService::instance().setSilentMode(true);
  RooMsgService::instance().setStreamStatus(1,false);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::MsgTopic::Plotting);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::MsgTopic::InputArguments);
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

  bool isHbb = false;
  // std::vector<std::string> histFolders = {"btag_DeepBoosted_H4qvsQCD"};
  std::vector<std::string> histFolders = {"btag_DeepBoosted_H4qvsQCD", "btag_DeepBoosted_H4qvsQCDptdep", "btag_DeepBoosted_H4qvsQCDp2", "btag_DeepBoosted_H4qvsQCDp02", "btag_DeepBoosted_H4qvsQCDpt1000", "btag_DeepBoosted_H4qvsQCDpt1000p2", "btag_DeepBoosted_H4qvsQCDpt1000p02",
  "btag_DeepBoosted_H4qvsQCDptdep_x3", "btag_DeepBoosted_H4qvsQCDptdep_x2x3", "btag_DeepBoosted_H4qvsQCDptdep_x1x3", "btag_DeepBoosted_H4qvsQCDmassdep_x3",
  "btag_DeepBoosted_H4qvsQCDmassdep2_x3", "btag_DeepBoosted_H4qvsQCDmassdep_x2x3", "btag_DeepBoosted_H4qvsQCDmassdep_x1x3", "btag_DeepBoosted_H4qvsQCDmassdep_x1x2" };
  if (isHbb) histFolders = {"btag_DeepBoosted_HbbvsQCD", "btag_DeepBoosted_probHbb", "tau21" };
  std::vector<std::string> collections = {"Puppi"};
  std::vector<std::string> channels = {"muonchannel", "electronchannel"};
  std::vector<std::string> years = {"2016", "2017", "2018", "RunII"};



  std::unique_ptr<CreateRooWorkspace> roo;

  if (argc>1) {
    std::string histFolder, channel, collection, year;
    for (int i = 1; i < argc; i++) {
      if (std::find(collections.begin(), collections.end(), argv[i]) != collections.end() ) collection = argv[i];
      if (std::find(histFolders.begin(), histFolders.end(), argv[i]) != histFolders.end() ) histFolder = argv[i];
      if (std::find(channels.begin(), channels.end(), argv[i]) != channels.end() ) channel = argv[i];
      if (std::find(years.begin(), years.end(), argv[i]) != years.end() ) year = argv[i];
    }
    std::cout << histFolder << " " << channel << " " << collection << " " << year << '\n';
    roo.reset(new CreateRooWorkspace(year,collection, channel, histFolder));
    roo->Process();

  } else {
    std::cout << "more " << '\n';
    for (std::string year: years) {
      for (std::string collection: collections) {
        for (std::string channel: channels) {
          for (std::string histFolder: histFolders) {
            roo.reset(new CreateRooWorkspace(year,collection, channel, histFolder));
            roo->Process();
          }
        }
      }
    }
  }

  return 0;

}
