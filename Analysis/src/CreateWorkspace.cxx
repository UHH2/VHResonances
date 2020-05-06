#include "CreateWorkspace.hpp"
#include "TSpline.h"

RooRealVar* ReDoVar(RooRealVar& var, TString name) { return new RooRealVar(var,var.getTitle().ReplaceAll("bg",name));}

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
    if (x<xmin || x>xmax || pull<=0) continue;
    chi2+=pull*pull;
    nbins++;
  }
}


void CreateWorkspace(std::string studies, std::string histFolder, std::string channel, std::string collection, std::string year, bool isHbb, bool fitCR, bool doObs) {

  std::cout << "****************************************" << std::endl;
  std::cout << "             CreateWorkspace            " << std::endl;
  std::cout << "----------------------------------------" << std::endl;
  std::cout << "histFolder" << std::string(15-10,' ') << histFolder << '\n';
  std::cout << "channel"    << std::string(15-7, ' ') << channel    << '\n';
  std::cout << "collection" << std::string(15-10,' ') << collection << '\n';
  std::cout << "year"       << std::string(15-4, ' ') << year       << '\n';
  std::cout << "****************************************\n" << std::endl;

  /**************************************************
  *             PLOTTING VARIABLES                  *
  **************************************************/
  writeExtraText = true;       // if extra text
  extraText  = "Work in progress" ;//"Preliminary";
  lumi_13TeV  = TString::Format("%.1f fb^{-1}", lumi_map.at(year).at("lumi_fb"));

  std::map<TString, int> Colors;

  Colors["Landau"]    = kRed+1;
  Colors["Flatte"]    = kRed+1;
  Colors["LN"]        = kRed+1;
  Colors["LE"]        = kBlue+1;
  Colors["GE"]        = kGreen+1;
  Colors["BkgPdf3p"]  = kGreen+1;
  Colors["BkgPdf4p"]  = kRed+1;

  Colors["NO"]        = kBlue+1;
  Colors["CB"]        = kOrange+1;
  Colors["Exp_1"]      = kRed;
  Colors["Exp_2"]      = kMagenta+1;
  Colors["Exp_3"]      = kGreen+3;
  Colors["Exp_4"]      = kGreen+1;
  Colors["Exp_5"]      = kYellow+1;
  Colors["Exp_6"]      = kBlack;

  Colors["bkg_pred"]  = kBlue+1;
  Colors["data"]      = kRed+1;
  Colors["main_bkg_CR"] = kGreen+1;
  Colors["main_bkg_SR"] = kOrange+1;
  Colors["main_bkg"]  = kAzure+7;
  Colors["extra_bkg"] = kOrange+1;
  Colors["DY_CR"]     = kGreen-2;
  Colors["DY_SR"]     = kOrange-2;
  Colors["DY"]        = kOrange-2;
  Colors["TTbar"]     = kOrange+10;
  Colors["VV"]        = kGreen+2;
  Colors["WW"]        = kGreen;
  Colors["WZ"]        = kGreen+3;
  Colors["ZZ"]        = kGreen-10;

  /**************************************************
  *               VARIABLE DEFINITIONS              *
  **************************************************/

  std::string Module    = "SignalRegion";
  std::string Histtype  = "ZprimeCandidate";
  // std::string HistName  = "Zprime_mass_rebin_full";
  // std::string HistName  = "Zprime_mass";
  // std::string HistName  = "Zprime_mass_rebin1";
  std::string HistName  = "Zprime_mass_rebin2";
  int rebin = 30;
  bool dorebin = true;
  bool doBinWidth = false;
  int bin = 50;
  int bin2 = 100;
  std::vector<double> bins_Zprime_rebin;
  for (double i = 0; i < 1500; i+=bin) bins_Zprime_rebin.push_back(i);
  for (double i = 1500+bin; i <= 4000; i+=bin2) bins_Zprime_rebin.push_back(i);
  bins_Zprime_rebin.push_back(4000);


  bool debug = false;
  // debug = true;
  bool doCheckPlots = true;
  bool doBkgPlots = true;
  // bool doCheckPlots = false;
  // bool doBkgPlots = false;
  bool doSignalPlots = true;
  bool doFtest = false;
  bool doPlotRatio = true;
  // bool doPlotRatio = false;

  // std::string plotting_mode = "pdf";
  std::string plotting_mode = "eps";
  std::string plotting_mode_2 = "pdf";

  // TString nameYaxis = "Events/GeV";
  TString nameYaxis = doBinWidth? "Events/bin" :"Events";
  TString nameXaxis = "m(Z')";
  TString nameRatioaxis = doPlotRatio?"Hist/Fit": "Pull";

  TString extra_text = doFtest? "_Ftest_": "";

  double x_lo     = 200;
  double x_hi     = 10000;
  double plot_lo  = 200;
  double plot_hi  = 4000;
  double plot_ylo = 1.1*1e-03;
  double plot_yhi = 1e07;
  // double x_lo     = 750;
  // double x_hi     = 2700;
  // double plot_lo  = 750;
  // double plot_hi  = 2700;
  // double fit_lo_turnOn = 580; // not used
  // double fit_hi_tails = 3560; // not used
  // double fit_lo   = 580;
  // double fit_hi   = 3560;
  double fit_lo   = 600;
  double fit_hi   = 3000;
  double pull_lo  = 600;
  double pull_hi  = 3000;


  // std::vector<TString> BkgNames = {"DY", "TTbar", "WZ", "WW", "ZZ"}; //TODO
  std::vector<TString> BkgNames = {"DY", "TTbar", "WZ","ZZ"};
  std::string BkgName   = "DY";
  std::string dataName  = channel.substr(0, channel.find("channel")); dataName[0] = toupper(dataName[0]) ; dataName = "DATA_Single"+dataName;

  /*
  # ##        #######     ###    ########     ##     ## ####  ######  ########  #######   ######   ########     ###    ##     ##  ######
  # ##       ##     ##   ## ##   ##     ##    ##     ##  ##  ##    ##    ##    ##     ## ##    ##  ##     ##   ## ##   ###   ### ##    ##
  # ##       ##     ##  ##   ##  ##     ##    ##     ##  ##  ##          ##    ##     ## ##        ##     ##  ##   ##  #### #### ##
  # ##       ##     ## ##     ## ##     ##    #########  ##   ######     ##    ##     ## ##   #### ########  ##     ## ## ### ##  ######
  # ##       ##     ## ######### ##     ##    ##     ##  ##        ##    ##    ##     ## ##    ##  ##   ##   ######### ##     ##       ##
  # ##       ##     ## ##     ## ##     ##    ##     ##  ##  ##    ##    ##    ##     ## ##    ##  ##    ##  ##     ## ##     ## ##    ##
  # ########  #######  ##     ## ########     ##     ## ####  ######     ##     #######   ######   ##     ## ##     ## ##     ##  ######
  */

  std::string user            = std::getenv("USER");
  std::string Path_ANALYSIS   = std::getenv("CMSSW_BASE"); Path_ANALYSIS += "/src/UHH2/HiggsToWWTagger/";
  std::string Path_NFS        = "/nfs/dust/cms/user/"+user+"/";
  std::string Path_SFRAME     = Path_NFS+"sframe_all/";
  std::string Path_STORAGE    = Path_NFS+"WorkingArea/File/Analysis/";
  std::string Path_SPlotter   = Path_NFS+"WorkingArea/SFramePlotter/";

  std::string AnalysisPath    = Path_ANALYSIS+"Analysis/";
  std::string StoragePath     = Path_STORAGE+year+"/";
  std::string PrefixrootFile  = "uhh2.AnalysisModuleRunner.";

  std::string syst = "nominal";
  std::string unique_name_complete = year+"_"+collection+"_"+channel+"_"+histFolder;
  std::string unique_name = "_"+year+"_"+collection+"_"+channel;
  std::string filepath    = StoragePath+Module+"/"+collection+"/"+channel+"/"+syst+"/";
  std::string workingDir  = AnalysisPath+(isHbb? "Limits/Hbb_": "Limits/")+studies+"/"+year+"/"+collection+"/"+channel+"/"+histFolder+"/";
  gSystem->Exec(("mkdir -p "+workingDir+"/datacards").c_str());

  std::string dataFileName = PrefixrootFile+"DATA."+dataName+"_"+year+"_noTree.root";
  std::string bkgFileName  = PrefixrootFile+"MC.MC_"+BkgName+"_"+year+"_noTree.root";
  std::map<TString, TString> bkgFileName_map;
  for (auto bkg: BkgNames) bkgFileName_map[bkg] = PrefixrootFile+"MC.MC_"+bkg+"_"+year+"_noTree.root";

  // std::string SGname = Histtype+"_HWW"+histFolder+"_SR/"+HistName; // WE DON'T WANT THE BR INVOLVED
  std::string SGname = Histtype+"_"+histFolder+"_SR/"+HistName;
  std::string SRname = Histtype+"_"+histFolder+"_SR/"+HistName;
  std::string CRname = Histtype+"_"+histFolder+"_CR/"+HistName;

  // TFile file_WS("WS.root","RECREATE");
  std::ofstream DataCard;
  DataCard.open (workingDir+"datacards/OutputFit_"+histFolder+".txt");
  RooWorkspace *ws=new RooWorkspace((channel+"_"+year).c_str());

  TFile *f_data = new TFile((filepath+dataFileName).c_str());
  std::map<TString, TFile *> f_bkg_map;
  for (auto [bkg, name]: bkgFileName_map) f_bkg_map[bkg] = new TFile(filepath+name);

  if (debug) {
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "Opening\t" << filepath+dataFileName << '\n';
    for (auto x : bkgFileName_map) std::cout << "Opening\t" << filepath+x.second << '\n';
    std::cout << "folder: SGname " << SGname << '\n';
    std::cout << "folder: SRname " << SRname << '\n';
    std::cout << "folder: CRname " << CRname << '\n';
    std::cout << "----------------------------------------" << std::endl;
  }


  std::map<TString, TH1F *> histo_map;
  for (auto [bkg, file]: f_bkg_map) {
    histo_map[bkg+"_CR"] = (TH1F*)((TH1F*)file->Get((CRname).c_str()))->Clone(bkg+" CR");
    histo_map[bkg+"_SR"] = (TH1F*)((TH1F*)file->Get((SRname).c_str()))->Clone(bkg+" SR");
  }

  if (debug) for (auto x : histo_map) std::cout << "histo " << x.first << std::string(10-x.first.Length(),' ') << x.second->Integral() << '\n';

  histo_map["main_bkg_CR"] = (TH1F*)((TH1F*)f_bkg_map[BkgName]->Get((CRname).c_str()))->Clone(("main_bkg_CR: "+BkgName).c_str());
  histo_map["main_bkg_SR"] = (TH1F*)((TH1F*)f_bkg_map[BkgName]->Get((SRname).c_str()))->Clone(("main_bkg_SR: "+BkgName).c_str());
  histo_map["extra_bkg_CR"] = (TH1F*)((TH1F*)f_bkg_map["TTbar"]->Get((CRname).c_str()))->Clone("extra_bkg_CR: TTbar");
  histo_map["extra_bkg_SR"] = (TH1F*)((TH1F*)f_bkg_map["TTbar"]->Get((SRname).c_str()))->Clone("extra_bkg_SR: TTbar");
  histo_map["VV_CR"] = (TH1F*)((TH1F*)f_bkg_map["WZ"]->Get((CRname).c_str()))->Clone("extra_bkg_CR: VV");
  histo_map["VV_SR"] = (TH1F*)((TH1F*)f_bkg_map["WZ"]->Get((SRname).c_str()))->Clone("extra_bkg_SR: VV");

  for (auto bkg: BkgNames) {
    histo_map[bkg+"_ratio"]= (TH1F*)(histo_map[bkg+"_SR"])->Clone(bkg+"_ratio");
    for (std::string reg: {"CR","SR"}) {
      if (bkg==BkgName) continue;
      histo_map["main_bkg_"+reg]->Add(histo_map[bkg+"_"+reg]);
      histo_map["main_bkg_"+reg]->SetName((std::string(histo_map["main_bkg_"+reg]->GetName())+"+"+bkg));
      if (bkg=="TTbar") continue;
      histo_map["extra_bkg_"+reg]->Add(histo_map[bkg+"_"+reg]);
      histo_map["extra_bkg_"+reg]->SetName((std::string(histo_map["main_bkg_"+reg]->GetName())+"+"+bkg));
      if (bkg=="WZ") continue;
      histo_map["VV_"+reg]->Add(histo_map[bkg+"_"+reg]);
    }
  }


  if (doObs) histo_map["data"] = (TH1F*)((TH1F*)f_data->Get((SRname).c_str()))->Clone("data: data SR");
  else histo_map["data"] = (TH1F*)((TH1F*)f_bkg_map[BkgName]->Get((SRname).c_str()))->Clone(("h_data_fake: "+BkgName+" SR").c_str());
  histo_map["norm"] = (TH1F*)((TH1F*)f_data->Get((SRname).c_str()))->Clone("data: data SR");

  if (fitCR) histo_map["bkg_pred"]  = (TH1F*)((TH1F*)f_data->Get((CRname).c_str()))->Clone("bkg_pred: data CR");
  else histo_map["bkg_pred"]  = (TH1F*)((TH1F*)f_bkg_map[BkgName]->Get((SRname).c_str()))->Clone(("bkg_pred: "+BkgName+" SR").c_str());


  for (TString name : {"main", "extra"}) {
    histo_map[name+"_bkg_ratio"]= (TH1F*)(histo_map[name+"_bkg_SR"])->Clone((name+"_bkg_ratio"));
  }

  if (debug){
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "CR hist_name:\t" << CRname << '\n';
    std::cout << "SR hist_name:\t" << SRname << '\n';
    std::cout << "fitCR:\t" << fitCR << "\t" << histo_map["bkg_pred"]->GetName() << '\n';
    std::cout << "doObs:\t" << doObs << "\t" << histo_map["data"]->GetName() << '\n';
    std::cout << "----------------------------------------" << std::endl;
  }

  std::vector<TString> Modes = {"bkg_pred", "data", "main_bkg_CR", "main_bkg_SR", "DY_CR", "DY_SR"};

  // int normFactor = histo_map["DY_CR"]->Integral()/histo_map["DY_SR"]->Integral(); //TODO

  for (auto [mode, h]: histo_map) {

    if (debug) std::cout << "pre " << mode << "\t" << CalculateIntegral(h,fit_lo,fit_hi,doBinWidth) << '\n';

    if (dorebin) {
      if (rebin) h->Rebin(rebin);
      else h = dynamic_cast<TH1F*>(h->Rebin(bins_Zprime_rebin.size()-1, h->GetName(), &bins_Zprime_rebin[0]));

      h->Scale(1,doBinWidth?"width":"");
    }
    if (debug) std::cout << "rebin " << mode << "\t" << CalculateIntegral(h,fit_lo,fit_hi,doBinWidth) << '\n';

    // Removing bins with low stat TODO
    if (mode=="DY_SR") {
      for (int i = 0; i < h->GetNbinsX()+1; i++) {
        if (h->GetBinContent(i)<2*1e-04) { h->SetBinContent(i,0); h->SetBinError(i,0); }
        // if (h->GetBinCenter(i)>1280 && h->GetBinCenter(i)<1320) { h->SetBinContent(i,0); h->SetBinError(i,0); }
      }
    }

    // if (mode=="DY_CR") {
    //   for (size_t i = 0; i < h->GetNbinsX()+1; i++) {
    //     h->SetBinContent(i,h->GetBinContent(i)/normFactor);
    //     h->SetBinError(i,h->GetBinError(i)/TMath::Sqrt(normFactor));
    //   }
    //   h->Scale(normFactor);
    // }
  }
  // histo_map["DY_SR"]->Scale(2); //TODO

  for (auto [name, h]: histo_map) {
    if (name.Contains("ratio")) {
      h->Divide(histo_map[TString(h->GetName()).ReplaceAll("SR","CR")]);
    }
  }

  /*
  # ########     ###    ########    ###              #######  ########   ######
  # ##     ##   ## ##      ##      ## ##            ##     ## ##     ## ##    ##
  # ##     ##  ##   ##     ##     ##   ##           ##     ## ##     ## ##
  # ##     ## ##     ##    ##    ##     ##          ##     ## ########   ######
  # ##     ## #########    ##    #########          ##     ## ##     ##       ##
  # ##     ## ##     ##    ##    ##     ##          ##     ## ##     ## ##    ##
  # ########  ##     ##    ##    ##     ##           #######  ########   ######
  */


  RooRealVar* x_var = new RooRealVar("x_var", "m_{Zprime} (GeV)", x_lo, x_hi);
  if (dorebin && rebin!=0) {
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
  std::unique_ptr<RooPlot> plotter;
  plotter.reset(x_var->frame(plot_lo,plot_hi));

  // Normalize h_MC_SR to h_Data_SR pretending it's data but has shape of bkg_pred in SR TODO
  double nEventsSR  = CalculateIntegral(histo_map["norm"],fit_lo,fit_hi,doBinWidth);
  if (!doObs) {
    histo_map["data"]->Scale(nEventsSR/CalculateIntegral(histo_map["data"],fit_lo,fit_hi,doBinWidth));
    // histo_map["data"]->Scale(1./CalculateFractionArea(histo_map["data"],fit_lo,fit_hi, x_lo, x_hi,doBinWidth));
    histo_map["data"]->SetName((std::string(histo_map["data"]->GetName())+"*nEventsSR").c_str());
  }

  // For the datacard
  DataCard << "=== RooFit data fit result to be entered in datacard === " << std::endl;
  DataCard << BkgName+" Background number of events = " << nEventsSR << std::endl;

  // histo_map["bkg_pred"]->Scale(nEventsSR/CalculateIntegral(histo_map["bkg_pred"],fit_lo,fit_hi,doBinWidth));//TODO
  RooDataHist* data_obs = new RooDataHist("data_obs", histo_map["data"]->GetName(), RooArgList(*x_var), histo_map["data"]);//TODO check it's imported correctly
  ws->import(*data_obs);

  if (debug) {
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "data_obs\t"  << histo_map["data"]->GetName() << '\n';
    std::cout << "nEventsSR\t" << nEventsSR << '\n';
    std::cout << "norm_data\t" << CalculateIntegral(histo_map["data"],fit_lo,fit_hi,doBinWidth) << '\n';
    std::cout << "----------------------------------------" << std::endl;
  }

  std::map<TString, std::map<TString, RooAbsPdf*> > Fits_map;
  std::map<TString, std::map<TString, bool> > doFits_map;

  for (auto x : Colors) { for (auto mode : Modes) doFits_map[mode][x.first]  = false; }

  // for (auto x : { "NO", "CB", "Exp_1", "Exp_2", "Exp_3", "Exp_4", "Exp_5", "Exp_6"}) { for (auto mode : Modes) doFits_map[mode][x]  = true; }
  for (auto x : {"CB", "Exp_1", "Exp_2", "Exp_3", "Exp_4"}) { for (auto mode : Modes) doFits_map[mode][x]  = true; }
  // for (auto x : { "NO", "CB", "Exp_3"}) { for (auto mode : Modes) doFits_map[mode][x]  = true; }
  // BkgPdf4p, BkgPdf3p TODO
  // for (auto x : { "NO", "CB", "Exp_4"}) { for (auto mode : Modes) doFits_map[mode][x]  = true; }
  // for (auto x : { "CB"}) { for (auto mode : Modes) doFits_map[mode][x]  = true; }

  if (debug) {
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "Fitting Modes:\t\t";
    for (auto x : doFits_map) std::cout << x.first << '\t';
    std::cout << "\nFitting Functions:\t";
    for (auto x : doFits_map["bkg_pred"]) if (x.second) std::cout << x.first << '\t';
    std::cout << '\n';
    std::cout << "----------------------------------------" << std::endl;
  }

  double pull_min = GetRange(histo_map["bkg_pred"],pull_lo);
  double pull_max = GetRange(histo_map["bkg_pred"],pull_hi);
  // double fit_min_turnOn = GetRange(histo_map["bkg_pred"],fit_lo_turnOn); // not used
  // double fit_max_tails = GetRange(histo_map["bkg_pred"],fit_hi_tails); // not used
  // double fit_min = GetRange(histo_map["bkg_pred"],fit_lo);
  // double fit_max = GetRange(histo_map["bkg_pred"],fit_hi);

  std::map<TString, double> fit_min, fit_max;
  for (auto [mode,h] : histo_map) {
    fit_min[mode] = GetRange(h,fit_lo);
    fit_max[mode] = GetRange(h,fit_hi);
  }

  /*
  #  ######  ##     ## ########  ######  ##    ##    ########  ##        #######  ########  ######
  # ##    ## ##     ## ##       ##    ## ##   ##     ##     ## ##       ##     ##    ##    ##    ##
  # ##       ##     ## ##       ##       ##  ##      ##     ## ##       ##     ##    ##    ##
  # ##       ######### ######   ##       #####       ########  ##       ##     ##    ##     ######
  # ##       ##     ## ##       ##       ##  ##      ##        ##       ##     ##    ##          ##
  # ##    ## ##     ## ##       ##    ## ##   ##     ##        ##       ##     ##    ##    ##    ##
  #  ######  ##     ## ########  ######  ##    ##    ##        ########  #######     ##     ######
  */



  if (doCheckPlots) {
    TCanvas* c_CRvsSR = tdrCanvas(("bkg_pred vs data"+unique_name_complete).c_str(), plot_lo, plot_hi, plot_ylo, plot_yhi, nameXaxis, nameYaxis);
    c_CRvsSR->SetLogy(1);
    TLegend *leg_CRvsSR = tdrLeg(0.50,0.70,0.78,0.9, 0.03);
    for (auto mode : Modes) {
      int color = Colors[mode];
      histo_map[mode]->SetMarkerSize(0.5);
      tdrDraw(histo_map[mode], "", kFullDotLarge, color, kSolid, color, 3004, color);
      leg_CRvsSR->AddEntry(histo_map[mode], Form("%s", histo_map[mode]->GetName()) ,"lp");
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
    TLegend *leg_ratios = tdrLeg(0.50,0.70,0.78,0.9, 0.03);
    for (auto [hname,h]: histo_map) {
      if (!hname.Contains("SR")) continue;
      TString name = TString(hname).ReplaceAll("_SR","");
      if (name=="WW") continue;
      if (name=="WZ") continue;
      if (name=="ZZ") continue;
      int color = Colors[name];
      TGraphAsymmErrors* ratio = new TGraphAsymmErrors();
      TGraphAsymmErrors* ratio2 = new TGraphAsymmErrors();
      TH1F* h_SR = (TH1F*)h->Clone();
      TH1F* h_CR = (TH1F*)histo_map[name+"_CR"]->Clone();
      h_SR->Rebin(5);
      h_CR->Rebin(5);
      ratio->Divide(h_SR, h_CR, "pois");
      ratio2->Divide(h_SR, h_CR, "pois");
      // TCanvas* c_ = tdrCanvas(hname, plot_lo, plot_hi, plot_ylo, 1e2, nameXaxis, nameYaxis);
      // c_->SetLogy(1);
      TF1 *f1 = new TF1("f1","pol0",fit_min[hname],fit_max[hname]);
      // TF1 *f2 = new TF1("f2","[0]+TMath::Log10([1]+[2]*x)",500,plot_hi);
      // TF1 *f4 = new TF1("f3","[0]-TMath::Exp([1]+[2]*x)",500,plot_hi);
      // TF1 *f4 = new TF1("f4","[0]-TMath::Exp([1]+TMath::Exp([2]+[3]*x))",500,plot_hi);
      f1->SetParameter(0,0.01);
      f1->SetLineColor(color);
      f1->SetLineWidth(3);
      ratio->Fit(f1,"RS0Q");
      // (TVirtualFitter::GetFitter())->GetConfidenceIntervals(f1,0.68);
      if (BkgName == name) {
        DataCard << name+" rateParam " <<  f1->GetParameter(0) << std::endl;
        // Next number should be ~ 1
        // DataCard << name+" rateParam " <<  f1->GetParameter(0)*CalculateIntegral(histo_map[name+"_CR"],fit_lo,fit_hi,doBinWidth)/CalculateIntegral(histo_map[name+"_SR"],fit_lo,fit_hi,doBinWidth) << std::endl;
      }
      c_ratios->cd(1);
      tdrDraw(ratio, "P0", kFullDotLarge, color, kSolid, color, 3004, color);
      // f4->Draw("same");
      leg_ratios->AddEntry(ratio, name+std::string(name=="TTbar"?8:(name.Length()==2?13:3),' ')+TString::Format("#chi^{2}/n.d.f.=%.2f p-value=%.2f", f1->GetChisquare()/f1->GetNDF(), f1->GetProb()) ,"lp");

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
      tdrDraw(ratio2, "P0", kFullDotLarge, color, kSolid, color, 3004, color);
    }
    c_ratios->cd(1);
    leg_ratios->Draw("same");

    c_ratios->SaveAs((workingDir+"Plot_Ratios_"+histFolder+"."+plotting_mode).c_str());
    if (plotting_mode_2!="") c_ratios->SaveAs((workingDir+"Plot_Ratios_"+histFolder+"."+plotting_mode_2).c_str());
  }

  /*
  # ########     ###     ######  ##    ##  ######   ########   #######  ##     ## ##    ## ########     ######## #### ########
  # ##     ##   ## ##   ##    ## ##   ##  ##    ##  ##     ## ##     ## ##     ## ###   ## ##     ##    ##        ##     ##
  # ##     ##  ##   ##  ##       ##  ##   ##        ##     ## ##     ## ##     ## ####  ## ##     ##    ##        ##     ##
  # ########  ##     ## ##       #####    ##   #### ########  ##     ## ##     ## ## ## ## ##     ##    ######    ##     ##
  # ##     ## ######### ##       ##  ##   ##    ##  ##   ##   ##     ## ##     ## ##  #### ##     ##    ##        ##     ##
  # ##     ## ##     ## ##    ## ##   ##  ##    ##  ##    ##  ##     ## ##     ## ##   ### ##     ##    ##        ##     ##
  # ########  ##     ##  ######  ##    ##  ######   ##     ##  #######   #######  ##    ## ########     ##       ####    ##
  */


  std::cout << "********************" << '\n';
  std::cout << "*  Background FIT  *" << '\n';
  std::cout << "********************" << '\n';

  std::map<TString, RooDataHist *> rooHist_map;

  if (dorebin && rebin!=0) {
    for (auto mode: Modes) rooHist_map[mode] = new RooDataHist(mode, mode, RooArgList(*x_var), RooFit::Import(*histo_map[mode],doBinWidth));
  } else {
    for (auto mode: Modes) rooHist_map[mode] = new RooDataHist(mode, mode, RooArgList(*x_var), histo_map[mode]);
  }

  // TODO?
  // for (auto [mode,h]: rooHist_map) {
  //   // Store RooDataHist format for bkg_pred and data et all
  //   ws->import(*h);
  //   // Store TH1F format for bkg_pred and data et all
  //   ws->import(*histo_map[mode],Form("h_%s",+h->GetName()));
  //   // Store RooHistPdf format for bkg_pred and data et all
  //   ws->import(*(new RooHistPdf( Form("hp_%s",+h->GetName()),"",RooArgSet(*x_var),*h)));
  // }



  if (debug) {
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "Fitting range:\t fit_lo\tfit_hi" << '\n';
    // std::cout << "sub_fit" << std::string(10,' ') << fit_min_turnOn << "\t" << fit_max_tails << std::endl; // not used
    for (auto mode : Modes) std::cout << mode << std::string(17-mode.Length(),' ') << fit_min[mode] << "\t" << fit_max[mode] << std::endl;
    std::cout << "----------------------------------------" << std::endl;
  }

  std::map<TString, std::map<TString, RooFitResult*> > FitRes_map;

  RooRealVar* bg_CB_p1 = new RooRealVar("bg_CB_p1", "bg_CB_p1", +1,   0,   10);
  RooRealVar* bg_CB_p2 = new RooRealVar("bg_CB_p2", "bg_CB_p2", +10,  0,   100);
  RooRealVar* bg_CB_p3 = new RooRealVar("bg_CB_p3", "bg_CB_p3", +500, 100, 1000);
  RooRealVar* bg_CB_p4 = new RooRealVar("bg_CB_p4", "bg_CB_p4", +100, 10,  500);
  for (auto mode: Modes) Fits_map[mode]["CB"] = new RevCrystalBall(mode+unique_name+"_CB", mode+unique_name+"_CB", *x_var, *ReDoVar(*bg_CB_p1,mode+unique_name), *ReDoVar(*bg_CB_p2,mode+unique_name), *ReDoVar(*bg_CB_p3,mode+unique_name), *ReDoVar(*bg_CB_p4,mode+unique_name));

  RooRealVar* bg_Novo_p1 = new RooRealVar("bg_Novo_p1", "bg_Novo_p1",  500, 200,  1000.);
  RooRealVar* bg_Novo_p2 = new RooRealVar("bg_Novo_p2", "bg_Novo_p2",  100,  10,  200.);
  RooRealVar* bg_Novo_p3 = new RooRealVar("bg_Novo_p3", "bg_Novo_p3", -0.6,  -2,  2. );
  for (auto mode: Modes) Fits_map[mode]["NO"] = new RooNovosibirsk(mode+unique_name+"_NO",mode+unique_name+"_NO",*x_var, *ReDoVar(*bg_Novo_p1,mode+unique_name), *ReDoVar(*bg_Novo_p2,mode+unique_name), *ReDoVar(*bg_Novo_p3,mode+unique_name));

  RooRealVar* bg_Exp_1_p1 = new RooRealVar("bg_Exp_1_p1", "bg_Exp_1_p1", -4., -100, 100);
  for (auto mode: Modes) Fits_map[mode]["Exp_1"] = new PolinomialExponent_1p(mode+unique_name+"_Exp_1",mode+unique_name+"_Exp_1",*x_var, *ReDoVar(*bg_Exp_1_p1,mode+unique_name));

  RooRealVar* bg_Exp_2_p1 = new RooRealVar("bg_Exp_2_p1", "bg_Exp_2_p1", -4., -100, 100);
  RooRealVar* bg_Exp_2_p2 = new RooRealVar("bg_Exp_2_p2", "bg_Exp_2_p2", 0.4, -100, 100);
  for (auto mode: Modes) Fits_map[mode]["Exp_2"] = new PolinomialExponent_2p(mode+unique_name+"_Exp_2",mode+unique_name+"_Exp_2",*x_var, *ReDoVar(*bg_Exp_2_p1,mode+unique_name), *ReDoVar(*bg_Exp_2_p2,mode+unique_name));

  RooRealVar* bg_Exp_3_p1 = new RooRealVar("bg_Exp_3_p1", "bg_Exp_3_p1", -4., -100, 100);
  RooRealVar* bg_Exp_3_p2 = new RooRealVar("bg_Exp_3_p2", "bg_Exp_3_p2", 0.4, -100, 100);
  RooRealVar* bg_Exp_3_p3 = new RooRealVar("bg_Exp_3_p3", "bg_Exp_3_p3", -0.1, -100, 100);
  for (auto mode: Modes) Fits_map[mode]["Exp_3"] = new PolinomialExponent_3p(mode+unique_name+"_Exp_3",mode+unique_name+"_Exp_3",*x_var, *ReDoVar(*bg_Exp_3_p1,mode+unique_name), *ReDoVar(*bg_Exp_3_p2,mode+unique_name), *ReDoVar(*bg_Exp_3_p3,mode+unique_name));

  RooRealVar* bg_Exp_4_p1 = new RooRealVar("bg_Exp_4_p1", "bg_Exp_4_p1", 4., -100, 100);
  RooRealVar* bg_Exp_4_p2 = new RooRealVar("bg_Exp_4_p2", "bg_Exp_4_p2", -10, -100,  10);
  RooRealVar* bg_Exp_4_p3 = new RooRealVar("bg_Exp_4_p3", "bg_Exp_4_p3", 4, -100, 100);
  RooRealVar* bg_Exp_4_p4 = new RooRealVar("bg_Exp_4_p4", "bg_Exp_4_p4", -0.4, -100, 100);
  for (auto mode: Modes) Fits_map[mode]["Exp_4"] = new PolinomialExponent_4p(mode+unique_name+"_Exp_4",mode+unique_name+"_Exp_4",*x_var, *ReDoVar(*bg_Exp_4_p1,mode+unique_name), *ReDoVar(*bg_Exp_4_p2,mode+unique_name), *ReDoVar(*bg_Exp_4_p3,mode+unique_name), *ReDoVar(*bg_Exp_4_p4,mode+unique_name));

  // RooRealVar* bg_Exp_5_p1 = new RooRealVar("bg_Exp_5_p1", "bg_Exp_5_p1", -8.46, -9.00, -8.00);
  // RooRealVar* bg_Exp_5_p2 = new RooRealVar("bg_Exp_5_p2", "bg_Exp_5_p2",  7.75,  7.00,  9.00);
  // RooRealVar* bg_Exp_5_p3 = new RooRealVar("bg_Exp_5_p3", "bg_Exp_5_p3", -7.45, -8.00, -6.00);
  // RooRealVar* bg_Exp_5_p4 = new RooRealVar("bg_Exp_5_p4", "bg_Exp_5_p4",  3.35,  3.00,  4.00);
  // RooRealVar* bg_Exp_5_p5 = new RooRealVar("bg_Exp_5_p5", "bg_Exp_5_p5", -5.37, -6.00, -4.00);
  // for (auto mode: Modes) Fits_map[mode]["Exp_5"] = new PolinomialExponent_5p(mode+unique_name+"_Exp_5",mode+unique_name+"_Exp_5",*x_var, *ReDoVar(*bg_Exp_5_p1,mode+unique_name), *ReDoVar(*bg_Exp_5_p2,mode+unique_name), *ReDoVar(*bg_Exp_5_p3,mode+unique_name), *ReDoVar(*bg_Exp_5_p4,mode+unique_name), *ReDoVar(*bg_Exp_5_p5,mode+unique_name));

  // RooRealVar* bg_Exp_6_p1 = new RooRealVar("bg_Exp_6_p1", "bg_Exp_6_p1", -3.68, -4.0, -3.0 );
  // RooRealVar* bg_Exp_6_p2 = new RooRealVar("bg_Exp_6_p2", "bg_Exp_6_p2", -1.57, -2.0, -1.0 );
  // RooRealVar* bg_Exp_6_p3 = new RooRealVar("bg_Exp_6_p3", "bg_Exp_6_p3", +1.56, +1.0, +2.0 );
  // RooRealVar* bg_Exp_6_p4 = new RooRealVar("bg_Exp_6_p4", "bg_Exp_6_p4", -1.14, -1.5, -1.0 );
  // RooRealVar* bg_Exp_6_p5 = new RooRealVar("bg_Exp_6_p5", "bg_Exp_6_p5", +5.29, +5.0, +5.5 );
  // RooRealVar* bg_Exp_6_p6 = new RooRealVar("bg_Exp_6_p6", "bg_Exp_6_p6", -7.05, -9.5, -9.0 );
  // for (auto mode: Modes) Fits_map[mode]["Exp_6"] = new PolinomialExponent_6p(mode+unique_name+"_Exp_6",mode+unique_name+"_Exp_6",*x_var, *ReDoVar(*bg_Exp_6_p1,mode+unique_name), *ReDoVar(*bg_Exp_6_p2,mode+unique_name), *ReDoVar(*bg_Exp_6_p3,mode+unique_name), *ReDoVar(*bg_Exp_6_p4,mode+unique_name), *ReDoVar(*bg_Exp_6_p5,mode+unique_name), *ReDoVar(*bg_Exp_6_p6,mode+unique_name));

  // RooRealVar* bg_BkgPdf4p_p1 = new RooRealVar("bg_BkgPdf4p_p1", "bg_BkgPdf4p_p1", +2.0, +2.0, +4.0);
  // RooRealVar* bg_BkgPdf4p_p2 = new RooRealVar("bg_BkgPdf4p_p2", "bg_BkgPdf4p_p2", +4.2, +3.0, +5.0);
  // RooRealVar* bg_BkgPdf4p_p3 = new RooRealVar("bg_BkgPdf4p_p3", "bg_BkgPdf4p_p3", +1.4, +0.5, +2.0);
  // RooRealVar* bg_BkgPdf4p_p4 = new RooRealVar("bg_BkgPdf4p_p4", "bg_BkgPdf4p_p4", -1.0, -2.0, -0.5);
  // for (auto mode: Modes) Fits_map[mode]["BkgPdf4p"] = new BkgPdf4p(mode+unique_name+"_BkgPdf4p",mode+unique_name+"_BkgPdf4p",*x_var, *ReDoVar(*bg_BkgPdf4p_p1,mode+unique_name), *ReDoVar(*bg_BkgPdf4p_p2,mode+unique_name), *ReDoVar(*bg_BkgPdf4p_p3,mode+unique_name), *ReDoVar(*bg_BkgPdf4p_p4,mode+unique_name));

  // RooRealVar* bg_BkgPdf3p_p1 = new RooRealVar("bg_BkgPdf3p_p1", "bg_BkgPdf3p_p1", -10., -11., -9.0);
  // RooRealVar* bg_BkgPdf3p_p2 = new RooRealVar("bg_BkgPdf3p_p2", "bg_BkgPdf3p_p2", +7.0, +6.0, +8.0);
  // RooRealVar* bg_BkgPdf3p_p3 = new RooRealVar("bg_BkgPdf3p_p3", "bg_BkgPdf3p_p3", +3.7, +3.0, +4.0);

  // RooRealVar* bg_BkgPdf3p_p1 = new RooRealVar("bg_BkgPdf3p_p1", "bg_BkgPdf3p_p1", -13.2, -1000, 1000);
  // RooRealVar* bg_BkgPdf3p_p2 = new RooRealVar("bg_BkgPdf3p_p2", "bg_BkgPdf3p_p2", +9.1, -1000, 1000);
  // RooRealVar* bg_BkgPdf3p_p3 = new RooRealVar("bg_BkgPdf3p_p3", "bg_BkgPdf3p_p3", +2.5, -100, 100);

  // RooRealVar* bg_BkgPdf3p_p1 = new RooRealVar("bg_BkgPdf3p_p1", "bg_BkgPdf3p_p1", -4.3, -5.3, -3.3);
  // RooRealVar* bg_BkgPdf3p_p2 = new RooRealVar("bg_BkgPdf3p_p2", "bg_BkgPdf3p_p2", +5.7, +4.7, +7.7);
  // RooRealVar* bg_BkgPdf3p_p3 = new RooRealVar("bg_BkgPdf3p_p3", "bg_BkgPdf3p_p3", +2.5, +1.5, +3.5);
  // for (auto mode: Modes) Fits_map[mode]["BkgPdf3p"] = new BkgPdf3p(mode+unique_name+"_BkgPdf3p",mode+unique_name+"_BkgPdf3p",*x_var, *ReDoVar(*bg_BkgPdf3p_p1,mode+unique_name), *ReDoVar(*bg_BkgPdf3p_p2,mode+unique_name), *ReDoVar(*bg_BkgPdf3p_p3,mode+unique_name));

  std::ofstream output;
  if (doFtest) output.open(workingDir+"datacards/output_"+year+"_"+histFolder+".txt");

  if (doBkgPlots) {
    for (TString mode: Modes) {
      // if (mode!="bkg_pred" && mode!="DY_CR" && mode!="DY_SR") continue;
      // if (mode!="bkg_pred" && mode!="DY_CR") continue;
      if (mode!="bkg_pred") continue;
      std::map<TString, RooHist*> hpull;
      std::map<TString, RooHist*> hratio;
      std::map<TString, double> chi2_map;

      TCanvas* c_bg = tdrDiCanvas("Events"+mode, plot_lo, plot_hi, plot_ylo, plot_yhi, doPlotRatio?0.8:-6, doPlotRatio?1.2:6, nameXaxis, nameYaxis, nameRatioaxis);
      c_bg->cd(1)->SetLogy(1);
      // plotter = x_var->frame(plot_lo,plot_hi);
      plotter.reset(x_var->frame(plot_lo,plot_hi));

      rooHist_map[mode]->plotOn(plotter.get(),RooFit::DataError(RooAbsData::Poisson));
      // rooHist_map[mode]->plotOn(plotter.get(),RooFit::DataError(RooAbsData::SumW2));
      // rooHist_map[mode]->plotOn(plotter.get());

      for (auto const& [model,dofit] : doFits_map[mode] ) {
        if (dofit) {
          int color = Colors[model];
          std::cout << "PDF FIT " << mode << "\t" << model << '\n';
          FitRes_map[mode][model] = Fits_map[mode][model]->fitTo(*rooHist_map[mode], RooFit::Range(fit_min[mode], fit_max[mode]), RooFit::SumW2Error(kTRUE), RooFit::Save(), RooFit::Verbose(kFALSE), RooFit::PrintEvalErrors(-1));
          // if (!doFtest && model.Contains("Exp") && model!="Exp_2" && model!="Exp_3" && model!="Exp_4" ) continue;
          Fits_map[mode][model]->plotOn(plotter.get(), RooFit::LineColor(color), RooFit::Range(pull_min, pull_max, kFALSE));
          // Fits_map[mode][model]->plotOn(plotter, RooFit::LineColor(color), RooFit::Range(fit_min[mode], fit_max[mode], kFALSE));

          // CREATE RATIO BETWEEN CURVE AND HIST
          double xstart,xstop,y ;
          RooCurve* curve = (RooCurve*) plotter->findObject(0,RooCurve::Class()) ;
          RooHist* hist_ = (RooHist*) plotter->findObject(0,RooHist::Class()) ;
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
          TString name = model+std::string(7-model.Length(), ' ' )+" #chi^2/n.d.f. = "+TString::Format("%.1f",chi2_map[model])+" p-value = "+TString::Format("%.2f",pv)+"\%";
          // TString name = model+std::string(5-model.Length(), ' ' )+"np="+TString::Format("%d",npf);
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

      TLegend *leg_bg = tdrLeg(0.40,0.55,0.7,0.85, 0.032);
      tdrHeader(leg_bg, histo_map[mode]->GetName(), 12, 0.04, 42, kBlack, true);


      RooPlot* Pullplotter = x_var->frame(plot_lo,plot_hi);
      for (auto x : doPlotRatio?hratio:hpull) {
        TString model = x.first;
        if (model.Contains("Exp_") || model.Contains("BkgPdf")) {
          int deg = atoi(new char(model.Contains("Exp_")?model[4]: model[6]));
          TString old_model = x.first;
          old_model.ReplaceAll(std::to_string(deg),std::to_string(deg-1));
          if (chi2_map.find(old_model)!= chi2_map.end()) {
            double FTest = DoFTest(chi2_map[old_model],chi2_map[model],deg-1,deg,hpull[model]->GetN());
            if (doFtest) output << mode << " " << model << " FTest " << FTest << " chi2: " << chi2_map[model] << '\n';
            x.second->SetName(x.second->GetName()+TString::Format(" F-test = %.2f",FTest));
          }
        }
        leg_bg->AddEntry(x.second, x.second->GetName() ,"l");
        Pullplotter->addPlotable(x.second,"P same");
      }
      leg_bg->Draw("same");

      c_bg->cd(2);
      Pullplotter->SetNdivisions(505,"Y");
      Pullplotter->Draw("same");
      TLine *line=new TLine(plot_lo, 0, plot_hi, 0);
      line->SetLineWidth(2);
      line->Draw("same");

      c_bg->SaveAs(workingDir+"Fit_Bg_all_"+mode+"_"+histFolder+extra_text+"."+plotting_mode);
      if (plotting_mode_2!="") c_bg->SaveAs(workingDir+"Fit_Bg_all_"+mode+"_"+histFolder+extra_text+"."+plotting_mode_2);


      if (doCheckPlots) {
        for (auto [model,pull] : hpull) {
          TCanvas* c_gauss = tdrCanvas("gaus_"+mode+"_"+model, -6, 6, 0.00001,100, nameXaxis, nameYaxis);
          TH1F* gauss = new TH1F(model, model,13, -6,6);
          TF1 *f1 = new TF1("f1","gaus",-5,5);
          for (int i = 0; i < pull->GetN(); i++) {
            double x,y;
            pull->GetPoint(i,x,y);
            gauss->Fill(y);
          }
          gauss->Fit(f1,"RMQ");
          tdrDraw(gauss, "Hist", kSolid, kBlack, kSolid, kBlack, 3000, kBlack);
          TLegend *leg_gauss = tdrLeg(0.70,0.60,0.9,0.85, 0.035);
          tdrHeader(leg_gauss, mode+" "+model, 12, 0.04, 42, kBlack, true);
          c_gauss->SaveAs(workingDir+"Gauss_Bg_all_"+mode+"_"+model+"_"+histFolder+extra_text+"."+plotting_mode);
          if (plotting_mode_2!="") c_gauss->SaveAs(workingDir+"Gauss_Bg_all_"+mode+"_"+model+"_"+histFolder+extra_text+"."+plotting_mode_2);
        }
      }


    }
  }

  if (doFtest) output.close();

  /*
  # ######## ########     ###    ##    ##  ######  ######## ######## ########     ######## ##     ## ##    ##  ######  ######## ####  #######  ##    ##  ######
  #    ##    ##     ##   ## ##   ###   ## ##    ## ##       ##       ##     ##    ##       ##     ## ###   ## ##    ##    ##     ##  ##     ## ###   ## ##    ##
  #    ##    ##     ##  ##   ##  ####  ## ##       ##       ##       ##     ##    ##       ##     ## ####  ## ##          ##     ##  ##     ## ####  ## ##
  #    ##    ########  ##     ## ## ## ##  ######  ######   ######   ########     ######   ##     ## ## ## ## ##          ##     ##  ##     ## ## ## ##  ######
  #    ##    ##   ##   ######### ##  ####       ## ##       ##       ##   ##      ##       ##     ## ##  #### ##          ##     ##  ##     ## ##  ####       ##
  #    ##    ##    ##  ##     ## ##   ### ##    ## ##       ##       ##    ##     ##       ##     ## ##   ### ##    ##    ##     ##  ##     ## ##   ### ##    ##
  #    ##    ##     ## ##     ## ##    ##  ######  ##       ######## ##     ##    ##        #######  ##    ##  ######     ##    ####  #######  ##    ##  ######
  */

  TString mode_CR = BkgName+"_CR";
  TString mode_SR = BkgName+"_SR";
  TCanvas* c_TF = tdrCanvas("TF", plot_lo, plot_hi, 1e-03, 1e00, nameXaxis, "A.U");
  c_TF->SetLogy(1);
  TLegend *leg_TF = tdrLeg(0.40,0.60,0.7,0.85, 0.035);
  tdrHeader(leg_TF, "Transfer Fuctions", 12, 0.04, 42, kBlack, true);
  gStyle->SetOptFit(kFALSE);

  TGraphAsymmErrors* ratio_hist = new TGraphAsymmErrors();
  ratio_hist->Divide(histo_map[mode_SR], histo_map[mode_CR], "pois");
  TF1 *f1 = new TF1("TF","pol0",fit_min[mode_CR],fit_max[mode_CR]);
  ratio_hist->Fit(f1,"RQ");
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
  tdrDraw(ratio_hist, "P0", kFullDotLarge, kBlack, kSolid, kBlack, 3004, kBlack);
  leg_TF->AddEntry(ratio_hist, Form("%s/%s", histo_map[mode_SR]->GetName(),histo_map[mode_CR]->GetName()) ,"lep");
  c_TF->SaveAs(workingDir+"TFs_"+histFolder+extra_text+"."+plotting_mode);
  if (plotting_mode_2!="") c_TF->SaveAs(workingDir+"TFs_"+histFolder+extra_text+"."+plotting_mode_2);




  /*
  # ########   #######   ######  ##     ## ########  ######  ##    ## ########  ##        #######  ########  ######
  # ##     ## ##     ## ##    ## ##     ## ##       ##    ## ##   ##  ##     ## ##       ##     ##    ##    ##    ##
  # ##     ## ##     ## ##       ##     ## ##       ##       ##  ##   ##     ## ##       ##     ##    ##    ##
  # ##     ## ##     ## ##       ######### ######   ##       #####    ########  ##       ##     ##    ##     ######
  # ##     ## ##     ## ##       ##     ## ##       ##       ##  ##   ##        ##       ##     ##    ##          ##
  # ##     ## ##     ## ##    ## ##     ## ##       ##    ## ##   ##  ##        ##       ##     ##    ##    ##    ##
  # ########   #######   ######  ##     ## ########  ######  ##    ## ##        ########  #######     ##     ######
  */



  if (doCheckPlots) {
    TCanvas* c_data_obs = tdrCanvas("data_obs", plot_lo, plot_hi, plot_ylo, 1e03, nameXaxis, nameYaxis);
    c_data_obs->SetLogy(1);
    // plotter = x_var->frame(plot_lo,plot_hi);
    plotter.reset(x_var->frame(plot_lo,plot_hi));
    data_obs->plotOn(plotter.get(), RooFit::LineColor(kBlack), RooFit::MarkerColor(kBlack));
    plotter->Draw("same");
    TLegend *leg_data_obs = tdrLeg(0.50,0.70,0.9,0.9, 0.035);
    leg_data_obs->AddEntry(data_obs,  Form("%s: %s", data_obs->GetName(), data_obs->GetTitle()) ,"lp");
    leg_data_obs->Draw("same");

    c_data_obs->SaveAs((workingDir+"Plot_data_obs_"+histFolder+"."+plotting_mode).c_str());
    if (plotting_mode_2!="") c_data_obs->SaveAs((workingDir+"Plot_data_obs_"+histFolder+"."+plotting_mode_2).c_str());

    // TODO Do Transfer functions
    TCanvas* c_inputs = tdrCanvas("inputs DataCard", plot_lo, plot_hi, plot_ylo, plot_yhi, nameXaxis, nameYaxis);
    c_inputs->SetLogy(1);
    // histo_map["bkg_pred"]->SetLineColor(kBlue+1);
    // histo_map["data"]->SetLineColor(kRed+1);
    // histo_map["main_bkg_SR"]->SetLineColor(kGreen+1);
    histo_map["bkg_pred"]->Draw("same");
    histo_map["data"]->Draw("same");
    histo_map["main_bkg_SR"]->Draw("same");
    TLegend *leg_CRvsSR = tdrLeg(0.50,0.70,0.78,0.9, 0.03);
    leg_CRvsSR->AddEntry(histo_map["bkg_pred"], Form("dofit: %s", histo_map["bkg_pred"]->GetName()) ,"lp");
    leg_CRvsSR->AddEntry(histo_map["data"], Form("extract obs. Limits: %s", histo_map["data"]->GetName()) ,"lp");
    leg_CRvsSR->AddEntry(histo_map["main_bkg_SR"], Form("extract exp. Limits: %s", histo_map["main_bkg_SR"]->GetName()) ,"lp");
    leg_CRvsSR->Draw("same");
    c_inputs->SaveAs((workingDir+"Plot_Inputs_"+histFolder+"."+plotting_mode).c_str());
    if (plotting_mode_2!="") c_inputs->SaveAs((workingDir+"Plot_Inputs_"+histFolder+"."+plotting_mode_2).c_str());

  }

  /*
  #  ######  ########  #######  ########  ########    ##      ##  #######  ########  ##    ##  ######  ########     ###     ######  ########
  # ##    ##    ##    ##     ## ##     ## ##          ##  ##  ## ##     ## ##     ## ##   ##  ##    ## ##     ##   ## ##   ##    ## ##
  # ##          ##    ##     ## ##     ## ##          ##  ##  ## ##     ## ##     ## ##  ##   ##       ##     ##  ##   ##  ##       ##
  #  ######     ##    ##     ## ########  ######      ##  ##  ## ##     ## ########  #####     ######  ########  ##     ## ##       ######
  #       ##    ##    ##     ## ##   ##   ##          ##  ##  ## ##     ## ##   ##   ##  ##         ## ##        ######### ##       ##
  # ##    ##    ##    ##     ## ##    ##  ##          ##  ##  ## ##     ## ##    ##  ##   ##  ##    ## ##        ##     ## ##    ## ##
  #  ######     ##     #######  ##     ## ########     ###  ###   #######  ##     ## ##    ##  ######  ##        ##     ##  ######  ########
  */


  // ws->import(*(new RooHistPdf(SgName.c_str(),"",*x_var,*Sg_Hist))); //TODO decide weather to go for shaped analysis or not

  for (auto mode: Modes) {
    for (auto const& [model,dofit] : doFits_map[mode] ) {
      if (!dofit || (model.Contains("Exp") && model!="Exp_2" && model!="Exp_3" && model!="Exp_4") ) continue;
      ws->import(*Fits_map[mode][model], RooFit::Silence());
      x_var->setRange("fitting", fit_lo, fit_hi);
      double i_fit = ((RooAbsReal*)((RooAbsPdf*)Fits_map[mode][model])->createIntegral(*x_var, RooFit::NormSet(*x_var), RooFit::Range("fitting")))->getVal();
      x_var->setRange("total", x_lo, x_hi);
      double i_tot = ((RooAbsReal*)((RooAbsPdf*)Fits_map[mode][model])->createIntegral(*x_var, RooFit::NormSet(*x_var), RooFit::Range("total")))->getVal();
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

      if (debug) std::cout << mode+unique_name+"_"+model << " integral " << nEventsSR*i_tot/i_fit << " " << CalculateIntegral(histo_map[mode],fit_lo,fit_hi,doBinWidth)*i_tot/i_fit << " " << CalculateIntegral(histo_map[mode],fit_lo,fit_hi,doBinWidth)/CalculateFractionArea(histo_map[mode],fit_lo,fit_hi, x_lo, x_hi,doBinWidth) << " " << rooHist_map[mode]->sum(!doBinWidth) << std::endl;
      DataCard  << mode+unique_name+"_"+model << " integral " << nEventsSR*i_tot/i_fit << " " << CalculateIntegral(histo_map[mode],fit_lo,fit_hi,doBinWidth)*i_tot/i_fit << " " << CalculateIntegral(histo_map[mode],fit_lo,fit_hi,doBinWidth)/CalculateFractionArea(histo_map[mode],fit_lo,fit_hi, x_lo, x_hi,doBinWidth) << " " << rooHist_map[mode]->sum(!doBinWidth) << std::endl;


      RooArgSet* model_params = Fits_map[mode][model]->getParameters(*x_var) ;
      TString name; int from = 0;
      while (TString(model_params->contentsString()).Tokenize(name, from, ",")) {
        RooRealVar* par = dynamic_cast<RooRealVar*>(model_params->find(name));
        DataCard << Form("%s%sparam %.3f %.3f", name.Data(), std::string(60-name.Length(),' ').c_str(), par->getVal(), par->getError()) <<std::endl;
      }
    }
  }

  /*
  #  ######  ####  ######   ##    ##    ###    ##          ######## #### ########
  # ##    ##  ##  ##    ##  ###   ##   ## ##   ##          ##        ##     ##
  # ##        ##  ##        ####  ##  ##   ##  ##          ##        ##     ##
  #  ######   ##  ##   #### ## ## ## ##     ## ##          ######    ##     ##
  #       ##  ##  ##    ##  ##  #### ######### ##          ##        ##     ##
  # ##    ##  ##  ##    ##  ##   ### ##     ## ##          ##        ##     ##
  #  ######  ####  ######   ##    ## ##     ## ########    ##       ####    ##
  */




  std::cout << "****************" << '\n';
  std::cout << "*  SIGNAL FIT  *" << '\n';
  std::cout << "****************" << '\n';

  std::ofstream SignalProperties;
  SignalProperties.open(workingDir+"datacards/SignalProperties_"+year+"_"+histFolder+".txt");

  // std::vector<int> masses = {600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2500};
  std::vector<int> masses = {600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500, 6000, 7000, 8000}; // TODO take if from constants.hpp
  // std::vector<int> masses = {1000};
  // std::vector<int> masses = {3000, 3500, 4000, 4500, 5000};
  std::map<int, int> colors = {{600, kRed}, {800, kGreen}, {1000, kViolet}, {1200, kBlue}, {1400, kBlack}, {1600, kOrange}, {1800, kAzure}, {2000, kSpring}, {2500, kPink},
  {3000, kRed}, {3500, kGreen}, {4000, kViolet}, {4500, kBlue}, {5000, kBlack}, {5500, kOrange}, {6000, kAzure}, {7000, kSpring}, {8000, kPink}};


  // double xsec_ref_ = 0.1; // this is to mantain the signal strenght close to 1; (remember to multiply for this Normalization when plotting)
  double xsec_ref_ = (xsec_ref.find(histFolder) != xsec_ref.end())? xsec_ref.at(histFolder): xsec_ref.at("default_value"); // default_value = 1
  DataCard << "xsec_ref " << " " << xsec_ref_ <<std::endl;


  RooPlot* Pullplotter_all = x_var->frame(x_lo,x_hi);
  RooPlot* plotter_all = x_var->frame(x_lo,x_hi);
  // TCanvas* c_sg_all = tdrDiCanvas("signal", plot_lo, x_hi, 2*1e-03, 40, -6, 6, nameXaxis, nameYaxis, "Pull");
  TCanvas* c_sg_all = tdrDiCanvas("signal", plot_lo, x_hi, 1*1e-03, 1, -6, 6, nameXaxis, nameYaxis, "Pull");
  TPaveText* pave;

  std::map<std::string, std::vector<double>> parameters;
  for (std::string name : {"Masses","nevents", "nevents_err", "mean","sigma","alpha","k","mean_err","sigma_err","alpha_err","k_err","chi2","pvalue"}) parameters[name] = std::vector<double>(masses.size(), 0);
  int i_mass = -1;

  for (int & mass : masses) {
    i_mass++;
    // if (mass!=1000) continue;
    parameters["Masses"].at(i_mass) = (double)mass;
    double rangeLo = mass*(1-5./28.);
    double rangeHi = mass*(1+3.5/28.);
    // if (mass<=1000) {
    //   rangeLo = mass*(1-4./28.);
    // }
    double plotLo = mass*(1-10./28.);
    double plotHi = mass*(1+10./28.);
    double ymax = (50-mass*1./100)*1.5;
    if (mass<1200) ymax = (-20+mass*6./100)*1.4;
    if (mass>3500) ymax = (30-mass*3./1000)*1.3;

    if (mass==8000 && channel=="electronchannel") ymax /= 2.;
    if (mass==8000 && year!="RunII") rangeLo -= 100.;
    if (mass==8000) rangeHi -= 100.;
    if (mass==6000) rangeLo -= 500.;
    if (mass==6000) rangeHi -= 200.;
    if (mass==5500) rangeLo -= 400.;
    if (mass==5500) rangeHi -= 400.;
    if (mass==5000) rangeLo -= 300.;
    if (mass==4500) rangeLo -= 500.;
    if (mass==4000) rangeLo -= 300.;
    if (mass==4000) rangeHi -= 100.;
    if (mass==3500) rangeLo -= 200.;
    if (mass==3500) rangeHi -= 100.;
    // if (mass==3000) rangeLo -= 200.;
    if (mass==3000) rangeHi -= 100.;
    if (mass==2000) rangeLo -= 100.;
    if (mass==2000) rangeHi -= 50.;
    if (mass==1800) rangeLo -= 100.;
    if (mass==1600) rangeLo -= 100.;
    if (mass==1400) rangeLo -= 100.;
    if (mass==1200) rangeLo -= 100.;
    if (mass==1000) rangeLo -= 50.;
    if (mass==1000) rangeHi += 50.;
    if (mass==800) rangeLo -= 50.;
    if (mass==800) rangeHi += 50.;
    if (mass==600) rangeLo -= 50.;
    if (mass==600) rangeHi += 50.;

    ymax *= lumi_map.at(year).at("lumi_fb")/lumi_map.at("RunII").at("lumi_fb");
    if (year=="2017") ymax *= 1.1;

    if (year=="2018" && mass==7000) rangeLo -= 300;
    if (year=="2018" && mass==7000) rangeHi -= 300;

    ymax *= xsec_ref_/0.1;

    if (debug) std::cout << "fit mass: " << mass << " between " << rangeLo << " and " << rangeHi << " plotting between " <<  plotLo << " and " << plotHi << " ymax " << ymax << '\n';

    std::string SgName = "M"+std::to_string(mass);
    std::string fnameSignal = "MC_ZprimeToZH_"+SgName;
    // std::string fnameSignal = "MC_ZprimeToZHToWW_"+SgName;
    // if (isHbb) fnameSignal = "MC_ZprimeToZHTobb_"+SgName;

    TFile *f_signal=new TFile((filepath+PrefixrootFile+"MC."+fnameSignal+"_"+year+"_noTree.root").c_str());
    TH1F *h_signal=(TH1F*)f_signal->Get(SGname.c_str());

    if (dorebin) {
      if (rebin) h_signal->Rebin(rebin);
      else h_signal = dynamic_cast<TH1F*>(h_signal->Rebin(bins_Zprime_rebin.size()-1, h_signal->GetName(), &bins_Zprime_rebin[0]));
      h_signal->Scale(1,doBinWidth?"width":"");
    }
    h_signal->Scale(xsec_ref_);

    RooRealVar* sg_p0 =new RooRealVar(("sg"+unique_name+"_p0").c_str(), ("sg"+unique_name+"_p0").c_str(), mass, rangeLo, rangeHi);
    RooRealVar* sg_p1 =new RooRealVar(("sg"+unique_name+"_p1").c_str(), ("sg"+unique_name+"_p1").c_str(), (mass>=3000)? 150: 40, (mass>=3000)? 120: 10., (mass>=3000)? 400:110.);

    RooRealVar* sg_p2 =new RooRealVar(("sg"+unique_name+"_p2").c_str(), ("sg"+unique_name+"_p2").c_str(), 1, -10, 10);
    RooRealVar* sg_p3 =new RooRealVar(("sg"+unique_name+"_p3").c_str(), ("sg"+unique_name+"_p3").c_str(), 0.1, -10, 10);


    // bool isCB = false;
    // RooGaussian* ModelSg = new RooGaussian(SgName.c_str(), "Signal Prediction", *x_var, *sg_p0, *sg_p1);
    bool isCB = true;
    RooCBShape* ModelSg = new RooCBShape(SgName.c_str(), "Signal Prediction", *x_var, *sg_p0, *sg_p1, *sg_p2, *sg_p3);
    RooDataHist* Sg_Hist;// = new RooDataHist(SgName.c_str(), SgName.c_str(), RooArgList(*x_var), h_signal);

    if (dorebin && rebin!=0) Sg_Hist = new RooDataHist(SgName.c_str(), SgName.c_str(), RooArgList(*x_var), RooFit::Import(*h_signal,doBinWidth));
    else Sg_Hist = new RooDataHist(SgName.c_str(), SgName.c_str(), RooArgList(*x_var), h_signal);

    RooFitResult *r_sg=ModelSg->fitTo(*Sg_Hist, RooFit::Range(GetRange(h_signal,rangeLo), GetRange(h_signal,rangeHi)), RooFit::SumW2Error(kTRUE), RooFit::Save(), RooFit::Verbose(kFALSE), RooFit::PrintEvalErrors(-1));

    if (doSignalPlots) {
      RooPlot* Pullplotter = x_var->frame(plotLo,plotHi);

      plotter.reset(x_var->frame(plotLo,plotHi));
      TCanvas* c_sg = tdrDiCanvas(SgName.c_str(), plotLo,plotHi, 2*1e-03, ymax, -6, 6, nameXaxis, nameYaxis, "Pull");
      // c_sg->cd(1)->SetLogy(1);
      Sg_Hist->plotOn(plotter.get());
      ModelSg->plotOn(plotter.get(), RooFit::VisualizeError(*r_sg, 1), RooFit::FillColor(kRed-7), RooFit::FillStyle(3001));
      ModelSg->plotOn(plotter.get(), RooFit::VisualizeError(*r_sg, 2), RooFit::FillColor(kRed-9), RooFit::FillStyle(3001));
      ModelSg->plotOn(plotter.get(), RooFit::LineColor(kRed));
      Sg_Hist->plotOn(plotter.get(), RooFit::LineColor(kBlack), RooFit::MarkerColor(kBlack));
      RooHist* hpull = plotter->pullHist();
      int npf = r_sg->floatParsFinal().getSize();
      // int ndf = hpull->GetN()-npf;
      int ndf;
      double chi2;
      CalculateChiSquare(chi2, ndf, hpull, rangeLo,rangeHi); //TODO do the same for bkg plotter->chiSquare(npf);
      ndf -=npf;
      double pv = TMath::Prob(chi2*ndf, ndf)*100;
      SignalProperties << SgName+" signal number of events = " << CalculateIntegral(h_signal,rangeLo-100,rangeHi+100,doBinWidth) <<""<<std::endl;//TODO
      SignalProperties << " chi2-pvalue "  << chi2 << " " << pv <<std::endl;
      SignalProperties << sg_p0->GetName() << "   param   "<<sg_p0->getVal()<<" "<<sg_p0->getError()<<std::endl;
      SignalProperties << sg_p1->GetName() << "   param   "<<sg_p1->getVal()<<" "<<sg_p1->getError()<<std::endl;
      SignalProperties << sg_p2->GetName() << "   param   "<<sg_p2->getVal()<<" "<<sg_p2->getError()<<std::endl;
      SignalProperties << sg_p3->GetName() << "   param   "<<sg_p3->getVal()<<" "<<sg_p3->getError()<<std::endl;

      pave = new TPaveText(0.7,0.6,0.8,0.8,"NDC");
      pave->SetBorderSize(0); pave->SetTextSize(0.03); pave->SetLineColor(1); pave->SetLineStyle(1);
      pave->SetLineWidth(2); pave->SetFillColor(0); pave->SetFillStyle(0);
      pave->AddText(TString::Format("Fit range = [%.0f,%.0f]", GetRange(h_signal,rangeLo),GetRange(h_signal,rangeHi)));
      pave->AddText(TString::Format("#chi2 = %.1f", chi2));
      pave->AddText(TString::Format("p-value = %.2f", pv));
      pave->AddText(TString::Format("#mu = %2.3f +- %2.3f", sg_p0->getVal(),sg_p0->getError()));
      pave->AddText(TString::Format("#sigma = %2.3f +- %2.3f", sg_p1->getVal(),sg_p1->getError()));
      pave->AddText(TString::Format("#alpha = %2.3f +- %2.3f", sg_p2->getVal(),sg_p2->getError()));
      pave->AddText(TString::Format("k = %2.3f +- %2.3f", sg_p3->getVal(),sg_p3->getError()));
      // parameters["nevents"].at(i_mass)     = CalculateIntegral(h_signal,rangeLo-100,rangeHi+100,doBinWidth);//TODO
      parameters["nevents"].at(i_mass)     = CalculateIntegral(h_signal,x_lo,x_hi,doBinWidth);//TODO what range?
      parameters["nevents_err"].at(i_mass) = TMath::Sqrt(parameters["nevents"].at(i_mass));
      parameters["mean"].at(i_mass)        = sg_p0->getVal();
      parameters["mean_err"].at(i_mass)    = sg_p0->getError();
      parameters["sigma"].at(i_mass)       = sg_p1->getVal();
      parameters["sigma_err"].at(i_mass)   = sg_p1->getError();
      parameters["alpha"].at(i_mass)       = sg_p2->getVal();
      parameters["alpha_err"].at(i_mass)   = sg_p2->getError();
      parameters["k"].at(i_mass)           = sg_p3->getVal();
      parameters["k_err"].at(i_mass)       = sg_p3->getError();
      parameters["chi2"].at(i_mass)        = chi2;
      parameters["pvalue"].at(i_mass)      = pv;

      ModelSg->plotOn(plotter_all, RooFit::LineColor(kRed));
      c_sg->cd(1);
      plotter->Draw("same");
      h_signal->Draw("same");
      pave->Draw("same");
      c_sg->cd(2);
      Pullplotter->addPlotable(hpull,"P same");
      Pullplotter->SetNdivisions(505,"Y");
      Pullplotter->Draw("same");

      TLine *line=new TLine(plotLo, 0, plotHi, 0);
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



      // file_WS.cd();
      // h_signal->SetNameTitle(("H"+std::to_string(mass)).c_str(),("H"+std::to_string(mass)).c_str());
      // h_signal->Write();
      ws->import(*ModelSg, RooFit::Silence());
    }

    // file_WS.cd();
    // for (auto [mode, h]: histo_map) {
    //   h->SetNameTitle(mode,mode);
    //   h->Write();
    //   DataCard << mode+" integral = " << h->Integral()<<""<<std::endl;
    // }
    //
    // data_obs->Write();
    // file_WS.Close();

    //TODO fix the norm of signal
    // DataCard << SgName+" signal number of events = " << h_signal->GetSumOfWeights()<<""<<std::endl;
    DataCard << SgName+" signal number of events = " << CalculateIntegral(h_signal,rangeLo-100,rangeHi+100,doBinWidth) <<""<<std::endl;
    // DataCard << SgName+" signal number of events = " << h_signal->Integral(h_signal->FindBin(rangeLo-100),h_signal->FindBin(rangeHi+100), doBinWidth?"width":"")<<""<<std::endl;

    DataCard << SgName << " " << sg_p0->GetName() << "   param   "<<sg_p0->getVal()<<" "<<sg_p0->getError()<<std::endl;
    DataCard << SgName << " " << sg_p1->GetName() << "   param   "<<sg_p1->getVal()<<" "<<sg_p1->getError()<<std::endl;
    if (isCB) DataCard << SgName << " " << sg_p2->GetName() << "   param   "<<sg_p2->getVal()<<" "<<sg_p2->getError()<<std::endl;
    if (isCB) DataCard << SgName << " " << sg_p3->GetName() << "   param   "<<sg_p3->getVal()<<" "<<sg_p3->getError()<<std::endl;

  }

  if (doSignalPlots) {
    c_sg_all->SaveAs((workingDir+"Fit_Sg_all_"+histFolder+"."+plotting_mode).c_str());
    if (plotting_mode_2!="") c_sg_all->SaveAs((workingDir+"Fit_Sg_all_"+histFolder+"."+plotting_mode_2).c_str());

    gStyle->SetOptFit(kFALSE);
    TCanvas* c_sg_par;
    TGraphErrors* gr_par;
    TF1* fit_par;
    TSpline3* fit_par_spline;
    TPaveText* pave_par;
    std::vector<double> dummy(parameters["Masses"].size(),50);

    for (auto x: parameters) {
      if (x.first=="Masses" || x.first.find("_err") != std::string::npos) continue;
      double y_max = (x.first=="mean")? 9000: ((x.first=="sigma")? 500:((x.first=="pvalue")? 100:((x.first=="nevents")? 400*lumi_map.at(year).at("lumi_fb")/lumi_map.at("RunII").at("lumi_fb"):5)));
      c_sg_par = tdrCanvas(("c_sg_"+x.first).c_str(), plot_lo, x_hi, 2*1e-03, y_max, nameXaxis, (x.first).c_str());
      if (x.first=="chi2" || x.first=="pvalue") gr_par = new TGraphErrors(parameters["Masses"].size(), &(parameters["Masses"][0]), &(x.second[0]), &(dummy[0]), &(dummy[0]));
      else gr_par = new TGraphErrors(parameters["Masses"].size(), &(parameters["Masses"][0]), &(x.second[0]), &(dummy[0]), &(parameters[x.first+"_err"][0]));
      fit_par = new TF1((x.first).c_str(),(x.first=="nevents")? "pol3":((x.first=="sigma")?"pol2":"pol1"),500, 8100);
      fit_par_spline = new TSpline3(("spline_"+x.first).c_str(), &(parameters["Masses"][0]), &(x.second[0]), parameters["Masses"].size(), "b2e2");
      fit_par_spline->SetLineColor(kOrange+1);
      fit_par_spline->Draw("same");
      // if (x.first=="nevents") fit_par->SetParameters();
      fit_par->SetLineColor(kRed+1); fit_par->SetLineWidth(3);
      gr_par->Fit(fit_par,"RSQ");
      // pave_par = new TPaveText(0.55,0.7,0.82,0.9,"NDC");
      pave_par = new TPaveText(0.7,0.6,0.8,0.8,"NDC");
      pave_par->SetBorderSize(0); pave_par->SetTextSize(0.03); pave_par->SetLineColor(1); pave_par->SetLineStyle(1);
      pave_par->SetLineWidth(2); pave_par->SetFillColor(0); pave_par->SetFillStyle(0);
      for (int i = 0; i < fit_par->GetNumberFreeParameters(); i++) {
        if (fit_par->GetParameter(i)<1e-03) pave_par->AddText(TString::Format("p%d = (%2.3f +- %2.3f)*10^{-3}", i,fit_par->GetParameter(i)*1e03,fit_par->GetParError(i)*1e03));
        else pave_par->AddText(TString::Format("p%d = %2.3f +- %2.3f", i,fit_par->GetParameter(i),fit_par->GetParError(i)));
      }
      pave_par->AddText(TString::Format("#chi2/n.d.f. = %.1f", fit_par->GetChisquare()/fit_par->GetNDF()));
      pave_par->AddText(TString::Format("p-value = %.2f", fit_par->GetProb()));
      gStyle->SetOptFit(kFALSE);
      tdrDraw(gr_par, "P", kFullDotLarge, kRed+1, kSolid, kRed+1, 3004, kRed+1);
      pave_par->Draw("same");
      c_sg_par->SaveAs((workingDir+"Fit_Sg_"+x.first+"_"+histFolder+"."+plotting_mode).c_str());
      if (plotting_mode_2!="") c_sg_par->SaveAs((workingDir+"Fit_Sg_"+x.first+"_"+histFolder+"."+plotting_mode_2).c_str());
    }

  }




  if (debug) ws->Print();
  ws->writeToFile((workingDir+"/datacards/ws_"+histFolder+".root").c_str());
  if (debug) std::cout << "SAVING " << workingDir+"/datacards/ws_"+histFolder+".root" << '\n';

  DataCard.close();
  std::cout << "End of CreateWorkspace" << '\n';

  return;
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
  std::vector<std::string> histFolders = {"btag_DeepBoosted_H4qvsQCD", "btag_DeepBoosted_H4qvsQCDp2", "btag_DeepBoosted_H4qvsQCDp02", "btag_DeepBoosted_H4qvsQCDpt1000", "btag_DeepBoosted_H4qvsQCDpt1000p2", "btag_DeepBoosted_H4qvsQCDpt1000p02"};
  if (isHbb) histFolders = {"btag_DeepBoosted_HbbvsQCD", "btag_DeepBoosted_probHbb", "tau21" };
  std::vector<std::string> collections = {"Puppi"};
  std::vector<std::string> channels = {"muonchannel", "electronchannel"};
  std::vector<std::string> years = {"2016", "2017", "2018", "RunII"};
  // std::vector<std::string> channels = {"muonchannel"};
  // std::vector<std::string> years = {"RunII"};

  std::string studies = "nominal";
  if (argc>1) {
    std::string histFolder, channel, collection, year;
    for (int i = 1; i < argc; i++) {
      if (std::find(collections.begin(), collections.end(), argv[i]) != collections.end() ) collection = argv[i];
      if (std::find(histFolders.begin(), histFolders.end(), argv[i]) != histFolders.end() ) histFolder = argv[i];
      if (std::find(channels.begin(), channels.end(), argv[i]) != channels.end() ) channel = argv[i];
      if (std::find(years.begin(), years.end(), argv[i]) != years.end() ) year = argv[i];
    }
    std::cout << histFolder << " " << channel << " " << collection << " " << year << '\n';
    CreateWorkspace(studies,histFolder, channel, collection, year, isHbb);

  } else {
    std::cout << "more " << '\n';
    for (std::string year: years) {
      for (std::string collection: collections) {
        for (std::string channel: channels) {
          for (std::string histFolder: histFolders) {
            CreateWorkspace(studies,histFolder, channel, collection, year, isHbb);
          }
        }
      }
    }
  }

  return 0;

}
