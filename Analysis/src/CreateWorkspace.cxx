#include "CreateWorkspace.hpp"
#include "TSpline.h"

double GetRange(TH1F* h, double x) { return h->GetXaxis()->GetBinLowEdge(h->FindBin(x)); }

double CalculateIntegral(TH1F* h, double min, double max, bool doBinWidth) { return h->Integral(h->FindBin(min),h->FindBin(max),doBinWidth?"width":"");}

double CalculateFractionAreaPDF(RooAbsPdf* PDF, RooRealVar x_var, double fit_lo, double fit_hi) {
  x_var.setRange("fitting", fit_lo, fit_hi);
  return ((RooAbsReal*)(PDF)->createIntegral(x_var, RooFit::NormSet(x_var), RooFit::Range("fitting")))->getVal();
}

double DoFTest(double chi2_1, double chi2_2, double npar_1, double npar_2, double n) {return 1. - TMath::FDistI(((chi2_1-chi2_2)/(npar_2-npar_1))/(chi2_2/(n-npar_2-1)), npar_2-npar_1, n-npar_2); }

void CalculateChiSquare(double& chi2, int& nbins, RooHist* hpull, double xmin, double xmax) {
  // This method should take only the fitted points. plotter->chiSquare(npf) is on the full range.
  // diff in chi2 wrt to sum(hpull) are due to avarage(pull) vs interpolate(chiSquare) values in calculations
  chi2 = 0; nbins=0;
  for(int i=0 ; i<hpull->GetN() ; i++) {
    double x,pull;
    hpull->GetPoint(i,x,pull);
    if (x<xmin || x>xmax) continue;
    chi2+=pull*pull;
    nbins++;
  }
}


void CreateRooWorkspace::CalculateSignalFittingRange(double mass, double& rangeLo, double& rangeHi, double& plotLo, double& plotHi, double& ymax) {


  // rangeLo = mass*(1-5./28.);
  // rangeHi = mass*(1+3.5/28.);
  // rangeLo = 0.776*mass-56;
  // rangeHi = 1.085*mass+46;
  //
  // if (mass==1200) rangeLo = 790;
  // if (mass==1200) rangeHi = 1390;

  if (channel=="invisiblechannel"){
    rangeLo = 0.7742*mass-113;
    rangeHi = 1.0857*mass+46;

    if (mass==1000) { rangeLo = 1000; rangeHi = 1300;}
    if (mass==1200) { rangeLo = 1000; rangeHi = 1700;}
    if (mass==1400) { rangeLo = 1100; rangeHi = 1900;} // 1000- 1500
    if (mass==1600) { rangeLo = 1000; rangeHi = 2100;}
    if (mass==1800) { rangeLo = 1000; rangeHi = 2500;}
    // if (mass==1800) { rangeLo = myMin; rangeHi = myMax;}
    if (mass==2000) { rangeLo = 1200; rangeHi = 2600;} // 1300-2300 1500-2800
    // if (mass==2000) { rangeLo = myMin; rangeHi = myMax;}
    if (mass==2500) { rangeLo = 1000; rangeHi = 3400;}
    if (mass==3000) { rangeLo = 2100; rangeHi = 3300;} // 21-22-2300
    if (mass==3500) { rangeLo = 1100; rangeHi = 3900;} // 11-1200  38-3900
    if (mass==4000) { rangeLo = 1800; rangeHi = 4800;}
    if (mass==4500) { rangeLo = 1300; rangeHi = 4900;}
    // if (mass==4500) { rangeLo = myMin; rangeHi = myMax;}
    if (mass==5000) { rangeLo = 1400; rangeHi = 5900;} //6000
    // if (mass==5000) { rangeLo = myMin; rangeHi = myMax;}
    if (mass==5500) { rangeLo = 2500; rangeHi = 6000;}
    if (mass==6000) { rangeLo = 2800; rangeHi = 7000;}
    if (mass==7000) { rangeLo = 3000; rangeHi = 7600;}
    if (mass==8000) { rangeLo = 4000; rangeHi = 8800;}

    plotLo = mass*(1-10./28.);
    plotHi = mass*(1+10./28.);

    // Temporary changes
    plotLo = mass*(1-20./28.);
    plotHi = mass*(1+20./28.);

    if (FindInString("ZH",histFolder)) ymax = 5;
    else ymax = 8;

  } else {
    bool isEle = FindInString("ele",channel);
    rangeLo = 0.7742*mass-113;
    rangeHi = 1.0857*mass+46;

    if (mass==1000) { rangeLo = 700;  rangeHi = 1400;}
    if (mass==1200) { rangeLo = 700;  rangeHi = 1600;}

    if (mass==1400) { rangeLo = 1000; rangeHi = isEle? 1600: 1800;} //OK-ish
    if (mass==1600) { rangeLo = 1000; rangeHi = 1900;}//OK-ish TODO
    if (mass==1800) { rangeLo = 1100; rangeHi = 2000;}
    if (mass==2000) { rangeLo = 1200; rangeHi = 2200;}//OK //14-15-16 2200
    // if (mass==2000) { rangeLo = myMin; rangeHi = myMax;}
    if (mass==2500) { rangeLo = 1800; rangeHi = 2700;}//OK
    if (mass==3000) { rangeLo = 2100; rangeHi = 3300;}//OK-ish
    if (mass==3500) { rangeLo = 2400; rangeHi = 3900;}
    // if (mass==4000) { rangeLo = isEle? 2900: 2800; rangeHi = isEle? 4400: 4400;}//2900-4400
    if (mass==4000) { rangeLo = 3000; rangeHi = isEle?4400:4500;}//36-45 31-43 ele 30-42 38-43
    if (mass==4500) { rangeLo = 3400; rangeHi = 4900;}//OK
    if (mass==5000) { rangeLo = isEle? 4300: 3300; rangeHi = isEle? 5300: 5400;}//OK
    // if (mass==4000) { rangeLo = myMin; rangeHi = myMax;}
    if (mass==5500) { rangeLo = 4400; rangeHi = 6000;}
    if (mass==6000) { rangeLo = 4600; rangeHi = 6600;}
    if (mass==7000) { rangeLo = 5700; rangeHi = 7600;}
    if (mass==8000) { rangeLo = 6800; rangeHi = 8800;}

    plotLo = mass*(1-10./28.);
    plotHi = mass*(1+10./28.);

    ymax = 1.0;
    if (FindInString("100",HistName)) ymax = 2.5;

    if (FindInString("ZH",histFolder)) ymax = 2.5;
    else ymax = 5;

  }

  ymax *= lumi_map.at(year).at("lumi_fb")/lumi_map.at("RunII").at("lumi_fb");
  if (year=="2017") ymax *= 1.1;

  std::string hname = GetSgName(mass);
  rangeLo = GetRange(histo_map[hname].get(), rangeLo);
  rangeHi = GetRange(histo_map[hname].get(), rangeHi);
  plotLo  = GetRange(histo_map[hname].get(), plotLo);
  plotHi  = GetRange(histo_map[hname].get(), plotHi);

}

std::string GetSgName(int mass, std::string syst) {
  if (syst=="nominal") return "M"+std::to_string(mass);
  else return "M"+std::to_string(mass)+syst; // Be careful if you change it. It's needed as input to combine!!
}


// CreateRooWorkspace::CreateRooWorkspace(std::string year_, std::string collection_, std::string channel_, std::string histFolder_) : year(year_), collection(collection_), channel(channel_), histFolder(histFolder_) {
CreateRooWorkspace::CreateRooWorkspace(std::string year_, std::string collection_, std::string channel_, std::string histFolder_, std::string min_, std::string max_): year(year_), collection(collection_), channel(channel_), histFolder(histFolder_) {
  /*
  &  &&&&&&   &&&&&&&  &&    &&  &&&&&&  &&&&&&&& &&&&&&&&  &&     &&  &&&&&&  &&&&&&&&  &&&&&&&  &&&&&&&&
  & &&    && &&     && &&&   && &&    &&    &&    &&     && &&     && &&    &&    &&    &&     && &&     &&
  & &&       &&     && &&&&  && &&          &&    &&     && &&     && &&          &&    &&     && &&     &&
  & &&       &&     && && && &&  &&&&&&     &&    &&&&&&&&  &&     && &&          &&    &&     && &&&&&&&&
  & &&       &&     && &&  &&&&       &&    &&    &&   &&   &&     && &&          &&    &&     && &&   &&
  & &&    && &&     && &&   &&& &&    &&    &&    &&    &&  &&     && &&    &&    &&    &&     && &&    &&
  &  &&&&&&   &&&&&&&  &&    &&  &&&&&&     &&    &&     &&  &&&&&&&   &&&&&&     &&     &&&&&&&  &&     &&
  */


  std::cout << "****************************************" << std::endl;
  std::cout << "             CreateWorkspace            " << std::endl;
  std::cout << "----------------------------------------" << std::endl;
  std::cout << "year"       << std::string(15-4, ' ') << year       << '\n';
  std::cout << "channel"    << std::string(15-7, ' ') << channel    << '\n';
  std::cout << "collection" << std::string(15-10,' ') << collection << '\n';
  std::cout << "histFolder" << std::string(15-10,' ') << histFolder << '\n';
  std::cout << "myMin" << std::string(15-10,' ') << min_ << '\n';
  std::cout << "myMax" << std::string(15-10,' ') << max_ << '\n';
  std::cout << "****************************************\n" << std::endl;

  myMin = (min_=="")? 0 : stoi(min_);
  myMax = (max_=="")? 0 : stoi(max_);
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

  if (myMin>0) studies = studies+"_M_"+std::to_string(myMin);
  if (myMax>0) studies = studies+"_M_"+std::to_string(myMax);

  workingDir      = Path_ANALYSIS+"Analysis/"+(isHbb? "Limits/Hbb_": "Limits/")+studies+"/"+year+"/"+collection+"/"+channel+"/"+histFolder+"/";
  gSystem->Exec(("mkdir -p "+workingDir+"/datacards").c_str());

  //Set up list of syst, including variations, for the correct channel
  std::set_union(SystematicsScale.begin(), SystematicsScale.end(), SystematicsShape.begin(), SystematicsShape.end(), std::back_inserter(SystematicsAll));
  if (FindInVector(SystNames, "all")) {
    SystNames.erase(SystNames.begin()+FindInVector(SystNames,"all"));
    for (std::string syst: SystematicsAll) {
      if (!FindInString("muon",channel) && FindInString("tracking",syst)) continue;
      if (!FindInString("muon",channel) && FindInString("isolation",syst)) continue;
      if (!FindInString("muon",channel) && FindInString("MuonScale",syst)) continue;
      // Remove some of these when running with more systematics for the invisiblechannel
      if (FindInString("invisible",channel) && FindInString("pu",syst)) continue;
      if (FindInString("invisible",channel) && FindInString("btag",syst)) continue;
      if (FindInString("invisible",channel) && FindInString("prefiring",syst)) continue;
      if (FindInString("invisible",channel) && FindInString("id",syst)) continue;
      if (FindInString("invisible",channel) && FindInString("trigger",syst)) continue;
      if (FindInString("invisible",channel) && FindInString("reco",syst)) continue;
      // Be careful if you change it. It's needed as input to combine!!
      for (std::string var: {"Up","Down"}) {
        SystNames.push_back(syst+var);
      }
    }
  }
  if (debug) for (auto & syst : SystNames) std::cout << syst<< std::endl; //TODO print nicely


  //Select hist name
  Module    = "SignalRegion";
  Histtype  = "ZprimeCandidate";
  std::string massText = "mass";
  if (channel=="invisiblechannel") massText="mass_transversal";
  // HistName  = "Zprime_"+massText+"_rebin_full";
  // HistName  = "Zprime_"+massText+"";
  // HistName  = "Zprime_"+massText+"_rebin1";
  // HistName  = "Zprime_"+massText+"_rebin2";
  // HistName  = "Zprime_"+massText;
  // HistName  = "Zprime_"+massText+"_rebin30";
  HistName  = "Zprime_"+massText+"_rebin100";

  fit_lo_SR = ranges.at("SR").at(channel).at("fit_lo");
  fit_hi_SR = ranges.at("SR").at(channel).at("fit_hi");
  fit_lo_CR = ranges.at("CR").at(channel).at("fit_lo");
  fit_hi_CR = ranges.at("CR").at(channel).at("fit_hi");

  unique_name = "_"+channel+"_"+year;
  unique_name = TString(unique_name).ReplaceAll("channel","");
  unique_name_complete = unique_name.substr(1)+"_"+histFolder;
  filepath    = Path_STORAGE+year+"/"+Module+"/"+collection+"/"+channel+"/nominal/";

  dataName  = channel.substr(0, channel.find("channel")); dataName[0] = toupper(dataName[0]) ; dataName = "DATA_Single"+dataName;
  if (channel=="invisiblechannel") dataName = "DATA_MET";
  dataFileName = PrefixrootFile+"DATA."+dataName+"_"+year+"_noTree.root";

  x_var.reset(new RooRealVar("x_var", "m_{Zprime} (GeV)", x_lo, x_hi));
  ws.reset(new RooWorkspace((unique_name.substr(1)).c_str()));
  DataCard.open (workingDir+"datacards/OutputFit_"+histFolder+".txt");
  DataCard << "=== RooFit data fit result to be entered in datacard === " << std::endl;
  if (doFtest) output.open(workingDir+"datacards/output_"+year+"_"+histFolder+".txt");
  SignalProperties.open(workingDir+"datacards/SignalProperties_"+year+"_"+histFolder+".txt");

  nameXaxis = "m(Z')";
  if (channel=="invisiblechannel") nameXaxis = "m_{T}(Z')";
  nameYaxis = doBinWidth? "Events/bin" :"Events";
  nameRatioaxis = doPlotRatio?"Hist/Fit": "Pull";

  // rebin = 30;
  rebin = 0;
  bin = 50;
  bin2 = 100;
  for (double i = 0; i < 1500; i+=bin) bins_Zprime_rebin.push_back(i);
  for (double i = 1500+bin; i <= 4000; i+=bin2) bins_Zprime_rebin.push_back(i);
  bins_Zprime_rebin.push_back(4000);

  for (std::string name : NameSgPars) SgPars[name] = std::vector<double>(MyMassPoints.size(), 0);

  // xsec_ref_ = 0.1; // this is to mantain the signal strenght close to 1; (remember to multiply for this Normalization when plotting)
  xsec_ref_ = (xsec_ref.find(histFolder) != xsec_ref.end())? xsec_ref.at(histFolder): xsec_ref.at("default_value"); // default_value = 1
  DataCard << "xsec_ref " << " " << xsec_ref_ <<std::endl;

  if(channel!="invisiblechannel") {
    Modes.erase(Modes.begin()+distance(Modes.begin(), std::find(Modes.begin(), Modes.end(), "WJets_SR")));
    Modes.erase(Modes.begin()+distance(Modes.begin(), std::find(Modes.begin(), Modes.end(), "WJets_CR")));
  }

  for (auto x : Colors) { for (auto mode : Modes) doFits_map[mode][x.first]  = false; }

  // for (auto x : { "NO", "CB", "Exp_1", "Exp_2", "Exp_3", "Exp_4", "Exp_5", "Exp_6"}) { for (auto mode : Modes) doFits_map[mode][x]  = true; }
  // for (auto x : {"CB", "Exp_1", "Exp_2", "Exp_3" }) { for (auto mode : Modes) doFits_map[mode][x]  = true; }
  // for (auto x : {"Exp_2", "Exp_3" }) { for (auto mode : Modes) doFits_map[mode][x]  = true; }
  // for (auto x : { "NO", "CB", "Exp_2", "Exp_3"}) { for (auto mode : Modes) doFits_map[mode][x]  = true; }
  // for (auto x : { "Exp_2", "Exp_3", "Exp_4"}) { for (auto mode : Modes) doFits_map[mode][x]  = true; }
  for (auto x : { "Exp_2"}) { for (auto mode : Modes) doFits_map[mode][x]  = true; }
  // BkgPdf4p, BkgPdf3p TODO
  // for (auto x : { "NO", "CB", "Exp_4"}) { for (auto mode : Modes) doFits_map[mode][x]  = true; }
  // for (auto x : { "CB"}) { for (auto mode : Modes) doFits_map[mode][x]  = true; }

};




void CreateRooWorkspace::LoadFiles() {
  /*
  & &&        &&&&&&&     &&&    &&&&&&&&     &&     && &&&&  &&&&&&  &&&&&&&&  &&&&&&&   &&&&&&   &&&&&&&&     &&&    &&     &&  &&&&&&
  & &&       &&     &&   && &&   &&     &&    &&     &&  &&  &&    &&    &&    &&     && &&    &&  &&     &&   && &&   &&&   &&& &&    &&
  & &&       &&     &&  &&   &&  &&     &&    &&     &&  &&  &&          &&    &&     && &&        &&     &&  &&   &&  &&&& &&&& &&
  & &&       &&     && &&     && &&     &&    &&&&&&&&&  &&   &&&&&&     &&    &&     && &&   &&&& &&&&&&&&  &&     && && &&& &&  &&&&&&
  & &&       &&     && &&&&&&&&& &&     &&    &&     &&  &&        &&    &&    &&     && &&    &&  &&   &&   &&&&&&&&& &&     &&       &&
  & &&       &&     && &&     && &&     &&    &&     &&  &&  &&    &&    &&    &&     && &&    &&  &&    &&  &&     && &&     && &&    &&
  & &&&&&&&&  &&&&&&&  &&     && &&&&&&&&     &&     && &&&&  &&&&&&     &&     &&&&&&&   &&&&&&   &&     && &&     && &&     &&  &&&&&&
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
    if (channel!="invisiblechannel" && bkg=="WJets") continue;
    std::string fname = filepath+PrefixrootFile+"MC.MC_"+bkg+"_"+year+"_noTree.root";
    if (debug) std::cout << "Opening\t" << fname << '\n';
    f_map[bkg].reset(new TFile(fname.c_str()));

    for (std::string reg: {"CR","SR"}) {
      std::string hname = reg=="SR"? SRname: CRname;
      histo_map[bkg+"_"+reg].reset((TH1F*)((TH1F*)f_map[bkg]->Get((hname).c_str()))->Clone((bkg+" "+reg).c_str()));
      //Create sum of all MC samples
      if (histo_map.count("MC_"+reg)) histo_map["MC_"+reg]->Add(histo_map[bkg+"_"+reg].get());
      else histo_map["MC_"+reg].reset((TH1F*)((TH1F*)f_map[bkg]->Get((hname).c_str()))->Clone(("MC_"+reg).c_str()));
    }
  }

  for (auto syst: SystNames) {
    for (const int & mass : MyMassPoints) {

      // std::string SgName = "M"+std::to_string(mass);
      std::string fname  = filepath+PrefixrootFile+"MC.MC_ZprimeToZH_"+(channel=="invisiblechannel"?"inv_" : "")+"M"+std::to_string(mass)+"_"+year+"_noTree.root";
      std::string hmname = GetSgName(mass, syst);
      std::string hname  = SRname;
      if (!isNominalFolder(syst)) fname = TString(fname).ReplaceAll("nominal",syst).ReplaceAll("Up","_up").ReplaceAll("Down","_down"); // Be careful if you change it. It's needed as input to combine!!
      else if (!isNominalSyst(syst)) hname = TString(hname).ReplaceAll(Histtype,Histtype+"_"+syst).ReplaceAll("Up","_up").ReplaceAll("Down","_down"); // Be careful if you change it. It's needed as input to combine!!

      if (FindInString("murmuf",syst)) continue;
      if (FindInString("NNPDF",syst)) continue;
      if (debug) std::cout << "Opening\t" << fname << " " << hmname << " " << hname << '\n';

      f_map[hmname].reset(new TFile(fname.c_str()));
      histo_map[hmname].reset((TH1F*)((TH1F*)f_map[hmname]->Get((hname).c_str()))->Clone((hmname+" SR").c_str()));
    }
  }

  // Calculate QCD scale variations
  if (FindInVector(SystNames, "murmufUp") && FindInVector(SystNames, "murmufDown")) {
    for (const int & mass : MyMassPoints) {
      std::string syst = "murmuf";
      std::string fname  = filepath+PrefixrootFile+"MC.MC_ZprimeToZH_"+(channel=="invisiblechannel"?"inv_" : "")+"M"+std::to_string(mass)+"_"+year+"_noTree.root";
      std::string hname  = SRname;
      std::string hmname = GetSgName(mass, syst);
      std::string hmname_Up = GetSgName(mass, syst+"Up");
      std::string hmname_Down = GetSgName(mass, syst+"Down");
      if (!isNominalFolder(syst)) fname = TString(fname).ReplaceAll("nominal",syst).ReplaceAll("Up","_up").ReplaceAll("Down","_down"); // Be careful if you change it. It's needed as input to combine!!

      f_map[hmname].reset(new TFile(fname.c_str()));

      histo_map[hmname_Up].reset((TH1F*)((TH1F*)f_map[hmname]->Get((hname).c_str()))->Clone((hmname_Up+" SR").c_str()));
      histo_map[hmname_Down].reset((TH1F*)((TH1F*)f_map[hmname]->Get((hname).c_str()))->Clone((hmname_Down+" SR").c_str()));

      std::unordered_map<std::string, std::unique_ptr<TH1F> > histos_temp;
      for (auto var: Var_murmuf) {
        std::string hname_ = TString(hname).ReplaceAll(Histtype,Histtype+"_"+syst+"_"+var).Data();
        histos_temp[var].reset((TH1F*)((TH1F*)f_map[hmname]->Get((hname_).c_str()))->Clone((hname_+" SR").c_str()));
      }

      for (int bin = 0; bin < histo_map[hmname_Up]->GetNbinsX()+1; bin++) {
        std::vector<double> y_vals;
        for (const auto& x : histos_temp) y_vals.push_back(x.second->GetBinContent(bin));
        histo_map[hmname_Up]->SetBinContent(bin, *max_element(y_vals.begin(), y_vals.end()));
        histo_map[hmname_Down]->SetBinContent(bin, *min_element(y_vals.begin(), y_vals.end()));
      }
    }
  }

  // Calculate PDF variations
  if (FindInVector(SystNames, "NNPDFUp") && FindInVector(SystNames, "NNPDFDown")) {
    for (const int & mass : MyMassPoints) {
      std::string syst = "NNPDF";
      std::string fname  = filepath+PrefixrootFile+"MC.MC_ZprimeToZH_"+(channel=="invisiblechannel"?"inv_" : "")+"M"+std::to_string(mass)+"_"+year+"_noTree.root";
      std::string hname  = SRname;
      std::string hmname = GetSgName(mass, syst);
      std::string hmname_Up = GetSgName(mass, syst+"Up");
      std::string hmname_Down = GetSgName(mass, syst+"Down");
      if (!isNominalFolder(syst)) fname = TString(fname).ReplaceAll("nominal",syst).ReplaceAll("Up","_up").ReplaceAll("Down","_down"); // Be careful if you change it. It's needed as input to combine!!

      f_map[hmname].reset(new TFile(fname.c_str()));

      histo_map[hmname_Up].reset((TH1F*)((TH1F*)f_map[hmname]->Get((hname).c_str()))->Clone((hmname_Up+" SR").c_str()));
      histo_map[hmname_Down].reset((TH1F*)((TH1F*)f_map[hmname]->Get((hname).c_str()))->Clone((hmname_Down+" SR").c_str()));

      std::unordered_map<std::string, std::unique_ptr<TH1F> > histos_temp;
      for (int var = 0; var < PDF_variations ; var++) {
        std::string hname_ = TString(hname).ReplaceAll(Histtype,Histtype+"_"+syst+"_"+std::to_string(var)).Data();
        histos_temp[std::to_string(var)].reset((TH1F*)((TH1F*)f_map[hmname]->Get((hname_).c_str()))->Clone((hname_+" SR").c_str()));
        //TODO  normalize?
      }

      for (int bin = 0; bin < histo_map[hmname_Up]->GetNbinsX()+1; bin++) {
        std::vector<double> y_vals;
        for (const auto& x : histos_temp) y_vals.push_back(x.second->GetBinContent(bin));
        double mean = 0.; double stdev=0;
        for (const auto& x : y_vals) mean += x;
        mean /= PDF_variations;
        for (const auto& x : y_vals) stdev += (x-mean)*(x-mean);
        stdev = TMath::Sqrt(stdev/PDF_variations);
        histo_map[hmname_Up]->SetBinContent(bin,   histo_map[hmname_Up]->GetBinContent(bin) + stdev );
        histo_map[hmname_Down]->SetBinContent(bin, histo_map[hmname_Down]->GetBinContent(bin) - stdev );
      }
    }
  }



  histo_map["data"].reset((TH1F*)((TH1F*)f_map["data"]->Get((SRname).c_str()))->Clone("data SR"));
  histo_map["data_CR"].reset((TH1F*)((TH1F*)f_map["data"]->Get((CRname).c_str()))->Clone("data CR"));

  histo_map["MC_fake_CR"].reset((TH1F*)(histo_map["MC_CR"]->Clone("MC fake CR")));

  TRandom3 *r3=new TRandom3();

  for (int i = 0; i < histo_map["MC_fake_CR"]->GetNbinsX()+1; i++) {
    double start_val = histo_map["MC_fake_CR"]->GetBinContent(i);
    double new_val =0 , new_err=0;

    new_val = r3->PoissonD(start_val);
    new_err = TMath::Sqrt(new_val);

    if (new_val<2) {
      start_val *= 1000;
      new_val = r3->PoissonD(start_val);
      new_val /= 1000;
      new_err = TMath::Sqrt(new_val);
    }

    histo_map["MC_fake_CR"]->SetBinContent(i,new_val); histo_map["MC_fake_CR"]->SetBinError(i,new_err);


  }


  if (fitMC) {
    //TODO are the the cases we want?
    // bkg_pred = DY_SR for lepton channel
    // bkg_pred = DY_SR + WJets_SR for invisiblechannel
    histo_map["bkg_pred"].reset(((TH1F*)histo_map["DY_SR"]->Clone("bkg_pred: DY SR")));
    if (channel=="invisiblechannel"){
      histo_map["bkg_pred"]->Add(histo_map["WJets_SR"].get());
      histo_map["bkg_pred"]->SetName("bkg_pred: DY+WJets SR");
    }
  } else {
    histo_map["bkg_pred"].reset(((TH1F*) histo_map["data"]->Clone("bkg_pred: data SR") ));
  }

  // Need to make the histogram available after closing the files.
  for (const auto& x : histo_map) x.second->SetDirectory(0);

  if (debug){
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "CR hist_name:\t" << CRname << '\n';
    std::cout << "SR hist_name:\t" << SRname << '\n';
    std::cout << "fitMC:\t"        << BoolToString(fitMC) << '\n';
    std::cout << "doObs:\t"        << BoolToString(doObs) << '\n';
    std::cout << "bkg_pred on: "   << histo_map["bkg_pred"]->GetName() << '\n';
    std::cout << "data on: "       << histo_map["data"]->GetName() << '\n';
    std::cout << "----------------------------------------" << std::endl;
  }

  if (debug) for (const auto& x : histo_map) std::cout << "histo " << x.first << " " << std::string(20-x.first.size(),' ') << x.second->Integral() << '\n';

}

void CreateRooWorkspace::PrepocessHistos() {
  /*
  & &&&&&&&&  &&&&&&&&  &&&&&&&& &&&&&&&&   &&&&&&&   &&&&&&  &&&&&&&&  &&&&&&   &&&&&&     &&     && &&&&  &&&&&&  &&&&&&&&  &&&&&&&   &&&&&&
  & &&     && &&     && &&       &&     && &&     && &&    && &&       &&    && &&    &&    &&     &&  &&  &&    &&    &&    &&     && &&    &&
  & &&     && &&     && &&       &&     && &&     && &&       &&       &&       &&          &&     &&  &&  &&          &&    &&     && &&
  & &&&&&&&&  &&&&&&&&  &&&&&&   &&&&&&&&  &&     && &&       &&&&&&    &&&&&&   &&&&&&     &&&&&&&&&  &&   &&&&&&     &&    &&     &&  &&&&&&
  & &&        &&   &&   &&       &&        &&     && &&       &&             &&       &&    &&     &&  &&        &&    &&    &&     &&       &&
  & &&        &&    &&  &&       &&        &&     && &&    && &&       &&    && &&    &&    &&     &&  &&  &&    &&    &&    &&     && &&    &&
  & &&        &&     && &&&&&&&& &&         &&&&&&&   &&&&&&  &&&&&&&&  &&&&&&   &&&&&&     &&     && &&&&  &&&&&&     &&     &&&&&&&   &&&&&&
  */

  for (const auto& x : histo_map) {
    std::string mode = x.first;
    if (debug) std::cout << "pre " << mode << "\t" << CalculateIntegral(x.second.get(),fit_lo_SR,fit_hi_SR,doBinWidth) << '\n';
    if (dorebin) {
      if (rebin) histo_map[mode]->Rebin(rebin);
      else histo_map[mode].reset(dynamic_cast<TH1F*>(histo_map[mode]->Rebin(bins_Zprime_rebin.size()-1, x.second->GetName(), &bins_Zprime_rebin[0])));

      histo_map[mode]->Scale(1,doBinWidth?"width":"");
    }
    if (debug) std::cout << "rebin " << mode << "\t" << CalculateIntegral(x.second.get(),fit_lo_SR,fit_hi_SR,doBinWidth) << '\n';

    // Removing bins with low stats
    if (mode=="DY_SR" || mode=="MC_SR" || (mode=="bkg_pred" && fitMC)) {//TODO Check for the invisiblechannel
      for (int i = 0; i < x.second->GetNbinsX()+1; i++) {
        if (x.second->GetBinContent(i)<2*1e-03) { histo_map[mode]->SetBinContent(i,0); histo_map[mode]->SetBinError(i,0); }
      }
    }
  }

  //Create ratio of SR/CR for all samples. Do it only now otherwise you might have binning problems
  for (const auto& x : histo_map) {
    if (x.first.find("_SR") != std::string::npos) {
      std::string newName = TString(x.first).ReplaceAll("SR","ratio").Data();
      histo_map[newName].reset((TH1F*)(x.second->Clone(newName.c_str())));
      histo_map[newName]->Divide(histo_map[TString(x.first).ReplaceAll("SR","CR").Data()].get());
    }
  }


  // Normalize signal to arbitraty xsec.
  for (auto syst: SystNames) { for (const int & mass : MyMassPoints) { histo_map[GetSgName(mass,syst)]->Scale(xsec_ref_); } }

}

void CreateRooWorkspace::NormaliseData() {

  /*
  & &&&&&&&&     &&&    &&&&&&&&    &&&              &&&&&&&  &&&&&&&&   &&&&&&
  & &&     &&   && &&      &&      && &&            &&     && &&     && &&    &&
  & &&     &&  &&   &&     &&     &&   &&           &&     && &&     && &&
  & &&     && &&     &&    &&    &&     &&          &&     && &&&&&&&&   &&&&&&
  & &&     && &&&&&&&&&    &&    &&&&&&&&&          &&     && &&     &&       &&
  & &&     && &&     &&    &&    &&     &&          &&     && &&     && &&    &&
  & &&&&&&&&  &&     &&    &&    &&     &&           &&&&&&&  &&&&&&&&   &&&&&&
  */

  // Normalize h_MC_SR to h_Data_SR pretending it's data but has shape of bkg_pred in SR

  nEventsSR  = CalculateIntegral(histo_map["data"].get(),fit_lo_SR,fit_hi_SR,doBinWidth);

  if (fitMC) {
    // Normalization is taken such that it matches in the fitting range
    histo_map["bkg_pred"]->Scale(nEventsSR/CalculateIntegral(histo_map["bkg_pred"].get(),fit_lo_SR,fit_hi_SR,doBinWidth));
    histo_map["bkg_pred"]->SetName((std::string(histo_map["bkg_pred"]->GetName())+"*nEventsSR").c_str());
  }

  if (debug) {
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "data_obs\t"  << histo_map["data"]->GetName() << '\n';
    std::cout << "bkg_pred\t"  << histo_map["bkg_pred"]->GetName() << '\n';
    std::cout << "nEventsSR\t" << nEventsSR << '\n';
    std::cout << "norm_data\t" << CalculateIntegral(histo_map["bkg_pred"].get(),fit_lo_SR,fit_hi_SR,doBinWidth) << '\n';
    std::cout << "----------------------------------------" << std::endl;
  }


  histo_map["data_obs"].reset((TH1F*)(histo_map["data"]->Clone("data_obs")));
  for (int i = 0; i < histo_map["data_obs"]->GetNbinsX()+1; i++) {
    if (histo_map["data_obs"]->GetXaxis()->GetBinUpEdge(i)<fit_lo_SR || histo_map["data_obs"]->GetXaxis()->GetBinLowEdge(i)>fit_hi_SR){
      histo_map["data_obs"]->SetBinContent(i,0); histo_map["data_obs"]->SetBinError(i,0);

    }
  }

  data_obs.reset(new RooDataHist("data_obs", histo_map["data_obs"]->GetName(), RooArgList(*x_var), histo_map["data_obs"].get()));
  DataCard << " Background number of events = " << nEventsSR << " " << CalculateIntegral(histo_map["data"].get(),0,10000,doBinWidth) << std::endl;
}


void CreateRooWorkspace::DoRebin() {
  /*
  & &&&&&&&&   &&&&&&&     &&&&&&&&  &&&&&&&& &&&&&&&&  &&&& &&    &&
  & &&     && &&     &&    &&     && &&       &&     &&  &&  &&&   &&
  & &&     && &&     &&    &&     && &&       &&     &&  &&  &&&&  &&
  & &&     && &&     &&    &&&&&&&&  &&&&&&   &&&&&&&&   &&  && && &&
  & &&     && &&     &&    &&   &&   &&       &&     &&  &&  &&  &&&&
  & &&     && &&     &&    &&    &&  &&       &&     &&  &&  &&   &&&
  & &&&&&&&&   &&&&&&&     &&     && &&&&&&&& &&&&&&&&  &&&& &&    &&
  */


  if (! dorebin || rebin==0) return;

  // RooBinning binning(bins_Zprime_rebin.size()-1, &bins_Zprime_rebin[0], "rebin");
  int Nbins = 0;
  double xmin = 1e6;
  double xmax = 0;
  for (int i=0; i<histo_map["bkg_pred"]->GetNbinsX()+1; ++i){
    if (histo_map["bkg_pred"]->GetXaxis()->GetBinLowEdge(i)>=fit_lo_SR){
      xmin = std::min(xmin, histo_map["bkg_pred"]->GetXaxis()->GetBinLowEdge(i));
      if (histo_map["bkg_pred"]->GetXaxis()->GetBinUpEdge(i)<=fit_hi_SR){
        xmax = std::max(xmax, histo_map["bkg_pred"]->GetXaxis()->GetBinUpEdge(i));
        ++Nbins;
      }
    }
  }
  RooBinning binningFitting(Nbins, xmin, xmax, "rebin");
  if (debug) std::cout << "binning fitting " << Nbins << " " << xmin << " " << xmax << binningFitting << '\n';
  x_var->setBinning(binningFitting);//TODO also when no rebin?

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
  & &&&& &&    && &&&& &&&&&&&& &&&&    &&&    &&       &&&& &&&&&&&& &&&&&&&&    &&&&&&&&  &&&&&&&&  &&&&&&&&  &&&&&&
  &  &&  &&&   &&  &&     &&     &&    && &&   &&        &&       &&  &&          &&     && &&     && &&       &&    &&
  &  &&  &&&&  &&  &&     &&     &&   &&   &&  &&        &&      &&   &&          &&     && &&     && &&       &&
  &  &&  && && &&  &&     &&     &&  &&     && &&        &&     &&    &&&&&&      &&&&&&&&  &&     && &&&&&&    &&&&&&
  &  &&  &&  &&&&  &&     &&     &&  &&&&&&&&& &&        &&    &&     &&          &&        &&     && &&             &&
  &  &&  &&   &&&  &&     &&     &&  &&     && &&        &&   &&      &&          &&        &&     && &&       &&    &&
  & &&&& &&    && &&&&    &&    &&&& &&     && &&&&&&&& &&&& &&&&&&&& &&&&&&&&    &&        &&&&&&&&  &&        &&&&&&
  */

  std::string FitName;

  for (auto syst: SystNames) {
    for (const int & mass : MyMassPoints) {
      FitName = GetSgName(mass,syst);
      std::string ParName = "_"+FitName+unique_name;

      if (channel=="invisiblechannel"){
        fitPars[FitName].emplace_back(new RooRealVar(("sg_p0"+ParName).c_str(), ("sg_p0_"+ParName).c_str(), mass, mass*0.7, mass*1.05)); // mu
        fitPars[FitName].emplace_back(new RooRealVar(("sg_p1"+ParName).c_str(), ("sg_p1_"+ParName).c_str(), 80,10,600));                 // sigma - width
        fitPars[FitName].emplace_back(new RooRealVar(("sg_p2"+ParName).c_str(), ("sg_p2_"+ParName).c_str(), 0.3, -2, 10));
        fitPars[FitName].emplace_back(new RooRealVar(("sg_p3"+ParName).c_str(), ("sg_p3_"+ParName).c_str(), 40, -50, 150));
      } else {
        fitPars[FitName].emplace_back(new RooRealVar(("sg_p0"+ParName).c_str(), ("sg_p0_"+ParName).c_str(), mass, mass*0.7, mass*1.3));
        fitPars[FitName].emplace_back(new RooRealVar(("sg_p1"+ParName).c_str(), ("sg_p1_"+ParName).c_str(), mass*0.02+20., 10., 500));
        fitPars[FitName].emplace_back(new RooRealVar(("sg_p2"+ParName).c_str(), ("sg_p2_"+ParName).c_str(), 0.7, 0, 10));
        fitPars[FitName].emplace_back(new RooRealVar(("sg_p3"+ParName).c_str(), ("sg_p3_"+ParName).c_str(), 5.0, 0, 20));
      }

      Fits_map[FitName]["Gauss"].reset(new RooGaussian(FitName.c_str(), "Signal Prediction", *x_var, *fitPars[FitName][0], *fitPars[FitName][1]));
      Fits_map[FitName]["CB"].reset(new RooCBShape(FitName.c_str(), "Signal Prediction", *x_var, *fitPars[FitName][0], *fitPars[FitName][1], *fitPars[FitName][2], *fitPars[FitName][3]));
      // Fits_map[FitName]["CB_fit"].reset(new RooCBShape(FitName.c_str(), "Signal Prediction", *x_var_sig, *fitPars[(FitName"_fit").c_str()][0], *fitPars[(FitName"_fit").c_str()][1], *fitPars[(FitName"_fit").c_str()][2], *fitPars[(FitName"_fit").c_str()][3]));
    }
  }


  for (auto mode: Modes) {

    // Old functions
    FitName = mode+"_NO";
    fitPars[FitName].emplace_back(new RooRealVar((mode+"_Novo_p1"+unique_name).c_str(), (mode+"_Novo_p1"+unique_name).c_str(),  500, 200,  1000.));
    fitPars[FitName].emplace_back(new RooRealVar((mode+"_Novo_p2"+unique_name).c_str(), (mode+"_Novo_p2"+unique_name).c_str(),  100,  10,  200.));
    fitPars[FitName].emplace_back(new RooRealVar((mode+"_Novo_p3"+unique_name).c_str(), (mode+"_Novo_p3"+unique_name).c_str(), -0.6,  -2,  2. ));
    Fits_map[mode]["NO"].reset(new RooNovosibirsk((FitName+unique_name).c_str(), (FitName+unique_name).c_str(),*x_var, *fitPars[FitName][0], *fitPars[FitName][1], *fitPars[FitName][2]));

    //TODO To be put in the old functions
    FitName = mode+"_CB";
    fitPars[FitName].emplace_back(new RooRealVar((mode+"_CB_p1"+unique_name).c_str(), (mode+"_CB_p1"+unique_name).c_str(), +1,   0,   10));
    fitPars[FitName].emplace_back(new RooRealVar((mode+"_CB_p2"+unique_name).c_str(), (mode+"_CB_p2"+unique_name).c_str(), +10,  0,   100));
    fitPars[FitName].emplace_back(new RooRealVar((mode+"_CB_p3"+unique_name).c_str(), (mode+"_CB_p3"+unique_name).c_str(), +500, 100, 1000));
    fitPars[FitName].emplace_back(new RooRealVar((mode+"_CB_p4"+unique_name).c_str(), (mode+"_CB_p4"+unique_name).c_str(), +100, 10,  500));
    Fits_map[mode]["CB"].reset(new RevCrystalBall((FitName+unique_name).c_str(), (FitName+unique_name).c_str(), *x_var, *fitPars[FitName][0], *fitPars[FitName][1], *fitPars[FitName][2], *fitPars[FitName][3]));

    FitName = mode+"_Exp_1";
    fitPars[FitName].emplace_back(new RooRealVar((mode+"_Exp_1_p1"+unique_name).c_str(), (mode+"_Exp_1_p1"+unique_name).c_str(), 0., -100, 100));
    Fits_map[mode]["Exp_1"].reset(new PolinomialExponent_1p((FitName+unique_name).c_str(), (FitName+unique_name).c_str(),*x_var, *fitPars[FitName][0]));

    if (channel=="invisiblechannel"){
      FitName = mode+"_Exp_2";
      fitPars[FitName].emplace_back(new RooRealVar((mode+"_Exp_2_p1"+unique_name).c_str(), (mode+"_Exp_2_p1"+unique_name).c_str(), -40., -100, 100));
      fitPars[FitName].emplace_back(new RooRealVar((mode+"_Exp_2_p2"+unique_name).c_str(), (mode+"_Exp_2_p2"+unique_name).c_str(), -4, -100, 100));
      Fits_map[mode]["Exp_2"].reset(new PolinomialExponent_2p((FitName+unique_name).c_str(), (FitName+unique_name).c_str(),*x_var, *fitPars[FitName][0], *fitPars[FitName][1]));

      FitName = mode+"_Exp_3";
      fitPars[FitName].emplace_back(new RooRealVar((mode+"_Exp_3_p1"+unique_name).c_str(), (mode+"_Exp_3_p1"+unique_name).c_str(), -40., -100, 100));
      fitPars[FitName].emplace_back(new RooRealVar((mode+"_Exp_3_p2"+unique_name).c_str(), (mode+"_Exp_3_p2"+unique_name).c_str(), -10, -100, 100));
      fitPars[FitName].emplace_back(new RooRealVar((mode+"_Exp_3_p3"+unique_name).c_str(), (mode+"_Exp_3_p3"+unique_name).c_str(), 10, -100, 100));
      Fits_map[mode]["Exp_3"].reset(new PolinomialExponent_3p((FitName+unique_name).c_str(), (FitName+unique_name).c_str(),*x_var, *fitPars[FitName][0], *fitPars[FitName][1], *fitPars[FitName][2]));
    } else {
      FitName = mode+"_Exp_2";
      fitPars[FitName].emplace_back(new RooRealVar((mode+"_Exp_2_p1"+unique_name).c_str(), (mode+"_Exp_2_p1"+unique_name).c_str(), 0., -4, 4));
      fitPars[FitName].emplace_back(new RooRealVar((mode+"_Exp_2_p2"+unique_name).c_str(), (mode+"_Exp_2_p2"+unique_name).c_str(), 0., -10, 10));
      Fits_map[mode]["Exp_2"].reset(new PolinomialExponent_2p((FitName+unique_name).c_str(), (FitName+unique_name).c_str(),*x_var, *fitPars[FitName][0], *fitPars[FitName][1]));

      FitName = mode+"_Exp_3";
      fitPars[FitName].emplace_back(new RooRealVar((mode+"_Exp_3_p1"+unique_name).c_str(), (mode+"_Exp_3_p1"+unique_name).c_str(), 0., -4, 4));
      fitPars[FitName].emplace_back(new RooRealVar((mode+"_Exp_3_p2"+unique_name).c_str(), (mode+"_Exp_3_p2"+unique_name).c_str(), 0.4, -100, 100));
      fitPars[FitName].emplace_back(new RooRealVar((mode+"_Exp_3_p3"+unique_name).c_str(), (mode+"_Exp_3_p3"+unique_name).c_str(), -0.1, -100, 100));
      Fits_map[mode]["Exp_3"].reset(new PolinomialExponent_3p((FitName+unique_name).c_str(), (FitName+unique_name).c_str(),*x_var, *fitPars[FitName][0], *fitPars[FitName][1], *fitPars[FitName][2]));
    }
    FitName = mode+"_Exp_4";
    fitPars[FitName].emplace_back(new RooRealVar((mode+"_Exp_4_p1"+unique_name).c_str(), (mode+"_Exp_4_p1"+unique_name).c_str(), 4., -100, 100));
    fitPars[FitName].emplace_back(new RooRealVar((mode+"_Exp_4_p2"+unique_name).c_str(), (mode+"_Exp_4_p2"+unique_name).c_str(), -10, -100,  10));
    fitPars[FitName].emplace_back(new RooRealVar((mode+"_Exp_4_p3"+unique_name).c_str(), (mode+"_Exp_4_p3"+unique_name).c_str(), 4, -100, 100));
    fitPars[FitName].emplace_back(new RooRealVar((mode+"_Exp_4_p4"+unique_name).c_str(), (mode+"_Exp_4_p4"+unique_name).c_str(), -0.4, -100, 100));
    Fits_map[mode]["Exp_4"].reset(new PolinomialExponent_4p((FitName+unique_name).c_str(), (FitName+unique_name).c_str(),*x_var, *fitPars[FitName][0], *fitPars[FitName][1], *fitPars[FitName][2], *fitPars[FitName][3]));

    FitName = mode+"_Exp_5";
    fitPars[FitName].emplace_back(new RooRealVar((mode+"_Exp_5_p1"+unique_name).c_str(), (mode+"_Exp_5_p1"+unique_name).c_str(), -8.46, -9.00, -8.00));
    fitPars[FitName].emplace_back(new RooRealVar((mode+"_Exp_5_p2"+unique_name).c_str(), (mode+"_Exp_5_p2"+unique_name).c_str(),  7.75,  7.00,  9.00));
    fitPars[FitName].emplace_back(new RooRealVar((mode+"_Exp_5_p3"+unique_name).c_str(), (mode+"_Exp_5_p3"+unique_name).c_str(), -7.45, -8.00, -6.00));
    fitPars[FitName].emplace_back(new RooRealVar((mode+"_Exp_5_p4"+unique_name).c_str(), (mode+"_Exp_5_p4"+unique_name).c_str(),  3.35,  3.00,  4.00));
    fitPars[FitName].emplace_back(new RooRealVar((mode+"_Exp_5_p5"+unique_name).c_str(), (mode+"_Exp_5_p5"+unique_name).c_str(), -5.37, -6.00, -4.00));
    Fits_map[mode]["Exp_5"].reset(new PolinomialExponent_5p((FitName+unique_name).c_str(), (FitName+unique_name).c_str(),*x_var, *fitPars[FitName][0], *fitPars[FitName][1], *fitPars[FitName][2], *fitPars[FitName][3], *fitPars[FitName][4]));

    FitName = mode+"_Exp_6";
    fitPars[FitName].emplace_back(new RooRealVar((mode+"_Exp_6_p1"+unique_name).c_str(), (mode+"_Exp_6_p1"+unique_name).c_str(), -3.68, -4.0, -3.0 ));
    fitPars[FitName].emplace_back(new RooRealVar((mode+"_Exp_6_p2"+unique_name).c_str(), (mode+"_Exp_6_p2"+unique_name).c_str(), -1.57, -2.0, -1.0 ));
    fitPars[FitName].emplace_back(new RooRealVar((mode+"_Exp_6_p3"+unique_name).c_str(), (mode+"_Exp_6_p3"+unique_name).c_str(), +1.56, +1.0, +2.0 ));
    fitPars[FitName].emplace_back(new RooRealVar((mode+"_Exp_6_p4"+unique_name).c_str(), (mode+"_Exp_6_p4"+unique_name).c_str(), -1.14, -1.5, -1.0 ));
    fitPars[FitName].emplace_back(new RooRealVar((mode+"_Exp_6_p5"+unique_name).c_str(), (mode+"_Exp_6_p5"+unique_name).c_str(), +5.29, +5.0, +5.5 ));
    fitPars[FitName].emplace_back(new RooRealVar((mode+"_Exp_6_p6"+unique_name).c_str(), (mode+"_Exp_6_p6"+unique_name).c_str(), -7.05, -9.5, -9.0 ));
    Fits_map[mode]["Exp_6"].reset(new PolinomialExponent_6p((FitName+unique_name).c_str(), (FitName+unique_name).c_str(),*x_var, *fitPars[FitName][0], *fitPars[FitName][1], *fitPars[FitName][2], *fitPars[FitName][3], *fitPars[FitName][4], *fitPars[FitName][5]));

    // Not fitting functions
    // FitName = mode+"_BkgPdf4p";
    // fitPars[FitName].emplace_back(new RooRealVar((mode+"_BkgPdf4p_p1"+unique_name).c_str(), (mode+"_BkgPdf4p_p1"+unique_name).c_str(), +2.0, +2.0, +4.0));
    // fitPars[FitName].emplace_back(new RooRealVar((mode+"_BkgPdf4p_p2"+unique_name).c_str(), (mode+"_BkgPdf4p_p2"+unique_name).c_str(), +4.2, +3.0, +5.0));
    // fitPars[FitName].emplace_back(new RooRealVar((mode+"_BkgPdf4p_p3"+unique_name).c_str(), (mode+"_BkgPdf4p_p3"+unique_name).c_str(), +1.4, +0.5, +2.0));
    // fitPars[FitName].emplace_back(new RooRealVar((mode+"_BkgPdf4p_p4"+unique_name).c_str(), (mode+"_BkgPdf4p_p4"+unique_name).c_str(), -1.0, -2.0, -0.5));
    // Fits_map[mode]["BkgPdf4p"].reset(new BkgPdf4p((FitName+unique_name).c_str(), (FitName+unique_name).c_str(),*x_var, *fitPars[FitName][0], *fitPars[FitName][1], *fitPars[FitName][2], *fitPars[FitName][3]));
    //
    // FitName = mode+"_BkgPdf3p";
    // fitPars[FitName].emplace_back(new RooRealVar((mode+"_BkgPdf3p_p1"+unique_name).c_str(), (mode+"_BkgPdf3p_p1"+unique_name).c_str(), -10., -11., -9.0));
    // fitPars[FitName].emplace_back(new RooRealVar((mode+"_BkgPdf3p_p2"+unique_name).c_str(), (mode+"_BkgPdf3p_p2"+unique_name).c_str(), +7.0, +6.0, +8.0));
    // fitPars[FitName].emplace_back(new RooRealVar((mode+"_BkgPdf3p_p3"+unique_name).c_str(), (mode+"_BkgPdf3p_p3"+unique_name).c_str(), +3.7, +3.0, +4.0));

    // fitPars[FitName].emplace_back(new RooRealVar((mode+"_BkgPdf3p_p1"+unique_name).c_str(), (mode+"_BkgPdf3p_p1"+unique_name).c_str(), -13.2, -1000, 1000));
    // fitPars[FitName].emplace_back(new RooRealVar((mode+"_BkgPdf3p_p2"+unique_name).c_str(), (mode+"_BkgPdf3p_p2"+unique_name).c_str(), +9.1, -1000, 1000));
    // fitPars[FitName].emplace_back(new RooRealVar((mode+"_BkgPdf3p_p3"+unique_name).c_str(), (mode+"_BkgPdf3p_p3"+unique_name).c_str(), +2.5, -100, 100));
    //
    // fitPars[FitName].emplace_back(new RooRealVar((mode+"_BkgPdf3p_p1"+unique_name).c_str(), (mode+"_BkgPdf3p_p1"+unique_name).c_str(), -4.3, -5.3, -3.3));
    // fitPars[FitName].emplace_back(new RooRealVar((mode+"_BkgPdf3p_p2"+unique_name).c_str(), (mode+"_BkgPdf3p_p2"+unique_name).c_str(), +5.7, +4.7, +7.7));
    // fitPars[FitName].emplace_back(new RooRealVar((mode+"_BkgPdf3p_p3"+unique_name).c_str(), (mode+"_BkgPdf3p_p3"+unique_name).c_str(), +2.5, +1.5, +3.5));
    // Fits_map[mode]["BkgPdf3p"].reset(new BkgPdf3p((FitName+unique_name).c_str(), (FitName+unique_name).c_str(),*x_var, *fitPars[FitName][0], *fitPars[FitName][1], *fitPars[FitName][2]));


  }
}


void CreateRooWorkspace::CreateRooDataHist() {

  /*
  &  &&&&&&  &&&&&&&&  &&&&&&&&    &&&    &&&&&&&& &&&&&&&&             &&&&&&&&   &&&&&&&   &&&&&&&  &&&&&&&&     &&&    &&&&&&&&    &&&    &&     && &&&&  &&&&&&  &&&&&&&&
  & &&    && &&     && &&         && &&      &&    &&                   &&     && &&     && &&     && &&     &&   && &&      &&      && &&   &&     &&  &&  &&    &&    &&
  & &&       &&     && &&        &&   &&     &&    &&                   &&     && &&     && &&     && &&     &&  &&   &&     &&     &&   &&  &&     &&  &&  &&          &&
  & &&       &&&&&&&&  &&&&&&   &&     &&    &&    &&&&&&               &&&&&&&&  &&     && &&     && &&     && &&     &&    &&    &&     && &&&&&&&&&  &&   &&&&&&     &&
  & &&       &&   &&   &&       &&&&&&&&&    &&    &&                   &&   &&   &&     && &&     && &&     && &&&&&&&&&    &&    &&&&&&&&& &&     &&  &&        &&    &&
  & &&    && &&    &&  &&       &&     &&    &&    &&                   &&    &&  &&     && &&     && &&     && &&     &&    &&    &&     && &&     &&  &&  &&    &&    &&
  &  &&&&&&  &&     && &&&&&&&& &&     &&    &&    &&&&&&&&             &&     &&  &&&&&&&   &&&&&&&  &&&&&&&&  &&     &&    &&    &&     && &&     && &&&&  &&&&&&     &&
  */


  if (dorebin && rebin!=0) {
    for (auto mode: Modes) rooHist_map[mode].reset(new RooDataHist(mode.c_str(), mode.c_str(), RooArgList(*x_var), RooFit::Import(*histo_map[mode].get(),doBinWidth)));
    for (auto syst: SystNames) {
      for (const int & mass : MyMassPoints) {
        std::string hname = GetSgName(mass,syst);
        rooHist_map[hname].reset(new RooDataHist(hname.c_str(), hname.c_str(), RooArgList(*x_var), RooFit::Import(*histo_map[hname].get(),doBinWidth)));
      }
    }
  } else {
    for (auto mode: Modes) {
      rooHist_map[mode].reset(new RooDataHist(mode.c_str(), mode.c_str(), RooArgList(*x_var), histo_map[mode].get()));
    }
    for (auto syst: SystNames) {
      for (const int & mass : MyMassPoints) {
        std::string hname = GetSgName(mass,syst);
        rooHist_map[hname].reset(new RooDataHist(hname.c_str(), hname.c_str(), RooArgList(*x_var), histo_map[hname].get()));
      }
    }
  }

}

void CreateRooWorkspace::DoFits() {
  /*
  & &&&&&&&&   &&&&&&&              &&&&&&&& &&&& &&&&&&&&  &&&&&&
  & &&     && &&     &&             &&        &&     &&    &&    &&
  & &&     && &&     &&             &&        &&     &&    &&
  & &&     && &&     &&             &&&&&&    &&     &&     &&&&&&
  & &&     && &&     &&             &&        &&     &&          &&
  & &&     && &&     &&             &&        &&     &&    &&    &&
  & &&&&&&&&   &&&&&&&              &&       &&&&    &&     &&&&&&
  */

  for (const auto& x : histo_map) {
    fit_min[x.first] = GetRange(x.second.get(),fit_lo_SR);
    fit_max[x.first] = GetRange(x.second.get(),fit_hi_SR);
    // if ("DY_SR"==x.first) fit_min[x.first] = GetRange(x.second.get(),fit_SR);
    // if ("DY_CR"==x.first) fit_min[x.first] = GetRange(x.second.get(),fit_lo_SR);
    // if ("DY_CR"==x.first) fit_max[x.first] = GetRange(x.second.get(),fit_max_CR);
    if (FindInString("CR", x.first)) fit_min[x.first] = GetRange(x.second.get(),1300);//TODO
    if (FindInString("CR", x.first)) fit_max[x.first] = GetRange(x.second.get(),6000);
    if (FindInString("data_CR", x.first)) fit_min[x.first] = GetRange(x.second.get(),fit_lo_CR);
    if (FindInString("data_CR", x.first)) fit_max[x.first] = GetRange(x.second.get(),fit_hi_CR);

    std::string mass = x.first;
    if (mass.compare(0,1,"M")==0 && mass.compare(0,2,"MC")!=0) {
      CalculateSignalFittingRange(std::stoi(mass.substr(1, mass.find("_")-1)), fit_min[x.first], fit_max[x.first], plot_min[x.first], plot_max[x.first], y_max[x.first]); // TODO
      // Normalization into the fitting range. Need to be corrected for the effected area of the PDF later.
      nEventsSignal[x.first] = CalculateIntegral(histo_map[x.first].get(),fit_min[x.first],fit_max[x.first],doBinWidth);//TODO
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
    std::cout << "Fitting range:\t fit_min\t fit_max" << '\n';
    // std::cout << "sub_fit" << std::string(10,' ') << fit_min_turnOn << "\t" << fit_max_tails << std::endl; // not used
    for (auto mode : Modes) std::cout << mode << std::string(17-mode.size(),' ') << fit_min[mode] << "\t" << fit_max[mode] << std::endl;
    std::cout << "----------------------------------------" << std::endl;
  }

  for (auto mode: Modes) {
    for (auto const& [model,dofit] : doFits_map[mode] ) {
      if (!dofit) continue;
      if (mode=="data") {
        FitRes_map[mode][model].reset(Fits_map[mode][model]->fitTo(*rooHist_map[mode], RooFit::Range(fit_min[mode], fit_max[mode]), RooFit::SumW2Error(kFALSE), RooFit::Minimizer("Minuit2"), RooFit::Save(), RooFit::Verbose(kFALSE), RooFit::PrintEvalErrors(-1)));
      } else {
        FitRes_map[mode][model].reset(Fits_map[mode][model]->fitTo(*rooHist_map[mode], RooFit::Range(fit_min[mode], fit_max[mode]), RooFit::SumW2Error(kTRUE), RooFit::Save(), RooFit::Verbose(kFALSE), RooFit::PrintEvalErrors(-1)));
      }
    }
  }

  std::cout << "********************" << '\n';
  std::cout << "*    Signal Fits   *" << '\n';
  std::cout << "********************" << '\n';

  for (auto syst: SystNames) {
    for (const int & mass : MyMassPoints) {
      std::string hname = GetSgName(mass,syst);
      if (debug) std::cout << "fit mass:" << syst << "\t" << mass << "\t" << fit_min[hname] << " -- " << fit_max[hname] << '\n';
      FitRes_map[hname][FitSignal].reset(Fits_map[hname][FitSignal]->fitTo(*rooHist_map[hname], RooFit::Range(fit_min[hname], fit_max[hname]), RooFit::SumW2Error(kTRUE), RooFit::Save(), RooFit::Verbose(kFALSE), RooFit::PrintEvalErrors(-1)));
    }
  }

}


void CreateRooWorkspace::ImportToWorkspace() {
  /*
  & &&&& &&     && &&&&&&&&   &&&&&&&  &&&&&&&&  &&&&&&&& &&&&&&&&  &&&&&&&  &&      &&  &&&&&&&  &&&&&&&&  &&    &&  &&&&&&  &&&&&&&&     &&&     &&&&&&  &&&&&&&&
  &  &&  &&&   &&& &&     && &&     && &&     &&    &&       &&    &&     && &&  &&  && &&     && &&     && &&   &&  &&    && &&     &&   && &&   &&    && &&
  &  &&  &&&& &&&& &&     && &&     && &&     &&    &&       &&    &&     && &&  &&  && &&     && &&     && &&  &&   &&       &&     &&  &&   &&  &&       &&
  &  &&  && &&& && &&&&&&&&  &&     && &&&&&&&&     &&       &&    &&     && &&  &&  && &&     && &&&&&&&&  &&&&&     &&&&&&  &&&&&&&&  &&     && &&       &&&&&&
  &  &&  &&     && &&        &&     && &&   &&      &&       &&    &&     && &&  &&  && &&     && &&   &&   &&  &&         && &&        &&&&&&&&& &&       &&
  &  &&  &&     && &&        &&     && &&    &&     &&       &&    &&     && &&  &&  && &&     && &&    &&  &&   &&  &&    && &&        &&     && &&    && &&
  & &&&& &&     && &&         &&&&&&&  &&     &&    &&       &&     &&&&&&&   &&&  &&&   &&&&&&&  &&     && &&    &&  &&&&&&  &&        &&     &&  &&&&&&  &&&&&&&&
  */

  // Reduce the x range to avoid strange fitting behaviour: very delicate!
  x_var.reset(new RooRealVar("x_var", "m_{Zprime} (GeV)", x_lo_short, x_hi_short));
  std::string FitName;
  for (auto syst: SystNames) { for (const int & mass : MyMassPoints) {
    FitName = GetSgName(mass,syst);
    Fits_map[FitName]["CB"].reset(new RooCBShape(FitName.c_str(), "Signal Prediction", *x_var, *fitPars[FitName][0], *fitPars[FitName][1], *fitPars[FitName][2], *fitPars[FitName][3]));
    ws->import(*Fits_map[FitName][FitSignal], RooFit::Silence());
  }}
  for (auto mode: Modes) { for (auto const& [model,dofit] : doFits_map[mode] ) {
    if (dofit) {
      FitName = mode+"_"+model;
      if (model=="NO") Fits_map[mode]["NO"].reset(new RooNovosibirsk((FitName+unique_name).c_str(), (FitName+unique_name).c_str(),*x_var, *fitPars[FitName][0], *fitPars[FitName][1], *fitPars[FitName][2]));
      if (model=="CB") Fits_map[mode]["CB"].reset(new RevCrystalBall((FitName+unique_name).c_str(), (FitName+unique_name).c_str(), *x_var, *fitPars[FitName][0], *fitPars[FitName][1], *fitPars[FitName][2], *fitPars[FitName][3]));
      if (model=="Exp_2") Fits_map[mode]["Exp_2"].reset(new PolinomialExponent_2p((FitName+unique_name).c_str(), (FitName+unique_name).c_str(),*x_var, *fitPars[FitName][0], *fitPars[FitName][1]));
      if (model=="Exp_3") Fits_map[mode]["Exp_3"].reset(new PolinomialExponent_3p((FitName+unique_name).c_str(), (FitName+unique_name).c_str(),*x_var, *fitPars[FitName][0], *fitPars[FitName][1], *fitPars[FitName][2]));
      if (model=="Exp_4") Fits_map[mode]["Exp_4"].reset(new PolinomialExponent_4p((FitName+unique_name).c_str(), (FitName+unique_name).c_str(),*x_var, *fitPars[FitName][0], *fitPars[FitName][1], *fitPars[FitName][2], *fitPars[FitName][3]));

      ws->import(*Fits_map[mode][model], RooFit::Silence());
    }
  }}
  data_obs.reset(new RooDataHist("data_obs", histo_map["data_obs"]->GetName(), RooArgList(*x_var), histo_map["data_obs"].get()));

  ws->import(*data_obs.get());

  if (debug) ws->Print();
}


void CreateRooWorkspace::DoPlots() {
  /*
  & &&&&&&&&   &&&&&&&  &&&&&&&&  &&        &&&&&&&  &&&&&&&&  &&&&&&
  & &&     && &&     && &&     && &&       &&     &&    &&    &&    &&
  & &&     && &&     && &&     && &&       &&     &&    &&    &&
  & &&     && &&     && &&&&&&&&  &&       &&     &&    &&     &&&&&&
  & &&     && &&     && &&        &&       &&     &&    &&          &&
  & &&     && &&     && &&        &&       &&     &&    &&    &&    &&
  & &&&&&&&&   &&&&&&&  &&        &&&&&&&&  &&&&&&&     &&     &&&&&&
  */

  PlotBkgFit();
  for (auto syst: SystNames) PlotSignals(syst);
  PlotSgPars();
}


void CreateRooWorkspace::PlotBkgFit() {
  plotter.reset(x_var->frame(plot_lo,plot_hi));

  for (auto mode: Modes) {
    // if (mode!="bkg_pred" && mode!="DY_CR" && mode!="DY_SR" && mode!="WJets_CR" && mode!="WJets_SR") continue;
    // if (mode!="bkg_pred" && mode!="DY_CR") continue;
    // if (mode=="data") continue;

    std::unordered_map<std::string, RooHist*> hpull;
    std::unordered_map<std::string, RooHist*> hratio;
    std::unordered_map<std::string, double> chi2_map;

    TCanvas* c_bg = tdrDiCanvas(("Events"+mode).c_str(), plot_lo, plot_hi, plot_ylo, plot_yhi, doPlotRatio?0.8:-6, doPlotRatio?1.2:6, nameXaxis, nameYaxis, nameRatioaxis);
    c_bg->cd(1)->SetLogy(1);
    // plotter = x_var->frame(plot_lo,plot_hi);
    plotter.reset(x_var->frame(plot_lo,plot_hi));

    // if (mode=="DY_SR") rooHist_map["data"]->plotOn(plotter.get(), RooFit::LineColor(kRed+1), RooFit::DataError(RooAbsData::Poisson));
    if (mode=="data") rooHist_map[mode]->plotOn(plotter.get(),RooFit::DataError(RooAbsData::Poisson));
    else rooHist_map[mode]->plotOn(plotter.get(),RooFit::DataError(RooAbsData::SumW2), RooFit::Name(mode.c_str()));
    // rooHist_map[mode]->plotOn(plotter.get(),RooFit::DataError(RooAbsData::SumW2));
    // rooHist_map[mode]->plotOn(plotter.get());

    std::unique_ptr<TLegend>leg_bg(tdrLeg(0.40,0.55,0.7,0.85, 0.032));
    tdrHeader(leg_bg.get(), histo_map[mode]->GetName(), 12, 0.04, 42, kBlack, true);
    if (strcmp(mode.c_str(),"bkg_pred")==0) leg_bg->AddEntry(plotter.get()->findObject(mode.c_str()),"MC","p");

    std::unique_ptr<RooCurve> curve;
    std::unique_ptr<RooHist> hist_;
    for (auto const& [model,dofit] : doFits_map[mode] ) {
      if (dofit) {
        int color = Colors[model];
        // if (!doFtest && model.Contains("Exp") && model!="Exp_2" && model!="Exp_3" && model!="Exp_4" ) continue;
        Fits_map[mode][model]->plotOn(plotter.get(), RooFit::LineColor(color), RooFit::Range(fit_min[mode], fit_max[mode], kFALSE));
        // Fits_map[mode][model]->plotOn(plotter, RooFit::LineColor(color), RooFit::Range(fit_min[mode], fit_max[mode], kFALSE));

        // CREATE RATIO BETWEEN CURVE AND HIST
        double xstart, xstop, y;
        curve.reset((RooCurve*) ( (RooCurve*) plotter->findObject(0,RooCurve::Class()))->Clone(""));
        hist_.reset((RooHist*) ( (RooHist*) plotter->findObject(0,RooHist::Class()))->Clone(""));
        curve->GetPoint(0,xstart,y);
        curve->GetPoint(curve->GetN()-1,xstop,y);
        hratio[model] = new RooHist();
        for(int i=0 ; i<hist_->GetN() ; i++) {
          double x,point;
          hist_->GetPoint(i,x,point) ;
          if (x<xstart || x>xstop) continue ;
          double yy = point/ curve->interpolate(x);
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

        std::string isCR = FindInString("CR", mode) ? "CR" : "SR";

        ranges.at("CR").at(channel).at("show_lo");
        ranges.at("CR").at(channel).at("show_hi");

        double show_lo = ranges.at(isCR).at(channel).at("show_lo");
        double show_hi = ranges.at(isCR).at(channel).at("show_hi");
        if (FindInString("data", mode) && FindInString("CR", mode)) {
          show_hi = ranges.at("CR").at(channel).at("show_hi_data");
        }
        Fits_map[mode][model]->plotOn(plotter.get(), RooFit::LineColor(color), RooFit::Range(show_lo, show_hi, kFALSE), RooFit::FillColor(color), RooFit::VisualizeError(*FitRes_map[mode][model], 0.7), RooFit::FillStyle(3001), RooFit::VLines());
        Fits_map[mode][model]->plotOn(plotter.get(), RooFit::LineColor(color), RooFit::Range(show_lo, show_hi, kFALSE));
        hpull[model] = plotter->pullHist();
        hratio[model]->SetName(name);
        hratio[model]->SetMarkerColor(color);
        hratio[model]->SetLineColor(color);

        hpull[model]->SetName(name);
        hpull[model]->SetMarkerColor(color);
        hpull[model]->SetLineColor(color);

        // Remove points with meaningless pull
        for(int i=0 ; i<hpull[model]->GetN() ; i++) {
          double x,point;
          hpull[model]->GetPoint(i,x,point);
          if (hpull[model]->GetErrorYlow(i)==0) hpull[model]->SetPoint(i,x,-10);
        }
      }
    }

    plotter->Draw("same");

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

void CreateRooWorkspace::PlotSignals(std::string syst) {

  RooPlot* Pullplotter_all = x_var->frame(RooFit::Range(x_lo-10,x_hi+10), RooFit::Name("sgALL_pull"));
  std::unique_ptr<RooPlot> plotter_all(x_var->frame(RooFit::Range(x_lo,x_hi), RooFit::Name("signal_all")));
  std::unique_ptr<TCanvas> c_sg_all(tdrDiCanvas("signal", plot_lo, x_hi, 1*1e-03, 1, -6, 6, nameXaxis, nameYaxis, "Pull"));

  int i_mass = -1;
  for (const int & mass : MyMassPoints) {
    i_mass++;
    std::string SgName = GetSgName(mass,syst);

    RooPlot* Pullplotter = x_var->frame(RooFit::Range(plot_min[SgName],plot_max[SgName]),RooFit::Name(SgName.c_str()));

    if (debug) std::cout << "plotting mass:" << mass << "\t" << plot_min[SgName] << " -- " << plot_max[SgName] << "\tymax " << y_max[SgName] << '\n';

    std::unique_ptr<TLegend> leg_sig(tdrLeg(0.7,0.60,0.9,0.78, 0.05));

    if (isNominalSyst(syst)) {
      std::unique_ptr<RooPlot> plotter_syst(x_var->frame(RooFit::Range(plot_min[SgName],plot_max[SgName]), RooFit::Name("signal_syst")));
      std::unique_ptr<TLegend> leg_syst(tdrLeg(0.75,0.50,0.9,0.7, 0.035));
      TCanvas* c_sg_syst = tdrCanvas(SgName.c_str(), plot_min[SgName],plot_max[SgName], 2*1e-03, y_max[SgName],nameXaxis, nameYaxis);
      rooHist_map[SgName]->plotOn(plotter_syst.get(), RooFit::LineColor(kWhite), RooFit::MarkerColor(kWhite));
      for (auto syst: SystNames) {
        Fits_map[GetSgName(mass,syst)][FitSignal]->plotOn(plotter_syst.get(), RooFit::LineColor(Colors[syst]));
        TLine* line = new TLine(); line->SetLineColor(Colors[syst]);
        leg_syst->AddEntry(line,  syst.c_str() ,"l");
      }
      Fits_map[SgName][FitSignal]->plotOn(plotter_syst.get(), RooFit::VisualizeError(*FitRes_map[SgName][FitSignal], 1), RooFit::FillColor(kRed-7), RooFit::FillStyle(3001));
      Fits_map[SgName][FitSignal]->plotOn(plotter_syst.get(), RooFit::VisualizeError(*FitRes_map[SgName][FitSignal], 2), RooFit::FillColor(kRed-9), RooFit::FillStyle(3001));
      Fits_map[SgName][FitSignal]->plotOn(plotter_syst.get(), RooFit::LineColor(Colors["nominal"]), RooFit::Name((SgName+FitSignal+"fit").c_str()));
      leg_syst->Draw("same");
      plotter_syst->Draw("same");
      c_sg_syst->SaveAs(workingDir+"Fit_Sg"+TString(SgName).ReplaceAll("nominal","syst")+"_"+histFolder+extra_text+"."+plotting_mode);
      if (plotting_mode_2!="") c_sg_syst->SaveAs(workingDir+"Fit_Sg"+SgName+"_syst_"+histFolder+extra_text+"."+plotting_mode_2);
    }

    plotter.reset(x_var->frame(plot_min[SgName],plot_max[SgName]));
    TCanvas* c_sg = tdrDiCanvas(SgName.c_str(), plot_min[SgName],plot_max[SgName], 2*1e-03, y_max[SgName], -6, 6, nameXaxis, nameYaxis, "Pull");
    rooHist_map[SgName]->plotOn(plotter.get(), RooFit::LineColor(kBlack), RooFit::MarkerColor(kBlack));
    Fits_map[SgName][FitSignal]->plotOn(plotter.get(), RooFit::VisualizeError(*FitRes_map[SgName][FitSignal], 1), RooFit::FillColor(kRed-7), RooFit::FillStyle(3001));
    Fits_map[SgName][FitSignal]->plotOn(plotter.get(), RooFit::VisualizeError(*FitRes_map[SgName][FitSignal], 2), RooFit::FillColor(kRed-9), RooFit::FillStyle(3001));
    Fits_map[SgName][FitSignal]->plotOn(plotter.get(), RooFit::LineColor(Colors["nominal"]));
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
    SignalProperties << SgName+" signal number of events = " << nEventsSignal[SgName] << std::endl;
    SignalProperties << " chi2-ndf-pvalue "  << chi2 << " " << ndf << " " << pv << std::endl;
    SignalProperties << fitPars[SgName][0]->GetName() << "   param   " << fitPars[SgName][0]->getVal() << " " << fitPars[SgName][0]->getError() << std::endl;
    SignalProperties << fitPars[SgName][1]->GetName() << "   param   " << fitPars[SgName][1]->getVal() << " " << fitPars[SgName][1]->getError() << std::endl;
    SignalProperties << fitPars[SgName][2]->GetName() << "   param   " << fitPars[SgName][2]->getVal() << " " << fitPars[SgName][2]->getError() << std::endl;
    SignalProperties << fitPars[SgName][3]->GetName() << "   param   " << fitPars[SgName][3]->getVal() << " " << fitPars[SgName][3]->getError() << std::endl;


    nEventsSignal[SgName] /= CalculateFractionAreaPDF(Fits_map[SgName][FitSignal].get(), *x_var.get(), fit_min[SgName], fit_max[SgName]);
    SignalProperties << SgName+" signal number of events corrected = " << nEventsSignal[SgName] << std::endl;

    std::unique_ptr<TPaveText> pave(new TPaveText(0.7,0.3,0.8,0.6,"NDC"));
    pave->SetBorderSize(0); pave->SetTextSize(0.03); pave->SetLineColor(1); pave->SetLineStyle(1);
    pave->SetLineWidth(2); pave->SetFillColor(0); pave->SetFillStyle(0);
    pave->AddText(TString::Format("Fit range = [%.0f,%.0f]", fit_min[SgName],fit_max[SgName]));
    pave->AddText(TString::Format("Integral%.1f", nEventsSignal[SgName]));
    pave->AddText(TString::Format("#chi^{2}/n.d.f. = %.1f", chi2));
    pave->AddText(TString::Format("p-value = %.1f", pv));
    pave->AddText(TString::Format("#mu = %2.3f +- %2.3f", fitPars[SgName][0]->getVal(),fitPars[SgName][0]->getError()));
    pave->AddText(TString::Format("#sigma = %2.3f +- %2.3f", fitPars[SgName][1]->getVal(),fitPars[SgName][1]->getError()));
    pave->AddText(TString::Format("#alpha = %2.3f +- %2.3f", fitPars[SgName][2]->getVal(),fitPars[SgName][2]->getError()));
    pave->AddText(TString::Format("k = %2.3f +- %2.3f", fitPars[SgName][3]->getVal(),fitPars[SgName][3]->getError()));
    SgPars["Masses"].at(i_mass)       = (double)mass;
    SgPars["nevents"].at(i_mass)      = nEventsSignal[SgName];
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

    Fits_map[SgName][FitSignal]->plotOn(plotter_all.get(), RooFit::LineColor(kRed), RooFit::Name((SgName+FitSignal).c_str()));
    c_sg->cd(1);

    leg_sig->AddEntry(plotter_all.get()->findObject((SgName+FitSignal).c_str()),"Signal sample","p");
    leg_sig->AddEntry(plotter_all.get()->findObject((SgName+FitSignal).c_str()),"Fit","l");
    leg_sig->Draw("same");

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
  &  &&&&&&  &&     && &&&&&&&&  &&&&&&  &&    &&    &&&&&&&&  &&        &&&&&&&  &&&&&&&&  &&&&&&
  & &&    && &&     && &&       &&    && &&   &&     &&     && &&       &&     &&    &&    &&    &&
  & &&       &&     && &&       &&       &&  &&      &&     && &&       &&     &&    &&    &&
  & &&       &&&&&&&&& &&&&&&   &&       &&&&&       &&&&&&&&  &&       &&     &&    &&     &&&&&&
  & &&       &&     && &&       &&       &&  &&      &&        &&       &&     &&    &&          &&
  & &&    && &&     && &&       &&    && &&   &&     &&        &&       &&     &&    &&    &&    &&
  &  &&&&&&  &&     && &&&&&&&&  &&&&&&  &&    &&    &&        &&&&&&&&  &&&&&&&     &&     &&&&&&
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
  // histo_map["MC_SR"]->SetLineColor(kGreen+1);
  histo_map["bkg_pred"]->Draw("same");
  histo_map["data"]->Draw("same");
  histo_map["MC_SR"]->Draw("same");
  std::unique_ptr<TLegend> leg_inputs(tdrLeg(0.50,0.70,0.78,0.9, 0.03));
  leg_inputs->AddEntry(histo_map["bkg_pred"].get(), Form("dofit: %s", histo_map["bkg_pred"]->GetName()) ,"lp");
  leg_inputs->AddEntry(histo_map["data"].get(), Form("extract obs. Limits: %s", histo_map["data"]->GetName()) ,"lp");
  leg_inputs->AddEntry(histo_map["MC_SR"].get(), Form("extract exp. Limits: %s", histo_map["MC_SR"]->GetName()) ,"lp");
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
      // DataCard << name+" rateParam " <<  f1->GetParameter(0)*CalculateIntegral(histo_map[name+"_CR"],fit_lo_SR,fit_hi_SR,doBinWidth)/CalculateIntegral(histo_map[name+"_SR"],fit_lo_SR,fit_hi_SR,doBinWidth) << std::endl;
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
  & &&&& &&    && &&&&&&&&  &&     && &&&&&&&& &&&&&&&&     &&&    &&&&&&&&    &&&     &&&&&&     &&&    &&&&&&&&  &&&&&&&&   &&&&&&
  &  &&  &&&   && &&     && &&     &&    &&    &&     &&   && &&      &&      && &&   &&    &&   && &&   &&     && &&     && &&    &&
  &  &&  &&&&  && &&     && &&     &&    &&    &&     &&  &&   &&     &&     &&   &&  &&        &&   &&  &&     && &&     && &&
  &  &&  && && && &&&&&&&&  &&     &&    &&    &&     && &&     &&    &&    &&     && &&       &&     && &&&&&&&&  &&     &&  &&&&&&
  &  &&  &&  &&&& &&        &&     &&    &&    &&     && &&&&&&&&&    &&    &&&&&&&&& &&       &&&&&&&&& &&   &&   &&     &&       &&
  &  &&  &&   &&& &&        &&     &&    &&    &&     && &&     &&    &&    &&     && &&    && &&     && &&    &&  &&     && &&    &&
  & &&&& &&    && &&         &&&&&&&     &&    &&&&&&&&  &&     &&    &&    &&     &&  &&&&&&  &&     && &&     && &&&&&&&&   &&&&&&
  */


  for (auto mode: Modes) {
    for (auto const& [model,dofit] : doFits_map[mode] ) {
      if (!dofit || (model.find("Exp")!= std::string::npos && model!="Exp_2" && model!="Exp_3" && model!="Exp_4") ) continue;
      // double frac = CalculateFractionAreaPDF(Fits_map[mode][model].get(), *x_var.get(), fit_lo_SR, fit_hi_SR);
      double frac = CalculateFractionAreaPDF(Fits_map[mode][model].get(), *x_var.get(), fit_min[mode], fit_max[mode]);
      if (debug) std::cout << mode+"_"+model+unique_name << " integral " << nEventsSR/frac  << " " << frac << " " << CalculateFractionAreaPDF(Fits_map[mode][model].get(), *x_var.get(), fit_min[mode], fit_max[mode]) << std::endl;
      DataCard  << mode+"_"+model+unique_name << " integral " << nEventsSR/frac  << " " << frac << std::endl;


      std::unique_ptr<RooArgSet> model_params(Fits_map[mode][model]->getParameters(*x_var));
      TString name; int from = 0;
      while (TString(model_params->contentsString()).Tokenize(name, from, ",")) {
        RooRealVar* par = dynamic_cast<RooRealVar*>(model_params->find(name));
        DataCard << Form("%s%sparam %.3f %.3f", name.Data(), std::string(60-name.Length(),' ').c_str(), par->getVal(), par->getError()) <<std::endl;
      }
    }
  }

  for (auto syst: SystNames) {
    for (const int & mass : MyMassPoints) {
      std::string SgName = GetSgName(mass,syst);
      DataCard << SgName+" signal number of events = " << nEventsSignal[SgName] <<""<<std::endl;
      for (unsigned int i = 0; i < fitPars[SgName].size(); i++) {
        DataCard << SgName << " " << fitPars[SgName][i]->GetName() << "   param   " <<fitPars[SgName][i]->getVal() << " " << fitPars[SgName][i]->getError() <<  std::endl;
      }
    }
  }

}




void CreateRooWorkspace::Process() {
  /*
  & &&&&&&&&  &&&&&&&&   &&&&&&&   &&&&&&  &&&&&&&&  &&&&&&   &&&&&&
  & &&     && &&     && &&     && &&    && &&       &&    && &&    &&
  & &&     && &&     && &&     && &&       &&       &&       &&
  & &&&&&&&&  &&&&&&&&  &&     && &&       &&&&&&    &&&&&&   &&&&&&
  & &&        &&   &&   &&     && &&       &&             &&       &&
  & &&        &&    &&  &&     && &&    && &&       &&    && &&    &&
  & &&        &&     &&  &&&&&&&   &&&&&&  &&&&&&&&  &&&&&&   &&&&&&
  */


  LoadFiles();//OK
  PrepocessHistos();//OK
  NormaliseData(); //OK
  DoRebin(); //TODO put it in the correct place
  InitializePDFs(); //OK
  CreateRooDataHist(); //OK
  DoFits(); //OK
  DoPlots();
  ImportToWorkspace(); //OK
  InputDatacards();
};


int main(int argc, char** argv){
  /*
  & &&     &&    &&&    &&&& &&    &&
  & &&&   &&&   && &&    &&  &&&   &&
  & &&&& &&&&  &&   &&   &&  &&&&  &&
  & && &&& && &&     &&  &&  && && &&
  & &&     && &&&&&&&&&  &&  &&  &&&&
  & &&     && &&     &&  &&  &&   &&&
  & &&     && &&     && &&&& &&    &&
  */

  gErrorIgnoreLevel = kFatal;
  RooMsgService::instance().setSilentMode(true);
  RooMsgService::instance().setStreamStatus(1,false);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::MsgTopic::Plotting);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::MsgTopic::InputArguments);
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

  bool isHbb = false;
  std::vector<std::string> histFolders = { "DeepAk8_H4qvsQCD_massdep", "DeepAk8_ZHccvsQCD_MD",
  "DeepAk8_HccvsQCD_MD", "DeepAk8_H4qvsQCD_MD", "DeepAk8_H4qvsQCD_massdep_HccvsQCD_MD",
  "DeepAk8_H4qvsQCD", "DeepAk8_HccvsQCD", "DeepAk8_ZHccvsQCD", "DeepAk8_H4qvsQCD_massdep_HccvsQCD",
  "DeepAk8_H4qvsQCD_massdep_ZHccvsQCD", "DeepAk8_H4qvsQCD_massdep_ZHccvsQCD_MD", "tau42",
  "DeepAk8_HccvsQCD2", "DeepAk8_ZHccvsQCD_MD2"};

  if (isHbb) histFolders = {"tau21" };
  std::vector<std::string> collections = {"Puppi"};
  std::vector<std::string> channels = {"muonchannel", "electronchannel", "invisiblechannel"};
  std::vector<std::string> years = {"2016", "2017", "2018", "RunII"};


  std::unique_ptr<CreateRooWorkspace> roo;


  if (argc>1) {
    std::string histFolder, channel, collection, year, min="", max="";
    for (int i = 1; i < argc; i++) {
      if (std::find(collections.begin(), collections.end(), argv[i]) != collections.end() ) collection = argv[i];
      if (std::find(histFolders.begin(), histFolders.end(), argv[i]) != histFolders.end() ) histFolder = argv[i];
      if (std::find(channels.begin(), channels.end(), argv[i]) != channels.end() ) channel = argv[i];
      if (std::find(years.begin(), years.end(), argv[i]) != years.end() ) year = argv[i];
    }
    if (argc>5) min = argv[5];
    if (argc>6) max = argv[6];

    std::cout << argc << " " << year << " " << collection << " " << channel << " " << histFolder << " " << min << " " << max << '\n';
    roo.reset(new CreateRooWorkspace(year,collection, channel, histFolder, min, max));
    roo->Process();

  } else {
    std::cout << "more " << '\n';
    for (std::string year: years) {
      for (std::string collection: collections) {
        for (std::string channel: channels) {
          for (std::string histFolder: histFolders) {
            // roo.reset(new CreateRooWorkspace(year,collection, channel, histFolder));
            roo->Process();
          }
        }
      }
    }
  }

  return 0;

}
