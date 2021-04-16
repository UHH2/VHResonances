#include "Plotter.hpp"


std::string Split(std::string word, int index, std::string split="/") {
  TObjArray *tx = TString(word).Tokenize(split);
  return ((TObjString *)(tx->At(index)))->String().Data();
}

SampleInfo::SampleInfo(std::string uniqueName, std::string fileName, double weight, Type type, std::string legName, int color, int linestyle) {
  if (false) PrintGreen("--> SampleInfo");

  SetUniqueName(uniqueName);
  SetLegName(legName);
  SetFileName(fileName);
  SetColor(color);
  SetLinestyle(linestyle);
  SetWeight(weight);
  switch(type) {
    case 0  : SetIsMC(); break;
    case 1  : SetIsSignal(); break;
    case 2  : SetIsData(); break;
  };

  if (false) {
    PrintInfo("UniqueName",GetUniqueName());
    PrintInfo("LegName",GetLegName());
    PrintInfo("FileName",GetFileName());
    PrintInfo("Color",GetColor());
    PrintInfo("Linestyle",GetLinestyle());
    PrintInfo("Weight",GetWeight());
    PrintInfo("IsMC",GetIsMC());
    PrintInfo("IsSignal",GetIsSignal());
    PrintInfo("IsData",GetIsData());
    PrintInfo("IsToStack",GetIsToStack());
  }
}

void Plotter::AddSample(std::string uniqueName, std::string fileName, double weight, Type type, std::string legName, int color, int linestyle) {
  Samples[uniqueName].reset(new SampleInfo(uniqueName, fileName+"_"+year+"_noTree.root", weight, type, legName, color, linestyle));
  histnames.push_back(uniqueName);
  if (Samples[uniqueName]->GetIsToStack()) histstacknames.push_back(uniqueName);
}

void Plotter::SetEnv() {
  if (debug) PrintGreen("--> SetEnv");

  user            = std::getenv("USER");
  Path_ANALYSIS   = std::getenv("CMSSW_BASE"); Path_ANALYSIS += "/src/UHH2/VHResonances/";
  Path_STORAGE    = "/nfs/dust/cms/user/"+user+"/WorkingArea/File/Analysis/";
  PrefixrootFile  = "uhh2.AnalysisModuleRunner.";
  relPath = year+"/"+module+"/"+collections+"/"+channel+"/nominal/";
  inputdir = Path_STORAGE+relPath;
  outputdir = Path_ANALYSIS+"Analysis/PlotterSteer/plots/";
  gSystem->Exec(("mkdir -p "+outputdir).c_str());

  if (debug) {
    PrintInfo("user",user);
    PrintInfo("Path_ANALYSIS",Path_ANALYSIS);
    PrintInfo("Path_STORAGE",Path_STORAGE);
    PrintInfo("PrefixrootFile",PrefixrootFile);
    PrintInfo("inputdir",inputdir);
    PrintInfo("outputdir",outputdir);
  }

  if(GetChannel()!="muonchannel" && FindInVector(systematics,"MuonScale")) systematics.erase(systematics.begin()+distance(systematics.begin(), std::find(systematics.begin(), systematics.end(), "MuonScale")));
  if (module=="SignalRegion") {
    for (const auto& syst: SystematicsScale) {
      if (!FindInString("muon",GetChannel()) && FindInString("tracking",syst)) continue;
      if (!FindInString("muon",GetChannel()) && FindInString("isolation",syst)) continue;
      if (FindInString("invisible",GetChannel()) && FindInString("id",syst)) continue;
      if (FindInString("invisible",GetChannel()) && FindInString("trigger",syst)) continue;
      if (FindInString("invisible",GetChannel()) && FindInString("reco",syst)) continue;
      systematics.push_back(syst);
    }
  }

  if (!plotSyst) systematics.clear();

  WIP();
  SetYear(TString::Format("%.1f fb^{-1}", lumi_map.at(year).at("lumi_fb")));

  lumi_unc = lumi_map.at(year).at("uncertainty");

  stack.reset(new THStack(("stack"+relPath+hname).c_str(), "" ));
  stack_unc.reset(new THStack(("stack_unc"+relPath+hname).c_str(), "" ));
  for (const auto& syst: systematics){
    if (debug) PrintInfo("syst",syst);
    for (const auto& var: {"up","down"}){
      stacksts[syst+var].reset(new THStack(("stack"+relPath+syst+var).c_str(), "" ));
    }
  }
}


void Plotter::LoadHists(){
  if (debug) PrintGreen("--> LoadHists");

  for (const auto& sample: histnames){
    if (debug) PrintGreen(inputdir+PrefixrootFile+Samples[sample]->GetFileName());
    if (debug) PrintGreen(hname);
    std::unique_ptr<TFile> f_; f_.reset(new TFile((inputdir+PrefixrootFile+Samples[sample]->GetFileName()).c_str()));
    histos[sample].reset((TH1D*)f_->Get(hname.c_str()));
    histos[sample]->SetDirectory(0);
    histos[sample]->Scale(Samples[sample]->GetWeight());
    histos[sample]->Rebin(GetRebin());
    f_->Close();
    if(binwidth==0) binwidth = histos[sample]->GetBinWidth(1);
    if(nbins==1) nbins = histos[sample]->GetNbinsX();
    else if(nbins != histos[sample]->GetNbinsX()) {PrintGreen(nbins);PrintGreen(histos[sample]->GetNbinsX()); throw std::runtime_error("Hists have different binning. Be carreful.");}
    if (!Samples[sample]->GetIsToStack()) continue;
    for (const auto& syst: systematics){
      if(syst=="lumi") {
        for (const auto& var: {"up","down"}){
          uncertainties[sample][syst+var].reset((TH1D*)(histos[sample]->Clone((sample+"_"+syst+var).c_str())));
          uncertainties[sample][syst+var]->SetDirectory(0);
          uncertainties[sample][syst+var]->Scale(Samples[sample]->GetWeight());
        }
        uncertainties[sample][syst+"up"]->Scale(1+lumi_unc/100);
        uncertainties[sample][syst+"down"]->Scale(1-lumi_unc/100);
      } else if(syst=="murmuf") {
        TString fname_ = inputdir+PrefixrootFile+Samples[sample]->GetFileName();
        if (debug) PrintGreen(syst);
        if (debug) PrintGreen(fname_);
        if (debug) PrintGreen(hname);
        std::unique_ptr<TFile> f_sys; f_sys.reset(new TFile(fname_));
        uncertainties[sample][syst+"up"].reset((TH1D*)((TH1D*)f_sys->Get(hname.c_str()))->Clone((syst+"up").c_str()));
        uncertainties[sample][syst+"down"].reset((TH1D*)((TH1D*)f_sys->Get(hname.c_str()))->Clone((syst+"down").c_str()));
        uncertainties[sample][syst+"up"]->SetDirectory(0);
        uncertainties[sample][syst+"down"]->SetDirectory(0);
        uncertainties[sample][syst+"up"]->Scale(Samples[sample]->GetWeight());
        uncertainties[sample][syst+"down"]->Scale(Samples[sample]->GetWeight());
        uncertainties[sample][syst+"up"]->Rebin(GetRebin());
        uncertainties[sample][syst+"down"]->Rebin(GetRebin());

        std::unordered_map<std::string, std::unique_ptr<TH1F> > histos_temp;
        for (auto var: Var_murmuf) {
          TString hname_ = TString(hname).ReplaceAll(Split(hname,0,"_"), Split(hname,0,"_")+"_"+syst+"_"+var);
          histos_temp[var].reset((TH1F*)f_sys->Get(hname_));
          histos_temp[var]->SetDirectory(0);
          histos_temp[var]->Scale(Samples[sample]->GetWeight());
          histos_temp[var]->Rebin(GetRebin());
        }

        for (int bin = 0; bin < uncertainties[sample][syst+"up"]->GetNbinsX()+1; bin++) {
          std::vector<double> y_vals;
          for (const auto& x : histos_temp) y_vals.push_back(x.second->GetBinContent(bin));
          uncertainties[sample][syst+"up"]->SetBinContent(bin, *max_element(y_vals.begin(), y_vals.end()));
          uncertainties[sample][syst+"down"]->SetBinContent(bin, *min_element(y_vals.begin(), y_vals.end()));
        }
        f_sys->Close();
      } else if(syst=="NNPDF") {
        TString fname_ = inputdir+PrefixrootFile+Samples[sample]->GetFileName();
        if (debug) PrintGreen(syst);
        if (debug) PrintGreen(fname_);
        if (debug) PrintGreen(hname);
        std::unique_ptr<TFile> f_sys; f_sys.reset(new TFile(fname_));
        uncertainties[sample][syst+"up"].reset((TH1D*)((TH1D*)f_sys->Get(hname.c_str()))->Clone((syst+"up").c_str()));
        uncertainties[sample][syst+"down"].reset((TH1D*)((TH1D*)f_sys->Get(hname.c_str()))->Clone((syst+"down").c_str()));
        uncertainties[sample][syst+"up"]->SetDirectory(0);
        uncertainties[sample][syst+"down"]->SetDirectory(0);
        uncertainties[sample][syst+"up"]->Scale(Samples[sample]->GetWeight());
        uncertainties[sample][syst+"down"]->Scale(Samples[sample]->GetWeight());
        uncertainties[sample][syst+"up"]->Rebin(GetRebin());
        uncertainties[sample][syst+"down"]->Rebin(GetRebin());

        std::unordered_map<std::string, std::unique_ptr<TH1F> > histos_temp;
        for (int var = 0; var < PDF_variations ; var++) {
          TString hname_ = TString(hname).ReplaceAll(Split(hname,0,"_"), Split(hname,0,"_")+"_"+syst+"_"+std::to_string(var));
          histos_temp[std::to_string(var)].reset((TH1F*)f_sys->Get(hname_));
          histos_temp[std::to_string(var)]->SetDirectory(0);
          histos_temp[std::to_string(var)]->Scale(Samples[sample]->GetWeight());
          histos_temp[std::to_string(var)]->Rebin(GetRebin());
          //TODO  normalize?
        }

        for (int bin = 0; bin < uncertainties[sample][syst+"up"]->GetNbinsX()+1; bin++) {
          std::vector<double> y_vals;
          for (const auto& x : histos_temp) y_vals.push_back(x.second->GetBinContent(bin));
          double mean = 0.; double stdev=0;
          for (const auto& x : y_vals) mean += x;
          mean /= PDF_variations;
          for (const auto& x : y_vals) stdev += (x-mean)*(x-mean);
          stdev = TMath::Sqrt(stdev/PDF_variations);
          uncertainties[sample][syst+"up"]->SetBinContent(bin, uncertainties[sample][syst+"up"]->GetBinContent(bin) + stdev );
          uncertainties[sample][syst+"down"]->SetBinContent(bin, uncertainties[sample][syst+"down"]->GetBinContent(bin) - stdev );
        }
        f_sys->Close();
      } else {
        for (const auto& var: {"up","down"}){
          TString fname_ = inputdir+PrefixrootFile+Samples[sample]->GetFileName();
          TString hname_ = hname;
          if (!isNominalFolder(syst)) fname_ = fname_.ReplaceAll("nominal",syst+"_"+var);
          else if (!isNominalSyst(syst)) hname_ = hname_.ReplaceAll(Split(hname,0,"_"), Split(hname,0,"_")+"_"+syst+"_"+var);
          if (debug) PrintGreen(syst);
          if (debug) PrintGreen(fname_);
          if (debug) PrintGreen(hname_);
          std::unique_ptr<TFile> f_sys; f_sys.reset(new TFile(fname_));
          uncertainties[sample][syst+var].reset((TH1D*)f_sys->Get(hname_));
          uncertainties[sample][syst+var]->SetDirectory(0);
          uncertainties[sample][syst+var]->Scale(Samples[sample]->GetWeight());
          uncertainties[sample][syst+var]->Rebin(GetRebin());
          f_sys->Close();
        }
      }
    }
  }
}


void Plotter::FindRanges(){
  if (debug) PrintGreen("--> FindRanges");

  for (const auto& sample: histnames){

    for(int i =1; i<nbins; i++) {
      bool reset = histos[sample]->GetBinContent(i)<=0;
      reset += (IsSetXmin() && histos[sample]->GetXaxis()->GetBinUpEdge(i)<=xmin);
      reset += (IsSetXmax() && histos[sample]->GetXaxis()->GetBinLowEdge(i)>=xmax);
      reset += (IsSetYmin() && histos[sample]->GetBinContent(i)<=ymin);
      if(reset) {
        histos[sample]->SetBinContent(i,0);
        histos[sample]->SetBinError(i,0);
      } else {
        if (!IsSetXmin()) xmin = std::min(xmin, histos[sample]->GetXaxis()->GetBinLowEdge(i));
        if (!IsSetXmax()) xmax = std::max(xmin, histos[sample]->GetXaxis()->GetBinUpEdge(i));
        if (!IsSetYmin()) ymin = std::min(ymin, histos[sample]->GetBinContent(i));
        if (!IsSetYmax()) ymax = std::max(ymax, histos[sample]->GetBinContent(i));
      }
    }
  }

  if (!FindInString("sum",hname)) xmin -= 1.5*binwidth;
  if (!FindInString("sum",hname)) xmax += 1.5*binwidth;
  if (!IsSetYmin()) ymin *= 0.90;
  if (!IsSetYmax()) ymax *= 120;
  if (!IsSetYmin()) ymin = std::max(ymin,0.0011);
  // ymax *= 1.20; //ok for non log
  if (debug) {
    PrintInfo("xmin", xmin);
    PrintInfo("xmax", xmax);
    PrintInfo("ymin", ymin);
    PrintInfo("ymax", ymax);
  }
}


void Plotter::MakeStackHist(){
  if (debug) PrintGreen("--> MakeStackHist");
  std::reverse(histstacknames.begin(),histstacknames.end());
  for (const auto& sample: histstacknames) {
    if (!Samples[sample]->GetIsToStack()) continue;
    histos[sample]->SetFillColor(Samples[sample]->GetColor());
    histos[sample]->SetLineWidth(0);
    stack->Add(histos[sample].get());
    stack_unc->Add(histos[sample].get());
    for (const auto& systematic: systematics){
      for (const auto& var: {"up","down"}){
        stacksts[systematic+var]->Add(uncertainties[sample][systematic+var].get());
      }
    }
  }
  h_err.reset(new TH1F(*(TH1F*)(stack->GetStack()->Last())));
  if (!plotSyst) return;
  h_err_syst.reset(new TH1F(*(TH1F*)(stack_unc->GetStack()->Last())));
  for (int bin = 0; bin < h_err_syst->GetNbinsX()+1; bin++) {
    double err_up = h_err_syst->GetBinError(bin)*h_err_syst->GetBinError(bin);
    double err_down = h_err_syst->GetBinError(bin)*h_err_syst->GetBinError(bin);
    for (const auto& x : stacksts) {
      TH1F* h_ = new TH1F(*(TH1F*)(x.second->GetStack()->Last()));
      double err = TMath::Abs(h_->GetBinContent(bin)-h_err_syst->GetBinContent(bin));
      if (FindInString("down",x.first)) err_down += err*err;
      else err_up += err*err;
    }
    h_err_syst->SetBinError(bin, (err_up>err_down)? TMath::Sqrt(err_up) : TMath::Sqrt(err_down));
  }
}


void Plotter::MakePlot(){
  if (debug) PrintGreen("--> MakePlot");

  TCanvas* canvas = tdrDiCanvas(("stack"+relPath+hname).c_str(), xmin, xmax, ymin, ymax, 0.3, 1.7, xTitle, "Events", "Data/Pred.", false);
  TLegend* leg = tdrLeg(0.6, 0.9 - 0.035*histnames.size(), 0.9, 0.9, 0.035, 42);

  canvas->cd(1)->SetLogy(1);
  stack->Draw("hist same");

  // draw uncertainty of stack
  tdrDraw(h_err.get(), "E2", 0, kGray+1, 0, kGray+1, 3005, kGray+1);
  if(plotSyst) tdrDraw(h_err_syst.get(), "E2", 0, kGray+2, 0, kGray+2, 3005, kGray+2);

  std::string name_hist_ratio;
  for (const auto& sample: histnames){
    if (Samples[sample]->GetIsToStack()) {
      leg->AddEntry(histos[sample].get(), (Samples[sample]->GetLegName()).c_str(),"f");
      continue;
    }
    if (Samples[sample]->GetIsData()) name_hist_ratio = sample;
    std::string opt = Samples[sample]->GetIsData()? "P E": "hist";
    histos[sample]->SetLineWidth(3);
    // histos[sample]->SetMarkerSize(1.15);
    tdrDraw(histos[sample].get(), opt, kFullDotLarge,  Samples[sample]->GetColor(), Samples[sample]->GetLinestyle(), Samples[sample]->GetColor(), 0, Samples[sample]->GetColor());
    opt = Samples[sample]->GetIsData()? "lp": "l";
    leg->AddEntry(histos[sample].get(), (Samples[sample]->GetLegName()).c_str(),opt.c_str());
  }

  fixOverlay();
  canvas->cd(2);

  TLegend* leg_ratio = tdrLeg(0.17, 0.75, 0.49, 0.89, 0.08);
  leg_ratio->SetNColumns(2);

  TH1F* h_ratio  = (TH1F*)(histos[name_hist_ratio]->Clone("ratio"));
  TH1F* h_ratiostat = (TH1F*)(stack->GetStack()->Last()->Clone("ratiostat"));
  std::unordered_map<std::string, std::unique_ptr<TH1F>> h_ratiosyst;

  for (const auto& systematic: systematics){
    for (const auto& var: {"up","down"}){
      h_ratiosyst[systematic+var].reset((TH1F*)(stacksts[systematic+var]->GetStack()->Last()->Clone(("ratiostat"+systematic+var).c_str())));
    }
  }

  TGraphAsymmErrors* h_ratiotot = new TGraphAsymmErrors();

  for(int j=1; j<nbins+1; j++){
    double x_center = h_ratio->GetXaxis()->GetBinCenter(j);
    double val_n = h_ratio->GetBinContent(j);
    double err_n = h_ratio->GetBinError(j);
    double val_d = h_ratiostat->GetBinContent(j);
    double err_d = h_ratiostat->GetBinError(j);
    if(val_d > 0){
      h_ratio->SetBinContent(j, val_n/val_d);
      h_ratio->SetBinError(j, err_n/val_d);
      h_ratiostat->SetBinContent(j, 1.);
      h_ratiostat->SetBinError(j, err_d/val_d);
      double systup = 0;
      double systdown = 0;
      for (const auto& systematic: systematics){
        double x;
        x = fabs(h_ratiosyst[systematic+"up"]->GetBinContent(j)-val_d);
        if (!isnan(x)) systup += x*x;
        x = fabs(h_ratiosyst[systematic+"down"]->GetBinContent(j)-val_d);
        if (!isnan(x)) systdown += x*x;
      }
      h_ratiotot->SetPoint(j, x_center, 1);
      double x_up = h_ratio->GetXaxis()->GetBinUpEdge(j)-x_center;
      double x_low = x_center - h_ratio->GetXaxis()->GetBinLowEdge(j);
      double y_up = TMath::Sqrt(systup+err_d*err_d)/val_d;
      double y_low = TMath::Sqrt(systdown+err_d*err_d)/val_d;
      h_ratiotot->SetPointError(j, x_low, x_up, y_low, y_up);

    }
    else{
      h_ratio->SetBinContent(j, 0.);
      h_ratio->SetBinError(j, 0.);
      h_ratiostat->SetBinContent(j, 0.);
      h_ratiostat->SetBinError(j, 0.);
      h_ratiotot->SetPoint(j, x_center, 0);
      h_ratiotot->SetPointError(j, x_center - h_ratio->GetXaxis()->GetBinLowEdge(j), h_ratio->GetXaxis()->GetBinUpEdge(j)-x_center, 0,0);
    }
  }

  if (plotSyst) tdrDraw(h_ratiotot, "E2", 0, kGray+2, 0, kGray+2, 1000, kGray+2);
  tdrDraw(h_ratiostat, "E2", 0, kGray+1, 0, kGray+1, 1000, kGray+1);
  tdrDraw(h_ratio, "E", 20, kBlack, 0, kBlack, 1000, kBlack);
  leg_ratio->AddEntry(h_ratiostat,"Stat.","f");
  if (plotSyst) leg_ratio->AddEntry(h_ratiotot,"Stat. #oplus Syst.","f");
  fixOverlay();

  canvas->SaveAs((outputdir+year+"_"+module+"_"+collections+"_"+channel+"_nominal_"+Split(hname,0)+"_"+Split(hname,1)+(plotSyst?"_syst":"")+".pdf").c_str());

}


Plotter::Plotter(std::string module_, std::string hname_, std::string year_, std::string channel_, std::string collections_) {
  if (debug) PrintGreen("--> Start");
  module = module_;
  hname = hname_;
  year = year_;
  channel = channel_;
  collections = collections_;

  if (debug) {
    PrintInfo("module",module);
    PrintInfo("hname",hname);
    PrintInfo("year",year);
    PrintInfo("channel",channel);
    PrintInfo("collections",collections);
  }
  SetEnv();
}

void Plotter::Process(){
  if (year!="RunII" && plotSyst) return;
  if (debug) PrintGreen("--> Process");
  LoadHists();
  FindRanges();
  MakeStackHist();
  MakePlot();
}

Plotter::~Plotter(){ if (debug) PrintGreen("--> End");}


/*
&  &&&&&&   &&&&&&&& &&    && &&&&&&&& &&&&&&&&     &&&    &&           &&&&&&  &&&&&&&& &&&&&&&& &&&&&&&& &&&& &&    &&  &&&&&&    &&&&&&
& &&    &&  &&       &&&   && &&       &&     &&   && &&   &&          &&    && &&          &&       &&     &&  &&&   && &&    &&  &&    &&
& &&        &&       &&&&  && &&       &&     &&  &&   &&  &&          &&       &&          &&       &&     &&  &&&&  && &&        &&
& &&   &&&& &&&&&&   && && && &&&&&&   &&&&&&&&  &&     && &&           &&&&&&  &&&&&&      &&       &&     &&  && && && &&   &&&&  &&&&&&
& &&    &&  &&       &&  &&&& &&       &&   &&   &&&&&&&&& &&                && &&          &&       &&     &&  &&  &&&& &&    &&        &&
& &&    &&  &&       &&   &&& &&       &&    &&  &&     && &&          &&    && &&          &&       &&     &&  &&   &&& &&    &&  &&    &&
&  &&&&&&   &&&&&&&& &&    && &&&&&&&& &&     && &&     && &&&&&&&&     &&&&&&  &&&&&&&&    &&       &&    &&&& &&    &&  &&&&&&    &&&&&&
*/

void SetData(Plotter* plotter){
  std::string ch = plotter->GetChannel();
  std::string name;
  if (ch=="muonchannel") name = "SingleMuon";
  if (ch=="electronchannel") name = "SingleElectron";
  if (ch=="invisiblechannel") name = "MET";
  plotter->AddSample("Data",        "DATA.DATA_"+name,   1.0,   typeData,   "Data",     kBlack,    kSolid);
}

void SetMC(Plotter* plotter){
  plotter->AddSample("DY",          "MC.MC_DY",               1.0,   typeMC,     "DY",       kOrange-2, kSolid);
  if (plotter->GetChannel()=="invisiblechannel") {
    plotter->AddSample("WJets",     "MC.MC_WJets",            1.0,   typeMC,     "WJets",    kAzure+7,  kSolid);
  }
  plotter->AddSample("VV",          "MC.MC_VV",               1.0,   typeMC,     "VV",       kGreen+2,  kSolid);
  plotter->AddSample("TTbar",       "MC.MC_TTbar",            1.0,   typeMC,     "TTbar",    kRed+1,    kSolid);
}

void SetSignal(Plotter* plotter){
  std::string name = "ZprimeToZH";
  if (plotter->GetChannel()=="invisiblechannel") name = "ZprimeToZH_inv";
  plotter->AddSample("Zprime2TeV",  "MC.MC_"+name+"_M2000", 0.001, typeSignal, "Z' 2 TeV", kBlack,    kSolid);
  plotter->AddSample("Zprime3TeV",  "MC.MC_"+name+"_M3000", 0.001, typeSignal, "Z' 3 TeV", kBlack,    kDashed);
  plotter->AddSample("Zprime4TeV",  "MC.MC_"+name+"_M4000", 0.001, typeSignal, "Z' 4 TeV", kBlack,    kDotted);
}



void PlotDistribution(std::string ch, std::string year, std::string module, std::string hname, double xmin, double xmax, double ymin, double ymax, std::string xtitle) {
  std::unique_ptr<Plotter> plotter;
  plotter.reset(new Plotter(module, hname, year, ch));
  SetData(plotter.get());
  SetMC(plotter.get());
  SetSignal(plotter.get());
  plotter->SetXRange(xmin,xmax);
  plotter->SetYRange(ymin, ymax);
  plotter->SetXTitle(xtitle);
  if(FindInString("_pt", hname)) plotter->SetRebin(5);
  plotter->Process();
}


void PlotZprimeMass(std::string ch="muonchannel", std::string year="RunII", std::string module="Selection") {
  std::string hname="ZprimeCandidate_ScaleFactors/Zprime_mass_rebin100";
  if (module=="SignalRegion") hname = TString(hname).ReplaceAll("ScaleFactors","Selection").Data();
  bool isInv = ch=="invisiblechannel";
  if (isInv) hname = TString(hname).ReplaceAll("mass", "mass_transversal");
  double ymax = isInv? 7.01e6 :2.01e5;
  std::string xtitle = isInv? "M_{T}(Z') [GeV]" : "M(Z') [GeV]";
  PlotDistribution(ch, year, module, hname, 700, 5000, 1.5e-02, ymax, xtitle);
}


void PlotCount(std::string ch="muonchannel", std::string year="RunII", std::string module="Selection") {
  std::string hname="ZprimeCandidate_ScaleFactors/sum_event_weights";
  if (module=="SignalRegion") hname = TString(hname).ReplaceAll("ScaleFactors","Selection").Data();
  bool isInv = ch=="invisiblechannel";
  double ymax = isInv? 7.01e12 :7.01e9;
  PlotDistribution(ch, year, module, hname, 0.45, 1.55, 1.5, ymax, "Counting Experiment");
}

void PlotTagger(std::string ch="muonchannel", std::string year="RunII", std::string module="Selection", std::string var="ZHccvsQCD", bool isMD=true) {
  std::string hname= "ZprimeCandidate_ScaleFactors/H_btag_DeepBoosted_"+var;
  if (isMD) hname = TString(hname).ReplaceAll("DeepBoosted", "MassDecorrelatedDeepBoosted").Data();
  if (module=="SignalRegion") hname = TString(hname).ReplaceAll("ScaleFactors","Selection").Data();
  bool isInv = ch=="invisiblechannel";
  double ymax = isInv? 7.01e12 :7.01e9;
  PlotDistribution(ch, year, module, hname, 0, 1, 1.5e-02, ymax, var);
}


void PlotPt(std::string ch="muonchannel", std::string year="RunII", std::string module="Selection", std::string hname_="Z_pt", std::string xname="p_{T} [GeV]") {
  std::string hname="ZprimeCandidate_ScaleFactors/"+hname_;
  if (module=="SignalRegion") hname = TString(hname).ReplaceAll("ScaleFactors","Selection").Data();
  bool isInv = ch=="invisiblechannel";
  double ymax = isInv? 4.01e7 :8.01e5;
  PlotDistribution(ch, year, module, hname, 200, 2500, 1.5e-02, ymax, xname);
}

void PlotEtaPhi(std::string ch="muonchannel", std::string year="RunII", std::string module="Selection", std::string hname_="Z_eta", std::string xname="#eta") {
  std::string hname="ZprimeCandidate_ScaleFactors/"+hname_;
  if (module=="SignalRegion") hname = TString(hname).ReplaceAll("ScaleFactors","Selection").Data();
  bool isInv = ch=="invisiblechannel";
  double ymax = isInv? 7.01e8 :2.01e8;
  PlotDistribution(ch, year, module, hname, -3, 3, 1.5e-02, ymax, xname);
}

void PlotHDistributions(std::string  ch="muonchannel", std::string year="RunII", std::string module="Selection") {
  std::string extra = "jet";
  PlotPt(ch,  year, module, "H_pt",  "p_{T}^{"+extra+"} [GeV]");
  PlotEtaPhi(ch, year, module, "H_phi", "#phi^{"+extra+"}" );
  PlotEtaPhi(ch, year, module, "H_eta", "#eta^{"+extra+"}" );

  if (module=="Selection") PlotTagger(ch, year, module, "ZHccvsQCD", true);
  PlotTagger(ch, year, module, "H4qvsQCD", false);
}

void PlotZDistributions(std::string  ch="muonchannel", std::string year="RunII", std::string module="Selection") {
  bool isInv = ch=="invisiblechannel";
  bool isEle = ch=="electronchannel";
  std::string extra = isInv? "miss": (isEle? "ee": "#mu#mu");
  PlotPt(ch,     year, module, "Z_pt" , "p_{T}^{"+extra+"} [GeV]");
  PlotEtaPhi(ch, year, module, "Z_phi", "#phi^{"+extra+"}" );
  PlotEtaPhi(ch, year, module, "Z_eta", "#eta^{"+extra+"}" );
}





int main(){
  gErrorIgnoreLevel = kError;
  for (const auto& ch: {"muonchannel","electronchannel","invisiblechannel"}) {
    for (const auto& year: {"2016","2017","2018","RunII"}) {
      for (const auto& module: {"Selection","SignalRegion"}) {
        PrintGreen<std::string>({ch,year,module});
        PlotZprimeMass(ch,year,module);
        PlotCount(ch,year,module);
        PlotHDistributions(ch,year,module);
        PlotZDistributions(ch,year,module);
      }
    }
  }
  return 0;
}
