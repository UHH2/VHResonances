#include "CalculateSignalEfficiencies.hpp"
#include <math.h>


/*
Given X as genering channel (e.g. ZeeHbb), then one can define:
- TOT = Sum_X X
- BR_X = X/TOT
Given a cut C, such that X->X_c, then one can define:
- tot efficiency as e_TOT = Sum_X X_c/TOT
- channel efficiency as e_x = X_c/TOT
- relative channel efficiency as e_x^rel = X_c/X
- normalized channel efficiency as e_x^norm = e_x^rel*BR_X = e_x
Then the following behaviour ar to be expected.
- e_TOT = sum_X e_x
- e_TOT = sum_X e_x^rel*BR_X = sum_X e_x^norm

For sake of simplicity we calculate only channel efficiency, since the relative
channel efficiency can be derived only with a BR.

We do not run over "Zelse" (aka tauchannel) because we are almost 100% not sensitive to it.
We run over Helse, since ~1/3 of the final signal is from there.
Zlep involves only muon and electron. In very good approximation can be considered as Zll,
since Ztau is not present (valid after the request of no leptons of the other type).

*/

double EfficiencyError(double count, double N)  { double eff = count/N; return TMath::Sqrt(eff*(1-eff)/N); }

bool isCSRegion(std::string tag) {return (tag.find("SR")!=std::string::npos || tag.find("CR")!=std::string::npos); }

//assume Z->ll only
const std::unordered_map<std::string,double> BRs ={
  {"Zmumu",       0.33}, {"Zee",         0.33}, {"Zelse",       0.33},
  {"HtoWW",       0.10}, {"Htobb",       0.58}, {"Helse",       0.32},
  {"ZmumuHtoWW",  0.03}, {"ZmumuHtobb",  0.19}, {"ZmumuHelse",  0.11},
  {"ZeeHtoWW",    0.03}, {"ZeeHtobb",    0.19}, {"ZeeHelse",    0.11},
  {"ZelseHtoWW",  0.03}, {"ZelseHtobb",  0.19}, {"ZelseHelse",  0.11},
};

void CalculateSignalEfficiencies(std::string histFolder) {

  std::string user            = std::getenv("USER");
  std::string Path_NFS        = "/nfs/dust/cms/user/"+user+"/";
  std::string Path_STORAGE    = Path_NFS+"WorkingArea/File/Analysis/";
  std::string Path_ANALYSIS   = std::getenv("CMSSW_BASE"); Path_ANALYSIS += "/src/UHH2/VHResonances/Analysis/";

  std::string outdir = Path_ANALYSIS+"/SignalEfficiencies/";
  std::string prefix = "uhh2.AnalysisModuleRunner.MC.";
  std::string syst = "nominal";

  std::vector<std::string> years =  {"2016", "2017", "2018", "RunII"};
  // std::vector<std::string> years =  {"RunII"};
  // std::vector<std::string> collections =  {"Puppi", "CHS", "HOTVR"};
  std::vector<std::string> collections =  {"Puppi"};
  std::vector<std::string> channels =  {"muonchannel", "electronchannel", "leptonchannel", "invisiblechannel" };
  // std::vector<std::string> channels =  {"muonchannel"};
  std::vector<std::string> decaymodes =  {"bb", "WW", "else", "Inc" };

  double x_min = 800;
  double x_max = 5200;
  double y_min = 0.01;
  double y_max = 10;
  TString x_name = "M_{Z'} (GeV)";
  TString y_name = "Selection efficiency";

  writeExtraText = true;       // if extra text
  extraText  = "Simulation";
  extraText2 = "Work in progress";
  lumi_13TeV  = "";
  // ForThesis();

  std::map<std::string, mypair_I > Cuts;

  // Cuts.insert(std::pair<std::string, mypair_I>("0_nocuts",                              mypair_I("weights",kBlack))); // only as control
  Cuts.insert(std::pair<std::string, mypair_I>("1_Triggers",                            mypair_I("Trigger",                       kRed+1)));
  Cuts.insert(std::pair<std::string, mypair_I>("2_Lepton selection",                    mypair_I("NLeptonSel",                    kOrange+1)));
  // Cuts.insert(std::pair<std::string, mypair_I>("3_#Delta R(ll) < 1.0",                  mypair_I("DeltaRDiLepton",                kOrange+1)));
  // Cuts.insert(std::pair<std::string, mypair_I>("4_#Delta#phi(ll,#it{jet}) #geq #pi/2",  mypair_I("JetDiLeptonPhiAngular",         kGreen+1)));
  // Cuts.insert(std::pair<std::string, mypair_I>("5_preselection",                        mypair_I("Preselection",                  kBlack))); // only as control
  // Cuts.insert(std::pair<std::string, mypair_I>("5_preselection",                        mypair_I("Preselection",                  kBlack))); // only as control
  Cuts.insert(std::pair<std::string, mypair_I>("5_Angular cuts",                           mypair_I("MuonScale",                  kOrange-2)));
  // Cuts.insert(std::pair<std::string, mypair_I>("7_Z' Reconstruction",                   mypair_I("ZprimeReco",                    kGreen+2))); // only as control
  Cuts.insert(std::pair<std::string, mypair_I>("7_b tag veto",                          mypair_I("ZprimeSelection",               kBlue+1)));
  // Cuts.insert(std::pair<std::string, mypair_I>("8_#it{p}_{T}^{#it{ll}}/m_{Z'} #geq 0.2",mypair_I("PTMassCut",                     kBlue+1))); //Very similar to ZprimeSelection
  // Cuts.insert(std::pair<std::string, mypair_I>("11_FullSelection",                       mypair_I("ExtraCleaning",                 kOrange-1))); // TODO The SFs here play a role

  Cuts.insert(std::pair<std::string, mypair_I>("99_H4qvsQCD",                           mypair_I("DeepAk8_H4qvsQCD_massdep_SR",              kGreen+3)));
  // Cuts.insert(std::pair<std::string, mypair_I>("99_HccvsQCD",                           mypair_I("DeepAk8_HccvsQCD2_SR",                  kGreen+3)));
  Cuts.insert(std::pair<std::string, mypair_I>("99_ZHccvsQCD",                         mypair_I("DeepAk8_ZHccvsQCD_MD2_SR",              kAzure+1)));
  // Cuts.insert(std::pair<std::string, mypair_I>("99_Hcc_MD",                             mypair_I("DeepAk8_HccvsQCD_MD_SR",                   kRed+1)));
  // Cuts.insert(std::pair<std::string, mypair_I>("99_H4q_MD",                             mypair_I("DeepAk8_H4qvsQCD_MD_SR",                   kAzure-7)));
  // Cuts.insert(std::pair<std::string, mypair_I>("99_H4q_PD_Hcc_MD",                      mypair_I("DeepAk8_H4qvsQCD_massdep_HccvsQCD_MD_SR",  kOrange+1)));
  // Cuts.insert(std::pair<std::string, mypair_I>("99_H4q",                                mypair_I("DeepAk8_H4qvsQCD_SR",                      kCyan+3)));
  // Cuts.insert(std::pair<std::string, mypair_I>("99_Hcc",                                mypair_I("DeepAk8_HccvsQCD_SR",                      kViolet+1)));
  // Cuts.insert(std::pair<std::string, mypair_I>("99_ZHcc",                               mypair_I("DeepAk8_ZHccvsQCD_SR",                     kCyan+1)));
  // Cuts.insert(std::pair<std::string, mypair_I>("99_H4q_PD_Hcc",                         mypair_I("DeepAk8_H4qvsQCD_massdep_HccvsQCD_SR",     kOrange-1)));
  // Cuts.insert(std::pair<std::string, mypair_I>("99_H4q_PD_ZHcc",                        mypair_I("DeepAk8_H4qvsQCD_massdep_ZHccvsQCD_SR",    kMagenta+2)));
  // Cuts.insert(std::pair<std::string, mypair_I>("99_H4q_PD_ZHcc_MD",                     mypair_I("DeepAk8_H4qvsQCD_massdep_ZHccvsQCD_MD_SR", kGreen+1)));


  Cuts.insert(std::pair<std::string, mypair_I>("20_tau42_SR",                           mypair_I("tau42_SR",                      kOrange+1)));
  Cuts.insert(std::pair<std::string, mypair_I>("21_tau42_CR",                           mypair_I("tau42_CR",                      kOrange-2)));


  std::map<std::string, TGraph* > Plot_ComparisonFinal;
  std::map<std::string, std::vector<double> > SignalEfficiencies;
  std::vector<std::string> order_norm;
  std::vector<double> dummy(MyMassPoints.size(),0.001);

  for (std::string year: years) {
    for (std::string collection: collections) {
      for (std::string channel: channels) {
        for (std::string decaymode: decaymodes) {
          std::string namePlot = decaymode+"_"+collection+"_"+channel+"_"+year+"_"+histFolder;
          std::cout << namePlot << std::endl;
          TString namePlotShort = "Hto"+namePlot; namePlotShort.ReplaceAll("Puppi_","").ReplaceAll("muon","#mu-").ReplaceAll("electron","e-").ReplaceAll("lepton","l-").ReplaceAll("_"," ");

          bool isInc = decaymode=="Inc";
          bool isLeptonChannel = channel=="leptonchannel";
          bool isInvisibleChannel = channel=="invisiblechannel";

          TString hname = "sum_event_weights";

          if (!isInvisibleChannel){
            if (!isInc || !isLeptonChannel) hname += "_";
            if (!isLeptonChannel) hname = hname+"Z"+ ((channel=="muonchannel") ? "mumu" : "ee" );
            if (!isInc) hname += "Hto"+decaymode;
            if (decaymode=="else") hname.ReplaceAll("Hto","H");
          }else{
            if (!isInc) hname += "_Hto"+decaymode;
            if (decaymode=="else") hname.ReplaceAll("Hto","H");
          }

          std::string PresectionStorePath = Path_STORAGE+year+"/Preselection/"+collection+"/"+channel+"/"+syst+"/";
          std::string SectionStorePath    = Path_STORAGE+year+"/Selection/"+collection+"/"+channel+"/"+syst+"/";
          std::string CSRStorePath        = Path_STORAGE+year+"/SignalRegion/"+collection+"/"+channel+"/"+syst+"/";
          // Reset values
          SignalEfficiencies.clear();
          order_norm.clear();

          for (std::pair<std::string, mypair_I> element : Cuts) {
            if (!isCSRegion(element.first)) order_norm.push_back(element.first);
            SignalEfficiencies[element.first] = std::vector<double>(MyMassPoints.size(), 0);
          }

          // To Load the leptonhisto we loop over muon and electron channels
          // Maybe nont too smart, but it works
          int loop = (isLeptonChannel)?2:1;

          std::string additionalText = "";
          if (channel=="invisiblechannel"){additionalText = "_inv";}

          for (int lep = 0; lep < loop ; lep++) {

            int Mass_index = 0;
            for (int MassPoint : MyMassPoints) {
              std::string MassName  = std::to_string((int)MassPoint);
              TString fn_presel = PresectionStorePath+prefix+"MC_ZprimeToZH"+additionalText+"_M"+MassName+"_"+year+"_noTree.root";
              TString fn_sel    = SectionStorePath+prefix+"MC_ZprimeToZH"+additionalText+"_M"+MassName+"_"+year+"_noTree.root";
              TString fn_csr    = CSRStorePath+prefix+"MC_ZprimeToZH"+additionalText+"_M"+MassName+"_"+year+"_noTree.root";

              if (isLeptonChannel) {
                fn_presel.ReplaceAll("leptonchannel","electronchannel");
                fn_sel.ReplaceAll("leptonchannel","electronchannel");
                fn_csr.ReplaceAll("leptonchannel","electronchannel");
                if (lep==1) {
                  fn_presel.ReplaceAll("Zee","Zmumu");
                  fn_sel.ReplaceAll("Zee","Zmumu");
                  fn_csr.ReplaceAll("Zee","Zmumu");
                  fn_presel.ReplaceAll("electronchannel","muonchannel");
                  fn_sel.ReplaceAll("electronchannel","muonchannel");
                  fn_csr.ReplaceAll("electronchannel","muonchannel");
                }
              }

              TFile *file_presel  = new TFile(fn_presel);
              TFile *file_sel     = new TFile(fn_sel);
              TFile *file_csr     = new TFile(fn_csr);

              // Load histograms
              for (std::pair<std::string, mypair_I> element : Cuts) {
                std::string tag = element.first;
                std::string cut = element.second.first;

                bool isCSR = isCSRegion(tag); // Maybe not too interesting in all the plots produced. Think about it first.
                TH1F* h_tot = (TH1F*)file_presel->Get("ZprimeCandidate_weights/sum_event_weights");
                if (isCSR) h_tot = (TH1F*)file_csr->Get("ZprimeCandidate_Selection/sum_event_weights");
                double tot_event = h_tot->GetBinContent(1);

                // Try to load hist from both Preselection and selection file.
                TH1F* h_preselection = (TH1F*)file_presel->Get("ZprimeCandidate_"+cut+"/"+hname);
                TH1F* h_selection = (TH1F*)file_sel->Get("ZprimeCandidate_"+cut+"/"+hname);
                if (isCSR) h_selection = (TH1F*)file_csr->Get("ZprimeCandidate_"+cut+"/"+hname);

                //Decide which is the right one.
                TH1F* h_sel;
                if (h_preselection) h_sel = h_preselection;
                else if (h_selection) h_sel = h_selection;
                else h_sel = (TH1F*)file_csr->Get("ZprimeCandidate_"+cut+"/"+hname); //In case the hist in the SignalRegion file (eg. 9_Hcc)

                double sel_event = h_sel->GetBinContent(1);

                if (!(h_preselection) && ! (h_selection)) {//Remove effect of Scale Factors
                  TH1F* h_SFs = (TH1F*)file_sel->Get("ZprimeCandidate_ScaleFactors/"+hname);
                  TH1F* h_noSFs = (TH1F*)file_sel->Get("ZprimeCandidate_PTMassCut/"+hname);
                  sel_event *= h_noSFs->GetBinContent(1)/h_SFs->GetBinContent(1);
                }
                SignalEfficiencies[tag][Mass_index] += sel_event/tot_event;
              }
              file_presel->Close(); file_sel->Close(); file_csr->Close();
              Mass_index++;
            }
          }

          // Plot efficiencies
          lumi_13TeV  = TString::Format("%.1f fb^{-1}", lumi_map.at(year).at("lumi_fb"));

          y_max = 10;
          TCanvas* canv_eff = tdrCanvas(("canv_eff"+namePlot).c_str(), x_min, x_max, y_min, y_max, x_name, y_name);
          canv_eff->SetLogy(1);
          TCanvas* canv_eff_SR = tdrCanvas(("canv_eff_SR"+namePlot).c_str(), x_min, x_max, y_min, y_max, x_name, y_name);
          canv_eff_SR->SetLogy(1);
          // TLegend *leg_eff = tdrLeg(0.50,0.68,0.89,0.89, 0.030, 42, kBlack);
          TLegend *leg_eff = tdrLeg(0.40,0.68,0.89,0.89, 0.030, 42, kBlack);
          leg_eff->SetNColumns(3);
          // TLegend *leg_eff_SR = tdrLeg(0.40,0.70,0.89,0.89, 0.030, 42, kBlack);
          TLegend *leg_eff_SR = tdrLeg(0.40,0.68,0.89,0.89, 0.030, 42, kBlack);
          leg_eff_SR->SetNColumns(3);

          for (std::pair<std::string, mypair_I> element : Cuts) {
            std::string tag = element.first;
            if (tag.find("CR")!=std::string::npos) continue; // TODO plot also CR
            int color = element.second.second;
            bool isCSR = isCSRegion(tag);
            TGraph* gr_eff = new TGraphErrors(MyMassPoints.size(), &(MyMassPoints[0]), &(SignalEfficiencies[tag][0]));
            gr_eff->SetLineWidth(2);
            if (isCSR) canv_eff_SR->cd();
            else canv_eff->cd();
            tdrDraw(gr_eff, "lp", kFullDotLarge, color, kSolid, color, 1000, color);
            if (isCSR) leg_eff_SR->AddEntry(gr_eff, tag.substr(tag.find("_")+1).c_str(),"lp");
            else leg_eff->AddEntry(gr_eff, tag.substr(tag.find("_")+1).c_str(),"lp");

            //Save graph for further usage. We care only about the SR eff.
            if(tag=="99_ZHccvsQCD") Plot_ComparisonFinal[namePlot] = gr_eff;
          }

          canv_eff->cd();
          leg_eff->Draw("same");
          canv_eff->SaveAs((outdir+"Eff_HTo"+namePlot+".pdf").c_str());
          leg_eff->Delete();

          canv_eff_SR->cd();
          leg_eff_SR->Draw("same");
          canv_eff_SR->SaveAs((outdir+"Eff_SR_HTo"+namePlot+".pdf").c_str());
          leg_eff_SR->Delete();

          // Create eff normalized wrt the previous step
          // NB. This is nothing related to e_x^norm
          y_max = 1.5;
          TCanvas* canv_eff_norm = tdrCanvas(("canv_eff_norm"+namePlot).c_str(), x_min, x_max, y_min, y_max, x_name, y_name);
          // TLegend *leg_eff_norm = tdrLeg(0.50,0.68,0.89,0.89, 0.030, 42, kBlack);
          TLegend *leg_eff_norm = tdrLeg(0.40,0.68,0.89,0.89, 0.030, 42, kBlack);
          leg_eff_norm->SetNColumns(3);

          for (unsigned int i = 0; i < order_norm.size(); i++) {
            std::string tag = order_norm[i];
            std::string tag_norm = order_norm[std::max((int)(i-1),0)];
            int color = Cuts[tag].second;
            std::vector<double> eff_norm(MyMassPoints.size(), 0);
            for (unsigned int m = 0; m < MyMassPoints.size(); m++) eff_norm[m] = SignalEfficiencies[tag][m]/SignalEfficiencies[tag_norm][m];
            TGraph* gr_eff = new TGraphErrors(MyMassPoints.size(), &(MyMassPoints[0]), &(eff_norm[0]));
            gr_eff->SetLineWidth(2);
            tdrDraw(gr_eff, "lp", kFullDotLarge, color, kSolid, color, 1000, color);
            leg_eff_norm->AddEntry(gr_eff, tag.substr(tag.find("_")+1).c_str(),"lp");
          }

          leg_eff_norm->Draw("same");
          canv_eff_norm->SaveAs((outdir+"Eff_norm_HTo"+namePlot+".pdf").c_str());
          leg_eff_norm->Delete();

          // Create eff relative e_x^rel
          y_max = 1.5;
          TCanvas* canv_eff_rel = tdrCanvas(("canv_rel_eff"+namePlot).c_str(), x_min, x_max, y_min, y_max, x_name, y_name);
          // canv_eff_rel->SetLogy(1);
          // TLegend *leg_eff_rel = tdrLeg(0.50,0.68,0.89,0.89, 0.030, 42, kBlack);
          TLegend *leg_eff_rel = tdrLeg(0.40,0.68,0.89,0.89, 0.030, 42, kBlack);
          leg_eff_rel->SetNColumns(2);

          for (std::pair<std::string, mypair_I> element : Cuts) {
            std::string tag = element.first;
            if (isCSRegion(tag)) continue; //Interested only in the selection part.
            if (tag.find("CR")!=std::string::npos) continue; // TODO plot also CR
            int color = element.second.second;
            std::vector<double> eff_rel(MyMassPoints.size(), 0);
            TString BRname = hname;
            BRname.ReplaceAll("sum_event_weights_","");
            double BR = (BRname!="sum_event_weights")? BRs.at(BRname.Data()) : 1;
            for (unsigned int m = 0; m < MyMassPoints.size(); m++) {
              eff_rel[m] = SignalEfficiencies[tag][m]/BR;
              // This happens only because of the finite approximation of the BR. Force it to 1.
              if (eff_rel[m] > 1) eff_rel[m] = round(eff_rel[m]);
            }
            TGraph* gr_eff = new TGraphErrors(MyMassPoints.size(), &(MyMassPoints[0]), &(eff_rel[0]));
            gr_eff->SetLineWidth(2);
            canv_eff_rel->cd();
            tdrDraw(gr_eff, "lp", kFullDotLarge, color, kSolid, color, 1000, color);
            leg_eff_rel->AddEntry(gr_eff, tag.substr(tag.find("_")+1).c_str(),"lp");
          }

          canv_eff_rel->cd();
          leg_eff_rel->Draw("same");
          canv_eff_rel->SaveAs((outdir+"Eff_rel_HTo"+namePlot+".pdf").c_str());
          leg_eff_rel->Delete();

        }
      }
    }
  }
  // Now all the eff are stored. One can plot extra comparisons

  TCanvas* canv_ComparisonFinal;
  y_min = 5*1e-04;
  y_max = 3.;
  TLegend* leg_ComparisonFinal;
  int color, lineStyle, markerSyle=kFullDotLarge;
  std::string namePlot, nameLeg;
  if (FindInVector(years, "RunII")>0) {

    lumi_13TeV  = TString::Format("%.1f fb^{-1}", lumi_map.at("RunII").at("lumi_fb"));
    canv_ComparisonFinal = tdrCanvas("canv_ComparisonFinal_RunII", x_min, x_max, y_min, y_max, x_name, y_name);
    canv_ComparisonFinal->SetLogy(1);
    leg_ComparisonFinal = tdrLeg(0.40,0.68,0.89,0.89, 0.030, 42, kBlack);
    leg_ComparisonFinal->SetNColumns(3);

    // Plot invisible channel
    if (std::count(channels.begin(), channels.end(), "invisiblechannel")>0){

      color = kBlack+3;  markerSyle=kFullTriangleDown; lineStyle = kSolid;  namePlot = "WW_Puppi_invisiblechannel_RunII_"+histFolder;     nameLeg = "Z#nu#nuHWW";       tdrDraw(Plot_ComparisonFinal[namePlot], "lp", markerSyle, color, lineStyle, color, 1000, color); leg_ComparisonFinal->AddEntry(Plot_ComparisonFinal[namePlot], nameLeg.c_str(),"lp");
      color = kBlack+1;  markerSyle=kFullSquare;       lineStyle = kDashed; namePlot = "bb_Puppi_invisiblechannel_RunII_"+histFolder;     nameLeg = "Z#nu#nuHbb";       tdrDraw(Plot_ComparisonFinal[namePlot], "lp", markerSyle, color, lineStyle, color, 1000, color); leg_ComparisonFinal->AddEntry(Plot_ComparisonFinal[namePlot], nameLeg.c_str(),"lp");
      color = kBlack+2;  markerSyle=kFullTriangleUp;   lineStyle = kDotted; namePlot = "else_Puppi_invisiblechannel_RunII_"+histFolder;   nameLeg = "Z#nu#nuHelse";     tdrDraw(Plot_ComparisonFinal[namePlot], "lp", markerSyle, color, lineStyle, color, 1000, color); leg_ComparisonFinal->AddEntry(Plot_ComparisonFinal[namePlot], nameLeg.c_str(),"lp");
    }

    if (channels.size()!=1 || channels[0] != "invisiblechannel"){ // Plot the other channels as well, if given

      color = kRed+1;    markerSyle=kFullDotLarge;     lineStyle = kSolid;  namePlot = "WW_Puppi_muonchannel_RunII_"+histFolder;       nameLeg = "Z#mu#muHWW";   tdrDraw(Plot_ComparisonFinal[namePlot], "lp", markerSyle, color, lineStyle, color, 1000, color); leg_ComparisonFinal->AddEntry(Plot_ComparisonFinal[namePlot], nameLeg.c_str(),"lp");
      color = kOrange+1; markerSyle=kFullCircle;       lineStyle = kDashed; namePlot = "bb_Puppi_muonchannel_RunII_"+histFolder;       nameLeg = "Z#mu#muHbb";   tdrDraw(Plot_ComparisonFinal[namePlot], "lp", markerSyle, color, lineStyle, color, 1000, color); leg_ComparisonFinal->AddEntry(Plot_ComparisonFinal[namePlot], nameLeg.c_str(),"lp");
      color = kOrange-2; markerSyle=kOpenCircle;       lineStyle = kDotted; namePlot = "else_Puppi_muonchannel_RunII_"+histFolder;     nameLeg = "Z#mu#muHelse"; tdrDraw(Plot_ComparisonFinal[namePlot], "lp", markerSyle, color, lineStyle, color, 1000, color); leg_ComparisonFinal->AddEntry(Plot_ComparisonFinal[namePlot], nameLeg.c_str(),"lp");
      color = kBlue+1;   markerSyle=kFullDiamond;      lineStyle = kSolid;  namePlot = "WW_Puppi_electronchannel_RunII_"+histFolder;   nameLeg = "ZeeHWW";       tdrDraw(Plot_ComparisonFinal[namePlot], "lp", markerSyle, color, lineStyle, color, 1000, color); leg_ComparisonFinal->AddEntry(Plot_ComparisonFinal[namePlot], nameLeg.c_str(),"lp");
      color = kViolet+1; markerSyle=kFullCross;        lineStyle = kDashed; namePlot = "bb_Puppi_electronchannel_RunII_"+histFolder;   nameLeg = "ZeeHbb";       tdrDraw(Plot_ComparisonFinal[namePlot], "lp", markerSyle, color, lineStyle, color, 1000, color); leg_ComparisonFinal->AddEntry(Plot_ComparisonFinal[namePlot], nameLeg.c_str(),"lp");
      color = kViolet;   markerSyle=kOpenSquare;       lineStyle = kDotted; namePlot = "else_Puppi_electronchannel_RunII_"+histFolder; nameLeg = "ZeeHelse";     tdrDraw(Plot_ComparisonFinal[namePlot], "lp", markerSyle, color, lineStyle, color, 1000, color); leg_ComparisonFinal->AddEntry(Plot_ComparisonFinal[namePlot], nameLeg.c_str(),"lp");
      color = kGreen+3;  markerSyle=kFullTriangleDown; lineStyle = kSolid;  namePlot = "WW_Puppi_leptonchannel_RunII_"+histFolder;     nameLeg = "ZllHWW";       tdrDraw(Plot_ComparisonFinal[namePlot], "lp", markerSyle, color, lineStyle, color, 1000, color); leg_ComparisonFinal->AddEntry(Plot_ComparisonFinal[namePlot], nameLeg.c_str(),"lp");
      color = kGreen+1;  markerSyle=kFullSquare;       lineStyle = kDashed; namePlot = "bb_Puppi_leptonchannel_RunII_"+histFolder;     nameLeg = "ZllHbb";       tdrDraw(Plot_ComparisonFinal[namePlot], "lp", markerSyle, color, lineStyle, color, 1000, color); leg_ComparisonFinal->AddEntry(Plot_ComparisonFinal[namePlot], nameLeg.c_str(),"lp");
      color = kGreen+2;  markerSyle=kFullTriangleUp;   lineStyle = kDotted; namePlot = "else_Puppi_leptonchannel_RunII_"+histFolder;   nameLeg = "ZllHelse";     tdrDraw(Plot_ComparisonFinal[namePlot], "lp", markerSyle, color, lineStyle, color, 1000, color); leg_ComparisonFinal->AddEntry(Plot_ComparisonFinal[namePlot], nameLeg.c_str(),"lp");
    }
    leg_ComparisonFinal->Draw("same");
    canv_ComparisonFinal->SaveAs((outdir+"Eff_ComparisonFinal_RunII_"+histFolder+".pdf").c_str());
  }
  lumi_13TeV  = TString::Format("%.1f fb^{-1}", lumi_map.at("RunII").at("lumi_fb"));
  canv_ComparisonFinal = tdrCanvas("canv_ComparisonFinal_Inc", x_min, x_max, y_min, y_max, x_name, y_name);
  canv_ComparisonFinal->SetLogy(1);
  leg_ComparisonFinal = tdrLeg(0.30,0.68,0.89,0.89, 0.025, 42, kBlack);
  leg_ComparisonFinal->SetNColumns(3);

  for (std::string year: years) {

    if (year=="2016")  {color = kGreen+2;  markerSyle=kFullTriangleDown;}
    if (year=="2017")  {color = kAzure+1;  markerSyle=kFullTriangleUp;}
    if (year=="2018")  {color = kRed+1;    markerSyle=kFullSquare;}
    if (year=="RunII") {color = kOrange+1; markerSyle=kFullDotLarge;}
    for (std::string channel: channels) {
      namePlot = "Inc_Puppi_"+channel+"_"+year+"_"+histFolder;
      nameLeg = year+std::string(6-year.size(), ' ' );
      if (channel=="muonchannel") {     lineStyle = kDashed; nameLeg += "Z#mu#muHInc";}
      if (channel=="electronchannel") { lineStyle = kDotted; nameLeg += "ZeeHInc";}
      if (channel=="leptonchannel") {   lineStyle = kSolid; nameLeg += "ZllHInc";}
      if (channel=="invisiblechannel"){ lineStyle = kSolid; nameLeg += "Z#nu#nuHInc";}
      tdrDraw(Plot_ComparisonFinal[namePlot], "lp", markerSyle, color, lineStyle, color, 1000, color); leg_ComparisonFinal->AddEntry(Plot_ComparisonFinal[namePlot], nameLeg.c_str(),"lp");
    }
  }
  leg_ComparisonFinal->Draw("same");
  canv_ComparisonFinal->SaveAs((outdir+"Eff_ComparisonFinal_Inc_"+histFolder+".pdf").c_str());



  if (FindInVector(years, "RunII")>0) {

    lumi_13TeV  = TString::Format("%.1f fb^{-1}", lumi_map.at("RunII").at("lumi_fb"));
    canv_ComparisonFinal = tdrCanvas("canv_ComparisonFinal_RunII_Higgs_", x_min, x_max, y_min, y_max, x_name, y_name);
    canv_ComparisonFinal->SetLogy(1);
    leg_ComparisonFinal = tdrLeg(0.40,0.68,0.89,0.89, 0.030, 42, kBlack);
    leg_ComparisonFinal->SetNColumns(3);

    // Plot invisible channel
    if (std::count(channels.begin(), channels.end(), "invisiblechannel")>0){
      color = kBlack;    markerSyle=kFullTriangleUp;   lineStyle = kSolid;  namePlot = "else_Puppi_invisiblechannel_RunII_"+histFolder;  nameLeg = "Z#nu#nuHelse";    tdrDraw(Plot_ComparisonFinal[namePlot], "lp", markerSyle, color, lineStyle, color, 1000, color); leg_ComparisonFinal->AddEntry(Plot_ComparisonFinal[namePlot], nameLeg.c_str(),"lp");
      color = kBlack-1;  markerSyle=kFullDotLarge;     lineStyle = kSolid;  namePlot = "Inc_Puppi_invisiblechannel_RunII_"+histFolder;   nameLeg = "Z#nu#nuHInc";     tdrDraw(Plot_ComparisonFinal[namePlot], "lp", markerSyle, color, lineStyle, color, 1000, color); leg_ComparisonFinal->AddEntry(Plot_ComparisonFinal[namePlot], nameLeg.c_str(),"lp");
    }
    if (channels.size()!=1 || channels[0] != "invisiblechannel"){ // Plot the other channels as well, if given

      color = kGreen+3;  markerSyle=kFullTriangleDown; lineStyle = kSolid;  namePlot = "WW_Puppi_leptonchannel_RunII_"+histFolder;    nameLeg = "ZllHWW";      tdrDraw(Plot_ComparisonFinal[namePlot], "lp", markerSyle, color, lineStyle, color, 1000, color); leg_ComparisonFinal->AddEntry(Plot_ComparisonFinal[namePlot], nameLeg.c_str(),"lp");
      color = kGreen+1;  markerSyle=kFullSquare;       lineStyle = kSolid;  namePlot = "bb_Puppi_leptonchannel_RunII_"+histFolder;    nameLeg = "ZllHbb";      tdrDraw(Plot_ComparisonFinal[namePlot], "lp", markerSyle, color, lineStyle, color, 1000, color); leg_ComparisonFinal->AddEntry(Plot_ComparisonFinal[namePlot], nameLeg.c_str(),"lp");
      color = kGreen+2;  markerSyle=kFullTriangleUp;   lineStyle = kSolid;  namePlot = "else_Puppi_leptonchannel_RunII_"+histFolder;  nameLeg = "ZllHelse";    tdrDraw(Plot_ComparisonFinal[namePlot], "lp", markerSyle, color, lineStyle, color, 1000, color); leg_ComparisonFinal->AddEntry(Plot_ComparisonFinal[namePlot], nameLeg.c_str(),"lp");
      color = kViolet+1; markerSyle=kOpenSquare;       lineStyle = kDotted; namePlot = "Inc_Puppi_muonchannel_RunII_"+histFolder;     nameLeg = "Z#mu#muHInc"; tdrDraw(Plot_ComparisonFinal[namePlot], "lp", markerSyle, color, lineStyle, color, 1000, color); leg_ComparisonFinal->AddEntry(Plot_ComparisonFinal[namePlot], nameLeg.c_str(),"lp");
      color = kAzure+1;  markerSyle=kFullCross;        lineStyle = kDotted; namePlot = "Inc_Puppi_electronchannel_RunII_"+histFolder; nameLeg = "ZeeHInc";     tdrDraw(Plot_ComparisonFinal[namePlot], "lp", markerSyle, color, lineStyle, color, 1000, color); leg_ComparisonFinal->AddEntry(Plot_ComparisonFinal[namePlot], nameLeg.c_str(),"lp");
      color = kRed+1;    markerSyle=kFullDotLarge;     lineStyle = kSolid;  namePlot = "Inc_Puppi_leptonchannel_RunII_"+histFolder;   nameLeg = "ZllHInc";     tdrDraw(Plot_ComparisonFinal[namePlot], "lp", markerSyle, color, lineStyle, color, 1000, color); leg_ComparisonFinal->AddEntry(Plot_ComparisonFinal[namePlot], nameLeg.c_str(),"lp");
    }
    leg_ComparisonFinal->Draw("same");
    canv_ComparisonFinal->SaveAs((outdir+"Eff_ComparisonFinal_RunII_Higgs_"+histFolder+".pdf").c_str());
  }
}

//TODO import this func from Utils
int FindInVector(const std::vector<std::string>& vec, const std::string& el) {
  int index = -1;
  // Find given element in vector
  auto it = std::find(vec.begin(), vec.end(), el);
  if (it != vec.end()) index = distance(vec.begin(), it);
  return index;
}


int main(int argc, char** argv){
  gSystem->Exec("mkdir -p ./SignalEfficiencies");

  std::vector<std::string> histFolders = { "DeepAk8_H4qvsQCD_massdep", "DeepAk8_ZHccvsQCD_MD",
  "DeepAk8_HccvsQCD_MD", "DeepAk8_H4qvsQCD_MD", "DeepAk8_H4qvsQCD_massdep_HccvsQCD_MD",
  "DeepAk8_H4qvsQCD", "DeepAk8_HccvsQCD", "DeepAk8_ZHccvsQCD", "DeepAk8_H4qvsQCD_massdep_HccvsQCD",
  "DeepAk8_H4qvsQCD_massdep_ZHccvsQCD", "DeepAk8_H4qvsQCD_massdep_ZHccvsQCD_MD",
  "DeepAk8_HccvsQCD2", "DeepAk8_ZHccvsQCD_MD2", "tau42"};

  if (argc>1) {
    std::string histFolder;
    for (int i = 1; i < argc; i++) {
      if (std::find(histFolders.begin(), histFolders.end(), argv[i]) != histFolders.end() ) histFolder = argv[i];
    }
    std::cout << histFolder << '\n';
    CalculateSignalEfficiencies(histFolder);

  } else {
    std::cout << "more " << '\n';
    for (std::string histFolder: histFolders) CalculateSignalEfficiencies(histFolder);
  }
  return 0;
}
