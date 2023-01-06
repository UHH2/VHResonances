#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Utils.h"
#include "UHH2/common/include/JetIds.h"
#include "UHH2/VHResonances/include/HiggsToWWHists.h"

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace std;
using namespace uhh2;

HiggsToWWHists::HiggsToWWHists(Context& ctx, const string& dname, const string& condMatch_, const string & condMatchStatus_): HistsBase(ctx, dname), condMatch(condMatch_), condMatchStatus(condMatchStatus_) {

  isInvisible = string2bool(ctx.get("invisiblechannel"));

  h_ZprimeCandidates = ctx.get_handle<vector<ZprimeCandidate>>("ZprimeCandidate");
  h_HDecay = ctx.get_handle<float>("HDecay");
  h_ZDecay = ctx.get_handle<float>("ZDecay");
  h_ZprimeDecay = ctx.get_handle<float>("ZprimeDecay");

  massType = "m"; massPlotName = "mass";
  if (isInvisible){
    massType =  "m_{T}";
    massPlotName = "mass_transversal";
  }

  // book all histograms here
  book_TH1F("sum_event_weights",                "counting experiment", 1, 0.5, 1.5);
  book_TH1F("sum_event_weights_HtoWW",          "counting experiment", 1, 0.5, 1.5);
  book_TH1F("sum_event_weights_Htobb",          "counting experiment", 1, 0.5, 1.5);
  book_TH1F("sum_event_weights_Htocc",          "counting experiment", 1, 0.5, 1.5);
  book_TH1F("sum_event_weights_Htogg",          "counting experiment", 1, 0.5, 1.5);
  book_TH1F("sum_event_weights_Htotautau",      "counting experiment", 1, 0.5, 1.5);
  book_TH1F("sum_event_weights_Helse",          "counting experiment", 1, 0.5, 1.5);
  book_TH1F("sum_event_weights_Zee",            "counting experiment", 1, 0.5, 1.5);
  book_TH1F("sum_event_weights_Zmumu",          "counting experiment", 1, 0.5, 1.5);
  book_TH1F("sum_event_weights_Zelse",          "counting experiment", 1, 0.5, 1.5);
  book_TH1F("sum_event_weights_ZeeHtoWW",       "counting experiment", 1, 0.5, 1.5);
  book_TH1F("sum_event_weights_ZeeHtobb",       "counting experiment", 1, 0.5, 1.5);
  book_TH1F("sum_event_weights_ZeeHtocc",       "counting experiment", 1, 0.5, 1.5);
  book_TH1F("sum_event_weights_ZeeHtogg",       "counting experiment", 1, 0.5, 1.5);
  book_TH1F("sum_event_weights_ZeeHtotautau",   "counting experiment", 1, 0.5, 1.5);
  book_TH1F("sum_event_weights_ZeeHelse",       "counting experiment", 1, 0.5, 1.5);
  book_TH1F("sum_event_weights_ZmumuHtoWW",     "counting experiment", 1, 0.5, 1.5);
  book_TH1F("sum_event_weights_ZmumuHtobb",     "counting experiment", 1, 0.5, 1.5);
  book_TH1F("sum_event_weights_ZmumuHtocc",     "counting experiment", 1, 0.5, 1.5);
  book_TH1F("sum_event_weights_ZmumuHtogg",     "counting experiment", 1, 0.5, 1.5);
  book_TH1F("sum_event_weights_ZmumuHtotautau", "counting experiment", 1, 0.5, 1.5);
  book_TH1F("sum_event_weights_ZmumuHelse",     "counting experiment", 1, 0.5, 1.5);
  book_TH1F("sum_event_weights_ZelseHtoWW",     "counting experiment", 1, 0.5, 1.5);
  book_TH1F("sum_event_weights_ZelseHtobb",     "counting experiment", 1, 0.5, 1.5);
  book_TH1F("sum_event_weights_ZelseHtocc",     "counting experiment", 1, 0.5, 1.5);
  book_TH1F("sum_event_weights_ZelseHtogg",     "counting experiment", 1, 0.5, 1.5);
  book_TH1F("sum_event_weights_ZelseHtotautau", "counting experiment", 1, 0.5, 1.5);
  book_TH1F("sum_event_weights_ZelseHelse",     "counting experiment", 1, 0.5, 1.5);
  book_TH1F("Zprime_number",                    "number of Zprime",    6, -.5, 5.5);

  // Zprime reconstruction

  book_TH1F("Zprime_"+massPlotName+"_rebin1",  massType + "^{Zprime} [GeV/c^{2}]", 9900, 0, 9900);
  book_TH1F("Zprime_"+massPlotName+"_rebin10", massType + "^{Zprime} [GeV/c^{2}]", 990,  0, 9900);
  book_TH1F("Zprime_"+massPlotName+"_rebin30", massType + "^{Zprime} [GeV/c^{2}]", 330,  0, 9900);
  book_TH1F("Zprime_"+massPlotName+"_rebin100",massType + "^{Zprime} [GeV/c^{2}]", 100,  0,10000);

  for (const string & name: {"Zprime","Z","H"}) {
    if (name=="Zprime")   book_TH1F(name+"_"+massPlotName, massType + "^{"+name+"} [GeV/c^{2}]", 40 ,700, 5200);
    else if (name=="Z")   book_TH1F(name+"_mass", "m^{"+name+"} [GeV/c^{2}]",  21, 80.5, 101.5);
    else if ((name=="H")) book_TH1F(name+"_mass",  "m^{"+name+"} [GeV/c^{2}]", 40,  0.,  200.);
    book_TH1F(name+"_pt",     "p_{T}^{"   +name+"} [GeV/c]",500,   0, 5000);
    book_TH1F(name+"_energy", "energy^{"  +name+"} [GeV]",   40, 200, 2200);
    book_TH1F(name+"_eta",    "#eta^{"    +name+"}",        100,  -5,    5);
    book_TH1F(name+"_phi",    "#phi^{"    +name+"}",         50,  -5,    5);
  }
  book_TH1F("l1_pt",          "p_{T}^{l1} [GeV/c]",         500,   0, 5000);
  book_TH1F("l2_pt",          "p_{T}^{l2} [GeV/c]",         500,   0, 5000);

  book_TH1F("delta_eta_H_Z", "#Delta#eta(H,Z)", 50, 0.0, 5.0);
  book_TH1F("delta_phi_H_Z", "#Delta#phi(H,Z)", 35, 0.0, 3.5);
  book_TH1F("delta_R_H_Z",   "#Delta R(H,Z)",   35, 1.5, 5.0);
  book_TH1F("delta_R_ll",    "#Delta R(l,l)",   15, 0.0, 1.5);
  book_TH2F("delta_R_llvsZprime"+massPlotName, ";"+massType+"^{Zprime} [GeV/c^{2}];#Delta#R(l,l)", 330, 0, 9900, 15,  0, 1.5);
  book_TH2F("PtZvsZprime"+massPlotName,        ";"+massType+"^{Zprime} [GeV/c^{2}];#Delta#R(l,l)", 330, 0, 9900, 50,200, 2200);

  book_TH1F("delta_R_subjets",       "#Delta R(sj_1,sj_2)",   20, 0.0, 1.0);
  book_TH1F("delta_R_subjets_Hcc",   "#Delta R(sj_1,sj_2)",   20, 0.0, 1.0);
  book_TH1F("delta_R_subjets_H4q",   "#Delta R(sj_1,sj_2)",   20, 0.0, 1.0);
  book_TH1F("delta_R_subjets_Hbb",   "#Delta R(sj_1,sj_2)",   20, 0.0, 1.0);
  book_TH1F("delta_R_subjets_Helse", "#Delta R(sj_1,sj_2)",   20, 0.0, 1.0);

  book_TH2F("etaphi_H",  ";#eta^{H};#phi^{H}", 50,-2.4, 2.4, 60,-3, 3);
  book_TH2F("etaphi_ll", ";#eta^{lep};#phi^{lep}", 50,-2.4, 2.4, 60,-3, 3);

  for (std::string & disc : discriminators) {
    if (FindInString("tau", disc)) {
      book_TH1F("H_"+disc,"#"+disc+"^{H}",30, -0.01, 1.01);
      book_TH2F("Zprime"+massPlotName+"vs"+disc, ";"+massType+"^{Zprime} [GeV/c^{2}];"+disc, 330, 0, 9900, 30, -0.01, 1.01 );
    }
    else if (FindInString("ParticleNet_mass", disc)) {
      book_TH1F("H_"+disc, disc+"^{H}",40,  0.,  200.);
      book_TH2F("Zprime"+massPlotName+"vs"+disc, ";"+massType+"^{Zprime} [GeV/c^{2}];"+disc, 330, 0, 9900, 40,  0.,  200.);
    }
    else if (FindInString("btag", disc)) {
      book_TH1F("H_"+disc, disc+"^{H}",30, -0.01, 1.01);
      book_TH2F("Zprime"+massPlotName+"vs"+disc, ";"+massType+"^{Zprime} [GeV/c^{2}];"+disc, 330, 0, 9900, 30, -0.01, 1.01 );
    }
    else if (FindInString("chi2", disc)) book_TH1F("H_"+disc, disc+"^{H}",35,0,70);
    else book_TH1F("H_"+disc, disc+"^{H} [GeV/c^{2}]",40,  0.,  200.);
  }

  book_TH1F("H_chitot","#chi_{TOT}",35,0,70);
  book_TH1F("H_chitot1","#chi_{TOT}",35,0,70);
  book_TH1F("H_chitot2","#chi_{TOT}",35,0,70);
  book_TH1F("H_chi_H", "#chi^{H}",35,0,70);
  book_TH1F("H_chi_Z", "#chi^{Z}",35,0,70);

  book_TH1F("H_Match","Match^{H}",19, 0, 19);
  book_TH1F("H_MatchingStatus","MatchingStatus^{H}",10, 0, 10);
  book_TH2F("H_MatchvsH_MatchingStatus",";Match^{H};MatchingStatus^{H}",19, 0, 19, 10, 0, 10);
  for (int i=1;i<20;i++) {
    H1("H_Match")->GetXaxis()->SetBinLabel(i,MatchingToString(i-1).c_str());
    H2("H_MatchvsH_MatchingStatus")->GetXaxis()->SetBinLabel(i,MatchingToString(i-1).c_str());
  }
  for (int i=1;i<11;i++) {
    H1("H_MatchingStatus")->GetXaxis()->SetBinLabel(i,MatchingStatusToString(i-1).c_str());
    H2("H_MatchvsH_MatchingStatus")->GetYaxis()->SetBinLabel(i,MatchingStatusToString(i-1).c_str());
  }

  book_TH1F("btags_DeepCSV","btags_DeepCSV", 4, 0, 4);
  book_TH2F("nsubjet_btags_DeepCSV",";DeepCSV^{WP,H}_{subjet1};DeepCSV^{WP,H}_{subjet2}",4, 0, 4, 4, 0, 4);
  H1("btags_DeepCSV")->GetXaxis()->SetBinLabel(1,"no b-tag");
  H1("btags_DeepCSV")->GetXaxis()->SetBinLabel(2,"loose");
  H1("btags_DeepCSV")->GetXaxis()->SetBinLabel(3,"medium");
  H1("btags_DeepCSV")->GetXaxis()->SetBinLabel(4,"tight");
  H2("nsubjet_btags_DeepCSV")->GetXaxis()->SetBinLabel(1,"no b-tag");
  H2("nsubjet_btags_DeepCSV")->GetXaxis()->SetBinLabel(2,"loose");
  H2("nsubjet_btags_DeepCSV")->GetXaxis()->SetBinLabel(3,"medium");
  H2("nsubjet_btags_DeepCSV")->GetXaxis()->SetBinLabel(4,"tight");
  H2("nsubjet_btags_DeepCSV")->GetYaxis()->SetBinLabel(1,"no b-tag");
  H2("nsubjet_btags_DeepCSV")->GetYaxis()->SetBinLabel(2,"loose");
  H2("nsubjet_btags_DeepCSV")->GetYaxis()->SetBinLabel(3,"medium");
  H2("nsubjet_btags_DeepCSV")->GetYaxis()->SetBinLabel(4,"tight");

  if (!isInvisible){
    book_TH1F("Zprime_ptinv"+massPlotName, "p_{T}^{ll}/"+massType+"(jet,ll)", 10, 0, 1);
  } else {
    book_TH1F("Zprime_ptinv"+massPlotName, "E_{T}^{miss}/"+massType+"(Z')", 10, 0, 1);
  }

  for (const std::string & disc : {"Hcc", "HccMD"}) {
    for (const std::string & sel : {"tot", "Twp", "Mwp", "Lwp"}) {
      auto name = disc+"_ptvsmatch_"+sel;
      book_TH2F(name,";p_{T}^{H} [GeV];MatchingStatus^{H}", 7, 0, 7, 19, 0, 19);
      H2(name)->GetXaxis()->SetBinLabel(1,"200");
      H2(name)->GetXaxis()->SetBinLabel(2,"250");
      H2(name)->GetXaxis()->SetBinLabel(3,"300");
      H2(name)->GetXaxis()->SetBinLabel(4,"350");
      H2(name)->GetXaxis()->SetBinLabel(5,"400");
      H2(name)->GetXaxis()->SetBinLabel(6,"450");
      H2(name)->GetXaxis()->SetBinLabel(7,"500");
      for (int i=1;i<20;i++) H2(name)->GetYaxis()->SetBinLabel(i,MatchingToString(i-1).c_str());
    }
  }

  book_TH1F("HT_event", "HT_event",   50, 0, 3000);
  book_TH1F("ST_event", "ST_event",   50, 0, 3000);
  book_TH1F("HT_Zprime", "HT_Zprime", 50, 0, 3000);
  book_TH1F("ST_Zprime", "ST_Zprime", 50, 0, 3000);

  book_TH2F("ST_ZprimevsZprime"+massPlotName, ";ST_Zprime;"+massType+"^{Zprime} [GeV/c^{2}]", 50, 0, 3000, 330, 0, 9900);
  book_TH2F("ST_ZprimevsDeepBoosted",         ";ST_Zprime;btag_DeepBoosted_H4qvsQCD",         50, 0, 3000,  30, -0.01, 1.01);

  for (std::string & disc : discriminators_subjets) {
    book_TH1F("H_"+disc+"_subjet",       disc+"^{subjet}",                     30, -0.01, 1.01);
    book_TH1F("H_"+disc+"_subjet1",      disc+"^{subjet1}",                    30, -0.01, 1.01);
    book_TH1F("H_"+disc+"_subjet2",      disc+"^{subjet2}",                    30, -0.01, 1.01);
    book_TH1F("H_"+disc+"_subjet21",     disc+"^{subjet2/subjet1}",            30, -0.01, 1.01);
    book_TH2F("H_"+disc+"_subjet12", ";"+disc+"^{subjet1};"+disc+"^{subjet2}", 30, -0.01, 1.01, 30, -0.01, 1.01);
  }

  for (std::string & disc : discriminators_Extra) {
    book_TH1F("H_"+disc, disc+"^{H}",30, -0.01, 1.01);
    book_TH2F("Zprime"+massPlotName+"vs"+disc, ";"+massType+"^{Zprime} [GeV/c^{2}];"+disc, 330, 0, 9900, 30, -0.01, 1.01 );
  }
  book_TH2F("H_btag_DeepBoosted_HbbvsHcc_2D", ";DeepBoosted_Hbb;DeepBoosted_Hcc", 30, -0.01, 1.01, 30, -0.01, 1.01 );

}


void HiggsToWWHists::fill(const Event & event){

  auto weight = event.weight;
  fill_H1("sum_event_weights", 1., weight);

  if (event.is_valid(h_HDecay) && event.is_valid(h_ZprimeDecay)) {
    ZprimeDecay HDec = static_cast<ZprimeDecay>(int(event.get(h_HDecay)));
    ZprimeDecay ZDec = static_cast<ZprimeDecay>(int(event.get(h_ZDecay)));
    ZprimeDecay ZprimeDec = static_cast<ZprimeDecay>(int(event.get(h_ZprimeDecay)));

    if(HDec==HWW)     fill_H1("sum_event_weights_HtoWW",     1., weight);
    if(HDec==Hbb)     fill_H1("sum_event_weights_Htobb",     1., weight);
    if(HDec==Hcc)     fill_H1("sum_event_weights_Htocc",     1., weight);
    if(HDec==Hgg)     fill_H1("sum_event_weights_Htogg",     1., weight);
    if(HDec==Htautau) fill_H1("sum_event_weights_Htotautau", 1., weight);
    if(HDec==Helse)   fill_H1("sum_event_weights_Helse",     1., weight);
    if(ZDec==Zee)     fill_H1("sum_event_weights_Zee",       1., weight);
    if(ZDec==Zmumu)   fill_H1("sum_event_weights_Zmumu",     1., weight);
    if(ZDec==Zelse)   fill_H1("sum_event_weights_Zelse",     1., weight);

    if(ZprimeDec==ZeeHWW)       fill_H1("sum_event_weights_ZeeHtoWW",       1., weight);
    if(ZprimeDec==ZeeHbb)       fill_H1("sum_event_weights_ZeeHtobb",       1., weight);
    if(ZprimeDec==ZeeHcc)       fill_H1("sum_event_weights_ZeeHtocc",       1., weight);
    if(ZprimeDec==ZeeHgg)       fill_H1("sum_event_weights_ZeeHtogg",       1., weight);
    if(ZprimeDec==ZeeHtautau)   fill_H1("sum_event_weights_ZeeHtotautau",   1., weight);
    if(ZprimeDec==ZeeHelse)     fill_H1("sum_event_weights_ZeeHelse",       1., weight);
    if(ZprimeDec==ZmumuHWW)     fill_H1("sum_event_weights_ZmumuHtoWW",     1., weight);
    if(ZprimeDec==ZmumuHbb)     fill_H1("sum_event_weights_ZmumuHtobb",     1., weight);
    if(ZprimeDec==ZmumuHcc)     fill_H1("sum_event_weights_ZmumuHtocc",     1., weight);
    if(ZprimeDec==ZmumuHgg)     fill_H1("sum_event_weights_ZmumuHtogg",     1., weight);
    if(ZprimeDec==ZmumuHtautau) fill_H1("sum_event_weights_ZmumuHtotautau", 1., weight);
    if(ZprimeDec==ZmumuHelse)   fill_H1("sum_event_weights_ZmumuHelse",     1., weight);
    if(ZprimeDec==ZelseHWW)     fill_H1("sum_event_weights_ZelseHtoWW",     1., weight);
    if(ZprimeDec==ZelseHbb)     fill_H1("sum_event_weights_ZelseHtobb",     1., weight);
    if(ZprimeDec==ZelseHcc)     fill_H1("sum_event_weights_ZelseHtocc",     1., weight);
    if(ZprimeDec==ZelseHgg)     fill_H1("sum_event_weights_ZelseHtogg",     1., weight);
    if(ZprimeDec==ZelseHtautau) fill_H1("sum_event_weights_ZelseHtotautau", 1., weight);
    if(ZprimeDec==ZelseHelse)   fill_H1("sum_event_weights_ZelseHelse",     1., weight);
  }

  if (! event.is_valid(h_ZprimeCandidates)) return;

  vector<ZprimeCandidate> ZprimeCandidates = event.get(h_ZprimeCandidates);

  // fill the histograms.

  fill_H1("Zprime_number",  ZprimeCandidates.size(), weight);

  std::string match, matchstatus;

  for (const ZprimeCandidate & cand: ZprimeCandidates) {
    match = cand.has_discriminator("Match")? MatchingToString(FloatToMatching(cand.discriminator("Match"))) : MatchingToString(0);
    matchstatus = cand.has_discriminator("MatchingStatus")? MatchingStatusToString(FloatToMatching(cand.discriminator("MatchingStatus"))) : MatchingStatusToString(0);

    if (condMatch!="") {
      if (condMatch=="else") {
        if (match!="HWWMatch" && match!="HbbMatch" && match!="HZZMatch") continue;
      } else if (condMatch!=match) continue;
    }
    if (condMatchStatus!="" && condMatchStatus!=matchstatus) continue;

    fill_H1("Zprime_"+massPlotName,             cand.Zprime_mass(), weight);
    fill_H1("Zprime_"+massPlotName+"_rebin1",   cand.Zprime_mass(), weight);
    fill_H1("Zprime_"+massPlotName+"_rebin10",  cand.Zprime_mass(), weight);
    fill_H1("Zprime_"+massPlotName+"_rebin30",  cand.Zprime_mass(), weight);
    fill_H1("Zprime_"+massPlotName+"_rebin100", cand.Zprime_mass(), weight);

    fill_H1("Zprime_pt",     cand.pt(),         weight);
    fill_H1("Zprime_energy", cand.energy(),     weight);
    fill_H1("Zprime_eta",    cand.eta(),        weight);
    fill_H1("Zprime_phi",    cand.phi(),        weight);

    fill_H1("Z_mass",        cand.Z().v4().M(), weight);
    fill_H1("Z_pt",          cand.Z().pt(),     weight);
    fill_H1("Z_energy",      cand.Z().energy(), weight);
    fill_H1("Z_eta",         cand.Z().eta(),    weight);
    fill_H1("Z_phi",         cand.Z().phi(),    weight);

    fill_H1("H_mass",        cand.H().v4().M(), weight);
    fill_H1("H_pt",          cand.H().pt(),     weight);
    fill_H1("H_energy",      cand.H().energy(), weight);
    fill_H1("H_eta",         cand.H().eta(),    weight);
    fill_H1("H_phi",         cand.H().phi(),    weight);

    double delta_eta_H_Z = 0;
    double delta_phi_H_Z = 0;
    double delta_R_H_Z   = 0;
    double delta_R_ll = 0;

    delta_phi_H_Z = fabs(deltaPhi(cand.H(), cand.Z()));
    if (!isInvisible) {
      fill_H1("l1_pt", cand.leptons()[0].pt(), weight);
      fill_H1("l2_pt", cand.leptons()[1].pt(), weight);
      delta_eta_H_Z = fabs(cand.H().eta() - cand.Z().eta());
      delta_R_H_Z = deltaR(cand.H(), cand.Z());
      delta_R_ll = deltaR(cand.leptons()[0], cand.leptons()[1]);
      H2("etaphi_ll")->Fill(cand.leptons()[0].eta(), cand.leptons()[0].phi(), weight);
      H2("etaphi_ll")->Fill(cand.leptons()[1].eta(), cand.leptons()[1].phi(), weight);
    }
    H2("etaphi_H")->Fill(cand.H().eta(), cand.H().phi(), weight);

    fill_H1("delta_eta_H_Z", delta_eta_H_Z, weight);
    fill_H1("delta_phi_H_Z", delta_phi_H_Z, weight);
    fill_H1("delta_R_H_Z",   delta_R_H_Z,   weight);
    fill_H1("delta_R_ll",    delta_R_ll,    weight);

    H2("delta_R_llvsZprime"+massPlotName)->Fill(cand.Zprime_mass(), delta_R_ll, weight);

    H2("PtZvsZprime"+massPlotName)->Fill(cand.Zprime_mass(), cand.Z().pt(), weight);

    double delta_R_sj = (cand.H().subjets().size()>=2)? deltaR(cand.H().subjets().at(0),cand.H().subjets().at(1)) : -1;
    fill_H1("delta_R_subjets", delta_R_sj, weight);
    if (match=="HbbMatch")      fill_H1("delta_R_subjets_Hbb", delta_R_sj, weight);
    else if (match=="HWWMatch") fill_H1("delta_R_subjets_H4q", delta_R_sj, weight);
    else if (match=="HccMatch") fill_H1("delta_R_subjets_Hcc", delta_R_sj, weight);
    else fill_H1("delta_R_subjets_Helse", delta_R_sj, weight);

    double chi_H = (cand.H().softdropmass()-HMASS)/HWIDTH;
    double chi_Z = (cand.Z().v4().M()-ZMASS)/ZWIDTH;
    if (isInvisible){ chi_Z=0; } // For the invisible channel, Z is not taken into account for the chi calculations.

    fill_H1("H_chitot", chi_H +chi_Z , weight);
    fill_H1("H_chitot1",TMath::Power(chi_H,2) +TMath::Power(chi_Z,2) , weight);
    fill_H1("H_chitot2",TMath::Sqrt(TMath::Power(chi_H,2) +TMath::Power(chi_Z,2)) , weight);
    fill_H1("H_chi_H", chi_H, weight);
    fill_H1("H_chi_Z", chi_Z, weight);

    if (cand.discriminator("btag_DeepCSV_tight"))  H1("btags_DeepCSV")->Fill("tight",  weight);
    else if (cand.discriminator("btag_DeepCSV_medium")) H1("btags_DeepCSV")->Fill("medium", weight);
    else if (cand.discriminator("btag_DeepCSV_loose"))  H1("btags_DeepCSV")->Fill("loose",  weight);
    else H1("btags_DeepCSV")->Fill("no b-tag",  weight);

    std::vector<std::string> index_tag(2,"no b-tag");
    for (int sj_ind = 0; sj_ind <  (int)cand.discriminator("subjets"); sj_ind++) {
      if (sj_ind>1) continue;
      if ((bool)cand.discriminator("btag_DeepCSV_tight_subjet_"+std::to_string(sj_ind))) index_tag[sj_ind] = "tight";
      else if ((bool)cand.discriminator("btag_DeepCSV_medium_subjet_"+std::to_string(sj_ind))) index_tag[sj_ind] = "medium";
      else if ((bool)cand.discriminator("btag_DeepCSV_loose_subjet_"+std::to_string(sj_ind))) index_tag[sj_ind] = "loose";
      else index_tag[sj_ind] = "no b-tag";
    }
    H2("nsubjet_btags_DeepCSV")->Fill(index_tag[0].c_str(), index_tag[1].c_str(), weight);

    fill_H1("Zprime_ptinv"+massPlotName,  cand.Z().pt()/cand.Zprime_mass(), weight);

    H1("H_Match")->Fill(match.c_str(), weight);
    H1("H_MatchingStatus")->Fill(matchstatus.c_str(), weight);
    H2("H_MatchvsH_MatchingStatus")->Fill(match.c_str(), matchstatus.c_str(), weight);

    double HccvsQCD    = cand.discriminator("btag_DeepBoosted_HccvsQCD");
    double HccvsQCD_MD = cand.discriminator("btag_DeepBoosted_HccvsQCD_MD");
    double H_pt = cand.H().pt();

    std::string H_pt_str = "";
    if (H_pt>500) H_pt_str = "500";
    else if (H_pt>450) H_pt_str = "450";
    else if (H_pt>400) H_pt_str = "400";
    else if (H_pt>350) H_pt_str = "350";
    else if (H_pt>300) H_pt_str = "300";
    else if (H_pt>250) H_pt_str = "250";
    else if (H_pt>200) H_pt_str = "200";

    // Values hard-coded taken from Loukas
    H2("Hcc_ptvsmatch_tot")->Fill(H_pt_str.c_str(), match.c_str(), weight);
    H2("HccMD_ptvsmatch_tot")->Fill(H_pt_str.c_str(), match.c_str(), weight);
    if (HccvsQCD>0.90)      H2("Hcc_ptvsmatch_Twp")->Fill(H_pt_str.c_str(), match.c_str(), weight);
    else if (HccvsQCD>0.80) H2("Hcc_ptvsmatch_Mwp")->Fill(H_pt_str.c_str(), match.c_str(), weight);
    else if (HccvsQCD>0.70) H2("Hcc_ptvsmatch_Lwp")->Fill(H_pt_str.c_str(), match.c_str(), weight);

    if (HccvsQCD_MD>0.90)      H2("HccMD_ptvsmatch_Twp")->Fill(H_pt_str.c_str(), match.c_str(), weight);
    else if (HccvsQCD_MD>0.80) H2("HccMD_ptvsmatch_Mwp")->Fill(H_pt_str.c_str(), match.c_str(), weight);
    else if (HccvsQCD_MD>0.70) H2("HccMD_ptvsmatch_Lwp")->Fill(H_pt_str.c_str(), match.c_str(), weight);

    double HT = cand.H().pt();
    double ST = cand.H().pt()+cand.Z().pt();
    fill_H1("HT_Zprime", HT, weight);
    fill_H1("ST_Zprime", ST, weight);
    H2("ST_ZprimevsZprime"+massPlotName)->Fill( ST, cand.Zprime_mass(), weight);
    H2("ST_ZprimevsDeepBoosted")->Fill(ST, cand.has_discriminator("btag_DeepBoosted_H4qvsQCD")? cand.discriminator("btag_DeepBoosted_H4qvsQCD"): 9999, weight);

    for (std::string & disc : discriminators) {
      double val = 9999;
      if (cand.has_discriminator(disc)) val = cand.discriminator(disc);
      fill_H1("H_"+disc, val, weight);
      if (FindInString("chi2", disc) || FindInString("SDmass", disc)) continue;
      H2("Zprime"+massPlotName+"vs"+disc)->Fill(cand.Zprime_mass(), val, weight);
    }

    int nsubjet = cand.H().subjets().size();
    for (std::string & disc : discriminators_subjets) {
      Jet subjet1, subjet2;
      if (nsubjet>0) subjet1 = cand.H().subjets().at(0);
      if (nsubjet>1) subjet2 = cand.H().subjets().at(1);
      double sub1=9999, sub2=0;
      if (disc=="btag_DeepJet") {                    sub1 = nsubjet>0 ? subjet1.btag_DeepJet() :                     9999; sub2 = nsubjet>1 ? subjet2.btag_DeepJet() : 0; }
      if (disc=="btag_DeepCSV") {                    sub1 = nsubjet>0 ? subjet1.btag_DeepCSV() :                     9999; sub2 = nsubjet>1 ? subjet2.btag_DeepCSV() : 0; }
      if (disc=="btag_DeepFlavour_bb") {             sub1 = nsubjet>0 ? subjet1.btag_DeepFlavour_bb() :              9999; sub2 = nsubjet>1 ? subjet2.btag_DeepFlavour_bb() : 0; }
      if (disc=="btag_DeepFlavour_b") {              sub1 = nsubjet>0 ? subjet1.btag_DeepFlavour_b() :               9999; sub2 = nsubjet>1 ? subjet2.btag_DeepFlavour_b() : 0; }
      if (disc=="btag_DeepFlavour_lepb") {           sub1 = nsubjet>0 ? subjet1.btag_DeepFlavour_lepb() :            9999; sub2 = nsubjet>1 ? subjet2.btag_DeepFlavour_lepb() : 0; }
      if (disc=="btag_DeepFlavour_uds") {            sub1 = nsubjet>0 ? subjet1.btag_DeepFlavour_uds() :             9999; sub2 = nsubjet>1 ? subjet2.btag_DeepFlavour_uds() : 0; }
      if (disc=="btag_DeepFlavour_g") {              sub1 = nsubjet>0 ? subjet1.btag_DeepFlavour_g() :               9999; sub2 = nsubjet>1 ? subjet2.btag_DeepFlavour_g() : 0; }
      if (disc=="btag_DeepFlavour_c") {              sub1 = nsubjet>0 ? subjet1.btag_DeepFlavour_c() :               9999; sub2 = nsubjet>1 ? subjet2.btag_DeepFlavour_c() : 0; }
      fill_H1("H_"+disc+"_subjet", (sub1>sub2)?sub1:sub2, weight);
      fill_H1("H_"+disc+"_subjet1", sub1, weight);
      fill_H1("H_"+disc+"_subjet2", sub2, weight);
      fill_H1("H_"+disc+"_subjet21", sub2/(sub1+sub2), weight);
      H2("H_"+disc+"_subjet12")->Fill(sub1, sub2, weight);
    }

    for (std::string & disc : discriminators_Extra) {
      double val=0;
      if (disc=="btag_DeepBoosted_TvsQCD")   val = cand.H().btag_DeepBoosted_TvsQCD();
      if (disc=="btag_DeepBoosted_WvsQCD")   val = cand.H().btag_DeepBoosted_WvsQCD();
      if (disc=="btag_DeepBoosted_ZvsQCD")   val = cand.H().btag_DeepBoosted_ZvsQCD();
      if (disc=="btag_DeepBoosted_HbbvsQCD") val = cand.H().btag_DeepBoosted_HbbvsQCD();
      if (disc=="btag_MassDecorrelatedDeepBoosted_TvsQCD")    val = cand.H().btag_MassDecorrelatedDeepBoosted_TvsQCD();
      if (disc=="btag_MassDecorrelatedDeepBoosted_WvsQCD")    val = cand.H().btag_MassDecorrelatedDeepBoosted_WvsQCD();
      if (disc=="btag_MassDecorrelatedDeepBoosted_ZHbbvsQCD") val = cand.H().btag_MassDecorrelatedDeepBoosted_ZHbbvsQCD();
      if (disc=="btag_MassDecorrelatedDeepBoosted_ZvsQCD")    val = cand.H().btag_MassDecorrelatedDeepBoosted_ZvsQCD();
      if (disc=="btag_MassDecorrelatedDeepBoosted_ZbbvsQCD")  val = cand.H().btag_MassDecorrelatedDeepBoosted_ZbbvsQCD();
      if (disc=="btag_MassDecorrelatedDeepBoosted_HbbvsQCD")  val = cand.H().btag_MassDecorrelatedDeepBoosted_HbbvsQCD();
      if (disc=="btag_DeepBoosted_HbbvsHcc") val = cand.H().btag_DeepBoosted_probHbb()/(cand.H().btag_DeepBoosted_probHcc()+cand.H().btag_DeepBoosted_probHbb());
      if (disc=="btag_DeepBoosted_HvsQCD") val = (cand.H().btag_DeepBoosted_probHcc()+cand.H().btag_DeepBoosted_probHbb())/(cand.H().btag_DeepBoosted_probHcc()+cand.H().btag_DeepBoosted_probHbb()+GetQCD(cand.H(), false));
      if (disc=="btag_DeepBoosted_ZbbvsQCD") val = cand.H().btag_DeepBoosted_ZbbvsQCD();
      if (disc=="btag_BoostedDoubleSecondaryVertexAK8")  val = cand.H().btag_BoostedDoubleSecondaryVertexAK8();
      if (disc=="btag_BoostedDoubleSecondaryVertexCA15") val = cand.H().btag_BoostedDoubleSecondaryVertexCA15();
      if (disc=="btag_DeepDoubleBvLJet_probHbb") val = cand.H().btag_DeepDoubleBvLJet_probHbb();
      if (disc=="btag_DeepDoubleBvLJet_probQCD") val = cand.H().btag_DeepDoubleBvLJet_probQCD();
      if (disc=="btag_DeepDoubleCvBJet_probHbb") val = cand.H().btag_DeepDoubleCvBJet_probHbb();
      if (disc=="btag_DeepDoubleCvBJet_probHcc") val = cand.H().btag_DeepDoubleCvBJet_probHcc();
      if (disc=="btag_DeepDoubleCvLJet_probHcc") val = cand.H().btag_DeepDoubleCvLJet_probHcc();
      if (disc=="btag_DeepDoubleCvLJet_probQCD") val = cand.H().btag_DeepDoubleCvLJet_probQCD();
      if (disc=="btag_MassIndependentDeepDoubleBvLJet_probHbb") val = cand.H().btag_MassIndependentDeepDoubleBvLJet_probHbb();
      if (disc=="btag_MassIndependentDeepDoubleBvLJet_probQCD") val = cand.H().btag_MassIndependentDeepDoubleBvLJet_probQCD();
      if (disc=="btag_MassIndependentDeepDoubleCvBJet_probHbb") val = cand.H().btag_MassIndependentDeepDoubleCvBJet_probHbb();
      if (disc=="btag_MassIndependentDeepDoubleCvBJet_probHcc") val = cand.H().btag_MassIndependentDeepDoubleCvBJet_probHcc();
      if (disc=="btag_MassIndependentDeepDoubleCvLJet_probHcc") val = cand.H().btag_MassIndependentDeepDoubleCvLJet_probHcc();
      if (disc=="btag_MassIndependentDeepDoubleCvLJet_probQCD") val = cand.H().btag_MassIndependentDeepDoubleCvLJet_probQCD();
      if (disc=="btag_DeepBoosted_probQCDb")      val = cand.H().btag_DeepBoosted_probQCDb();
      if (disc=="btag_DeepBoosted_probQCDbb")     val = cand.H().btag_DeepBoosted_probQCDbb();
      if (disc=="btag_DeepBoosted_probQCDc")      val = cand.H().btag_DeepBoosted_probQCDc();
      if (disc=="btag_DeepBoosted_probQCDcc")     val = cand.H().btag_DeepBoosted_probQCDcc();
      if (disc=="btag_DeepBoosted_probQCDothers") val = cand.H().btag_DeepBoosted_probQCDothers();
      if (disc=="btag_DeepBoosted_probTbqq")      val = cand.H().btag_DeepBoosted_probTbqq();
      if (disc=="btag_DeepBoosted_probTbcq")      val = cand.H().btag_DeepBoosted_probTbcq();
      if (disc=="btag_DeepBoosted_probTbq")       val = cand.H().btag_DeepBoosted_probTbq();
      if (disc=="btag_DeepBoosted_probTbc")       val = cand.H().btag_DeepBoosted_probTbc();
      if (disc=="btag_DeepBoosted_probWqq")       val = cand.H().btag_DeepBoosted_probWqq();
      if (disc=="btag_DeepBoosted_probWcq")       val = cand.H().btag_DeepBoosted_probWcq();
      if (disc=="btag_DeepBoosted_probZcc")       val = cand.H().btag_DeepBoosted_probZcc();
      if (disc=="btag_DeepBoosted_probZqq")       val = cand.H().btag_DeepBoosted_probZqq();
      if (disc=="btag_DeepBoosted_probZbb")       val = cand.H().btag_DeepBoosted_probZbb();
      if (disc=="btag_DeepBoosted_probHbb")       val = cand.H().btag_DeepBoosted_probHbb();
      if (disc=="btag_DeepBoosted_probHcc")       val = cand.H().btag_DeepBoosted_probHcc();
      if (disc=="btag_DeepBoosted_probHqqqq")       val = cand.H().btag_DeepBoosted_probHqqqq();
      if (disc=="btag_DeepBoosted_raw_score_qcd") val = cand.H().btag_DeepBoosted_raw_score_qcd();
      if (disc=="btag_DeepBoosted_raw_score_top") val = cand.H().btag_DeepBoosted_raw_score_top();
      if (disc=="btag_DeepBoosted_raw_score_w")   val = cand.H().btag_DeepBoosted_raw_score_w();
      if (disc=="btag_DeepBoosted_raw_score_z")   val = cand.H().btag_DeepBoosted_raw_score_z();
      if (disc=="btag_DeepBoosted_raw_score_h")   val = cand.H().btag_DeepBoosted_raw_score_h();
      if (disc=="btag_MassDecorrelatedDeepBoosted_bbvsLight")     val = cand.H().btag_MassDecorrelatedDeepBoosted_bbvsLight();
      if (disc=="btag_MassDecorrelatedDeepBoosted_ccvsLight")     val = cand.H().btag_MassDecorrelatedDeepBoosted_ccvsLight();
      if (disc=="btag_MassDecorrelatedDeepBoosted_probHbb")       val = cand.H().btag_MassDecorrelatedDeepBoosted_probHbb();
      if (disc=="btag_MassDecorrelatedDeepBoosted_probQCDc")      val = cand.H().btag_MassDecorrelatedDeepBoosted_probQCDc();
      if (disc=="btag_MassDecorrelatedDeepBoosted_probQCDbb")     val = cand.H().btag_MassDecorrelatedDeepBoosted_probQCDbb();
      if (disc=="btag_MassDecorrelatedDeepBoosted_probTbqq")      val = cand.H().btag_MassDecorrelatedDeepBoosted_probTbqq();
      if (disc=="btag_MassDecorrelatedDeepBoosted_probTbcq")      val = cand.H().btag_MassDecorrelatedDeepBoosted_probTbcq();
      if (disc=="btag_MassDecorrelatedDeepBoosted_probTbq")       val = cand.H().btag_MassDecorrelatedDeepBoosted_probTbq();
      if (disc=="btag_MassDecorrelatedDeepBoosted_probQCDothers") val = cand.H().btag_MassDecorrelatedDeepBoosted_probQCDothers();
      if (disc=="btag_MassDecorrelatedDeepBoosted_probQCDb")      val = cand.H().btag_MassDecorrelatedDeepBoosted_probQCDb();
      if (disc=="btag_MassDecorrelatedDeepBoosted_probTbc")       val = cand.H().btag_MassDecorrelatedDeepBoosted_probTbc();
      if (disc=="btag_MassDecorrelatedDeepBoosted_probWqq")       val = cand.H().btag_MassDecorrelatedDeepBoosted_probWqq();
      if (disc=="btag_MassDecorrelatedDeepBoosted_probQCDcc")     val = cand.H().btag_MassDecorrelatedDeepBoosted_probQCDcc();
      if (disc=="btag_MassDecorrelatedDeepBoosted_probHcc")       val = cand.H().btag_MassDecorrelatedDeepBoosted_probHcc();
      if (disc=="btag_MassDecorrelatedDeepBoosted_probZcc")       val = cand.H().btag_MassDecorrelatedDeepBoosted_probZcc();
      if (disc=="btag_MassDecorrelatedDeepBoosted_proWcq")        val = cand.H().btag_MassDecorrelatedDeepBoosted_proWcq();
      if (disc=="btag_MassDecorrelatedDeepBoosted_probZqq")       val = cand.H().btag_MassDecorrelatedDeepBoosted_probZqq();
      if (disc=="btag_MassDecorrelatedDeepBoosted_probHqqqq")     val = cand.H().btag_MassDecorrelatedDeepBoosted_probHqqqq();
      if (disc=="btag_MassDecorrelatedDeepBoosted_probZbb")       val = cand.H().btag_MassDecorrelatedDeepBoosted_probZbb();

      if (disc=="btag_ParticleNet_mass") {
        val = cand.discriminator(disc);
        fill_H1("H_"+disc, val, weight);
      } else {
        fill_H1("H_"+disc, val, weight);
        H2("Zprime"+massPlotName+"vs"+disc)->Fill(cand.Zprime_mass(), val, weight);
      }
    }
    H2("H_btag_DeepBoosted_HbbvsHcc_2D")->Fill(cand.H().btag_DeepBoosted_probHbb(), cand.H().btag_DeepBoosted_probHcc(), weight);

  }

  double HT = 0;
  // vector<std::vector<TopJet>> topjets = event.get(h_topjets);
  for (const auto & jet : *event.toppuppijets ) HT += jet.pt();
  fill_H1("HT_event", HT, weight);
  double ST = HT;
  for (const auto & electron : *event.electrons ) ST += electron.pt();
  for (const auto & muon : *event.muons ) ST += muon.pt();
  fill_H1("ST_event", ST, weight);
}

HiggsToWWHists::~HiggsToWWHists(){}


/*
#  ######  #### ##     ## ########  ##       ########    ########  #######  ########      ######  ##    ##  ######  ########
# ##    ##  ##  ###   ### ##     ## ##       ##          ##       ##     ## ##     ##    ##    ##  ##  ##  ##    ##    ##
# ##        ##  #### #### ##     ## ##       ##          ##       ##     ## ##     ##    ##         ####   ##          ##
#  ######   ##  ## ### ## ########  ##       ######      ######   ##     ## ########      ######     ##     ######     ##
#       ##  ##  ##     ## ##        ##       ##          ##       ##     ## ##   ##            ##    ##          ##    ##
# ##    ##  ##  ##     ## ##        ##       ##          ##       ##     ## ##    ##     ##    ##    ##    ##    ##    ##
#  ######  #### ##     ## ##        ######## ########    ##        #######  ##     ##     ######     ##     ######     ##
*/



HiggsToWWHistsSlim::HiggsToWWHistsSlim(Context& ctx, const string& dname, const string& condMatch_, const string & condMatchStatus_): HistsBase(ctx, dname), condMatch(condMatch_), condMatchStatus(condMatchStatus_) {

  isInvisible = string2bool(ctx.get("invisiblechannel"));
  h_ZprimeCandidates = ctx.get_handle<vector<ZprimeCandidate>>("ZprimeCandidate");

  massType = "m"; massPlotName = "mass";
  if (isInvisible){
    massType =  "m_T";
    massPlotName = "mass_transversal";
  }
  // book all histograms here
  book_TH1F("sum_event_weights",            "counting experiment", 1, 0.5, 1.5);

  for (const string & name: {"Zprime","Z","H"}) {
    if (name=="Zprime")   book_TH1F(name+"_"+massPlotName, massType + "^"+name+" [GeV]", 40 ,700, 4700);
    else if (name=="Z")   book_TH1F(name+"_mass", "m^"+name+" [GeV]",  21, 80.5, 101.5);
    else if ((name=="H")) book_TH1F(name+"_mass",  "m^"+name+" [GeV]", 40,  0.,  200.);
    book_TH1F(name+"_pt",      "p_{T}^"   +name+" [GeV]",  500,  0,  5000);
    book_TH1F(name+"_energy",  "energy^"  +name+" [GeV]",  40,  200, 2200);
    book_TH1F(name+"_eta",     "#eta"     +name,           100,  -5,    5);
    book_TH1F(name+"_phi",     "#phi"     +name,           50,   -5,    5);
  }

  book_TH1F("H_btag_MassDecorrelatedDeepBoosted_ZHccvsQCD", "ZHccvsQCD", 30, -0.01, 1.01);
  book_TH1F("H_btag_DeepBoosted_H4qvsQCD", "H4qvsQCD", 30, -0.01, 1.01);
  book_TH1F("Zprime_"+massPlotName+"_rebin1",  massType + "^{Zprime} [GeV/c^{2}]", 9900, 0, 9900);
  book_TH1F("Zprime_"+massPlotName+"_rebin10", massType + "^{Zprime} [GeV/c^{2}]", 990,  0, 9900);
  book_TH1F("Zprime_"+massPlotName+"_rebin30", massType + "^{Zprime} [GeV/c^{2}]", 330,  0, 9900);
  book_TH1F("Zprime_"+massPlotName+"_rebin100",massType + "^{Zprime} [GeV/c^{2}]", 100,  0,10000);

}


void HiggsToWWHistsSlim::fill(const Event & event){

  auto weight = event.weight;
  fill_H1("sum_event_weights", 1., weight);
  if (! event.is_valid(h_ZprimeCandidates)) return;

  vector<ZprimeCandidate> ZprimeCandidates = event.get(h_ZprimeCandidates);

  // fill the histograms.
  std::string match, matchstatus;

  for (const ZprimeCandidate & cand: ZprimeCandidates) {
    match = cand.has_discriminator("Match")? MatchingToString(FloatToMatching(cand.discriminator("Match"))) : MatchingToString(0);
    matchstatus = cand.has_discriminator("MatchingStatus")? MatchingStatusToString(FloatToMatching(cand.discriminator("MatchingStatus"))) : MatchingStatusToString(0);

    if (condMatch!="") {
      if (condMatch=="else") {
        if (match!="HWWMatch" && match!="HbbMatch" && match!="HZZMatch") continue;
      } else if (condMatch!=match) continue;
    }
    if (condMatchStatus!="" && condMatchStatus!=matchstatus) continue;
    fill_H1("Zprime_"+massPlotName+"_rebin1",   cand.Zprime_mass(), weight);
    fill_H1("Zprime_"+massPlotName+"_rebin10",  cand.Zprime_mass(), weight);
    fill_H1("Zprime_"+massPlotName+"_rebin30",  cand.Zprime_mass(), weight);
    fill_H1("Zprime_"+massPlotName+"_rebin100", cand.Zprime_mass(), weight);

    fill_H1("Zprime_pt",     cand.pt(),         weight);
    fill_H1("Zprime_energy", cand.energy(),     weight);
    fill_H1("Zprime_eta",    cand.eta(),        weight);
    fill_H1("Zprime_phi",    cand.phi(),        weight);

    fill_H1("Z_mass",        cand.Z().v4().M(), weight);
    fill_H1("Z_pt",          cand.Z().pt(),     weight);
    fill_H1("Z_energy",      cand.Z().energy(), weight);
    fill_H1("Z_eta",         cand.Z().eta(),    weight);
    fill_H1("Z_phi",         cand.Z().phi(),    weight);

    fill_H1("H_mass",        cand.H().v4().M(), weight);
    fill_H1("H_pt",          cand.H().pt(),     weight);
    fill_H1("H_energy",      cand.H().energy(), weight);
    fill_H1("H_eta",         cand.H().eta(),    weight);
    fill_H1("H_phi",         cand.H().phi(),    weight);

    fill_H1("H_btag_MassDecorrelatedDeepBoosted_ZHccvsQCD", cand.H().btag_MassDecorrelatedDeepBoosted_ZHccvsQCD(), weight);
    fill_H1("H_btag_DeepBoosted_H4qvsQCD", cand.H().btag_DeepBoosted_H4qvsQCD(), weight);

  }

}

HiggsToWWHistsSlim::~HiggsToWWHistsSlim(){}
