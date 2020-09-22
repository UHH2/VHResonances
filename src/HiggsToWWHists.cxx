#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Utils.h"
#include "UHH2/common/include/JetIds.h"
#include "UHH2/VHResonances/include/HiggsToWWHists.h"

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace std;
using namespace uhh2;


#define MYJETTAGNNLOOP(func)\
func(btag_BoostedDoubleSecondaryVertexAK8)\
func(btag_BoostedDoubleSecondaryVertexCA15)\
func(btag_DeepDoubleBvLJet_probHbb)\
func(btag_DeepDoubleBvLJet_probQCD)\
func(btag_DeepDoubleCvBJet_probHbb)\
func(btag_DeepDoubleCvBJet_probHcc)\
func(btag_DeepDoubleCvLJet_probHcc)\
func(btag_DeepDoubleCvLJet_probQCD)\
func(btag_MassIndependentDeepDoubleBvLJet_probHbb)\
func(btag_MassIndependentDeepDoubleBvLJet_probQCD)\
func(btag_MassIndependentDeepDoubleCvBJet_probHbb)\
func(btag_MassIndependentDeepDoubleCvBJet_probHcc)\
func(btag_MassIndependentDeepDoubleCvLJet_probHcc)\
func(btag_MassIndependentDeepDoubleCvLJet_probQCD)\
func(btag_DeepBoosted_TvsQCD)\
func(btag_DeepBoosted_WvsQCD)\
func(btag_DeepBoosted_ZvsQCD)\
func(btag_DeepBoosted_ZbbvsQCD)\
func(btag_DeepBoosted_HbbvsQCD)\
func(btag_MassDecorrelatedDeepBoosted_TvsQCD)\
func(btag_MassDecorrelatedDeepBoosted_ZHccvsQCD)\
func(btag_MassDecorrelatedDeepBoosted_WvsQCD)\
func(btag_MassDecorrelatedDeepBoosted_ZHbbvsQCD)\
func(btag_MassDecorrelatedDeepBoosted_ZvsQCD)\
func(btag_MassDecorrelatedDeepBoosted_ZbbvsQCD)\
func(btag_MassDecorrelatedDeepBoosted_HbbvsQCD)\
func(btag_DeepBoosted_probQCDb)\
func(btag_DeepBoosted_probQCDbb)\
func(btag_DeepBoosted_probQCDc)\
func(btag_DeepBoosted_probQCDcc)\
func(btag_DeepBoosted_probQCDothers)\
func(btag_DeepBoosted_probTbqq)\
func(btag_DeepBoosted_probTbcq)\
func(btag_DeepBoosted_probTbq)\
func(btag_DeepBoosted_probTbc)\
func(btag_DeepBoosted_probWqq)\
func(btag_DeepBoosted_probWcq)\
func(btag_DeepBoosted_probZcc)\
func(btag_DeepBoosted_probZqq)\
func(btag_DeepBoosted_probZbb)\
func(btag_DeepBoosted_probHbb)\
func(btag_DeepBoosted_probHcc)\
func(btag_DeepBoosted_raw_score_qcd)\
func(btag_DeepBoosted_raw_score_top)\
func(btag_DeepBoosted_raw_score_w)\
func(btag_DeepBoosted_raw_score_z)\
func(btag_MassDecorrelatedDeepBoosted_bbvsLight)\
func(btag_MassDecorrelatedDeepBoosted_ccvsLight)\
func(btag_MassDecorrelatedDeepBoosted_probHbb)\
func(btag_MassDecorrelatedDeepBoosted_probQCDc)\
func(btag_MassDecorrelatedDeepBoosted_probQCDbb)\
func(btag_MassDecorrelatedDeepBoosted_probTbqq)\
func(btag_MassDecorrelatedDeepBoosted_probTbcq)\
func(btag_MassDecorrelatedDeepBoosted_probTbq)\
func(btag_MassDecorrelatedDeepBoosted_probQCDothers)\
func(btag_MassDecorrelatedDeepBoosted_probQCDb)\
func(btag_MassDecorrelatedDeepBoosted_probTbc)\
func(btag_MassDecorrelatedDeepBoosted_probWqq)\
func(btag_MassDecorrelatedDeepBoosted_probQCDcc)\
func(btag_MassDecorrelatedDeepBoosted_probHcc)\
func(btag_MassDecorrelatedDeepBoosted_probZcc)\
func(btag_MassDecorrelatedDeepBoosted_proWcq)\
func(btag_MassDecorrelatedDeepBoosted_probZqq)\
func(btag_MassDecorrelatedDeepBoosted_probZbb)\



#define JETSNNBOOK(mytag)\
book_TH1F("H_"+MyString(#mytag), MyString(#mytag)+"^{H}",101,0,1.01);\
book_TH1F("H_"+MyString(#mytag)+"_rebin", MyString(#mytag)+"^{H}",30,0,1.02);\
book_TH2F("Zprime"+massPlotName+"vs"+MyString(#mytag), ";"+massType+"^{Zprime} [GeV/c^{2}];"+MyString(#mytag), 330, 0, 9900, 30,0,1.02 );\

#define JETSNNFILL(mytag)\
fill_H1("H_"+MyString(#mytag), cand.H().mytag(), weight);\
fill_H1("H_"+MyString(#mytag)+"_rebin", cand.H().mytag(), weight);\
H2("Zprime"+massPlotName+"vs"+MyString(#mytag))->Fill(cand.Zprime_mass(), cand.H().mytag(), weight);\


#define MYSUBJETTAGLOOP(func)\
func(btag_combinedSecondaryVertex)\
func(btag_combinedSecondaryVertexMVA)\
func(btag_DeepCSV)\
func(btag_DeepFlavour_bb)\
func(btag_DeepFlavour_b)\
func(btag_DeepFlavour_lepb)\
func(btag_DeepFlavour_uds)\
func(btag_DeepFlavour_g)\
func(btag_DeepFlavour_c)\
func(btag_DeepJet)\

#define SUBJETSBOOK(mytag)\
isLong = MyString(#mytag).find("combinedSecondaryVertex")!=std::string::npos;\
book_TH1F("H_"+MyString(#mytag)+"_subjet",   MyString(#mytag)+"^{subjet}",  isLong? 202: 101, isLong? -1.01: -0.01, isLong? 1.01: 1.01);\
book_TH1F("H_"+MyString(#mytag)+"_subjet1",  MyString(#mytag)+"^{subjet1}", isLong? 202: 101, isLong? -1.01: -0.01, isLong? 1.01: 1.01);\
book_TH1F("H_"+MyString(#mytag)+"_subjet2",  MyString(#mytag)+"^{subjet2}", isLong? 202: 101, isLong? -1.01: -0.01, isLong? 1.01: 1.01);\
book_TH1F("H_"+MyString(#mytag)+"_subjet21", MyString(#mytag)+"^{subjet2/subjet1}", isLong? 202: 101, isLong? -1.01: -0.01, isLong? 1.01: 1.01);\
book_TH2F("H_"+MyString(#mytag)+"_subjet12", ";"+MyString(#mytag)+"^{subjet1}"+MyString(#mytag)+"^{subjet2}", isLong? 202: 101, isLong? -1.01: -0.01, isLong? 1.01: 1.01, isLong? 202: 101, isLong? -1.01: -0.01, isLong? 1.01: 1.01);\

#define SUBJETSFILL(mytag)\
sub1 = cand.H().subjets().size()>0 ? cand.H().subjets().at(0).mytag() : 9999;\
sub2 = cand.H().subjets().size()>1 ? cand.H().subjets().at(1).mytag() : 0;\
fill_H1("H_"+MyString(#mytag)+"_subjet", (sub1>sub2)?sub1:sub2, weight);\
fill_H1("H_"+MyString(#mytag)+"_subjet1", sub1, weight);\
fill_H1("H_"+MyString(#mytag)+"_subjet2", sub2, weight);\
fill_H1("H_"+MyString(#mytag)+"_subjet21", sub2/(sub1+sub2), weight);\
H2("H_"+MyString(#mytag)+"_subjet12")->Fill(sub1, sub2, weight);\


HiggsToWWHists::HiggsToWWHists(Context& ctx, const string& dname, const string& condMatch_, const string & condMatchStatus_): HistsBase(ctx, dname), condMatch(condMatch_), condMatchStatus(condMatchStatus_) {

  h_ZprimeCandidates = ctx.get_handle<vector<ZprimeCandidate>>("ZprimeCandidate");
  h_HDecay = ctx.get_handle<float>("HDecay");
  h_ZDecay = ctx.get_handle<float>("ZDecay");
  h_ZprimeDecay = ctx.get_handle<float>("ZprimeDecay");

  massType = "m"; massPlotName = "mass";
  if (string2bool(ctx.get("invisiblechannel"))){
    massType =  "m_T";
    massPlotName = "mass_transversal";
  }

  // book all histograms here
  book_TH1F("sum_event_weights", "counting experiment", 1, 0.5, 1.5);
  book_TH1F("sum_event_weights_HtoWW", "counting experiment", 1, 0.5, 1.5);
  book_TH1F("sum_event_weights_Htobb", "counting experiment", 1, 0.5, 1.5);
  book_TH1F("sum_event_weights_Zee", "counting experiment", 1, 0.5, 1.5);
  book_TH1F("sum_event_weights_Zmumu", "counting experiment", 1, 0.5, 1.5);
  book_TH1F("sum_event_weights_Zelse", "counting experiment", 1, 0.5, 1.5);
  book_TH1F("sum_event_weights_Helse", "counting experiment", 1, 0.5, 1.5);
  book_TH1F("sum_event_weights_ZeeHtoWW", "counting experiment", 1, 0.5, 1.5);
  book_TH1F("sum_event_weights_ZeeHtobb", "counting experiment", 1, 0.5, 1.5);
  book_TH1F("sum_event_weights_ZeeHelse", "counting experiment", 1, 0.5, 1.5);
  book_TH1F("sum_event_weights_ZmumuHtoWW", "counting experiment", 1, 0.5, 1.5);
  book_TH1F("sum_event_weights_ZmumuHtobb", "counting experiment", 1, 0.5, 1.5);
  book_TH1F("sum_event_weights_ZmumuHelse", "counting experiment", 1, 0.5, 1.5);
  book_TH1F("sum_event_weights_ZelseHtoWW", "counting experiment", 1, 0.5, 1.5);
  book_TH1F("sum_event_weights_ZelseHtobb", "counting experiment", 1, 0.5, 1.5);
  book_TH1F("sum_event_weights_ZelseHelse", "counting experiment", 1, 0.5, 1.5);
  book_TH1F("Zprime_number",  "number of Zprime",                21, -.5, 20.5);

  // Zprime reconstruction
  vector<float> bins_Zprime_rebin_full;
  for (float i = 0; i <= 10000; i+=100) bins_Zprime_rebin_full.push_back(i);
  vector<float> bins_Zprime_rebin1;
  for (float i = 0; i <= 2000; i+=100) bins_Zprime_rebin1.push_back(i);
  for (float i = 2500; i <= 3000; i+=500) bins_Zprime_rebin1.push_back(i);
  bins_Zprime_rebin1.push_back(10000);

  // book_TH1F("Zprime_"+massPlotName+"_rebin_full", massType + "^{Zprime} [GeV/c^{2}]", bins_Zprime_rebin_full.size()-1, &bins_Zprime_rebin_full[0]);
  // book_TH1F("Zprime_"+massPlotName+"_rebin1", massType + "^{Zprime} [GeV/c^{2}]", bins_Zprime_rebin1.size()-1, &bins_Zprime_rebin1[0]);
  // book_TH1F("Zprime_"+massPlotName+"_rebin2", massType + "^{Zprime} [GeV/c^{2}]", 10000, 0, 10000);
  book_TH1F("Zprime_"+massPlotName+"_rebin30",massType + "^{Zprime} [GeV/c^{2}]", 330, 0, 9900);

  for (const string & name: {"Zprime","Z","H"}) {
    // bool check = (name.compare("Zprime")) == 0;
    // book_TH1F(name+"_mass",    "m^"       +name+" [GeV/c^{2}]", check?300:40, check?0:70,check? 3000:110);
    if (name=="Zprime")   book_TH1F(name+"_"+massPlotName, massType + "^"+name+" [GeV/c^{2}]", 300,  0, 3000);
    else if (name=="Z")   book_TH1F(name+"_mass", "m^"+name+" [GeV/c^{2}]", 40,  70, 110);
    else if ((name=="H")) book_TH1F(name+"_mass", "m^"+name+" [GeV/c^{2}]", 10,   0, 200);
    if ((name=="H")) book_TH1F(name+"_mass1GeV",  "m^"+name+" [GeV/c^{2}]", 200,  0, 200);
    book_TH1F(name+"_pt",      "p_{T}^"   +name+" [GeV]",       bins_Zprime_rebin_full.size()-1, &bins_Zprime_rebin_full[0]);
    // book_TH1F(name+"_energy",  "energy^"  +name+" [GeV]",       bins_Zprime_rebin1.size()-1, &bins_Zprime_rebin1[0]);
    // book_TH1F(name+"_energy",  "energy^"  +name+" [GeV]",       1000,0,10000);
    // book_TH1F(name+"_eta",     "#eta"     +name,                100,-5,5);
    // book_TH1F(name+"_phi",     "#phi"     +name,                50,-M_PI,M_PI);
  }

  book_TH1F("delta_phi_H_Z","#Delta#phi(H,Z)",50,0,M_PI);

  for (std::string & disc : discriminators) {
    if (FindInString("tau", disc)) {
      book_TH1F("H_"+disc,"#"+disc+"^{H}",101,0,1.01);
      book_TH1F("H_"+disc+"_rebin","#"+disc+"^{H}",30,0,1.02);
      book_TH2F("Zprime"+massPlotName+"vs"+disc, ";"+massType+"^{Zprime} [GeV/c^{2}];"+disc, 330, 0, 9900, 30,0,1.02 );
    }
    else if (FindInString("btag", disc) || FindInString("NN", disc) || FindInString("DCL", disc) ) {
      book_TH1F("H_"+disc, disc+"^{H}",101,0,1.01);
      book_TH1F("H_"+disc+"_rebin", disc+"^{H}",30,0,1.02);
      book_TH2F("Zprime"+massPlotName+"vs"+disc, ";"+massType+"^{Zprime} [GeV/c^{2}];"+disc, 330, 0, 9900, 30,0,1.02 );
    }
    else if (FindInString("chi2", disc)) book_TH1F("H_"+disc, disc+"^{H}",35,0,70);
    else book_TH1F("H_"+disc, disc+"^{H} [GeV/c^{2}]",100,0,300);
  }

  book_TH1F("H_Match","Match^{H}",17, 0, 17);
  book_TH1F("H_MatchingStatus","MatchingStatus^{H}",10, 0, 10);
  book_TH2F("H_MatchvsH_MatchingStatus",";Match^{H};MatchingStatus^{H}",17, 0, 17, 10, 0, 10);
  for (int i=1;i<18;i++) {
    H1("H_Match")->GetXaxis()->SetBinLabel(i,MatchingToString(i-1).c_str());
    H2("H_MatchvsH_MatchingStatus")->GetXaxis()->SetBinLabel(i,MatchingToString(i-1).c_str());
  }
  for (int i=1;i<11;i++) {
    H1("H_MatchingStatus")->GetXaxis()->SetBinLabel(i,MatchingStatusToString(i-1).c_str());
    H2("H_MatchvsH_MatchingStatus")->GetYaxis()->SetBinLabel(i,MatchingStatusToString(i-1).c_str());
  }

  book_TH1F("H_NN_HWWvsQCD","NN_{HWWvsQCD}^{H}",100,0,1.01);
  book_TH1F("H_NN_HWWvsQCD_1","NN_{HWWvsQCD_1}^{H}",100,0,1.01);
  book_TH1F("H_NN_HWWvsQCD_2","NN_{HWWvsQCD_2}^{H}",100,0,1.01);

  book_TH1F("btags_DeepCSV","btags_DeepCSV", 4, 0, 4);
  book_TH2F("nsubjet_btags_DeepCSV",";DeepCSV^{WP,H}_{subjet1};DeepCSV^{WP,H}_{subjet2}",4, 0, 4, 4, 0, 4);
  // book_TH2F("nsubjet_btags_DeepCSV",";nsubjet^{H};btags_DeepCSV^{H}",41, -.5, 40.5, 4, 0, 4);
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

  book_TH1F("Zprime_ptinv"+massPlotName,  "p_{T}^{ll}/"+massType+"(jet,ll)",         200, 0, 2);
  book_TH1F("HT_event", "HT_event", 300, 0, 3000);
  book_TH1F("ST_event", "ST_event", 300, 0, 3000);
  book_TH1F("HT_Zprime", "HT_Zprime", 300, 0, 3000);
  book_TH1F("ST_Zprime", "ST_Zprime", 300, 0, 3000);

  book_TH2F("ST_ZprimevsZprime"+massPlotName,      ";ST_Zprime;"+massType+"^{Zprime} [GeV/c^{2}]", 300, 0, 3000,330, 0, 9900);
  book_TH2F("ST_ZprimevsDeepBoosted",              ";ST_Zprime;btag_DeepBoosted_H4qvsQCD", 300, 0, 3000, 30,0,1.02 );

  bool isLong;
  MYSUBJETTAGLOOP(SUBJETSBOOK)
  // MYJETTAGNNLOOP(JETSNNBOOK)
  for (std::string & disc : discriminators_Extra) {
    book_TH1F("H_"+disc, disc+"^{H}",101,0,1.01);
    book_TH1F("H_"+disc+"_rebin", disc+"^{H}",30,0,1.02);
    book_TH2F("Zprime"+massPlotName+"vs"+disc, ";"+massType+"^{Zprime} [GeV/c^{2}];"+disc, 330, 0, 9900, 30,0,1.02 );
  }
  book_TH2F("H_btag_DeepBoosted_HbbvsHcc2", ";DeepBoosted_Hbb; DeepBoosted_Hcc", 101,0,1.01, 101,0,1.01 );

}


void HiggsToWWHists::fill(const Event & event){

  auto weight = event.weight;
  fill_H1("sum_event_weights", 1., weight);

  if (event.is_valid(h_HDecay) && event.is_valid(h_ZprimeDecay)) {
    ZprimeDecay HDec = static_cast<ZprimeDecay>(int(event.get(h_HDecay)));
    ZprimeDecay ZDec = static_cast<ZprimeDecay>(int(event.get(h_ZDecay)));
    ZprimeDecay ZprimeDec = static_cast<ZprimeDecay>(int(event.get(h_ZprimeDecay)));

    if(HDec==HWW)   fill_H1("sum_event_weights_HtoWW", 1., weight);
    if(HDec==Hbb)   fill_H1("sum_event_weights_Htobb", 1., weight);
    if(HDec==Helse) fill_H1("sum_event_weights_Helse", 1., weight);
    if(ZDec==Zee)   fill_H1("sum_event_weights_Zee", 1., weight);
    if(ZDec==Zmumu) fill_H1("sum_event_weights_Zmumu", 1., weight);
    if(ZDec==Zelse) fill_H1("sum_event_weights_Zelse", 1., weight);

    if(ZprimeDec==ZeeHWW)     fill_H1("sum_event_weights_ZeeHtoWW",   1., weight);
    if(ZprimeDec==ZeeHbb)     fill_H1("sum_event_weights_ZeeHtobb",   1., weight);
    if(ZprimeDec==ZeeHelse)   fill_H1("sum_event_weights_ZeeHelse",   1., weight);
    if(ZprimeDec==ZmumuHWW)   fill_H1("sum_event_weights_ZmumuHtoWW", 1., weight);
    if(ZprimeDec==ZmumuHbb)   fill_H1("sum_event_weights_ZmumuHtobb", 1., weight);
    if(ZprimeDec==ZmumuHelse) fill_H1("sum_event_weights_ZmumuHelse", 1., weight);
    if(ZprimeDec==ZelseHWW)   fill_H1("sum_event_weights_ZelseHtoWW", 1., weight);
    if(ZprimeDec==ZelseHbb)   fill_H1("sum_event_weights_ZelseHtobb", 1., weight);
    if(ZprimeDec==ZelseHelse) fill_H1("sum_event_weights_ZelseHelse", 1., weight);
  }

  if (! event.is_valid(h_ZprimeCandidates)) return;

  vector<ZprimeCandidate> ZprimeCandidates = event.get(h_ZprimeCandidates);

  // fill the histograms.

  fill_H1("Zprime_number",  ZprimeCandidates.size(),weight);

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

    fill_H1("Zprime_"+massPlotName,               cand.Zprime_mass(),weight);
    // fill_H1("Zprime_"+massPlotName+"_rebin_full", cand.Zprime_mass(),weight);
    // fill_H1("Zprime_"+massPlotName+"_rebin1",     cand.Zprime_mass(),weight);
    // fill_H1("Zprime_"+massPlotName+"_rebin2",     cand.Zprime_mass(),weight);
    fill_H1("Zprime_"+massPlotName+"_rebin30",    cand.Zprime_mass(),weight);
    fill_H1("Zprime_pt",              cand.pt(),weight);
    // fill_H1("Zprime_energy",          cand.energy(),weight);
    // fill_H1("Zprime_eta",             cand.eta(),weight);
    // fill_H1("Zprime_phi",             cand.phi(),weight);

    fill_H1("Z_mass",     cand.Z().v4().M(),weight);
    fill_H1("Z_pt",       cand.Z().pt(),weight);
    // fill_H1("Z_energy",   cand.Z().energy(),weight);
    // fill_H1("Z_eta",      cand.Z().eta(),weight);
    // fill_H1("Z_phi",      cand.Z().phi(),weight);

    fill_H1("H_mass",     cand.H().v4().M(),weight);
    fill_H1("H_mass1GeV", cand.H().v4().M(),weight);
    fill_H1("H_pt",       cand.H().pt(),weight);
    // fill_H1("H_energy",   cand.H().energy(),weight);
    // fill_H1("H_eta",      cand.H().eta(),weight);
    // fill_H1("H_phi",      cand.H().phi(),weight);

    double delta_phi_H_Z = fabs(deltaPhi(cand.H(), cand.Z()));
    fill_H1("delta_phi_H_Z", delta_phi_H_Z , weight);

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

    fill_H1("Zprime_ptinv"+massPlotName,  cand.Z().pt()/cand.Zprime_mass(),weight);

    for (std::string & disc : discriminators) {
      fill_H1("H_"+disc, cand.has_discriminator(disc)? cand.discriminator(disc): 9999 ,weight);
      if (FindInString("chi2", disc) || FindInString("SDmass", disc)) continue;
      fill_H1("H_"+disc+"_rebin", cand.has_discriminator(disc)? cand.discriminator(disc): 9999 ,weight);
      H2("Zprime"+massPlotName+"vs"+disc)->Fill(cand.Zprime_mass(), cand.has_discriminator(disc)? cand.discriminator(disc): 9999, weight);
    }

    fill_H1("H_NN_HWWvsQCD", (cand.has_discriminator("NN_HWW")&&cand.has_discriminator("NN_QCD"))? cand.discriminator("NN_HWW")/(cand.discriminator("NN_HWW")+cand.discriminator("NN_QCD")): 9999 ,weight);
    fill_H1("H_NN_HWWvsQCD_1", (cand.has_discriminator("NN_HWW_1")&&cand.has_discriminator("NN_QCD_1"))? cand.discriminator("NN_HWW_1")/(cand.discriminator("NN_HWW_1")+cand.discriminator("NN_QCD")): 9999 ,weight);
    fill_H1("H_NN_HWWvsQCD_2", (cand.has_discriminator("NN_HWW_2")&&cand.has_discriminator("NN_QCD_2"))? cand.discriminator("NN_HWW_2")/(cand.discriminator("NN_HWW_2")+cand.discriminator("NN_QCD")): 9999 ,weight);
    H1("H_Match")->Fill(match.c_str(), weight);
    H1("H_MatchingStatus")->Fill(matchstatus.c_str(), weight);
    H2("H_MatchvsH_MatchingStatus")->Fill(match.c_str(), matchstatus.c_str(), weight);

    double HT = cand.H().pt();
    double ST = cand.H().pt()+cand.Z().pt();
    fill_H1("HT_Zprime", HT, weight);
    fill_H1("ST_Zprime", ST, weight);
    H2("ST_ZprimevsZprime"+massPlotName)->Fill( ST, cand.Zprime_mass(), weight);
    H2("ST_ZprimevsDeepBoosted")->Fill(ST, cand.has_discriminator("btag_DeepBoosted_H4qvsQCD")? cand.discriminator("btag_DeepBoosted_H4qvsQCD"): 9999, weight);
    double sub1, sub2;
    MYSUBJETTAGLOOP(SUBJETSFILL)
    // MYJETTAGNNLOOP(JETSNNFILL)

    for (std::string & disc : discriminators_Extra) {
      double val=0;
      if (disc=="btag_BoostedDoubleSecondaryVertexAK8") val = cand.H().btag_BoostedDoubleSecondaryVertexAK8();
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
      if (disc=="btag_DeepBoosted_TvsQCD") val = cand.H().btag_DeepBoosted_TvsQCD();
      if (disc=="btag_DeepBoosted_WvsQCD") val = cand.H().btag_DeepBoosted_WvsQCD();
      if (disc=="btag_DeepBoosted_ZvsQCD") val = cand.H().btag_DeepBoosted_ZvsQCD();
      if (disc=="btag_DeepBoosted_ZbbvsQCD") val = cand.H().btag_DeepBoosted_ZbbvsQCD();
      if (disc=="btag_DeepBoosted_HbbvsQCD") val = cand.H().btag_DeepBoosted_HbbvsQCD();
      if (disc=="btag_MassDecorrelatedDeepBoosted_TvsQCD") val = cand.H().btag_MassDecorrelatedDeepBoosted_TvsQCD();
      if (disc=="btag_MassDecorrelatedDeepBoosted_ZHccvsQCD") val = cand.H().btag_MassDecorrelatedDeepBoosted_ZHccvsQCD();
      if (disc=="btag_MassDecorrelatedDeepBoosted_WvsQCD") val = cand.H().btag_MassDecorrelatedDeepBoosted_WvsQCD();
      if (disc=="btag_MassDecorrelatedDeepBoosted_ZHbbvsQCD") val = cand.H().btag_MassDecorrelatedDeepBoosted_ZHbbvsQCD();
      if (disc=="btag_MassDecorrelatedDeepBoosted_ZvsQCD") val = cand.H().btag_MassDecorrelatedDeepBoosted_ZvsQCD();
      if (disc=="btag_MassDecorrelatedDeepBoosted_ZbbvsQCD") val = cand.H().btag_MassDecorrelatedDeepBoosted_ZbbvsQCD();
      if (disc=="btag_MassDecorrelatedDeepBoosted_HbbvsQCD") val = cand.H().btag_MassDecorrelatedDeepBoosted_HbbvsQCD();
      if (disc=="btag_DeepBoosted_probQCDb") val = cand.H().btag_DeepBoosted_probQCDb();
      if (disc=="btag_DeepBoosted_probQCDbb") val = cand.H().btag_DeepBoosted_probQCDbb();
      if (disc=="btag_DeepBoosted_probQCDc") val = cand.H().btag_DeepBoosted_probQCDc();
      if (disc=="btag_DeepBoosted_probQCDcc") val = cand.H().btag_DeepBoosted_probQCDcc();
      if (disc=="btag_DeepBoosted_probQCDothers") val = cand.H().btag_DeepBoosted_probQCDothers();
      if (disc=="btag_DeepBoosted_probTbqq") val = cand.H().btag_DeepBoosted_probTbqq();
      if (disc=="btag_DeepBoosted_probTbcq") val = cand.H().btag_DeepBoosted_probTbcq();
      if (disc=="btag_DeepBoosted_probTbq") val = cand.H().btag_DeepBoosted_probTbq();
      if (disc=="btag_DeepBoosted_probTbc") val = cand.H().btag_DeepBoosted_probTbc();
      if (disc=="btag_DeepBoosted_probWqq") val = cand.H().btag_DeepBoosted_probWqq();
      if (disc=="btag_DeepBoosted_probWcq") val = cand.H().btag_DeepBoosted_probWcq();
      if (disc=="btag_DeepBoosted_probZcc") val = cand.H().btag_DeepBoosted_probZcc();
      if (disc=="btag_DeepBoosted_probZqq") val = cand.H().btag_DeepBoosted_probZqq();
      if (disc=="btag_DeepBoosted_probZbb") val = cand.H().btag_DeepBoosted_probZbb();
      if (disc=="btag_DeepBoosted_probHbb") val = cand.H().btag_DeepBoosted_probHbb();
      if (disc=="btag_DeepBoosted_probHcc") val = cand.H().btag_DeepBoosted_probHcc();
      if (disc=="btag_DeepBoosted_raw_score_qcd") val = cand.H().btag_DeepBoosted_raw_score_qcd();
      if (disc=="btag_DeepBoosted_raw_score_top") val = cand.H().btag_DeepBoosted_raw_score_top();
      if (disc=="btag_DeepBoosted_raw_score_w") val = cand.H().btag_DeepBoosted_raw_score_w();
      if (disc=="btag_DeepBoosted_raw_score_z") val = cand.H().btag_DeepBoosted_raw_score_z();
      if (disc=="btag_MassDecorrelatedDeepBoosted_bbvsLight") val = cand.H().btag_MassDecorrelatedDeepBoosted_bbvsLight();
      if (disc=="btag_MassDecorrelatedDeepBoosted_ccvsLight") val = cand.H().btag_MassDecorrelatedDeepBoosted_ccvsLight();
      if (disc=="btag_MassDecorrelatedDeepBoosted_probHbb") val = cand.H().btag_MassDecorrelatedDeepBoosted_probHbb();
      if (disc=="btag_MassDecorrelatedDeepBoosted_probQCDc") val = cand.H().btag_MassDecorrelatedDeepBoosted_probQCDc();
      if (disc=="btag_MassDecorrelatedDeepBoosted_probQCDbb") val = cand.H().btag_MassDecorrelatedDeepBoosted_probQCDbb();
      if (disc=="btag_MassDecorrelatedDeepBoosted_probTbqq") val = cand.H().btag_MassDecorrelatedDeepBoosted_probTbqq();
      if (disc=="btag_MassDecorrelatedDeepBoosted_probTbcq") val = cand.H().btag_MassDecorrelatedDeepBoosted_probTbcq();
      if (disc=="btag_MassDecorrelatedDeepBoosted_probTbq") val = cand.H().btag_MassDecorrelatedDeepBoosted_probTbq();
      if (disc=="btag_MassDecorrelatedDeepBoosted_probQCDothers") val = cand.H().btag_MassDecorrelatedDeepBoosted_probQCDothers();
      if (disc=="btag_MassDecorrelatedDeepBoosted_probQCDb") val = cand.H().btag_MassDecorrelatedDeepBoosted_probQCDb();
      if (disc=="btag_MassDecorrelatedDeepBoosted_probTbc") val = cand.H().btag_MassDecorrelatedDeepBoosted_probTbc();
      if (disc=="btag_MassDecorrelatedDeepBoosted_probWqq") val = cand.H().btag_MassDecorrelatedDeepBoosted_probWqq();
      if (disc=="btag_MassDecorrelatedDeepBoosted_probQCDcc") val = cand.H().btag_MassDecorrelatedDeepBoosted_probQCDcc();
      if (disc=="btag_MassDecorrelatedDeepBoosted_probHcc") val = cand.H().btag_MassDecorrelatedDeepBoosted_probHcc();
      if (disc=="btag_MassDecorrelatedDeepBoosted_probZcc") val = cand.H().btag_MassDecorrelatedDeepBoosted_probZcc();
      if (disc=="btag_MassDecorrelatedDeepBoosted_proWcq") val = cand.H().btag_MassDecorrelatedDeepBoosted_proWcq();
      if (disc=="btag_MassDecorrelatedDeepBoosted_probZqq") val = cand.H().btag_MassDecorrelatedDeepBoosted_probZqq();
      if (disc=="btag_MassDecorrelatedDeepBoosted_probZbb") val = cand.H().btag_MassDecorrelatedDeepBoosted_probZbb();

      if (disc=="btag_DeepBoosted_HbbvsHcc") val = cand.H().btag_DeepBoosted_probHbb()/(cand.H().btag_DeepBoosted_probHcc()+cand.H().btag_DeepBoosted_probHbb());
      if (disc=="btag_DeepBoosted_HvsQCD") val = (cand.H().btag_DeepBoosted_probHcc()+cand.H().btag_DeepBoosted_probHbb())/(cand.H().btag_DeepBoosted_probQCDb()+cand.H().btag_DeepBoosted_probQCDbb()+cand.H().btag_DeepBoosted_probQCDc()+cand.H().btag_DeepBoosted_probQCDcc()+cand.H().btag_DeepBoosted_probQCDothers());

      fill_H1("H_"+disc, val, weight);
      fill_H1("H_"+disc+"_rebin", val, weight);
      H2("Zprime"+massPlotName+"vs"+disc)->Fill(cand.Zprime_mass(), val, weight);
    }
    H2("H_btag_DeepBoosted_HbbvsHcc2")->Fill(cand.H().btag_DeepBoosted_probHbb(), cand.H().btag_DeepBoosted_probHcc(), weight);

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

  h_ZprimeCandidates = ctx.get_handle<vector<ZprimeCandidate>>("ZprimeCandidate");

  massType = "m"; massPlotName = "mass";
  if (string2bool(ctx.get("invisiblechannel"))){
    massType =  "m_T";
    massPlotName = "mass_transversal";
  }
  // book all histograms here
  book_TH1F("Zprime_"+massPlotName+"_rebin30",massType + "^{Zprime} [GeV/c^{2}]", 330, 0, 9900);

}


void HiggsToWWHistsSlim::fill(const Event & event){

  auto weight = event.weight;
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
    fill_H1("Zprime_"+massPlotName+"_rebin30",    cand.Zprime_mass(),weight);

  }

}

HiggsToWWHistsSlim::~HiggsToWWHistsSlim(){}
