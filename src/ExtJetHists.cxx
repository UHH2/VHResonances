#include "UHH2/core/include/Event.h"

#include "UHH2/common/include/Utils.h"

#include "UHH2/VHResonances/include/ExtJetHists.h"
#include "UHH2/VHResonances/include/Utils.hpp"
#include "UHH2/VHResonances/include/HistsBase.hpp"

#include "TH1F.h"
#include <type_traits>
#include <iostream>

using namespace std;
using namespace uhh2;

#define MYTAGLOOP(func)\
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
func(btag_DeepBoosted_H4qvsQCD)\
func(btag_MassDecorrelatedDeepBoosted_TvsQCD)\
func(btag_MassDecorrelatedDeepBoosted_ZHccvsQCD)\
func(btag_MassDecorrelatedDeepBoosted_WvsQCD)\
func(btag_MassDecorrelatedDeepBoosted_ZHbbvsQCD)\
func(btag_MassDecorrelatedDeepBoosted_ZvsQCD)\
func(btag_MassDecorrelatedDeepBoosted_ZbbvsQCD)\
func(btag_MassDecorrelatedDeepBoosted_HbbvsQCD)\
func(btag_MassDecorrelatedDeepBoosted_H4qvsQCD)\
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
func(btag_DeepBoosted_probHqqqq)\
func(btag_DeepBoosted_raw_score_qcd)\
func(btag_DeepBoosted_raw_score_top)\
func(btag_DeepBoosted_raw_score_w)\
func(btag_DeepBoosted_raw_score_z)\
func(btag_DeepBoosted_raw_score_h)\
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
func(btag_MassDecorrelatedDeepBoosted_probHqqqq)\
func(btag_MassDecorrelatedDeepBoosted_probZbb)\

const std::string MyString(const std::string & tag) {return tag;};

#define MYTAGBOOK(mytag)\
isLong = MyString(#mytag).find("BoostedDoubleSecondary")!=std::string::npos;\
book_TH1F(MyString(#mytag)+histSuffix,MyString(#mytag)+"^{"+axisSuffix+"}", isLong? 202: 101, isLong? -1.01: -0.01, isLong? 1.01: 1.01);\

#define MYTAGFILL(mytag)\
fill_H1(MyString(#mytag)+histSuffix, jet.mytag(),weight);\

/*
█  ██████  ██████  ███    ██ ███████ ████████ ██████  ██    ██  ██████ ████████  ██████  ██████
█ ██      ██    ██ ████   ██ ██         ██    ██   ██ ██    ██ ██         ██    ██    ██ ██   ██
█ ██      ██    ██ ██ ██  ██ ███████    ██    ██████  ██    ██ ██         ██    ██    ██ ██████
█ ██      ██    ██ ██  ██ ██      ██    ██    ██   ██ ██    ██ ██         ██    ██    ██ ██   ██
█  ██████  ██████  ██   ████ ███████    ██    ██   ██  ██████   ██████    ██     ██████  ██   ██
*/

ExtJetHists::ExtJetHists(Context & ctx, const string & dname, const string & collection_, const unsigned int NumberOfPlottedJets_): HistsBase(ctx, dname), NumberOfPlottedJets(NumberOfPlottedJets_), collection(collection_){
  if (collection.find("top") != string::npos || collection.find("Top") != string::npos || collection.find("hotvrPuppi") != string::npos) isTop = true;
  h_jets = ctx.get_handle<vector<Jet>>(collection);
  h_topjets = ctx.get_handle<vector<TopJet>>(collection);
  h_image = ctx.get_handle<std::vector<tensorflow::Tensor>>("Images");

  book_TH1F("sum_event_weights", "counting experiment", 1, 0.5, 1.5);
  book_TH1F("weights","weights",4000, -2, 2);
  book_TH1F("number","number of jets",21, -.5, 20.5);
  book_TH1F("btags_DeepCSV","DeepCSV^{WP}", 4, 0, 4);
  H1("btags_DeepCSV")->GetXaxis()->SetBinLabel(1,"no b-tag");
  H1("btags_DeepCSV")->GetXaxis()->SetBinLabel(2,"loose");
  H1("btags_DeepCSV")->GetXaxis()->SetBinLabel(3,"medium");
  H1("btags_DeepCSV")->GetXaxis()->SetBinLabel(4,"tight");
  Btag_map["DeepCSV_loose"]   = BTag(BTag::DEEPCSV, BTag::WP_LOOSE);
  Btag_map["DeepCSV_medium"]  = BTag(BTag::DEEPCSV, BTag::WP_MEDIUM);
  Btag_map["DeepCSV_tight"]   = BTag(BTag::DEEPCSV, BTag::WP_TIGHT);

  book_jetHist("_jet", "jet",20,5000,isTop);
  vector<double> minPt {20,20,20,20};
  vector<double> maxPt {5000,1000,500,350};
  for(unsigned int i =0; i<NumberOfPlottedJets; i++){
    if(i<4) book_jetHist("_"+to_string(i+1),axis_suffix[i],minPt[i],maxPt[i],isTop);
    else book_jetHist("_"+to_string(i+1), to_string(i+1)+"-th jet",20,500,isTop);
  }
  book_TH1F("deltaRmin_1", "#Delta R_{min}(1st jet,nearest jet)", 40, 0, 8.0);
  book_TH1F("deltaRmin_2", "#Delta R_{min}(2nd jet,nearest jet)", 40, 0, 8.0);


}

void ExtJetHists::fill(const Event & event){
  if (event.is_valid(h_topjets))    fill_internal(event, event.get(h_topjets));
  else if (event.is_valid(h_jets))  fill_internal(event, event.get(h_jets));
}


/*
█ ██████   ██████   ██████  ██   ██      █████  ███    ██ ██████      ███████ ██ ██      ██
█ ██   ██ ██    ██ ██    ██ ██  ██      ██   ██ ████   ██ ██   ██     ██      ██ ██      ██
█ ██████  ██    ██ ██    ██ █████       ███████ ██ ██  ██ ██   ██     █████   ██ ██      ██
█ ██   ██ ██    ██ ██    ██ ██  ██      ██   ██ ██  ██ ██ ██   ██     ██      ██ ██      ██
█ ██████   ██████   ██████  ██   ██     ██   ██ ██   ████ ██████      ██      ██ ███████ ███████
*/


template <typename T>
void ExtJetHists::fill_internal(const Event & event, vector<T> jets){
  auto weight = event.weight;
  fill_H1("sum_event_weights", 1., weight);
  fill_H1("weights", weight, 1);
  fill_H1("number", jets.size(), weight);
  std::map<std::string, int> btag_map = {{"no b-tag",0},{"loose",0},{"medium",0},{"tight",0}};
  for(unsigned int i = 0; i <jets.size(); i++){
    T& jet = jets[i];
    if (Btag_map["DeepCSV_tight"](jet, event)) btag_map["tight"]+=1;
    else if (Btag_map["DeepCSV_medium"](jet, event)) btag_map["medium"]+=1;
    else if (Btag_map["DeepCSV_loose"](jet, event)) btag_map["loose"]+=1;
    else btag_map["no b-tag"]+=1;
    fill_jetHist(event, "_jet",jet);
    if(i<4) fill_jetHist(event, "_"+to_string(i+1),jet);
    else if (i<NumberOfPlottedJets) fill_jetHist(event, "_"+to_string(i+1),jet);
    if(i < 2){
      auto next_jet = closestParticle(jet, jets);
      auto drmin = next_jet ? deltaR(jet, *next_jet) : numeric_limits<float>::infinity();
      if(i==0) fill_H1("deltaRmin_1", drmin, weight);
      else fill_H1("deltaRmin_2", drmin, weight);
    }
  }

  for (auto [tag,val]: btag_map) H1("btags_DeepCSV")->Fill(tag.c_str(), val*weight);

}

void ExtJetHists::book_jetHist(const string & histSuffix, const string & axisSuffix, double minPt, double maxPt, bool isTop){
  book_TH1F("mass"+histSuffix,"M "+axisSuffix+" [GeV/c^{2}]",100,0,300);
  book_TH1F("mT"+histSuffix,"m_{T} "+axisSuffix+" [GeV/c^{2}]",100,0,1000);
  book_TH1F("pt"+histSuffix,"p_{T} "+axisSuffix+" [GeV]",50,minPt,maxPt);
  book_TH1F("eta"+histSuffix,"#eta "+axisSuffix,100,-5,5);
  book_TH1F("phi"+histSuffix,"#phi "+axisSuffix,50,-M_PI,M_PI);
  book_TH1F("csv"+histSuffix,"csv-disriminator "+axisSuffix,50,0,1);
  book_TH1F("DeepCSV"+histSuffix,"DeepCSV-disriminator "+axisSuffix,50,0,1);
  book_TH1F("flavor"+histSuffix,"flavor "+axisSuffix,200,-100,100);
  book_TH1F("deltaR_lepton"+histSuffix,"#Delta R(lepton,"+axisSuffix+")",40, 0, 8.0);
  book_TH1F("deltaphi_lepton"+histSuffix,"#Delta#phi(lepton,"+axisSuffix+")",50,0,M_PI);
  book_TH1F("deltaR_muon1"+histSuffix,"#Delta R(muon1,"+axisSuffix+")",40, 0, 8.0);
  book_TH1F("deltaphi_muon1"+histSuffix,"#Delta#phi(muon1,"+axisSuffix+")",50,0,M_PI);
  book_TH1F("deltaR_muon2"+histSuffix,"#Delta R(muon2,"+axisSuffix+")",40, 0, 8.0);
  book_TH1F("deltaphi_muon2"+histSuffix,"#Delta#phi(muon2,"+axisSuffix+")",50,0,M_PI);
  book_TH1F("deltaR_ele1"+histSuffix,"#Delta R(ele1,"+axisSuffix+")",40, 0, 8.0);
  book_TH1F("deltaphi_ele1"+histSuffix,"#Delta#phi(ele1,"+axisSuffix+")",50,0,M_PI);
  book_TH1F("deltaR_ele2"+histSuffix,"#Delta R(ele2,"+axisSuffix+")",40, 0, 8.0);
  book_TH1F("deltaphi_ele2"+histSuffix,"#Delta#phi(ele2,"+axisSuffix+")",50,0,M_PI);
  book_TH1F("deltaphi_jet_MET"+histSuffix,"#Delta#phi(E_{T}^{miss},"+axisSuffix+")",50,0,M_PI);
  book_TH1F("ptMET_ptJet_sine"+histSuffix,"2 p_{T}(MET)*p_{T}(jet)*sin(#Delta#phi)/(p_{T}(MET) + p_{T}(Jet)), scalar, ("+axisSuffix+")",40, 0, 8.0);
  book_TH1F("ptMET_ptJet_sine_vect"+histSuffix,"2 p_{T}(MET)*p_{T}(jet)*sin(#Delta#phi)/(p_{T}(MET) + p_{T}(Jet)), vectorial, ("+axisSuffix+")",40, 0, 8.0);
  book_TH1F("Zprime_inv_M_T"+histSuffix, "m_T of Zprime [GeV/c^{2}] ("+axisSuffix+")", 300,  0, 3000);
  book_TH1F("jetArea"+histSuffix,"jetArea^{"+axisSuffix+"}",150,0,15);
  book_TH1F("jetArea_pt200_300"+histSuffix,"jetArea^{"+axisSuffix+",pt(200,300)}",150,0,15);
  book_TH1F("jetArea_pt300_400"+histSuffix,"jetArea^{"+axisSuffix+",pt(300,400)}",150,0,15);
  book_TH2F("jetAreavspt"+histSuffix,";jetArea^{"+axisSuffix+"};p_{T} "+axisSuffix+" [GeV]",150,0,15,50,minPt,maxPt);
  if (isTop) {
    book_TH1F("SDmass"+histSuffix,"SDmass^{"+axisSuffix+"} [GeV/c^{2}]",100,0,300);
    book_TH1F("SDmassvsinvMass"+histSuffix,"SDmass/invMass^{"+axisSuffix+"}",100,-1,2);
    book_TH1F("tau1"+histSuffix,"#tau1^{"+axisSuffix+"}",20,0,1.01);
    book_TH1F("tau2"+histSuffix,"#tau2^{"+axisSuffix+"}",20,0,1.01);
    book_TH1F("tau3"+histSuffix,"#tau3^{"+axisSuffix+"}",20,0,1.01);
    book_TH1F("tau4"+histSuffix,"#tau4^{"+axisSuffix+"}",20,0,1.01);
    book_TH1F("tau1_groomed"+histSuffix,"#tau1_groomed^{"+axisSuffix+"}",20,0,1.01);
    book_TH1F("tau2_groomed"+histSuffix,"#tau2_groomed^{"+axisSuffix+"}",20,0,1.01);
    book_TH1F("tau3_groomed"+histSuffix,"#tau3_groomed^{"+axisSuffix+"}",20,0,1.01);
    book_TH1F("tau4_groomed"+histSuffix,"#tau4_groomed^{"+axisSuffix+"}",20,0,1.01);
    book_TH1F("tau21"+histSuffix,"#tau21^{"+axisSuffix+"}",20,0,1.01);
    book_TH1F("tau31"+histSuffix,"#tau31^{"+axisSuffix+"}",20,0,1.01);
    book_TH1F("tau41"+histSuffix,"#tau41^{"+axisSuffix+"}",20,0,1.01);
    book_TH1F("tau32"+histSuffix,"#tau32^{"+axisSuffix+"}",20,0,1.01);
    book_TH1F("tau42"+histSuffix,"#tau42^{"+axisSuffix+"}",20,0,1.01);
    book_TH1F("tau43"+histSuffix,"#tau43^{"+axisSuffix+"}",20,0,1.01);
    book_TH1F("nsubjet"+histSuffix,"nsubjet^{"+axisSuffix+"}",41, -.5, 40.5);
    // book_TH2F("nsubjet_btags_DeepCSV"+histSuffix,";nsubjet^{"+axisSuffix+"};btags_DeepCSV^{"+axisSuffix+"}",41, -.5, 40.5, 4, 0, 4);
    book_TH2F("nsubjet_btags_DeepCSV"+histSuffix,";DeepCSV^{WP,"+axisSuffix+"}_{subjet1};DeepCSV^{WP,"+axisSuffix+"}_{subjet2}",4, 0, 4, 4, 0, 4);
    H2("nsubjet_btags_DeepCSV"+histSuffix)->GetXaxis()->SetBinLabel(1,"no b-tag");
    H2("nsubjet_btags_DeepCSV"+histSuffix)->GetXaxis()->SetBinLabel(2,"loose");
    H2("nsubjet_btags_DeepCSV"+histSuffix)->GetXaxis()->SetBinLabel(3,"medium");
    H2("nsubjet_btags_DeepCSV"+histSuffix)->GetXaxis()->SetBinLabel(4,"tight");
    H2("nsubjet_btags_DeepCSV"+histSuffix)->GetYaxis()->SetBinLabel(1,"no b-tag");
    H2("nsubjet_btags_DeepCSV"+histSuffix)->GetYaxis()->SetBinLabel(2,"loose");
    H2("nsubjet_btags_DeepCSV"+histSuffix)->GetYaxis()->SetBinLabel(3,"medium");
    H2("nsubjet_btags_DeepCSV"+histSuffix)->GetYaxis()->SetBinLabel(4,"tight");
    book_TH1F("Match"+histSuffix,"Match^{"+axisSuffix+"}",17, 0, 17);
    book_TH1F("MatchingStatus"+histSuffix,"MatchingStatus^{"+axisSuffix+"}",10, 0, 10);
    book_TH2F("MatchvsMatchingStatus"+histSuffix,";Match^{"+axisSuffix+"};MatchingStatus^{"+axisSuffix+"}",17, 0, 17, 10, 0, 10);
    for (int i=1;i<18;i++) {
      H1("Match"+histSuffix)->GetXaxis()->SetBinLabel(i,MatchingToString(i-1).c_str());
      H2("MatchvsMatchingStatus"+histSuffix)->GetXaxis()->SetBinLabel(i,MatchingToString(i-1).c_str());
    }
    for (int i=1;i<11;i++) {
      H1("MatchingStatus"+histSuffix)->GetXaxis()->SetBinLabel(i,MatchingStatusToString(i-1).c_str());
      H2("MatchvsMatchingStatus"+histSuffix)->GetYaxis()->SetBinLabel(i,MatchingStatusToString(i-1).c_str());
    }
    book_TH1F("NN_HWW"+histSuffix,"NN_HWW^{"+axisSuffix+"}",100, 0, 1.01);
    book_TH1F("NN_Hbb"+histSuffix,"NN_Hbb^{"+axisSuffix+"}",100, 0, 1.01);
    book_TH1F("NN_QCD"+histSuffix,"NN_QCD^{"+axisSuffix+"}",100, 0, 1.01);
    book_TH1F("NN_Top"+histSuffix,"NN_Top^{"+axisSuffix+"}",100, 0, 1.01);
    book_TH1F("NN_W"+histSuffix,"NN_W^{"+axisSuffix+"}",100, 0, 1.01);
    book_TH1F("NN_Z"+histSuffix,"NN_Z^{"+axisSuffix+"}",100, 0, 1.01);
    book_TH1F("NN_HWW_1"+histSuffix,"NN_HWW_1^{"+axisSuffix+"}",100, 0, 1.01);
    book_TH1F("NN_Hbb_1"+histSuffix,"NN_Hbb_1^{"+axisSuffix+"}",100, 0, 1.01);
    book_TH1F("NN_QCD_1"+histSuffix,"NN_QCD_1^{"+axisSuffix+"}",100, 0, 1.01);
    book_TH1F("NN_Top_1"+histSuffix,"NN_Top_1^{"+axisSuffix+"}",100, 0, 1.01);
    book_TH1F("NN_W_1"+histSuffix,"NN_W_1^{"+axisSuffix+"}",100, 0, 1.01);
    book_TH1F("NN_Z_1"+histSuffix,"NN_Z_1^{"+axisSuffix+"}",100, 0, 1.01);
    book_TH1F("NN_HWW_2"+histSuffix,"NN_HWW_2^{"+axisSuffix+"}",100, 0, 1.01);
    book_TH1F("NN_Hbb_2"+histSuffix,"NN_Hbb_2^{"+axisSuffix+"}",100, 0, 1.01);
    book_TH1F("NN_QCD_2"+histSuffix,"NN_QCD_2^{"+axisSuffix+"}",100, 0, 1.01);
    book_TH1F("NN_Top_2"+histSuffix,"NN_Top_2^{"+axisSuffix+"}",100, 0, 1.01);
    book_TH1F("NN_W_2"+histSuffix,"NN_W_2^{"+axisSuffix+"}",100, 0, 1.01);
    book_TH1F("NN_Z_2"+histSuffix,"NN_Z_2^{"+axisSuffix+"}",100, 0, 1.01);
    book_TH1F("CNN_HWW"+histSuffix,"CNN_HWW^{"+axisSuffix+"}",100, 0, 1.01);
    book_TH1F("CNN_QCD"+histSuffix,"CNN_QCD^{"+axisSuffix+"}",100, 0, 1.01);
    book_TH2F("NN_HWW_2vsCNN_HWW"+histSuffix,";NN_HWW_2^{"+axisSuffix+"}"+";CNN_HWW^{"+axisSuffix+"}",100, 0, 1.01,100, 0, 1.01);
    book_TH2F("imagept"+histSuffix,";#eta;#phi}",40, -1, 1.01, 40, 1, 1.01);
    book_TH2F("imageCH"+histSuffix,";#eta;#phi}",40, -1, 1.01, 40, 1, 1.01);
    book_TH2F("imageNH"+histSuffix,";#eta;#phi}",40, -1, 1.01, 40, 1, 1.01);
    bool isLong;
    MYTAGLOOP(MYTAGBOOK)
  }
  book_TH2F("jetptvsdeltaphi_dilep"+histSuffix,";#Delta#phi(lepton,"+axisSuffix+"; p_{T} "+axisSuffix+" [GeV]",50, 0, M_PI,50,minPt,maxPt);
  return;
}

template <>
void ExtJetHists::fill_jetHist<Jet>(const Event & event, const string& histSuffix, const Jet& jet){
  auto weight = event.weight;
  fill_H1("mass"+histSuffix, jet.v4().M(), weight);
  fill_H1("mT"+histSuffix, jet.v4().Mt(), weight);
  fill_H1("pt"+histSuffix, jet.pt(), weight);
  fill_H1("eta"+histSuffix, jet.eta(), weight);
  fill_H1("phi"+histSuffix, jet.phi(), weight);
  fill_H1("csv"+histSuffix, jet.btag_combinedSecondaryVertex(), weight);
  fill_H1("DeepCSV"+histSuffix, jet.btag_DeepCSV(), weight);
  fill_H1("flavor"+histSuffix, jet.pdgId(), weight);
  for (const auto & muon : *event.muons ) {
    fill_H1("deltaR_lepton"+histSuffix, deltaR(jet, muon), weight);
    fill_H1("deltaphi_lepton"+histSuffix, deltaPhi(jet, muon), weight);
    for (const auto & muon2 : *event.muons ) {
      if (muon2==muon) continue;
      auto diLep = muon.v4() + muon2.v4();
      auto Dphi = deltaPhi(diLep, jet);
      fill_H2("jetptvsdeltaphi_dilep"+histSuffix, Dphi, jet.pt(), weight);
    }
  }
  for (const auto & electron : *event.electrons ) {
    fill_H1("deltaR_lepton"+histSuffix, deltaR(jet, electron), weight);
    fill_H1("deltaphi_lepton"+histSuffix, deltaPhi(jet, electron), weight);
    for (const auto & electron2 : *event.electrons ) {
      if (electron2==electron) continue;
      auto diLep = electron.v4() + electron2.v4();
      auto Dphi = deltaPhi(diLep, jet);
      fill_H2("jetptvsdeltaphi_dilep"+histSuffix, Dphi, jet.pt(), weight);
    }
  }
  double delta_phi_jet_met = deltaPhi(jet, *event.met);
  fill_H1("deltaphi_jet_MET"+histSuffix, delta_phi_jet_met , weight);
  fill_H1("ptMET_ptJet_sine"+histSuffix, 2*event.met->pt()*jet.pt()*sin(delta_phi_jet_met)/(event.met->pt()+jet.pt()),weight);
  fill_H1("ptMET_ptJet_sine_vect"+histSuffix, 2*event.met->pt()*jet.pt()*sin(delta_phi_jet_met)/((event.met->v4()+jet.v4()).pt()),weight);
  fill_H1("Zprime_inv_M_T"+histSuffix, sqrt(2*event.met->pt() * jet.pt() * (1-cos(delta_phi_jet_met))), weight);
  fill_H1("jetArea"+histSuffix, jet.jetArea(), weight);
  if (jet.pt()>200 || jet.pt()<300) fill_H1("jetArea_pt200_300"+histSuffix,jet.jetArea(), weight);
  if (jet.pt()>300 || jet.pt()<400) fill_H1("jetArea_pt300_400"+histSuffix,jet.jetArea(), weight);
  fill_H2("jetAreavspt"+histSuffix,jet.jetArea(), jet.pt(), weight);
  if (event.muons) {
    if (event.muons->size() > 0) {
      fill_H1("deltaR_muon1"+histSuffix, deltaR(jet, (*event.muons)[0]), weight);
      fill_H1("deltaphi_muon1"+histSuffix, deltaPhi(jet, (*event.muons)[0]), weight);
    }
    if (event.muons->size() > 1) {
      fill_H1("deltaR_muon2"+histSuffix, deltaR(jet, (*event.muons)[1]), weight);
      fill_H1("deltaphi_muon2"+histSuffix, deltaPhi(jet, (*event.muons)[1]), weight);
    }
  }
  if (event.electrons) {
    if (event.electrons->size() > 0) {
      fill_H1("deltaR_muon1"+histSuffix, deltaR(jet, (*event.electrons)[0]), weight);
      fill_H1("deltaphi_muon1"+histSuffix, deltaPhi(jet, (*event.electrons)[0]), weight);
    }
    if (event.electrons->size() > 1) {
      fill_H1("deltaR_muon2"+histSuffix, deltaR(jet, (*event.electrons)[1]), weight);
      fill_H1("deltaphi_muon2"+histSuffix, deltaPhi(jet, (*event.electrons)[1]), weight);
    }
  }
  return;
}

template <>
void ExtJetHists::fill_jetHist<TopJet>(const Event & event, const string& histSuffix, const TopJet& jet){
  auto weight = event.weight;
  fill_jetHist<Jet>(event, histSuffix,jet);
  fill_H1("SDmass"+histSuffix, jet.softdropmass(), weight);
  fill_H1("SDmassvsinvMass"+histSuffix, (jet.softdropmass()/jet.v4().M()<2) ? jet.softdropmass()/jet.v4().M() : 1.8, weight);
  fill_H1("tau1"+histSuffix, jet.tau1(), weight);
  fill_H1("tau2"+histSuffix, jet.tau2(), weight);
  fill_H1("tau3"+histSuffix, jet.tau3(), weight);
  fill_H1("tau4"+histSuffix, jet.tau4(), weight);
  fill_H1("tau1_groomed"+histSuffix, jet.tau1_groomed(), weight);
  fill_H1("tau2_groomed"+histSuffix, jet.tau2_groomed(), weight);
  fill_H1("tau3_groomed"+histSuffix, jet.tau3_groomed(), weight);
  fill_H1("tau4_groomed"+histSuffix, jet.tau4_groomed(), weight);
  fill_H1("tau21"+histSuffix, (jet.tau1()!=0) ? (jet.tau2()/jet.tau1()) : -1, weight);
  fill_H1("tau31"+histSuffix, (jet.tau1()!=0) ? (jet.tau3()/jet.tau1()) : -1, weight);
  fill_H1("tau41"+histSuffix, (jet.tau1()!=0) ? (jet.tau4()/jet.tau1()) : -1, weight);
  fill_H1("tau32"+histSuffix, (jet.tau2()!=0) ? (jet.tau3()/jet.tau2()) : -1, weight);
  fill_H1("tau42"+histSuffix, (jet.tau2()!=0) ? (jet.tau4()/jet.tau2()) : -1, weight);
  fill_H1("tau43"+histSuffix, (jet.tau3()!=0) ? (jet.tau4()/jet.tau3()) : -1, weight);
  fill_H1("nsubjet"+histSuffix, jet.subjets().size(), weight);

  // int tight = 0; int medium = 0; int loose = 0; int notag = 0;
  // for(auto subjet: jet.subjets()){
  //   if (Btag_map["DeepCSV_tight"](subjet, event)) tight+=1;
  //   else if (Btag_map["DeepCSV_medium"](subjet, event)) medium+=1;
  //   else if (Btag_map["DeepCSV_loose"](subjet, event)) loose+=1;
  //   else notag+=1;
  // }
  // fill_H2("nsubjet_btags_DeepCSV"+histSuffix,jet.subjets().size(), 1, notag*weight);
  // fill_H2("nsubjet_btags_DeepCSV"+histSuffix,jet.subjets().size(), 2, loose*weight);
  // fill_H2("nsubjet_btags_DeepCSV"+histSuffix,jet.subjets().size(), 3, medium*weight);
  // fill_H2("nsubjet_btags_DeepCSV"+histSuffix,jet.subjets().size(), 4, tight*weight);


  std::vector<std::string> index_tag(2,"no b-tag"); int sj_ind = 0;
  for(auto subjet: jet.subjets()){
    if (sj_ind>1) continue;
    if (Btag_map["DeepCSV_tight"](subjet, event)) index_tag[sj_ind] = "tight";
    else if (Btag_map["DeepCSV_medium"](subjet, event)) index_tag[sj_ind] = "medium";
    else if (Btag_map["DeepCSV_loose"](subjet, event)) index_tag[sj_ind] = "loose";
    else index_tag[sj_ind] = "no b-tag";
    sj_ind++;
  }

  H2("nsubjet_btags_DeepCSV"+histSuffix)->Fill(index_tag[0].c_str(), index_tag[1].c_str(), weight);


  std::string match = jet.has_tag(TopJet::Matching)? MatchingToString(FloatToMatching(jet)) : MatchingToString(0);
  std::string matchstatus = jet.has_tag(TopJet::MatchingStatus)? MatchingStatusToString(FloatToMatchingStatus(jet)) : MatchingStatusToString(0);
  H1("Match"+histSuffix)->Fill(match.c_str(), weight);
  H1("MatchingStatus"+histSuffix)->Fill(matchstatus.c_str(), weight);
  H2("MatchvsMatchingStatus"+histSuffix)->Fill(match.c_str(),matchstatus.c_str(), weight);
  fill_H1("NN_HWW"+histSuffix,    jet.has_tag(TopJet::NN_HWW)   ? jet.get_tag(TopJet::NN_HWW)   : 9999, weight);
  fill_H1("NN_Hbb"+histSuffix,    jet.has_tag(TopJet::NN_Hbb)   ? jet.get_tag(TopJet::NN_Hbb)   : 9999, weight);
  fill_H1("NN_QCD"+histSuffix,    jet.has_tag(TopJet::NN_QCD)   ? jet.get_tag(TopJet::NN_QCD)   : 9999, weight);
  fill_H1("NN_Top"+histSuffix,    jet.has_tag(TopJet::NN_Top)   ? jet.get_tag(TopJet::NN_Top)   : 9999, weight);
  fill_H1("NN_W"+histSuffix,      jet.has_tag(TopJet::NN_W)     ? jet.get_tag(TopJet::NN_W)     : 9999, weight);
  fill_H1("NN_Z"+histSuffix,      jet.has_tag(TopJet::NN_Z)     ? jet.get_tag(TopJet::NN_Z)     : 9999, weight);
  fill_H1("NN_HWW_1"+histSuffix,  jet.has_tag(TopJet::NN_HWW_1) ? jet.get_tag(TopJet::NN_HWW_1) : 9999, weight);
  fill_H1("NN_Hbb_1"+histSuffix,  jet.has_tag(TopJet::NN_Hbb_1) ? jet.get_tag(TopJet::NN_Hbb_1) : 9999, weight);
  fill_H1("NN_QCD_1"+histSuffix,  jet.has_tag(TopJet::NN_QCD_1) ? jet.get_tag(TopJet::NN_QCD_1) : 9999, weight);
  fill_H1("NN_Top_1"+histSuffix,  jet.has_tag(TopJet::NN_Top_1) ? jet.get_tag(TopJet::NN_Top_1) : 9999, weight);
  fill_H1("NN_W_1"+histSuffix,    jet.has_tag(TopJet::NN_W_1)   ? jet.get_tag(TopJet::NN_W_1)   : 9999, weight);
  fill_H1("NN_Z_1"+histSuffix,    jet.has_tag(TopJet::NN_Z_1)   ? jet.get_tag(TopJet::NN_Z_1)   : 9999, weight);
  fill_H1("NN_HWW_2"+histSuffix,  jet.has_tag(TopJet::NN_HWW_2) ? jet.get_tag(TopJet::NN_HWW_2) : 9999, weight);
  fill_H1("NN_Hbb_2"+histSuffix,  jet.has_tag(TopJet::NN_Hbb_2) ? jet.get_tag(TopJet::NN_Hbb_2) : 9999, weight);
  fill_H1("NN_QCD_2"+histSuffix,  jet.has_tag(TopJet::NN_QCD_2) ? jet.get_tag(TopJet::NN_QCD_2) : 9999, weight);
  fill_H1("NN_Top_2"+histSuffix,  jet.has_tag(TopJet::NN_Top_2) ? jet.get_tag(TopJet::NN_Top_2) : 9999, weight);
  fill_H1("NN_W_2"+histSuffix,    jet.has_tag(TopJet::NN_W_2)   ? jet.get_tag(TopJet::NN_W_2)   : 9999, weight);
  fill_H1("NN_Z_2"+histSuffix,    jet.has_tag(TopJet::NN_Z_2)   ? jet.get_tag(TopJet::NN_Z_2)   : 9999, weight);
  fill_H1("CNN_HWW"+histSuffix,   jet.has_tag(TopJet::CNN_HWW)  ? jet.get_tag(TopJet::CNN_HWW)  : 9999, weight);
  fill_H1("CNN_QCD"+histSuffix,   jet.has_tag(TopJet::CNN_QCD)  ? jet.get_tag(TopJet::CNN_QCD)  : 9999, weight);
  fill_H2("NN_HWW_2vsCNN_HWW"+histSuffix, jet.has_tag(TopJet::NN_HWW_2)  ? jet.get_tag(TopJet::NN_HWW_2)  : 9999, jet.has_tag(TopJet::CNN_HWW)  ? jet.get_tag(TopJet::CNN_HWW)  : 9999, weight);
  if (event.is_valid(h_image)) {
    int index = std::distance( event.get(h_topjets).begin(), std::find(event.get(h_topjets).begin(), event.get(h_topjets).end(), jet));
    auto image = event.get(h_image)[index];
    for (int i = 0; i < 40; i++) {
      for (int j = 0; j < 40; j++) {
        H2("imagept"+histSuffix)->SetBinContent(i+1,j+1, H2("imagept"+histSuffix)->GetBinContent(i+1,j+1)+image.tensor<float, 4>()(0,0,i,j)*weight);
        H2("imageCH"+histSuffix)->SetBinContent(i+1,j+1, H2("imageCH"+histSuffix)->GetBinContent(i+1,j+1)+image.tensor<float, 4>()(0,1,i,j)*weight);
        H2("imageNH"+histSuffix)->SetBinContent(i+1,j+1, H2("imageNH"+histSuffix)->GetBinContent(i+1,j+1)+image.tensor<float, 4>()(0,2,i,j)*weight);
      }
    }
  }
  MYTAGLOOP(MYTAGFILL)

  return;
}
