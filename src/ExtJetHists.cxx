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
  book_TH1F("mass"+histSuffix,"M^{"+axisSuffix+"} [GeV/c^{2}]",100,0,300);
  book_TH1F("mT"+histSuffix,"m_{T}^{"+axisSuffix+"} [GeV/c^{2}]",100,0,1000);
  book_TH1F("pt"+histSuffix,"p_{T}^{"+axisSuffix+"} [GeV/c]",50,minPt,maxPt);
  book_TH1F("eta"+histSuffix,"#eta^{"+axisSuffix+"}",100,-5,5);
  book_TH1F("phi"+histSuffix,"#phi^{"+axisSuffix+"}",50,-M_PI,M_PI);
  book_TH2F("etaphi"+histSuffix,";#eta^"+axisSuffix+";#phi^"+axisSuffix, 500, -5, 5, 500,-3, 3);
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
  book_TH1F("jetCHF"+histSuffix,"jetCHF^{"+axisSuffix+"}",100,0,1.);
  book_TH1F("jetNEF"+histSuffix,"jetNEF^{"+axisSuffix+"}",100,0,1.);
  book_TH1F("jetNHF"+histSuffix,"jetNHF^{"+axisSuffix+"}",100,0,1.);
  book_TH1F("jetCEF"+histSuffix,"jetCEF^{"+axisSuffix+"}",100,0,1.);
  book_TH1F("jet_muonEnergyFraction"+histSuffix,"jet_muonEnergyFraction^{"+axisSuffix+"}",100,0,1.);
  book_TH1F("jet_photonEnergyFraction"+histSuffix,"jet_photonEnergyFraction^{"+axisSuffix+"}",100,0,1.);
  book_TH1F("jet_chargedMultiplicity"+histSuffix,"jet_chargedMultiplicity^{"+axisSuffix+"}",100,0,1.);
  book_TH1F("jet_neutralMultiplicity"+histSuffix,"jet_neutralMultiplicity^{"+axisSuffix+"}",100,0,1.);
  book_TH1F("jet_muonMultiplicity"+histSuffix,"jet_muonMultiplicity^{"+axisSuffix+"}",100,0,1.);
  book_TH1F("jet_electronMultiplicity"+histSuffix,"jet_electronMultiplicity^{"+axisSuffix+"}",100,0,1.);
  book_TH1F("jet_photonMultiplicity"+histSuffix,"jet_photonMultiplicity^{"+axisSuffix+"}",100,0,1.);
  book_TH1F("jet_puppiMultiplicity"+histSuffix,"jet_puppiMultiplicity^{"+axisSuffix+"}",100,0,1.);
  book_TH1F("jet_neutralPuppiMultiplicity"+histSuffix,"jet_neutralPuppiMultiplicity^{"+axisSuffix+"}",100,0,1.);
  book_TH1F("jet_neutralHadronPuppiMultiplicity"+histSuffix,"jet_neutralHadronPuppiMultiplicity^{"+axisSuffix+"}",100,0,1.);
  book_TH1F("jet_photonPuppiMultiplicity"+histSuffix,"jet_photonPuppiMultiplicity^{"+axisSuffix+"}",100,0,1.);
  book_TH1F("jet_HFHadronPuppiMultiplicity"+histSuffix,"jet_HFHadronPuppiMultiplicity^{"+axisSuffix+"}",100,0,1.);
  book_TH1F("jet_HFEMPuppiMultiplicity"+histSuffix,"jet_HFEMPuppiMultiplicity^{"+axisSuffix+"}",100,0,1.);
  book_TH1F("jetArea"+histSuffix,"jetArea^{"+axisSuffix+"}",150,0,15);
  book_TH1F("jetArea_pt200_300"+histSuffix,"jetArea^{"+axisSuffix+",pt(200,300)}",150,0,15);
  book_TH1F("jetArea_pt300_400"+histSuffix,"jetArea^{"+axisSuffix+",pt(300,400)}",150,0,15);
  book_TH2F("jetAreavspt"+histSuffix,";jetArea^{"+axisSuffix+"};p_{T} "+axisSuffix+" [GeV]",150,0,15,50,minPt,maxPt);
  book_TH2F("jetmassvspt"+histSuffix,";M^{"+axisSuffix+"};p_{T} "+axisSuffix+" [GeV]",100,0,300,50,minPt,maxPt);
  if (isTop) {
    book_TH1F("SDmass"+histSuffix,"SDmass^{"+axisSuffix+"} [GeV/c^{2}]",100,0,300);
    book_TH1F("SDmassvsinvMass"+histSuffix,"SDmass/invMass^{"+axisSuffix+"}",100,-1,2);
    book_TH2F("jetSDvspt"+histSuffix,";SDmass^{"+axisSuffix+"};p_{T} "+axisSuffix+" [GeV]",100,0,300,50,minPt,maxPt);
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
    book_TH1F("Match"+histSuffix,"Match^{"+axisSuffix+"}",19, 0, 19);
    book_TH1F("MatchingStatus"+histSuffix,"MatchingStatus^{"+axisSuffix+"}",10, 0, 10);
    book_TH2F("MatchvsMatchingStatus"+histSuffix,";Match^{"+axisSuffix+"};MatchingStatus^{"+axisSuffix+"}",19, 0, 19, 10, 0, 10);
    for (int i=1;i<20;i++) {
      H1("Match"+histSuffix)->GetXaxis()->SetBinLabel(i,MatchingToString(i-1).c_str());
      H2("MatchvsMatchingStatus"+histSuffix)->GetXaxis()->SetBinLabel(i,MatchingToString(i-1).c_str());
    }
    for (int i=1;i<11;i++) {
      H1("MatchingStatus"+histSuffix)->GetXaxis()->SetBinLabel(i,MatchingStatusToString(i-1).c_str());
      H2("MatchvsMatchingStatus"+histSuffix)->GetYaxis()->SetBinLabel(i,MatchingStatusToString(i-1).c_str());
    }

    for (std::string & disc : discriminators) {
      if (FindInString("tau", disc))       book_TH1F(disc+histSuffix,"#"+disc+"^{"+axisSuffix+"}", 30, -0.01, 1.01);
      else if (FindInString("btag", disc)) book_TH1F(disc+histSuffix, disc+"^{"+axisSuffix+"}", 30, -0.01, 1.01);
      else if (FindInString("chi2", disc)) book_TH1F(disc+histSuffix, disc+"^{"+axisSuffix+"}", 35, 0, 70);
      else book_TH1F(disc+histSuffix, disc+"^{"+axisSuffix+"}[GeV/c^{2}]", 40, 0., 200.);
    }

    for (std::string & disc : discriminators_subjets) {
      book_TH1F(disc+"_subjet"+histSuffix,   disc+"_{subjet}^{"+axisSuffix+"}",          30, -0.01, 1.01);
      book_TH1F(disc+"_subjet1"+histSuffix,  disc+"_{subjet1}^{"+axisSuffix+"}",         30, -0.01, 1.01);
      book_TH1F(disc+"_subjet2"+histSuffix,  disc+"_{subjet2}^{"+axisSuffix+"}",         30, -0.01, 1.01);
      book_TH1F(disc+"_subjet21"+histSuffix, disc+"_{subjet2/subjet1}^{"+axisSuffix+"}", 30, -0.01, 1.01);
      book_TH2F(disc+"_subjet12"+histSuffix, ";"+disc+"_{subjet1}^{"+axisSuffix+"};"+disc+"^{subjet2}^{"+axisSuffix+"}", 30, -0.01, 1.01, 30, -0.01, 1.01);
    }

    for (std::string & disc : discriminators_Extra) {
      bool isLong = FindInString("BoostedDoubleSecondary", disc);
      book_TH1F(disc+histSuffix, disc+"^{"+axisSuffix+"}", isLong? 60: 30, isLong? -1.01: -0.01, isLong? 1.01: 1.01);
    }

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
  fill_H2("etaphi"+histSuffix, jet.eta(), jet.phi(), weight);
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
  fill_H1("jetNEF"+histSuffix,jet.neutralEmEnergyFraction(),weight);
  fill_H1("jetNHF"+histSuffix,jet.neutralHadronEnergyFraction(),weight);
  fill_H1("jetCEF"+histSuffix,jet.chargedEmEnergyFraction(),weight);
  fill_H1("jetCHF"+histSuffix,jet.chargedHadronEnergyFraction(),weight);
  fill_H1("jet_muonEnergyFraction"+histSuffix,jet.muonEnergyFraction(),weight);
  fill_H1("jet_photonEnergyFraction"+histSuffix,jet.photonEnergyFraction(),weight);
  fill_H1("jet_chargedMultiplicity"+histSuffix,jet.chargedMultiplicity(),weight);
  fill_H1("jet_neutralMultiplicity"+histSuffix,jet.neutralMultiplicity(),weight);
  fill_H1("jet_muonMultiplicity"+histSuffix,jet.muonMultiplicity(),weight);
  fill_H1("jet_electronMultiplicity"+histSuffix,jet.electronMultiplicity(),weight);
  fill_H1("jet_photonMultiplicity"+histSuffix,jet.photonMultiplicity(),weight);
  fill_H1("jet_puppiMultiplicity"+histSuffix,jet.puppiMultiplicity(),weight);
  fill_H1("jet_neutralPuppiMultiplicity"+histSuffix,jet.neutralPuppiMultiplicity(),weight);
  fill_H1("jet_neutralHadronPuppiMultiplicity"+histSuffix,jet.neutralHadronPuppiMultiplicity(),weight);
  fill_H1("jet_photonPuppiMultiplicity"+histSuffix,jet.photonPuppiMultiplicity(),weight);
  fill_H1("jet_HFHadronPuppiMultiplicity"+histSuffix,jet.HFHadronPuppiMultiplicity(),weight);
  fill_H1("jet_HFEMPuppiMultiplicity"+histSuffix,jet.HFEMPuppiMultiplicity(),weight);
  fill_H1("jetArea"+histSuffix, jet.jetArea(), weight);
  if (jet.pt()>200 || jet.pt()<300) fill_H1("jetArea_pt200_300"+histSuffix,jet.jetArea(), weight);
  if (jet.pt()>300 || jet.pt()<400) fill_H1("jetArea_pt300_400"+histSuffix,jet.jetArea(), weight);
  fill_H2("jetAreavspt"+histSuffix, jet.jetArea(), jet.pt(), weight);
  fill_H2("jetmassvspt"+histSuffix, jet.v4().M(), jet.pt(), weight);
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
  fill_H2("jetSDvspt"+histSuffix, jet.softdropmass(), jet.pt(), weight);
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


  for (std::string & disc : discriminators) {
    double val=0;
    if (disc=="btag_MassDecorrelatedDeepBoosted_HccvsQCD") val = GetHccvsQCD(jet, true);
    if (disc=="btag_DeepBoosted_HccvsQCD") val = GetHccvsQCD(jet, false);
    if (disc=="btag_MassDecorrelatedDeepBoosted_ZHccvsQCD") val = GetZHccvsQCD(jet, true);
    if (disc=="btag_DeepBoosted_ZHccvsQCD") val = GetZHccvsQCD(jet, false);
    if (disc=="btag_MassDecorrelatedDeepBoosted_H4qvsQCD") val = jet.btag_MassDecorrelatedDeepBoosted_H4qvsQCD();
    if (disc=="btag_DeepBoosted_H4qvsQCD") val = jet.btag_DeepBoosted_H4qvsQCD();
    fill_H1(disc+histSuffix, val, weight);
  }

  int nsubjet = jet.subjets().size();
  for (std::string & disc : discriminators_subjets) {
    Jet subjet1, subjet2;
    if (nsubjet>0) subjet1 = jet.subjets().at(0);
    if (nsubjet>1) subjet2 = jet.subjets().at(1);
    double sub1=9999, sub2=0;
    if (disc=="btag_DeepJet") {                    sub1 = nsubjet>0 ? subjet1.btag_DeepJet() :                     9999; sub2 = nsubjet>1 ? subjet2.btag_DeepJet() : 0; }
    if (disc=="btag_DeepCSV") {                    sub1 = nsubjet>0 ? subjet1.btag_DeepCSV() :                     9999; sub2 = nsubjet>1 ? subjet2.btag_DeepCSV() : 0; }
    if (disc=="btag_DeepFlavour_bb") {             sub1 = nsubjet>0 ? subjet1.btag_DeepFlavour_bb() :              9999; sub2 = nsubjet>1 ? subjet2.btag_DeepFlavour_bb() : 0; }
    if (disc=="btag_DeepFlavour_b") {              sub1 = nsubjet>0 ? subjet1.btag_DeepFlavour_b() :               9999; sub2 = nsubjet>1 ? subjet2.btag_DeepFlavour_b() : 0; }
    if (disc=="btag_DeepFlavour_lepb") {           sub1 = nsubjet>0 ? subjet1.btag_DeepFlavour_lepb() :            9999; sub2 = nsubjet>1 ? subjet2.btag_DeepFlavour_lepb() : 0; }
    if (disc=="btag_DeepFlavour_uds") {            sub1 = nsubjet>0 ? subjet1.btag_DeepFlavour_uds() :             9999; sub2 = nsubjet>1 ? subjet2.btag_DeepFlavour_uds() : 0; }
    if (disc=="btag_DeepFlavour_g") {              sub1 = nsubjet>0 ? subjet1.btag_DeepFlavour_g() :               9999; sub2 = nsubjet>1 ? subjet2.btag_DeepFlavour_g() : 0; }
    if (disc=="btag_DeepFlavour_c") {              sub1 = nsubjet>0 ? subjet1.btag_DeepFlavour_c() :               9999; sub2 = nsubjet>1 ? subjet2.btag_DeepFlavour_c() : 0; }
    fill_H1(disc+"_subjet"+histSuffix, (sub1>sub2)?sub1:sub2, weight);
    fill_H1(disc+"_subjet1"+histSuffix, sub1, weight);
    fill_H1(disc+"_subjet2"+histSuffix, sub2, weight);
    fill_H1(disc+"_subjet21"+histSuffix, sub2/(sub1+sub2), weight);
    H2(disc+"_subjet12"+histSuffix)->Fill(sub1, sub2, weight);
  }

  for (std::string & disc : discriminators_Extra) {
    double val=0;
    if (disc=="btag_DeepBoosted_TvsQCD") val = jet.btag_DeepBoosted_TvsQCD();
    if (disc=="btag_DeepBoosted_WvsQCD") val = jet.btag_DeepBoosted_WvsQCD();
    if (disc=="btag_DeepBoosted_ZvsQCD") val = jet.btag_DeepBoosted_ZvsQCD();
    if (disc=="btag_DeepBoosted_HbbvsQCD") val = jet.btag_DeepBoosted_HbbvsQCD();
    if (disc=="btag_MassDecorrelatedDeepBoosted_TvsQCD") val = jet.btag_MassDecorrelatedDeepBoosted_TvsQCD();
    if (disc=="btag_MassDecorrelatedDeepBoosted_WvsQCD") val = jet.btag_MassDecorrelatedDeepBoosted_WvsQCD();
    if (disc=="btag_MassDecorrelatedDeepBoosted_ZHbbvsQCD") val = jet.btag_MassDecorrelatedDeepBoosted_ZHbbvsQCD();
    if (disc=="btag_MassDecorrelatedDeepBoosted_ZvsQCD") val = jet.btag_MassDecorrelatedDeepBoosted_ZvsQCD();
    if (disc=="btag_MassDecorrelatedDeepBoosted_ZbbvsQCD") val = jet.btag_MassDecorrelatedDeepBoosted_ZbbvsQCD();
    if (disc=="btag_MassDecorrelatedDeepBoosted_HbbvsQCD") val = jet.btag_MassDecorrelatedDeepBoosted_HbbvsQCD();
    if (disc=="btag_DeepBoosted_HbbvsHcc") val = jet.btag_DeepBoosted_probHbb()/(jet.btag_DeepBoosted_probHcc()+jet.btag_DeepBoosted_probHbb());
    if (disc=="btag_DeepBoosted_HvsQCD") val = (jet.btag_DeepBoosted_probHcc()+jet.btag_DeepBoosted_probHbb())/(jet.btag_DeepBoosted_probHcc()+jet.btag_DeepBoosted_probHbb()+GetQCD(jet, false));
    if (disc=="btag_DeepBoosted_ZbbvsQCD") val = jet.btag_DeepBoosted_ZbbvsQCD();
    if (disc=="btag_BoostedDoubleSecondaryVertexAK8") val = jet.btag_BoostedDoubleSecondaryVertexAK8();
    if (disc=="btag_BoostedDoubleSecondaryVertexCA15") val = jet.btag_BoostedDoubleSecondaryVertexCA15();
    if (disc=="btag_DeepDoubleBvLJet_probHbb") val = jet.btag_DeepDoubleBvLJet_probHbb();
    if (disc=="btag_DeepDoubleBvLJet_probQCD") val = jet.btag_DeepDoubleBvLJet_probQCD();
    if (disc=="btag_DeepDoubleCvBJet_probHbb") val = jet.btag_DeepDoubleCvBJet_probHbb();
    if (disc=="btag_DeepDoubleCvBJet_probHcc") val = jet.btag_DeepDoubleCvBJet_probHcc();
    if (disc=="btag_DeepDoubleCvLJet_probHcc") val = jet.btag_DeepDoubleCvLJet_probHcc();
    if (disc=="btag_DeepDoubleCvLJet_probQCD") val = jet.btag_DeepDoubleCvLJet_probQCD();
    if (disc=="btag_MassIndependentDeepDoubleBvLJet_probHbb") val = jet.btag_MassIndependentDeepDoubleBvLJet_probHbb();
    if (disc=="btag_MassIndependentDeepDoubleBvLJet_probQCD") val = jet.btag_MassIndependentDeepDoubleBvLJet_probQCD();
    if (disc=="btag_MassIndependentDeepDoubleCvBJet_probHbb") val = jet.btag_MassIndependentDeepDoubleCvBJet_probHbb();
    if (disc=="btag_MassIndependentDeepDoubleCvBJet_probHcc") val = jet.btag_MassIndependentDeepDoubleCvBJet_probHcc();
    if (disc=="btag_MassIndependentDeepDoubleCvLJet_probHcc") val = jet.btag_MassIndependentDeepDoubleCvLJet_probHcc();
    if (disc=="btag_MassIndependentDeepDoubleCvLJet_probQCD") val = jet.btag_MassIndependentDeepDoubleCvLJet_probQCD();
    if (disc=="btag_DeepBoosted_probQCDb") val = jet.btag_DeepBoosted_probQCDb();
    if (disc=="btag_DeepBoosted_probQCDbb") val = jet.btag_DeepBoosted_probQCDbb();
    if (disc=="btag_DeepBoosted_probQCDc") val = jet.btag_DeepBoosted_probQCDc();
    if (disc=="btag_DeepBoosted_probQCDcc") val = jet.btag_DeepBoosted_probQCDcc();
    if (disc=="btag_DeepBoosted_probQCDothers") val = jet.btag_DeepBoosted_probQCDothers();
    if (disc=="btag_DeepBoosted_probTbqq") val = jet.btag_DeepBoosted_probTbqq();
    if (disc=="btag_DeepBoosted_probTbcq") val = jet.btag_DeepBoosted_probTbcq();
    if (disc=="btag_DeepBoosted_probTbq") val = jet.btag_DeepBoosted_probTbq();
    if (disc=="btag_DeepBoosted_probTbc") val = jet.btag_DeepBoosted_probTbc();
    if (disc=="btag_DeepBoosted_probWqq") val = jet.btag_DeepBoosted_probWqq();
    if (disc=="btag_DeepBoosted_probWcq") val = jet.btag_DeepBoosted_probWcq();
    if (disc=="btag_DeepBoosted_probZcc") val = jet.btag_DeepBoosted_probZcc();
    if (disc=="btag_DeepBoosted_probZqq") val = jet.btag_DeepBoosted_probZqq();
    if (disc=="btag_DeepBoosted_probZbb") val = jet.btag_DeepBoosted_probZbb();
    if (disc=="btag_DeepBoosted_probHbb") val = jet.btag_DeepBoosted_probHbb();
    if (disc=="btag_DeepBoosted_probHcc") val = jet.btag_DeepBoosted_probHcc();
    if (disc=="btag_DeepBoosted_probHqqqq") val = jet.btag_DeepBoosted_probHqqqq();
    if (disc=="btag_DeepBoosted_raw_score_qcd") val = jet.btag_DeepBoosted_raw_score_qcd();
    if (disc=="btag_DeepBoosted_raw_score_top") val = jet.btag_DeepBoosted_raw_score_top();
    if (disc=="btag_DeepBoosted_raw_score_w") val = jet.btag_DeepBoosted_raw_score_w();
    if (disc=="btag_DeepBoosted_raw_score_z") val = jet.btag_DeepBoosted_raw_score_z();
    if (disc=="btag_DeepBoosted_raw_score_h") val = jet.btag_DeepBoosted_raw_score_h();
    if (disc=="btag_MassDecorrelatedDeepBoosted_bbvsLight") val = jet.btag_MassDecorrelatedDeepBoosted_bbvsLight();
    if (disc=="btag_MassDecorrelatedDeepBoosted_ccvsLight") val = jet.btag_MassDecorrelatedDeepBoosted_ccvsLight();
    if (disc=="btag_MassDecorrelatedDeepBoosted_probHbb") val = jet.btag_MassDecorrelatedDeepBoosted_probHbb();
    if (disc=="btag_MassDecorrelatedDeepBoosted_probQCDc") val = jet.btag_MassDecorrelatedDeepBoosted_probQCDc();
    if (disc=="btag_MassDecorrelatedDeepBoosted_probQCDbb") val = jet.btag_MassDecorrelatedDeepBoosted_probQCDbb();
    if (disc=="btag_MassDecorrelatedDeepBoosted_probTbqq") val = jet.btag_MassDecorrelatedDeepBoosted_probTbqq();
    if (disc=="btag_MassDecorrelatedDeepBoosted_probTbcq") val = jet.btag_MassDecorrelatedDeepBoosted_probTbcq();
    if (disc=="btag_MassDecorrelatedDeepBoosted_probTbq") val = jet.btag_MassDecorrelatedDeepBoosted_probTbq();
    if (disc=="btag_MassDecorrelatedDeepBoosted_probQCDothers") val = jet.btag_MassDecorrelatedDeepBoosted_probQCDothers();
    if (disc=="btag_MassDecorrelatedDeepBoosted_probQCDb") val = jet.btag_MassDecorrelatedDeepBoosted_probQCDb();
    if (disc=="btag_MassDecorrelatedDeepBoosted_probTbc") val = jet.btag_MassDecorrelatedDeepBoosted_probTbc();
    if (disc=="btag_MassDecorrelatedDeepBoosted_probWqq") val = jet.btag_MassDecorrelatedDeepBoosted_probWqq();
    if (disc=="btag_MassDecorrelatedDeepBoosted_probQCDcc") val = jet.btag_MassDecorrelatedDeepBoosted_probQCDcc();
    if (disc=="btag_MassDecorrelatedDeepBoosted_probHcc") val = jet.btag_MassDecorrelatedDeepBoosted_probHcc();
    if (disc=="btag_MassDecorrelatedDeepBoosted_probZcc") val = jet.btag_MassDecorrelatedDeepBoosted_probZcc();
    if (disc=="btag_MassDecorrelatedDeepBoosted_proWcq") val = jet.btag_MassDecorrelatedDeepBoosted_proWcq();
    if (disc=="btag_MassDecorrelatedDeepBoosted_probZqq") val = jet.btag_MassDecorrelatedDeepBoosted_probZqq();
    if (disc=="btag_MassDecorrelatedDeepBoosted_probHqqqq") val = jet.btag_MassDecorrelatedDeepBoosted_probHqqqq();
    if (disc=="btag_MassDecorrelatedDeepBoosted_probZbb") val = jet.btag_MassDecorrelatedDeepBoosted_probZbb();

    fill_H1(disc+histSuffix, val, weight);
  }

  #define MYTAGFILL(mytag)\
  fill_H1(disc+histSuffix, jet.mytag(),weight);\

  return;
}
