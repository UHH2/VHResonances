#include "UHH2/core/include/Event.h"
#include "UHH2/VHResonances/include/DiLeptonHists.h"

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace std;
using namespace uhh2;

DiLeptonHists::DiLeptonHists(Context & ctx, const string & dname, const string & jetcollection,const string & topjetcollection): HistsBase(ctx, dname){

  if (jetcollection != "") h_jets = ctx.get_handle<vector<Jet>>(jetcollection);
  if (topjetcollection != "") h_topjets = ctx.get_handle<vector<TopJet>>(topjetcollection);

  // book all histograms here
  for (const string & lepton: {"Electron", "Muon"}) {
    string dilepton = (lepton=="Electron") ? "ee" : "#mu#mu";
    book_TH1F("di"+lepton+"_number",  "number of "+dilepton,                11, -.5, 10.5);
    book_TH1F("di"+lepton+"_charge",  "charge of "+dilepton,                5,-2.5,2.5);
    book_TH1F("di"+lepton+"_mass",    "m^"        +dilepton+" [GeV/c^{2}]", 40,70,110);
    book_TH1F("di"+lepton+"_pt",      "p_{T}^"    +dilepton+" [GeV]",       30,0,1500);
    book_TH1F("di"+lepton+"_eta",     "#eta "     +dilepton,                100,-5,5);
    book_TH1F("di"+lepton+"_phi",     "#phi "     +dilepton,                50,-M_PI,M_PI);
    book_TH1F("di"+lepton+"_DR12",    "#Delta R(" +dilepton+")",            100, 0, M_PI);
    book_TH1F("di"+lepton+"_jet_Dphi","#Delta#phi("+dilepton+",jet)",       100, 0, M_PI);
    book_TH1F("di"+lepton+"_jet_Deta","#Delta#eta("+dilepton+",jet)",       100, 0, 3*M_PI);
    book_TH1F("di"+lepton+"_jet_DR",  "#Delta R("+dilepton+",jet)",         100, 0, 3*M_PI);
    book_TH2F("di"+lepton+"_pt1_pt2", ";PT1;pT2",                           100, 0, 500, 100, 0, 500.);
  }

}


void DiLeptonHists::fill(const Event & event){

  vector<Jet> jets;
  if (event.is_valid(h_topjets)) jets.assign((&event.get(h_topjets))->begin(), (&event.get(h_topjets))->end());
  else if (event.is_valid(h_jets)) jets.assign((&event.get(h_jets))->begin(), (&event.get(h_jets))->end());
  else throw logic_error("DiLeptonHists: Impossible case");

  // fill the histograms.
  auto weight = event.weight;

  for (const string & lepton: {"Electron", "Muon"}) {
    vector<Particle> leptons;
    if (lepton == "Muon") leptons.assign(event.muons->begin(), event.muons->end());
    else if (lepton == "Electron") leptons.assign(event.electrons->begin(), event.electrons->end());
    else throw logic_error("DiLeptonHists: Impossible case");

    fill_H1("di"+lepton+"_number",  leptons.size(),weight);
    if (leptons.size() < min_leptons) continue;
    const auto& lep1 = leptons.at(0);
    const auto& lep2 = leptons.at(1);
    auto dilep = lep1.v4() + lep2.v4();
    fill_H1("di"+lepton+"_charge",  lep1.charge()+lep2.charge(),weight);
    fill_H1("di"+lepton+"_mass",    dilep.M(),weight);
    fill_H1("di"+lepton+"_pt",      dilep.Pt(),weight);
    fill_H1("di"+lepton+"_eta",     dilep.Eta(),weight);
    fill_H1("di"+lepton+"_phi",     dilep.Phi(),weight);
    fill_H1("di"+lepton+"_DR12",    uhh2::deltaR(lep1, lep2),weight);
    if (jets.size()>0) fill_H1("di"+lepton+"_jet_Dphi", deltaPhi(dilep, jets.at(0)), weight);
    if (jets.size()>0) fill_H1("di"+lepton+"_jet_Deta", dilep.Eta() - jets.at(0).eta(), weight);
    if (jets.size()>0) fill_H1("di"+lepton+"_jet_DR", uhh2::deltaR(dilep, jets.at(0)), weight);
    fill_H2("di"+lepton+"_pt1_pt2", lep1.pt(), lep2.pt(),weight);
  }

}

DiLeptonHists::~DiLeptonHists(){}
