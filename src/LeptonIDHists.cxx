#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/VHResonances/include/LeptonIDHists.h"
#include "UHH2/VHResonances/include/Utils.hpp"

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace std;
using namespace uhh2;

LeptonIDHists::LeptonIDHists(Context & ctx, const string & dname): HistsBase(ctx, dname) {


  muonchannel     = string2bool(ctx.get("muonchannel"));
  electronchannel = string2bool(ctx.get("electronchannel"));

  if (muonchannel && electronchannel) throw std::runtime_error("In LeptonIDHists.cxx: Choose exactly one lepton channel.");

  if (electronchannel) {
    eleIds["ele_ID_kincut"]                   = PtEtaSCCut(min_lepton_pt, min_lepton_eta);
    eleIds["ele_ID_MVA_loose_iso"]            = AndId<Electron>(ElectronID_MVA_Fall17_loose_iso,    PtEtaSCCut(min_lepton_pt, min_lepton_eta));
    eleIds["ele_ID_MVA_loose_noIso"]          = AndId<Electron>(ElectronID_MVA_Fall17_loose_noIso,  PtEtaSCCut(min_lepton_pt, min_lepton_eta));
    eleIds["ele_ID_HEEP"]                     = AndId<Electron>(ElectronID_HEEP_RunII_25ns,         PtEtaSCCut(min_lepton_pt, min_lepton_eta));
    eleIds["ele_ID_veto"]                     = AndId<Electron>(ElectronID_Fall17_veto,             PtEtaSCCut(min_lepton_pt, min_lepton_eta));
    eleIds["ele_ID_veto_noIso"]               = AndId<Electron>(ElectronID_Fall17_veto_noIso,       PtEtaSCCut(min_lepton_pt, min_lepton_eta));
    eleIds["ele_ID_loose"]                    = AndId<Electron>(ElectronID_Fall17_loose,            PtEtaSCCut(min_lepton_pt, min_lepton_eta));
    eleIds["ele_ID_loose_noIso"]              = AndId<Electron>(ElectronID_Fall17_loose_noIso,      PtEtaSCCut(min_lepton_pt, min_lepton_eta));
    eleIds["ele_ID_medium"]                   = AndId<Electron>(ElectronID_Fall17_medium,           PtEtaSCCut(min_lepton_pt, min_lepton_eta));
    eleIds["ele_ID_medium_noIso"]             = AndId<Electron>(ElectronID_Fall17_medium_noIso,     PtEtaSCCut(min_lepton_pt, min_lepton_eta));
    eleIds["ele_ID_tight"]                    = AndId<Electron>(ElectronID_Fall17_tight,            PtEtaSCCut(min_lepton_pt, min_lepton_eta));
    eleIds["ele_ID_tight_noIso"]              = AndId<Electron>(ElectronID_Fall17_tight_noIso,      PtEtaSCCut(min_lepton_pt, min_lepton_eta));
    for (const auto& id :eleIds) leptoncollections.push_back(id.first);
  }

  if (muonchannel) {
    muoIds["muon_ID_kincut"]                  = PtEtaCut(min_lepton_pt, min_lepton_eta);
    muoIds["muon_ID_CutBasedIdLoose"]         = AndId<Muon>(MuonID(Muon::CutBasedIdLoose),          PtEtaCut(min_lepton_pt, min_lepton_eta));
    muoIds["muon_ID_CutBasedIdMedium"]        = AndId<Muon>(MuonID(Muon::CutBasedIdMedium),         PtEtaCut(min_lepton_pt, min_lepton_eta));
    muoIds["muon_ID_CutBasedIdMediumPrompt"]  = AndId<Muon>(MuonID(Muon::CutBasedIdMediumPrompt),   PtEtaCut(min_lepton_pt, min_lepton_eta));
    muoIds["muon_ID_CutBasedIdTight"]         = AndId<Muon>(MuonID(Muon::CutBasedIdTight),          PtEtaCut(min_lepton_pt, min_lepton_eta));
    muoIds["muon_ID_CutBasedIdGlobalHighPt"]  = AndId<Muon>(MuonID(Muon::CutBasedIdGlobalHighPt),   PtEtaCut(min_lepton_pt, min_lepton_eta));
    muoIds["muon_ID_CutBasedIdTrkHighPt"]     = AndId<Muon>(MuonID(Muon::CutBasedIdTrkHighPt),      PtEtaCut(min_lepton_pt, min_lepton_eta));
    muoIds["muon_ID_MvaLoose"]                = AndId<Muon>(MuonID(Muon::MvaLoose),                 PtEtaCut(min_lepton_pt, min_lepton_eta));
    muoIds["muon_ID_MvaMedium"]               = AndId<Muon>(MuonID(Muon::MvaMedium),                PtEtaCut(min_lepton_pt, min_lepton_eta));
    muoIds["muon_ID_MvaTight"]                = AndId<Muon>(MuonID(Muon::MvaTight),                 PtEtaCut(min_lepton_pt, min_lepton_eta));
    muoIds["muon_ID_TkIsoLoose"]              = AndId<Muon>(MuonID(Muon::TkIsoLoose),               PtEtaCut(min_lepton_pt, min_lepton_eta));
    muoIds["muon_ID_TkIsoTight"]              = AndId<Muon>(MuonID(Muon::TkIsoTight),               PtEtaCut(min_lepton_pt, min_lepton_eta));
    for (const auto& id :muoIds) leptoncollections.push_back(id.first);
  }

  string name = muonchannel? "#mu" : "e";
  for (const std::string& lepcoll : leptoncollections) {

    if (muonchannel)     { h_muon[lepcoll] = ctx.get_handle<std::vector<Muon>>(lepcoll); }
    if (electronchannel) { h_ele[lepcoll]  = ctx.get_handle<std::vector<Electron>>(lepcoll); }

    for (const std::string& cut : {"", "_sel", "_match", "_2IDs", "_sel_2IDs", "_match_2IDs"}) {
      book_TH1F(lepcoll+"_pt12"+cut, "p_{T}^{"+name+"_1,"+name+"_2} [GeV]",   50, 0.0, 5000);
      book_TH1F(lepcoll+"_DR12"+cut, "#Delta R("+name+"_1,"+name+"_2)",      100, 0.0, 5);
    }
  }
}


void LeptonIDHists::fill(const Event & event){

  // fill the histograms.
  auto weight = event.weight;

  int gen_ind1 =-1, gen_ind2 = -1;

  for( auto & gp : *event.genparticles) {
    if (abs(gp.pdgId())!=ParticleID::ZPrime) continue;

    const GenParticle* d1 = gp.daughter(event.genparticles, 1);
    const GenParticle* d2 = gp.daughter(event.genparticles, 2);
    int d1_ID = abs(d1->pdgId());
    int d2_ID = abs(d2->pdgId());
    if (!DoubleDecay(d1_ID, d2_ID, Decay::ZH)) continue;
    int index = (d1_ID==ParticleID::Z)? 1: 2;
    const GenParticle* sd1 = gp.daughter(event.genparticles, index)->daughter(event.genparticles, 1);
    const GenParticle* sd2 = gp.daughter(event.genparticles, index)->daughter(event.genparticles, 2);

    int sd1_ID = sd1->pdgId();
    int sd2_ID = sd2->pdgId();
    if (!DoubleDecay(abs(sd1_ID), abs(sd2_ID), (muonchannel? Decay::mumu : Decay::ee))) continue;
    gen_ind1 = sd1->index();
    gen_ind2 = sd2->index();
  }

  if (gen_ind1 == -1 || gen_ind2 == -1) return;

  GenParticle gen_lep1 = event.genparticles->at(gen_ind1);
  GenParticle gen_lep2 = event.genparticles->at(gen_ind2);

  vector<Particle> leptons;
  if (muonchannel) leptons.assign(event.muons->begin(), event.muons->end());
  else if (electronchannel) leptons.assign(event.electrons->begin(), event.electrons->end());
  else throw logic_error("LeptonIDHists: Impossible case");

  int lepton_size = leptons.size();
  if (lepton_size < 2) return;

  for (unsigned int i = 0; i < leptons.size(); i++) {
    Particle lep1 = leptons.at(i);
    for (unsigned int j = i+1; j < leptons.size(); j++) {
      Particle lep2 = leptons.at(j);
      double deltaR11 = uhh2::deltaR(lep1, gen_lep1);
      double deltaR22 = uhh2::deltaR(lep2, gen_lep2);
      double deltaR12 = uhh2::deltaR(lep1, gen_lep2);
      double deltaR21 = uhh2::deltaR(lep2, gen_lep1);
      bool gen_match = (deltaR11<0.1 && deltaR22<0.1) || (deltaR12<0.1 && deltaR21<0.1);
      auto DR = uhh2::deltaR(lep1, lep2);
      auto diLep = lep1.v4() + lep2.v4();
      double diLepPt = diLep.Pt();
      bool sel = (fabs(diLep.M() - ZMASS) < ZWIDTH) && (lep1.charge() + lep2.charge() == 0) && diLep.pt()>min_dilep_pt && min_DR_dilep < DR && DR < max_DR_dilep;
      bool goodID1=false, goodID2=false, TrkhighID1 = false, TrkhighID2=false;
      for (const std::string& lepcoll : leptoncollections) {
        if (muonchannel) {
          goodID1 = muoIds[lepcoll](event.muons->at(i), event);
          goodID2 = muoIds[lepcoll](event.muons->at(j), event);
          TrkhighID1 = muoIds["muon_ID_CutBasedIdTrkHighPt"](event.muons->at(i), event);
          TrkhighID2 = muoIds["muon_ID_CutBasedIdTrkHighPt"](event.muons->at(j), event);
        }
        if (electronchannel) {
          goodID1 = eleIds[lepcoll](event.electrons->at(i), event);
          goodID2 = eleIds[lepcoll](event.electrons->at(j), event);
        }
        if (!(goodID1 || goodID2)) continue;
        if (goodID1 && goodID2) {
          fill_H1(lepcoll+"_pt12", diLepPt, weight);
          fill_H1(lepcoll+"_DR12", DR, weight);
          if (sel) {
            fill_H1(lepcoll+"_pt12_sel", diLepPt, weight);
            fill_H1(lepcoll+"_DR12_sel", DR, weight);
          }
          if (gen_match) {
            fill_H1(lepcoll+"_pt12_match", diLepPt, weight);
            fill_H1(lepcoll+"_DR12_match", DR, weight);
          }
        }

        if ((TrkhighID1 && TrkhighID2) && (goodID1 || goodID2)) {
          fill_H1(lepcoll+"_pt12_2IDs", diLepPt, weight);
          fill_H1(lepcoll+"_DR12_2IDs", DR, weight);
          if (sel) {
            fill_H1(lepcoll+"_pt12_sel_2IDs", diLepPt, weight);
            fill_H1(lepcoll+"_DR12_sel_2IDs", DR, weight);
          }
          if (gen_match) {
            fill_H1(lepcoll+"_pt12_match_2IDs", diLepPt, weight);
            fill_H1(lepcoll+"_DR12_match_2IDs", DR, weight);
          }
        }

      }
    }
  }
}


LeptonIDHists::~LeptonIDHists(){}
