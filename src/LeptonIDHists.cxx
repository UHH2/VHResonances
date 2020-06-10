#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/VHResonances/include/LeptonIDHists.h"

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace std;
using namespace uhh2;

LeptonIDHists::LeptonIDHists(Context & ctx, const string & dname): HistsBase(ctx, dname) {

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

  for (const auto& id :eleIds) leptoncollections.push_back(id.first);
  for (const auto& id :muoIds) leptoncollections.push_back(id.first);

  for (const std::string& lepcoll : leptoncollections) {

    bool isMuon = lepcoll.find("muon") != std::string::npos;
    bool isEle  = lepcoll.find("ele") != std::string::npos;
    if (!isMuon && !isEle) throw logic_error("LeptonIDHists: Impossible case");

    string name;
    if (isMuon) { name = "#mu"; h_muon[lepcoll] = ctx.get_handle<std::vector<Muon>>(lepcoll); }
    if (isEle)  { name = "e";   h_ele[lepcoll]  = ctx.get_handle<std::vector<Electron>>(lepcoll); }

    book_TH1F(lepcoll+"_number",        "number of "+name,                      11,   -0.5, 10.5);
    book_TH1F(lepcoll+"_pt12",          "p_{T}^{"+name+"_1,"+name+"_2} [GeV]",  50,    0.0, 5000);
    book_TH1F(lepcoll+"_DR12",          "#Delta R("+name+"_1,"+name+"_2)",      100,   0.0, 4);
    book_TH1F(lepcoll+"_pt12_sel",      "p_{T}^{"+name+"_1,"+name+"_2} [GeV]",  50,    0.0, 5000);
    book_TH1F(lepcoll+"_DR12_sel",      "#Delta R("+name+"_1,"+name+"_2)",      100,   0.0, 4);
    book_TH1F(lepcoll+"_pt12_2IDs",     "p_{T}^{"+name+"_1,"+name+"_2} [GeV]",  50,    0.0, 5000);
    book_TH1F(lepcoll+"_DR12_2IDs",     "#Delta R("+name+"_1,"+name+"_2)",      100,   0.0, 4);
    book_TH1F(lepcoll+"_pt12_sel_2IDs", "p_{T}^{"+name+"_1,"+name+"_2} [GeV]",  50,    0.0, 5000);
    book_TH1F(lepcoll+"_DR12_sel_2IDs", "#Delta R("+name+"_1,"+name+"_2)",      100,   0.0, 4);
  }
}


void LeptonIDHists::fill(const Event & event){

  // fill the histograms.
  auto weight = event.weight;

  for (const std::string& lepcoll : leptoncollections) {

    bool isMuon = lepcoll.find("muon") != std::string::npos;
    bool isEle  = lepcoll.find("ele") != std::string::npos;
    if (!isMuon && !isEle) throw logic_error("LeptonIDHists: Impossible case");

    vector<Particle> leptons;
    if (isMuon) leptons.assign(event.muons->begin(), event.muons->end());
    else if (isEle) leptons.assign(event.electrons->begin(), event.electrons->end());
    else throw logic_error("ZprimeCandidateReconstruction: Impossible case");

    int lepton_size = leptons.size();
    fill_H1(lepcoll+"_number",  lepton_size, weight);
    if (lepton_size < 2) continue;
    for (unsigned int i = 0; i < leptons.size(); i++) {
      Particle lep1 = leptons.at(i);
      for (unsigned int j = i+1; j < leptons.size(); j++) {
        Particle lep2 = leptons.at(j);
        bool goodID1=false, goodID2=false, TrkhighID1 = false, TrkhighID2=false;
        if (isMuon) {
          goodID1 = muoIds[lepcoll](event.muons->at(i), event);
          goodID2 = muoIds[lepcoll](event.muons->at(j), event);
          TrkhighID1 = muoIds["muon_ID_CutBasedIdTrkHighPt"](event.muons->at(i), event);
          TrkhighID2 = muoIds["muon_ID_CutBasedIdTrkHighPt"](event.muons->at(j), event);
        }
        if (isEle) {
          goodID1 = eleIds[lepcoll](event.electrons->at(i), event);
          goodID2 = eleIds[lepcoll](event.electrons->at(j), event);
        }
        if (!(goodID1 || goodID2)) continue;
        auto DR = uhh2::deltaR(lep1, lep2);
        auto diLep = lep1.v4() + lep2.v4();
        if (goodID1 && goodID2){
          if (i==0 && j==1) {
            fill_H1(lepcoll+"_pt12", diLep.Pt(), weight);
            fill_H1(lepcoll+"_DR12", DR, weight);
          }
          if( (fabs(diLep.M() - ZMASS) < ZWIDTH) && (lep1.charge() + lep2.charge() == 0) && diLep.pt()>min_dilep_pt && min_DR_dilep < DR && DR < max_DR_dilep) {
            fill_H1(lepcoll+"_pt12_sel", diLep.Pt(), weight);
            fill_H1(lepcoll+"_DR12_sel", uhh2::deltaR(lep1, lep2), weight);
          }
        }
        if ((goodID1 && TrkhighID2) || (goodID2 && TrkhighID1)) {
          if (i==0 && j==1) {
            fill_H1(lepcoll+"_pt12_2IDs", diLep.Pt(), weight);
            fill_H1(lepcoll+"_DR12_2IDs", DR, weight);
          }
          if( (fabs(diLep.M() - ZMASS) < ZWIDTH) && (lep1.charge() + lep2.charge() == 0) && diLep.pt()>min_dilep_pt && min_DR_dilep < DR && DR < max_DR_dilep) {
            fill_H1(lepcoll+"_pt12_sel_2IDs", diLep.Pt(), weight);
            fill_H1(lepcoll+"_DR12_sel_2IDs", uhh2::deltaR(lep1, lep2), weight);
          }
        }
      }
    }
  }
}

LeptonIDHists::~LeptonIDHists(){}
