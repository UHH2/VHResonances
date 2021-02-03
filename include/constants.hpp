#pragma once

#include <string>
#include <unordered_map>

#include "TMath.h"

#include "UHH2/common/include/JetIds.h"

// DEFINE GLOBAL VARIABLES

const float ZMASS  = 91.0;
const float ZWIDTH = 10.0;

const float HMASS  = 120.0;
const float HWIDTH = 15.0;

const int min_leptons = 2;
const float min_jet_pt = 30.0;
const float min_topjet_pt = 200.0;
// const float min_lepton_pt = 30.0;
const float min_lepton_pt = 52.0;
const float min_lepton_eta = 2.4;
// const float min_jet_dilep_delta_phi = 2.7;
const float min_MET_pt = 250.0;

const float min_dilep_pt = 200;
const float min_DR_dilep = 0.0;
const float max_DR_dilep = 1.0;
const float min_jet_dilep_delta_phi = 2.0;
const float max_jet_dilep_delta_phi = M_PI;
const float min_Dphi_AK8jet_MET = 2.0;
const float min_Dphi_AK4jet_MET = 0.5;

const float min_Z_pt_ZH_mass = 0.2;
const float min_Z_pt_ZH_mass_invisible = 0.4;
const float min_ZH_mass = 700;

const BTag::algo BTag_algo = BTag::DEEPCSV;
const BTag::wp BTag_wp = BTag::WP_LOOSE;

const double TaggerThr = 0.8;


// Theory Corrections
const bool do_EWK = true;
const bool do_QCD_EWK = false;
const bool do_QCD_NLO  = true;
const bool do_QCD_NNLO = false;

// Taken from :
// https://twiki.cern.ch/twiki/bin/view/CMS/EgHLTRunIISummary
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/B2GMuon
// https://twiki.cern.ch/twiki/bin/view/CMS/MuonHLT2017

const std::unordered_map<std::string, std::map<std::string, std::pair<int, int>>>
Trigger_run_validity = {
  { "2016", {
    { "HLT_Mu50_v*",                              std::pair(272760, 284044) }, // 272007
    { "HLT_TkMu50_v*",                            std::pair(274954, 284044) },
    { "HLT_Ele27_WPTight_Gsf_v*",                 std::pair(272760, 284044) }, // 273158
    { "HLT_Ele115_CaloIdVT_GsfTrkIdT_v*",         std::pair(272760, 284044) }, // 273158, 284044
    { "HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v*", std::pair(274954, 284044) },                    // most restrictive numbers
    { "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v*", std::pair(274954, 284044) },                    // most restrictive numbers
    { "HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight_v*", std::pair(274954, 284044) }, // most restrictive numbers
    { "HLT_PFMET110_PFMHT110_IDTight_v*",         std::pair(274954, 284044) },                    // most restrictive numbers
    { "HLT_PFMET120_PFMHT120_IDTight_v*",         std::pair(274954, 284044) },                    // most restrictive numbers
    { "HLT_PFMET170_NotCleaned_v*",               std::pair(274954, 284044) },                    // most restrictive numbers
    { "HLT_PFMET170_HBHECleaned*",                std::pair(274954, 284044) },                    // most restrictive numbers
    { "HLT_Photon175_v*",                         std::pair(272760, 284044) }, // 273158, 284044
    { "HLT_AK8PFJet40_v*",                        std::pair(272760, 284044) },
    { "HLT_AK8PFJet60_v*",                        std::pair(272760, 284044) },
    { "HLT_AK8PFJet80_v*",                        std::pair(272760, 284044) },
    { "HLT_AK8PFJet140_v*",                       std::pair(272760, 284044) },
    { "HLT_AK8PFJet200_v*",                       std::pair(272760, 284044) },
    { "HLT_AK8PFJet260_v*",                       std::pair(272760, 284044) },
    { "HLT_AK8PFJet320_v*",                       std::pair(272760, 284044) },
    { "HLT_AK8PFJet400_v*",                       std::pair(272760, 284044) },
    { "HLT_AK8PFJet450_v*",                       std::pair(272760, 284044) },
    { "HLT_AK8PFJet500_v*",                       std::pair(272760, 284044) },
  }},
  { "2017", {
    { "HLT_Mu50_v*",                              std::pair(296070, 306460) }, // 297020, 306462
    { "HLT_OldMu100_v*",                          std::pair(299368, 306460) }, //306462 // 299370
    { "HLT_TkMu100_v*",                           std::pair(299368, 306460) }, // 306462
    { "HLT_Ele35_WPTight_Gsf_v*",                 std::pair(296070, 306460)}, // 297050, 306460
    { "HLT_Ele115_CaloIdVT_GsfTrkIdT_v*",         std::pair(299368, 306460)},
    { "HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v*", std::pair(299368, 306460) },                    // most restrictive numbers
    { "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v*", std::pair(299368, 306460) },                    // most restrictive numbers
    { "HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v*", std::pair(299368, 306460) },                    // most restrictive numbers
    { "HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v*", std::pair(299368, 306460) },                    // most restrictive numbers
    { "HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight_v*", std::pair(299368, 306460) }, // most restrictive numbers
    { "HLT_PFMET110_PFMHT110_IDTight_v*",         std::pair(299368, 306460) },                    // most restrictive numbers
    { "HLT_PFMET120_PFMHT120_IDTight_v*",         std::pair(299368, 306460) },                    // most restrictive numbers
    { "HLT_PFMET130_PFMHT130_IDTight_v*",         std::pair(299368, 306460) },                    // most restrictive numbers
    { "HLT_PFMET140_PFMHT140_IDTight_v*",         std::pair(299368, 306460) },                    // most restrictive numbers
    { "HLT_PFMETTypeOne110_PFMHT110_IDTight_v*",  std::pair(299368, 306460) },                    // most restrictive numbers
    { "HLT_PFMETTypeOne120_PFMHT120_IDTight_v*",  std::pair(299368, 306460) },                    // most restrictive numbers
    { "HLT_PFMETTypeOne130_PFMHT130_IDTight_v*",  std::pair(299368, 306460) },                    // most restrictive numbers
    { "HLT_PFMETTypeOne140_PFMHT140_IDTight_v*",  std::pair(299368, 306460) },                    // most restrictive numbers
    { "HLT_PFMET200_NotCleaned_v*",               std::pair(299368, 306460) },                    // most restrictive numbers
    { "HLT_PFMET200_HBHECleaned_v*",              std::pair(299368, 306460) },                    // most restrictive numbers
    { "HLT_PFMET200_HBHE_BeamHaloCleaned_v*",     std::pair(299368, 306460) },                    // most restrictive numbers
    { "HLT_PFMET250_HBHECleaned*",                std::pair(299368, 306460) },                    // most restrictive numbers
    { "HLT_Photon200_v*",                         std::pair(299368, 306460) }, // 273158, 284044
    { "HLT_AK8PFJet40_v*",                        std::pair(296070, 306460) },
    { "HLT_AK8PFJet60_v*",                        std::pair(296070, 306460) },
    { "HLT_AK8PFJet80_v*",                        std::pair(296070, 306460) },
    { "HLT_AK8PFJet140_v*",                       std::pair(296070, 306460) },
    { "HLT_AK8PFJet200_v*",                       std::pair(296070, 306460) },
    { "HLT_AK8PFJet260_v*",                       std::pair(296070, 306460) },
    { "HLT_AK8PFJet320_v*",                       std::pair(296070, 306460) },
    { "HLT_AK8PFJet400_v*",                       std::pair(296070, 306460) },
    { "HLT_AK8PFJet450_v*",                       std::pair(296070, 306460) },
    { "HLT_AK8PFJet500_v*",                       std::pair(296070, 306460) },
  }},
  { "2018", {
    { "HLT_Mu50_v*",                              std::pair(315252, 325175) }, // 315252, 325273
    { "HLT_OldMu100_v*",                          std::pair(315252, 325175) }, // 315252, 325273
    { "HLT_TkMu100_v*",                           std::pair(315252, 325175) }, // 315252, 325273
    { "HLT_Ele32_WPTight_Gsf_v*",                 std::pair(315252, 325175) }, // 315257, 325172
    { "HLT_Ele115_CaloIdVT_GsfTrkIdT_v*",         std::pair(315252, 325175)}, // 315257,325172
    { "HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v*", std::pair(315252, 325175) },                    // most restrictive numbers
    { "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v*", std::pair(315252, 325175) },                    // most restrictive numbers
    { "HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v*", std::pair(315252, 325175) },                    // most restrictive numbers
    { "HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v*", std::pair(315252, 325175) },                    // most restrictive numbers
    { "HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight_v*", std::pair(315252, 325175) }, // most restrictive numbers
    { "HLT_PFMET110_PFMHT110_IDTight_v*",         std::pair(315252, 325175) },                    // most restrictive numbers
    { "HLT_PFMET120_PFMHT120_IDTight_v*",         std::pair(315252, 325175) },                    // most restrictive numbers
    { "HLT_PFMET130_PFMHT130_IDTight_v*",         std::pair(315252, 325175) },                    // most restrictive numbers
    { "HLT_PFMET140_PFMHT140_IDTight_v*",         std::pair(315252, 325175) },                    // most restrictive numbers
    { "HLT_PFMETTypeOne110_PFMHT110_IDTight_v*",  std::pair(315252, 325175) },                    // most restrictive numbers
    { "HLT_PFMETTypeOne120_PFMHT120_IDTight_v*",  std::pair(315252, 325175) },                    // most restrictive numbers
    { "HLT_PFMETTypeOne130_PFMHT130_IDTight_v*",  std::pair(315252, 325175) },                    // most restrictive numbers
    { "HLT_PFMETTypeOne140_PFMHT140_IDTight_v*",  std::pair(315252, 325175) },                    // most restrictive numbers
    { "HLT_PFMET200_NotCleaned_v*",               std::pair(315252, 325175) },                    // most restrictive numbers
    { "HLT_PFMET200_HBHECleaned_v*",              std::pair(315252, 325175) },                    // most restrictive numbers
    { "HLT_PFMET200_HBHE_BeamHaloCleaned_v*",     std::pair(315252, 325175) },                    // most restrictive numbers
    { "HLT_PFMET250_HBHECleaned*",                std::pair(315252, 325175) },                    // most restrictive numbers
    { "HLT_Photon200_v*",                         std::pair(315252, 325175) }, // 273158, 325175
    { "HLT_PFHT180_v*",                           std::pair(315252, 325175) },
    { "HLT_PFHT250_v*",                           std::pair(315252, 325175) },
    { "HLT_PFHT350_v*",                           std::pair(315252, 325175) },
    { "HLT_PFHT370_v*",                           std::pair(315252, 325175) },
    { "HLT_PFHT430_v*",                           std::pair(315252, 325175) },
    { "HLT_PFHT510_v*",                           std::pair(315252, 325175) },
    { "HLT_PFHT590_v*",                           std::pair(315252, 325175) },
    { "HLT_PFHT680_v*",                           std::pair(315252, 325175) },
    { "HLT_PFHT780_v*",                           std::pair(315252, 325175) },
    { "HLT_PFHT890_v*",                           std::pair(315252, 325175) },
    { "HLT_PFHT1050_v*",                           std::pair(315252, 325175) },
    { "HLT_AK8PFJet40_v*",                        std::pair(315252, 325175) },
    { "HLT_AK8PFJet60_v*",                        std::pair(315252, 325175) },
    { "HLT_AK8PFJet80_v*",                        std::pair(315252, 325175) },
    { "HLT_AK8PFJet140_v*",                       std::pair(315252, 325175) },
    { "HLT_AK8PFJet200_v*",                       std::pair(315252, 325175) },
    { "HLT_AK8PFJet260_v*",                       std::pair(315252, 325175) },
    { "HLT_AK8PFJet320_v*",                       std::pair(315252, 325175) },
    { "HLT_AK8PFJet400_v*",                       std::pair(315252, 325175) },
    { "HLT_AK8PFJet450_v*",                       std::pair(315252, 325175) },
    { "HLT_AK8PFJet500_v*",                       std::pair(315252, 325175) },
  }},
};

const std::unordered_map<std::string, std::map<std::string, float>>
lumi_map = {
  { "2016", {
    { "lumiPlot",     36},
    { "lumi_fb",      35.9},
    { "lumi_pb",      35920},
    { "uncertainty",  2.5},
  }},
  { "2017", {
    { "lumiPlot",     41},
    { "lumi_fb",      41.5},
    { "lumi_pb",      41530},
    { "uncertainty",  2.3},
  }},
  { "2018", {
    { "lumiPlot",     59},
    { "lumi_fb",      59.7},
    { "lumi_pb",      59740},
    { "uncertainty",  2.5},
  }},
  { "RunII", {
    { "lumiPlot",     137},
    { "lumi_fb",      137.19},
    { "lumi_pb",      137190},
    { "uncertainty",  1.8},
  }},
};


const std::unordered_map<std::string, std::unordered_map<std::string, std::pair<std::string, std::string > > >
ScaleFactors_map = {//std::pair(filename,histname)
  { "2016", {
    { "Muon_Tracking",            std::pair("Muon_Tracking_SF_2016", "Tracking_SF")},
    { "Muon_Reconstruction",      std::pair("Muon_Reconstruction_SF_2016", "Reconstruction_SF")},
    { "Muon_HighPtID",            std::pair("Muon_ID_SF_2016_RunBCDEFGH", "NUM_HighPtID_DEN_genTracks_eta_pair_newTuneP_probe_pt")},
    { "Muon_TrkHighPtID",         std::pair("Muon_ID_SF_2016_TrkHighPt", "scalefactor")},
    { "Muon_Isolation",           std::pair("Muon_Isolation_SF_2016_RunBCDEFGH", "NUM_TightRelIso_DEN_TightIDandIPCut_eta_pt")},
    { "Muon_Trigger",             std::pair("Muon_Trigger_SF_2016_RunBCDEFGH", "Mu50_OR_TkMu50_PtEtaBins/abseta_pt_ratio")},
    { "Electron_LooseID",         std::pair("Electron_ID_SF_2016_loose", "EGamma_SF2D")},
    { "Electron_Reconstruction",  std::pair("Electron_Reconstruction_SF_2016", "EGamma_SF2D")},
    { "Electron_Trigger",         std::pair("Electron_Trigger_SF_2016", "SF_TH2F")},
    { "Jet_Tagger",               std::pair("Tagger_SF_2016", "SF_Var_FlavX_2016")},
  }},
  { "2017", {
    { "Muon_Tracking",            std::pair("Muon_Tracking_SF_2017", "Tracking_SF")},
    { "Muon_Reconstruction",      std::pair("Muon_Reconstruction_SF_2017", "Reconstruction_SF")},
    { "Muon_HighPtID",            std::pair("Muon_ID_SF_2017_RunBCDEF", "NUM_HighPtID_DEN_genTracks_pair_newTuneP_probe_pt_abseta")},
    { "Muon_TrkHighPtID",         std::pair("Muon_ID_SF_2017_RunBCDEF", "NUM_TrkHighPtID_DEN_genTracks_pair_newTuneP_probe_pt_abseta")},
    { "Muon_Isolation",           std::pair("Muon_Isolation_SF_2017_RunBCDEF", "NUM_TightRelTkIso_DEN_TrkHighPtID_pair_newTuneP_probe_pt_abseta")},
    { "Muon_Trigger",             std::pair("Muon_Trigger_SF_2017_RunBCDEF", "Mu50_PtEtaBins/abseta_pt_ratio")},
    { "Electron_LooseID",         std::pair("Electron_ID_SF_2017_loose", "EGamma_SF2D")},
    { "Electron_Reconstruction",  std::pair("Electron_Reconstruction_SF_2017", "EGamma_SF2D")},
    { "Electron_Trigger",         std::pair("Electron_Trigger_SF_2017", "SF_TH2F")},
    { "Jet_Tagger",               std::pair("Tagger_SF_2017", "SF_Var_FlavX_2017")},
  }},
  { "2018", {
    { "Muon_Tracking",            std::pair("Muon_Tracking_SF_2018", "Tracking_SF")},
    { "Muon_Reconstruction",      std::pair("Muon_Reconstruction_SF_2018", "Reconstruction_SF")},
    { "Muon_HighPtID",            std::pair("Muon_ID_SF_2018_RunABCD", "NUM_HighPtID_DEN_TrackerMuons_pair_newTuneP_probe_pt_abseta")},
    { "Muon_TrkHighPtID",         std::pair("Muon_ID_SF_2018_RunABCD", "NUM_TrkHighPtID_DEN_TrackerMuons_pair_newTuneP_probe_pt_abseta")},
    { "Muon_Isolation",           std::pair("Muon_Isolation_SF_2018_RunABCD", "NUM_TightRelTkIso_DEN_TrkHighPtID_pair_newTuneP_probe_pt_abseta")},
    { "Muon_Trigger",             std::pair("Muon_Trigger_SF_2018", "Mu50_OR_OldMu100_OR_TkMu100_PtEtaBins/abseta_pt_ratio")},
    { "Electron_LooseID",         std::pair("Electron_ID_SF_2018_loose", "EGamma_SF2D")},
    { "Electron_Reconstruction",  std::pair("Electron_Reconstruction_SF_2018", "EGamma_SF2D")},
    { "Electron_Trigger",         std::pair("Electron_Trigger_SF_2018", "SF_TH2F")},
    { "Jet_Tagger",               std::pair("Tagger_SF_2017", "SF_Var_FlavX_2017")},
  }},
};


const std::unordered_map<std::string, float> xsec_ref = { {"btag_DeepBoosted_H4qvsQCD", 0.001}, {"btag_DeepBoosted_H4qvsQCDptdep", 0.001}, {"default_value", 0.001}};

const std::vector<double> MassPoints = {600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500, 6000, 7000, 8000};


inline double PtToMass(double pt) { return pt*2;}

inline double MassToPt(double mass) { return mass/2;}

inline double PtToMass2(double pt) { return (pt-22.5)/0.36;}

inline double MassToPt2(double mass) { return mass*0.36 +22.5;}

const double MassDependentCut_old_value = -100;
inline double MassDependentCut_old(double pt) { double x = PtToMass(pt); return 5.03e-03+1.7e07*TMath::Power(x,-3);}
const double MassDependentCut_value = -10;
inline double MassDependentCut(double pt) { double x = PtToMass(pt); return 0.00238601739564+15407844.1061*TMath::Power(x,-3);}
const double MassDependentCut_cc_value = -40; // For Hcc
inline double MassDependentCut_cc(double pt) { double x = PtToMass(pt); return 0.110034415361+236887594.133*TMath::Power(x,-3);}

const double MassDependentCut2_cc_value = -20; // For Hcc
inline double MassDependentCut2_cc(double pt) { double x = PtToMass(pt); return 0.155+115000000*TMath::Power(x,-3);}
const double MassDependentCut3_cc_value = -30; // For Hcc
inline double MassDependentCut3_cc(double pt) { double x = PtToMass(pt); return 0.385+151000000*TMath::Power(x,-3);}

inline const char* BoolToString(bool b) { return b ? "true" : "false";}

enum taggers {NN_IsHiggs=1000, NN_IsQCD, NN_IsTop};

enum ParticleID { H=25, W=24, Z=23, t=6, b=5, c=4, s=3, d=2, u=1, g=21, ZPrime=9000001};
enum Matching { unknown=-1, noMatch=0, gluonMatch=1, qMatch=2, HMatch=3, HbbMatch=4, HWWMatch=5, HccMatch=6, HqqMatch=7, HggMatch=8, HtautauMatch=9, HZZMatch=10, topMatch=11, tWbMatch=12, WMatch=13, WqqMatch=14, WllMatch=15, ZMatch=16, ZqqMatch=17, ZllMatch=18};
enum MatchingStatus { Unknown=-1, NotMatched=0, MotherMatched=1, DaughterMatched=2, Hadronic=3, Hadronic1=4, Hadronic2=5, Hadronic3=6, SemiLep=7, FullLep=8, SemiMatched=9 };
enum Decay {nodecay, leptonic, semileptonic, hadronic, gluon, light, cc, bb, Wb, WW, ZZ, ee, mumu, tautau, ll, nunu, HWWfullLep, HWWfullHad, HWWsemiLep, ZH};
enum ZprimeDecay {nomatch=0, Zee=1, Zmumu=2, Zelse=3, HWW=10, Hbb=20, Hcc=30, Hgg=40, Helse=50, ZeeHWW=11, ZmumuHWW=12, ZelseHWW=13, ZeeHbb=21, ZmumuHbb=22, ZelseHbb=23, ZeeHcc=31, ZmumuHcc=32, ZelseHcc=33, ZeeHgg=41, ZmumuHgg=42, ZelseHgg=43, ZeeHelse=51, ZmumuHelse=52, ZelseHelse=53};

inline std::string DecayToString(const int & tagname) {
  if(tagname == nodecay)      return "nodecay";
  if(tagname == leptonic)     return "leptonic";
  if(tagname == semileptonic) return "semileptonic";
  if(tagname == hadronic)     return "hadronic";
  if(tagname == gluon)        return "gluon";
  if(tagname == light)        return "light";
  if(tagname == cc)           return "cc";
  if(tagname == bb)           return "bb";
  if(tagname == Wb)           return "Wb";
  if(tagname == WW)           return "WW";
  if(tagname == ZZ)           return "ZZ";
  if(tagname == ee)           return "ee";
  if(tagname == mumu)         return "mumu";
  if(tagname == tautau)       return "tautau";
  if(tagname == ll)           return "ll";
  if(tagname == nunu)         return "nunu";
  if(tagname == HWWfullLep)   return "HWWfullLep";
  if(tagname == HWWfullHad)   return "HWWfullHad";
  if(tagname == HWWsemiLep)   return "HWWsemiLep";
  return "unknown";
}

inline Matching FloatToMatching(const float & tagname_) { return static_cast<Matching>(int(tagname_)); }

inline MatchingStatus FloatToMatchingStatus(const float & tagname_) { return static_cast<MatchingStatus>(int(tagname_)); }

inline Matching StringToMatching(const std::string & tagname) {
  if(tagname == "noMatch")      return noMatch;
  if(tagname == "gluonMatch")   return gluonMatch;
  if(tagname == "qMatch")       return qMatch;
  if(tagname == "HMatch")       return HMatch;
  if(tagname == "HbbMatch")     return HbbMatch;
  if(tagname == "HWWMatch")     return HWWMatch;
  if(tagname == "HqqMatch")     return HqqMatch;
  if(tagname == "HggMatch")     return HggMatch;
  if(tagname == "HccMatch")     return HccMatch;
  if(tagname == "HtautauMatch") return HtautauMatch;
  if(tagname == "HZZMatch")     return HZZMatch;
  if(tagname == "topMatch")     return topMatch;
  if(tagname == "tWbMatch")     return tWbMatch;
  if(tagname == "WMatch")       return WMatch;
  if(tagname == "WqqMatch")     return WqqMatch;
  if(tagname == "WllMatch")     return WllMatch;
  if(tagname == "ZMatch")       return ZMatch;
  if(tagname == "ZqqMatch")     return ZqqMatch;
  if(tagname == "ZllMatch")     return ZllMatch;
  return unknown;
}

inline std::string MatchingToString(const float & tagname_) {
  Matching tagname = FloatToMatching(tagname_);
  if(tagname == noMatch)      return "noMatch";
  if(tagname == gluonMatch)   return "gluonMatch";
  if(tagname == qMatch)       return "qMatch";
  if(tagname == HMatch)       return "HMatch";
  if(tagname == HbbMatch)     return "HbbMatch";
  if(tagname == HWWMatch)     return "HWWMatch";
  if(tagname == HqqMatch)     return "HqqMatch";
  if(tagname == HccMatch)     return "HccMatch";
  if(tagname == HggMatch)     return "HggMatch";
  if(tagname == HtautauMatch) return "HtautauMatch";
  if(tagname == HZZMatch)     return "HZZMatch";
  if(tagname == topMatch)     return "topMatch";
  if(tagname == tWbMatch)     return "tWbMatch";
  if(tagname == WMatch)       return "WMatch";
  if(tagname == WqqMatch)     return "WqqMatch";
  if(tagname == WllMatch)     return "WllMatch";
  if(tagname == ZMatch)       return "ZMatch";
  if(tagname == ZqqMatch)     return "ZqqMatch";
  if(tagname == ZllMatch)     return "ZllMatch";
  return "unknown";
}


inline std::string MatchingStatusToString(const float & tagname_) {
  MatchingStatus tagname = FloatToMatchingStatus(tagname_);
  if(tagname == NotMatched)       return "NotMatched";
  if(tagname == MotherMatched)    return "MotherMatched";
  if(tagname == DaughterMatched)  return "DaughterMatched";
  if(tagname == Hadronic)         return "Hadronic";
  if(tagname == Hadronic1)        return "Hadronic1";
  if(tagname == Hadronic2)        return "Hadronic2";
  if(tagname == Hadronic3)        return "Hadronic3";
  if(tagname == SemiLep)          return "SemiLep";
  if(tagname == FullLep)          return "FullLep";
  if(tagname == SemiMatched)     return "SemiMatched";
  return "Unknown";
}

inline MatchingStatus StringToMatchingStatus(const std::string & tagname) {
  if(tagname == "NotMatched")       return NotMatched;
  if(tagname == "MotherMatched")    return MotherMatched;
  if(tagname == "DaughterMatched")  return DaughterMatched;
  if(tagname == "Hadronic")         return Hadronic;
  if(tagname == "Hadronic1")        return Hadronic1;
  if(tagname == "Hadronic2")        return Hadronic2;
  if(tagname == "Hadronic3")        return Hadronic3;
  if(tagname == "SemiLep")          return SemiLep;
  if(tagname == "FullLep")          return FullLep;
  if(tagname == "SemiMatched")      return SemiMatched;
  return Unknown;
}


inline std::string ZprimeDecayToString(const int & tagname) {
  if(tagname == nomatch)     return "nomatch";
  if(tagname == Zee)         return "Zee";
  if(tagname == Zmumu)       return "Zmumu";
  if(tagname == Zelse)       return "Zelse";
  if(tagname == HWW)         return "HWW";
  if(tagname == Hbb)         return "Hbb";
  if(tagname == Hcc)         return "Hcc";
  if(tagname == Hgg)         return "Hgg";
  if(tagname == Helse)       return "Helse";
  if(tagname == ZeeHWW)      return "ZeeHWW";
  if(tagname == ZmumuHWW)    return "ZmumuHWW";
  if(tagname == ZelseHWW)    return "ZelseHWW";
  if(tagname == ZeeHbb)      return "ZeeHbb";
  if(tagname == ZmumuHbb)    return "ZmumuHbb";
  if(tagname == ZelseHbb)    return "ZelseHbb";
  if(tagname == ZeeHcc)      return "ZeeHcc";
  if(tagname == ZmumuHcc)    return "ZmumuHcc";
  if(tagname == ZelseHcc)    return "ZelseHcc";
  if(tagname == ZeeHgg)      return "ZeeHgg";
  if(tagname == ZmumuHgg)    return "ZmumuHgg";
  if(tagname == ZelseHgg)    return "ZelseHgg";
  if(tagname == ZeeHelse)    return "ZeeHelse";
  if(tagname == ZmumuHelse)  return "ZmumuHelse";
  if(tagname == ZelseHelse)  return "ZelseHelse";
  return "unknown";
}
