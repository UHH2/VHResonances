#pragma once

#include <string>
#include <unordered_map>

#include "TMath.h"

// DEFINE GLOBAL VARIABLES

const float ZMASS  = 91.0;
const float ZWIDTH = 10.0;

const int min_leptons = 2;
const float min_jet_pt = 30.0;
const float min_topjet_pt = 200.0;
// const float min_lepton_pt = 30.0;
const float min_lepton_pt = 52.0;
const float min_lepton_eta = 2.4;
// const float max_lepton_iso = 0.15;
// const float min_jet_dilep_delta_phi = 2.7;

const float min_jet_dilep_pt = 150;
const float min_DR_dilep = 0.0;
const float max_DR_dilep = 1.0;
const float min_jet_dilep_delta_phi = M_PI/2;
const float max_jet_dilep_delta_phi = M_PI;

const float min_Z_pt_ZH_mass = 0.2;

// Taken from :
// https://twiki.cern.ch/twiki/bin/view/CMS/EgHLTRunIISummary
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/B2GMuon
// https://twiki.cern.ch/twiki/bin/view/CMS/MuonHLT2017

const std::unordered_map<std::string, std::map<std::string, std::pair<int, int>>>
Trigger_run_validity = {
  { "2016", {
    { "HLT_Mu50_v*", std::pair(272760, 284044) }, // 272007
    { "HLT_TkMu50_v*", std::pair(274954, 284044) },
    { "HLT_Ele27_WPTight_Gsf_v*", std::pair(272760, 284044) }, // 273158
    { "HLT_Ele115_CaloIdVT_GsfTrkIdT_v*", std::pair(272760, 284044) }, // 273158, 284044
  }},
  { "2017", {
    { "HLT_Mu50_v*", std::pair(296070, 306460) }, // 297020, 306462
    { "HLT_OldMu100_v*", std::pair(299368, 306460) }, //306462 // 299370
    { "HLT_TkMu100_v*", std::pair(299368, 306460) }, // 306462
    { "HLT_Ele35_WPTight_Gsf_v*", std::pair(296070, 306460)}, // 297050, 306460
    { "HLT_Ele115_CaloIdVT_GsfTrkIdT_v*", std::pair(299368, 306460)},
  }},
  { "2018", {
    { "HLT_Mu50_v*", std::pair(315252, 325175) }, // 315252, 325273
    { "HLT_OldMu100_v*", std::pair(315252, 325175) }, // 315252, 325273
    { "HLT_TkMu100_v*", std::pair(315252, 325175) }, // 315252, 325273
    { "HLT_Ele32_WPTight_Gsf_v*", std::pair(315252, 325175) }, // 315257, 325172
    { "HLT_Ele115_CaloIdVT_GsfTrkIdT_v*", std::pair(315252, 325175)}, // 315257,325172
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

const std::unordered_map<std::string, float> xsec_ref = { {"btag_DeepBoosted_H4qvsQCD", 0.001}, {"btag_DeepBoosted_H4qvsQCDptdep", 0.001}, {"default_value", 0.001}};

const std::vector<double> MassPoints = {600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500, 6000, 7000, 8000};


inline double PtToMass(double pt) { return pt*2;}

inline double MassToPt(double mass) { return mass/2;}

inline double PtToMass2(double pt) { return (pt-22.5)/0.36;}

inline double MassToPt2(double mass) { return mass*0.36 +22.5;}




// const double PtDepentdentCut_value_x3 = -1000; // Fit on pt exp(x0+x3)
// inline double PtDepentdentCut_x3(double pt)   { return 5.65e-03+1.90e+07*TMath::Power(PtToMass(pt),-3);}
// const double PtDepentdentCut_value_x2x3 = -2000; // Fit on pt exp(x0+x2+x3)
// inline double PtDepentdentCut_x2x3(double pt) { return 9.31e-03-1.46e+04*TMath::Power(PtToMass(pt),-2)+2.60e+07*TMath::Power(PtToMass(pt),-3);}
// const double PtDepentdentCut_value_x1x3 = -3000; // Fit on pt exp(x0+x1+x3)
// inline double PtDepentdentCut_x1x3(double pt) { return 1.25e-02-1.31e+01*TMath::Power(PtToMass(pt),-1)+2.15e+07*TMath::Power(PtToMass(pt),-3);}


// const double MassDepentdentCut_value_x3 = -100; // Fit on Mass exp(x0+x3). pt = m/2
// inline double MassDepentdentCut_x3(double pt) { double x = PtToMass(pt); return 5.08e-03+2.72e07*TMath::Power(x,-3);}
// const double MassDepentdentCut2_value_x3 = -150; // Fit on Mass exp(x0+x3). Fit on pt vs m
// inline double MassDepentdentCut2_x3(double pt) { double x = PtToMass(pt); return 5.08e-03+2.72e07*TMath::Power(x,-3);}
// const double MassDepentdentCut_value_x2x3 = -200; // Fit on Mass exp(x0+x2+x3). pt = m/2
// inline double MassDepentdentCut_x2x3(double pt) { double x = PtToMass(pt); return 4.06e-03+1.15e+04*TMath::Power(x,-2)+1.83e+07*TMath::Power(x,-3);}
// const double MassDepentdentCut_value_x1x3 = -300; // Fit on Mass exp(x0+x1+x3). pt = m/2
// inline double MassDepentdentCut_x1x3(double pt) { double x = PtToMass(pt); return 3.32e-03+5.83*TMath::Power(x,-1)+2.41e+07*TMath::Power(x,-3);}
// const double MassDepentdentCut_value_x1x2 = -400; // Fit on Mass exp(x0+x1+x2). pt = m/2
// inline double MassDepentdentCut_x1x2(double pt) { double x = PtToMass(pt); return 6.72e-03-1.97e+01*TMath::Power(x,-1)+4.83e+04*TMath::Power(x,-2);}


const double PtDepentdentCut_value_x3 = -1000; // Fit on pt exp(x0+x3)
inline double PtDepentdentCut_x3(double pt)   { return 5.8e-03+1.9e+07*TMath::Power(PtToMass(pt),-3);}
const double PtDepentdentCut_value_x2x3 = -2000; // Fit on pt exp(x0+x2+x3)
inline double PtDepentdentCut_x2x3(double pt) { return 8.6e-03-1.3e+04*TMath::Power(PtToMass(pt),-2)+2.5e+07*TMath::Power(PtToMass(pt),-3);}
const double PtDepentdentCut_value_x1x3 = -3000; // Fit on pt exp(x0+x1+x3)
inline double PtDepentdentCut_x1x3(double pt) { return 1.1e-02-1.1e+01*TMath::Power(PtToMass(pt),-1)+2.1e+07*TMath::Power(PtToMass(pt),-3);}
const double PtDepentdentCut_value_x1x2 = -4000; // Fit on pt exp(x0+x1+x2)
inline double PtDepentdentCut_x1x2(double pt) { return 2.4e-02-6.7e+01*TMath::Power(PtToMass(pt),-1)+6.7e+04*TMath::Power(PtToMass(pt),-2);}
const double PtDepentdentCut_value_x2 = -5000; // Fit on pt exp(x0+x2)
inline double PtDepentdentCut_x2(double pt) { return 2.3e-03+2.3e04*TMath::Power(PtToMass(pt),-2);}
const double PtDepentdentCut_value_sqrt = -6000; // Fit on pt [0]+[2]*pow(x,[1])
inline double PtDepentdentCut_sqrt(double pt) { return 3.2e-03+7.8e+05*TMath::Power(PtToMass(pt),-2.5e+00);}


const double MassDepentdentCut_value_x3 = -100; // Fit on Mass exp(x0+x3). pt = m/2
inline double MassDepentdentCut_x3(double pt) { double x = PtToMass(pt); return 5.03e-03+2.7e07*TMath::Power(x,-3);}
const double MassDepentdentCut2_value_x3 = -150; // Fit on Mass exp(x0+x3). Fit on pt vs m
inline double MassDepentdentCut2_x3(double pt) { double x = PtToMass(pt); return 5.03e-03+2.7e07*TMath::Power(x,-3);}
const double MassDepentdentCut_value_x2x3 = -200; // Fit on Mass exp(x0+x2+x3). pt = m/2
inline double MassDepentdentCut_x2x3(double pt) { double x = PtToMass(pt); return 4.3e-03+9.7e+03*TMath::Power(x,-2)+2.0e+07*TMath::Power(x,-3);}
const double MassDepentdentCut_value_x1x3 = -300; // Fit on Mass exp(x0+x1+x3). pt = m/2
inline double MassDepentdentCut_x1x3(double pt) { double x = PtToMass(pt); return 3.8e-03+4.8e+00*TMath::Power(x,-1)+2.5e+07*TMath::Power(x,-3);}
const double MassDepentdentCut_value_x1x2 = -400; // Fit on Mass exp(x0+x1+x2). pt = m/2
inline double MassDepentdentCut_x1x2(double pt) { double x = PtToMass(pt); return 7.8e-03-2.3e+01*TMath::Power(x,-1)+5.0e+04*TMath::Power(x,-2);}
const double MassDepentdentCut_value_x2 = -500; // Fit on Mass exp(x0+x2). pt = m/2
inline double MassDepentdentCut_x2(double pt) { double x = PtToMass(pt); return 2.3e-03+2.3e04*TMath::Power(x,-2);}
const double MassDepentdentCut_value_sqrt = -600; // Fit on Mass [0]+[2]*pow(x,[1]). pt = m/2
inline double MassDepentdentCut_sqrt(double pt) { double x = PtToMass(pt); return 2.5e-03+5.6e+04*TMath::Power(x,-2.1e+00);}

inline const char * const BoolToString(bool b) { return b ? "true" : "false";}

enum taggers {NN_IsHiggs=1000, NN_IsQCD, NN_IsTop};

enum ParticleID { H=25, W=24, Z=23, t=6, b=5, c=4, s=3, d=2, u=1, g=21, ZPrime=9000001};
enum Matching { unknown=-1, noMatch=0, gluonMatch=1, qMatch=2, HMatch=3, HbbMatch=4, HWWMatch=5, HqqMatch=6, HtautauMatch=7, HZZMatch=8, topMatch=9, tWbMatch=10, WMatch=11, WqqMatch=12, WllMatch=13, ZMatch=14, ZqqMatch=15, ZllMatch=16};
enum MatchingStatus { Unknown=-1, NotMatched=0, MotherMatched=1, DaughterMatched=2, Hadronic=3, Hadronic1=4, Hadronic2=5, Hadronic3=6, SemiLep=7, FullLep=8, SemiMatched=9 };
enum Decay {nodecay, leptonic, semileptonic, hadronic, gluon, light, bb, Wb, WW, ZZ, ee, mumu, tautau, ll, nunu, HWWfullLep, HWWfullHad, HWWsemiLep, ZH};
enum ZprimeDecay {nomatch=0, Zee=1, Zmumu=2, Zelse=3, HWW=10, Hbb=20, Helse=30, ZeeHWW=11, ZmumuHWW=12, ZelseHWW=13, ZeeHbb=21, ZmumuHbb=22, ZelseHbb=23, ZeeHelse=31, ZmumuHelse=32, ZelseHelse=33};

inline std::string DecayToString(const int & tagname) {
  if(tagname == nodecay)      return "nodecay";
  if(tagname == leptonic)     return "leptonic";
  if(tagname == semileptonic) return "semileptonic";
  if(tagname == hadronic)     return "hadronic";
  if(tagname == gluon)        return "gluon";
  if(tagname == light)        return "light";
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
