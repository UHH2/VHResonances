#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/TopJet.h"
#include "UHH2/HiggsToWWTagger/include/constants.hpp"
#include "UHH2/HiggsToWWTagger/include/Utils.hpp"


Matching FloatToMatching(const TopJet & jet) { return static_cast<Matching>(int(jet.get_tag(TopJet::Matching))); }
MatchingStatus FloatToMatchingStatus(const TopJet & jet) { return static_cast<MatchingStatus>(int(jet.get_tag(TopJet::MatchingStatus))); }

bool XOR( bool a, bool b) { return (!a&&b) || (!b && a); };

bool isLeptonic(int pdgId) { return (fabs(pdgId)>= 11 && fabs(pdgId)<=18) ? true : false; }
bool isHadronic(int pdgId) { return (fabs(pdgId) <= 5) ? true : false; }

bool DobleDecay(int pdgId1, int pdgId2, Decay decay) {
  if (decay==nodecay) return true;
  if ( isLeptonic(pdgId1) && isLeptonic(pdgId2) && decay==leptonic)     return true;
  if ( isLeptonic(pdgId1) && isHadronic(pdgId2) && decay==semileptonic) return true;
  if ( isHadronic(pdgId1) && isLeptonic(pdgId2) && decay==semileptonic) return true;
  if ( isHadronic(pdgId1) && isHadronic(pdgId2) && decay==hadronic)     return true;
  if ( fabs(pdgId1) == 21 && fabs(pdgId2) == 21 && decay==gluon)        return true;
  if ( fabs(pdgId1) <=  4 && fabs(pdgId2) <=  4 && decay==light)        return true;
  if ( fabs(pdgId1) ==  5 && fabs(pdgId2) ==  5 && decay==bb)           return true;
  if ( fabs(pdgId1) == 23 && fabs(pdgId2) == 23 && decay==ZZ)           return true;
  if ( fabs(pdgId1) == 24 && fabs(pdgId2) == 24 && decay==WW)           return true;
  if ( fabs(pdgId1) == 24 && fabs(pdgId2) ==  5 && decay==Wb)           return true;
  if ( fabs(pdgId1) ==  5 && fabs(pdgId2) == 24 && decay==Wb)           return true;
  if ( fabs(pdgId1) == 11 && fabs(pdgId2) == 11 && decay==ee)           return true;
  if ( fabs(pdgId1) == 13 && fabs(pdgId2) == 13 && decay==mumu)         return true;
  if ( fabs(pdgId1) == 15 && fabs(pdgId2) == 15 && decay==tautau)       return true;
  if ( fabs(pdgId1) == 11 && fabs(pdgId2) == 11 && decay==ll)           return true;
  if ( fabs(pdgId1) == 13 && fabs(pdgId2) == 13 && decay==ll)           return true;
  if ( fabs(pdgId1) == 12 && fabs(pdgId2) == 12 && decay==nunu)         return true;
  if ( fabs(pdgId1) == 14 && fabs(pdgId2) == 14 && decay==nunu)         return true;
  if ( fabs(pdgId1) == 15 && fabs(pdgId2) == 16 && decay==nunu)         return true;
  if ( fabs(pdgId1) == 23 && fabs(pdgId2) == 25 && decay==ZH)           return true;
  if ( fabs(pdgId1) == 25 && fabs(pdgId2) == 23 && decay==ZH)           return true;
  return false;
}

bool FindInString(const std::string& search, const std::string& str) {return str.find(search)!=std::string::npos ;}
