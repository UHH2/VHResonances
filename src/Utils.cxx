#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/TopJet.h"
#include "UHH2/VHResonances/include/constants.hpp"
#include "UHH2/VHResonances/include/Utils.hpp"


Matching FloatToMatching(const TopJet & jet) { return static_cast<Matching>(int(jet.get_tag(TopJet::Matching))); }
MatchingStatus FloatToMatchingStatus(const TopJet & jet) { return static_cast<MatchingStatus>(int(jet.get_tag(TopJet::MatchingStatus))); }

bool XOR( bool a, bool b) { return (!a&&b) || (!b && a); };

bool isLeptonic(int pdgId) { return (fabs(pdgId)>= 11 && fabs(pdgId)<=18) ? true : false; }
bool isHadronic(int pdgId) { return (fabs(pdgId) <= 5) ? true : false; }

bool DoubleDecay(int pdgId1, int pdgId2, Decay decay) {
  if (decay==nodecay) return true;
  if ( isLeptonic(pdgId1) && isLeptonic(pdgId2) && decay==leptonic)     return true;
  if ( isLeptonic(pdgId1) && isHadronic(pdgId2) && decay==semileptonic) return true;
  if ( isHadronic(pdgId1) && isLeptonic(pdgId2) && decay==semileptonic) return true;
  if ( isHadronic(pdgId1) && isHadronic(pdgId2) && decay==hadronic)     return true;
  if ( fabs(pdgId1) == 21 && fabs(pdgId2) == 21 && decay==gluon)        return true;
  if ( fabs(pdgId1) <=  3 && fabs(pdgId2) <=  3 && decay==light)        return true;
  if ( fabs(pdgId1) ==  4 && fabs(pdgId2) ==  4 && decay==cc)           return true;
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

int FindInVector(const std::vector<std::string>& vec, const std::string& el) {
  int index = -1;
  // Find given element in vector
  auto it = std::find(vec.begin(), vec.end(), el);
  if (it != vec.end()) index = distance(vec.begin(), it);
  return index;
}

const std::string MyString(const std::string & tag) {return tag;};// Needed just to deal with preprocessing macros


double GetQCD(const TopJet& j, bool isPN, bool isMD) {
  if (isMD) {
    if (isPN) {
      return (j.btag_MassDecorrelatedParticleNetJetTags_probQCDbb()+j.btag_MassDecorrelatedParticleNetJetTags_probQCDcc()+j.btag_MassDecorrelatedParticleNetJetTags_probQCDb()+j.btag_MassDecorrelatedParticleNetJetTags_probQCDc()+j.btag_MassDecorrelatedParticleNetJetTags_probQCDothers());
    } else {
      return (j.btag_MassDecorrelatedDeepBoosted_probQCDb()+j.btag_MassDecorrelatedDeepBoosted_probQCDbb()+j.btag_MassDecorrelatedDeepBoosted_probQCDc()+j.btag_MassDecorrelatedDeepBoosted_probQCDcc()+j.btag_MassDecorrelatedDeepBoosted_probQCDothers());
    }
  }
  else {
    if (isPN) {
      return (j.btag_ParticleNetJetTags_probQCDbb()+j.btag_ParticleNetJetTags_probQCDcc()+j.btag_ParticleNetJetTags_probQCDb()+j.btag_ParticleNetJetTags_probQCDc()+j.btag_ParticleNetJetTags_probQCDothers());
    } else {
      return (j.btag_DeepBoosted_probQCDb()+j.btag_DeepBoosted_probQCDbb()+j.btag_DeepBoosted_probQCDc()+j.btag_DeepBoosted_probQCDcc()+j.btag_DeepBoosted_probQCDothers());
    }
  }
}

double GetHccvsQCD(const TopJet& j, bool isPN, bool isMD) {
  double Hcc = 0;
  if (isMD) Hcc = isPN? j.btag_MassDecorrelatedParticleNetJetTags_probXcc(): j.btag_MassDecorrelatedDeepBoosted_probHcc();
  else Hcc = isPN? j.btag_ParticleNetJetTags_probHcc() : j.btag_DeepBoosted_probHcc();
  return Hcc/(Hcc+GetQCD(j,isPN,isMD));
}

double GetZHccvsQCD(const TopJet& j, bool isPN, bool isMD) {
  double ZHcc = 0;
  if (isMD) ZHcc = isPN? j.btag_MassDecorrelatedParticleNetJetTags_probXcc() : j.btag_MassDecorrelatedDeepBoosted_probZcc()+j.btag_MassDecorrelatedDeepBoosted_probHcc();
  else ZHcc = isPN? (j.btag_MassDecorrelatedParticleNetJetTags_probXcc()+j.btag_ParticleNetJetTags_probZcc()) : (j.btag_DeepBoosted_probZcc()+j.btag_DeepBoosted_probHcc());
  return ZHcc/(ZHcc+GetQCD(j,isPN,isMD));
}
