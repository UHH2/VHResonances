#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Utils.h"
#include "UHH2/common/include/Utils.h"

#include "UHH2/VHResonances/include/constants.hpp"
#include "UHH2/VHResonances/include/Utils.hpp"
#include "UHH2/VHResonances/include/GenLevelJetMatch.h"

#include <stdexcept>

using namespace std;
using namespace uhh2;

bool GenLevelJetMatch::MatchGenPart(const Event & event, TopJet& jet, double Dr) {

  jet.set_tag(TopJet::Matching, (float)noMatch);
  jet.set_tag(TopJet::MatchingStatus, (float)NotMatched);

  for (auto& case_: {H,t,W,Z,b,c,s,d,u,g}) {

    if ( FloatToMatching(jet)!=noMatch && FloatToMatching(jet)!=qMatch && FloatToMatchingStatus(jet)!=MotherMatched ) continue;

    std::vector<GenParticle> vec(event.genparticles->size());
    auto it = std::copy_if ((*event.genparticles).begin(), (*event.genparticles).end(), vec.begin(), [case_](GenParticle i){ return fabs(i.pdgId())==case_;} );
    vec.resize(std::distance(vec.begin(),it));

    for (auto& gp: vec ) {
      if (deltaR(gp,jet)>Dr) continue;

      if      (fabs(gp.pdgId()) == H) jet.set_tag(TopJet::Matching, (float)HMatch);
      else if (fabs(gp.pdgId()) == W) jet.set_tag(TopJet::Matching, (float)WMatch);
      else if (fabs(gp.pdgId()) == Z) jet.set_tag(TopJet::Matching, (float)ZMatch);
      else if (fabs(gp.pdgId()) == t) jet.set_tag(TopJet::Matching, (float)topMatch);
      else if (fabs(gp.pdgId()) == b) jet.set_tag(TopJet::Matching, (float)qMatch);
      else if (fabs(gp.pdgId()) == c) jet.set_tag(TopJet::Matching, (float)qMatch);
      else if (fabs(gp.pdgId()) == s) jet.set_tag(TopJet::Matching, (float)qMatch);
      else if (fabs(gp.pdgId()) == d) jet.set_tag(TopJet::Matching, (float)qMatch);
      else if (fabs(gp.pdgId()) == u) jet.set_tag(TopJet::Matching, (float)qMatch);
      else if (fabs(gp.pdgId()) == g) jet.set_tag(TopJet::Matching, (float)qMatch); // TODO gluonMatch
      else continue;

      jet.set_tag(TopJet::MatchingStatus, (float)MotherMatched);
      if (FloatToMatching(jet)==qMatch) continue;

      const GenParticle* d1 = gp.daughter(event.genparticles, 1);
      const GenParticle* d2 = gp.daughter(event.genparticles, 2);
      if (!d1 || !d2) continue;
      if (! (deltaR(*d1,jet)<Dr && deltaR(*d2,jet)<Dr)) continue;
      int ID1 = d1->pdgId();
      int ID2 = d2->pdgId();

      if      (fabs(gp.pdgId())== H && DobleDecay(ID1, ID2,WW))       jet.set_tag(TopJet::Matching, (float)HWWMatch);
      else if (fabs(gp.pdgId())== H && DobleDecay(ID1, ID2,bb))       jet.set_tag(TopJet::Matching, (float)HbbMatch);
      else if (fabs(gp.pdgId())== H && DobleDecay(ID1, ID2,light))    jet.set_tag(TopJet::Matching, (float)HqqMatch);
      else if (fabs(gp.pdgId())== H && DobleDecay(ID1, ID2,gluon))    jet.set_tag(TopJet::Matching, (float)HqqMatch);
      else if (fabs(gp.pdgId())== H && DobleDecay(ID1, ID2,ZZ))       jet.set_tag(TopJet::Matching, (float)HZZMatch);
      else if (fabs(gp.pdgId())== H && DobleDecay(ID1, ID2,tautau))   jet.set_tag(TopJet::Matching, (float)HtautauMatch);
      else if (fabs(gp.pdgId())== t && DobleDecay(ID1, ID2,Wb))       jet.set_tag(TopJet::Matching, (float)tWbMatch);
      else if (fabs(gp.pdgId())== W && DobleDecay(ID1, ID2,hadronic)) jet.set_tag(TopJet::Matching, (float)WqqMatch);
      else if (fabs(gp.pdgId())== W && DobleDecay(ID1, ID2,leptonic)) jet.set_tag(TopJet::Matching, (float)WllMatch);
      else if (fabs(gp.pdgId())== Z && DobleDecay(ID1, ID2,hadronic)) jet.set_tag(TopJet::Matching, (float)ZqqMatch);
      else if (fabs(gp.pdgId())== Z && DobleDecay(ID1, ID2,leptonic)) jet.set_tag(TopJet::Matching, (float)ZllMatch);
      else continue;

      jet.set_tag(TopJet::MatchingStatus, (float)DaughterMatched);
      if (!DobleDecay(ID1, ID2,WW) && !DobleDecay(ID1, ID2,Wb) && !DobleDecay(ID1, ID2,ZZ)) continue;

      std::vector<Decay> Daughters = {nodecay,nodecay}; int nInside = 0;

      for (int i = 0; i <= 1; i++) {
        const GenParticle* dau = gp.daughter(event.genparticles, i+1);
        if (fabs(dau->pdgId())==b) {Daughters[i] = hadronic; nInside++; continue;}
        const GenParticle* sd1 = dau->daughter(event.genparticles, 1);
        const GenParticle* sd2 = dau->daughter(event.genparticles, 2);
        if (!sd1 || !sd2) continue;
        if (DobleDecay(sd1->pdgId(), sd2->pdgId(),hadronic)) Daughters[i] = hadronic;
        if (DobleDecay(sd1->pdgId(), sd2->pdgId(),leptonic)) Daughters[i] = leptonic;
        if (deltaR(*sd1,jet)<Dr) nInside++;
        if (deltaR(*sd2,jet)<Dr) nInside++;
      }

      if      (Daughters[0]==leptonic && Daughters[1]==leptonic)                jet.set_tag(TopJet::MatchingStatus, (float)FullLep);
      else if (Daughters[0]==leptonic && Daughters[1]==hadronic)                jet.set_tag(TopJet::MatchingStatus, (float)SemiLep);
      else if (Daughters[0]==hadronic && Daughters[1]==leptonic)                jet.set_tag(TopJet::MatchingStatus, (float)SemiLep);
      else if (Daughters[0]==hadronic && Daughters[1]==hadronic && nInside==4)  jet.set_tag(TopJet::MatchingStatus, (float)Hadronic);
      else if (Daughters[0]==hadronic && Daughters[1]==hadronic && nInside==3)  jet.set_tag(TopJet::MatchingStatus, (float)Hadronic3);
      else if (Daughters[0]==hadronic && Daughters[1]==hadronic && nInside==2)  jet.set_tag(TopJet::MatchingStatus, (float)Hadronic2);
      else if (Daughters[0]==hadronic && Daughters[1]==hadronic && nInside==1)  jet.set_tag(TopJet::MatchingStatus, (float)Hadronic1);
      else if (nInside!=0)                                                      jet.set_tag(TopJet::MatchingStatus, (float)SemiMatched);

    }

  }

  return true;
}


GenLevelJetMatch::GenLevelJetMatch(float Dr_): Dr(Dr_) {

}

GenLevelJetMatch::GenLevelJetMatch(Context & ctx, const string& jetCollection) {
  Dr = string2double(ctx.get("MatchingRadius","0.8"));
  h_topjets_ = ctx.get_handle<std::vector<TopJet>>(jetCollection);
}

bool GenLevelJetMatch::process(Event & event){
  if (event.isRealData) return true;

  std::vector<TopJet>* topjets(0);
  if (event.is_valid(h_topjets_)) topjets = &event.get(h_topjets_);

  for (auto& jet : *topjets) MatchGenPart(event, jet, Dr);

  return true;
}
