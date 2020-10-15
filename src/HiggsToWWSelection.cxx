#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/VHResonances/include/HiggsToWWSelection.h"

#include <stdexcept>
#include <set>

using namespace std;
using namespace uhh2;

JetDiLeptonPhiAngularSelection::JetDiLeptonPhiAngularSelection (float pt_min_, float phi_min_, float phi_max_, float min_Dphi_jet_MET_, TString lepton_, const Event::Handle<vector<TopJet> > & topjetcollection_): pt_min(pt_min_), phi_min(phi_min_), phi_max(phi_max_), min_Dphi_jet_MET(min_Dphi_jet_MET_), lepton(lepton_), topjetcollection(topjetcollection_){}

bool JetDiLeptonPhiAngularSelection::passes(const Event& event){

  const auto & jets = event.get(topjetcollection);

  vector<Particle> leptons;

  if (lepton == "muons") leptons.assign(event.muons->begin(), event.muons->end());
  else if (lepton == "electrons") leptons.assign(event.electrons->begin(), event.electrons->end());
  else if (lepton != "invisible") throw logic_error("JetDiLeptonPhiAngularSelection: Impossible case");

  if(jets.size() < 1 ) throw std::runtime_error("JetDiLeptonPhiAngularSelection::JetDiLeptonPhiAngularSelection -- unexpected number of jets");

  if (lepton == "invisible"){
    // invisblechannel
    for(const auto & jet: jets){
      double Dphi = deltaPhi(jet, *event.met);
      if (Dphi > min_Dphi_jet_MET) return true;
    }
    return false;
  }

  // leptonchannel
  if(leptons.size() < min_leptons ) throw std::runtime_error("JetDiLeptonPhiAngularSelection::JetDiLeptonPhiAngularSelection -- unexpected number of leptons");
  for (unsigned int i = 0; i < leptons.size(); i++) {
    Particle lep1 = leptons.at(i);
    for (unsigned int j = i+1; j < leptons.size(); j++) {
      Particle lep2 = leptons.at(j);
      auto diLep = lep1.v4() + lep2.v4();
      if( (fabs(diLep.M() - ZMASS) < ZWIDTH) && (lep1.charge() + lep2.charge() == 0) && diLep.pt()>pt_min ) {
        for(const auto & jet: jets){
          auto Dphi = deltaPhi(diLep, jet);
          if( phi_min < Dphi  && Dphi< phi_max) return true;
        }
      }
    }
  }
  return false;
}




DeltaRDiLepton::DeltaRDiLepton (float DR_min_, float DR_max_, TString lepton_): DR_min(DR_min_), DR_max(DR_max_), lepton(lepton_) { };

bool DeltaRDiLepton::passes(const Event& event){

  vector<Particle> leptons;

  if (lepton == "muons") leptons.assign(event.muons->begin(), event.muons->end());
  else if (lepton == "electrons") leptons.assign(event.electrons->begin(), event.electrons->end());
  else throw logic_error("DeltaRDiLepton: Impossible case");

  if(leptons.size() < min_leptons ) throw std::runtime_error("DeltaRDiLepton::DeltaRDiLepton -- unexpected number of leptons");
  for (unsigned int i = 0; i < leptons.size(); i++) {
    Particle lep1 = leptons.at(i);
    for (unsigned int j = i+1; j < leptons.size(); j++) {
      Particle lep2 = leptons.at(j);
      auto DR = uhh2::deltaR(lep1, lep2);
      if (DR< DR_max &&  DR> DR_min) return true;
    }
  }
  return false;
}

/* Remove events if the number of jets with Delta Phi(jet, MET) > minDeltaPhi is outside the interval [minJets, maxJets].
Can pass -1 as maxJets to require all jets to pass.
*/
DeltaPhiJetMETCut::DeltaPhiJetMETCut(Context& ctx, string jetCollection_, float minDeltaPhi_, int minJets_, int maxJets_): jetCollection(jetCollection_), minDeltaPhi(minDeltaPhi_), minJets(minJets_), maxJets(maxJets_){
    h_jets = ctx.get_handle<vector<Jet>>(jetCollection);
    h_topjets = ctx.get_handle<vector<TopJet>>(jetCollection);
};

bool DeltaPhiJetMETCut::passes(const Event& event){
  vector<Jet> jets;
  if (event.is_valid(h_topjets)) jets.assign((&event.get(h_topjets))->begin(), (&event.get(h_topjets))->end());
  else if (event.is_valid(h_jets)) jets.assign((&event.get(h_jets))->begin(), (&event.get(h_jets))->end());
  else throw logic_error("DeltaPhiJetMETCut: No valid jet collection given.");

  // Count how many jets fullfill Delta Phi(jet, MET) > minDeltaPhi
  int jetsCounter = 0;
  // Count how many jets have failed Delta Phi(jet, MET) > minDeltaPhi
  int jetsCounterFailing = 0;

  for(const auto & jet: jets){
    double Dphi = fabs(deltaPhi(jet, *event.met));
    if (Dphi > minDeltaPhi){jetsCounter ++;}
    else {jetsCounterFailing ++;}
  }

  // Pass if the number of jets fulfilling the criteria lies between minJets and maxJets.
  if (minJets <= jetsCounter && jetsCounter <=maxJets) return true;

  // Pass if maxJets is -1 (require all jets to pass) and number of jets fulfilling the criteria is at least minJets.
  if (minJets <= jetsCounter && maxJets == -1 && jetsCounterFailing==0) return true;

  return false;
}


ZprimeCandidateID::ZprimeCandidateID (const Event::Handle<vector<ZprimeCandidate> > & h_ZprimeCandidates_ ): h_ZprimeCandidates(h_ZprimeCandidates_){}

bool ZprimeCandidateID::operator()(const ZprimeCandidate& cand, const uhh2::Event & event) const {

  double min_chi2 = 10000;
  for(const auto & cand_: event.get(h_ZprimeCandidates)) min_chi2 = std::min(min_chi2, cand_.discriminator("chi2"));

  return cand.Zprime_mass()>=min_ZH_mass && cand.discriminator("btag_DeepCSV_loose")==0 && min_chi2==cand.discriminator("chi2");
}



PTMassCut::PTMassCut (float cut_min_, const Event::Handle<vector<ZprimeCandidate> > & h_ZprimeCandidates_ ): cut_min(cut_min_), h_ZprimeCandidates(h_ZprimeCandidates_){}

bool PTMassCut::passes(const Event& event){

  const auto & ZprimeCandidates = event.get(h_ZprimeCandidates);

  for(const auto & cand: ZprimeCandidates){
    auto Z = cand.Z();
    auto H = cand.H();
    auto ZPrime = Z.v4()+H.v4();
    float cut = Z.pt()/ZPrime.M();
    if( cut > cut_min) return true;
  }

  return false;
}



TaggerCut::TaggerCut (float cut_min_, float cut_max_, float pt_min_, string tagger_, const Event::Handle<vector<ZprimeCandidate> > & h_ZprimeCandidates_): cut_min(cut_min_), cut_max(cut_max_), pt_min(pt_min_), tagger(tagger_), h_ZprimeCandidates(h_ZprimeCandidates_){}

bool TaggerCut::passes(const Event& event){

  const std::vector<ZprimeCandidate> & ZprimeCandidates = event.get(h_ZprimeCandidates);
  if(ZprimeCandidates.size()!=1 ) throw std::runtime_error("TaggerCut::TaggerCut -- unexpected number of ZprimeCandidates");

  auto jet = ZprimeCandidates.at(0).H();
  if (jet.pt()>pt_min && pt_min>0) return true; //TODO pt_dep cut for all the cases?
  float tag = 0;
  if (tagger.find("CNN_")!=std::string::npos)    tag = jet.has_tag(jet.tagname2tag(tagger))? jet.get_tag(jet.tagname2tag(tagger)) : -1;
  else if (tagger.find("MC_")!=std::string::npos)tag = jet.get_tag(jet.tagname2tag(tagger));
  else if (tagger=="btag_DeepBoosted_H4qvsQCD")  tag = jet.btag_DeepBoosted_H4qvsQCD();
  else if (tagger=="btag_DeepBoosted_HbbvsQCD")  tag = jet.btag_DeepBoosted_HbbvsQCD();
  else if (tagger=="btag_DeepBoosted_probHbb")   tag = jet.btag_DeepBoosted_probHbb();
  else if (tagger=="btag_DeepBoosted_probHqqqq") tag = jet.btag_DeepBoosted_probHqqqq();
  else if (tagger=="tau21") tag = (jet.tau1()!=0) ? (jet.tau2()/jet.tau1()) : -1;
  else if (tagger=="tau31") tag = (jet.tau1()!=0) ? (jet.tau3()/jet.tau1()) : -1;
  else if (tagger=="tau41") tag = (jet.tau1()!=0) ? (jet.tau4()/jet.tau1()) : -1;
  else if (tagger=="tau32") tag = (jet.tau2()!=0) ? (jet.tau3()/jet.tau2()) : -1;
  else if (tagger=="tau42") tag = (jet.tau2()!=0) ? (jet.tau4()/jet.tau2()) : -1;
  else if (tagger=="tau43") tag = (jet.tau3()!=0) ? (jet.tau4()/jet.tau3()) : -1;
  else if (tagger=="tau21_groomed") tag = (jet.tau1_groomed()!=0 && jet.tau2_groomed()>0) ? (jet.tau2_groomed()/jet.tau1_groomed()) : -1;
  else if (tagger=="tau31_groomed") tag = (jet.tau1_groomed()!=0 && jet.tau3_groomed()>0) ? (jet.tau3_groomed()/jet.tau1_groomed()) : -1;
  else if (tagger=="tau41_groomed") tag = (jet.tau1_groomed()!=0 && jet.tau4_groomed()>0) ? (jet.tau4_groomed()/jet.tau1_groomed()) : -1;
  else if (tagger=="tau32_groomed") tag = (jet.tau2_groomed()!=0 && jet.tau3_groomed()>0) ? (jet.tau3_groomed()/jet.tau2_groomed()) : -1;
  else if (tagger=="tau42_groomed") tag = (jet.tau2_groomed()!=0 && jet.tau4_groomed()>0) ? (jet.tau4_groomed()/jet.tau2_groomed()) : -1;
  else if (tagger=="tau43_groomed") tag = (jet.tau3_groomed()!=0 && jet.tau4_groomed()>0) ? (jet.tau4_groomed()/jet.tau3_groomed()) : -1;
  else throw std::runtime_error("TaggerCut::passes -- tagger not defined!");

  if (pt_min == MassDepentdentCut_value) return (cut_min<cut_max)? tag>MassDepentdentCut(jet.pt()): tag<=MassDepentdentCut(jet.pt()); //TODO not too nice

  return (cut_min < tag  && tag < cut_max);

}


MultiBTagSubJetID::MultiBTagSubJetID(JetId SubjetId, int nsubjets): m_SubjetId(SubjetId), m_nsubjets(nsubjets){}

bool MultiBTagSubJetID::operator()(const TopJet & topjet, const uhh2::Event & event) const {

  int pass_SubjetId = 0;
  for(auto subjet : topjet.subjets()){ if(m_SubjetId(subjet,event)) pass_SubjetId++; }

  return pass_SubjetId >= m_nsubjets;
}



BlindDataSelection::BlindDataSelection(Context& ctx) {
  is_BlindSelection = ctx.get("dataset_type") == "true";
  is_mc = ctx.get("dataset_type") == "MC";
  h_is_Blind = ctx.get_handle<bool>("is_Blind");
}
bool BlindDataSelection::passes(const Event & event){

  if(is_mc || !is_BlindSelection) return true;
  return event.get(h_is_Blind);
}
