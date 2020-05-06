#include "UHH2/core/include/Event.h"

#include "UHH2/common/include/Utils.h"

#include "UHH2/VHResonances/include/HistsBase.hpp"
#include "UHH2/VHResonances/include/ExtJetHists.h"
#include "UHH2/VHResonances/include/ExtJetHistsWithConditions.h"
#include "UHH2/VHResonances/include/GenLevelJetMatch.h"
#include "UHH2/VHResonances/include/Utils.hpp"

using namespace std;
using namespace uhh2;


/*
█ ██████   ██████   ██████  ██   ██      █████  ███    ██ ██████      ███████ ██ ██      ██
█ ██   ██ ██    ██ ██    ██ ██  ██      ██   ██ ████   ██ ██   ██     ██      ██ ██      ██
█ ██████  ██    ██ ██    ██ █████       ███████ ██ ██  ██ ██   ██     █████   ██ ██      ██
█ ██   ██ ██    ██ ██    ██ ██  ██      ██   ██ ██  ██ ██ ██   ██     ██      ██ ██      ██
█ ██████   ██████   ██████  ██   ██     ██   ██ ██   ████ ██████      ██      ██ ███████ ███████
*/


void ExtJetHistsWithConditions::fill(const Event & event){
  // if (event.is_valid(h_jets)) fill_internal(event, event.get(h_jets));
  // else
  if (event.is_valid(h_topjets)) fill_internal(event, event.get(h_topjets));
}

template <typename T>
void ExtJetHistsWithConditions::fill_internal(const Event & event, vector<T> jets){
  auto weight = event.weight;
  fill_H1("sum_event_weights", 1., weight);
  fill_H1("weights", weight, 1);
  fill_H1("number", jets.size(), weight);
  for(unsigned int i = 0; i <jets.size(); i++){
    T& jet = jets[i];
    if (condMatch!=""       && !jet.has_tag(TopJet::Matching)) continue;
    if (condMatchStatus!="" && !jet.has_tag(TopJet::MatchingStatus)) continue;
    if (condMatch!=""       && FloatToMatching(jet)!=StringToMatching(condMatch) ) continue;
    if (condMatchStatus!="" && FloatToMatchingStatus(jet)!=StringToMatchingStatus(condMatchStatus) ) continue;

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
}
