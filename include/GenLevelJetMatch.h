#pragma once

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/HiggsToWWTagger/include/constants.hpp"

class GenLevelJetMatch: public uhh2::AnalysisModule {
public:
  GenLevelJetMatch(uhh2::Context & ctx, const std::string & jetCollection="topjets");
  GenLevelJetMatch(float Dr_);
  virtual bool process(uhh2::Event & event) override;

  bool MatchGenPart(const uhh2::Event & event, TopJet& jet, double Dr);

private:
  float Dr;
  uhh2::Event::Handle<std::vector<TopJet>> h_topjets_;
};
