#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/LorentzVector.h"

#include "UHH2/VHResonances/include/HistsBase.hpp"
#include "UHH2/VHResonances/include/ExtJetHists.h"

#include <vector>
#include <string>


class ExtJetHistsWithConditions: public ExtJetHists {
public:
  ExtJetHistsWithConditions(uhh2::Context& ctx, const std::string& dname, const std::string& collection, std::string condMatch_="", std::string condMatchStatus_="", const unsigned int NumberOfPlottedJets=4): ExtJetHists(ctx, dname, collection, NumberOfPlottedJets), condMatch(condMatch_), condMatchStatus(condMatchStatus_) {};
  virtual void fill(const uhh2::Event&) override;
protected:
  template <typename T>
  void fill_internal(const uhh2::Event&, std::vector<T>);
  std::string condMatch;
  std::string condMatchStatus;
};
