#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/VHResonances/include/HistsBase.hpp"
#include "UHH2/VHResonances/include/constants.hpp"

class DiLeptonHists: public HistsBase {
public:

  DiLeptonHists(uhh2::Context&, const std::string&, const std::string & jetcollection ="", const std::string & topjetcollection="");
  virtual void fill(const uhh2::Event&) override;
  virtual ~DiLeptonHists();

private:
  const std::vector<std::string> leptonNames = {"Electron", "Muon"};
  uhh2::Event::Handle<std::vector<Jet> > h_jets;
  uhh2::Event::Handle<std::vector<TopJet> > h_topjets;

};
