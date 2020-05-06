#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/HiggsToWWTagger/include/HistsBase.hpp"
#include "UHH2/HiggsToWWTagger/include/constants.hpp"

class DiLeptonHists: public HistsBase {
public:

	DiLeptonHists(uhh2::Context&, const std::string&, const std::string & jetcollection ="", const std::string & topjetcollection="");
  virtual void fill(const uhh2::Event&) override;
	virtual ~DiLeptonHists();

private:
	uhh2::Event::Handle<std::vector<Jet> > h_jets;
  uhh2::Event::Handle<std::vector<TopJet> > h_topjets;

};
