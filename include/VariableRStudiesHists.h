#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/LorentzVector.h"

#include "UHH2/VHResonances/include/constants.hpp"
#include "UHH2/VHResonances/include/HistsBase.hpp"
#include "TString.h"

#include <vector>
#include <string>


class VariableRStudiesHists: public HistsBase {
public:
  VariableRStudiesHists(uhh2::Context& , const std::string& dname, const std::string& collection_);
  virtual void fill(const uhh2::Event&) override;
protected:
  void book_ParticleHist(const std::string &, const std::string &, double, double);
  void fill_ParticleHist(const uhh2::Event &, const std::string &, const GenParticle &);
  template <typename T>
  void fill_internal(const uhh2::Event&, std::vector<T>);

  bool is_mc;
  std::string collection;

  uhh2::Event::Handle<std::vector<TopJet> > h_topjets;
};
