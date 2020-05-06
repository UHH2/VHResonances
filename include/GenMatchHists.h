#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/LorentzVector.h"

#include "UHH2/VHResonances/include/constants.hpp"
#include "UHH2/VHResonances/include/HistsBase.hpp"
#include "TString.h"

#include <vector>
#include <string>


class GenMatchHists: public HistsBase {
public:
  GenMatchHists(uhh2::Context&, const std::string&);
  virtual void fill(const uhh2::Event&) override;
protected:
  void book_ParticleHist(const std::string &, const std::string &, double, double);
  void book_FamilyTree(const std::string&);
  void fill_ParticleHist(const uhh2::Event &, const std::string &, const GenParticle &);
  void fill_FamilyTree(const uhh2::Event &, const std::string &, const GenParticle &);

  bool is_mc;
};
