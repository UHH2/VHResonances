#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/LorentzVector.h"
#include "UHH2/common/include/JetIds.h"

#include "UHH2/VHResonances/include/constants.hpp"
#include "UHH2/VHResonances/include/HistsBase.hpp"
#include "TString.h"

#include <vector>
#include <string>


// save diagnostic state
#pragma GCC diagnostic push
// turn off the specific warning. Can also use "-Wall"
#pragma GCC diagnostic ignored "-Wignored-qualifiers"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"
#pragma GCC diagnostic pop


class ExtJetHists: public HistsBase {
public:
  ExtJetHists(uhh2::Context&, const std::string&, const std::string& ,const unsigned int NumberOfPlottedJets=4);
  virtual void fill(const uhh2::Event&) override;

protected:
  unsigned int NumberOfPlottedJets;
  std::vector<std::string> axis_suffix {"first jet","second jet","third jet","fourth jet"};
  bool isTop = false;
  void book_jetHist(const std::string&, const std::string&, double, double, bool);
  template <typename T>
  void fill_jetHist(const uhh2::Event&,const std::string&, const T&);
  template <typename T>
  void fill_internal(const uhh2::Event&, std::vector<T>);

  std::string collection;
  uhh2::Event::Handle<std::vector<Jet> > h_jets;
  uhh2::Event::Handle<std::vector<TopJet> > h_topjets;
  uhh2::Event::Handle<std::vector<tensorflow::Tensor> > h_image;
  std::map<TString, JetId> Btag_map;
};
