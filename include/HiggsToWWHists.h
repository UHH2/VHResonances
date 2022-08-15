#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/VHResonances/include/HistsBase.hpp"
#include "UHH2/VHResonances/include/Utils.hpp"
#include "UHH2/VHResonances/include/HiggsToWWSelection.h"

class HiggsToWWHists: public HistsBase {
public:

  HiggsToWWHists(uhh2::Context& ctx, const std::string& dname, const std::string& condMatch_="", const std::string& condMatchStatus_="");
  virtual void fill(const uhh2::Event&) override;
  virtual ~HiggsToWWHists();

private:
  uhh2::Event::Handle< std::vector<ZprimeCandidate> > h_ZprimeCandidates;
  uhh2::Event::Handle<float> h_HDecay;
  uhh2::Event::Handle<float> h_ZDecay;
  uhh2::Event::Handle<float> h_ZprimeDecay;

  bool isInvisible;
  std::string condMatch;
  std::string condMatchStatus;
  std::string massType;
  std::string massPlotName;

  std::vector<std::string> discriminators = {"chi2", "SDmass", "btag_ParticleNet_mass",
  "tau1", "tau2", "tau3", "tau4", "tau21", "tau31", "tau41", "tau32", "tau42", "tau43",
  "btag_DeepBoosted_HccvsQCD_MD", "btag_DeepBoosted_HccvsQCD", "btag_DeepBoosted_ZHccvsQCD_MD", "btag_DeepBoosted_ZHccvsQCD",
  "btag_DeepBoosted_H4qvsQCD_MD", "btag_DeepBoosted_H4qvsQCD",
  "btag_ParticleNet_HccvsQCD","btag_ParticleNet_HccvsQCD_MD","btag_ParticleNet_ZHccvsQCD",
  "btag_ParticleNet_ZHccvsQCD_MD","btag_ParticleNet_H4qvsQCD","btag_ParticleNet_HbbvsQCD",
  "btag_ParticleNet_HccvsQCD_2","btag_ParticleNet_XbbvsQCD_MD","btag_ParticleNet_XccvsQCD_MD"};

  std::vector<std::string> discriminators_subjets =  {"btag_DeepJet",
  "btag_DeepCSV", "btag_DeepFlavour_bb", "btag_DeepFlavour_b", "btag_DeepFlavour_lepb", "btag_DeepFlavour_uds", "btag_DeepFlavour_g", "btag_DeepFlavour_c"};

  std::vector<std::string> discriminators_Extra = {};

};


class HiggsToWWHistsSlim: public HistsBase {
public:

  HiggsToWWHistsSlim(uhh2::Context& ctx, const std::string& dname, const std::string& condMatch_="", const std::string& condMatchStatus_="");
  virtual void fill(const uhh2::Event&) override;
  virtual ~HiggsToWWHistsSlim();

private:
  uhh2::Event::Handle< std::vector<ZprimeCandidate> > h_ZprimeCandidates;
  bool isInvisible;
  std::string condMatch;
  std::string condMatchStatus;
  std::string massType;
  std::string massPlotName;
};
