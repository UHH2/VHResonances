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

  std::string condMatch;
  std::string condMatchStatus;
  std::string massType;
  std::string massPlotName;

  std::vector<std::string> discriminators= {"chi2", "SDmass", "tau1", "tau2", "tau3", "tau4", "tau21", "tau31", "tau41", "tau32", "tau42", "tau43",
  "NN_HWW", "NN_Hbb", "NN_QCD", "NN_Top", "NN_W", "NN_Z", "NN_HWW_1", "NN_Hbb_1", "NN_QCD_1", "NN_Top_1", "NN_W_1", "NN_Z_1", "NN_HWW_2", "NN_Hbb_2", "NN_QCD_2", "NN_Top_2", "NN_W_2", "NN_Z_2",
  "btag_DeepBoosted_H4qvsQCD", "btag_MassDecorrelatedDeepBoosted_H4qvsQCD", "btag_DeepBoosted_probHqqqq", "btag_DeepBoosted_raw_score_h", "btag_MassDecorrelatedDeepBoosted_probHqqqq",
  "CNN_HWW", "CNN_QCD", "CNN_Top", "CNN_W", "CNN_Z", "DCL_HWW", "DCL_QCD", "DCL_Top", "DCL_W", "DCL_Z"};
};


class HiggsToWWHistsSlim: public HistsBase {
public:

  HiggsToWWHistsSlim(uhh2::Context& ctx, const std::string& dname, const std::string& condMatch_="", const std::string& condMatchStatus_="");
  virtual void fill(const uhh2::Event&) override;
  virtual ~HiggsToWWHistsSlim();

private:
  uhh2::Event::Handle< std::vector<ZprimeCandidate> > h_ZprimeCandidates;
  std::string condMatch;
  std::string condMatchStatus;
  std::string massType;
  std::string massPlotName;
};
