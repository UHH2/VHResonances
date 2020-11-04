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

  std::vector<std::string> discriminators = {"chi2", "SDmass", "tau1", "tau2", "tau3", "tau4", "tau21", "tau31", "tau41", "tau32", "tau42", "tau43",
  "btag_DeepBoosted_H4qvsQCD", "btag_MassDecorrelatedDeepBoosted_H4qvsQCD", "btag_DeepBoosted_probHqqqq",
  "btag_DeepBoosted_raw_score_h", "btag_MassDecorrelatedDeepBoosted_probHqqqq" };


  std::vector<std::string> discriminators_Extra = {"btag_BoostedDoubleSecondaryVertexAK8", "btag_BoostedDoubleSecondaryVertexCA15",
  "btag_DeepDoubleBvLJet_probHbb", "btag_DeepDoubleBvLJet_probQCD", "btag_DeepDoubleCvBJet_probHbb", "btag_DeepDoubleCvBJet_probHcc", "btag_DeepDoubleCvLJet_probHcc","btag_DeepDoubleCvLJet_probQCD",
  "btag_DeepBoosted_WvsQCD", "btag_DeepBoosted_ZvsQCD", "btag_DeepBoosted_ZbbvsQCD", "btag_DeepBoosted_HbbvsQCD",
  "btag_DeepBoosted_probQCDb", "btag_DeepBoosted_probQCDbb", "btag_DeepBoosted_probQCDc", "btag_DeepBoosted_probQCDcc", "btag_DeepBoosted_probQCDothers",
  "btag_DeepBoosted_probWqq", "btag_DeepBoosted_probWcq", "btag_DeepBoosted_probZcc", "btag_DeepBoosted_probZqq", "btag_DeepBoosted_probZbb",
  "btag_DeepBoosted_probHbb", "btag_DeepBoosted_probHcc", "btag_DeepBoosted_raw_score_qcd", "btag_DeepBoosted_raw_score_w", "btag_DeepBoosted_raw_score_z",
  "btag_DeepBoosted_HbbvsHcc", "btag_DeepBoosted_HccvsQCD", "btag_DeepBoosted_HvsQCD" };

  std::vector<std::string> discriminators_subjets =  {"btag_combinedSecondaryVertex", "btag_combinedSecondaryVertexMVA", "btag_DeepJet",
  "btag_DeepCSV", "btag_DeepFlavour_bb", "btag_DeepFlavour_b", "btag_DeepFlavour_lepb", "btag_DeepFlavour_uds", "btag_DeepFlavour_g", "btag_DeepFlavour_c"};

  //
  // std::vector<std::string> discriminators_MyNN = { "NN_HWW", "NN_Hbb", "NN_QCD", "NN_Top", "NN_W", "NN_Z",
  // "NN_HWW_1", "NN_Hbb_1", "NN_QCD_1", "NN_Top_1", "NN_W_1", "NN_Z_1",
  // "NN_HWW_2", "NN_Hbb_2", "NN_QCD_2", "NN_Top_2", "NN_W_2", "NN_Z_2",
  // "CNN_HWW", "CNN_QCD", "CNN_Top", "CNN_W", "CNN_Z", "DCL_HWW", "DCL_QCD", "DCL_Top", "DCL_W", "DCL_Z",
  // };
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
