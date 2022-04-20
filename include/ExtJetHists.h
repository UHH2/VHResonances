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

// #include <eigen3/Eigen/Dense>
// #include <UHH2/Eigen/Eigen/Dense>



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


  std::vector<std::string> discriminators = {"btag_MassDecorrelatedDeepBoosted_HccvsQCD", "btag_DeepBoosted_HccvsQCD",
  "btag_MassDecorrelatedDeepBoosted_ZHccvsQCD", "btag_DeepBoosted_ZHccvsQCD", "btag_MassDecorrelatedDeepBoosted_H4qvsQCD", "btag_DeepBoosted_H4qvsQCD"};

  std::vector<std::string> discriminators_subjets =  {"btag_DeepJet", "btag_DeepCSV",
  "btag_DeepFlavour_bb", "btag_DeepFlavour_b", "btag_DeepFlavour_lepb", "btag_DeepFlavour_uds", "btag_DeepFlavour_g", "btag_DeepFlavour_c"};

  std::vector<std::string> discriminators_Extra = {"btag_DeepBoosted_WvsQCD", "btag_DeepBoosted_ZvsQCD", "btag_DeepBoosted_HvsQCD",
  "btag_MassDecorrelatedDeepBoosted_WvsQCD", "btag_MassDecorrelatedDeepBoosted_ZvsQCD", "btag_MassDecorrelatedDeepBoosted_ZHbbvsQCD",
  "btag_DeepBoosted_HbbvsQCD", "btag_MassDecorrelatedDeepBoosted_HbbvsQCD", "btag_MassDecorrelatedDeepBoosted_bbvsLight", "btag_MassDecorrelatedDeepBoosted_ccvsLight"};

  // std::vector<std::string> discriminators_Extra = {
  //   "btag_DeepBoosted_TvsQCD", "", "btag_MassDecorrelatedDeepBoosted_TvsQCD", "btag_DeepBoosted_HbbvsHcc",
  //   "btag_BoostedDoubleSecondaryVertexAK8", "btag_BoostedDoubleSecondaryVertexCA15", "btag_DeepDoubleBvLJet_probHbb", "btag_DeepDoubleBvLJet_probQCD",
  //   "btag_DeepDoubleCvBJet_probHbb", "btag_DeepDoubleCvBJet_probHcc", "btag_DeepDoubleCvLJet_probHcc", "btag_DeepDoubleCvLJet_probQCD",
  //   "btag_MassIndependentDeepDoubleBvLJet_probHbb", "btag_MassIndependentDeepDoubleBvLJet_probQCD", "btag_MassIndependentDeepDoubleCvBJet_probHbb",
  //   "btag_MassIndependentDeepDoubleCvBJet_probHcc", "btag_MassIndependentDeepDoubleCvLJet_probHcc", "btag_MassIndependentDeepDoubleCvLJet_probQCD",
  //   "btag_DeepBoosted_probQCDb", "btag_DeepBoosted_probQCDbb", "btag_DeepBoosted_probQCDc", "btag_DeepBoosted_probQCDcc", "btag_DeepBoosted_probQCDothers",
  //   "btag_DeepBoosted_probTbqq", "btag_DeepBoosted_probTbcq", "btag_DeepBoosted_probTbq", "btag_DeepBoosted_probTbc", "btag_DeepBoosted_probWqq",
  //   "btag_DeepBoosted_raw_score_top", "btag_MassDecorrelatedDeepBoosted_probQCDc", "btag_MassDecorrelatedDeepBoosted_probQCDbb", "btag_MassDecorrelatedDeepBoosted_probTbqq",
  //   "btag_MassDecorrelatedDeepBoosted_probTbcq", "btag_MassDecorrelatedDeepBoosted_probTbq", "btag_MassDecorrelatedDeepBoosted_probQCDothers",
  //   "btag_MassDecorrelatedDeepBoosted_probQCDb", "btag_MassDecorrelatedDeepBoosted_probTbc", "btag_MassDecorrelatedDeepBoosted_probQCDcc",
  //   "btag_MassDecorrelatedDeepBoosted_probHcc", "btag_MassDecorrelatedDeepBoosted_probZcc", "btag_MassDecorrelatedDeepBoosted_proWcq",
  //   "btag_MassDecorrelatedDeepBoosted_probHbb", "btag_MassDecorrelatedDeepBoosted_probZbb", "btag_MassDecorrelatedDeepBoosted_probHqqqq",
  //   "btag_MassDecorrelatedDeepBoosted_probWqq", "btag_MassDecorrelatedDeepBoosted_probZqq", "btag_MassDecorrelatedDeepBoosted_ZbbvsQCD",
  //   "btag_DeepBoosted_ZbbvsQCD", "btag_DeepBoosted_probWcq", "btag_DeepBoosted_probZcc", "btag_DeepBoosted_probZqq",
  //   "btag_DeepBoosted_probZbb", "btag_DeepBoosted_probHbb", "btag_DeepBoosted_probHcc", "btag_DeepBoosted_probHqqqq",
  //   "btag_DeepBoosted_raw_score_qcd", "btag_DeepBoosted_raw_score_w", "btag_DeepBoosted_raw_score_z", "btag_DeepBoosted_raw_score_h",
  // };

};
