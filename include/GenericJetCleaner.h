#pragma once

#include <iostream>
#include <memory>
#include <unordered_map>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Utils.h"

#include "UHH2/common/include/JetCorrections.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/Utils.h"


class GenericJetCleaner: public uhh2::AnalysisModule {

public:
  explicit GenericJetCleaner(uhh2::Context & ctx, const std::string& jetLabel_, const bool& isTopJet_, const JetId& jetID_, const TopJetId& TopjetID_, const MuonId& muoID_, const ElectronId& eleID_);
  virtual bool process(uhh2::Event & event) override;
  virtual ~GenericJetCleaner();

  void PrintInputs();

protected:
  std::unordered_map<std::string, std::vector<std::string> > JEC_corr;
  // std::unordered_map<std::string, std::unique_ptr<JetCorrector> > JetCorrCHS;
  std::unordered_map<std::string, std::unique_ptr<GenericJetCorrector> > JetCorr;
  std::unordered_map<std::string, std::unique_ptr<GenericTopJetCorrector> > TopJetCorr;
  std::unordered_map<std::string, std::unique_ptr<JetLeptonCleaner_by_KEYmatching> > JLC;
  std::unordered_map<std::string, std::unique_ptr<GenericSubJetCorrector> > TopJetSubjetCorr;
  std::unique_ptr<GenericJetResolutionSmearer> JetResolutionSmearer;
  std::unique_ptr<AnalysisModule> SoftDropMassUpdate;
  std::unique_ptr<JetCleaner> Jetcleaner;
  std::unique_ptr<TopJetCleaner> TopJetcleaner;

  std::unordered_map<std::string, std::string> strings;
  std::unordered_map<std::string, bool> bools;

  uhh2::Event::Handle<std::vector<Jet> > h_jets;

  JetId jetID;
  TopJetId TopjetID;

  MuonId muoIDcleaning;
  ElectronId eleIDcleaning;

  std::unordered_map<std::string, std::vector<std::string>> runs = { {"2016", runPeriods2016}, {"2017", runPeriods2017}, {"2018", runPeriods2018}};

};
