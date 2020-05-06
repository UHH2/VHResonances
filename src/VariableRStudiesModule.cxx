#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Selection.h"
#include "UHH2/core/include/Utils.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/common/include/EventHists.h"
#include "UHH2/common/include/ElectronHists.h"
#include "UHH2/common/include/MuonHists.h"
#include "UHH2/common/include/MuonIds.h"
#include "UHH2/common/include/ElectronIds.h"
#include "UHH2/common/include/ObjectIdUtils.h"
#include "UHH2/common/include/JetIds.h"
#include "UHH2/common/include/JetCorrections.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/CommonModules.h"
#include "UHH2/common/include/JetCorrections.h"
#include "UHH2/common/include/MCWeight.h"
#include "UHH2/common/include/EventVariables.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/EventVariables.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/AdditionalSelections.h"
#include "UHH2/common/include/LumiSelection.h"
#include "UHH2/common/include/TriggerSelection.h"
#include "UHH2/common/include/Utils.h"

#include "UHH2/HiggsToWWTagger/include/ModuleBase.h"
#include "UHH2/HiggsToWWTagger/include/constants.hpp"
#include "UHH2/HiggsToWWTagger/include/ExtJetHists.h"
#include "UHH2/HiggsToWWTagger/include/ExtJetHistsWithConditions.h"
#include "UHH2/HiggsToWWTagger/include/GenMatchHists.h"
#include "UHH2/HiggsToWWTagger/include/VariableRStudiesHists.h"
#include "UHH2/HiggsToWWTagger/include/GenLevelJetMatch.h"
#include "UHH2/HiggsToWWTagger/include/PrintingFunctions.hpp"

#include "UHH2/HiggsToWWTagger/include/GenericJetCleaner.h"

using namespace std;

/*
█ ██████  ███████ ███████ ██ ███    ██ ██ ████████ ██  ██████  ███    ██
█ ██   ██ ██      ██      ██ ████   ██ ██    ██    ██ ██    ██ ████   ██
█ ██   ██ █████   █████   ██ ██ ██  ██ ██    ██    ██ ██    ██ ██ ██  ██
█ ██   ██ ██      ██      ██ ██  ██ ██ ██    ██    ██ ██    ██ ██  ██ ██
█ ██████  ███████ ██      ██ ██   ████ ██    ██    ██  ██████  ██   ████
*/

class VariableRStudiesModule: public ModuleBASE {

public:

  explicit VariableRStudiesModule(uhh2::Context&);
  virtual bool process(uhh2::Event&) override;
  void book_histograms(uhh2::Context&);
  void fill_histograms(uhh2::Event&, string);

protected:

  std::vector<std::string> histogram_tags = { "VariableRStudies"};

};


void VariableRStudiesModule::book_histograms(uhh2::Context& ctx){
  for(const auto & tag : histogram_tags){
    string mytag;
    mytag = "gen_"            + tag; book_HFolder(mytag, new GenMatchHists(ctx,mytag));
    mytag = "gen_matchPuppi_" + tag; book_HFolder(mytag, new VariableRStudiesHists(ctx,mytag,"toppuppijets"));
    mytag = "gen_matchCHS_"   + tag; book_HFolder(mytag, new VariableRStudiesHists(ctx,mytag,"topjets"));

  }
}

void VariableRStudiesModule::fill_histograms(uhh2::Event& event, string tag){
  std::vector<string> mytags = {"gen_", "gen_matchPuppi_", "gen_matchCHS_"};
  for (auto& mytag : mytags) {
    HFolder(mytag+ tag)->fill(event);
  }
}

/*
█  ██████  ██████  ███    ██ ███████ ████████ ██████  ██    ██  ██████ ████████  ██████  ██████
█ ██      ██    ██ ████   ██ ██         ██    ██   ██ ██    ██ ██         ██    ██    ██ ██   ██
█ ██      ██    ██ ██ ██  ██ ███████    ██    ██████  ██    ██ ██         ██    ██    ██ ██████
█ ██      ██    ██ ██  ██ ██      ██    ██    ██   ██ ██    ██ ██         ██    ██    ██ ██   ██
█  ██████  ██████  ██   ████ ███████    ██    ██   ██  ██████   ██████    ██     ██████  ██   ██
*/

VariableRStudiesModule::VariableRStudiesModule(uhh2::Context& ctx){

  book_histograms(ctx);


}


/*
█ ██████  ██████   ██████   ██████ ███████ ███████ ███████
█ ██   ██ ██   ██ ██    ██ ██      ██      ██      ██
█ ██████  ██████  ██    ██ ██      █████   ███████ ███████
█ ██      ██   ██ ██    ██ ██      ██           ██      ██
█ ██      ██   ██  ██████   ██████ ███████ ███████ ███████
*/

bool VariableRStudiesModule::process(uhh2::Event& event) {

  fill_histograms(event, "VariableRStudies");

  return true;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the VariableRStudiesModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(VariableRStudiesModule)
