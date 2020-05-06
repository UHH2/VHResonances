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
#include "UHH2/HiggsToWWTagger/include/GenLevelJetMatch.h"
#include "UHH2/HiggsToWWTagger/include/HiggsToWWSelections.h"
#include "UHH2/HiggsToWWTagger/include/DiLeptonHists.h"
#include "UHH2/HiggsToWWTagger/include/HiggsToWWModules.h"
#include "UHH2/HiggsToWWTagger/include/HiggsToWWHists.h"

using namespace std;

/*
█ ██████  ███████ ███████ ██ ███    ██ ██ ████████ ██  ██████  ███    ██
█ ██   ██ ██      ██      ██ ████   ██ ██    ██    ██ ██    ██ ████   ██
█ ██   ██ █████   █████   ██ ██ ██  ██ ██    ██    ██ ██    ██ ██ ██  ██
█ ██   ██ ██      ██      ██ ██  ██ ██ ██    ██    ██ ██    ██ ██  ██ ██
█ ██████  ███████ ██      ██ ██   ████ ██    ██    ██  ██████  ██   ████
*/

class SignalRegionModule: public ModuleBASE {

public:

  explicit SignalRegionModule(uhh2::Context&);
  virtual bool process(uhh2::Event&) override;
  void book_histograms(uhh2::Context&);
  void fill_histograms(uhh2::Event&, string);
  void book_handles(uhh2::Context&);
  void export_weights(uhh2::Event&);
  void PrintInputs();

protected:

  // Define variables

  std::vector<std::string> histogram_tags = {"Selection", "FSSelection",
  "btag_DeepBoosted_H4qvsQCD_SR",           "btag_DeepBoosted_H4qvsQCD_CR",
  "btag_DeepBoosted_H4qvsQCDp2_SR",         "btag_DeepBoosted_H4qvsQCDp2_CR",
  "btag_DeepBoosted_H4qvsQCDp02_SR",        "btag_DeepBoosted_H4qvsQCDp02_CR",
  "btag_DeepBoosted_H4qvsQCDpt1000_SR",     "btag_DeepBoosted_H4qvsQCDpt1000_CR",
  "btag_DeepBoosted_H4qvsQCDpt1000p2_SR",   "btag_DeepBoosted_H4qvsQCDpt1000p2_CR",
  "btag_DeepBoosted_H4qvsQCDpt1000p02_SR",  "btag_DeepBoosted_H4qvsQCDpt1000p02_CR",
  "btag_DeepBoosted_H4qvsQCDpt1500_SR",     "btag_DeepBoosted_H4qvsQCDpt1500_CR",
  "btag_DeepBoosted_H4qvsQCDpt1500p2_SR",   "btag_DeepBoosted_H4qvsQCDpt1500p2_CR",
  "btag_DeepBoosted_H4qvsQCDpt1500p02_SR",  "btag_DeepBoosted_H4qvsQCDpt1500p02_CR",

  "btag_DeepBoosted_HbbvsQCD_SR",   "btag_DeepBoosted_HbbvsQCD_CR",
  // "btag_DeepBoosted_probHbb_SR",    "btag_DeepBoosted_probHbb_CR",
  // "btag_DeepBoosted_probHqqqq_SR",  "btag_DeepBoosted_probHqqqq_CR",
  //"NN_SR", "NN_CR", "NN_1_SR", "NN_1_CR", "NN_2_SR", "NN_2_CR",  "CNN_SR", "CNN_CR",
  // "tau21_SR", "tau21_CR", "tau31_SR", "tau31_CR", "tau41_SR", "tau41_CR", "tau32_SR", "tau32_CR", "tau42_SR", "tau42_CR", "tau43_SR", "tau43_CR",
  "tau42_SR", "tau42_CR"};

  std::vector<std::string> weight_tags = {"weight_lumi", "weight_GLP", "weight_pu", "weight_pu_up", "weight_pu_down", "HDecay", "ZDecay", "ZprimeDecay"};


  std::unordered_map<std::string, std::string> MS;
  std::unordered_map<std::string, bool> MB;

  Event::Handle<std::vector<TopJet> > h_topjets;
  Event::Handle<std::vector<ZprimeCandidate> > h_ZprimeCandidates;

  // Define selections

  // std::unique_ptr<Selection> NN_SR_selection, NN_CR_selection;
  // std::unique_ptr<Selection> NN_SR_selection_1, NN_CR_selection_1;
  // std::unique_ptr<Selection> NN_SR_selection_2, NN_CR_selection_2;
  // std::unique_ptr<Selection> CNN_SR_selection, CNN_CR_selection;


  std::unique_ptr<Selection> btag_DeepBoosted_H4qvsQCD_SR_selection,  btag_DeepBoosted_H4qvsQCD_CR_selection;

  std::unique_ptr<Selection> btag_DeepBoosted_H4qvsQCDp2_SR_selection,        btag_DeepBoosted_H4qvsQCDp2_CR_selection;
  std::unique_ptr<Selection> btag_DeepBoosted_H4qvsQCDp02_SR_selection,       btag_DeepBoosted_H4qvsQCDp02_CR_selection;
  std::unique_ptr<Selection> btag_DeepBoosted_H4qvsQCDpt1000_SR_selection,    btag_DeepBoosted_H4qvsQCDpt1000_CR_selection;
  std::unique_ptr<Selection> btag_DeepBoosted_H4qvsQCDpt1000p2_SR_selection,  btag_DeepBoosted_H4qvsQCDpt1000p2_CR_selection;
  std::unique_ptr<Selection> btag_DeepBoosted_H4qvsQCDpt1000p02_SR_selection, btag_DeepBoosted_H4qvsQCDpt1000p02_CR_selection;
  std::unique_ptr<Selection> btag_DeepBoosted_H4qvsQCDpt1500_SR_selection,    btag_DeepBoosted_H4qvsQCDpt1500_CR_selection;
  std::unique_ptr<Selection> btag_DeepBoosted_H4qvsQCDpt1500p2_SR_selection,  btag_DeepBoosted_H4qvsQCDpt1500p2_CR_selection;
  std::unique_ptr<Selection> btag_DeepBoosted_H4qvsQCDpt1500p02_SR_selection, btag_DeepBoosted_H4qvsQCDpt1500p02_CR_selection;

  std::unique_ptr<Selection> btag_DeepBoosted_HbbvsQCD_SR_selection,  btag_DeepBoosted_HbbvsQCD_CR_selection;
  // std::unique_ptr<Selection> btag_DeepBoosted_probHbb_SR_selection,   btag_DeepBoosted_probHbb_CR_selection;
  // std::unique_ptr<Selection> btag_DeepBoosted_probHqqqq_SR_selection, btag_DeepBoosted_probHqqqq_CR_selection;

  // std::unique_ptr<Selection> tau21_SR_selection, tau21_CR_selection;
  // std::unique_ptr<Selection> tau31_SR_selection, tau31_CR_selection;
  // std::unique_ptr<Selection> tau41_SR_selection, tau41_CR_selection;
  // std::unique_ptr<Selection> tau32_SR_selection, tau32_CR_selection;
  std::unique_ptr<Selection> tau42_SR_selection, tau42_CR_selection;
  // std::unique_ptr<Selection> tau43_SR_selection, tau43_CR_selection;

};


void SignalRegionModule::PrintInputs() {
  std::cout << "****************************************" << std::endl;
  std::cout << "           SignalRegionModule           " << std::endl;
  std::cout << "----------------------------------------" << std::endl;
  for (auto x : MS) std::cout << x.first << std::string( 18-x.first.size(), ' ' ) << x.second << '\n';
  for (auto x : MB)   std::cout << x.first << std::string( 18-x.first.size(), ' ' ) << (x.second? "true" : "false") << '\n';
  std::cout << "****************************************\n" << std::endl;
}

void SignalRegionModule::book_handles(uhh2::Context& ctx) {
  string tag;
  for(const auto & tag : weight_tags) {
    if (!MB["is_mc"] && tag.find("weight_pu")!=std::string::npos) continue;
    book_WFolder(tag+"_in", new Event::Handle< float >, ctx.declare_event_output< float >(tag));
    book_WFolder(tag+"_out", new Event::Handle< float >, ctx.declare_event_input< float >(tag));
  }
}


void SignalRegionModule::export_weights(uhh2::Event& event) {
  for(const auto & tag : weight_tags) {
    if (!MB["is_mc"] && tag.find("weight_pu")!=std::string::npos) continue;
    event.set(WFolder(tag+"_out"), event.get(WFolder(tag+"_in")));
  }
}


void SignalRegionModule::book_histograms(uhh2::Context& ctx){
  for(const auto & tag : histogram_tags){
    string mytag;
    mytag = "nTopJet_"                + tag; book_HFolder(mytag, new ExtJetHists(ctx,mytag, MS["topjetLabel"] ));
    mytag = "ZprimeCandidate_"        + tag; book_HFolder(mytag, new HiggsToWWHists(ctx,mytag));
    mytag = "ZprimeCandidate_HWW"     + tag; book_HFolder(mytag, new HiggsToWWHists(ctx,mytag, "HWWMatch"));
    mytag = "ZprimeCandidate_Hbb"     + tag; book_HFolder(mytag, new HiggsToWWHists(ctx,mytag, "HbbMatch"));
    mytag = "ZprimeCandidate_HZZ"     + tag; book_HFolder(mytag, new HiggsToWWHists(ctx,mytag, "HZZMatch"));
    mytag = "ZprimeCandidate_else"    + tag; book_HFolder(mytag, new HiggsToWWHists(ctx,mytag, "else"));
    mytag = "ZprimeCandidate_PUup"    + tag; book_HFolder(mytag, new HiggsToWWHists(ctx,mytag));
    mytag = "ZprimeCandidate_PUdown"  + tag; book_HFolder(mytag, new HiggsToWWHists(ctx,mytag));
  }
}

void SignalRegionModule::fill_histograms(uhh2::Event& event, string tag){
  std::vector<string> mytags = {"nTopJet_", "ZprimeCandidate_", "ZprimeCandidate_HWW", "ZprimeCandidate_Hbb", "ZprimeCandidate_HZZ", "ZprimeCandidate_else"};
  for (auto& mytag : mytags) HFolder(mytag+ tag)->fill(event);
  string mytag;
  float save_weight = event.weight;
  if (!event.isRealData) event.weight = save_weight*event.get(WFolder("weight_pu_up_in"))/event.get(WFolder("weight_pu_in"));
  mytag = "ZprimeCandidate_PUup";   HFolder(mytag+ tag)->fill(event);
  if (!event.isRealData) event.weight = save_weight*event.get(WFolder("weight_pu_down_in"))/event.get(WFolder("weight_pu_in"));
  mytag = "ZprimeCandidate_PUdown"; HFolder(mytag+ tag)->fill(event);
  event.weight = save_weight;
}

/*
█  ██████  ██████  ███    ██ ███████ ████████ ██████  ██    ██  ██████ ████████  ██████  ██████
█ ██      ██    ██ ████   ██ ██         ██    ██   ██ ██    ██ ██         ██    ██    ██ ██   ██
█ ██      ██    ██ ██ ██  ██ ███████    ██    ██████  ██    ██ ██         ██    ██    ██ ██████
█ ██      ██    ██ ██  ██ ██      ██    ██    ██   ██ ██    ██ ██         ██    ██    ██ ██   ██
█  ██████  ██████  ██   ████ ███████    ██    ██   ██  ██████   ██████    ██     ██████  ██   ██
*/

SignalRegionModule::SignalRegionModule(uhh2::Context& ctx){

  // Set up variables

  MS["dataset_version"] = ctx.get("dataset_version");
  MB["is_mc"]           = ctx.get("dataset_type") == "MC";
  MB["isPuppi"]         = string2bool(ctx.get("isPuppi"));
  MB["isCHS"]           = string2bool(ctx.get("isCHS"));
  MB["isHOTVR"]         = string2bool(ctx.get("isHOTVR"));
  MB["muonchannel"]     = string2bool(ctx.get("muonchannel"));
  MB["electronchannel"] = string2bool(ctx.get("electronchannel"));

  MB["is_ZH"]   = MS["dataset_version"].find("MC_ZprimeToZH")!=std::string::npos;
  MB["is_bb"]   = MS["dataset_version"].find("Tobb")!=std::string::npos;
  MB["is_WW"]   = MS["dataset_version"].find("ToWW")!=std::string::npos;
  MB["is_else"] = MS["dataset_version"].find("_extra")!=std::string::npos;


  if (MB["isPuppi"] == MB["isCHS"] && MB["isPuppi"] == MB["isHOTVR"]) throw std::runtime_error("In SignalRegionModule.cxx: Choose exactly one jet collection.");

  MS["topjetLabel"] = MB["isCHS"]? "topjets": (MB["isPuppi"]? "toppuppijets": (MB["isHOTVR"]? "hotvrPuppi": ""));

  // Set up histograms:

  book_histograms(ctx);
  book_handles(ctx);
  PrintInputs();

  // Set up selections

  h_topjets = ctx.get_handle<std::vector<TopJet>>(MS["topjetLabel"]);
  h_ZprimeCandidates = ctx.declare_event_input<std::vector<ZprimeCandidate>>("ZprimeCandidate");
  // TODO get_handle? declare_event_output

  // NN_SR_selection.reset(new TaggerCut(0.2, 1.0, "MC_HWW", h_ZprimeCandidates));
  // NN_CR_selection.reset(new TaggerCut(0.0, 0.2, "MC_HWW", h_ZprimeCandidates));
  // NN_SR_selection_1.reset(new TaggerCut(0.8, 1.0, "MC_HWW_1", h_ZprimeCandidates));
  // NN_CR_selection_1.reset(new TaggerCut(0.0, 0.2, "MC_HWW_1", h_ZprimeCandidates));
  // NN_SR_selection_2.reset(new TaggerCut(0.8, 1.0, "MC_HWW_2", h_ZprimeCandidates));
  // NN_CR_selection_2.reset(new TaggerCut(0.0, 0.2, "MC_HWW_2", h_ZprimeCandidates));
  // CNN_SR_selection.reset(new TaggerCut(0.8, 1.0, "CNN_HWW", h_ZprimeCandidates));
  // CNN_CR_selection.reset(new TaggerCut(0.0, 0.2, "CNN_HWW", h_ZprimeCandidates));


  btag_DeepBoosted_H4qvsQCD_SR_selection.reset(new TaggerCut(0.8, 1.0, -1, "btag_DeepBoosted_H4qvsQCD", h_ZprimeCandidates));
  btag_DeepBoosted_H4qvsQCD_CR_selection.reset(new TaggerCut(0.0, 0.2, -1, "btag_DeepBoosted_H4qvsQCD", h_ZprimeCandidates));



  btag_DeepBoosted_H4qvsQCDp2_SR_selection.reset(new TaggerCut(0.2, 1.0, -1, "btag_DeepBoosted_H4qvsQCD", h_ZprimeCandidates));
  btag_DeepBoosted_H4qvsQCDp2_CR_selection.reset(new TaggerCut(0.0, 0.2, -1, "btag_DeepBoosted_H4qvsQCD", h_ZprimeCandidates));
  btag_DeepBoosted_H4qvsQCDp02_SR_selection.reset(new TaggerCut(0.02, 1.0, -1, "btag_DeepBoosted_H4qvsQCD", h_ZprimeCandidates));
  btag_DeepBoosted_H4qvsQCDp02_CR_selection.reset(new TaggerCut(0.0, 0.02, -1, "btag_DeepBoosted_H4qvsQCD", h_ZprimeCandidates));
  btag_DeepBoosted_H4qvsQCDpt1000_SR_selection.reset(new TaggerCut(0.8, 1.0, 1000, "btag_DeepBoosted_H4qvsQCD", h_ZprimeCandidates));
  btag_DeepBoosted_H4qvsQCDpt1000_CR_selection.reset(new TaggerCut(0.0, 0.8, 1000, "btag_DeepBoosted_H4qvsQCD", h_ZprimeCandidates));
  btag_DeepBoosted_H4qvsQCDpt1000p2_SR_selection.reset(new TaggerCut(0.2, 1.0, 1000, "btag_DeepBoosted_H4qvsQCD", h_ZprimeCandidates));
  btag_DeepBoosted_H4qvsQCDpt1000p2_CR_selection.reset(new TaggerCut(0.0, 0.2, 1000, "btag_DeepBoosted_H4qvsQCD", h_ZprimeCandidates));
  btag_DeepBoosted_H4qvsQCDpt1000p02_SR_selection.reset(new TaggerCut(0.02, 1.0, 1000, "btag_DeepBoosted_H4qvsQCD", h_ZprimeCandidates));
  btag_DeepBoosted_H4qvsQCDpt1000p02_CR_selection.reset(new TaggerCut(0.0, 0.02, 1000, "btag_DeepBoosted_H4qvsQCD", h_ZprimeCandidates));

  btag_DeepBoosted_H4qvsQCDpt1500_SR_selection.reset(new TaggerCut(0.8, 1.0, 1500, "btag_DeepBoosted_H4qvsQCD", h_ZprimeCandidates));
  btag_DeepBoosted_H4qvsQCDpt1500_CR_selection.reset(new TaggerCut(0.0, 0.8, 1500, "btag_DeepBoosted_H4qvsQCD", h_ZprimeCandidates));
  btag_DeepBoosted_H4qvsQCDpt1500p2_SR_selection.reset(new TaggerCut(0.2, 1.0, 1500, "btag_DeepBoosted_H4qvsQCD", h_ZprimeCandidates));
  btag_DeepBoosted_H4qvsQCDpt1500p2_CR_selection.reset(new TaggerCut(0.0, 0.2, 1500, "btag_DeepBoosted_H4qvsQCD", h_ZprimeCandidates));
  btag_DeepBoosted_H4qvsQCDpt1500p02_SR_selection.reset(new TaggerCut(0.02, 1.0, 1500, "btag_DeepBoosted_H4qvsQCD", h_ZprimeCandidates));
  btag_DeepBoosted_H4qvsQCDpt1500p02_CR_selection.reset(new TaggerCut(0.0, 0.02, 1500, "btag_DeepBoosted_H4qvsQCD", h_ZprimeCandidates));


  btag_DeepBoosted_HbbvsQCD_SR_selection.reset(new TaggerCut(0.8, 1.0, -1, "btag_DeepBoosted_HbbvsQCD", h_ZprimeCandidates));
  btag_DeepBoosted_HbbvsQCD_CR_selection.reset(new TaggerCut(0.0, 0.2, -1, "btag_DeepBoosted_HbbvsQCD", h_ZprimeCandidates));

  // btag_DeepBoosted_probHbb_SR_selection.reset(new TaggerCut(0.8, 1.0, -1, "btag_DeepBoosted_probHbb", h_ZprimeCandidates));
  // btag_DeepBoosted_probHbb_CR_selection.reset(new TaggerCut(0.0, 0.2, -1, "btag_DeepBoosted_probHbb", h_ZprimeCandidates));
  //
  // btag_DeepBoosted_probHqqqq_SR_selection.reset(new TaggerCut(0.8, 1.0, -1, "btag_DeepBoosted_probHqqqq", h_ZprimeCandidates));
  // btag_DeepBoosted_probHqqqq_CR_selection.reset(new TaggerCut(0.0, 0.2, -1, "btag_DeepBoosted_probHqqqq", h_ZprimeCandidates));

  std::string groom = MB["isHOTVR"]? "_groomed": "";

  // tau21_SR_selection.reset(new TaggerCut(0.0, 0.6, -1, "tau21"+groom, h_ZprimeCandidates));
  // tau21_CR_selection.reset(new TaggerCut(0.6, 1.0, -1, "tau21"+groom, h_ZprimeCandidates));
  //
  // tau31_SR_selection.reset(new TaggerCut(0.0, 0.4, -1, "tau31"+groom, h_ZprimeCandidates));
  // tau31_CR_selection.reset(new TaggerCut(0.4, 1.0, -1, "tau31"+groom, h_ZprimeCandidates));

  // tau41_SR_selection.reset(new TaggerCut(0.0, 0.3, -1, "tau41"+groom, h_ZprimeCandidates));
  // tau41_CR_selection.reset(new TaggerCut(0.3, 1.0, -1, "tau41"+groom, h_ZprimeCandidates));
  //
  // tau32_SR_selection.reset(new TaggerCut(0.0, 0.6, -1, "tau32"+groom, h_ZprimeCandidates));
  // tau32_CR_selection.reset(new TaggerCut(0.6, 1.0, -1, "tau32"+groom, h_ZprimeCandidates));

  tau42_SR_selection.reset(new TaggerCut(0.0, 0.5, -1, "tau42"+groom, h_ZprimeCandidates));
  tau42_CR_selection.reset(new TaggerCut(0.5, 1.0, -1, "tau42"+groom, h_ZprimeCandidates));
  //
  // tau43_SR_selection.reset(new TaggerCut(0.0, 0.8, -1, "tau43"+groom, h_ZprimeCandidates));
  // tau43_CR_selection.reset(new TaggerCut(0.8, 1.0, -1, "tau43"+groom, h_ZprimeCandidates));
  //

}


/*
█ ██████  ██████   ██████   ██████ ███████ ███████ ███████
█ ██   ██ ██   ██ ██    ██ ██      ██      ██      ██
█ ██████  ██████  ██    ██ ██      █████   ███████ ███████
█ ██      ██   ██ ██    ██ ██      ██           ██      ██
█ ██      ██   ██  ██████   ██████ ███████ ███████ ███████
*/

bool SignalRegionModule::process(uhh2::Event& event) {

  event.weight = event.get(WFolder("weight_GLP_in"));

  fill_histograms(event, "Selection");

  if (event.get(h_ZprimeCandidates).size()!=1) return false; // TODO!!!
  // if (event.get(h_ZprimeCandidates).size()!=1 && event.get(h_ZprimeCandidates)[0].discriminator("SDmass")<50) return false; // TODO!!!
  //
  //
  // if (event.get(h_ZprimeCandidates).size()!=1)std::cout << "new Event" << '\n';
  // for (auto cand : event.get(h_ZprimeCandidates)) {
  //   if (event.get(h_ZprimeCandidates).size()!=1) {
  //     std::cout << "\t " << cand.discriminator("chi2") << " " << cand.Zprime_mass() << " " << cand.H().v4().M() << " " << cand.discriminator("SDmass") << " " << MatchingStatusToString(cand.discriminator("MatchingStatus")) << " " << MatchingToString(cand.discriminator("Match"));
  //     // std::cout << " " << cand.discriminator("tau1") << " " << cand.discriminator("tau2") << " " << cand.discriminator("tau3") << " " << cand.discriminator("tau4");
  //     // std::cout << " " << cand.Z().pt() << " " <<  cand.leptons()[0].pt() << " " << cand.leptons()[1].pt();
  //     // std::cout << " " << cand.Z().pt() << " " << cand.Z().eta() << " " <<  cand.Z().phi();
  //     std::cout << "\n\t\t" << cand.H().pt() << " " << cand.H().eta() << " " <<  cand.H().phi() ;
  //     std::cout << "\n\t\t" << uhh2::deltaPhi(cand.Z(),cand.H()) << " " << uhh2::deltaR(cand.Z(),cand.H()) << " " << uhh2::deltaPhi(cand.leptons()[0],cand.H()) << " " << uhh2::deltaR(cand.leptons()[0],cand.H()) << " " << uhh2::deltaPhi(cand.leptons()[1],cand.H()) << " " << uhh2::deltaR(cand.leptons()[1],cand.H());
  //     std::cout << "\n\t\t" << uhh2::deltaPhi(event.get(h_ZprimeCandidates)[0].H(),cand.H()) << " " << uhh2::deltaR(event.get(h_ZprimeCandidates)[0].H(),cand.H()) << " " << fabs(event.get(h_ZprimeCandidates)[0].H().eta() - cand.H().eta());
  //     std::cout << '\n';
  //   }
  // }
  //
  // if (event.get(h_ZprimeCandidates).size()!=1) {
  //   export_weights(event);
  //
  //   return false;
  // }

  // if (MB["is_ZH"]){
  //   ZprimeDecay ZDecay = static_cast<ZprimeDecay>(int(event.get(WFolder("ZDecay_in"))));
  //   ZprimeDecay HDecay = static_cast<ZprimeDecay>(int(event.get(WFolder("HDecay_in"))));
  //   if (MB["muonchannel"]     && ZDecay != Zmumu) return false;
  //   if (MB["electronchannel"] && ZDecay != Zee)   return false;
  //   if (MB["is_bb"]           && HDecay != Hbb)   return false;
  //   if (MB["is_WW"]           && HDecay != HWW)   return false;
  //   if (MB["is_else"]         && HDecay != Helse) return false;
  // }

  // fill_histograms(event, "FSSelection");

  // if(NN_SR_selection->passes(event))    fill_histograms(event, "NN_SR");
  // if(NN_CR_selection->passes(event))    fill_histograms(event, "NN_CR");
  // if(NN_SR_selection_1->passes(event))  fill_histograms(event, "NN_1_SR");
  // if(NN_CR_selection_1->passes(event))  fill_histograms(event, "NN_1_CR");
  // if(NN_SR_selection_2->passes(event))  fill_histograms(event, "NN_2_SR");
  // if(NN_CR_selection_2->passes(event))  fill_histograms(event, "NN_2_CR");
  // if(CNN_SR_selection->passes(event))  fill_histograms(event, "CNN_SR");
  // if(CNN_CR_selection->passes(event))  fill_histograms(event, "CNN_CR");

  if(btag_DeepBoosted_H4qvsQCD_SR_selection->passes(event))   fill_histograms(event, "btag_DeepBoosted_H4qvsQCD_SR");
  if(btag_DeepBoosted_H4qvsQCD_CR_selection->passes(event))   fill_histograms(event, "btag_DeepBoosted_H4qvsQCD_CR");

  if(btag_DeepBoosted_H4qvsQCDp2_SR_selection->passes(event))   fill_histograms(event, "btag_DeepBoosted_H4qvsQCDp2_SR");
  if(btag_DeepBoosted_H4qvsQCDp2_CR_selection->passes(event))   fill_histograms(event, "btag_DeepBoosted_H4qvsQCDp2_CR");
  if(btag_DeepBoosted_H4qvsQCDp02_SR_selection->passes(event))   fill_histograms(event, "btag_DeepBoosted_H4qvsQCDp02_SR");
  if(btag_DeepBoosted_H4qvsQCDp02_CR_selection->passes(event))   fill_histograms(event, "btag_DeepBoosted_H4qvsQCDp02_CR");
  if(btag_DeepBoosted_H4qvsQCDpt1000_SR_selection->passes(event))   fill_histograms(event, "btag_DeepBoosted_H4qvsQCDpt1000_SR");
  if(btag_DeepBoosted_H4qvsQCDpt1000_CR_selection->passes(event))   fill_histograms(event, "btag_DeepBoosted_H4qvsQCDpt1000_CR");
  if(btag_DeepBoosted_H4qvsQCDpt1000p2_SR_selection->passes(event))   fill_histograms(event, "btag_DeepBoosted_H4qvsQCDpt1000p2_SR");
  if(btag_DeepBoosted_H4qvsQCDpt1000p2_CR_selection->passes(event))   fill_histograms(event, "btag_DeepBoosted_H4qvsQCDpt1000p2_CR");
  if(btag_DeepBoosted_H4qvsQCDpt1000p02_SR_selection->passes(event))   fill_histograms(event, "btag_DeepBoosted_H4qvsQCDpt1000p02_SR");
  if(btag_DeepBoosted_H4qvsQCDpt1000p02_CR_selection->passes(event))   fill_histograms(event, "btag_DeepBoosted_H4qvsQCDpt1000p02_CR");

  if(btag_DeepBoosted_H4qvsQCDpt1500_SR_selection->passes(event))   fill_histograms(event, "btag_DeepBoosted_H4qvsQCDpt1500_SR");
  if(btag_DeepBoosted_H4qvsQCDpt1500_CR_selection->passes(event))   fill_histograms(event, "btag_DeepBoosted_H4qvsQCDpt1500_CR");
  if(btag_DeepBoosted_H4qvsQCDpt1500p2_SR_selection->passes(event))   fill_histograms(event, "btag_DeepBoosted_H4qvsQCDpt1500p2_SR");
  if(btag_DeepBoosted_H4qvsQCDpt1500p2_CR_selection->passes(event))   fill_histograms(event, "btag_DeepBoosted_H4qvsQCDpt1500p2_CR");
  if(btag_DeepBoosted_H4qvsQCDpt1500p02_SR_selection->passes(event))   fill_histograms(event, "btag_DeepBoosted_H4qvsQCDpt1500p02_SR");
  if(btag_DeepBoosted_H4qvsQCDpt1500p02_CR_selection->passes(event))   fill_histograms(event, "btag_DeepBoosted_H4qvsQCDpt1500p02_CR");


  if(btag_DeepBoosted_HbbvsQCD_SR_selection->passes(event))   fill_histograms(event, "btag_DeepBoosted_HbbvsQCD_SR");
  if(btag_DeepBoosted_HbbvsQCD_CR_selection->passes(event))   fill_histograms(event, "btag_DeepBoosted_HbbvsQCD_CR");

  // if(btag_DeepBoosted_probHbb_SR_selection->passes(event))    fill_histograms(event, "btag_DeepBoosted_probHbb_SR");
  // if(btag_DeepBoosted_probHbb_CR_selection->passes(event))    fill_histograms(event, "btag_DeepBoosted_probHbb_CR");
  //
  // if(btag_DeepBoosted_probHqqqq_SR_selection->passes(event))  fill_histograms(event, "btag_DeepBoosted_probHqqqq_SR");
  // if(btag_DeepBoosted_probHqqqq_CR_selection->passes(event))  fill_histograms(event, "btag_DeepBoosted_probHqqqq_CR");

  // if(tau21_SR_selection->passes(event)) fill_histograms(event, "tau21_SR");
  // if(tau21_CR_selection->passes(event)) fill_histograms(event, "tau21_CR");

  // if(tau31_SR_selection->passes(event)) fill_histograms(event, "tau31_SR");
  // if(tau31_CR_selection->passes(event)) fill_histograms(event, "tau31_CR");

  // if(tau41_SR_selection->passes(event)) fill_histograms(event, "tau41_SR");
  // if(tau41_CR_selection->passes(event)) fill_histograms(event, "tau41_CR");

  // if(tau32_SR_selection->passes(event)) fill_histograms(event, "tau32_SR");
  // if(tau32_CR_selection->passes(event)) fill_histograms(event, "tau32_CR");

  if(tau42_SR_selection->passes(event)) fill_histograms(event, "tau42_SR");
  if(tau42_CR_selection->passes(event)) fill_histograms(event, "tau42_CR");

  // if(tau43_SR_selection->passes(event)) fill_histograms(event, "tau43_SR");
  // if(tau43_CR_selection->passes(event)) fill_histograms(event, "tau43_CR");

  export_weights(event);

  return true;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the SignalRegionModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(SignalRegionModule)
