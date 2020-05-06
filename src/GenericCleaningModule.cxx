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
#include "UHH2/HiggsToWWTagger/include/PrintingFunctions.hpp"
#include "UHH2/HiggsToWWTagger/include/GenericJetCleaner.h"
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

class GenericCleaningModule: public ModuleBASE {

public:

  explicit GenericCleaningModule(uhh2::Context&);
  virtual bool process(uhh2::Event&) override;
  void book_histograms(uhh2::Context&);
  void fill_histograms(uhh2::Event&, string);
  void book_handles(uhh2::Context&);
  void PrintInputs();

protected:

  // Define variables
  std::vector<std::string> histogram_tags = { "nocuts", "NBoostedJetCHS", "NBoostedJetPuppi", "NBoostedJetHOTVR"};

  std::unordered_map<std::string, std::string> strings;
  std::unordered_map<std::string, bool> bools;

  Event::Handle<std::vector<Jet> > h_jetsCHS, h_jetsPuppi;
  Event::Handle<std::vector<TopJet> > h_topjetsCHS, h_topjetsPuppi, h_topjetsHOTVR;

  // Define common modules
  std::unique_ptr<uhh2::Selection> lumi_selection;
  std::vector<std::unique_ptr<AnalysisModule>> weightsmodules, modules;
  std::unique_ptr<uhh2::AndSelection> metfilters_selection;
  std::unique_ptr<GenericJetCleaner> GJCCHS, GJCTopCHS, GJCPuppi, GJCTopPuppi, GJCTopHOTVR;

  // Define selections
  std::shared_ptr<Selection> NBoostedJetSelCHS, NBoostedJetSelPuppi, NBoostedJetSelHOTVR;

};

void GenericCleaningModule::book_handles(uhh2::Context& ctx) {
  string tag;
  tag = "weight_lumi";  book_WFolder(tag, new Event::Handle< float >, ctx.declare_event_output< float >(tag));
  tag = "weight_GLP";   book_WFolder(tag, new Event::Handle< float >, ctx.declare_event_output< float >(tag));
  tag = "weight_pu";    book_WFolder(tag, new Event::Handle< float >, bools["is_mc"]? ctx.get_handle< float >(tag) : ctx.declare_event_output< float >(tag));
}

void GenericCleaningModule::PrintInputs() {
  std::cout << "****************************************" << std::endl;
  std::cout << "          GenericCleaningModule         " << std::endl;
  std::cout << "----------------------------------------" << std::endl;
  for (auto x : strings) std::cout << x.first << std::string( 18-x.first.size(), ' ' ) << x.second << '\n';
  for (auto x : bools)   std::cout << x.first << std::string( 18-x.first.size(), ' ' ) << (x.second? "true" : "false") << '\n';
  std::cout << "****************************************\n" << std::endl;
}

void GenericCleaningModule::book_histograms(uhh2::Context& ctx){
  for(const auto & tag : histogram_tags){
    string mytag;
    mytag = "event_"        + tag; book_HFolder(mytag, new EventHists(ctx,mytag));
    mytag = "gen_"          + tag; book_HFolder(mytag, new GenMatchHists(ctx,mytag));
    mytag = "nTopJetCHS_"   + tag; book_HFolder(mytag, new ExtJetHists(ctx,mytag,strings["topjetLabelCHS"]));
    mytag = "nTopJetPuppi_" + tag; book_HFolder(mytag, new ExtJetHists(ctx,mytag,strings["topjetLabelPuppi"]));
    mytag = "nTopJetHOTVR_" + tag; book_HFolder(mytag, new ExtJetHists(ctx,mytag,strings["topjetLabelHOTVR"]));
  }
}

void GenericCleaningModule::fill_histograms(uhh2::Event& event, string tag){
  std::vector<string> mytags = {"event_", "gen_", "nTopJetCHS_", "nTopJetPuppi_", "nTopJetHOTVR_"};
  for (auto& mytag : mytags) HFolder(mytag+ tag)->fill(event);
}

/*
█  ██████  ██████  ███    ██ ███████ ████████ ██████  ██    ██  ██████ ████████  ██████  ██████
█ ██      ██    ██ ████   ██ ██         ██    ██   ██ ██    ██ ██         ██    ██    ██ ██   ██
█ ██      ██    ██ ██ ██  ██ ███████    ██    ██████  ██    ██ ██         ██    ██    ██ ██████
█ ██      ██    ██ ██  ██ ██      ██    ██    ██   ██ ██    ██ ██         ██    ██    ██ ██   ██
█  ██████  ██████  ██   ████ ███████    ██    ██   ██  ██████   ██████    ██     ██████  ██   ██
*/

GenericCleaningModule::GenericCleaningModule(uhh2::Context& ctx){

  // Set up variables
  strings["year"]           = ctx.get("year");
  strings["SysType_PU"]     = ctx.get("SysType_PU");
  bools["is_mc"]            = ctx.get("dataset_type") == "MC";
  bools["lumisel"]          = string2bool(ctx.get("lumisel"));
  bools["mclumiweight"]     = string2bool(ctx.get("mclumiweight"));
  bools["mcpileupreweight"] = string2bool(ctx.get("mcpileupreweight"));
  bools["isPuppi"]          = string2bool(ctx.get("isPuppi"));
  bools["isCHS"]            = string2bool(ctx.get("isCHS"));
  bools["isHOTVR"]          = string2bool(ctx.get("isHOTVR"));
  bools["eleid"]            = string2bool(ctx.get("eleid"));
  bools["muid"]             = string2bool(ctx.get("muid"));
  bools["tauid"]            = string2bool(ctx.get("tauid"));
  bools["metfilters"]       = string2bool(ctx.get("metfilters"));

  strings["jetLabelCHS"] = "jets";
  strings["jetLabelPuppi"] = "jetsAk4Puppi";
  strings["topjetLabelCHS"] = "topjets";
  strings["topjetLabelPuppi"] = "toppuppijets";
  strings["topjetLabelHOTVR"] = "hotvrPuppi";

  h_jetsCHS   = ctx.get_handle<std::vector<Jet>>(strings["jetLabelCHS"]);
  h_jetsPuppi = ctx.get_handle<std::vector<Jet>>(strings["jetLabelPuppi"]);
  h_topjetsCHS   = ctx.get_handle<std::vector<TopJet>>(strings["topjetLabelCHS"]);
  h_topjetsPuppi = ctx.get_handle<std::vector<TopJet>>(strings["topjetLabelPuppi"]);
  h_topjetsHOTVR = ctx.get_handle<std::vector<TopJet>>(strings["topjetLabelHOTVR"]);

  // Set up histograms:

  book_histograms(ctx);
  book_handles(ctx);
  PrintInputs();

  // Set up selections
  const MuonId muoId(AndId<Muon> (MuonID(Muon::CutBasedIdTrkHighPt), PtEtaCut(30, min_lepton_eta)));
  const ElectronId eleId(AndId<Electron>(ElectronID_MVA_Fall17_loose_noIso, PtEtaSCCut(min_lepton_pt, min_lepton_eta)));

  const JetId jetIdCHS(AndId<Jet>   (JetPFID(JetPFID::WP_TIGHT_CHS), PtEtaCut(30, min_lepton_eta),  NoLeptonInJet("all", eleId, muoId, bools["isHOTVR"]? 0.8: -1)));
  const JetId jetIdPuppi(AndId<Jet> (JetPFID(JetPFID::WP_TIGHT_PUPPI),PtEtaCut(30, min_lepton_eta), NoLeptonInJet("all", eleId, muoId, bools["isHOTVR"]? 0.8: -1)));
  const TopJetId topjetIdCHS(AndId<TopJet>   (JetPFID(JetPFID::WP_TIGHT_CHS),   PtEtaCut(200, min_lepton_eta), NoLeptonInJet("all", eleId, muoId, bools["isHOTVR"]? 0.8: -1)));
  const TopJetId topjetIdPuppi(AndId<TopJet> (JetPFID(JetPFID::WP_TIGHT_PUPPI), PtEtaCut(200, min_lepton_eta), NoLeptonInJet("all", eleId, muoId, bools["isHOTVR"]? 0.8: -1)));

  PrimaryVertexId pvid = StandardPrimaryVertexId(); // TODO
  modules.emplace_back(new PrimaryVertexCleaner(pvid));
  if(bools["is_mc"]) {
    // ctx.declare_event_input<std::vector<Particle> >(ctx.get("TopJetCollectionGEN"), "topjetsGEN"); THIS IS AK8 but not softdropmass
    if(bools["mclumiweight"])  weightsmodules.emplace_back(new MCLumiWeight(ctx));
    if(bools["mcpileupreweight"]) weightsmodules.emplace_back(new MCPileupReweight(ctx,strings["SysType_PU"]));
  } else {
    if(bools["lumisel"]) lumi_selection.reset(new LumiSelection(ctx));
  }

  if(bools["eleid"]) modules.emplace_back(new ElectronCleaner(eleId));
  if(bools["muid"])  modules.emplace_back(new MuonCleaner(muoId));

  weightsmodules.emplace_back(new GenLevelJetMatch(ctx,strings["topjetLabelCHS"]));
  weightsmodules.emplace_back(new GenLevelJetMatch(ctx,strings["topjetLabelPuppi"]));
  weightsmodules.emplace_back(new GenLevelJetMatch(ctx,strings["topjetLabelHOTVR"]));
  weightsmodules.emplace_back(new FinalStateMatching(ctx));

  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2
  metfilters_selection.reset(new AndSelection(ctx, "metfilters"));
  metfilters_selection->add<TriggerSelection>("goodVertices", "Flag_goodVertices");
  metfilters_selection->add<NPVSelection>("1 good PV",1,-1,pvid); /* Not a metfilter. Used to select 1 good PV */
  metfilters_selection->add<TriggerSelection>("globalSuperTightHalo2016Filter", "Flag_globalSuperTightHalo2016Filter");
  metfilters_selection->add<TriggerSelection>("HBHENoiseFilter", "Flag_HBHENoiseFilter");
  metfilters_selection->add<TriggerSelection>("HBHENoiseIsoFilter", "Flag_HBHENoiseIsoFilter");
  metfilters_selection->add<TriggerSelection>("EcalDeadCellTriggerPrimitiveFilter", "Flag_EcalDeadCellTriggerPrimitiveFilter");
  if (strings["year"] != "2016") metfilters_selection->add<EcalBadCalibSelection>("EcalBadCalibSelection"); /*TODO check 2016*/ // Use this instead of Flag_ecalBadCalibFilter, uses ecalBadCalibReducedMINIAODFilter in ntuple_generator
  if (strings["year"] != "2016") metfilters_selection->add<TriggerSelection>("BadPFMuonFilter", "Flag_BadPFMuonFilter"); /*TODO check 2016, maybe Extra_BadPFMuonFilter */
  if (!bools["is_mc"]) metfilters_selection->add<TriggerSelection>("eeBadScFilter", "Flag_eeBadScFilter"); /* TODO Not recommended for MC, but do check */
  /* metfilters_selection->add<TriggerSelection>("BadChargedCandidateFilter", "Flag_BadChargedCandidateFilter"); TODO Not recommended, under review.Separate module in ntuple_generator for 2016v2*/


  GJCCHS.reset(new      GenericJetCleaner(ctx, strings["jetLabelCHS"],       false, jetIdCHS,   topjetIdCHS,   muoId,  eleId));
  GJCTopCHS.reset(new   GenericJetCleaner(ctx, strings["topjetLabelCHS"],    true,  jetIdCHS,   topjetIdCHS,   muoId,  eleId));
  GJCPuppi.reset(new    GenericJetCleaner(ctx, strings["jetLabelPuppi"],     false, jetIdPuppi, topjetIdPuppi, muoId,  eleId));
  GJCTopPuppi.reset(new GenericJetCleaner(ctx, strings["topjetLabelPuppi"],  true,  jetIdPuppi, topjetIdPuppi, muoId,  eleId));
  GJCTopHOTVR.reset(new GenericJetCleaner(ctx, strings["topjetLabelHOTVR"],  true,  jetIdPuppi, topjetIdPuppi, muoId,  eleId));

  NBoostedJetSelCHS.reset(new NTopJetSelection(1,-1,topjetIdCHS,h_topjetsCHS));
  NBoostedJetSelPuppi.reset(new NTopJetSelection(1,-1,topjetIdPuppi,h_topjetsPuppi));
  NBoostedJetSelHOTVR.reset(new NTopJetSelection(1,-1,topjetIdPuppi,h_topjetsHOTVR));

  for (auto x: {"offlineSlimmedPrimaryVertices", "slimmedElectronsUSER", "slimmedMuonsUSER", "genInfo", "GenParticles", "slimmedGenJets",
  "genjetsAk8SubstructureSoftDrop", "jetsAk4CHS", "jetsAk4Puppi", "slimmedMETs", "slimmedMETsPuppi", "triggerNames", "triggerPrescales",
  "triggerPrescalesL1max", "triggerPrescalesL1min", "triggerResults", "beamspot_x0", "beamspot_y0", "beamspot_z0", "isRealData",
  "luminosityBlock", "passEcalBadCalib", "prefiringWeight", "prefiringWeightDown", "prefiringWeightUp", "rho", "run", "year", "hotvrGen"}) {
    ctx.undeclare_event_output(x);
  }

}


/*
█ ██████  ██████   ██████   ██████ ███████ ███████ ███████
█ ██   ██ ██   ██ ██    ██ ██      ██      ██      ██
█ ██████  ██████  ██    ██ ██      █████   ███████ ███████
█ ██      ██   ██ ██    ██ ██      ██           ██      ██
█ ██      ██   ██  ██████   ██████ ███████ ███████ ███████
*/

bool GenericCleaningModule::process(uhh2::Event& event) {

  if ((event.year).find(strings["year"])==std::string::npos) throw std::runtime_error("In GenericCleaningModule.cxx: You are running on "+event.year+" sample with a "+strings["year"]+" year config file. Fix this.");

  auto weight_gen = event.weight;
  fill_histograms(event, "nocuts");

  if(event.isRealData && bools["lumisel"]) if(!lumi_selection->passes(event)) return false;

  //  MCLumiWeight, MCPileupReweight, GenLevelJetMatch, FinalStateMatching
  for(auto & m : weightsmodules) m->process(event);

  double weight_pu = 1;
  if (bools["is_mc"]) weight_pu = event.get(WFolder("weight_pu"));
  else event.set(WFolder("weight_pu"), weight_pu);
  event.set(WFolder("weight_lumi"), event.weight/(weight_gen*weight_pu));
  event.set(WFolder("weight_GLP"), event.weight);

  // PrimaryVertexCleaner, ElectronCleaner, MuonCleaner
  for(auto & m : modules) m->process(event);

  if(bools["metfilters"]) if(!metfilters_selection->passes(event)) return false;

  GJCCHS->process(event);
  GJCTopCHS->process(event);
  GJCPuppi->process(event);
  GJCTopPuppi->process(event);
  GJCTopHOTVR->process(event);

  sort_by_pt<Muon>(*event.muons);
  sort_by_pt<Electron>(*event.electrons);
  sort_by_pt<Jet>(event.get(h_jetsCHS));
  sort_by_pt<Jet>(event.get(h_jetsPuppi));
  sort_by_pt<TopJet>(event.get(h_topjetsCHS));
  sort_by_pt<TopJet>(event.get(h_topjetsPuppi));
  sort_by_pt<TopJet>(event.get(h_topjetsHOTVR));

  bool passBoostCHS   = NBoostedJetSelCHS->passes(event);
  bool passBoostPuppi = NBoostedJetSelPuppi->passes(event);
  bool passBoostHOTVR = NBoostedJetSelHOTVR->passes(event);

  if (passBoostCHS) fill_histograms(event, "NBoostedJetCHS");
  if (passBoostPuppi) fill_histograms(event, "NBoostedJetPuppi");
  if (passBoostHOTVR) fill_histograms(event, "NBoostedJetHOTVR");

  if(!passBoostCHS && !passBoostPuppi && !passBoostHOTVR) return false;

  return true;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the GenericCleaningModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(GenericCleaningModule)
