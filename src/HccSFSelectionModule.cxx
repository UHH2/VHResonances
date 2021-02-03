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
#include "UHH2/common/include/JetHists.h"
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
#include "UHH2/common/include/CollectionProducer.h"
#include "UHH2/common/include/DetectorCleaning.h"
#include "UHH2/common/include/PDFWeights.h"

#include "UHH2/VHResonances/include/ModuleBase.h"
#include "UHH2/VHResonances/include/constants.hpp"
#include "UHH2/VHResonances/include/ExtJetHists.h"
#include "UHH2/VHResonances/include/ExtJetHistsWithConditions.h"
#include "UHH2/VHResonances/include/GenMatchHists.h"
#include "UHH2/VHResonances/include/GenLevelJetMatch.h"
#include "UHH2/VHResonances/include/PrintingFunctions.hpp"
#include "UHH2/VHResonances/include/GenericJetCleaner.h"
#include "UHH2/VHResonances/include/HiggsToWWSelection.h"
#include "UHH2/VHResonances/include/HiggsToWWModules.h"
#include "UHH2/VHResonances/include/HiggsToWWHists.h"
#include "UHH2/VHResonances/include/DiLeptonHists.h"

using namespace std;

/*
█ ██████  ███████ ███████ ██ ███    ██ ██ ████████ ██  ██████  ███    ██
█ ██   ██ ██      ██      ██ ████   ██ ██    ██    ██ ██    ██ ████   ██
█ ██   ██ █████   █████   ██ ██ ██  ██ ██    ██    ██ ██    ██ ██ ██  ██
█ ██   ██ ██      ██      ██ ██  ██ ██ ██    ██    ██ ██    ██ ██  ██ ██
█ ██████  ███████ ██      ██ ██   ████ ██    ██    ██  ██████  ██   ████
*/

class HccSFSelectionModule: public ModuleBASE {

public:

  explicit HccSFSelectionModule(uhh2::Context&);
  virtual bool process(uhh2::Event&) override;
  void book_histograms(uhh2::Context&);
  void fill_histograms(uhh2::Event&, string);
  void book_handles(uhh2::Context&);
  void PrintInputs();

protected:

  // Define variables
  std::string NameModule = "PreselectionModule";
  std::vector<std::string> histogram_tags = { "nocuts", "weights", "Trigger", "HEM", "cleaned", "NBoostedJet", "NLeptonSel", "LepInTopJet", "LepInJet"};

  std::unordered_map<std::string, std::string> MS;
  std::unordered_map<std::string, bool> MB;

  Event::Handle<std::vector<Jet> > h_jets;
  Event::Handle<std::vector<TopJet> > h_topjets;
  // Define common modules
  std::unique_ptr<uhh2::Selection> lumi_selection;
  std::vector<std::unique_ptr<AnalysisModule>> weightsmodules, modules;
  std::unique_ptr<uhh2::AndSelection> metfilters_selection;
  std::unique_ptr<GenericJetCleaner> GJC, GTJC;
  // std::unique_ptr<AnalysisModule> TTbarGenmodule;

  std::unique_ptr<TopJetCleaner> TopJetcleaner;
  std::unique_ptr<JetCleaner> Jetcleaner;

  // Define selections
  std::unordered_map<std::string, std::unique_ptr<Selection>> Trigger_selection;
  std::shared_ptr<Selection> NBoostedJetSel, CleanedJetSel, CleanedTopJetSel;
  std::shared_ptr<VetoSelection> VetoLeptonSel;
  std::shared_ptr<Selection> NoLeptonSel, NLeptonSel;
  std::unique_ptr<Selection> HEMEventCleaner_Selection;


};

void HccSFSelectionModule::book_handles(uhh2::Context& ctx) {
  string tag;
  tag = "weight_lumi";  book_WFolder(tag, new Event::Handle< float >, ctx.declare_event_output< float >(tag));
  tag = "weight_GLP";   book_WFolder(tag, new Event::Handle< float >, ctx.declare_event_output< float >(tag));
  tag = "weight_pu";    book_WFolder(tag, new Event::Handle< float >, MB["is_mc"]? ctx.get_handle< float >(tag) : ctx.declare_event_output< float >(tag));
}

void HccSFSelectionModule::PrintInputs() {
  std::cout << "****************************************" << std::endl;
  std::cout << "             "+NameModule+"             " << std::endl;
  std::cout << "----------------------------------------" << std::endl;
  for (auto x : MS) std::cout << x.first << std::string( 18-x.first.size(), ' ' ) << x.second << '\n';
  for (auto x : MB) std::cout << x.first << std::string( 18-x.first.size(), ' ' ) << BoolToString(x.second) << '\n';
  std::cout << "****************************************\n" << std::endl;
}

void HccSFSelectionModule::book_histograms(uhh2::Context& ctx) {
  for(const auto & tag : histogram_tags){
    string mytag;
    mytag = "event_"    + tag; book_HFolder(mytag, new EventHists(ctx,mytag));
    mytag = "gen_"      + tag; book_HFolder(mytag, new GenMatchHists(ctx,mytag));
    mytag = "nTopJet_"  + tag; book_HFolder(mytag, new ExtJetHists(ctx,mytag, MS["topjetLabel"]));
    mytag = "nJet_"     + tag; book_HFolder(mytag, new ExtJetHists(ctx,mytag, MS["jetLabel"]));
    mytag = "ele_"      + tag; book_HFolder(mytag, new ElectronHists(ctx,mytag, MS["topjetLabel"]));
    mytag = "muon_"     + tag; book_HFolder(mytag, new MuonHists(ctx,mytag, MS["topjetLabel"]));
    mytag = "diLepton_" + tag; book_HFolder(mytag, new DiLeptonHists(ctx,mytag, "", MS["topjetLabel"]));
    mytag = "BTagEff_"  + tag; book_HFolder(mytag, new BTagMCEfficiencyHists(ctx, mytag,BTag(BTag::DEEPCSV, BTag::WP_LOOSE), MS["topjetLabel"]));
    mytag = "Lumi_"     + tag; book_HFolder(mytag, new LuminosityHists(ctx, mytag));
  }
}

void HccSFSelectionModule::fill_histograms(uhh2::Event& event, string tag){
  std::vector<string> mytags = {"event_", "gen_", "nTopJet_", "nJet_", "ele_", "muon_", "diLepton_", "BTagEff_", "Lumi_"};
  for (auto& mytag : mytags) HFolder(mytag+ tag)->fill(event);
}

/*
█  ██████  ██████  ███    ██ ███████ ████████ ██████  ██    ██  ██████ ████████  ██████  ██████
█ ██      ██    ██ ████   ██ ██         ██    ██   ██ ██    ██ ██         ██    ██    ██ ██   ██
█ ██      ██    ██ ██ ██  ██ ███████    ██    ██████  ██    ██ ██         ██    ██    ██ ██████
█ ██      ██    ██ ██  ██ ██      ██    ██    ██   ██ ██    ██ ██         ██    ██    ██ ██   ██
█  ██████  ██████  ██   ████ ███████    ██    ██   ██  ██████   ██████    ██     ██████  ██   ██
*/

HccSFSelectionModule::HccSFSelectionModule(uhh2::Context& ctx){

  // Set up variables
  MS["year"]              = ctx.get("year");
  MS["dataset_version"]   = ctx.get("dataset_version");
  MB["is_mc"]             = ctx.get("dataset_type") == "MC";
  MB["isPuppi"]           = string2bool(ctx.get("isPuppi"));
  MB["isCHS"]             = string2bool(ctx.get("isCHS"));
  MB["isHOTVR"]           = string2bool(ctx.get("isHOTVR"));
  MB["charmchannel"]      = string2bool(ctx.get("charmchannel"));
  MS["SysType_PU"]        = ctx.get("SysType_PU");
  MB["lumisel"]           = string2bool(ctx.get("lumisel"));
  MB["mclumiweight"]      = string2bool(ctx.get("mclumiweight"));
  MB["mcpileupreweight"]  = string2bool(ctx.get("mcpileupreweight"));
  MB["eleid"]             = string2bool(ctx.get("eleid"));
  MB["muid"]              = string2bool(ctx.get("muid"));
  MB["tauid"]             = string2bool(ctx.get("tauid"));
  MB["metfilters"]        = string2bool(ctx.get("metfilters"));

  if ((MB["isPuppi"] && MB["isCHS"]) || (MB["isPuppi"] && MB["isHOTVR"]) || (MB["isCHS"] && MB["isHOTVR"]) ) throw std::runtime_error("In "+NameModule+".cxx: Choose exactly one jet collection.");
  if (!MB["charmchannel"]) throw std::runtime_error("In "+NameModule+".cxx: Only charm channel intended.");

  MS["jetLabel"]    = MB["isCHS"]? "jets":    (MB["isPuppi"]? "jetsAk4Puppi": (MB["isHOTVR"]? "jetsAk4Puppi": ""));
  MS["topjetLabel"] = MB["isCHS"]? "topjets": (MB["isPuppi"]? "toppuppijets": (MB["isHOTVR"]? "hotvrPuppi": ""));

  JetPFID::wp JETwp = MB["isCHS"]? JetPFID::WP_TIGHT_CHS: (MB["isPuppi"]? JetPFID::WP_TIGHT_PUPPI : (MB["isHOTVR"]? JetPFID::WP_TIGHT_PUPPI : JetPFID::WP_LOOSE_CHS));

  h_jets = ctx.get_handle<std::vector<Jet>>(MS["jetLabel"]);
  h_topjets = ctx.get_handle<std::vector<TopJet>>(MS["topjetLabel"]);

  // Set up histograms:

  book_histograms(ctx);
  book_handles(ctx);
  PrintInputs();

  // Set up selections

  const MuonId muoId = AndId<Muon>(PtEtaCut(8, min_lepton_eta), PtEtaCut(8, min_lepton_eta));
  const ElectronId eleId = AndId<Electron>(PtEtaSCCut(8, min_lepton_eta), PtEtaSCCut(8, min_lepton_eta));

  const JetId jetId = AndId<Jet> (JetPFID(JETwp), PtEtaCut(min_jet_pt, min_lepton_eta));
  const TopJetId topjetId = AndId<TopJet> (JetPFID(JETwp), PtEtaCut(min_topjet_pt, min_lepton_eta));
  // const TopJetId topjetId = AndId<TopJet> (JetPFID(JETwp), PtEtaCut(min_topjet_pt, min_lepton_eta), NoLeptonInJet("all", eleId, muoId, MB["isHOTVR"]? 0.8: -1));
  const TopJetId topjetIdSoftMuon = AndId<TopJet> (LeptonInJet("muon", eleId, muoId, 0.8), LeptonInJet("muon", eleId, muoId, 0.8));
  const JetId jetIdSoftMuon = AndId<Jet> (LeptonInJet("muon", eleId, muoId, 0.4), LeptonInJet("muon", eleId, muoId, 0.4));
  // TODO

  PrimaryVertexId pvid = StandardPrimaryVertexId(); // TODO
  modules.emplace_back(new PrimaryVertexCleaner(pvid));
  if(MB["is_mc"]) {
    if(MB["mclumiweight"])  weightsmodules.emplace_back(new MCLumiWeight(ctx));
    if(MB["mcpileupreweight"]) weightsmodules.emplace_back(new MCPileupReweight(ctx,MS["SysType_PU"]));
  } else {
    if(MB["lumisel"]) lumi_selection.reset(new LumiSelection(ctx));
  }

  if(MB["eleid"]) modules.emplace_back(new ElectronCleaner(eleId));
  if(MB["muid"])  modules.emplace_back(new MuonCleaner(muoId));

  weightsmodules.emplace_back(new GenLevelJetMatch(ctx,MS["topjetLabel"]));
  weightsmodules.emplace_back(new FinalStateMatching(ctx));
  weightsmodules.emplace_back(new NLOCorrections(ctx));

  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2
  metfilters_selection.reset(new AndSelection(ctx, "metfilters"));
  metfilters_selection->add<TriggerSelection>("goodVertices", "Flag_goodVertices");
  metfilters_selection->add<NPVSelection>("1 good PV",1,-1,pvid); /* Not a metfilter. Used to select 1 good PV */
  metfilters_selection->add<TriggerSelection>("globalSuperTightHalo2016Filter", "Flag_globalSuperTightHalo2016Filter");
  metfilters_selection->add<TriggerSelection>("HBHENoiseFilter", "Flag_HBHENoiseFilter");
  metfilters_selection->add<TriggerSelection>("HBHENoiseIsoFilter", "Flag_HBHENoiseIsoFilter");
  metfilters_selection->add<TriggerSelection>("EcalDeadCellTriggerPrimitiveFilter", "Flag_EcalDeadCellTriggerPrimitiveFilter");
  metfilters_selection->add<TriggerSelection>("BadPFMuonFilter", "Flag_BadPFMuonFilter");
  if (MS["year"] != "2016") metfilters_selection->add<EcalBadCalibSelection>("EcalBadCalibSelection"); /*TODO check 2016*/ // Use this instead of Flag_ecalBadCalibFilter, uses ecalBadCalibReducedMINIAODFilter in ntuple_generator
  if (!MB["is_mc"]) metfilters_selection->add<TriggerSelection>("eeBadScFilter", "Flag_eeBadScFilter"); /* TODO Not recommended for MC, but do check */
  /* metfilters_selection->add<TriggerSelection>("BadChargedCandidateFilter", "Flag_BadChargedCandidateFilter"); TODO Not recommended, under review.*/

  //Quick fix for Detector issues
  // HEMEventCleaner_Selection.reset(new HEMCleanerSelection(ctx, MS["jetLabel"], true, true, true));
  HEMEventCleaner_Selection.reset(new AndSelection(ctx)); // HEM important for inv channel only. DiLep selection reduces prob drastically

  GJC.reset( new GenericJetCleaner(ctx, MS["jetLabel"],    false, jetId, topjetId, muoId, eleId));
  GTJC.reset(new GenericJetCleaner(ctx, MS["topjetLabel"], true,  jetId, topjetId, muoId, eleId));

  // PDFReweight_module.reset(new PDFReweight(ctx));

  for (auto& t : Trigger_run_validity.at(MS["year"])) {
    //if (MB["charmchannel"] && !FindInString("AK8PFJet", t.first) ) continue;
    if (MB["charmchannel"] && !FindInString("PFHT", t.first) ) continue;
    Trigger_selection[t.first].reset(new TriggerSelection( t.first ));
  }

  NBoostedJetSel.reset(new NTopJetSelection(2,-1,topjetId,h_topjets));

  CleanedJetSel.reset(new NJetSelection(1,-1,jetIdSoftMuon,h_jets));
  CleanedTopJetSel.reset(new NTopJetSelection(2,-1,topjetIdSoftMuon,h_topjets));


  NLeptonSel.reset(new NMuonSelection(1));
  // NoLeptonSel.reset(new NElectronSelection(1));
  // VetoLeptonSel.reset(new VetoSelection(NoLeptonSel));

  // if (MB["do_TTgen"]) TTbarGenmodule.reset(new TTbarGenProducer(ctx));

  TopJetcleaner.reset(new TopJetCleaner(ctx, topjetIdSoftMuon, MS["topjetLabel"]));
  Jetcleaner.reset(new JetCleaner(ctx, jetIdSoftMuon, MS["jetLabel"]));

}


/*
█ ██████  ██████   ██████   ██████ ███████ ███████ ███████
█ ██   ██ ██   ██ ██    ██ ██      ██      ██      ██
█ ██████  ██████  ██    ██ ██      █████   ███████ ███████
█ ██      ██   ██ ██    ██ ██      ██           ██      ██
█ ██      ██   ██  ██████   ██████ ███████ ███████ ███████
*/

bool HccSFSelectionModule::process(uhh2::Event& event) {

  if ((event.year).find(MS["year"])==std::string::npos) throw std::runtime_error("In "+NameModule+".cxx: You are running on "+event.year+" sample with a "+MS["year"]+" year config file. Fix this.");

  auto weight_gen = event.weight;
  // if (MB["do_TTgen"]) TTbarGenmodule->process(event);
  fill_histograms(event, "nocuts");

  if(event.isRealData && MB["lumisel"]) if(!lumi_selection->passes(event)) return false;

  //  MCLumiWeight, MCPileupReweight, GenLevelJetMatch, FinalStateMatching, NLO corrections
  for(auto & m : weightsmodules) m->process(event);

  double weight_pu = 1;
  if (MB["is_mc"]) weight_pu = event.get(WFolder("weight_pu"));
  else event.set(WFolder("weight_pu"), weight_pu);
  event.set(WFolder("weight_lumi"), event.weight/(weight_gen*weight_pu));
  event.set(WFolder("weight_GLP"), event.weight);

  fill_histograms(event, "weights");

  bool pass_triggers_OR = false;

  for (auto& el : Trigger_selection) {
    if (event.isRealData && (event.run < Trigger_run_validity.at(MS["year"]).at(el.first).first || event.run > Trigger_run_validity.at(MS["year"]).at(el.first).second) ) continue;
    pass_triggers_OR += el.second->passes(event);
    if (pass_triggers_OR) break;
  }
  if (!pass_triggers_OR) return false;

  fill_histograms(event, "Trigger");

  //Effective in 2018 only.
  // Here the assumption is that it should check for cleaned jets to avoid overlap with leptons (Tight WP recommanded).
  if(!HEMEventCleaner_Selection->passes(event)) return false;
  fill_histograms(event, "HEM");

  // PrimaryVertexCleaner, ElectronCleaner, MuonCleaner
  for(auto & m : modules) m->process(event);

  GJC->process(event);
  GTJC->process(event);

  if(MB["metfilters"]) if(!metfilters_selection->passes(event)) return false;

  sort_by_pt<Muon>(*event.muons);
  sort_by_pt<Electron>(*event.electrons);
  sort_by_pt<Jet>(event.get(h_jets));
  sort_by_pt<TopJet>(event.get(h_topjets));

  fill_histograms(event, "cleaned");
  //
  // if(!VetoLeptonSel->passes(event)) return false;
  // fill_histograms(event, "Veto");

  if(!NBoostedJetSel->passes(event)) return false;
  fill_histograms(event, "NBoostedJet");

  if(!NLeptonSel->passes(event)) return false;
  fill_histograms(event, "NLeptonSel");

  TopJetcleaner->process(event);
  if(!CleanedTopJetSel->passes(event)) return false;
  fill_histograms(event, "LepInTopJet");

  Jetcleaner->process(event);
  if(!CleanedJetSel->passes(event)) return false;
  fill_histograms(event, "LepInJet");

  return true;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the HccSFSelectionModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(HccSFSelectionModule)
