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

#include "UHH2/VHResonances/include/ModuleBase.h"
#include "UHH2/VHResonances/include/constants.hpp"
#include "UHH2/VHResonances/include/ExtJetHists.h"
#include "UHH2/VHResonances/include/ExtJetHistsWithConditions.h"
#include "UHH2/VHResonances/include/GenMatchHists.h"
#include "UHH2/VHResonances/include/GenLevelJetMatch.h"
#include "UHH2/VHResonances/include/PrintingFunctions.hpp"
#include "UHH2/VHResonances/include/GenericJetCleaner.h"
#include "UHH2/VHResonances/include/HiggsToWWSelections.h"
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

class PreselectionModule: public ModuleBASE {

public:

  explicit PreselectionModule(uhh2::Context&);
  virtual bool process(uhh2::Event&) override;
  void book_histograms(uhh2::Context&);
  void fill_histograms(uhh2::Event&, string);
  void book_handles(uhh2::Context&);
  void PrintInputs();

protected:

  // Define variables
  std::string NameModule = "PreselectionModule";
  std::vector<std::string> histogram_tags = { "nocuts", "weights", "cleaned", "Veto", "NLeptonSel", "NBoostedJet", "DeltaRDiLepton", "JetDiLeptonPhiAngular"};

  std::unordered_map<std::string, std::string> MS;
  std::unordered_map<std::string, bool> MB;

  Event::Handle<std::vector<Jet> > h_jets;
  Event::Handle<std::vector<TopJet> > h_topjets;

  // Define common modules
  std::unique_ptr<uhh2::Selection> lumi_selection;
  std::vector<std::unique_ptr<AnalysisModule>> weightsmodules, modules;
  std::unique_ptr<uhh2::AndSelection> metfilters_selection;
  std::unique_ptr<GenericJetCleaner> GJC, GTJC;

  // Define selections
  std::shared_ptr<Selection> NBoostedJetSel;
  std::shared_ptr<VetoSelection> VetoLeptonSel;
  std::shared_ptr<Selection> NoLeptonSel, NLeptonSel, DeltaRDiLepton_selection, JetDiLeptonPhiAngularSel;

};

void PreselectionModule::book_handles(uhh2::Context& ctx) {
  string tag;
  tag = "weight_lumi";  book_WFolder(tag, new Event::Handle< float >, ctx.declare_event_output< float >(tag));
  tag = "weight_GLP";   book_WFolder(tag, new Event::Handle< float >, ctx.declare_event_output< float >(tag));
  tag = "weight_pu";    book_WFolder(tag, new Event::Handle< float >, MB["is_mc"]? ctx.get_handle< float >(tag) : ctx.declare_event_output< float >(tag));
}

void PreselectionModule::PrintInputs() {
  std::cout << "****************************************" << std::endl;
  std::cout << "             "+NameModule+"             " << std::endl;
  std::cout << "----------------------------------------" << std::endl;
  for (auto x : MS) std::cout << x.first << std::string( 18-x.first.size(), ' ' ) << x.second << '\n';
  for (auto x : MB) std::cout << x.first << std::string( 18-x.first.size(), ' ' ) << BoolToString(x.second) << '\n';
  std::cout << "****************************************\n" << std::endl;
}

void PreselectionModule::book_histograms(uhh2::Context& ctx) {
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
    mytag = "ZprimeCandidate_" + tag; book_HFolder(mytag, new HiggsToWWHists(ctx,mytag));
  }
}

void PreselectionModule::fill_histograms(uhh2::Event& event, string tag){
  std::vector<string> mytags = {"event_", "gen_", "nTopJet_", "nJet_", "ele_", "muon_", "diLepton_", "BTagEff_", "ZprimeCandidate_"};
  for (auto& mytag : mytags) HFolder(mytag+ tag)->fill(event);
}

/*
█  ██████  ██████  ███    ██ ███████ ████████ ██████  ██    ██  ██████ ████████  ██████  ██████
█ ██      ██    ██ ████   ██ ██         ██    ██   ██ ██    ██ ██         ██    ██    ██ ██   ██
█ ██      ██    ██ ██ ██  ██ ███████    ██    ██████  ██    ██ ██         ██    ██    ██ ██████
█ ██      ██    ██ ██  ██ ██      ██    ██    ██   ██ ██    ██ ██         ██    ██    ██ ██   ██
█  ██████  ██████  ██   ████ ███████    ██    ██   ██  ██████   ██████    ██     ██████  ██   ██
*/

PreselectionModule::PreselectionModule(uhh2::Context& ctx){

  // Set up variables
  MS["year"]              = ctx.get("year");
  MB["is_mc"]             = ctx.get("dataset_type") == "MC";
  MB["isPuppi"]           = string2bool(ctx.get("isPuppi"));
  MB["isCHS"]             = string2bool(ctx.get("isCHS"));
  MB["isHOTVR"]           = string2bool(ctx.get("isHOTVR"));
  MB["muonchannel"]       = string2bool(ctx.get("muonchannel"));
  MB["electronchannel"]   = string2bool(ctx.get("electronchannel"));
  MB["invisiblechannel"]   = string2bool(ctx.get("invisiblechannel"));
  MS["SysType_PU"]        = ctx.get("SysType_PU");
  MB["lumisel"]           = string2bool(ctx.get("lumisel"));
  MB["mclumiweight"]      = string2bool(ctx.get("mclumiweight"));
  MB["mcpileupreweight"]  = string2bool(ctx.get("mcpileupreweight"));
  MB["eleid"]             = string2bool(ctx.get("eleid"));
  MB["muid"]              = string2bool(ctx.get("muid"));
  MB["tauid"]             = string2bool(ctx.get("tauid"));
  MB["metfilters"]        = string2bool(ctx.get("metfilters"));

  if ((MB["isPuppi"] && MB["isCHS"]) || (MB["isPuppi"] && MB["isHOTVR"]) || (MB["isCHS"] && MB["isHOTVR"]) ) throw std::runtime_error("In PreselectionModule.cxx: Choose exactly one jet collection.");
  if ((MB["muonchannel"] && MB["electronchannel"]) || (MB["muonchannel"] && MB["invisiblechannel"]) || (MB["electronchannel"] && MB["invisiblechannel"])) throw std::runtime_error("In PreselectionModule.cxx: Choose exactly one lepton channel.");

  MS["leptons"] = MB["muonchannel"]? "muons": (MB["electronchannel"]? "electrons": "");

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

  const MuonId muoId = AndId<Muon>(MuonID(Muon::CutBasedIdTrkHighPt), PtEtaCut(min_lepton_pt, min_lepton_eta));
  // const ElectronId eleId = AndId<Electron>(ElectronID_Fall17_loose, PtEtaSCCut(min_lepton_pt, min_lepton_eta));
  const ElectronId eleId = AndId<Electron>(ElectronTagID(Electron::cutBasedElectronID_Fall17_94X_V2_loose), PtEtaSCCut(min_lepton_pt, min_lepton_eta));
  // cutBasedElectronID_Summer16_80X_V1_loose

  const JetId jetId = AndId<Jet> (JetPFID(JETwp), PtEtaCut(min_jet_pt, min_lepton_eta));
  const TopJetId topjetId = AndId<TopJet> (JetPFID(JETwp), PtEtaCut(min_topjet_pt, min_lepton_eta), NoLeptonInJet("all", eleId, muoId, MB["isHOTVR"]? 0.8: -1));

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
  if (MS["year"] != "2016") metfilters_selection->add<EcalBadCalibSelection>("EcalBadCalibSelection"); /*TODO check 2016*/ // Use this instead of Flag_ecalBadCalibFilter, uses ecalBadCalibReducedMINIAODFilter in ntuple_generator
  if (MS["year"] != "2016") metfilters_selection->add<TriggerSelection>("BadPFMuonFilter", "Flag_BadPFMuonFilter"); /*TODO check 2016, maybe Extra_BadPFMuonFilter */
  if (!MB["is_mc"]) metfilters_selection->add<TriggerSelection>("eeBadScFilter", "Flag_eeBadScFilter"); /* TODO Not recommended for MC, but do check */
  /* metfilters_selection->add<TriggerSelection>("BadChargedCandidateFilter", "Flag_BadChargedCandidateFilter"); TODO Not recommended, under review.Separate module in ntuple_generator for 2016v2*/

  GJC.reset( new GenericJetCleaner(ctx, MS["jetLabel"],    false, jetId, topjetId, muoId, eleId));
  GTJC.reset(new GenericJetCleaner(ctx, MS["topjetLabel"], true,  jetId, topjetId, muoId, eleId));

  NBoostedJetSel.reset(new NTopJetSelection(1,-1,topjetId,h_topjets));

  if (MB["muonchannel"]) {
    NLeptonSel.reset(new NMuonSelection(min_leptons)); // min_leptons=2
    NoLeptonSel.reset(new NElectronSelection(1));
  }
  if (MB["electronchannel"]) {
    NLeptonSel.reset(new NElectronSelection(min_leptons)); // min_leptons=2
    NoLeptonSel.reset(new NMuonSelection(1));
  }
  if (MB["invisiblechannel"]) {
    NLeptonSel.reset(new NElectronSelection(0,0));
    NoLeptonSel.reset(new NMuonSelection(1));
  }

  DeltaRDiLepton_selection.reset(new DeltaRDiLepton(min_DR_dilep, max_DR_dilep, MS["leptons"]));
  JetDiLeptonPhiAngularSel.reset(new JetDiLeptonPhiAngularSelection(min_dilep_pt, min_jet_dilep_delta_phi, max_jet_dilep_delta_phi, MS["leptons"], h_topjets));
  VetoLeptonSel.reset(new VetoSelection(NoLeptonSel));

  if (MB["invisiblechannel"]) {
      DeltaRDiLepton_selection.reset(new AndSelection(ctx));
      JetDiLeptonPhiAngularSel.reset(new AndSelection(ctx));
  }
}


/*
█ ██████  ██████   ██████   ██████ ███████ ███████ ███████
█ ██   ██ ██   ██ ██    ██ ██      ██      ██      ██
█ ██████  ██████  ██    ██ ██      █████   ███████ ███████
█ ██      ██   ██ ██    ██ ██      ██           ██      ██
█ ██      ██   ██  ██████   ██████ ███████ ███████ ███████
*/

bool PreselectionModule::process(uhh2::Event& event) {

  if ((event.year).find(MS["year"])==std::string::npos) throw std::runtime_error("In "+NameModule+".cxx: You are running on "+event.year+" sample with a "+MS["year"]+" year config file. Fix this.");

  auto weight_gen = event.weight;
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

  if(!VetoLeptonSel->passes(event)) return false;
  fill_histograms(event, "Veto");

  if(!NLeptonSel->passes(event)) return false;
  fill_histograms(event, "NLeptonSel");

  if(!NBoostedJetSel->passes(event)) return false;
  fill_histograms(event, "NBoostedJet");

  if(!DeltaRDiLepton_selection->passes(event)) return false;
  fill_histograms(event, "DeltaRDiLepton");

  if(!JetDiLeptonPhiAngularSel->passes(event)) return false;
  fill_histograms(event, "JetDiLeptonPhiAngular");

  return true;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the PreselectionModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(PreselectionModule)
