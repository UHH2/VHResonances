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

class SelectionModule: public ModuleBASE {

public:

  explicit SelectionModule(uhh2::Context&);
  virtual bool process(uhh2::Event&) override;
  void book_histograms(uhh2::Context&);
  void fill_histograms(uhh2::Event&, string);
  void book_handles(uhh2::Context&);
  void export_weights(uhh2::Event&);
  void PrintInputs();

protected:

  // Define variables
  std::string NameModule = "SelectionModule";
  std::vector<std::string> histogram_tags = {"Preselection", "Trigger", "MuonScale", "ZprimeReco", "ZprimeSelection", "PTMassCut", "ScaleFactors"};
  std::vector<std::string> weight_tags = {"weight_lumi", "weight_GLP", "weight_pu", "weight_pu_up", "weight_pu_down", "HDecay", "ZDecay", "ZprimeDecay", "weight_btag","weight_btag_up", "weight_btag_down"};

  std::unordered_map<std::string, std::string> MS;
  std::unordered_map<std::string, bool> MB;

  Event::Handle<std::vector<Jet> > h_jets;
  Event::Handle<std::vector<TopJet> > h_topjets;
  Event::Handle<std::vector<ZprimeCandidate> > h_ZprimeCandidates;

  Event::Handle<std::vector<tensorflow::Tensor> > h_image;

  // Define selections

  std::unique_ptr<Selection> PTMassCut_selection;
  std::unique_ptr<AnalysisModule> ZprimeCandidateReconstruction_module;
  std::unique_ptr<AnalysisModule> CollectionProducer_module;
  std::unordered_map<std::string, std::unique_ptr<Selection>> Trigger_selection;
  std::unordered_map<std::string, std::unique_ptr<AnalysisModule>> ScaleFactors_module;
  std::unique_ptr<AnalysisModule> MuonScaleVariations_module;

};


void SelectionModule::PrintInputs() {
  std::cout << "****************************************" << std::endl;
  std::cout << "             "+NameModule+"             " << std::endl;
  std::cout << "----------------------------------------" << std::endl;
  for (auto x : MS) std::cout << x.first << std::string( 18-x.first.size(), ' ' ) << x.second << '\n';
  for (auto x : MB) std::cout << x.first << std::string( 18-x.first.size(), ' ' ) << BoolToString(x.second) << '\n';
  std::cout << "****************************************\n" << std::endl;
}

void SelectionModule::book_handles(uhh2::Context& ctx) {
  for(const auto & tag : weight_tags) {
    if (!MB["is_mc"] && tag.find("weight_pu")!=std::string::npos) continue;
    book_WFolder(tag+"_in",  new Event::Handle< float >, (tag.find("btag")!=std::string::npos)? ctx.get_handle< float >(tag) : ctx.declare_event_input< float >(tag));
    book_WFolder(tag+"_out", new Event::Handle< float >, ctx.declare_event_output< float >(tag));
  }
}


void SelectionModule::export_weights(uhh2::Event& event) {
  for(const auto & tag : weight_tags) {
    if (!MB["is_mc"] && tag.find("weight_pu")!=std::string::npos) continue;
    event.set(WFolder(tag+"_out"), event.get(WFolder(tag+"_in")));
  }
  event.set(WFolder("weight_GLP_out"), event.weight);
}


void SelectionModule::book_histograms(uhh2::Context& ctx){
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

void SelectionModule::fill_histograms(uhh2::Event& event, string tag){
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

SelectionModule::SelectionModule(uhh2::Context& ctx){

  // Set up variables
  MS["year"]              = ctx.get("year");
  MB["is_mc"]             = ctx.get("dataset_type") == "MC";
  MB["isPuppi"]           = string2bool(ctx.get("isPuppi"));
  MB["isCHS"]             = string2bool(ctx.get("isCHS"));
  MB["isHOTVR"]           = string2bool(ctx.get("isHOTVR"));
  MB["muonchannel"]       = string2bool(ctx.get("muonchannel"));
  MB["electronchannel"]   = string2bool(ctx.get("electronchannel"));
  MB["invisiblechannel"]   = string2bool(ctx.get("invisiblechannel"));

  if ((MB["isPuppi"] && MB["isCHS"]) || (MB["isPuppi"] && MB["isHOTVR"]) || (MB["isCHS"] && MB["isHOTVR"]) ) throw std::runtime_error("In SelectionModule.cxx: Choose exactly one jet collection.");
  if ((MB["muonchannel"] && MB["electronchannel"]) || (MB["muonchannel"] && MB["invisiblechannel"]) || (MB["electronchannel"] && MB["invisiblechannel"])) throw std::runtime_error("In SelectionModule.cxx: Choose exactly one lepton channel.");

  MS["jetLabel"]    = MB["isCHS"]? "jets":    (MB["isPuppi"]? "jetsAk4Puppi": (MB["isHOTVR"]? "jetsAk4Puppi": ""));
  MS["topjetLabel"] = MB["isCHS"]? "topjets": (MB["isPuppi"]? "toppuppijets": (MB["isHOTVR"]? "hotvrPuppi": ""));

  // Set up histograms:

  book_histograms(ctx);
  book_handles(ctx);
  PrintInputs();

  // Set up selections

  h_jets = ctx.get_handle<std::vector<Jet>>(MS["jetLabel"]);
  h_topjets = ctx.get_handle<std::vector<TopJet>>(MS["topjetLabel"]);
  h_ZprimeCandidates = ctx.get_handle<std::vector<ZprimeCandidate>>("ZprimeCandidate");
  h_image = ctx.get_handle<std::vector<tensorflow::Tensor>>("Images");

  MS["leptons"] = MB["muonchannel"]? "muons": (MB["electronchannel"]? "electrons": (MB["invisiblechannel"]? "invisible": ""));
  for (auto& t : Trigger_run_validity.at(MS["year"])) {
    if (MB["muonchannel"] && t.first.find("Mu")==std::string::npos ) continue;
    if (MB["muonchannel"] &&  t.first.find("NoMu")!=std::string::npos ) continue;
    if (MB["electronchannel"] && t.first.find("Ele")==std::string::npos ) continue;
    if (MB["invisiblechannel"] && t.first.find("MET")==std::string::npos ) continue;
    Trigger_selection[t.first].reset(new TriggerSelection( t.first ));
  }

  //Scale factors
  MuonScaleVariations_module.reset(new MuonScaleVariations(ctx));

  ScaleFactors_module["BTag"].reset(new MCBTagScaleFactor(ctx, BTag_algo, BTag_wp, MS["topjetLabel"], "nominal", "lt"));
  ScaleFactors_module["SFs"].reset(new ScaleFactorsManager(ctx, h_ZprimeCandidates));

  ZprimeCandidateReconstruction_module.reset(new ZprimeCandidateReconstruction(ctx, min_dilep_pt, min_DR_dilep, max_DR_dilep, min_jet_dilep_delta_phi, max_jet_dilep_delta_phi, MS["leptons"], MS["topjetLabel"]));
  CollectionProducer_module.reset(new CollectionProducer<ZprimeCandidate>( ctx, "ZprimeCandidate", "ZprimeCandidate", (ZprimeCandidate_ID)ZprimeCandidateID(h_ZprimeCandidates)));
  PTMassCut_selection.reset(new PTMassCut(min_Z_pt_ZH_mass, h_ZprimeCandidates));

}


/*
█ ██████  ██████   ██████   ██████ ███████ ███████ ███████
█ ██   ██ ██   ██ ██    ██ ██      ██      ██      ██
█ ██████  ██████  ██    ██ ██      █████   ███████ ███████
█ ██      ██   ██ ██    ██ ██      ██           ██      ██
█ ██      ██   ██  ██████   ██████ ███████ ███████ ███████
*/

bool SelectionModule::process(uhh2::Event& event) {

  if ((event.year).find(MS["year"])==std::string::npos) throw std::runtime_error("In "+NameModule+".cxx: You are running on "+event.year+" sample with a "+MS["year"]+" year config file. Fix this.");

  event.weight = event.get(WFolder("weight_GLP_in"));

  fill_histograms(event, "Preselection");

  bool pass_triggers_OR = false;

  for (auto& el : Trigger_selection) {
    if (event.isRealData && (event.run < Trigger_run_validity.at(MS["year"]).at(el.first).first || event.run > Trigger_run_validity.at(MS["year"]).at(el.first).second) ) continue;
    pass_triggers_OR += el.second->passes(event);
    if (pass_triggers_OR) break;
  }
  if (!pass_triggers_OR) return false;

  fill_histograms(event, "Trigger");

  MuonScaleVariations_module->process(event);
  fill_histograms(event, "MuonScale");

  // cut delta Phi between MET and all TopJets at min_Dphi_MET (2.0)
  if (MB["invisiblechannel"]){
    const auto & jets = event.get(h_topjets);
    for(const auto & jet: jets){
      double Dphi = fabs(deltaPhi(jet, *event.met));
      if (Dphi < min_Dphi_MET) return false;
    }
  }

  ZprimeCandidateReconstruction_module->process(event);
  fill_histograms(event, "ZprimeReco");

  CollectionProducer_module->process(event);
  if(event.get(h_ZprimeCandidates).size()<1) return false;
  if(event.get(h_ZprimeCandidates).size()>1) {
    for(const auto & cand: event.get(h_ZprimeCandidates)) {
      if (cand.discriminator("SDmass") > cand.H().v4().M()) return false;
      if (cand.discriminator("SDmass") < 60) return false;
    }
  }
  fill_histograms(event, "ZprimeSelection");

  if(!PTMassCut_selection->passes(event)) return false;
  fill_histograms(event, "PTMassCut");

  for (auto& el : ScaleFactors_module) el.second->process(event);
  fill_histograms(event, "ScaleFactors");

  export_weights(event);

  return true;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the SelectionModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(SelectionModule)
