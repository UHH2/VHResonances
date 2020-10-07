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

#include "UHH2/VHResonances/include/ModuleBase.h"
#include "UHH2/VHResonances/include/constants.hpp"
#include "UHH2/VHResonances/include/ExtJetHists.h"
#include "UHH2/VHResonances/include/ExtJetHistsWithConditions.h"
#include "UHH2/VHResonances/include/GenMatchHists.h"
#include "UHH2/VHResonances/include/GenLevelJetMatch.h"
#include "UHH2/VHResonances/include/HiggsToWWSelection.h"
#include "UHH2/VHResonances/include/DiLeptonHists.h"
#include "UHH2/VHResonances/include/HiggsToWWModules.h"
#include "UHH2/VHResonances/include/HiggsToWWHists.h"

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
  std::string GetSystName(const std::string& syst, const std::string& var);

protected:

  // Define variables

  std::vector<std::string> histogram_tags = {"Selection",
  "btag_DeepBoosted_H4qvsQCD_SR",           "btag_DeepBoosted_H4qvsQCD_CR",
  "btag_DeepBoosted_H4qvsQCDmassdep_SR", "btag_DeepBoosted_H4qvsQCDmassdep_CR",
  "btag_DeepBoosted_H4qvsQCDmassdep_2_SR", "btag_DeepBoosted_H4qvsQCDmassdep_2_CR",
  "btag_DeepBoosted_H4qvsQCDmassdep_3_SR", "btag_DeepBoosted_H4qvsQCDmassdep_3_CR",
  "btag_DeepBoosted_H4qvsQCDmassdep_bb_SR", "btag_DeepBoosted_H4qvsQCDmassdep_bb_CR",
  "btag_DeepBoosted_H4qvsQCDmassdep_cc_SR", "btag_DeepBoosted_H4qvsQCDmassdep_cc_CR",
  "btag_DeepBoosted_H4qvsQCDmassdep_gg_SR", "btag_DeepBoosted_H4qvsQCDmassdep_gg_CR",
  "btag_DeepBoosted_H4qvsQCDmassdep_cc_2_SR", "btag_DeepBoosted_H4qvsQCDmassdep_cc_2_CR",
  "btag_DeepBoosted_H4qvsQCDmassdep_gg_2_SR", "btag_DeepBoosted_H4qvsQCDmassdep_gg_2_CR",
  "btag_DeepBoosted_H4qvsQCDmassdep_cc_3_SR", "btag_DeepBoosted_H4qvsQCDmassdep_cc_3_CR",
  "btag_DeepBoosted_H4qvsQCDmassdep_gg_3_SR", "btag_DeepBoosted_H4qvsQCDmassdep_gg_3_CR",
  "tau42_SR", "tau42_CR"};

  std::vector<std::string> weight_tags = {"weight_lumi", "weight_GLP", "HDecay", "ZDecay", "ZprimeDecay"};
  std::vector<std::string> Systematics =  {"pu", "btag", "prefiring", "id", "tracking", "trigger", "reco"};
  std::vector<std::string> Variations = {"", "up", "down"};


  std::unordered_map<std::string, std::string> MS;
  std::unordered_map<std::string, bool> MB;

  Event::Handle<std::vector<TopJet> > h_topjets;
  Event::Handle<std::vector<ZprimeCandidate> > h_ZprimeCandidates;

  // Define selections

  std::unique_ptr<Selection> btag_DeepBoosted_H4qvsQCD_SR_selection,            btag_DeepBoosted_H4qvsQCD_CR_selection;
  std::unique_ptr<Selection> btag_DeepBoosted_H4qvsQCDmassdep_SR_selection,  btag_DeepBoosted_H4qvsQCDmassdep_CR_selection;
  std::unique_ptr<Selection> tau42_SR_selection, tau42_CR_selection;
};


void SignalRegionModule::PrintInputs() {
  std::cout << "****************************************" << std::endl;
  std::cout << "           SignalRegionModule           " << std::endl;
  std::cout << "----------------------------------------" << std::endl;
  for (auto x : MS) std::cout << x.first << std::string( 18-x.first.size(), ' ' ) << x.second << '\n';
  for (auto x : MB) std::cout << x.first << std::string( 18-x.first.size(), ' ' ) << BoolToString(x.second) << '\n';
  std::cout << "****************************************\n" << std::endl;
}

std::string SignalRegionModule::GetSystName(const std::string& syst, const std::string& var) {
  string tag = "weight_";
  if (FindInString("id", syst) || FindInString("tracking", syst) || FindInString("trigger", syst) || FindInString("reco", syst)) {
    tag += MB["muonchannel"]? "sfmu_": (MB["electronchannel"]? "sfelec_": "");
  }
  tag += syst;
  if (var!="") tag += "_"+var;
  if (FindInString("prefiring", syst)) {
    std::string Var = var;
    Var[0] = toupper(Var[0]);
    tag = "prefiringWeight"+Var;
  }
  return tag;
}

void SignalRegionModule::book_handles(uhh2::Context& ctx) {
  for(const auto & tag : weight_tags) {
    book_WFolder(tag+"_in",  new Event::Handle< float >, ctx.declare_event_input< float >(tag));
    book_WFolder(tag+"_out", new Event::Handle< float >, ctx.declare_event_output< float >(tag));
  }
  if (!MB["is_mc"]) return;

  for(const auto & syst : Systematics) {
    for(const auto & var : Variations) {
      string tag = GetSystName(syst, var);
      book_WFolder(tag+"_in",  new Event::Handle< float >, ctx.declare_event_input< float >(tag));
      book_WFolder(tag+"_out", new Event::Handle< float >, ctx.declare_event_output< float >(tag));
    }
  }


}


void SignalRegionModule::export_weights(uhh2::Event& event) {
  for(const auto & tag : weight_tags) event.set(WFolder(tag+"_out"), event.get(WFolder(tag+"_in")));
  if (!MB["is_mc"]) return;

  for(const auto & syst : Systematics) {
    for(const auto & var : Variations) {
      string tag = GetSystName(syst, var);
      event.set(WFolder(tag+"_out"), event.get(WFolder(tag+"_in")));
    }
  }
}


void SignalRegionModule::book_histograms(uhh2::Context& ctx){
  for(const auto & tag : histogram_tags){
    string mytag;
    mytag = "diLepton_"           + tag; book_HFolder(mytag, new DiLeptonHists(ctx,mytag, "", MS["topjetLabel"]));
    mytag = "ZprimeCandidate_"    + tag; book_HFolder(mytag, new HiggsToWWHists(ctx,mytag));
    mytag = "nTopJet_"            + tag; book_HFolder(mytag, new ExtJetHists(ctx,mytag, MS["topjetLabel"] ));

    for (std::string& syst :  Systematics){
      for (std::string& var : Variations){
        if (var=="") continue;
        mytag = "ZprimeCandidate_"+syst+"_"+var+"_"+tag; book_HFolder(mytag, new HiggsToWWHistsSlim(ctx,mytag));
      }
    }

  }
}

void SignalRegionModule::fill_histograms(uhh2::Event& event, string tag){
  std::vector<string> mytags = {"nTopJet_", "ZprimeCandidate_","diLepton_"};
  for (auto& mytag : mytags) HFolder(mytag+ tag)->fill(event);
  string mytag;
  float save_weight = event.weight;
  for (std::string& syst :  Systematics){
    for (std::string& var : Variations){
      if (var=="") continue;
      if (!event.isRealData) event.weight = save_weight*event.get(WFolder(GetSystName(syst, var)+"_in"))/event.get(WFolder(GetSystName(syst, "")+"_in"));
      mytag = "ZprimeCandidate_"+syst+"_"+var+"_"+tag; HFolder(mytag)->fill(event);
    }
  }

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

  MS["dataset_version"]  = ctx.get("dataset_version");
  MB["is_mc"]            = ctx.get("dataset_type") == "MC";
  MB["isPuppi"]          = string2bool(ctx.get("isPuppi"));
  MB["isCHS"]            = string2bool(ctx.get("isCHS"));
  MB["isHOTVR"]          = string2bool(ctx.get("isHOTVR"));
  MB["muonchannel"]      = string2bool(ctx.get("muonchannel"));
  MB["electronchannel"]  = string2bool(ctx.get("electronchannel"));
  MB["invisiblechannel"] = string2bool(ctx.get("invisiblechannel"));

  MB["is_ZH"]   = MS["dataset_version"].find("MC_ZprimeToZH")!=std::string::npos;
  MB["is_bb"]   = MS["dataset_version"].find("Tobb")!=std::string::npos;
  MB["is_WW"]   = MS["dataset_version"].find("ToWW")!=std::string::npos;
  MB["is_else"] = MS["dataset_version"].find("_extra")!=std::string::npos;

  if ((MB["isPuppi"] && MB["isCHS"]) || (MB["isPuppi"] && MB["isHOTVR"]) || (MB["isCHS"] && MB["isHOTVR"]) ) throw std::runtime_error("In SelectionModule.cxx: Choose exactly one jet collection.");
  if ((MB["muonchannel"] && MB["electronchannel"]) || (MB["muonchannel"] && MB["invisiblechannel"]) || (MB["electronchannel"] && MB["invisiblechannel"])) throw std::runtime_error("In SelectionModule.cxx: Choose exactly one lepton channel.");

  MS["topjetLabel"] = MB["isCHS"]? "topjets": (MB["isPuppi"]? "toppuppijets": (MB["isHOTVR"]? "hotvrPuppi": ""));

  if(MB["electronchannel"]) Systematics.erase(Systematics.begin()+FindInVector(Systematics,"tracking"));

  // Set up histograms:
  book_histograms(ctx);
  book_handles(ctx);
  PrintInputs();

  // Set up selections
  h_topjets = ctx.get_handle<std::vector<TopJet>>(MS["topjetLabel"]);
  h_ZprimeCandidates = ctx.declare_event_input<std::vector<ZprimeCandidate>>("ZprimeCandidate");

  btag_DeepBoosted_H4qvsQCD_SR_selection.reset(new TaggerCut(0.8, 1.0, -1, "btag_DeepBoosted_H4qvsQCD", h_ZprimeCandidates));
  btag_DeepBoosted_H4qvsQCD_CR_selection.reset(new TaggerCut(0.0, 0.2, -1, "btag_DeepBoosted_H4qvsQCD", h_ZprimeCandidates));
  btag_DeepBoosted_H4qvsQCDmassdep_SR_selection.reset(new TaggerCut(0, 1,  MassDepentdentCut_value, "btag_DeepBoosted_H4qvsQCD", h_ZprimeCandidates));
  btag_DeepBoosted_H4qvsQCDmassdep_CR_selection.reset(new TaggerCut(0, -1, MassDepentdentCut_value, "btag_DeepBoosted_H4qvsQCD", h_ZprimeCandidates));

  std::string groom = MB["isHOTVR"]? "_groomed": "";
  tau42_SR_selection.reset(new TaggerCut(0.0, 0.5, -1, "tau42"+groom, h_ZprimeCandidates));
  tau42_CR_selection.reset(new TaggerCut(0.5, 1.0, -1, "tau42"+groom, h_ZprimeCandidates));

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

  if(btag_DeepBoosted_H4qvsQCD_SR_selection->passes(event))   fill_histograms(event, "btag_DeepBoosted_H4qvsQCD_SR");
  if(btag_DeepBoosted_H4qvsQCD_CR_selection->passes(event))   fill_histograms(event, "btag_DeepBoosted_H4qvsQCD_CR");

  ZprimeCandidate cand = event.get(h_ZprimeCandidates)[0];

  double delta_R_ll = deltaR(cand.leptons()[0], cand.leptons()[1]);
  double jet_pt = cand.H().pt();
  double ptgcut = (-5e-04*jet_pt*2);
  double ptccut = (-3.5e-04*jet_pt*2+0.7);
  double cc_score = cand.H().btag_DeepBoosted_probHcc();
  double bb_score = cand.H().btag_DeepBoosted_probHbb();
  double qcd_score = cand.H().btag_DeepBoosted_probQCDb()+cand.H().btag_DeepBoosted_probQCDbb()+cand.H().btag_DeepBoosted_probQCDc()+cand.H().btag_DeepBoosted_probQCDcc()+cand.H().btag_DeepBoosted_probQCDothers();
  double Hcc = cc_score/(cc_score+qcd_score);
  double Hbb = bb_score/(bb_score+qcd_score);
  double sj1_g = -1; double sj2_g = -1;
  if (cand.H().subjets().size()>0) sj1_g = cand.H().subjets()[0].btag_DeepFlavour_g();
  if (cand.H().subjets().size()>1) sj2_g = cand.H().subjets()[1].btag_DeepFlavour_g();

  if(btag_DeepBoosted_H4qvsQCDmassdep_SR_selection->passes(event)) {
    fill_histograms(event, "btag_DeepBoosted_H4qvsQCDmassdep_3_SR");
    if (delta_R_ll>0.4) return false;
    fill_histograms(event, "btag_DeepBoosted_H4qvsQCDmassdep_2_SR");
    if (cand.H().pt()<350 || cand.Z().pt()<350 ) return false;
    fill_histograms(event, "btag_DeepBoosted_H4qvsQCDmassdep_SR");
    if (Hbb>0.8)                    fill_histograms(event, "btag_DeepBoosted_H4qvsQCDmassdep_bb_SR");
    if (Hcc>0.8)                    fill_histograms(event, "btag_DeepBoosted_H4qvsQCDmassdep_cc_SR");
    if (sj1_g<0.1 || sj2_g<0.1)     fill_histograms(event, "btag_DeepBoosted_H4qvsQCDmassdep_gg_SR");
    if (Hcc>0.1)                    fill_histograms(event, "btag_DeepBoosted_H4qvsQCDmassdep_cc_2_SR");
    if (sj1_g<0.8 || sj2_g<0.8)     fill_histograms(event, "btag_DeepBoosted_H4qvsQCDmassdep_gg_2_SR");
    if (Hcc>ptccut)                 fill_histograms(event, "btag_DeepBoosted_H4qvsQCDmassdep_cc_3_SR");
    if (sj1_g<ptgcut||sj2_g<ptgcut) fill_histograms(event, "btag_DeepBoosted_H4qvsQCDmassdep_gg_3_SR");
  }
  if(btag_DeepBoosted_H4qvsQCDmassdep_CR_selection->passes(event)) {
    fill_histograms(event, "btag_DeepBoosted_H4qvsQCDmassdep_CR");
    if (delta_R_ll>0.4) return false;
    fill_histograms(event, "btag_DeepBoosted_H4qvsQCDmassdep_2_CR");
    if (cand.H().pt()<350 || cand.Z().pt()<350 ) return false;
    fill_histograms(event, "btag_DeepBoosted_H4qvsQCDmassdep_3_CR");
    if (!(Hbb>0.8))                   fill_histograms(event, "btag_DeepBoosted_H4qvsQCDmassdep_bb_CR");
    if (!(Hcc>0.8))                   fill_histograms(event, "btag_DeepBoosted_H4qvsQCDmassdep_cc_CR");
    if (!(sj1_g<0.1 || sj2_g<0.1))    fill_histograms(event, "btag_DeepBoosted_H4qvsQCDmassdep_gg_CR");
    if (!(Hcc>0.1))                   fill_histograms(event, "btag_DeepBoosted_H4qvsQCDmassdep_cc_2_CR");
    if (!(sj1_g<0.8 || sj2_g<0.8))    fill_histograms(event, "btag_DeepBoosted_H4qvsQCDmassdep_gg_2_CR");
    if (!(Hcc>ptccut))                fill_histograms(event, "btag_DeepBoosted_H4qvsQCDmassdep_cc_3_CR");
    if (!(sj1_g<ptgcut||sj2_g<ptgcut))fill_histograms(event, "btag_DeepBoosted_H4qvsQCDmassdep_gg_3_CR");
  }

  if(tau42_SR_selection->passes(event)) fill_histograms(event, "tau42_SR");
  if(tau42_CR_selection->passes(event)) fill_histograms(event, "tau42_CR");

  export_weights(event);

  return true;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the SignalRegionModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(SignalRegionModule)
