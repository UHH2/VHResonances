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

  std::string NameModule = "SignalRegionModule";
  std::vector<std::string> histogram_tags = {"Selection", "ExtraCleaning",
  "DeepAk8_ZHccvsQCD_MD_SR",                  "DeepAk8_ZHccvsQCD_MD_CR",
  "DeepAk8_ZHccvsQCD_MD2_SR",                 "DeepAk8_ZHccvsQCD_MD2_CR",
  "DeepAk8_HccvsQCD_SR",                      "DeepAk8_HccvsQCD_CR",
  "DeepAk8_H4qvsQCD_SR",                      "DeepAk8_H4qvsQCD_CR",
  "DeepAk8_H4qvsQCD_massdep_SR",              "DeepAk8_H4qvsQCD_massdep_CR",
  "DeepAk8_H4qvsQCD_massdep_HccvsQCD_SR",     "DeepAk8_H4qvsQCD_massdep_HccvsQCD_CR",
  "PN_ZHccvsQCD_MD_SR",                       "PN_ZHccvsQCD_MD_CR",
  "PN_ZHccvsQCD_MD2_SR",                      "PN_ZHccvsQCD_MD2_CR",
  "PN_HccvsQCD_SR",                           "PN_HccvsQCD_CR",
  "PN_H4qvsQCD_SR",                           "PN_H4qvsQCD_CR",
  "tau42_SR", "tau42_CR","tau21_SR", "tau21_CR"};

  std::vector<std::string> weight_tags = {"weight_lumi", "weight_GLP", "HDecay", "ZDecay", "ZprimeDecay"};
  // std::vector<std::string> Systematics = {"pu", "btag", "prefiring", "id", "isolation", "tracking", "trigger", "reco", "taggerSF", "murmuf", "NNPDF"};
  std::vector<std::string> Systematics = {"pu", "prefiring", "id", "isolation", "trigger", "reco", "taggerSF", "murmuf", "NNPDF"};

  std::vector<std::string> Var_murmuf = {"upup", "upnone", "noneup", "nonedown", "downnone", "downdown"};
  std::vector<std::string> Variations = {"", "up", "down"};
  int PDF_variations = 100; int pdfindex_shift = 9;
  std::unordered_map<std::string, std::string> PFDs = {
    {"NNPDF", "NNPDF31_lo_as_0130"},
    // {"NNPDF", "NNPDF31_nnlo_as_0118_nf_4_mc_hessian"},
    // To keep as memo!
    // {"NNPDF31_lo_as_0130", "NNPDF31_lo_as_0130"}, used for PDF reweight in invisiblechannel channel
    // {"PDF4LHC15_nnlo_100", "PDF4LHC15_nnlo_100"},
    // {"NNPDF31_nnlo_as_0118_mc_hessian_pdfas", "NNPDF31_nnlo_as_0118_mc_hessian_pdfas"},// same results as NNPDF31_nnlo_as_0118_nf_4_mc_hessian
    // {"NNPDF31_nnlo_hessian_pdfas", "NNPDF31_nnlo_hessian_pdfas"},//This is the default stored in ntuples for the lepton channel
    // {"PDF4LHC15_nlo_mc_pdfas", "PDF4LHC15_nlo_mc_pdfas"},
    // {"PDF4LHC15_nlo_mc", "PDF4LHC15_nlo_mc"},
    // {"PDF4LHC15_nnlo_mc_pdfas", "PDF4LHC15_nnlo_mc_pdfas"},
    // {"PDF4LHC15_nnlo_mc", "PDF4LHC15_nnlo_mc"},
  };

  std::unordered_map<std::string, std::string> MS;
  std::unordered_map<std::string, bool> MB;

  Event::Handle<std::vector<TopJet> > h_topjets;
  Event::Handle<std::vector<ZprimeCandidate> > h_ZprimeCandidates;

  std::unordered_map<std::string, std::vector<double> > PDF_weights;
  std::unordered_map<std::string, std::unique_ptr<PDFWeights> > m_pdfweights;

  // Define selections

  std::unique_ptr<Selection> DeepAk8_H4qvsQCD_massdep_SR_selection;
  std::unique_ptr<Selection> tau42_SR_selection, tau21_SR_selection;

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
  if (FindInString("id", syst) || FindInString("isolation", syst) || FindInString("tracking", syst) || FindInString("trigger", syst) || FindInString("reco", syst)) {
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
    if (FindInString("PDF", syst)) continue;
    for(const auto & var : FindInString("murmuf", syst) ? Var_murmuf: Variations) {
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
    if (FindInString("PDF", syst)) continue;
    for(const auto & var : FindInString("murmuf", syst) ? Var_murmuf: Variations) {
      string tag = GetSystName(syst, var);
      event.set(WFolder(tag+"_out"), event.get(WFolder(tag+"_in")));
    }
  }
}


void SignalRegionModule::book_histograms(uhh2::Context& ctx){
  for(const auto & tag : histogram_tags){
    string mytag;
    mytag = "event_"              + tag; book_HFolder(mytag, new EventHists(ctx,mytag));
    mytag = "diLepton_"           + tag; book_HFolder(mytag, new DiLeptonHists(ctx,mytag, "", MS["topjetLabel"]));
    mytag = "ZprimeCandidate_"    + tag; book_HFolder(mytag, new HiggsToWWHists(ctx,mytag));
    mytag = "nTopJet_"            + tag; book_HFolder(mytag, new ExtJetHists(ctx,mytag, MS["topjetLabel"] ));
    if (FindInString("_CR", tag)) continue;
    if (!FindInString("ZHccvsQCD", tag)) continue;
    for (std::string& syst :  Systematics){
      if (FindInString("PDF", syst)) {
        for(int i=0; i<PDF_variations; i++){
          std::string var = to_string(i);
          mytag = "ZprimeCandidate_"+syst+"_"+var+"_"+tag; book_HFolder(mytag, new HiggsToWWHistsSlim(ctx,mytag));
        }
      } else {
        for (std::string& var : FindInString("murmuf", syst) ? Var_murmuf: Variations){
          if (var=="") continue;
          mytag = "ZprimeCandidate_"+syst+"_"+var+"_"+tag; book_HFolder(mytag, new HiggsToWWHistsSlim(ctx,mytag));
        }
      }
    }
  }
}

void SignalRegionModule::fill_histograms(uhh2::Event& event, string tag) {
  std::vector<string> mytags = {"event_", "nTopJet_", "ZprimeCandidate_","diLepton_"};
  for (auto& mytag : mytags) HFolder(mytag+ tag)->fill(event);
  string mytag;
  float save_weight = event.weight;
  if (FindInString("_CR", tag)) return;
  if (!FindInString("ZHccvsQCD", tag)) return;

  if (MB["is_ZH"]) {for (std::pair<std::string, std::string> pdf : PFDs) PDF_weights[pdf.first] = m_pdfweights[pdf.first]->GetWeightList(event);}

  for (std::string& syst :  Systematics){
    if (FindInString("PDF", syst)) {
      for(int i=0; i<PDF_variations; i++){
        if (MB["is_ZH"]) {
          double weightFactor = (syst=="PDF")? event.genInfo->systweights().at(i+pdfindex_shift)/event.genInfo->systweights().at(pdfindex_shift) : PDF_weights[syst][i];
          event.weight = save_weight*weightFactor;
        } else {event.weight = save_weight;}
        std::string var = to_string(i);
        mytag = "ZprimeCandidate_"+syst+"_"+var+"_"+tag; HFolder(mytag)->fill(event);
      }
    } else {
      for (std::string& var : FindInString("murmuf", syst) ? Var_murmuf: Variations){
        if (var=="") continue;
        if (!event.isRealData) {
          if (FindInString("murmuf", syst)) event.weight = save_weight*event.get(WFolder(GetSystName(syst, var)+"_in"));
          else event.weight = save_weight*event.get(WFolder(GetSystName(syst, var)+"_in"))/event.get(WFolder(GetSystName(syst, "")+"_in"));
        }
        mytag = "ZprimeCandidate_"+syst+"_"+var+"_"+tag; HFolder(mytag)->fill(event);
      }
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
  MS["year"]              = ctx.get("year");
  MS["dataset_version"]   = ctx.get("dataset_version");
  MB["is_mc"]             = ctx.get("dataset_type") == "MC";
  MB["isPuppi"]           = string2bool(ctx.get("isPuppi"));
  MB["isCHS"]             = string2bool(ctx.get("isCHS"));
  MB["isHOTVR"]           = string2bool(ctx.get("isHOTVR"));
  MB["muonchannel"]       = string2bool(ctx.get("muonchannel"));
  MB["electronchannel"]   = string2bool(ctx.get("electronchannel"));
  MB["invisiblechannel"]  = string2bool(ctx.get("invisiblechannel"));

  MB["is_ZH"]   = MS["dataset_version"].find("MC_ZprimeToZH")!=std::string::npos;
  MB["is_bb"]   = MS["dataset_version"].find("Tobb")!=std::string::npos;
  MB["is_WW"]   = MS["dataset_version"].find("ToWW")!=std::string::npos;
  MB["is_else"] = MS["dataset_version"].find("_extra")!=std::string::npos;

  if ((MB["isPuppi"] && MB["isCHS"]) || (MB["isPuppi"] && MB["isHOTVR"]) || (MB["isCHS"] && MB["isHOTVR"]) ) throw std::runtime_error("In "+NameModule+".cxx: Choose exactly one jet collection.");
  if ((MB["muonchannel"] && MB["electronchannel"]) || (MB["muonchannel"] && MB["invisiblechannel"]) || (MB["electronchannel"] && MB["invisiblechannel"])) throw std::runtime_error("In "+NameModule+".cxx: Choose exactly one lepton channel.");

  MS["topjetLabel"] = MB["isCHS"]? "topjets": (MB["isPuppi"]? "toppuppijets": (MB["isHOTVR"]? "hotvrPuppi": ""));

  // if(!MB["muonchannel"]) Systematics.erase(Systematics.begin()+FindInVector(Systematics,"tracking"));
  if(!MB["muonchannel"]) Systematics.erase(Systematics.begin()+FindInVector(Systematics,"isolation"));
  if(MB["invisiblechannel"]) Systematics.erase(Systematics.begin()+FindInVector(Systematics,"id"));
  if(MB["invisiblechannel"]) Systematics.erase(Systematics.begin()+FindInVector(Systematics,"trigger"));
  if(MB["invisiblechannel"]) Systematics.erase(Systematics.begin()+FindInVector(Systematics,"reco"));

  h_topjets = ctx.get_handle<std::vector<TopJet>>(MS["topjetLabel"]);
  h_ZprimeCandidates = ctx.declare_event_input<std::vector<ZprimeCandidate>>("ZprimeCandidate");

  // Set up histograms:
  book_histograms(ctx);
  book_handles(ctx);
  PrintInputs();

  // Set up selections
  for (std::pair<std::string, std::string> pdf : PFDs) m_pdfweights[pdf.first].reset(new PDFWeights(pdf.second));


  // TODO
  DeepAk8_H4qvsQCD_massdep_SR_selection.reset(new TaggerCut(0, 1,  MassDependentCut_value, "btag_DeepBoosted_H4qvsQCD", h_ZprimeCandidates));

  std::string groom = MB["isHOTVR"]? "_groomed": "";
  tau42_SR_selection.reset(new TaggerCut(0.0, 0.5, -1, "tau42"+groom, h_ZprimeCandidates));
  tau21_SR_selection.reset(new TaggerCut(0.0, 0.4, -1, "tau21"+groom, h_ZprimeCandidates));

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

  if (event.get(h_ZprimeCandidates).size()!=1) return false;
  // if (event.get(h_ZprimeCandidates).size()!=1 && event.get(h_ZprimeCandidates)[0].discriminator("SDmass")<50) return false; // TODO!!!

  ZprimeCandidate cand = event.get(h_ZprimeCandidates)[0];
  if (cand.Zprime_mass()<800) return false;
  if(!MB["invisiblechannel"]){ if(deltaR(cand.leptons()[0], cand.leptons()[1])>0.45) return false;}
  fill_histograms(event, "ExtraCleaning");

  bool ZHccvsQCD_MD_pass = cand.discriminator("btag_DeepBoosted_ZHccvsQCD_MD")>TaggerThr;
  bool HccvsQCD_pass     = cand.discriminator("btag_DeepBoosted_HccvsQCD")>TaggerThr;
  bool H4qvsQCD_pass     = cand.discriminator("btag_DeepBoosted_H4qvsQCD")>TaggerThr;

  bool PN_ZHccvsQCD_MD_pass = cand.discriminator("btag_ParticleNet_ZHccvsQCD")>TaggerThr;
  bool PN_HccvsQCD_pass     = cand.discriminator("btag_ParticleNet_HccvsQCD")>TaggerThr;
  bool PN_H4qvsQCD_pass     = cand.discriminator("btag_ParticleNet_H4qvsQCD")>TaggerThr;

  bool H4qvsQCD_massdep_pass  = DeepAk8_H4qvsQCD_massdep_SR_selection->passes(event);

  if(ZHccvsQCD_MD_pass) fill_histograms(event, "DeepAk8_ZHccvsQCD_MD_SR");
  else fill_histograms(event, "DeepAk8_ZHccvsQCD_MD_CR");

  if(ZHccvsQCD_MD_pass && cand.H().softdropmass()>30) fill_histograms(event, "DeepAk8_ZHccvsQCD_MD2_SR");
  else fill_histograms(event, "DeepAk8_ZHccvsQCD_MD2_CR");

  if(PN_ZHccvsQCD_MD_pass) fill_histograms(event, "PN_ZHccvsQCD_MD_SR");
  else fill_histograms(event, "PN_ZHccvsQCD_MD_CR");

  if(PN_ZHccvsQCD_MD_pass && cand.H().softdropmass()>30) fill_histograms(event, "PN_ZHccvsQCD_MD2_SR");
  else fill_histograms(event, "PN_ZHccvsQCD_MD2_CR");

  if(HccvsQCD_pass) fill_histograms(event, "DeepAk8_HccvsQCD_SR");
  else fill_histograms(event, "DeepAk8_HccvsQCD_CR");

  if(PN_HccvsQCD_pass) fill_histograms(event, "PN_HccvsQCD_SR");
  else fill_histograms(event, "PN_HccvsQCD_CR");

  if(H4qvsQCD_pass) fill_histograms(event, "DeepAk8_H4qvsQCD_SR");
  else fill_histograms(event, "DeepAk8_H4qvsQCD_CR");

  if(PN_H4qvsQCD_pass) fill_histograms(event, "PN_H4qvsQCD_SR");
  else fill_histograms(event, "PN_H4qvsQCD_CR");

  if(H4qvsQCD_massdep_pass) {
    fill_histograms(event, "DeepAk8_H4qvsQCD_massdep_SR");
    if (HccvsQCD_pass)     fill_histograms(event, "DeepAk8_H4qvsQCD_massdep_HccvsQCD_SR");
  } else {
    fill_histograms(event, "DeepAk8_H4qvsQCD_massdep_CR");
    if (!HccvsQCD_pass)     fill_histograms(event, "DeepAk8_H4qvsQCD_massdep_HccvsQCD_CR");
  }

  if(tau42_SR_selection->passes(event)) fill_histograms(event, "tau42_SR");
  else                                  fill_histograms(event, "tau42_CR");
  if(tau21_SR_selection->passes(event)) fill_histograms(event, "tau21_SR");
  else                                  fill_histograms(event, "tau21_CR");

  export_weights(event);

  return true;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the SignalRegionModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(SignalRegionModule)
