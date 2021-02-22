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

class PDFReweightModule: public ModuleBASE {

public:

  explicit PDFReweightModule(uhh2::Context&);
  virtual bool process(uhh2::Event&) override;
  void book_histograms(uhh2::Context&);
  void fill_histograms(uhh2::Event&, string);
  void book_handles(uhh2::Context&);
  void PrintInputs();

protected:

  // Define variables
  std::string NameModule = "PDFReweightModule";
  std::vector<std::string> histogram_tags = { "nocuts", "weights"};
  std::vector<std::string> Systematics = {"PDF", "NNPDF", "NNPDF31_lo_as_0130", "PDF4LHC15_nnlo_100"};

  int PDF_variations = 100; int pdfindex_shift = 9;
  std::unordered_map<std::string, std::string> PFDs = {
    {"NNPDF", "NNPDF31_nnlo_as_0118_nf_4_mc_hessian"},
    {"NNPDF31_lo_as_0130", "NNPDF31_lo_as_0130"},
    {"PDF4LHC15_nnlo_100", "PDF4LHC15_nnlo_100"},
  };

  std::unordered_map<std::string, std::string> MS;
  std::unordered_map<std::string, bool> MB;

  // Define common modules
  std::unique_ptr<AnalysisModule> PDFReweight_module;
  std::unordered_map<std::string, std::vector<double> > PDF_weights;
  std::unordered_map<std::string, std::unique_ptr<PDFWeights> > m_pdfweights;

};


void PDFReweightModule::PrintInputs() {
  std::cout << "****************************************" << std::endl;
  std::cout << "             "+NameModule+"             " << std::endl;
  std::cout << "----------------------------------------" << std::endl;
  for (auto x : MS) std::cout << x.first << std::string( 18-x.first.size(), ' ' ) << x.second << '\n';
  for (auto x : MB) std::cout << x.first << std::string( 18-x.first.size(), ' ' ) << BoolToString(x.second) << '\n';
  std::cout << "****************************************\n" << std::endl;
}

void PDFReweightModule::book_histograms(uhh2::Context& ctx) {
  for(const auto & tag : histogram_tags){
    string mytag;
    mytag = "event_"    + tag; book_HFolder(mytag, new EventHists(ctx,mytag));
    mytag = "gen_"      + tag; book_HFolder(mytag, new GenMatchHists(ctx,mytag));
    mytag = "nTopJet_"  + tag; book_HFolder(mytag, new ExtJetHists(ctx,mytag, MS["topjetLabel"]));
    mytag = "nJet_"     + tag; book_HFolder(mytag, new ExtJetHists(ctx,mytag, MS["jetLabel"]));
    for (std::string& syst :  Systematics){
      for(int i=0; i<PDF_variations; i++){
        std::string var = to_string(i);
        mytag = "ZprimeCandidate_"+syst+"_"+var+"_"+tag; book_HFolder(mytag, new HiggsToWWHistsSlim(ctx,mytag));
      }
    }
  }
}

void PDFReweightModule::fill_histograms(uhh2::Event& event, string tag){
  std::vector<string> mytags = {"event_", "gen_", "nTopJet_", "nJet_"};
  for (auto& mytag : mytags) HFolder(mytag+ tag)->fill(event);
  string mytag;
  float save_weight = event.weight;

  for (std::pair<std::string, std::string> pdf : PFDs) PDF_weights[pdf.first] = m_pdfweights[pdf.first]->GetWeightList(event);

  for (std::string& syst :  Systematics){
    for(int i=0; i<PDF_variations; i++){
      double weightFactor = (syst=="PDF")? event.genInfo->systweights().at(i+pdfindex_shift)/event.genInfo->systweights().at(pdfindex_shift) : PDF_weights[syst][i];
      event.weight = save_weight*weightFactor;
      std::string var = to_string(i);
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

PDFReweightModule::PDFReweightModule(uhh2::Context& ctx){

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

  if ((MB["isPuppi"] && MB["isCHS"]) || (MB["isPuppi"] && MB["isHOTVR"]) || (MB["isCHS"] && MB["isHOTVR"]) ) throw std::runtime_error("In "+NameModule+".cxx: Choose exactly one jet collection.");
  if ((MB["muonchannel"] && MB["electronchannel"]) || (MB["muonchannel"] && MB["invisiblechannel"]) || (MB["electronchannel"] && MB["invisiblechannel"])) throw std::runtime_error("In "+NameModule+".cxx: Choose exactly one lepton channel.");

  if (!FindInString("MC_ZprimeToZH", MS["dataset_version"])) throw std::runtime_error("In "+NameModule+".cxx: Not meant to run on not signal samples.");

  MS["leptons"] = MB["muonchannel"]? "muons": (MB["electronchannel"]? "electrons": (MB["invisiblechannel"]? "invisible": ""));

  MS["jetLabel"]    = MB["isCHS"]? "jets":    (MB["isPuppi"]? "jetsAk4Puppi": (MB["isHOTVR"]? "jetsAk4Puppi": ""));
  MS["topjetLabel"] = MB["isCHS"]? "topjets": (MB["isPuppi"]? "toppuppijets": (MB["isHOTVR"]? "hotvrPuppi": ""));

  // Set up histograms:
  book_histograms(ctx);
  PrintInputs();

  // Set up selections
  PDFReweight_module.reset(new PDFReweight(ctx));
  for (std::pair<std::string, std::string> pdf : PFDs) m_pdfweights[pdf.first].reset(new PDFWeights(pdf.second));
}


/*
█ ██████  ██████   ██████   ██████ ███████ ███████ ███████
█ ██   ██ ██   ██ ██    ██ ██      ██      ██      ██
█ ██████  ██████  ██    ██ ██      █████   ███████ ███████
█ ██      ██   ██ ██    ██ ██      ██           ██      ██
█ ██      ██   ██  ██████   ██████ ███████ ███████ ███████
*/

bool PDFReweightModule::process(uhh2::Event& event) {

  if ((event.year).find(MS["year"])==std::string::npos) throw std::runtime_error("In "+NameModule+".cxx: You are running on "+event.year+" sample with a "+MS["year"]+" year config file. Fix this.");

  fill_histograms(event, "nocuts");

  if (MB["invisiblechannel"] && FindInString("MC_ZprimeToZH_inv", MS["dataset_version"])) PDFReweight_module->process(event);

  fill_histograms(event, "weights");

  return true;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the PDFReweightModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(PDFReweightModule)
