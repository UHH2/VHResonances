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

class HEMIssueStudyModule: public ModuleBASE {

public:

  explicit HEMIssueStudyModule(uhh2::Context&);
  virtual bool process(uhh2::Event&) override;
  void book_histograms(uhh2::Context&);
  void fill_histograms(uhh2::Event&, string);
  void book_handles(uhh2::Context&);
  void PrintInputs();

protected:

  // Define variables
  std::string NameModule = "HEMIssueStudyModule";
  std::vector<std::string> histogram_tags = { "nocuts", "weights", "HEM", "cleaned", "JetDiLeptonPhiAngular", "QCDRejection", "ScaleFactors", "ExtraCleaning", "DeepAk8_ZHccvsQCD_MD_SR"};

  std::unordered_map<std::string, std::string> MS;
  std::unordered_map<std::string, bool> MB;

  Event::Handle<std::vector<Jet> > h_jets;
  Event::Handle<std::vector<TopJet> > h_topjets;
  Event::Handle<std::vector<ZprimeCandidate> > h_ZprimeCandidates;
  // Define common modules
  std::unique_ptr<uhh2::Selection> lumi_selection;
  std::vector<std::unique_ptr<AnalysisModule>> weightsmodules, modules;
  std::unique_ptr<uhh2::AndSelection> metfilters_selection;
  std::unique_ptr<GenericJetCleaner> GJC, GTJC;

  // Define selections
  std::unordered_map<std::string, std::unique_ptr<Selection>> Trigger_selection;
  std::shared_ptr<Selection> NBoostedJetSel;
  std::shared_ptr<VetoSelection> VetoLeptonSel;
  std::shared_ptr<Selection> NoLeptonSel, NLeptonSel, DeltaRDiLepton_selection, JetDiLeptonPhiAngularSel;
  std::unique_ptr<Selection> HEMEventCleaner_Selection;

  std::unique_ptr<Selection> PTMassCut_selection, DeltaPhiJetMETCut_TopJets_selection, DeltaPhiJetMETCut_Jets_selection;
  std::unique_ptr<AnalysisModule> ZprimeCandidateReconstruction_module;
  std::unique_ptr<AnalysisModule> CollectionProducer_module;
  std::unique_ptr<AnalysisModule> MCScaleVariation_module;
  std::unordered_map<std::string, std::unique_ptr<AnalysisModule>> ScaleFactors_module;
  std::unique_ptr<AnalysisModule> MuonScaleVariations_module;


};

void HEMIssueStudyModule::book_handles(uhh2::Context& ctx) {
  string tag;
  tag = "weight_lumi";  book_WFolder(tag, new Event::Handle< float >, ctx.declare_event_output< float >(tag));
  tag = "weight_GLP";   book_WFolder(tag, new Event::Handle< float >, ctx.declare_event_output< float >(tag));
  tag = "weight_pu";    book_WFolder(tag, new Event::Handle< float >, MB["is_mc"]? ctx.get_handle< float >(tag) : ctx.declare_event_output< float >(tag));
}

void HEMIssueStudyModule::PrintInputs() {
  std::cout << "****************************************" << std::endl;
  std::cout << "             "+NameModule+"             " << std::endl;
  std::cout << "----------------------------------------" << std::endl;
  for (auto x : MS) std::cout << x.first << std::string( 18-x.first.size(), ' ' ) << x.second << '\n';
  for (auto x : MB) std::cout << x.first << std::string( 18-x.first.size(), ' ' ) << BoolToString(x.second) << '\n';
  std::cout << "****************************************\n" << std::endl;
}

void HEMIssueStudyModule::book_histograms(uhh2::Context& ctx) {
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
    mytag = "ZprimeCandidate_" + tag; book_HFolder(mytag, new HiggsToWWHists(ctx,mytag));
  }
}

void HEMIssueStudyModule::fill_histograms(uhh2::Event& event, string tag){
  std::vector<string> mytags = {"event_", "gen_", "nTopJet_", "nJet_", "ele_", "muon_", "diLepton_", "BTagEff_", "Lumi_", "ZprimeCandidate_"};
  for (auto& mytag : mytags) HFolder(mytag+ tag)->fill(event);
}

/*
█  ██████  ██████  ███    ██ ███████ ████████ ██████  ██    ██  ██████ ████████  ██████  ██████
█ ██      ██    ██ ████   ██ ██         ██    ██   ██ ██    ██ ██         ██    ██    ██ ██   ██
█ ██      ██    ██ ██ ██  ██ ███████    ██    ██████  ██    ██ ██         ██    ██    ██ ██████
█ ██      ██    ██ ██  ██ ██      ██    ██    ██   ██ ██    ██ ██         ██    ██    ██ ██   ██
█  ██████  ██████  ██   ████ ███████    ██    ██   ██  ██████   ██████    ██     ██████  ██   ██
*/

HEMIssueStudyModule::HEMIssueStudyModule(uhh2::Context& ctx){

  // Set up variables
  MS["year"]              = ctx.get("year");
  MS["dataset_version"]   = ctx.get("dataset_version");
  MB["is_SingleElectron"] = FindInString("SingleElectron", MS["dataset_version"]);
  MB["is_SinglePhoton"]   = FindInString("SinglePhoton",   MS["dataset_version"]);
  MB["is_mc"]             = ctx.get("dataset_type") == "MC";
  MB["isPuppi"]           = string2bool(ctx.get("isPuppi"));
  MB["isCHS"]             = string2bool(ctx.get("isCHS"));
  MB["isHOTVR"]           = string2bool(ctx.get("isHOTVR"));
  MB["muonchannel"]       = string2bool(ctx.get("muonchannel"));
  MB["electronchannel"]   = string2bool(ctx.get("electronchannel"));
  MB["invisiblechannel"]  = string2bool(ctx.get("invisiblechannel"));
  MS["SysType_PU"]        = ctx.get("SysType_PU");
  MB["lumisel"]           = string2bool(ctx.get("lumisel"));
  MB["mclumiweight"]      = string2bool(ctx.get("mclumiweight"));
  MB["mcpileupreweight"]  = string2bool(ctx.get("mcpileupreweight"));
  MB["eleid"]             = string2bool(ctx.get("eleid"));
  MB["muid"]              = string2bool(ctx.get("muid"));
  MB["tauid"]             = string2bool(ctx.get("tauid"));
  MB["metfilters"]        = string2bool(ctx.get("metfilters"));

  MB["doHEM"]             = string2bool(ctx.get("doHEM"));

  if ((MB["isPuppi"] && MB["isCHS"]) || (MB["isPuppi"] && MB["isHOTVR"]) || (MB["isCHS"] && MB["isHOTVR"]) ) throw std::runtime_error("In "+NameModule+".cxx: Choose exactly one jet collection.");
  if ((MB["muonchannel"] && MB["electronchannel"]) || (MB["muonchannel"] && MB["invisiblechannel"]) || (MB["electronchannel"] && MB["invisiblechannel"])) throw std::runtime_error("In "+NameModule+".cxx: Choose exactly one lepton channel.");

  MS["leptons"] = MB["muonchannel"]? "muons": (MB["electronchannel"]? "electrons": (MB["invisiblechannel"]? "invisible": ""));

  MS["jetLabel"]    = MB["isCHS"]? "jets":    (MB["isPuppi"]? "jetsAk4Puppi": (MB["isHOTVR"]? "jetsAk4Puppi": ""));
  MS["topjetLabel"] = MB["isCHS"]? "topjets": (MB["isPuppi"]? "toppuppijets": (MB["isHOTVR"]? "hotvrPuppi": ""));

  JetPFID::wp JETwp = MB["isCHS"]? JetPFID::WP_TIGHT_CHS: (MB["isPuppi"]? JetPFID::WP_TIGHT_PUPPI : (MB["isHOTVR"]? JetPFID::WP_TIGHT_PUPPI : JetPFID::WP_LOOSE_CHS));

  h_jets = ctx.get_handle<std::vector<Jet>>(MS["jetLabel"]);
  h_topjets = ctx.get_handle<std::vector<TopJet>>(MS["topjetLabel"]);
  h_ZprimeCandidates = ctx.get_handle<std::vector<ZprimeCandidate>>("ZprimeCandidate");

  // Set up histograms:

  book_histograms(ctx);
  book_handles(ctx);
  PrintInputs();

  // Set up selections

  const MuonId muoId = AndId<Muon>(MuonID(Muon::CutBasedIdTrkHighPt), PtEtaCut(min_lepton_pt, min_lepton_eta), MuonID(Muon::PFIsoTight));
  const ElectronId eleId = AndId<Electron>(ElectronTagID(Electron::cutBasedElectronID_Fall17_94X_V2_loose), PtEtaSCCut(min_lepton_pt, min_lepton_eta));

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
  metfilters_selection->add<TriggerSelection>("BadPFMuonFilter", "Flag_BadPFMuonFilter");
  if (MS["year"] != "2016") metfilters_selection->add<EcalBadCalibSelection>("EcalBadCalibSelection"); /*TODO check 2016*/ // Use this instead of Flag_ecalBadCalibFilter, uses ecalBadCalibReducedMINIAODFilter in ntuple_generator
  if (!MB["is_mc"]) metfilters_selection->add<TriggerSelection>("eeBadScFilter", "Flag_eeBadScFilter"); /* TODO Not recommended for MC, but do check */
  /* metfilters_selection->add<TriggerSelection>("BadChargedCandidateFilter", "Flag_BadChargedCandidateFilter"); TODO Not recommended, under review.*/

  //Quick fix for Detector issues
  HEMEventCleaner_Selection.reset(new HEMCleanerSelection(ctx, MS["jetLabel"], MS["topjetLabel"], true, false));

  GJC.reset( new GenericJetCleaner(ctx, MS["jetLabel"],    false, jetId, topjetId, muoId, eleId));
  GTJC.reset(new GenericJetCleaner(ctx, MS["topjetLabel"], true,  jetId, topjetId, muoId, eleId));

  for (auto& t : Trigger_run_validity.at(MS["year"])) {
    if (MB["muonchannel"] && !FindInString("Mu", t.first) ) continue;
    if (MB["muonchannel"] && FindInString("NoMu", t.first) ) continue;
    if (MB["electronchannel"] && !FindInString("Ele", t.first) && !FindInString("Pho", t.first) ) continue;
    if (MB["invisiblechannel"] && !FindInString("MET", t.first) ) continue;
    Trigger_selection[t.first].reset(new TriggerSelection( t.first ));
  }

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
  JetDiLeptonPhiAngularSel.reset(new JetDiLeptonPhiAngularSelection(min_dilep_pt, min_jet_dilep_delta_phi, max_jet_dilep_delta_phi, min_Dphi_AK8jet_MET, MS["leptons"], h_topjets));
  VetoLeptonSel.reset(new VetoSelection(NoLeptonSel));

  if (MB["invisiblechannel"]) {
    // do not run the DeltaRDiLepton_selection on the invisiblechannel
    DeltaRDiLepton_selection.reset(new AndSelection(ctx));
  }

  //Scale factors
  MuonScaleVariations_module.reset(new MuonScaleVariations(ctx));

  MCScaleVariation_module.reset(new MCScaleVariation(ctx));

  ScaleFactors_module["BTag"].reset(new MCBTagScaleFactor(ctx, BTag_algo, BTag_wp, MS["topjetLabel"], "nominal", "lt"));
  ScaleFactors_module["SFs"].reset(new ScaleFactorsManager(ctx, h_ZprimeCandidates));

  ZprimeCandidateReconstruction_module.reset(new ZprimeCandidateReconstruction(ctx, min_dilep_pt, min_DR_dilep, max_DR_dilep, min_jet_dilep_delta_phi, max_jet_dilep_delta_phi, MS["leptons"], MS["topjetLabel"]));
  CollectionProducer_module.reset(new CollectionProducer<ZprimeCandidate>( ctx, "ZprimeCandidate", "ZprimeCandidate", (ZprimeCandidate_ID)ZprimeCandidateID(h_ZprimeCandidates)));

  float min_Z_pt_ZH_mass_cut = min_Z_pt_ZH_mass;
  if (MS["leptons"]=="invisible"){min_Z_pt_ZH_mass_cut = min_Z_pt_ZH_mass_invisible;}
  PTMassCut_selection.reset(new PTMassCut(min_Z_pt_ZH_mass_cut, h_ZprimeCandidates, MS["leptons"]));

  // Delta Phi cut between MET and all TopJets at 2.0
  DeltaPhiJetMETCut_TopJets_selection.reset(new DeltaPhiJetMETCut(ctx, MS["topjetLabel"], min_Dphi_AK8jet_MET, 0, -1));

  // Delta Phi cut between MET and all Jets at 0.5 (QCD rejection)
  DeltaPhiJetMETCut_Jets_selection.reset(new DeltaPhiJetMETCut(ctx, MS["jetLabel"], min_Dphi_AK4jet_MET, 0, -1));

}


/*
█ ██████  ██████   ██████   ██████ ███████ ███████ ███████
█ ██   ██ ██   ██ ██    ██ ██      ██      ██      ██
█ ██████  ██████  ██    ██ ██      █████   ███████ ███████
█ ██      ██   ██ ██    ██ ██      ██           ██      ██
█ ██      ██   ██  ██████   ██████ ███████ ███████ ███████
*/

bool HEMIssueStudyModule::process(uhh2::Event& event) {

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

  bool pass_triggers_OR = false;
  bool pass_Ele_triggers_Photon_Dataset = false;

  for (auto& el : Trigger_selection) {
    if (event.isRealData && (event.run < Trigger_run_validity.at(MS["year"]).at(el.first).first || event.run > Trigger_run_validity.at(MS["year"]).at(el.first).second) ) continue;
    bool pass = el.second->passes(event);
    // For 2016 and 2017 the SinglePhoton and SingleElectron datasets are separete.
    // To avoid double counting, we consider eleTriggers in the SingleElectron
    // and we veto eleTriggers in the SinglePhoton
    if (MS["year"]!="2018") {
      if (MB["is_SingleElectron"] && !FindInString("Ele", el.first)) continue;
      if (MB["is_SinglePhoton"] && FindInString("Ele", el.first) && pass) {pass_Ele_triggers_Photon_Dataset = true; break;}
    }
    pass_triggers_OR += pass;
    pass_triggers_OR += el.second->passes(event);
    if (pass_triggers_OR && !MB["is_SinglePhoton"]) break;
  }
  if (!pass_triggers_OR || pass_Ele_triggers_Photon_Dataset) return false;

  //Effective in 2018 only.
  // Here the assumption is that it should check for cleaned jets to avoid overlap with leptons (Tight WP recommanded).
  if(MB["doHEM"] && !HEMEventCleaner_Selection->passes(event)) return false;
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

  if(!VetoLeptonSel->passes(event)) return false;

  if(!NLeptonSel->passes(event)) return false;

  if(!NBoostedJetSel->passes(event)) return false;

  if (MB["invisiblechannel"]) {
    if(event.met->pt()<min_MET_pt) return false;
  }

  if(!DeltaRDiLepton_selection->passes(event)) return false;

  if(!JetDiLeptonPhiAngularSel->passes(event)) return false;
  fill_histograms(event, "JetDiLeptonPhiAngular");

  if (MB["invisiblechannel"]){
    // QCD rejection, cut at Delta Phi between all jets and MET at min_Dphi_AK4jet_MET (0.5).
    if (!DeltaPhiJetMETCut_Jets_selection->passes(event)) return false;
    // Cut delta Phi between MET and all TopJets at min_Dphi_AK8jet_MET (2.0)
    if (!DeltaPhiJetMETCut_TopJets_selection->passes(event)) return false;
  }
  fill_histograms(event, "QCDRejection");

  MuonScaleVariations_module->process(event);

  ZprimeCandidateReconstruction_module->process(event);

  CollectionProducer_module->process(event);
  if(event.get(h_ZprimeCandidates).size()<1) return false;
  if(event.get(h_ZprimeCandidates).size()>1) {
    for(const auto & cand: event.get(h_ZprimeCandidates)) {
      if (cand.discriminator("SDmass") > cand.H().v4().M()) return false;
      if (cand.discriminator("SDmass") < 60) return false;
    }
  }

  if(!PTMassCut_selection->passes(event)) return false;

  MCScaleVariation_module->process(event);
  for (auto& el : ScaleFactors_module) el.second->process(event);
  fill_histograms(event, "ScaleFactors");

  if (event.get(h_ZprimeCandidates).size()!=1) return false;
  // if (event.get(h_ZprimeCandidates).size()!=1 && event.get(h_ZprimeCandidates)[0].discriminator("SDmass")<50) return false; // TODO!!!

  ZprimeCandidate cand = event.get(h_ZprimeCandidates)[0];
  if (cand.Zprime_mass()<800) return false;
  if(!MB["invisiblechannel"]){ if(deltaR(cand.leptons()[0], cand.leptons()[1])>0.45) return false;}
  fill_histograms(event, "ExtraCleaning");

  bool ZHccvsQCD_MD_pass = cand.discriminator("btag_DeepBoosted_ZHccvsQCD_MD")>TaggerThr;
  if(ZHccvsQCD_MD_pass) fill_histograms(event, "DeepAk8_ZHccvsQCD_MD_SR");

  return true;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the HEMIssueStudyModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(HEMIssueStudyModule)
