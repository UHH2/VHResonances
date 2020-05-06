#include "UHH2/HiggsToWWTagger/include/GenericJetCleaner.h"
#include "UHH2/HiggsToWWTagger/include/HiggsToWWModules.h"
#include <boost/algorithm/string.hpp>

using namespace std;
using namespace uhh2;

void GenericJetCleaner::PrintInputs() {
  std::cout << "****************************************" << std::endl;
  std::cout << "           GenericJetCleaner            " << std::endl;
  std::cout << "----------------------------------------" << std::endl;
  for (auto x : strings) std::cout << x.first << std::string( 18-x.first.size(), ' ' ) << x.second << '\n';
  for (auto x : bools)   std::cout << x.first << std::string( 18-x.first.size(), ' ' ) << (x.second? "true" : "false") << '\n';
  std::cout << "****************************************\n" << std::endl;
}

GenericJetCleaner::GenericJetCleaner(Context & ctx, const string& jetLabel_, const bool& isTopJet_, const JetId& jetID_, const TopJetId& TopjetID_, const MuonId& muoID_, const ElectronId& eleID_):
jetID(jetID_), TopjetID(TopjetID_), muoIDcleaning(muoID_), eleIDcleaning(eleID_) {

  strings["jetLabel"] = jetLabel_;
  bools["isTopJet"]   = isTopJet_;

  strings["year"]   = ctx.get("year");
  strings["JEC"]    = ctx.get("JEC_Version");
  strings["JER"]    = ctx.get("JER_Version");
  bools["is_mc"]    = ctx.get("dataset_type") == "MC";
  bools["isPuppi"]  = string2bool(ctx.get("isPuppi"));
  bools["isCHS"]    = string2bool(ctx.get("isCHS"));
  bools["isHOTVR"]  = string2bool(ctx.get("isHOTVR"));
  bools["isJLC"]    = string2bool(ctx.get("jlc", "false"));
  bools["isJEC"]    = string2bool(ctx.get("jec", "true"));
  bools["isTopJEC"] = string2bool(ctx.get("topjec", "true"));
  bools["isJER"]    = string2bool(ctx.get("jersmear", "true"));
  bools["isTopJER"] = string2bool(ctx.get("topjersmear", "true"));
  bools["isMet"]    = string2bool(ctx.get("do_metcorrection", "false"));
  bools["jetid"]    = string2bool(ctx.get("jetid", "true"));
  bools["topjetid"] = string2bool(ctx.get("topjetid", "true"));

  strings["genjetLabel"]   = bools["isTopJet"]? (bools["isHOTVR"]? "hotvrGen" : "gentopjets") : "genjets"; // TODO make it more general for user
  strings["jet_coll"]      = bools["isTopJet"]? "AK8" : "AK4"; strings["jet_coll"] += bools["isCHS"]? "PFchs" : "PFPuppi"; // Assume correction for HOTVR=Puppi
  bools["do_softdropcorr"] = bools["isTopJEC"] && bools["isTopJet"];

  PrintInputs();

  // Must be added manually. Don't remove
  runs[strings["year"]].push_back("MC");

  string jecTag = strings["JEC"].substr(0,strings["JEC"].find("_V"));
  string jecVer = strings["JEC"].substr(strings["JEC"].find("_V")+2,strings["JEC"].size()-strings["JEC"].find("_V")-2);
  string sfFilename  = "JRDatabase/textFiles/"+strings["JER"]+"_MC/"+strings["JER"]+"_MC_SF_"+strings["jet_coll"]+".txt";
  string resFilename = "JRDatabase/textFiles/"+strings["JER"]+"_MC/"+strings["JER"]+"_MC_PtResolution_"+strings["jet_coll"]+".txt";

  // Define containers

  for (const std::string & run : runs[strings["year"]]) {
    JEC_corr[run] = (run=="MC")? JERFiles::JECFilesMC(jecTag, jecVer, strings["jet_coll"]) : JERFiles::JECFilesDATA(jecTag, jecVer, strings["jet_coll"], run);

    if ( bools["isJLC"] && !bools["isTopJet"] && !bools["isHOTVR"]) { // TODO JLC for HOTVR
      JLC[run].reset(new JetLeptonCleaner_by_KEYmatching(ctx, JEC_corr[run], strings["jetLabel"]));
      if ( muoIDcleaning ) JLC[run]->set_muon_id(muoIDcleaning);
      if ( eleIDcleaning ) JLC[run]->set_electron_id(eleIDcleaning);
    }

    if ( bools["isJEC"]    && !bools["isTopJet"]                      ) JetCorr[run].reset(new GenericJetCorrector(ctx, JEC_corr[run],strings["jetLabel"]));
    if ( bools["isTopJEC"] &&  bools["isTopJet"] && !bools["isHOTVR"] ) TopJetCorr[run].reset(new GenericTopJetCorrector(ctx, JEC_corr[run],strings["jetLabel"]));
    if ( bools["isTopJEC"] &&  bools["isTopJet"]                      ) TopJetSubjetCorr[run].reset(new GenericSubJetCorrector(ctx, JEC_corr[run],strings["jetLabel"]));
    if ( bools["isTopJEC"] &&  bools["isTopJet"] &&  bools["isHOTVR"] ) {
      if (run!="MC") TopJetSubjetCorr[run]->set_HOTVR(bools["isHOTVR"]);
      else TopJetSubjetCorr[run]->set_doJER(ctx, sfFilename, resFilename, strings["genjetLabel"], bools["isHOTVR"]);
    }

  }

  if ( bools["isJER"]           && !(bools["isHOTVR"] && bools["isTopJet"]) ) JetResolutionSmearer.reset(new GenericJetResolutionSmearer(ctx, strings["jetLabel"], strings["genjetLabel"], sfFilename, resFilename));

  if ( bools["do_softdropcorr"] && !bools["isHOTVR"]  ) SoftDropMassUpdate.reset(new SDMassCalculator(ctx,strings["jetLabel"]));
  if ( bools["jetid"]           && !bools["isTopJet"] ) Jetcleaner.reset(new JetCleaner(ctx, jetID, strings["jetLabel"]));
  if ( bools["topjetid"]        &&  bools["isTopJet"] ) TopJetcleaner.reset(new TopJetCleaner(ctx, TopjetID, strings["jetLabel"]));


}

bool GenericJetCleaner::process(Event & event) {

  std::unordered_map<std::string, bool > apply_run;

  for (const std::string & run : runs[strings["year"]]) apply_run[run] = false;

  bool apply_all = false;
  if (!bools["is_mc"]) {
    for (const std::string & run : runs[strings["year"]]) {
      if (run=="MC") continue;
      if (run_number_map.at(strings["year"]).at(run).first <= event.run && event.run <= run_number_map.at(strings["year"]).at(run).second) apply_run[run] = true;
      apply_all+=apply_run[run];
    }
  } else {apply_run["MC"] = true; apply_all+=apply_run["MC"];}
  if (apply_all != 1) throw std::runtime_error("In GenericJetCleaner.cxx: Sum of apply_all when applying JECs is not == 1. Fix this.");


  for (const std::string run : runs[strings["year"]]) {
    if (apply_run[run]){
      if ( bools["isJLC"]    && !bools["isTopJet"] && !bools["isHOTVR"] ) JLC[run]->process(event);
      if ( bools["isJEC"]    && !bools["isTopJet"]                      ) JetCorr[run]->process(event);
      // if ( !bools["isTopJet"] && bools["isMet"]                         ) JetCorr[run]->correct_met(event, bools["isCHS"]); // TODO class does't has the function
      if ( bools["isTopJEC"] &&  bools["isTopJet"] && !bools["isHOTVR"] ) TopJetCorr[run]->process(event);
      if ( bools["isTopJEC"] &&  bools["isTopJet"]                      ) TopJetSubjetCorr[run]->process(event);
    }
  }

  if ( bools["isJER"]           && !(bools["isHOTVR"] && bools["isTopJet"]) ) JetResolutionSmearer->process(event);
  if ( bools["do_softdropcorr"] && !bools["isHOTVR"]                        ) SoftDropMassUpdate->process(event);
  if ( bools["jetid"]           && !bools["isTopJet"]                       ) Jetcleaner->process(event);
  if ( bools["topjetid"]        &&  bools["isTopJet"]                       ) TopJetcleaner->process(event);

  return true;
}

GenericJetCleaner::~GenericJetCleaner(){}
