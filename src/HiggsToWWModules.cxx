#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Utils.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/common/include/MuonIds.h"
#include "UHH2/VHResonances/include/constants.hpp"
#include "UHH2/VHResonances/include/HiggsToWWModules.h"
#include <UHH2/VHResonances/include/ZprimeCandidate.h>

#include <stdexcept>
#include <set>

using namespace std;
using namespace uhh2;


FinalStateMatching::FinalStateMatching(Context & ctx) {
  h_ZDecay_ = ctx.declare_event_output<float>("ZDecay");
  h_HDecay_ = ctx.declare_event_output<float>("HDecay");
  h_ZprimeDecay_ = ctx.declare_event_output<float>("ZprimeDecay");
  skipMatching = !((ctx.get("dataset_type") == "MC" ) && (ctx.get("dataset_version").find("Zprime") != std::string::npos));
  GenParticles_printer.reset(new GenParticlesPrinter(ctx));
}

bool FinalStateMatching::process(Event &event){

  if (skipMatching) {
    event.set(h_ZDecay_, nomatch);
    event.set(h_HDecay_, nomatch);
    event.set(h_ZprimeDecay_, nomatch);
    return true;
  }

  std::map<std::string, int> ParticleFlavorMap;
  std::map<std::string, ZprimeDecay> ParticleDecayMap = {{"Z",nomatch}, {"H",nomatch}};

  for(const auto & gp: *event.genparticles){
    if (fabs(gp.pdgId())!=ZPrime) continue;
    const GenParticle* d1 = gp.daughter(event.genparticles, 1);
    const GenParticle* d2 = gp.daughter(event.genparticles, 2);
    ParticleFlavorMap["Z"] = 1;
    ParticleFlavorMap["H"] = 2;
    if (!d1 || !d2) continue;
    if (DoubleDecay(d1->pdgId(), d2->pdgId(),ZH)) {
      if (d2->pdgId()!=25) {
        ParticleFlavorMap["Z"] = 2;
        ParticleFlavorMap["H"] = 1;
      }
    }
    else throw std::runtime_error("In FinalStateSelection: Zprime is decaying into "+std::to_string(d1->pdgId())+" and "+std::to_string(d2->pdgId()));


    const GenParticle* Z_dau1 = gp.daughter(event.genparticles, ParticleFlavorMap["Z"])->daughter(event.genparticles, 1);
    const GenParticle* Z_dau2 = gp.daughter(event.genparticles, ParticleFlavorMap["Z"])->daughter(event.genparticles, 2);
    if (Z_dau1 && Z_dau2) {
      if (DoubleDecay(Z_dau1->pdgId(), Z_dau2->pdgId(),ee))        ParticleDecayMap["Z"] = Zee;
      else if (DoubleDecay(Z_dau1->pdgId(), Z_dau2->pdgId(),mumu)) ParticleDecayMap["Z"] = Zmumu;
      else ParticleDecayMap["Z"] = Zelse;
    }

    const GenParticle* H_dau1 = gp.daughter(event.genparticles, ParticleFlavorMap["H"])->daughter(event.genparticles, 1);
    const GenParticle* H_dau2 = gp.daughter(event.genparticles, ParticleFlavorMap["H"])->daughter(event.genparticles, 2);
    if (H_dau1 && H_dau2) {
      if (DoubleDecay(H_dau1->pdgId(), H_dau2->pdgId(),bb)) ParticleDecayMap["H"] = Hbb;
      else if (DoubleDecay(H_dau1->pdgId(), H_dau2->pdgId(),WW)) {
        const GenParticle* W_dau1_1 = H_dau1->daughter(event.genparticles, 1);
        const GenParticle* W_dau1_2 = H_dau1->daughter(event.genparticles, 2);
        const GenParticle* W_dau2_1 = H_dau2->daughter(event.genparticles, 1);
        const GenParticle* W_dau2_2 = H_dau2->daughter(event.genparticles, 2);
        auto ID1_1 = W_dau1_1->pdgId();
        auto ID1_2 = W_dau1_2->pdgId();
        auto ID2_1 = W_dau2_1->pdgId();
        auto ID2_2 = W_dau2_2->pdgId();
        if (W_dau1_1 && W_dau1_2 && W_dau2_1 && W_dau2_2 && DoubleDecay(ID1_1, ID1_2,hadronic) && DoubleDecay(ID2_1, ID2_2,hadronic) ){
          if (DoubleDecay(ID1_1, ID1_2,hadronic) && DoubleDecay(ID2_1, ID2_2,hadronic)) {
            ParticleDecayMap["H"] = HWW;
          } else ParticleDecayMap["H"] = Helse;
        } else ParticleDecayMap["H"] = Helse;
      }
      else ParticleDecayMap["H"] = Helse;
    }
    // else {
    //   std::cout << H_dau1 << " " << H_dau2 << " " << gp.daughter(event.genparticles, ParticleFlavorMap["H"])->status() << " " << gp.daughter(event.genparticles, ParticleFlavorMap["H"])->index() << " " << gp.daughter(event.genparticles, ParticleFlavorMap["H"])->mother1() << " " << gp.daughter(event.genparticles, ParticleFlavorMap["H"])->mother2() << " " << gp.daughter(event.genparticles, ParticleFlavorMap["H"])->daughter1() << " " << gp.daughter(event.genparticles, ParticleFlavorMap["H"])->daughter2() << " " << gp.daughter(event.genparticles, ParticleFlavorMap["H"])->spin() << '\n';
    //   GenParticles_printer->process(event);
    //
    // }

    event.set(h_ZDecay_, (float)ParticleDecayMap["Z"]);
    event.set(h_HDecay_, (float)ParticleDecayMap["H"]);
    event.set(h_ZprimeDecay_, (float)(ParticleDecayMap["H"]+ParticleDecayMap["Z"]));
  }


  return true;
}





#define MYTAGSETDISCR(mytag)\
candidate.set_discriminators(#mytag,   jet.mytag());\


ZprimeCandidateReconstruction::ZprimeCandidateReconstruction(Context& ctx, float pt_min_, float DR_min_, float DR_max_, float phi_min_, float phi_max_, const string& lepton_, const string& topjetcollection_) : pt_min(pt_min_), DR_min(DR_min_), DR_max(DR_max_), phi_min(phi_min_), phi_max(phi_max_), lepton(lepton_), topjetcollection(topjetcollection_){

  h_topjets = ctx.get_handle<vector<TopJet>>(topjetcollection_);
  h_ZprimeCandidates_ = ctx.declare_event_output<vector<ZprimeCandidate>>("ZprimeCandidate");

}


bool ZprimeCandidateReconstruction::process(Event& event){
  assert(event.muons || event.electrons);

  const auto & jets = event.get(h_topjets);

  vector<Particle> leptons;
  if (lepton == "muons") leptons.assign(event.muons->begin(), event.muons->end());
  else if (lepton == "electrons") leptons.assign(event.electrons->begin(), event.electrons->end());
  else throw logic_error("ZprimeCandidateReconstruction: Impossible case");

  // Declare output

  std::vector<ZprimeCandidate> candidates;

  if(jets.size() < 1 ) throw std::runtime_error("ZprimeCandidateReconstruction::ZprimeCandidateReconstruction -- unexpected number of jets");
  if(leptons.size() < min_leptons ) throw std::runtime_error("ZprimeCandidateReconstruction::ZprimeCandidateReconstruction -- unexpected number of leptons");

  std::map<TString, JetId> Btag_map;
  Btag_map["DeepCSV_loose"]   = BTag(BTag::DEEPCSV, BTag::WP_LOOSE);
  Btag_map["DeepCSV_medium"]  = BTag(BTag::DEEPCSV, BTag::WP_MEDIUM);
  Btag_map["DeepCSV_tight"]   = BTag(BTag::DEEPCSV, BTag::WP_TIGHT);

  for (unsigned int i = 0; i < leptons.size(); i++) {
    Particle lep1 = leptons.at(i);
    for (unsigned int j = i+1; j < leptons.size(); j++) {
      Particle lep2 = leptons.at(j);
      if (lepton == "muons") {
        bool TrkhighID1, TrkhighID2;
        TrkhighID1 = MuonID(Muon::CutBasedIdGlobalHighPt)(event.muons->at(i), event);
        TrkhighID2 = MuonID(Muon::CutBasedIdGlobalHighPt)(event.muons->at(j), event);
        if (!TrkhighID1 && !TrkhighID2) continue;
      }
      auto DR = uhh2::deltaR(lep1, lep2);
      auto diLep = lep1.v4() + lep2.v4();
      if( (fabs(diLep.M() - ZMASS) < ZWIDTH) && (lep1.charge() + lep2.charge() == 0) && diLep.pt()>pt_min && DR_min < DR && DR < DR_max) {
        for(const auto & jet: jets){
          auto Dphi = deltaPhi(diLep, jet);
          if( phi_min < Dphi  && Dphi< phi_max){
            auto ZplusJet = diLep + jet.v4();
            double chi2 = TMath::Power(((jet.softdropmass()-HMASS)/HWIDTH),2)+ TMath::Power(((diLep.M()-ZMASS)/ZWIDTH),2);
            ZprimeCandidate candidate;
            candidate.set_Zprime(ZplusJet);
            candidate.set_Z(diLep);
            candidate.set_H(jet);
            candidate.set_jets_leptonic({lep1,lep2});
            candidate.set_discriminators("chi2",   chi2);
            candidate.set_discriminators("subjets",jet.subjets().size());
            candidate.set_discriminators("btag_DeepCSV_loose",  (double)MultiBTagSubJetID(Btag_map["DeepCSV_loose"])(jet, event));
            candidate.set_discriminators("btag_DeepCSV_medium", (double)MultiBTagSubJetID(Btag_map["DeepCSV_medium"])(jet, event));
            candidate.set_discriminators("btag_DeepCSV_tight",  (double)MultiBTagSubJetID(Btag_map["DeepCSV_tight"])(jet, event));
            int sj_ind = 0;
            for(auto subjet : jet.subjets()) {
              candidate.set_discriminators("btag_DeepCSV_loose_subjet_"+std::to_string(sj_ind),  (double)Btag_map["DeepCSV_loose"](subjet,event));
              candidate.set_discriminators("btag_DeepCSV_medium_subjet_"+std::to_string(sj_ind), (double)Btag_map["DeepCSV_medium"](subjet,event));
              candidate.set_discriminators("btag_DeepCSV_tight_subjet_"+std::to_string(sj_ind),  (double)Btag_map["DeepCSV_tight"](subjet,event));
              sj_ind++;
            }
            candidate.set_discriminators("SDmass",    jet.softdropmass());
            candidate.set_discriminators("tau1",      jet.tau1());
            candidate.set_discriminators("tau2",      jet.tau2());
            candidate.set_discriminators("tau3",      jet.tau3());
            candidate.set_discriminators("tau4",      jet.tau4());
            candidate.set_discriminators("tau21",     (jet.tau1()!=0) ? (jet.tau2()/jet.tau1()) : -1);
            candidate.set_discriminators("tau31",     (jet.tau1()!=0) ? (jet.tau3()/jet.tau1()) : -1);
            candidate.set_discriminators("tau41",     (jet.tau1()!=0) ? (jet.tau4()/jet.tau1()) : -1);
            candidate.set_discriminators("tau32",     (jet.tau2()!=0) ? (jet.tau3()/jet.tau2()) : -1);
            candidate.set_discriminators("tau42",     (jet.tau2()!=0) ? (jet.tau4()/jet.tau2()) : -1);
            candidate.set_discriminators("tau43",     (jet.tau3()!=0) ? (jet.tau4()/jet.tau3()) : -1);
            candidate.set_discriminators("NN_HWW",    jet.has_tag(TopJet::NN_HWW)   ? jet.get_tag(TopJet::NN_HWW)   : 9999);
            candidate.set_discriminators("NN_Hbb",    jet.has_tag(TopJet::NN_Hbb)   ? jet.get_tag(TopJet::NN_Hbb)   : 9999);
            candidate.set_discriminators("NN_QCD",    jet.has_tag(TopJet::NN_QCD)   ? jet.get_tag(TopJet::NN_QCD)   : 9999);
            candidate.set_discriminators("NN_Top",    jet.has_tag(TopJet::NN_Top)   ? jet.get_tag(TopJet::NN_Top)   : 9999);
            candidate.set_discriminators("NN_W",      jet.has_tag(TopJet::NN_W)     ? jet.get_tag(TopJet::NN_W)     : 9999);
            candidate.set_discriminators("NN_Z",      jet.has_tag(TopJet::NN_Z)     ? jet.get_tag(TopJet::NN_Z)     : 9999);
            candidate.set_discriminators("NN_HWW_1",  jet.has_tag(TopJet::NN_HWW_1) ? jet.get_tag(TopJet::NN_HWW_1) : 9999);
            candidate.set_discriminators("NN_Hbb_1",  jet.has_tag(TopJet::NN_Hbb_1) ? jet.get_tag(TopJet::NN_Hbb_1) : 9999);
            candidate.set_discriminators("NN_QCD_1",  jet.has_tag(TopJet::NN_QCD_1) ? jet.get_tag(TopJet::NN_QCD_1) : 9999);
            candidate.set_discriminators("NN_Top_1",  jet.has_tag(TopJet::NN_Top_1) ? jet.get_tag(TopJet::NN_Top_1) : 9999);
            candidate.set_discriminators("NN_W_1",    jet.has_tag(TopJet::NN_W_1)   ? jet.get_tag(TopJet::NN_W_1)   : 9999);
            candidate.set_discriminators("NN_Z_1",    jet.has_tag(TopJet::NN_Z_1)   ? jet.get_tag(TopJet::NN_Z_1)   : 9999);
            candidate.set_discriminators("NN_HWW_2",  jet.has_tag(TopJet::NN_HWW_2) ? jet.get_tag(TopJet::NN_HWW_2) : 9999);
            candidate.set_discriminators("NN_Hbb_2",  jet.has_tag(TopJet::NN_Hbb_2) ? jet.get_tag(TopJet::NN_Hbb_2) : 9999);
            candidate.set_discriminators("NN_QCD_2",  jet.has_tag(TopJet::NN_QCD_2) ? jet.get_tag(TopJet::NN_QCD_2) : 9999);
            candidate.set_discriminators("NN_Top_2",  jet.has_tag(TopJet::NN_Top_2) ? jet.get_tag(TopJet::NN_Top_2) : 9999);
            candidate.set_discriminators("NN_W_2",    jet.has_tag(TopJet::NN_W_2)   ? jet.get_tag(TopJet::NN_W_2)   : 9999);
            candidate.set_discriminators("NN_Z_2",    jet.has_tag(TopJet::NN_Z_2)   ? jet.get_tag(TopJet::NN_Z_2)   : 9999);
            candidate.set_discriminators("CNN_HWW",   jet.has_tag(TopJet::CNN_HWW)  ? jet.get_tag(TopJet::CNN_HWW)  : 9999);
            candidate.set_discriminators("CNN_Hbb",   jet.has_tag(TopJet::CNN_Hbb)  ? jet.get_tag(TopJet::CNN_Hbb)  : 9999);
            candidate.set_discriminators("CNN_QCD",   jet.has_tag(TopJet::CNN_QCD)  ? jet.get_tag(TopJet::CNN_QCD)  : 9999);
            candidate.set_discriminators("CNN_Top",   jet.has_tag(TopJet::CNN_Top)  ? jet.get_tag(TopJet::CNN_Top)  : 9999);
            candidate.set_discriminators("CNN_W",     jet.has_tag(TopJet::CNN_W)    ? jet.get_tag(TopJet::CNN_W)    : 9999);
            candidate.set_discriminators("CNN_Z",     jet.has_tag(TopJet::CNN_Z)    ? jet.get_tag(TopJet::CNN_Z)    : 9999);
            candidate.set_discriminators("DCL_HWW",   jet.has_tag(TopJet::DCL_HWW)  ? jet.get_tag(TopJet::DCL_HWW)  : 9999);
            candidate.set_discriminators("DCL_Hbb",   jet.has_tag(TopJet::DCL_Hbb)  ? jet.get_tag(TopJet::DCL_Hbb)  : 9999);
            candidate.set_discriminators("DCL_QCD",   jet.has_tag(TopJet::DCL_QCD)  ? jet.get_tag(TopJet::DCL_QCD)  : 9999);
            candidate.set_discriminators("DCL_Top",   jet.has_tag(TopJet::DCL_Top)  ? jet.get_tag(TopJet::DCL_Top)  : 9999);
            candidate.set_discriminators("DCL_W",     jet.has_tag(TopJet::DCL_W)    ? jet.get_tag(TopJet::DCL_W)    : 9999);
            candidate.set_discriminators("DCL_Z",     jet.has_tag(TopJet::DCL_Z)    ? jet.get_tag(TopJet::DCL_Z)    : 9999);

            candidate.set_discriminators("Match",           jet.has_tag(TopJet::Matching)       ? jet.get_tag(TopJet::Matching)       : 0);
            candidate.set_discriminators("MatchingStatus",  jet.has_tag(TopJet::MatchingStatus) ? jet.get_tag(TopJet::MatchingStatus) : 0);
            MYTAGSETDISCR(btag_DeepBoosted_H4qvsQCD)
            MYTAGSETDISCR(btag_MassDecorrelatedDeepBoosted_H4qvsQCD)
            MYTAGSETDISCR(btag_DeepBoosted_probHqqqq)
            MYTAGSETDISCR(btag_DeepBoosted_raw_score_h)
            MYTAGSETDISCR(btag_MassDecorrelatedDeepBoosted_probHqqqq)
            candidates.emplace_back(candidate);
          }
        }
      }
    }
  }

  // Set the handle with all candidates
  event.set(h_ZprimeCandidates_, candidates);

  return true;
}



SDMassCalculator::SDMassCalculator(uhh2::Context & ctx, const std::string & jetCollName) {
  h_topjets_ = ctx.get_handle<std::vector<TopJet>>(jetCollName);
}

bool SDMassCalculator::process(uhh2::Event & event) {
  std::vector<TopJet>* topjets(0);

  if (event.is_valid(h_topjets_)) topjets = &event.get(h_topjets_);
  else throw std::runtime_error("SDMassCalculator::process -- invalid handle to topjets");

  for (auto & jet : *topjets) jet.set_softdropmass(calcSDmass(jet));
  return true;
}

float SDMassCalculator::calcSDmass(const TopJet & jet) {
  // Calculate uncorrected SD mass from subjets
  LorentzVector puppi_softdrop;
  for (auto & subjet : jet.subjets()) puppi_softdrop += subjet.v4();
  return inv_mass_safe(puppi_softdrop);
}



BlindData::BlindData(Context& ctx, const Event::Handle<vector<ZprimeCandidate> > & h_ZprimeCandidates): h_ZprimeCandidates_(h_ZprimeCandidates){
  h_is_Blind = ctx.get_handle<bool>("is_Blind");
}
bool BlindData::process(Event & event){

  bool is_Blind = false;
  const auto & ZprimeCandidates = event.get(h_ZprimeCandidates_);
  for(auto & cand: ZprimeCandidates){
    float disc = cand.discriminator("is_Blind");
    if (disc==0) is_Blind=true;
  }
  event.set(h_is_Blind, is_Blind);
  return true;
}



/*
#  ######   ######     ###    ##       ########    ########    ###     ######  ########  #######  ########   ######
# ##    ## ##    ##   ## ##   ##       ##          ##         ## ##   ##    ##    ##    ##     ## ##     ## ##    ##
# ##       ##        ##   ##  ##       ##          ##        ##   ##  ##          ##    ##     ## ##     ## ##
#  ######  ##       ##     ## ##       ######      ######   ##     ## ##          ##    ##     ## ########   ######
#       ## ##       ######### ##       ##          ##       ######### ##          ##    ##     ## ##   ##         ##
# ##    ## ##    ## ##     ## ##       ##          ##       ##     ## ##    ##    ##    ##     ## ##    ##  ##    ##
#  ######   ######  ##     ## ######## ########    ##       ##     ##  ######     ##     #######  ##     ##  ######
*/

// Generic Class for Applying SFs
void ScaleFactorsFromHistos::LoadHisto(TFile* file, std::string name, std::string hname) {
  histos[name].reset((TH1F*)file->Get(hname.c_str()));
  histos[name]->SetDirectory(0);
};

double ScaleFactorsFromHistos::Evaluator(std::string hname, double var) {
  // invalid cases
  if (var == uhh2::infinity) return 1.0;

  int firstBin = 1;
  int lastBin  = histos[hname]->GetNbinsX();
  double h_min = histos[hname]->GetBinCenter(firstBin)-histos[hname]->GetBinError(firstBin);
  double h_max = histos[hname]->GetBinCenter(lastBin)+histos[hname]->GetBinError(lastBin);
  double var_for_eval = var;
  var_for_eval = (var_for_eval > h_min) ? var_for_eval : h_min+0.001;
  var_for_eval = (var_for_eval < h_max) ? var_for_eval : h_max-0.001;
  return histos[hname]->GetBinContent(histos[hname]->FindBin(var_for_eval));
};


/*
# ######## ##     ## ########  #######  ########  ##    ##     ######  ########
#    ##    ##     ## ##       ##     ## ##     ##  ##  ##     ##    ## ##
#    ##    ##     ## ##       ##     ## ##     ##   ####      ##       ##
#    ##    ######### ######   ##     ## ########     ##        ######  ######
#    ##    ##     ## ##       ##     ## ##   ##      ##             ## ##
#    ##    ##     ## ##       ##     ## ##    ##     ##       ##    ## ##
#    ##    ##     ## ########  #######  ##     ##    ##        ######  ##
*/



// Apply Theory weights
NLOCorrections::NLOCorrections(uhh2::Context& ctx) {

  std::string dataset_version = ctx.get("dataset_version");
  // Corrections for 2017 and 2018 are the same. 2016 is different
  is2017 = ! FindInString("2016", ctx.get("year"));

  //TODO it's arbitrary.
  is_Wjets  = FindInString("MC_WJets",dataset_version);
  is_DY     = FindInString("MC_DY",dataset_version);
  is_Znn    = FindInString("MC_Znunu",dataset_version);
  is_Zjets  = is_DY || is_Znn;

  std::string folder_ = ctx.get("NLOCorrections"); //TODO better name
  for (const std::string& proc: {"w","z"}) {
    TFile* file_ = new TFile((folder_+"merged_kfactors_"+proc+"jets.root").c_str());
    for (const std::string& corr: {"ewk","qcd","qcd_ewk"}) LoadHisto(file_, proc+"_"+corr, "kfactor_monojet_"+corr);
    file_->Close();
  }
  for (const std::string& proc: {"dy","znn"}) {
    TFile* file_ = new TFile((folder_+"kfac_"+proc+"_filter.root").c_str());
    LoadHisto(file_, proc+"_qcd_2017", "kfac_"+proc+"_filter");
    file_->Close();
  }
  TFile* file_ = new TFile((folder_+"2017_gen_v_pt_qcd_sf.root").c_str());
  LoadHisto(file_, "w_qcd_2017", "wjet_dress_inclusive");
  file_->Close();
  file_ = new TFile((folder_+"lindert_qcd_nnlo_sf.root").c_str());
  for (const std::string& proc: {"eej", "evj", "vvj"}) LoadHisto(file_, proc+"_qcd_nnlo", proc);
  file_->Close();
}


double NLOCorrections::GetPartonObjectPt(uhh2::Event& event, ParticleID objID) {
  for(const auto & gp : *event.genparticles) {if (gp.pdgId()==objID) return gp.pt(); }
  return uhh2::infinity;
};


bool NLOCorrections::process(uhh2::Event& event){
  // Sample dependant corrections
  if (!is_Wjets && !is_Zjets && event.isRealData) return true;

  double objpt = uhh2::infinity, theory_weight = 1.0;
  std::string process = "";

  // TODO only DY or also Znn?
  if (is_DY) theory_weight /= default_DY_kFactor;

  if (is_Zjets) objpt = GetPartonObjectPt(event,ParticleID::Z);
  if (is_Wjets) objpt = GetPartonObjectPt(event,ParticleID::W);

  if (is_Zjets) process = "z";
  if (is_Wjets) process = "w";

  if (do_QCD_EWK) theory_weight *= Evaluator(process+"_qcd_ewk",objpt);
  else {
    if (do_EWK) theory_weight *= Evaluator(process+"_ewk",objpt);
    if (do_QCD_NLO) {
      if (is2017) {
        if (is_DY)  process = "dy";
        if (is_Znn) process = "znn";
      }
      theory_weight *= Evaluator(process+"_qcd"+(is2017?"_2017":""),objpt);
    }
  }

  if (do_QCD_NNLO) {
    if (is_DY)    process = "eej";
    if (is_Znn)   process = "vvj";
    if (is_Wjets) process = "evj";
    theory_weight *= Evaluator(process+"_qcd_nnlo",objpt);
  }

  event.weight *= theory_weight;

  return true;
}
