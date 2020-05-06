#include "UHH2/core/include/Event.h"

#include "UHH2/common/include/Utils.h"

#include "UHH2/VHResonances/include/HistsBase.hpp"
#include "UHH2/VHResonances/include/VariableRStudiesHists.h"
#include "UHH2/VHResonances/include/GenLevelJetMatch.h"
#include "UHH2/VHResonances/include/Utils.hpp"

#include "TH1F.h"
#include <type_traits>
#include <iostream>

using namespace std;
using namespace uhh2;


/*
█  ██████  ██████  ███    ██ ███████ ████████ ██████  ██    ██  ██████ ████████  ██████  ██████
█ ██      ██    ██ ████   ██ ██         ██    ██   ██ ██    ██ ██         ██    ██    ██ ██   ██
█ ██      ██    ██ ██ ██  ██ ███████    ██    ██████  ██    ██ ██         ██    ██    ██ ██████
█ ██      ██    ██ ██  ██ ██      ██    ██    ██   ██ ██    ██ ██         ██    ██    ██ ██   ██
█  ██████  ██████  ██   ████ ███████    ██    ██   ██  ██████   ██████    ██     ██████  ██   ██
*/


VariableRStudiesHists::VariableRStudiesHists(uhh2::Context& ctx, const std::string& dname, const std::string& collection_): HistsBase(ctx, dname), collection(collection_) {

  is_mc = ctx.get("dataset_type") == "MC";
  h_topjets = ctx.get_handle<std::vector<TopJet>>(collection);

  vector<string> ParticleNames = {"Z", "Higgs", "W", "top"};
  for (auto& name : ParticleNames) {
    book_ParticleHist(name, name, 20, 1500);
  }

}

void VariableRStudiesHists::fill(const Event & event){
  if(!is_mc) return;
  assert(event.genparticles);
  string name;
  for( auto & gp : *event.genparticles){
    switch (abs(gp.pdgId())) {
      case ParticleID::Z  : name = "Z"; break;
      case ParticleID::H  : name = "Higgs"; break;
      case ParticleID::W  : name = "W"; break;
      case ParticleID::t  : name = "top"; break;
      default: continue;
    }
    fill_ParticleHist(event, name, gp);
  }
}


/*
█ ██████   ██████   ██████  ██   ██      █████  ███    ██ ██████      ███████ ██ ██      ██
█ ██   ██ ██    ██ ██    ██ ██  ██      ██   ██ ████   ██ ██   ██     ██      ██ ██      ██
█ ██████  ██    ██ ██    ██ █████       ███████ ██ ██  ██ ██   ██     █████   ██ ██      ██
█ ██   ██ ██    ██ ██    ██ ██  ██      ██   ██ ██  ██ ██ ██   ██     ██      ██ ██      ██
█ ██████   ██████   ██████  ██   ██     ██   ██ ██   ████ ██████      ██      ██ ███████ ███████
*/

void VariableRStudiesHists::book_ParticleHist(const string & histSuffix, const string & axisSuffix, double minPt, double maxPt){
  book_TH1F("N_"      +histSuffix, "Number of " +axisSuffix, 21, -0.5, 20.5);
  book_TH1F("pt_"     +histSuffix, "p_{T} "     +axisSuffix, 50,minPt,maxPt);
  book_TH1F("mass_"   +histSuffix, "M "         +axisSuffix, 100, 0, 200);
  book_TH1F("flavor_" +histSuffix, "flavor "    +axisSuffix, 101, -50.5, 50.5);
  for (auto mass : {"", "_OS"}) {
    book_TH2F("DeltaR_max_vspt_"+histSuffix+mass, 3000, 0, 3000, 300, 0., 3.);
    for (auto centre : {"_centre", "_jet"}) {
      for (auto dau : {"1", "2", "3", "4"}) {
        if ((histSuffix=="Z"|| histSuffix=="W") and atoi(dau)>2) continue;
        if (histSuffix=="top" and atoi(dau)>3) continue;
        for (int part = 1; part <= atoi(dau); part++) {
          book_TH2F("DeltaR_vspt_"+histSuffix+"_dau"+dau+"_part"+to_string(part)+centre+mass, 3000, 0, 3000, 300, 0., 3.);
        }
      }
    }
  }
  return;
}

void VariableRStudiesHists::fill_ParticleHist(const Event & event, const string & histSuffix, const GenParticle& gp) {
  if(!event.is_valid(h_topjets) || event.get(h_topjets).size()==0) return;
  auto jet = closestParticle(gp, event.get(h_topjets));

  std::vector<const GenParticle*> daughters;
  std::vector<const GenParticle*> genPart;

  daughters.push_back(gp.daughter(event.genparticles, 1));
  daughters.push_back(gp.daughter(event.genparticles, 2));

  for (auto & dau : daughters) {
    if (!dau) continue;
    if (fabs(dau->pdgId())==24 || fabs(dau->pdgId())==23) {
      const GenParticle* sub1 = dau->daughter(event.genparticles, 1);
      const GenParticle* sub2 = dau->daughter(event.genparticles, 2);
      if (sub1 && fabs(sub1->pdgId())<6) genPart.push_back(sub1);
      if (sub2 && fabs(sub2->pdgId())<6) genPart.push_back(sub2);
    } else if ( fabs( dau->pdgId())<6 )  genPart.push_back(dau);
  }


  double deltaR_max=-1;
  double x_cm = 0;
  double y_cm = 0;
  double norm = 0;

  for (unsigned int i = 0; i < genPart.size(); i++) {
    const GenParticle* gp1 = genPart.at(i);
    x_cm += gp1->eta()*gp1->pt();
    y_cm += gp1->phi()*gp1->pt();
    norm += gp1->pt();
    for (unsigned int j = i+1; j < genPart.size(); j++) {
      const GenParticle* gp2 = genPart.at(j);
      deltaR_max = std::max(deltaR(*gp1,*gp2),deltaR_max);
    }
  }

  std::vector<double> DeltaRs, DeltaRs_jet;
  Particle* baricentrum = new Particle();
  baricentrum->set_eta(x_cm/norm);
  baricentrum->set_phi(y_cm/norm);

  // for (unsigned int i = 0; i < genPart.size(); i++) {
  //   const GenParticle* gp1 = genPart.at(i);
  //   DeltaRs.push_back(deltaR(*baricentrum,*gp1));
  //   DeltaRs_jet.push_back(deltaR(*jet,*gp1));
  // }

  for (auto& gp: genPart) {
    DeltaRs.push_back(deltaR(*baricentrum,*gp));
    DeltaRs_jet.push_back(deltaR(*jet,*gp));
  }

  std::sort(DeltaRs.begin(),     DeltaRs.end(),     [](double a, double b) {return a < b;});
  std::sort(DeltaRs_jet.begin(), DeltaRs_jet.end(), [](double a, double b) {return a < b;});


  auto weight = event.weight;
  fill_H1("N_"      +histSuffix, 1,           weight);
  fill_H1("pt_"     +histSuffix, gp.pt(),     weight);
  fill_H1("mass_"   +histSuffix, gp.v4().M(), weight);
  fill_H1("flavor_" +histSuffix, gp.pdgId(),  weight);

  string mass = (gp.v4().M()>70) ? "" : "_OS";
  fill_H2("DeltaR_max_vspt_"+histSuffix+mass, gp.pt(), deltaR_max, weight);

  for (size_t part = 0; part < DeltaRs.size(); part++) {
    TString name = "DeltaR_vspt_"+histSuffix+"_dau"; name += DeltaRs.size(); name +="_part"; name += part+1; name += "_centre"+mass;
    fill_H2((string)name, gp.pt(), DeltaRs[part], weight);
    name = "DeltaR_vspt_"+histSuffix+"_dau"; name += DeltaRs.size(); name +="_part"; name += part+1; name += "_jet"+mass;
    fill_H2((string)name, gp.pt(), DeltaRs_jet[part], weight);
  }

  return;
}
