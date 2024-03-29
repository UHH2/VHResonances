#include "UHH2/core/include/Event.h"

#include "UHH2/common/include/Utils.h"

#include "UHH2/VHResonances/include/HistsBase.hpp"
#include "UHH2/VHResonances/include/GenMatchHists.h"

#include "TH1F.h"
#include "TH2F.h"
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


GenMatchHists::GenMatchHists(Context & ctx, const string & dname): HistsBase(ctx, dname) {

  is_mc = ctx.get("dataset_type") == "MC";

  vector<string> ParticleNames = {"Z", "Higgs", "Wplus", "Wminus", "t", "tbar", "ZPrime"};
  for (auto& name : ParticleNames) {
    book_ParticleHist(name, name, 20, 1500);
    book_FamilyTree(name);
  }

  book_TH1F("MET", "missing E_{T}", 200, 0, 2000);
  book_TH1F("METphi", "missing E_{T} #phi", 200,-5,5);
  book_TH1F("qScale", "qScale", 100, 0, 3000);
  book_TH1F("PU_pT_hat_max", "PU p_T hat max", 100, 0, 3000);
  book_TH1F("binningValues", "binning values", 100, 0, 3000);
  book_TH1F("HT", "H_T", 100, 0, 3000);
  book_TH2F("MET_METphi", ";missing E_{T};missing E_{T} #phi", 200,0, 2000, 200,-5,5);
}

void GenMatchHists::fill(const Event & event){
  if(!is_mc) return;
  assert(event.genparticles);

  fill_H1("MET", event.met->pt(), event.weight);
  fill_H1("METphi", event.met->phi(), event.weight);
  fill_H1("qScale", event.genInfo->qScale(), event.weight);
  fill_H1("PU_pT_hat_max", event.genInfo->PU_pT_hat_max(), event.weight);

  if (event.genInfo->binningValues().size() != 0){
    fill_H1("binningValues", event.genInfo->binningValues()[0], event.weight);
  }

  double ht=0;
  constexpr const int invalid_daughter = (unsigned short)(-1);

  for(const auto & gp : *event.genparticles){
    if(gp.daughter1() != invalid_daughter || gp.daughter2() != invalid_daughter) continue;
    // if we are here, it means we have a final state particle.
    // Add to HT in case it is a parton (quark -- including b but not top -- or gluon).
    // Note that the exact HT definition depends on the madgraph configuration, but this
    // should cover the most common case.
    int id = abs(gp.pdgId());
    if((id >= 1 && id <= 5) || (id == 21)){
      ht += gp.pt();
    }
  }

  fill_H1("HT", ht, event.weight);
  H2("MET_METphi")->Fill(event.met->pt(), event.met->phi(), event.weight);

  string name;
  for( auto & gp : *event.genparticles){
    switch (gp.pdgId()) {
      case ParticleID::Z      : name = "Z";      break;
      case ParticleID::H      : name = "Higgs";  break;
      case ParticleID::W      : name = "Wplus";  break;
      case -ParticleID::W     : name = "Wminus"; break;
      case ParticleID::t      : name = "t";      break;
      case -ParticleID::t     : name = "tbar";   break;
      case ParticleID::ZPrime : name = "ZPrime"; break;

      default: continue;
    }
    fill_ParticleHist(event, name, gp);
    fill_FamilyTree(event, name, gp);
  }
}


/*
█ ██████   ██████   ██████  ██   ██      █████  ███    ██ ██████      ███████ ██ ██      ██
█ ██   ██ ██    ██ ██    ██ ██  ██      ██   ██ ████   ██ ██   ██     ██      ██ ██      ██
█ ██████  ██    ██ ██    ██ █████       ███████ ██ ██  ██ ██   ██     █████   ██ ██      ██
█ ██   ██ ██    ██ ██    ██ ██  ██      ██   ██ ██  ██ ██ ██   ██     ██      ██ ██      ██
█ ██████   ██████   ██████  ██   ██     ██   ██ ██   ████ ██████      ██      ██ ███████ ███████
*/

void GenMatchHists::book_ParticleHist(const string & histSuffix, const string & axisSuffix, double minPt, double maxPt){
  book_TH1F("N_"      +histSuffix, "Number of " +axisSuffix, 21, -0.5, 20.5);
  book_TH1F("pt_"     +histSuffix, "p_{T} "     +axisSuffix, 50,minPt,maxPt);
  book_TH1F("eta_"    +histSuffix, "#eta "      +axisSuffix, 100,-5,5);
  book_TH1F("phi_"    +histSuffix, "#phi "      +axisSuffix, 50,-M_PI,M_PI);
  book_TH1F("mass_"   +histSuffix, "M "         +axisSuffix, (histSuffix=="ZPrime")? 1000 : 100, 0, (histSuffix=="ZPrime")? 10000 : 200);
  // book_TH1F("mass_"   +histSuffix, "M "         +axisSuffix, 100, 0, 200);
  book_TH1F("flavor_" +histSuffix, "flavor "    +axisSuffix, 101, -50.5, 50.5);
  book_TH2F("DeltaRmaxvspt_"+histSuffix, 100, 0, 1000, 50, 0., 2.);
  book_TH2F("DeltaR1vspt_"  +histSuffix, 100, 0, 1000, 50, 0., 2.);
  book_TH2F("DeltaR2vspt_"  +histSuffix, 100, 0, 1000, 50, 0., 2.);
  book_TH2F("DeltaR3vspt_"  +histSuffix, 100, 0, 1000, 50, 0., 2.);
  book_TH2F("DeltaR4vspt_"  +histSuffix, 100, 0, 1000, 50, 0., 2.);
  book_TH2F("DeltaRmaxvspt_"+histSuffix+"_OS", 100, 0, 1000, 50, 0., 2.);
  book_TH2F("DeltaR1vspt_"  +histSuffix+"_OS", 100, 0, 1000, 50, 0., 2.);
  book_TH2F("DeltaR2vspt_"  +histSuffix+"_OS", 100, 0, 1000, 50, 0., 2.);
  book_TH2F("DeltaR3vspt_"  +histSuffix+"_OS", 100, 0, 1000, 50, 0., 2.);
  book_TH2F("DeltaR4vspt_"  +histSuffix+"_OS", 100, 0, 1000, 50, 0., 2.);
  return;
}

void GenMatchHists::fill_ParticleHist(const Event & event, const string & histSuffix, const GenParticle& gp) {
  auto weight = event.weight;
  fill_H1("N_"      +histSuffix, 1,           weight);
  fill_H1("pt_"     +histSuffix, gp.pt(),     weight);
  fill_H1("eta_"    +histSuffix, gp.eta(),    weight);
  fill_H1("phi_"    +histSuffix, gp.phi(),    weight);
  fill_H1("mass_"   +histSuffix, gp.v4().M(), weight);
  fill_H1("flavor_" +histSuffix, gp.pdgId(),  weight);


  std::vector<const GenParticle*> genPart;
  std::vector<const GenParticle*> daughters;

  daughters.push_back(gp.daughter(event.genparticles, 1));
  daughters.push_back(gp.daughter(event.genparticles, 2));

  for (auto & dau : daughters) {
    if (!dau) continue;
    if (fabs(dau->pdgId())==24 || fabs(dau->pdgId())==23) {
      const GenParticle* sub1 = dau->daughter(event.genparticles, 1);
      const GenParticle* sub2 = dau->daughter(event.genparticles, 2);
      if (sub1 && fabs(sub1->pdgId())<6 && fabs(sub1->pdgId())>0) genPart.push_back(sub1);
      if (sub2 && fabs(sub2->pdgId())<6 && fabs(sub2->pdgId())>0) genPart.push_back(sub2);
    } else if (fabs( dau->pdgId())<6 && fabs( dau->pdgId())>0) genPart.push_back(dau);
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

  std::vector<double> DeltaRs;
  Particle * baricentrum = new Particle();
  baricentrum->set_eta(x_cm/norm);
  baricentrum->set_phi(y_cm/norm);

  for (unsigned int i = 0; i < genPart.size(); i++) {
    const GenParticle* gp1 = genPart.at(i);
    DeltaRs.push_back(deltaR(*baricentrum,*gp1));
  }

  std::sort(DeltaRs.begin(), DeltaRs.end(), [](double a, double b) {return a < b;});


  if (gp.v4().M()>70) fill_H2("DeltaRmaxvspt_" +histSuffix, gp.pt(), deltaR_max, weight);
  else fill_H2("DeltaRmaxvspt_" +histSuffix+"_OS", gp.pt(), deltaR_max, weight);

  for (size_t i = 0; i < DeltaRs.size(); i++) {
    TString name = "DeltaR"; name += i+1; name += "vspt_" +histSuffix; if (gp.v4().M()<70) name+="_OS";
    fill_H2((string)name, gp.pt(), DeltaRs[i], weight);
  }
  return;
}


void GenMatchHists::book_FamilyTree(const string & ParticleName) {
  book_ParticleHist("mother1"   +ParticleName, "mother1"  +ParticleName, 20,1500);
  book_ParticleHist("mother2"   +ParticleName, "mother2"  +ParticleName, 20,1500);
  book_ParticleHist("daughter1" +ParticleName, "daughter1"+ParticleName, 20,1500);
  book_ParticleHist("daughter2" +ParticleName, "daughter2"+ParticleName, 20,1500);
  book_TH1F("deltaR_daughters"  +ParticleName, "#Delta R(daughters) "+ParticleName, 100, 0, M_PI);
  book_TH1F("deltaR_daughter1"  +ParticleName, "#Delta R(daughter1) "+ParticleName, 100, 0, M_PI);
  book_TH1F("deltaR_daughter2"  +ParticleName, "#Delta R(daughter2) "+ParticleName, 100, 0, M_PI);
}

void GenMatchHists::fill_FamilyTree(const Event & event, const string & ParticleName, const GenParticle& gp) {
  if (gp.mother(event.genparticles, 1))   fill_ParticleHist(event, "mother1"   +ParticleName, *gp.mother(event.genparticles, 1));
  if (gp.mother(event.genparticles, 2))   fill_ParticleHist(event, "mother2"   +ParticleName, *gp.mother(event.genparticles, 2));
  if (gp.daughter(event.genparticles, 1)) fill_ParticleHist(event, "daughter1" +ParticleName, *gp.daughter(event.genparticles, 1));
  if (gp.daughter(event.genparticles, 2)) fill_ParticleHist(event, "daughter2" +ParticleName, *gp.daughter(event.genparticles, 2));
  if (gp.daughter(event.genparticles, 1)) fill_H1("deltaR_daughter1" +ParticleName, deltaR(gp,*gp.daughter(event.genparticles, 1)), event.weight);
  if (gp.daughter(event.genparticles, 2)) fill_H1("deltaR_daughter2" +ParticleName, deltaR(gp,*gp.daughter(event.genparticles, 2)), event.weight);
  if (gp.daughter(event.genparticles, 1) && gp.daughter(event.genparticles, 2)) fill_H1("deltaR_daughters" +ParticleName, deltaR(*gp.daughter(event.genparticles, 1),*gp.daughter(event.genparticles, 2)), event.weight);
}
