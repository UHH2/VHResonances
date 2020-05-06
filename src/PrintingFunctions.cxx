#include "UHH2/VHResonances/include/PrintingFunctions.hpp"



using namespace std;
using namespace uhh2;

void PrintGenPartsTable(Event& event, TString title) {
  cout << title << '\n';
  TableOutput to({"id", "ind", "d1", "d2", "mo1", "mo2", "pt", "eta", "phi"});
  vector<string> v;
  for(const auto & gp : *event.genparticles){
    v.clear();
    v.push_back(to_string(gp.pdgId())); v.push_back(to_string(gp.index()));
    v.push_back(to_string(gp.daughter1())); v.push_back(to_string(gp.daughter2()));
    v.push_back(to_string(gp.mother1())); v.push_back(to_string(gp.mother2()));
    v.push_back(double2string(gp.pt(), 2)); v.push_back(double2string(gp.eta(), 2)); v.push_back(double2string(gp.phi(), 2));
    to.add_row(v);
  }
  to.print(cout);
}


void PrintLeptonTable(Event& event, TString lepton, TString title) {
  vector<Particle> leptons;

  if (lepton == "muons") leptons.assign(event.muons->begin(), event.muons->end());
  else if (lepton == "electrons") leptons.assign(event.electrons->begin(), event.electrons->end());
  else throw logic_error("PrintLeptonTable: Impossible case");

  cout << lepton << "\t" << title << '\n';
  TableOutput to({"pt", "eta", "phi", "E", "mass", "charge"});
  vector<string> v;
  for(const auto & lep : leptons){
    v.clear();
    v.push_back(double2string(lep.pt(), 2)); v.push_back(double2string(lep.eta(), 2)); v.push_back(double2string(lep.phi(), 2));
    v.push_back(double2string(lep.energy(), 2)); v.push_back(double2string(lep.v4().M(), 2)); v.push_back(double2string(lep.charge(), 2));
    to.add_row(v);
  }
  to.print(cout);
}




void PrintJetTable(std::vector<TopJet> & Jets, TString title) {
  // vector<Jet> jets;

  // jets.assign(Jets->begin(), Jets->end());

  // if (collection == "topjets") jets.assign(event.topjets->begin(), event.topjets->end());
  // else if (collection == "jets") jets.assign(event.jets->begin(), event.jets->end());
  // else throw logic_error("PrintLeptonTable: Impossible case");

  cout << title << '\n';
  TableOutput to({"pt", "eta", "phi", "E", "mass",
  "totEnFr", "muonEnFr", "charEmEnFr", "phoEnFr", "neutEmEnFr", "neutHadEnFr", "charHadEnFr",
  "muonMult", "eleMult", "phoMult", "charMult", "neutMult", "HadMult", "nOfDaugh",
  "puppiMult", "neutPuppiMult", "neutHadPuppiMult", "phoPuppiMult", "HFHadPuppiMult", "HFEMPuppiMult",});
  // "btag_CSV", "btag_CSVMVA", "btag_DeepCSV", "btag_BDSVAK8", "btag_BDSVCA15", "JEC_f_raw", "isMatched"});

  vector<string> v;
  for(const auto & jet : Jets){
    // if (abs(jet.muonEnergyFraction()+jet.chargedEmEnergyFraction()+jet.neutralEmEnergyFraction()+jet.neutralHadronEnergyFraction()+jet.chargedHadronEnergyFraction()-1)<0.05) continue;
    v.clear();
    v.push_back(double2string(jet.pt(), 2)); v.push_back(double2string(jet.eta(), 2)); v.push_back(double2string(jet.phi(), 2)); v.push_back(double2string(jet.energy(), 2)); v.push_back(double2string(jet.v4().M(), 2));
    v.push_back(double2string(jet.muonEnergyFraction()+jet.chargedEmEnergyFraction()+jet.neutralEmEnergyFraction()+jet.neutralHadronEnergyFraction()+jet.chargedHadronEnergyFraction(), 2));
    v.push_back(double2string(jet.muonEnergyFraction(), 2)); v.push_back(double2string(jet.chargedEmEnergyFraction(), 2)); v.push_back(double2string(jet.photonEnergyFraction(), 2));
    v.push_back(double2string(jet.neutralEmEnergyFraction(), 2)); v.push_back(double2string(jet.neutralHadronEnergyFraction(), 2)); v.push_back(double2string(jet.chargedHadronEnergyFraction(), 2));
    v.push_back(to_string(jet.muonMultiplicity())); v.push_back(to_string(jet.electronMultiplicity())); v.push_back(to_string(jet.photonMultiplicity()));
    v.push_back(to_string(jet.chargedMultiplicity())); v.push_back(to_string(jet.neutralMultiplicity())); v.push_back(to_string(jet.chargedMultiplicity()+jet.neutralMultiplicity())); v.push_back(to_string(jet.numberOfDaughters()));

    v.push_back(to_string(jet.puppiMultiplicity())); v.push_back(to_string(jet.neutralPuppiMultiplicity())); v.push_back(to_string(jet.neutralHadronPuppiMultiplicity()));
    v.push_back(to_string(jet.photonPuppiMultiplicity())); v.push_back(to_string(jet.HFHadronPuppiMultiplicity())); v.push_back(to_string(jet.HFEMPuppiMultiplicity()));

    // v.push_back(double2string(jet.btag_combinedSecondaryVertex(), 2)); v.push_back(double2string(jet.btag_combinedSecondaryVertexMVA(), 2)); v.push_back(double2string(jet.btag_DeepCSV(), 2));
    // v.push_back(double2string(jet.btag_BoostedDoubleSecondaryVertexAK8(), 2)); v.push_back(double2string(jet.btag_BoostedDoubleSecondaryVertexCA15(), 2)); v.push_back(double2string(jet.JEC_factor_raw(), 2));
    // if (jet.has_tag(TopJet::isMatched)) v.push_back(double2string((float)jet.get_tag(TopJet::isMatched), 2));
    // else {v.push_back("NF");}
    to.add_row(v);
  }
  to.print(cout);
}
