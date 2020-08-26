#pragma once

#include <string>
#include <iostream>
#include <memory>
#include <unordered_map>


#include "UHH2/core/include/Particle.h"
#include "UHH2/core/include/TopJet.h"

class ZprimeCandidate : public Particle {

public:

  ZprimeCandidate(){};

  // Getter
  double Zprime_mass() const{return m_Zprime_mass;}
  Particle Z() const{return m_Z;}
  TopJet H() const{return m_H;}
  std::vector<Particle> leptons() const{return m_leptons;}

  double discriminator(const std::string& key) const {
    auto it = m_discriminators.find(key);
    if(it == m_discriminators.end()) throw std::runtime_error("ZprimeCandidate::discriminator(): discriminator with key '" + key + "' not set");
    return it->second;
  }

  bool has_discriminator(const std::string& key) const {return m_discriminators.find(key) != m_discriminators.end();}

  // Setters
  void set_Zprime(double pt, double eta, double phi, double energy) {
    this->set_charge(0);
    this->set_pt(pt);
    this->set_eta(eta);
    this->set_phi(phi);
    this->set_energy(energy);
    m_Zprime_mass = this->v4().M();
  }
  void set_Zprime(LorentzVector candidate) {
    this->set_charge(0);
    this->set_pt(candidate.Pt());
    this->set_eta(candidate.Eta());
    this->set_phi(candidate.Phi());
    this->set_energy(candidate.energy());
    m_Zprime_mass = this->v4().M();
  }
  void set_Z(LorentzVector candidate) {
    m_Z.set_charge(0);
    m_Z.set_pt(candidate.Pt());
    m_Z.set_eta(candidate.Eta());
    m_Z.set_phi(candidate.Phi());
    m_Z.set_energy(candidate.energy());
  }
  void set_H(LorentzVector candidate) {
    m_H.set_charge(0);
    m_H.set_pt(candidate.Pt());
    m_H.set_eta(candidate.Eta());
    m_H.set_phi(candidate.Phi());
    m_H.set_energy(candidate.energy());
  }
  void set_Z(Particle x) {m_Z=x;}
  void set_H(TopJet x) {m_H=x;}
  void set_jets_leptonic(std::vector<Particle> x) {m_leptons=x;}
  void set_Zprime_MT (double x) {m_Zprime_mass = x;}
  void set_discriminators(const std::string& key, double discr) { m_discriminators[key] = discr;}

private:

  double m_Zprime_mass;
  Particle m_Z;
  TopJet m_H;
  std::vector<Particle> m_leptons;
  std::unordered_map<std::string, double > m_discriminators;

};
