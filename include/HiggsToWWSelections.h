#pragma once

#include <cmath>
#include "UHH2/core/include/fwd.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Selection.h"
#include "UHH2/common/include/JetIds.h"
#include "UHH2/VHResonances/include/constants.hpp"

class JetDiLeptonPhiAngularSelection: public uhh2::Selection {
public:
  JetDiLeptonPhiAngularSelection(float pt_min, float phi_min, float phi_max, TString lepton, const uhh2::Event::Handle<std::vector<TopJet> > & topjetcollection );
  virtual bool passes(const uhh2::Event& event) override;
private:
  float pt_min, phi_min, phi_max;
  TString lepton;
  uhh2::Event::Handle<std::vector<TopJet> > topjetcollection;
};



class DeltaRDiLepton: public uhh2::Selection {
public:
  DeltaRDiLepton(float DR_min, float DR_max, TString lepton);
  virtual bool passes(const uhh2::Event& event) override;
private:
  float DR_min, DR_max;
  TString lepton;
};



class DeltaPhiMET: public uhh2::Selection {
public:
  DeltaPhiMET(float Dphi_min_, bool is_invisiblechannel_, const uhh2::Event::Handle<std::vector<TopJet> > & topjetcollection_);
  virtual bool passes(const uhh2::Event& event) override;
private:
  float Dphi_min;
  bool is_invisiblechannel;
  uhh2::Event::Handle<std::vector<TopJet> > topjetcollection;
};



typedef std::function<bool (const ZprimeCandidate &, const uhh2::Event &)> ZprimeCandidate_ID;

class ZprimeCandidateID {
public:

  explicit ZprimeCandidateID(const uhh2::Event::Handle<std::vector<ZprimeCandidate> > & h_ZprimeCandidates);
  bool operator()(const ZprimeCandidate& cand, const uhh2::Event & event) const;
private:

  uhh2::Event::Handle< std::vector<ZprimeCandidate> > h_ZprimeCandidates;

};


class PTMassCut: public uhh2::Selection {
public:
  PTMassCut(float cut_min, const uhh2::Event::Handle<std::vector<ZprimeCandidate> > & h_ZprimeCandidates);
  virtual bool passes(const uhh2::Event& event) override;
private:
  float cut_min;
  uhh2::Event::Handle< std::vector<ZprimeCandidate> > h_ZprimeCandidates;

};



class TaggerCut: public uhh2::Selection {
public:
  TaggerCut(float cut_min, float cut_max, float pt_min, std::string tagger, const uhh2::Event::Handle<std::vector<ZprimeCandidate> > & topjetcollection);
  virtual bool passes(const uhh2::Event& event) override;
private:
  float cut_min, cut_max, pt_min;
  std::string tagger;
  uhh2::Event::Handle<std::vector<ZprimeCandidate> > h_ZprimeCandidates;
};


class MultiBTagSubJetID {
public:

  explicit MultiBTagSubJetID(JetId SubjetId, int nsubjets =2);
  bool operator()(const TopJet & topjet, const uhh2::Event & event) const;

private:
  JetId m_SubjetId;
  int m_nsubjets;
};




class BlindDataSelection: public uhh2::Selection {
public:
  BlindDataSelection(uhh2::Context& ctx);
  virtual bool passes(const uhh2::Event& event) override;
private:
  uhh2::Event::Handle<bool> h_is_Blind;
  bool is_BlindSelection, is_mc;
};
