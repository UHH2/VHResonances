#pragma once

#include <UHH2/core/include/AnalysisModule.h>
#include <UHH2/core/include/Event.h>

#include "UHH2/common/include/JetIds.h"
#include "UHH2/common/include/PrintingModules.h"
#include "UHH2/common/include/MCWeight.h"

#include "UHH2/VHResonances/include/Utils.hpp"
#include "UHH2/VHResonances/include/ModuleBase.h"
#include "UHH2/VHResonances/include/ZprimeCandidate.h"
#include "UHH2/VHResonances/include/HiggsToWWSelections.h"
#include "UHH2/VHResonances/include/GeneralizedEndpoint.hpp"

class FinalStateMatching: public uhh2::AnalysisModule {
public:
  explicit FinalStateMatching(uhh2::Context & ctx);

  virtual bool process(uhh2::Event & event) override;

private:
  uhh2::Event::Handle<float> h_ZDecay_;
  uhh2::Event::Handle<float> h_HDecay_;
  uhh2::Event::Handle<float> h_ZprimeDecay_;
  bool skipMatching=false;
  std::unique_ptr<uhh2::AnalysisModule> GenParticles_printer;

};



class ZprimeCandidateReconstruction : public uhh2::AnalysisModule {

public:
  explicit ZprimeCandidateReconstruction(uhh2::Context& ctx, float pt_min, float DR_min, float DR_max, float phi_min, float phi_max, const std::string& lepton, const std::string& topjetcollection);
  virtual bool process(uhh2::Event&) override;
  virtual void setDiscriminators(uhh2::Event&, ZprimeCandidate& candiate,
    ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<double> > ZplusJet,
    ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<double> > diLep,
    Particle lep1, Particle lep2, bool MuonID1, bool MuonID2, TopJet jet,
    unsigned int i, unsigned int j, std::map<TString, JetId> Btag_map);

private:

  uhh2::Event::Handle< std::vector<TopJet> > h_topjets;
  uhh2::Event::Handle< std::vector<ZprimeCandidate> > h_ZprimeCandidates_;

  float pt_min, DR_min, DR_max, phi_min, phi_max;
  std::string lepton;
  std::string topjetcollection;

};

class SDMassCalculator: public uhh2::AnalysisModule {
public:
  explicit SDMassCalculator(uhh2::Context & ctx, const std::string & jetCollName="topjets");
  virtual ~SDMassCalculator() {};
  virtual bool process(uhh2::Event&) override;
  float calcSDmass(const TopJet & jet);
private:
  uhh2::Event::Handle<std::vector<TopJet>> h_topjets_;
  uhh2::Event::Handle<std::vector<GenTopJet>> h_gentopjets_;
};


class BlindData : public uhh2::AnalysisModule {

public:
  explicit BlindData(uhh2::Context& ctx, const uhh2::Event::Handle<std::vector<ZprimeCandidate> >& h_ZprimeCandidates);
  virtual bool process(uhh2::Event&) override;

private:

  uhh2::Event::Handle< std::vector<ZprimeCandidate> > h_ZprimeCandidates_;
  uhh2::Event::Handle<bool> h_is_Blind;

};

// Generic Class for Applying SFs
class ScaleFactorsFromHistos : public uhh2::AnalysisModule {

public:
  void LoadHisto(TFile* file, std::string name, std::string hname);
  double Evaluator(std::string hname, double var);

protected:
  std::unordered_map<std::string, std::unique_ptr<TH1F> > histos;

};

// Apply Theory weights //TODO make multiple inheritance
class NLOCorrections : public ScaleFactorsFromHistos {

public:
  explicit NLOCorrections(uhh2::Context& ctx);
  virtual bool process(uhh2::Event&) override;
  double GetPartonObjectPt(uhh2::Event& event, ParticleID objID);

private:
  bool is_Wjets, is_Zjets, is_DY, is_Znn, is2017;

};




// Generic Class for Applying SFs
class ScaleFactorsManager : public uhh2::AnalysisModule {

public:
  explicit ScaleFactorsManager(uhh2::Context & ctx, const uhh2::Event::Handle<std::vector<ZprimeCandidate> >&);
  virtual bool process(uhh2::Event&) override;

protected:
  uhh2::Event::Handle< std::vector<ZprimeCandidate> > h_ZprimeCandidates_;
  std::unordered_map<std::string, std::unique_ptr<MCMuonScaleFactor> > SFs_muo;
  std::unordered_map<std::string, std::unique_ptr<MCElecScaleFactor> > SFs_ele;
  bool muonchannel;
  bool electronchannel;

};



// Generic Class for Applying Muon Scale Variations
class MuonScaleVariations : public uhh2::AnalysisModule {

public:
  explicit MuonScaleVariations(uhh2::Context & ctx);
  virtual bool process(uhh2::Event&) override;

private:
  int mode;
  std::unique_ptr<GeneralizedEndpoint> GE;

};
