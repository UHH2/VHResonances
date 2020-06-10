#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/common/include/ObjectIdUtils.h"
#include "UHH2/common/include/ElectronIds.h"
#include "UHH2/common/include/MuonIds.h"
#include "UHH2/VHResonances/include/HistsBase.hpp"
#include "UHH2/VHResonances/include/constants.hpp"

class LeptonIDHists: public HistsBase {
public:

	LeptonIDHists(uhh2::Context&, const std::string&);
	virtual void fill(const uhh2::Event&) override;
	virtual ~LeptonIDHists();

private:
	std::unordered_map<std::string, uhh2::Event::Handle<std::vector<Muon> > > h_muon;
	std::unordered_map<std::string, uhh2::Event::Handle<std::vector<Electron> > > h_ele;

	std::unordered_map<std::string, ElectronId> eleIds;
  std::unordered_map<std::string, MuonId> muoIds;
	std::vector<std::string> leptoncollections;

};
