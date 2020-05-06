#pragma once

#include <stdexcept>

#include "UHH2/core/include/Utils.h"
#include "UHH2/core/include/Event.h"


void PrintGenPartsTable(uhh2::Event& event, TString title="");
void PrintLeptonTable(uhh2::Event& event, TString lepton, TString title="");
void PrintJetTable(std::vector<TopJet> & Jets, TString title="");
