#pragma once

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/TopJet.h"
#include "UHH2/VHResonances/include/constants.hpp"

Matching FloatToMatching(const TopJet & jet);
MatchingStatus FloatToMatchingStatus(const TopJet & jet);

bool isLeptonic(int pdgId);
bool isHadronic(int pdgId);
bool DobleDecay(int pdgId1, int pdgId2, Decay decay);

bool FindInString(const std::string& search, const std::string& str);
