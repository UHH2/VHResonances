#include <iostream>
#include <TROOT.h>
#include "TSystem.h"
#include <stdlib.h>
#include <stdio.h>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TFile.h"
#include "TLatex.h"
#include "UHH2/VHResonances/include/constants.hpp"
#include "UHH2/VHResonances/Analysis/include/tdrstyle_all.hpp"

typedef std::pair<std::string, int> mypair_I;

double EfficiencyError(double count, double N);
bool isCSRegion(std::string tag);
void CalculateSignalEfficiencies();
