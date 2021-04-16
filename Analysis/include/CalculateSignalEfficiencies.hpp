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

const std::vector<double> MyMassPoints = {1000, 1200, 1400, 1600, 1800, 2000, 2500, 3000, 3500, 4000, 4500, 5000};

double EfficiencyError(double count, double N);
bool isCSRegion(std::string tag);
void CalculateSignalEfficiencies(std::string histFolder ="DeepAk8_ZHccvsQCD_MD");


int FindInVector(const std::vector<std::string>& vec, const std::string& el);
