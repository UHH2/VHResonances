#pragma once
#include <iostream>
#include <TROOT.h>
#include "TSystem.h"
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TFile.h"
#include "TLatex.h"
#include "UHH2/core/include/Utils.h"
#include "UHH2/HiggsToWWTagger/include/constants.hpp"
#include "UHH2/HiggsToWWTagger/Analysis/include/tdrstyle_all.hpp"

// double plot_hi = 2800;
extern const double plot_lo = 300;
extern const double plot_hi = 8200;
extern const unsigned int nPoints=MassPoints.size();


extern const std::vector<double> theo_xsec     = {1.089e+04,3.415e+03,1.505e+03,7.783e+02,4.387e+02,2.620e+02,1.626e+02,1.038e+02,3.640e+01,1.371e+01,5.330e+00,2.133e+00,8.652e-01,3.604e-01,1.535e-01,6.639e-02,0,0};
extern const std::vector<double> theo_xsec_err = {6.314e+00,4.186e+01,9.512e-01,4.371e-01,2.517e-01,1.494e-01,9.351e-02,5.563e-02,2.002e-02,7.165e-03,2.522e-03,9.392e-04,4.011e-04,1.515e-04,6.865e-05,2.653e-05,0,0};

extern const std::vector<double> MassPoints_Hbb = { 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000, 3100, 3200, 3300, 3400, 3500, 3600, 3700, 3800, 3900, 4000, 4100, 4200, 4300, 4400, 4500, 4600, 4700, 4800, 4900, 5000};
extern const std::vector<double> expectedHbb    = { 18.675, 10.9125, 7.3125, 5.5125, 4.2937, 3.4781, 2.8594, 2.4281, 2.0719, 1.7906, 1.5609, 1.3781, 1.2328, 1.1156, 1.0125, 0.9281, 0.8578, 0.8016, 0.7547, 0.7125, 0.6797, 0.6539, 0.6328, 0.6141, 0.6, 0.5906, 0.5813, 0.5766, 0.5813, 0.5859, 0.5953, 0.6094, 0.6234, 0.6492, 0.6797, 0.7172, 0.7594, 0.8156, 0.8203, 0.8297, 0.8438, 0.8578, 0.8766};
extern const std::vector<double> expectedHbb0b  = {14.325, 15.8625, 17.85, 20.325, 17.625, 14.625, 11.9625, 9.7312, 7.9313, 6.3938, 5.2875, 4.4437, 3.7875, 3.2625, 2.8312, 2.4844, 2.1844, 1.9312, 1.7062, 1.5141, 1.35, 1.2047, 1.0828, 0.9797, 0.893, 0.8203, 0.757, 0.7031, 0.6609, 0.6281, 0.5977, 0.5719, 0.5484, 0.5344, 0.5238, 0.518, 0.5156, 0.5203, 0.4945, 0.4734, 0.457, 0.4453, 0.4336};

void PlotLimits(bool doObs = false, bool isHbb = false);
