#pragma once

#include "TROOT.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TPad.h"
#include "TLatex.h"
#include "TLine.h"
#include "TBox.h"
#include "TASImage.h"
#include "TPaveLabel.h"
#include "TLegendEntry.h"
#include "TFrame.h"
#include "TCanvas.h"
#include "TH1D.h"

#include <iostream>
#include <algorithm>

//////////////////////////////////////////
// New CMS Style from 2014              //
// https://ghm.web.cern.ch/ghm/plots/   //
// Merged all macros into one           //
//////////////////////////////////////////

const std::string red("\x1b[0;31m");
const std::string green("\x1b[0;32m");
const std::string yellow("\x1b[0;33m");
const std::string cyan("\x1b[0;36m");
const std::string magenta("\x1b[0;35m");
const std::string reset("\x1b[0m");

extern float cmsTextFont;
extern TString cmsText;
extern bool writeExtraText;
extern TString extraText;
extern TString extraText2;
extern float extraTextFont;
extern float lumiTextSize;
extern float lumiTextOffset;
extern float cmsTextSize;
extern float cmsTextOffset;
extern float relPosX;
extern float relPosY;
extern float relExtraDY;
extern float extraOverCmsTextSize;
extern TString lumi_13TeV;
extern TString lumi_8TeV;
extern TString lumi_7TeV;
extern TString lumi_sqrtS;
extern bool drawLogo;
extern bool kSquare;
extern bool kRectangular;



void CMSOff();
void ExtraTextOff();

void ForThesis();
void WIP(bool isSimulation=false);
void SetSimulation();
void SetYear(TString year);

// tdrGrid: Turns the grid lines on (true) or off (false)
void tdrGrid(bool gridOn);
// fixOverlay: Redraws the axis
void fixOverlay();
void setTDRStyle();
void CMS_lumi( TPad* pad, int iPeriod=4, int iPosX=11 );

////////////////////
// General Macro	//
////////////////////

TCanvas* tdrCanvas(const char* canvName, double x_min, double x_max, double y_min, double y_max, const char* nameXaxis, const char* nameYaxis, bool square = kRectangular, int iPeriod = 4, int iPos = 11);
void tdrCanvasSetAxes(TCanvas *canv, double x_min, double x_max, double y_min, double y_max);
TCanvas* tdrDiCanvas(const char* canvName, double x_min, double x_max, double y_min, double y_max, double y_min2, double y_max2,  const char* nameXaxis, const char* nameYaxis, const char* nameYaxis2, bool square = kRectangular, int iPeriod = 4, int iPos = 11);
void setNewTitle(TString name = "new title");
void tdrDraw(TH1* h, std::string opt, int marker=kFullCircle, int mcolor = kBlack, int lstyle=kSolid, int lcolor=-1, int fstyle=1001, int fcolor=kYellow+1);
void tdrDraw(TGraph* g, std::string opt, int marker=kFullCircle, int mcolor = kBlack, int lstyle=kSolid, int lcolor=-1, int fstyle=1001, int fcolor=kYellow+1);
TLegend *tdrLeg(double x1, double y1, double x2, double y2, double textSize=0.045, int textFont=42, int textColor=kBlack);
void tdrHeader(TLegend* leg, TString legTitle, int textAlign = 12, double textSize = 0.04, int textFont = 42, int textColor = kBlack, bool isToRemove = true);
