#pragma once
#include <iostream>
#include <map>

// Taken from:
// https://github.com/cms-analysis/SUSYBSMAnalysis-Zprime2muAnalysis/blob/mini-AOD/src/GeneralizedEndpoint.h
class GeneralizedEndpoint {
 public:
   GeneralizedEndpoint();
   virtual ~GeneralizedEndpoint();
   float GeneralizedEndpointPt(float MuonPt, int MuonCharge, float MuonEta, float MuonPhi, int Mode, int verbose=0);
 private:
   std::map<int,std::map<int,float> > _Correction;
   std::map<int,std::map<int,float> > _CorrectionError;

 };
