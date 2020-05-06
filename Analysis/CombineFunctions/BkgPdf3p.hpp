/*****************************************************************************
* Project: RooFit                                                           *
*                                                                           *
* This code was autogenerated by RooClassFactory                            *
*****************************************************************************/

#pragma once

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"

class BkgPdf3p : public RooAbsPdf {
public:
  BkgPdf3p() {} ;
  BkgPdf3p(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _p0, RooAbsReal& _p1, RooAbsReal& _p2);
  BkgPdf3p(const BkgPdf3p& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new BkgPdf3p(*this,newname); }
  inline virtual ~BkgPdf3p() { }

  /*
  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;
  */

protected:

  RooRealProxy x;
  RooRealProxy p0;
  RooRealProxy p1;
  RooRealProxy p2;

  Double_t evaluate() const ;

private:

  ClassDef(BkgPdf3p,2)  // PDF of dijet function with 4 free parameters
};
