#pragma once

#include <iostream>
#include <memory>
#include <unordered_map>

#include <UHH2/core/include/AnalysisModule.h>
#include <UHH2/core/include/Event.h>
#include <UHH2/core/include/Hists.h>
#include "UHH2/core/include/Selection.h"

class ModuleBASE: public uhh2::AnalysisModule {

public:
  explicit ModuleBASE() {}
  explicit ModuleBASE(uhh2::Context&) {}
  virtual bool process(uhh2::Event&);

protected:
  typedef std::unordered_map<std::string, std::unique_ptr<uhh2::Hists> >::iterator hfolder_itr;
  std::unordered_map<std::string, std::unique_ptr<uhh2::Hists> > hfolder_;

  std::unordered_map<std::string, std::shared_ptr<uhh2::Event::Handle< float > > > WFolder_;
  std::unordered_map<std::string, std::shared_ptr<uhh2::Event::Handle< double > > > WWFolder_;

  template<typename F> void book_HFolder(const std::string&, F*);
  template<typename F=uhh2::Hists> F* HFolder(const std::string&);

  template<typename F> void book_WFolder(const std::string&, F*, F);
  template<typename F=uhh2::Event::Handle<float>> F WFolder(const std::string&);

  template<typename F> void book_WWFolder(const std::string&, F*, F);
  template<typename F=uhh2::Event::Handle<double>> F WWFolder(const std::string&);
};


/*
##     ## ####  ######  ########  #######   ######   ########     ###    ##     ##  ######
##     ##  ##  ##    ##    ##    ##     ## ##    ##  ##     ##   ## ##   ###   ### ##    ##
##     ##  ##  ##          ##    ##     ## ##        ##     ##  ##   ##  #### #### ##
#########  ##   ######     ##    ##     ## ##   #### ########  ##     ## ## ### ##  ######
##     ##  ##        ##    ##    ##     ## ##    ##  ##   ##   ######### ##     ##       ##
##     ##  ##  ##    ##    ##    ##     ## ##    ##  ##    ##  ##     ## ##     ## ##    ##
##     ## ####  ######     ##     #######   ######   ##     ## ##     ## ##     ##  ######
*/



template<typename F>
void ModuleBASE::book_HFolder(const std::string& dname_, F* f_ptr){

  if(hfolder_.find(dname_) != hfolder_.end()) throw std::runtime_error("ModuleBASE::book_HFolder -- already existing folder key: "+dname_);
  else hfolder_[dname_].reset(f_ptr);
  return;
}

template<typename F>
F* ModuleBASE::HFolder(const std::string& key_){
  F* f(0);
  if(hfolder_.find(key_) != hfolder_.end()) f = static_cast<F*>(hfolder_[key_].get());
  else throw std::runtime_error("ModuleBASE::HFolder -- folder key not found: "+key_);
  if(!f) throw std::runtime_error("ModuleBASE::HFolder -- null pointer to folder: "+key_);
  return f;
}


/*
##     ##    ###    ##    ## ########  ##       ########  ######
##     ##   ## ##   ###   ## ##     ## ##       ##       ##    ##
##     ##  ##   ##  ####  ## ##     ## ##       ##       ##
######### ##     ## ## ## ## ##     ## ##       ######    ######
##     ## ######### ##  #### ##     ## ##       ##             ##
##     ## ##     ## ##   ### ##     ## ##       ##       ##    ##
##     ## ##     ## ##    ## ########  ######## ########  ######
*/



template<typename F>
void ModuleBASE::book_WFolder(const std::string& dname_, F* f_ptr, F f_){

  if(WFolder_.find(dname_) != WFolder_.end()) throw std::runtime_error("ModuleBASE::book_WFolder -- already existing folder key: "+dname_);
  else { WFolder_[dname_].reset(f_ptr); *WFolder_[dname_] = f_;}
  return;
}

template<typename F>
F ModuleBASE::WFolder(const std::string& key_){
  F* f(0);
  if(WFolder_.find(key_) != WFolder_.end()) f = static_cast<F*>(WFolder_[key_].get());
  else throw std::runtime_error("ModuleBASE::WFolder -- folder key not found: "+key_);
  if(!f) throw std::runtime_error("ModuleBASE::WFolder -- null pointer to folder: "+key_);
  return *f;
}
