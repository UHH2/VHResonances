#include <UHH2/VHResonances/include/HistsBase.hpp>

using namespace std;
using namespace uhh2;

void HistsBase::book_TH1F(const string& name_, const int xnbins_, const double xmin_, const double xmax_){
  if(h1.find(name_) == h1.end()) h1[name_] = book<TH1F>(name_.c_str(), name_.c_str(), xnbins_, xmin_, xmax_);
  else throw runtime_error("HistsBase::book_TH1F -- histogram key already existing: "+name_);
  return;
}

void HistsBase::book_TH1F(const string& name_,const string& title_, const int xnbins_, const double xmin_, const double xmax_){
  if(h1.find(name_) == h1.end()) h1[name_] = book<TH1F>(name_.c_str(), title_.c_str(), xnbins_, xmin_, xmax_);
  else throw runtime_error("HistsBase::book_TH1F -- histogram key already existing: "+name_);
  return;
}

void HistsBase::book_TH1F(const string& name_, const string& title_, const int xnbins_, const double* xbin_ls_){
  if(h1.find(name_) == h1.end()) h1[name_] = book<TH1F>(name_.c_str(), title_.c_str(), xnbins_, xbin_ls_);
  else throw runtime_error("HistsBase::book_TH1F -- histogram key already existing: "+name_);
  return;
}

void HistsBase::book_TH1F(const string& name_, const string& title_, const int xnbins_, const float* xbin_ls_){
  if(h1.find(name_) == h1.end()) h1[name_] = book<TH1F>(name_.c_str(), title_.c_str(), xnbins_, xbin_ls_);
  else throw runtime_error("HistsBase::book_TH1F -- histogram key already existing: "+name_);
  return;
}


void HistsBase::book_TH2F(const string& name_, const int xnbins_, const double xmin_, const double xmax_, const int ynbins_, const double ymin_, const double ymax_){
  if(h2.find(name_) == h2.end()) h2[name_] = book<TH2F>(name_.c_str(), name_.c_str(), xnbins_, xmin_, xmax_, ynbins_, ymin_, ymax_);
  else throw runtime_error("HistsBase::book_TH2F -- histogram key already existing: "+name_);
  return;
}

void HistsBase::book_TH2F(const string& name_, const string& title_, const int xnbins_, const double xmin_, const double xmax_, const int ynbins_, const double ymin_, const double ymax_){
  if(h2.find(name_) == h2.end()) h2[name_] = book<TH2F>(name_.c_str(), title_.c_str(), xnbins_, xmin_, xmax_, ynbins_, ymin_, ymax_);
  else throw runtime_error("HistsBase::book_TH2F -- histogram key already existing: "+name_);
  return;
}

void HistsBase::book_TH2F(const string& name_, const int xnbins_, const double* xbin_ls_, const int ynbins_, const double* ybin_ls_){
  if(h2.find(name_) == h2.end()) h2[name_] = book<TH2F>(name_.c_str(), name_.c_str(), xnbins_, xbin_ls_, ynbins_, ybin_ls_);
  else throw runtime_error("HistsBase::book_TH2F -- histogram key already existing: "+name_);
  return;
}

TH1F* HistsBase::H1(const string& key_){
  TH1F* h(0);
  if(h1.find(key_) != h1.end()) h = h1[key_];
  else throw runtime_error("HistsBase::H1 -- histogram key not found: "+key_);
  return h;
}

TH2F* HistsBase::H2(const string& key_){
  TH2F* h(0);
  if(h2.find(key_) != h2.end()) h = h2[key_];
  else throw runtime_error("HistsBase::H2 -- histogram key not found: "+key_);
  return h;
}

void HistsBase::fill_H1(const string& key_, double content, double weight){
  TH1F* h = H1(key_);
  h->Fill(content, weight);
  return;
}

void HistsBase::fill_H2(const string& key_, double content1, double content2, double weight){
  TH2F* h = H2(key_);
  h->Fill(content1, content2, weight);
  return;
}
