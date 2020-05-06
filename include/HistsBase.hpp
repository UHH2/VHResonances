#pragma once

#include <string>
#include <unordered_map>

#include <TH1F.h>
#include <TH2F.h>

#include <UHH2/core/include/Hists.h>
#include <UHH2/core/include/Event.h>

class HistsBase : public uhh2::Hists {

public:
  explicit HistsBase(uhh2::Context& ctx, const std::string& dirname): uhh2::Hists(ctx, dirname) {}

protected:

  void book_TH1F(const std::string&, const int, const double, const double);
  void book_TH1F(const std::string&, const std::string&, const int, const double, const double);
  void book_TH1F(const std::string&, const std::string&, const int, const double*);
  void book_TH1F(const std::string&, const std::string&, const int, const float*);

  void book_TH2F(const std::string&, const int, const double, const double, const int, const double, const double);
  void book_TH2F(const std::string&, const std::string&, const int, const double, const double, const int, const double, const double);
  void book_TH2F(const std::string&, const int, const double*, const int, const double*);

  TH1F* H1(const std::string&);
  TH2F* H2(const std::string&);

  void fill_H1(const std::string&, double, double);
  void fill_H2(const std::string&, double, double, double);

  std::unordered_map<std::string, TH1F*> h1;
  std::unordered_map<std::string, TH2F*> h2;
};
