#ifndef SIG_PLOTTER_H
#define SIG_PLOTTER_H

#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <utility>
#include <map>
#include <iomanip>
#include "AtlasStyle.h"
#include "TMultiGraph.h"
#include "TGraph.h"
#include "TLine.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TFile.h"
#include "TDirectoryFile.h"
#include "TH1D.h"
using namespace std;

class sig_plotter{
 public:
  //worker functions shared by tmva_li_train
  sig_plotter(bool ks);

  void do_plots(const string&sigdir,const string&pdir)const;

  //tr_opt 0 (false) use MET; 1 for inv_ll variables; 2 (default) for inv_bb; 3 (both); else std
  bool is_bkg(TString bg,int which=-1)const;
  TString li_tag(int which)const;
  TString bg_tag(int bg)const;
  TString bg_label(int bg)const;
  std::pair<TString,int> var_type(const TString& var)const;
  std::vector<TString> sorted_vars(const std::vector<TString>&vars)const;
  std::pair<std::vector<TString>,std::vector<TString>> li_vars(int which)const;
  int n_bg_cases()const{return 2;}
  int n_li_cases()const{return 6;}
  int n_std_min()const{return -1;}

 protected:
  TString ith_var_tag(int vt,int iv)const;
  TH1D* bdt_dist(const string&file,int sbg)const;
  void print_bdts(TCanvas*c,const string&ddir,const string&odir,int vt,int bg,int iv=-1)const;
  pair<double,double> optimal_sig(TH1D*sig,TH1D*big)const;
  pair<double,double> sig_from_files(const string&ddir,int vt,int bg,int iv)const;
  double sig_from_log(const string&fname)const;
  void auto_limits(vector<TH1D*> loki,double zero=0.)const;
  TStyle* set_style()const;
  TLatex plot_latex()const;
 private:
  string bdt_dir;
  bool keep_spectators;
};

#endif
