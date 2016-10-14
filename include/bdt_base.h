#ifndef BDT_BASE_H
#define BDT_BASE_H

#include <dirent.h>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TColor.h"
#include "TSystemDirectory.h"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <utility>
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
#include "TH2D.h"
#include "TFitResultPtr.h"
#include "HistoTransform.h"
using namespace std;

const bool always_tag=true;

class bdt_base{
 public:
  bdt_base(const string&in,const string&out,const string&tag,int nt=2,int nj=2,int eff_cut=77,int ratio=20,bool cts=false,int lptv=1,bool use_btvs=true,bool wmll=true,bool met=false);
  pair<double,double>cut_sig(TH1D*sig,TH1D*big,double cut)const;
  pair<double,double>cut_sig(pair<TH1D*,TH1D*>sbg,double cut)const{return cut_sig(sbg.first,sbg.second,cut);}
  pair<double,double> cumulat_sig(TH1D*sig,TH1D*big)const;
  pair<double,double> cumulat_sig(pair<TH1D*,TH1D*>sbg)const{return cumulat_sig(sbg.first,sbg.second);}
  pair<double,double> optimal_sig(TH1D*sig,TH1D*big,bool simple=false)const;
  pair<double,double> optimal_sig(pair<TH1D*,TH1D*>sbg,bool simple=false)const{return optimal_sig(sbg.first,sbg.second,simple);}
  int which_sam(TString sam)const;
  TString sam_tag(int sam)const;
  TString sam_label(int sam)const;
  int sam_min()const{return 0;}
  int n_sam()const{return 7;}
  void print_info_plot_2d(const vector<vector<double> >&info,const vector<TString>&xlabels,const vector<TString>&ylabels,const TString&name,const TString&pdir,const TString&title)const;
  void set_tag(const string&tag){identifier=tag;}
  TString filename(const vector<TString>&vars,int bg=-1)const{return tmva_out_dir+"TMVA"+sam_tag(bg==0?-1:bg)+ftag(vars)+".root";}

 protected:
  TString ftag(const vector<TString>&vars)const;

  bool is_single_bg(int bg)const{return bg>1&&bg<n_sam();}
  bool process_sample(int bg,int sam)const;
  pair<TH1D*,TH1D*> transformation_DF(TH1D*sig,TH1D*big,bool d,bool husk=true,double zs=4.5,double zb=4.5)const;
  pair<TH1D*,TH1D*> transformation_DF(pair<TH1D*,TH1D*>sbg,bool d,bool husk=true,double zs=4.5,double zb=4.5)const{return transformation_DF(sbg.first,sbg.second,d,husk,zs,zb);}

  TString these_cuts()const{return cuts(ntag,njet,beff_cut,charm_ratio,cts_mv2,loptv,widemll,metcut);}
  TString cuts(int nt,int nj,int cut,int ratio,bool cts,bool lptv,bool wmll,bool met)const;
  TString btag_str()const{return "-"+btag_var()+cut_str();}
  TString btag_label()const{return btag_var()+" "+cut_label();}
  TString btag_var(int ratio,bool cts)const;
  TString cut_str(int cut,int ratio)const;
  TString cut_label(int cut,int ratio)const;
  TString region_str(int nt,int nj,bool lptv)const;
  TString region_label(int nt,int nj,bool lptv)const;
  TString tdf_tstr(int tDF)const;
  string tdf_str(int tDF)const{return string(tdf_tstr(tDF));}
  double mv2c_cut(int cut,int ratio)const;

  //default versions using class members
  TString btag_var()const{return btag_var(charm_ratio,cts_mv2);}
  TString cut_str()const{return cut_str(beff_cut,charm_ratio);}
  TString cut_label()const{return cut_label(beff_cut,charm_ratio);}
  TString region_str()const{return region_str(ntag,njet,loptv);}
  TString region_label()const{return region_label(ntag,njet,loptv);}
  double mv2c_cut()const{return mv2c_cut(beff_cut,charm_ratio);}
  
  //scut work
  string check_dir(const string&dir)const;//ensures strings end with a "/"; less legwork later
  string check_isdir(const string&dir)const;//check_dir and aborts if directory doesn't exist
  string check_mkdir(const string&dir)const;//check_dir and makes directory if it doesn't exist
  TString check_dir(const TString&dir)const{return TString(check_dir(string(dir.Data())));}//ensures strings end with a "/"; less legwork later
  TString check_mkdir(const TString&dir)const{return TString(check_mkdir(string(dir.Data())));}//ensures strings end with a "/"; less legwork later
  bool file_exists(const string&name)const;
  int color(int which)const;
  TStyle* set_style()const;
  TLatex plot_latex()const;
  void auto_limits(vector<TH1D*> loki,double zero=0.)const;

  //mutators--for modifying cuts
  void set_njet(int nj){njet=nj;}
  void set_ntag(int nt){ntag=nt;}
  void set_ptv(bool pt){loptv=pt;}
  void set_cut(int cut){beff_cut=cut;}
  void set_ratio(int ratio){charm_ratio=ratio;}
  void set_cts(bool cts){cts_mv2=cts;}

  //mutators for in and out directories (btagging directory structure looks different because of truth tagging stats)
  void set_in(const string&in){tmva_in_dir=check_isdir(in);}
  void set_out(const string&out){tmva_out_dir=check_mkdir(out);}

  bool debug;
  //tmva in and out
  string tmva_in_dir,tmva_out_dir,identifier;
  //event selection flags
  int ntag,njet,beff_cut,charm_ratio;
  bool cts_mv2;
  int loptv;
  bool use_btagvs,widemll,metcut;//use_btagvs is whether we use b-tagging variables--not possible for certain iterations of truth tagging

  //configuration parameters
  map<int,map<int,float> > btag_ratio_cuts;int def_cratio,def_beff_cut;
};

#endif
