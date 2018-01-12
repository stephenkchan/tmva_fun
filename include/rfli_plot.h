#ifndef RFLI_PLOT_H
#define RFLI_PLOT_H

#include "rfli_rank.h"
#include "bdt_testing.h"
//#include "TMVA/correlations.h"//"correlations.C"

class rfli_plot: public rfli_rank{
 public:
  rfli_plot(const string&in,const string&out,int ptv=0,int nj=2,bool btvs=false,int whic=1);
  vector<TString> all_plots(const TString&pdir,int tDF=0);
  pair<double,double> single_plots_log_check(int whic,int nj,int ptv,const TString&pdir,int tDF=0,bool split=true);
  pair<double,double> single_plots(int whic,int nj,int ptv,const TString&pdir,int tDF=0,bool split=true);

 protected:
  bool splitptv;//were events split into pTV categories for training? (yes gives you in principle better separation but lower stats (esp hi pTV))
  void quad_sum_error(vector<double>&sigs,vector<double>&errors)const;
  string var_log(int tDF=0)const{return log(tag()+tdf_str(tDF));}

  int n_li_cases()const{return 6;}
  int n_std_min()const{return -1;}
};

#endif
