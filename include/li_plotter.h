#ifndef LI_PLOTTER_H
#define LI_PLOTTER_H

#include "li_ranking.h"
#include "bdt_testing.h"
//#include "TMVA/correlations.h"//"correlations.C"

class li_plotter: public li_ranking{
 public:
  li_plotter(const string&in,const string&out,int lptv=0,int nj=2,bool btvs=true);
  vector<TString> all_plots(const TString&pdir,int tDF=0);
  pair<TString,double> single_plots(int whic,int nj,bool lpt,const TString&pdir,int tDF=0,bool split=false);

 protected:
  bool splitptv;//were events split into pTV categories for training? (yes gives you in principle better separation but lower stats (esp hi pTV))
  string var_log(int tDF=0)const{return log(tag()+tdf_str(tDF));}

  int n_li_cases()const{return 6;}
  int n_std_min()const{return -1;}
};

#endif
