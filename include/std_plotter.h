#ifndef STD_PLOTTER_H
#define STD_PLOTTER_H

#include "std_ranking.h"
#include "bdt_testing.h"
//#include "TMVA/correlations.h"//"correlations.C"

class std_plotter: public std_ranking{
 public:
  std_plotter(const string&in,const string&out);
/*   ~std_plotter(); */
  vector<TString> all_plots(const TString&pdir,int tDF=0);
  pair<TString,double> single_plots(int nt,int nj,bool ptv,const TString&pdir,int tDF);


 protected:
  TCanvas*c;
};

#endif
