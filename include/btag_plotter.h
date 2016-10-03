#ifndef BTAG_PLOTTER_H
#define BTAG_PLOTTER_H

#include "btag_ranking.h"
#include "bdt_testing.h"
//#include "TMVA/correlations.h"//"correlations.C"

class btag_plotter: public btag_ranking{
 public:
  btag_plotter(const string&top);
/*   ~btag_plotter(); */
  vector<TString> all_plots(const TString&pdir,int tDF=0);
  pair<TString,double> single_plots(int cut,int ratio,bool cts,int nj,bool lptv,const TString&pdir,int tDF);


 protected:
  TCanvas*c;
};

#endif
