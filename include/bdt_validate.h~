#ifndef BDT_VALIDATE_H
#define BDT_VALIDATE_H

#include "bdt_trainer.h"

class bdt_validate: public bdt_trainer{
 public:
  bdt_validate(const string&in,const string&out,const string&tag,const vector<TString>&vs,bool overwrite=false,int nt=2,int nj=2,int eff_cut=70,int ratio=20,bool cts=true,bool lptv=true,bool wmll=true,bool met=false);
  bdt_validate(const vector<TString>&vs,bool overwrite,const bdt_base&base);
  bdt_validate(const bdt_trainer&train);
  pair<double,double> val_sig(int bg=-1,int trans_DF=0){return cumulat_sig(get_validation_bdt_dists(bg,trans_DF));}//optimal_sig(get_validation_bdt_dists(bg,trans_DF));}
  string get_in_dir()const{return tmva_in_dir;}
  string get_out_dir()const{return tmva_out_dir;}
  pair<TH1D*,TH1D*> get_bdt_dists()const{return start;}
  pair<TH1D*,TH1D*> get_validation_bdt_dists(int bg=-1,int trans_DF=0);
  void print_validation_bdts(TCanvas *c,const TString&pdir,int bg=-1,int trans_DF=0,const string&tag=""){print_bdt_plots(c,get_validation_bdt_dists(bg,trans_DF),pdir,"val-"+(tag==""?identifier:tag)+sam_tag(bg)+(trans_DF>0&&identifier.find(tdf_str(trans_DF).c_str())==string::npos?tdf_str(trans_DF):""),sam_label(bg),true);}
  
};
#endif
