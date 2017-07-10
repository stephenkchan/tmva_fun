#ifndef BDT_TESTING_H
#define BDT_TESTING_H

#include "bdt_trainer.h"
#include "bdt_validate.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

class bdt_testing: public bdt_trainer{
 public:
  bdt_testing(const string&in,const string&out,const string&tag,const vector<TString>&vs,bool overwrite=false,int nt=2,int nj=2,int eff_cut=70,int ratio=20,bool cts=true,int lptv=2,bool btv=true,bool wmll=true,bool met=false);
  bdt_testing(const vector<TString>&vs,bool overwrite,const bdt_base&base);
  bdt_testing(const bdt_trainer&train);
  pair<double,double> test_sig(int bg=-1,int trans_DF=0){return optimal_sig(make_testing_bdt_dists(bg,trans_DF));}
  string get_in_dir()const{return tmva_in_dir;}
  string get_out_dir()const{return tmva_out_dir;}
  pair<TH1D*,TH1D*> get_bdt_dists()const{return start;}
  pair<TH1D*,TH1D*> make_testing_bdt_dists(int bg=-1,int trans_DF=0,bool split=false)const;
  pair<double,double> print_test_bdts(TCanvas *c,const TString&pdir,int bg=-1,int trans_F=0,bool split=true){return print_bdt_plots(c,make_testing_bdt_dists(bg,trans_F,split),pdir,"test-"+identifier+sam_tag(bg),sam_label(bg)+" "+region_label(),false);}//bdt_validate((bdt_trainer)*this).val_sig(bg,trans_F).first);}
  
 protected:
  bool passes_cut(const map<TString,pair<float,TBranch*> >& vars,int ev,int nj)const;
  bool passes_cut(TTree*t,int ev,const TString& cut)const;//weird event numbers, so deprecated--promote when you update MC releases
};
#endif
