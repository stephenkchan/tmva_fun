#ifndef BDT_TRAINER_H
#define BDT_TRAINER_H

#include "TTree.h"
#include "TMVA/Factory.h"
#include "TMVA/MethodBDT.h"
#include "TMVA/Tools.h"
#include "bdt_base.h"
using namespace std;
using namespace TMVA;

class bdt_trainer : public bdt_base{
 public:
  bdt_trainer(const string&in,const string&out,const string&tag,const vector<TString>&vs,bool overwrite,int nt=2,int nj=2,int eff_cut=70,int ratio=20,bool cts=true,int lptv=1,bool btv=true,bool wmll=true,bool met=false,bool do_validation=true,int seed=0,bool _isdb=false);
  bdt_trainer(const vector<TString>&vs,bool overwrite,const bdt_base&base,bool do_validation=true,int seed=0,bool _isdb=false);
  ~bdt_trainer();
  void train(int bg=-99)const{train(bg,validate,kseed);}
  pair<double,double> significance();
  pair<double,double> print_bdt_plots(TCanvas *c,pair<TH1D*,TH1D*>iris,const TString&pdir,const TString&tag,const TString&bglabel,bool funky_style=false,double cut=-99999);
  pair<double,double> print_bdt_plots(TCanvas *c,const TString&pdir,int bg=-1){return print_bdt_plots(c,start,pdir,"train-"+identifier+sam_tag(bg),sam_label(bg)+" "+region_label(),true);}
  void bdt_dists(bool overwrite,bool do_validation,int seed,int bg=-99);
  pair<double,double> evts_S_B(int bg)const;

 protected:

  //
  bool is_signal(int sam)const;
  void train(int bg,bool do_validation,int seed)const;
  void silence_tree(TTree*tree)const;
  TString root_file(int bg,bool do_validation,int seed)const{return (is_standard(do_validation,seed)?filename(vars,bg):production_filename(bg));}
  TString root_file(int bg=-1)const{return root_file(bg,validate,kseed);}
  void initialize(bool overwrite,bool do_validation,int seed,bool db);
  bool is_standard(bool do_validation,int seed)const{return do_validation && seed==0;}
  TString xml_tag(bool do_validation,int seed,bool split=true)const{return (is_standard(do_validation,seed)?"TMVAClassification"+ftag(vars,split):production_tag());}
  TString xml_file(bool do_validation,int seed,bool split=true)const{return xml_tag(do_validation,seed,split)+"_BDT.weights.xml";}
  TString xml_file(bool split)const{return xml_file(validate,kseed,split);}
  TString production_tag(int bg=-1)const{ return Form("Harvard_13TeV_llbb_%s_%iof2",identifier.c_str(),kseed%2==1);}
  TString production_filename(int bg=-1)const{return Form("%s%s.root",tmva_out_dir.c_str(),production_tag(bg).Data());}

  //file names should look like:
  //TMVAClassification_BDT_LiverpoolBmham_8TeV_llbb_2tag2jet_ptv0_120_ZllH125_0of2_root5.34.05_v3.0.xml;

  //variables to train, ordered
  vector<TString> vars;

  //production type: if validation is shut off, then events are split into two groups, not three
  //this means no holdout method, just TMVA's testing and training---we do this to harmonize
  //and boost statistics for final results; the kseed is the seed for training, hopefully k-folds into
  //even/odd
  bool validate;int kseed;bool isdb;

  //configuration parameters
  vector<TFile*>opfiles;

  //the bdt training distributions
  pair<TH1D*,TH1D*> start;
};

#endif
