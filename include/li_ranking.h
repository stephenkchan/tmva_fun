#ifndef LI_RANKING_H
#define LI_RANKING_H
#include "bdt_ranker.h"

class li_ranking : public bdt_ranker{
 public:
  li_ranking(const string&in,const string&out,bool lptv=false,int nj=2,bool btv=true);
  void all_ranks(bool overwrite=false,int tDF=0);
  void single_rank(int whic,int nj,bool lpt,bool overwrite=false,int tDF=0);
  void finish_rank(int whic,int nj,bool lpt,bool overwrite=false,int tDF=0);

 protected:
  void set_vars(int whic){ranked=vars(whic).first;unranked=vars(whic).second;}
  pair<vector<TString>,vector<TString>> vars(int whic)const;
  string var_tag(int whic)const;
  string tag()const{return var_tag(which)+string(region_str());}
  void set_pars(int whic,int nj,bool lpt){set_which(whic);set_njet(nj);set_ptv(lpt);set_vars(whic);}
  void set_which(int whic){which=whic;}

  int n_li_cases()const{return 6;}
  int n_std_min()const{return -1;}
  vector<int>jets()const{return{2,3};}
  vector<int>ptvs()const{return{1,0};}//it's technically a bool; true is lo pTV
  vector<int>cases()const{return{-2,-1,0,1,2,3,4};}
  int which;
};

#endif
