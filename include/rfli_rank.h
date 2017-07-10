#ifndef RFLI_RANK_H
#define RFLI_RANK_H
#include "bdt_ranker.h"

//fuck backwards compatability and all that invisible shit; 2017-04-12 onwards only
class rfli_rank : public bdt_ranker{
 public:
  rfli_rank(const string&in,const string&out,int lptv=0,int nj=2,bool btv=true,int whic=1);
  void all_ranks(bool overwrite=false,int tDF=0);
  void single_rank(int whic,int nj,int lpt,bool overwrite=false,int tDF=0);
  void finish_rank(int whic,int nj,int lpt,bool overwrite=false,int tDF=0);
  void set_tftag(int tF){set_tag(tag()+tdf_str(tF));}
  string tag()const{return var_tag(which)+string(region_str(true));}

 protected:
  void set_vars(int whic){ranked=vars(whic).first;unranked=vars(whic).second;}
  pair<vector<TString>,vector<TString>> vars(int whic)const;
  string var_tag(int whic)const;
  void set_pars(int whic,int nj,int lpt){set_which(whic);set_njet(nj);set_ptv(lpt);}
  void set_which(int whic){which=whic;set_vars(whic);}

  vector<int>jets()const{return{2,3};}
  vector<int>ptvs()const{return{1,2};}//here, it's "reasonable" bins; lo ptv isn't useful evidently
  vector<int>cases()const{return{0,1,3,8};}//cut,std,li-met,rf-sel
  int which;
};

#endif
