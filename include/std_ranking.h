#ifndef STD_RANKING_H
#define STD_RANKING_H
#include "bdt_ranker.h"

class std_ranking : public bdt_ranker{
 public:
  std_ranking(const string&in,const string&out,int nt=2,int nj=2,bool ptv=true);
  void all_ranks(bool overwrite=false,int tDF=0);
  void single_rank(bool overwrite=false,int tDF=0);
  void finish_rank(bool overwrite=false,int tDF=0);
  void single_rank(int nt,int nj,bool ptv,bool overwrite,int tDF);

 protected:
  vector<TString> std_unranked(int nj)const;
  vector<int>tags()const{return {2};}
  vector<int>jets()const{return {2,3};}
  vector<int>ptvs()const{return {1,0};}
};

#endif
