#ifndef BTAG_RANKING_H
#define BTAG_RANKING_H
#include "bdt_ranker.h"

class btag_ranking : public bdt_ranker{
 public:
  btag_ranking(const string&top,int cut=70,int ratio=20,bool cts=true,int nj=2,bool lptv=true);
  void all_ranks(bool overwrite=false,int tDF=0);
  void single_rank(int cut,int ratio,bool cts,bool overwrite=false,int tDF=0);
  void single_rank(bool overwrite=false,int tDF=0);
  void finish_rank(bool overwrite=false,int tDF=0);

 protected:
  void set_btag(int cut,int ratio,bool cts);
  string btag_dir(int cut,bool isout)const{return btag_dir(tmva_top_dir,cut,isout);}
  string btag_dir(const string& dir,int cut,bool isout)const;
  vector<TString> std_unranked_btag(int ratio,bool cts,int nj,int cut)const;
  vector<int>cuts()const{return {60,70,77,85};}
  vector<int>ratios()const{return {20};}//{0,10,20};}
  vector<int>cts()const{return {0};}//,1};}
  vector<int>jets()const{return {2,3};}
  vector<int>ptvs()const{return {1,0};}
  string tmva_top_dir;
};

#endif
