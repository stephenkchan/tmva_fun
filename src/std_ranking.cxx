#include "std_ranking.h"

std_ranking::std_ranking(const string&in,const string&out,int nt,int nj,bool ptv):bdt_ranker(in,out,{"mBB","dRBB"},std_unranked(nj),nt,nj,70,20,true,ptv){}

void std_ranking::all_ranks(bool overwrite,int tDF){
  for(auto nt:tags()){
    for(auto nj:jets()){
      for(auto ptv:ptvs())single_rank(nt,nj,ptv,overwrite,tDF);
    }
  }
}

void std_ranking::single_rank(bool overwrite,int tDF){
  cout<<"Begin std ranking for cuts region: "<<region_label()<<endl;
  do_ranking(string(region_str()),overwrite,tDF);
}

void std_ranking::finish_rank(bool overwrite,int tDF){
  cout<<"Begin std ranking for cuts region: "<<region_label()<<endl;
  finish_ranking(string(region_str()),overwrite,tDF);
}

void std_ranking::single_rank(int nt,int nj,bool ptv,bool overwrite,int tDF){
  set_njet(nj);set_ntag(nt);set_ptv(ptv);
  cout<<nt<<" "<<nj<<" "<<ptv<<"; "<<ntag<<" "<<njet<<" "<<loptv<<": Begin std ranking for cuts region: "<<region_label()<<endl;
  do_ranking(string(region_str()),overwrite,tDF);
}

vector<TString> std_ranking::std_unranked(int nj)const{
  vector<TString>rankme={"dEtaBB","pTV","pTB1","pTB2","dPhiVBB","dEtaVBB","mLL","MET","MV2c20B1","MV2c20B2"};
  if(nj>=3){ rankme.push_back("mBBJ");rankme.push_back("pTJ3");}
  return rankme;
}
