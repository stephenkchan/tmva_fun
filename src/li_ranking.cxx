#include "li_ranking.h"

li_ranking::li_ranking(const string&in,const string&out,bool lptv,int nj,bool btvs):bdt_ranker(in,out,vector<TString>(),vector<TString>(),2,nj,77,20,false,lptv,btvs){}

void li_ranking::all_ranks(bool overwrite,int tDF){
  for(auto nj:jets()){
    for(auto lpt:ptvs()){
      for(auto whic:cases())single_rank(whic,nj,lpt,overwrite,tDF);
    }
  }
}

void li_ranking::single_rank(int whic,int nj,bool lpt,bool overwrite,int tDF){
  set_pars(whic,nj,lpt);
  cout<<"Begin training for type: "<<tag()<<endl;
  do_ranking(tag(),overwrite,tDF);
}
void li_ranking::finish_rank(int whic,int nj,bool lpt,bool overwrite,int tDF){
  set_pars(whic,nj,lpt);
  cout<<"Begin training for type: "<<tag()<<endl;
  finish_ranking(tag(),overwrite,tDF);
}

pair<vector<TString>,vector<TString> > li_ranking::vars(int whic)const{
  //training variables
  vector<TString> base=(whic<0?vector<TString>({"mBB","dRBB"}):vector<TString>({"j0_j1","angle_j0_j1"})),unranked;
  if(use_btagvs){
    //the cuts are lowest (highest) efficiencies so that we don't feed TMVA distributions of just one value
    //for pseudocontinuous tagging for tagged (untagged) jets; constant variables cause TMVA to crash
    TString bv=btag_var();
    if(beff_cut>60)unranked={bv+"B1",bv+"B2"};
    if(njet>2&&beff_cut<85)unranked.push_back(bv+"J3");
  }
  if(whic==777)return{{"mBB","dRBB"},unranked};
  string date=tmva_in_dir.substr(1+tmva_in_dir.find("/201"),tmva_in_dir.find("/",tmva_in_dir.find("/201")+1)-tmva_in_dir.find("/201"));
  if(date<"20160308")base[1]="dRBB";//angle_j0_j1 not implemented before then
  //"cut-based" default
  if(whic<n_std_min()||whic>=n_li_cases())return pair<vector<TString>,vector<TString> >(base,vector<TString>());
  vector<TString> train_std={"dEtaBB","pTV","dPhiVBB","dEtaVBB","mLL","MET"};
  if(use_btagvs||always_tag){
    train_std.push_back("pTB1");train_std.push_back("pTB2");
    if(njet>2)unranked.push_back("pTJ3");
  }
  vector<TString> train_li={"j0_l0", "j0_l1", "j1_l0", "j1_l1", "l0_l1", "angle_bbz_bbll","gamma_ZHz"};
  vector<TString> li_ll={"inv_ll", "angle_inv_ll"},li_bb={"inv_bb", "angle_inv_bb"},li_lb0={"inv_lb0","angle_inv_lb0"};
  vector<TString> std_met={"MET"};

  //add li base variables for whic>=0, else the standard variables
  if(whic>=0) unranked.insert(unranked.end(),train_li.begin(),train_li.end());
  else unranked.insert(unranked.end(),train_std.begin(),train_std.end());

  //MET candidates, if any; the max case takes all the invisibles
  if(whic==1) unranked.insert(unranked.end(),std_met.begin(),std_met.end());
  if(whic==2||whic==n_li_cases()-1) unranked.insert(unranked.end(),li_ll.begin(),li_ll.end());//ll inv
  if(whic==3||whic==n_li_cases()-1) unranked.insert(unranked.end(),li_bb.begin(),li_bb.end());//bb inv
  if(whic==4||whic==n_li_cases()-1) unranked.insert(unranked.end(),li_lb0.begin(),li_lb0.end());//lb0 inv

  return pair<vector<TString>,vector<TString> >(base,unranked);
}

string li_ranking::var_tag(int whic)const{
  string varset="cut";
  if(whic==0)varset="li-only";
  else if(whic==1)varset="li-met";
  else if(whic==2)varset="li-ll";
  else if(whic==3)varset="li-bb";
  else if(whic==4)varset="li-lb0";
  else if(whic==777)varset="test";
  else if(whic==n_li_cases()-1)varset="li-all";
  else if(whic>=n_std_min()&&whic<0)varset="std";
  varset+=(use_btagvs?"-pct":"-nopct");
  return varset;
}


