#include "rfli_rank.h"

rfli_rank::rfli_rank(const string&in,const string&out,int ptv,int nj,bool btvs,int whic)
  :bdt_ranker(in,out,vector<TString>(),vector<TString>(),2,nj,70,10,false,ptv,btvs),which(whic){set_vars(whic);}

void rfli_rank::all_ranks(bool overwrite,int tDF){
  for(auto nj:jets()){
    for(auto lpt:ptvs()){
      for(auto whic:cases())single_rank(whic,nj,lpt,overwrite,tDF);
    }
  }
}

void rfli_rank::single_rank(int whic,int nj,int ptv,bool overwrite,int tDF){
  set_pars(whic,nj,ptv);
  cout<<"Begin training for type: "<<tag()<<endl;
  do_ranking(tag(),overwrite,tDF);
}
void rfli_rank::finish_rank(int whic,int nj,int lpt,bool overwrite,int tDF){
  set_pars(whic,nj,lpt);
  cout<<"Begin training for type: "<<tag()<<endl;
  finish_ranking(tag(),overwrite,tDF);
}

pair<vector<TString>,vector<TString> > rfli_rank::vars(int whic)const{
  //training variables
  //let the default be "cut-based"
  vector<TString> base={"mBB","dRBB"},unranked;
  if(whic<=0)return{base,unranked};//cut-based
  else if(whic==2||whic==3)base=vector<TString>{"j0_j1"};
  else if(whic<9&&whic>3)base=vector<TString>{"MH"};
  else if(whic==99)return pair<vector<TString>,vector<TString> >(vector<TString>{"MH"},vector<TString>{"MCM","MZ"});
  else if(whic!=1){
    cerr<<"You have specified a bad which "<<which<<endl;
    exit(9);
  }if(use_btagvs){
    //the cuts are lowest (highest) efficiencies so that we don't feed TMVA distributions of just one value
    //for pseudocontinuous tagging for tagged (untagged) jets; constant variables cause TMVA to crash
    TString bv=btag_var();
    if(beff_cut>60)unranked={bv+"B1",bv+"B2"};
    if(njet>2&&beff_cut<85)unranked.push_back(bv+"J3");
  }
  vector<TString> std={"pTV","dPhiVBB","dEtaVBB","mLL","MET"};
//   if(use_btagvs||always_tag){
  std.push_back("pTB1");std.push_back("pTB2");
  if(njet>2){
    std.push_back("pTJ3");
    std.push_back("mBBJ");
  }
  vector<TString> li={"j0_l0", "j0_l1", "j1_l0", "j1_l1", "l0_l1", "angle_bbz_bbll","angle_bb_z","gamma_ZHz"},li_met={"MET"},rf={"MCM","MZ","cosCM","cosZ","cosH","Rpt","Rpz"},rf_met={"Rmet","dphiCMMet"},rfdp={"dphiLABCM","dphiCMZ","dphiCMH"},rf_select={"Rmet","dphiCMH"};
  
  if(whic==1)unranked.insert(unranked.end(),std.begin(),std.end());
  else if(whic<=3){//LI's
    unranked.insert(unranked.end(),li.begin(),li.end());
    if(whic==3)unranked.insert(unranked.end(),li_met.begin(),li_met.end());
  }
  else if(whic<9){//RF's
    unranked.insert(unranked.end(),rf.begin(),rf.end());
    if(whic==5||whic==7)unranked.insert(unranked.end(),rf_met.begin(),rf_met.end());//include rf met
    if(whic==6||whic==7)unranked.insert(unranked.end(),rfdp.begin(),rfdp.end());//include rf dphi's
    if(which==8)unranked.insert(unranked.end(),rf_select.begin(),rf_select.end());//include rf special
  }
  return pair<vector<TString>,vector<TString> >(base,unranked);
}

string rfli_rank::var_tag(int whic)const{
  string varset="cut";
  if(whic==1)varset="std";
  else if(whic==2)varset="li";
  else if(whic==3)varset="li-met";
  else if(whic==4)varset="rf";
  else if(whic==5)varset="rf-met";
  else if(whic==6)varset="rfdp";
  else if(whic==7)varset="rfdp-met";
  else if(whic==8)varset="rf-sel";
  else if(whic==99)varset="test";
  varset+=(use_btagvs?"-pct":"");
  return varset;
}


