#include "btag_ranking.h"

btag_ranking::btag_ranking(const string&top,int cut,int ratio,bool cts,int nj,bool lptv):bdt_ranker(btag_dir(top,cut,false),btag_dir(top,cut,true),{"mBB","dRBB"},std_unranked_btag(ratio,cts,nj,cut),2,nj,cut,ratio,cts,lptv),tmva_top_dir(check_isdir(top)){}

string btag_ranking::btag_dir(const string& dir,int cut,bool isout)const{
  ostringstream dirs;dirs<<check_dir(dir)<<"eff_"<<cut<<"/";
  if(isout)dirs<<"output/";
  return dirs.str();
}
void btag_ranking::all_ranks(bool overwrite,int tDF){
  for(auto cut:cuts()){
    for(auto ratio:ratios()){
      for(auto ct:cts())single_rank(cut,ratio,ct,overwrite,tDF);
    }
  }
}

void btag_ranking::single_rank(int cut,int ratio,bool cts,bool overwrite,int tDF){
  set_btag(cut,ratio,cts);
  cout<<"Begin btag ranks for: "<<btag_label()<<endl;
  do_ranking(string(btag_str()+region_str()),overwrite,tDF);
}
void btag_ranking::single_rank(bool overwrite,int tDF){
  cout<<"Begin btag ranks for: "<<btag_label()<<endl;
  do_ranking(string(btag_str()+region_str()),overwrite,tDF);
}

void btag_ranking::finish_rank(bool overwrite,int tDF){
  cout<<"Begin btag ranks for: "<<btag_label()<<endl;
  finish_ranking(string(btag_str()+region_str()),overwrite,tDF);
}

void btag_ranking::set_btag(int cut,int ratio,bool cts){
  set_cut(cut);set_ratio(ratio);set_cts(cts);
  set_in(btag_dir(cut,false));set_out(btag_dir(cut,true));
  unranked=std_unranked_btag(ratio,cts,njet,cut);
}

vector<TString> btag_ranking::std_unranked_btag(int ratio,bool cts,int nj,int cut)const{
  vector<TString>rankme={"dEtaBB","pTV","pTB1","pTB2","dPhiVBB","dEtaVBB","mLL","MET"};
  //the cuts are lowest (highest) efficiencies so that we don't feed TMVA distributions of just one value
  //for pseudocontinuous tagging for tagged (untagged) jets; constant variables cause TMVA to crash
  if(nj>=3){  
    rankme.push_back("mBBJ");rankme.push_back("pTJ3");
    if(cut<85)rankme.push_back(btag_var(ratio,cts)+"J3");
  }
  if(cut>60){rankme.push_back(btag_var(ratio,cts)+"B1");rankme.push_back(btag_var(ratio,cts)+"B2");}
  return rankme;
}
